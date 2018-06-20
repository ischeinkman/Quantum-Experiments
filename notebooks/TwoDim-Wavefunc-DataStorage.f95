! Module for storing data from a 2D wavefunction simulation 
! The data is filtered so only the max_points highest absolute values 
! are stored; in addition, only non-zero points are stored.

module Wavefunc_DataStorage 

    use iso_fortran_env, only : REAL64
    use omp_lib
    implicit none

    public :: DataStore, TimeNode, DataPoint, store_data, custom_abs, verify_storage

    ! A single data point, in a time slice, with 
    ! spacial coordinates and the complex output
    type DataPoint 
        Real(kind=Real64) :: x, y
        Complex(kind=Real64) :: value
    end type DataPoint

    ! A full time slice
    type TimeNode
        Real(kind=Real64) :: t
        type(DataPoint), dimension(:), allocatable :: data
        type(TimeNode), pointer :: next
    end type TimeNode

    ! The linked list datastore
    type DataStore
        type(TimeNode), pointer :: head, tail
        Integer size
    end type DataStore

    contains 

    elemental pure function DataPoint_Check(dp1, dp2) result(are_eq) 
        type(DataPoint), intent(in) :: dp1, dp2
        LOGICAL :: are_eq
        are_eq = dp1%x == dp2%x .and. dp1%y == dp2%y .and. dp1%value == dp2%value
    end function

    subroutine verify_storage(storage, retval)
        type(DataStore), intent(in) :: storage
        LOGICAL, intent(out) ::retval

        Integer :: acc
        type(TimeNode), pointer :: data_idx
        Integer :: node_shape(1), pt_idx

        Print *, "START THE VERIFICATION!"
        Print *, "......................."
        Print *, "Checking size. "
        if(storage%size .eq. 0) then 
            retval = .TRUE.
            return
        end if
        Print *, "Storage is claiming size of: ", storage%size
        if(.not. ASSOCIATED(storage%head)) then
            retval = .FALSE.
            return
        end if
        Print *, "Storage is claiming a head."
        Print *, "Head is: ", storage%head%t
        data_idx => storage%head

        retval = .TRUE.
        Print *, "Starting verification loop..."
        acc = 0
        do while(ASSOCIATED(data_idx))
            
            node_shape = shape(data_idx%data)
            Print *, "  Node shape: ", node_shape

            acc = acc + node_shape(1)

            do pt_idx = 1, node_shape(1)
                retval = retval .and. data_idx%data(pt_idx)%value .ne. 0
                if(.not. retval) then
                    Print *, "Got value of ", data_idx%data(pt_idx), "; now failing."
                    return 
                end if
            end do 
            data_idx => data_idx%next
        end do 

        retval = (acc .eq. storage%size)

        if(.not. retval) then 
            Print *, "Got acc of ", acc, " vs size of ", storage%size
        end if 

    end subroutine

    ! Private method to append a single node to the end of the linked list
    subroutine append_node(storage, node)
        type(DataStore), intent(inout) :: storage
        type(TimeNode), target, intent(in) :: node

        Integer, dimension(1) :: node_shape

        node_shape = shape(node%data)

        ! If the list is empty 
        if (.not. ASSOCIATED(storage%head) .or. storage%size .eq. 0) then
            allocate(storage%head)
            storage%head = node
            storage%tail => storage%head
            storage%size = node_shape(1)
            return
        end if

        ! Otherwise
        allocate(storage%tail%next)
        storage%tail%next = node
        storage%tail => storage%tail%next
        NULLIFY(storage%tail%next)
        storage%size = storage%size + node_shape(1)
    end subroutine

    ! Store data from a full, uncompressed universe time slice into 
    ! a data store
    subroutine store_data(storage, data, t, d, xmin, ymin, max_points, smallest_val)

        ! Parameters:

        type (DataStore), intent(inout) :: storage ! The data store
        Complex(kind=Real64), dimension (:, :), intent(in) :: data ! The universe to store
        Real(kind=Real64), intent(in) :: t ! The current time value of the slice
        Real(kind=Real64), intent(in) ::  d, xmin, ymin ! Values to convert from an index to spacial coordinate
        Integer, intent(in) :: max_points ! The maximum number of points to store from the slice
        Real(kind=Real64), intent(in) :: smallest_val ! The smallest allowed abs 

        ! Local vars: 
        Integer, dimension(2) :: data_shape ! The dimensions of the universe
        Integer :: xind, yind ! Indices to iterate through the slice
        Integer :: acc ! The total number of nonzero data points
        type(DataPoint), dimension(:), allocatable :: buffer ! A buffer to hold all the datapoints
        Real(kind=Real64) :: x, y ! Temp vars to hold the converted spacial coord
        type(TimeNode), target :: retval ! The object to be added into the store
        Complex(kind=real64) :: tempvar

        ! Initialize the new timenode
        retval%t = t
        NULLIFY(retval%next)
        

        ! Get the total number of points we have, as well as the first max_points
        ! points just in case we can easily store then 

        data_shape = SHAPE(data)
        allocate(buffer(max_points))
        acc = 0
        !$OMP parallel do default(shared) private(xind, yind, x, y, tempvar)
        do xind = 1, data_shape(1)
            do yind = 1, data_shape(2)
                tempvar = data(xind, yind)
                if(tempvar .eq. 0) then
                    cycle
                end if
                if(custom_abs(tempvar) .ge. smallest_val) then

                    x = (xind - 1) * d + xmin
                    y = (yind - 1) * d + ymin

                    !$OMP critical(build_buffer)
                    acc = acc + 1
                    if (acc .le. max_points) then 
                        buffer(acc) = DataPoint(x, y, tempvar)
                    end if 
                    !$OMP end critical(build_buffer)

                end if
            end do
        end do
        !$OMP end parallel do

        ! If we are in the correct range, immediately store the new data
        if(acc .le. max_points) then 
            allocate(retval%data(acc))
            retval%data(1:acc) = buffer(1:acc)
            deallocate(buffer)
            call append_node(storage, retval)
            return 
        end if

        ! Otherwise, reallocate the buffer and collect ALL the points
        deallocate(buffer)
        allocate(buffer(acc))

        do xind=1, data_shape(1)

            ! We can break when we know we have all the points
            if (acc .eq. 0) then
                exit 
            end if

            !$OMP parallel do default(shared) private(yind, x, y, tempvar)
            do yind = 1, data_shape(2)

                if(acc .eq. 0) then 
                    cycle 
                end if

                tempvar = data(xind, yind)
                if(tempvar .eq. 0) then
                    cycle
                end if

                if(custom_abs(tempvar) .ge. smallest_val) then 
                    x = (xind - 1) * d + xmin
                    y = (yind - 1) * d + ymin

                    !$OMP critical 
                    buffer(acc) = DataPoint(x, y, tempvar)
                    acc = acc - 1;
                    !$OMP end critical
                end if

                ! The inner break
                if (acc .eq. 0) then 
                    cycle 
                end if 
            end do
            !$OMP end parallel do
        end do

        ! Sort the buffer to get the right number of points
        call parallel_partial_sort(buffer, max_points)

        allocate(retval%data(max_points))
        retval%data(1:max_points) = buffer(1:max_points)
        deallocate(buffer)
        call append_node(storage, retval)

    end subroutine

    ! Quickselect algorithm 
    ! Sorts data such that the items between 1 and max_points are smaller
    ! than those from max_points to the end of the array
    ! Source: https://en.wikipedia.org/wiki/Quickselect
    subroutine partial_sort(data, max_points)

        type(DataPoint), dimension (:), intent(inout) :: data
        Integer, intent(in) :: max_points

        Integer :: left, right, pivot, partition_idx
        Integer, dimension(1) :: data_shape
        type(DataPoint) :: pivotValue, tmp

        left = 1
        ! Start at the end
        data_shape = shape(data)
        right = data_shape(1)

        if(right .le. max_points) then 
            return 
        end if

        do
            ! FUTURE: Get a better pivot algorithm
            pivot = (right + left) /2 

            ! Partition subroutine

            pivotValue = data(pivot)

            ! Push the pivot to the end so we don't sort it
            tmp = data(pivot) 
            data(pivot) = data(right)
            data(right) = tmp

            ! Go through the list, storing where the border is in pivot
            pivot = left
            do partition_idx = left, right - 1
                if (custom_abs(data(partition_idx)%value) .lt. custom_abs(pivotValue%value)) then
                    tmp = data(pivot)
                    data(pivot) = data(partition_idx)
                    data(partition_idx) = tmp
                    pivot = pivot + 1
                end if
            end do
            
            ! Put the current pivot value back
            tmp = data(pivot) 
            data(pivot) = data(right)
            data(right) = tmp
            
            ! Fin partion

            ! Figure out the next subarray to check
            if (pivot .eq. max_points) then
                return 
            else if(pivot .gt. max_points) then 
                right = pivot - 1
            else 
                left = pivot + 1
            end if
        end do

    end subroutine

    ! A custom complex absolute value function for dealing with small numbers
    elemental function custom_abs(inpval) result(retval)
        Complex(kind=REAL64), intent(in) :: inpval
        Real(kind=REAL64) :: retval 
        
        Real(kind=REAL64) :: realpart, imgpart
        Real(kind=REAL64) :: expn


        ! If standard abs works just use that
        retval = abs(inpval)
        if (retval .ne. 0 .or. inpval .eq. 0) then
            return
        end if

        realpart = REAL(inpval)
        imgpart = AIMAG(inpval)

        ! If one of the parts is 0 then we can use the standard, not square-rooting
        ! real number abs
        if (realpart .eq. 0) then
            retval = abs(imgpart)
            return 
        else if (imgpart .eq. 0) then
            retval = abs(realpart)
            return
        end if

        ! We fix the rounding by either removing the 10 factor if we can,
        ! or just multiplying and then dividing by a big number if we cant
        expn = 1D1 ** (-1 * LOG10(realpart))

        if (expn .ne. expn .or. expn .eq. 0) then 
            expn = 1D1**200
        end if

        retval = abs(inpval * expn) /expn

    end function    
    
    ! Quickselect algorithm 
    ! Sorts data such that the items between 1 and max_points are smaller
    ! than those from max_points to the end of the array
    ! Source: https://en.wikipedia.org/wiki/Quickselect
    subroutine parallel_partial_sort(arr, max_points)

        type(DataPoint), dimension (:), intent(inout) :: arr
        Integer, intent(in) :: max_points

        Integer, parameter :: MIN_PAR_SIZE = 100


        ! Shared variables
        ! ----------------
        type(DataPoint), dimension(:), allocatable :: start_buffer
        Integer :: data_idx
        Integer :: thread_count
        Integer :: buffer_idx
        Integer :: rept_count
        Integer :: left_to_find
        Integer, dimension(1) :: data_shape

        ! Semishared variables
        ! --------------------
        type(DataPoint), dimension(:, :), allocatable :: buffers
        Integer, dimension(:), allocatable :: buffer_sizes
        Real(kind=Real64), dimension(:), allocatable :: borders
        
        ! Private variables
        ! ------------------
        Integer :: idx
        Integer :: thread_num
        Real(Kind=Real64) :: left_border, right_border
        Real(Kind=Real64) :: tmpval

        data_shape = shape(arr)
        if(data_shape(1) .le. max_points) then 
            return 
        end if

        allocate(start_buffer(data_shape(1)))
        start_buffer(:) = arr(:)

        data_idx = 1

        !$OMP PARALLEL shared(thread_count)
        if (omp_get_thread_num() .eq. 0) then 
            thread_count = omp_get_num_threads()
        end if
        !$OMP END PARALLEL

        allocate(buffer_sizes(thread_count))
        allocate(borders(thread_count - 1))
        allocate(buffers(thread_count, data_shape(1)))
        rept_count = 0

        do while(data_idx .le. max_points)
            ! Set each border to a custom_abs(value) from the data
            call getPivotValues(start_buffer, data_shape(1), borders)
            rept_count = rept_count +1
            call sort(borders)
            do buffer_idx = 2, thread_count - 2
                if(any(borders(1:buffer_idx - 1) .eq. borders(buffer_idx)) .or. & 
                   any(borders(buffer_idx + 1 : thread_count - 1) .eq. borders(buffer_idx)) &
                ) then
                    Print *, "DUP BORDER: ", borders
                    Print *, arr(-1)
                end if
            end do

            ! Reset the buffers
            !allocate(buffers(thread_count, data_shape(1)))
            buffer_sizes(:) = 0


            ! Parallelly partition the buffer into the different regions
            !$OMP parallel default(shared) private(tmpval, thread_num, idx, left_border, right_border)
            thread_num = omp_get_thread_num() +1
                
            ! Set up this thread's borders
            if(thread_num .eq. 1) then 
                left_border = -1D99
                right_border = borders(1)
            else if(thread_num .eq. thread_count) then 
                left_border = borders(thread_count - 1)
                right_border = 1D99
            else 
                left_border = borders(thread_num - 1)
                right_border = borders(thread_num)
            end if

            ! Construct this thread's partition
            if(left_border .ne. right_border) then 
                do idx = 1 , data_shape(1)
                    tmpval = custom_abs(start_buffer(idx)%value)
                    if(tmpval .lt. right_border .and. tmpval .ge. left_border) then
                        buffer_sizes(thread_num) = buffer_sizes(thread_num) + 1
                        buffers(thread_num, buffer_sizes(thread_num)) = start_buffer(idx)
                    end if
                end do
            end if
            !$OMP end parallel

            do buffer_idx = 1, thread_count
                if(buffer_sizes(buffer_idx) .lt. 1) then 
                    cycle 
                end if 

                left_to_find = max_points - data_idx + 1
                if(buffer_sizes(buffer_idx) .eq. left_to_find) then 
                    arr(data_idx:max_points) = buffers(buffer_idx, 1:left_to_find)
                    data_idx = max_points
                    exit 
                else if(buffer_sizes(buffer_idx) .gt. left_to_find) then
                    start_buffer(1:buffer_sizes(buffer_idx)) = buffers(buffer_idx, 1:buffer_sizes(buffer_idx))
                    data_shape(1) = buffer_sizes(buffer_idx)
                    EXIT
                else 
                    arr(data_idx:buffer_sizes(buffer_idx) + data_idx - 1) = buffers(buffer_idx, 1:buffer_sizes(buffer_idx))
                    data_idx = data_idx + buffer_sizes(buffer_idx)
                end if
            end do

            if(data_idx .eq. max_points) then
                exit
            else if (data_idx .gt. max_points) then 
                Print *, "SOMEHOW OVERFILLED ON PAR QSELECT."
                Print *, arr(-1)
            else if(data_shape(1) .le. thread_count * 2 .or. data_shape(1) .le. MIN_PAR_SIZE) then 
                call partial_sort(start_buffer(1:data_shape(1)), max_points - data_idx + 1)
                arr(data_idx:max_points) = start_buffer(1:max_points - data_idx + 1)
                EXIT
            end if

            !deallocate(buffers)

        end do

        deallocate(start_buffer)
        deallocate(buffers)
        deallocate(buffer_sizes)
        deallocate(borders)
    end subroutine

    ! Gets distinct random pivot values if possible
    subroutine getPivotValues(arr, arr_size, pivots) 
        type(DataPoint), dimension (:), intent(in) :: arr
        Integer, intent(in) :: arr_size
        Real(kind=Real64), dimension (:), intent(out) :: pivots
        
        Integer, parameter :: MAX_RETRIES = 10

        Integer :: pivot_size, idx, retry_num

        pivot_size = size(pivots)

        if (pivot_size .ge. arr_size) then
            pivots(1:arr_size) = custom_abs(arr%value)
            pivots(arr_size+1 : pivot_size) = custom_abs(arr(arr_size)%value)
            return 
        end if
        
        call RANDOM_NUMBER(pivots)

        pivots = 1 + INT(pivots * (arr_size - 1))

        pivots = custom_abs(arr(INT(pivots))%value)
        
        do idx = 1, pivot_size - 1
            do retry_num = 0, MAX_RETRIES
                if (pivots(idx) .ne. 0 .and. & 
                    .not. any(pivots(1:idx-1) == pivots(idx)) .and. &
                    .not. any(pivots(idx+1 :) == pivots(idx)) &
                ) then 
                    EXIT
                end if
                call RANDOM_NUMBER(pivots(idx))
                pivots(idx) = 1 + INT(pivots(idx) * (arr_size - 1))
                pivots(idx) = custom_abs(arr(INT(pivots(idx)))%value)
            end do
        end do
    end subroutine

    ! Sorts small amounts of data. 
    subroutine sort(data)

        Real(kind=Real64), dimension (:), intent(inout) :: data

        Integer :: i, min_idx
        Real(kind=Real64) :: tmp
        
        do i =1, size(data) - 1
            min_idx = MINLOC(data(i:), 1) + i - 1
            if(data(i) .gt. data(min_idx)) then 
                tmp = data(i)
                data(i) = data(min_idx)
                data(min_idx) = tmp
            end if
        end do 
    end subroutine

    subroutine verify_parallel_psort(full_dt, max_points)
        type(DataPoint), dimension(:), intent(in) :: full_dt
        Integer, intent(in) :: max_points
        
        type(DataPoint), dimension(:), allocatable :: seq_buffer, par_buffer
        Integer, dimension(1) :: data_shape
        Integer :: seq_idx, par_idx


        data_shape = shape(full_dt)

        ! Get the parallel data

        allocate(par_buffer(data_shape(1)))
        par_buffer = full_dt
        call parallel_partial_sort(par_buffer, max_points)
        
        ! Get the seq data
        
        allocate(seq_buffer(data_shape(1)))
        seq_buffer = full_dt
        call partial_sort(seq_buffer, max_points)

        ! Verify

        if (max_points .le. 10) then
            Print *, "Got seq of:"
            Print *, seq_buffer
            Print *, "\nGot par of:"
            Print *, par_buffer
        end if
       
        if(minval(custom_abs(seq_buffer(1:max_points)%value)) .ne. minval(custom_abs(par_buffer(1:max_points)%value))) then
            Print *, "Inequal min. "
            Print *, seq_buffer(-2)
        else if(maxval(custom_abs(seq_buffer(1:max_points)%value)) .ne. maxval(custom_abs(par_buffer(1:max_points)%value))) then
            Print *, "Inequal max: ", & 
                     "Seq: ", maxval(custom_abs(seq_buffer(1:max_points)%value)), &
                     "Par: ", maxval(custom_abs(par_buffer(1:max_points)%value))
            Print *, "Full seq:"
            Print *, abs(seq_buffer(1:max_points)%value)
            Print *, "Full par:"
            Print *, abs(par_buffer(1:max_points)%value)
            Print *, seq_buffer(-2)
        end if

        do seq_idx = 1, max_points
            if(.not. any(DataPoint_Check(par_buffer, seq_buffer(seq_idx)))) then 
                Print *, "Par is missing: ", seq_idx, " => ", seq_buffer(seq_idx)
                !Print *, seq_buffer(-2)
            end if
        end do
        Print *, "Check for missing pars. "
        do par_idx = 1, max_points
            if(.not. any(DataPoint_Check(seq_buffer, par_buffer(par_idx)))) then 
                Print *, "Seq is missing: ", par_idx, " => ", par_buffer(par_idx)
                Print *, par_buffer(-2)
            end if
        end do
        Print *, "Checked for missing seqs."


        Print *, "Verification complete."

    end subroutine

end module