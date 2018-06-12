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
        ! $OMP parallel do default(shared) private(xind, yind, x, y, tempvar)
        do xind = 1, data_shape(1)
            do yind = 1, data_shape(2)
                tempvar = data(xind, yind)
                if(tempvar .eq. 0) then
                    cycle
                end if
                if(custom_abs(tempvar) .ge. smallest_val) then

                    x = (xind - 1) * d + xmin
                    y = (yind - 1) * d + ymin

                    ! $OMP critical(build_buffer)
                    acc = acc + 1
                    if (acc .le. max_points) then 
                        buffer(acc) = DataPoint(x, y, tempvar)
                    end if 
                    ! $OMP end critical(build_buffer)

                end if
            end do
        end do
        ! $OMP end parallel do

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

        ! $OMP parallel do default(shared) private(xind, yind, x, y, tempvar)
        do xind=1, data_shape(1)

            ! We can break when we know we have all the points
            if (acc .eq. 0) then
                exit 
            end if

            do yind = 1, data_shape(2)

                tempvar = data(xind, yind)
                if(tempvar .eq. 0) then
                    cycle
                end if

                if(custom_abs(tempvar) .ge. smallest_val) then 
                    x = (xind - 1) * d + xmin
                    y = (yind - 1) * d + ymin

                    ! $OMP critical 
                    buffer(acc) = DataPoint(x, y, tempvar)
                    acc = acc - 1;
                    ! $OMP end critical
                end if

                ! The inner break
                if (acc .eq. 0) then 
                    exit 
                end if 
            end do
        end do
        ! $OMP end parallel do

        ! Sort the buffer to get the right number of points
        call partial_sort(buffer, max_points)

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

end module