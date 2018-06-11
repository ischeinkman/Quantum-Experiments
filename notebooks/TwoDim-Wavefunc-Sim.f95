
module Wave_Packet_Simulation

    use iso_fortran_env, only : REAL64
    use omp_lib
    implicit none
    !private 
    !public::start_wave_packet, first_iteration_at_point, iterate_psi_at_point, &
    !        calc_psi_coeff, calc_psiOff_coeff, calc_psiV_coeff, normalize_matrix, &
    !        custom_abs, DataPoint

    type DataPoint
        REAL(kind=REAL64) :: t, x, y
        COMPLEX(kind=REAL64) :: value
    end type DataPoint

    contains

    ! Simple wave packet centered at (0,0) with momentum px, py and width sigma
    pure function start_wave_packet(x, y, px, py, sigma, hbar, pi) result (psi_at_point)
        REAL(kind=REAL64), intent(in) :: x, y, px, py, sigma, hbar, pi
        COMPLEX(kind=REAL64) :: psi_at_point
        
        REAL(kind=REAL64) :: coeff_val
        COMPLEX(kind=REAL64) :: exp_val

        !coeff_val = sqrt(sigma *  sqrt(pi))
        !exp_val = (px * x + py * y) * (0, 1.0)/hbar - (x * x + y * y)/(2 * sigma * sigma) 

        !psi_at_point = coeff_val * exp(exp_val)

        if (x .eq. 0 .and. y .eq. 0) then 
            psi_at_point = 1
        else 
            psi_at_point = 0
        end if

        
    end function start_wave_packet


    ! The first step iteration, since we don't have access to t-s at that point
    pure function first_iteration_at_point(xind, yind, s, d, hbar, m, curV, curPsi) result(psi_next)
        Integer, intent(in) :: xind, yind
        Real(kind=REAL64), intent(in) :: s, d, hbar, m, curV
        Complex(kind=REAL64), intent(in), dimension(:, :) :: curPsi

        Complex(kind=REAL64) :: psi_next

        Complex(kind=REAL64) :: offPsi, offPsi_coeff, psiV_coeff, psi_coeff, curPsiPoint
        Integer, dimension(2) :: realitySize


        psi_coeff = 1 - (0, -2) *hbar * s/(m * d * d)
        psiV_coeff = (0, -1) * s /hbar
        offPsi_coeff = (0, 0.5) * hbar * s/(m * d * d)
        
        curPsiPoint = curPsi(xind, yind)
        
        ! We default the value of Psi to 0 at the borders, since we are assuming
        ! a particle in a finite box.
        realitySize = shape(curPsi)
        if (xind .gt. 1 .and. xind .lt. realitySize(1)) then 
            offPsi = curPsi(xind + 1, yind) + curPsi(xind - 1, yind)
        else if(xind .lt. realitySize(1)) then 
            offPsi = curPsi(xind + 1, yind) + 0
        else 
            offPsi = curPsi(xind - 1, yind) + 0
        end if
        
        ! Same for y
        if (yind .gt. 1 .and. yind .lt. realitySize(2)) then 
            offPsi = offPsi + curPsi(xind, yind + 1) + curPsi(xind, yind - 1)
        else if(yind .lt. realitySize(2)) then 
            offPsi = offPsi + curPsi(xind, yind + 1) + 0
        else 
            offPsi = offPsi + curPsi(xind, yind - 1) + 0
        end if

        psi_next = curPsiPoint * psi_coeff + curPsiPoint * curV * psiV_coeff + offPsi * offPsi_coeff

    end function first_iteration_at_point

    ! Get the coefficient for the \Psi(x, y, t) term
    pure function calc_psi_coeff(m, hbar, s, d) result(retval)
        REAL(kind=REAL64), intent(in) :: m, hbar, s, d
        COMPLEX(kind=REAL64) :: retval

        retval = (-4.0, 1.0) * s * hbar / (m * d * d)
        
    end function calc_psi_coeff

    ! Get the coefficient for the \V(x, y, t) * \Psi(x, y, t) term
    pure function calc_psiV_coeff(m, hbar, s, d) result(retval)
        REAL(kind=REAL64), intent(in) :: m, hbar, s, d
        COMPLEX(kind=REAL64) :: retval

        retval = (-2.0, 1.0) * s/hbar
        
    end function calc_psiV_coeff

    ! Get the coefficient for the gradient components summation term
    pure function calc_psiOff_coeff(m, hbar, s, d) result(retval)
        REAL(kind=REAL64), intent(in) :: m, hbar, s, d
        COMPLEX(kind=REAL64) :: retval

        retval = hbar * s * (0.0, 1.0) / (m * d * d)
    end function calc_psiOff_coeff

    ! General function to iterate the wave at a point
    pure function iterate_psi_at_point(xind, yind, psi_coeff, psiV_coeff, psiOff_coeff, curV, psi_prev, psi_cur) result(psi_next)
        INTEGER, intent(in) :: xind, yind
        REAL(kind=REAL64), intent(in) :: curV
        COMPLEX(kind=REAL64), intent(in) :: psi_coeff, psiV_coeff, psiOff_coeff
        Complex(kind=REAL64), dimension(:, :), intent(in) :: psi_prev, psi_cur
        COMPLEX(kind=REAL64) :: psi_next

        COMPLEX(kind=REAL64) :: offPsi
        Integer, dimension(2) :: realitySize
        
        ! We default the value of Psi to 0 at the borders, since we are assuming
        ! a particle in a finite box.
        realitySize = shape(psi_cur)
        if (xind .gt. 1 .and. xind .lt. realitySize(1)) then 
            offPsi = psi_cur(xind + 1, yind) + psi_cur(xind - 1, yind)
        else if(xind .lt. realitySize(1)) then 
            offPsi = psi_cur(xind + 1, yind) + 0
        else 
            offPsi = psi_cur(xind - 1, yind) + 0
        end if
        
        ! Same for y
        if (yind .gt. 1 .and. yind .lt. realitySize(2)) then 
            offPsi = offPsi + psi_cur(xind, yind + 1) + psi_cur(xind, yind - 1)
        else if(yind .lt. realitySize(2)) then 
            offPsi = offPsi + psi_cur(xind, yind + 1) + 0
        else 
            offPsi = offPsi + psi_cur(xind, yind - 1) + 0
        end if

        psi_next = psi_prev(xind, yind) + (psi_coeff + curV * psiV_coeff) * psi_cur(xind, yind) + psiOff_coeff * offPsi

    end function iterate_psi_at_point

    ! Normalizes the wave function data, so sum(abs(matrix)**2) = 1.
    subroutine normalize_matrix(matrix)
        COMPLEX(kind=REAL64), intent(inout), dimension(:, :) :: matrix

        INTEGER, dimension(2) :: dims
        REAL(kind=REAL64) :: sum_squares
        INTEGER :: xind, yind

        dims = SHAPE(matrix)
        sum_squares = 0.0

        ! Get the normalization coefficient at that time
        ! $OMP parallel do default(shared) private(xind, yind) reduction (+:sum_squares)
        do xind = 1, dims(1)
            do concurrent (yind = 1:dims(2))
                if (matrix(xind, yind) .ne. 0) then 
                    sum_squares = sum_squares + 1D100 * custom_abs(matrix(xind, yind))**2
                end if 
            end do
        end do
        ! $OMP end parallel do 
        
        sum_squares = sum_squares /1D100

        ! Short Circuit on a bad universe
        if (sum_squares .eq. 0) then
            Print *, "Bad universe. " 
            return 
        else if(abs(sum_squares - 1) .lt. 0.0001) then
            Print *, "SPICY! Sumsq = ", sum_squares
        end if
        
        ! Get the normalization coefficient at that time

        ! $OMP parallel do default(shared) private(xind, yind) shared(sum_squares, matrix)
        do xind = 1, dims(1)
            do concurrent (yind = 1:dims(2))
                if (matrix(xind, yind) .ne. 0) then 
                    matrix(xind, yind) = matrix(xind, yind) / sqrt(sum_squares)
                end if 
            end do
        end do
        ! $OMP end parallel do

    end subroutine 

    ! A custom complex absolute value function for dealing with small numbers
    elemental function custom_abs(inpval) result(retval)
        Complex(kind=REAL64), intent(in) :: inpval
        Real(kind=REAL64) :: retval 
        
        Real(kind=REAL64) :: realpart, imgpart
        Real(kind=REAL64) :: expn


        retval = abs(inpval)
        if (retval .ne. 0 .or. inpval .eq. 0) then
            return
        end if

        realpart = REAL(inpval)
        imgpart = AIMAG(inpval)

        if (realpart .eq. 0) then
            retval = abs(imgpart)
            return 
        else if (imgpart .eq. 0) then
            retval = abs(realpart)
            return
        end if

        expn = 1D1 ** (-1 * LOG10(realpart))

        if (expn .ne. expn .or. expn .eq. 0) then 
            expn = 1D1**200
        end if

        retval = abs(inpval * expn) /expn

    end function    

    ! Store the truncated universe data
    subroutine store_data(storage, data, t, d, xmin, ymin, max_points)
        type(DataPoint), dimension(0:), intent(inout) :: storage
        REAL(kind=REAL64), intent(in) :: t, d, xmin, ymin
        INTEGER, intent(in) :: max_points
        COMPLEX(kind=REAL64), dimension(:, :), intent(in) :: data

        Integer, dimension(1) :: storage_dims
        Integer, dimension(2) :: data_dims 
        Integer :: storage_offset, stidx, xidx, yidx
        REAL(kind=REAL64) :: x, y
        
        storage_dims = shape(storage)
        data_dims = shape(data)

        storage_offset = NINT(storage(0)%x)

        stidx = storage_offset 

        !$OMP Parallel do default(shared) private(xidx, yidx, x, y)
        do xidx = 1, data_dims(1)

            ! Is the array filled?
            if (storage(0)%t .ne. 0) then 
                cycle
            end if

            ! Did it get filled but we errored? 
            if (storage(0)%x .ge. storage_dims(1) -1) then 
                !$OMP critical (set_end_flag_header)
                storage(0)%t = 1
                !$OMP end critical (set_end_flag_header)
                cycle
            end if

            do yidx = 1, data_dims(2)

                ! Is the array filled?
                if (storage(0)%t .ne. 0) then 
                    exit
                end if

                ! Did it get filled but we errored? 
                if (storage(0)%x .ge. storage_dims(1) -1) then 
                    !$OMP critical (set_end_flag_header)
                    storage(0)%t = 1
                    !$OMP end critical (set_end_flag_header)
                    exit
                end if

                ! Store data
                if(data(xidx, yidx) .ne. 0) then

                    
                    ! Construct the positional points
                    x = d * (xidx - 1) + xmin
                    y = d * (yidx - 1) + ymin

                    ! Did we hit the per-limit cap?
                    if (max_points .gt. 0 .and. stidx - storage_offset .ge. max_points) then
                        !$OMP Critical (swap_min_lock)
                        call swap_min(storage,storage_offset, stidx, DataPoint(t, x, y, data(xidx, yidx)))
                        !$OMP End Critical (swap_min_lock)
                        cycle
                    end if
                    
                    !$OMP Critical (store_data)

                    ! Store it
                    storage(stidx) = DataPoint(t, x, y, data(xidx, yidx))
                    stidx = stidx + 1

                    ! Check if we reached the limits
                    if ( stidx + 1 .ge. storage_dims(1)) then
                        storage(0)%t = 1
                    else if ( max_points .gt. 0 .and. stidx - storage_offset .ge. max_points) then
                        Print *, "Hit max points at time: ", t
                        Print *, "Now switching to min swap between points ", storage_offset, " and ", stidx
                    end if

                    !$OMP End critical (store_data)
                end if
            end do
        end do
        !$OMP end parallel do

        storage(0)%x = stidx 
    end subroutine 

    subroutine swap_min(subarray, begin, end, ndata)
        type(DataPoint), dimension(:), intent(inout) :: subarray
        Integer, intent(in) :: begin, end
        type(DataPoint), intent(in) :: ndata

        Integer :: minIdx, thread_min_idx, dataidx, per_thread_count
        Integer :: num_threads, cur_thread, start_idx, end_idx
        Integer, allocatable, dimension(:) :: min_per_thread

        if (ndata%value  .eq. 0) then 
            return 
        end if

        !$OMP Parallel default(shared) 
        if (omp_get_thread_num() .eq. 0) then
            num_threads = omp_get_num_threads()
        end if
        !$OMP end parallel
        
        allocate(min_per_thread(num_threads))
        per_thread_count = (end - begin)/num_threads


        !$OMP Parallel default(shared) private(dataidx, cur_thread, start_idx, end_idx, thread_min_idx)
        cur_thread = omp_get_thread_num()

        start_idx = begin + cur_thread * num_threads
        if(cur_thread .eq. num_threads - 1) then
            end_idx = end
        else 
            end_idx = begin + (cur_thread + 1) * num_threads
        end if

        thread_min_idx = start_idx

        do dataidx = start_idx + 1, end_idx
            if (custom_abs(subarray(thread_min_idx)%value) .gt. custom_abs(subarray(dataidx)%value)) then
                thread_min_idx = dataidx
            end if
        end do

        min_per_thread(cur_thread + 1) = thread_min_idx
        !$OMP END Parallel

        minIdx = min_per_thread(1)

        do cur_thread = 2, num_threads 
            if (custom_abs(subarray(minIdx)%value) .lt. custom_abs(subarray(min_per_thread(cur_thread))%value)) then
                minIdx = min_per_thread(cur_thread)
            end if
        end do 

        if(custom_abs(subarray(minIdx)%value) .lt. custom_abs(ndata%value)) then
            subarray(minIdx) = ndata
        end if

        deallocate(min_per_thread)

    end subroutine

end module Wave_Packet_Simulation

program main

    use Wave_Packet_Simulation
    
    implicit none

    ! Constants of reality 
    REAL(kind=REAL64), parameter :: PI  = 4.D0 * DATAN(1.D0) ! Gurantees max accuracy
    REAL(kind=REAL64), parameter :: hbar = 1.054571800D-34!J * s

    ! Constants of this simulation

    ! Data storage variables
    Integer, parameter :: timesize = 2000, points_per_time = 100000, storage_sz = points_per_time * timesize
    Integer, parameter :: max_points = 400000
    Integer, parameter :: xsize = 1000, ysize = 1000 ! Universe size
    COMPLEX(kind=REAL64), dimension(:, :), allocatable, target :: universe_a, universe_b, universe_c ! The closest time slices
    COMPLEX(kind=REAL64), dimension(:, :), pointer :: prev_data, curr_data, next_data, tmp_ptr ! Pointers to each slice to more efficiently rotate

    ! Universe bounds
    REAL(kind=REAL64), parameter :: xmin = -1.0, xmax = 1.0, ymin = -1.0, ymax = 1.0
    REAL(kind=REAL64), parameter :: tmax = 5.0


    type(DataPoint), dimension(:), allocatable :: storage


    REAL(kind=REAL64), parameter :: s = (tmax)/timesize ! time step, seconds
    REAL(kind=REAL64), parameter :: d = MAX((xmax - xmin)/xsize, (ymax -ymin)/ysize) ! distance step, meters; equal to equal to (xmax - xmin)/xsize

    REAL(kind=REAL64), parameter :: px = 1.0, py = 0 ! initial momentum, m/s
    REAL(kind=REAL64), parameter :: sigma = 4D-4 ! the width in meters?

    REAL(kind=REAL64), parameter :: m = 9.10938356D-31 ! Mass of particle (electron here) in kg

    Integer, parameter :: fd = 1 ! Output CSV file descriptor



    ! Local variables
    Integer :: xind, yind, tind ! Iteration indices
    REAL(kind=REAL64) :: x, y, t ! Temp variables for converting from the index
    REAL(kind=REAL64) :: curV ! Temp variable for the potential at a point
    COMPLEX(kind=REAL64) :: psi_coeff, psiV_coeff, psi_off_coeff ! Coefficients for the iterator

    allocate(storage(0:storage_sz))
    allocate(universe_a(xsize, ysize))
    allocate(universe_b(xsize, ysize))
    allocate(universe_c(xsize, ysize))
    
    Print *, "Hello!"
    Print *, "Going from t=0 to t=", timesize * s, " with steps of ", s, & 
             " over a region of [",xmin, ", ", xmin + d*xsize, "] x [", ymin, ", ", ymin + d*ysize, &
             "] with steps of ", d

    storage(0:) = DataPoint(0, 0, 0, (0, 0))
    storage(0) = DataPoint(0, 1, 0, (0, 0))
    
    tmp_ptr => universe_a
    prev_data => universe_a ! The universe at time t+s
    curr_data => universe_b ! The universe at time t
    next_data => universe_c ! The universe at time t+s


    !$omp parallel do default(shared) private(x, y, xind, yind)
    do xind = 1,xsize
        x = (xind - 1) * d + xmin
        do concurrent (yind = 1:ysize)
            y = (yind - 1) * d + ymin
            prev_data(xind, yind) = start_wave_packet(x, y, px, py, sigma, hbar, pi)
        end do
    end do
    !$omp end parallel do

    call normalize_matrix(prev_data)
    call store_data(storage, prev_data, 0D0, d, xmin, ymin, max_points)

    Print *, "Finished t = ", 0, " step"
    
    !$omp parallel do default(shared) private(x, y, xind, yind, curV)
    do xind = 1,xsize
        x = (xind - 1) * d + xmin
        do concurrent (yind = 1:ysize)
            y = (yind - 1) * d + ymin
            curV = V(x, y, s)
            curr_data(xind, yind) = first_iteration_at_point(xind, yind, s, d, hbar, m, curV, prev_data)
        end do
    end do
    !$omp end parallel do

    call normalize_matrix(curr_data)

    call store_data(storage, curr_data, s, d, xmin, ymin, max_points)

    Print *, "Finished t = ", 1, " step"

    psi_coeff = calc_psi_coeff(m, hbar, s, d)
    psiV_coeff = calc_psiV_coeff(m, hbar, s, d)
    psi_off_coeff = calc_psiOff_coeff(m, hbar, s, d)
    
    tind = 3
    do while(NINT(storage(0)%t) .eq. 0 .and. t .lt. tmax)
        t = (tind - 1) * s

        !$omp parallel do default(shared) private(x, y, xind, yind, curV)
        do xind = 1,xsize
            x = (xind - 1) * d + xmin
            do concurrent (yind = 1:ysize)
                y = (yind - 1) * d + ymin
                curV = V(x, y, t)
                next_data(xind, yind) = iterate_psi_at_point( &
                    xind, yind, psi_coeff, psiV_coeff, psi_off_coeff, curV,  &
                    prev_data, curr_data &
                )
            end do
        end do
        !$omp end parallel do

        call normalize_matrix(next_data)
        call store_data(storage, next_data, t, d, xmin, ymin, max_points)
        tmp_ptr => prev_data
        prev_data => curr_data
        curr_data => next_data
        next_data => tmp_ptr
        Print *, "Finished t = ", tind -1, " step"
        tind = tind + 1
    end do 

    ! Write the data to a disk
    ! To save space we only write the non-zero data, not the entire 
    ! universe at every time point
    
    open(fd, file="Simulationdata.csv")

    do xind = 1, storage_sz
        if(storage(xind)%value .eq. 0) then
            exit 
        end if 
        Write (fd, *) storage(xind)%t, ", ", &
                      storage(xind)%x, ", ", &
                      storage(xind)%y, ", ", &
                      real(storage(xind)%value), ",", &
                      aimag(storage(xind)%value)
    end do

    deallocate(universe_a)
    deallocate(universe_b)
    deallocate(universe_c)
    deallocate(storage)

    contains

    ! The potential function 
    pure function V(x, y, t) result(retval)
        REAL(kind=REAL64), intent(in) :: x, y, t
        REAL(kind=REAL64) :: retval

        ! BIG potential at the edges of the box
        if (t .ge. 0 .and. (x .gt. 1 .or. y .gt. 1 .or. x .lt. -1 .or. y .lt. -1)) then
            retval = 9D99
        else

            ! Right now we have no potential inside
            retval = 0
        end if

    end function V



end program main