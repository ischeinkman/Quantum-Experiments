
module Wave_Packet_Simulation

    use iso_fortran_env, only : REAL64
    use Wavefunc_DataStorage
    use omp_lib
    implicit none

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

end module Wave_Packet_Simulation

program main

    use Wave_Packet_Simulation
    
    implicit none

    ! Constants of reality 
    REAL(kind=REAL64), parameter :: PI  = 4.D0 * DATAN(1.D0) ! Gurantees max accuracy
    REAL(kind=REAL64), parameter :: hbar = 1.054571800D-34!J * s

    ! Constants of this simulation

    ! Data storage variables
    type(DataStore) ::  storage
    Integer, parameter :: timesize = 200, points_per_time = 100000, storage_sz = points_per_time * timesize
    Integer, parameter :: max_points = 400000
    Integer, parameter :: xsize = 1000, ysize = 1000 ! Universe size
    COMPLEX(kind=REAL64), dimension(:, :), allocatable, target :: universe_a, universe_b, universe_c ! The closest time slices
    COMPLEX(kind=REAL64), dimension(:, :), pointer :: prev_data, curr_data, next_data, tmp_ptr ! Pointers to each slice to more efficiently rotate

    ! Universe bounds
    REAL(kind=REAL64), parameter :: xmin = -1.0, xmax = 1.0, ymin = -1.0, ymax = 1.0
    REAL(kind=REAL64), parameter :: tmax = 1.0
    REAL(kind=REAL64), parameter :: min_percent = 1D-6


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
    type(TimeNode), pointer :: data_idx
    
    storage%size = 0
    Nullify(storage%head)
    Nullify(storage%tail)

    allocate(universe_a(xsize, ysize))
    allocate(universe_b(xsize, ysize))
    allocate(universe_c(xsize, ysize))
    
    Print *, "Hello!"
    Print *, "Going from t=0 to t=", timesize * s, " with steps of ", s, & 
             " over a region of [",xmin, ", ", xmin + d*xsize, "] x [", ymin, ", ", ymin + d*ysize, &
             "] with steps of ", d

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
    call store_data(storage, prev_data, 0D0, d, xmin, ymin, max_points, min_percent)

    Print *, "Finished t = ", 0, " step with ", storage%size, " total points."
    
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

    call store_data(storage, curr_data, s, d, xmin, ymin, max_points, min_percent)

    Print *, "Finished t = ", 1, " step with ", storage%size, " total points."

    psi_coeff = calc_psi_coeff(m, hbar, s, d)
    psiV_coeff = calc_psiV_coeff(m, hbar, s, d)
    psi_off_coeff = calc_psiOff_coeff(m, hbar, s, d)
    
    tind = 3
    t = s
    do while(storage%size .le. storage_sz .and. t .lt. tmax)
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
        call store_data(storage, next_data, t, d, xmin, ymin, max_points, min_percent)
        tmp_ptr => prev_data
        prev_data => curr_data
        curr_data => next_data
        next_data => tmp_ptr
        if(MOD(tind-1, 100) .eq. 0) then 
            Print *, "Finished t = ", tind-1, " step with ", storage%size, " total data points."
        end if
        tind = tind + 1
    end do 

    ! Write the data to a disk
    ! To save space we only write the non-zero data, not the entire 
    ! universe at every time point
    
    open(fd, file="Simulationdata.csv")



    data_idx => storage%head
    do while(ASSOCIATED(data_idx))

        do xind = 1, size(data_idx%data)
            Write (fd, *) data_idx%t, ", ", &
                        data_idx%data(xind)%x, ", ", &
                        data_idx%data(xind)%y, ", ", &
                        real(data_idx%data(xind)%value), ",", &
                        aimag(data_idx%data(xind)%value)
        end do 
        data_idx => data_idx%next 
    end do

    deallocate(universe_a)
    deallocate(universe_b)
    deallocate(universe_c)

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