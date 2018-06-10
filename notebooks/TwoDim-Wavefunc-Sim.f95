
program main
implicit none


    ! Constants of reality 
    REAL*8, parameter :: PI  = 4.D0 * DATAN(1.D0) ! Gurantees max accuracy
    REAL*8, parameter :: hbar = 1.054571800D-34!J * s

    ! Constants of this simulation
    Integer, parameter :: timesize = 200, xsize = 1000, ysize = 1000 ! Universe size
    COMPLEX * 16, dimension(timesize, xsize, ysize) :: data ! The universe itself

    REAL*8, parameter :: s = 20D-4 ! time step, seconds
    REAL*8, parameter :: d = 20D-4 ! distance step, meters; equal to equal to (xmax - xmin)/xsize

    REAL*8, parameter :: px = 1.0, py = 0 ! initial momentum, m/s
    REAL*8, parameter :: sigma = 4D-4 ! the width in meters?

    REAL*8, parameter :: m = 9.10938356D-31 ! Mass of particle (electron here) in kg

    Integer, parameter :: fd = 1 ! Output CSV file descriptor



    ! Local variables
    Integer :: xind, yind, tind ! Iteration indices
    REAL*8 :: x, y, t ! Temp variables for converting from the index
    REAL*8 :: curV ! Temp variable for the potential at a point
    COMPLEX*16 :: psi_coeff, psiV_coeff, psi_off_coeff ! Coefficients for the iterator

    do xind = 1,xsize
        x = (xind - 1) * d - 1.0
        !$omp parallel do default(shared) private(y, yind)
        do yind = 1,ysize
            y = (yind - 1) * d - 1.0
            data(1, xind, yind) = start_wave_packet(x, y, px, py, sigma, hbar, pi)
        end do
        !$omp end parallel do
    end do

    call normalize_matrix(data(1, :, :))

    Print *, "Finished t = ", 0, " step"
    
    do xind = 1,xsize
        x = (xind - 1) * d - 1.0
        !$omp parallel do default(shared) private(y, yind, curV)
        do yind = 1,ysize
            y = (yind - 1) * d - 1.0
            curV = V(x, y, s)
            data(2, xind, yind) = first_iteration_at_point(xind, yind, s, d, hbar, m, curV, data(1, :, :))
        end do
        !$omp end parallel do
    end do

    call normalize_matrix(data(2, :, :))

    Print *, "Finished t = ", 1, " step"

    psi_coeff = calc_psi_coeff(m, hbar, s, d)
    psiV_coeff = calc_psiV_coeff(m, hbar, s, d)
    psi_off_coeff = calc_psiOff_coeff(m, hbar, s, d)
    
    do tind = 3, timesize
        t = (tind - 1) * s
        do xind = 1,xsize
            x = (xind - 1) * d - 1.0
            !$omp parallel do default(shared) private(y, yind, curV)
            do yind = 1,ysize
                y = (yind - 1) * d - 1.0
                curV = V(x, y, t)
                data(tind, xind, yind) = iterate_psi_at_point( &
                    xind, yind, psi_coeff, psiV_coeff, psi_off_coeff, curV,  &
                    data(tind - 2, : , :), data(tind - 1, : , :) &
                )
            end do
            !$omp end parallel do
        end do
        call normalize_matrix(data(tind, :, :))
        Print *, "Finished t = ", tind -1, " step"
    end do 

    ! Write the data to a disk
    ! To save space we only write the non-zero data, not the entire 
    ! universe at every time point
    
    open(fd, file="Simulationdata.csv")

    do tind = 1, timesize
        do xind = 1,xsize
            do yind = 1,ysize
                if(data(tind, xind, yind) .ne. 0) then
                    x = (xind - 1) * d - 1.0
                    y = (yind - 1) * d - 1.0
                    t = (tind - 1) * s 
                    Write (fd, *) t, ", ", x, ", ", y, ", ", real(data(tind, xind, yind)), ",",  imag(data(tind, xind, yind))
                end if
            end do
        end do
    end do


contains


    ! Simple wave packet centered at (0,0) with momentum px, py and width sigma
    pure function start_wave_packet(x, y, px, py, sigma, hbar, pi) result (psi_at_point)
        REAL*8, intent(in) :: x, y, px, py, sigma, hbar, pi
        COMPLEX*16 :: psi_at_point
        
        REAL * 8 :: coeff_val
        COMPLEX*16 :: exp_val

        !coeff_val = sqrt(sigma *  sqrt(pi))
        !exp_val = (px * x + py * y) * (0, 1.0)/hbar - (x * x + y * y)/(2 * sigma * sigma) 

        !psi_at_point = coeff_val * exp(exp_val)

        if (x .eq. 0 .and. y .eq. 0) then 
            psi_at_point = 1
        else 
            psi_at_point = 0
        end if

        
    end function start_wave_packet


    ! The potential function 
    pure function V(x, y, t) result(retval)
        REAL*8, intent(in) :: x, y, t
        REAL*8 :: retval

        ! BIG potential at the edges of the box
        if (t .ge. 0 .and. (x .gt. 1 .or. y .gt. 1 .or. x .lt. -1 .or. y .lt. -1)) then
            retval = 9D99
        else

            ! Right now we have no potential inside
            retval = 0
        end if

    end function V

    ! The first step iteration, since we don't have access to t-s at that point
    pure function first_iteration_at_point(xind, yind, s, d, hbar, m, curV, curPsi) result(psi_next)
        Integer, intent(in) :: xind, yind
        Real*8, intent(in) :: s, d, hbar, m, curV
        Complex*16, intent(in), dimension(:, :) :: curPsi

        Complex*16 :: psi_next

        Complex*16 :: offPsi, offPsi_coeff, psiV_coeff, psi_coeff, curPsiPoint
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
        REAL*8, intent(in) :: m, hbar, s, d
        COMPLEX*16 :: retval

        retval = (-4.0, 1.0) * s * hbar / (m * d * d)
        
    end function calc_psi_coeff

    ! Get the coefficient for the \V(x, y, t) * \Psi(x, y, t) term
    pure function calc_psiV_coeff(m, hbar, s, d) result(retval)
        REAL*8, intent(in) :: m, hbar, s, d
        COMPLEX*16 :: retval

        retval = (-2.0, 1.0) * s/hbar
        
    end function calc_psiV_coeff

    ! Get the coefficient for the gradient components summation term
    pure function calc_psiOff_coeff(m, hbar, s, d) result(retval)
        REAL*8, intent(in) :: m, hbar, s, d
        COMPLEX*16 :: retval

        retval = hbar * s * (0.0, 1.0) / (m * d * d)
    end function calc_psiOff_coeff

    ! General function to iterate the wave at a point
    pure function iterate_psi_at_point(xind, yind, psi_coeff, psiV_coeff, psiOff_coeff, curV, psi_prev, psi_cur) result(psi_next)
        INTEGER, intent(in) :: xind, yind
        REAL*8, intent(in) :: curV
        COMPLEX*16, intent(in) :: psi_coeff, psiV_coeff, psiOff_coeff, psi_prev(:, :), psi_cur(:, :)
        COMPLEX*16 :: psi_next

        COMPLEX * 16 :: offPsi
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

    subroutine normalize_matrix(matrix)
        COMPLEX*16, intent(inout), dimension(:, :) :: matrix

        INTEGER, dimension(2) :: dims
        REAL*8 :: sum_squares, old_sum

        dims = SHAPE(matrix)
        sum_squares = 0.0

        ! Get the normalization coefficient at that time
        do xind = 1, dims(1)
            do yind = 1, dims(2)
                old_sum = sum_squares
                sum_squares = sum_squares + 1D100 * custom_abs(matrix(xind, yind))**2
            end do
        end do
        
        sum_squares = sum_squares /1D100

        ! Short Circuit on a bad universe
        if (sum_squares .eq. 0) then
            Print *, "Bad universe. " 
            return 
        else if(abs(sum_squares - 1) .lt. 0.0001) then
            Print *, "SPICY! Sumsq = ", sum_squares
        end if
        
        ! Get the normalization coefficient at that time
        do xind = 1, dims(1)
            do yind = 1, dims(2)
                matrix(xind, yind) = matrix(xind, yind) / sqrt(sum_squares)
            end do
        end do

    end subroutine 

    pure function custom_abs(inpval) result(retval)
        Complex*16, intent(in) :: inpval
        Real*8 :: retval 
        
        Real*8 :: realpart, imgpart
        Real*8 :: expn


        retval = abs(inpval)
        if (retval .ne. 0 .or. inpval .eq. 0) then
            return
        end if

        realpart = REAL(inpval)
        imgpart = IMAG(inpval)

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
end program main