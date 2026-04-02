module physical_geodesy
    use iso_fortran_env
    implicit none
    !> Exponente para la base B^IND
    integer, parameter :: IND = 960

    !> Parametros básicos para los polinomios asociados
    !>  de legendre
    real(kind=8), parameter :: PI = 3.141592653589793_8
    real(kind=8), parameter :: ROOT3 = 1.7320508075688772_8

    !> Parametros del método Fukushima
    real(kind=8), parameter :: BIG = 2.0_8**IND
    real(kind=8), parameter :: BIGI = 2.0_8**(-IND)
    real(kind=8), parameter :: BIGS = 2.0_8 ** (IND/2)
    real(kind=8), parameter :: BIGSI = 2.0_8 ** (-IND/2)

contains

    subroutine x2f(x, ix, P)
        real(kind=8), intent(in) :: x
        integer(kind=4), intent(in) :: ix
        real(kind=8), intent(out) :: P

        if (ix == 0) then
            P = x 
        else if (ix > 0 ) then
            P = x * BIG ** ix
        else 
            P = x * BIGI ** ix 
        end if
    end subroutine x2f

    subroutine xnorm(x, ix)
        real(kind=8), intent(inout) :: x  
        integer(kind=4), intent(inout) :: ix

        if (abs(x) >= BIGS .and. x/= 0.0_8) then
            x = x * BIGI; ix = ix + 1
        else if (abs(x) <= BIGSI .and. x/=0.0_8) then
            x = x * BIG; ix = ix -1
        end if
    end subroutine xnorm

    subroutine xl2sum(f,g,x,ix,y,iy,z,iz)
        real(kind=8), intent(in) :: f,g,x,y
        integer(kind=4), intent(in) :: ix,iy
        real(kind=8), intent(out) :: z
        integer(kind=4), intent(out) :: iz
        integer(kind=4) :: id

        id = ix - iy

        if (id == 0) then
            z = f * x + g * y; iz = ix
        else if (id == 1) then
            z = f * x + g * (y*BIGI); iz = ix
        else if (id == -1) then
            z = g * y + f * (x*BIGI); iz = iy
        else if (id > 1) then
            z = f * x; iz = ix
        else 
            z = g * y; iz = iy
        end if
        call xnorm(z, iz)
    end subroutine

    subroutine FCM_AFL_FUKU(Nmax, teta, P)
        integer(kind=4), intent(in) :: Nmax
        real(kind=8), intent(in) :: teta
        real(kind=8), allocatable, intent(out) :: P(:,:)

        integer(kind=4) :: n, m, ix, ix1, ix2, iz
        real(kind=8) :: a, b, ct, st, x, x1, x2, z

        st = sin(teta * PI/180.0_8)
        ct = cos(teta * PI/180.0_8)

        allocate(P(0:Nmax,0:Nmax))
        P = 0.0_8

        call x2f(1.0_8, 0, P(0,0))

        do m = 0, Nmax
            if (m == 0) then
                x = 1.0_8; ix = 0
                call x2f(x, ix, P(0,0))
            else if (m == 1) then
                x = ROOT3 * ct; ix = 0
                call x2f(x, ix, P(1,1))
            else
                x = sqrt((2.0_8 * m + 1) / (2.0_8 * m)) &
                    * ct * x
                call xnorm(x,ix)
                call x2f(x, ix, P(m,m))
            end if

            x1 = x; ix1 = ix
            if (m < Nmax) then
                x1 = sqrt(2.0_8 * m + 3) * st * x1
                call xnorm(x1, ix1)
                call x2f(x1, ix1, P(m+1,m))
            end if

            x2 = x; ix2 = ix
            do n = m+2, Nmax
                a = sqrt(((2.0_8*n - 1)*(2.0_8*n + 1)) / &
                    ((n-m)*(n+m)))
                b = sqrt(((2.0_8*n + 1)*(n+m-1)*(n-m-1)) / &
                    ((n-m)*(n+m)*(2.0_8*n - 3)))

                call xl2sum(a*st,-b,x1,ix1,x2,ix2, z, iz)
                call x2f(z,iz,P(n,m))
                x2 = x1 ; ix2 = ix1
                x1 = z; ix1 = iz
            end do
        end do
    end subroutine FCM_AFL_FUKU

    pure subroutine grid(ini, fin, steps, vec)
        real(kind=8), intent(in) :: ini, fin
        integer(kind=4), intent(in) :: steps
        real(kind=8), intent(out) :: vec(:)

        integer :: n
        real(kind=8) :: dist

        dist = (fin - ini)/real(steps, kind=8)
        do concurrent(n = 0 : steps)
            vec(n+1) = real(ini + n * dist,kind=8)
        end do
    end subroutine grid
end module physical_geodesy