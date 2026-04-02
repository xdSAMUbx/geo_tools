program main
    use iso_fortran_env
    use physical_geodesy
    implicit none
    integer(kind=4) :: n, m, Nmax = 10
    real(kind=8) :: lat
    real(kind=8), allocatable :: P(:,:)

    lat = 10.0_8
    call FCM_AFL_FUKU(Nmax,lat, P)

    do m = 0, Nmax
        do n = m, Nmax
            write(*,'("P (",I3,",",I3," ) =",ES15.5)') n, m, P(n,m)
        end do
    end do
end program main