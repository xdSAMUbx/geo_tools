program main
    use iso_fortran_env
    use physical_geodesy
    implicit none
    integer(kind=4) :: n, m, j, steps = 120, Nmax = 5
    real(kind=8) :: lat_ini, lat_fin, ini, fin
    real(kind=8), allocatable :: P(:,:), pru(:)


    lat_ini = 86.0_8; lat_fin = 85.0_8
    allocate(pru(0:steps), P(0:Nmax,0:Nmax))

    call cpu_time(ini)
    !call grid(lat_ini,lat_fin,steps,pru)
    !do j = 0, steps
    call FCM_AFL_FUKU(Nmax, lat_ini, P)
        !write(*,*) "El valor es: ", pru(j), "Con su AFL", P(Nmax,Nmax)
    !end do
    call cpu_time(fin)
    do m = 0, Nmax
        do n = m, Nmax
            write(*,*) "P(",n,",",m,") = ", P(n,m)
        end do
    end do
end program main