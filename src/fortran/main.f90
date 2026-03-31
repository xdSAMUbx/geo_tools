program main
    use iso_fortran_env
    use prueba
    implicit none
    integer :: n = 2
    call prueba_1(n)
    write(*,*) n
end program main