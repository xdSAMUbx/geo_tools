module prueba
    use iso_fortran_env
    implicit none
contains
    pure subroutine prueba_1(n)
        integer, intent(inout) :: n
        n = 2 * n
    end subroutine prueba_1
end module prueba