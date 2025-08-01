module mpix_harmonize_wrapper
    use iso_c_binding
    implicit none

interface

    subroutine mpix_harmonize_f2c(comm, flag) bind(C, name='MPIX_Harmonize_f2c')
      import
      integer(C_INT), VALUE, INTENT(IN) :: comm
      integer(C_INT), INTENT(OUT) :: flag
    end subroutine mpix_harmonize_f2c

end interface

contains
    subroutine mpix_harmonize_f(comm, flag)
      integer, INTENT(IN) :: comm
      integer, OPTIONAL, INTENT(OUT) :: flag
      integer :: c_flag

      call mpix_harmonize_f2c(comm, c_flag)
      if (present(flag)) flag = c_flag

    end subroutine mpix_harmonize_f

end module mpix_harmonize_wrapper
