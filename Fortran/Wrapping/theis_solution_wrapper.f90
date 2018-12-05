module theis_solution_wrapper

    use iso_c_binding, only: c_double, c_int
    use theis_solution, only: analyticaltheis

    implicit none

    contains

        subroutine c_theis_solution(k, phi, P0, Q0, h, rho, nu, C, r, nObservations, time, p) bind(c)
            
            real(c_double),intent(in) :: k !
            real(c_double),intent(in) :: phi
            real(c_double),intent(in) :: P0 !
            real(c_double),intent(in) :: Q0 !qm
            real(c_double),intent(in) :: h !
            real(c_double),intent(in) :: rho
            real(c_double),intent(in) :: nu
            real(c_double),intent(in) :: C
            real(c_double),intent(in) :: r
            integer(c_int),intent(in) :: nObservations
            real(c_double),intent(in) :: time(nObservations)
            real(c_double),intent(out) :: p(nObservations)
            integer(c_int) :: i
            real(c_double) :: pressure
            real(c_double) :: t
            
            call analyticaltheis(k, phi, P0, Q0, h, rho, nu, C, r, nObservations, time, p)
        
        end subroutine c_theis_solution

end module theis_solution_wrapper