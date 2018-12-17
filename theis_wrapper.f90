module theis_wrapper

  use iso_c_binding, only: c_double, c_int
  use theis_main, only: theis

  implicit none

  contains

  subroutine c_theis(k,nu,phi,rho,c,b,Q0,P0,r,numData,time,pressure) bind(c)
    real(c_double), intent(in) :: k,nu,phi,rho,c,b,Q0,P0,r
    integer(c_int), intent(in) :: numData
    real(c_double), dimension(numData), intent(in) :: time
    real(c_double), dimension(numData), intent(out) :: pressure
    call theis(k,nu,phi,rho,c,b,Q0,P0,r,numData,time,pressure)
  end subroutine c_theis

end module theis_wrapper

! gfortran -o only_theis.exe main.f90 theis_main.o variable_types.o problem_data.o variable_parameters.o models.o theis_solution.o utility_functions.o
