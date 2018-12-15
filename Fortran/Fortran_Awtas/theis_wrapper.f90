module theis_wrapper

  use iso_c_binding, only: c_double, c_int
  use theis_main, only: Theis

  implicit none

  contains

  subroutine c_Theis(k,nu,phi,rho,c,b,Q0,P0,r,t0,dt,t1,numData,pressure) bind(c)
    real(c_double), intent(in) :: k,nu,phi,rho,c,b,Q0,P0,r,t0,dt,t1
    integer(c_int), intent(in) :: numData
    ! real(DP), intent(in) :: t(numData)
    real(c_double), dimension(numData), intent(out) :: pressure
    call Theis(k,nu,phi,rho,c,b,Q0,P0,r,t0,dt,t1,numData,pressure)
  end subroutine c_Theis

end module theis_wrapper

! gfortran -o only_theis.exe main.f90 theis_main.o variable_types.o problem_data.o variable_parameters.o models.o theis_solution.o utility_functions.o
