module radial1d_wrapper

  use iso_c_binding, only: c_double, c_int
  use radial1d_main, only: radial1d

  implicit none

  contains

  subroutine c_radial1d(phi,k,Pressure0,X0,rw,thick,CR,COND,RHOR,&
                        COMP,ConstRate,distFromWell,numData,time,pressure) bind(c)
    real(c_double), intent(in) :: phi,k,Pressure0,X0,rw,thick,CR,COND,RHOR,COMP,ConstRate,distFromWell
    integer(c_int), intent(in) :: numData
    real(c_double), dimension(numData), intent(in) :: time
    real(c_double), dimension(numData), intent(out) :: pressure
    call radial1d(phi,k,Pressure0,X0,rw,thick,CR,COND,RHOR,&
                  COMP,ConstRate,distFromWell,numData,time,pressure)
  end subroutine c_radial1d

end module radial1d_wrapper

! gfortran -o only_theis.exe main.f90 theis_main.o variable_types.o problem_data.o variable_parameters.o models.o theis_solution.o utility_functions.o
