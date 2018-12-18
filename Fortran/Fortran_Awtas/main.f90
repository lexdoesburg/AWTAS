! gfortran -o main.exe main.f90

! gfortran -o main.exe main.f90 variable_types.o problem_data.o variable_parameters.o models.o noise.o theis_solution.o utility_functions.o homogeneousporous.o thermodynamics.o numericalsimulator1d_routines.o NumericalSimulator1D.o gammafunction.o matrixsolvers.o modelprogress.o theis_main.o

! gfortran -o only_theis.exe main.f90 theis_main.o variable_types.o problem_data.o variable_parameters.o models.o theis_solution.o utility_functions.o

! variable_types > problem_data > utility_functions > theis_solution > variable_parameters > models > theis_main

! program Main
!   use theis_main
!   use radial1d_main
!   use variable_types
! 	use problem_data
!
!   implicit none
!   real(DP) :: k,nu,phi,rho,c,b,Q0,P0,r,t0,dt,t1
!   integer(I4B) :: numData
!   real(DP), allocatable :: pressure(:)
!   k=0.000000000001
!   nu=0.0001111
!   phi=0.1
!   rho=813.37
!   c=0.001303
!   b=100
!   Q0=-0.005
!   P0=3600000
!   r=0.05
!   t0=0
!   dt=200
!   t1=54000
! 	numData=271
!   allocate(pressure(numData))
!
!   call Theis(k,nu,phi,rho,c,b,Q0,P0,r,t0,dt,t1,numData,pressure)
! 	open(unit=1, file='tester1.txt',status='replace')
! 	write (1,'(e30.25)') pressure
! 	close(1)
! end program Main

! gfortran -o radial_main.exe main.f90 variable_types.o problem_data.o variable_parameters.o models.o noise.o theis_solution.o utility_functions.o homogeneousporous.o thermodynamics.o numericalsimulator1d_routines.o NumericalSimulator1D.o gammafunction.o matrixsolvers.o modelprogress.o radial1d_main.o

program Main
  use radial1d_main
  use variable_types

  implicit none

  real(DP) :: k,phi,Pressure0,X0,rw,thick,CR,COND,RHOR,COMP,ConstRate,distFromWell,t0,dt,t1
  integer(I4B) :: numData
  real(DP), allocatable :: pressure(:)
  k=0.0000000000001
  phi=0.1
  Pressure0=4000000
  X0=200
  rw=0.1
  thick=100
  CR=1000
  COND=2.5
  RHOR=2500
  COMP=0
  ConstRate=5
  distFromWell=0.1 ! distance of observation point from action well
  t0=0
  dt=200
  t1=54000
	numData=271
  allocate(pressure(numData))
  call radial1d(phi,k,Pressure0,X0,rw,thick,CR,COND,RHOR,&
                      COMP,ConstRate,distFromWell,numData,t0,dt,t1,pressure)
	open(unit=1, file='radial_tester2.txt',status='replace')
	write (1,'(e30.25)') pressure
	close(1)
end program Main
