! gfortran -o radial.exe radial1d_executable.f90 variable_types.o problem_data.o variable_parameters.o models.o noise.o theis_solution.o utility_functions.o homogeneousporous.o thermodynamics.o numericalsimulator1d_routines.o NumericalSimulator1D.o gammafunction.o matrixsolvers.o modelprogress.o radial1d_main.o

program Radial1d_executable
  use radial1d_main
  use variable_types

  implicit none
  real(DP) :: k,phi,Pressure0,X0,rw,thick,CR,COND,RHOR,COMP,ConstRate,distFromWell
  integer(I4B) :: numData
  real(DP), allocatable :: time(:)
  real(DP), allocatable :: pressure(:)
  character(len=60) :: inFile
  character(len=60) :: outFile
  integer(I4B) :: in,out
  print  *, "Radial1d_executable: Get command line arguments"
  call get_command_argument(1,inFile)
  call get_command_argument(2,outFile)

  in=10
  out=11

  print  *, "Radial1d_executable: Open and read input file"
  open(unit=in,file=inFile,status='old')
  read(in,*) k,phi,Pressure0,X0,rw,thick,CR,COND,RHOR,COMP,ConstRate,distFromWell,numData
  allocate(time(numData))  ! further on you only have 21 elements
  read(in,*) time          ! so, read them in
  close(in)


  allocate(pressure(numData))
  print  *, "Radial1d_executable: Call radial1d subroutine from radial1d_main.f90"
  call radial1d(phi,k,Pressure0,X0,rw,thick,CR,COND,RHOR,&
                      COMP,ConstRate,distFromWell,numData,time,pressure)
  print  *, "Radial1d_executable: Open and write output file"
	open(unit=out, file=outFile, status='replace')
	write (out,'(e30.25)') pressure
	close(out)

  print  *, "Radial1d_executable: Deallocate pressure and time arrays"
  deallocate(pressure)
  deallocate(time)
end program Radial1d_executable
