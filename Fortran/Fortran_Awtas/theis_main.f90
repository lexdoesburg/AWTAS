! ! For theis solution
! ! gfortran -o pseudodata pseudodata.f90 variable_types.o problem_data.o variable_parameters.o models.o noise.o theis_solution.o utility_functions.o
!
! ! For homogeneous porous simulator
! ! gfortran -o theis_test.exe theis_main.f90 variable_types.o problem_data.o variable_parameters.o models.o noise.o theis_solution.o utility_functions.o homogeneousporous.o thermodynamics.o numericalsimulator1d_routines.o NumericalSimulator1D.o gammafunction.o matrixsolvers.o modelprogress.o
! program Theis
!
! ! Generates artificial test data by adding random noise to the results
! ! from an AWTAS model.
!
! ! Assume 1 pump, and 1 observation point per observation well.
!
! use variable_types
! use problem_data
! use variable_parameters
! use models
! use noise
!
! implicit none
!
! character(len=60) :: inFile
! character(len=60) :: outFile
!
! ! Global variables:
! integer(I4B) :: out,dat
! integer(I4B) :: MaxTotalNData
! integer(I4B) :: ObsPointNo,i,DataNo
! real(DP) :: t,t0,dt,t1
! real(DP), allocatable :: variable(:)
! real(DP) :: k,nu,phi,rho,c,b,Q0,P0,r
! type(datapoint),allocatable      :: ReadData(:)
! integer(I4B) :: lower,upper
!
! external updatemodelprogress
!
! call get_command_argument(1,inFile)
! call get_command_argument(2,outFile)
!
! dat=10
! out=11
!
! MaxTotalNData=5000
!
! open (unit=out,file=outFile,status='replace')
! NObsPoints=1
! ModelType=0
! ! Problem dimensions:
! NPumps=1
! MaxNPumpTimes=1000
!
! NVariables=2
! NFixedParameters=7
! NReservoirConditions=2
!
! allocate(Pump(NPumps),ObsPoint(NObsPoints))
! allocate(ReadData(MaxTotalNData))
! allocate(PumpData(NPumps,MaxNPumpTimes))
! allocate(FixedParameter(NFixedParameters),ReservoirCondition(NReservoirConditions))
! allocate(variable(NVariables))
!
! Pump(1)%NData=1
! Pump(1)%Scheme=1
!
! open(unit=dat,file=inFile,status='old')
! read(dat,*) k,nu,phi,rho,c,b,Q0,P0,r,t0,dt,t1
! variable(1)=k
! variable(2)=phi
! ReservoirCondition(1)=P0
! FixedParameter(1)=0.0_DP
! FixedParameter(2)=b
! FixedParameter(3)=0.1
! FixedParameter(4)=nu
! FixedParameter(5)=rho
! FixedParameter(6)=c
! PumpData(1,1)%rate=Q0
! ObsPoint(1)%Position%x(1)=r
! close(dat)
!
! DataNo=0
! TotalNData=0
! ObsPoint(1)%Property=1
! ObsPoint(1)%DataOffset=0.0
! ObsPoint(1)%Position%x(2)=0.0  ! Not used.
! ObsPoint(1)%Position%x(3)=0.0  ! Not used.
! ObsPoint(1)%Error=0.0
! ObsPoint(1)%DataIndex=1
! ObsPoint(1)%NData=0
!
! do t=t0,t1,dt
!   DataNo=DataNo+1
!   if (DataNo<=MaxTotalNData) then
!     ReadData(DataNo)%time=t
!     ObsPoint(1)%NData=ObsPoint(1)%NData+1
!   else
!     write (*,'(a,i5)') 'Too many data- can handle only',MaxTotalNData
!     stop
!   end if
! end do
! TotalNData=TotalNData+ObsPoint(1)%NData
! allocate(TestData(TotalNData))
! TestData=ReadData(1:TotalNData)
! ! Assign ObsPoint errors to TestData:
! ObsPoint%Weight=1.0_dp
! TestData%Weight=1.0_dp
! lower=ObsPoint(1)%DataIndex
! upper=lower+ObsPoint(1)%NData-1
! TestData(lower:upper)%error=ObsPoint(1)%error
! ! Generate modelled values:
! TestData%ModelledValue=model(variable,updatemodelprogress)
! ! Write file:
! write (out,'(10(e15.5,1x))') (variable(i),i=1,NVariables)
! write (out,'(i15)') NObsPoints
! write (out,*)
! DataNo=0
! write (out,'(i15,1x,i15)') ObsPointNo,ObsPoint(1)%Property
! write (out,'(e15.5,1x,e15.5)') ObsPoint(1)%Position%x(1),ObsPoint(1)%DataOffset
! write (out,'(i15,1x,e15.5)') ObsPoint(1)%NData,ObsPoint(1)%Error
! do i=1,ObsPoint(1)%NData
!   DataNo=DataNo+1
!   write (out,'(e30.25,1x,e30.25)') TestData(DataNo)%time,TestData(DataNo)%ModelledValue
! end do
!
! close(out)
!
! stop
!
! end program Theis
!
! !.....................................................................................
!
! subroutine updatemodelprogress(nummodelruns,timestepsize, &
!   progressfraction, analysisstopped)
! ! This is a dummy routine to pass into the fast solver.
!   use variable_types
!   integer(I4B), intent(in) :: nummodelruns
!   real(DP), intent(in) :: timestepsize,progressfraction
!   logical(LGT), intent(inout) :: analysisstopped
!
!   return
! end subroutine updatemodelprogress


module theis_main

  ! use iso_c_binding
  use variable_types
  use problem_data
  use variable_parameters
  use models

  implicit none

  contains

  subroutine Theis(k,nu,phi,rho,c,b,Q0,P0,r,t0,dt,t1,numData,pressure)
    ! Arguments
    real(DP), intent(in) :: k,nu,phi,rho,c,b,Q0,P0,r,t0,dt,t1
    integer(I4B), intent(in) :: numData
    ! real(DP), intent(in) :: t(numData)
    real(DP), dimension(numData), intent(out) :: pressure
    ! Locals
    integer(I4B) :: MaxTotalNData
    integer(I4B) :: ObsPointNo,i,DataNo
    real(DP), allocatable :: variable(:)
    type(datapoint),allocatable :: ReadData(:)
    integer(I4B) :: lower,upper,t

    ! external updatemodelprogress

    MaxTotalNData=5000

    ! open (unit=out,file=outFile,status='replace')
    NObsPoints=1
    ModelType=0
    ! Problem dimensions:
    NPumps=1
    MaxNPumpTimes=1000

    NVariables=2
    NFixedParameters=7
    NReservoirConditions=2

    allocate(Pump(NPumps),ObsPoint(NObsPoints))
    allocate(ReadData(MaxTotalNData))
    allocate(PumpData(NPumps,MaxNPumpTimes))
    allocate(FixedParameter(NFixedParameters),ReservoirCondition(NReservoirConditions))
    allocate(variable(NVariables))

    Pump(1)%NData=1
    Pump(1)%Scheme=1

    variable(1)=k
    variable(2)=phi
    ReservoirCondition(1)=P0
    FixedParameter(1)=0.0_DP
    FixedParameter(2)=b
    FixedParameter(3)=0.1
    FixedParameter(4)=nu
    FixedParameter(5)=rho
    FixedParameter(6)=c
    PumpData(1,1)%rate=Q0
    ObsPoint(1)%Position%x(1)=r
    ! close(dat)

    DataNo=0
    TotalNData=0
    ObsPoint(1)%Property=1
    ObsPoint(1)%DataOffset=0.0
    ObsPoint(1)%Position%x(2)=0.0  ! Not used.
    ObsPoint(1)%Position%x(3)=0.0  ! Not used.
    ObsPoint(1)%Error=0.0
    ObsPoint(1)%DataIndex=1
    ObsPoint(1)%NData=0

    do t=t0,t1,dt
      DataNo=DataNo+1
      if (DataNo<=MaxTotalNData) then
        ReadData(DataNo)%time=t
        ObsPoint(1)%NData=ObsPoint(1)%NData+1
      else
        write (*,'(a,i5)') 'Too many data- can handle only',MaxTotalNData
        stop
      end if
    end do
    TotalNData=TotalNData+ObsPoint(1)%NData
    allocate(TestData(TotalNData))
    TestData=ReadData(1:TotalNData)
    ! Assign ObsPoint errors to TestData:
    ObsPoint%Weight=1.0_dp
    TestData%Weight=1.0_dp
    lower=ObsPoint(1)%DataIndex
    upper=lower+ObsPoint(1)%NData-1
    TestData(lower:upper)%error=ObsPoint(1)%error
    ! Generate modelled values:
    TestData%ModelledValue=model(variable,updatemodelprogress)
    pressure = TestData%ModelledValue
    return
  end subroutine Theis

  subroutine updatemodelprogress(nummodelruns,timestepsize, &
    progressfraction, analysisstopped)
    use variable_types
    integer(I4B), intent(in) :: nummodelruns
    real(DP), intent(in) :: timestepsize,progressfraction
    logical(LGT), intent(inout) :: analysisstopped
    return
  end subroutine updatemodelprogress

end module theis_main

! gfortran -o only_theis.exe theis_main.f90 variable_types.o problem_data.o variable_parameters.o models.o theis_solution.o utility_functions.o
