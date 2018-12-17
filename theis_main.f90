! ! For theis solution
! ! gfortran -o pseudodata pseudodata.f90 variable_types.o problem_data.o variable_parameters.o models.o noise.o theis_solution.o utility_functions.o
!
! ! For homogeneous porous simulator
! ! gfortran -o theis_test.exe theis_main.f90 variable_types.o problem_data.o variable_parameters.o models.o noise.o theis_solution.o utility_functions.o homogeneousporous.o thermodynamics.o numericalsimulator1d_routines.o NumericalSimulator1D.o gammafunction.o matrixsolvers.o modelprogress.o

module theis_main

  ! use iso_c_binding
  use variable_types
  use problem_data
  use variable_parameters
  use models

  implicit none

  contains

  subroutine theis(k,nu,phi,rho,c,b,Q0,P0,r,numData,time,pressure)
    ! Arguments
    real(DP), intent(in) :: k,nu,phi,rho,c,b,Q0,P0,r
    integer(I4B), intent(in) :: numData
    real(DP), dimension(numData), intent(in) :: time
    real(DP), dimension(numData), intent(out) :: pressure
    ! Locals
    integer(I4B) :: MaxTotalNData
    integer(I4B) :: ObsPointNo,i,DataNo
    real(DP), allocatable :: variable(:)
    type(datapoint),allocatable :: ReadData(:)
    integer(I4B) :: lower,upper

    ! external updatemodelprogress

    MaxTotalNData=10000

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
    ObsPoint(1)%Position%x(2)=0.0  ! Not used.
    ObsPoint(1)%Position%x(3)=0.0  ! Not used.
    
    DataNo=0
    TotalNData=0
    ObsPoint(1)%Property=1
    ObsPoint(1)%DataOffset=0.0
    ObsPoint(1)%Error=0.0
    ObsPoint(1)%DataIndex=1
    ObsPoint(1)%NData=0

    do i=1,numData
      DataNo=DataNo+1
      if (DataNo<=MaxTotalNData) then
        ReadData(DataNo)%time=time(i)
        ObsPoint(1)%NData=ObsPoint(1)%NData+1
      else
        ! write (*,'(a,i5)') 'Too many data- can handle only',MaxTotalNData
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
    ! call DestroyProblemDataArrays
    deallocate(Pump)
    deallocate(ObsPoint)
    deallocate(ReadData)
    deallocate(PumpData)
    deallocate(FixedParameter)
    deallocate(ReservoirCondition)
    deallocate(variable)
    deallocate(TestData)
    return
  end subroutine theis

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
