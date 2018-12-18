
! For homogeneous porous simulator
! gfortran -o pseudodata pseudodata.f90 variable_types.o problem_data.o variable_parameters.o models.o noise.o theis_solution.o utility_functions.o homogeneousporous.o thermodynamics.o numericalsimulator1d_routines.o NumericalSimulator1D.o gammafunction.o matrixsolvers.o modelprogress.o
module radial1d_main

  use variable_types
  use problem_data
  use variable_parameters
  use models

  implicit none

  contains

  subroutine radial1d(phi,k,Pressure0,X0,rw,thick,CR,COND,RHOR,&
                      COMP,ConstRate,distFromWell,numData,t0,dt,t1,pressure)
                
    ! Arguments
    real(DP), intent(in) :: phi,k,Pressure0,X0,rw,thick,CR,COND,RHOR,COMP,ConstRate,distFromWell
    real(DP), intent(in) :: t0,dt,t1
    integer(I4B), intent(in) :: numData
    real(DP), dimension(numData), intent(out) :: pressure

    ! Local variables:
    integer(I4B) :: MaxTotalNData
    integer(I4B) :: ObsPointNo,i,DataNo,t
    real(DP), allocatable :: variable(:)
    type(datapoint), allocatable :: ReadData(:)
    integer(I4B) :: lower,upper

    ! external updatemodelprogress

    MaxTotalNData=10000
    NObsPoints=1
    ModelType=1
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
    allocate(PumpSchemeParams(1,1)) ! added by  lex

    Pump(1)%NData=1
    Pump(1)%Scheme=1 ! Constant rate of flow

    variable(1)=k
    variable(2)=phi
    ReservoirCondition(1)=Pressure0
    ReservoirCondition(2)=X0
    FixedParameter(1)=0.0_DP
    FixedParameter(2)=thick
    FixedParameter(3)=rw
    FixedParameter(4)=CR
    FixedParameter(5)=COND
    FixedParameter(6)=RHOR
    FixedParameter(7)=COMP
    ! StepRate deleted
    Pump(1)%Enthalpy=0
    PumpSchemeParams(1,1)=ConstRate ! lex added
	  Pump(1)%StepFlows=.true.


    ! Pump(1)%NData=numData

    DataNo=0
    TotalNData=0

    ObsPoint(1)%Property = 1 ! Property (1:pressure, 2:temperature, 3: enthalpy)
    ObsPoint(1)%DataOffset=0.0
    ObsPoint(1)%Position%x(1)=distFromWell ! Observation point distance from action well
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

    ! Deallocate arrays
    deallocate(Pump)
    deallocate(ObsPoint)
    deallocate(ReadData)
    deallocate(PumpData)
    deallocate(FixedParameter)
    deallocate(ReservoirCondition)
    deallocate(variable)
    deallocate(PumpSchemeParams)
    deallocate(TestData)

    return
  end subroutine radial1d

  !.....................................................................................

  subroutine updatemodelprogress(nummodelruns,timestepsize, &
    progressfraction, analysisstopped)
  ! This is a dummy routine to pass into the fast solver.
    use variable_types
    integer(I4B), intent(in) :: nummodelruns
    real(DP), intent(in) :: timestepsize,progressfraction
    logical(LGT), intent(inout) :: analysisstopped

    return
  end subroutine updatemodelprogress

end module radial1d_main
