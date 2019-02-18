module radial1d_main

  use variable_types
  use problem_data
  use HomogeneousPorousSimulator

  implicit none

  contains
!-------10--------10--------10--------10--------10--------10--------10--------80
  subroutine radial1d(phi,k,Pressure0,X0,rw,thick,CR,COND,RHOR,&
                      COMP,ConstRate,distFromWell,numData,time,pressure)
    ! Arguments
    real(DP), intent(in) :: phi,k,Pressure0,X0,rw,thick,CR,COND,RHOR,COMP,ConstRate,distFromWell
    ! real(DP), intent(in) :: t0,dt,t1
    integer(I4B), intent(in) :: numData
    real(DP), dimension(numData), intent(in) :: time
    real(DP), dimension(numData), intent(out) :: pressure

    ! Local variables:
    integer(I4B) :: MaxTotalNData
    integer(I4B) :: ObsPointNo,i,DataNo,t
    real(DP), allocatable :: variable(:)
    type(datapoint), allocatable :: ReadData(:)
    integer(I4B) :: lower,upper

    ! external updatemodelprogress
    print  *, "radial1d_main: Assigning helper values"
    MaxTotalNData=10000
    NObsPoints=1
    ModelType=1
    ! Problem dimensions:
    NPumps=1
    MaxNPumpTimes=1000

    NVariables=2
    NFixedParameters=7
    NReservoirConditions=2

    print  *, "radial1d_main: Allocating arrays"
    allocate(Pump(NPumps),ObsPoint(NObsPoints))
    allocate(ReadData(MaxTotalNData))
    allocate(PumpData(NPumps,MaxNPumpTimes))
    allocate(FixedParameter(NFixedParameters),ReservoirCondition(NReservoirConditions))
    allocate(variable(NVariables))
    allocate(PumpSchemeParams(1,1)) ! added by  lex

    print  *, "radial1d_main: Assigning variables/parameters"
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
    print  *, "radial1d_main: Assigning observation point data"
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
    print  *, "radial1d_main: Allocating time in do"
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

    print  *, "radial1d_main: Allocate TestData"
    TotalNData=TotalNData+ObsPoint(1)%NData
    allocate(TestData(TotalNData))
    TestData=ReadData(1:TotalNData)

    ! Assign ObsPoint errors to TestData:
    print  *, "radial1d_main: Assign obspoint errors"
    ObsPoint%Weight=1.0_dp
    TestData%Weight=1.0_dp
    lower=ObsPoint(1)%DataIndex
    upper=lower+ObsPoint(1)%NData-1
    TestData(lower:upper)%error=ObsPoint(1)%error

    ! Generate modelled values:
    print  *, "radial1d_main: Getting modelled values (call model through model.f90)"
    TestData%ModelledValue=homogeneousporous(variable,updatemodelprogress)
    pressure = TestData%ModelledValue

    ! Deallocate arrays
    print  *, "radial1d_main: Deallocating arrays"
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
    ! print  *, "radial1d_main: Inside updatemodelprogress"
    return
  end subroutine updatemodelprogress

end module radial1d_main
