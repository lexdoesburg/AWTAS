program run_radial1d
  !
  ! This program is used to test that the wrapper actually works as expected before moving on to use cython to wrap the code into python.
  !
    use variable_types
    ! use call_radial1d
    use radial1d_wrapper

    implicit none

    real(DP)     :: Porosity, Permeability, LayerThickness, ActionWellRadius
    real(DP)     :: RockSpecificHeat, RockHeatConductivity, RockDensity, RockCompressibility
    real(DP)     :: InitialPressure, InitialX
    integer(I4B) :: InjectionWell, NumPumpTimes, NumObservationPoints, TotalNumData, PumpingScheme
    real(DP)     :: InjectionEnthalpy, MassFlowrate, FlowDuration
    real(DP), allocatable :: PumpTime(:), PumpRate(:), Time(:), ObsPointRadialLocation(:), pressure(:)
    integer(I4B), allocatable :: ObsPointNumData(:), ObsPointProperty(:)
    integer(I4B) :: Deliverability
    real(DP)     :: ProductionIndex, CutoffPressure
    real(DP), allocatable :: ModelledValue(:)

    integer(I4B) :: i
    real(DP) :: start_time, end_time

    Porosity = 0.1
    Permeability = 1d-14
    LayerThickness = 100
    ActionWellRadius = 0.0
    RockSpecificHeat = 1000
    RockHeatConductivity = 2.5
    RockDensity = 2500
    RockCompressibility = 0.0
    InitialPressure = 40d5
    InitialX = 200
    InjectionWell = 0
    NumPumpTimes = 3 ! 1
    NumObservationPoints = 2 ! 
    TotalNumData = 43 !500
    PumpingScheme = 4
    InjectionEnthalpy = 0
    MassFlowrate = -8
    FlowDuration = 20000
    Deliverability = 0
    ProductionIndex = 0
    CutoffPressure = 0

    allocate(PumpTime(NumPumpTimes))
    allocate(PumpRate(NumPumpTimes))
    allocate(Time(TotalNumData*2))
    allocate(pressure(TotalNumData*2))
    allocate(ModelledValue(TotalNumData*2))
    allocate(ObsPointRadialLocation(NumObservationPoints))
    allocate(ObsPointNumData(NumObservationPoints))
    allocate(ObsPointProperty(NumObservationPoints))

    PumpTime(1) = 0.0
    PumpRate(1) = -8
    PumpTime(2) = 10000
    PumpRate(2) = 0
    PumpTime(3) = 40000
    PumpRate(3) = -5
    ObsPointRadialLocation(1) = 0.01   
    ObsPointNumData(1) = TotalNumData
    ObsPointProperty(1) = 1
    ObsPointRadialLocation(2) = 0.01
    ObsPointNumData(2) = TotalNumData
    ObsPointProperty(2) = 1


    open (unit=11,file="fortran_src/time_data.txt",status='old')
    do i = 1,TotalNumData
      read(11,*) Time(i),pressure(i)
    end do
    close(11)
    open (unit=11,file="fortran_src/time_data.txt",status='old')
    do i = TotalNumData+1,TotalNumData*2
      read(11,*) Time(i),pressure(i)
    end do
    close(11)

    call cpu_time(start_time)
    ! ! call radial1d
    ! call Radial1D(Porosity, Permeability, LayerThickness, ActionWellRadius,&
    !                 RockSpecificHeat, RockHeatConductivity, RockDensity, RockCompressibility,&
    !                 InitialPressure, InitialX, InjectionWell, InjectionEnthalpy,NumPumpTimes,&
    !                 NumObservationPoints, TotalNumData*2, PumpingScheme, MassFlowrate, FlowDuration,&
    !                 PumpTime, PumpRate, Time, ObsPointRadialLocation, ObsPointNumData,&
    !                 ObsPointProperty,Deliverability, ProductionIndex, CutoffPressure,&
    !                 ModelledValue)

    ! ! radial1d wrapper
    ! call c_Radial1D(Porosity, Permeability, LayerThickness,&
    !                 ActionWellRadius, RockSpecificHeat, RockHeatConductivity, RockDensity,&
    !                 RockCompressibility, InitialPressure, InitialX, InjectionWell,&
    !                 InjectionEnthalpy, NumPumpTimes, NumObservationPoints, TotalNumData,&
    !                 PumpingScheme, MassFlowrate, FlowDuration, PumpTime, PumpRate, &
    !                 Time, ObsPointRadialLocation, ObsPointNumData, ObsPointProperty,&
    !                 Deliverability, ProductionIndex, CutoffPressure,&
    !                 ModelledValue)
    call c_radial1d(Porosity, Permeability, LayerThickness,&
    ActionWellRadius, RockSpecificHeat, RockHeatConductivity, RockDensity,&
    RockCompressibility, InitialPressure, InitialX, InjectionWell,&
    InjectionEnthalpy, NumPumpTimes, NumObservationPoints, TotalNumData,&
    PumpingScheme, MassFlowrate, FlowDuration,&
    Time,&
    Deliverability, ProductionIndex, CutoffPressure,&
    ModelledValue)

    call cpu_time(end_time)
    
    print '("Time elapsed = ",f6.3," seconds.")', end_time-start_time

    open(unit=1, file='radial_test2.txt',status='replace')
    do i = 1, TotalNumData*2
      write (1,'(e30.25","e30.25)') Time(i),ModelledValue(i)
    end do
    close(1)

end program run_radial1d