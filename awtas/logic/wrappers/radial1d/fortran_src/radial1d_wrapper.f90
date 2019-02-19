module radial1d_wrapper

  use iso_c_binding, only: c_double, c_int
  use call_radial1d, only: radial1d

  implicit none

  contains

  subroutine c_radial1d(Porosity, Permeability, LayerThickness,&
      ActionWellRadius, RockSpecificHeat, RockHeatConductivity, RockDensity,&
      RockCompressibility, InitialPressure, InitialX, InjectionWell,&
      InjectionEnthalpy, NumPumpTimes, NumObservationPoints, TotalNumData,&
      PumpingScheme, MassFlowrate, FlowDuration, PumpTime, PumpRate, &
      Time, ObsPointRadialLocation, ObsPointNumData, ObsPointProperty,&
      Deliverability, ProductionIndex, CutoffPressure,&
      ModelledValue) bind(c)

  ! subroutine c_radial1d(Porosity, Permeability, LayerThickness,&
  !   ActionWellRadius, RockSpecificHeat, RockHeatConductivity, RockDensity,&
  !   RockCompressibility, InitialPressure, InitialX, InjectionWell,&
  !   InjectionEnthalpy, NumPumpTimes, NumObservationPoints, TotalNumData,&
  !   PumpingScheme, MassFlowrate, FlowDuration,&
  !   Time,&
  !   Deliverability, ProductionIndex, CutoffPressure,&
  !   ModelledValue) bind(c)

    ! Input arguments
    real(c_double), intent(in) :: Porosity, Permeability ! Variables
    real(c_double), intent(in) :: LayerThickness, ActionWellRadius,& 
                                  RockSpecificHeat, RockHeatConductivity,&
                                  RockDensity, RockCompressibility ! Fixed Parameters
    real(c_double), intent(in) :: InitialPressure, InitialX ! Initial reservoir conditions
    integer(c_int), intent(in) :: InjectionWell ! True if the well is an injection well, false if production
    real(c_double), intent(in) :: InjectionEnthalpy
    integer(c_int), intent(in) :: NumPumpTimes
    integer(c_int), intent(in) :: NumObservationPoints
    integer(c_int), intent(in) :: TotalNumData
    integer(c_int), intent(in) :: PumpingScheme
    real(c_double), intent(in) :: MassFlowrate, FlowDuration ! Flow parameters for step flow and constant flow
    real(c_double), intent(in), dimension(NumPumpTimes) :: PumpTime, PumpRate ! Flow parameters for measured flows
    real(c_double), intent(in), dimension(TotalNumData) :: Time
    real(c_double), intent(in), dimension(NumObservationPoints) :: ObsPointRadialLocation
    integer(c_int), intent(in), dimension(NumObservationPoints) :: ObsPointNumData,&
                                                                   ObsPointProperty ! Observation point info arrays
    integer(c_int), intent(in) :: Deliverability
    real(c_double), intent(in) :: ProductionIndex, CutoffPressure

    ! Outputs
    ! integer(c_int), intent(out) :: ExecutionFlag ! Flag which signals whether the program executed properly or not. Flag = 1 if everything ran fine, 2 if missing input, 3 if simulator did not run as expected
    real(c_double), intent(out), dimension(TotalNumData) :: ModelledValue

    ! ! Locals
    ! logical(4) :: InjectionWell_logical
    ! logical(4) :: Deliverability_logical
    ! integer(c_int) :: i
    ! if (InjectionWell == 0) then
    !   InjectionWell_logical = .False.
    ! else
    !   InjectionWell_logical = .True.
    ! end if

    ! if (Deliverability == 0) then
    !   Deliverability_logical = .False.
    ! else
    !   Deliverability_logical = .True.
    ! end if
    ! print *, Porosity
    ! print *, Permeability
    ! print *, LayerThickness
    ! print *, ActionWellRadius
    ! print *, RockSpecificHeat
    ! print *, RockHeatConductivity
    ! print *, RockDensity
    ! print *, RockCompressibility
    ! print *, InitialPressure
    ! print *, InitialX
    ! print *, InjectionWell
    ! print *, InjectionEnthalpy
    ! print *, NumPumpTimes
    ! print *, NumObservationPoints
    ! print *, TotalNumData
    ! print *, PumpingScheme
    ! print *, MassFlowrate
    ! print *, FlowDuration
    ! do i=1,NumPumpTimes
    !   print *, PumpTime(i)
    ! end do
    ! do i=1,NumPumpTimes
    !   print *, PumpRate(i)
    ! end do
    ! do i=1,NumObservationPoints
    !   print *, ObsPointRadialLocation(i)
    ! end do
    ! do i=1,NumObservationPoints
    !   print *, ObsPointNumData(i)
    ! end do
    ! do i=1,NumObservationPoints
    !   print *, ObsPointProperty(i)
    ! end do
    ! print *, Deliverability
    ! print *, ProductionIndex
    ! print *, CutoffPressure
    ! print *, 'Inside fortran wrapper'
    
    ! print *, 'Printing time'
    ! do i=1,TotalNumData
    !   print *, Time(i)
    ! end do
    ! print *, 'Finished printing time'

    ! print *, 'Going into fortran caller'
    ! real(c_double), dimension(NumPumpTimes) :: PumpRate, PumpTime
    ! real(c_double), dimension(NumObservationPoints) :: ObsPointRadialLocation
    ! integer(c_int), dimension(NumObservationPoints) :: ObsPointNumData, ObsPointProperty
  
    ! PumpRate(1) = 0
    ! PumpTime(1) = 0
    ! ObsPointRadialLocation(1) = 0.01
    ! ObsPointNumData(1) = TotalNumData
    ! ObsPointProperty(1) = 1
    ! print *, Deliverability
    ! print *, ProductionIndex
    ! print *, CutoffPressure

    call radial1d(Porosity, Permeability, LayerThickness,&
        ActionWellRadius, RockSpecificHeat, RockHeatConductivity, RockDensity,&
        RockCompressibility, InitialPressure, InitialX, InjectionWell,&
        InjectionEnthalpy, NumPumpTimes, NumObservationPoints, TotalNumData,&
        PumpingScheme, MassFlowrate, FlowDuration, PumpTime, PumpRate, &
        Time, ObsPointRadialLocation, ObsPointNumData, ObsPointProperty,&
        Deliverability, ProductionIndex, CutoffPressure,&
        ModelledValue)

  end subroutine c_radial1d

end module radial1d_wrapper
