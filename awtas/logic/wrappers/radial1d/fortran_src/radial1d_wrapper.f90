module radial1d_wrapper

  use iso_c_binding, only: c_double, c_int
  use call_radial1d, only: radial1d

  implicit none

  contains

  subroutine c_radial1d(Porosity, Permeability, LayerThickness, ActionWellRadius,&
    RockSpecificHeat, RockHeatConductivity, RockDensity, RockCompressibility,&
    InitialPressure, InitialX, InjectionWell, InjectionEnthalpy,NumPumpTimes,&
    NumObservationPoints, TotalNumData, PumpingScheme,&
    PumpTime, PumpRate, Time, ObsPointRadialLocation, ObsPointNumData,&
    ObsPointProperty,Deliverability, ProductionIndex, CutoffPressure,&
    NumGridBlocks, NumConstantGridBlocks, ConstantGridBlockSize,&
    GridBlockGrowthFactor, ModelledValue, StatusFlag) bind(c)

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
    ! real(c_double), intent(in) :: MassFlowrate, FlowDuration ! Flow parameters for step flow and constant flow
    real(c_double), intent(in), dimension(NumPumpTimes) :: PumpTime, PumpRate ! Flow parameters for measured flows
    real(c_double), intent(in), dimension(TotalNumData) :: Time
    real(c_double), intent(in), dimension(NumObservationPoints) :: ObsPointRadialLocation
    integer(c_int), intent(in), dimension(NumObservationPoints) :: ObsPointNumData,&
                                                                   ObsPointProperty ! Observation point info arrays
    integer(c_int), intent(in) :: Deliverability
    real(c_double), intent(in) :: ProductionIndex, CutoffPressure
    integer(c_int), intent(in) :: NumGridBlocks, NumConstantGridBlocks
    real(c_double), intent(in) :: ConstantGridBlockSize, GridBlockGrowthFactor


    ! Outputs
    integer(c_int), intent(out) :: StatusFlag ! Flag which signals whether the program executed properly or not. Flag = 0 if successful run, 1 if thermodynamic errors, 2 if failure due to too many timestep reductions.
    real(c_double), intent(out), dimension(TotalNumData) :: ModelledValue

    ! ! Locals
    ! logical(4) :: InjectionWell_logical
    ! logical(4) :: Deliverability_logical
    integer(c_int) :: i
    real(c_double) :: start_time1, end_time, start_time2
    call cpu_time(start_time1)
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
    
    ! ! Check the time array was imported through cython correctly
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

    call cpu_time(start_time2)
    call radial1d(Porosity, Permeability, LayerThickness, ActionWellRadius,&
          RockSpecificHeat, RockHeatConductivity, RockDensity, RockCompressibility,&
          InitialPressure, InitialX, InjectionWell, InjectionEnthalpy,NumPumpTimes,&
          NumObservationPoints, TotalNumData, PumpingScheme,&
          PumpTime, PumpRate, Time, ObsPointRadialLocation, ObsPointNumData,&
          ObsPointProperty,Deliverability, ProductionIndex, CutoffPressure,&
          NumGridBlocks, NumConstantGridBlocks, ConstantGridBlockSize, GridBlockGrowthFactor,&
          ModelledValue, StatusFlag)
    call cpu_time(end_time)
    print *, 'Overall radial1d function call time: ', end_time-start_time2
    print *, 'Overall radial1d_wrapper (inside) time: ', end_time-start_time1
  end subroutine c_radial1d

end module radial1d_wrapper
