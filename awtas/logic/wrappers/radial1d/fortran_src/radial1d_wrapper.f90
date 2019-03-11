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

    call radial1d(Porosity, Permeability, LayerThickness, ActionWellRadius,&
          RockSpecificHeat, RockHeatConductivity, RockDensity, RockCompressibility,&
          InitialPressure, InitialX, InjectionWell, InjectionEnthalpy,NumPumpTimes,&
          NumObservationPoints, TotalNumData, PumpingScheme,&
          PumpTime, PumpRate, Time, ObsPointRadialLocation, ObsPointNumData,&
          ObsPointProperty,Deliverability, ProductionIndex, CutoffPressure,&
          NumGridBlocks, NumConstantGridBlocks, ConstantGridBlockSize, GridBlockGrowthFactor,&
          ModelledValue, StatusFlag)

  end subroutine c_radial1d

end module radial1d_wrapper
