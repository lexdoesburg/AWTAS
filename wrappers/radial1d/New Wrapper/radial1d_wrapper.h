extern void c_radial1d(double* Porosity, double* Permeability, double* LayerThickness,
                        double* ActionWellRadius, double* RockSpecificHeat, double* RockHeatConductivity, double* RockDensity,
                        double* RockCompressibility, double* InitialPressure, double* InitialX, int* InjectionWell,
                        double* InjectionEnthalpy, int* NumPumpTimes, int* NumObservationPoints, int* TotalNumData,
                        int* PumpingScheme, double* MassFlowrate, double* FlowDuration, double* PumpTime, double* PumpRate,
                        double* Time, double* ObsPointRadialLocation, int* ObsPointNumData, int* ObsPointProperty,
                        int* Deliverability, double* ProductionIndex, double* CutoffPressure, double* ModelledValue);
// extern void c_radial1d(double* Porosity, double* Permeability, double* LayerThickness,
//                         double* ActionWellRadius, double* RockSpecificHeat, double* RockHeatConductivity, double* RockDensity,
//                         double* RockCompressibility, double* InitialPressure, double* InitialX, int* InjectionWell,
//                         double* InjectionEnthalpy, int* NumPumpTimes, int* NumObservationPoints, int* TotalNumData,
//                         int* PumpingScheme, double* MassFlowrate, double* FlowDuration,
//                         double* Time,
//                         int* Deliverability, double* ProductionIndex, double* CutoffPressure, double* ModelledValue);