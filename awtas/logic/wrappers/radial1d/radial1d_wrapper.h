extern void c_radial1d(double* Porosity, double* Permeability, double* LayerThickness,
                        double* ActionWellRadius, double* RockSpecificHeat, double* RockHeatConductivity, double* RockDensity,
                        double* RockCompressibility, double* InitialPressure, double* InitialX, int* InjectionWell,
                        double* InjectionEnthalpy, int* NumPumpTimes, int* NumObservationPoints, int* TotalNumData,
                        int* PumpingScheme, double* PumpTime, double* PumpRate,
                        double* Time, double* ObsPointRadialLocation, int* ObsPointNumData, int* ObsPointProperty,
                        int* Deliverability, double* ProductionIndex, double* CutoffPressure, int* NumGridBlocks,
                        int* NumConstantGridBlocks, double* ConstantGridBlockSize, double* GridBlockGrowthFactor,
                        double* ModelledValue, int* StatusFlag);
                        