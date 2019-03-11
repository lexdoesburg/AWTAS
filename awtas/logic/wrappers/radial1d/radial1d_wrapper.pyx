from numpy cimport ndarray
from numpy import empty

cdef extern from "radial1d_wrapper.h":
    extern void c_radial1d(double* Porosity, double* Permeability, double* LayerThickness,
        double* ActionWellRadius, double* RockSpecificHeat, double* RockHeatConductivity, double* RockDensity,
        double* RockCompressibility, double* InitialPressure, double* InitialX, int* InjectionWell,
        double* InjectionEnthalpy, int* NumPumpTimes, int* NumObservationPoints, int* TotalNumData,
        int* PumpingScheme, double* PumpTime, double* PumpRate,
        double* Time, double* ObsPointRadialLocation, int* ObsPointNumData, int* ObsPointProperty,
        int* Deliverability, double* ProductionIndex, double* CutoffPressure, int* NumGridBlocks,
        int* NumConstantGridBlocks, double* ConstantGridBlockSize, double* GridBlockGrowthFactor,
        double* ModelledValue, int* StatusFlag)

# # Using cython memoryviews for input arrays
# def radial1d(double Porosity, double Permeability, double LayerThickness,
#         double ActionWellRadius, double RockSpecificHeat, double RockHeatConductivity, double RockDensity,
#         double RockCompressibility, double InitialPressure, double InitialX, int InjectionWell,
#         double InjectionEnthalpy, int NumPumpTimes, int NumObservationPoints, int TotalNumData,
#         int PumpingScheme, double MassFlowrate, double FlowDuration, double[:] PumpTime, double[:] PumpRate,
#         double[:] Time, double[:] ObsPointRadialLocation, int[:] ObsPointNumData, int[:] ObsPointProperty,
#         int Deliverability, double ProductionIndex, double CutoffPressure):

# Using cimported numpy ndarrays for input arrays
def radial1d(double Porosity, double Permeability, double LayerThickness,
        double ActionWellRadius, double RockSpecificHeat, double RockHeatConductivity, double RockDensity,
        double RockCompressibility, double InitialPressure, double InitialX, int InjectionWell,
        double InjectionEnthalpy, int NumPumpTimes, int NumObservationPoints, int TotalNumData,
        int PumpingScheme, ndarray[double] PumpTime, ndarray[double] PumpRate,
        ndarray[double] Time, ndarray[double] ObsPointRadialLocation, ndarray[int] ObsPointNumData, ndarray[int] ObsPointProperty,
        int Deliverability = 0, double ProductionIndex = 0.0, double CutoffPressure = 0.0, int NumGridBlocks=100, int NumConstantGridBlocks=20,
        double ConstantGridBlockSize=0.01, double GridBlockGrowthFactor=1.2):
    """
    Cython wrapper of the fortran homogeneous porous simulator.
    """
    cdef int StatusFlag
    cdef ndarray[double] ModelledValue = empty(TotalNumData)

    c_radial1d(&Porosity, &Permeability, &LayerThickness,
                &ActionWellRadius, &RockSpecificHeat, &RockHeatConductivity, &RockDensity,
                &RockCompressibility, &InitialPressure, &InitialX, &InjectionWell,
                &InjectionEnthalpy, &NumPumpTimes, &NumObservationPoints, &TotalNumData,
                &PumpingScheme, &PumpTime[0], &PumpRate[0],
                &Time[0], &ObsPointRadialLocation[0], &ObsPointNumData[0], &ObsPointProperty[0],
                &Deliverability, &ProductionIndex, &CutoffPressure, &NumGridBlocks, &NumConstantGridBlocks,
                &ConstantGridBlockSize, &GridBlockGrowthFactor, &ModelledValue[0], &StatusFlag)

    return ModelledValue, StatusFlag