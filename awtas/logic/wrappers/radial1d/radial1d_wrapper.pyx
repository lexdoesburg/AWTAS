from numpy cimport ndarray
from numpy import empty

cdef extern from "radial1d_wrapper.h":
    extern void c_radial1d(double* Porosity, double* Permeability, double* LayerThickness,
        double* ActionWellRadius, double* RockSpecificHeat, double* RockHeatConductivity, double* RockDensity,
        double* RockCompressibility, double* InitialPressure, double* InitialX, int* InjectionWell,
        double* InjectionEnthalpy, int* NumPumpTimes, int* NumObservationPoints, int* TotalNumData,
        int* PumpingScheme, double* MassFlowrate, double* FlowDuration, double* PumpTime, double* PumpRate,
        double* Time, double* ObsPointRadialLocation, int* ObsPointNumData, int* ObsPointProperty,
        int* Deliverability, double* ProductionIndex, double* CutoffPressure, double* ModelledValue)

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
        int PumpingScheme, double MassFlowrate, double FlowDuration, ndarray[double] PumpTime, ndarray[double] PumpRate,
        ndarray[double] Time, ndarray[double] ObsPointRadialLocation, ndarray[int] ObsPointNumData, ndarray[int] ObsPointProperty,
        int Deliverability, double ProductionIndex, double CutoffPressure):

    import time
   
    cdef ndarray[double] ModelledValue = empty(TotalNumData)

    start_time = time.time()
    c_radial1d(&Porosity, &Permeability, &LayerThickness,
                &ActionWellRadius, &RockSpecificHeat, &RockHeatConductivity, &RockDensity,
                &RockCompressibility, &InitialPressure, &InitialX, &InjectionWell,
                &InjectionEnthalpy, &NumPumpTimes, &NumObservationPoints, &TotalNumData,
                &PumpingScheme, &MassFlowrate, &FlowDuration, &PumpTime[0], &PumpRate[0],
                &Time[0], &ObsPointRadialLocation[0], &ObsPointNumData[0], &ObsPointProperty[0],
                &Deliverability, &ProductionIndex, &CutoffPressure, &ModelledValue[0])
    end_time = time.time()
    print('Radial1d wrapper time (from cython): {}'.format(end_time-start_time))
    return ModelledValue