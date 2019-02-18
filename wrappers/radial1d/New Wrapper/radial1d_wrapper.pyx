from numpy cimport ndarray
# from libcpp cimport bool
from numpy import empty

# cdef extern from "radial1d_wrapper.h":
#     extern void c_radial1d(double* Porosity, double* Permeability, double* LayerThickness,
#         double* ActionWellRadius, double* RockSpecificHeat, double* RockHeatConductivity, double* RockDensity,
#         double* RockCompressibility, double* InitialPressure, double* InitialX, int* InjectionWell,
#         double* InjectionEnthalpy, int* NumPumpTimes, int* NumObservationPoints, int* TotalNumData,
#         int* PumpingScheme, double* MassFlowrate, double* FlowDuration, double* PumpTime, double* PumpRate,
#         double* Time, double* ObsPointRadialLocation, int* ObsPointNumData, int* ObsPointProperty,
#         int* Deliverability, double* ProductionIndex, double* CutoffPressure, double* ModelledValue)

# def radial1d(double Porosity, double Permeability, double LayerThickness,
#         double ActionWellRadius, double RockSpecificHeat, double RockHeatConductivity, double RockDensity,
#         double RockCompressibility, double InitialPressure, double InitialX, int InjectionWell,
#         double InjectionEnthalpy, int NumPumpTimes, int NumObservationPoints, int TotalNumData,
#         int PumpingScheme, double MassFlowrate, double FlowDuration, ndarray[double] PumpTime, ndarray[double] PumpRate,
#         ndarray[double] Time, ndarray[double] ObsPointRadialLocation, ndarray[int] ObsPointNumData, ndarray[int] ObsPointProperty,
#         int Deliverability, double ProductionIndex, double CutoffPressure):

   
#     cdef ndarray[double] ModelledValue = empty(TotalNumData)
#     print('Printing time from cython')
#     print(Time)
#     print('Finished printing time from cython')
#     print('Going into fortran function')
#     c_radial1d(&Porosity, &Permeability, &LayerThickness,
#                 &ActionWellRadius, &RockSpecificHeat, &RockHeatConductivity, &RockDensity,
#                 &RockCompressibility, &InitialPressure, &InitialX, &InjectionWell,
#                 &InjectionEnthalpy, &NumPumpTimes, &NumObservationPoints, &TotalNumData,
#                 &PumpingScheme, &MassFlowrate, &FlowDuration, &PumpTime[0], &PumpRate[0],
#                 &Time[0], &ObsPointRadialLocation[0], &ObsPointNumData[0], &ObsPointProperty[0],
#                 &Deliverability, &ProductionIndex, &CutoffPressure, &ModelledValue[0])

#     return ModelledValue

# cdef extern from "radial1d_wrapper.h":
#     extern void c_radial1d(double* Porosity, double* Permeability, double* LayerThickness,
#         double* ActionWellRadius, double* RockSpecificHeat, double* RockHeatConductivity, double* RockDensity,
#         double* RockCompressibility, double* InitialPressure, double* InitialX, int* InjectionWell,
#         double* InjectionEnthalpy, int* NumPumpTimes, int* NumObservationPoints, int* TotalNumData,
#         int* PumpingScheme, double* MassFlowrate, double* FlowDuration,
#         double* Time,
#         int* Deliverability, double* ProductionIndex, double* CutoffPressure, double* ModelledValue)

# def radial1d(double Porosity, double Permeability, double LayerThickness,
#         double ActionWellRadius, double RockSpecificHeat, double RockHeatConductivity, double RockDensity,
#         double RockCompressibility, double InitialPressure, double InitialX, int InjectionWell,
#         double InjectionEnthalpy, int NumPumpTimes, int NumObservationPoints, int TotalNumData,
#         int PumpingScheme, double MassFlowrate, double FlowDuration,
#         ndarray[double] Time,
#         int Deliverability, double ProductionIndex, double CutoffPressure):

   
#     cdef ndarray[double] ModelledValue = empty(TotalNumData)
#     print('Printing time from cython')
#     print(Time)
#     print('Finished printing time from cython')
#     print('Going into fortran function')
#     # c_radial1d(&Porosity, &Permeability, &LayerThickness,
#     #             &ActionWellRadius, &RockSpecificHeat, &RockHeatConductivity, &RockDensity,
#     #             &RockCompressibility, &InitialPressure, &InitialX, &InjectionWell,
#     #             &InjectionEnthalpy, &NumPumpTimes, &NumObservationPoints, &TotalNumData,
#     #             &PumpingScheme, &MassFlowrate, &FlowDuration, &PumpTime[0], &PumpRate[0],
#     #             &Time[0], &ObsPointRadialLocation[0], &ObsPointNumData[0], &ObsPointProperty[0],
#     #             &Deliverability, &ProductionIndex, &CutoffPressure, &ModelledValue[0])
#     c_radial1d(&Porosity, &Permeability, &LayerThickness,
#             &ActionWellRadius, &RockSpecificHeat, &RockHeatConductivity, &RockDensity,
#             &RockCompressibility, &InitialPressure, &InitialX, &InjectionWell,
#             &InjectionEnthalpy, &NumPumpTimes, &NumObservationPoints, &TotalNumData,
#             &PumpingScheme, &MassFlowrate, &FlowDuration,
#             &Time[0],
#             &Deliverability, &ProductionIndex, &CutoffPressure, &ModelledValue[0])

#     return ModelledValue

cdef extern from "radial1d_wrapper.h":
    extern void c_radial1d(double* Porosity, double* Permeability, double* LayerThickness,
        double* ActionWellRadius, double* RockSpecificHeat, double* RockHeatConductivity, double* RockDensity,
        double* RockCompressibility, double* InitialPressure, double* InitialX, int* InjectionWell,
        double* InjectionEnthalpy, int* NumPumpTimes, int* NumObservationPoints, int* TotalNumData,
        int* PumpingScheme, double* MassFlowrate, double* FlowDuration, double* PumpTime, double* PumpRate,
        double* Time, double* ObsPointRadialLocation, int* ObsPointNumData, int* ObsPointProperty,
        int* Deliverability, double* ProductionIndex, double* CutoffPressure, double* ModelledValue)

def radial1d(double Porosity, double Permeability, double LayerThickness,
        double ActionWellRadius, double RockSpecificHeat, double RockHeatConductivity, double RockDensity,
        double RockCompressibility, double InitialPressure, double InitialX, int InjectionWell,
        double InjectionEnthalpy, int NumPumpTimes, int NumObservationPoints, int TotalNumData,
        int PumpingScheme, double MassFlowrate, double FlowDuration, double[:] PumpTime, double[:] PumpRate,
        double[:] Time, double[:] ObsPointRadialLocation, int[:] ObsPointNumData, int[:] ObsPointProperty,
        int Deliverability, double ProductionIndex, double CutoffPressure):

   
    cdef ndarray[double] ModelledValue = empty(TotalNumData)
    print('Printing time from cython')
    print(Time)
    print('Finished printing time from cython')
    print('Going into fortran function')
    print('FlowRate = {}'.format(MassFlowrate))
    c_radial1d(&Porosity, &Permeability, &LayerThickness,
                &ActionWellRadius, &RockSpecificHeat, &RockHeatConductivity, &RockDensity,
                &RockCompressibility, &InitialPressure, &InitialX, &InjectionWell,
                &InjectionEnthalpy, &NumPumpTimes, &NumObservationPoints, &TotalNumData,
                &PumpingScheme, &MassFlowrate, &FlowDuration, &PumpTime[0], &PumpRate[0],
                &Time[0], &ObsPointRadialLocation[0], &ObsPointNumData[0], &ObsPointProperty[0],
                &Deliverability, &ProductionIndex, &CutoffPressure, &ModelledValue[0])

    return ModelledValue