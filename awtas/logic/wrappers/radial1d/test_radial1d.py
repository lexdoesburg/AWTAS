import numpy as np
import time
import radial1d_wrapper

print('inside')
print(radial1d_wrapper.__doc__)

# # double Porosity
# # double Permeability
# # double LayerThickness
# # double ActionWellRadius
# # double RockSpecificHeat
# # double RockHeatConductivity
# # double RockDensity
# # double RockCompressibility
# # double InitialPressure
# # double InitialX
# # bint InjectionWell
# # double InjectionEnthalpy
# # int NumPumpTimes
# # int NumObservationPoints
# # int TotalNumData
# # int PumpingScheme
# # double MassFlowrate
# # double FlowDuration
# # ndarray[double] PumpTime
# # ndarray[double] PumpRate
# # ndarray[double] Time
# # ndarray[double] ObsPointRadialLocation
# # ndarray[int] ObsPointNumData
# # ndarray[int] ObsPointProperty
# # bint Deliverability
# # double ProductionIndex
# # double CutoffPressure


time1, expected_pressure1 = np.genfromtxt('Pwell.dat', delimiter=',').T

time2, expected_pressure2 = np.genfromtxt('Pwell2.dat', delimiter=',').T
expected_pressure1 = expected_pressure1 * 1e5
expected_pressure2 = expected_pressure2 * 1e5

# time1 = np.linspace(0,86400,200)
# time2 = np.linspace(0,86400,200)
# time1 = np.array([1.00000001e-07,2.00000002e-07,4.00000005e-07,8.00000009e-07
# ,1.60000002e-06,3.20000004e-06,6.40000007e-06,1.28000001e-05
# ,2.56000003e-05,5.12000006e-05,1.01199999e-04,2.01200004e-04
# ,4.01199999e-04,8.01199989e-04,1.60119997e-03,3.20119993e-03
# ,6.40119985e-03,1.28012002e-02,2.56012008e-02,5.12011983e-02
# ,1.01201199e-01,2.01201200e-01,4.01201189e-01,8.01201224e-01
# ,1.60120118e+00,3.20120120e+00,6.40120125e+00,1.28012009e+01
# ,2.56012020e+01,5.12012024e+01,1.01201202e+02,2.01201202e+02
# ,4.01201202e+02,8.01201172e+02,1.60120117e+03,3.20120117e+03
# ,6.40120117e+03,1.28012012e+04,2.16000000e+04,3.91975977e+04
# ,4.23000000e+04,6.39000000e+04,8.64000000e+04])

# time2 = np.array([1.00000001e-07,2.00000002e-07,4.00000005e-07,8.00000009e-07
# ,1.60000002e-06,3.20000004e-06,6.40000007e-06,1.28000001e-05
# ,2.56000003e-05,5.12000006e-05,1.01199999e-04,2.01200004e-04
# ,4.01199999e-04,8.01199989e-04,1.60119997e-03,3.20119993e-03
# ,6.40119985e-03,1.28012002e-02,2.56012008e-02,5.12011983e-02
# ,1.01201199e-01,2.01201200e-01,4.01201189e-01,8.01201224e-01
# ,1.60120118e+00,3.20120120e+00,6.40120125e+00,1.28012009e+01
# ,2.56012020e+01,5.12012024e+01,1.01201202e+02,2.01201202e+02
# ,3.01201202e+02,5.01201202e+02,9.01201172e+02,1.30120117e+03
# ,2.10120117e+03,3.70120117e+03,5.30120117e+03,8.50120117e+03
# ,1.17012012e+04,1.81012012e+04,2.16000000e+04,3.44000000e+04
# ,4.23000000e+04,6.39000000e+04,8.64000000e+04])

time1 = time1.copy()
time2 = time2.copy()

initial_pressure = 40e5
initial_temp1 = 200.0
initial_sv2 = 0.9
const_flowrate1 = -8.0
const_flowrate2 = -2.0
total_data1 = len(time1)
total_data2 = len(time2)

phi = 0.1
k = 1e-14
layer_thickness = 100.0
well_radius = 0.0
rock_specific_heat = 1000.0
rock_heat_conductivity = 2.5
rock_density = 2500.0
rock_compressibility = 0.0
injection_well = 0
injection_enthalpy = 0.0
num_pump_times = 1
num_observation_points = 1
pumping_scheme = 1 # constant rate 1, measured 0
flow_duration = 0.0
pump_times1 = np.zeros(1, dtype=float)
pump_rates1 = np.zeros(1, dtype=float)
pump_times2 = np.zeros(1, dtype=float)
pump_rates2 = np.zeros(1, dtype=float)
obs_point_locations = np.ndarray(1, dtype=float)
obs_point_locations[0] = 0.01 - 1e-5 # There appears to be some loss of precision occuring due to types
obs_point_num_data1 = np.ndarray(1, dtype=np.int32)
obs_point_num_data1[0] = len(time1)
obs_point_num_data2 = np.ndarray(1, dtype=np.int32)
obs_point_num_data2[0] = len(time2)
obs_point_property = np.ones(1, dtype=np.int32)
obs_point_property[0] = 1
deliverability = 0
production_index = 0.0
cutoff_pressure = 0.0
pump_times1[0] = 0.0
pump_rates1[0] = -8
pump_times2[0] = 0.0
pump_rates2[0] = -2

# pump_times = np.array([0, 10000, 30000, 50000], dtype=float)
# pump_rates = np.array([-1, 0, -2, -0.5], dtype=float)

# num_blocks = 100
# num_constant_blocks = 20
# print(num_constant_blocks)
# constant_block_size = 0.01
# block_growth_factor = 1.2

print('Calling radial1d first time')
start_time1 = time.time()
pressure1, flag1 = radial1d_wrapper.radial1d(phi, k, layer_thickness, well_radius, rock_specific_heat, rock_heat_conductivity, rock_density, rock_compressibility, initial_pressure,
                        initial_temp1, injection_well, injection_enthalpy, num_pump_times, num_observation_points, total_data1, pumping_scheme,
                        pump_times1, pump_rates1, time1, obs_point_locations, obs_point_num_data1, obs_point_property, deliverability, production_index, cutoff_pressure)
                        # num_blocks, num_constant_blocks, constant_block_size, block_growth_factor)

end_time1 = time.time()
print('Calling radial1d second time')
start_time2 = time.time()
pressure2, flag2 = radial1d_wrapper.radial1d(phi, k, layer_thickness, well_radius, rock_specific_heat, rock_heat_conductivity, rock_density, rock_compressibility, initial_pressure,
                        initial_sv2, injection_well, injection_enthalpy, num_pump_times, num_observation_points, total_data2, pumping_scheme,
                        pump_times2, pump_rates2, time2, obs_point_locations, obs_point_num_data2, obs_point_property, deliverability, production_index, cutoff_pressure)
                        # num_blocks, num_constant_blocks, constant_block_size, block_growth_factor)
end_time2 = time.time()


# print('Calling radial1d first time')
# start_time1 = time.time()
# pressure1 = radial1d_wrapper.radial1d(phi, k, layer_thickness, well_radius, rock_specific_heat, rock_heat_conductivity, rock_density, rock_compressibility, initial_pressure,
#                         initial_temp1, injection_well, injection_enthalpy, num_pump_times, num_observation_points, total_data1, pumping_scheme, const_flowrate1, flow_duration,
#                         time1, deliverability, production_index, cutoff_pressure)
# end_time1 = time.time()
# print('Calling radial1d second time')
# pressure2 = radial1d_wrapper.radial1d(phi, k, layer_thickness, well_radius, rock_specific_heat, rock_heat_conductivity, rock_density, rock_compressibility, initial_pressure,
#                         initial_sv2, injection_well, injection_enthalpy, num_pump_times, num_observation_points, total_data2, pumping_scheme, const_flowrate2, flow_duration,
#                         time2, deliverability, production_index, cutoff_pressure)
# end_time2 = time.time()

print('Model 1 runtime = {}'.format(end_time1-start_time1))
print('Model 2 runtime = {}'.format(end_time2-start_time2))

print('Flag 1 = {} Flag 2 = {}'.format(flag1, flag2))


import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=1, ncols=2)
# Model 1
ax[0].plot(time1, expected_pressure1, 'k-', label='TOUGH2')
ax[0].plot(time1, pressure1, 'r--', label='AWTAS')
ax[0].legend(loc='best')
ax[0].set_title('Model 1')
ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('Pressure (Pa)')
# Model 2
ax[1].plot(time2, expected_pressure2, 'k-', label='TOUGH2')
ax[1].plot(time2, pressure2, 'r--', label='AWTAS')
ax[1].legend(loc='best')
ax[1].set_title('Model 2')
ax[1].set_xlabel('Time (s)')
ax[1].set_ylabel('Pressure (Pa)')

# time_eclipse1, pressure_eclipse1 = np.genfromtxt('Output_Model1.txt', delimiter=',').T
# time_eclipse2, pressure_eclipse2 = np.genfromtxt('Output_Model2.txt', delimiter=',').T
# ax[0].plot(time_eclipse1, pressure_eclipse1, 'g--', label='Eclipse')
# ax[1].plot(time_eclipse2, pressure_eclipse2, 'g--', label='Eclipse')
# ax[0].legend(loc='best')
# ax[1].legend(loc='best')

# # -O3 - with deliverability and recharge removed
# Model 1 runtime = 0.05554509162902832
# Model 2 runtime = 0.11745405197143555

# # Np ndarrays -O3
# Model 1 runtime = 0.08332204818725586
# Model 2 runtime = 0.19403386116027832

# # Np ndarrays -Ofast
# Model 1 runtime = 0.08932924270629883
# Model 2 runtime = 0.1617128849029541

# # Cython memory views -O3
# Model 1 runtime = 0.08358311653137207
# Model 2 runtime = 0.17027878761291504

# # Cython memory views -Ofast
# Model 1 runtime = 0.07443404197692871
# Model 2 runtime = 0.14401578903198242
plt.show()

