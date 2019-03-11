import numpy as np
import matplotlib.pyplot as plt
import time as time_module

import sys
import os
import inspect
current_dir = os.path.dirname(os.path.realpath(inspect.getfile(inspect.currentframe()))) # https://stackoverflow.com/questions/714063/importing-modules-from-parent-folder/11158224#11158224
awtas_main_dir = os.path.dirname(current_dir)
sys.path.insert(0, awtas_main_dir)

# To run without changing sys path use 'python -m tests.test_model' in command line while in main AWTAS directory

from awtas.logic import model
from awtas.logic import data as data_class

"""
This file was used to manually test the different forward models and inverse modelling. Excuse the mess.
NOTE: Can't confirm if all of this code will still run as there have been many updates to the model and data classes since last used.
"""

def setup_model(model_type, parameters, variables, time, test_num, noise=True, sd=1e-5, generate_datafile=False):
    model_type = model_type.lower()
    if model_type == 'theis':
        test_model = model.Theis_Solution()
    # elif model_type == 'theis fortran':
    #     test_model = model.Theis_Solution_Fortran()
    elif model_type == 'radial1d':
        test_model = model.Radial_1D()
    # elif model_type == 'test':
    #     test_model = model.Test_Model()
    test_model.generate_data(variables, parameters, time, noise=noise, sd=sd, save_file=generate_datafile, filename='{}_test{}.dat'.format(type(test_model).__name__, test_num))
    return test_model

def plot_solution(model, dual_plot=False):
    if dual_plot:
        f, ax = plt.subplots(nrows=1, ncols=2)
        f.suptitle('Observed Data vs Approximated Data for {} Model'.format(model.__class__.__name__))

        # Plot pressure vs log time graph
        ax[0].semilogx(np.log(model.data.time),model.data.observation,'kx',label="Observed Data")
        ax[0].semilogx(np.log(model.data.time),model.data.approximation,'r-',label="Approximated Data")
        ax[0].set_xlabel("Log Time (log(s))")
        ax[0].set_ylabel("Pressure (Pa)")
        ax[0].legend(loc="best")

        # Plot pressure vs time graph
        ax[1].plot(model.data.time,model.data.observation,'kx',label="Observed Data")
        ax[1].plot(model.data.time,model.data.approximation,'r-',label="Approximated Data")
        ax[1].set_xlabel("Time (s)")
        ax[1].set_ylabel("Pressure (Pa)")
    else:
        plt.semilogx(np.log(model.data.time), model.data.observation,"kx",label="Observed Data")
        plt.semilogx(np.log(model.data.time), model.data.approximation,"r-",label="Approximated Data")
        plt.title('Observed Data vs Approximated Data for {} Model'.format(model.__class__.__name__))
        plt.xlabel("Log Time (log(s))")
        plt.ylabel("Pressure (Pa)")
        plt.legend(loc="best")
    plt.show()

def print_optimal_parameters(model, variables, test_num):
    variable_names = ['porosity', 'permeability']
    optimal_parameters = model.find_model_parameters()
    for i in range(len(variables)):
        if i < len(variable_names):
            print('Test {} optimal {}: {}'.format(test_num, variable_names[i], optimal_parameters[i]))
        else:
            print('Test {} optimal parameter {}: {}'.format(test_num, i+1, optimal_parameters[i]))
        print('    Error: {}'.format(abs(variables[i]-optimal_parameters[i])))


# # -----------------------------------------------------------------------------------
# # Test Case 1: Parameters from page 11 AWTAS
# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.005 # m^3/s
# k = 1e-12 # m^2
# phi = 0.1
# rho = 813.37 # Water at 240 degrees celsius
# nu = 0.0001111 # Water at 240 degrees celsius
# C = 0.001303 # Water at 240 degrees celsius
# # t = 15 # days
# # t = t * 24 * 60 * 60
# t = 54000 # secs

# print(t)
# time = np.linspace(0, t, num=270)
# parameters = [p0, qm, h, rho, nu, C, r]
# test1_variables = [phi, k]

# # Build model
# # theis_test_1 = setup_model('theis', parameters, test1_variables, time, test_num=1, sd=1e-4)
# theis_test_1 = setup_model('theis', parameters, test1_variables, time, test_num=1, sd=150, noise=True, generate_datafile=True)
# theis_test_5 = setup_model('theis', parameters, test1_variables, time, test_num=5, sd=300, generate_datafile=True)
# theis_test_5.data.set_error(300)
# # -----------------------------------------------------------------------------------

# # -----------------------------------------------------------------------------------
# # Test 2 parameters
# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.015 # m^3/s
# k = 4e-12 # m^2
# phi = 0.01
# rho = 527.59 # Water at 360 degrees celsius
# nu = 0.0000603 # Water at 360 degrees celsius
# C = 0.03748 # Water at 360 degrees celsius
# t = 54000 # seconds
# time = np.linspace(0, t, num=100)
# parameters = [p0, qm, h, rho, nu, C, r]
# test2_variables = [phi, k]

# # Build model
# # theis_test_2 = setup_model('theis', parameters, test2_variables, time, test_num=2, sd=2e-5)
# theis_test_2 = setup_model('theis', parameters, test2_variables, time, test_num=2, sd=300, generate_datafile=True)

# # -----------------------------------------------------------------------------------

# # -----------------------------------------------------------------------------------
# # Test 3 parameters
# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.015 # m^3/s
# k = 1e-14 # m^2
# phi = 0.01
# rho = 527.59 # Water at 360 degrees celsius
# nu = 0.0000603 # Water at 360 degrees celsius
# C = 0.03748 # Water at 360 degrees celsius
# t = 54000 # seconds
# time = np.linspace(0, t, num=100)
# parameters = [p0, qm, h, rho, nu, C, r]
# test3_variables = [phi, k]

# # Build model
# # theis_test_3 = setup_model('theis', parameters, test3_variables, time, test_num=3, sd=1e-10)
# theis_test_3 = setup_model('theis', parameters, test3_variables, time, test_num=3, sd=0.0005, generate_datafile=True)

# # -----------------------------------------------------------------------------------

# # -----------------------------------------------------------------------------------
# # Test Case 1: Parameters from page 11 AWTAS
# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.005 # m^3/s
# k = 1e-12 # m^2
# phi = 0.1
# rho = 813.37 # Water at 240 degrees celsius
# nu = 0.0001111 # Water at 240 degrees celsius
# C = 0.001303 # Water at 240 degrees celsius
# t = 0.5 # days
# t = t * 24 * 60 * 60
# print(t)
# time = np.linspace(0, t, num=100)
# parameters = [p0, qm, h, rho, nu, C, r]
# test4_variables = [phi, k]

# # Build model
# # theis_test_4 = setup_model('theis', parameters, test4_variables, time, test_num=4, sd=3e-5)
# theis_test_4 = setup_model('theis', parameters, test4_variables, time, test_num=4, sd=90, generate_datafile=True)

# # -----------------------------------------------------------------------------------

# # -----------------------------------------------------------------------------------
# print_optimal_parameters(theis_test_1, test1_variables, 1)
# plot_solution(theis_test_1)

# print_optimal_parameters(theis_test_2, test2_variables, 2)
# plot_solution(theis_test_2)

# print_optimal_parameters(theis_test_3, test3_variables, 3)
# plot_solution(theis_test_3)

# print_optimal_parameters(theis_test_4, test4_variables, 3)
# plot_solution(theis_test_4)

# # ----------- Old find_model_parameters function-------------------------------------

# optimal_parameters = theis_test_1.find_model_parameters(verbose=True)
# plot_solution(theis_test_1, dual_plot=True)

# optimal_parameters = theis_test_2.find_model_parameters(verbose=True)
# plot_solution(theis_test_2, dual_plot=True)

# optimal_parameters = theis_test_3.find_model_parameters(verbose=True)
# plot_solution(theis_test_3, dual_plot=True)

# optimal_parameters = theis_test_4.find_model_parameters(verbose=True)
# plot_solution(theis_test_4, dual_plot=True)

# optimal_parameters = theis_test_5.find_model_parameters(verbose=True)
# plot_solution(theis_test_5, dual_plot=True)

# # ----------- New find_model_parameters function-------------------------------------

# optimal_parameters = theis_test_1.find_model_parameters2(verbose=True)
# plot_solution(theis_test_1, dual_plot=True)

# optimal_parameters = theis_test_2.find_model_parameters2(verbose=True)
# plot_solution(theis_test_2, dual_plot=True)

# optimal_parameters = theis_test_3.find_model_parameters2(verbose=True)
# plot_solution(theis_test_3, dual_plot=True)

# optimal_parameters = theis_test_4.find_model_parameters2(verbose=True)
# plot_solution(theis_test_4, dual_plot=True)

# optimal_parameters = theis_test_5.find_model_parameters2(verbose=True)
# plot_solution(theis_test_5, dual_plot=True)

# # -----------------------------------------------------------------------------------
# plt.plot(theis_test_1.data.time, theis_test_1.model([0.27, 6.95e-12]), 'k--', label='Initial')
# plt.plot(theis_test_1.data.time, theis_test_1.model(optimal_parameters), 'r-', label='Approximation')
# plt.plot(theis_test_1.data.time, theis_test_1.data.observation, 'kx', label='Data')
# plt.legend(loc='best')
# plt.show()

# # -----------------------------------------------------------------------------------

# time,pressure = np.genfromtxt('theis_test1.txt', delimiter=' ', skip_header=0).T
# print(time)
# print(pressure)
# # plt.semilogx(np.log(time),pressure,'kx',label='fortran')
# # plt.semilogx(np.log(theis_test_1.data.time), theis_test_1.data.observation, 'ro',label='python')
# plt.plot(time,pressure,'kx',label='fortran')
# plt.plot(theis_test_1.data.time, theis_test_1.data.observation, 'rx',label='python')
# plt.legend(loc='best')
# plt.show()


# # -----------------------------------------------------------------------------------
# # Test difficult function
# p1 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.005 # m^3/s
# k = 1e-12 # m^2
# phi = 0.1
# rho = 813.37 # Water at 240 degrees celsius
# nu = 0.0001111 # Water at 240 degrees celsius
# C = 0.001303 # Water at 240 degrees celsius
# parameters = [p1, qm, h, rho, nu, C, r]

# phi = 0.05
# k = 1.79e-13
# p0 = 3.1256e6
# x0 = 0.17
# variables = [phi, k, p0, x0]

# t = 54000 # seconds
# time = np.linspace(0, t, num=100)
# # data = data_class.Data()
# # data.set_time(time)

# # Build the model and data
# test_model = model.Test_Model()
# test_model.generate_data(variables, parameters, time,noise=True,sd=5e-8)
# parameters = test_model.find_model_parameters()
# print(parameters)
# plot_solution(test_model, dual_plot=True)

# # # ----------------------------------------------------------------------------------
# # Test Radial1d Model

# k=1e-13
# phi=0.1
# p0=40e5
# X0=200
# rw=0.1
# thick=100
# CR=1000
# COND=2.5
# RHOR=2500
# COMP=0
# ConstRate=5
# distFromWell=0.1 # distance of observation point from action well
# t = 54000
# num_data = 271

# time = np.linspace(0, t, num=num_data)
# parameters = [p0, X0, rw, thick, CR, COND, RHOR, COMP, ConstRate, distFromWell]
# radial1_variables = [phi, k]

# # # Build model
# # # theis_test_4 = setup_model('theis', parameters, test4_variables, time, test_num=4, sd=3e-5)
# # radial_test_1 = setup_model('radial1d', parameters, radial1_variables, time, test_num=1, sd=15000,generate_datafile=True)

# time,pressure = np.genfromtxt('Radial_1D_test1.dat', delimiter=',', skip_header=6).T

# radial_test_1_data = data_class.Data(model_type='radial1d')
# radial_test_1_data.observation = pressure
# radial_test_1_data.error = 15000
# radial_test_1_data.set_known_parameters(parameters)
# radial_test_1_data.set_time(time)
# radial_test_1 = model.Radial_1D(radial_test_1_data)

# # radial_test_1.data.observation = pressure
# # radial_test_1.data.time = time
# # radial_test_1.data.parameters = parameters
# # radial_test_1.data.error = 15000

# # print(radial_test_1.data.parameters)
# # print(radial_test_1.data.time)

# # # # pressure = radial_test_1.model([0.1, 1e-13])
# # # from radial1d_wrapper import radial1d
# # # pressure2 = radial1d(0.09293594872248502, 1.0036337151556582e-13, p0, X0, rw, thick, CR, COND, RHOR, COMP, ConstRate, distFromWell, num_data, np.linspace(0, t, num=num_data))
# # # # print(time)
# # # # print(pressure)
# # # plt.plot(time,pressure,'kx')
# # # plt.plot(time,pressure2,'r-')
# # # plt.show()

# # print('Start')
# start = time_module.time()
# optimal_parameters = radial_test_1.find_model_parameters2(verbose=True)
# end = time_module.time()
# print(radial_test_1.data.approximation)

# print('Optimal phi = {} Optimal k = {}'.format(optimal_parameters[0],optimal_parameters[1]))
# # # plt.plot(radial_test_1.data.time, radial_test_1.data.approximation)
# plot_solution(radial_test_1, dual_plot=True)
# print('Time to find solution = ', end-start)
# # radial_data = data_class.create_data('radial1d', filename=None, time=None, observation=None, parameters=parameters, error=None)
# # print(radial_data.parameters)
# # radial_data.generate_datafile('testing_radial1d_data.txt')

# # new_radial = data_class.create_data('radial1d')
# # new_radial.read_file('testing_radial1d_data.txt')
# # print(new_radial.parameters)

# # # -----------------------------------------------------------------------------------



# # plt.plot(t, data,"k-",label="Synthetic Data (Theis Solution W/Noise)")
# start = time_module.time()
# plt.semilogx(np.log(time), theis_model.data.observation,"kx",label="Synthetic Data (Theis Solution W/Noise)")
# plt.semilogx(np.log(time), theis_model.data.approximation,"r-",label="Approximated Curve")
# # plt.plot(t, p, "g-", label="Theis Analytic Solution")
# plt.title("Observed Data vs Fitted Curve")
# plt.xlabel("Log Time (s)")
# plt.ylabel("Pressure (Pa)")
# plt.legend(loc="best")
# end = time_module.time()
# print('Time elapsed = {}'.format(end - start))
# plt.show()

# # Test 2 parameters
# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.015 # m^3/s
# k = 4e-12 # m^2
# phi = 0.01
# rho = 527.59 # Water at 360 degrees celsius
# nu = 0.0000603 # Water at 360 degrees celsius
# C = 0.03748 # Water at 360 degrees celsius
# t = 54000 # seconds
# time = np.linspace(0, t, num=100)
# parameters = [p0, qm, h, rho, nu, C, r]

# variables = [phi, k]
# # Build the model and data
# theis_model = model.Theis_Solution()
# theis_model.generate_data(variables, parameters, time, noise=True, sd=2e-5, save_file=True, filename="theis_test2_datafile_v1_1.dat")
# # theis_model.data.generate_datafile("theis_test1_datafile2.txt")
# opt_phi, opt_k = theis_model.find_model_parameters()
# print(opt_phi)
# print(opt_k)

# # plt.plot(t, data,"k-",label="Synthetic Data (Theis Solution W/Noise)")
# start = time_module.time()
# plt.semilogx(np.log(time), theis_model.data.observation,"kx",label="Synthetic Data (Theis Solution W/Noise)")
# plt.semilogx(np.log(time), theis_model.data.approximation,"r-",label="Approximated Curve")
# # plt.plot(t, p, "g-", label="Theis Analytic Solution")
# plt.title("Observed Data vs Fitted Curve")
# plt.xlabel("Log Time (s)")
# plt.ylabel("Pressure (Pa)")
# plt.legend(loc="best")
# end = time_module.time()
# print('Time elapsed = {}'.format(end - start))
# plt.show()


# # Test 3 parameters
# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.015 # m^3/s
# k = 1e-14 # m^2
# phi = 0.01
# rho = 527.59 # Water at 360 degrees celsius
# nu = 0.0000603 # Water at 360 degrees celsius
# C = 0.03748 # Water at 360 degrees celsius
# t = 54000 # seconds
# time = np.linspace(0, t, num=100)
# parameters = [p0, qm, h, rho, nu, C, r]

# variables = [phi, k]
# # Build the model and data
# theis_model = model.Theis_Solution()
# theis_model.generate_data(variables, parameters, time, noise=True, sd=1e-10, save_file=True, filename="theis_test3_datafile_v1_1.dat")
# # theis_model.data.generate_datafile("theis_test1_datafile2.txt")
# opt_phi, opt_k = theis_model.find_model_parameters()
# print(opt_phi)
# print(opt_k)

# # plt.plot(t, data,"k-",label="Synthetic Data (Theis Solution W/Noise)")
# start = time_module.time()
# plt.semilogx(np.log(time), theis_model.data.observation,"kx",label="Synthetic Data (Theis Solution W/Noise)")
# plt.semilogx(np.log(time), theis_model.data.approximation,"r-",label="Approximated Curve")
# # plt.plot(t, p, "g-", label="Theis Analytic Solution")
# plt.title("Observed Data vs Fitted Curve")
# plt.xlabel("Log Time (s)")
# plt.ylabel("Pressure (Pa)")
# plt.legend(loc="best")
# end = time_module.time()
# print('Time elapsed = {}'.format(end - start))
# plt.show()



# if __name__ == '__main__':
#     plot_solution(theis_test_1)

# # -----------------------------------------------------------------------------------
# # Try two-phase production SKG9D case
# # -----------------------------------------------------------------------------------
# # from awtas.logic.wrappers.radial1d import radial1d_wrapper

# # time_days, expected_pressure = np.genfromtxt('SKG9D_press.dat', delimiter=', ').T 
# # time_days, expected_enthalpy = np.genfromtxt('SKG9D_enth.dat', delimiter=', ').T
# # time_secs_ent, actual_enthalpy = np.genfromtxt('two_phase_ent.txt', delimiter=',').T
# # actual_enthalpy = actual_enthalpy/1000
# # time_secs = time_days * 86400
# # # expected_pressure = expected_pressure * 1e5

# # flow_times_secs, flow_rates = np.genfromtxt('SKG9D_flow.dat', delimiter=', ').T 
# # flow_times_days = flow_times_secs / 86400

# # time_secs1, modelled_value = np.genfromtxt('two_phase.txt', delimiter=',').T
# # time_days1 = time_secs1/86400

# # k = 7.27e-16
# # phi = 0.091

# # initial_pressure = 138.86733939886042 * 1e5
# # # initial_temp = 336
# # initial_sv = 0.1
# # layer_thickness = 200
# # rock_specific_heat = 1000
# # rock_heat_conductivity = 2.5
# # rock_density = 2500
# # rock_compressibility = 0.0
# # well_radius = 0.2
# # injection_well = 0
# # injection_enthalpy = 0.0
# # # time = time_secs[:90]
# # time = time_secs

# # total_data = len(time)
# # num_observation_points = 1
# # num_pump_times = len(flow_rates)
# # pumping_scheme = 0
# # mass_flowrate = 0
# # flow_duration = 0
# # pump_times = flow_times_secs.copy()
# # pump_rates = flow_rates.copy()
# # obs_point_locations = np.array([well_radius], dtype=float)
# # obs_point_num_data = np.array([total_data], dtype=np.int32)
# # obs_point_property = np.array([1], dtype=np.int32)
# # obs_point_property2 = np.array([3], dtype=np.int32)
# # deliverability = 0
# # production_index = 0
# # cutoff_pressure = 0
# # num_blocks = 100
# # num_constant_blocks = 10
# # constant_block_size = 0.1
# # block_growth_factor = 1.2

# # time_days2 = time/86400
# # from time import time as timer
# # start1 = timer()
# # modelled_pressure, flag1 = radial1d_wrapper.radial1d(phi, k, layer_thickness, well_radius, rock_specific_heat, rock_heat_conductivity, rock_density, rock_compressibility, initial_pressure,
# #                                 initial_sv, injection_well, injection_enthalpy, num_pump_times, num_observation_points, total_data, pumping_scheme, mass_flowrate, flow_duration,
# #                                 pump_times, pump_rates, time, obs_point_locations, obs_point_num_data, obs_point_property, deliverability, production_index, cutoff_pressure,
# #                                 num_blocks, num_constant_blocks, constant_block_size, block_growth_factor)
# # end1 = timer()

# # start2 = timer()
# # modelled_enthalpy, flag2 = radial1d_wrapper.radial1d(phi, k, layer_thickness, well_radius, rock_specific_heat, rock_heat_conductivity, rock_density, rock_compressibility, initial_pressure,
# #                                 initial_sv, injection_well, injection_enthalpy, num_pump_times, num_observation_points, total_data, pumping_scheme, mass_flowrate, flow_duration,
# #                                 pump_times, pump_rates, time, obs_point_locations, obs_point_num_data, obs_point_property2, deliverability, production_index, cutoff_pressure,
# #                                 num_blocks, num_constant_blocks, constant_block_size, block_growth_factor)
# # end2= timer()
# # print('Time elapsed 1 = {}'.format(end1-start1))
# # print('Time elapsed 2 = {}'.format(end2-start2))

# # f, ax1 = plt.subplots(1,2)
# # # ax2 = ax1[0].twinx()
# # # ax3 = ax1[1].twinx()
# # ax1[0].plot(time_days, expected_pressure, 'k-', label='Measured pressure')
# # ax1[0].plot(time_days2, modelled_pressure/1e5, 'r--', label='Simulated pressure')
# # ax1[0].legend(loc='best')
# # # ax1[0].plot(time_days1, modelled_value/1e5, 'r--') # from txt file

# # # ax2.plot(flow_times_days, flow_rates, 'b-')
# # measured_enth = ax1[1].plot(time_days, expected_enthalpy, 'k-', label='Measured enthalpy')
# # # ax1[1].plot(time_days1, actual_enthalpy, 'r--')
# # sim_enth = ax1[1].plot(time_days2, modelled_enthalpy/1000, 'r--', label='Simulated enthalpy')
# # # flows = ax3.plot(flow_times_days, flow_rates, 'b--', label='Flow rates')
# # # lines = measured_enth + sim_enth + flows
# # # labels = [l.get_label() for l in lines]
# # # ax1[1].legend(lines, labels, loc='best')
# # ax1[1].legend(loc='best')
# # # ax1[0].set_xlim([0,97])
# # # ax1[1].set_xlim([0,97])

# # # ax2.set_xlim([0,97])
# # ax1[0].set_ylim([-5,130])
# # ax1[1].set_ylim([1500,2400])
# # # ax2.set_ylim([-6,0])
# # # ax3.set_ylim([-6,0])

# # plt.show()

# #  Overall radial1d function call time:    10.134601000000000     
# #  Overall radial1d_wrapper (inside) time:    10.134609000000001     
# # Radial1d wrapper time (from cython): 10.221587896347046
# #  Inside homogeneous porous
# #  Execution:            2
# #  Time elapsed inside homogeneous porous = 
# #    10.081318999999999     
# #  Overall radial1d function call time:    10.081356000000001     
# #  Overall radial1d_wrapper (inside) time:    10.081359000000001     
# # Radial1d wrapper time (from cython): 10.172600746154785
# # Time elapsed 1 = 10.222235918045044
# # Time elapsed 2 = 10.17296290397644


# time_days_pres, expected_pressure = np.genfromtxt('SKG9D_press.dat', delimiter=', ').T 
# time_days_ent, expected_enthalpy = np.genfromtxt('SKG9D_enth.dat', delimiter=', ').T
# time_secs_ent, actual_enthalpy = np.genfromtxt('two_phase_ent.txt', delimiter=',').T
# expected_enthalpy = expected_enthalpy * 1000
# time_secs_pres = time_days_pres * 86400
# time_secs_ent = time_days_ent * 86400
# time_pressure = time_secs_pres.copy()
# time_enthalpy = time_secs_ent.copy()

# expected_pressure = expected_pressure * 1e5

# flow_times_secs, flow_rates = np.genfromtxt('SKG9D_flow.dat', delimiter=', ').T 
# flow_times_days = flow_times_secs / 86400

# time = time_secs_pres.copy()
# pump_rates = flow_rates.copy()
# pump_times = flow_times_secs.copy()

# k = 7.27e-16
# phi = 0.091

# initial_pressure = 138.86733939886042 * 1e5
# initial_x = 0.1
# layer_thickness = 200
# rock_specific_heat = 1000
# rock_heat_conductivity = 2.5
# rock_density = 2500
# rock_compressibility = 0.0
# action_well_radius = 0.2

# parameters = [initial_pressure, initial_x, action_well_radius, layer_thickness, rock_specific_heat, rock_heat_conductivity, rock_density, rock_compressibility]
# variables = [phi, k]

# # Set up observation points
# radial_location = [0.099, 0.099]
# property = ['Pressure', 'Enthalpy']
# num_data = [len(time_pressure), len(time_enthalpy)]
# times = [time_pressure, time_enthalpy]
# observations = [expected_pressure, expected_enthalpy]
# observation_points = data_class.ObservationPoints(radial_location, property, num_data, times, observations)

# # # observation_points = data_class.ObservationPoints(radial_location=0.099, property='pressure', num_data=len(time), times=time, observations=expected_pressure)
# # print(observation_points.num_data)
# # print(observation_points.num_observation_points)
# # print(observation_points.property)
# # print(observation_points.radial_location)
# # print(observation_points.times)
# # print(observation_points.observations)

# # time = time[:75]

# pump_info = data_class.Pump(pumping_scheme='Measured Flows', flow_rates=pump_rates, flow_times=pump_times, deliverability=False, production_index=0.0, cutoff_pressure=0.0)
# grid_info = data_class.Grid(num_blocks=100, num_constant_blocks=10, constant_block_size=0.1, block_growth_factor=1.2)

# data = data_class.Data('radial1d', parameters_list=parameters, pump_info=pump_info, observation_points=observation_points, grid_info=grid_info)

# # print(data.time)
# # print(data.observation)

# new_data = data_class.Data('radial1d', 'radial1d_testcase1.txt')

# # radial1d_model = model.Radial_1D(data=data)
# radial1d_model = model.Radial_1D(data=new_data)



# # optimal_parameters = radial1d_model.find_model_parameters2(verbose=True, initial_guess=[0.091, 8e-16], single_run=False)
# # # print('Optimal params = ', optimal_parameters)
# # # optimal_parameters = [0.04580637995823384, 8.473963146776042e-16]
# optimal_parameters = [phi, k]
# # radial1d_model.data.time = time_secs.copy()
# modelled_value = radial1d_model.model(optimal_parameters)

# radial1d_model.data.observation_points.store_modelled_values(modelled_value)
# print(radial1d_model.data.observation_points.modelled_values)
# print(len(radial1d_model.data.observation_points.modelled_values[0]))
# print(len(radial1d_model.data.observation_points.modelled_values[1]))

# plt.plot(radial1d_model.data.observation_points.times[0]/86400, radial1d_model.data.observation_points.observations[0]/1e5, 'k-', label='Observed Pressure')
# plt.plot(radial1d_model.data.observation_points.times[0]/86400, radial1d_model.data.observation_points.modelled_values[0]/1e5, 'r--', label='Modelled Pressure')
# plt.plot(radial1d_model.data.observation_points.times[1]/86400, radial1d_model.data.observation_points.observations[1]/1e5, 'g-', label='Observed Enthalpy')
# plt.plot(radial1d_model.data.observation_points.times[1]/86400, radial1d_model.data.observation_points.modelled_values[1]/1e5, 'b--', label='Modelled Enthalpy')
# # plt.plot(data.time/86400, modelled_value/1e5, 'r--', label='Modelled value')
# plt.legend(loc='best')
# plt.show()
# # xd = []
# # print(type(xd))
# # print(type(data.observation_points.num_data))

# # data.generate_datafile('radial1d_testcase1.txt', variables={'Porosity' : phi, 'Permeability' : k})


# # print('Fixed parameters')
# # print(new_data.fixed_parameters)
# # print(new_data.reservoir_conditions)

# # print('Grid info')
# # print(new_data.grid_info.num_blocks)
# # print(new_data.grid_info.num_constant_blocks)
# # print(new_data.grid_info.constant_block_size)
# # print(new_data.grid_info.block_growth_factor)

# # print('Pump info')
# # print(new_data.pump_info.num_pump_times)
# # print(new_data.pump_info.injection_well)
# # print(new_data.pump_info.injection_enthalpy)
# # print(new_data.pump_info.deliverability)
# # print(new_data.pump_info.production_index)
# # print(new_data.pump_info.cutoff_pressure)
# # print(new_data.pump_info.flow_times)
# # print(new_data.pump_info.flow_rates)

# # print('Observation point info')
# # print(new_data.observation_points.num_observation_points)
# # print(new_data.observation_points.radial_location)
# # print(new_data.observation_points.num_data)
# # print(new_data.observation_points.property)
# # print(new_data.observation_points.times)
# # print(new_data.observation_points.observations)

