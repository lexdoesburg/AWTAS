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


def setup_model(model_type, parameters, variables, time, test_num, noise=True, sd=1e-5, generate_datafile=False):
    model_type = model_type.lower()
    if model_type == 'theis':
        test_model = model.Theis_Solution()
    elif model_type == 'theis fortran':
        test_model = model.Theis_Solution_Fortran()
    elif model_type == 'radial1d':
        test_model = model.Radial_1D()
    elif model_type == 'test':
        test_model = model.Test_Model()
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


# -----------------------------------------------------------------------------------
# Test Case 1: Parameters from page 11 AWTAS
p0 = 3.6e6 # Pa
h = 100 # m
r = 0.05 # m
qm = -0.005 # m^3/s
k = 1e-12 # m^2
phi = 0.1
rho = 813.37 # Water at 240 degrees celsius
nu = 0.0001111 # Water at 240 degrees celsius
C = 0.001303 # Water at 240 degrees celsius
# t = 15 # days
# t = t * 24 * 60 * 60
t = 54000 # secs

print(t)
time = np.linspace(0, t, num=270)
parameters = [p0, qm, h, rho, nu, C, r]
test1_variables = [phi, k]

# Build model
# theis_test_1 = setup_model('theis', parameters, test1_variables, time, test_num=1, sd=1e-4)
theis_test_1 = setup_model('theis', parameters, test1_variables, time, test_num=1, sd=150, noise=True)
theis_test_5 = setup_model('theis', parameters, test1_variables, time, test_num=1, sd=300)
theis_test_5.data.set_error(300)
# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
# Test 2 parameters
p0 = 3.6e6 # Pa
h = 100 # m
r = 0.05 # m
qm = -0.015 # m^3/s
k = 4e-12 # m^2
phi = 0.01
rho = 527.59 # Water at 360 degrees celsius
nu = 0.0000603 # Water at 360 degrees celsius
C = 0.03748 # Water at 360 degrees celsius
t = 54000 # seconds
time = np.linspace(0, t, num=100)
parameters = [p0, qm, h, rho, nu, C, r]
test2_variables = [phi, k]

# Build model
# theis_test_2 = setup_model('theis', parameters, test2_variables, time, test_num=2, sd=2e-5)
theis_test_2 = setup_model('theis', parameters, test2_variables, time, test_num=2, sd=300)

# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
# Test 3 parameters
p0 = 3.6e6 # Pa
h = 100 # m
r = 0.05 # m
qm = -0.015 # m^3/s
k = 1e-14 # m^2
phi = 0.01
rho = 527.59 # Water at 360 degrees celsius
nu = 0.0000603 # Water at 360 degrees celsius
C = 0.03748 # Water at 360 degrees celsius
t = 54000 # seconds
time = np.linspace(0, t, num=100)
parameters = [p0, qm, h, rho, nu, C, r]
test3_variables = [phi, k]

# Build model
# theis_test_3 = setup_model('theis', parameters, test3_variables, time, test_num=3, sd=1e-10)
theis_test_3 = setup_model('theis', parameters, test3_variables, time, test_num=3, sd=0.0005)

# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
# Test Case 1: Parameters from page 11 AWTAS
p0 = 3.6e6 # Pa
h = 100 # m
r = 0.05 # m
qm = -0.005 # m^3/s
k = 1e-12 # m^2
phi = 0.1
rho = 813.37 # Water at 240 degrees celsius
nu = 0.0001111 # Water at 240 degrees celsius
C = 0.001303 # Water at 240 degrees celsius
t = 0.5 # days
t = t * 24 * 60 * 60
print(t)
time = np.linspace(0, t, num=100)
parameters = [p0, qm, h, rho, nu, C, r]
test4_variables = [phi, k]

# Build model
# theis_test_4 = setup_model('theis', parameters, test4_variables, time, test_num=4, sd=3e-5)
theis_test_4 = setup_model('theis', parameters, test4_variables, time, test_num=4, sd=90)

# -----------------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------------
# Try two-phase production SKG9D case
# -----------------------------------------------------------------------------------
from awtas.logic.wrappers.radial1d import radial1d_wrapper

time_days, expected_pressure = np.genfromtxt('SKG9D_press.dat', delimiter=', ').T 
time_days, expected_enthalpy = np.genfromtxt('SKG9D_enth.dat', delimiter=', ').T
time_secs_ent, actual_enthalpy = np.genfromtxt('two_phase_ent.txt', delimiter=',').T
actual_enthalpy = actual_enthalpy/1000
time_secs = time_days * 86400
# # print(time_secs)
# time_secs = np.array([0.,86400.,172800.,259200.,345600.,432000.,518400.
# ,604800.,691200.,777600.,864000.,950400.,1036800.,1123200.
# ,1209600.,1296000.,1382400.,1468800.,1555200.,1641600.,1728000.
# ,1814400.,1900800.,1987200.,2073600.,2160000.,2246400.,2332800.
# ,2419200.,2505600.,2592000.,2678400.,2764800.,2851200.,2937600.
# ,3024000.,3110400.,3196800.,3283200.,3369600.,3456000.,3542400.
# ,3628800.,3715200.,3801600.,3888000.,3974400.,4060800.,4147200.
# ,4492800.,4579200.,4665600.,4752000.,4924800.,5097600.,5184000.
# ,5270400.,5356800.,5443200.,5529600.,5616000.,5702400.,5788800.
# ,5875200.,5961600.,6048000.,6134400.,6220800.,6307200.,6393600.
# ,6480000.,6566400.,6652800.,6739200.,6825600.,6912000.,6998400.
# ,7084800.,7171200.,7257600.,7344000.,7516800.,7603200.,7689600.
# ,7776000.,7862400.,7948800.,8121600.,8208000.,8294400.,8380800.
# ,9158400.,9244800.,9417600.,9504000.,9590400.,9763200.,9849600.
# ,9936000.,10022400.,10108800.,10195200.,10281600.,10368000.,10454400.
# ,10540800.,10627200.,10713600.,10800000.,10886400.,10972800.,11059200.
# ,11145600.,11232000.,11318400.,11404800.,11491200.,11577600.,11664000.
# ,11750400.,11836800.,11923200.,12009600.])
# expected_pressure = expected_pressure * 1e5

flow_times_secs, flow_rates = np.genfromtxt('SKG9D_flow.dat', delimiter=', ').T 
# flow_times_secs = np.array([0.,43200.,86400.,130000.,173000.,216000.,259000.
# ,302000.,346000.,389000.,432000.,475000.,518000.,562000.
# ,605000.,648000.,691000.,734000.,778000.,821000.,864000.
# ,907000.,950000.,994000.,1040000.,1080000.,1120000.,1170000.
# ,1210000.,1250000.,1300000.,1340000.,1380000.,1430000.,1470000.
# ,1510000.,1560000.,1600000.,1640000.,1680000.,1730000.,1770000.
# ,1810000.,1860000.,1900000.,1940000.,1990000.,2030000.,2070000.
# ,2120000.,2160000.,2200000.,2250000.,2290000.,2330000.,2380000.
# ,2420000.,2460000.,2510000.,2550000.,2590000.,2640000.,2680000.
# ,2720000.,2760000.,2810000.,2850000.,2890000.,2940000.,2980000.
# ,3020000.,3070000.,3110000.,3150000.,3200000.,3240000.,3280000.
# ,3330000.,3370000.,3410000.,3460000.,3500000.,3540000.,3590000.
# ,3630000.,3670000.,3720000.,3760000.,3800000.,3840000.,3890000.
# ,3930000.,3970000.,4020000.,4060000.,4100000.,4150000.,4190000.
# ,4230000.,4280000.,4320000.,4360000.,4410000.,4450000.,4490000.
# ,4540000.,4580000.,4620000.,4670000.,4710000.,4750000.,4800000.
# ,4840000.,4880000.,4920000.,4970000.,5010000.,5050000.,5100000.
# ,5140000.,5180000.,5230000.,5270000.,5310000.,5360000.,5400000.
# ,5440000.,5490000.,5530000.,5570000.,5620000.,5660000.,5700000.
# ,5750000.,5790000.,5830000.,5880000.,5920000.,5960000.,6000000.
# ,6050000.,6090000.,6130000.,6180000.,6220000.,6260000.,6310000.
# ,6350000.,6390000.,6440000.,6480000.,6520000.,6570000.,6610000.
# ,6650000.,6700000.,6740000.,6780000.,6830000.,6870000.,6910000.
# ,6960000.,7000000.,7040000.,7080000.,7130000.,7170000.,7210000.
# ,7260000.,7300000.,7340000.,7390000.,7430000.,7470000.,7520000.
# ,7560000.,7600000.,7650000.,7690000.,7730000.,7780000.,7820000.
# ,7860000.,7910000.,7950000.,7990000.,8040000.,8080000.,8120000.
# ,8160000.,8210000.,8250000.,8290000.,8340000.,8380000.,8420000.
# ,8470000.,8510000.,8550000.,8600000.,8640000.,8680000.,8730000.
# ,8770000.,8810000.,8860000.,8900000.,8940000.,8990000.,9030000.
# ,9070000.,9120000.,9160000.,9200000.,9240000.,9290000.,9330000.
# ,9370000.,9420000.,9460000.,9500000.,9550000.,9590000.,9630000.
# ,9680000.,9720000.,9760000.,9810000.,9850000.,9890000.,9940000.
# ,9980000.,10000000.,10100000.,10100000.,10200000.,10200000.,10200000.
# ,10300000.,10300000.,10400000.,10400000.,10500000.,10500000.,10500000.
# ,10600000.,10600000.,10700000.,10700000.,10800000.,10800000.,10800000.
# ,10900000.,10900000.,11000000.,11000000.,11100000.,11100000.,11100000.
# ,11200000.,11200000.,11300000.,11300000.,11400000.,11400000.,11400000.
# ,11500000.,11500000.,11600000.,11600000.,11700000.,11700000.,11800000.
# ,11800000.,11800000.,11900000.,11900000.,12000000.,12000000.])

# flow_rates = np.array([-5.82,-5.14,-5.19,-5.25,-4.69,-4.14,-3.89,-3.64,-3.57
# ,-3.5 ,-3.64,-3.78,-3.72,-3.67,-3.44,-3.22,-3.25,-3.28
# ,-3.71,-4.14,-3.85,-3.56,-3.51,-3.47,-3.58,-3.69,-3.56
# ,-3.42,-3.76,-4.11,-3.78,-3.44,-3.61,-3.78,-3.53,-3.28
# ,-3.31,-3.33,-3.29,-3.25,-3.31,-3.36,-3.29,-3.22,-3.26
# ,-3.31,-3.53,-3.75,-3.46,-3.17,-3.29,-3.42,-3.42,-3.42
# ,-3.47,-3.53,-3.32,-3.11,-3.11,-3.11,-3.39,-3.67,-3.15
# ,-2.64,-2.75,-2.86,-2.83,-2.81,-2.81,-2.81,-3.1, -3.39
# ,-3.18,-2.97,-3.07,-3.17,-3.13,-3.08,-3.04,-3.,-2.96
# ,-2.92,-2.96,-3.,-3.07,-3.14,-3.1, -3.06,-2.99,-2.92
# ,-2.86,-2.81,-3.01,-3.22,-3.13,-3.03,-1.51, 0.,0.
# ,0.,0.,0.,-1.99,-3.97,-3.88,-3.78,-3.88,-3.97
# ,-3.78,-3.58,-1.79, 0.,-2.15,-4.31,-4.22,-4.14,-4.08
# ,-4.03,-4.04,-4.06,-4.07,-4.08,-4.17,-4.25,-4.35,-4.44
# ,-4.53,-4.61,-4.5, -4.39,-4.4, -4.42,-4.32,-4.22,-4.28
# ,-4.33,-4.28,-4.22,-4.39,-4.56,-4.5, -4.44,-4.5, -4.56
# ,-4.36,-4.17,-4.39,-4.61,-4.5, -4.39,-4.28,-4.17,-4.17
# ,-4.17,-4.08,-4.,-4.06,-4.11,-3.88,-3.64,-3.61,-3.58
# ,-3.15,-2.72,-2.93,-3.14,-3.04,-2.94,-2.72,-2.5, -2.51
# ,-2.53,-2.51,-2.5, -2.17,-1.83,-1.74,-1.64,-1.67,-1.69
# ,-1.47,-1.25,-1.54,-1.83,-1.78,-1.72,-1.69,-1.67,-1.81
# ,-1.94,-1.9, -1.86,-1.83,-1.81,-0.903,0.,0.,0.
# ,0.,0.,0.,0.,0.,0.,0.,0.,0.
# ,0.,0.,0.,-3.68,-7.36,-9.44,-11.5,-11.5,-11.6
# ,-11.5,-11.5,-11.3,-11.2,-11.1,-10.9,-10.7,-10.6,-10.4
# ,-10.2,-10.3,-10.4,-10.3,-10.1,-10.,-9.86,-9.96,-10.1
# ,-10.4,-10.7, -9.79,-8.89,-8.85,-8.81,-8.78,-8.75,-8.89
# ,-9.03,-9.03,-9.03,-9.24,-9.44,-9.51,-9.58,-9.33,-9.08
# ,-9.15,-9.22,-9.33,-9.44,-9.38,-9.31,-9.18,-9.06,-9.04
# ,-9.03,-8.96,-8.89,-9.1, -9.31,-9.1, -8.89,-8.99,-9.08
# ,-9.15,-9.22,-9.26,-9.31,-9.38,-9.44,-9.24,-9.03,-9.03])

# print(flow_times_secs)
# print(flow_rates)
flow_times_days = flow_times_secs / 86400

time_secs1, modelled_value = np.genfromtxt('two_phase.txt', delimiter=',').T
time_days1 = time_secs1/86400

k = 7.27e-16
phi = 0.091


# from copy import copy
initial_pressure = 138.86733939886042 * 1e5
# initial_temp = 336
initial_sv = 0.1
layer_thickness = 200
rock_specific_heat = 1000
rock_heat_conductivity = 2.5
rock_density = 2500
rock_compressibility = 0.0
well_radius = 0.0
injection_well = 0
injection_enthalpy = 0.0
# time = time_secs[:90]
time = time_secs

total_data = len(time)
num_observation_points = 1
num_pump_times = len(flow_rates)
pumping_scheme = 0
mass_flowrate = 0
flow_duration = 0
pump_times = flow_times_secs.copy()
pump_rates = flow_rates.copy()
obs_point_locations = np.array([0.01], dtype=float)
obs_point_num_data = np.array([total_data], dtype=np.int32)
obs_point_property = np.array([1], dtype=np.int32)
obs_point_property2 = np.array([3], dtype=np.int32)
deliverability = 0
production_index = 0
cutoff_pressure = 0


print(time)
time_days2 = time/86400
from time import time as timer
start1 = timer()
modelled_pressure = radial1d_wrapper.radial1d(phi, k, layer_thickness, well_radius, rock_specific_heat, rock_heat_conductivity, rock_density, rock_compressibility, initial_pressure,
                                initial_sv, injection_well, injection_enthalpy, num_pump_times, num_observation_points, total_data, pumping_scheme, mass_flowrate, flow_duration,
                                pump_times, pump_rates, time, obs_point_locations, obs_point_num_data, obs_point_property, deliverability, production_index, cutoff_pressure)
end1 = timer()


start2 = timer()
modelled_enthalpy = radial1d_wrapper.radial1d(phi, k, layer_thickness, well_radius, rock_specific_heat, rock_heat_conductivity, rock_density, rock_compressibility, initial_pressure,
                                initial_sv, injection_well, injection_enthalpy, num_pump_times, num_observation_points, total_data, pumping_scheme, mass_flowrate, flow_duration,
                                pump_times, pump_rates, time, obs_point_locations, obs_point_num_data, obs_point_property2, deliverability, production_index, cutoff_pressure)
end2= timer()
print('Time elapsed 1 = {}'.format(end1-start1))
print('Time elapsed 2 = {}'.format(end2-start2))

f, ax1 = plt.subplots(1,2)
# ax2 = ax1[0].twinx()
# ax3 = ax1[1].twinx()
ax1[0].plot(time_days, expected_pressure, 'k-', label='Measured pressure')
ax1[0].plot(time_days2, modelled_pressure/1e5, 'r--', label='Simulated pressure')
ax1[0].legend(loc='best')
# ax1[0].plot(time_days1, modelled_value/1e5, 'r--') # from txt file

# ax2.plot(flow_times_days, flow_rates, 'b-')
measured_enth = ax1[1].plot(time_days, expected_enthalpy, 'k-', label='Measured enthalpy')
# ax1[1].plot(time_days1, actual_enthalpy, 'r--')
sim_enth = ax1[1].plot(time_days2, modelled_enthalpy/1000, 'r--', label='Simulated enthalpy')
# flows = ax3.plot(flow_times_days, flow_rates, 'b--', label='Flow rates')
# lines = measured_enth + sim_enth + flows
# labels = [l.get_label() for l in lines]
# ax1[1].legend(lines, labels, loc='best')
ax1[1].legend(loc='best')
# ax1[0].set_xlim([0,97])
# ax1[1].set_xlim([0,97])

# ax2.set_xlim([0,97])
ax1[0].set_ylim([-5,130])
ax1[1].set_ylim([1500,2400])
# ax2.set_ylim([-6,0])
# ax3.set_ylim([-6,0])

plt.show()

#  Overall radial1d function call time:    10.134601000000000     
#  Overall radial1d_wrapper (inside) time:    10.134609000000001     
# Radial1d wrapper time (from cython): 10.221587896347046
#  Inside homogeneous porous
#  Execution:            2
#  Time elapsed inside homogeneous porous = 
#    10.081318999999999     
#  Overall radial1d function call time:    10.081356000000001     
#  Overall radial1d_wrapper (inside) time:    10.081359000000001     
# Radial1d wrapper time (from cython): 10.172600746154785
# Time elapsed 1 = 10.222235918045044
# Time elapsed 2 = 10.17296290397644