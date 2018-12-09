import numpy as np
import matplotlib.pyplot as plt
import time as time_module

import model
import data as data_class

def setup_model(model_type, parameters, variables, time, test_num, noise=True, sd=1e-5, generate_datafile=True):
    model_type = model_type.lower()
    if model_type == 'theis':
        test_model = model.Theis_Solution()
    elif model_type == 'theis fortran':
        test_model = model.Theis_Solution_Fortran()
    elif model_type == 'radial 1d':
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
t = 15 # days
t = t * 24 * 60 * 60
print(t)
time = np.linspace(0, t, num=1000)
parameters = [p0, qm, h, rho, nu, C, r]
test1_variables = [phi, k]

# Build model
# theis_test_1 = setup_model('theis', parameters, test1_variables, time, test_num=1, sd=1e-4)
theis_test_1 = setup_model('theis', parameters, test1_variables, time, test_num=1, sd=500)

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
theis_test_2 = setup_model('theis', parameters, test2_variables, time, test_num=2, sd=20)

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
theis_test_3 = setup_model('theis', parameters, test3_variables, time, test_num=3, sd=0.0001)

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
theis_test_4 = setup_model('theis', parameters, test4_variables, time, test_num=4, sd=30)

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

optimal_parameters = theis_test_1.find_model_parameters()
plot_solution(theis_test_1)

optimal_parameters = theis_test_2.find_model_parameters()
plot_solution(theis_test_2)

optimal_parameters = theis_test_3.find_model_parameters()
plot_solution(theis_test_3)

optimal_parameters = theis_test_4.find_model_parameters()
plot_solution(theis_test_4)
# # -----------------------------------------------------------------------------------

# # ----------------------------------------------------------------------------------
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
# # ----------------------------------------------------------------------------------





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