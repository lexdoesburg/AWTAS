import os

import numpy as np
import matplotlib.pyplot as plt
import time as time_module

import model2 as model
import data2 as data_class


def transf0(theta):
    """transform parameter vector to a new base adapted for automated calibration

    Args:
        theta (list): parameter vector in original base
    """
    xlim = [0.01, 0.2]	# initial phi range
    ylim = [1e-16, 1e-12]	# initial k range
    return [(theta[0]-xlim[0])/(xlim[1]-xlim[0]), (theta[1]-ylim[0])/(ylim[1]-ylim[0])]

# Import data from text file
data = data_class.Data()
data.read_file('SKG9D_press.txt')
# data.read_file('SKG9D_test1.dat')
# data.set_error(1.2)
# Create model object
new_model = model.SKG9D(data)

# new_model.generate_data([0.1, 2.6e-15], data.time, parameters=None, noise = False, sd = 300/1e5, save_file=True, filename='SKG9D_test1.dat')
# plt.plot(data.time, new_model.data.observation,'r-')
# new_model.generate_data([0.1, 2.6e-15], data.time, parameters=None, noise = True, sd = 1.2, save_file=True, filename='SKG9D_test1.dat')
# plt.plot(data.time, new_model.data.observation,'kx')
# plt.show()

# Check data was read correctly
# print(data.time)
# print(data.observation)
# print('length', len(data.observation))

# Initialise
# 103.07, 0.2332 optimal
# initial_parameters = np.array([102., .2, 0.082, 2.76e-15])
# X0 = new_model.transf0(initial_parameters)
# r0 = new_model.wellbore_obj(X0)
# print(r0)

# Find solution using model
start = time_module.clock()
guess = np.array([0.04, .6e-15])
# guess = None
parameters = new_model.find_model_parameters2(guess, verbose=True, single_run=False)

end = time_module.clock()
print('---------\nFinding Parameters Time elapsed = {}'.format(end - start))
print('---------\n---------\n{}\n---------\n---------\n'.format(parameters))

with open('simulation_parameters.txt', 'a') as file:
    # file.write('Initial Pressure (Bar): {}\n'.format(parameters[0]))
    # file.write('Initial Steam Mass Fraction: {}\n'.format(parameters[1]))
    file.write('Porosity: {}\n'.format(parameters[0]))
    file.write('Permeability: {}\n'.format(parameters[1]))
    file.write('Time for simulation: {}\n'.format(end-start))

# initial_parameters = np.array([102.50286809, 0.22538374])
# initial_parameters = np.array([103.07, 0.2332,0.082,2.76e-15])

# optimal_parameters = [0.015755192744286974, 3.703236263357016e-15]
# optimal_parameters = [0.01570910057030036, 3.703079600761992e-15]
# optimal_parameters = [1.04999999e-01, 2.58974857e-15]
# optimal_parameters = [0.10499999948594617, 2.602715815877674e-15]
# optimal_parameters = [0.07192956166621617, 2.9056658976041045e-15]
# optimal_parameters = [1.12315854e-01, 2.53425890e-15]
# optimal_parameters = [5.36353651e-02, 3.14698023e-15]
optimal_parameters = parameters
c = new_model.model(optimal_parameters)
# c = new_model.model(transf0(optimal_parameters))
c_chi = new_model.chi_squared()
print(c_chi)

p = new_model.model([0.082, 2.76e-15])
# p = new_model.model(transf0([0.082, 2.76e-15]))
# p = new_model.model(transf0([0.1, 2.6e-15]))
chi_squared = new_model.chi_squared()
print(chi_squared)
# # print(p)

start = time_module.time()
# plt.semilogx(np.log(data.time[:len(new_model.data.approximation)]), c, "r-",label="Approximation")
# plt.semilogx(np.log(data.time), p, "g-", label="Optimal")
# plt.semilogx(np.log(data.time), data.observation,"kx",label="Observed Data")
# plt.xlabel("Log Time (Days)")
plt.plot(data.time, data.observation,"kx",label="Observed Data")
plt.plot(data.time[:len(c)], c, "r-",label="Approximation")
plt.plot(data.time, p, "k--", label="Optimal")
plt.xlabel("Time (Days)")
plt.title("Observed Data vs Fitted Curve")
plt.ylabel("Pressure (Bar)")
plt.legend(loc="best")
end = time_module.time()
# print('Time elapsed = {}'.format(end - start))
plt.show()

