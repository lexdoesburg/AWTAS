import os

import numpy as np
import matplotlib.pyplot as plt
import time as time_module

import model2 as model
import data2 as data_class


# Import data from text file
data = data_class.Data()
data.read_file('SKG9D_press.txt')
# Create model object
new_model = model.SKG9D(data)
# Check data was read correctly
print(data.time)
print(data.observation)
print('length', len(data.observation))
# Initialise
# 103.07, 0.2332 optimal
# initial_parameters = np.array([102., .2, 0.082, 2.76e-15])
# X0 = new_model.transf0(initial_parameters)
# r0 = new_model.wellbore_obj(X0)
# print(r0)

# # Find solution using model
# start = time_module.clock()
# guess = np.array([0.04, 2e-15])
# parameters = new_model.find_model_parameters2(guess, verbose=True)
# end = time_module.clock()
# print('---------\nFinding Parameters Time elapsed = {}'.format(end - start))
# print('---------\n---------\n{}\n---------\n---------\n'.format(parameters))

# with open('simulation_parameters.txt', 'a') as file:
#     # file.write('Initial Pressure (Bar): {}\n'.format(parameters[0]))
#     # file.write('Initial Steam Mass Fraction: {}\n'.format(parameters[1]))
#     file.write('Porosity: {}\n'.format(parameters[0]))
#     file.write('Permeability: {}\n'.format(parameters[1]))
#     file.write('Time for simulation: {}\n'.format(end-start))

# initial_parameters = np.array([102.50286809, 0.22538374])
# initial_parameters = np.array([103.07, 0.2332,0.082,2.76e-15])

optimal_parameters = [0.145, 2.5e-15]
# optimal_parameters = parameters

c = new_model.model(optimal_parameters)
c_chi = new_model.chi_squared()
print(c_chi)

p = new_model.model([0.082, 2.76e-15])
chi_squared = new_model.chi_squared()
print(chi_squared)
# print(p)

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




# start = time_module.time()
# plt.semilogx(np.log(data.time), new_model.data.observation,"kx",label="Observation Data")
# # plt.semilogx(np.log(time), theis_model.data.approximation,"r-",label="Approximated Curve")
# # plt.plot(t, p, "g-", label="Theis Analytic Solution")
# plt.title("Observed Data vs Fitted Curve")
# plt.xlabel("Log Time (Days)")
# plt.ylabel("Pressure (Bar)")
# plt.legend(loc="best")
# end = time_module.time()
# print('Time elapsed = {}'.format(end - start))
# plt.show()