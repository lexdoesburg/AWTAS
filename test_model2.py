import os

import numpy as np
import matplotlib.pyplot as plt
import time as time_module

import model2 as model
import data as data_class


# Import data from text file
data = data_class.Data()
data.read_file('SKG9D_press.txt')
# Create model object
new_model = model.SKG9D(data)
# Check data was read correctly
print(data.time)
print(data.observation)
print(len(data.observation))
# Initialise
# 103.07, 0.2332 optimal
# initial_parameters = np.array([102., .2, 0.082, 2.76e-15])
# X0 = new_model.transf0(initial_parameters)
# r0 = new_model.wellbore_obj(X0)
# print(r0)

initial_parameters = None
# Find solution using model
start = time_module.clock()
parameters = new_model.find_model_parameters(initial_parameters=initial_parameters)
end = time_module.clock()
print('---------\nFinding Parameters Time elapsed = {}'.format(end - start))
print('---------\n---------\n{}\n---------\n---------\n'.format(parameters))

with open('simulation_parameters.txt', 'w') as file:
    file.write('Initial Pressure (Bar): {}\n'.format(parameters[0]))
    file.write('Initial Steam Mass Fraction: {}\n'.format(parameters[1]))
    file.write('Porosity: {}\n'.format(parameters[2]))
    file.write('Permeability: {}\n'.format(parameters[3]))
    file.write('Time for simulation: {}\n'.format(end-start))

# initial_parameters = np.array([102.50286809, 0.22538374])
# initial_parameters = np.array([103.07, 0.2332,0.082,2.76e-15])
# p = new_model.model(initial_parameters)
# t[0] = 0
# print(t)
# print(p)

start = time_module.time()
# plt.semilogx(np.log(data.time), new_model.data.approximation,"rx",label="Model Data")
# plt.semilogx(np.log(data.time), data.observation,"kx",label="Observed Data")
# parameters = [100, 0.2, 0.075, 2e-15]
# approx = new_model.model(parameters)
plt.plot(data.time[:len(data.observation)], new_model.data.approximation, "r-",label="Model Data")
plt.plot(data.time, data.observation,"kx",label="Observed Data")
# plt.plot(t, p, "g-", label="Theis Analytic Solution")
plt.title("Observed Data vs Fitted Curve")
plt.xlabel("Log Time (Days)")
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