import numpy as np
import pytest

import sys
import os
import inspect
current_dir = os.path.dirname(os.path.realpath(inspect.getfile(inspect.currentframe()))) # https://stackoverflow.com/questions/714063/importing-modules-from-parent-folder/11158224#11158224
awtas_main_dir = os.path.dirname(current_dir)
sys.path.insert(0, awtas_main_dir)
# To run without changing sys path use 'python -m tests.test_datastructure' in command line while in main AWTAS directory
from awtas.logic import model
from awtas.logic import data as datastructure

# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.005 # m^3/s
# k = 1e-12 # m^2
# phi = 0.1
# rho = 813.37 # Water at 240 degrees celsius
# nu = 0.0001111 # Water at 240 degrees celsius
# C = 0.001303 # Water at 240 degrees celsius
# # t = 43200 # seconds
# time = np.linspace(0, 43200, num=100)
# parameters = [p0, qm, h, rho, nu, C, r]

# # Test 1 - Test reading file
# # print("Test 1")
# # data_1 = datastructure.Data('testdata.txt', parameters=parameters)

# # print(data_1.time)
# # print(data_1.observation)
# # print(data_1.parameters) 
# # print(data_1.approximation) # None
# # print(data_1.phi) # None
# # print(data_1.k) # None

# # # Test 2 - Set Parameters
# # print("\nTest 2")
# # data_2 = datastructure.Data()
# # print(data_2.time)
# # print(data_2.observation)
# # print(data_2.parameters) 
# # print(data_2.approximation) # None
# # print(data_2.phi) # None
# # print(data_2.k) # None

# # data_2.set_known_parameters(parameters)
# # print(data_2.parameters)

# # data_1.generate_datafile('test_generation.txt')
# # data_3 = datastructure.Data('test_generation.txt')
# # print(data_3.time)
# # print(data_3.observation)
# # print(data_3.parameters)

# # data = datastructure.Data()
# # def generate_parameters_dictionary(parameters, model_type):
# #     return data.create_parameter_dictionary(parameters, model_type)

# # def test_create_parameter_dictionary():
# #     print(generate_parameters_dictionary(None,None))

# # test_create_parameter_dictionary()

# data = datastructure.Data(model_type='Theis')

# print(data.parameters)
# data.fill_parameter_dictionary(parameters)
# print(data.parameters)


# time, pressure = np.genfromtxt('Theis_Solution_test1.dat', delimiter=',', skip_header=4).T
# data.time = time
# data.observation = pressure
# data.generate_datafile('xd123.txt')

# # ----
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

# data = datastructure.Data(model_type='Radial1d')

# print(data.parameters)
# data.fill_parameter_dictionary(parameters)
# print(data.parameters)


# time, pressure = np.genfromtxt('Theis_Solution_test1.dat', delimiter=',', skip_header=4).T
# data.time = time
# data.observation = pressure
# data.generate_datafile('xd12345.txt')


# new_data = datastructure.Data(model_type='Radial1d',filename='xd12345.txt')
# print()
# print()
# print(new_data.parameters)
# print(new_data.time)
# print(new_data.observation)
