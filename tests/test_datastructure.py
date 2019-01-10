import numpy as np
import pytest
import data as datastructure
import model

# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.005 # m^3/s
# k = 1e-12 # m^2
# phi = 0.1
# rho = 813.37 # Water at 240 degrees celsius
# nu = 0.0001111 # Water at 240 degrees celsius
# C = 0.001303 # Water at 240 degrees celsius
# t = 43200 # seconds
# time = np.linspace(0, 43200, num=100)
# parameters = [p0, qm, h, rho, nu, C, r]

# Test 1 - Test reading file
# print("Test 1")
# data_1 = datastructure.Data('testdata.txt', parameters=parameters)

# print(data_1.time)
# print(data_1.observation)
# print(data_1.parameters) 
# print(data_1.approximation) # None
# print(data_1.phi) # None
# print(data_1.k) # None

# # Test 2 - Set Parameters
# print("\nTest 2")
# data_2 = datastructure.Data()
# print(data_2.time)
# print(data_2.observation)
# print(data_2.parameters) 
# print(data_2.approximation) # None
# print(data_2.phi) # None
# print(data_2.k) # None

# data_2.set_known_parameters(parameters)
# print(data_2.parameters)

# data_1.generate_datafile('test_generation.txt')
# data_3 = datastructure.Data('test_generation.txt')
# print(data_3.time)
# print(data_3.observation)
# print(data_3.parameters)

data = datastructure.Data()
def generate_parameters_dictionary(parameters, model_type):
    return data.create_parameter_dictionary(parameters, model_type)

def test_create_parameter_dictionary():
    print(generate_parameters_dictionary(None,None))

test_create_parameter_dictionary()
