import pytest
import numpy as np
import time as time_module
import model
import data as datastructure

def get_expected_output(filename):
    time, pressure = np.genfromtxt(filename, delimiter=',', skip_header=12).T
    _, values = np.genfromtxt(filename, delimiter=',', skip_footer=len(time)).T
    parameters = values[:-2]
    variables = values[-2:]
    return variables, parameters, time, pressure

# def get_model_parameters(filename):
#     _, values = np.genfromtxt(filename, delimiter=',').T
#     parameters = values[:-2]
#     variables = values[-2:]
#     return variables, parameters

def setup_model(time, parameters):
    theis_data = datastructure.Data(time=time, parameters_list=parameters, model_type='theis')
    theis_model = model.Theis_Solution(data=theis_data)
    return theis_model

def get_actual_output(model, variables):
    pressure = model.model(variables)
    return pressure

def test_theis_output1():
    variables, parameters, time, expected_pressure = get_expected_output('theis_test1_output.txt')
    theis_model = setup_model(time, parameters)
    actual_pressure = get_actual_output(theis_model, variables)
    assert actual_pressure == expected_pressure

