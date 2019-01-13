import pytest
import numpy as np
import model
import data

def get_expected_output(test_file):
    with open(test_file, 'r') as file:
        file.readline() # Skip header
        parameters_string = file.readline()
        parameters = [float(value) for value in parameters_string.split(',')]
        # print(parameters)
        file.readline() # Skip more metadata
        variables_string = file.readline()
        variables = [float(value) for value in variables_string.split(',')]
        # print(variables)
    time, expected_pressure = np.genfromtxt(test_file, delimiter=',', skip_header=5).T
    return variables, parameters, time, expected_pressure

def get_actual_output(variables, parameters, time):
    theis_data = data.Data(time=time, parameters=parameters)
    theis_model = model.Theis_Solution(theis_data)
    actual_pressure = theis_model.model(variables)
    return actual_pressure

@pytest.mark.parametrize('test_file', ['theis_testcase1.txt'])

def test_theis_output(test_file):
    variables, parameters, time, expected_pressure = get_expected_output(test_file)
    actual_pressure = get_actual_output(variables, parameters, time)
    # print(actual_pressure - expected_pressure)
    assert np.allclose(actual_pressure, expected_pressure)