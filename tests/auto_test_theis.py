import pytest
import numpy as np
import os

from awtas.logic import model
from awtas.logic import data

# NOTE: Tests fail since data structure was updated - need to update tests.

#----------------------------------
#   To run this test: Type in command line interface
#           python -m pytest -v -s <path to file>
#----------------------------------

test_file_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_datafiles')
theis_testcases = ['theis_testcase1.txt', 'theis_testcase2.txt']

def get_expected_output(testfile):
    with open(testfile, 'r') as file:
        file.readline() # Skip header
        parameters_string = file.readline()
        parameters = [float(value) for value in parameters_string.split(',')]
        # print(parameters)
        file.readline() # Skip more metadata
        variables_string = file.readline()
        variables = [float(value) for value in variables_string.split(',')]
        # print(variables)
    time, expected_pressure = np.genfromtxt(testfile, delimiter=',', skip_header=5).T
    return variables, parameters, time, expected_pressure

def get_actual_output(variables, parameters, time, model_type):
    test_data = data.Data(time=time, parameters_list=parameters, model_type=model_type)
    if model_type == 'theis':
        test_model = model.Theis_Solution(test_data)
    elif model_type == 'radial1d':
        # test_model = model.Radial_1D(test_data)
        pass
    actual_pressure = test_model.model(variables)
    return actual_pressure

@pytest.mark.parametrize('theis_testfile', [os.path.join(test_file_dir, testcase) for testcase in theis_testcases])

def test_theis_output(theis_testfile):
    variables, parameters, time, expected_pressure = get_expected_output(theis_testfile)
    actual_pressure = get_actual_output(variables, parameters, time, model_type='theis')
    # print(actual_pressure - expected_pressure)
    assert np.allclose(actual_pressure, expected_pressure)

# @pytest.mark.parametrize('radial1d_testfile', ['radial1d_testcase1.txt'])

# def test_radial1d_output(radial1d_testfile):
#     variables, parameters, time, expected_pressure = get_expected_output(radial1d_testfile)
#     actual_pressure = get_actual_output(variables, parameters, time, model_type='radial1d')
#     print(actual_pressure - expected_pressure)
#     assert np.allclose(actual_pressure, expected_pressure)