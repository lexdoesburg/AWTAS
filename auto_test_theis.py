import pytest
import numpy as np
import model
import data

#----------------------------------
#   To run this test: Type in command line interface
#       [Windows]   python -m pytest -v -s H:\\Summer_Project\\AWTAS\\auto_test_theis.py
#       [Mac]       python -m pytest -v -s /Users/lexdoesburg/Documents/Uni2018/Summer_Research/Summer_Project/AWTAS/auto_test_theis.py
#----------------------------------

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

@pytest.mark.parametrize('theis_testfile', ['theis_testcase1.txt', 'theis_testcase2.txt'])

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