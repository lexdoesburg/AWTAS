#!/usr/bin/env python

# <examples/doc_parameters_basic.py>
import numpy as np
import model
import data as data_class
from lmfit import Minimizer, Parameters, report_fit, Model, minimize
import time
from scipy.special import exp1
# create data to be fitted
data = data_class.Data()
data.read_file('Theis_Solution_test1.dat')
t = data.time
model_obj = model.Theis_Solution()

# define objective function: returns the array to be minimized
def fcn2min(params, t, data):
    """Model a decaying sine wave and subtract data."""
    model = test_model(params, t)
    return model - data


def test_model(parameters, t):
    phi = parameters['phi']
    k = parameters['k']
    p0 = parameters['p0']
    qm = parameters['qm']
    rho = parameters['rho']
    nu = parameters['nu']
    h = parameters['h']
    C = parameters['C']
    r = parameters['r']

    D = k/(nu*phi*rho*C)
    with np.errstate(divide="ignore", invalid="ignore"): # Hides 'RuntimeWarning: invalid value encountered in divide' if t[0] == 0.
        # p = p0 + (qm/(r*np.pi*k*(h/nu)))*exp1((r**2)/(4*D*self.data.time)) 
        p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*t)) # Same as fortran output
    if t[0] <= 1e-7: # Check if initial reading is at time 0
        p[0] = p0
    return p


# ---------------------
#   Appears to fail for very small changes in pressure
# ---------------------


# create a set of Parameters
params = Parameters()
# Test 1
params.add_many(('phi', None, True, 0.01, 0.2),
                ('k', None, True, 1e-16, 1e-12),
                ('p0', 3.6e6, False),
                ('qm', -0.005, False),
                ('h', 100, False),
                ('rho', 813.37, False),
                ('nu', 0.0001111, False),
                ('C', 0.001303, False),
                ('r', 0.05, False))

# # Test 3
# params.add_many(('phi', None, True, 0.01, 0.2),
#                 ('k', None, True, 1e-16, 1e-12),
#                 ('p0', 3.6e6, False),
#                 ('qm', -0.015, False),
#                 ('h', 100, False),
#                 ('rho', 527.59, False),
#                 ('nu', 0.0000603, False),
#                 ('C', 0.03748, False),
#                 ('r', 0.05, False))


# do fit, here with leastsq model
start = time.clock()
# fit_kws={'diag':(1/t.mean(), 1./data.observation.mean())}
result = minimize(fcn2min, params, method='leastsq', args=(t,), kws={'data':data.observation}, diag=(1/t.mean(), 1./data.observation.mean()))
end = time.clock()
print('Time to find solution: ',end-start)
print(report_fit(result))

final = data.observation + result.residual

# try to plot results
try:
    import matplotlib.pyplot as plt
    plt.plot(t, data.observation, 'k+')
    plt.plot(t, final, 'r')
    plt.show()
except ImportError:
    pass
# <end of examples/doc_parameters_basic.py>