import numpy as np
from scipy.special import exp1 # pylint: disable-msg=E0611
from scipy.optimize import curve_fit
from scipy.optimize import leastsq

def theis_solution(p0, qm, k, h, phi, rho, nu, C, r, t):
    """
    Calculate and return the pressure for the given input.

    Inputs: p0 - Initial pressure
            qm - Mass flowrate (constant, -ve for production)
            k - Permeability
            h - Thickness
            phi - Porosity
            rho - Density
            nu - Kinematic viscosity
            C - Compressibility
            r - Radius of the well
            t - Time
    
    Output: p - pressure at a given well radius and time
    """
    # D = (k/nu)/(phi*rho*C) 
    D = k/(nu*phi*rho*C)
    p = p0 + (qm/(r*np.pi*k*(h/nu)))*exp1((r**2)/(4*D*t)) # Double check exponential integral is correct.
    # p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*t)) # Double check exponential integral is correct.
    # if t[0] == 0:
        # p[0] = 0
    return p
    
# def theis_function(t, phi, k):
#     pass

def theis_residual(parameters, p0, qm, h, rho, nu, C, r, t, data):
    phi, k = parameters
    return data - theis_solution(p0, qm, k, h, phi, rho, nu, C, r, t)

def generate_data(phi, k, n, time, p0, qm, h, rho, nu, C, r, noise = False, sd = 2.5e-4):
    """
    Generate approximated data using the Theis solution for a guess of porosity and permeability.
    """
    t = np.linspace(1, time, num=n)
    # t = np.linspace(0, time, num=n) # this is the line we want
    p = theis_solution(p0, qm, k, h, phi, rho, nu, C, r, t)
    if noise:
        np.random.seed(0)
        p += p*sd*np.random.randn(p.shape[0])
    return p


# def chi_squared(x, data, time, p0, qm, h, rho, nu, C, r):
#     """
#     Non-linear function to be minimised.

#     Inputs: x - 1D array of porosity and permeability
#             data - Known data for a geothermal well-test
#             approximation - Modelled approximation to the data
    
#     Output: chi_squared - Measure of discrepency between data and model response.
#     """
#     n = len(data)
#     approximation = generate_data(x[0], x[1], n, time, p0, qm, h, rho, nu, C, r, noise=True)
#     sd = np.std(data) # Standard deviation of known data
#     return sum(((data-approximation)/sd)**2)


def find_model_parameters(data, p0, qm, h, rho, nu, C, r, t, phi=0.1, k=10e-13):
    """
    Find the model parameters porosity and permeability.

    Inputs: data - Observed data for a geothermal well-test
            phi - Intial guess of porosity (default is arbitrary - from AWTAS page 11)
            k - Initial guess of permeability (default is arbitrary - from AWTAS page 11)
    
    Output: phi - Estimated value of porosity
            k - Estimated value of permeability
    """
    initial_parameters = np.array([phi, k])
    optimal_parameters, flag = leastsq(theis_residual, initial_parameters, args=(p0, qm, h, rho, nu, C, r, t, data))
    phi, k = optimal_parameters
    return phi, k

# Use forward model to determine approximation data.
# Perform non-linear optimisation to get variable parameters phi and permeability.