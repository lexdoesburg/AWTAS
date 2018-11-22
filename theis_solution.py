import numpy as np
from scipy.special import exp1 # pylint: disable-msg=E0611
from scipy.optimize import curve_fit
# import Numdifftools as nd

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
    x = (r**2)/(4*D*t)

    # Alternate approximation to the E1(u)
    gamma = np.euler_gamma
    wu = -gamma - np.log((r**2)/(4*D*t))
    p1 = p0 + (qm/(r*np.pi*k*(h/nu)))*wu
    # Using scipy exp1
    e1 = exp1((r**2)/(4*D*t))
    p = p0 + (qm/(r*np.pi*k*(h/nu)))*exp1((r**2)/(4*D*t)) # Double check exponential integral is correct.
    # p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*t)) # Double check exponential integral is correct.
    # if t[0] == 0:
        # p[0] = 0
    return p
    
def theis_function():
    pass

def generate_data(phi, k, n, time, p0, qm, h, rho, nu, C, r):
    """
    Generate approximated data using the Theis solution for a guess of porosity and permeability.
    """
    t = np.linspace(0, time, num=n)
    return theis_solution(p0, qm, k, h, phi, rho, nu, C, r, t)


def chi_squared(x, data, time, p0, qm, h, rho, nu, C, r):
    """
    Non-linear function to be minimised.

    Inputs: x - 1D array of porosity and permeability
            data - Known data for a geothermal well-test
            approximation - Modelled approximation to the data
    
    Output: chi_squared - Measure of discrepency between data and model response.
    """
    n = len(data)
    approximation = generate_data(x[0], x[1], n, time, p0, qm, h, rho, nu, C, r)
    sd = np.std(data) # Standard deviation of known data
    return sum((data-approximation)/sd)**2

# Use forward model to determine approximation data.
# Perform non-linear optimisation to get variable parameters phi and permeability.

def find_model_parameters(data, phi=0.1, k=10e-13):
    """
    Find the model parameters porosity and permeability.

    Inputs: data - Known data for a geothermal well-test
            phi - Intial guess of porosity
            k - Initial guess of permeability
    
    Output: phi - Estimated value of porosity
            k - Estimated value of permeability
    """
    pass