import numpy as np
from scipy.special import exp1 # pylint: disable-msg=E0611
from scipy.optimize import newton
# import Numdifftools as nd

def theis_solution(p0, q, k, h, phi, rho, nu, C, r, t):
    """
    Calculate and return the pressure for the given input.

    Inputs: p0 - Initial pressure
            q - Mass flowrate (constant, -ve for production)
            k - Permeability
            h - Thickness
            phi - Porosity
            rho - Density
            nu - Kinematic viscosity
            C - Compressibility
            r - Radius of the well
            t - Time
    
    Output: p - pressure
    """
    D = (k/nu)/(phi*rho*C) #
    p = p0 + (q/(r*np.pi*k*(h/nu)))*exp1((r**2)/4*D*t) # Double check exponential integral is correct.
    return p
    

def find_model_parameters(data, phi, k):
    """
    Find the model parameters phi and permeability.

    Inputs: data - Known data for a geothermal well-test
            phi - Porosity
            k - Permeability
    
    Output: phi - Porosity
            k - Permeability
    """
    pass


def chi_squared(data, approximation):
    """
    Non-linear function to be minimised.

    Inputs: data - Known data for a geothermal well-test
            approximation - Modelled approximation to the data
    
    Output: chi_squared - Measure of discrepency between data and model response.
    """
    assert(len(data) == len(approximation), "Number of data points are not equal.") 
    sd = np.std(data) # Standard deviation of known data
    return sum((data-approximation)/sd)**2

# Use forward model to determine approximation data.
# Perform non-linear optimisation to get variable parameters phi and permeability.
