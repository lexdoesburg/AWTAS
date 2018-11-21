import numpy as np
from scipy.special import exp1 # pylint: disable-msg=E0611
from scipy.optimize import newton
# import Numdifftools as nd

def theis_solution(p0, q, k, h, porosity, density, v, C, r, t):
    """
    Calculate and return the pressure for the given input.

    Inputs: p0 - Initial pressure
            q - Mass flowrate (constant, -ve for production)
            k - Permeability
            h - Thickness
            porosity - Porosity
            density - Density
            v - Kinematic viscosity
            C - Compressibility
            r - Radius of the well
            t - Time
    
    Output: p - pressure
    """
    D = (k/v)/(porosity*density*C) #
    p = p0 + (q/(r*np.pi*k*(h/v)))*exp1((r**2)/4*D*t) # Double check exponential integral is correct.
    return p
    

def find_model_parameters(data, porosity, k):
    """
    Find the model parameters porosity and permeability.

    Inputs: data - Known data for a geothermal well-test
            porosity - Porosity
            k - Permeability
    
    Output: porosity - Porosity
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

