import numpy as np
from scipy.special import exp1 # pylint: disable-msg=E0611
# from scipy.optimize import curve_fit
from scipy.optimize import leastsq

def theis_solution(p0, qm, k, h, phi, rho, nu, C, r, t):
    """
    Calculate and return the pressure for the given input using the analytical Theis solution.

    Inputs: p0 - Initial pressure
            qm - Mass flowrate (constant, -ve for production)
            k - Permeability
            h - Thickness
            phi - Porosity
            rho - Density
            nu - Kinematic viscosity
            C - Compressibility
            r - Radius of the well
            t - 1D array of measurement times
    
    Output: p(t) - pressure at a given well radius
    """
    D = k/(nu*phi*rho*C)
    # u = np.divide(r**2,4*D*t, out=np.zeros_like(t), where = t != 0)
    with np.errstate(divide="ignore", invalid="ignore"): # Hides 'RuntimeWarning: invalid value encountered in divide' if t[0] == 0.
        p = p0 + (qm/(r*np.pi*k*(h/nu)))*exp1((r**2)/(4*D*t)) # Double check exponential integral is correct
        # p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*t)) # Double check exponential integral is correct
    if t[0] <= 1e-7: # Check if initial reading is at time 0
        p[0] = p0
    return p
    

def theis_residual(parameters, p0, qm, h, rho, nu, C, r, t, data):
    """
    Calculate the residual of the Theis solution (difference between observed data and estimated data)
    """
    phi, k = parameters
    return data - theis_solution(p0, qm, k, h, phi, rho, nu, C, r, t)

def generate_data(phi, k, n, time, p0, qm, h, rho, nu, C, r, noise = False, sd = 2.5e-4, save_file=False):
    """
    Generate approximated data using the Theis solution for a guess of porosity and permeability.
    """
    t = np.linspace(0, time, num=n)
    p = theis_solution(p0, qm, k, h, phi, rho, nu, C, r, t)
    if noise:
        np.random.seed(0) # Set random seed to 0 for consistancy in testing
        p += p*sd*np.random.randn(p.shape[0])

    if save_file:
        generate_datafile("testdata.txt", t, p)
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


# def theis_function(t, phi, k):
#     pass


def find_model_parameters(data, p0, qm, h, rho, nu, C, r, t, phi=0.1, k=1e-14, curve_fit=False):
    """
    Find the model parameters porosity and permeability.

    Inputs: data - Observed data for a geothermal well-test
            p0 - Initial pressure
            qm - Mass flowrate (constant, -ve for production)
            k - Permeability
            h - Thickness
            phi - Porosity
            rho - Density
            nu - Kinematic viscosity
            C - Compressibility
            r - Radius of the well
            t - 1D array of measurement times
            phi - Intial guess of porosity (default is arbitrary - from AWTAS page 11)
            k - Initial guess of permeability (default is arbitrary - from AWTAS page 11)
    
    Output: phi - Estimated value of porosity
            k - Estimated value of permeability
    """

    """
    NOTES:
        - A suitable range for porosity is 0.01 - 0.2, permeability is 1e-16—1e-12
        - Leastsq is highly sensitive to permeability (with porosity = 0.1): Doesn't find optimal permeability for 1e-15—1e-16 and for 1e-9 and upwards.
            - Also sensitive to porosity (with k = 1e-14): Doesn't find either parameter for porosity > 0.125 and any value < 0.1
            - Actual value of porosity = 0.1, permeability = 1e-12. (with permeability guess of 1e-13, porosity guess can be less accurate).
    """

    if curve_fit:
        # was going to possibly use scipy.optimize.curve_fit if no initial parameters were given
        pass
    else:
        initial_parameters = np.array([phi, k])
        optimal_parameters, flag = leastsq(theis_residual, initial_parameters, args=(p0, qm, h, rho, nu, C, r, t, data))
        phi, k = optimal_parameters
    return phi, k

# Use forward model to determine approximation data.
# Perform non-linear optimisation to get variable parameters phi and permeability.

def generate_datafile(filename, time, measurement):
    with open(filename, 'w') as file:
        file.write('Log Time(s),Pressure(Pa)\n')
        for i in range(len(time)):
            if time[i] >= 1e-7:
                file.write('{},{}\n'.format(time[i], measurement[i]))
            else:
                file.write('{},{}\n'.format(0, measurement[i]))


def read_data(filename):
    x_data, y_data = np.genfromtxt(filename, delimiter=',', skip_header=1).T
    # print(x_data)
    # print(y_data)
    return x_data, y_data