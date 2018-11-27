import numpy as np
from scipy.special import exp1
# from scipy.optimize import curve_fit
from scipy.optimize import leastsq

import data as data_class

class Model():
    def __init__(self, data=None):
        """
        Initialise the model
        """
        self.data = data # Data structure containing observed pressure measurements, time of measurements and well parameters

    def model(self, phi, k):
        """
        Model function.
        """
        raise NotImplementedError('Model function has not been implemented.')

    def residual_function(self, parameters):
        """
        Calculate the residual of the model (difference between observed data and estimated data)
        """
        phi, k = parameters
        return self.data.observation - self.model(phi, k)

    def __chi_squared(self):
        sd = np.std(self.data.observation)
        chi_squared = np.sum((self.data.observation-self.data.approximation)/sd)**2
        return chi_squared

    def find_model_parameters(self, phi=None, k=None, curve_fit=False):
        """
        Find the model parameters porosity and permeability.

        Inputs: phi (float) - Intial guess of porosity (default is arbitrary - from AWTAS page 11)
                k (float) - Initial guess of permeability (default is arbitrary - from AWTAS page 11)
                curve_fit (boolean)
        
        Output: phi (float) - Estimated value of porosity
                k (float)- Estimated value of permeability

        NOTES:
        - A suitable range for porosity is 0.01 - 0.2, permeability is 1e-16—1e-12
        - Leastsq is highly sensitive to permeability (with porosity = 0.1): Doesn't find optimal permeability for 1e-15—1e-16 and for 1e-9 and upwards.
            - Also sensitive to porosity (with k = 1e-14): Doesn't find either parameter for porosity > 0.125 and any value < 0.1
            - Actual value of porosity = 0.1, permeability = 1e-12. (with permeability guess of 1e-13, porosity guess can be less accurate).
        - Look at using either curve_fit or least_squares instead.
        """

        if curve_fit:
            # was going to possibly use scipy.optimize.curve_fit if no initial parameters were given
            pass
        else:
            calls = 0
            # broken = False
            estimates = np.ndarray(shape=(3,25))
            i = 0
            for k in [1e-16, 1e-15, 1e-14, 1e-13, 1e-12]:
                for phi in [0.2, 0.15, 0.1, 0.05, 0.01]:
                    initial_parameters = np.array([phi, k])
                    optimal_parameters, flag = leastsq(self.residual_function, initial_parameters) # if flag is 1 - found a good soln
                    self.data.set_approximation(self.model(optimal_parameters[0], optimal_parameters[1])) # Store the approximated data in the data structure
                    chi_squared = self.__chi_squared()
                    estimates[:, i] = optimal_parameters[0], optimal_parameters[1], chi_squared
                    print('Phi: {} k: {} Chi squared: {}'.format(phi, k, chi_squared))
                    i += 1
                    calls += 1
                #     if chi_squared <= 20:
                #         broken = True
                #         break
                # if broken:
                #     break
            index = np.argmin(estimates[2])
            print('index = {} chi squared = {}'.format(index, estimates[2, index]))
            print('Function called: {} times'.format(calls))

            phi, k = estimates[:2, index]
            print('Phi {}, k {}'.format(phi, k))
            # phi, k = optimal_parameters
        self.data.set_unknown_parameters(phi, k) # Store phi and k in data structure
        self.data.set_approximation(self.model(phi, k)) # Store the approximated data in the data structure
        return phi, k

    def generate_data(self, phi, k, time, parameters, noise = False, sd = 2.5e-4, save_file=False, filename="example_datafile.txt"):
        """
        Generate approximated data using the Theis solution for a guess of porosity and permeability.
        """
        # filename="datafile_{}.txt".format(datetime.now().strftime("%d-%M-%Y_%H:%M") # Argument
        self.data = data_class.Data()
        self.data.set_time(time)
        self.data.set_known_parameters(parameters)
        p = self.model(phi, k)
        if noise:
            np.random.seed(0) # Set random seed to 0 for consistency in testing
            p += p*sd*np.random.randn(p.shape[0])

        if save_file:
            self.__generate_datafile(filename, p)
        self.data.set_observation(p)
        return self.data
    
    def __generate_datafile(self, filename, measurement):
        with open(filename, 'w') as file:
            file.write('Time(s),Pressure(Pa)\n')
            for i in range(len(self.data.time)):
                if self.data.time[i] >= 1e-7:
                    file.write('{},{}\n'.format(self.data.time[i], measurement[i]))
                else:
                    file.write('{},{}\n'.format(0, measurement[i]))


class Theis_Solution(Model):
    def model(self, phi, k):
        """
        Calculate and return the pressure for the given input using the analytical Theis solution.

        Inputs: phi - Porosity
                k - Permeability
                Makes use of self.data for remaining parameters
        
        Output: p(t) - pressure at a given time (or time array).
        """
        p0, qm, h, rho, nu, C, r = self.data.parameters # Unpack the well parameters
        D = k/(nu*phi*rho*C)
        with np.errstate(divide="ignore", invalid="ignore"): # Hides 'RuntimeWarning: invalid value encountered in divide' if t[0] == 0.
            p = p0 + (qm/(r*np.pi*k*(h/nu)))*exp1((r**2)/(4*D*self.data.time)) # TODO: Double check exponential integral is correct
            # p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*t)) # TODO: Double check exponential integral is correct
        if self.data.time[0] <= 1e-7: # Check if initial reading is at time 0
            p[0] = p0
        return p