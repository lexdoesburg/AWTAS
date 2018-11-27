import numpy as np
from scipy.special import exp1
from scipy.optimize import curve_fit, leastsq, least_squares
# from scipy.optimize import leastsq

import data as data_class

class Model():
    def __init__(self, data=None):
        """
        Initialise the model
        """
        self.data = data # Data structure containing observed pressure measurements, time of measurements and well parameters

    def model(self, time, phi, k):
        """
        Model function.
        """
        raise NotImplementedError('Model function has not been implemented.')

    def residual_function(self, parameters):
        """
        Calculate the residual of the model (difference between observed data and estimated data)
        """
        phi, k = parameters
        return self.data.observation - self.model(self.data.time, phi, k)
    
    def weighted_error(self):
        sd = np.std(self.data.observation)
        return np.sum((self.data.observation-self.data.approximation)/sd)**2

    def find_model_parameters(self, phi=10, k=1e-20, curvefit=False, leastsquares=False, bounded=False):
        initial_parameters = np.array([phi, k])
        if curvefit:
            if bounded:
                popt, pcov = curve_fit(self.model, self.data.time, self.data.observation, initial_parameters, bounds=([0.005, 1e-17], [0.205, 1e-11]))
            else:
                popt, pcov = curve_fit(self.model, self.data.time, self.data.observation, initial_parameters)
            phi, k = popt
        elif leastsquares:
            if bounded:
                result = least_squares(self.residual_function, initial_parameters, bounds=([0.005, 1e-17], [0.205, 1e-11]))
            else:
                result = least_squares(self.residual_function, initial_parameters)
            phi, k = result.x
            # print('Cost = {}'.format(result.cost))
        else:
            optimal_parameters, flag = leastsq(self.residual_function, initial_parameters) # if flag is 1 - found a good soln
            phi, k = optimal_parameters
        self.data.set_unknown_parameters(phi, k) # Store phi and k in data structure
        self.data.set_approximation(self.model(self.data.time, phi, k)) # Store the approximated data in the data structure
        return phi, k

    def generate_data(self, phi, k, time, parameters, noise = False, sd = 2.5e-4, save_file=False, filename="example_datafile.txt"):
        """
        Generate approximated data using the Theis solution for a guess of porosity and permeability.
        """
        # filename="datafile_{}.txt".format(datetime.now().strftime("%d-%M-%Y_%H:%M") # Argument
        self.data = data_class.Data()
        self.data.set_time(time)
        self.data.set_known_parameters(parameters)
        p = self.model(self.data.time, phi, k)
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
    def model(self, time, phi, k):
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
            p = p0 + (qm/(r*np.pi*k*(h/nu)))*exp1((r**2)/(4*D*time)) # TODO: Double check exponential integral is correct
            # p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*t)) # TODO: Double check exponential integral is correct
        if time[0] <= 1e-7: # Check if initial reading is at time 0
            p[0] = p0
        return p