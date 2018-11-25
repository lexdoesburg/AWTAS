import numpy as np
from scipy.special import exp1 # pylint: disable-msg=E0611
# from scipy.optimize import curve_fit
from scipy.optimize import leastsq

from datetime import datetime

import data as data_class

class Model():
    def __init__(self, data=None):
    # def __init__(self, p0, qm, h, rho, nu, C, r, t, data=None):
        # self.p0 = p0 # Initial pressure
        # self.qm = qm # Mass flowrate (constant, -ve for production)
        # self.h = h # Thickness
        # self.rho = rho # Density
        # self.nu = nu # Kinematic viscosity
        # self.C = C # Compressibility
        # self.r = r # Radius
        # self.t = t # 1D array of measurement times
        self.data = data # Observed pressure/time measurements

    def model(self, phi, k):
        pass

    def residual_function(self, parameters):
        """
        Calculate the residual of the model (difference between observed data and estimated data)
        """
        phi, k = parameters
        return self.data.observation - self.model(phi, k)

    def find_model_parameters(self, phi=0.1, k=1e-14, curve_fit=False):
        if curve_fit:
            # was going to possibly use scipy.optimize.curve_fit if no initial parameters were given
            pass
        else:
            initial_parameters = np.array([phi, k])
            optimal_parameters, flag = leastsq(self.residual_function, initial_parameters)
            phi, k = optimal_parameters
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
        p0, qm, h, rho, nu, C, r = self.data.parameters
        D = k/(nu*phi*rho*C)
        # u = np.divide(r**2,4*D*t, out=np.zeros_like(t), where = t != 0)
        with np.errstate(divide="ignore", invalid="ignore"): # Hides 'RuntimeWarning: invalid value encountered in divide' if t[0] == 0.
            p = p0 + (qm/(r*np.pi*k*(h/nu)))*exp1((r**2)/(4*D*self.data.time)) # Double check exponential integral is correct
            # p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*t)) # Double check exponential integral is correct
        if self.data.time[0] <= 1e-7: # Check if initial reading is at time 0
            p[0] = p0
        # self.data.set_approximation(p)
        return p