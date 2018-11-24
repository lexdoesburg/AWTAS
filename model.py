import numpy as np
from scipy.special import exp1 # pylint: disable-msg=E0611
# from scipy.optimize import curve_fit
from scipy.optimize import leastsq

from datetime import datetime

class Model():
    def __init__(self, p0, qm, h, rho, nu, C, r, t, data=None):
        self.p0 = p0 # Initial pressure
        self.qm = qm # Mass flowrate (constant, -ve for production)
        self.h = h # Thickness
        self.rho = rho # Density
        self.nu = nu # Kinematic viscosity
        self.C = C # Compressibility
        self.r = r # Radius
        self.t = t # 1D array of measurement times
        self.data = data # Observed pressure/time measurements

    def model(self, phi, k):
        pass

    def residual_function(self, parameters):
        """
        Calculate the residual of the model (difference between observed data and estimated data)
        """
        phi, k = parameters
        return self.data - self.model(phi, k)

    def find_model_parameters(self, phi=0.1, k=1e-14, curve_fit=False):
        if curve_fit:
            # was going to possibly use scipy.optimize.curve_fit if no initial parameters were given
            pass
        else:
            initial_parameters = np.array([phi, k])
            optimal_parameters, flag = leastsq(self.residual_function, initial_parameters)
            phi, k = optimal_parameters
        return phi, k

    def generate_data(self, phi, k, noise = False, sd = 2.5e-4, save_file=False, filename="example_datafile.txt"):
        """
        Generate approximated data using the Theis solution for a guess of porosity and permeability.
        """
        # filename="datafile_{}.txt".format(datetime.now().strftime("%d-%M-%Y_%H:%M") # Argument
        p = self.model(phi, k)
        if noise:
            np.random.seed(0) # Set random seed to 0 for consistancy in testing
            p += p*sd*np.random.randn(p.shape[0])

        if save_file:
            self.__generate_datafile(filename, p)
        self.data = p
        return self.data
    
    def __generate_datafile(self, filename, measurement):
        with open(filename, 'w') as file:
            file.write('Time(s),Pressure(Pa)\n')
            for i in range(len(self.t)):
                if self.t[i] >= 1e-7:
                    file.write('{},{}\n'.format(self.t[i], measurement[i]))
                else:
                    file.write('{},{}\n'.format(0, measurement[i]))


class Theis_Solution(Model):
    def model(self, phi, k):
        D = k/(self.nu*phi*self.rho*self.C)
        # u = np.divide(r**2,4*D*t, out=np.zeros_like(t), where = t != 0)
        with np.errstate(divide="ignore", invalid="ignore"): # Hides 'RuntimeWarning: invalid value encountered in divide' if t[0] == 0.
            p = self.p0 + (self.qm/(self.r*np.pi*k*(self.h/self.nu)))*exp1((self.r**2)/(4*D*self.t)) # Double check exponential integral is correct
            # p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*t)) # Double check exponential integral is correct
        if self.t[0] <= 1e-7: # Check if initial reading is at time 0
            p[0] = self.p0
        return p