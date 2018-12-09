import numpy as np
from scipy.special import exp1
from scipy.optimize import curve_fit
from scipy.optimize import leastsq

import lmfit

import data as data_class

class Model():
    def __init__(self, data=None):
        """
        Initialise the model
        """
        self.data = data # Data structure containing observed pressure measurements, time of measurements and well parameters
        self.calls = 0

    def model(self, variables):
        """
        Model function.

        Input: variables (list or np.array): List or numpy array of variable parameters for the model.

        Output: modelled pressures at reading times.
        """
        raise NotImplementedError('Model function has not been implemented.')

    def residual_function(self, variables, error=None):
        """
        Calculate the residual of the model (difference between observed data and estimated data)
        """
        # phi, k = parameters
        self.calls += 1
        if not error:
            residual = self.data.observation - self.model(variables)
        else:
            residual = self.__chi_squared(error)
        return residual

    def __chi_squared(self, error=None):
        # sd = np.std(self.data.observation)
        if error:
            chi_squared = ((self.data.observation-self.data.approximation)/error)**2
            # chi_squared = np.sum(((self.data.observation-self.data.approximation)/error)**2)
        else:
            chi_squared = np.sum((self.data.observation-self.data.approximation))**2
        return chi_squared

    def find_model_params_test(self):
        pass

    def find_model_parameters(self, variables=None, curve_fit=False):
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

        # if curve_fit:
        #     # was going to possibly use scipy.optimize.curve_fit if no initial parameters were given
        #     pass
        # else:
        calls = 0
        # broken = False
        estimates = np.ndarray(shape=(3,25))
        i = 0
        for k in [1e-16, 1e-15, 1e-14, 1e-13, 1e-12]:
            for phi in [0.2, 0.15, 0.1, 0.05, 0.01]:
                initial_parameters = np.array([phi, k])
                optimal_parameters, flag = leastsq(self.residual_function, initial_parameters) # if flag is 1 - found a good soln                
                # optimal_parameters, cov_x, infodict, mesg, flag = leastsq(self.residual_function, initial_parameters, full_output=1) # if flag is 1 - found a good soln
                self.data.set_approximation(self.model(optimal_parameters)) # Store the approximated data in the data structure
                chi_squared = self.__chi_squared()
                # print('Nfev = ', infodict['nfev'])
                # print('Chi squared: {} Function evaluated at output: {}'.format(chi_squared, np.sum(infodict['fvec'])))
                estimates[:, i] = optimal_parameters[0], optimal_parameters[1], chi_squared                
                # estimates[:, i] = optimal_parameters[0], optimal_parameters[1], abs(np.sum(infodict['fvec']))
                # print('Phi: {} k: {} Chi squared: {}'.format(phi, k, chi_squared))
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

        optimal_parameters = estimates[:2, index]
        print('Optimal phi: {} Optimal k: {}'.format(optimal_parameters[0], optimal_parameters[1]))
        # print('Phi {}, k {}'.format(phi, k))
        # phi, k = optimal_parameters
        self.data.set_unknown_parameters(optimal_parameters) # Store phi and k in data structure
        self.data.set_approximation(self.model(optimal_parameters)) # Store the approximated data in the data structure
        print('Calls = ', self.calls)
        self.calls = 0
        return optimal_parameters

    def generate_data(self, variables, parameters, time, noise = False, sd = 2.5e-4, save_file=False, filename="{}_testdata.dat".format('theis_soln')):
        """
        Generate approximated data using the Theis solution for a guess of porosity and permeability.
        """
        # filename="datafile_{}.txt".format(datetime.now().strftime("%d-%M-%Y_%H:%M") # Argument
        self.data = data_class.Data()
        self.data.set_time(time)
        self.data.set_known_parameters(parameters)
        p = self.model(variables)
        if noise:
            np.random.seed(0) # Set random seed to 0 for consistency in testing
            p += sd*np.random.randn(p.shape[0])
            print('SD = {}'.format(np.std(p)))
        self.data.set_observation(p)
        if save_file:
            self.data.generate_datafile(filename, variables=variables)
        # return self.data
    
    # def __generate_datafile(self, filename, measurement):
    #     with open(filename, 'w') as file:
    #         file.write('Time(s),Pressure(Pa)\n')
    #         for i in range(len(self.data.time)):
    #             if self.data.time[i] >= 1e-7:
    #                 file.write('{},{}\n'.format(self.data.time[i], measurement[i]))
    #             else:
    #                 file.write('{},{}\n'.format(0, measurement[i]))


class Theis_Solution(Model):
    def model(self, variables):
        """
        Calculate and return the pressure for the given input using the analytical Theis solution.

        Inputs: phi - Porosity
                k - Permeability
                Makes use of self.data for remaining parameters
        
        Output: p(t) - pressure at a given time (or time array).
        """
        phi, k = variables
        p0, qm, h, rho, nu, C, r = self.data.parameters # Unpack the well parameters
        D = k/(nu*phi*rho*C)
        with np.errstate(divide="ignore", invalid="ignore"): # Hides 'RuntimeWarning: invalid value encountered in divide' if t[0] == 0.
            # p = p0 + (qm/(r*np.pi*k*(h/nu)))*exp1((r**2)/(4*D*self.data.time)) 
            p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*self.data.time)) # Same as fortran output
        if self.data.time[0] <= 1e-7: # Check if initial reading is at time 0
            p[0] = p0
        return p

# import theis_solution_fortran as ts

class Theis_Solution_Fortran(Model):
    def model(self, parameters):
        p0, qm, h, rho, nu, C, r = self.data.parameters
        phi, k = parameters
        # num_observations = len(self.data.time)
        # p = ts.theis_solution(k, phi, p0, qm, h, rho, nu, C, r, num_observations, self.data.time)
        # return p
        pass

# import NumericalSimulator1D as radial_1D

class Radial_1D(Model):
    def model(self, parameters):
        # Unpack the parameters which we are trying to find

        # Get relevant data from data structure
        
        # Run the model to get result
        # p = radial_1D.NumericalSolution1D(inputs)
        # return p
        pass

class Test_Model(Model):
    def model(self, parameters):
        phi,k,p0,x0 = parameters
        p1, qm, h, rho, nu, C, r = self.data.parameters
        G = qm*h*2.34/(nu*C)
        H = phi*k/h*rho
        I = x0/C + k*r
        p = p0 + G*H*I*self.data.time
        return p
    
    def curve_fit_func(self, time, p0, x0, k, phi):
        # parameters = [phi, k, p0, x0]
        p1, qm, h, rho, nu, C, r = self.data.parameters
        G = qm*h*2.34/(nu*C)
        H = phi*k/h*rho
        I = x0/C + k*r
        p = p0 + G*H*I*self.data.time
        return p


    def find_model_parameters(self, variables=None, curvefit=False):
        if curvefit:
            optimal_parameters, pcov = curve_fit(self.curve_fit_func, self.data.observation, self.data.time, p0=variables, bounds=(0,[0.2,1e-12,5e6,1]))
        else:
            calls = 0
            estimates = np.ndarray(shape=(5,625))
            i = 0
            if variables:
                optimal_parameters, flag = leastsq(self.residual_function, variables)
            else:
                for p0 in np.linspace(1e6,4e6,5):
                    for x0 in np.linspace(0.05,0.3,5):
                        for k in [1e-16, 1e-15, 1e-14, 1e-13, 1e-12]:
                            for phi in [0.2, 0.15, 0.1, 0.05, 0.01]:
                                initial_parameters = np.array([phi, k, p0, x0])
                                optimal_parameters, flag = leastsq(self.residual_function, initial_parameters) # if flag is 1 - found a good soln
                                self.data.set_approximation(self.model(optimal_parameters)) # Store the approximated data in the data structure
                                chi_squared = abs(self.__chi_squared())
                                estimates[:, i] = optimal_parameters[0], optimal_parameters[1], optimal_parameters[2], optimal_parameters[3], chi_squared
                                # estimates[:, i] = optimal_parameters, chi_squared
                                # print('Phi: {} k: {} Chi squared: {}'.format(phi, k, chi_squared))
                                i += 1
                                calls += 1
                #     if chi_squared <= 20:
                #         broken = True
                #         break
                # if broken:
                #     break
            index = np.argmin(estimates[4])
            print('index = {} chi squared = {}'.format(index, estimates[4, index]))
            print('Function called: {} times'.format(calls))

            optimal_parameters = estimates[:4, index]
            # print('Phi {}, k {}'.format(phi, k))
            print('Calls = ', self.calls)
            self.calls = 0
        # phi, k = optimal_parameters
        # self.data.set_unknown_parameters(phi, k, ) # Store phi and k in data structure
        self.data.set_approximation(self.model(optimal_parameters)) # Store the approximated data in the data structure
        return optimal_parameters
    
    def __chi_squared(self):
        sd = np.std(self.data.observation)
        chi_squared = np.sum((self.data.observation-self.data.approximation)/sd)**2
        return chi_squared