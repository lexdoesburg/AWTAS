import numpy as np
from scipy.special import exp1
from scipy.optimize import curve_fit
from scipy.optimize import leastsq

import lmfit

import data2 as data_class

import time

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
        phi, k = variables
        self.calls += 1
        start = time.clock()
        approximation = self.model(variables)
        if len(approximation) != len(self.data.observation):
            num_to_append = len(self.data.observation)-len(approximation)
            values = np.zeros(num_to_append)
            # values.fill(np.inf)
            approximation = np.append(approximation, values)
            print('Appending {} zeros to the approximation'.format(num_to_append))
        if error is None:
            # Calling without estimated error   
            # print("Call without estimated error")
            residual = self.data.observation - approximation
        else:
            # Calling with estimated error   
            # print("Call with estimated error")   
            residual = (self.data.observation - approximation)/error
        end = time.clock()
        print('Calling residual ({}): Phi: {} k: {} Chi-squared: {} Time Elapsed: {}'.format(self.calls, phi, k, self.__chi_squared(error), end-start))
        return residual

    def __chi_squared(self, error=None):
        # sd = np.std(self.data.observation)
        if len(self.data.approximation) != len(self.data.observation):
            num_to_append = len(self.data.observation)-len(self.data.approximation)
            values = np.zeros(num_to_append)
            # values.fill(np.inf)
            self.data.approximation = np.append(self.data.approximation, values)
            print('Appending {} zeros to the approximation'.format(num_to_append))

        if error is not None:
            chi_squared = np.sum((self.data.observation-self.data.approximation)/error)**2
            # chi_squared = ((self.data.observation-self.data.approximation)/error)**2
            # chi_squared = np.sum(((self.data.observation-self.data.approximation)/error)**2)
        else:
            chi_squared = np.sum(self.data.observation-self.data.approximation)**2
        return chi_squared

    def func(self, t, phi, k):
        variables = phi, k
        return self.model(variables)


    def find_model_params_test(self):
        xData = self.data.time
        yData = self.data.observation
        initial_guess = [0.105, 1e-14]
        popt, pcov = curve_fit(self.func, xData, yData, p0=initial_guess, diag=(1./xData.mean(),1./yData.mean()) )
        print(popt)
        return popt

    def generate_initial_guess(self, initial_guess=None):
        k_range = np.array([1e-16, 1e-12]) # Range of feasible permeabilities
        phi_range = (0.01, 0.2) # Range of feasible porosities
        
        # def parameters_in_range(parameters):
        #     phi, k = parameters
        #     range_k = k_range/1e-16
        #     k = k/1e-16
        #     eps = 1e-7
        #     if phi_range[0] <= phi <= phi_range[1] and k-range_k[0] >= eps and range_k[1]-k >= eps:
        #         return True
        #     else:
        #         return False

        if initial_guess is not None:
            # print('Initial guess supplied')
            phi, k = initial_guess
            k_magnitude = np.floor(np.log10(k))
            print('k magnitude = ', k_magnitude)
            if k_magnitude-1 < -16:
                k_guess = [k_range[0], k, 10**(k_magnitude+1)]
            elif k_magnitude+1 > -12:
                k_guess = [10**(k_magnitude-1), k, k_range[1]]
            else:
                k_guess = [10**(k_magnitude-1), k, 10**(k_magnitude+1)]
            
            phi_search_range = 0.04 # Will search 0.07 above and below the given guess
            if phi - phi_search_range < 0.01:
                phi_guess = [phi_range[0], phi, phi+phi_search_range]
            elif phi + phi_search_range > 0.2:
                phi_guess = [phi-phi_search_range, phi, phi_range[1]]
            else:
                phi_guess = [phi-phi_search_range, phi, phi+phi_search_range]
        else:
            # print('No initial guess supplied')
            k_guess = [1e-16, 1e-14, 1e-12]
            phi_guess = [0.01, 0.105, 0.2]
            # phi = 0.105

        best_estimate = np.empty(3) # Format, phi, k, chi_sq
        best_estimate.fill(np.inf)

        for k in k_guess:
            for phi in phi_guess:
                initial_parameters = [phi, k]
                optimal_parameters, flag = leastsq(self.residual_function, initial_parameters, args=(self.data.error), diag=(1./1e-2,1./1e-16)) # if flag is 1 - found a good soln
                chi_squared = self.__chi_squared(self.data.error)
                print('Phi: {} k: {} Chi squared: {}'.format(phi,k,chi_squared))
                # chi_squared_magnitude = np.floor(np.log10(chi_squared))
                # truth_test = parameters_in_range(optimal_parameters)
                # print('Phi: {} k: {} Chi squared: {} Truth test: {}'.format(phi,k,chi_squared,truth_test))
                if chi_squared < best_estimate[2]:
                    best_estimate[0], best_estimate[1], best_estimate[2] = optimal_parameters[0], optimal_parameters[1], chi_squared
                # print('Current Best Estimate: {}'.format(best_estimate))
                with open('simulation_parameters.txt', 'a') as file:
                    file.write('Current estimate: {}\n'.format(best_estimate))
                    
        return best_estimate[:-1]

    def find_model_parameters2(self, initial_guess=None, verbose=False):
        start_time = time.clock()
        initial_parameters = self.generate_initial_guess(initial_guess)
        optimal_parameters, flag = leastsq(self.residual_function, initial_parameters, args=(self.data.error), diag=(1./1e-2,1./1e-16)) # if flag is 1 - found a good soln
        chi_squared = self.__chi_squared(self.data.error)
        self.data.set_unknown_parameters(optimal_parameters)
        end_time = time.clock()
        if verbose:
            print('Total model calls: {} Total time spent: {}'.format(self.calls, (end_time-start_time)))
            print('Initial Guess: {}'.format(initial_parameters))
            print('Optimal Phi: {} Optimal k: {}'.format(optimal_parameters[0], optimal_parameters[1]))
        self.calls = 0
        return optimal_parameters

    def find_model_parameters(self, variables=None, curve_fit=False, verbose=False):
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
        start_time = time.clock()
        # if curve_fit:
        #     # was going to possibly use scipy.optimize.curve_fit if no initial parameters were given
        #     pass
        # else:
        calls = 0
        # broken = False
        estimates = np.ndarray(shape=(3,25))
        i = 0
        # Could run this for the 3-5 values of k with a fixed phi and then take the best solution as starting point
        for k in [1e-16, 1e-15, 1e-14, 1e-13, 1e-12]:
            for phi in [0.2, 0.15, 0.1, 0.05, 0.01]:
                initial_parameters = np.array([phi, k])
                # optimal_parameters, flag = leastsq(self.residual_function, initial_parameters) # if flag is 1 - found a good soln                
                # if self.data.error:
                #     # Calling with estimated error
                #     # print("Call with estimated error")
                optimal_parameters, flag = leastsq(self.residual_function, initial_parameters, args=(self.data.error)) # if flag is 1 - found a good soln
                # else: 
                    # Calling without estimated error   
                    # print("Call without estimated error")            
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
        optimal_parameters = estimates[:2, index]
        self.data.set_unknown_parameters(optimal_parameters) # Store phi and k in data structure
        self.data.set_approximation(self.model(optimal_parameters)) # Store the approximated data in the data structure

        end_time = time.clock()
        if verbose:
            print('Total model calls: {} Total time spent: {}'.format(self.calls, (end_time-start_time)))
            print('Optimal Phi: {} Optimal k: {}'.format(optimal_parameters[0], optimal_parameters[1]))
        self.calls = 0
        return optimal_parameters

    def generate_data(self, variables, parameters, time, noise = False, sd = 150, save_file=False, filename="{}_testdata.dat".format('theis_soln')):
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
            # magnitude = 10**(np.floor(np.log10(np.average(p))))
            noise = np.random.randn(p.shape[0])
            # noise = (sd*(noise/np.std(noise)))*magnitude
            noise = sd * (noise/np.std(noise))
            p += noise
            # p += sd*np.random.randn(p.shape[0])
            print('SD = {}'.format(np.std(noise)))
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
    def model(self, parameters):
        """
        Calculate and return the pressure for the given input using the analytical Theis solution.

        Inputs: phi - Porosity
                k - Permeability
                Makes use of self.data for remaining parameters
        
        Output: p(t) - pressure at a given time (or time array).
        """
        phi, k = parameters
        p0, qm, h, rho, nu, C, r = self.data.parameters # Unpack the well parameters
        D = k/(nu*phi*rho*C)
        with np.errstate(divide="ignore", invalid="ignore"): # Hides 'RuntimeWarning: invalid value encountered in divide' if t[0] == 0.
            p = p0 + (qm/(r*np.pi*k*(h/nu)))*exp1((r**2)/(4*D*self.data.time)) # TODO: Double check exponential integral is correct
            # p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*t)) # TODO: Double check exponential integral is correct
        if self.data.time[0] <= 1e-7: # Check if initial reading is at time 0
            p[0] = p0
        return p

import os
from t2data import *
from t2listing import *
import csv

class SKG9D(Model):
    def model(self, parameters, verbose=False):
        """
        """
        # Parameter space

        # Parameters to determine
        p0 = 103.07e5 # Initial reservoir pressure (Pa)
        x0 = 0.2332		# initial steam mass fraction
        # porosity = .082		# initial porosity
        # permeability = 2.76e-15	# initial permeability (m2)
        # porosity = parameters[2]
        # permeability = parameters[3]
        porosity, permeability = parameters
        
        start_datfile = time.time()
        # write values in the AUTOUGH2 dat file
        dat = t2data('SKG9D.DAT')
        dat.grid.rocktype['IGNIM'].porosity = porosity
        dat.grid.rocktype['IGNIM'].permeability = [permeability, permeability, permeability]
        dat.parameter['default_incons'] = [p0, x0]
        dat.write('SKG9D_1.DAT')
        end_datfile = time.time()
        datfile_time = end_datfile-start_datfile

        start_simulator = time.time()
        # run AUTOUGH2
        dat.run(simulator='AUTOUGH2_5.exe', silent = True)
        end_simulator = time.time()
        simulator_time = end_simulator-start_simulator

        start_csv = time.time()
        # time vector indicated when result must be returned
        t_data = []
        with open('SKG9D_press.csv', 'r') as csvfile:
            reader=csv.reader(csvfile, delimiter=',')
            for row in reader:								
                t_data.append(float(row[0]))
        end_csv = time.time()
        csv_time = end_csv-start_csv
        
        try:
            start_listing = time.time()
            # read LISTING file
            lst = t2listing('SKG9D_1.LISTING', skip_tables = ['connection'])
            [(th, h), (tp, p)] = lst.history([('g', ('  A 1', 'WEL 1'), 'Enthalpy'), ('e', '  A 1', 'Pressure')])
            th = np.array(th)*(1/24.)*(1/3600.)	# convert time from seconds to days
            p = np.array(p)*1.e-5				# convert pressure from Pa to bars
            # h = np.array(h)*1.e-3				# convert enthalpy from J/kg to kJ/kg
            list_index = [(np.abs(th-t)).argmin() for t in t_data if t<=th[-1]]	# index list when result must be returned
            # list_t = th[list_index]
            # list_h = h[list_index]
            list_p = p[list_index]
            end_listing = time.time()
            listing_time = end_listing-start_listing
        except:
            # print('Inside exception')
            list_p = [0]*123

        start_deletion = time.time()
        # Delete negative pressures
        while list_p[-1] < 0.:
            # list_t = list_t[:-2]
            # list_h = list_h[:-2]
            list_p = list_p[:-2]
        self.data.set_approximation(list_p)
        end_deletion = time.time()
        deletion_time = end_deletion-start_deletion

        if verbose:
            print('AUTOUGH Model Time:')
            print('    Datfile time = {}'.format(datfile_time))
            print('    Simulation time = {}'.format(simulator_time))
            print('    Csv time = {}'.format(csv_time))
            print('    Listing time = {}'.format(listing_time))
            print('    Deletion time = {}'.format(deletion_time))

        return list_p
    
    def chi_squared(self, error=None):
        # sd = np.std(self.data.observation)
        if len(self.data.approximation) != len(self.data.observation):
            num_to_append = len(self.data.observation)-len(self.data.approximation)
            values = np.zeros(num_to_append)
            # values.fill(np.inf)
            self.data.approximation = np.append(self.data.approximation, values)
            print('Appending {} zeros to the approximation'.format(num_to_append))

        if error is not None:
            chi_squared = np.sum((self.data.observation-self.data.approximation)/error)**2
            # chi_squared = ((self.data.observation-self.data.approximation)/error)**2
            # chi_squared = np.sum(((self.data.observation-self.data.approximation)/error)**2)
        else:
            chi_squared = np.sum(self.data.observation-self.data.approximation)**2
        return chi_squared