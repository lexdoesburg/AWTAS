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
        if phi >= 0.22 or phi <= 0.0095:
            self.data.approximation = [9999999999]*len(self.data.observation)
            return [9999999999]*len(self.data.observation)
        elif k >= 1.2e-12 or k <= 0.8e-16:
            self.data.approximation = [9999999999]*len(self.data.observation)
            return [9999999999]*len(self.data.observation)

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
        # phi, k = self.transf1(variables)
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
        else:
            chi_squared = np.sum(self.data.observation-self.data.approximation)**2
        return chi_squared

    def generate_initial_guess(self, initial_guess=None, single_run=False):
        k_range = np.array([1e-16, 1e-12]) # Range of feasible permeabilities
        phi_range = (0.01, 0.2) # Range of feasible porosities
        
        if initial_guess is not None and not single_run:
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
        elif initial_guess is not None and single_run:
            phi_guess = [initial_guess[0]]
            k_guess = [initial_guess[1]]
        else:
            # print('No initial guess supplied')
            k_guess = [1e-16, 1e-15, 1e-14, 1e-13, 1e-12]
            phi_guess = [0.01, 0.105, 0.2]
            # phi = 0.105

        best_estimate = np.empty(3) # Format, phi, k, chi_sq
        best_estimate.fill(np.inf)
        iterat = 0
        for k in k_guess:
            for phi in phi_guess:
                iterat += 1
                initial_parameters = [phi, k]
                # initial_parameters = self.transf0(initial_parameters)
                # optimal_parameters, flag = leastsq(self.residual_function, initial_parameters, args=(self.data.error), diag=(1./1e-2,1./1e-15), epsfcn=1.7) # if flag is 1 - found a good soln
                optimal_parameters, flag = leastsq(self.residual_function, initial_parameters, args=(self.data.error), epsfcn=1.7) # if flag is 1 - found a good soln
                # optimal_parameters = self.transf1(optimal_parameters)
                chi_squared = self.__chi_squared(self.data.error)
                print('Iter {}: Phi: {} k: {} Chi squared: {}'.format(iterat, phi,k,chi_squared))
                # chi_squared_magnitude = np.floor(np.log10(chi_squared))
                # truth_test = parameters_in_range(optimal_parameters)
                # print('Phi: {} k: {} Chi squared: {} Truth test: {}'.format(phi,k,chi_squared,truth_test))
                if chi_squared < best_estimate[2]:
                    best_estimate[0], best_estimate[1], best_estimate[2] = optimal_parameters[0], optimal_parameters[1], chi_squared
                # print('Current Best Estimate: {}'.format(best_estimate))
                with open('simulation_parameters.txt', 'a') as file:
                    file.write('Current estimate: {}\n'.format(best_estimate))
                    
        return best_estimate[:-1]

    def find_model_parameters2(self, initial_guess=None, verbose=False, single_run=False):
        start_time = time.clock()
        # initial_guess = self.transf0(initial_guess)
        # print('Initial guess = {}'.format(initial_guess))
        initial_parameters = self.generate_initial_guess(initial_guess, single_run)
        if not single_run:
            # optimal_parameters, flag = leastsq(self.residual_function, initial_parameters, args=(self.data.error), diag=(1./1e-2,1./1e-15), epsfcn=1.7) # if flag is 1 - found a good soln            
            optimal_parameters, flag = leastsq(self.residual_function, initial_parameters, args=(self.data.error),epsfcn=1.7) # if flag is 1 - found a good soln
        else:
            optimal_parameters = initial_parameters
        chi_squared = self.__chi_squared(self.data.error)
        self.data.set_unknown_parameters(optimal_parameters)
        end_time = time.clock()
        if verbose:
            print('Total model calls: {} Total time spent: {}'.format(self.calls, (end_time-start_time)))
            print('Initial Guess: {}'.format(initial_parameters))
            print('Optimal Phi: {} Optimal k: {} Chi-squared: {}'.format(optimal_parameters[0], optimal_parameters[1], chi_squared))
        self.calls = 0
        return optimal_parameters

    def transf0(self, theta):
        """transform parameter vector to a new base adapted for automated calibration

        Args:
            theta (list): parameter vector in original base
        """
        xlim = [0.01, 0.2]	# initial phi range
        ylim = [1e-16, 1e-12]	# initial k range
        return [(theta[0]-xlim[0])/(xlim[1]-xlim[0]), (theta[1]-ylim[0])/(ylim[1]-ylim[0])]


    def transf1(self, X):
        """transform parameter vector to their original base

        Args:
            X (list): parameter vector in adapted base
        """
        xlim = [0.01, 0.2]	# initial phi range
        ylim = [1e-16, 1e-12]	# initial k range
        return [X[0]*(xlim[1]-xlim[0])+xlim[0], X[1]*(ylim[1]-ylim[0])+ylim[0]]

   
    def generate_data(self, variables, time, parameters=None, noise = False, sd = 150, save_file=False, filename="{}_testdata.dat".format('theis_soln')):
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

import os
from t2data import *
from t2listing import *
import csv

class SKG9D(Model):
    # parameter space
        
    def transf1(self, X):
        """transform parameter vector to their original base

        Args:
            X (list): parameter vector in adapted base
        """
        xlim = [0.01, 0.2]	# initial phi range
        ylim = [1e-16, 1e-12]	# initial k range
        return [X[0]*(xlim[1]-xlim[0])+xlim[0], X[1]*(ylim[1]-ylim[0])+ylim[0]]


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
        # porosity, permeability = self.transf1(parameters)
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

        if list_p == []:
            list_p = [0]*123

        start_deletion = time.time()
        # Delete negative pressures
        while list_p[-1] < 0.:
            # list_t = list_t[:-2]
            # list_h = list_h[:-2]
            list_p = list_p[:-1]
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