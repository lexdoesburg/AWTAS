import numpy as np
from scipy.special import exp1
# from scipy.optimize import curve_fit
from scipy.optimize import leastsq

import data as data_class
import time

class Model():
    def __init__(self, data=None):
        """
        Initialise the model
        """
        self.function_calls = 0
        self.data = data # Data structure containing observed pressure measurements, time of measurements and well parameters

    def model(self, parameters):
        """
        Model function.
        """
        raise NotImplementedError('Model function has not been implemented.')

    def residual_function(self, parameters):
        """
        Calculate the residual of the model (difference between observed data and estimated data)
        """
        # phi, k = parameters
        self.function_calls += 1
        print('--------\nCalled residual function {} time(s)\n--------'.format(self.function_calls))
        print('Parameters = {}'.format(parameters))
        model_data = self.model(parameters)
        if len(self.data.observation) != len(model_data):
            print('Model finished early - return large residual of 100000\n--------')
            return [100000]*len(self.data.observation)
        else:
            print('Current chi_squared = {}\nCurrent sum of residuals = {}\n--------'.format(self.__chi_squared(), np.sum(self.data.observation - model_data)))

        return self.data.observation - model_data

    def __chi_squared(self):
        sd = np.std(self.data.observation)
        chi_squared = np.sum((self.data.observation-self.data.approximation)/sd)**2
        return chi_squared

    def find_model_parameters(self, phi=None, k=None, curve_fit=False, initial_parameters=None):
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
            # estimates = np.ndarray(shape=(3,25))
            # i = 0
            # for k in [1e-16, 1e-15, 1e-14, 1e-13, 1e-12]:
            #     for phi in [0.2, 0.15, 0.1, 0.05, 0.01]:
                    # initial_parameters = np.array([phi, k])
            for p0 in np.linspace(100, 105, 3):
                for x0 in np.linspace(0.225, 0.25, 3):
                    # initial_parameters = np.array([p0, x0, 0.08, 2.8e-15])
                    calls += 1
                    print('---Calling Optimisation Function {}---'.format(calls))
                    initial_parameters = np.array([p0, x0])
                    optimal_parameters, flag = leastsq(self.residual_function, initial_parameters) # if flag is 1 - found a good soln
                    print('Current Optimal Parameters = {}'.format(optimal_parameters))
                    # self.data.set_approximation(self.model(optimal_parameters)) # Store the approximated data in the data structure
                    # chi_squared = self.__chi_squared()
                    # estimates[:, i] = optimal_parameters[0], optimal_parameters[1], chi_squared
                    # print('Phi: {} k: {} Chi squared: {}'.format(phi, k, chi_squared))
                    # i += 1
                #     if chi_squared <= 20:
                #         broken = True
                #         break
                # if broken:
                #     break
            # index = np.argmin(estimates[2])
            # print('index = {} chi squared = {}'.format(index, estimates[2, index]))
            # print('Function called: {} times'.format(calls))

            # parameters = estimates[:2, index]
            # phi, k = parameters
            # print('Phi {}, k {}'.format(parameters[0], parameters[1]))
            # phi, k = optimal_parameters
        # self.data.set_unknown_parameters(phi, k) # Store phi and k in data structure
        self.data.set_approximation(self.model(optimal_parameters)) # Store the approximated data in the data structure
        print('Function calls = {}'.format(self.function_calls))
        self.function_calls = 0
        return optimal_parameters

    def generate_data(self, phi, k, time, parameters, noise = False, sd = 2.5e-4, save_file=False, filename="example_datafile.txt"):
        """
        Generate approximated data using the Theis solution for a guess of porosity and permeability.
        """
        # filename="datafile_{}.txt".format(datetime.now().strftime("%d-%M-%Y_%H:%M") # Argument
        self.data = data_class.Data()
        self.data.set_time(time)
        self.data.set_known_parameters(parameters)
        p = self.model([phi, k])
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
    def model(self, parameters):
        """
        """
        # Parameter space

        # Parameters to determine
        p0 = parameters[0]*1.e5 # Initial reservoir pressure (Pa)
        x0 = parameters[1]		# initial steam mass fraction
        porosity = .082		# initial porosity
        permeability = 2.76e-15	# initial permeability (m2)
        # porosity = parameters[2]
        # permeability = parameters[3]

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

        start_deletion = time.time()
        # Delete negative pressures
        while list_p[-1] < 0.:
            # list_t = list_t[:-2]
            # list_h = list_h[:-2]
            list_p = list_p[:-2]
        self.data.set_approximation(list_p)
        end_deletion = time.time()
        deletion_time = end_deletion-start_deletion

        print('AUTOUGH Model Time:')
        print('    Datfile time = {}'.format(datfile_time))
        print('    Simulation time = {}'.format(simulator_time))
        print('    Csv time = {}'.format(csv_time))
        print('    Listing time = {}'.format(listing_time))
        print('    Deletion time = {}'.format(deletion_time))
        return list_p

    def transf0(self, theta):
        """transform parameter vector to a new base adapted for automated calibration

        Args:
            theta (list): parameter vector in original base
        """
        xlim = [90., 125.]	# initial pressure range
        ylim = [0.05, .5]	# inital steam fraction range
        return [(theta[0]-xlim[0])/(xlim[1]-xlim[0]), (theta[1]-ylim[0])/(ylim[1]-ylim[0])]
	
    def transf1(self, X):
        """transform parameter vector to their original base

        Args:
            X (list): parameter vector in adapted base
        """
        xlim = [90., 125.]	# initial pressure range
        ylim = [0.05, .5]	# inital steam fraction range
        return [X[0]*(xlim[1]-xlim[0])+xlim[0], X[1]*(ylim[1]-ylim[0])+ylim[0]]

    def wellbore_obj(self, X):
        """Returns the objective function associated to the wellbore model

        Args:
            X (list): parameter vector in a base adapted for automated calibration
        """
        
        # runs the AUTOUGH2 model
        list_t, list_h, list_p = self.model(self.transf1(X))
        
        # load production data
        t_data = []						# initiate time vector
        h_data = []						# initiate enthalpy vector
        with open('SKG9D_enth.csv', 'r') as csvfile:	# open enthalpy production history csv file
            reader=csv.reader(csvfile, delimiter=',')	# open reader
            for row in reader:						# go through each row
                t_data.append(float(row[0]))			# save time
                h_data.append(float(row[1]))			# save enthalpy
        p_data = []						# initiate pressure vector
        with open('SKG9D_press.csv', 'r') as csvfile:	# open pressure production history csv file
            reader=csv.reader(csvfile, delimiter=',')	# open reader
            for row in reader:						# go through each row								
                p_data.append(float(row[1]))			# save pressure
        
        # calculate the objective function
        obj = 0.	# initiate
        w = .046	# enthalpy weight
        for i in range(len(list_t)):	# runs through time list
            obj += ((list_h[i]-h_data[i])*w)**2		# add weighted sqared enthalpy difference
        for i in range(len(list_p)):	# runs through time list
            obj += (list_p[i]-p_data[i])**2			# add sqared pressure difference
        
        # add penalty for models not running to the end
        n_diff = len(t_data)-len(list_t)	# numbre of time dots not reached by the model
        obj += n_diff*800.					# penalty
        
        # reduce objective function value for practical reasons
        obj = obj*1.e-5
        return obj