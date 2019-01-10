import numpy as np
from scipy.special import exp1
from scipy.optimize import curve_fit
from scipy.optimize import leastsq

# import lmfit

import data as data_class

import time

class Model():
    def __init__(self, data=None, model_type=None):
        """
        Initialise the model
        """
        self.data = data # Data structure containing observed pressure measurements, time of measurements and well parameters
        self.model_type = model_type
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
        if error is None:
            # Calling without estimated error   
            # print("Call without estimated error")   
            residual = self.data.observation - self.model(variables)
        else:
            # Calling with estimated error   
            # print("Call with estimated error")   
            residual = (self.data.observation - self.model(variables))/error
        return residual

    def __chi_squared(self, error=None):
        # sd = np.std(self.data.observation)
        if error is not None:
            chi_squared = np.sum((self.data.observation-self.data.approximation)/error)**2
            # chi_squared = ((self.data.observation-self.data.approximation)/error)**2
            # chi_squared = np.sum(((self.data.observation-self.data.approximation)/error)**2)
        else:
            chi_squared = np.sum(self.data.observation-self.data.approximation)**2
        return chi_squared

    # def func(self, t, phi, k):
    #     variables = phi, k
    #     return self.model(variables)


    # def find_model_params_test(self):
    #     xData = self.data.time
    #     yData = self.data.observation
    #     initial_guess = [0.105, 1e-14]
    #     popt, pcov = curve_fit(self.func, xData, yData, p0=initial_guess, diag=(1./xData.mean(),1./yData.mean()) )
    #     print(popt)
    #     return popt

    def generate_initial_guess(self, initial_guess=None, single_run=False):
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
            
            phi_search_range = 0.07 # Will search 0.07 above and below the given guess
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
            # Guess' for Theis solution
            # k_guess = [1e-16, 1e-15, 1e-14, 1e-13, 1e-12]
            # phi_guess = np.linspace(0.01, 0.2, 10)
            # Guess' for Radial1d
            k_guess = [1e-15, 1e-14, 1e-13, 1e-12]
            phi_guess = np.linspace(0.01, 0.14, 10)
            # phi = 0.105

        best_estimate = np.empty(3) # Format, phi, k, chi_sq
        best_estimate.fill(np.inf)

        for k in k_guess:
            for phi in phi_guess:
                initial_parameters = [phi, k]
                # optimal_parameters, flag = leastsq(self.residual_function, initial_parameters, args=(self.data.error), epsfcn=None) # if flag is 1 - found a good soln
                self.data.approximation = self.model(initial_parameters)
                chi_squared = self.__chi_squared(self.data.error)
                print('Phi: {} k: {} Chi squared: {}'.format(phi,k,chi_squared))
                # chi_squared_magnitude = np.floor(np.log10(chi_squared))
                # truth_test = parameters_in_range(optimal_parameters)
                # print('Phi: {} k: {} Chi squared: {} Truth test: {}'.format(phi,k,chi_squared,truth_test))
                if chi_squared < best_estimate[2]:
                    # best_estimate[0], best_estimate[1], best_estimate[2] = optimal_parameters[0], optimal_parameters[1], chi_squared
                    best_estimate[0], best_estimate[1], best_estimate[2] = phi, k, chi_squared
                    
                # print('Current Best Estimate: {}'.format(best_estimate))
                    
        return best_estimate[:-1]

    def find_model_parameters2(self, initial_guess=None, verbose=False, single_run=False):
        start_time = time.clock()
        initial_parameters = self.generate_initial_guess(initial_guess, single_run=single_run)
        if not single_run:
            optimal_parameters, flag = leastsq(self.residual_function, initial_parameters, args=(self.data.error), epsfcn=None) # if flag is 1 - found a good soln
        else:
            optimal_parameters = initial_parameters
        self.data.set_unknown_parameters(optimal_parameters)
        self.data.set_approximation(self.model(optimal_parameters)) # Store the approximated data in the data structure
        end_time = time.clock()
        if verbose:
            chi_squared = self.__chi_squared(self.data.error)
            print('Total model calls: {} Total time spent: {}'.format(self.calls, (end_time-start_time)))
            print('Initial Guess: {}'.format(initial_parameters))
            print('Optimal Phi: {} Optimal k: {} Chi-Squared: {}'.format(optimal_parameters[0], optimal_parameters[1], chi_squared))
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
        self.data = data_class.Data(model_type=self.model_type)
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


class Theis_Solution(Model):
    def __init__(self, data=None):
        """
        Initialise the model
        """
        super().__init__(data=data, model_type='theis')

    def model(self, variables):
        """
        Calculate and return the pressure for the given input using the analytical Theis solution.

        Inputs: phi - Porosity
                k - Permeability
                Makes use of self.data for remaining parameters
        
        Output: p(t) - pressure at a given time (or time array).
        """
        phi, k = variables
        # p0, qm, h, rho, nu, C, r = self.data.parameters # Unpack the well parameters
        # Unpack the well parameters
        p0 = self.data.parameters['Initial Pressure']['Value']
        qm = self.data.parameters['Mass Flowrate']['Value']
        h = self.data.parameters['Layer Thickness']['Value']
        rho = self.data.parameters['Density']['Value']
        nu = self.data.parameters['Kinematic Viscosity']['Value']
        C = self.data.parameters['Compressibility']['Value']
        r = self.data.parameters['Radius']['Value']

        D = k/(nu*phi*rho*C)
        with np.errstate(divide="ignore", invalid="ignore"): # Hides 'RuntimeWarning: invalid value encountered in divide' if t[0] == 0.
            # p = p0 + (qm/(r*np.pi*k*(h/nu)))*exp1((r**2)/(4*D*self.data.time)) 
            p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*self.data.time)) # Same as fortran output
        if self.data.time[0] <= 1e-7: # Check if initial reading is at time 0
            p[0] = p0
        self.data.approximation = p # Set the approximation for quick estimation of chi_sq
        return p

    # def test_model(self, params, t):
    #     phi = params['phi']
    #     k = params['k']
    #     p0 = params['p0']
    #     qm = params['qm']
    #     rho = params['rho']
    #     nu = params['nu']
    #     h = params['h']
    #     C = params['C']
    #     r = params['r']

    #     D = k/(nu*phi*rho*C)
    #     with np.errstate(divide="ignore", invalid="ignore"): # Hides 'RuntimeWarning: invalid value encountered in divide' if t[0] == 0.
    #         # p = p0 + (qm/(r*np.pi*k*(h/nu)))*exp1((r**2)/(4*D*self.data.time)) 
    #         p = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*t)) # Same as fortran output
    #     if t[0] <= 1e-7: # Check if initial reading is at time 0
    #         p[0] = p0
    #     return p

# import theis_wrapper as theis_fortran

# class Theis_Solution_Fortran(Model):
#     def model(self, parameters):
#         p0, qm, h, rho, nu, C, r = self.data.parameters
#         phi, k = parameters
#         num_observations = len(self.data.time)
#         p = theis_fortran.theis(k, nu, phi, rho, C, h, qm, p0, r, num_observations, self.data.time)
#         return p

from radial1d_wrapper import radial1d

class Radial_1D(Model):
    def __init__(self, data=None):
        """
        Initialise the model
        """
        super().__init__(data=data, model_type='radial1d')

    def model(self, parameters):
        # Unpack the parameters which we are trying to find
        phi, k = parameters
        print('Phi = {} k = {}'.format(phi,k))
        # Get relevant data from data structure
        # p0, X0, rw, thick, CR, COND, RHOR, COMP, ConstRate, distFromWell = self.data.parameters
        p0 = self.data.parameters['Initial Pressure']['Value']
        X0 = self.data.parameters['Initial Temperature']['Value']
        rw = self.data.parameters['Action Well Radius']['Value']
        thick = self.data.parameters['Layer Thickness']['Value']
        CR = self.data.parameters['Rock Specific Heat']['Value']
        COND = self.data.parameters['Rock Conductivity']['Value']
        RHOR = self.data.parameters['Rock Density']['Value']
        COMP = self.data.parameters['Rock Compressibility']['Value']
        ConstRate = self.data.parameters['Mass Flowrate']['Value']
        distFromWell = self.data.parameters['Observation Point Distance']['Value']
        # Run the model to get result
        print('Calling radial1d')
        pressure = radial1d(phi, k, p0, X0, rw, thick, CR, COND, RHOR, COMP, ConstRate, distFromWell, 271, np.linspace(0, 54000, 271))
        print('Returning pressure')
        return pressure