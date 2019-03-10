import numpy as np
from scipy.special import exp1
from scipy.optimize import leastsq
# from scipy.optimize import curve_fit
# import lmfit
import awtas.logic.data as data_class


class Model():
    """
    Parent class which defines the inverse modelling functions for the use of determining unknown
    well parameters in well test analysis. Sub classes which inherit from this class will define
    their own forward model which will be used by the inverse modelling code.
    """
    def __init__(self, data=None, model_type=None):
        """
        Initialise the model object.
        """
        self.data = data
        self.model_type = model_type
        self.calls = 0


    def model(self, variables):
        """
        The forward model - calculates the modelled approximation of the well response for a given
        set of input variables.

        Inputs: 
            variables (iterable) - List or numpy array of variable parameters for the model.
        Outputs: 
            modelled_value (np.array) - Modelled approximation of well response at the times specified 
                                        in the data structure.
        """
        raise NotImplementedError('Model function has not been implemented.')


    def residual_function(self, variables, error=None):
        """
        Calculates the residual of the model for a given set of input variables.

        Inputs:
            variables (iterable) - List or numpy array of variable parameters for the model.
        Outputs:
            residual (np.array) - Difference between the observed measurement and the modelled approximation.
        """
        if error is not None:
            # Calculate without estimated error (error is unknown)
            modelled_value = self.model(variables)  
            residual = self.data.observation - modelled_value
        else:
            # Calculate with estimated error (error is known) - residual is weighted by the errors
            modelled_value = self.model(variables)
            residual = (self.data.observation - modelled_value)/error
        # Increment the number of times the forward model has been called
        self.calls += 1
        return residual


    def __chi_squared(self, modelled_value):
        """
        Calculate a measure of the least-squares error, weighted by the standard deviation of 
        any known errors in the observed data.

        Inputs:
            modelled_value (np.array) - The modelled approximation of the well response.
        Outputs:
            chi_squared (float) - Measure of (weighted) least-squares error.
        """
        if self.data.error is not None:
            chi_squared = np.sum(((self.data.observation-modelled_value)/self.data.error)**2)
        else:
            chi_squared = np.sum((self.data.observation-modelled_value)**2)
        return chi_squared


    def generate_initial_guess(self, verbose=False):
        """
        Calls the forward model at discrete points within the feasible range of variable combinations.
        The combination of variables which gives a well response that most closely fits the observed
        data is then returned as the initial guess to be supplied to the inverse modelling code.

        Outputs:
            variables (iterable) - List or numpy array of an initial guess of the variable parameters for
                                   the model.
        """
        # k_range = (1e-16, 1e-12) # Range of feasible permeabilities
        # phi_range = (0.01, 0.2) # Range of feasible porosities
        # if initial_guess is not None:
        #     # Initial guess supplied
        #     phi, k = initial_guess
        #     k_magnitude = np.floor(np.log10(k))
        #     print('k magnitude = ', k_magnitude)
        #     if k_magnitude-1 < -16:
        #         k_guess = [k_range[0], k, 10**(k_magnitude+1)]
        #     elif k_magnitude+1 > -12:
        #         k_guess = [10**(k_magnitude-1), k, k_range[1]]
        #     else:
        #         k_guess = [10**(k_magnitude-1), k, 10**(k_magnitude+1)]
        #     phi_search_range = 0.07 # Will search 0.07 above and below the given guess
        #     if phi - phi_search_range < 0.01:
        #         phi_guess = [phi_range[0], phi, phi+phi_search_range]
        #     elif phi + phi_search_range > 0.2:
        #         phi_guess = [phi-phi_search_range, phi, phi_range[1]]
        #     else:
        #         phi_guess = [phi-phi_search_range, phi, phi+phi_search_range]
        # else:
        #     # No initial guess supplied
        #     k_guess = [1e-16, 1e-15, 1e-14, 1e-13, 1e-12]
        #     phi_guess = np.linspace(0.01, 0.2, 10)

        # Range of variable values to call model over
        k_guess = [1e-16, 1e-15, 1e-14, 1e-13, 1e-12]
        phi_guess = np.linspace(0.01, 0.2, 10)

        # Array to store the best variable combination and its chi_squared value
        best_estimate = np.empty(3) # Format = [phi, k, chi_squared]
        best_estimate.fill(np.inf)

        for k in k_guess:
            for phi in phi_guess:
                initial_guess = [phi, k]
                modelled_value = self.model(initial_guess)
                self.calls += 1
                chi_squared = self.__chi_squared(modelled_value)
                if verbose:
                    print('Current Phi: {}, Current k: {}, Chi-Squared: {}'.format(phi,k,chi_squared))
                # If the new chi_squared value is less than the current best estimate, update the best estimate
                #TODO: Change the condition for best estimate being updated to the chi_squared value which is closest to the total number of observed data points.
                if chi_squared < best_estimate[2]:
                    best_estimate[0], best_estimate[1], best_estimate[2] = phi, k, chi_squared

        # Return the initial guess of variable values                 
        return best_estimate[:-1]


    def find_model_parameters(self, initial_guess=None, verbose=False):
        """
        The inverse model - calls the forward model repeatedly using non-linear least-squares optimisation
        to determine the combination of unknown variable values that give the well response that best fits
        to the observed data.

        Inputs:
            initial_guess (iterable, optional) - A list or np.array of an initial guess of the unknown
                                                 variable values. If this is not supplied, an initial guess
                                                 will be estimated at the cost of extra forward model calls.
            verbose (bool) - If true prints information while the inverse model runs.
        Outpus:
            optimal_parameters (list) - The optimal values of the unknown variables as found by the inverse
                                        model.
        """
        # If initial guess is not supplied then find a suitable initial guess
        if not initial_guess:
            initial_guess = self.generate_initial_guess(verbose)
        
        # Call the inverse modelling code - non-linear least-squares optimisation
        optimal_parameters, flag = leastsq(self.residual_function, initial_guess, args=(self.data.error), epsfcn=None) # if flag is 1 - found a good soln
        self.data.set_unknown_parameters(optimal_parameters) # Store the optimal parameters in the data structure
        self.data.set_approximation(self.model(optimal_parameters)) # Store the approximated data in the data structure
        self.calls += 1

        # Print information about the inverse modelling
        if verbose:
            chi_squared = self.__chi_squared(self.data.approximation)
            print('Initial Phi: {}, Initial k: {}'.format(initial_guess[0], initial_guess[1]))
            print('Optimal Phi: {}, Optimal k: {}'.format(optimal_parameters[0], optimal_parameters[1]))
            print('Total model calls: {}, Chi-Squared: {}'.format(self.calls, chi_squared))
        
        # Reset the number of forward model calls
        self.calls = 0
        return optimal_parameters


    # Old inverse modelling function
    # def find_model_parameters(self, variables=None, curve_fit=False, verbose=False):
        # """
        # Find the model parameters porosity and permeability.

        # Inputs: phi (float) - Intial guess of porosity
        #         k (float) - Initial guess of permeability
        #         curve_fit (boolean)
        
        # Output: phi (float) - Estimated value of porosity
        #         k (float)- Estimated value of permeability
        # """
        # start_time = time.clock()
        # # if curve_fit:
        # #     # was going to possibly use scipy.optimize.curve_fit if no initial parameters were given
        # #     pass
        # # else:
        # calls = 0
        # # broken = False
        # estimates = np.ndarray(shape=(3,25))
        # i = 0
        # # Could run this for the 3-5 values of k with a fixed phi and then take the best solution as starting point
        # for k in [1e-16, 1e-15, 1e-14, 1e-13, 1e-12]:
        #     for phi in [0.2, 0.15, 0.1, 0.05, 0.01]:
        #         initial_parameters = np.array([phi, k])
        #         # optimal_parameters, flag = leastsq(self.residual_function, initial_parameters) # if flag is 1 - found a good soln                
        #         # if self.data.error:
        #         #     # Calling with estimated error
        #         #     # print("Call with estimated error")
        #         optimal_parameters, flag = leastsq(self.residual_function, initial_parameters, args=(self.data.error)) # if flag is 1 - found a good soln
        #         # else: 
        #             # Calling without estimated error   
        #             # print("Call without estimated error")            
        #             # optimal_parameters, cov_x, infodict, mesg, flag = leastsq(self.residual_function, initial_parameters, full_output=1) # if flag is 1 - found a good soln
        #         self.data.set_approximation(self.model(optimal_parameters)) # Store the approximated data in the data structure
        #         chi_squared = self.__chi_squared(self.data.approximation)
        #         # print('Nfev = ', infodict['nfev'])
        #         # print('Chi squared: {} Function evaluated at output: {}'.format(chi_squared, np.sum(infodict['fvec'])))
        #         estimates[:, i] = optimal_parameters[0], optimal_parameters[1], chi_squared                
        #         # estimates[:, i] = optimal_parameters[0], optimal_parameters[1], abs(np.sum(infodict['fvec']))
        #         # print('Phi: {} k: {} Chi squared: {}'.format(phi, k, chi_squared))
        #         i += 1
        #         calls += 1
        #     #     if chi_squared <= 20:
        #     #         broken = True
        #     #         break
        #     # if broken:
        #     #     break
        # index = np.argmin(estimates[2])
        # optimal_parameters = estimates[:2, index]
        # self.data.set_unknown_parameters(optimal_parameters) # Store phi and k in data structure
        # self.data.set_approximation(self.model(optimal_parameters)) # Store the approximated data in the data structure
        # end_time = time.clock()
        # if verbose:
        #     print('Total model calls: {} Total time spent: {}'.format(self.calls, (end_time-start_time)))
        #     print('Optimal Phi: {} Optimal k: {}'.format(optimal_parameters[0], optimal_parameters[1]))
        # self.calls = 0
        # return optimal_parameters


    def generate_data(self, variables, parameters, time, noise = False, sd = 150., save_file=False, filename="testdata.dat"):
        """
        TODO: This function should be a part of data.py not model.py and should be rewritten

        Generate data using the forward model for given input of variables, and specified data.
        Note: Hasn't been tested since data.py was updated so may not work for theis solution and
              definitely will not work with homogeneous porous (radial1d) model.

        Inputs:
            variables (iterable) - Variable values to be used.
            parameters (iterable) - Well parameters to be used.
            time (iterable) - Times at which modelled response should be calculated.
            noise (bool) - If true noise will be added to the modelled response.
            sd (float) - Standard deviation of the noise/errors to be added.
            save_file (bool) - If true a text file will be written.
            filename (string) - Name of the file to be written.
        """
        # Create new data structure
        self.data = data_class.Data(model_type=self.model_type)
        self.data.set_time(time)
        self.data.set_known_parameters(parameters)

        # Calculate the modelled response
        modelled_value = self.model(variables)

        # Add noise to the data
        if noise:
            self.data.set_error(sd)
            np.random.seed(0) # Set random seed to 0 for consistency in testing
            # magnitude = 10**(np.floor(np.log10(np.average(p))))
            noise = np.random.randn(modelled_value.shape[0])
            # noise = (sd*(noise/np.std(noise)))*magnitude
            noise = sd * (noise/np.std(noise))
            modelled_value += noise
            # modelled_value += sd*np.random.randn(p.shape[0])
        
        # Save the data to the data structure
        self.data.set_observation(modelled_value)

        # Write a file
        if save_file:
            variables_dict = {'Porosity' : variables[0], 'Permeability' : variables[1]}
            self.data.generate_datafile(filename, variables=variables_dict)


class Theis_Solution(Model):
    """
    Analytical Theis solution model.
    """
    def __init__(self, data=None):
        """
        Initialise the model.
        """
        super().__init__(data=data, model_type='theis')

    def model(self, variables):
        """
        The forward model - calculates the modelled approximation of the well pressure response for a
        given set of input variables.

        Inputs: 
            variables (iterable) - List or numpy array of variable parameters for the model. This
                                   model uses the iterable format of [porosity, permeability].
        Outputs: 
            pressure (np.array) - Modelled approximation of well pressure response at the times specified 
                                  in the data structure.
        """
        # Unpack the variables
        phi, k = variables

        # Unpack the well parameters from the data structure
        p0 = self.data.reservoir_conditions['Initial Pressure']['Value']
        qm = self.data.fixed_parameters['Mass Flowrate']['Value']
        h = self.data.fixed_parameters['Layer Thickness']['Value']
        rho = self.data.fixed_parameters['Density']['Value']
        nu = self.data.fixed_parameters['Kinematic Viscosity']['Value']
        C = self.data.fixed_parameters['Compressibility']['Value']
        r = self.data.fixed_parameters['Action Well Radius']['Value']

        # Perform the calculations
        D = k/(nu*phi*rho*C)
        with np.errstate(divide="ignore", invalid="ignore"): # Hides 'RuntimeWarning: invalid value encountered in divide' if t[0] == 0.
            pressure = p0 + ((qm*nu)/(4*np.pi*k*h))*exp1((r**2)/(4*D*self.data.time))
        if self.data.time[0] <= 1e-7: # Check if initial reading is at time 0
            pressure[0] = p0

        return pressure

# from theis_wrapper import theis
# class Theis_Solution_Fortran(Theis_Solution):
    # def model(self, variables):
        # # p0, qm, h, rho, nu, C, r = self.data.parameters
        # p0 = self.data.parameters['Initial Pressure']['Value']
        # qm = self.data.parameters['Mass Flowrate']['Value']
        # h = self.data.parameters['Layer Thickness']['Value']
        # rho = self.data.parameters['Density']['Value']
        # nu = self.data.parameters['Kinematic Viscosity']['Value']
        # C = self.data.parameters['Compressibility']['Value']
        # r = self.data.parameters['Radius']['Value']
        # phi, k = variables
        # num_observations = len(self.data.time)
        # p = theis(k, nu, phi, rho, C, h, qm, p0, r, num_observations, self.data.time)
        # return p

from awtas.logic.wrappers.radial1d import radial1d_wrapper

class Radial_1D(Model):
    """
    Homogeneous porous, single layer, radial 1d numerical model. This model is a wrapper of Fortran code
    which performs the calculations.
    """
    def __init__(self, data=None):
        """
        Initialise the model.
        """
        super().__init__(data=data, model_type='radial1d')


    def model(self, variables):
        """
        The forward model - calculates the modelled approximation of the well response for a given
        set of input variables.

        Inputs: 
            variables (iterable) - List or numpy array of variable parameters for the model. This
                                   model uses the iterable format of [porosity, permeability].
        Outputs: 
            modelled_value (np.array) - Modelled approximation of well response at the times specified 
                                        in the data structure.
        """
        # Unpack the variables
        phi, k = variables

        # Unpack the well parameters from the data structure
        initial_pressure = self.data.reservoir_conditions['Initial Pressure']['Value']
        initial_x = self.data.reservoir_conditions['Initial X']['Value']
        well_radius = self.data.fixed_parameters['Action Well Radius']['Value']
        layer_thickness = self.data.fixed_parameters['Layer Thickness']['Value']
        rock_specific_heat = self.data.fixed_parameters['Rock Specific Heat']['Value']
        rock_heat_conductivity = self.data.fixed_parameters['Rock Heat Conductivity']['Value']
        rock_density = self.data.fixed_parameters['Rock Density']['Value']
        rock_compressibility = self.data.fixed_parameters['Rock Compressibility']['Value']
        pump_info = self.data.pump_info
        observation_points = self.data.observation_points
        grid_info = self.data.grid_info

        # Run the model to get the response
        modelled_value, execution_flag = radial1d_wrapper.radial1d(phi, k, layer_thickness, well_radius, rock_specific_heat, rock_heat_conductivity, rock_density, rock_compressibility,
            initial_pressure, initial_x, pump_info.injection_well, pump_info.injection_enthalpy, pump_info.num_pump_times, observation_points.num_observation_points, self.data.total_num_data, pump_info.pumping_scheme,
            pump_info.flow_times, pump_info.flow_rates, self.data.time, observation_points.radial_location, observation_points.num_data, observation_points.property, pump_info.deliverability, pump_info.production_index, pump_info.cutoff_pressure,
            grid_info.num_blocks, grid_info.num_constant_blocks, grid_info.constant_block_size, grid_info.block_growth_factor)

        return modelled_value


def create_model(model_type, data=None):
    """
    Creates an instance of a model derived class based on the model type.

    Inputs:
        model_type (string) - "theis" produces a model instance for the analytical theis solution,
                              "radial1d" produces a model instance for the homogeneous porous numerical solution.
        data (data instance) - Data structure of the same model type observed measurements, times of measurements
                               and the relevant well parameters for that model type.
    Outputs:
        model (model instance) - Model instance with functions to do forward and inverse modelling.
    """
    model_type = model_type.lower()
    if model_type == 'theis':
        model = Theis_Solution(data)
    # elif model_type == 'theis_fortran':
    #     # model = Theis_Solution_Fortran(data)
    #     pass
    elif model_type == 'radial1d':
        model = Radial_1D(data)
    else:
        raise ValueError('Error that model type does not exist.')
    return model