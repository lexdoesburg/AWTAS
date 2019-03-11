import numpy as np

# Dictionary defining the parameters used in each model - the key is the model_type (update when new models are added)
_model_parameters = {
    'theis' : # Well parameters required for the analytical theis solution model
        {
            'Reservoir Conditions' : ['Initial Pressure'],
            'Fixed Parameters' : ['Mass Flowrate', 'Layer Thickness', 'Density', 'Kinematic Viscosity', 'Compressibility', 'Action Well Radius'],
            'Variables' : ['Porosity', 'Permeability']
        },
    'radial1d' : # Well parameters required for the numerical homogeneous porous model
        {
            'Reservoir Conditions' : ['Initial Pressure', 'Initial X'],
            'Fixed Parameters' : ['Layer Thickness', 'Action Well Radius', 'Rock Specific Heat', 'Rock Heat Conductivity', 'Rock Density', 'Rock Compressibility'],
            'Variables' : ['Porosity', 'Permeability']
        }
    }

# Default units for each parameter. Can look to implement changing units via combo box in GUI later.
_default_parameter_units = {
    # The units are formatted to look better when displayed in the GUI
    'Reservoir Conditions' : 
        {
        'Initial Pressure' : {'Units':'Pa'},
        'Initial X' : {'Units': {'Initial Temperature':'°C', # Initial X can be either temperature or vapour saturation
                                 'Initial Vapour Saturation':'Dimensionless'}}
        },
    'Fixed Parameters' : 
        {
        'Mass Flowrate' : {'Units':'kg/s'},
        'Layer Thickness' : {'Units':'m'},
        'Density' : {'Units':'kg/m<sup>3</sup>'},
        'Kinematic Viscosity' : {'Units':'m<sup>2</sup>/s'},
        'Compressibility' : {'Units':'1/Pa'},
        'Action Well Radius' : {'Units':'m'},
        'Rock Specific Heat' : {'Units':'J/kgK'},
        'Rock Heat Conductivity' : {'Units':'W/mK'},
        'Rock Density' : {'Units':'kg/m<sup>3</sup>'},
        'Rock Compressibility' : {'Units':'1/Pa'}
        },
    'Variables' :
        {
        'Porosity' : {'Units':'Dimensionless'},
        'Permeability' : {'Units':'m<sup>2</sup>'}
        }    
    }

# Dictionary defining the upper and lower bounds of different parameters.
parameter_bounds = {
    'Porosity' : [0.01, 0.2],
    'Permeability' : [1e-16, 1e-12]
    }

class Data():
    """
    This is a class which defines the data structure which contains all of the relevant information for a given well
    model.
    """
    def __init__(self, model_type, filename=None, time=None, observation=None, parameters_list=None, pump_info=None, observation_points=None, grid_info=None, error=1):
        """
        Initialise the data structure.

        Note (#TODO:): The theis solution model uses an older form of the data structure, and therefore there is a lot of code that is 
                       specific to the theis solution. This will eventually updated so all model types use the same form of the data structure.

        Inputs:
            model_type (string) - Defines the model type the data will be used for: either "theis" for analytical theis solution model
                                  or "radial1d" for 1d radial homogeneous porous model.
            filename (string) - Filename of a datafile to read the data from.
            time (iterable) - Array of times where observations were recorded and where modelled values should be calculated.
            observation (iterable) - Array of observation values.
            parameters_list (iterable) - Array of known well parameters. Note the order is specified for different model types.
                                         See set_known_parameters docstring below for the ordering for each model type.
            pump_info (Pump object instance) - Object which stores the pump information for the model. (See Pump class below for more details).
            observation_points (ObservationPoints object instance) - Object which stores the observation point information. (See ObservationPoints
                                                                     class below for more details).
            grid_info (Grid object instance) - Object which stores the grid generation information. (See Grid class below for more details).
            error (float) - Value defining the known standard deviation of errors from the observed data reading (only used for Theis solution model)
        """
        # Check the model type is valid
        self.model_type = model_type.lower()
        assert self.model_type_valid(), 'Input model type is invalid try one of {} instead.'.format(_model_parameters.keys())
        self.approximation = None # array of approximated pressure data

        # Initialise the parameter dictionaries which will hold the fixed well parameters
        self.reservoir_conditions, self.fixed_parameters, self.variables = self._initialise_parameter_dictionaries()

        self.pump_info = pump_info
        self.observation_points = observation_points
        self.grid_info = grid_info
        self.initial_x = None # String that states if initial x is temperature or vapour saturation (for use within the GUI)
        
        # Initialise the time, observations and errors.
        if self.model_type is not 'theis' and self.observation_points:
            self.time = np.concatenate(self.observation_points.times)
            self.observation = np.concatenate(self.observation_points.observations)
            self._update_errors()
        else:
            self.time = time # array of time data
            self.observation = observation # array of observed pressure data
            self.error = error # Estimated error in the readings # Move the error to observation points (each could have different errors)

        # Initialise the total number of data
        if model_type == 'radial1d' and observation_points:
            self.total_num_data = sum(observation_points.num_data)
        elif model_type == 'theis' and time:
            self.total_num_data = len(time)
        else:
            self.total_num_data = 0

        # If a list of parameters was supplied fill the parameter dictionaries
        if parameters_list:
            self._fill_parameter_dictionaries(parameters_list)
        
        # If a filename is supplied read the data in from that file
        if filename:
            self.read_file(filename)
    

    def _initialise_parameter_dictionaries(self):
        """
        Helper function used during initialisation to build dictionaries that will hold the fixed well parameter data.

        The dictionary is of the form {'Well parameter' : {'Value' : <some value>, 'Units' : <units>}}
        The value of a given parameter can be found by: 
            value = parameter_dictionary['Well Parameter']['Value'], and similarly for the parameter's units.
        """
        # Initialise the dictionaries
        reservoir_conditions = {}
        fixed_parameters = {}
        variables = {}

        # Fill the dictionaries with only the relevant model parameters as defined in the _model_parameters dictionary
        for reservoir_condition in _model_parameters[self.model_type]['Reservoir Conditions']:
            reservoir_conditions[reservoir_condition] = {'Value':None, 'Units':_default_parameter_units['Reservoir Conditions'][reservoir_condition]['Units']}
        for parameter in _model_parameters[self.model_type]['Fixed Parameters']:
            fixed_parameters[parameter] = {'Value':None, 'Units':_default_parameter_units['Fixed Parameters'][parameter]['Units']}
        for variable in _model_parameters[self.model_type]['Variables']:
            variables[variable] = {'Value':None, 'Units':_default_parameter_units['Variables'][variable]['Units']}

        return reservoir_conditions, fixed_parameters, variables


    def _fill_parameter_dictionaries(self, parameters_list): 
        """
        Helper function used to fill the parameter dictionaries with values.

        Inputs:
            parameters_list (iterable) - Array of known well parameters. Note the order is specified for different model types.
                                         See set_known_parameters docstring below for the ordering for each model type.
        """       
        if parameters_list:
            if self.model_type == 'theis':
                self.reservoir_conditions['Initial Pressure']['Value'] = parameters_list[0]
                self.fixed_parameters['Mass Flowrate']['Value'] = parameters_list[1]
                self.fixed_parameters['Layer Thickness']['Value'] = parameters_list[2]
                self.fixed_parameters['Density']['Value'] = parameters_list[3]
                self.fixed_parameters['Kinematic Viscosity']['Value'] = parameters_list[4]
                self.fixed_parameters['Compressibility']['Value'] = parameters_list[5]
                self.fixed_parameters['Action Well Radius']['Value'] = parameters_list[6]
            elif self.model_type == 'radial1d':
                self.reservoir_conditions['Initial Pressure']['Value'] = parameters_list[0]
                self.reservoir_conditions['Initial X']['Value'] = parameters_list[1]
                self.fixed_parameters['Action Well Radius']['Value'] = parameters_list[2]
                self.fixed_parameters['Layer Thickness']['Value'] = parameters_list[3]
                self.fixed_parameters['Rock Specific Heat']['Value'] = parameters_list[4]
                self.fixed_parameters['Rock Heat Conductivity']['Value'] = parameters_list[5]
                self.fixed_parameters['Rock Density']['Value'] = parameters_list[6]
                self.fixed_parameters['Rock Compressibility']['Value'] = parameters_list[7]
                self._set_initial_x()
               

    def _update_errors(self):
        """
        Helper function used to set the error array from each different observation point standard deviation of observed
        error
        """
        errors = []
        # Concatenate the error from each observation point into a single list
        for i in range(self.observation_points.num_observation_points):
            obs_point_error = [self.observation_points.errors[i]] * self.observation_points.num_data[i]
            errors.extend(obs_point_error)

        # Convert the error list into a numpy array
        self.error = np.fromiter(errors, dtype=float)


    def _set_initial_x(self):
        """
        Helper function used to change the label of the 'Initial X' parameter depending on if it is an initial vapour saturation
        or an initial temperature. This is used by the GUI to show the proper label instead of 'Initial X'.
        """
        if self.model_type == 'radial1d':
            if self.reservoir_conditions['Initial X']['Value'] < 1.0:
                self.initial_x = 'Initial Vapour Saturation'
            else:
                self.initial_x = 'Initial Temperature'


    def read_file(self, filename):
        """
        Populates the data structure using a data file. The format of the data files is strict. Examples of the file format
        for each model can be found in /awtas/gui/example_datafiles.

        Note (#TODO:): Currently there is a different format for the theis solution model as compared to the numerical
                       homogeneous porous (radial1d) model. This will be updated such that all models use the same format.
        
        Inputs:
            filename (string) - the filename or path to a file to be read.
        """
        if self.model_type == 'theis':
            # Run the old code specifically for the analytical theis solution
            with open(filename, 'r') as file:
                # Read the error value
                file.readline() # Skip first line
                info_values = file.readline()
                info_values = info_values.rstrip().split(',') # Convert string to list of strings
                try:
                    self.error = float(info_values[0])
                    if self.error <= 1e-6:
                        self.error = 1
                except ValueError:
                    # If the error can't be converted to a float then don't store any error value
                    self.error = 1
                file.readline() # Skip blank line

                # Read parameter values
                parameter_labels = file.readline()
                parameter_labels = parameter_labels.rstrip().split(',') # Convert string to list of strings
                parameter_values = file.readline() 
                parameter_values = [float(value) for value in parameter_values.split(',')] # Convert string to list of floats
                for parameter_name, value in zip(parameter_labels, parameter_values):
                    if parameter_name in self.reservoir_conditions.keys():
                        self.reservoir_conditions[parameter_name]['Value'] = value # This works only if the file uses the same parameter names as the dictionary key.
                    elif parameter_name in self.fixed_parameters.keys():
                        self.fixed_parameters[parameter_name]['Value'] = value # This works only if the file uses the same parameter names as the dictionary key.
                    else:
                        raise ValueError('Parameter name "{}" not recognised from input data file.'.format(parameter_name))
            
            # Read the time and observation values
            self.time, self.observation = np.genfromtxt(filename, delimiter=',', skip_header=7).T
        else: 
            # Run the new code for homogeneous porous (radial1d) and later models

            # Initialise lists for the creation of the observation points object
            obs_point_properties = []
            obs_point_locations = []
            obs_point_num_data = []
            obs_point_times = []
            obs_point_observations = []
            obs_point_errors = []

            # Start reading file
            with open(filename, 'r') as file:
                # for counter, line in enumerate(file):
                for line in file:
                    # print('Line {}: {}'.format(counter, line))
                    if 'KNOWN WELL PARAMETERS' in line:
                        # Read the known well parameters into the parameter dictionaries
                        labels = file.readline()
                        labels = labels.rstrip().split(',') # Convert string to list of strings
                        values = file.readline()
                        values = [float(value) for value in values.split(',')] # Convert string to list of floats
                        for label, value in zip(labels, values):
                            if label in self.reservoir_conditions.keys():
                                self.reservoir_conditions[label]['Value'] = value # This works only if the file uses the same parameter names as the dictionary key.
                            elif label in self.fixed_parameters.keys():
                                self.fixed_parameters[label]['Value'] = value # This works only if the file uses the same parameter names as the dictionary key.
                            else:
                                raise ValueError('Parameter name "{}" not recognised from input data file. Try one of {} or {} instead.'.format(label, self.reservoir_conditions.keys(), self.fixed_parameters.keys()))
                        self._set_initial_x() # Update the initial x label
                    elif 'GRID INFORMATION' in line:
                        # Read the grid information
                        labels = file.readline()
                        labels = labels.rstrip().split(',') # Convert string to list of strings
                        values = file.readline()
                        values = values.rstrip().split(',') # Convert string to list of strings
                        for label, value in zip(labels, values):
                            if label == 'Number of Grid Blocks':
                                num_blocks = np.int32(value)
                            elif label == 'Number of Constant Sized Blocks':
                                num_constant_blocks = np.int32(value)
                            elif label == 'Constant Block Size':
                                constant_block_size = float(value)
                            elif label == 'Block Growth Factor':
                                block_growth_factor = float(value)

                        # Update the grid info
                        self.grid_info = Grid(num_blocks=num_blocks, num_constant_blocks=num_constant_blocks, constant_block_size=constant_block_size, block_growth_factor=block_growth_factor)
                    elif 'PUMP INFORMATION' in line:
                        # Read the pump information
                        labels = file.readline()
                        labels = labels.rstrip().split(',') # Convert string to list of strings
                        values = file.readline()
                        values = values.rstrip().split(',') # Convert string to list of strings
                        injection_enthalpy = None
                        for label, value in zip(labels, values):
                            if label == 'Pumping Scheme':
                                pumping_scheme = value
                            elif label == 'Number of Pump Times':
                                num_pump_times = int(value)
                            elif label == 'Deliverability':
                                if value.lower() in ['false', 'no']:
                                    deliverability = False
                                elif value.lower() in ['true', 'yes']:
                                    deliverability = True
                                else:
                                    deliverability = bool(int(value))
                            elif label == 'Production Index':
                                production_index = float(value)
                            elif label == 'Cutoff Pressure':
                                cutoff_pressure = float(value)
                            elif label == 'Injection Enthalpy':
                                injection_enthalpy = float(value)

                        # Skip 2 lines
                        file.readline()
                        file.readline()

                        # Initialise lists to store the pump time and flow rates
                        pump_times = []
                        pump_rates = []

                        # Read the pump time and flow rates
                        for i in range(num_pump_times):
                            flow_info = file.readline()
                            flow_info = [float(value) for value in flow_info.split(',')]
                            pump_times.append(flow_info[0])
                            pump_rates.append(flow_info[1])
                        
                        # Update the pump information
                        self.pump_info = Pump(pumping_scheme=pumping_scheme,flow_rates=pump_rates, flow_times=pump_times, deliverability=deliverability, production_index=production_index, cutoff_pressure=cutoff_pressure, injection_enthalpy=injection_enthalpy)
                    elif 'OBSERVATION POINT' in line:
                        # Read the information about the observation point
                        error_supplied = False # True if an error is specified
                        labels = file.readline()
                        labels = labels.rstrip().split(',') # Convert string to list of strings
                        values = file.readline()
                        values = values.rstrip().split(',') # Convert string to list of strings
                        for label, value in zip(labels, values):
                            if label == 'Property Observed':
                                obs_point_properties.append(value)
                            elif label == 'Radial Location [m]':
                                obs_point_locations.append(float(value))
                            elif label == 'Number of Observations':
                                num_data = int(value)
                                obs_point_num_data.append(num_data)
                            elif label == 'Error Standard Deviation':
                                error_supplied = True
                                error = float(value)
                                obs_point_errors.append(error)
                        if not error_supplied:
                            # If no errors supplied set the error to 1, essentially meaning the error is unknown
                            obs_point_errors.append(1)
                        
                        # Skip 2 lines
                        file.readline()
                        file.readline()

                        # Read the time and observation values
                        times = np.ndarray(shape=num_data, dtype=float)
                        observations = np.ndarray(shape=num_data, dtype=float)
                        for i in range(num_data):
                            observation_info = file.readline()
                            observation_info = [float(value) for value in observation_info.split(',')]
                            times[i] = observation_info[0]
                            observations[i] = observation_info[1]
                        obs_point_times.append(times)
                        obs_point_observations.append(observations)

            # Update the remaining items
            self.observation_points = ObservationPoints(radial_location=obs_point_locations, property=obs_point_properties, num_data=obs_point_num_data, times=obs_point_times, observations=obs_point_observations, errors=obs_point_errors)                
            self.time = np.concatenate(self.observation_points.times)
            self.observation = np.concatenate(self.observation_points.observations)
            self.total_num_data = len(self.time)
            self._update_errors()


    def set_known_parameters(self, parameters):
        """
        Sets the known well parameters.
        
        Inputs:
            parameters (iterable) - Array of parameters in a specified order based on the model type.
                Analytical Theis Model Ordering:
                    parameters[0] = intial pressure
                    parameters[1] = mass flowrate
                    parameters[2] = thickness
                    parameters[3] = density
                    parameters[4] = kinematic viscosity
                    parameters[5] = compressibility
                    parameters[6] = action well radius

                Radial 1D Model Ordering:
                    parameters[0] = intial pressure
                    parameters[1] = initial temperature/vapour saturation
                    parameters[2] = action well radius
                    parameters[3] = layer thickness
                    parameters[4] = rock specific heat
                    parameters[5] = rock heat conductivity
                    parameters[6] = rock density    
                    parameters[7] = rock compressibility    
        """
        self._fill_parameter_dictionaries(parameters)


    def set_unknown_parameters(self, variables):
        """
        Sets the unknown well parameters.

        Inputs:
            variables (iterable) - Array of variables in a specified order.
                Ordering:
                variables[0] = porosity
                variables[1] = permeability
        """
        self.variables['Porosity']['Value'] = variables[0]
        self.variables['Permeability']['Value'] = variables[1]


    def set_time(self, time):
        """
        Set the time array.

        Inputs:
            time (iterable) - Array of time values (floats).
        """
        # Update the total number of data
        self.total_num_data = len(time)

        # Set the time
        self.time = time
        if self.model_type == 'radial1d':
            self.observation_points._store_times(time)


    def set_observation(self, observation):
        """
        Set the observation array.

        Inputs:
            observation (iterable) - Array of observation values (floats).
        """
        self.observation = observation
        if self.model_type == 'radial1d':
            self.observation_points._store_observation(observation)


    def set_approximation(self, approximation):
        """
        Set the modelled approximation array.

        Inputs:
            approximation (iterable) - Array of modelled values (floats).
        """
        self.approximation = approximation
        if self.model_type == 'radial1d':
            self.observation_points.store_modelled_values(approximation)


    def set_error(self, error):
        """
        Set the error.

        Note: Currently only suitable for theis solution model.

        Inputs:
            error (float) - Standard deviation of errors of observation data
        """
        assert self.model_type == 'theis', 'The set_error function is not yet working for {} model.'.format(self.model_type)
        self.error = error
    

    def get_well_type(self):
        """
        Determines the type of well the data corresponds to (for use within the GUI).
        """
        if self.model_type == 'theis':
            flow_rate = self.fixed_parameters['Mass Flowrate']['Value']
            if flow_rate > 0:
                well_type = 'Injection'
            elif flow_rate < 0:
                well_type = 'Production'
        else:
            if self.pump_info.injection_well == 1:
                well_type = 'Injection'
            else:
                well_type = 'Production' 
        return well_type


    def create_pump_info(self, pumping_scheme, flow_rates, flow_times, deliverability, production_index, cutoff_pressure, injection_enthalpy=None):
        """
        Sets the pump info.

        Inputs:
            pumping_scheme (string) - Defines the type of flow data, either 'Measured Flows', 'Constant Flow' or 'Step Flows'
            flow_rates (float/iterable) - Pump flow rates
            flow_times (float/iterable) - Times at which the flows occur. Note: There is a special case for the pumping scheme
                                          'Step Flows' where the flow_time is the duration of the flow, not the time the flow begins.
            deliverability (bool) - True if the pump is on deliverability.
            production_index (float) - Production index for deliverability.
            cutoff_pressure (float) -  Cuttoff pressure for deliverability.
            injection_enthalpy (float, optional) - If the well is an injection well then supply injection enthalpy.
        """
        self.pump_info = Pump(pumping_scheme, flow_rates, flow_times, deliverability, production_index, cutoff_pressure, injection_enthalpy=injection_enthalpy)
    

    def create_observation_points(self, radial_location, property, num_data, times, observations, errors=None):
        """
        Sets the observation point info.

        Inputs:
            radial_location (float/list) - The radial location of the observation point in meters. Or a list of radial locations
                                           for a set of observation points.
            property (string/list) - String defining the property observed, either 'Deliverability', 'Pressure', 'Temperature'
                                         or 'Enthalpy'. Or a list of properties for a set of observation points.
            num_data (integer/list) - Total number of observation measurements. Or a list of numbers of data for a set of observation
                                      points.
            times (np.ndarray/list) - Numpy array of times. Or a list of numpy arrays for a set of observation points.
            observations (np.ndarray/list) - Numpy array of observations. Or a list of numpy arrays for a set of observation points.
            errors (float/int/list, optional) - Standard deviation of observation errors. Or a list of errors for a set of observation 
                                                points. If not specified, error is assumed unknown.
        """
        self.observation_points = ObservationPoints(radial_location, property, num_data, times, observations, errors)


    def create_grid_info(self, num_blocks, num_constant_blocks, constant_block_size, block_growth_factor):
        """
        Sets the grid info.

        Inputs:
            num_blocks (integer) - The number of grid blocks to be generated.
            num_constant_blocks (integer) - The number of grid blocks which will have a fixed size.
            constant_block_size (float) - The fixed size which constant blocks will have (meters).
            block_growth_factor (float) - The factor by which non-constant blocks will grow by.
        """
        self.grid_info = Grid(num_blocks, num_constant_blocks, constant_block_size, block_growth_factor)


    def model_type_valid(self):
        """
        Check if the model type is valid.
        """
        if self.model_type in _model_parameters.keys():
            return True
        else:
            return False


    def add_noise(self, sd):
        """
        Adds noise of a specified standard deviation to the observation data.

        Inputs:
            sd (float) - Standard deviation of noise added.
        """
        if self.observation is not None:
            np.random.seed(0) # Set random seed to 0 for consistency in testing
            noise = np.random.randn(self.observation.shape[0])
            noise = sd * (noise/np.std(noise))
            self.observation += noise
            self.observation_points._store_observation(self.observation)


    def write_output_file(self, filename):
        """
        Writes an output file containing the estimated variable values and the approximated model
        response.

        Inputs:
            filename (string) - File name or path.
        """
        with open(filename, 'w') as file:
            # First write the variable values
            file.write('ESTIMATED VARIABLE VALUES\n')
            variable_labels = ''
            variable_values = ''
            for i, (variable, info) in enumerate(self.variables.items()):
                if i != len(self.variables.keys()) - 1:
                    variable_labels += '{},'.format(variable)
                    variable_values += '{},'.format(info['Value'])
                else:
                    variable_labels += '{}\n'.format(variable)
                    variable_values += '{}\n'.format(info['Value'])
            file.write(variable_labels)
            file.write(variable_values)
            file.write('\n')
            
            # Next write the modelled approximation
            if self.model_type == 'theis':
                file.write('MODELLED APPROXIMATION\n')
                file.write('Time [s],Approximated Value\n')
                for i in range(len(self.time)):
                    file.write('{},{}\n'.format(self.time[i], self.approximation[i]))
            else:
                for i in range(self.observation_points.num_observation_points):
                    file.write('OBSERVATION POINT {} MODELLED APPROXIMATION\n'.format(i+1))
                    file.write('Time [s],Approximated Value\n')
                    for j in range(self.observation_points.num_data[i]):
                        file.write('{},{}\n'.format(self.observation_points.times[i][j],
                                                    self.observation_points.modelled_values[i][j]))
                    if i != self.observation_points.num_observation_points - 1: # Write a new line between each observation point
                        file.write('\n')


    def generate_datafile(self, filename, variables=None):
        """
        Generates a data file using the current data stored in the data structure.

        Inputs:
            filename (string) - File name or path.
            variables (dictionary, optional) - A dictionary that contains known values of the 'unknown'
                                               variables. Useful for checking how well the unknown variables
                                               have been fit.
        """
        if self.model_type == 'theis':
            # Run old code for theis solution
            with open(filename, 'w') as file:
                # Write known errors and variable values
                info_labels = 'Error Standard Deviation,'
                if self.error:
                    info_values = '{},'.format(self.error)
                else:
                    info_values = 'Unknown,'
                if variables:
                    for i, (variable, value) in enumerate(variables.items()):
                        if i != len(variables.keys()) - 1:
                            info_labels += '{},'.format(variable)
                            info_values += '{},'.format(value)
                        else:
                            info_labels += '{}\n'.format(variable)
                            info_values += '{}\n'.format(value)
                else:
                    info_labels += 'Variable Values\n'
                    info_values += 'Unknown\n'
                file.write(info_labels)
                file.write(info_values)
                file.write('\n')

                # Write known well parameters
                parameter_labels = ''
                parameter_values = ''
                for parameter in _model_parameters[self.model_type]['Reservoir Conditions']:
                    parameter_labels += '{},'.format(parameter)
                    parameter_values += '{},'.format(self.reservoir_conditions[parameter]['Value'])
                for parameter in _model_parameters[self.model_type]['Fixed Parameters']:
                    if parameter != _model_parameters[self.model_type]['Fixed Parameters'][-1]:
                        parameter_labels += '{},'.format(parameter)
                        parameter_values += '{},'.format(self.fixed_parameters[parameter]['Value'])
                    else:
                        parameter_labels += '{}\n'.format(parameter)
                        parameter_values += '{}\n'.format(self.fixed_parameters[parameter]['Value'])
                file.write(parameter_labels)
                file.write(parameter_values)
                file.write('\n')
                
                # Write the time of observation and observed pressure readings:
                file.write('Time [s],Pressure Observation [Pa]\n')
                for i in range(len(self.time)):
                    file.write('{},{}\n'.format(self.time[i], self.observation[i]))
        else: 
            # New code for homogeneous porous (radial 1d) and later models
            with open(filename, 'w') as file:
                # If known variable values are specified write them first
                if variables:
                    file.write('KNOWN VARIABLE VALUES\n')
                    variable_labels = ''
                    variable_values = ''
                    for i, (variable, value) in enumerate(variables.items()):
                        if i != len(variables.keys()) - 1:
                            variable_labels += '{},'.format(variable)
                            variable_values += '{},'.format(value)
                        else:
                            variable_labels += '{}\n'.format(variable)
                            variable_values += '{}\n'.format(value)
                    file.write(variable_labels)
                    file.write(variable_values)
                    file.write('\n')

                # Write known well parameters
                file.write('KNOWN WELL PARAMETERS\n')
                parameter_labels = ''
                parameter_values = ''
                for parameter in _model_parameters[self.model_type]['Reservoir Conditions']:
                    parameter_labels += '{},'.format(parameter)
                    parameter_values += '{},'.format(self.reservoir_conditions[parameter]['Value'])
                for parameter in _model_parameters[self.model_type]['Fixed Parameters']:
                    if parameter != _model_parameters[self.model_type]['Fixed Parameters'][-1]:
                        parameter_labels += '{},'.format(parameter)
                        parameter_values += '{},'.format(self.fixed_parameters[parameter]['Value'])
                    else:
                        parameter_labels += '{}\n'.format(parameter)
                        parameter_values += '{}\n'.format(self.fixed_parameters[parameter]['Value'])
                file.write(parameter_labels)
                file.write(parameter_values)
                file.write('\n')

                # Write grid information if grid required for model
                if self.grid_info:
                    file.write('GRID INFORMATION\n')
                    grid_labels = 'Number of Grid Blocks,Number of Constant Sized Blocks,Constant Block Size,Block Growth Factor\n'
                    grid_values = '{},{},{},{}\n'.format(self.grid_info.num_blocks, self.grid_info.num_constant_blocks,
                                                        self.grid_info.constant_block_size, self.grid_info.block_growth_factor)
                    file.write(grid_labels)
                    file.write(grid_values)
                    file.write('\n')

                # Write pump information
                file.write('PUMP INFORMATION\n')
                # Hard-coded - could think of a better solution
                pump_info_labels = 'Pumping Scheme,Number of Pump Times,Deliverability,Production Index,Cutoff Pressure'
                pump_info_values = '{},{},{},{},{}'.format(_pump_scheme_from_flag[self.pump_info.pumping_scheme], self.pump_info.num_pump_times,
                                                        self.pump_info.deliverability, self.pump_info.production_index,
                                                        self.pump_info.cutoff_pressure)
                if self.pump_info.injection_well == 1:
                    pump_info_labels += ',Injection Enthalpy\n'
                    pump_info_values += ',{}\n'.format(self.pump_info.injection_enthalpy)
                else:
                    pump_info_labels += '\n'
                    pump_info_values += '\n'
                file.write(pump_info_labels)
                file.write(pump_info_values)
                file.write('\n')

                # Write the pump times and pump rates
                file.write('Time [s],Mass Flowrate [kg/s]\n')
                for i in range(self.pump_info.num_pump_times):
                    file.write('{},{}\n'.format(self.pump_info.flow_times[i], self.pump_info.flow_rates[i]))
                file.write('\n')

                # Write observation point information
                for i in range(self.observation_points.num_observation_points):
                    file.write('OBSERVATION POINT {}\n'.format(i+1))
                    # Hard-coded - could think of a better solution
                    observation_point_labels = 'Property Observed,Radial Location [m],Number of Observations'
                    observation_point_values = '{},{},{}'.format(_observation_property_from_flag[self.observation_points.property[i]],
                                                                self.observation_points.radial_location[i],
                                                                self.observation_points.num_data[i])
                    if self.observation_points.errors[i] <= 1.0 + 1e-6:
                        observation_point_labels += '\n'
                        observation_point_values += '\n'
                    else:       
                        observation_point_labels += ',Error Standard Deviation\n'
                        observation_point_values += ',{}\n'.format(self.observation_points.errors[i])
                            
                    file.write(observation_point_labels)
                    file.write(observation_point_values)
                    file.write('\n')

                    # Write all observations
                    file.write('Time [s],Observation\n')
                    for j in range(self.observation_points.num_data[i]):
                        file.write('{},{}\n'.format(self.observation_points.times[i][j], self.observation_points.observations[i][j]))
                    if i != self.observation_points.num_observation_points - 1: # Write a new line between each observation point
                        file.write('\n')


#----------------------------------------------------------------------------------------------------------
# Observation point information class 
#----------------------------------------------------------------------------------------------------------

# Type of measurement data that was recorded at the observation point (mapping readable names to the fortran equivalent flags)
_observation_properties = {
    'Deliverability' : 0,
    'Pressure' : 1,
    'Temperature' : 2,
    'Enthalpy' : 3
    }

# Opposite ordering of the above dictionary to get meaningful name from fortran flag
_observation_property_from_flag = {
    0 : 'Deliverability',
    1 : 'Pressure',
    2 : 'Temperature',
    3 : 'Enthalpy'
    }

class ObservationPoints:
    """
    Observation point class which holds the information about all observation points for a given well. 
    This class is necessary because the fortran model wrappers require the input data to be of specific types,
    and that can be controlled here.
    """
    """
    TODO:
    Observation point class isn't very intuitive. Would be better to have a single class and each observation point object 
    is a single instance of that class. Currently the class is for multiple observation points.
    This way was the quickest way to work best with the fortran code currently - can be improved later to use single object instance of
    the class for each observation point and then combine the required information from each obs point to be passed to fortran.
    """
    def __init__(self, radial_location=[], property=[], num_data=[], times=[], observations=[], errors=None):
        """
        Initialise the observation points information class.

        Inputs:
            radial_location (float/list) - The radial location of the observation point in meters. Or a list of radial locations
                                           for a set of observation points.
            property (string/list) - String defining the property observed, either 'Deliverability', 'Pressure', 'Temperature'
                                         or 'Enthalpy'. Or a list of properties for a set of observation points.
            num_data (integer/list) - Total number of observation measurements. Or a list of numbers of data for a set of observation
                                      points.
            times (np.ndarray/list) - Numpy array of times. Or a list of numpy arrays for a set of observation points.
            observations (np.ndarray/list) - Numpy array of observations. Or a list of numpy arrays for a set of observation points.
            errors (float/int/list, optional) - Standard deviation of observation errors. Or a list of errors for a set of observation
                                                points. If not specified, error is assumed unknown.
        """
        # If inputs are not iterable, convert to iterables
        if type(radial_location) is float:
            radial_location = [radial_location]
        if type(property) is str:
            property = [property]
        if type(num_data) is int:
            num_data = [num_data]
        if type(times) is np.ndarray:
            times = [times]
        if type(observations) is np.ndarray:
            observations = [observations]
        if type(errors) in [int, float]:
            errors = [errors] * len(num_data)
        elif errors is None:
            errors = [1] * len(num_data)

        # Assert the information lists are the same size 
        assert len(radial_location) == len(property) == len(num_data) == len(times) == len(observations) == len(errors), 'Trying to add multiple observation points but insufficent information supplied.'
        
        # Set the propertys to the correct fortran flag based on the input string
        if property:
            self.property = [_observation_properties[p] for p in property]
        else:
            self.property = property

        # Finish initialising the relevant information, using the correct data types.        
        self.num_observation_points = len(radial_location)
        self.radial_location = np.fromiter(radial_location, dtype=float)
        self.num_data = np.fromiter(num_data, dtype=np.int32)
        self.property = np.fromiter(self.property, dtype=np.int32)
        self.times = times
        self.observations = observations
        self.errors = errors
        self.modelled_values = []

        # Units will be useful when displaying observation point information in the GUI. 
        # TODO: Display observation point information in GUI.
        self.units = {
            'Radial location' : 'm'
        }

    
    def add_observation_points(self, radial_location, property, num_data, times, observations, errors=None):
        """
        Add observation point(s). Either takes lists each in identical order, or individual values (for one observation point).

        Inputs:
            radial_location (float/list) - The radial location of the observation point in meters. Or a list of radial locations
                                           for a set of observation points.
            property (string/list) - String defining the property observed, either 'Deliverability', 'Pressure', 'Temperature'
                                         or 'Enthalpy'. Or a list of properties for a set of observation points.
            num_data (integer/list) - Total number of observation measurements. Or a list of numbers of data for a set of observation
                                      points.
            times (np.ndarray/list) - Numpy array of times. Or a list of numpy arrays for a set of observation points.
            observations (np.ndarray/list) - Numpy array of observations. Or a list of numpy arrays for a set of observation points.
            errors (float/int/list, optional) - Standard deviation of observation errors. Or a list of errors for a set of observation
                                                points. If not specified, error is assumed unknown.
        """
        if type(radial_location) is list and type(property) is list and type(num_data) is list and type(times) is list and type(observations) is list:
            # Multiple observation points being added
        
            assert len(radial_location) == len(property) == len(num_data) == len(times) == len(observations), 'Trying to add multiple observation points but insufficent information supplied.'
            if errors is not None:
                if type(errors) is list:
                    # Assert the number of errors supplied is correct
                    assert len(errors) == len(radial_location)
                    self.errors.extend(errors)
                else:
                    errors = [errors] * len(radial_location)
                    self.errors.extend(errors)
            else:
                errors = [1] * len(radial_location)
                self.errors.extend(errors)
            
            # Add new observation information to existing arrays.
            np.append(self.radial_location, radial_location)
            np.append(self.property, [_observation_properties[property_type] for property_type in property])
            np.append(self.num_data, num_data)
            self.times.extend(times)
            self.observations.extend(observations)
        elif type(radial_location) in [float, int] and type(property) is str and type(num_data) is int:
            # Single observation point being added
                    
            if errors is not None:
                if type(errors) in [int, float]:
                    self.errors.append(errors)
            else:
                errors = 1
                self.errors.append(errors)
            np.append(self.radial_location, radial_location)
            np.append(self.property, _observation_properties[property])
            np.append(self.num_data, num_data)
            assert len(times) == num_data == len(observations), 'Given {} observation times and {} observation values when expecting {}.'.format(len(times), len(observations), num_data)
            self.times.append(times)
            self.observations.append(observations)
            
    
    # The next 3 functions are all very similar - could possibly be moved to a single function.
    def store_modelled_values(self, modelled_values):
        """
        Sets the modelled values for each observation point.

        Inputs:
            modelled_values (iterable) - Array of modelled values.
        """
        current_index = 0
        for i in range(self.num_observation_points):
            next_index = current_index + self.num_data[i]
            self.modelled_values.append(np.fromiter(modelled_values[current_index:next_index], dtype=float))
            current_index = next_index


    def _store_observation(self, observation):
        """
        Sets the observations for each observation point.

        Inputs:
            observation (iterable) - Array of observations.
        """
        current_index = 0
        self.observations = []
        for i in range(self.num_observation_points):
            next_index = current_index + self.num_data[i]
            self.observations.append(np.fromiter(observation[current_index:next_index], dtype=float))
            current_index = next_index


    def _store_times(self, time):
        """
        Sets the times for each observation point.

        Inputs:
            time (iterable) - Array of times.
        """
        current_index = 0
        self.times = []
        for i in range(self.num_observation_points):
            next_index = current_index + self.num_data[i]
            self.observations.append(np.fromiter(time[current_index:next_index], dtype=float))
            current_index = next_index


#----------------------------------------------------------------------------------------------------------
# Pump information class
#----------------------------------------------------------------------------------------------------------

# Type of flow that is occurring (mapping readable names to the fortran equivalent flags)
_pump_schemes = {
    'Measured Flows' : 0, # Flowrate changes at different times
    'Constant Flow' : 1,  # Flowrate constant through the whole time period
    'Step Flows' : 4      # Flow only active for a certain amount of time
    }

# Opposite ordering of the above dictionary to get meaningful name from fortran flag
_pump_scheme_from_flag = {
    0 : 'Measured Flows', # Flowrate changes at different times
    1 : 'Constant Flow',  # Flowrate constant through the whole time period
    4 : 'Step Flows'      # Flow only active for a certain amount of time
    }

class Pump:
    """
    Pump class which holds the information about the pump for a given well. 
    """
    def __init__(self, pumping_scheme, flow_rates, flow_times, deliverability, production_index, cutoff_pressure, injection_enthalpy=None):
        """
        Initialise the pump object.

        Inputs:
            pumping_scheme (string) - Defines the type of flow data, either 'Measured Flows', 'Constant Flow' or 'Step Flows'
            flow_rates (float/iterable) - Pump flow rates
            flow_times (float/iterable) - Times at which the flows occur. Note: There is a special case for the pumping scheme
                                          'Step Flows' where the flow_time is the duration of the flow, not the time the flow begins.
            deliverability (bool) - True if the pump is on deliverability.
            production_index (float) - Production index for deliverability.
            cutoff_pressure (float) -  Cuttoff pressure for deliverability.
            injection_enthalpy (float, optional) - If the well is an injection well then supply injection enthalpy.
        """
        # Set the pumping scheme to the correct fortran flag based on the input string
        self.pumping_scheme = _pump_schemes[pumping_scheme]

        # If input rates and times are not iterable, convert to iterables
        if type(flow_rates) in [float, int]:
            flow_rates = [flow_rates]
        if type(flow_times) in [float, int]:
            flow_times = [flow_times]

        self.flow_rates = np.fromiter(flow_rates, dtype=float)
        self.flow_times = np.fromiter(flow_times, dtype=float)
        assert len(flow_rates) == len(flow_times), 'Flow rate and flow time array lengths do not match.'
        self.num_pump_times = len(flow_rates)

        # If a pump can only be production or injection we can infer the pump type based on the sign of the flow rates (-ve for production, +ve for injection)
        try:
            min_flow = np.min(self.flow_rates[np.nonzero(self.flow_rates)])
        except ValueError:
            raise ValueError('No meaningful flow rate specified (all flow rates are 0).')

        if min_flow > 0: # Injection
            self.injection_well = 1
            if injection_enthalpy is None:
                raise ValueError('No injection enthalpy specified for injection well.')
            self.injection_enthalpy = injection_enthalpy
        else: # Production
            self.injection_well = 0
            self.injection_enthalpy = 0.0

        if deliverability:
            self.deliverability = 1
            self.production_index = production_index
            self.cutoff_pressure = cutoff_pressure
        else:
            self.deliverability = 0
            self.production_index = 0.0
            self.cutoff_pressure = 0.0

        # Units will be useful when displaying pump information in the GUI. 
        # TODO: Display pump information in GUI.
        self.units = {
            'Mass flux' : 'kg/s',
            'Flow time' : 's',
            'Injection Enthalpy' : 'J/kg',
            'Cutoff Pressure' : 'Pa',
            'Production Index' : 'TODO:'
        }


#----------------------------------------------------------------------------------------------------------
# Grid information class 
#----------------------------------------------------------------------------------------------------------
class Grid:
    """
    Holds information used to generate a radial 1d grid. Could be updated to use existing grid
    information (volumes, areas etc).
    """
    def __init__(self, num_blocks=100, num_constant_blocks=20, constant_block_size=0.01, block_growth_factor=1.2):
        """
        Inialise the grid info class

        Inputs:
            num_blocks (integer) - The number of grid blocks to be generated.
            num_constant_blocks (integer) - The number of grid blocks which will have a fixed size.
            constant_block_size (float) - The fixed size which constant blocks will have (meters).
            block_growth_factor (float) - The factor by which non-constant blocks will grow by.
        """
        assert num_blocks >= num_constant_blocks, 'Number of constant blocks cannot exceed the total number of blocks.'
        self.num_blocks = num_blocks
        self.num_constant_blocks = num_constant_blocks
        self.constant_block_size = constant_block_size
        self.block_growth_factor = block_growth_factor