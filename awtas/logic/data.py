import numpy as np


class Grid:
    """
    Holds information about the grid. Could be updated to use existing grids (volumes, areas etc).
    """
    def __init__(self, num_blocks=100, num_constant_blocks=20, constant_block_size=0.01, block_growth_factor=1.2):
        assert num_blocks >= num_constant_blocks, 'Number of constant blocks cannot exceed the total number of blocks.'
        self.num_blocks = num_blocks
        self.num_constant_blocks = num_constant_blocks
        self.constant_block_size = constant_block_size
        self.block_growth_factor = block_growth_factor


# Type of measurement data that was recorded at the observation point (mapping readable names to the fortran equivalent flags)
_observation_properties = {
    'Deliverability' : 0,
    'Pressure' : 1,
    'Temperature' : 2,
    'Enthalpy' : 3
    }

_observation_property_from_flag = {
    # Opposite ordering of the above dictionary to get meaningful name from fortran flag
    0 : 'Deliverability',
    1 : 'Pressure',
    2 : 'Temperature',
    3 : 'Enthalpy'
}

"""
Observation point class isn't highly intuitive. Would be better to have a class and each observation point object 
is a single observation point. Currently the class is for multiple observation points.
This way was the quickest way to work best with the fortran code currently - can be improved later to use single object instance of
the class for each observation point and then combine the required information from each obs point to be passed to fortran.
"""
class ObservationPoints:
    def __init__(self, radial_location=[], property=[], num_data=[], times=[], observations=[]):
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
        assert len(radial_location) == len(property) == len(num_data) == len(times) == len(observations), 'Trying to add multiple observation points but insufficent information supplied.'
        
        if property:
            self.property = [_observation_properties[p] for p in property]
        else:
            self.property = property

        self.num_observation_points = len(radial_location)
        self.radial_location = np.fromiter(radial_location, dtype=float)
        self.num_data = np.fromiter(num_data, dtype=np.int32)
        self.property = np.fromiter(self.property, dtype=np.int32)
        self.times = times
        self.observations = observations
        self.modelled_values = []

        self.units = {
            'Radial location' : 'm'
        }

    
    def add_observation_points(self, radial_location, property, num_data, times, observations):
        """
        Add observation point(s). Either takes three lists each in identical order, or three individual values.
        """
        if type(radial_location) is list and type(property) is list and type(num_data) is list and type(times) is list and type(observations) is list:
            assert len(radial_location) == len(property) == len(num_data) == len(times) == len(observations), 'Trying to add multiple observation points but insufficent information supplied.'
            np.append(self.radial_location, radial_location)
            np.append(self.property, [_observation_properties[property_type] for property_type in property])
            np.append(self.num_data, num_data)
            self.times.extend(times)
            self.observations.extend(observations)
        elif type(radial_location) is float and type(property) is str and type(num_data) is int:
            np.append(self.radial_location, radial_location)
            np.append(self.property, _observation_properties[property])
            np.append(self.num_data, num_data)
            assert len(times) == num_data == len(observations), 'Given {} observation times and {} observation values when expecting {}.'.format(len(times), len(observations), num_data)
            self.times.append(times)
            self.observations.append(observations)
    
    def store_modelled_values(self, modelled_values):
        current_index = 0
        for i in range(self.num_observation_points):
            next_index = current_index + self.num_data[i]
            self.modelled_values.append(np.fromiter(modelled_values[current_index:next_index], dtype=float))
            current_index = next_index


# Type of flow that is occurring (mapping readable names to the fortran equivalent flags)
_pump_schemes = {
    'Measured Flows' : 0, # Flowrate changes at different times
    'Constant Flow' : 1,  # Flowrate constant through the whole time period
    'Step Flows' : 4      # Flow only active for a certain amount of time
    }

_pump_scheme_from_flag = {
    # Opposite ordering of the above dictionary to get meaningful name from fortran flag
    0 : 'Measured Flows', # Flowrate changes at different times
    1 : 'Constant Flow',  # Flowrate constant through the whole time period
    4 : 'Step Flows'      # Flow only active for a certain amount of time
    }


class Pump:
    def __init__(self, pumping_scheme, flow_rates, flow_times, deliverability, production_index, cutoff_pressure, injection_well=None, injection_enthalpy=None):
        self.pumping_scheme = _pump_schemes[pumping_scheme]

        # If input rates and times are not iterable, convert to iterables
        if type(flow_rates) is float:
            flow_rates = [flow_rates]
        if type(flow_times) is float:
            flow_times = [flow_times]

        self.flow_rates = np.fromiter(flow_rates, dtype=float)
        self.flow_times = np.fromiter(flow_times, dtype=float)
        assert len(flow_rates) == len(flow_times), 'Flow rate and flow time array lengths do not match.'
        self.num_pump_times = len(flow_rates)

        # if injection_well: 
        #     self.injection_well = 1
        #     self.injection_enthalpy = injection_enthalpy
        # else:
        #     self.injection_well = 0
        #     self.injection_enthalpy = 0.0

        # If a pump can only be production or injection we can infer the pump type based on the sign of the flow rates (-ve for production, +ve for injection)
        try:
            min_flow = np.min(self.flow_rates[np.nonzero(self.flow_rates)])
        except ValueError:
            raise ValueError('No meaningful flowrate specified (all flowrates are 0).')

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

        self.units = {
            'Mass flux' : 'kg/s',
            'Flow time' : 's',
            'Injection Enthalpy' : 'J',
            'Cutoff Pressure' : 'Pa',
            'Production Index' : 'TODO'
        }


# Parameters used in each model - the key is the model_type (update when new models are added)
_model_parameters = {
    'theis' : 
        {
            'Reservoir Conditions' : ['Initial Pressure'],
            'Fixed Parameters' : ['Mass Flowrate', 'Layer Thickness', 'Density', 'Kinematic Viscosity', 'Compressibility', 'Action Well Radius'],
            'Variables' : ['Porosity', 'Permeability']
        },
    'radial1d' :
        {
            'Reservoir Conditions' : ['Initial Pressure', 'Initial X'],
            'Fixed Parameters' : ['Layer Thickness', 'Action Well Radius', 'Rock Specific Heat', 'Rock Heat Conductivity', 'Rock Density', 'Rock Compressibility'],
            'Variables' : ['Porosity', 'Permeability']
        }
    }


# Default units for each parameter. Can look to implement changing units via combo box later.
_default_parameter_units = {
    # The units are formatted to look better when displayed in the GUI
    'Reservoir Conditions' : 
        {
        'Initial Pressure' : {'Units':'Pa'},
        'Initial X' : {'Units': {'Initial Temperature':'°C', 
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


class Data():
    """
    This is a class which defines the data structure which contains (eventually extended for use in all model types).
    """
    def __init__(self, model_type, filename=None, time=None, observation=None, parameters_list=None, pump_info=None, observation_points=None, grid_info=None, error=None):
        self.model_type = model_type.lower()
        assert self.model_type_valid(), 'Input model type is invalid try one of {} instead.'.format(_model_parameters.keys())
        self.approximation = None # array of approximated pressure data
        self.reservoir_conditions, self.fixed_parameters, self.variables = self._initialise_parameter_dictionaries()
        self.pump_info = pump_info
        self.observation_points = observation_points
        self.grid_info = grid_info
        self.initial_x = None # String that states if initial x is temperature or vapour saturation (for use within the GUI)

        if self.model_type is not 'theis' and self.observation_points:
            # if self.observation_points.num_observation_points == 1:
            #     self.time = self.observation_points.times
            #     self.observation = self.observation_points.observations
            # else: # Otherwise multiple observation points.
            self.time = np.concatenate(self.observation_points.times)
            self.observation = np.concatenate(self.observation_points.observations)
        else:
            self.time = time # array of time data
            self.observation = observation # array of observed pressure data

        if model_type == 'radial1d' and observation_points:
            self.total_num_data = sum(observation_points.num_data)
        elif model_type == 'theis' and time:
            self.total_num_data = len(time)
        else:
            self.total_num_data = 0

        if parameters_list:
            self._fill_parameter_dictionaries(parameters_list)
        
        self.error = error # Estimated error in the readings # Move the error to observation points (each could have different errors

        # If a filename is given read the data in
        if filename:
            self.read_file(filename)
    

    def _initialise_parameter_dictionaries(self):
        reservoir_conditions = {}
        fixed_parameters = {}
        variables = {}
        for reservoir_condition in _model_parameters[self.model_type]['Reservoir Conditions']:
            reservoir_conditions[reservoir_condition] = {'Value':None, 'Units':_default_parameter_units['Reservoir Conditions'][reservoir_condition]['Units']}
        for parameter in _model_parameters[self.model_type]['Fixed Parameters']:
            fixed_parameters[parameter] = {'Value':None, 'Units':_default_parameter_units['Fixed Parameters'][parameter]['Units']}
        for variable in _model_parameters[self.model_type]['Variables']:
            variables[variable] = {'Value':None, 'Units':_default_parameter_units['Variables'][variable]['Units']}
        return reservoir_conditions, fixed_parameters, variables


    def _fill_parameter_dictionaries(self, parameters_list):        
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
                if self.reservoir_conditions['Initial X']['Value'] < 1.0:
                    self.initial_x = 'Initial Vapour Saturation'
                else:
                    self.initial_x = 'Initial Temperature'
               

    def read_file(self, filename):
        # TODO Update theis solution code to work with the improved code used for radial1d
        if self.model_type == 'theis':
            # Old code
            with open(filename, 'r') as file:
                # Read error
                file.readline() # Skip first line
                info_values = file.readline()
                info_values = info_values.rstrip().split(',') # Convert string to list of strings
                try:
                    self.error = float(info_values[0])
                    if self.error <= 1e-6:
                        self.error = None
                except ValueError:
                    # If the error can't be converted to a float then don't store any error value
                    self.error = None
                print('Read error value = {}'.format(self.error))
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
            # New code for homogeneous porous (radial1d) and later models
            obs_point_properties = []
            obs_point_locations = []
            obs_point_num_data = []
            obs_point_times = []
            obs_point_observations = []
            with open(filename, 'r') as file:
                for counter, line in enumerate(file):
                # for line in file:
                    print('Line {}: {}'.format(counter, line))
                    if 'KNOWN WELL PARAMETERS' in line:
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
                    elif 'GRID INFORMATION' in line:
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
                        self.grid_info = Grid(num_blocks=num_blocks, num_constant_blocks=num_constant_blocks, constant_block_size=constant_block_size, block_growth_factor=block_growth_factor)
                    elif 'PUMP INFORMATION' in line:
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
                        file.readline()
                        file.readline()
                        pump_times = []
                        pump_rates = []
                        for i in range(num_pump_times):
                            flow_info = file.readline()
                            flow_info = [float(value) for value in flow_info.split(',')]
                            pump_times.append(flow_info[0])
                            pump_rates.append(flow_info[1])
                        self.pump_info = Pump(pumping_scheme=pumping_scheme,flow_rates=pump_rates, flow_times=pump_times, deliverability=deliverability, production_index=production_index, cutoff_pressure=cutoff_pressure, injection_enthalpy=injection_enthalpy)
                    elif 'OBSERVATION POINT' in line:
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
                        file.readline()
                        file.readline()
                        times = np.ndarray(shape=num_data, dtype=float)
                        observations = np.ndarray(shape=num_data, dtype=float)
                        for i in range(num_data):
                            observation_info = file.readline()
                            observation_info = [float(value) for value in observation_info.split(',')]
                            times[i] = observation_info[0]
                            observations[i] = observation_info[1]
                        obs_point_times.append(times)
                        obs_point_observations.append(observations)

            self.observation_points = ObservationPoints(radial_location=obs_point_locations, property=obs_point_properties, num_data=obs_point_num_data, times=obs_point_times, observations=obs_point_observations)                
            self.time = np.concatenate(self.observation_points.times)
            self.observation = np.concatenate(self.observation_points.observations)
            self.total_num_data = len(self.time)

    def set_known_parameters(self, parameters):
        """
        parameters = list of parameters

        Analytical Theis Model List Ordering:
            parameters[0] = intial pressure
            parameters[1] = mass flowrate
            parameters[2] = thickness
            parameters[3] = density
            parameters[4] = kinematic viscosity
            parameters[5] = compressibility
            parameters[6] = action well radius

        Radial 1D Model List Ordering:
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
        self.variables = variables


    def set_time(self, time):
        if self.observation and len(self.observation) == len(time):
            self.total_num_data = len(time)
        else:
            print('Number of recorded observations is not equal to the number of recorded times.')
        self.time = time


    def set_observation(self, observation):
        self.observation = observation


    def set_approximation(self, approximation):
        self.approximation = approximation


    def set_error(self, error):
        self.error = error
    

    def create_pump_info(self, pumping_scheme, flow_rates, flow_times, deliverability, production_index, cutoff_pressure, injection_well, injection_enthalpy):
        self.pump_info = Pump(pumping_scheme, flow_rates, flow_times, deliverability, production_index, cutoff_pressure, injection_well, injection_enthalpy)
    

    def create_observation_points(self, radial_location, property, num_data, times, observations):
        self.observation_points = ObservationPoints(radial_location, property, num_data, times, observations)


    def create_grid_info(self, num_blocks, num_constant_blocks, constant_block_size, block_growth_factor):
        self.grid_info = Grid(num_blocks, num_constant_blocks, constant_block_size, block_growth_factor)


    def model_type_valid(self):
        if self.model_type in _model_parameters.keys():
            return True
        else:
            return False

    def generate_datafile(self, filename, variables=None):
        # TODO Update theis solution code to work with the improved code used for radial1d model
        if self.model_type == 'theis':
            # Old code
            with open(filename, 'w') as file:
                # Write known errors and variable values
                info_labels = 'SD of Error,'
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
                # if variables:
                #     for i, (variable, value) in enumerate(variables.items()):

                #         if i != len(variables) - 1:
                #             file.write('{} : {}, '.format(variable, value))
                #         else:
                #             file.write('{} : {}'.format(variable, value))
                #             if self.error:
                #                 file.write(', Known standard deviation of error : {} [Pa]\n'.format(self.error))
                #             else:
                #                 file.write('\n')
                # else:
                #     if self.error:
                #         file.write('Variable values unknown, Known standard deviation of error : {} [Pa]\n'.format(self.error))
                #     else:
                #         file.write('\n')
                
                # Write the time of observation and observed pressure readings:
                file.write('Time [s],Pressure Observation [Pa]\n')
                for i in range(len(self.time)):
                    file.write('{},{}\n'.format(self.time[i], self.observation[i]))
                # for parameter in _model_parameters[self.model_type]:
                #     if parameter != _model_parameters[self.model_type][-1]:
                #         file.write('{},'.format(parameter))
                #     else:
                #         file.write('{}\n'.format(parameter))
                # for parameter in _model_parameters[self.model_type]:
                #     if parameter != _model_parameters[self.model_type][-1]:
                #         file.write('{},'.format(self.parameters[parameter]['Value']))
                #     else:
                #         file.write('{}\n'.format(self.parameters[parameter]['Value']))
                #         if variables:
                #             file.write('Known Porosity : {}, Known Permeability : {}, '.format(variables[0], variables[1]))
                #             if self.error:
                #                 file.write('Standard Deviation of Errors : {}\n'.format(self.error))
                #             else:
                #                 file.write('Standard Deviation of Errors : Unknown\n')
                #         else:
                #             file.write('\n')
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
                    observation_point_labels = 'Property Observed,Radial Location [m],Number of Observations\n'
                    observation_point_values = '{},{},{}\n'.format(_observation_property_from_flag[self.observation_points.property[i]],
                                                                self.observation_points.radial_location[i],
                                                                self.observation_points.num_data[i])
                    file.write(observation_point_labels)
                    file.write(observation_point_values)
                    file.write('\n')

                    # Write all observations
                    file.write('Time [s],Observation\n')
                    for j in range(self.observation_points.num_data[i]):
                        file.write('{},{}\n'.format(self.observation_points.times[i][j], self.observation_points.observations[i][j]))
                    if i != self.observation_points.num_observation_points - 1: # Write a new line between each observation point
                        file.write('\n')
