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
_observation_point_properties = {
    'deliverability' : 0,
    'pressure' : 1,
    'temperature' : 2,
    'enthalpy' : 3
    }


class ObservationPoints:
    def __init__(self, radial_location=[], property=[], num_data=[]):
        if type(radial_location) is float:
            radial_location = [radial_location]
        if type(property) is str:
            property = [property]
        if type(num_data) is int:
            num_data = [num_data]
        assert len(radial_location) == len(property) == len(num_data), 'Trying to add multiple observation points but insufficent information supplied.'
        
        if property:
            self.property = [_observation_point_properties[p] for p in property]
        else:
            self.property = property

        self.num_observation_points = len(radial_location)
        self.radial_location = np.fromiter(radial_location, dtype=float)
        self.num_data = np.fromiter(num_data, dtype=np.int32)
        self.property = np.fromiter(self.property, dtype=np.int32)
        self.units = {
            'Radial location' : 'm'
        }

    
    def add_observation_points(self, radial_location, property, num_data):
        """
        Add observation point(s). Either takes three lists each in identical order, or three individual values.
        """
        if type(radial_location) is list and type(property) is list and type(num_data) is list:
            assert len(radial_location) == len(property) == len(num_data), 'Trying to add multiple observation points but insufficent information supplied.'
            np.append(self.radial_location, radial_location)
            np.append(self.property, [_observation_point_properties[property_type] for property_type in property])
            np.append(self.num_data, num_data)
        elif type(radial_location) is float and type(property) is str and type(num_data) is int:
            np.append(self.radial_location, radial_location)
            np.append(self.property, _observation_point_properties[property])
            np.append(self.num_data, num_data)


# Type of flow that is occurring (mapping readable names to the fortran equivalent flags)
_pump_schemes = {
    'measured' : 0, # Measured flows (flowrate changes at different times)
    'constant' : 1, # Constant flow (flowrate constant through the whole time period)
    'step' : 4      # Step flow (flow only active for a certain amount of time)
    }


class Pump:
    def __init__(self, pumping_scheme, flow_rates, flow_times, injection_well, injection_enthalpy, deliverability, production_index, cutoff_pressure):
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

        if injection_well:
            self.injection_well = 1
        else:
            self.injection_well = 0

        self.injection_enthalpy = injection_enthalpy

        if deliverability:
            self.deliverability = 1
        else:
            self.deliverability = 0
        
        self.production_index = production_index
        self.cutoff_pressure = cutoff_pressure

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
    # The units are formatted so they will look better when displayed in the GUI
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
    def __init__(self, model_type, filename=None, time=None, observation=None, parameters_list=None, pump_info=None, observation_points=None, grid_info=Grid(), error=None):
        self.model_type = model_type.lower()
        self.time = time # array of time data
        self.observation = observation # array of observed pressure data
        self.approximation = None # array of approximated pressure data
        self.reservoir_conditions, self.fixed_parameters, self.variables = self._initialise_parameter_dictionaries()
        self.pump_info = pump_info
        self.observation_points = observation_points
        self.grid_info = grid_info
        self.initial_x = None # String that states if initial x is temperature or vapour saturation (for use within the GUI)

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
        # TODO Update this work with the new structure

        # self.time, self.observation = np.genfromtxt(self.filename, delimiter=',', skip_header=1).T
        with open(filename, 'r') as file:
            metadata = file.readline()
            metadata = metadata.rstrip().split(',') # Convert string to list of strings
            parameter_values = file.readline() 
            parameter_values = [float(value) for value in parameter_values.split(',')] # Convert string to list of floats
            for parameter_name, value in zip(metadata, parameter_values):
                if parameter_name in self.reservoir_conditions.keys():
                    self.reservoir_conditions[parameter_name]['Value'] = value # This works only if the file uses the same parameter names as the dictionary key.
                elif parameter_name in self.fixed_parameters.keys():
                    self.fixed_parameters[parameter_name]['Value'] = value # This works only if the file uses the same parameter names as the dictionary key.
                else:
                    raise ValueError('Parameter name "{}" not recognised from input data file.'.format(parameter_name))
        self.time, self.observation = np.genfromtxt(filename, delimiter=',', skip_header=4).T

    
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
    

    def create_pump_info(self, pumping_scheme, flow_rates, flow_times, injection_well, injection_enthalpy, deliverability, production_index, cutoff_pressure):
        self.pump_info = Pump(pumping_scheme, flow_rates, flow_times, injection_well, injection_enthalpy, deliverability, production_index, cutoff_pressure)
    

    def create_observation_points(self, radial_location, property, num_data):
        self.observation_points = ObservationPoints(radial_location, property, num_data)


    def create_grid_info(self, num_blocks, num_constant_blocks, constant_block_size, block_growth_factor):
        self.grid_info = Grid(num_blocks, num_constant_blocks, constant_block_size, block_growth_factor)


    def generate_datafile(self, filename, variables=None):
        # TODO Update this work with the new structure
        with open(filename, 'w') as file:
            # Write known well parameters first:
            # file.write('Initial Pressure,Mass Flowrate,Thickness,Density,Kinematic Viscosity,Compressibility,Radius\n')
            for parameter in _model_parameters[self.model_type]['Reservoir Conditions'] + _model_parameters[self.model_type]['Fixed Parameters']:
                if parameter != _model_parameters[self.model_type]['Fixed Parameters'][-1]:
                    file.write('{},'.format(parameter))
                else:
                    file.write('{}\n'.format(parameter))
            for parameter in _model_parameters[self.model_type]['Reservoir Conditions'] + _model_parameters[self.model_type]['Fixed Parameters']:
                if parameter != _model_parameters[self.model_type]['Fixed Parameters'][-1]:
                    if parameter in _model_parameters[self.model_type]['Reservoir Conditions']:
                        file.write('{},'.format(self.reservoir_conditions[parameter]['Value']))
                    else:
                        file.write('{},'.format(self.fixed_parameters[parameter]['Value']))
                else:
                    if parameter in _model_parameters[self.model_type]['Reservoir Conditions']:
                        file.write('{},'.format(self.reservoir_conditions[parameter]['Value']))
                    else:
                        file.write('{},'.format(self.fixed_parameters[parameter]['Value']))
                    if variables:
                        file.write('Known Porosity : {}, Known Permeability : {}, '.format(variables[0], variables[1]))
                        if self.error:
                            file.write('Standard Deviation of Errors : {}\n'.format(self.error))
                        else:
                            file.write('Standard Deviation of Errors : Unknown\n')
                    else:
                        file.write('\n')

            # for parameter in self.parameters:
            #     if parameter is not self.parameters[-1]:
            #         file.write('{},'.format(parameter)) # Comma-separated values
            #     else:
            #         file.write('{}\n'.format(parameter)) # Add a blank line separation.
            #         if variables:
            #             file.write('-------- Actual Porosity = {} --- Actual Permeability = {} --------\n'.format(variables[0], variables[1]))
            #         else:
            #             file.write('\n')

            # Write the time of observation and observed pressure readings:
            file.write('Time [s],Pressure Observation [Pa]\n')
            # j = len(self.parameters)
            for i in range(len(self.time)):
                # if j > 0:
                #     file.write('{}, {}, {}\n'.format(self.time[i], self.observation[i], self.parameters[i]))
                # else:
                file.write('{},{}\n'.format(self.time[i], self.observation[i]))
                # j -= 1