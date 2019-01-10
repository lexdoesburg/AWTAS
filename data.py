import numpy as np

def create_data(model_type, filename=None, time=None, observation=None, parameters=None, error=None):
    model_type = model_type.lower()
    if model_type == 'theis':
        data = Data(filename, time, observation, parameters, error)
    elif model_type == 'radial1d':
        data = Radial1d_Data(filename, time, observation, parameters, error)
    return data

class DataPoint():
    def __init__(self, time, observation, approximation=None, error=None, weight=None):
        self.time = time
        self.observation = observation
        self.approximation = approximation
        self.error = error
        self.weight = weight

class Data():
    """
    This is a class which defines the problem data structure (eventually extended for use in all model types).
    """
    __model_parameters = {
                         'theis' : ['Initial Pressure', 'Mass Flowrate', 'Layer Thickness', 'Density', 'Kinematic Viscosity', 'Compressibility', 'Radius'],
                         'radial1d' : ['Initial Pressure', 'Mass Flowrate', 'Layer Thickness', 'Initial Temperature', 'Initial Vapour Saturation', 'Action Well Radius', 'Rock Specific Heat', 'Rock Conductivity', 'Rock Density', 'Rock Compressibility', 'Observation Point Distance']
                         }
    
    __default_parameter_units = {
                                'Initial Pressure' : {'Units':'Pa', 'Models':['theis', 'radial1d']},
                                'Mass Flowrate' : {'Units':'Kg/s', 'Models':['theis', 'radial1d']},
                                'Layer Thickness' : {'Units':'m', 'Models':['theis', 'radial1d']},
                                'Density' : {'Units':'Kg/m<sup>3</sup>', 'Models':['theis']},
                                'Kinematic Viscosity' : {'Units':'m<sup>2</sup>/s', 'Models':['theis']},
                                'Compressibility' : {'Units':'1/Pa', 'Models':['theis']},
                                'Radius' : {'Units':'m', 'Models':['theis']},
                                'Initial Temperature' : {'Units':'TODO', 'Models':['radial1d']},
                                'Initial Vapour Pressure' : {'Units':'TODO', 'Models':['radial1d']},
                                'Action Well Radius' : {'Units':'TODO', 'Models':['radial1d']},
                                'Rock Specific Heat' : {'Units':'TODO', 'Models':['radial1d']},
                                'Rock Conductivity' : {'Units':'TODO', 'Models':['radial1d']},
                                'Rock Density' : {'Units':'TODO', 'Models':['radial1d']},
                                'Rock Compressibility' : {'Units':'TODO', 'Models':['radial1d']},
                                'Observation Point Distance' : {'Units':'TODO', 'Models':['radial1d']}
                                }

    def __init__(self, filename=None, time=None, observation=None, parameters_list=None, error=None, model_type='Theis'):
        self.model_type = model_type.lower()
        self.filename = filename # the filename of where the data was imported from
        self.time = time # array of time data
        self.observation = observation # array of observed pressure data
        self.approximation = None # array of approximated pressure data
        self.parameters = self._initialise_parameter_dictionary()
        if parameters_list:
            self.fill_parameter_dictionary(parameters_list)
        self.variables = [None]*2
        
        self.error = error # Estimated error in the readings

        # If a filename is given read the data in
        if filename:
            self.read_file()
    
    def _initialise_parameter_dictionary(self):
        # if hasattr(self, 'parameters'):
        #     raise ValueError('Parameter dictionary has already been initialised')
        # else:
        #     self.parameters = {'Initial Pressure' : {'Value':None, 'Units':'Pa'},
        #                     'Mass Flowrate' : {'Value':None, 'Units':'Kg/s'},
        #                     'Layer Thickness' : {'Value':None, 'Units':'m'},
        #                     'Density' : {'Value':None, 'Units':'Kg/m3'},
        #                     'Kinematic Viscosity' : {'Value':None, 'Units':'m2/s'},
        #                     'Compressibility' : {'Value':None, 'Units':'1/Pa'},
        #                     'Radius' : {'Value':None, 'Units':'m'},
        #                     'Initial Temperature' : {'Value':None, 'Units':'TODO'},
        #                     'Initial Vapour Saturation' : {'Value':None, 'Units':'TODO'},
        #                     'Action Well Radius' : {'Value':None, 'Units':'TODO'},
        #                     'Rock Specific Heat' : {'Value':None, 'Units':'TODO'},
        #                     'Rock Conductivity' : {'Value':None, 'Units':'TODO'},
        #                     'Rock Density' : {'Value':None, 'Units':'TODO'},
        #                     'Rock Compressibility' : {'Value':None, 'Units':'TODO'},
        #                     'Observation Point Distance' : {'Value':None, 'Units':'TODO'}}
        parameters = {}
        for parameter in self.__model_parameters[self.model_type]:
            parameters[parameter] = {'Value':None, 'Units':self.__default_parameter_units[parameter]['Units']}
        return parameters


    def fill_parameter_dictionary(self, parameters_list):        
        if parameters_list:
            if self.model_type == 'theis':
                self.parameters['Initial Pressure']['Value'] = parameters_list[0]
                self.parameters['Mass Flowrate']['Value'] = parameters_list[1]
                self.parameters['Layer Thickness']['Value'] = parameters_list[2]
                self.parameters['Density']['Value'] = parameters_list[3]
                self.parameters['Kinematic Viscosity']['Value'] = parameters_list[4]
                self.parameters['Compressibility']['Value'] = parameters_list[5]
                self.parameters['Radius']['Value'] = parameters_list[6]
            elif self.model_type == 'radial1d':
                self.parameters['Initial Pressure']['Value'] = parameters_list[0]
                self.parameters['Initial Temperature']['Value'] = parameters_list[1]
                self.parameters['Action Well Radius']['Value'] = parameters_list[2]
                self.parameters['Layer Thickness']['Value'] = parameters_list[3]
                self.parameters['Rock Specific Heat']['Value'] = parameters_list[4]
                self.parameters['Rock Conductivity']['Value'] = parameters_list[5]
                self.parameters['Rock Density']['Value'] = parameters_list[6]
                self.parameters['Rock Compressibility']['Value'] = parameters_list[7]
                self.parameters['Mass Flowrate']['Value'] = parameters_list[8]
                self.parameters['Observation Point Distance']['Value'] = parameters_list[9]


    def read_file(self, filename=None):
        if filename:
            self.filename = filename
        # self.time, self.observation = np.genfromtxt(self.filename, delimiter=',', skip_header=1).T
        with open(self.filename, 'r') as file:
            metadata = file.readline()
            metadata = metadata.rstrip().split(',') # Convert string to list of strings
            parameter_values = file.readline() 
            parameter_values = [float(value) for value in parameter_values.split(',')] # Convert string to list of floats
            for parameter_name, value in zip(metadata, parameter_values):
                parameter_name = parameter_name.lower()
                if parameter_name == "initial pressure" or parameter_name == "p0":
                    self.parameters[0] = value
                elif parameter_name == "mass flowrate" or parameter_name == "qm":
                    self.parameters[1] = value
                elif parameter_name == "thickness" or parameter_name == "h":
                    self.parameters[2] = value
                elif parameter_name == "density" or parameter_name == "rho":
                    self.parameters[3] = value
                elif parameter_name == "kinematic viscosity" or parameter_name == "nu":
                    self.parameters[4] = value
                elif parameter_name == "compressibility" or parameter_name == "c":
                    self.parameters[5] = value
                elif parameter_name == "radius" or parameter_name == "r":
                    self.parameters[6] = value
        
        self.time, self.observation = np.genfromtxt(self.filename, delimiter=',', skip_header=4).T

    
    def set_known_parameters(self, parameters):
        """
        parameters = list of parameters
        parameters[0] = intial pressure
        parameters[1] = mass flowrate
        parameters[2] = thickness
        parameters[3] = density
        parameters[4] = kinematic viscosity
        parameters[5] = compressibility
        parameters[6] = radius
        """
        # self.parameters = parameters
        self.fill_parameter_dictionary()

    def set_unknown_parameters(self, variables):
        # self.phi = phi
        # self.k = k
        self.variables = variables

    def set_time(self, time):
        self.time = time
    
    def set_observation(self, observation):
        self.observation = observation
    
    def set_approximation(self, approximation):
        self.approximation = approximation

    def set_error(self, error):
        self.error = error
    
    def generate_datafile(self, filename, variables=None):
        with open(filename, 'w') as file:
            # Write known well parameters first:
            file.write('Initial Pressure,Mass Flowrate,Thickness,Density,Kinematic Viscosity,Compressibility,Radius\n')
            for parameter in self.parameters:
                if parameter is not self.parameters[-1]:
                    file.write('{},'.format(parameter)) # Comma-separated values
                else:
                    file.write('{}\n'.format(parameter)) # Add a blank line separation.
                    if variables:
                        file.write('-------- Actual Porosity = {} --- Actual Permeability = {} --------\n'.format(variables[0], variables[1]))
                    else:
                        file.write('\n')
            # Write the time of observation and observed pressure readings:
            file.write('Time (s),Pressure Observation (Pa)\n')
            # j = len(self.parameters)
            for i in range(len(self.time)):
                # if j > 0:
                #     file.write('{}, {}, {}\n'.format(self.time[i], self.observation[i], self.parameters[i]))
                # else:
                file.write('{},{}\n'.format(self.time[i], self.observation[i]))
                # j -= 1


class Radial1d_Data():
    """
    This is a class which defines the problem data structure (eventually extended for use in all model types).
    """
    def __init__(self, filename=None, time=None, observation=None, parameters=None, error=None):
        self.parameter_names = ['Initial Pressure', 'Temperature', 'Action Well Radius', 'Layer Thickness', 'Rock Specific Heat', 'Rock Conductivity', 'Rock Density', 'Rock Compressibility', 'Flow Rate', 'Observation Point Distance']
        self.filename = filename # the filename of where the data was imported from
        self.time = time # array of time data
        self.observation = observation # array of observed pressure data
        self.approximation = None # array of approximated pressure data
        if parameters:
            self.parameters = parameters # list of the well parameters
        else:
            self.parameters = [None]*10
        # self.phi = None # estimated porosity (float)
        # self.k = None # estimated permeability (float)
        self.variables = [None]*2
        
        self.error = error # Estimated error in the readings

        # If a filename is given read the data in
        if filename:
            self.read_file()
    
    def read_file(self, filename=None):
        if filename:
            self.filename = filename
        # self.time, self.observation = np.genfromtxt(self.filename, delimiter=',', skip_header=1).T
        with open(self.filename, 'r') as file:
            metadata = file.readline()
            metadata = metadata.rstrip().split(',') # Convert string to list of strings
            parameter_values = file.readline() 
            parameter_values = [float(value) for value in parameter_values.split(',')] # Convert string to list of floats
            for parameter_name, value in zip(metadata, parameter_values):
                parameter_name = parameter_name.lower()
                if parameter_name == "initial pressure" or parameter_name == "p0":
                    self.parameters[0] = value
                elif parameter_name in ["temperature", "t", "vapour saturation"]:
                    self.parameters[1] = value
                elif parameter_name == "action well radius" or parameter_name == "r":
                    self.parameters[2] = value
                elif parameter_name == "layer thickness" or parameter_name == "h":
                    self.parameters[3] = value
                elif parameter_name == "rock specific heat" or parameter_name == "cr":
                    self.parameters[4] = value
                elif parameter_name == "rock conductivity" or parameter_name == "cond":
                    self.parameters[5] = value
                elif parameter_name == "rock density" or parameter_name == "rhor":
                    self.parameters[6] = value
                elif parameter_name == "rock compressibility" or parameter_name == "comp":
                    self.parameters[7] = value
                elif parameter_name == "flow rate" or parameter_name == "q":
                    self.parameters[8] = value
                elif parameter_name == "observation point distance" or parameter_name == "dist":
                    self.parameters[9] = value
        
        self.time, self.observation = np.genfromtxt(self.filename, delimiter=',', skip_header=4).T

    
    def set_known_parameters(self, parameters):
        """
        [p0, X0, rw, thick, CR, COND, RHOR, COMP, ConstRate, distFromWell]
        parameters = list of parameters
        parameters[0] = intial pressure
        parameters[1] = temperature / vapour saturation
        parameters[2] = action well radius
        parameters[3] = layer thickness
        parameters[4] = rock specific heat
        parameters[5] = rock heat conductivity
        parameters[6] = rock density
        parameters[7] = rock compressibility
        parameters[8] = flow rate
        parameters[9] = observation point distance from well
        """
        self.parameters = parameters

    def set_unknown_parameters(self, variables):
        # self.phi = phi
        # self.k = k
        self.variables = variables

    def set_time(self, time):
        self.time = time
    
    def set_observation(self, observation):
        self.observation = observation
    
    def set_approximation(self, approximation):
        self.approximation = approximation

    def set_error(self, error):
        self.error = error
    
    def generate_datafile(self, filename, variables=None):
        with open(filename, 'w') as file:
            # Write known well parameters first:
            file.write('Initial Pressure,Temperature,Action Well Radius,Layer Thickness,Rock Specific Heat,Rock Conductivity,Rock Density,Rock Compressibility,Flow Rate,Observation Point Distance\n')
            j = 0
            for parameter in self.parameters:
                if j < 9:
                    file.write('{},'.format(parameter)) # Comma-separated values
                else:
                    file.write('{}\n'.format(parameter)) # Add a blank line separation.
                    if variables:
                        file.write('-------- Actual Porosity = {} --- Actual Permeability = {} --------\n'.format(variables[0], variables[1]))
                    else:
                        file.write('\n')
                j = j+1
            # Write the time of observation and observed pressure readings:
            file.write('Time (s),Pressure Observation (Pa)\n')
            # j = len(self.parameters)
            if self.time != None and self.observation != None:
                for i in range(len(self.time)):
                    # if j > 0:
                    #     file.write('{}, {}, {}\n'.format(self.time[i], self.observation[i], self.parameters[i]))
                    # else:
                    file.write('{},{}\n'.format(self.time[i], self.observation[i]))
                    # j -= 1
            else:
                file.write('None,None')


class Theis_Data(Data):
    def __init__(self, filename=None, time=None, observation=None, parameters=[None]*7, error=None):
        self.filename = filename # the filename of where the data was imported from
        self.time = time # array of time data
        self.observation = observation # array of observed pressure data
        self.approximation = None # array of approximated pressure data

        p0_dict = {'Value':parameters[0], 'Units':'Pa'}
        qm_dict = {'Value':parameters[1], 'Units':'Kg/s'}
        h_dict = {'Value':parameters[2], 'Units':'m'}
        rho_dict = {'Value':parameters[3], 'Units':'Kg/m3'}
        nu_dict = {'Value':parameters[4], 'Units':'m2/s'}
        c_dict = {'Value':parameters[5], 'Units':'1/Pa'}
        r_dict = {'Value':parameters[6], 'Units':'m'}

        self.parameters = {'Initial Pressure':p0_dict, 'Mass Flowrate':qm_dict, 'Layer Thickness':h_dict, 'Density':rho_dict,
              'Kinematic Viscosity':nu_dict, 'Compressibility':c_dict, 'Radius':r_dict}

        self.variables = dict.fromkeys(['Initial Pressure'])        
        # self.variables = dict.fromkeys(['Porosity', 'Permeability'])
        self.variables = {'Porosity':{'Value':None, 'Units':''}, 'Permeability':{'Value':None, 'Units':'m2'}}

        self.error = error # Estimated error in the readings

        # If a filename is given read the data in
        if filename:
            self.read_file()

    def read_file(self, filename=None):
        if filename:
            self.filename = filename
        # self.time, self.observation = np.genfromtxt(self.filename, delimiter=',', skip_header=1).T
        with open(self.filename, 'r') as file:
            metadata = file.readline()
            metadata = metadata.rstrip().split(',') # Convert string to list of strings
            parameter_values = file.readline() 
            parameter_values = [float(value) for value in parameter_values.split(',')] # Convert string to list of floats
            for parameter_name, value in zip(metadata, parameter_values):
                parameter_name = parameter_name.lower()
                if parameter_name == "initial pressure" or parameter_name == "p0":
                    self.parameters['Initial Pressure']['Value'] = value
                elif parameter_name == "mass flowrate" or parameter_name == "qm":
                    self.parameters['Mass Flowrate']['Value'] = value
                elif parameter_name == "thickness" or parameter_name == "h":
                    self.parameters['Layer Thickness']['Value'] = value
                elif parameter_name == "density" or parameter_name == "rho":
                    self.parameters['Density']['Value'] = value
                elif parameter_name == "kinematic viscosity" or parameter_name == "nu":
                    self.parameters['Kinematic Viscosity']['Value'] = value
                elif parameter_name == "compressibility" or parameter_name == "c":
                    self.parameters['Compressibility']['Value'] = value
                elif parameter_name == "radius" or parameter_name == "r":
                    self.parameters['Radius']['Value'] = value
        
        self.time, self.observation = np.genfromtxt(self.filename, delimiter=',', skip_header=4).T

    
    def set_known_parameters(self, parameters):
        """
        parameters = list of parameters
        parameters[0] = intial pressure
        parameters[1] = mass flowrate
        parameters[2] = thickness
        parameters[3] = density
        parameters[4] = kinematic viscosity
        parameters[5] = compressibility
        parameters[6] = radius
        """
        self.parameters = parameters

    def set_unknown_parameters(self, variables):
        # self.phi = phi
        # self.k = k
        self.variables['Porosity']['Value'], self.variables['Permeability']['Value'] = variables

    def set_time(self, time):
        self.time = time
    
    def set_observation(self, observation):
        self.observation = observation
    
    def set_approximation(self, approximation):
        self.approximation = approximation

    def set_error(self, error):
        self.error = error
    
    def generate_datafile(self, filename, variables=None):
        with open(filename, 'w') as file:
            # Write known well parameters first:
            # file.write('Initial Pressure,Mass Flowrate,Thickness,Density,Kinematic Viscosity,Compressibility,Radius\n')
            i = 0
            for parameter in self.parameters.keys():
                if i < len(self.parameters)-1:
                    file.write('{},'.format(parameter)) # Comma-separated values
                else:
                    file.write('{}\n'.format(parameter)) # Add a blank line separation.
                i += 1
            i = 0
            for value in self.parameters.values():
                if i < len(self.parameters)-1:
                    file.write('{},'.format(value)) # Comma-separated values
                else:
                    file.write('{}\n'.format(value))
                    if variables:
                        file.write('-------- Actual Porosity = {} --- Actual Permeability = {} --------\n'.format(variables['Porosity']['Value'], variables['Permeability']['Value']))
                    else:
                        file.write('\n')
            # Write the time of observation and observed pressure readings:
            file.write('Time (s),Pressure Observation (Pa)\n')
            # j = len(self.parameters)
            for i in range(len(self.time)):
                # if j > 0:
                #     file.write('{}, {}, {}\n'.format(self.time[i], self.observation[i], self.parameters[i]))
                # else:
                file.write('{},{}\n'.format(self.time[i], self.observation[i]))
                # j -= 1