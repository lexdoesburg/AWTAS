import numpy as np
 
class Data():
    """
    This is a class which defines the problem data structure (eventually extended for use in all model types).
    """
    def __init__(self, filename=None, time=None, observation=None, parameters=None):
        self.filename = filename # the filename of where the data was imported from
        self.time = time # array of time data
        self.observation = observation # array of observed pressure data
        self.approximation = None # array of approximated pressure data
        if parameters:
            self.parameters = parameters # list of the well parameters
        else:
            self.parameters = [None]*7
        self.phi = None # estimated porosity (float)
        self.k = None # estimated permeability (float)
        
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
        self.parameters = parameters

    def set_unknown_parameters(self, phi, k):
        self.phi = phi
        self.k = k

    def set_time(self, time):
        self.time = time
    
    def set_observation(self, observation):
        self.observation = observation
    
    def set_approximation(self, approximation):
        self.approximation = approximation
    
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

