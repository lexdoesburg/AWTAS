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
        self.parameters = parameters # list of the well parameters
        self.phi = None # estimated porosity (float)
        self.k = None # estimated permeability (float)
        
        # If a filename is given read the data in
        if filename:
            self.read_file()
    
    def read_file(self, filename=None):
        if filename:
            self.filename = filename
        self.time, self.observation = np.genfromtxt(self.filename, delimiter=',', skip_header=1).T
    
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
    
    def generate_datafile(self, filename):
        with open(filename, 'w') as file:
            file.write('Time, Pressure Observation, Parameters\n')
            j = len(self.parameters)
            for i in range(len(self.time)):
                if j > 0:
                    file.write('{}, {}, {}\n'.format(self.time[i], self.observation[i], self.parameters[i]))
                else:
                    file.write('{}, {}\n'.format(self.time[i], self.observation[i]))
                j -= 1

