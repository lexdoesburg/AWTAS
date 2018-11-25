# import numpy as np
from numpy import genfromtext

class Data():
    def __init__(self, filename=None, time=None, observed=None, parameters=None):
        self.filename = filename
        self.time = time # list of time data
        self.observed = observed # list of pressure data
        self.approximation = None
        # self.p0 = None
        # self.qm = None
        # self.h = None
        # self.rho = None
        # self.nu = None
        # self.C = None
        # self.r = None
        self.parameters = parameters #list of the parameters
        self.phi = None
        self.k = None
        
        if filename:
            self.read_file()
    
    def read_file(self):
        if self.filename:
            self.time, self.observed = genfromtxt(self.filename, delimiter=',', skip_header=1).T
    
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
    
    def set_pressure(self, pressure):
        self.pressure = pressure
    
    def set_approximation(self, approximation):
        self.approximation = approximation

