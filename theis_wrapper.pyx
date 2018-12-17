from numpy cimport ndarray
from numpy import empty

cdef extern from "theis_wrapper.h":
    void c_theis(double* k, double* nu, double* phi, double* rho, double* c, double* b, double* Q0, double* P0, double* r, int* numData, double* time, double* pressure)

def theis(double k, double nu, double phi, double rho, double c, double b, double Q0, double P0, double r, int numData, ndarray[double] time):
    cdef ndarray[double] pressure = empty(numData)
    c_theis(&k, &nu, &phi, &rho, &c, &b, &Q0, &P0, &r, &numData, &time[0], &pressure[0])
    #c_theis(&k, &nu, &phi, &rho, &c, &b, &Q0, &P0, &r, &t0, &dt, &t1, &numData, <double*> pressure.data)
    return pressure