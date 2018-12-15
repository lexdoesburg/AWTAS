from numpy cimport ndarray
from numpy import empty

cdef extern from "theis_wrapper.h":
    void c_Theis(double* k, double* nu, double* phi, double* rho, double* c, double* b, double* Q0, double* P0, double* r, double* t0, double* dt, double* t1, int* numData, double* pressure)

def theis(double k, double nu, double phi, double rho, double c, double b, double Q0, double P0, double r, double t0, double dt, double t1, int numData):
    cdef ndarray[double] pressure = empty(numData)
    #c_Theis(&k, &nu, &phi, &rho, &c, &b, &Q0, &P0, &r, &t0, &dt, &t1, &numData, &pressure[0])
    c_Theis(&k, &nu, &phi, &rho, &c, &b, &Q0, &P0, &r, &t0, &dt, &t1, &numData, <double*> pressure.data)
    return pressure