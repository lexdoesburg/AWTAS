from numpy cimport ndarray
from numpy import empty

cdef extern from "radial1d_wrapper.h":
    void c_radial1d(double* phi, double* k, double* Pressure0, double* X0, double* rw, double* thick, double* CR, double* COND, double* RHOR, double* COMP, double* ConstRate, double* distFromWell, int* numData, double* time, double* pressure)

def radial1d(double phi, double k, double Pressure0, double X0, double rw, double thick, double CR, double COND, double RHOR, double COMP, double ConstRate, double distFromWell, int numData, ndarray[double] time):
    cdef ndarray[double] pressure = empty(numData)
    c_radial1d(&phi, &k, &Pressure0, &X0, &rw, &thick, &CR, &COND, &RHOR, &COMP, &ConstRate, &distFromWell, &numData, &time[0], &pressure[0])
    #c_theis(&k, &nu, &phi, &rho, &c, &b, &Q0, &P0, &r, &t0, &dt, &t1, &numData, <double*> pressure.data)
    return pressure