import numpy as np
import matplotlib.pyplot as plt
from radial1d_wrapper import radial1d

time, pressure_tough2 = np.genfromtxt('Pwell.dat', delimiter='   ').T
print(time)

time = [float(t) for t in time]
print(time)
time = np.asarray(time)

k=0.1e-13
phi=0.100000001
Pressure0=40e5
X0=200
rw=0.1
thick=100
CR=1000
COND=2.5
RHOR=2500
COMP=0
ConstRate=-5
distFromWell=0 # distance of observation point from action well
numData=len(time)

pressure = radial1d(phi, k, Pressure0, X0, rw, thick, CR, COND, RHOR, COMP, ConstRate, distFromWell, numData, time)

plt.plot(time, pressure_tough2, 'k-')
plt.plot(time, pressure/1e5, 'r--')
plt.show()
