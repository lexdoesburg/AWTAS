
# import subprocess
# import time

# program = '/Users/lexdoesburg/Documents/Uni2018/Summer_Research/Summer_Project/AWTAS/Fortran/Fortran_Awtas/theis_test'
# # arguments = ['wellfit_theis.dat','theis_testdata.txt']
# start = time.time()
# subprocess.call([program, 'wellfit_theis.dat', 'theis_testdata.txt'])
# end = time.time()

# print('Time elapsed: {}'.format(end-start))

# print('Called program')

import time

from radial1d_wrapper import radial1d
import numpy as np
import matplotlib.pyplot as plt

k=1.0036337151556582e-13
phi=0.09293594872248502
Pressure0=40e5
X0=200
rw=0.1
thick=100
CR=1000
COND=2.5
RHOR=2500
COMP=0
ConstRate=5
distFromWell=0.1 # distance of observation point from action well
numData=271
time_array=np.linspace(0,54000,numData)

# num_readings = []
# time_readings = []
# for i in range(100,10000,100):
#     print('Entering for loop for {} data points'.format(i))
#     time_array=np.linspace(0,54000,i)
#     numData = len(time_array)
#     start = time.time()
#     pressure = radial1d(phi, k, Pressure0, X0, rw, thick, CR, COND, RHOR, COMP, ConstRate, distFromWell, numData, time_array)
#     end = time.time()
#     print('Radial1D ',end-start)
#     time_readings.append(end-start)
#     num_readings.append(i)
#     with open('radial1d_scale1.txt', 'a') as file:
#         file.write('{},{}\n'.format(i,end-start))
# num_readings,time_readings = np.genfromtxt('radial1d_scale1.txt', delimiter=',').T
# plt.plot(num_readings, time_readings)
# plt.xlabel('Number of data')
# plt.ylabel('Time to evaluate')
# plt.show()

start = time.time()
pressure = radial1d(phi, k, Pressure0, X0, rw, thick, CR, COND, RHOR, COMP, ConstRate, distFromWell, numData, time_array)
end = time.time()
print(pressure)
print(len(pressure))
print('Time elapsed: {}'.format(end-start))

plt.plot(time_array/3600,pressure/1e5,'k-')
plt.xlabel('Time [Hours]')
plt.ylabel('Pressure [Bar]')
plt.title('Homogeneous Porous Simulator Solution')
plt.show()

# for i in range(1000):
#     pressure = radial1d(phi, k, Pressure0, X0, rw, thick, CR, COND, RHOR, COMP, ConstRate, distFromWell, numData, time_array)
#     print('Call {}'.format(i))
# print('Completed 1000 calls')


# # Testing the radial1d through an executable
# import subprocess
# import time
# import numpy as np
# import matplotlib.pyplot as plt

# program = '/Users/lexdoesburg/Documents/Uni2018/Summer_Research/Summer_Project/AWTAS/Fortran/Fortran_Awtas/radial.exe'
# # arguments = ['wellfit_theis.dat','theis_testdata.txt']
# infile = 'radial1d_exe.txt'
# outfile = 'radial1d_exe_output.txt'
# print('Calling fortran executable')
# start = time.time()
# subprocess.call([program, infile, outfile])
# pressure = np.genfromtxt(outfile).T
# end = time.time()
# print('Call time: {}'.format(end-start))
# print(pressure)


# outfile = 'radial1d_exe_output.txt'
# pressure = np.genfromtxt(outfile).T
# print(pressure)

# plt.plot(np.linspace(0,54000,271),pressure)
# plt.show()

