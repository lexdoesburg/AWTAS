
# import subprocess
# import time
# import model
# import data
# import numpy as np
# import matplotlib.pyplot as plt
# from theis_wrapper import theis
# from theis_solution import theis_solution

# program = '/Users/lexdoesburg/Documents/Uni2018/Summer_Research/Summer_Project/AWTAS/Fortran/Fortran_Awtas/theis_test.exe'
# # arguments = ['wellfit_theis.dat','theis_testdata.txt']
# fe_start = time.time()
# subprocess.call([program, 'wellfit_theis.dat', 'theis_testdata.txt'])
# fe_end = time.time()

# # Fortran module
# fm_start = time.time()

# P0 = 3.6e6 # Pa
# b = 100 # m
# r = 0.05 # m
# Q0 = -0.005 # m^3/s
# k = 1e-12 # m^2
# phi = 0.1
# rho = 813.37 # Water at 240 degrees celsius
# nu = 0.0001111 # Water at 240 degrees celsius
# c = 0.001303 # Water at 240 degrees celsius
# t0 = 0
# dt = 200
# t1 = 54000
# numData = 271

# time_array = np.linspace(0,54000,10000)
# # fm_pressure = theis(k, nu, phi, rho, c, b, Q0, P0, r, t0, dt, t1, numData)
# print('Entering for loop')
# for i in range(100,10000,250):    
#     time_array = np.linspace(0,54000,i)
#     start = time.time()
#     fm_pressure = theis(k, nu, phi, rho, c, b, Q0, P0, r, len(time_array), time_array)
#     end = time.time()
#     with open('fortran_theis_scaling.txt', 'a') as file:
#         file.write('{},{}\n'.format(i,end-start))
# print('Exiting for loop')

# fm_end = time.time()
# # print(fm_pressure)

# # Fortran module through python model
# ftm_start = time.time()
# theis_data = data.Data()
# theis_data.time = time_array

# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.005 # m^3/s
# k = 1e-12 # m^2
# phi = 0.1
# rho = 813.37 # Water at 240 degrees celsius
# nu = 0.0001111 # Water at 240 degrees celsius
# C = 0.001303 # Water at 240 degrees celsius

# parameters = [p0, qm, h, rho, nu, C, r]

# theis_data.parameters = parameters

# fortran_theis = model.Theis_Solution_Fortran(theis_data)
# ftm_pressure = fortran_theis.model([phi,k])
# ftm_end = time.time()

# # Python Module
# pm_start = time.time()

# theis_data = data.Data()
# theis_data.time = time_array

# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.005 # m^3/s
# k = 1e-12 # m^2
# phi = 0.1
# rho = 813.37 # Water at 240 degrees celsius
# nu = 0.0001111 # Water at 240 degrees celsius
# C = 0.001303 # Water at 240 degrees celsius

# parameters = [p0, qm, h, rho, nu, C, r]

# theis_data.parameters = parameters
# theis_model = model.Theis_Solution(theis_data)
# pm_pressure = theis_model.model([phi,k])
# pm_end = time.time()

# rp_start = time.time()
# rp_pressure = theis_solution(p0, qm, k, h, phi, rho, nu, C, r, time_array)
# rp_end = time.time()

# # print(pm_pressure-fm_pressure)
# print('Fortran Executable time elapsed: {}'.format(fe_end-fe_start))
# # print('Fortran Module time elapsed: {}'.format(fm_end-fm_start))
# print('Fortran Module through python model time elapsed: {}'.format(ftm_end-ftm_start))
# print('Python time elapsed: {}'.format(pm_end-pm_start))
# print('Raw python time elapsed: {}'.format(rp_end-rp_start))

# time,fe_pressure = np.genfromtxt('theis_testdata.txt', delimiter=' ', skip_header=6).T

# plt.plot(time,fe_pressure,'k-',label='Fortran Executable')
# # plt.plot(theis_data.time, fm_pressure,'g--',label='Fortran Module')
# plt.plot(theis_data.time, ftm_pressure,'g--',label='Fortran Module through python model')
# plt.plot(theis_data.time, pm_pressure,'r:',label='Python')
# plt.legend(loc='best')
# plt.show()


# # import subprocess
# # import time
# # import model
# # import data
# # import numpy as np
# # import matplotlib.pyplot as plt

# # program = '/Users/lexdoesburg/Documents/Uni2018/Summer_Research/Summer_Project/AWTAS/Fortran/Fortran_Awtas/xd.exe'
# # start = time.time()
# # subprocess.call([program])
# # end = time.time()

# # print('Fortran Time elapsed: {}'.format(end-start))

import numpy as np
import matplotlib.pyplot as plt
# time, pressure_1d = np.genfromtxt('1d_radial1.txt', delimiter=',', skip_header=0).T
# time2, pressure_1d2 = np.genfromtxt('1d_radial2.txt', delimiter=',', skip_header=0).T
# # print(time)
# # print(pressure_1d)
# plt.plot(time/3600,pressure_1d/1e5,'k-')
# plt.plot(time2/3600,pressure_1d2/1e5,'r--')
# plt.xlabel('Time [Hours]')
# plt.ylabel('Pressure [Bar]')
# plt.title('Homogeneous Porous Simulator Solution')
# plt.show()

time, pressure_1d = np.genfromtxt('radial_tester.txt', delimiter=',', skip_header=0).T
pressure_1d2 = np.genfromtxt('radial_tester2.txt', delimiter=',', skip_header=0).T
# print(time)
# print(pressure_1d)
print(pressure_1d-pressure_1d2)

# plt.plot(time/3600,pressure_1d/1e5,'k-')
plt.plot(time/3600,pressure_1d2/1e5,'k-')
plt.xlabel('Time [Hours]')
plt.ylabel('Pressure [Bar]')
plt.title('Homogeneous Porous Simulator Solution')
plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
# num, time_vals = np.genfromtxt('fortran_theis_scaling.txt', delimiter=',', skip_header=0).T
# plt.plot(num,time_vals)
# plt.show()
