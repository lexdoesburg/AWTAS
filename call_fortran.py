
import subprocess
import time
import model
import data
import numpy as np
import matplotlib.pyplot as plt

program = '/Users/lexdoesburg/Documents/Uni2018/Summer_Research/Summer_Project/AWTAS/Fortran/Fortran_Awtas/theis_test.exe'
# arguments = ['wellfit_theis.dat','theis_testdata.txt']
start = time.time()
subprocess.call([program, 'wellfit_theis.dat', 'theis_testdata.txt'])
end = time.time()

print('Fortran Time elapsed: {}'.format(end-start))

start = time.time()
theis_data = data.Data()
theis_data.time = np.linspace(0,54000,271)

p0 = 3.6e6 # Pa
h = 100 # m
r = 0.05 # m
qm = -0.005 # m^3/s
k = 1e-12 # m^2
phi = 0.1
rho = 813.37 # Water at 240 degrees celsius
nu = 0.0001111 # Water at 240 degrees celsius
C = 0.001303 # Water at 240 degrees celsius

parameters = [p0, qm, h, rho, nu, C, r]

theis_data.parameters = parameters
theis_model = model.Theis_Solution(theis_data)
pressure = theis_model.model([phi,k])
end = time.time()

print('Python Time elapsed: {}'.format(end-start))

time,pressure1 = np.genfromtxt('theis_testdata.txt', delimiter=' ', skip_header=6).T

plt.plot(time,pressure1,'k-',label='Fortran')
plt.plot(theis_data.time, pressure,'r--',label='Python')
plt.legend(loc='best')
plt.show()


# import subprocess
# import time
# import model
# import data
# import numpy as np
# import matplotlib.pyplot as plt

# program = '/Users/lexdoesburg/Documents/Uni2018/Summer_Research/Summer_Project/AWTAS/Fortran/Fortran_Awtas/xd.exe'
# start = time.time()
# subprocess.call([program])
# end = time.time()

# print('Fortran Time elapsed: {}'.format(end-start))