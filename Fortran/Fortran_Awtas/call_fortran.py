
import subprocess
import time

program = '/Users/lexdoesburg/Documents/Uni2018/Summer_Research/Summer_Project/AWTAS/Fortran/Fortran_Awtas/theis_test'
# arguments = ['wellfit_theis.dat','theis_testdata.txt']
start = time.time()
subprocess.call([program, 'wellfit_theis.dat', 'theis_testdata.txt'])
end = time.time()

print('Time elapsed: {}'.format(end-start))

print('Called program')