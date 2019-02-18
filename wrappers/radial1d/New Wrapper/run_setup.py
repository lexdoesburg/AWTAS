import os
import shutil
from glob import glob

if os.path.isdir('build'):
    shutil.rmtree('build')
if os.path.isfile('radial1d_wrapper.cpython-36m-darwin.so'):
    os.remove('radial1d_wrapper.cpython-36m-darwin.so')

print('python3 setup.py build_ext --inplace')
os.system('python3 setup.py build_ext --inplace')


# Remove the compiled object & mod files
# target_dir = os.path.join(os.getcwd(), 'fortran_src')
target_dir = os.getcwd()

for file_extension in ('*.o', '*.mod', '*.c'):
    for file in glob(os.path.join(target_dir, file_extension)):
        os.remove(file)