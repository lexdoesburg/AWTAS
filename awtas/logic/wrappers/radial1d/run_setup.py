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
current_dir = os.getcwd()

for file_extension in ('*.o', '*.mod', '*.c'):
    for file in glob(os.path.join(current_dir, file_extension)):
        os.remove(file)

if os.path.isdir('build'):
    shutil.rmtree('build')