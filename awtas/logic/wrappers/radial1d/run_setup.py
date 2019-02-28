import os
import shutil
from glob import glob

"""
This script should be run from the directory that it is in for no errors
"""

# Remove old wrapper file
if os.path.isfile('radial1d_wrapper.cpython-36m-darwin.so'):
    os.remove('radial1d_wrapper.cpython-36m-darwin.so')

# Build the new wrapper
print('python setup.py build_ext --inplace')
os.system('python setup.py build_ext --inplace')

# Remove the extra files produced when building the wrapper
current_dir = os.getcwd()
for file_extension in ('*.o', '*.mod', '*.c'):
    for file in glob(os.path.join(current_dir, file_extension)):
        os.remove(file)

if os.path.isdir('build'):
    shutil.rmtree('build')