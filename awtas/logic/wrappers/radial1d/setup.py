from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

from numpy import get_include
from glob import glob
import os

"""
This script should be run through run_setup.py and run_setup.py should be run from the directory it is in for no errors.

Alternatively this script can be run directly from the directory it is in using the following line in the command line interface:
      python setup.py build_ext --inplace

Code based on answer by IanH (https://stackoverflow.com/questions/22404060/fortran-cython-workflow)
"""

# Fortran source files
source_files = ['variable_types.f90',
                'utility_functions.f90',
                'problem_data.f90',
                'thermodynamics.for',
                'matrixsolvers.for',
                'gammafunction.for',
                'modelprogress.f90',
                'numericalsimulator1d_routines.for',
                'numericalsimulator1d.f90',
                'homogeneousporoussimulator.f90',
                'call_radial1d.f90',
                'radial1d_wrapper.f90']

# Gfortran compilation flags
# NOTE: Remove '-mmacosx...' if trying on windows.
compile_flags = ['-Ofast', '-fPIC', '-mmacosx-version-min=10.12'] # -Ofast can be interchanged for -O3 for slightly reduced performance.
flags = ' '.join(compile_flags)
file_prefixes = [file.split('.')[0] for file in source_files]
object_files = [file_prefix + '.o' for file_prefix in file_prefixes]

# compile the fortran modules
fortran_source_dir = os.path.join(os.getcwd(), 'fortran_src')

for i, file in enumerate(source_files):
    print('gfortran {source} -c -o {obj} {extra_flags}'.format(source=os.path.join(fortran_source_dir, file), obj=object_files[i], extra_flags=flags))
    os.system('gfortran {source} -c -o {obj} {extra_flags}'.format(source=os.path.join(fortran_source_dir, file), obj=object_files[i], extra_flags=flags))

ext_modules = [Extension(# Module name
                         'radial1d_wrapper',
                         # Cython source
                         ['radial1d_wrapper.pyx'],
                         # Extra compile arguments for gcc
                         extra_compile_args=['-Ofast', '-fPIC'],
                         # library_dirs=['/Users/lexdoesburg/Documents/Uni2018/Summer_Research/Summer_Project/AWTAS/wrappers/radial1d/'],
                         libraries=['gfortran'],
                         # Other arguments and files to link
                         # NOTE: Remove '+ [-mmacosx...]' if trying on windows.
                         extra_link_args=object_files + ['-mmacosx-version-min=10.12'])]

setup(name = 'radial1d_wrapper',
      cmdclass = {'build_ext': build_ext},
      include_dirs = [get_include()], # include numpy headers
      ext_modules = ext_modules)
