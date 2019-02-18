from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

from numpy import get_include
from glob import glob
import os

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
                # 'radial1d_main.f90',
                'call_radial1d.f90',
                'radial1d_wrapper.f90']

compile_flags = ['-Ofast', '-fPIC', '-mmacosx-version-min=10.12']
flags = ' '.join(compile_flags)
file_prefixes = [file.split('.')[0] for file in source_files]
object_files = [file_prefix + '.o' for file_prefix in file_prefixes]

# compile the fortran modules (mac)
for i, file in enumerate(source_files):
    print('gfortran {source} -c -o {obj} {extra_flags}'.format(source=file, obj=object_files[i], extra_flags=flags))
    os.system('gfortran {source} -c -o {obj} {extra_flags}'.format(source=file, obj=object_files[i], extra_flags=flags))

ext_modules = [Extension(# module name:
                         'radial1d_wrapper',
                         # source file:
                         ['radial1d_wrapper.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-Ofast', '-fPIC'],
                        #  library_dirs=['/Users/lexdoesburg/Documents/Uni2018/Summer_Research/Summer_Project/AWTAS/wrappers/radial1d/'],
                         libraries=['gfortran'],
                         # other files to link to
                         extra_link_args=object_files + ['-mmacosx-version-min=10.12'])]

setup(name = 'radial1d_wrapper',
      cmdclass = {'build_ext': build_ext},
      include_dirs = [get_include()], # include numpy headers
      ext_modules = ext_modules)


#-----------------------------------------------------------
# To build extension run the following line in command line
#-----------------------------------------------------------
# python setup.py build_ext --inplace

