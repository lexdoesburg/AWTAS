from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
system('gfortran variable_types.f90 -c -o variable_types.o -O3 -fPIC -mmacosx-version-min=10.9')
system('gfortran variable_parameters.f90 -c -o variable_parameters.o -O3 -fPIC -mmacosx-version-min=10.9')
system('gfortran utility_functions.f90 -c -o utility_functions.o -O3 -fPIC -mmacosx-version-min=10.9')
system('gfortran problem_data.f90 -c -o problem_data.o -O3 -fPIC -mmacosx-version-min=10.9')
system('gfortran theis_solution.f90 -c -o theis_solution.o -O3 -fPIC -mmacosx-version-min=10.9')
system('gfortran models.f90 -c -o models.o -O3 -fPIC -mmacosx-version-min=10.9')
system('gfortran theis_main.f90 -c -o theis_main.o -O3 -fPIC -mmacosx-version-min=10.9')
system('gfortran theis_wrapper.f90 -c -o theis_wrapper.o -O3 -fPIC -mmacosx-version-min=10.9')


ext_modules = [Extension(# module name:
                         'theis_wrapper',
                         # source file:
                         ['theis_wrapper.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=['variable_types.o', 'variable_parameters.o',
                                          'utility_functions.o', 'problem_data.o',
                                          'theis_solution.o', 'models.o',
                                          'theis_main.o', 'theis_wrapper.o'])]

setup(name = 'theis_wrapper',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)

# python3 setup.py build_ext --inplace