from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system, sep

# # compile the fortran modules without linking (mac)
# system('gfortran variable_types.f90 -c -o variable_types.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran variable_parameters.f90 -c -o variable_parameters.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran utility_functions.f90 -c -o utility_functions.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran problem_data.f90 -c -o problem_data.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran theis_solution.f90 -c -o theis_solution.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran models.f90 -c -o models.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran theis_main.f90 -c -o theis_main.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran theis_wrapper.f90 -c -o theis_wrapper.o -O3 -fPIC -mmacosx-version-min=10.9')

# compile the fortran modules without linking (windows)
system('gfortran variable_types.f90 -c -o variable_types.o -O3 -fPIC')
system('gfortran variable_parameters.f90 -c -o variable_parameters.o -O3 -fPIC')
system('gfortran utility_functions.f90 -c -o utility_functions.o -O3 -fPIC')
system('gfortran problem_data.f90 -c -o problem_data.o -O3 -fPIC')
system('gfortran theis_solution.f90 -c -o theis_solution.o -O3 -fPIC')
system('gfortran models.f90 -c -o models.o -O3 -fPIC')
system('gfortran theis_main.f90 -c -o theis_main.o -O3 -fPIC')
system('gfortran theis_wrapper.f90 -c -o theis_wrapper.o -O3 -fPIC')



ext_modules = [Extension(# module name:
                         'theis_wrapper',
                         # source file:
                         ['theis_wrapper.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                        #  library_dirs=[r'C:\Program Files\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\lib\gcc\x86_64-w64-mingw32\8.1.0\\'],
                         libraries=['gfortran'],
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

# python setup.py build_ext --inplace --compiler=mingw64
# python3 setup.py build_ext --inplace

# import sys
# print(sys.version)