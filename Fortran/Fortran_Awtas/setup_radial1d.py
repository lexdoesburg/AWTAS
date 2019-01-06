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
# system('gfortran homogeneousporous.f90 -c -o homogeneousporous.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran thermodynamics.for -c -o thermodynamics.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran numericalsimulator1d_routines.for -c -o numericalsimulator1d_routines.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran NumericalSimulator1D.f90 -c -o NumericalSimulator1D.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran gammafunction.for -c -o gammafunction.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran matrixsolvers.for -c -o matrixsolvers.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran modelprogress.f90 -c -o modelprogress.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran radial1d_main.f90 -c -o radial1d_main.o -O3 -fPIC -mmacosx-version-min=10.9')
# system('gfortran radial1d_wrapper.f90 -c -o radial1d_wrapper.o -O3 -fPIC -mmacosx-version-min=10.9')


ext_modules = [Extension(# module name:
                         'radial1d_wrapper',
                         # source file:
                         ['radial1d_wrapper.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                        #  library_dirs=[r'C:\Program Files\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\lib\gcc\x86_64-w64-mingw32\8.1.0\\'],
                         libraries=['gfortran'],
                         # other files to link to
                         extra_link_args=['variable_types.o', 'variable_parameters.o',
                                          'utility_functions.o', 'problem_data.o',
                                          'theis_solution.o', 'models.o',
                                          'theis_main.o', 'homogeneousporous.o',
                                          'thermodynamics.o', 'numericalsimulator1d_routines.o',
                                          'NumericalSimulator1D.o', 'gammafunction.o',
                                          'matrixsolvers.o', 'modelprogress.o',
                                          'radial1d_main.o', 'radial1d_wrapper.o'])]

setup(name = 'radial1d_wrapper',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)

# python setup_radial1d.py build_ext --inplace --compiler=mingw32
# python3 setup_radial1d.py build_ext --inplace

# import sys
# print(sys.version)

# from distutils.core import setup
# from Cython.Build import cythonize

# setup(name='Hello world app',
#       ext_modules=cythonize("hello.pyx"))
