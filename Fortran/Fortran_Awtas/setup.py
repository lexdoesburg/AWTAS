from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
fortran_mod_comp = 'gfortran -c theis_main.f90'
print(fortran_mod_comp)
system(fortran_mod_comp)
shared_obj_comp = 'gfortran -c theis_wrapper.f90'
print(shared_obj_comp)
system(shared_obj_comp)

ext_modules = [Extension(# module name:
                         'TheisSolution',
                         # source file:
                         ['theis_wrapper.pyx'],
                         # other compile args for gcc
                        #  extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=['theis_main.o', 'theis_wrapper.o', 'variable_types.o', 'problem_data.o', 'variable_parameters.o', 'models.o', 'noise.o', 'theis_solution.o', 'utility_functions.o', 'homogeneousporous.o', 'thermodynamics.o', 'numericalsimulator1d_routines.o', 'NumericalSimulator1D.o', 'gammafunction.o', 'matrixsolvers.o', 'modelprogress.o'])]

setup(name = 'TheisSolution',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)

# python setup.py build_ext --inplace