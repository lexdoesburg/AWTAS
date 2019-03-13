Python3.6 was used for development.

# External packages used (all included with Anaconda3)
* NumPy (model, data, testing)
* Matplotlib (testing, GUI)
* SciPy (model)
* PyQt5 (GUI)

It is therefore recommended to install Python through Anaconda3.

To run the code it is also required to install the gfortran Fortran compiler. This is because the code also contains a Fortran wrapper (found in `awtas/logic/wrappers/radial1d/*.so`). This wrapper was built on MacOS 10.12 so will only work on Mac systems with MacOS 10.12 and up. 

To run the code on Windows this Fortran wrapper will need to be rebuilt. To do this it will require a Windows machine of the target operating system with gfortran installed. The MacOS wrapper (the `.so` file in `awtas/logic/wrappers/radial1d/`) can be deleted. The wrapper can then be rebuilt by slightly changing and then running `awtas/logic/wrappers/radial1d/setup.py`. There are some notes in the file itself about what will need changing. I can’t guarantee that it will work first try and there may need to be some other changes made to get it working. Once the wrapper has been rebuilt, the remainder of the Python code should work as expected in Windows.

# Running the program on MacOS 10.12 or greater
* Download the code and then using the terminal change directory to the folder that `run_awtas.py` is in and then run this script using the command `python run_awtas.py`.

# Using the program
The program is currently very simple and the way to use the program is as follows:
* Change to the desired model type in the top right corner. Note the analytical Theis model was mostly used for early experimentation and does not have as much functionality as the Homogeneous Porous model (cannot use multiple flows or observation types other than pressure).
* Import a data file. There are example data files in `awtas/gui/example_datafiles/` that can be used to try out the program.
* Choose whether or not to supply an initial guess of the unknown variable values (if no guess is supplied then one will be estimated).
* Click the fit curve button.
* Once the model has been fit the estimated unknown variable values will be displayed.
* The model output and estimated variable values can be exported to a data file by going to the `File menu > Export Results`.

# Making your own data files manually
* These data files have a strict format and if a new model is to be tested then ensure that none of the wording is changed. With that said the "KNOWN VARIABLE VALUES" section is optional and can be deleted if wanted. Ensure that all data is separated with a single comma (`,`) and no spaces.
* The different type of pumping schemes available for the Homogeneous Porous model are:
    * "Measured Flows" - Flow rates and times must be specified.
    * "Step Flows" - The flow rate and the duration of the flow must be specified.
    * "Constant Flow" - The flow rate must be specified, the flow time should be left as "0.0".
    * There is an example of each of these in the example_datafiles
* The different types of observation properties available for the Homogeneous Porous model are:
    * "Pressure" - Measurements in Pascals
    * "Temperature" - Measurements in degrees celsius
    * "Enthalpy" - Measurements in J/kg
    * "Deliverability" - Note deliverability hasn't been tested
* Note: the observation point radial location should be slightly smaller than the distance you want to estimate at. For example if you want to estimate at 0.1m then the observation location should be 0.0999m or similar. This is a workaround and needs to be addressed in the code.
* Also note: if the first block of the grid is the location where the estimation is to occur then the observation location should simply be smaller than the action well radius.
* Finally, the units of each parameter can be found in `data.py` in the `_default_parameter_units` dictionary.