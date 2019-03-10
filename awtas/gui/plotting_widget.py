from PyQt5.QtWidgets import QWidget, QPushButton, QVBoxLayout, QFileDialog, QMessageBox, QGroupBox, QGridLayout, QLabel, QSpacerItem, QSizePolicy, QComboBox, QCheckBox, QLineEdit
from PyQt5.QtCore import pyqtSlot, Qt
from PyQt5.QtGui import QDoubleValidator

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.ticker import ScalarFormatter

import numpy as np

import awtas.logic.data as data_class
import awtas.logic.model as model

import time

class PlotWidget(QWidget):
    def __init__(self, parent=None):
        super(PlotWidget, self).__init__(parent)
        self.plotting_canvas = PlottingCanvas(self)
        self.plot_toolbar = NavigationToolbar(self.plotting_canvas, self)
        self.current_model = 'Analytical Theis'
        self.data = data_class.Data(model_type='theis') # Default model is Theis solution
        self.model = None
        self.data_imported = False
        self.parameters_imported = False
        self.initial_guess = False
        self.initial_variables = None
        self.prev_initial_variables = None
        self.build_UI()
            

    def build_UI(self):
        # Choose model type combobox
        self.model_type_combobox = QComboBox(self)
        self.model_type_combobox.setToolTip('Change the type of model used')
        self.model_type_combobox.addItems(['Analytical Theis','Homogeneous Porous'])
        self.model_type_combobox.activated[str].connect(self.change_model)
        self.model_type_label = QLabel('Model Type: ')
        self.model_type_layout = QGridLayout()
        self.model_type_layout.addWidget(self.model_type_label,0,0)
        self.model_type_layout.addWidget(self.model_type_combobox,0,1)

        # Import data button
        self.import_data_button = QPushButton('Import Data', self)
        self.import_data_button.setToolTip('Import observation and time data')
        self.import_data_button.clicked.connect(self.plot_data)

        # Fit model button
        self.fit_button = QPushButton('Fit Curve', self)
        self.fit_button.setToolTip('Fit the target variables to best match the observation data')
        self.fit_button.clicked.connect(self.fit_data)
        self.fit_button.setEnabled(False)

        # Group the two buttons together in a single layout
        self.button_layout = QVBoxLayout()
        self.button_layout.addLayout(self.model_type_layout)
        self.button_layout.addWidget(self.import_data_button)
        self.button_layout.addWidget(self.fit_button)

        # Initial guess check box
        self.initial_guess_check_box = QCheckBox('Provide initial guess for variables?')
        self.initial_guess_check_box.stateChanged.connect(self.show_initial_guess_widgets)

        # Initial guess group box
        self.initial_guess_groupbox, self.initial_guess_groupbox_layout = self.create_info_groupbox('Variable Initial Guess')
        self.initial_guess_value_labels = self.populate_info_groupbox_layout(self.initial_guess_groupbox_layout, 'Variables', values_editable=True)
        self.initial_guess_groupbox.setVisible(False)

        # # Update initial guess button
        # self.initial_guess_button = QPushButton('Update')
        # self.initial_guess_button.setToolTip('Updates the current initial guess to use the above input values')
        # self.initial_guess_button.clicked.connect(self.set_initial_guess)
        # self.initial_guess_button.setVisible(False)

        self.initial_guess_layout = QVBoxLayout()
        self.initial_guess_layout.addWidget(self.initial_guess_check_box)
        self.initial_guess_layout.addWidget(self.initial_guess_groupbox)
        # self.initial_guess_layout.addWidget(self.initial_guess_button)

        # Create parameter information groupboxes
        self.reservoir_conditions_groupbox, self.reservoir_conditions_groupbox_layout = self.create_info_groupbox('Reservoir Conditions')
        self.fixed_parameters_groupbox, self.fixed_parameters_groupbox_layout = self.create_info_groupbox('Fixed Parameters')
        self.variables_groupbox, self.variables_groupbox_layout = self.create_info_groupbox('Variable Estimations')

        # Populate parameter information groupboxes
        self.reservoir_condition_value_labels = self.populate_info_groupbox_layout(self.reservoir_conditions_groupbox_layout, 'Reservoir Conditions')
        self.fixed_parameter_value_labels = self.populate_info_groupbox_layout(self.fixed_parameters_groupbox_layout, 'Fixed Parameters')
        self.variable_value_labels = self.populate_info_groupbox_layout(self.variables_groupbox_layout, 'Variables')

        # Create the parameter infomation side bar
        self.parameter_sidebar = QVBoxLayout()
        self.parameter_sidebar.addLayout(self.initial_guess_layout)
        self.parameter_sidebar.addWidget(self.variables_groupbox)
        self.parameter_sidebar.addWidget(self.reservoir_conditions_groupbox)
        self.parameter_sidebar.addWidget(self.fixed_parameters_groupbox)
        self.parameter_sidebar.addSpacerItem(QSpacerItem(300, 0, hPolicy=QSizePolicy.Expanding, vPolicy=QSizePolicy.Expanding))

        # Create overall widget layout
        self.layout = QGridLayout()
        self.layout.addWidget(self.plotting_canvas, 1, 0)
        self.layout.addWidget(self.plot_toolbar, 0, 0)
        self.layout.addLayout(self.button_layout, 0, 1)
        self.layout.addLayout(self.parameter_sidebar, 1, 1)

        # Make the canvas at least 640x480 pixels
        self.layout.setColumnMinimumWidth(0, 640)
        self.layout.setRowMinimumHeight(1, 480)

        # Prioritise the canvas to stretch over other components
        self.layout.setColumnStretch(0, 1)

        # Set the widget layout
        self.setLayout(self.layout)

    def set_initial_guess(self):
        initial_porosity = float(self.initial_guess_value_labels['Porosity'].text())
        initial_permeability = float(self.initial_guess_value_labels['Permeability'].text())
        self.initial_variables = [initial_porosity, initial_permeability]
        print(self.initial_variables)


    def show_initial_guess_widgets(self, state):
        if state == Qt.Checked:
            self.initial_guess = True
            self.initial_guess_groupbox.setVisible(True)
            # self.initial_guess_button.setVisible(True)
        else:
            self.initial_guess = False
            self.initial_guess_groupbox.setVisible(False)
            # self.initial_guess_button.setVisible(False)


    def create_info_groupbox(self, parameter_type):
        """
        Creates a new groupbox for the purpose of displaying parameter information.
        """
        groupbox = QGroupBox(parameter_type)
        groupbox.setAlignment(4)
        groupbox_grid = QGridLayout()
        groupbox.setLayout(groupbox_grid)
        return groupbox, groupbox_grid


    def clear_layout(self, layout):
        """
        Recursively removes all widgets and sub-layouts from a layout.

        Code from https://stackoverflow.com/questions/22623151/python-how-to-unassign-layout-from-groupbox-in-pyqt answer by ekhumoro
        """
        if layout is not None:
            while layout.count():
                layout_item = layout.takeAt(0)
                widget = layout_item.widget()
                if widget is not None:
                    widget.deleteLater()
                else:
                    self.clear_layout(layout_item.layout())


    @pyqtSlot(str)
    def change_model(self, model_type):
        """
        Update the parameter sidebar and data structure when a different model type is selected using the model type combobox.
        """
        # If the model type is changing update GUI
        if model_type != self.current_model: 

            # Update current model
            self.current_model = model_type

            # Change data structure
            if model_type == 'Analytical Theis':
                self.data = data_class.Data(model_type='theis')
            elif model_type == 'Homogeneous Porous':
                self.data = data_class.Data(model_type='radial1d')

            # Clear layouts
            self.clear_layout(self.reservoir_conditions_groupbox_layout)
            self.clear_layout(self.fixed_parameters_groupbox_layout)
            self.clear_layout(self.variables_groupbox_layout)
            
            # Repopulate layouts
            self.reservoir_condition_value_labels = self.populate_info_groupbox_layout(self.reservoir_conditions_groupbox_layout, 'Reservoir Conditions')
            self.fixed_parameter_value_labels = self.populate_info_groupbox_layout(self.fixed_parameters_groupbox_layout, 'Fixed Parameters')
            self.variable_value_labels = self.populate_info_groupbox_layout(self.variables_groupbox_layout, 'Variables')
            

    def produce_output_file(self):
        filename, _ = QFileDialog.getSaveFileName(self, 'Save file as', '', 'Data Files (*.txt *.csv *.dat)')
        self.data.write_output_file(filename)


    @pyqtSlot()
    def plot_data(self):
        self.clear_all_parameters()
        try:
            import_successful = self.import_data_from_file()
        except ValueError:
            import_successful = False
            error_message = QMessageBox()
            error_message.setIcon(QMessageBox.Critical)
            error_message.setText('Import Error')
            error_message.setInformativeText('Data was not imported successfully. Please check that you are using the correct data file\
                and that the data file is formatted correctly.')
            error_message.setWindowTitle('Import Error')
            error_message.exec_()
        if import_successful:
            self.plotting_canvas.plot_observed_data(self.data)
            # Update the reservoir condition labels since the initial temperature/vapour saturation label will change depending on the newly imported data
            self.clear_layout(self.reservoir_conditions_groupbox_layout) # Clear layouts
            self.reservoir_condition_value_labels = self.populate_info_groupbox_layout(self.reservoir_conditions_groupbox_layout, 'Reservoir Conditions') # Repopulate layouts
            self.update_parameter_labels()
            self.data_imported = True
            self.parameters_imported = True
            if self.parameters_imported and self.data_imported:
                self.fit_button.setEnabled(True)
        

    def import_data_from_file(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'Load data file', '', 'Data Files (*.txt *.csv *.dat)')
        if filename:
            self.data.read_file(filename=filename)
            data_imported = True
        else:
            data_imported = False
        return data_imported


    @pyqtSlot()
    def fit_data(self):
        if self.data and self.data.fixed_parameters and self.data.reservoir_conditions:
            self.plotting_canvas.clear_fitted_lines()
            self.model = model.create_model(model_type=self.data.model_type, data=self.data)
            if self.initial_guess:
                self.set_initial_guess()
                self.model.find_model_parameters2(initial_guess=self.initial_variables, single_run=False)
            else:
                self.model.find_model_parameters2()

            self.plotting_canvas.plot_fit(self.data)
            for (parameter, item) in self.data.variables.items():
                if parameter == 'Porosity':
                    self.variable_value_labels[parameter].setText('{:.6f}'.format(item['Value']))
                elif parameter == 'Permeability':
                    self.variable_value_labels[parameter].setText('{:.6e}'.format(item['Value']))
                else:
                    self.variable_value_labels[parameter].setText('{}'.format(item['Value']))
        else:
            error_message = QMessageBox()
            error_message.setIcon(QMessageBox.Critical)
            error_message.setText('Error')
            error_message.setInformativeText('Please import model data first')
            error_message.setWindowTitle('Error')
            error_message.exec_()


    def populate_info_groupbox_layout(self, groupbox_layout, parameter_type, values_editable=False):
        parameters = {
            'Reservoir Conditions' : self.data.reservoir_conditions.items(),
            'Fixed Parameters' : self.data.fixed_parameters.items(),
            'Variables' : self.data.variables.items()
        }

        # info_value_labels = []
        info_value_labels = {}
        

        for i, (parameter, info) in enumerate(parameters[parameter_type]):
            units = info['Units']
            parameter_label = parameter
            # Special case as initial x can either be temperature or vapour saturation
            if parameter == 'Initial X':
                if self.data.initial_x:
                    parameter_label = self.data.initial_x
                    units = units[parameter_label]
                else:
                    parameter_label = 'Initial Temperature/Vapour Saturation'
                    units = 'Dimensionless'

            if units == 'Dimensionless':
                label = '{param}: '.format(param=parameter_label)
            else:
                label = '{param} [{unit}]: '.format(param=parameter_label, unit=units)

            groupbox_layout.addWidget(QLabel(label), i, 0)
            if values_editable:
                value_label = QLineEdit()
                lower_bound, upper_bound = data_class.parameter_bounds[parameter]
                if (lower_bound + upper_bound)/2 < 1e-4:
                    notation = 1
                else:
                    notation = 0
                value_validator = QDoubleValidator(lower_bound, upper_bound, 20)
                value_validator.setNotation(notation)
                value_label.setValidator(value_validator)
                value_label.setToolTip('Enter a value between {} and {}'.format(lower_bound, upper_bound))
            else:
                value_label = QLabel()
            # info_value_labels.append(value_label)
            info_value_labels[parameter] = value_label
            groupbox_layout.addWidget(value_label, i, 1)
        
        return info_value_labels
    

    def update_parameter_labels(self):
        # for i, info in enumerate(self.data.fixed_parameters.values()):
        #     self.fixed_parameter_value_labels[i].setText('{}'.format(info['Value']))
        # for i, info in enumerate(self.data.reservoir_conditions.values()):
        #     self.reservoir_condition_value_labels[i].setText('{:.2f}'.format(info['Value']))
        for (parameter, info) in self.data.fixed_parameters.items():
            self.fixed_parameter_value_labels[parameter].setText('{}'.format(info['Value']))
        for (parameter, info) in self.data.reservoir_conditions.items():
            self.reservoir_condition_value_labels[parameter].setText('{:.2f}'.format(info['Value']))
    

    def clear_all_parameters(self):
        if self.parameters_imported:
            for label in list(self.reservoir_condition_value_labels.values()) + list(self.fixed_parameter_value_labels.values()) + list(self.variable_value_labels.values()):
                label.setText('')

class PlottingCanvas(FigureCanvas):
    def __init__(self, parent=None):
        # self.observation_property_info = {
        #     0 : {'Label' : 'Deliverability [Units TODO:]', 'Scale Factor' : 1 },
        #     1 : {'Label' : 'Pressure [Bar]', 'Scale Factor' : 1e-5 },
        #     2 : {'Label' : 'Temperature [°C]', 'Scale Factor' : 1 },
        #     3 : {'Label' : 'Enthalpy [kJ/kg]', 'Scale Factor' : 1e-3 }
        # }
        self.observation_scale = 1e-5 # Only currently useful for converting pressure from Pa to Bar
        self.figure = Figure()
        self.axes = None
        FigureCanvas.__init__(self, self.figure)
        self.setParent(parent)
        self.fitted_lines = []
        self.draw()
    

    def plot_observed_data(self, data):
        # Time bottleneck in here is in setting the layout of the plot (either using tight layout or draw).
        if self.figure.get_axes():
            self.figure.clear()
            self.fitted_lines = []
        self.axes = self.figure.add_subplot(1,1,1)
        self.axes.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))

        time_scale, x_label = self.get_time_axis_scale(data, log_scale=False)
        # TODO: Update theis solution code so that it uses observation points
        # TODO: Allow plotting of only single observation points in the future
        if data.model_type == 'theis':
            # self.axes.semilogx(np.log(data.time*time_scale), data.observation*self.observation_scale, 'kx', label='Observed Data')
            self.axes.plot(data.time*time_scale, data.observation*self.observation_scale, 'kx', label='Observed Data')
        elif data.model_type == 'radial1d':
            if data.observation_points.num_observation_points == 1:
                # observation_scale = self.observation_property_info[data.observation_points.property[0]]['Scale Factor']
                # self.axes.semilogx(np.log(data.observation_points.times[0]*time_scale), data.observation_points.observations[0]*self.observation_scale, 'kx', label='Observed Data')
                self.axes.plot(data.observation_points.times[0]*time_scale, data.observation_points.observations[0]*self.observation_scale, 'kx', label='Observed Data')
            else:
                for i in range(data.observation_points.num_observation_points):
                    # observation_scale = self.observation_property_info[data.observation_points.property[i]]['Scale Factor']
                    # self.axes.semilogx(np.log(data.observation_points.times[i]*time_scale), data.observation_points.observations[i]*self.observation_scale, 'x', label='Observed Data {}'.format(i+1))
                    self.axes.plot(data.observation_points.times[i]*time_scale, data.observation_points.observations[i]*self.observation_scale, 'x', label='Observed Data {}'.format(i+1))

        title = '{} Well Observed Values'.format(data.get_well_type())
        self.axes.set_title(title)
        self.axes.set_xlabel(x_label) # TODO: Allow for either log or linear time plots to plotted
        self.axes.set_ylabel('Pressure [Bar]') # TODO: Change this label adaptively depending on the type of data plotted
        self.axes.legend(loc='best')
        self.figure.tight_layout(pad = 2)
        self.draw()
    

    def get_time_axis_scale(self, data, log_scale=False):
        if log_scale:
            if data.time[-1] >= 50*24*3600: # Time series goes to at least 50 days
                time_scale = 1/3600
                axis_label = 'Log Time [hours]'
            else:
                time_scale = 1
                axis_label = 'Log Time [seconds]'
        else:
            if data.time[-1] < 18000:
                time_scale = 1
                axis_label = 'Time [seconds]'
            elif 18000 <= data.time[-1] <= 432000:
                time_scale = 1/3600
                axis_label = 'Time [hours]'
            else:
                time_scale = 1/86400
                axis_label = 'Time [days]'
        return time_scale, axis_label


    def clear_fitted_lines(self):
        if self.fitted_lines:
            num_lines = len(self.fitted_lines)
            for i in range(num_lines):
                self.fitted_lines[i][0].remove()
            self.fitted_lines = []
            self.draw()


    def plot_fit(self, data):
        time_scale, _ = self.get_time_axis_scale(data, log_scale=False)
        if data.model_type == 'theis':
            # self.fitted_lines.append(self.axes.semilogx(np.log(data.time*time_scale), data.approximation*self.observation_scale, 'r-', label='Fitted Approximation'))
            self.fitted_lines.append(self.axes.plot(data.time*time_scale, data.approximation*self.observation_scale, 'r-', label='Fitted Approximation'))
        elif data.model_type == 'radial1d':
            if data.observation_points.num_observation_points == 1:
                # observation_scale = self.observation_property_info[data.observation_points.property[0]]['Scale Factor']
                # self.fitted_lines.append(self.axes.semilogx(np.log(data.observation_points.times[0]*time_scale), data.observation_points.modelled_values[0]*self.observation_scale, 'r-', label='Fitted Approximation'))
                self.fitted_lines.append(self.axes.plot(data.observation_points.times[0]*time_scale, data.observation_points.modelled_values[0]*self.observation_scale, 'r-', label='Fitted Approximation'))
            else:
                for i in range(data.observation_points.num_observation_points):
                    # observation_scale = self.observation_property_info[data.observation_points.property[i]]['Scale Factor']
                    # self.fitted_lines.append(self.axes.semilogx(np.log(data.observation_points.times[i]*time_scale), data.observation_points.modelled_values[i]*self.observation_scale, '-', label='Fitted Approximation {}'.format(i+1)))
                    self.fitted_lines.append(self.axes.plot(data.observation_points.times[i]*time_scale, data.observation_points.modelled_values[i]*self.observation_scale, '-', label='Fitted Approximation {}'.format(i+1)))
        
        title = '{} Well Observed vs Modelled Values'.format(data.get_well_type())
        self.axes.set_title(title)
        self.axes.legend(loc='best')
        self.draw()
