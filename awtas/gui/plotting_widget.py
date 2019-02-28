from PyQt5.QtWidgets import QWidget, QPushButton, QVBoxLayout, QFileDialog, QMessageBox, QGroupBox, QGridLayout, QLabel, QSpacerItem, QSizePolicy, QComboBox
from PyQt5.QtCore import pyqtSlot

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

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
        self.init_UI()
            

    def init_UI(self):
        # Choose model type combobox
        self.model_type_combobox = QComboBox(self)
        self.model_type_combobox.addItems(['Analytical Theis','Homogeneous Porous'])
        self.model_type_combobox.activated[str].connect(self.change_model)
        self.model_type_label = QLabel('Model Type: ')
        self.model_type_layout = QGridLayout()
        self.model_type_layout.addWidget(self.model_type_label,0,0)
        self.model_type_layout.addWidget(self.model_type_combobox,0,1)

        # Import data button
        self.import_data_button = QPushButton('Import Data', self)
        self.import_data_button.setToolTip('Import pressure measurements and time data.')
        self.import_data_button.clicked.connect(self.plot_data)

        # Fit model button
        self.fit_button = QPushButton('Fit Curve', self)
        self.fit_button.setToolTip('Fit the permeability and porosity to match the data')
        self.fit_button.clicked.connect(self.fit_data)
        self.fit_button.setEnabled(False)

        # Group the two buttons together in a single layout
        self.button_layout = QVBoxLayout()
        self.button_layout.addLayout(self.model_type_layout)
        self.button_layout.addWidget(self.import_data_button)
        self.button_layout.addWidget(self.fit_button)

        # Create parameter information groupboxes
        self.reservoir_conditions_groupbox, self.reservoir_conditions_groupbox_layout = self.create_info_groupbox('Reservoir Conditions')
        self.fixed_parameters_groupbox, self.fixed_parameters_groupbox_layout = self.create_info_groupbox('Fixed Parameters')
        self.variables_groupbox, self.variables_groupbox_layout = self.create_info_groupbox('Variables')

        # Populate parameter information groupboxes
        self.reservoir_condition_value_labels = self.populate_info_groupbox_layout(self.reservoir_conditions_groupbox_layout, 'Reservoir Conditions')
        self.fixed_parameter_value_labels = self.populate_info_groupbox_layout(self.fixed_parameters_groupbox_layout, 'Fixed Parameters')
        self.variable_value_labels = self.populate_info_groupbox_layout(self.variables_groupbox_layout, 'Variables')

        # Create the parameter infomation side bar
        self.parameter_sidebar = QVBoxLayout()        
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

        # Make the canvas at least 500x500 pixels
        # Canvas size is good at 500x500 windows. 700x500 mac.
        self.layout.setColumnMinimumWidth(0, 700)
        self.layout.setRowMinimumHeight(1, 500)

        # Prioritise the canvas to stretch over other components
        self.layout.setColumnStretch(0, 1)

        # Set the widget layout
        self.setLayout(self.layout)


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
        print(model_type)
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
            

    @pyqtSlot()
    def plot_data(self):
        self.clear_all_parameters()
        # filename, _ = QFileDialog.getOpenFileName(self, 'Load data file', "", '*.dat;*.txt')
        # if filename:
        #     start = time.time()
        #     if self.data:
        #         self.data.read_file(filename=filename)
        #     else:
        #         self.data = data_class.Data(filename=filename)
        #     end = time.time()
        #     print('Time elapsed = {}'.format(end - start))
        self.import_data_from_file()
            # self.clear_all_parameters()

        start = time.time()
        self.plotting_canvas.plot_observed_data(self.data)
        end = time.time()
        print('Plotting Observed Data Time Elapsed = {}'.format(end - start))
        
        self.update_parameter_labels()
        self.data_imported = True
        self.parameters_imported = True
        if self.parameters_imported and self.data_imported:
            self.fit_button.setEnabled(True)


    def import_data_from_file(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'Load data file', "", '*.dat;*.txt')
        if filename:
            start = time.time()
            self.data.read_file(filename=filename)
            end = time.time()
            print('Time elapsed = {}'.format(end - start))


    @pyqtSlot()
    def fit_data(self):
        if self.data and self.data.fixed_parameters and self.data.reservoir_conditions:
            self.model = model.Theis_Solution(self.data)
            self.model.find_model_parameters()
            self.plotting_canvas.plot_fit(self.data)
            for i, label in enumerate(self.variable_value_labels):
                if i != len(self.data.variables)-1:
                    label.setText('{:.6f}'.format(self.data.variables[i]))
                else:
                    label.setText('{:.6e}'.format(self.data.variables[i]))
        else:
            error_message = QMessageBox()
            error_message.setIcon(QMessageBox.Critical)
            error_message.setText('Error')
            error_message.setInformativeText('Please enter Theis solution model parameters first')
            error_message.setWindowTitle('Error')
            error_message.exec_()


    def populate_info_groupbox_layout(self, groupbox_layout, parameter_type):
        parameters = {
            'Reservoir Conditions' : self.data.reservoir_conditions.items(),
            'Fixed Parameters' : self.data.fixed_parameters.items(),
            'Variables' : self.data.variables.items()
        }

        info_value_labels = []

        for i, (parameter, info) in enumerate(parameters[parameter_type]):
            units = info['Units']
            
            # Special case as initial x can either be temperature or vapour saturation
            if parameter == 'Initial X':
                if self.data.initial_x:
                    parameter = self.data.initial_x
                    units = units[parameter]
                else:
                    parameter = 'Initial Temperature/Vapour Saturation'
                    units = 'Dimensionless'

            if units == 'Dimensionless':
                label = '{param}: '.format(param=parameter)
            else:
                label = '{param} [{unit}]: '.format(param=parameter, unit=units)

            groupbox_layout.addWidget(QLabel(label), i, 0)
            value_label = QLabel()
            info_value_labels.append(value_label)
            groupbox_layout.addWidget(value_label, i, 1)
        
        return info_value_labels
    

    # def create_parameters_groupbox(self):
    #     parameters_groupbox = QGroupBox("Known Parameters")
    #     parameters_groupbox.setAlignment(4)
    #     # parameters_groupbox.setCheckable(True)
    #     # parameters_groupbox.setChecked(False)

    #     parameters_grid = QGridLayout()
    #     for i, (parameter, info) in enumerate(self.data.parameters.items()):
    #         parameters_grid.addWidget(QLabel('{} [{}]: '.format(parameter, info['Units'])), i, 0)
    #         new_label = QLabel()
    #         self.parameter_value_labels.append(new_label)
    #         parameters_grid.addWidget(new_label, i, 1)
    #     parameters_groupbox.setLayout(parameters_grid)

    #     variables_groupbox = QGroupBox('Unknown Parameters')
    #     variables_groupbox.setAlignment(4)

    #     variables_grid = QGridLayout()
    #     variables_grid.addWidget(QLabel('Porosity: '), 0, 0)
    #     variables_grid.addWidget(QLabel('Permeability [m<sup>2</sup>]: '), 1, 0)
    #     for i in range(variables_grid.rowCount()):
    #         new_label = QLabel()
    #         self.variable_value_labels.append(new_label)
    #         variables_grid.addWidget(new_label, i, 1)
    #     variables_groupbox.setLayout(variables_grid)

    #     fullWidget = QVBoxLayout()
    #     fullWidget.addWidget(parameters_groupbox)
    #     fullWidget.addWidget(variables_groupbox)
    #     # Spaceritem is good at 220, 0 for windows. 300, 0 for mac.
    #     fullWidget.addSpacerItem(QSpacerItem(300, 0, hPolicy=QSizePolicy.Expanding, vPolicy=QSizePolicy.Expanding))
    #     return fullWidget


    def update_parameter_labels(self):
        for i, info in enumerate(self.data.fixed_parameters.values()):
            self.fixed_parameter_value_labels[i].setText(str(info['Value']))
        for i, info in enumerate(self.data.reservoir_conditions.values()):
            self.reservoir_condition_value_labels[i].setText(str(info['Value']))
    

    def clear_all_parameters(self):
        if self.parameters_imported:
            for label in self.reservoir_condition_value_labels + self.fixed_parameter_value_labels + self.variable_value_labels:
                label.setText('')


class PlottingCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.figure = Figure()
        self.axes = None
        FigureCanvas.__init__(self, self.figure)
        self.setParent(parent)
        # self.plot_toolbar = NavigationToolbar(self.plotting_canvas, self)
        self.fitted_lines = []
        self.draw()
    

    def plot_observed_data(self, data):
        # Time bottleneck here is in setting the layout of the plot (either using tight layout or draw).
        start = time.clock()
        if self.figure.get_axes():
            self.figure.clear()
            self.fitted_lines = []
        end = time.clock()
        print('time figure clear: {}'.format(end-start))

        start = time.clock()
        self.axes = self.figure.add_subplot(1,1,1)
        end = time.clock()
        print('time add axes: {}'.format(end-start))
        start = time.clock()
        self.axes.semilogx(np.log(data.time), data.observation, 'kx', label='Observed Data')
        end = time.clock()
        print('time plot data (inc log data): {}'.format(end-start))
        start = time.clock()
        self.axes.set_xlabel('Log Time (log (s))')
        self.axes.set_ylabel('Pressure (Pa)')
        self.axes.legend(loc='lower left')
        end = time.clock()
        print('time set labels: {}'.format(end-start))
        start = time.clock()
        self.figure.tight_layout(pad = 2)
        end = time.clock()
        print('time tight layout: {}'.format(end-start))
        start = time.clock()
        self.draw()
        end = time.clock()
        print('time draw: {}'.format(end-start))
    

    def plot_fit(self, data):
        if self.fitted_lines:
            self.axes.lines.remove(self.fitted_lines[0])
        self.fitted_lines = self.axes.semilogx(np.log(data.time), data.approximation, 'r-', label='Fitted Approximation')
        self.axes.legend(loc='lower left')
        self.draw()
