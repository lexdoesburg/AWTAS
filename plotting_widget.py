from PyQt5.QtWidgets import QWidget, QPushButton, QVBoxLayout, QFileDialog, QMessageBox, QGroupBox, QGridLayout, QLabel, QSpacerItem, QSizePolicy
from PyQt5.QtCore import pyqtSlot

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import numpy as np

import data as data_class
import model

import time

class PlotWidget(QWidget):
    def __init__(self, parent=None):
        super(PlotWidget, self).__init__(parent)
        self.plotting_canvas = PlottingCanvas(self)
        self.plot_toolbar = NavigationToolbar(self.plotting_canvas, self)
        self.data = data_class.Data()
        self.model = None
        self.data_imported = False
        self.parameters_imported = False
        self.parameter_value_labels = []
        self.variable_value_labels = []
        self.init_UI()
    

    def init_UI(self):
        # Import data button
        self.import_data_button = QPushButton('Import Data', self)
        self.import_data_button.setToolTip('Import pressure measurements and time data.')
        self.import_data_button.clicked.connect(self.plot_data)

        # Fit model button
        self.fit_button = QPushButton('Fit Curve', self)
        self.fit_button.setToolTip('Fit the permeability and porosity to match the data')
        self.fit_button.clicked.connect(self.fit_data)
        self.fit_button.setEnabled(False)

        button_layout = QVBoxLayout()
        button_layout.addWidget(self.import_data_button)
        button_layout.addWidget(self.fit_button)
        # Model Dropdown menu
        # models_label = QLabel(self)
        # models_label.setText("Select Model: ")
        # self.models_dropdown = QComboBox(self)
        # # models_dropdown.addItem("Select Model")
        # self.models_dropdown.addItem("Theis Solution")
        # # models_dropdown.activated[str].connect(self.select_model)

        # Define the layout
        self.layout = QGridLayout()
        self.layout.addWidget(self.plot_toolbar, 0, 0)
        self.layout.addWidget(self.plotting_canvas, 1, 0)
        # self.layout.addWidget(models_label, 0, 1)
        # self.layout.addWidget(self.models_dropdown, 0, 2)
        # self.layout.addWidget(self.import_data_button, 0, 1)
        # self.layout.addWidget(self.fit_button, 2, 1)
        self.layout.addLayout(button_layout, 0, 1)
        # layout.addLayout(self.parameters_layout(parameter_names=self.parameter_names, default_values=self.default_values), 1, 1)
        self.layout.addLayout(self.create_parameters_groupbox(),1,1)
        # Prioritise the canvas to stretch
        self.layout.setColumnStretch(0, 1)
        print('First row stretch: {} Second row stretch: {} Third row stretch: {}'.format(self.layout.rowStretch(0), self.layout.rowStretch(1), self.layout.rowStretch(2)))
        print('First column stretch: {} Second column stretch: {}'.format(self.layout.columnStretch(0), self.layout.columnStretch(1)))
        # Make the canvas at least 500x500 pixels
        self.layout.setColumnMinimumWidth(0, 500)
        self.layout.setRowMinimumHeight(1, 500)
        self.setLayout(self.layout)

    # @pyqtSlot()
    # def select_model(self, model_type):
    #     if model_type == "Theis Solution":
    #         self.model = model.Theis_Solution(self.data)

    @pyqtSlot()
    def plot_data(self):
        self.clear_all_parameters()
        filename, _ = QFileDialog.getOpenFileName(self, 'Load data file', "", '*.dat;*.txt')
        if filename:
            start = time.time()
            if self.data:
                self.data.read_file(filename=filename)
            else:
                self.data = data_class.Data(filename=filename)
            end = time.time()
            print('Time elapsed = {}'.format(end - start))

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

    @pyqtSlot()
    def fit_data(self):
        if self.data and self.data.parameters:
            self.model = model.Theis_Solution(self.data)
            self.model.find_model_parameters()
            self.plotting_canvas.plot_fit(self.data)
            for i, label in enumerate(self.variable_value_labels):
                if i != len(self.data.variables)-1:
                    label.setText('{:.6f}'.format(self.data.variables[i]))
                else:
                    label.setText('{:.6E}'.format(self.data.variables[i]))
        else:
            error_message = QMessageBox()
            error_message.setIcon(QMessageBox.Critical)
            error_message.setText('Error')
            error_message.setInformativeText('Please enter Theis solution model parameters first')
            error_message.setWindowTitle('Error')
            error_message.exec_()

    def create_parameters_groupbox(self):
        parameters_groupbox = QGroupBox("Known Parameters")
        parameters_groupbox.setAlignment(4)
        # parameters_groupbox.setCheckable(True)
        # parameters_groupbox.setChecked(False)

        parameters_grid = QGridLayout()
        parameters_grid.addWidget(QLabel('Initial Pressure (Pa): '), 0, 0)
        parameters_grid.addWidget(QLabel('Mass Flowrate (Kg/s) '), 1, 0)
        parameters_grid.addWidget(QLabel('Thickness (m): '), 2, 0)
        parameters_grid.addWidget(QLabel('Density (kg/m3): '), 3, 0)
        parameters_grid.addWidget(QLabel('Kinematic Viscosity (m2/s): '), 4, 0)
        parameters_grid.addWidget(QLabel('Compressibility (1/Pa): '), 5, 0)
        parameters_grid.addWidget(QLabel('Radius (m): '), 6, 0)

        # if self.data.parameters[0]:
        #     for i, parameter in enumerate(self.data.parameters):
        #         parameters_grid.addWidget(QLabel(str(parameter)), i, 1)
        # else:
        for i in range(parameters_grid.rowCount()):
            new_label = QLabel()
            self.parameter_value_labels.append(new_label)
            parameters_grid.addWidget(new_label, i, 1)
        
        parameters_groupbox.setLayout(parameters_grid)

        variables_groupbox = QGroupBox('Unknown Parameters')
        variables_groupbox.setAlignment(4)

        variables_grid = QGridLayout()
        variables_grid.addWidget(QLabel('Porosity: '), 0, 0)
        variables_grid.addWidget(QLabel('Permeability (m2): '), 1, 0)
        for i in range(variables_grid.rowCount()):
            new_label = QLabel()
            self.variable_value_labels.append(new_label)
            variables_grid.addWidget(new_label, i, 1)
        # variables_grid.setVerticalSpacing()
        variables_groupbox.setLayout(variables_grid)

        fullWidget = QVBoxLayout()
        fullWidget.addWidget(parameters_groupbox)
        fullWidget.addWidget(variables_groupbox)
        fullWidget.addSpacerItem(QSpacerItem(220, 0, hPolicy=QSizePolicy.Expanding, vPolicy=QSizePolicy.Expanding))
        return fullWidget

    def update_parameter_labels(self):
        for i, label in enumerate(self.parameter_value_labels):
            label.setText(str(self.data.parameters[i]))
    
    def clear_all_parameters(self):
        for label in self.parameter_value_labels + self.variable_value_labels:
            label.setText('')

    # @pyqtSlot()
    # def save_parameters(self):
    #     parameters = []
    #     parameters.append(float(self.p0_input.text()))
    #     parameters.append(float(self.qm_input.text()))
    #     parameters.append(float(self.h_input.text()))
    #     parameters.append(float(self.rho_input.text()))
    #     parameters.append(float(self.nu_input.text()))
    #     parameters.append(float(self.C_input.text()))
    #     parameters.append(float(self.r_input.text()))
    #     self.data.set_known_parameters(parameters)
    #     self.parameters_imported = True
    #     if self.parameters_imported and self.data_imported:
    #         self.fit_button.setEnabled(True)

    # @pyqtSlot()
    # def check_box(self, frame):
    #     if self.params_checkbox.isChecked:
    #         frame.show()
    #     else:
    #         frame.hide()

    # def parameters_layout(self):
    #     # phi_label = QLabel('Porosity: ')
    #     # k_label = QLabel('Permeability: ')
    #     # phi_input = QLineEdit()
    #     # k_input = QLineEdit()
    #     # frame_layout = QGridLayout()
    #     # frame_layout.addWidget(phi_label, 0, 0)
    #     # frame_layout.addWidget(phi_input, 0, 1)
    #     # frame_layout.addWidget(k_label, 1, 0)
    #     # frame_layout.addWidget(k_input, 1, 1)

    #     # frame = QFrame()
    #     # frame.hide()
    #     # frame.setLayout(frame_layout)
    #     # self.params_checkbox = QCheckBox('Enter Porosity and Permeability Estimates')
    #     # self.params_checkbox.stateChanged.connect(lambda: self.check_box(frame))

    #     p0_label = QLabel('Initial Pressure')
    #     self.p0_input = QLineEdit()
    #     self.p0_input.setText('3.6e6')
    #     qm_label = QLabel('Mass Flowrate')
    #     self.qm_input = QLineEdit()
    #     self.qm_input.setText('-0.005')
    #     h_label = QLabel('Thickness')
    #     self.h_input = QLineEdit()
    #     self.h_input.setText('100')
    #     rho_label = QLabel('Density')
    #     self.rho_input = QLineEdit()
    #     self.rho_input.setText('813.37')
    #     nu_label = QLabel('Kinematic Viscosity')
    #     self.nu_input = QLineEdit()
    #     self.nu_input.setText('0.0001111')
    #     C_label = QLabel('Compressibility')
    #     self.C_input = QLineEdit()
    #     self.C_input.setText('0.001303')
    #     r_label = QLabel('Radius')
    #     self.r_input = QLineEdit()
    #     self.r_input.setText('0.05')
        
    #     input_button = QPushButton('Save Parameters', self)
    #     input_button.setToolTip('Import the above entered parameters')
    #     input_button.clicked.connect(self.save_parameters)

    #     layout = QGridLayout()
    #     # layout.addWidget(self.params_checkbox, 0, 0)
    #     # layout.addWidget(frame, 1, 0)
    #     # layout.addLayout(k_widget, 1, 1)
    #     layout.addWidget(p0_label, 2, 0)
    #     layout.addWidget(self.p0_input, 2, 1)
    #     layout.addWidget(qm_label, 3, 0)
    #     layout.addWidget(self.qm_input, 3, 1)
    #     layout.addWidget(h_label, 4, 0)
    #     layout.addWidget(self.h_input, 4, 1)
    #     layout.addWidget(rho_label, 5, 0)
    #     layout.addWidget(self.rho_input, 5, 1)
    #     layout.addWidget(nu_label, 6, 0)
    #     layout.addWidget(self.nu_input, 6, 1)
    #     layout.addWidget(C_label, 7, 0)
    #     layout.addWidget(self.C_input, 7, 1)
    #     layout.addWidget(r_label, 8, 0)
    #     layout.addWidget(self.r_input, 8, 1)
    #     layout.addWidget(input_button, 9, 1)

    #     # layout = QGridLayout()
    #     # layout.setVerticalSpacing(1)
    #     # row = 0
    #     # for parameter, default_value in zip(parameter_names, default_values):
    #     #     label = QLabel(parameter)
    #     #     input_box = QLineEdit()
    #     #     if default_value:
    #     #         input_box.setText(str(default_value))
    #     #     layout.addWidget(label, row, 0)
    #     #     layout.addWidget(input_box, row, 1)
    #     #     row += 1
    #     # input_button = QPushButton('Save Parameters', self)
    #     # input_button.setToolTip('Import the above entered parameters')
    #     # input_button.clicked.connect(self.save_parameters)
    #     # layout.addWidget(input_button, row, 0)
    #     return layout

class PlottingCanvas(FigureCanvas):
    def __init__(self, parent=None):
        # self.figure = plt.figure()
        self.figure = Figure()
        self.axes = None
        FigureCanvas.__init__(self, self.figure)
        self.setParent(parent)
        # self.plot_toolbar = NavigationToolbar(self.plotting_canvas, self)
        self.fitted_lines = []
        self.draw()
    
    def plot_observed_data(self, data):
        # Bottleneck here is in getting the layout of the plot set properly (either tight layout or draw)
        start = time.clock()
        # self.figure.clear()
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