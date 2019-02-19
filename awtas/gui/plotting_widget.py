from PyQt5.QtWidgets import QWidget, QPushButton, QVBoxLayout, QFileDialog, QMessageBox, QGroupBox, QGridLayout, QLabel, QSpacerItem, QSizePolicy, QComboBox
from PyQt5.QtCore import pyqtSlot

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import numpy as np

import awtas.logic.data as data_class
import awtas.logic.model as model
# import data as data_class
# import model

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
        # Choose model type combobox
        self.model_type_combobox = QComboBox(self)
        self.model_type_combobox.addItems(['Analytical Theis','Homogeneous Porous'])
        self.model_type_combobox.activated[str].connect(self.model_selection)
        model_type_label = QLabel('Model Type: ')
        model_type_layout = QGridLayout()
        model_type_layout.addWidget(model_type_label,0,0)
        model_type_layout.addWidget(self.model_type_combobox,0,1)
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
        button_layout = QVBoxLayout()
        button_layout.addLayout(model_type_layout)
        button_layout.addWidget(self.import_data_button)
        button_layout.addWidget(self.fit_button)

        # Define the overall widget layout
        self.layout = QGridLayout()
        self.layout.addWidget(self.plotting_canvas, 1, 0)
        self.layout.addWidget(self.plot_toolbar, 0, 0)
        # self.layout.addWidget(models_label, 0, 1)
        # self.layout.addWidget(self.models_dropdown, 0, 2)
        # self.layout.addWidget(self.import_data_button, 0, 1)
        # self.layout.addWidget(self.fit_button, 2, 1)
        self.layout.addLayout(button_layout, 0, 1)
        # layout.addLayout(self.parameters_layout(parameter_names=self.parameter_names, default_values=self.default_values), 1, 1)
        self.parameter_sidebar = self.create_parameters_groupbox()
        self.layout.addLayout(self.parameter_sidebar, 1, 1)

        # Make the canvas at least 500x500 pixels
        # Canvas is good at 500x500 windows. 700x500 mac.
        self.layout.setColumnMinimumWidth(0, 700)
        self.layout.setRowMinimumHeight(1, 500)

        # Prioritise the canvas to stretch over other components
        self.layout.setColumnStretch(0, 1)
        print('First row stretch: {} Second row stretch: {} Third row stretch: {}'.format(self.layout.rowStretch(0), self.layout.rowStretch(1), self.layout.rowStretch(2)))
        print('First column stretch: {} Second column stretch: {}'.format(self.layout.columnStretch(0), self.layout.columnStretch(1)))

        # Set the widget layout
        self.setLayout(self.layout)

        # self.show()

    # @pyqtSlot()
    # def select_model(self, model_type):
    #     if model_type == "Theis Solution":
    #         self.model = model.Theis_Solution(self.data)

    @pyqtSlot(str)
    def model_selection(self, text):
        """
        Function to change what happens when a different model type is selected
        """
        print(text)
        if text == 'Analytical Theis':
            self.data = data_class.Data(model_type='theis')
        elif text == 'Homogeneous Porous':
            self.data = data_class.Data(model_type='radial1d')


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
        if self.data and self.data.parameters:
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

    def create_parameters_groupbox(self):
        parameters_groupbox = QGroupBox("Known Parameters")
        parameters_groupbox.setAlignment(4)
        # parameters_groupbox.setCheckable(True)
        # parameters_groupbox.setChecked(False)

        parameters_grid = QGridLayout()
        for i, (parameter, info) in enumerate(self.data.parameters.items()):
            parameters_grid.addWidget(QLabel('{} [{}]: '.format(parameter, info['Units'])), i, 0)
            new_label = QLabel()
            self.parameter_value_labels.append(new_label)
            parameters_grid.addWidget(new_label, i, 1)
        parameters_groupbox.setLayout(parameters_grid)

        variables_groupbox = QGroupBox('Unknown Parameters')
        variables_groupbox.setAlignment(4)

        variables_grid = QGridLayout()
        variables_grid.addWidget(QLabel('Porosity: '), 0, 0)
        variables_grid.addWidget(QLabel('Permeability [m<sup>2</sup>]: '), 1, 0)
        for i in range(variables_grid.rowCount()):
            new_label = QLabel()
            self.variable_value_labels.append(new_label)
            variables_grid.addWidget(new_label, i, 1)
        variables_groupbox.setLayout(variables_grid)

        fullWidget = QVBoxLayout()
        fullWidget.addWidget(parameters_groupbox)
        fullWidget.addWidget(variables_groupbox)
        # Spaceritem is good at 220, 0 for windows. 300, 0 for mac.
        fullWidget.addSpacerItem(QSpacerItem(300, 0, hPolicy=QSizePolicy.Expanding, vPolicy=QSizePolicy.Expanding))
        return fullWidget

    def update_parameter_labels(self):
        # for i, label in enumerate(self.parameter_value_labels):
        #     label.setText(str(self.data.parameters[i]))
        for i, info in enumerate(self.data.parameters.values()):
            self.parameter_value_labels[i].setText(str(info['Value']))
        # pass
    
    def clear_all_parameters(self):
        for label in self.parameter_value_labels + self.variable_value_labels:
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
