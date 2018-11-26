import sys

from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QFileDialog, QGridLayout, QAction, qApp, QMenuBar, QMenu, QMessageBox, QLineEdit, QLabel, QInputDialog, QComboBox, QDoubleSpinBox
from PyQt5.QtWidgets import QVBoxLayout, QMainWindow
from PyQt5.QtCore import pyqtSlot

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

import numpy as np
# import for testing performance
import time

# import theis_solution as ts
import model
import data as data_class

class Main_Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = 'AWTAS'
        self.central_widget = AWTAS_App()
        self.initUI()
 
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setCentralWidget(self.central_widget)
        mainMenu = self.menuBar() 
        fileMenu = mainMenu.addMenu('File')
        editMenu = mainMenu.addMenu('Edit')
        viewMenu = mainMenu.addMenu('View')
        searchMenu = mainMenu.addMenu('Search')
        toolsMenu = mainMenu.addMenu('Tools')
        helpMenu = mainMenu.addMenu('Help')
 
        exitButton = QAction('Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)
        self.show()


class AWTAS_App(QWidget):
    def __init__(self, parent=None):
        super(AWTAS_App, self).__init__(parent)
        self.title = 'AWTAS'
        # self.left = 0
        # self.top = 30
        # self.width = 640
        # self.height = 480
        self.figure = plt.figure()
        self.axes = None
        # self.plotting_canvas = FigureCanvas(self.figure)

        self.plotting_canvas = PlottingCanvas()
        self.plot_toolbar = NavigationToolbar(self.plotting_canvas, self)

        self.data = data_class.Data()
        self.model = None

        self.data_imported = False
        self.parameters_imported = False

        # self.parameter_names = ['Initial Pressure', 'Mass Flowrate', 'Thickness', 'Density', 'Kinematic Viscosity', 'Compressibility', 'Radius']
        # self.default_values = [3.6e6, -0.005, 100, 813.37, 0.0001111, 0.001303, 0.05]

        self.init_UI()
 

    def init_UI(self):
        self.setWindowTitle(self.title)
        # self.setGeometry(self.left, self.top, self.width, self.height)

        # ---------------------------------- 
        # # Define exit action
        # exit_action = QAction(QIcon('exit.png'), '&Exit', self)        
        # exit_action.setShortcut('Ctrl+Q')
        # exit_action.setStatusTip('Exit application')
        # exit_action.triggered.connect(qApp.quit)

        # # Import Data
        # import_data = QAction('Load Data', self)
        # import_data.setShortcut('Ctrl+L')
        # import_data.setStatusTip('Load observation data from .txt file')
        # import_data.triggered.connect(self.plot_data)

        # # Menu_bar
        # menu_bar = QMenuBar(self)
        # file_menu = menu_bar.addMenu('&File')
        # file_menu.addAction(exit_action)
        # file_menu.addAction(import_data)
        # ----------------------------------

        # Import data button
        self.import_data_button = QPushButton('Import Data', self)
        self.import_data_button.setToolTip('Import pressure measurements and time data.')
        self.import_data_button.clicked.connect(self.plot_data)

        # Fit model button
        self.fit_button = QPushButton('Fit Curve', self)
        self.fit_button.setToolTip('Fit the permeability and porosity to match the data')
        self.fit_button.clicked.connect(self.fit_data)
        self.fit_button.setEnabled(False)

        # Model Dropdown menu
        # models_label = QLabel(self)
        # models_label.setText("Select Model: ")
        # self.models_dropdown = QComboBox(self)
        # # models_dropdown.addItem("Select Model")
        # self.models_dropdown.addItem("Theis Solution")
        # # models_dropdown.activated[str].connect(self.select_model)

        # Define the layout
        self.layout = QGridLayout()
        self.layout.addWidget(self.plotting_canvas, 1, 0)
        self.layout.addWidget(self.plot_toolbar, 0, 0)
        # self.layout.addWidget(models_label, 0, 1)
        # self.layout.addWidget(self.models_dropdown, 0, 2)
        self.layout.addWidget(self.import_data_button, 0, 1)
        self.layout.addWidget(self.fit_button, 2, 1)
        # layout.addLayout(self.parameters_layout(parameter_names=self.parameter_names, default_values=self.default_values), 1, 1)
        self.layout.addLayout(self.parameters_layout(), 1, 1)
        # Prioritise the canvas to stretch
        self.layout.setColumnStretch(0, 1)
        print('First row stretch: {} Second row stretch: {} Third row stretch: {}'.format(self.layout.rowStretch(0), self.layout.rowStretch(1), self.layout.rowStretch(2)))
        print('First column stretch: {} Second column stretch: {}'.format(self.layout.columnStretch(0), self.layout.columnStretch(1)))
        # Make the canvas at least 500x500 pixels
        self.layout.setColumnMinimumWidth(0, 500)
        self.layout.setRowMinimumHeight(1, 500)
        self.setLayout(self.layout)

        # self.show()

    # @pyqtSlot()
    # def select_model(self, model_type):
    #     if model_type == "Theis Solution":
    #         self.model = model.Theis_Solution(self.data)

    @pyqtSlot()
    def plot_data(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'Load data file', "", '*.txt')
        if filename:
            start = time.time()
            if self.data:
                self.data.read_file(filename=filename)
            else:
                self.data = data_class.Data(filename=filename)
            end = time.time()
            print('Time elapsed = {}'.format(end - start))

            start = time.time()
            # self.figure.clear()
            # self.axes = self.figure.add_subplot(1,1,1)
            # self.axes.semilogx(np.log(self.data.time), self.data.observation, 'kx', label='Observed Data')
            # self.axes.set_xlabel('Log Time (s)')
            # self.axes.set_ylabel('Pressure (Pa)')
            # self.axes.legend(loc='best')
            self.plotting_canvas.plot_observed_data(self.data)
            end = time.time()
            print('Time elapsed = {}'.format(end - start))
            print('Drawing')
            # self.plotting_canvas.draw()
            # self.figure.tight_layout()
        self.data_imported = True
        if self.parameters_imported and self.data_imported:
            self.fit_button.setEnabled(True)

    @pyqtSlot()
    def fit_data(self):
        if self.data and self.data.parameters:
            self.model = model.Theis_Solution(self.data)
            self.model.find_model_parameters()
            self.plotting_canvas.plot_fit(self.data)
            parameters_label = QLabel('Porosity: {:.6f}\tPermeability: {:.6E}'.format(self.data.phi, self.data.k))
            self.layout.addWidget(parameters_label, 2, 0)
            # self.axes.semilogx(np.log(self.data.time), self.data.approximation, 'r-', label='Fitted Approximation')
            # self.plotting_canvas.draw()
            
        else:
            error_message = QMessageBox()
            error_message.setIcon(QMessageBox.Critical)
            error_message.setText('Error')
            error_message.setInformativeText('Please enter Theis solution model parameters first')
            error_message.setWindowTitle('Error')
            error_message.exec_()
            # error_message = QErrorMessage()
            # error_message.showMessage('Please import data and enter model parameters first')
            # error_message.exec_()
    
    @pyqtSlot()
    def save_parameters(self):
        parameters = []
        parameters.append(float(self.p0_input.text()))
        parameters.append(float(self.qm_input.text()))
        parameters.append(float(self.h_input.text()))
        parameters.append(float(self.rho_input.text()))
        parameters.append(float(self.nu_input.text()))
        parameters.append(float(self.C_input.text()))
        parameters.append(float(self.r_input.text()))
        self.data.set_known_parameters(parameters)
        self.parameters_imported = True
        if self.parameters_imported and self.data_imported:
            self.fit_button.setEnabled(True)

    def parameters_layout(self):
        p0_label = QLabel('Initial Pressure')
        self.p0_input = QLineEdit()
        self.p0_input.setText('3.6e6')
        qm_label = QLabel('Mass Flowrate')
        self.qm_input = QLineEdit()
        self.qm_input.setText('-0.005')
        h_label = QLabel('Thickness')
        self.h_input = QLineEdit()
        self.h_input.setText('100')
        rho_label = QLabel('Density')
        self.rho_input = QLineEdit()
        self.rho_input.setText('813.37')
        nu_label = QLabel('Kinematic Viscosity')
        self.nu_input = QLineEdit()
        self.nu_input.setText('0.0001111')
        C_label = QLabel('Compressibility')
        self.C_input = QLineEdit()
        self.C_input.setText('0.001303')
        r_label = QLabel('Radius')
        self.r_input = QLineEdit()
        self.r_input.setText('0.05')
        
        input_button = QPushButton('Save Parameters', self)
        input_button.setToolTip('Import the above entered parameters')
        input_button.clicked.connect(self.save_parameters)

        layout = QGridLayout()
        layout.addWidget(p0_label, 0, 0)
        layout.addWidget(self.p0_input, 0, 1)
        layout.addWidget(qm_label, 1, 0)
        layout.addWidget(self.qm_input, 1, 1)
        layout.addWidget(h_label, 2, 0)
        layout.addWidget(self.h_input, 2, 1)
        layout.addWidget(rho_label, 3, 0)
        layout.addWidget(self.rho_input, 3, 1)
        layout.addWidget(nu_label, 4, 0)
        layout.addWidget(self.nu_input, 4, 1)
        layout.addWidget(C_label, 5, 0)
        layout.addWidget(self.C_input, 5, 1)
        layout.addWidget(r_label, 6, 0)
        layout.addWidget(self.r_input, 6, 1)
        layout.addWidget(input_button, 7, 1)

        # layout = QGridLayout()
        # layout.setVerticalSpacing(1)
        # row = 0
        # for parameter, default_value in zip(parameter_names, default_values):
        #     label = QLabel(parameter)
        #     input_box = QLineEdit()
        #     if default_value:
        #         input_box.setText(str(default_value))
        #     layout.addWidget(label, row, 0)
        #     layout.addWidget(input_box, row, 1)
        #     row += 1
        # input_button = QPushButton('Save Parameters', self)
        # input_button.setToolTip('Import the above entered parameters')
        # input_button.clicked.connect(self.save_parameters)
        # layout.addWidget(input_button, row, 0)
        return layout

class PlottingCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.figure = plt.figure()
        self.axes = None
        FigureCanvas.__init__(self, self.figure)
        self.setParent(parent)
        # self.plot_toolbar = NavigationToolbar(self.plotting_canvas, self)
        self.fitted_lines = []
        self.draw()
    
    def plot_observed_data(self, data):
        self.figure.clear()
        self.fitted_lines = []
        self.axes = self.figure.add_subplot(1,1,1)
        self.axes.semilogx(np.log(data.time), data.observation, 'kx', label='Observed Data')
        self.axes.set_xlabel('Log Time (s)')
        self.axes.set_ylabel('Pressure (Pa)')
        self.axes.legend(loc='lower left')
        self.figure.tight_layout(pad = 2)
        self.draw()
    
    def plot_fit(self, data):
        if self.fitted_lines:
            self.axes.lines.remove(self.fitted_lines[0])
        self.fitted_lines = self.axes.semilogx(np.log(data.time), data.approximation, 'r-', label='Fitted Approximation')
        self.axes.legend(loc='lower left')
        self.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Main_Window()
    sys.exit(app.exec_())