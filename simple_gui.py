import sys

from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QFileDialog, QGridLayout, QAction, qApp, QMenuBar, QMenu, QMessageBox, QLineEdit, QLabel, QInputDialog, QComboBox, QVBoxLayout
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QIcon

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

import numpy as np
# import for testing performance
import time

# import theis_solution as ts
import model
import data as data_class

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


        self.data = None
        self.model = None

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
        import_data_button = QPushButton('Import Data', self)
        import_data_button.setToolTip('Import pressure measurements and time data.')
        import_data_button.clicked.connect(self.plot_data)

        # Fit model button
        fit_button = QPushButton('Fit Curve', self)
        fit_button.setToolTip('Fit the permeability and porosity to match the data')
        fit_button.clicked.connect(self.fit_data)

        # Model Dropdown menu
        models_label = QLabel(self)
        models_label.setText("Select Model: ")
        models_dropdown = QComboBox(self)
        # models_dropdown.addItem("Select Model")
        models_dropdown.addItem("Theis Solution")
        # models_dropdown.activated[str].connect(self.select_model)

        # Define the layout
        layout = QGridLayout()
        layout.addWidget(self.plotting_canvas, 1, 0)
        layout.addWidget(self.plot_toolbar, 0, 0)
        layout.addWidget(models_label, 0, 1)
        layout.addWidget(models_dropdown, 0, 2)
        layout.addWidget(import_data_button, 2, 1)
        layout.addWidget(fit_button, 2, 2)
        layout.addLayout(self.parameters_layout(), 1, 1)
        self.setLayout(layout)

        self.show()

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

    @pyqtSlot()
    def fit_data(self):
        if self.data and self.data.parameters:
            self.model = model.Theis_Solution(self.data)
            self.model.find_model_parameters()
            self.plotting_canvas.plot_fit(self.data)
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
    ex = AWTAS_App()
    sys.exit(app.exec_())