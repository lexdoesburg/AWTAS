import sys

from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QFileDialog, QGridLayout
from PyQt5.QtCore import pyqtSlot
# from PyQt5.QtGui import QIcon
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

import numpy as np

import theis_solution as ts
# import for testing performance
import time

class AWTAS_App(QWidget):
    def __init__(self, parent=None):
        super(AWTAS_App, self).__init__(parent)
        self.title = 'AWTAS'
        self.left = 0
        self.top = 30
        # self.scale = 2
        self.width = 640
        self.height = 480
        self.figure = plt.figure()
        self.plotting_canvas = FigureCanvas(self.figure)
        self.plot_toolbar = NavigationToolbar(self.plotting_canvas, self)
        self.init_UI()
 

    def init_UI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        import_data_button = QPushButton('Import Data', self)
        import_data_button.setToolTip('Import pressure measurements and time data.')
        import_data_button.clicked.connect(self.plot_data)
        
        fit_button = QPushButton('Fit Curve', self)
        fit_button.setToolTip('Fit the permeability and porosity to match the data')
        fit_button.clicked.connect(self.fit_data)

        layout = QGridLayout()
        layout.addWidget(self.plotting_canvas, 0, 0)
        layout.addWidget(self.plot_toolbar, 1, 0)
        layout.addWidget(import_data_button, 0, 1)
        layout.addWidget(fit_button, 0, 2)
        self.setLayout(layout)

        self.show()


    @pyqtSlot()
    def plot_data(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'Load data file', "", '*.txt')
        if filename:
            start = time.time()

            x_data, y_data = ts.read_data(filename)

            end = time.time()
            print('Time elapsed = {}'.format(end - start))

            start = time.time()
            self.figure.clear()
            ax = self.figure.add_subplot(1,1,1)
            ax.semilogx(np.log(x_data), y_data, 'kx', label='Observed Data')
            ax.set_xlabel('Log Time (Sec)')
            ax.set_ylabel('Pressure (Pa)')
            ax.legend(loc='best')

            end = time.time()
            print('Time elapsed = {}'.format(end - start))

            return x_data, y_data

    @pyqtSlot()
    def fit_data(self):
        pass

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = AWTAS_App()
    sys.exit(app.exec_())