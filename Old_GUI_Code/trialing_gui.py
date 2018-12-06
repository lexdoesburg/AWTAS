import sys
 
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QSizePolicy, QMessageBox, QWidget, QPushButton, QFileDialog
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
 
 
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import numpy as np

import theis_solution as ts

# import for testing performance
import time
 
class App(QWidget):
 
    def __init__(self):
        super().__init__()
        self.left = 0
        self.top = 30
        self.title = 'AWTAS'
        self.width = 640
        self.height = 400
        self.init_UI()
 
    def init_UI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
 
        plot_canvas = PlotCanvas(self, width=5, height=4)
        plot_canvas.move(0,0)
 
        load_data_button = QPushButton('Load Data', self)
        load_data_button.setToolTip('This s an example button')
        load_data_button.move(500,0)
        load_data_button.resize(140,100)
        load_data_button.clicked.connect(lambda: self.plot_data(plot_canvas))

        fit_data_button = QPushButton('Fit Data', self)
        fit_data_button.setToolTip('This s an example button')
        fit_data_button.move(500,100)
        fit_data_button.resize(140,100)
        # fit_data_button.clicked.connect(lambda: self.fit_data(m, x_data, y_data))
 
        self.show()
 

    @pyqtSlot()
    def plot_data(self, plot_canvas):
        filename, _ = QFileDialog.getOpenFileName(self, 'Load data file', "", '*.txt')
        if filename:
            start = time.time()

            x_data, y_data = ts.read_data(filename)

            end = time.time()
            print('Time elapsed = {}'.format(end - start))

            start = time.time()

            plot_canvas.plot_observed_data(x_data, y_data)

            end = time.time()
            print('Time elapsed = {}'.format(end - start))

            return x_data, y_data
    

    @pyqtSlot()
    def fit_data(self, plot_canvas, x_data, y_data):
        pass


class PlotCanvas(FigureCanvas):
 
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(1, 1, 1)
 
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
 
        FigureCanvas.setSizePolicy(self,
                QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.draw()
 
 
    def plot_observed_data(self, x_data, y_data):
        ax = self.figure.add_subplot(1, 1 , 1)
        ax.semilogx(np.log(x_data), y_data, 'kx', label='Observed Data')
        ax.set_xlabel('Log Time (s)')
        ax.set_ylabel('Pressure (Pa)')
        ax.legend(loc='best')
        # ax.plot(x_data, y_data, 'kx', label='Observed Data')
        self.draw()
 
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())