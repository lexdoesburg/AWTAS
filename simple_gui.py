import sys

from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QFileDialog
from PyQt5.QtCore import pyqtSlot
# from PyQt5.QtGui import QIcon

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import theis_solution as ts
 
class AWTAS_App(QWidget):
 
    def __init__(self):
        super().__init__()
        self.title = 'AWTAS'
        self.left = 0
        self.top = 30
        self.scale = 2
        self.width = 640*self.scale
        self.height = 480*self.scale
        self.init_UI()
 

    def init_UI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        import_data_button = QPushButton('Import Data', self)
        import_data_button.setToolTip('Import pressure measurements and time data.')
        import_data_button.move(540*self.scale,50)
        filename = import_data_button.clicked.connect(self.get_filename)

        self.show()


    @pyqtSlot()
    def get_filename(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'Load data file', "", '*.txt')
        if filename:
            x_data, y_data = ts.read_data(filename)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = AWTAS_App()
    sys.exit(app.exec_())