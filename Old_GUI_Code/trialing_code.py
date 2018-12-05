# #!/usr/bin/env python
# #-*- coding:utf-8 -*-

# import random

# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
# from matplotlib.figure import Figure

# from PyQt5 import QtCore #conda install pyqt
# from PyQt5 import QtWidgets

# class MatplotlibWidget(QtWidgets.QWidget):
#     def __init__(self, parent=None):
#         super(MatplotlibWidget, self).__init__(parent)

#         self.figure = Figure()
#         self.canvas = FigureCanvasQTAgg(self.figure)

#         self.axis = self.figure.add_subplot(111)

#         self.layoutVertical = QtWidgets.QVBoxLayout(self)#QVBoxLayout
#         self.layoutVertical.addWidget(self.canvas)

# class ThreadSample(QtCore.QThread):
#     newSample = QtCore.pyqtSignal(list)

#     def __init__(self, parent=None):
#         super(ThreadSample, self).__init__(parent)

#     def run(self):
#         randomSample = random.sample(range(0, 10), 10)

#         self.newSample.emit(randomSample)

# class MyWindow(QtWidgets.QWidget):
#     def __init__(self, parent=None):
#         super(MyWindow, self).__init__(parent)

#         self.pushButtonPlot = QtWidgets.QPushButton(self)
#         self.pushButtonPlot.setText("Plot")
#         self.pushButtonPlot.clicked.connect(self.on_pushButtonPlot_clicked)

#         self.matplotlibWidget = MatplotlibWidget(self)

#         self.layoutVertical = QtWidgets.QVBoxLayout(self)
#         self.layoutVertical.addWidget(self.pushButtonPlot)
#         self.layoutVertical.addWidget(self.matplotlibWidget)

#         self.threadSample = ThreadSample(self)
#         self.threadSample.newSample.connect(self.on_threadSample_newSample)
#         self.threadSample.finished.connect(self.on_threadSample_finished)

#     @QtCore.pyqtSlot()
#     def on_pushButtonPlot_clicked(self):
#         self.samples = 0
#         self.matplotlibWidget.axis.clear()
#         self.threadSample.start()

#     @QtCore.pyqtSlot(list)
#     def on_threadSample_newSample(self, sample):
#         self.matplotlibWidget.axis.plot(sample)
#         self.matplotlibWidget.canvas.draw()

#     @QtCore.pyqtSlot()
#     def on_threadSample_finished(self):
#         self.samples += 1
#         if self.samples <= 2:
#             self.threadSample.start()

# if __name__ == "__main__":
#     import sys

#     app = QtWidgets.QApplication(sys.argv)
#     app.setApplicationName('MyWindow')

#     main = MyWindow()
#     main.resize(666, 333)
#     main.show()

#     sys.exit(app.exec_())


import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QGridLayout

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import time
import random

class Window(QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        self.button = QPushButton('Plot')
        self.button.clicked.connect(self.plot)

        # set the layout
        layout = QGridLayout()
        layout.addWidget(self.toolbar, 0, 0)
        layout.addWidget(self.canvas, 1, 0)
        layout.addWidget(self.button, 0, 1)
        self.setLayout(layout)

    def plot(self):
        ''' plot some random stuff '''
        # random data
        data = [random.random() for i in range(1000)]

        # instead of ax.hold(False)
        self.figure.clear()

        # create an axis
        ax = self.figure.add_subplot(111)

        # discards the old graph
        # ax.hold(False) # deprecated, see above

        # plot data
        ax.plot(data, '*-')
        # refresh canvas
        start = time.clock()
        self.canvas.draw()
        end = time.clock()
        print('Time: ',(end-start))

if __name__ == '__main__':
    app = QApplication(sys.argv)

    main = Window()
    main.show()

    sys.exit(app.exec_())