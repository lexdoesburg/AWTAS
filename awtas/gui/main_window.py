import sys
from PyQt5.QtWidgets import QApplication, QMainWindow

import awtas.gui.tab_widget as tab_widget

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = 'AWTAS'
        # self.plot_widget = PlotWidget(self)
        self.tabs = tab_widget.Tabs(self)
        self.init_UI()
        
    def init_UI(self):
        self.setWindowTitle(self.title)
        # self.setCentralWidget(self.plot_widget)
        self.setCentralWidget(self.tabs)
        self.show()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MainWindow()
    sys.exit(app.exec_())