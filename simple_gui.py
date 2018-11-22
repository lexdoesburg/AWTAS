import sys
from PyQt5.QtWidgets import QApplication, QWidget, QDesktopWidget
from PyQt5.QtGui import QIcon

class App(QWidget):
 
    def __init__(self):
        super().__init__()
        self.title = 'AWTAS'
        self.left = 0
        self.top = 30
        self.scale = 2
        self.width = 640*self.scale
        self.height = 480*self.scale
        self.initUI()
 
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.show()
 
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())