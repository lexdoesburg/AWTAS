from PyQt5.QtWidgets import QWidget, QTabWidget, QPushButton, QVBoxLayout
from PyQt5.QtCore import pyqtSlot

# import plotting_widget1 as plotting_widget
import awtas.gui.plotting_widget as plotting_widget

class Tabs(QWidget):
    def __init__(self, parent=None):
        super(Tabs, self).__init__(parent)
        self.layout = QVBoxLayout()
        self.tabs = QTabWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.remove_tab)
        self.tabs.setMovable(True)
        self.main_tab = plotting_widget.PlotWidget(self.tabs)
        self.tabs.addTab(self.main_tab, 'Analysis 1')
        self.new_tab_button = QPushButton('New')
        self.new_tab_button.setToolTip('Open a new analysis.')
        self.new_tab_button.clicked.connect(self.add_new_tab)
        self.tabs.setCornerWidget(self.new_tab_button)
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
    
    @pyqtSlot()
    def add_new_tab(self):
        new_tab = plotting_widget.PlotWidget(self)
        self.tabs.addTab(new_tab, 'Analysis {}'.format(self.tabs.count()+1))
        self.tabs.setCurrentIndex(self.tabs.count()-1)
    
    @pyqtSlot(int)
    def remove_tab(self, index):
        widget = self.tabs.widget(index)
        if widget is not None:
            widget.deleteLater()
        self.tabs.removeTab(index)
        if self.tabs.count() == 0:
            self.tabs.addTab(plotting_widget.PlotWidget(self), 'Analysis 1')