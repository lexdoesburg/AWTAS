from PyQt5.QtWidgets import QApplication, QMainWindow, QAction
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot

import awtas.gui.tab_widget as tab_widget

class MainWindow(QMainWindow):
    """
    Main GUI window.
    """
    def __init__(self):
        super().__init__()
        self.title = 'AWTAS'
        self.tab_widget = tab_widget.Tabs(self)
        self.init_UI()
        

    def init_UI(self):
        """
        Initialise the GUI.
        """
        self.setWindowTitle(self.title)
        self.setCentralWidget(self.tab_widget)

        # Create menu bar and file menu
        menu_bar = self.menuBar()
        # menu_bar.setNativeMenuBar(False)
        file_menu = menu_bar.addMenu('&File')
        
        # Import data action
        import_data = QAction('&Import Data', self)
        import_data.setShortcut('Ctrl+L')
        import_data.setStatusTip('Import data from file')
        import_data.triggered.connect(self.import_data)

        # Export/save data action
        export_data = QAction('&Export Results', self)
        export_data.setShortcut('Ctrl+S')
        export_data.setStatusTip('Export results to file')
        export_data.triggered.connect(self.export_results)

        # Exit action adapted from http://zetcode.com/gui/pyqt5/menustoolbars/
        exit = QAction('&Exit', self)
        exit.setShortcut('Ctrl+Q')
        exit.setStatusTip('Exit application')
        exit.triggered.connect(self.close)

        # Add file bar actions
        file_menu.addAction(import_data)
        file_menu.addAction(export_data)
        file_menu.addAction(exit)
        
        self.showMaximized() # self.show()


    def get_current_plot_widget(self):
        """
        Returns the current plotting widget object from the active tab.
        """
        current_tab = self.tab_widget.tabs.currentIndex()
        current_plot_widget = self.tab_widget.tabs.widget(current_tab)
        return current_plot_widget


    @pyqtSlot()
    def import_data(self):
        """
        Imports data from file when file>import data or ctrl + L is pressed.
        """
        current_plot_widget = self.get_current_plot_widget()
        current_plot_widget.import_data_action()


    @pyqtSlot()
    def export_results(self):
        """
        Exports the results to a data file when file>export data or ctrl + S is pressed.
        """
        current_plot_widget = self.get_current_plot_widget()
        current_plot_widget.produce_output_file()