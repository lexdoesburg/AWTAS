if __name__ == '__main__':
    import sys
    from PyQt5.QtWidgets import QApplication
    
    import awtas.gui.main_window as main_window
    
    app = QApplication(sys.argv)
    ex = main_window.MainWindow()
    sys.exit(app.exec_())