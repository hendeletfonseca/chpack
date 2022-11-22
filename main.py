from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QGridLayout, QPushButton, QLineEdit, QSizePolicy
from PyQt5 import QtCore
import sys
from pyViewerApp import pyViewerApp

ORGANIZATION_NAME = 'LCC-IC-UFF'
ORGANIZATION_DOMAIN = ''
APPLICATION_NAME = 'PyViewer'
SETTINGS_TRAY = 'settings/tray'

def main():
    # To ensure that every time you call QSettings not enter the data of your application, 
    # which will be the settings, you can set them globally for all applications   
    QtCore.QCoreApplication.setApplicationName(ORGANIZATION_NAME)
    QtCore.QCoreApplication.setOrganizationDomain(ORGANIZATION_DOMAIN)
    QtCore.QCoreApplication.setApplicationName(APPLICATION_NAME)

    # create pyqt5 app
    app = QApplication(sys.argv)
    app.setStyle('Fusion')

    # create the instance of our Window
    mw = pyViewerApp(APPLICATION_NAME)

    # showing all the widgets
    mw.show()

    # start the app
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()