import sys
import os
import subprocess
from design import *
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
from PyQt5.QtGui import QPixmap

class App(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        super().setupUi(self)
        self.btnFile.clicked.connect(self.fileDialog)
        self.btnPytomoview.clicked.connect(self.openPytomoview)
        self.btnCh2pp.clicked.connect(self.openCh2pp)
        self.btnTempView.clicked.connect(self.openTempView)
        self.btnFlowView.clicked.connect(self.openFlowView)
        self.filePath = None

    def fileDialog(self):
        file, _ = QFileDialog.getOpenFileName(
            self.centralwidget,
            'Open file',
            '',
            #options=QFileDialog.DontUseNativeDialog
        )
        self.inputOpenFile.setText(file)
        self.filePath = file
        print(self.filePath)
        i = len(self.filePath) - 1
        while (self.filePath[i] != '.'):
            self.filePath = self.filePath[:i]
            i -= 1
        self.filePath = self.filePath[:i]

    def openPytomoview(self):
        subprocess.Popen(f"python -c 'import pytomoviewer as ptv; ptv.run()'", stdout=subprocess.PIPE, shell=True)
        
    def openCh2pp(self):
        if self.filePath:
            subprocess.Popen(f"julia src/ch2pp.jl {self.filePath}", stdout=subprocess.PIPE, shell=True)
        print("finish")

    def openTempView(self):
        subprocess.Popen(f"python pyview/temp_view.py", stdout=subprocess.PIPE, shell=True)

    def openFlowView(self):
        subprocess.Popen(f"python pyview/flow_view.py", stdout=subprocess.PIPE, shell=True)

if __name__ == '__main__':
    qt = QApplication(sys.argv)
    app = App()
    app.show()
    sys.exit(qt.exec_())