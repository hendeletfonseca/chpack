import sys, os, subprocess
from design import *
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QMessageBox
from PyQt5.QtGui import QPixmap
from pyview.temp_view import *
from pyview.flow_view import *

class App(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        super().setupUi(self)

        self.btnFile.clicked.connect(self.fileDialog)
        self.btnPytomoview.clicked.connect(self.openPytomoview)
        self.btnCh2pp.clicked.connect(self.runCh2pp)
        self.btnTempView.clicked.connect(self.runTempView)
        self.btnFlowView.clicked.connect(self.runFlowView)
        
        self.filePath = None
        self.images = []

        self.errorMessage = QMessageBox()

    def fileDialog(self):
        file, _ = QFileDialog.getOpenFileName(
            self.centralwidget,
            'Open file',
            '',
            #options=QFileDialog.DontUseNativeDialog
        )
        self.inputOpenFile.setText(file)

        if self.inputOpenFile.text() != "":
            self.filePath = file
            i = len(self.filePath) - 1
            while (self.filePath[i] != '.'):
                self.filePath = self.filePath[:i]
                i -= 1
            self.filePath = self.filePath[:i]
            print(self.filePath)
        

    def openPytomoview(self):
        subprocess.Popen(f'python -c "import pytomoviewer as ptv; ptv.run()"', stdout=subprocess.PIPE, shell=True)
        
    def runCh2pp(self):
        if self.filePath:
            #subprocess.Popen(f"julia src/ch2pp.jl {self.filePath}", stdout=subprocess.PIPE, shell=True)
            os.system(f"julia src/ch2pp.jl {self.filePath}")
        else:
            self.errorMessage.setIcon(QMessageBox.Critical)
            self.errorMessage.setText("Not Path")
            self.errorMessage.setInformativeText("Please select a file")
            self.errorMessage.setWindowTitle("Error")
            self.errorMessage.exec_()

    def runTempView(self):
        #subprocess.Popen(f"python pyview/temp_view.py", stdout=subprocess.PIPE, shell=True)
        if self.filePath != None:
            #os.system("python pyview/temp_view.py")
            temp_view(self.filePath)
            img = QPixmap("pyview/images/temp_rhs_0.png")
            self.images.append(img)
            #self.original_img = QPixmap("pyview/images/temp_rhs_0.png")
            self.labelImg.setPixmap(self.images[0])
        else:
            self.errorMessage.setIcon(QMessageBox.Critical)
            self.errorMessage.setText("Not Path")
            self.errorMessage.setInformativeText("Please select a file")
            self.errorMessage.setWindowTitle("Error")
            self.errorMessage.exec_()

    def runFlowView(self):
        #subprocess.Popen(f"python pyview/flow_view.py", stdout=subprocess.PIPE, shell=True)
        #os.system("python pyview/flow_view.py")
        if self.filePath != None:
            flow_view(self.filePath)
            img = QPixmap("pyview/images/flow_rhs_0.png")
            self.images.append(img)
            self.labelImg.setPixmap(self.images[1])
        else:
            self.errorMessage.setIcon(QMessageBox.Critical)
            self.errorMessage.setText("Not Path")
            self.errorMessage.setInformativeText("Please select a file")
            self.errorMessage.setWindowTitle("Error")
            self.errorMessage.exec_()

if __name__ == '__main__':
    qt = QApplication(sys.argv)
    app = App()
    app.show()
    sys.exit(qt.exec_())