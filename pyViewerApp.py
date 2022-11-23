import os
import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QGridLayout, QPushButton, QLineEdit, QSizePolicy
import subprocess

class pyViewerApp(QMainWindow):

    def __init__(self, APPLICATION_NAME, parent=None):
        super().__init__(parent)
        self.setWindowTitle(APPLICATION_NAME)
        self.setFixedSize(600, 600)
        self.cw = QWidget()
        self.grid = QGridLayout(self.cw)

        self.display = QLineEdit()
        self.grid.addWidget(self.display, 0,0,1,5)
        self.display.setStyleSheet(
            '* {background: white; color: #000; font-size: 30px; font-weight: bold;}'
        )
        #self.display.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)

        '''self.add_btn(QPushButton("Pytomoviewer"),1,0,1,1, 
            run_ptv(self),
            '* {font-size: 10px;}'
        )'''
        self.add_btn(QPushButton("Pytomoviewer"),1,0,1,1,
            self.run_ptv(),
            '* {font-size: 10px;}'
        )
        self.add_btn(QPushButton("ch2pp"),1,1,1,1,
            lambda: subprocess.Popen(f"julia src/ch2pp.jl inp/{self.display.text()}", stdout=subprocess.PIPE, shell=True)
            #lambda: os.system(f"julia src/ch2pp.jl inp/{self.display.text()}")
        )
        self.add_btn(QPushButton("temp_view"),1,2,1,1,
            lambda: subprocess.Popen(f"python pyview/temp_view.py", stdout=subprocess.PIPE, shell=True)
        )
        self.add_btn(QPushButton("flow_view"),1,3,1,1,
            lambda: subprocess.Popen(f"python pyview/flow_view.py", stdout=subprocess.PIPE, shell=True)
        )
        
        self.setCentralWidget(self.cw)

    def add_btn(self, btn, row, col, rowspan, colspan, function, style=None):
        self.grid.addWidget(btn, row, col, rowspan, colspan)
        btn.clicked.connect(function)

        if (style):
            btn.setStyleSheet(style)

        #btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)

    def run_ptv(self):
        
        self.error = False
        
        try:
            self.os.system("python -c 'import pytomoviewer as ptv; ptv.run()'")
        except:
            self.os.system("pip install setuptools && pip install -U scikit-image && pip install git+https://github.com/LCC-UFF/pytomoviewer")
            error = True
        
        if error:
            os.system("python -c 'import pytomoviewer as ptv; ptv.run()'")
        
        return