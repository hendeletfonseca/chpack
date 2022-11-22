import os
import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QGridLayout, QPushButton, QLineEdit, QSizePolicy

class pyViewerApp(QMainWindow):
    def __init__(self, APPLICATION_NAME, parent=None):
        super().__init__(parent)
        self.setWindowTitle(APPLICATION_NAME)
        self.setFixedSize(400, 400)
        self.cw = QWidget()
        self.grid = QGridLayout(self.cw)

        self.display = QLineEdit()
        self.grid.addWidget(self.display, 0,0,1,5)
        self.display.setStyleSheet(
            '* {background: white; color: #000; font-size: 30px; font-weight: bold;}'
        )
        #self.display.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)

        self.add_btn(QPushButton("Pytomoviewer"),1,0,1,1,
            lambda: os.system("cd pytomoviewer && python -c 'import pytomoviewer as ptv; ptv.run()'"),
            '* {font-size: 10px;}'
        )
        self.add_btn(QPushButton("ch2pp"),1,1,1,1,
            lambda: os.system(f"julia src/ch2pp.jl inp/{self.display.text()}")
        )
        self.add_btn(QPushButton("temp_view"),1,2,1,1,
            lambda: os.system(f"python pyview/temp_view.py")
        )
        self.add_btn(QPushButton("flow_view"),1,3,1,1,
            lambda: os.system(f"python pyview/flow_view.py")
        )
        
        self.setCentralWidget(self.cw)

    def add_btn(self, btn, row, col, rowspan, colspan, function, style=None):
        self.grid.addWidget(btn, row, col, rowspan, colspan)
        btn.clicked.connect(function)

        if (style):
            btn.setStyleSheet(style)

        #btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)