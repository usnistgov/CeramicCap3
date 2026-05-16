
from PyQt5.QtWidgets import (
    QVBoxLayout,
    QGridLayout,
    QWidget,
    QTabWidget,
)

class MyTabWidget(QWidget):
    def __init__(self, parent):
        super().__init__(parent)
        self._layout = QVBoxLayout(self)
        # Initialize tab screen
        self.master     = QTabWidget()

        self.master.resize(300, 200) 
        self.tablabels=['setup','scatter','raw','PSA','alpha(f)','eta','msg','last status','config']
        self.mytabs =[]
        for l in self.tablabels:
            self.mytabs.append(QWidget())
        for n in range(len(self.tablabels)):
            self.master.addTab(self.mytabs[n],\
                    self.tablabels[n]) 

        self._layout.addWidget(self.master)
        self.setLayout(self._layout)
        
        glayout =  QGridLayout()
        self.mytabs[0].setLayout(glayout)
        glayout.addWidget(parent.circuit_setup, 0, 0)
        glayout =  QGridLayout()
        self.mytabs[1].setLayout(glayout)
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.scatterplots[i,j],i,j)
        glayout =  QGridLayout()
        self.mytabs[2].setLayout(glayout)
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.rawplots[i,j],i,j)
        glayout =  QGridLayout()
        self.mytabs[3].setLayout(glayout)
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.psaplots[i,j],i,j)
        glayout =  QGridLayout()
        self.mytabs[4].setLayout(glayout)
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.alphafplots[i,j],i,j)

        glayout =  QGridLayout()
        self.mytabs[5].setLayout(glayout)
        for j in range(2):
            glayout.addWidget(parent.etaplots[0,j],0,j)

        glayout =  QGridLayout()
        self.mytabs[6].setLayout(glayout)
        glayout.addWidget(parent.output,0,0)

        glayout =  QGridLayout()
        self.mytabs[7].setLayout(glayout)
        glayout.addWidget(parent.mstatus,0,0)

        glayout =  QGridLayout()
        self.mytabs[8].setLayout(glayout)
        glayout.addWidget(parent.config_editor,0,0)

        

        self.master.currentChanged.connect(parent.replot)            
