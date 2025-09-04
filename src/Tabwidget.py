
from PyQt5.QtGui import QFont,QIcon
from PyQt5.QtWidgets import (
    QVBoxLayout,
    QGridLayout,
    QWidget,
    QTabWidget,
)

class MyTabWidget(QWidget): 
    def __init__(self, parent): 
        super(QWidget, self).__init__(parent)  
        self.layout = QVBoxLayout(self) 
        # Initialize tab screen 
        self.master     = QTabWidget() 

        self.master.resize(300, 200) 
        self.tablabels=['scatter','raw','alpha','circles','msg']
        self.mytabs =[]
        for l in self.tablabels:
            self.mytabs.append(QWidget())
        for n in range(len(self.tablabels)):
            self.master.addTab(self.mytabs[n],\
                    self.tablabels[n]) 

        self.layout.addWidget(self.master) 
        self.setLayout(self.layout) 
        
        glayout =  QGridLayout()
        self.mytabs[0].setLayout(glayout) 
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.scatterplots[i,j],i,j)    
        glayout =  QGridLayout()
        self.mytabs[1].setLayout(glayout) 
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.rawplots[i,j],i,j)                
        glayout =  QGridLayout()
        self.mytabs[2].setLayout(glayout) 
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.abplots[i,j],i,j)                
        glayout =  QGridLayout()
        self.mytabs[3].setLayout(glayout) 
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.ciplots[i,j],i,j)                

        glayout =  QGridLayout()
        self.mytabs[-1].setLayout(glayout) 
        glayout.addWidget(parent.output,0,0)    

        self.master.currentChanged.connect(parent.replot)            
