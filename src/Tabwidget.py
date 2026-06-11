
from PyQt5.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QWidget,
    QTabWidget,
    QPushButton,
)

class MyTabWidget(QWidget):
    def __init__(self, parent):
        super().__init__(parent)
        self._layout = QVBoxLayout(self)
        # Initialize tab screen
        self.master     = QTabWidget()

        self.master.resize(300, 200) 
        self.tablabels=['setup','scatter','raw','resid','PSA','alpha(f)','alpha(t)','V1bal','msg','last status','config']
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
        scatter_grid_widget = QWidget()
        scatter_grid = QGridLayout(scatter_grid_widget)
        scatter_grid.setContentsMargins(0, 0, 0, 0)
        for i in range(2):
            for j in range(2):
                scatter_grid.addWidget(parent.scatterplots[i, j], i, j)
        self.mytabs[1].setLayout(QGridLayout())
        self.mytabs[1].layout().addWidget(scatter_grid_widget, 0, 0)
        glayout =  QGridLayout()
        self.mytabs[2].setLayout(glayout)
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.rawplots[i,j],i,j)
        glayout =  QGridLayout()
        self.mytabs[3].setLayout(glayout)
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.residplots[i,j],i,j)
        psa_vbox = QVBoxLayout()
        self.mytabs[4].setLayout(psa_vbox)
        psa_btn_row = QHBoxLayout()
        parent.psa_resid_btn = QPushButton('Show: Residuals')
        parent.psa_drift_btn = QPushButton('Drift sub: Off')
        psa_btn_row.addWidget(parent.psa_resid_btn)
        psa_btn_row.addWidget(parent.psa_drift_btn)
        psa_btn_row.addStretch()
        psa_vbox.addLayout(psa_btn_row)
        psa_grid_widget = QWidget()
        psa_grid = QGridLayout(psa_grid_widget)
        psa_grid.setContentsMargins(0, 0, 0, 0)
        for i in range(2):
            for j in range(2):
                psa_grid.addWidget(parent.psaplots[i, j], i, j)
        psa_vbox.addWidget(psa_grid_widget)
        glayout =  QGridLayout()
        self.mytabs[5].setLayout(glayout)
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.alphafplots[i,j],i,j)

        glayout =  QGridLayout()
        self.mytabs[6].setLayout(glayout)
        for j in range(2):
            glayout.addWidget(parent.etaplots[0,j],0,j)

        glayout =  QGridLayout()
        self.mytabs[7].setLayout(glayout)
        for i in range(2):
            for j in range(2):
                glayout.addWidget(parent.balanceplots[i,j],i,j)

        glayout =  QGridLayout()
        self.mytabs[8].setLayout(glayout)
        glayout.addWidget(parent.output,0,0)

        glayout =  QGridLayout()
        self.mytabs[9].setLayout(glayout)
        glayout.addWidget(parent.mstatus,0,0)

        glayout =  QGridLayout()
        self.mytabs[10].setLayout(glayout)
        glayout.addWidget(parent.config_editor, 0, 0)

        glayout.setRowStretch(0, 1)

        

        self.master.currentChanged.connect(parent.replot)            
