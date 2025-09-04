
from PyQt5.QtWidgets import (
    QWidget,
    QDoubleSpinBox,
)


class VoltField(QDoubleSpinBox):
    def __init__(self, parent,val=10):
        super(QWidget, self).__init__(parent)  
        self.setMinimumWidth(100)
        self.setMinimum(-10.0)
        self.setMaximum(10.0)
        self.setSingleStep(0.001)
        self.setDecimals(4) 
        self.setValue(val)  # Set initial value here
        self.parent =parent
    def conn(self):
        self.valueChanged.connect(self.parent.newValues)
