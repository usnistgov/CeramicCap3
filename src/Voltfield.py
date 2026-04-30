
from PyQt5.QtWidgets import QDoubleSpinBox


class VoltField(QDoubleSpinBox):
    def __init__(self, parent, val=10):
        super().__init__(parent)
        self.setMinimumWidth(100)
        self.setMinimum(-10.0)
        self.setMaximum(10.0)
        self.setSingleStep(0.001)
        self.setDecimals(4)
        self.setValue(val)
        self.valueChanged.connect(parent.newValues)
