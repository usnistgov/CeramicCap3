import configparser
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QScrollArea,
    QLabel, QLineEdit, QFormLayout, QPushButton, QGroupBox
)
from PyQt5.QtCore import pyqtSignal


class ConfigEditor(QWidget):
    configSaved = pyqtSignal()

    def __init__(self, cfgpath, parent=None):
        super().__init__(parent)
        self.cfgpath = cfgpath
        self._fields = {}

        outer = QVBoxLayout(self)

        note = QLabel('Note: saving here removes comments from config.ini')
        note.setStyleSheet('color: gray; font-style: italic;')
        outer.addWidget(note)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        inner = QWidget()
        self._inner_layout = QVBoxLayout(inner)
        self._inner_layout.setSpacing(8)
        scroll.setWidget(inner)
        outer.addWidget(scroll)

        btn_row = QHBoxLayout()
        btn_save = QPushButton('Save')
        btn_reload = QPushButton('Reload')
        btn_row.addWidget(btn_save)
        btn_row.addWidget(btn_reload)
        btn_row.addStretch()
        outer.addLayout(btn_row)

        btn_save.clicked.connect(self.save)
        btn_reload.clicked.connect(self.reload)

        self._build()

    def _build(self):
        cp = configparser.ConfigParser()
        cp.read(self.cfgpath)
        for section in cp.sections():
            box = QGroupBox(f'[{section}]')
            form = QFormLayout(box)
            self._fields[section] = {}
            for key in cp[section]:
                le = QLineEdit(cp[section][key].replace('\n', ' ').strip())
                form.addRow(key, le)
                self._fields[section][key] = le
            self._inner_layout.addWidget(box)
        self._inner_layout.addStretch()

    def reload(self):
        cp = configparser.ConfigParser()
        cp.read(self.cfgpath)
        for section, keys in self._fields.items():
            for key, le in keys.items():
                if cp.has_option(section, key):
                    le.setText(cp[section][key].replace('\n', ' ').strip())

    def save(self):
        cp = configparser.ConfigParser()
        cp.read(self.cfgpath)
        for section, keys in self._fields.items():
            for key, le in keys.items():
                cp[section][key] = le.text()
        with open(self.cfgpath, 'w') as f:
            cp.write(f)
        self.configSaved.emit()
