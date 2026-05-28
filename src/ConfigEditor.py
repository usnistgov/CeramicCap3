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
        btn_reload = QPushButton('Reload')
        btn_restore = QPushButton('Restore original flist')
        btn_row.addWidget(btn_reload)
        btn_row.addWidget(btn_restore)
        btn_row.addStretch()
        outer.addLayout(btn_row)

        btn_reload.clicked.connect(self.reload)
        btn_restore.clicked.connect(self.restore_flist)

        self._build()

    _READONLY_KEYS = {('meas', 'olist')}

    def _build(self):
        cp = configparser.ConfigParser()
        cp.read(self.cfgpath)
        for section in cp.sections():
            box = QGroupBox(f'[{section}]')
            form = QFormLayout(box)
            self._fields[section] = {}
            for key in cp[section]:
                le = QLineEdit(cp[section][key].replace('\n', ' ').strip())
                if (section.lower(), key.lower()) in self._READONLY_KEYS:
                    le.setReadOnly(True)
                    le.setStyleSheet('color: gray; background: #f0f0f0;')
                form.addRow(key, le)
                self._fields[section][key] = le
            self._inner_layout.addWidget(box)
        self._inner_layout.addStretch()
        # Auto-sync fstart whenever freqlist is edited
        meas = self._fields.get('MEAS', self._fields.get('meas', {}))
        flist_le = next((le for k, le in meas.items() if k.lower() == 'freqlist'), None)
        if flist_le:
            flist_le.editingFinished.connect(self._sync_fstart_from_flist)

    def reload(self):
        cp = configparser.ConfigParser()
        cp.read(self.cfgpath)
        for section, keys in self._fields.items():
            for key, le in keys.items():
                if cp.has_option(section, key):
                    le.setText(cp[section][key].replace('\n', ' ').strip())

    def _sync_fstart_from_flist(self):
        meas = self._fields.get('MEAS', self._fields.get('meas', {}))
        flist_key  = next((k for k in meas if k.lower() == 'freqlist'), None)
        fstart_key = next((k for k in meas if k.lower() == 'fstart'),   None)
        if not flist_key or not fstart_key:
            return
        parts = [p.strip() for p in meas[flist_key].text().split(',') if p.strip()]
        if parts:
            try:
                meas[fstart_key].setText(f'{float(parts[0]):.1f}')
            except ValueError:
                pass

    def restore_flist(self):
        meas = self._fields.get('MEAS', self._fields.get('meas', {}))
        olist_key = next((k for k in meas if k.lower() == 'olist'), None)
        flist_key = next((k for k in meas if k.lower() == 'freqlist'), None)
        if olist_key and flist_key:
            meas[flist_key].setText(meas[olist_key].text())
            self._sync_fstart_from_flist()

    def save(self):
        cp = configparser.ConfigParser()
        cp.read(self.cfgpath)
        for section, keys in self._fields.items():
            for key, le in keys.items():
                cp[section][key] = le.text()
        with open(self.cfgpath, 'w') as f:
            cp.write(f)
        self.configSaved.emit()
