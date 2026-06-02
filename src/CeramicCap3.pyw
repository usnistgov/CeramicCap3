import os, sys, pathlib, threading
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analysis'))
import traceback
import ctypes
import datetime
import time
import mplwidget
import numpy as np
import R2FMath
import serial
import serial.tools.list_ports
import TZA
import CCConfig
from Meas3 import Meas, reset_instruments, open_instruments, close_instruments
from Tabwidget import MyTabWidget
import CustomData
import CircuitSetup
import ConfigEditor
import analysismodule

from PyQt5.QtCore import (
    QMutex,
    QThread,
    Qt,
    pyqtSignal)

from PyQt5.QtGui import QKeySequence

from PyQt5.QtWidgets import (
    QApplication,
    QInputDialog,
    QLabel,
    QMainWindow,
    QPushButton,
    QShortcut,
    QVBoxLayout,
    QHBoxLayout,
    QWidget,
    QStatusBar,
    QProgressBar,
    QTextEdit,
)


_logdir = ''

_FSAMP_RATES = [1562.5, 3125.0, 6250.0, 12500.0, 25000.0, 50000.0,
                100000.0, 200000.0, 400000.0, 800000.0]

def pick_fsamp(fsig, fsamp_max=800000.0):
    target = 8.0 * fsig
    for r in _FSAMP_RATES:
        if r > target and r <= fsamp_max:
            return r
    return fsamp_max


class MainWindow(QMainWindow):
    stopDVM = pyqtSignal()
    resetDone = pyqtSignal()

    def __init__(self, mutex):
        super().__init__()
        self.mutex = mutex
        self.t0 = time.time()
        self.Npts = 8
        self.ver = 'unknown'
        self.quit = False
        self.measuring = False
        self.loopfinished = False
        self.yyyymmdir = ''
        self.mytext = []
        self.mytextmaxlen = 1000
        self.tzaport = []
        for p   in serial.tools.list_ports.comports():
            if p.description.startswith('TZA/OPM500'):
                self.tzaport.append(p.device)
        self.tza2 = TZA.TZA(self.tzaport[0])
        self.config = CCConfig.CCC()
        self.thread = QThread()
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        self.V1 = 1.811+2.578j   # small modulated center (SOUR1)
        self.V2 = -9.9+0j        # large constant (SOUR2)
        self.dV = 0.1
        self.g2 = 100
        self.tza2.set_fgain(self.g2)
        self.warmup_count = 0
        self._rm = self._sg1 = self._dvm = None   # opened in startMeas, closed after sweep

        self.le_fixed_g2 = None
        self.C1_eff = None
        self.V1_learned = {}
        self.Z1_learned = {}
        self._freq_converged = False
        self._started_from_V1_learned = False
        self.fsigold = -1
        self.rData = CustomData.NPoints(1000, 800000)
        self.rSet = CustomData.ThreeChannels(1000, 800000, 2, [], [], [], 0, 0, 0, ts=-1)
        self.livePhasors = [[] for _ in range(4)]
        self.allData = CustomData.AllData()
        self.parseconfig()
        self.myprint(f"Welcome to CeramicCap v{self.ver}")

        self.scatterplots = np.empty((2, 2), dtype=object)
        for i in range(2):
            for j in range(2):
                self.scatterplots[i, j] = mplwidget.MplWidget(rightax=False)
                self.scatterplots[i, j].setfmt("%.5f", "%.6f")
        self.rawplots = np.empty((2, 2), dtype=object)
        for i in range(2):
            for j in range(2):
                self.rawplots[i, j] = mplwidget.MplWidget(rightax=False)

        self.residplots = np.empty((2, 2), dtype=object)
        for i in range(2):
            for j in range(2):
                self.residplots[i, j] = mplwidget.MplWidget(rightax=False)

        self.psaplots = np.empty((2, 2), dtype=object)
        for i in range(2):
            for j in range(2):
                self.psaplots[i, j] = mplwidget.MplWidget(rightax=False)

        self.alphafplots = np.empty((2, 2), dtype=object)
        for i in range(2):
            for j in range(2):
                self.alphafplots[i, j] = mplwidget.MplWidget(rightax=(i == 1))

        self.etaplots = np.empty((1, 2), dtype=object)
        for j in range(2):
            self.etaplots[0, j] = mplwidget.MplWidget(rightax=False)

        self.balanceplots = np.empty((2, 2), dtype=object)
        for i in range(2):
            for j in range(2):
                self.balanceplots[i, j] = mplwidget.MplWidget(rightax=False)
                self.balanceplots[i, j].setfmt('%.4f', '%.5f')
        self.meta_V1rb = []
        self.meta_eta3 = []

        self.output = QTextEdit(self)
        self.output.resize(540, 200)
        self.output.setReadOnly(True)

        self.mstatus = QTextEdit(self)
        self.mstatus.resize(540, 200)
        self.mstatus.setReadOnly(True)

        self.buStart = QPushButton("Start")
        self.buStop  = QPushButton("Stop")
        self.buQuit  = QPushButton("Quit")
        button_row = QHBoxLayout()
        button_row.addWidget(self.buStart)
        button_row.addWidget(self.buStop)
        button_row.addWidget(self.buQuit)
        button_row.addStretch()

        self.circuit_setup = CircuitSetup.CircuitSetupWidget(self.config.cfgpath)
        self.config_editor = ConfigEditor.ConfigEditor(self.config.cfgpath)
        self.config_editor.configSaved.connect(self.parseconfig)
        self.tabWidget = MyTabWidget(self)
        self.psa_show_resid = False
        self.psa_resid_btn.clicked.connect(self._toggle_psa_mode)
        self.scatter_show_eta = True
        self.scatter_mode_btn.clicked.connect(self._toggle_scatter_mode)

        main_layout = QVBoxLayout()
        main_layout.addLayout(button_row)
        main_layout.addWidget(self.tabWidget)
        central_widget.setLayout(main_layout)

        self.progressBar = QProgressBar()
        self.progressBar.setMaximumWidth(400)
        self.progressBar.setRange(0, self.Npts)
        self.progressBar.setValue(0)
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        widget = QWidget(self)
        widget.setLayout(QHBoxLayout())
        self.sblabel = QLabel("")
        self.etalabel = QLabel("")
        self.freqlabel = QLabel("—")
        self.freqprogresslabel = QLabel("")
        self.sweepProgressBar = QProgressBar()
        self.sweepProgressBar.setMaximumWidth(200)
        self.sweepProgressBar.setRange(0, 1)
        self.sweepProgressBar.setValue(0)
        self.sweepProgressBar.setFormat('%v / %m freq')
        self.gainlabel = QLabel("g2=—")
        widget.layout().addWidget(self.sblabel)
        widget.layout().addWidget(self.progressBar)
        widget.layout().addWidget(self.freqlabel)
        widget.layout().addWidget(self.sweepProgressBar)
        widget.layout().addWidget(self.freqprogresslabel)
        widget.layout().addWidget(self.etalabel)
        widget.layout().addWidget(self.gainlabel)
        self.statusBar.showMessage('Welcome to CeramicCap', 5000)
        self.statusBar.addPermanentWidget(widget)

        self.buStart.clicked.connect(self.startMeas)
        self.buStop.clicked.connect(self.stopMeas)
        self.buStop.setEnabled(False)
        self.buQuit.clicked.connect(self.finishUp)
        self.resetDone.connect(self._on_reset_done)
        self.buStart.setEnabled(False)

        QShortcut(QKeySequence('Space'), self, self._toggle_meas)
        QShortcut(QKeySequence('Q'),     self, self.finishUp)

        self.setWindowTitle('CeramicCap3')
        self.statusBar.showMessage('Initializing instruments...', 0)
        threading.Thread(target=self._reset_bg, daemon=True).start()

    def _reset_bg(self):
        reset_instruments()
        self.resetDone.emit()

    def _on_reset_done(self):
        self.buStart.setEnabled(True)
        self.statusBar.showMessage('Ready — press Start to begin  [Space] start/stop  [Q] quit', 0)
        self.myprint('Instruments reset and ready')

    def _toggle_meas(self):
        if self.measuring:
            self.stopMeas()
        elif self.buStart.isEnabled():
            self.startMeas()

    def parseconfig(self):
        self.config.read()
        old_C1  = getattr(self, 'C1',  None)
        old_C2  = getattr(self, 'C2',  None)
        old_SN1 = getattr(self, 'SN1', None)
        old_SN2 = getattr(self, 'SN2', None)
        self.C1 = self.config.C1
        self.C2 = self.config.C2
        self.SN1 = self.config.SN1
        self.SN2 = self.config.SN2
        if self.C1 != old_C1 or self.C2 != old_C2 or self.SN1 != old_SN1 or self.SN2 != old_SN2:
            self.C1_eff = None
            self.V1_learned = {}
            self.Z1_learned = {}

        self.flist = self.config.flist
        if not self.measuring:
            self.fsig = self.config.fstart
        self.dvfrac = self.config.dvfrac
        self.nrMeas  = self.config.nrMeas
        self.fsamp    = self.config.fsamp
        self.forcefmax = self.config.forcefmax
        self.datadir = self.config.datadir
        self.yyyymmdir   = self.config.logdir
        self.rawdatadir  = self.config.rawdatadir
        self.saverawdata = self.config.saverawdata
        self.eta3_limit = self.config.eta3_limit
        self.max_step_frac = self.config.max_step_frac
        self.decay = self.config.decay
        self.Npts = self.config.nellipse
        if hasattr(self, 'progressBar'):
            self.progressBar.setRange(0, self.Npts)
        self.fixg2 = self.config.fixg2
        self._SAT_THRESHOLD = self.config.sat_threshold
        self.fixed_g2 = R2FMath.mingainvalue2(self.flist, self.C1) if self.fixg2 else None
        if self.C1 >= 1e-5:   # 10 µF cap: TZA gain must be 1
            self.fixg2 = True
            self.fixed_g2 = 1
        if self.le_fixed_g2 is not None:
            self.le_fixed_g2.setText(str(int(self.fixed_g2)) if self.fixed_g2 is not None else '— (fixg2 off)')
        self.ver = self.config.version
        global _logdir
        _logdir = self.config.logdir

    def RBupdate(self):
        self.gainlabel.setText(f"g2={self.g2}")

    def laUpdate(self):
        pass

    def setf(self, f):
        self.statusBar.showMessage('Next frequency: {0} kHz'.format(f/1000), 5000)
        self.fsig = f
        self.freqlabel.setText(f'{f/1000:.3g} kHz')
        self.allData.deletekey(f)

    def fp(self):
        f = self.fsig
        if f in self.flist:
            ix = self.flist.index(f)
            if ix + 1 < len(self.flist):
                self.setf(self.flist[ix + 1])

    def fm(self):
        f = self.fsig
        if f in self.flist:
            ix = self.flist.index(f)
            if ix-1 >= 0:
                self.setf(self.flist[ix-1])
            else:
                self.setf(self.flist[-1])


    def stopMeas(self):
        if not self.measuring:
            return
        self.measuring = False
        self.myprint("Stop requested — waiting for current point to finish")
        self.statusBar.showMessage('Stopping...', 0)
        self.stopDVM.emit()

    def finishUp(self):
        if not self.measuring:
            self.circuit_setup.save_config()
            self.close()
            return
        self.statusBar.showMessage('Wait for graceful exit')
        self.myprint("Stopping — waiting for current point to finish")
        self.measuring = False
        self.quit = True
        self.stopDVM.emit()

    def closeEvent(self, event):
        if self.measuring and not self.loopfinished:
            self.statusBar.showMessage('Wait for graceful exit')
            self.myprint("Close event called, but loop not finished")
            if not self.quit:
                self.stopDVM.emit()
            self.quit = True
            event.ignore()
        else:
            self.circuit_setup.save_config()
            self.statusBar.showMessage('Quitting')
            self.myprint("Quitting gracefully")
            event.accept()

    def _reduce_gain(self, g):
        """Reduce gain by one decade, floored at 1. Returns new gain."""
        if g <= 1:
            return 1
        new_g = max(1, g // 10)
        return new_g

    def onNewData(self, MyData: CustomData.NPoints):
        old_rData = self.rData
        self.livePhasors = [[] for _ in range(4)]
        self.rData = MyData
        self.rData.Res['switch_pos'] = 0 if self.switch_normal else 1

        # Gain adjustment is only permitted during warmup
        max_v4 = self.rData.max_raw_amplitude(2)
        in_warmup = not self._freq_converged
        sat_threshold = min(self._SAT_THRESHOLD, 0.9 * self.v4_range)
        sat4 = in_warmup and (not self.fixg2) and max_v4 > sat_threshold
        if sat4:
            new_g2 = self._reduce_gain(self.g2)
            if new_g2 < self.g2:
                self._g2_ceiling = self.g2  # never boost back to this gain or above
                self.g2 = new_g2
                self.tza2.set_fgain(self.g2)
                self.myprint(f'r gain {self.g2}')
                self.statusBar.showMessage(f'Saturation: V4 -> g2={self.g2} — gain reduced', 4000)
                self.RBupdate()
                self.warmup_count = 0
                self.meta_V1rb = []
                self.meta_eta3 = []
                old_rData.strip_raw()
                return  # discard saturated measurement; next goagain() uses the reduced gain
            # Already at minimum gain — try increasing the DAQ range before proceeding
            _V4_RANGES = [3, 18]
            try:
                vi = _V4_RANGES.index(self.v4_range)
                if vi < len(_V4_RANGES) - 1:
                    self.v4_range = _V4_RANGES[vi + 1]
                    self.myprint(f'V4 saturated at g2=1 — DAQ range -> {self.v4_range} V')
                    self.warmup_count = 0
                    self.meta_V1rb = []
                    self.meta_eta3 = []
                    old_rData.strip_raw()
                    return
            except ValueError:
                pass
            self.myprint(f'V4 saturated at minimum gain and maximum DAQ range — proceeding')

        # Low-signal check: if V4 is tiny during warmup, boost gain one decade
        # Do not boost to a gain that previously caused saturation at this frequency
        low4 = (in_warmup and (not self.fixg2) and max_v4 < self._SAT_THRESHOLD / 20
                and self.g2 < 100000 and self.g2 * 10 < getattr(self, '_g2_ceiling', float('inf')))
        if low4:
            self.g2 = min(100000, self.g2 * 10)
            self.tza2.set_fgain(self.g2)
            self.myprint(f'i gain {self.g2}')
            self.statusBar.showMessage(f'Low signal: V4 -> g2={self.g2} — gain increased', 4000)
            self.RBupdate()
            self.warmup_count = 0
            self.meta_V1rb = []
            self.meta_eta3 = []
            old_rData.strip_raw()
            return  # discard low-signal measurement; next goagain() uses the boosted gain

        self.V1 = self.rData.Res['V1_balance']
        m   = self.rData.Res['eta3fit_slope']
        v4m = self.rData.Res['eta3_mean']
        self.meta_V1rb.append(self.rData.Res['V1rb_center'])
        self.meta_eta3.append(self.rData.Res['eta3_mean'])
        V1rb_best = self.calc_V1rb_best()
        _v1bo_clamped = True   # conservative default: don't update C1_eff unless converged
        if V1rb_best is not None:
            if np.isfinite(V1rb_best.real) and np.isfinite(V1rb_best.imag):
                V1_center_approx = complex(np.mean(self.meta_V1rb))
                step_mag = abs(V1rb_best - V1_center_approx)
                max_step_approx = self.max_step_frac * abs(V1_center_approx) + 0.05
                _v1bo_clamped = step_mag >= 0.9 * max_step_approx
                clamp_flag = ' [clamped]' if _v1bo_clamped else ''
                V1bo_str = f'  V1bo={V1rb_best:.5f}{clamp_flag}'
                self.V1 = V1rb_best
            else:
                V1bo_str = f'  V1bo=non-finite({V1rb_best})'
        else:
            V1bo_str = ''
        _converged = abs(v4m) < self.eta3_limit
        self.myprint(f"fit m={m:.2f}  c={v4m:.4f}  |c|={abs(v4m):.4f}/{self.eta3_limit:.4f}  V1b={self.rData.Res['V1_balance']:.5f}{V1bo_str}")
        if self.rData.goodData:
            self.warmup_count += 1
            if _converged:
                self._freq_converged = True
                self.allData.append(self.rData)
                self.V1_learned[self.fsig] = self.V1
                if not _v1bo_clamped:
                    self._update_C1_eff()
                c1eff_str = f'{self.C1_eff:.4e}' if self.C1_eff is not None else 'None'
                clamp_str = ' [clamped—not updated]' if _v1bo_clamped else ''
                self.myprint(f'  C1_eff={c1eff_str}{clamp_str}')
            elif self.warmup_count == 1 and self._started_from_V1_learned and abs(v4m) > self.eta3_limit * 20:
                # V1_learned is stale (first measurement already far from balance); discard and restart
                del self.V1_learned[self.fsig]
                if self.fsig in self.Z1_learned:
                    del self.Z1_learned[self.fsig]
                self._started_from_V1_learned = False
                self.V1 = self.calcVsmall(self.fsig, V2=self.V2)
                self.meta_V1rb = []
                self.meta_eta3 = []
                self.myprint(f'  V1_learned stale (|c|={abs(v4m):.3f})—reset to calcVsmall {self.V1:.5f}')
        else:
            self.myprint('  goodData=False — non-finite combined3, point discarded')
        cur_sw = 0 if self.switch_normal else 1
        sw_str = 'straight' if self.switch_normal else 'cross'
        if not _converged:
            self.sblabel.setText(f"warmup |c|={abs(v4m):.4f} (lim {self.eta3_limit:.4f}) [{sw_str}] at {self.fsig/1000:5.2f} kHz")
        else:
            self.sblabel.setText("{0}/{1} [{2}] at {3:5.2f} kHz".format(
                self.allData.countf_sw(self.fsig, cur_sw), self.nrMeas, sw_str, self.fsig/1000))
        if self.fsig == self.fsigold:
            if self.allData.countf_sw(self.fsig, cur_sw) == self.nrMeas:
                if self.switch_normal:      # sw0 complete — switch to sw1
                    self.switch_normal = False
                    self.warmup_count = 0
                    self.meta_V1rb = []
                    self.meta_eta3 = []
                    self.myprint(f"straight done at {self.fsig/1000:.4g} kHz — switching to cross")
                else:                       # sw1 complete — save and advance
                    self.saveData(self.fsig)
                    self.switch_normal = True
                    self.freqs_done += 1
                    self.freqprogresslabel.setText(f'{self.freqs_done}/{self.total_freqs} freq')
                    self.sweepProgressBar.setValue(self.freqs_done)
                    elapsed = time.time() - self.run_start_time
                    secs_remaining = (self.total_freqs - self.freqs_done) * elapsed / self.freqs_done
                    eta = datetime.datetime.now() + datetime.timedelta(seconds=secs_remaining)
                    self.etalabel.setText(f"ETA {eta.strftime('%H:%M')}")
                    if self.fsig == self.flist[-1]:
                        self.measuring = False
                    else:
                        self.fp()
        self.replot()
        old_rData.strip_raw()

    def calcVsmall(self, f, V2=-9.0, R=50):
        C2 = self.C2
        iw = 1j * f * 2 * np.pi
        I1 = V2 * iw * C2 / (1 + iw * C2 * R)
        if f in self.Z1_learned:
            Z_C1 = self.Z1_learned[f]
        elif self.C1_eff is not None:
            Z_C1 = 1.0 / (iw * self.C1_eff)
        else:
            Z_C1 = 1.0 / (iw * self.C1)
        V1 = -I1 * (R + Z_C1)
        return V1

    def _update_C1_eff(self, R=50):
        """Infer effective C1 impedance from converged V1 and V2; store real part as C1_eff and full complex Z per frequency."""
        try:
            iw = 1j * 2 * np.pi * self.fsig
            I2 = self.V2 * iw * self.C2 / (1 + iw * self.C2 * R)
            Z_cap = -self.V1 / I2 - R
            C1_new = float(np.real(1.0 / (iw * Z_cap)))
            if C1_new > 0 and 0.1 * self.C1 < C1_new < 10 * self.C1:
                self.C1_eff = C1_new
                self.Z1_learned[self.fsig] = Z_cap
        except Exception:
            pass

    def startMeas(self):
        desc, ok = QInputDialog.getText(self, 'Run description', 'Enter a one-line description for this run:')
        if not ok:
            return
        self.run_description = desc.strip()
        self.config_editor.save()       # triggers configSaved → parseconfig
        self.circuit_setup.save_config()
        self.measuring = True
        desc = self.run_description
        if desc:
            self.myprint('=' * 56)
            self.myprint(desc)
            self.myprint('=' * 56)
        self.myprint("Measuring C1={0:3.0e} C2={1:3.0e}".format(self.C1, self.C2))
        self.myprint("S/N: {0} {1}".format(self.SN1, self.SN2))
        self.buStart.setEnabled(False)
        self.buStop.setEnabled(True)
        self.buStop.setStyleSheet('background-color: #c0392b; color: white; font-weight: bold;')
        self.allData = CustomData.AllData()
        self.runDataDir = self.ensureDir()
        self.meta_V1rb = []
        self.meta_eta3 = []
        self.C1_eff = None
        self.fsig = next((f for f in self.flist if f >= self.config.fstart), self.flist[-1])
        self.fsigold = -1
        self.warmup_count = 0
        self.v4_range = 3
        self.switch_normal = True
        self._run_logged = False
        self.progressBar.setValue(0)
        self.run_start_time = time.time()
        self.freqs_done = 0
        try:
            start_idx = self.flist.index(self.fsig)
        except ValueError:
            start_idx = 0
        self.total_freqs = len(self.flist) - start_idx
        self.freqlabel.setText(f'{self.fsig/1000:.3g} kHz')
        self.freqprogresslabel.setText(f'0/{self.total_freqs} freq')
        self.sweepProgressBar.setRange(0, self.total_freqs)
        self.sweepProgressBar.setValue(0)
        self.etalabel.setText('')
        self.setWindowTitle(f'CeramicCap3 — {self.run_description}')
        self.myprint("Starting frequency sweep")
        self.statusBar.showMessage('Measuring...', 0)
        self._rm, self._sg1, self._dvm = open_instruments()
        self.goagain()

    def goagain(self):
        if self.quit:
            if self._rm is not None:
                close_instruments(self._rm, self._sg1, self._dvm)
                self._rm = self._sg1 = self._dvm = None
            self.loopfinished = True
            self.close()
        elif self.measuring:
            self.thread = QThread()
            cur_fsamp = self.fsamp if self.forcefmax else pick_fsamp(self.fsig, self.fsamp)
            cur_nsamp = min(int(cur_fsamp), 160000)  # cap at 160k samples to limit transfer time
            fit_cache = self.rData._fit_cache if self.fsigold == self.fsig else {}
            self.mydvm = Meas(self.mutex, self.Npts, self.rawdatadir, self.saverawdata, cur_fsamp, cur_nsamp, chunk_periods=10, max_nhars=self.config.max_nhars, v4_range=self.v4_range, switch_normal=self.switch_normal, fit_cache=fit_cache, sg1=self._sg1, dvm=self._dvm, rm=self._rm)
            self.mydvm.moveToThread(self.thread)
            self.tza2.set_fgain(self.g2)
            if self.fsigold == self.fsig:
                while np.abs(self.V1) > 9:
                    self.V1 = self.V1*0.9
                    self.V2 = self.V2*0.9
                self.mydvm.storeV(self.V1, self.V2, self.dV, self.fsig, self.g2)
            else:
                self.myprint(f"  {self.fsig/1000:.4g} kHz  fsamp={cur_fsamp/1000:.4g} kHz  ".center(64, '-'))
                tempgain2 = R2FMath.newgainvalue2(self.fsig, self.C1, dV= self.dV,Vmax=3)
                self.g2 = self.fixed_g2 if self.fixg2 else tempgain2
                C1_used = self.C1_eff if self.C1_eff is not None else self.C1
                self.myprint(f'getgain C1_used={C1_used:.3e} (nom {self.C1:.1e}) tempgain2={tempgain2} g2={self.g2}')
                self.config.setGains(self.g2)
                self.tza2.set_fgain(self.g2)
                self.RBupdate()
                self.warmup_count = 0
                self._freq_converged = False
                self._g2_ceiling = float('inf')  # reset per-frequency gain ceiling
                self.sblabel.setText(f"warming up at {self.fsig/1000:5.2f} kHz")
                self.meta_V1rb = []
                self.meta_eta3 = []
                self.V2 = -9.9+0j
                if self.fsig in self.V1_learned:
                    self.V1 = self.V1_learned[self.fsig]
                    self._started_from_V1_learned = True
                    self.myprint(f'  V1=learned {self.V1:.5f}  V2={self.V2:.5f}')
                else:
                    self._started_from_V1_learned = False
                    self.V1 = self.calcVsmall(self.fsig, V2=self.V2)
                    src = 'Z1learned' if self.fsig in self.Z1_learned else 'calcVsmall'
                    self.myprint(f'  V1={src} {self.V1:.5f}  V2={self.V2:.5f}')
                while np.abs(self.V1) > 9:
                    self.V2 *= 0.9
                    self.V1 *= 0.9
                self.dV = np.abs(self.V1) * self.dvfrac
                self.mydvm.storeV(self.V1, self.V2, self.dV, self.fsig, self.g2)
            self.laUpdate()
            self.mydvm.logMessage.connect(self.myprint)
            self.mydvm.dataReady.connect(self.onNewData)
            self.mydvm.dataSetReady.connect(self.onNewSet)
            self.mydvm.finished.connect(self.thread.quit)
            self.mydvm.finished.connect(self.mydvm.deleteLater)
            self.thread.finished.connect(self.goagain)
            self.stopDVM.connect(self.mydvm.stop, Qt.DirectConnection)
            self.thread.started.connect(self.mydvm.start)
            self.thread.finished.connect(self.thread.deleteLater)
            self.thread.start()
            self.fsigold = self.fsig
        else:
            if self._rm is not None:
                close_instruments(self._rm, self._sg1, self._dvm)
                self._rm = self._sg1 = self._dvm = None
            stopped = not self.measuring
            self.myprint("Stopped" if stopped else "Frequency sweep complete")
            self.statusBar.showMessage('Stopped — press Start to measure again' if stopped
                                       else 'Sweep complete — press Start to measure again', 0)
            self.buStart.setEnabled(True)
            self.buStop.setEnabled(False)
            self.buStop.setStyleSheet('')
            self.etalabel.setText(f"Done {datetime.datetime.now().strftime('%H:%M:%S')}")
            title_suffix = 'Stopped' if stopped else 'Complete'
            self.setWindowTitle(f'CeramicCap3 — {getattr(self, "run_description", "")} [{title_suffix}]')

    def onNewSet(self, MySet: CustomData.ThreeChannels):
        self.rSet = MySet
        self.progressBar.setValue(MySet.i+1)
        if MySet.ts > 0:
            for j in range(3):
                self.livePhasors[j].append(MySet.Data[j].Vc)
        currentTab = self.tabWidget.master.tabText(self.tabWidget.master.currentIndex())
        if currentTab == 'scatter':
            self.plotscatter()
        elif currentTab in ('raw', 'PSA'):
            self.replot()

    def myprint(self, newtext):
        n = datetime.datetime.now()
        t = time.time()-self.t0
        t1 = n.strftime('%H:%M:%S')+' >{0:8.1f} < : '.format(t)+newtext
        t2 = n.strftime('%m/%d')+' '+t1+'\n'
        fn = 'CAP'+n.strftime('%m%d%Y')+'.LOG'
        if self.yyyymmdir != '':
            logdir = os.path.join(self.yyyymmdir, n.strftime('%Y%m'))
            pathlib.Path(logdir).mkdir(parents=True, exist_ok=True)
            with open(os.path.join(logdir, fn), "a") as file:
                file.write(t2)
        self.mytext.append(t1)
        while len(self.mytext) > self.mytextmaxlen:
            self.mytext.pop(0)

    def replot(self):
        currentIndex = self.tabWidget.master.currentIndex()
        tat = self.tabWidget.master.tabText(currentIndex)
        if getattr(self, '_prev_tab', None) == 'config' and tat != 'config':
            self.config_editor.save()
        self._prev_tab = tat
        if tat == 'scatter':
            self.plotscatter()
        elif tat == 'msg':
            self.showOutput()
        elif tat == 'raw':
            self.plotraw()
        elif tat == 'resid':
            self.plotresid()
        elif tat == 'PSA':
            self.plotpsa()
        elif tat == 'alpha(f)':
            self.plotalphaf()
        elif tat == 'eta':
            self.ploteta()
        elif tat == 'V1bal':
            self.plotbalance()
        elif tat == 'last status':
            self.showLastStatus()
        elif tat == 'config':
            self.config_editor.reload()

    def ensureDir(self):
        yymm = datetime.datetime.now().strftime('%y%m')
        dd   = datetime.datetime.now().strftime('%d')
        bd0 = os.path.join(self.datadir, self.SN1, yymm,dd)
        pathlib.Path(bd0).mkdir(parents=True, exist_ok=True)
        return bd0

    def saveData(self, f):
        bd0 = self.runDataDir
        dt = datetime.datetime.fromtimestamp(self.run_start_time)
        valname = self.config.recapdir[self.C1] + '-' + self.config.recapdir[self.C2]
        tag = dt.strftime('%Y%m%d_%H%M')

        if not self._run_logged:
            log_line = '\t'.join([tag, self.config.recapdir[self.C1],
                                   self.config.recapdir[self.C2],
                                   self.SN1, self.SN2, self.run_description]) + '\n'
            with open(os.path.join(self.config.datadir, 'runlog.txt'), 'a', encoding='utf-8') as lf:
                lf.write(log_line)
            self._run_logged = True

        fn_conf = f'conf_{valname}_{tag}.ini'
        if not os.path.exists(os.path.join(bd0, fn_conf)):
            import shutil
            shutil.copy2(self.config.cfgpath, os.path.join(bd0, fn_conf))

        # VOLT6_: raw phasors, one row per measurement point
        df_volt = self.allData.getRawPhasors(f, self.run_start_time)
        volt_path = os.path.join(bd0, f'VOLT6_{valname}_{tag}.dat')
        if not os.path.exists(volt_path):
            with open(volt_path, 'w') as fh:
                fh.write(f'# {self.run_description}\n')
                fh.write('# ' + ' '.join(df_volt.columns) + '\n')
        df_volt.to_csv(volt_path, mode='a', sep=' ', index=False, header=False, float_format='%.8g')

        # CC6_: one row per ellipse sweep; expand complex columns to re_/im_ pairs
        df_cc_raw = self.allData.getEllipseData(f)
        df_cc = pd.DataFrame()
        for col in df_cc_raw.columns:
            if pd.api.types.is_complex_dtype(df_cc_raw[col]):
                df_cc[f're_{col}'] = df_cc_raw[col].values.real
                df_cc[f'im_{col}'] = df_cc_raw[col].values.imag
            else:
                df_cc[col] = df_cc_raw[col]
        cc_path = os.path.join(bd0, f'CC6_{valname}_{tag}.dat')
        if not os.path.exists(cc_path):
            with open(cc_path, 'w') as fh:
                fh.write(f'# {self.run_description}\n')
                fh.write('# ' + ' '.join(df_cc.columns) + '\n')
        df_cc.to_csv(cc_path, mode='a', sep=' ', index=False, header=False, float_format='%.8g')

    def plotscatter(self):
        for j in range(2):
            for i in range(2):
                self.scatterplots[i, j].canvas.ax1.cla()
        show_eta = getattr(self, 'scatter_show_eta', True)
        if self.livePhasors[0]:
            V1 = np.array(self.livePhasors[0])
            V2 = np.array(self.livePhasors[1])
            V3 = np.array(self.livePhasors[2])
            p2 = V2/V1 if show_eta else V2
            p3 = V3/V1 if show_eta else V3
            for sl, mk in [(slice(None, None, 2), '+'), (slice(1, None, 2), 'x')]:
                self.scatterplots[0, 0].canvas.ax1.plot(np.real(V1[sl]), np.imag(V1[sl]), 'r' + mk)
                self.scatterplots[0, 1].canvas.ax1.plot(np.real(p2[sl]), np.imag(p2[sl]), 'g' + mk)
                self.scatterplots[1, 0].canvas.ax1.plot(np.real(p3[sl]), np.imag(p3[sl]), 'b' + mk)
        elif self.rData.Res['ts'] > 0:
            V1 = self.rData.raw3[:, 0]
            V2 = self.rData.raw3[:, 1]
            V3 = self.rData.raw3[:, 2]
            p2 = V2/V1 if show_eta else V2
            p3 = V3/V1 if show_eta else V3
            for sl, mk in [(slice(None, None, 2), 'o'), (slice(1, None, 2), 's')]:
                self.scatterplots[0, 0].canvas.ax1.plot(np.real(V1[sl]), np.imag(V1[sl]), 'r' + mk)
                self.scatterplots[0, 1].canvas.ax1.plot(np.real(p2[sl]), np.imag(p2[sl]), 'g' + mk)
                self.scatterplots[1, 0].canvas.ax1.plot(np.real(p3[sl]), np.imag(p3[sl]), 'b' + mk)

        ax00 = self.scatterplots[0, 0].canvas.ax1
        ax01 = self.scatterplots[0, 1].canvas.ax1
        ax10 = self.scatterplots[1, 0].canvas.ax1
        ax00.set_xlabel('Re(V1) / V'); ax00.set_ylabel('Im(V1) / V')
        if show_eta:
            ax01.set_xlabel('Re(eta2)'); ax01.set_ylabel('Im(eta2)')
            ax10.set_xlabel('Re(eta3)'); ax10.set_ylabel('Im(eta3)')
        else:
            ax01.set_xlabel('Re(V2) / V'); ax01.set_ylabel('Im(V2) / V')
            ax10.set_xlabel('Re(V3) / V'); ax10.set_ylabel('Im(V3) / V')

        for j in range(2):
            for i in range(2):
                self.scatterplots[i, j].canvas.draw()

    def plotalphaf(self):
        if self.allData.count() == 0:
            return
        for i in range(2):
            for j in range(2):
                ax = self.alphafplots[i, j].canvas.ax1
                ax.set_xscale('linear')
                ax.cla()

        freqs = self.allData.freqs()
        ana = {uf: [(analysismodule.analyze_block(nd.V1m, nd.V2m, nd.V3m, nd.V3m),
                     nd.Res.get('switch_pos', 0))
                    for nd in self.allData.entries(uf)]
               for uf in freqs}

        # straight (sw_pos=0): open circle;  cross (sw_pos=1): filled square
        SW_MARKERS = {0: ('o', 'none'), 1: ('s', None)}

        def add_scalar(ax, key, color):
            for sw, (mk, mfc) in SW_MARKERS.items():
                f_pts, v_pts, f_bar, v_bar, e_bar = [], [], [], [], []
                mkw = {'ms': 4}
                if mfc is not None:
                    mkw['mfc'] = mfc
                for uf in freqs:
                    vals = np.array([r[key] for r, sp in ana[uf] if sp == sw], dtype=float)
                    if len(vals) == 0:
                        continue
                    if len(vals) > 2:
                        f_bar.append(uf)
                        v_bar.append(np.mean(vals))
                        e_bar.append(np.std(vals, ddof=1))
                    else:
                        f_pts.extend([uf] * len(vals))
                        v_pts.extend(vals.tolist())
                if f_pts:
                    ax.plot(f_pts, v_pts, color + mk, **mkw)
                if f_bar:
                    ax.errorbar(f_bar, v_bar, yerr=e_bar, fmt=color + mk, capsize=3, **mkw)

        # Clear bx1 on the bottom row (rightax=True widgets)
        for j in range(2):
            self.alphafplots[1, j].canvas.bx1.cla()

        # Top-left: α₃
        ax_al = self.alphafplots[0, 0].canvas.ax1
        add_scalar(ax_al, 'al_right', 'r')
        ax_al.set_ylabel('α₃  (Re ΔC/C)')

        # Top-right: D₃
        ax_d = self.alphafplots[0, 1].canvas.ax1
        add_scalar(ax_d, 'D_right', 'b')
        ax_d.set_ylabel('D₃  (Im ΔC/C)')

        # Bottom row: |Y₄₂Z₄| (left) and ∠Y₄₂Z₄ (right)
        ax_mag = self.alphafplots[1, 0].canvas.ax1
        ax_ang = self.alphafplots[1, 1].canvas.ax1
        for sw, (mk, mfc) in SW_MARKERS.items():
            f_pts, mag_pts, ang_pts = [], [], []
            f_bar, mag_bar, mag_std, ang_bar, ang_std = [], [], [], [], []
            mkw = {'ms': 4}
            if mfc is not None:
                mkw['mfc'] = mfc
            for uf in freqs:
                g = np.array([1.0 / r['g_right'] for r, sp in ana[uf] if sp == sw], dtype=complex)
                if len(g) == 0:
                    continue
                mags = np.abs(g)
                angs = np.angle(g, deg=True)
                if len(g) > 2:
                    f_bar.append(uf)
                    mag_bar.append(np.mean(mags))
                    mag_std.append(np.std(mags, ddof=1))
                    ang_bar.append(np.mean(angs))
                    ang_std.append(np.std(angs, ddof=1))
                else:
                    f_pts.extend([uf] * len(g))
                    mag_pts.extend(mags.tolist())
                    ang_pts.extend(angs.tolist())
            if f_pts:
                ax_mag.plot(f_pts, mag_pts, 'r' + mk, **mkw)
                ax_ang.plot(f_pts, ang_pts, 'b' + mk, **mkw)
            if f_bar:
                ax_mag.errorbar(f_bar, mag_bar, yerr=mag_std, fmt='r' + mk, capsize=3, **mkw)
                ax_ang.errorbar(f_bar, ang_bar, yerr=ang_std, fmt='b' + mk, capsize=3, **mkw)
        ax_mag.set_yscale('log')
        ax_mag.set_ylabel('|Y₄₂Z₄|')
        ax_ang.set_ylabel('∠Y₄₂Z₄ (°)')

        for i in range(2):
            for j in range(2):
                self.alphafplots[i, j].canvas.ax1.set_xscale('log')
                self.alphafplots[i, j].canvas.draw()

    def ploteta(self):
        if self.rData.Res['ts'] <= 0 or not hasattr(self.rData, 'combined3'):
            return
        for j in range(2):
            self.etaplots[0, j].canvas.ax1.cla()
        mul = 1e6
        offset = self.rData.Res['ratio']

        def split(v):
            return (np.real(v) - offset) * mul, np.imag(v) * mul

        c4x, c4y   = split(self.rData.combined3)
        Vz3x, Vz3y = split(self.rData.Res['Vz3'])
        ax4 = self.etaplots[0, 1].canvas.ax1
        ax4.plot(c4x, c4y, 'm+')
        ax4.plot(Vz3x, Vz3y, 'r*', markersize=10)
        ax4.set_xlabel('(Re(γ₃η₃ − η₂) − ratio) × 10⁶')
        ax4.set_ylabel('Im(γ₃η₃ − η₂) × 10⁶')
        for j in range(2):
            self.etaplots[0, j].canvas.draw()

    def calc_V1rb_best(self, decay=None, max_step_frac=None):
        if decay is None:
            decay = self.decay
        if max_step_frac is None:
            max_step_frac = self.max_step_frac
        if len(self.meta_V1rb) < 3:
            return None
        V1rb = np.array(self.meta_V1rb)
        eta3 = np.array(self.meta_eta3)
        N = len(V1rb)
        # Exponential weights: most recent point = 1, older points decay geometrically.
        # decay=0.85 means a point 10 measurements ago has weight ~0.20.
        weights = decay ** np.arange(N - 1, -1, -1)
        w = np.sqrt(weights)   # scale rows so lstsq minimises weighted sum of squares
        xr, xi = V1rb.real, V1rb.imag
        er, ei = eta3.real, eta3.imag
        xr_mean = np.dot(weights, xr) / weights.sum()
        xi_mean = np.dot(weights, xi) / weights.sum()
        V1_center = xr_mean + 1j * xi_mean
        dxr, dxi = xr - xr_mean, xi - xi_mean
        A = np.block([
            [np.column_stack([dxr, -dxi, np.ones(N), np.zeros(N)])],
            [np.column_stack([dxi,  dxr, np.zeros(N), np.ones(N)])]
        ])
        W = np.concatenate([w, w])
        params, _, _, _ = np.linalg.lstsq(A * W[:, np.newaxis], np.concatenate([er, ei]) * W, rcond=None)
        a, b, cr, ci = params
        m_meta = a + 1j * b
        c_meta = cr + 1j * ci
        if abs(m_meta) < 1e-12:
            return None
        result = V1_center - c_meta / m_meta
        if not np.isfinite(result.real) or not np.isfinite(result.imag) or abs(result) > 10.0:
            return None
        # Clamp step size: never move more than max_step_frac of |V1_center| from cluster
        # centre in one iteration. Prevents runaway when the fit extrapolates far outside
        # the measured V1rb cluster (small slope, poorly conditioned).
        step = result - V1_center
        max_step = max_step_frac * abs(V1_center) + 0.05
        if abs(step) > max_step:
            result = V1_center + step * (max_step / abs(step))
        return result

    def plotbalance(self):
        for i in range(2):
            for j in range(2):
                self.balanceplots[i, j].canvas.ax1.cla()
        if not self.meta_V1rb:
            for i in range(2):
                for j in range(2):
                    self.balanceplots[i, j].canvas.draw()
            return
        V1rb = np.array(self.meta_V1rb[-6:])
        eta3 = np.array(self.meta_eta3[-6:])
        xr, xi = V1rb.real, V1rb.imag
        er, ei = eta3.real, eta3.imag
        axs = [[self.balanceplots[i, j].canvas.ax1 for j in range(2)] for i in range(2)]
        axs[0][0].set_xlabel('Re(V1rb)')
        axs[0][0].set_ylabel('Re(η₃)')
        axs[0][1].set_xlabel('Im(V1rb)')
        axs[0][1].set_ylabel('Re(η₃)')
        axs[1][0].set_xlabel('Re(V1rb)')
        axs[1][0].set_ylabel('Im(η₃)')
        axs[1][1].set_xlabel('Im(V1rb)')
        axs[1][1].set_ylabel('Im(η₃)')
        import matplotlib
        n_shades = self.nrMeas + 5
        cmap = matplotlib.colormaps['Blues']
        pt_colors = [cmap(0.25 + 0.75 * k / (n_shades - 1)) for k in range(n_shades)]
        for k, (xrk, xik, erk, eik) in enumerate(zip(xr, xi, er, ei)):
            c = pt_colors[min(k, n_shades - 1)]
            axs[0][0].plot(xrk, erk, '.', color=c, markersize=8)
            axs[0][1].plot(xik, erk, '.', color=c, markersize=8)
            axs[1][0].plot(xrk, eik, '.', color=c, markersize=8)
            axs[1][1].plot(xik, eik, '.', color=c, markersize=8)
        for row in axs:
            for ax in row:
                ax.tick_params(axis='x', rotation=30)
        V1rb_best = self.calc_V1rb_best()
        if V1rb_best is not None:
            xr_mean, xi_mean = np.mean(xr), np.mean(xi)
            dxr, dxi = xr - xr_mean, xi - xi_mean
            N = len(xr)
            A = np.block([
                [np.column_stack([dxr, -dxi, np.ones(N), np.zeros(N)])],
                [np.column_stack([dxi,  dxr, np.zeros(N), np.ones(N)])]
            ])
            params, _, _, _ = np.linalg.lstsq(A, np.concatenate([er, ei]), rcond=None)
            a, b, cr, ci = params
            xr_best, xi_best = V1rb_best.real, V1rb_best.imag
            er_fit = a * (xr - xr_mean) - b * (xi - xi_mean) + cr
            ei_fit = b * (xr - xr_mean) + a * (xi - xi_mean) + ci
            idx_r = np.argsort(xr)
            idx_i = np.argsort(xi)
            axs[0][0].plot(xr[idx_r], er_fit[idx_r], 'b-')
            axs[0][1].plot(xi[idx_i], er_fit[idx_i], 'b-')
            axs[1][0].plot(xr[idx_r], ei_fit[idx_r], 'b-')
            axs[1][1].plot(xi[idx_i], ei_fit[idx_i], 'b-')
            for ax in [axs[0][0], axs[1][0]]:
                ax.axvline(xr_best, color='r', linestyle='--')
            for ax in [axs[0][1], axs[1][1]]:
                ax.axvline(xi_best, color='r', linestyle='--')
            axs[0][0].plot(xr_best, 0, 'r*', markersize=10)
            axs[0][1].plot(xi_best, 0, 'r*', markersize=10)
            axs[1][0].plot(xr_best, 0, 'r*', markersize=10)
            axs[1][1].plot(xi_best, 0, 'r*', markersize=10)
        for i in range(2):
            for j in range(2):
                self.balanceplots[i, j].canvas.draw()

    def _toggle_scatter_mode(self):
        self.scatter_show_eta = not self.scatter_show_eta
        self.scatter_mode_btn.setText('Show: eta' if self.scatter_show_eta else 'Show: raw')
        self.plotscatter()

    def _toggle_psa_mode(self):
        self.psa_show_resid = not self.psa_show_resid
        self.psa_resid_btn.setText('Show: Raw' if self.psa_show_resid else 'Show: Residuals')
        self.plotpsa()

    def plotpsa(self):
        if self.rSet.ts <= 0:
            return
        fsamp = self.rSet.fsamp
        fsig  = self.rSet.fsig
        colors = ['r', 'g', 'b', 'm']
        labels = ['Ch1', 'Ch2', 'Ch3', 'Ch4']
        for idx, (row, col) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
            ax = self.psaplots[row, col].canvas.ax1
            ax.cla()
            if idx >= len(self.rSet.Data):
                self.psaplots[row, col].canvas.draw()
                continue
            d = self.rSet.Data[idx]
            raw = np.asarray(d.data, dtype=float)
            if self.psa_show_resid and d.fv is not None:
                signal = raw - np.asarray(d.fv, dtype=float)
            else:
                signal = raw
            n = len(signal)
            w = np.hanning(n)
            amp = np.abs(np.fft.rfft(signal * w)) / (np.sum(w) / 2)
            freqs = np.fft.rfftfreq(n, d=1.0 / fsamp)
            ax.semilogy(freqs[1:], amp[1:], colors[idx] + '-', linewidth=0.5)
            # Arrow from local floor up to the signal-frequency bin
            sig_bin = np.argmin(np.abs(freqs[1:] - fsig))
            n_adj = 8
            lo = max(0, sig_bin - n_adj)
            hi = min(len(freqs) - 2, sig_bin + n_adj)
            adj = np.concatenate([amp[1 + lo : 1 + sig_bin],
                                   amp[2 + sig_bin : 2 + hi]])
            arrow_top = adj.min() if len(adj) > 0 else amp[1:].min()
            ax.annotate('', xy=(fsig, arrow_top), xytext=(fsig, amp[1:].min()),
                        arrowprops=dict(arrowstyle='->', color='k', lw=1.5,
                                        mutation_scale=14))
            ax.set_xscale('log')
            ax.set_xlabel('f (Hz)')
            suffix = ' (resid)' if self.psa_show_resid else ''
            ax.set_ylabel(labels[idx] + suffix)
        for row in range(2):
            for col in range(2):
                self.psaplots[row, col].canvas.draw()

    def plotraw(self):
        if self.rSet.ts <= 0:
            return
        for j in range(2):
            for i in range(2):
                self.rawplots[i, j].canvas.ax1.cla()
        t = np.arange(len(self.rSet.Data[0].data))
        n = len(t)
        Vc = self.rSet.Data[0].Vc
        phi = -np.angle(-1j*Vc) % (2*np.pi)
        period = self.rSet.fsamp / self.rSet.fsig
        be0 = phi / (2*np.pi) * period
        ncycles = round((n/2 - be0) / period)
        be = int(be0 + ncycles * period)
        en = be + int(period * 2)
        self.rawplots[0, 0].canvas.ax1.plot(t[be:en], self.rSet.Data[0].data[be:en], 'r.')
        self.rawplots[0, 1].canvas.ax1.plot(t[be:en], self.rSet.Data[1].data[be:en], 'g.')
        self.rawplots[1, 0].canvas.ax1.plot(t[be:en], self.rSet.Data[2].data[be:en], 'b.')
        self.rawplots[0, 0].canvas.ax1.plot(t[be:en], self.rSet.Data[0].fv[be:en], 'k-')
        self.rawplots[0, 1].canvas.ax1.plot(t[be:en], self.rSet.Data[1].fv[be:en], 'k-')
        self.rawplots[1, 0].canvas.ax1.plot(t[be:en], self.rSet.Data[2].fv[be:en], 'k-')
        for j in range(2):
            for i in range(2):
                self.rawplots[i, j].canvas.draw()

    def plotresid(self):
        if self.rSet.ts <= 0:
            return
        colors = ['r', 'g', 'b', 'm']
        labels = ['Ch1 (V1)', 'Ch2 (V2)', 'Ch3 (V3)', 'Ch4 (V4)']
        plots  = [(0, 0), (0, 1), (1, 0), (1, 1)]
        for idx, (row, col) in enumerate(plots):
            ax = self.residplots[row, col].canvas.ax1
            ax.cla()
            if idx < len(self.rSet.Data):
                d = self.rSet.Data[idx]
                if d.data is not None and d.fv is not None:
                    resid = d.data - d.fv
                    ax.plot(resid, colors[idx] + ',')
                    ax.axhline(0, color='k', linewidth=0.5)
                    ax.set_ylabel(f'residual {labels[idx]}')
                    ax.set_xlabel('sample')
            self.residplots[row, col].canvas.draw()

    def showOutput(self):
        self.output.setText('\n'.join(self.mytext))
        scrollbar = self.output.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

    def showLastStatus(self):
        if self.rData.Res['ts'] > 0:
            o = ''
            for k, v in self.rData.Res.items():
                if isinstance(v, complex):
                    if v.imag > 0:
                        o += '{0:30}:{1:8.3f} + i {2:8.3f}\n'.format(k, v.real, v.imag)
                    else:
                        o += '{0:30}:{1:8.3f} - i {2:8.3f}\n'.format(k, v.real, -v.imag)
                elif isinstance(v, float):
                    o += '{0:30}:{1:8.3f}\n'.format(k, v)
                elif isinstance(v, int):
                    o += '{0:30}:{1}\n'.format(k, v)
            self.mstatus.setText(o)
            scrollbar = self.mstatus.verticalScrollBar()
            scrollbar.setValue(scrollbar.maximum())


def excepthook(exc_type, exc_value, exc_tb):
    tb = "".join(traceback.format_exception(exc_type, exc_value, exc_tb))
    now = datetime.datetime.now()
    print("error catched at ", now)
    print("error message:\n", tb)
    errlog = os.path.join(_logdir or os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'R2F_errlog.dat')
    with open(errlog, 'w') as fi:
        fi.write("error catched!:\n")
        fi.write("error message:\n" + tb)
    app.quit()


if __name__ == '__main__':
    mutex = QMutex()
    myappid = 'mycompany.myproduct.subproduct.version'
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    app = QApplication(sys.argv)
    sys.excepthook = excepthook
    window = MainWindow(mutex)
    window.show()
    sys.exit(app.exec())
