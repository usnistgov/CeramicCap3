import os, sys, pathlib, threading
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analysis'))
import analysismodule
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
from Meas3 import Meas, reset_instruments
from Tabwidget import MyTabWidget
import CustomData
import CircuitSetup
import ConfigEditor

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
        self.tza1 = TZA.TZA(self.tzaport[1])
        self.tza2 = TZA.TZA(self.tzaport[0])
        self.config = CCConfig.CCC()
        self.thread = QThread()
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        self.V1 = 1.811+2.578j   # small modulated center (SOUR1)
        self.V2 = -9.9+0j        # large constant (SOUR2)
        self.dV = 0.1
        self.g1 = 10000
        self.g2 = 100
        self.tza1.set_fgain(self.g1)
        self.tza2.set_fgain(self.g2)
        self.warmup_count = 0

        self.le_fixed_g1 = None
        self.le_fixed_g2 = None
        self.C41_eff = None
        self.fsigold = -1
        self.rData = CustomData.NPoints(1000, 800000)
        self.rSet = CustomData.FourChannels(1000, 800000, 2, [], [], [], [], 0, 0, 0, ts=-1)
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
        self.meta_eta4 = []

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

        main_layout = QVBoxLayout()
        main_layout.addLayout(button_row)
        main_layout.addWidget(self.tabWidget)
        central_widget.setLayout(main_layout)

        self.progressBar = QProgressBar()
        self.progressBar.setMaximumWidth(400)
        self.progressBar.setRange(0, self.Npts*2)
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
        self.gainlabel = QLabel("g1=—  g2=—")
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
        self.C31 = self.config.C31
        self.C32 = self.config.C32
        old_C41  = getattr(self, 'C41',  None)
        old_C42  = getattr(self, 'C42',  None)
        old_SN41 = getattr(self, 'SN41', None)
        old_SN42 = getattr(self, 'SN42', None)
        self.C41 = self.config.C41
        self.C42 = self.config.C42
        desc = getattr(self, 'run_description', '').strip()
        if desc:
            self.myprint('=' * 56)
            self.myprint(desc)
            self.myprint('=' * 56)
        self.myprint("Measuring C31={0:3.0e} C32= {1:3.0e} & C41={2:3.0e} C42={3:3.0e}".format(self.C31, self.C32, self.C41, self.C42))

        self.SN31 = self.config.SN31
        self.SN32 = self.config.SN32
        self.SN41 = self.config.SN41
        self.SN42 = self.config.SN42
        if self.C41 != old_C41 or self.C42 != old_C42 or self.SN41 != old_SN41 or self.SN42 != old_SN42:
            self.C41_eff = None

        self.myprint("S/N: {0} {1} {2} {3}".format(self.SN31, self.SN32, self.SN41, self.SN42))

        self.flist = self.config.flist
        if not self.measuring:
            self.fsig = self.config.fstart
        self.dvfrac = self.config.dvfrac
        self.nrMeas  = self.config.nrMeas
        self.fsamp   = self.config.fsamp
        self.nsamp   = self.config.nsamp
        self.datadir = self.config.datadir
        self.yyyymmdir   = self.config.logdir
        self.rawdatadir  = self.config.rawdatadir
        self.saverawdata = self.config.saverawdata
        self.nwarmup = self.config.nwarmup
        self.fixg1 = self.config.fixg1
        self.fixg2 = self.config.fixg2
        self._SAT_THRESHOLD = self.config.sat_threshold
        self.switching = self.config.switching
        self.fixed_g1 = R2FMath.mingainvalue1(self.flist, self.C31, self.dV) if self.fixg1 else None
        self.fixed_g2 = R2FMath.mingainvalue2(self.flist, self.C41) if self.fixg2 else None
        if self.le_fixed_g1 is not None:
            self.le_fixed_g1.setText(str(int(self.fixed_g1)) if self.fixed_g1 is not None else '— (fixg1 off)')
            self.le_fixed_g2.setText(str(int(self.fixed_g2)) if self.fixed_g2 is not None else '— (fixg2 off)')
        self.ver = self.config.version
        global _logdir
        _logdir = self.config.logdir

    def RBupdate(self):
        self.gainlabel.setText(f"g1={self.g1}  g2={self.g2}")

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

        # Gain adjustment is only permitted during warmup
        max_v3 = self.rData.max_raw_amplitude(2)
        max_v4 = self.rData.max_raw_amplitude(3)
        in_warmup = self.warmup_count <= self.nwarmup
        sat3 = in_warmup and (not self.fixg1) and max_v3 > self._SAT_THRESHOLD
        sat4 = in_warmup and (not self.fixg2) and max_v4 > self._SAT_THRESHOLD
        if sat3 or sat4:
            myop=f'r gain {sat3=} {sat4=} {self.g1} {self.g2}'
            if sat3:
                self.g1 = self._reduce_gain(self.g1)
                self.tza1.set_fgain(self.g1)
            if sat4:
                self.g2 = self._reduce_gain(self.g2)
                self.tza2.set_fgain(self.g2)
            myop+=f'to {self.g1} {self.g2}'
            self.myprint(myop)
            parts = []
            if sat3: parts.append(f'V3 → g1={self.g1}')
            if sat4: parts.append(f'V4 → g2={self.g2}')
            self.statusBar.showMessage(f'⚠  Saturation: {", ".join(parts)} — gain reduced', 4000)
            self.RBupdate()
            self.warmup_count = 0
            self.meta_V1rb = []
            self.meta_eta4 = []
            old_rData.strip_raw()
            return  # discard saturated measurement; next goagain() uses the reduced gain

        # Low-signal check: if V3 or V4 is tiny during warmup, boost gain one decade
        low3 = in_warmup and (not self.fixg1) and max_v3 < self._SAT_THRESHOLD / 20 and self.g1 < 100000
        low4 = in_warmup and (not self.fixg2) and max_v4 < self._SAT_THRESHOLD / 20 and self.g2 < 100000
        if low3 or low4:
            myop = f'i gain {low3=} {low4=} {self.g1} {self.g2}'
            if low3:
                self.g1 = min(100000, self.g1 * 10)
                self.tza1.set_fgain(self.g1)
            if low4:
                self.g2 = min(100000, self.g2 * 10)
                self.tza2.set_fgain(self.g2)
            myop += f' to {self.g1} {self.g2}'
            self.myprint(myop)
            parts = []
            if low3: parts.append(f'V3 → g1={self.g1}')
            if low4: parts.append(f'V4 → g2={self.g2}')
            self.statusBar.showMessage(f'↑  Low signal: {", ".join(parts)} — gain increased', 4000)
            self.RBupdate()
            self.warmup_count = 0
            self.meta_V1rb = []
            self.meta_eta4 = []
            old_rData.strip_raw()
            return  # discard low-signal measurement; next goagain() uses the boosted gain

        self.V1 = self.rData.Res['V1_balance']
        m   = self.rData.Res['eta4fit_slope']
        v4m = self.rData.Res['eta4_mean']
        self.meta_V1rb.append(self.rData.Res['V1rb_center'])
        self.meta_eta4.append(self.rData.Res['eta4_mean'])
        V1rb_best = self.calc_V1rb_best()
        if V1rb_best is not None:
            if np.isfinite(V1rb_best.real) and np.isfinite(V1rb_best.imag):
                V1bo_str = f'  V1bo={V1rb_best:.5f}'
                self.V1 = V1rb_best
            else:
                V1bo_str = f'  V1bo=non-finite({V1rb_best})'
        else:
            V1bo_str = ''
        self.myprint(f"fit m={m:.2f}  c={v4m:.4f}  V1b={self.rData.Res['V1_balance']:.5f}{V1bo_str}")
        self._update_C41_eff()
        if self.rData.goodData:
            self.warmup_count += 1
            if self.warmup_count > self.nwarmup:
                self.allData.append(self.rData)
        if self.warmup_count <= self.nwarmup:
            self.sblabel.setText(f"warmup {self.warmup_count}/{self.nwarmup} at {self.fsig/1000:5.2f} kHz")
        else:
            self.sblabel.setText("{0}/{1} good measurements at {2:5.2f} kHz".format(
                self.allData.countf(self.fsig), self.nrMeas, self.fsig/1000))
        if self.fsig == self.fsigold:
            if self.allData.countf(self.fsig) == self.nrMeas:
                self.saveData(self.fsig)
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
        C42 = self.C42
        C41 = self.C41_eff if self.C41_eff is not None else self.C41
        iw = 1j * f * 2 * np.pi
        I1 = V2 * iw * C42 / (1 + iw * C42 * R)
        V1 = -I1 * (R + 1 / (iw * C41))
        return V1

    def _update_C41_eff(self, R=50):
        """Infer effective C41 from the current measured V1_balance and V2."""
        try:
            iw = 1j * 2 * np.pi * self.fsig
            I2 = self.V2 * iw * self.C42 / (1 + iw * self.C42 * R)
            Z_cap = -self.V1 / I2 - R
            C41_new = float(np.real(1.0 / (iw * Z_cap)))
            if C41_new > 0 and 0.1 * self.C41 < C41_new < 10 * self.C41:
                self.C41_eff = C41_new
        except Exception:
            pass

    def startMeas(self):
        desc, ok = QInputDialog.getText(self, 'Run description', 'Enter a one-line description for this run:')
        if not ok:
            return
        self.run_description = desc.strip()
        self.config_editor.save()
        self.circuit_setup.save_config()
        self.parseconfig()
        self.measuring = True
        self.buStart.setEnabled(False)
        self.buStop.setEnabled(True)
        self.buStop.setStyleSheet('background-color: #c0392b; color: white; font-weight: bold;')
        self.allData = CustomData.AllData()
        self.runDataDir = self.ensureDir()
        self.meta_V1rb = []
        self.meta_eta4 = []
        self.C41_eff = None
        self.fsig = next((f for f in self.flist if f >= self.config.fstart), self.flist[-1])
        self.fsigold = -1
        self.warmup_count = 0
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
        self.goagain()

    def goagain(self):
        if self.quit:
            self.loopfinished = True
            self.close()
        elif self.measuring:
            self.thread = QThread()
            self.mydvm = Meas(self.mutex, self.Npts, self.rawdatadir, self.saverawdata, self.fsamp, self.nsamp, switching=self.switching, chunk_periods=10, max_nhars=self.config.max_nhars)
            self.mydvm.moveToThread(self.thread)
            self.tza1.set_fgain(self.g1)
            self.tza2.set_fgain(self.g2)
            if self.fsigold == self.fsig:
                while np.abs(self.V1) > 9:
                    self.V1 = self.V1*0.9
                    self.V2 = self.V2*0.9
                self.mydvm.storeV(self.V1, self.V2, self.dV, self.fsig, self.g1, self.g2)
            else:
                self.myprint(f"  {self.fsig/1000:.4g} kHz  ".center(56, '-'))
                tempgain1 = R2FMath.newgainvalue1(self.fsig, self.C31, dV= self.dV,Vmax=3)
                tempgain2 = R2FMath.newgainvalue2(self.fsig, self.C41, dV= self.dV,Vmax=3)

                self.g1 = self.fixed_g1 if self.fixg1 else tempgain1
                self.g2 = self.fixed_g2 if self.fixg2 else tempgain2
                #self.g1 = 1
                #self.g2 = 1
                C41_used = self.C41_eff if self.C41_eff is not None else self.C41
                self.myprint(f'getgain {self.C31=:.1e} C41_used={C41_used:.3e} (nom {self.C41:.1e}) {tempgain1=} {tempgain2=} {self.g1} {self.g2}')
                self.config.setGains(self.g1, self.g2)
                self.tza1.set_fgain(self.g1)
                self.tza2.set_fgain(self.g2)
                self.RBupdate()
                self.warmup_count = 0
                self.sblabel.setText(f"warming up at {self.fsig/1000:5.2f} kHz")
                self.meta_V1rb = []
                self.meta_eta4 = []
                V1_balance_prev = self.V1   # measured balance from previous frequency
                V2_prev        = self.V2   # V2 that was used at previous frequency
                self.V2 = -9.9+0j
                V1_model_next = self.calcVsmall(self.fsig, V2=self.V2)
                if self.fsigold > 0:
                    # Apply a multiplicative correction derived from the previous
                    # frequency: ratio = measured / model(f_prev).  Multiplying
                    # model(f_next) by this ratio corrects for parasitics without
                    # assuming the balance is identical across a large frequency jump.
                    V1_model_prev = self.calcVsmall(self.fsigold, V2=V2_prev)
                    if abs(V1_model_prev) > 1e-10:
                        ratio = V1_balance_prev / V1_model_prev
                        self.V1 = V1_model_next * ratio
                    else:
                        self.V1 = V1_model_next
                else:
                    self.V1 = V1_model_next
                # Sanity-check V1 against the plain nominal model (C41, not C41_eff).
                # If ratio correction or a corrupted C41_eff produced a wildly wrong
                # V1, fall back to the nominal so V2 is never scaled to zero.
                iw_nom = 1j * self.fsig * 2 * np.pi
                I1_nom = self.V2 * iw_nom * self.C42 / (1 + iw_nom * self.C42 * 50)
                V1_nominal = -I1_nom * (50 + 1 / (iw_nom * self.C41))
                if abs(V1_nominal) > 1e-10 and abs(self.V1) > 3 * abs(V1_nominal):
                    self.myprint(f'V1={self.V1:.4f} too far from nominal {V1_nominal:.4f} — using nominal')
                    self.V1 = V1_nominal
                while np.abs(self.V1) > 9:
                    self.V2 *= 0.9
                    self.V1 *= 0.9
                self.dV = np.abs(self.V1) * self.dvfrac
                self.myprint(f'{self.V1=:.5f} {self.V2=:.5f} {self.dV=:.5f}')
                self.mydvm.storeV(self.V1, self.V2, self.dV, self.fsig, self.g1, self.g2)
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

    def onNewSet(self, MySet: CustomData.FourChannels):
        self.rSet = MySet
        self.progressBar.setValue(MySet.i+1)
        if MySet.ts > 0:
            for j in range(4):
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
        bd0 = os.path.join(self.datadir, self.SN41, yymm,dd)
        pathlib.Path(bd0).mkdir(parents=True, exist_ok=True)
        return bd0

    def saveData(self, f):
        bd0 = self.runDataDir
        dt = datetime.datetime.fromtimestamp(self.run_start_time)
        valname = self.config.recapdir[self.C41]+'-'+self.config.recapdir[self.C42]
        fn_conf = 'conf_'+valname+'_'+dt.strftime('%Y%m%d_%H%M')+'.ini'
        if not os.path.exists(os.path.join(bd0, fn_conf)):
            import shutil
            shutil.copy2(self.config.cfgpath, os.path.join(bd0, fn_conf))
        fn = 'CC3_'+valname+'_'+dt.strftime('%Y%m%d_%H%M')+'.dat'
        mykeys = ['fsig', 'ts', 'alpha3mean', 'beta3mean', 'alpha4mean', 'beta4mean',
                  'Vz3', 'Vz4', 'gamma3', 'gamma4', 'V1cReadback']
        rdict = self.allData.getkeys(f, mykeys)
        if not os.path.exists(os.path.join(bd0, fn)):
            with open(os.path.join(bd0, fn), "w") as file:
                file.write('# {}\n'.format(self.run_description))
                o = '# '
                for k, v in rdict.items():
                    if isinstance(v[0], complex):
                        o += 'Re({0}) Im({0}) '.format(k)
                    else:
                        o += '{0} '.format(k)
                o += '\n'
                file.write(o)
        with open(os.path.join(bd0, fn), "a") as file:
            k = list(rdict)[0]
            L = len(rdict[k])
            for n in range(L):
                o = ''
                for k, v in rdict.items():
                    a = v[n]
                    if isinstance(a, complex):
                        o += '{0:12.9f} {1:12.9f} '.format(a.real, a.imag)
                    else:
                        o += '{0:12.9f} '.format(a)
                o += '\n'
                file.write(o)
        C = self.allData.getRawPhasors(f, self.run_start_time)
        fn = 'VOLT_'+valname+'_'+dt.strftime('%Y%m%d_%H%M')+'.dat'
        if not os.path.exists(os.path.join(bd0, fn)):
            with open(os.path.join(bd0, fn), "w") as file:
                file.write('# {}\n'.format(self.run_description))
                file.write('# frequency/Hz t/s  reV1 imV1 reV2 imV2 reV3 imV3 reV4 imV4 reV1set imV1set reV2set imV2set gain1 gain2\n')
        with open(os.path.join(bd0, fn), "a") as file:
            for n in range(np.shape(C)[0]):
                o = '{0:6.0f} {1:6.1f} '.format(C[n, 0], C[n, 1])
                for k in range(2, np.shape(C)[1]):
                    o += '{0:15.8f} '.format(C[n, k])
                o += '\n'
                file.write(o)

    def plotscatter(self):
        for j in range(2):
            for i in range(2):
                self.scatterplots[i, j].canvas.ax1.cla()
        if self.livePhasors[0]:
            V1 = np.array(self.livePhasors[0])
            V2 = np.array(self.livePhasors[1])
            V3 = np.array(self.livePhasors[2])
            V4 = np.array(self.livePhasors[3])
            eta2 = V2/V1
            eta3 = V3/V1
            eta4 = V4/V1
            self.scatterplots[0, 0].canvas.ax1.plot(np.real(V1),  np.imag(V1), 'r+')
            self.scatterplots[0, 1].canvas.ax1.plot(np.real(eta2), np.imag(eta2), 'g+')
            self.scatterplots[1, 0].canvas.ax1.plot(np.real(eta3), np.imag(eta3), 'b+')
            self.scatterplots[1, 1].canvas.ax1.plot(np.real(eta4), np.imag(eta4), 'm+')
        elif self.rData.Res['ts'] > 0:
            V1 = self.rData.ave4[:, 0]
            V2 = self.rData.ave4[:, 1]
            V3 = self.rData.ave4[:, 2]
            V4 = self.rData.ave4[:, 3]
            eta2 = V2/V1
            eta3 = V3/V1
            eta4 = V4/V1
                #ell = R2FMath.ComplexEllipse.fit_from_cmplx_points(ratio)
            self.scatterplots[0, 0].canvas.ax1.plot(np.real(V1),  np.imag(V1), 'ro')
            self.scatterplots[0, 1].canvas.ax1.plot(np.real(eta2), np.imag(eta2), 'go')
            self.scatterplots[1, 0].canvas.ax1.plot(np.real(eta3), np.imag(eta3), 'bo')
            self.scatterplots[1, 1].canvas.ax1.plot(np.real(eta4), np.imag(eta4), 'mo')

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

        freqs = sorted(self.allData.mydict.keys())
        ana = {uf: [analysismodule.analyze_block(nd.V1m, nd.V2m, nd.V3m, nd.V4m)
                    for nd in self.allData.mydict[uf]]
               for uf in freqs}

        def add_scalar(ax, key, color):
            f_pts, v_pts, f_bar, v_bar, e_bar = [], [], [], [], []
            for uf in freqs:
                vals = np.array([r[key] for r in ana[uf]], dtype=float)
                if len(vals) > 2:
                    f_bar.append(uf)
                    v_bar.append(np.mean(vals))
                    e_bar.append(np.std(vals, ddof=1))
                else:
                    f_pts.extend([uf] * len(vals))
                    v_pts.extend(vals.tolist())
            if f_pts:
                ax.plot(f_pts, v_pts, color + 'o', ms=4)
            if f_bar:
                ax.errorbar(f_bar, v_bar, yerr=e_bar, fmt=color + 'o', capsize=3, ms=4)

        from matplotlib.lines import Line2D as _L2D
        for ax, a_key, d_key, side in [
            (self.alphafplots[0, 0].canvas.ax1, 'al_left',  'D_left',  '₃'),
            (self.alphafplots[0, 1].canvas.ax1, 'al_right', 'D_right', '₄'),
        ]:
            add_scalar(ax, a_key, 'r')
            add_scalar(ax, d_key, 'b')
            ax.set_ylabel(f'α{side} / D{side}')
            ax.legend(handles=[
                _L2D([0], [0], color='r', marker='o', ls='', label=f'α{side}  (Re ΔC/C)'),
                _L2D([0], [0], color='b', marker='o', ls='', label=f'D{side}  (Im ΔC/C)'),
            ], fontsize=8)

        from matplotlib.lines import Line2D
        for ax, bx, g_key, label in [
            (self.alphafplots[1, 0].canvas.ax1, self.alphafplots[1, 0].canvas.bx1, 'g_left',  'Y₃₂Z₃'),
            (self.alphafplots[1, 1].canvas.ax1, self.alphafplots[1, 1].canvas.bx1, 'g_right', 'Y₄₂Z₄'),
        ]:
            ax.set_xscale('linear')
            ax.cla()
            bx.cla()
            f_pts, mag_pts, ang_pts = [], [], []
            f_bar, mag_bar, mag_std, ang_bar, ang_std = [], [], [], [], []
            for uf in freqs:
                g = np.array([1.0 / r[g_key] for r in ana[uf]], dtype=complex)
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
                ax.plot(f_pts, mag_pts, 'ro', ms=4)
                bx.plot(f_pts, ang_pts, 'bo', ms=4)
            if f_bar:
                ax.errorbar(f_bar, mag_bar, yerr=mag_std, fmt='ro', capsize=3, ms=4)
                bx.errorbar(f_bar, ang_bar, yerr=ang_std, fmt='bo', capsize=3, ms=4)
            ax.set_yscale('log')
            ax.set_ylabel('|' + label + '|')
            bx.set_ylabel('∠' + label + ' (°)')
            ax.legend(handles=[
                Line2D([0], [0], color='r', marker='o', ls='', label='|' + label + '|'),
                Line2D([0], [0], color='b', marker='o', ls='', label='∠' + label),
            ], fontsize=8)

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

        c3x, c3y   = split(self.rData.combined3)
        c4x, c4y   = split(self.rData.combined4)
        Vz3x, Vz3y = split(self.rData.Res['Vz3'])
        Vz4x, Vz4y = split(self.rData.Res['Vz4'])
        ax3 = self.etaplots[0, 0].canvas.ax1
        ax4 = self.etaplots[0, 1].canvas.ax1
        ax3.plot(c3x, c3y, 'b+')
        ax3.plot(Vz3x, Vz3y, 'r*', markersize=10)
        ax3.set_xlabel('(Re(γ₃η₃ − η₂) − ratio) × 10⁶')
        ax3.set_ylabel('Im(γ₃η₃ − η₂) × 10⁶')
        ax4.plot(c4x, c4y, 'm+')
        ax4.plot(Vz4x, Vz4y, 'r*', markersize=10)
        ax4.set_xlabel('(Re(γ₄η₄ − η₂) − ratio) × 10⁶')
        ax4.set_ylabel('Im(γ₄η₄ − η₂) × 10⁶')
        for j in range(2):
            self.etaplots[0, j].canvas.draw()

    def calc_V1rb_best(self, decay=0.85):
        if len(self.meta_V1rb) < 3:
            return None
        V1rb = np.array(self.meta_V1rb)
        eta4 = np.array(self.meta_eta4)
        N = len(V1rb)
        # Exponential weights: most recent point = 1, older points decay geometrically.
        # decay=0.85 means a point 10 measurements ago has weight ~0.20.
        weights = decay ** np.arange(N - 1, -1, -1)
        w = np.sqrt(weights)   # scale rows so lstsq minimises weighted sum of squares
        xr, xi = V1rb.real, V1rb.imag
        er, ei = eta4.real, eta4.imag
        xr_mean = np.dot(weights, xr) / weights.sum()
        xi_mean = np.dot(weights, xi) / weights.sum()
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
        result = (xr_mean + 1j * xi_mean) - c_meta / m_meta
        if not np.isfinite(result.real) or not np.isfinite(result.imag) or abs(result) > 10.0:
            return None
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
        eta4 = np.array(self.meta_eta4[-6:])
        xr, xi = V1rb.real, V1rb.imag
        er, ei = eta4.real, eta4.imag
        axs = [[self.balanceplots[i, j].canvas.ax1 for j in range(2)] for i in range(2)]
        axs[0][0].set_xlabel('Re(V1rb)')
        axs[0][0].set_ylabel('Re(η₄)')
        axs[0][1].set_xlabel('Im(V1rb)')
        axs[0][1].set_ylabel('Re(η₄)')
        axs[1][0].set_xlabel('Re(V1rb)')
        axs[1][0].set_ylabel('Im(η₄)')
        axs[1][1].set_xlabel('Im(V1rb)')
        axs[1][1].set_ylabel('Im(η₄)')
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
        self.rawplots[1, 1].canvas.ax1.plot(t[be:en], self.rSet.Data[3].data[be:en], 'm.')
        self.rawplots[0, 0].canvas.ax1.plot(t[be:en], self.rSet.Data[0].fv[be:en], 'k-')
        self.rawplots[0, 1].canvas.ax1.plot(t[be:en], self.rSet.Data[1].fv[be:en], 'k-')
        self.rawplots[1, 0].canvas.ax1.plot(t[be:en], self.rSet.Data[2].fv[be:en], 'k-')
        self.rawplots[1, 1].canvas.ax1.plot(t[be:en], self.rSet.Data[3].fv[be:en], 'k-')
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
