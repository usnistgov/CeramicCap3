import os, sys, pathlib, threading
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

from PyQt5.QtCore import (
    QMutex,
    QThread,
    Qt,
    pyqtSignal)

from PyQt5.QtWidgets import (
    QApplication,
    QLabel,
    QMainWindow,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QWidget,
    QSpacerItem,
    QSizePolicy,
    QStatusBar,
    QProgressBar,
    QTextEdit,
    QRadioButton,
    QButtonGroup
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
        self.ver = 3.0
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
        self.firstgood = False

        self.fsigold = -1
        self.rData = CustomData.NPoints(1000, 800000)
        self.rSet = CustomData.FourChannels(1000, 800000, 2, [], [], [], [], 0, 0, 0, ts=-1)
        self.livePhasors = [[] for _ in range(4)]
        self.allData = CustomData.AllData()
        self.myprint(f"Welcome to Version {self.ver}")
        self.parseconfig()

        self.scatterplots = np.empty((2, 2), dtype=object)
        for i in range(2):
            for j in range(2):
                self.scatterplots[i, j] = mplwidget.MplWidget(rightax=False)
                self.scatterplots[i, j].setfmt("%.5f", "%.6f")
        self.rawplots = np.empty((2, 2), dtype=object)
        for i in range(2):
            for j in range(2):
                self.rawplots[i, j] = mplwidget.MplWidget(rightax=False)

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

        self.output = QTextEdit(self)
        self.output.resize(540, 200)
        self.output.setReadOnly(True)

        self.mstatus = QTextEdit(self)
        self.mstatus.resize(540, 200)
        self.mstatus.setReadOnly(True)

        glayout = QHBoxLayout()
        vlayout = QVBoxLayout()
        vlayout.addItem(QSpacerItem(60, 20, QSizePolicy.Fixed, QSizePolicy.Fixed))

        self.buStart = QPushButton("Start")
        self.buQuit = QPushButton("Quit")
        self.bufp = QPushButton("f++")
        self.bufm = QPushButton("f--")
        vlayout.addWidget(self.buStart)
        vlayout.addWidget(self.buQuit)
        vlayout.addWidget(self.bufp)
        vlayout.addWidget(self.bufm)
        vlayout.addWidget(QLabel('current f:'))
        self.laf = QLabel('{0:5.1f} kHz'.format(self.fsig/1000))
        vlayout.addWidget(self.laf)
        vlayout.addWidget(QLabel('Re(V1):'))
        self.lareV1 = QLabel("{0:8.3f} V".format(np.real(self.V1)))
        vlayout.addWidget(self.lareV1)
        vlayout.addWidget(QLabel('Im(V1):'))
        self.laimV1 = QLabel('{0:8.3f} V'.format(np.imag(self.V1)))
        vlayout.addWidget(self.laimV1)
        vlayout.addItem(QSpacerItem(60, 20, QSizePolicy.Fixed, QSizePolicy.Expanding))

        glayout.addLayout(vlayout)
        self.circuit_setup = CircuitSetup.CircuitSetupWidget(self.config.cfgpath)
        self.tabWidget = MyTabWidget(self)
        glayout.addWidget(self.tabWidget)
        vlayout = QVBoxLayout()
        self.rg1 = QButtonGroup()
        self.rg2 = QButtonGroup()

        self.rb1 = [QRadioButton(i) for i in self.tza1.gaindict.values()]
        vlayout.addWidget(QLabel('Gain left'))
        for rb in self.rb1:
            self.rg1.addButton(rb)
            vlayout.addWidget(rb)
            if rb.text() == self.tza1.hgain:
                rb.setChecked(True)

        self.rb2 = [QRadioButton(i) for i in self.tza2.gaindict.values()]
        vlayout.addWidget(QLabel('Gain right'))
        for rb in self.rb2:
            self.rg2.addButton(rb)
            vlayout.addWidget(rb)
            if rb.text() == self.tza2.hgain:
                rb.setChecked(True)
        vlayout.addItem(QSpacerItem(60, 20, QSizePolicy.Fixed, QSizePolicy.Expanding))
        glayout.addLayout(vlayout)
        self.rg1.buttonToggled.connect(self.rb1Toggled)
        self.rg2.buttonToggled.connect(self.rb2Toggled)

        central_widget.setLayout(glayout)

        self.progressBar = QProgressBar()
        self.progressBar.setMaximumWidth(400)
        self.progressBar.setRange(0, self.Npts*2)
        self.progressBar.setValue(0)
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        widget = QWidget(self)
        widget.setLayout(QHBoxLayout())
        self.sblabel = QLabel("")
        widget.layout().addWidget(self.sblabel)
        widget.layout().addWidget(self.progressBar)
        self.statusBar.showMessage('Welcome to CeramicCap', 5000)
        self.statusBar.addPermanentWidget(widget)

        self.buStart.clicked.connect(self.startMeas)
        self.buQuit.clicked.connect(self.finishUp)
        self.bufp.clicked.connect(self.fp)
        self.bufm.clicked.connect(self.fm)
        self.resetDone.connect(self._on_reset_done)
        self.buStart.setEnabled(False)
        self.statusBar.showMessage('Initializing instruments...', 0)
        threading.Thread(target=self._reset_bg, daemon=True).start()

    def _reset_bg(self):
        reset_instruments()
        self.resetDone.emit()

    def _on_reset_done(self):
        self.buStart.setEnabled(True)
        self.statusBar.showMessage('Ready — press Start to begin', 0)
        self.myprint('Instruments reset and ready')

    def parseconfig(self):
        self.config.read()
        self.C31 = self.config.C31
        self.C32 = self.config.C32
        self.C41 = self.config.C41
        self.C42 = self.config.C42
        self.myprint("Measuring C31={0:3.0e} C32= {1:3.0e} & C41={2:3.0e} C42={3:3.0e}".format(self.C31, self.C32, self.C41, self.C42))

        self.SN31 = self.config.SN31
        self.SN32 = self.config.SN32
        self.SN41 = self.config.SN41
        self.SN42 = self.config.SN42

        self.myprint("S/N: {0} {1} {2} {3}".format(self.SN31, self.SN32, self.SN41, self.SN42))

        self.fsig = self.config.fstart
        self.flist = self.config.flist
        self.dvfrac = self.config.dvfrac
        self.nrMeas = self.config.nrMeas
        self.datadir = self.config.datadir
        self.yyyymmdir   = self.config.logdir
        self.rawdatadir  = self.config.rawdatadir
        self.saverawdata = self.config.saverawdata
        global _logdir
        _logdir = self.config.logdir

    def RBupdate(self):
        for rb in self.rb1:
            if rb.text() == self.tza1.hgain:
                rb.setChecked(True)
        for rb in self.rb2:
            if rb.text() == self.tza2.hgain:
                rb.setChecked(True)

    def laUpdate(self):
        self.laf.setText('{0:8.1f} kHz'.format(self.fsig/1000))
        self.lareV1.setText('{0:8.3f} V'.format(np.real(self.V1)))
        self.laimV1.setText('{0:8.3f} V'.format(np.imag(self.V1)))

    def setf(self, f):
        self.statusBar.showMessage('Next frequency: {0} kHz'.format(f/1000), 5000)
        self.fsig = f
        self.allData.deletekey(f)

    def fp(self):
        f = self.fsig
        if f in self.flist:
            ix = self.flist.index(f)
            if ix+1 < len(self.flist):
                self.setf(self.flist[ix+1])
            else:
                self.setf(self.flist[0])

    def fm(self):
        f = self.fsig
        if f in self.flist:
            ix = self.flist.index(f)
            if ix-1 >= 0:
                self.setf(self.flist[ix-1])
            else:
                self.setf(self.flist[-1])

    def rb1Toggled(self, button, checked):
        if checked:
            self.tza1.set_hgain(button.text())
            self.g1 = self.tza1.text_to_gain(button.text())

    def rb2Toggled(self, button, checked):
        if checked:
            self.tza2.set_hgain(button.text())
            self.g2 = self.tza2.text_to_gain(button.text())

    def finishUp(self):
        if not self.measuring:
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
            self.statusBar.showMessage('Quitting')
            self.myprint("Quitting gracefully")
            event.accept()


    def onNewData(self, MyData: CustomData.NPoints):
        self.livePhasors = [[] for _ in range(4)]
        self.rData = MyData
        self.V1 = self.rData.Res['V1_balance']
        m = self.rData.Res['V4fit_slope']
        c = self.rData.Res['V4fit_intercept']
        self.myprint(f"V4 fit: slope={m:.5f}  intercept={c:.5f}  V1_balance={self.V1:.5f}")
        if self.rData.goodData:
            if self.firstgood:
                self.allData.append(self.rData)
            else:
                self.firstgood = True
        self.sblabel.setText("{0}/{1} good measurements at {2:5.2f} kHz".format(
            self.allData.countf(self.fsig), self.nrMeas, self.fsig/1000))
        if self.fsig == self.fsigold:
            if self.allData.countf(self.fsig) == self.nrMeas:
                self.saveData(self.fsig)
                if self.fsig == self.flist[-1]:
                    self.measuring = False
                else:
                    self.fp()
        self.replot()

    def calcVsmall(self, f, V2=-9.0, R=50):
        C42 = self.C42
        C41 = self.C41
        iw = 1j*f*2*np.pi
        I1 = V2*iw*C42/(1+iw*C42*R)
        V1 = -I1*(R+1/(iw*C41))
        return V1

    def startMeas(self):
        self.circuit_setup.save_config()
        self.parseconfig()
        self.measuring = True
        self.buStart.setEnabled(False)
        self.fsig = next((f for f in self.flist if f >= self.config.fstart), self.flist[-1])
        self.fsigold = -1
        self.firstgood = False
        self.progressBar.setValue(0)
        self.myprint("Starting frequency sweep")
        self.statusBar.showMessage('Measuring...', 0)
        self.goagain()

    def goagain(self):
        if self.quit:
            self.loopfinished = True
            self.close()
        elif self.measuring:
            self.thread = QThread()
            self.mydvm = Meas(self.mutex, self.Npts, self.rawdatadir, self.saverawdata)
            self.mydvm.moveToThread(self.thread)
            self.tza1.set_fgain(self.g1)
            self.tza2.set_fgain(self.g2)
            if self.fsigold == self.fsig:
                while np.abs(self.V1) > 9:
                    self.V1 = self.V1*0.9
                    self.V2 = self.V2*0.9
                self.mydvm.storeV(self.V1, self.V2, self.dV, self.fsig, self.g1, self.g2)
            else:
                self.g1 = R2FMath.newgainvalue1(self.fsig, self.C31, self.C41)
                #self.g1 = R2FMath.newgainvalue1(self.fsig, self.C31, 1e-6)
                self.g2 = R2FMath.newgainvalue2(self.fsig, self.C41)
                #self.g1 = 1
                #self.g2 = 1
                self.myprint('pick gains {0} {1}'.format(self.g1, self.g2))
                self.config.setGains(self.g1, self.g2)
                self.tza1.set_fgain(self.g1)
                self.tza2.set_fgain(self.g2)
                self.RBupdate()
                self.firstgood = False
                self.V2 = -9.9+0j
                self.V1 = self.calcVsmall(self.fsig, V2=self.V2)
                while np.abs(self.V1) > 9:
                    self.V2 = self.V2*0.9
                    self.V1 = self.calcVsmall(self.fsig, V2=self.V2)
                self.dV = np.abs(self.V1) * self.dvfrac
                self.myprint(f'{self.V1=} {self.V2=} {self.dV=}')
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
            self.myprint("Frequency sweep complete")
            self.statusBar.showMessage('Sweep complete — press Start to measure again', 0)
            self.buStart.setEnabled(True)

    def onNewSet(self, MySet: CustomData.FourChannels):
        self.rSet = MySet
        self.progressBar.setValue(MySet.i+1)
        if MySet.ts > 0:
            for j in range(4):
                self.livePhasors[j].append(MySet.Data[j].Vc)
        self.plotscatter()
        currentTab = self.tabWidget.master.tabText(self.tabWidget.master.currentIndex())
        if currentTab != 'scatter':
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
        if tat == 'scatter':
            self.plotscatter()
        elif tat == 'msg':
            self.showOutput()
        elif tat == 'raw':
            self.plotraw()
        elif tat == 'PSA':
            self.plotpsa()
        elif tat == 'alpha(f)':
            self.plotalphaf()
        elif tat == 'eta':
            self.ploteta()
        elif tat == 'last status':
            self.showLastStatus()

    def ensureDir(self):
        yymm = datetime.datetime.now().strftime('%y%m')
        dd   = datetime.datetime.now().strftime('%d')
        bd0 = os.path.join(self.datadir, self.SN41, yymm,dd)
        pathlib.Path(bd0).mkdir(parents=True, exist_ok=True)
        return bd0

    def saveData(self, f):
        bd0 = self.ensureDir()
        dt = datetime.datetime.fromtimestamp(self.t0)
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
        C = self.allData.getRawPhasors(f, self.t0)
        fn = 'VOLT_'+valname+'_'+dt.strftime('%Y%m%d_%H%M')+'.dat'
        if not os.path.exists(os.path.join(bd0, fn)):
            with open(os.path.join(bd0, fn), "w") as file:
                file.write('# frequency/Hz t/s  reV1 imV1 reV2 imV2 reV3 imV3 reV4 imV4 reV1set imV1set reV2set imV2set \n')
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
                self.alphafplots[i, j].canvas.ax1.cla()

        freqs = sorted(self.allData.mydict.keys())

        def add_scalar(ax, getter, color):
            f_pts, v_pts, f_bar, v_bar, e_bar = [], [], [], [], []
            for uf in freqs:
                vals = np.array([getter(nd) for nd in self.allData.mydict[uf]], dtype=float)
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

        add_scalar(self.alphafplots[0, 0].canvas.ax1, lambda nd: nd.Res['alpha3mean'], 'r')
        add_scalar(self.alphafplots[0, 1].canvas.ax1, lambda nd: nd.Res['alpha4mean'], 'b')

        from matplotlib.lines import Line2D
        for ax, bx, key, label in [
            (self.alphafplots[1, 0].canvas.ax1, self.alphafplots[1, 0].canvas.bx1, 'gamma3', 'Y₃₂Z₃'),
            (self.alphafplots[1, 1].canvas.ax1, self.alphafplots[1, 1].canvas.bx1, 'gamma4', 'Y₄₂Z₄'),
        ]:
            ax.cla()
            bx.cla()
            f_pts, mag_pts, ang_pts = [], [], []
            f_bar, mag_bar, mag_std, ang_bar, ang_std = [], [], [], [], []
            for uf in freqs:
                g = np.array([1.0 / nd.Res[key] for nd in self.allData.mydict[uf]], dtype=complex)
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
            data = np.asarray(self.rSet.Data[idx].data, dtype=float)
            n = len(data)
            w = np.hanning(n)
            amp = np.abs(np.fft.rfft(data * w)) / (np.sum(w) / 2)
            freqs = np.fft.rfftfreq(n, d=1.0 / fsamp)
            ax.semilogy(freqs[1:], amp[1:], colors[idx] + '-', linewidth=0.5)
            ax.axvline(x=fsig, color='k', linestyle='--', linewidth=0.8)
            ax.set_xscale('log')
            ax.set_xlabel('f (Hz)')
            ax.set_ylabel(labels[idx])
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
