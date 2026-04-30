import os,sys,pathlib
import traceback
import ctypes
import datetime
import time
from shutil import copyfile
import mplwidget
import numpy as np 
import time
import R2FMath
import serial
import serial.tools.list_ports
import TZA
import spectral3
import CCConfig
from Voltfield import VoltField
from Meas3 import Meas
from Tabwidget  import MyTabWidget
import CustomData 
import scipy.optimize
import TZA


from PyQt5.QtCore import (
    Qt,
    QSize,
    QMutex, 
    QObject, 
    QThread, 
    pyqtSignal, 
    pyqtSlot,
    QTimer)

from PyQt5.QtGui import QFont,QIcon
from PyQt5.QtWidgets import (
    QApplication,
    QAction,
    QCheckBox,
    QLabel,
    QMainWindow,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QWidget,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QSpacerItem,
    QSizePolicy,
    QStatusBar,
    QSpinBox,
    QDoubleSpinBox,
    QAbstractSpinBox,
    QProgressBar,
    QMessageBox,
    QDialog,
    QComboBox,
    QLineEdit,
    QTextEdit,
    QRadioButton,
    QButtonGroup

)



class MainWindow(QMainWindow):
    stopDVM  = pyqtSignal()

    def __init__(self,mutex):
        self.tzaport =[]
        for p in serial.tools.list_ports.comports():
            if p.description.startswith('TZA/OPM500'):
                self.tzaport.append(p.device)
        self.tza1 = TZA.TZA(self.tzaport[1])
        self.tza2 = TZA.TZA(self.tzaport[0])
        self.t0=time.time()
        self.config =CCConfig.CCC()
        self.Npts = 8  # This is the number of double points, i.e.  the numbe rof points in the circle
        self.ver = 3.0
        self.quit=False
        self.loopfinished=False
        self.yyyymmdir = r'C:\DATA\CERAMIC\202509'
        self.mytext =[]
        self.mytextmaxlen =1000
        super().__init__()
        self.mutex = mutex
        self.thread = QThread()
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        self.V2 = 1.811+2.578j
        self.V1 = -9.9+0j
        self.dV = 0.1
        self.g1 = 10000
        self.g2 = 100
        self.tza1.set_fgain(self.g1)
        self.tza2.set_fgain(self.g2)
        self.firstgood=False
        
        self.fsigold = -1
        self.rData = CustomData.NPoints(1000,800000)
        self.rSet  = CustomData.FourChannels(1000,800000,2,[],[],[],[],0,0,0,ts=-1)
        self.allData = CustomData.AllData()
        self.myprint(f"Welcome to Version {self.ver}")
        self.parseconfig()

        quit = QAction("Quit", self)
        quit.triggered.connect(self.finishUp) 

        self.scatterplots=np.empty((2,2),dtype=object)
        for i in range(2):
            for j in range(2):
                self.scatterplots[i,j] = mplwidget.MplWidget(rightax=False)
                self.scatterplots[i,j].setfmt("%.5f","%.6f")
        self.rawplots=np.empty((2,2),dtype=object)
        for i in range(2):
            for j in range(2):
                self.rawplots[i,j] = mplwidget.MplWidget(rightax=False)

        self.abplots=np.empty((2,2),dtype=object)
        for i in range(2):
            for j in range(2):
                self.abplots[i,j] = mplwidget.MplWidget(rightax=False)


        self.ciplots=np.empty((2,2),dtype=object)
        for i in range(2):
            for j in range(2):
                self.ciplots[i,j] = mplwidget.MplWidget(rightax=False)

        self.alphafplots=np.empty((1,2),dtype=object)
        for j in range(2):
            self.alphafplots[0,j] = mplwidget.MplWidget(rightax=False)


        self.output =  QTextEdit(self)
        self.output.resize(540, 200)
        self.output.setReadOnly(True)
       
        self.mstatus =  QTextEdit(self)
        self.mstatus.resize(540, 200)
        self.mstatus.setReadOnly(True)

        glayout = QHBoxLayout()
        vlayout= QVBoxLayout()
        vlayout.addItem(QSpacerItem(60,20,QSizePolicy.Fixed,QSizePolicy.Fixed))
        vlayout.addWidget(QLabel('dV'))
        self.VF_dV= VoltField(self,self.dV)
        vlayout.addWidget(self.VF_dV)

        vlayout.addItem(QSpacerItem(60,20,QSizePolicy.Fixed,QSizePolicy.Fixed))

        self.buQuit  =  QPushButton("Quit")
        self.bufp    =  QPushButton("f++")
        self.bufm    =  QPushButton("f--")
        vlayout.addWidget(self.buQuit)
        vlayout.addWidget(self.bufp)
        vlayout.addWidget(self.bufm)
        vlayout.addWidget(QLabel('current f:'))
        self.laf = QLabel('{0:5.1f} kHz'.format(self.fsig/1000))
        vlayout.addWidget(self.laf)
        vlayout.addWidget(QLabel('Re(V1):'))
        self.lareV1 = QLabel("{0:8.3f} V".format(np.real(self.V1)))
        vlayout.addWidget(self.lareV1)
        vlayout.addWidget(QLabel('Re(V1):'))
        self.laimV1 = QLabel('{0:8.3f} V'.format(np.imag(self.V1)))
        vlayout.addWidget(self.laimV1)
        vlayout.addItem(QSpacerItem(60,20,QSizePolicy.Fixed,QSizePolicy.Expanding))


        glayout.addLayout(vlayout)
        self.tabWidget = MyTabWidget(self)   
        glayout.addWidget(self.tabWidget)
        vlayout= QVBoxLayout()
        self.rg1 = QButtonGroup()
        self.rg2 = QButtonGroup()

        self.rb1=[QRadioButton(i) for i in self.tza1.gaindict.values()]
        vlayout.addWidget(QLabel('Gain left'))
        for rb in self.rb1:
            self.rg1.addButton(rb)
            vlayout.addWidget(rb)
            if rb.text()==self.tza1.hgain:
                rb.setChecked(True)

        self.rb2=[QRadioButton(i) for i in self.tza2.gaindict.values()]
        vlayout.addWidget(QLabel('Gain right'))
        for rb in self.rb2:
            self.rg2.addButton(rb)
            vlayout.addWidget(rb)
            if rb.text()==self.tza2.hgain:
                rb.setChecked(True)
        vlayout.addItem(QSpacerItem(60,20,QSizePolicy.Fixed,QSizePolicy.Expanding))
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
        self.sblabel  = QLabel("") 
        widget.layout().addWidget(self.sblabel)
        widget.layout().addWidget(self.progressBar)
        self.statusBar.showMessage('Welcome to CeramicCap',5000)
        self.statusBar.addPermanentWidget(widget) 

        self.buQuit.clicked.connect(self.finishUp)
        self.bufp.clicked.connect(self.fp)
        self.bufm.clicked.connect(self.fm)
        self.goagain()


    def parseconfig(self):
        self.config.read()
        self.C31 = self.config.C31
        self.C32 = self.config.C32
        self.C41 = self.config.C41
        self.C42 = self.config.C42
        self.myprint("Measuring C31={0:3.0e} C32= {1:3.0e} & C41={2:3.0e} C42={3:3.0e}".format(self.C31,self.C32,self.C41,self.C42))

        self.SN31 = self.config.SN31
        self.SN32 = self.config.SN32
        self.SN41 = self.config.SN41
        self.SN42 = self.config.SN42

        self.myprint("S/N: {0} {1} {2} {3}".format(self.SN31,self.SN32,self.SN41,self.SN42))

        self.fsig =self.config.fstart

        self.finelist = self.config.finelist
        self.coarselist = self.config.coarselist
        self.flist = self.config.flist
        self.nrMeas = self.config.nrMeas
        self.datadir = self.config.datadir

    def RBupdate(self):
        for rb in self.rb1:
            if rb.text()==self.tza1.hgain:
                rb.setChecked(True)
        for rb in self.rb2:
            if rb.text()==self.tza2.hgain:
                rb.setChecked(True)
     

    def laUpdate(self):
        self.laf.setText('{0:8.1f} kHz'.format(self.fsig/1000))
        self.lareV1.setText('{0:8.3f} V'.format(np.real(self.V1)))
        self.laimV1.setText('{0:8.3f} V'.format(np.imag(self.V1)))
    
    def on_thread_finished(self):
        pass


    def setf(self,f):
        self.statusBar.showMessage('Next frequency: {0} kHz'.format(f/1000),5000)
        self.fsig = f
        self.allData.deletekey(f)
    

    def fp(self):
        f = self.fsig
        if f in self.flist:
            ix=self.flist.index(f)       
            if ix+1<len(self.flist):
                self.setf(self.flist[ix+1])
            else:
                self.setf(self.flist[0])
    
    def fm(self):
        f = self.fsig
        if f in self.flist:
            ix=self.flist.index(f)        
            if ix-1>=0:
                self.setf(self.flist[ix-1])
            else:
                self.setf(self.flist[-1])

    def rb1Toggled(self,button, checked):
        if checked:
            self.g1=self.tza1.text_to_gain(button.text())

            #self.tza1.set_hgain(button.text())

    def rb2Toggled(self,button, checked):
        if checked:
            self.tza2.set_hgain(button.text())
            self.g2=self.tza2.text_to_gain(button.text())




    def finishUp(self):
        self.statusBar.showMessage('Wait for graceful exit')
        self.myprint("Wait to finish up measurement for graceful exit")
        self.quit=True


    def closeEvent(self, event):
        if self.loopfinished==False:
            self.statusBar.showMessage('Wait for graceful exit')
            self.myprint("Close event called, but loop not finished")
            self.quit=True            
            event.ignore()
        else:
            self.statusBar.showMessage('Quitting')
            self.myprint("Quitting gracefully")
            event.accept()

    def newValues(self):
        self.dV = self.VF_dV.value()

    def onNewData(self,MyData: CustomData.NPoints):
        self.rData = MyData             #recent data
        self.V2 = self.rData.Res['Vz4']
        if self.rData.goodData == True:
            if self.firstgood==True:
                self.allData.append(self.rData)
            else:
                self.firstgood=True     
        self.sblabel.setText("{0}/{1} good measurements at {2:5.2f} kHz".format(\
            self.allData.countf(self.fsig),self.nrMeas,self.fsig/1000))
        if self.fsig == self.fsigold:
            if self.allData.countf(self.fsig)==self.nrMeas:
                self.saveData(self.fsig)
                self.fp()
        self.replot()
 
    def calcVsmall(self,f,V1=-9.0,R=50):
        C42 =self.C42
        C41 =self.C41
        iw = 1j*f*2*np.pi
        I1 = V1*iw*C42/(1+iw*C42*R)
        V2=-I1*(R+1/(iw*C41))
        return V2





    def goagain(self):
        if not self.quit:
            self.thread = QThread()
            self.mydvm = Meas(self.mutex,self,self.Npts)
            self.mydvm.moveToThread(self.thread)
            self.tza1.set_fgain(self.g1)
            self.tza2.set_fgain(self.g2)
            if self.fsigold == self.fsig:
                while np.abs(self.V2)>9:
                    self.V1=self.V1*0.9
                    self.V2=self.V2*0.9
                self.mydvm.storeV(self.V1,self.V2,self.dV,self.fsig,self.g1,self.g2)    
            else:
                self.g1 = R2FMath.newgainvalue1(self.fsig,self.C31,self.C41)
                self.g2 = R2FMath.newgainvalue2(self.fsig,self.C41)
                self.myprint('pick gains {0} {1}'.format(self.g1,self.g2))
                self.tza1.set_fgain(self.g1)
                self.tza2.set_fgain(self.g2)
                self.RBupdate()
                self.firstgood= False
                self.V1 = -9.9+0j
                self.V2 = self.calcVsmall(self.fsig,V1=self.V1)
                while np.abs(self.V2)>9:
                    self.V1 = self.V1*0.9
                    self.V2 = self.calcVsmall(self.fsig,V1=self.V1)
                self.dV = np.abs(self.V2)/200
                #self.myprint('New f={0:8.0f} kHz V1={1:4.3f} + {2:4.3f}i  dV={3:4.3f}'.format(self.fsig/1000,np.real(self.V1),np.imag(self.V1),self.dV))
                self.mydvm.storeV(self.V1,self.V2,self.dV,self.fsig,self.g1,self.g2)        
            self.laUpdate()
            self.mydvm.dataReady.connect(self.onNewData)
            self.mydvm.dataSetReady.connect(self.onNewSet)
            self.mydvm.finished.connect(self.thread.quit)          # stop thread event loop
            self.mydvm.finished.connect(self.mydvm.deleteLater)   # delete worker in its thread
            self.thread.finished.connect(self.goagain)        
            self.stopDVM.connect(self.mydvm.stop)
            self.thread.started.connect(self.mydvm.start) 
            self.thread.finished.connect(self.thread.deleteLater)
            self.thread.start()  
            self.fsigold =self.fsig
        else:
            self.loopfinished =True
            self.close()


    def onNewSet(self,MySet: CustomData.FourChannels):
        self.rSet = MySet             #recent data
        self.progressBar.setValue(MySet.i+1)
        self.replot()

    def myprint(self,newtext):
        n = datetime.datetime.now()
        t = time.time()-self.t0
        t1 = n.strftime('%H:%M:%S')+' >{0:8.1f} < : '.format(t)+ newtext
        t2 = n.strftime('%m/%d')+' '+ t1 +'\n'
        fn = 'CAP'+n.strftime('%m%d%Y')+'.LOG'
        if self.yyyymmdir!='':
            with open(os.path.join(self.yyyymmdir,fn), "a") as file:
                file.write(t2)
        self.mytext.append(t1)
        while len(self.mytext)>self.mytextmaxlen:
            self.mytext.pop(0)


    def replot(self):
        currentIndex=self.tabWidget.master.currentIndex()
        tat = self.tabWidget.master.tabText(currentIndex)
        if tat=='scatter':
           self.plotscatter()
        if tat=='msg':
           self.showOutput()
        if tat=='raw':
           self.plotraw()
        if tat=='alpha(t)':
           self.plotalpha()
        if tat=='circles':
           self.plotcircles()
        if tat=='alpha(f)':
           self.plotalphaf()
        if tat=='last status':           
            self.showLastStatus()

    def ensureDir(self):
        bd = self.datadir
        bd0 = os.path.join(bd,self.SN41)
        long_path = pathlib.Path(bd0)
        long_path.mkdir(parents=True, exist_ok=True)
        return bd0



    def saveData(self,f):
        bd0=self.ensureDir()
        dt = datetime.datetime.fromtimestamp(self.t0)
        valname = self.config.recapdir[self.C41]+'-'+self.config.recapdir[self.C42]
        fn='CC3_'+valname+'_'+dt.strftime('%Y%m%d_%H%M')+'.dat'
        mykeys = ['fsig','ts','alpha3mean','beta3mean','alpha4mean','beta4mean',\
                  'V2cplxcenter','V2cplxradius','V2cplxangle','V3cplxcenter','V3cplxradius','V3cplxangle',\
                  'V4cplxcenter','V4cplxradius','V4cplxangle','V1setamp','V2setcplxcenter','V2setcplxradius','V2setcplxangle',
                  'Vz3','Vz4']
        rdict = self.allData.getkeys(self.fsig,mykeys)
        L=0
        if os.path.exists(os.path.join(bd0,fn))==False:
            with open(os.path.join(bd0,fn), "w") as file:
                o='# '
                for k,v in rdict.items():
                    L = len(v)
                    if isinstance(v[0], complex):
                        o+='Re({0}) Im({0}) '.format(k)
                    else:
                        o+='{0} '.format(k)
                o+='\n'
                file.write(o)
        with open(os.path.join(bd0,fn), "a") as file:
            k =list(rdict)[0]
            L = len(rdict[k])
            for n in range(L):
                o=''
                for k,v in rdict.items():    
                    a=v[n]
                    if  isinstance(a, complex):
                        o+='{0:12.9f} {1:12.9f} '.format(a.real,a.imag)
                    else:
                        o+='{0:12.9f} '.format(a)
                o+='\n'
                file.write(o)
                #file.write('# frequency/Hz t/s  a3 b3 a4 b4  x2 y2 r2  x3 y3 r3 x4 y4 r4 g1 g2 |V2| xV1set yV1set rV1set\n')
        #C=self.allData.getAveVolts(self.fsig,self.t0)
        C=self.allData.getRawPhasors(self.fsig,self.t0)
        fn='VOLT_'+valname+'_'+dt.strftime('%Y%m%d_%H%M')+'.dat'
        if os.path.exists(os.path.join(bd0,fn))==False:
            with open(os.path.join(bd0,fn), "w") as file:
                file.write('# frequency/Hz t/s  reV1 imV1 reV2 imV2 reV3 imV3 reV4 imV4 reV1set imV1set reV2set imV2set \n')
        with open(os.path.join(bd0,fn), "a") as file:
            for n in range(np.shape(C)[0]):
                o = '{0:6.0f} {1:6.1f} '.format(C[n,0],C[n,1])
                for k in range(2,np.shape(C)[1]):
                    o+='{0:15.8f} '.format(C[n,k])
                o+='\n'
                file.write(o)

    def plotscatter(self):
        if self.rData.Res['ts']<=0:
            return
        for j in range(2):
            for i in range(2):
                self.scatterplots[i,j].canvas.ax1.cla()    
        self.scatterplots[0,0].canvas.ax1.plot(np.real(self.rData.ave4[:,0]),\
                                                np.imag(self.rData.ave4[:,0]),'ro')
        self.scatterplots[0,1].canvas.ax1.plot(np.real(self.rData.ave4[:,1]),\
                                                np.imag(self.rData.ave4[:,1]),'go')
        self.scatterplots[1,0].canvas.ax1.plot(np.real(self.rData.ave4[:,2]),\
                                                np.imag(self.rData.ave4[:,2]),'bo')
        self.scatterplots[1,1].canvas.ax1.plot(np.real(self.rData.ave4[:,3]),\
                                                np.imag(self.rData.ave4[:,3]),'mo')
        
        self.scatterplots[0,1].canvas.ax1.scatter(**self.rData.Cir[1].plot_mycircle(),marker='.', s=1,cmap='Greens')
        self.scatterplots[1,0].canvas.ax1.scatter(**self.rData.Cir[2].plot_mycircle(),marker='.', s=1,cmap='Blues')
        self.scatterplots[1,1].canvas.ax1.scatter(**self.rData.Cir[3].plot_mycircle(),marker='.', s=1,cmap='Purples')

            #if self.cbaxes.isChecked():
            #    self.scatterplots[1,0].canvas.ax1.axhline()
            #    self.scatterplots[1,0].canvas.ax1.axvline()
            #    self.scatterplots[1,1].canvas.ax1.axhline()
            #    self.scatterplots[1,1].canvas.ax1.axvline()            
        for j in range(2):
            for i in range(2):
                 self.scatterplots[i,j].canvas.draw()

    def plotalpha(self):
        if self.allData.countf(self.fsig)==0:
            return
        for j in range(2):
            for i in range(2):
                self.abplots[i,j].canvas.ax1.cla()    
        mul =1e6
        mykeys = ['ts','alpha3mean','beta3mean','alpha4mean','beta4mean']
        rdict = self.allData.getkeys(self.fsig,mykeys)
        self.abplots[0,0].canvas.ax1.plot(rdict['ts']-self.t0,rdict['alpha3mean']*mul,'ro')
        self.abplots[0,1].canvas.ax1.plot(rdict['ts']-self.t0,rdict['beta3mean']*mul,'bo')

        self.abplots[1,0].canvas.ax1.plot(rdict['ts']-self.t0,rdict['alpha4mean']*mul,'ro')
        self.abplots[1,1].canvas.ax1.plot(rdict['ts']-self.t0,rdict['beta4mean']*mul,'bo')
        for j in range(2):
            for i in range(2):
                 self.abplots[i,j].canvas.draw()


    def plotalphaf(self):
        if self.allData.count()==0:
            return
        for j in range(2):
            self.alphafplots[0,j].canvas.ax1.cla() 

        mykeys = ['fsig','alpha3mean','beta3mean','alpha4mean','beta4mean']
        rdict = self.allData.getdictf(mykeys)

        self.alphafplots[0,0].canvas.ax1.plot(rdict['fsig'],rdict['alpha3mean'],'ro')
        self.alphafplots[0,1].canvas.ax1.plot(rdict['fsig'],rdict['alpha4mean'],'bo')
        self.alphafplots[0,0].canvas.ax1.set_xscale('log')
        self.alphafplots[0,1].canvas.ax1.set_xscale('log')
        for j in range(2):
            self.alphafplots[0,j].canvas.draw()


    def plotcircles(self):
        if self.allData.countf(self.fsig)==0:
            return
        for j in range(2):
            for i in range(2):
                self.ciplots[i,j].canvas.ax1.cla()    
        mykeys =['ts','V2cplxradius','V3cplxradius','V4cplxradius','V1setcplxradius']
        rdict = self.allData.getkeys(self.fsig,mykeys)
        self.ciplots[0,0].canvas.ax1.plot(rdict['ts']-self.t0,np.abs(rdict['V2cplxradius']),'ro')
        self.ciplots[0,1].canvas.ax1.plot(rdict['ts']-self.t0,np.abs(rdict['V1setcplxradius']),'bo')

        self.ciplots[1,0].canvas.ax1.plot(rdict['ts']-self.t0,np.abs(rdict['V3cplxradius']),'ro')
        self.ciplots[1,1].canvas.ax1.plot(rdict['ts']-self.t0,np.abs(rdict['V4cplxradius']),'bo')

        self.ciplots[0,0].canvas.ax1.set_ylabel('r(V2 meas)')
        self.ciplots[0,1].canvas.ax1.set_ylabel('r(V1 set)')
        self.ciplots[1,0].canvas.ax1.set_ylabel('r(V3 meas)')
        self.ciplots[1,1].canvas.ax1.set_ylabel('r(V4 meas)')

        for j in range(2):
            for i in range(2):
                 self.ciplots[i,j].canvas.draw()


    def plotraw(self):
        if self.rSet.ts<=0:
            return
        for j in range(2):
            for i in range(2):
                self.rawplots[i,j].canvas.ax1.cla()    

        self.t = np.arange(len(self.rSet.Data[0].data))
        Vc=self.rSet.Data[0].Vc
        phi = -np.angle(-1j*Vc)
        while phi<0:
            phi=phi+2*np.pi
        t = phi/(2*np.pi*self.rSet.fsig/self.rSet.fsamp)
        be1=int(t)
        be=be1
        en =int(self.rSet.fsamp/self.rSet.fsig*2)+be1
        #be=0
        #en=-1

        self.rawplots[0,0].canvas.ax1.plot(self.t[be:en],self.rSet.Data[0].data[be:en],'r.')
        self.rawplots[0,1].canvas.ax1.plot(self.t[be:en],self.rSet.Data[1].data[be:en],'g.')
        self.rawplots[1,0].canvas.ax1.plot(self.t[be:en],self.rSet.Data[2].data[be:en],'b.')
        self.rawplots[1,1].canvas.ax1.plot(self.t[be:en],self.rSet.Data[3].data[be:en],'m.') 
        self.rawplots[0,0].canvas.ax1.plot(self.t[be:en],self.rSet.Data[0].fv[be:en],'k-')
        self.rawplots[0,1].canvas.ax1.plot(self.t[be:en],self.rSet.Data[1].fv[be:en],'k-')
        self.rawplots[1,0].canvas.ax1.plot(self.t[be:en],self.rSet.Data[2].fv[be:en],'k-')
        self.rawplots[1,1].canvas.ax1.plot(self.t[be:en],self.rSet.Data[3].fv[be:en],'k-')
        for j in range(2):
            for i in range(2):
                 self.rawplots[i,j].canvas.draw()

    def showOutput(self):
        o=''
        for n, i in enumerate(self.mytext):
            o+='{1}\n'.format(n,i)
        self.output.setText(o)
        scrollbar = self.output.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

    def showLastStatus(self):
        if self.rData.Res['ts']>0:
            o=''
            for k,v in self.rData.Res.items():
                if isinstance(v, complex):
                    if v.imag>0:
                        o+='{0:30}:{1:8.3f} + i {2:8.3f}\n'.format(k,v.real,v.imag)
                    else:
                        o+='{0:30}:{1:8.3f} - i {2:8.3f}\n'.format(k,v.real,-v.imag)
                elif isinstance(v, float):
                    o+='{0:30}:{1:8.3f}\n'.format(k,v)
                elif isinstance(v, int):
                    o+='{0:30}:{1}\n'.format(k,v)
            self.mstatus.setText(o)
            scrollbar = self.mstatus.verticalScrollBar()
            scrollbar.setValue(scrollbar.maximum())

       
def excepthook(exc_type, exc_value, exc_tb):
    tb = "".join(traceback.format_exception(exc_type, exc_value, exc_tb))
    now = datetime.datetime.now()
    print("error catched at ",now)
    print("error message:\n", tb)
    fi=open('R2F_errlog.dat','w')
    fi.write("error catched!:\n")
    fi.write("error message:\n"+ tb)
    fi.close()
    app.quit()
    #QApplication.quit()
    # or QtWidgets.QApplication.exit(0)

if __name__=='__main__':   
    mutex = QMutex()
    myappid = 'mycompany.myproduct.subproduct.version' # arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    app = QApplication(sys.argv)
    sys.excepthook = excepthook
    
    
    window = MainWindow(mutex)
    window.show()
    app.exec()
    sys.exit()
