import os,sys
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
from rightwindow import RightWindowWidget
from BridgeConfig import BridgeConfig
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
        self.ver = 3.0
        self.quit=False
        self.yyyymmdir = r'C:\DATA\CERAMIC\202509'
        self.mytext =[]
        self.mytextmaxlen =1000
        super().__init__()
        self.mutex = mutex
        self.thread = QThread()
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        self.V1 = 1.811+2.578j
        self.V2 = -9.9+0j
        self.dV = 0.1
        self.g1 = 1
        self.g2 = 1
        self.tza1.set_fgain(self.g1)
        self.tza2.set_fgain(self.g2)
        self.fsig =float(self.config.meas['FSTART'])
        self.fsigold = self.fsig
        self.rData = CustomData.EightPoints(1000,800000)
        self.rSet  = CustomData.FourChannels(1000,800000,2,[],[],[],[],0,0,0,ts=-1)
        self.allData = CustomData.AllData()

        self.finelist = [float(i) for i in self.config.meas['FINELIST'].split(',')]
        self.coarselist = [float(i) for i in self.config.meas['COARSELIST'].split(',')]

        if self.config.meas['USEFINE']==True:
            self.flist = self.finelist
        else:
            self.flist = self.coarselist
         

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


        self.output =  QTextEdit(self)
        self.output.resize(540, 200)
        self.output.setReadOnly(True)
       
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
        self.laf = QLabel('{0:8.0f} kHz'.format(self.fsig/1000))
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
        glayout.addLayout(vlayout)        
        self.rg1.buttonToggled.connect(self.rb1Toggled)
        self.rg2.buttonToggled.connect(self.rb2Toggled)

        #self.rightwindow = RightWindowWidget(self)     
        #glayout.addWidget(self.rightwindow)

        #glayout.addWidget(self.rightwindow)
        central_widget.setLayout(glayout)  
        self.myprint(f"Welcome to Version {self.ver}")
        
        self.progressBar = QProgressBar()
        self.progressBar.setMaximumWidth(400)
        self.progressBar.setRange(0, 8)
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


    def laUpdate(self):
        self.laf.setText('{0:8.0f} kHz'.format(self.fsig/1000))
        self.lareV1.setText('{0:8.3f} V'.format(np.real(self.V1)))
        self.laimV1.setText('{0:8.3f} V'.format(np.imag(self.V1)))
    
    def on_thread_finished(self):
        pass


    def setf(self,f):
        self.statusBar.showMessage('Next frequency: {0} kHz'.format(f/1000),5000)
        self.fsig = f
        self.myprint('setf in main program {0}'.format(self.fsig/1000))


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
            print(f"rb1 {button.text()=} {self.g1=}")

            #self.tza1.set_hgain(button.text())

    def rb2Toggled(self,button, checked):
        if checked:
            self.tza2.set_hgain(button.text())
            self.g2=self.tza2.text_to_gain(button.text())
            print(f"rb2 {button.text()=} {self.g2=}")




    def finishUp(self):
        self.quit=True

    def closeEvent(self, event):
        self.myprint("Quitting gracefully")

    def newValues(self):
        self.dV = self.VF_dV.value()

    def onNewData(self,MyData: CustomData.EightPoints):
        self.rData = MyData             #recent data
        self.V1 = self.rData.Vz3
        if self.rData.goodData == True:
            self.allData.append(self.rData)
        if self.fsig == self.fsigold:
            if self.allData.countf(self.fsig)==10:
                self.fp()
        self.replot()

    def goagain(self):
        if not self.quit:
            self.thread = QThread()
            self.mydvm = Meas(self.mutex,self)
            self.mydvm.moveToThread(self.thread)
            self.tza1.set_fgain(self.g1)
            self.tza2.set_fgain(self.g2)
            self.mydvm.storeV(self.V1,self.V2,self.dV,self.fsig)    
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
            self.close()


    def onNewSet(self,MySet: CustomData.FourChannels):
        self.rSet = MySet             #recent data
        self.progressBar.setValue(MySet.i+1)
        self.replot()

    def myprint(self,newtext):
        n = datetime.datetime.now()
        t = time.time()-self.t0
        t1 = n.strftime('%H:%M:%S')+' >{0:8.1f} < : '.format(t)+ newtext
        t2 = n.strftime('%m/%d')+ t1 +'\n'
        if self.yyyymmdir!='':
            with open(os.path.join(self.yyyymmdir,'log.dat'), "a") as file:
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
           self.showoutput()
        if tat=='raw':
           self.plotraw()
        if tat=='alpha':
           self.plotalpha()
        if tat=='circles':
           self.plotcircles()

    def plotscatter(self):
        if self.rData.ts<=0:
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
        
        self.scatterplots[0,1].canvas.ax1.plot(*self.rData.Cir[1].plot_mycircle(),'g-')
        self.scatterplots[1,0].canvas.ax1.plot(*self.rData.Cir[2].plot_mycircle(),'b-')
        self.scatterplots[1,1].canvas.ax1.plot(*self.rData.Cir[3].plot_mycircle(),'m-')

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

        t,a3 = self.allData.geta3(self.fsig,self.t0)
        t,b3 = self.allData.getb3(self.fsig,self.t0)
        t,a4 = self.allData.geta4(self.fsig,self.t0)
        t,b4 = self.allData.getb4(self.fsig,self.t0)                        

        self.abplots[0,0].canvas.ax1.plot(t,a3*1e6,'ro')
        self.abplots[0,1].canvas.ax1.plot(t,a4*1e6,'bo')

        self.abplots[1,0].canvas.ax1.plot(t,b3*1e6,'ro')
        self.abplots[1,1].canvas.ax1.plot(t,b4*1e6,'bo')
        for j in range(2):
            for i in range(2):
                 self.abplots[i,j].canvas.draw()


    def plotcircles(self):
        if self.allData.countf(self.fsig)==0:
            return
        for j in range(2):
            for i in range(2):
                self.ciplots[i,j].canvas.ax1.cla()    

        t, par = self.allData.getCircles(self.fsig,self.t0)

        self.ciplots[0,0].canvas.ax1.plot(t,par[:,0],'ro')
        self.ciplots[0,1].canvas.ax1.plot(t,par[:,2],'bo')

        self.ciplots[1,0].canvas.ax1.plot(t,par[:,1],'ro')
        self.ciplots[1,1].canvas.ax1.plot(t,par[:,3],'bo')
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
        phi = np.angle(Vc)
        if phi<0: phi=phi+2*np.pi
        t = phi/(2*np.pi)*self.rSet.fsamp/self.rSet.fsig
        be1=int(t)
        be=be1
        en =int(self.rSet.fsamp/self.rSet.fsig*2)+be1

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

    def showoutput(self):
        o=''
        for n, i in enumerate(self.mytext):
            o+='{1}\n'.format(n,i)
        self.output.setText(o)
        scrollbar = self.output.verticalScrollBar()
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
