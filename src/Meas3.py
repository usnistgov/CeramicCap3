
from PyQt5.QtCore import (
QObject, 
pyqtSignal,
pyqtSlot,
QTimer)

import time
import numpy as np 
import pyvisa
import R2FMath
import CustomData

class Meas(QObject):
    finished = pyqtSignal()
    dataReady  = pyqtSignal(CustomData.NPoints)
    dataSetReady  = pyqtSignal(CustomData.FourChannels)

    def __init__(self,mutex,parent,NDpts):
        super(QObject, self).__init__()
        self.NDpts = NDpts  #4 for a circle double points means in both switch positions
        self.Npts = 2*NDpts #8 for a circle double points means in both switch positions
        self.par =parent
        self.isidle =True
        self.mutex       = mutex
        self.co=-1
        self.Nhars=2
        self.runt = time.time()
        self.fsig = 1000
        self.fsamp = 800000 
        self.maxco=1200*2
        self.rm = pyvisa.ResourceManager()
        sg1_name = 'USB0::0x0957::0x2C07::MY57801011::0::INSTR' #keysight
        #sg1_name = 'USB0::0x1AB1::0x0642::DG1ZA192202024::0::INSTR'#rigol
        self.sg1 = self.rm.open_resource(sg1_name)
        dvm_name ='USB0::0x2A8D::0x8601::MY59000912::0::INSTR'
        self.dvm = self.rm.open_resource(dvm_name)
        #self.dvm.write('*RST')
        self.oldf=-1
        self.precmd()       

    def precmd(self):
        self.sg1.write('UNIT:ANGL DEG')
        self.sg1.write('VOLT:UNIT VPP')
        self.sg1.write('VOLT:OFFS 0')
        self.sg1.write('SOUR1:FUNC SIN')
        self.sg1.write('SOUR2:FUNC SIN')
        self.dvm.timeout=25000
        self.dvm.write('FORM3 REAL')
        self.dvm.write('ACQ3:VOLT 10,(@101:104)')
        self.dvm.write('SAMP3:RATE {0:8.2f},(@101:104)'.format(self.fsamp))
        self.dvm.write('SAMP3:COUN {0},(@101:104)'.format(int(self.fsamp/10)))
        self.dvm.write('INP3:COUP AC,(@101:104)')
        self.dvm.write('TRIG3:SOUR BUS,(@101:104)')

    def write1dbg(self,ostr,debug=False):
        self.sg1.write(ostr)
        if debug:
            print(ostr)
        time.sleep(0.001)
 
    def storeV(self,V1c,V2c,dV1,fsig,g1,g2):
        if self.isidle:
            self.V1c=V1c
            self.V2c=V2c
            self.dV1 = dV1
            self.fsig = fsig
            self.g1 = g1
            self.g2 = g2
        
    def prepForMeas(self):
        if self.co==0:
            self.par.myprint("V1= {0:8.3f}  V2={1:8.3f} dV1={2:8.3f}  f={3:5.1f} kHz".format(self.V1c,self.V2c,self.dV1,self.fsig/1000))
        if self.co%2==0:
            self.dvm.write('ROUT:OPEN (@211,248)')   
            self.dvm.write('ROUT:CLOS (@218,241)') # 8->1 1->4
        else:
            self.dvm.write('ROUT:OPEN (@218,241)')   
            self.dvm.write('ROUT:CLOS (@211,248)') # 1->1 8->4
        ang = (self.co//2)/self.NDpts*2*np.pi
        a = self.dV1*1.2
        b = self.dV1*0.8
        theta = 30/180*np.pi
        self.V1 = self.V1c+np.exp(1j*theta)*(a*np.cos(ang)+1j*b*np.sin(ang))
        #if self.co%8==0 or self.co%8==1:
            #self.V1 = self.V1c+self.dV1
        #elif self.co%8==2 or self.co%8==3:
            #self.V1 = self.V1c-self.dV1
        #elif self.co%8==4 or self.co%8==5:
            #self.V1 = self.V1c+1j*self.dV1
        #else:
            #self.V1 = self.V1c-1j*self.dV1
        self.V2 = self.V2c
        V1amp   = np.abs(self.V1)
        V2amp   = np.abs(self.V2)
        if V1amp>=10.0:
            print('Error V1>10 V -- rescaling')
            self.V1 = 9/np.abs(self.V1)*self.V1
            V1amp   = np.abs(self.V1)
        elif V2amp>10.0:
            print('Error V2>10 V -- rescaling')
            self.V2 = 9.9/np.abs(self.V2)*self.V2
            V2amp   = np.abs(self.V2)
        V1phase = np.angle(self.V1)/np.pi*180
        V2phase = np.angle(self.V2)/np.pi*180
        self.write1dbg('SOUR1:VOLT {0:8.4f}'.format(V1amp))
        self.write1dbg('SOUR2:VOLT {0:8.4f}'.format(V2amp))
        self.write1dbg('SOUR1:FREQ {0:8.4f}'.format(self.fsig))
        self.write1dbg('SOUR2:FREQ {0:8.4f}'.format(self.fsig))
        self.write1dbg('PHAS:SYNC')
        self.write1dbg('SOUR1:PHASE {0:8.4f}'.format(V1phase))
        self.write1dbg('SOUR2:PHASE {0:8.4f}'.format(V2phase))
        #self.par.myprint(f"cmd: Amp:{V1amp:.3f} Phase:{V1phase:.3f}")
        V1=float(self.sg1.query('SOUR1:VOLT?'))
        V2=float(self.sg1.query('SOUR2:VOLT?'))
        phase1=float(self.sg1.query('SOUR1:PHASE?'))
        phase2=float(self.sg1.query('SOUR2:PHASE?'))
        self.V1rb = V1 * np.exp(1j*phase1/180*np.pi)
        self.V2rb = V2 * np.exp(1j*phase2/180*np.pi)
        ret=self.sg1.query('SYSTem:ERRor?')
        retval = int(ret.split(',')[0])
        if retval!=0:
            self.par.myprint(f'Error: {ret}')


        
    def sendtrig(self):
        self.dvm.write('INIT3 (@101:104)')# %%
        self.dvm.write('*TRG')
       
    
    @pyqtSlot()
    def start(self):
        self.isidle=False
        start = time.time()
        self.rawN = CustomData.NPoints(self.fsig,self.fsamp,g1=self.g1,g2=self.g2,N=self.Npts)
        for self.co in range(self.Npts):
            self.prepForMeas()
            self.sendtrig()
            time.sleep(1)
            V=self.getvals()
        stop = time.time()
        self.isidle=True
        self.rawN.calc()
        self.dataReady.emit(self.rawN)
        self.finished.emit()

    @pyqtSlot()
    def stop(self):
        print('Good bye')


    def readCh(self,ch):
        """ ch = 1,2,3,4 """
        if ch not in (1,2,3,4):
            print('Channel not valid, you sent: {0}'.format(ch))
        ostr ='FETCH3? (@10{0})'.format(ch)
        values = self.dvm.query_binary_values(ostr, \
                datatype='f', is_big_endian=True)
        return values  
      
    def getvals(self):
        if self.co%2==0:
            ch1=self.dvm.query_binary_values('FETCH3? (@101)',  datatype='f', is_big_endian=True)
            ch2=self.dvm.query_binary_values('FETCH3? (@102)',  datatype='f', is_big_endian=True)
        else:
            ch2=self.dvm.query_binary_values('FETCH3? (@101)',  datatype='f', is_big_endian=True)
            ch1=self.dvm.query_binary_values('FETCH3? (@102)',  datatype='f', is_big_endian=True)
        ch3=self.dvm.query_binary_values('FETCH3? (@103)',  datatype='f', is_big_endian=True)
        ch4=self.dvm.query_binary_values('FETCH3? (@104)',  datatype='f', is_big_endian=True)
        self.rawN.setPoint(self.co,ch1,ch2,ch3,ch4,self.V1rb,self.V2rb,time.time())
        self.dataSetReady.emit(self.rawN.Data[self.co])
