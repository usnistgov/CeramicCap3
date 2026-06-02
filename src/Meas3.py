
from PyQt5.QtCore import (
    QObject,
    pyqtSignal,
    pyqtSlot)

import time
import threading
import numpy as np
import pyvisa
import CustomData
import datetime
import os


SG1_ADDR = 'USB0::0x0957::0x2C07::MY57801011::0::INSTR'
DVM_ADDR = 'USB0::0x2A8D::0x8601::MY59000912::0::INSTR'


def reset_instruments():
    """Reset signal generator and scanner at startup, then enable both SG outputs."""
    rm = pyvisa.ResourceManager()
    try:
        sg1 = rm.open_resource(SG1_ADDR)
        sg1.timeout = 15000
        sg1.write('*RST')
        sg1.query('*OPC?')
        sg1.write('OUTP1 ON')
        sg1.write('OUTP2 ON')
        sg1.close()
    except Exception:
        pass
    try:
        dvm = rm.open_resource(DVM_ADDR)
        dvm.timeout = 15000
        dvm.write('*RST')
        dvm.query('*OPC?')
        dvm.close()
    except Exception:
        pass
    rm.close()


def open_instruments():
    """Open and return (rm, sg1, dvm) for a measurement sweep."""
    rm = pyvisa.ResourceManager()
    sg1 = rm.open_resource(SG1_ADDR)
    sg1.timeout = 25000
    dvm = rm.open_resource(DVM_ADDR)
    dvm.timeout = 25000
    return rm, sg1, dvm


def close_instruments(rm, sg1, dvm):
    """Close sweep instruments; swallows errors so the caller is never blocked."""
    for res in (sg1, dvm):
        try:
            res.close()
        except Exception:
            pass
    try:
        rm.close()
    except Exception:
        pass


class Meas(QObject):
    finished = pyqtSignal()
    dataReady = pyqtSignal(CustomData.NPoints)
    dataSetReady = pyqtSignal(CustomData.ThreeChannels)
    logMessage = pyqtSignal(str)

    def __init__(self, mutex, NDpts, rawdatadir=r'c:\RAWDATA', saverawdata=False, fsamp=800000, nsamp=80000, chunk_periods=0, max_nhars=10, v4_range=3, switch_normal=True, fit_cache=None, sg1=None, dvm=None, rm=None):
        self.bdraw = rawdatadir
        self.saverawdata = saverawdata
        self.chunk_periods = chunk_periods
        self.Nhars = max_nhars
        super().__init__()
        self.NDpts = NDpts
        self.Npts = NDpts
        self.mutex = mutex
        self.isidle = True
        self._stop = False
        self._stop_event = threading.Event()
        self.runt = time.time()
        self.fsig = 1000
        self.fsamp = fsamp
        self.nsamp = nsamp
        self.V1_setpoints = None
        self.v4_range = v4_range        # DAQ voltage range for ch103 (TZA); normally 10 V, 100 V if needed
        self.switch_normal = switch_normal
        self.initial_fit_cache = fit_cache if fit_cache is not None else {}
        if sg1 is not None and dvm is not None:
            # Caller owns the connection; Meas must not close it.
            self._owns_instruments = False
            self.rm  = rm
            self.sg1 = sg1
            self.dvm = dvm
        else:
            self._owns_instruments = True
            self.rm  = pyvisa.ResourceManager()
            self.sg1 = self.rm.open_resource(SG1_ADDR)
            self.dvm = self.rm.open_resource(DVM_ADDR)
        self.precmd()


    def switch(self,switch_normal=True):
        self.switch_normal =switch_normal
        if self.switch_normal:
            self.dvm.write('ROUT:OPEN (@218,241)')
            self.dvm.write('ROUT:CLOS (@211,248)')        
        else:
            self.dvm.write('ROUT:OPEN (@211,248)')
            self.dvm.write('ROUT:CLOS (@218,241)')        



    def precmd(self):
        self.sg1.write('UNIT:ANGL DEG')
        self.sg1.write('VOLT:UNIT VPP')
        self.sg1.write('VOLT:OFFS 0')
        self.sg1.write('SOUR1:FUNC SIN')
        self.sg1.write('SOUR2:FUNC SIN')
        self.dvm.timeout = 25000

        self.dvm.write('FORM3 REAL')
        #self.dvm.write('ACQ3:VOLT 3,DIFF,AC,TIME,(@101)')    # V1 small signal (~1 V)
        self.dvm.write('ACQ3:VOLT 18,DIFF,AC,TIME,(@101:103)')   # V2 large signal (~9.9 V)
        #self.dvm.write(f'ACQ3:VOLT {self.v4_range},DIFF,AC,TIME,(@103)')  # TZA
        self.dvm.write('SAMP3:RATE {0:8.2f},(@101:103)'.format(self.fsamp))
        self.dvm.write('SAMP3:COUN {0},(@101:103)'.format(self.nsamp))
        #self.dvm.write('INP3:COUP AC,(@101:103)')
        self.dvm.write('TRIG3:SOUR BUS,(@101:103)')

    def write1dbg(self, ostr, debug=False):
        self.sg1.write(ostr)
        if debug:
            print(ostr)
        time.sleep(0.001)

    def storeV(self, V1c, V2c, dV1, fsig, g2):
        if self.isidle:
            self.V1c = V1c   # small modulated center
            self.V2c = V2c   # large constant
            self.dV1 = dV1   # modulation radius for V1
            self.fsig = fsig
            self.g2 = g2

    def storeV1pts(self, V1pts):
        if self.isidle:
            self.V1_setpoints = np.asarray(V1pts, dtype=complex)

    def prepForMeas(self):
        if self.co == 0:
            self.logMessage.emit("V1= {0:8.3f}  V2={1:8.3f} dV1={2:8.3f}  f={3:5.1f} kHz".format(
                self.V1c, self.V2c, self.dV1, self.fsig/1000))
        if self.V1_setpoints is not None:
            self.V1 = self.V1_setpoints[self.co]
        else:
            ang = self.co/self.NDpts*2*np.pi
            a = self.dV1*0.8
            b = self.dV1*1.2
            theta = 30/180*np.pi
            self.V1 = self.V1c + np.exp(1j*theta)*(a*np.cos(ang)+1j*b*np.sin(ang))
        self.V2 = self.V2c   # V2 is the large constant voltage
        V1amp = np.abs(self.V1)
        V2amp = np.abs(self.V2)
        if V2amp >= 10.0:
            self.logMessage.emit("Error V2>10 V -- rescaling")
            self.V2 = 9/np.abs(self.V2)*self.V2
            V2amp = np.abs(self.V2)
        elif V1amp > 10.0:
            self.logMessage.emit("Error V1>10 V -- rescaling")
            self.V1 = 9.9/np.abs(self.V1)*self.V1
            V1amp = np.abs(self.V1)
        V1phase = np.angle(self.V1)/np.pi*180
        V2phase = np.angle(self.V2)/np.pi*180
        self.write1dbg('SOUR1:VOLT {0:8.4f}'.format(V1amp))
        self.write1dbg('SOUR2:VOLT {0:8.4f}'.format(V2amp))
        self.write1dbg('SOUR1:FREQ {0:8.4f}'.format(self.fsig))
        self.write1dbg('SOUR2:FREQ {0:8.4f}'.format(self.fsig))
        self.write1dbg('PHAS:SYNC')
        self.write1dbg('SOUR1:PHASE {0:8.4f}'.format(V1phase))
        self.write1dbg('SOUR2:PHASE {0:8.4f}'.format(V2phase))

 
        ch1amp   = float(self.sg1.query('SOUR1:VOLT?'))
        ch1phase = float(self.sg1.query('SOUR1:PHASE?'))
        ch2amp   = float(self.sg1.query('SOUR2:VOLT?'))
        ch2phase = float(self.sg1.query('SOUR2:PHASE?'))
 
        self.V1rb = ch1amp * np.exp(1j*ch1phase/180*np.pi)  # V1 on SOUR1
        self.V2rb = ch2amp * np.exp(1j*ch2phase/180*np.pi)  # V2 on SOUR2
        ret = self.sg1.query('SYSTem:ERRor?')
        retval = int(ret.split(',')[0])
        if retval != 0:
            self.logMessage.emit(f'Error: {ret}')

    def sendtrig(self):
        self.dvm.write('INIT3 (@101:103)')
        self.dvm.write('*TRG')

    @pyqtSlot()
    def start(self):
        self._stop = False
        self._stop_event.clear()
        self.isidle = False
        self.switch(self.switch_normal)
        self._stop_event.wait(1.0)      # relay settle time
        if self._stop:
            self.isidle = True
            self.finished.emit()
            return
        self.rawN = CustomData.NPoints(self.fsig, self.fsamp, Nhars=self.Nhars, g2=self.g2,
                                       N=self.Npts, chunk_periods=self.chunk_periods,
                                       fit_cache=self.initial_fit_cache)
        for self.co in range(self.Npts):
            if self._stop:
                break
            self.prepForMeas()
            if self._stop:          # stop pressed during SG setup — skip this point
                break
            self.sendtrig()
            acq_time = self.nsamp / self.fsamp
            self._stop_event.wait(acq_time + 0.5)  # interruptible; +0.5 s margin for trigger latency
            try:
                self.getvals()
            except Exception as e:
                self.logMessage.emit(f"Point {self.co} failed ({type(e).__name__}: {e}) — discarding set")
            if self._stop:
                break
        self.isidle = True
        if not self._stop and self.rawN.is_complete:
            self.rawN.calc()
            self.dataReady.emit(self.rawN)
        elif not self._stop:
            self.logMessage.emit(
                f'set incomplete: {self.rawN._n_filled}/{self.rawN.N} pts acquired — skipping')
        self.finished.emit()
        if self._owns_instruments:
            for res in (self.sg1, self.dvm):
                try:
                    res.close()
                except Exception:
                    pass
            try:
                self.rm.close()
            except Exception:
                pass

    @pyqtSlot()
    def stop(self):
        self._stop = True
        self._stop_event.set()  # wake the wait() immediately; ABOR sent by the worker

    def _fetch(self, chan):
        """Fetch one digitizer channel; logs if the transfer takes unexpectedly long."""
        t0 = time.time()
        data = self.dvm.query_binary_values(f'FETCH3? (@{chan})', datatype='f', is_big_endian=True)
        dt = time.time() - t0
        if dt > 2.0:
            self.logMessage.emit(f'  FETCH @{chan} slow: {dt:.2f}s')
        return data

    def getvals(self):
        if self.switch_normal:
            ch1 = self._fetch(101)
            ch2 = self._fetch(102)
        else:
            ch2 = self._fetch(101)
            ch1 = self._fetch(102)

        ch3 = self._fetch(103)
        if self.saverawdata:
            now = datetime.datetime.now()
            timestamp = now.strftime("%Y%m%d_%H%M%S")
            bd = os.path.join(self.bdraw, now.strftime("%Y"), now.strftime("%m"), now.strftime("%d"))
            os.makedirs(bd, exist_ok=True)
            fn = os.path.join(bd, f"ceramic_raw_{timestamp}.npz")
            np.savez_compressed(fn, ch1=ch1, ch2=ch2, ch3=ch3)
        self.rawN.setPoint(self.co, ch1, ch2, ch3, self.V1rb, self.V2rb, time.time())
        self.dataSetReady.emit(self.rawN.Data[self.co])
