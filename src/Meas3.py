
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


class Meas(QObject):
    finished = pyqtSignal()
    dataReady = pyqtSignal(CustomData.NPoints)
    dataSetReady = pyqtSignal(CustomData.FourChannels)
    logMessage = pyqtSignal(str)

    def __init__(self, mutex, NDpts, rawdatadir=r'c:\RAWDATA', saverawdata=False, fsamp=800000, nsamp=80000, switching=True, chunk_periods=0, max_nhars=10, v4_range=10):
        self.bdraw = rawdatadir
        self.saverawdata = saverawdata
        self.switching = switching
        self.chunk_periods = chunk_periods
        self.Nhars = max_nhars
        super().__init__()
        self.NDpts = NDpts
        self.Npts = 2*NDpts
        self.mutex = mutex
        self.isidle = True
        self._stop = False
        self._stop_event = threading.Event()
        self.runt = time.time()
        self.fsig = 1000
        self.fsamp = fsamp
        self.nsamp = nsamp
        self.V1_setpoints = None
        self._relay_pos = -1   # tracks last commanded relay position to avoid redundant switches
        self.ch3_range = 0.3            # current DVM range for ch103; bumped on overflow
        self._ch3_range_pending = False # set when range changed; triggers reconfigure in prepForMeas
        self.v4_range = v4_range        # DAQ voltage range for ch104 (V4); normally 10 V, 100 V if needed
        self.rm = pyvisa.ResourceManager()
        self.sg1 = self.rm.open_resource(SG1_ADDR)
        self.dvm = self.rm.open_resource(DVM_ADDR)
        self.precmd()

    def precmd(self):
        self.sg1.write('UNIT:ANGL DEG')
        self.sg1.write('VOLT:UNIT VPP')
        self.sg1.write('VOLT:OFFS 0')
        self.sg1.write('SOUR1:FUNC SIN')
        self.sg1.write('SOUR2:FUNC SIN')
        self.dvm.timeout = 25000
        self.dvm.write('FORM3 REAL')
        self.dvm.write('ACQ3:VOLT 10,(@101,102)')
        self.dvm.write(f'ACQ3:VOLT {self.v4_range},(@104)')
        self.dvm.write('ACQ3:VOLT 0.3,(@103)')
        self.dvm.write('SAMP3:RATE {0:8.2f},(@101:104)'.format(self.fsamp))
        self.dvm.write('SAMP3:COUN {0},(@101:104)'.format(self.nsamp))
        self.dvm.write('INP3:COUP AC,(@101:104)')
        self.dvm.write('TRIG3:SOUR BUS,(@101:104)')

    def write1dbg(self, ostr, debug=False):
        self.sg1.write(ostr)
        if debug:
            print(ostr)
        time.sleep(0.001)

    def storeV(self, V1c, V2c, dV1, fsig, g1, g2):
        if self.isidle:
            self.V1c = V1c   # small modulated center
            self.V2c = V2c   # large constant
            self.dV1 = dV1   # modulation radius for V1
            self.fsig = fsig
            self.g1 = g1
            self.g2 = g2

    def storeV1pts(self, V1pts):
        if self.isidle:
            self.V1_setpoints = np.asarray(V1pts, dtype=complex)

    def prepForMeas(self):
        if self.co == 0:
            self.logMessage.emit("V1= {0:8.3f}  V2={1:8.3f} dV1={2:8.3f}  f={3:5.1f} kHz".format(
                self.V1c, self.V2c, self.dV1, self.fsig/1000))
        if self._ch3_range_pending:
            self.dvm.write(f'ACQ3:VOLT {self.ch3_range},(@103)')
            self.dvm.write('SAMP3:RATE {0:8.2f},(@103)'.format(self.fsamp))
            self.dvm.write('SAMP3:COUN {0},(@103)'.format(self.nsamp))
            self._ch3_range_pending = False
        new_pos = (self.co % 2) if self.switching else 0
        if new_pos != self._relay_pos:
            if new_pos == 0:
                self.dvm.write('ROUT:OPEN (@211,248)')
                self.dvm.write('ROUT:CLOS (@218,241)')
            else:
                self.dvm.write('ROUT:OPEN (@218,241)')
                self.dvm.write('ROUT:CLOS (@211,248)')
            self._relay_pos = new_pos
            self._stop_event.wait(1.0)  # let relay contacts settle before driving signal
        if self.V1_setpoints is not None:
            self.V1 = self.V1_setpoints[self.co//2]
        else:
            ang = (self.co//2)/self.NDpts*2*np.pi
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
        self.dvm.write('INIT3 (@101:104)')
        self.dvm.write('*TRG')

    @pyqtSlot()
    def start(self):
        self._stop = False
        self._stop_event.clear()
        self.isidle = False
        self.rawN = CustomData.NPoints(self.fsig, self.fsamp, Nhars=self.Nhars, g1=self.g1, g2=self.g2,
                                       N=self.Npts, chunk_periods=self.chunk_periods)
        for self.co in range(self.Npts):
            if self._stop:
                break
            self.prepForMeas()
            if self._stop:          # stop pressed during SG setup — skip this point
                break
            self.sendtrig()
            self._stop_event.wait(1.0)  # interruptible — wakes immediately on stop()
            try:
                self.getvals()
            except Exception as e:
                self.logMessage.emit(f"Point {self.co} failed ({type(e).__name__}: {e}) — discarding set")
            if self._stop:
                break
        self.isidle = True
        if not self._stop and min(self.rawN.ats) > 0:
            self.rawN.calc()
            self.dataReady.emit(self.rawN)
        self.finished.emit()
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

    _CH3_RANGES = [0.3, 3.0, 10.0]

    def getvals(self):
        relay_pos = (self.co % 2) if self.switching else 0
        if relay_pos == 0:
            ch2 = self.dvm.query_binary_values('FETCH3? (@101)', datatype='f', is_big_endian=True)
            ch1 = self.dvm.query_binary_values('FETCH3? (@102)', datatype='f', is_big_endian=True)
        else:
            ch1 = self.dvm.query_binary_values('FETCH3? (@101)', datatype='f', is_big_endian=True)
            ch2 = self.dvm.query_binary_values('FETCH3? (@102)', datatype='f', is_big_endian=True)
        ch3 = self.dvm.query_binary_values('FETCH3? (@103)', datatype='f', is_big_endian=True)
        ch4 = self.dvm.query_binary_values('FETCH3? (@104)', datatype='f', is_big_endian=True)

        # Auto-range ch3: if overflowed, bump range so next point acquires at the larger range
        ch3_arr = np.asarray(ch3)
        if np.max(np.abs(ch3_arr)) > 0.9 * self.ch3_range:
            idx = self._CH3_RANGES.index(self.ch3_range)
            if idx < len(self._CH3_RANGES) - 1:
                self.ch3_range = self._CH3_RANGES[idx + 1]
                self._ch3_range_pending = True
                self.logMessage.emit(f'ch3 overflow -> bumping range to {self.ch3_range} V')
        if self.saverawdata:
            now = datetime.datetime.now()
            timestamp = now.strftime("%Y%m%d_%H%M%S")
            bd = os.path.join(self.bdraw, now.strftime("%Y"), now.strftime("%m"), now.strftime("%d"))
            os.makedirs(bd, exist_ok=True)
            fn = os.path.join(bd, f"ceramic_raw_{timestamp}.npz")
            np.savez_compressed(fn, ch1=ch1, ch2=ch2, ch3=ch3, ch4=ch4)
        self.rawN.setPoint(self.co, ch1, ch2, ch3, ch4, self.V1rb, self.V2rb, time.time())
        self.dataSetReady.emit(self.rawN.Data[self.co])
