import serial
"""
Module for TZA amplifier
"""

class TZA:
    def __init__(self,comport):
        self.poldict ={0: 'not inverted',1:'inverted'}
        self.bwdict ={1: '10 kHz',2:'1 kHz',3: '100 Hz',4:'10 Hz'}
        self.gaindict ={1: 'x1',2:'x10',3: 'x100',4:'x1,000',5: 'x10,000',6:'x100,000'}
        self.gainfac ={1: 1,2:10,3:100,4:1000,5: 10000,6:100000}
        self.comport = comport
        self.is_port_open = False
        self.init_device()
        


    def open_port(self):
        if self.is_port_open == False:
            self.ser = serial.Serial(self.comport, baudrate=115200, bytesize=8, parity='N', stopbits=1, timeout=1000, xonxoff=0, rtscts=0)
            resp=self.get_resp(b'$U') #start command
            self.is_port_open=True

    def close_port(self):
        self.ser.close()
        self.is_port_open=False

    def get_settings(self):
        self.open_port()
        resp=self.get_resp(b'$F') #get polarity 0= not inverted 1 = inverted
        self.pol  = int(chr(resp[1]))
        self.hpol = self.poldict[self.pol]
        resp=self.get_resp(b'V?')  #Gain 1=x1, 2=x1e1,3 =x1e2, 4=x1e3,5=x1e4,6=x1e5
        self.gain = int(chr(resp[1]))
        self.hgain = self.gaindict[self.gain]
        resp=self.get_resp(b'B?')
        self.bw = int(chr(resp[1]))
        self.hbw = self.bwdict[self.bw]
        self.close_port()

    def set_gain(self,gain):
        if gain>=1 and gain<=6:
            self.open_port()
            self.get_resp(f"V{gain}".encode('utf-8'))
            self.get_settings()
            self.close_port()

    def set_hgain(self,hgain):
        for k,v in self.gaindict.items():
            if v == hgain:
                self.set_gain(k)

    def set_fgain(self,fgain):
        for k,v in self.gainfac.items():
            if v == fgain:
                self.set_gain(k)

    def text_to_gain(self,txt):
        for k,v in self.gaindict.items():
            if v == txt:
                return (self.gainfac[k])


    def set_bw(self,bw):
        if bw>=1 and bw<=4:
            self.open_port()
            self.get_resp(f"B{bw}".encode('utf-8'))
            self.get_settings()
            self.close_port()


    def init_device(self):
        self.open_port()
        resp=self.get_resp(b'$N') #set to non inverter
        resp=self.get_resp(b'B1') #ser bandwith 1,2,3,4 1=10 kHz, 2- 1 kHz, 3 = 100 Hz, 4 = 10Hz
        self.get_settings()
        self.close_port()

    def get_resp(self,cmd):
        self.ser.write(cmd)    
        ret =self.ser.read_until(b'\r')#
        return ret



    
        