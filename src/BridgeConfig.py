import numpy as np
import os

class BridgeConfig():
    def __init__(self):
        self.cfg ={'TOPLEFT':'1000pF','BOTTOMLEFT':'100pF','TOPRIGHT':'1000pF','BOTTOMRIGHT':'100pF'}
        self.dict_100nF_10nF = self.readfile('SN2519J00896_100n_10n.dat')

    def readfile(self,fn):
        bd =r"C:\python\CeramicCap2" #os.getcwd()
        fp=os.path.join(bd,'data',fn)
        mydict ={}
        with open(fp,'r') as file:
            for lines in file:
                flines =[ float(i) for i in lines.split()]
                mydict[flines[0]]=flines[1]+1j*flines[2]
        return mydict


    def getC41(self):
        k ='TOPRIGHT'
        if    self.cfg[k]=='100pF':
            return 1e-8
        elif self.cfg[k]=='1000pF':
            return 1e-9
        elif self.cfg[k]=='10nF':
            return 1e-8
        elif self.cfg[k]=='100nF':
            return 1e-7
        elif self.cfg[k]=='1uF':
            return 1e-6
        elif self.cfg[k]=='10uF':
            return 1e-5
        else:   raise Exception( f"cfg[{k}] invalid")

    def getC42(self):
        k ='BOTTOMRIGHT'
        if    self.cfg[k]=='100pF':
            return 1e-8
        elif self.cfg[k]=='1000pF':
            return 1e-9
        elif self.cfg[k]=='10nF':
            return 1e-8
        elif self.cfg[k]=='100nF':
            return 1e-7
        elif self.cfg[k]=='1uF':
            return 1e-6
        elif self.cfg[k]=='10uF':
            return 1e-5
        else:   raise Exception( f"cfg[{k}] invalid")


    def calcI2(self,f, R=50,V2=-9.0):
        C42 =self.getC42()
        iw = 1j*f*2*np.pi
        I2 = V2*iw*C42/(1+iw*C42*R+1)
        return I2

    def calcVsmall(self,f,R=50,V2=-9.0):
        C42 =self.getC42()
        C41 =self.getC41()
        iw = 1j*f*2*np.pi
        I2 = V2*iw*C42/(1+iw*C42*R+1)
        V1=-I2*(R+1/(iw*C41))
        return V1

    def setCfg(self,key,value):
        self.cfg[key]=value
        self.printCfg()

    def printCfg(self):
        for k,v in self.cfg.items():
            print(f"{k}: {v}")


    def getgain1(self,f):
        if self.cfg['TOPRIGHT']=='1000pF' and self.cfg['BOTTOMRIGHT']=='100pF':
            if f<30000: return 'x1,000'
            else: return 'x100'
        elif self.cfg['TOPRIGHT']=='100nF' and self.cfg['BOTTOMRIGHT']=='10nF':
            if f<30000: return 'x100'
            else: return 'x100'
        else: return 'x100'
           
    def getgain2(self,f):
        if self.cfg['TOPRIGHT']=='1000pF' and self.cfg['BOTTOMRIGHT']=='100pF':
            if f<30000: return 'x1,000'
            else: return 'x100'
        elif self.cfg['TOPRIGHT']=='10nF' and self.cfg['BOTTOMRIGHT']=='1000pF':
            if f<30000: return 'x100'
            else: return 'x10'
        elif self.cfg['TOPRIGHT']=='100nF' and self.cfg['BOTTOMRIGHT']=='10nF':
            if f<30000: return 'x1'
            else: return 'x1'
        elif self.cfg['TOPRIGHT']=='1uF' and self.cfg['BOTTOMRIGHT']=='100nF':
            if f<30000: return 'x1'
            else: return 'x1'
        elif self.cfg['TOPRIGHT']=='10uF' and self.cfg['BOTTOMRIGHT']=='1uF':
            if f<30000: return 'x1'
            else: return 'x1'

    def getdV(self,f):
        if self.cfg['TOPRIGHT']=='1000pF' and self.cfg['BOTTOMRIGHT']=='100pF':
            if f<30000: return 0.01
            else: return 0.01
        elif self.cfg['TOPRIGHT']=='10nF' and self.cfg['BOTTOMRIGHT']=='1000pF':
            if f<30000: return 0.01
            else: return 0.01
        elif self.cfg['TOPRIGHT']=='100nF' and self.cfg['BOTTOMRIGHT']=='10nF':
            if f<30000: return 0.01
            else: return 0.01
        elif self.cfg['TOPRIGHT']=='1uF' and self.cfg['BOTTOMRIGHT']=='100nF':
            if f<30000: return 0.01
            else: return 0.01
        elif self.cfg['TOPRIGHT']=='10uF' and self.cfg['BOTTOMRIGHT']=='1uF':
            Vsmall = self.calcVsmall(f)
            return 0.02*np.abs(Vsmall)

          

    def obtainfromdict(self,f,mydict):
        fl = list(mydict)
        ix = np.argmin(np.abs(np.array(fl)-f))
        ret = mydict[fl[ix]]
        print(f'obtainfromdict {f=} {ret=}')
        return ret

    def getV2(self,f):
        return -9.0+1j*0
        if self.cfg['TOPRIGHT']=='1000pF' and self.cfg['BOTTOMRIGHT']=='100pF':
            if f<30000:
                return -9.9+1j*0
            else:
                return -9.9+1j*0
        elif self.cfg['TOPRIGHT']=='10uF' and self.cfg['BOTTOMRIGHT']=='1uF':
            if f<=30000:
                return -0.99+1j*0
            elif f==100000:
                return -0.099+1j*0
            else:
                return -0.99+1j*0
        else:
            if f<30000:
                return -9.9+1j*0
            else:
                return -0.99+1j*0


    def getV1(self,f):
        return self.calcVsmall(f)
        if self.cfg['TOPRIGHT']=='1000pF' and self.cfg['BOTTOMRIGHT']=='100pF':
            if f<30000: return 0.99+1j*0.01
            else:  return 0.99+1j*0.01
        elif self.cfg['TOPRIGHT']=='10nF' and self.cfg['BOTTOMRIGHT']=='1000pF':
            if f<30000: return 0.99+1j*0.01
            else:  return 0.09+1j*0.001
        elif self.cfg['TOPRIGHT']=='100nF' and self.cfg['BOTTOMRIGHT']=='10nF':
            return self.obtainfromdict(f,self.dict_100nF_10nF)
        elif self.cfg['TOPRIGHT']=='1uF' and self.cfg['BOTTOMRIGHT']=='100nF': 
            myV1dict ={1000:1+0.26j,3000:1+0.26j,10000:2+2.5j,17000:3+3.7j,30000: 0.53+0.407j,52000:0.725+0.36j,100000: 0.99+0.1j}
            return self.obtainfromdict(f,myV1dict)
        elif self.cfg['TOPRIGHT']=='10uF' and self.cfg['BOTTOMRIGHT']=='1uF': 
            myV1dict ={1000: 0.9338+0.138j,3000:0.9338+0.138j,10000:0.9338+0.138j,17000:0.9736+0.158j,\
                       30000:1.00+0.09j,52000:0.99+0.04j,100000:0.1+0.005*1j}
            return self.obtainfromdict(f,myV1dict)