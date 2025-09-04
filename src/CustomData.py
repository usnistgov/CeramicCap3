import R2FMath
import numpy as np

class MyData:
    def __init__(self, ave4,ts):
        self.ave4 = ave4
        self.ts = ts

    def __str__(self):
        return f"TS: {self.ts}"


class SampleData:
    def __init__(self,fsig,fsamp,data,Nhars=2):
        self.data  = np.array(data)
        self.fsig  = fsig
        self.fsamp = fsamp
        self.f0 = self.fsig/self.fsamp
        self.Nhars = Nhars

    def setfmin(self,fmin):
        self.fmin = fmin

    def findf(self):
        self.fmin = R2FMath.get_f(self.data,self.f0)
        return self.fmin
        
    def fit(self):
        self.Vc, self.fv,self.c2 = R2FMath.fit_sine_cplx(self.data,self.fmin,self.Nhars)


class FourChannels:
    def __init__(self,fsig,fsamp,Nhars,ch1,ch2,ch3,ch4,V1c,V2c,i,ts):
        self.ts    = ts
        self.fsig  = fsig
        self.fsamp = fsamp
        self.i=i
        if ts<0:
            return
        self.Data =[]
        self.V1c=V1c
        self.V2c=V2c
        self.Data.append(SampleData(fsig,fsamp,ch1,Nhars))
        self.Data.append(SampleData(fsig,fsamp,ch2,Nhars))
        self.Data.append(SampleData(fsig,fsamp,ch3,Nhars))
        self.Data.append(SampleData(fsig,fsamp,ch4,Nhars))
        fmin = self.Data[0].findf()
        for i in range(4):
            self.Data[i].setfmin(fmin)
            self.Data[i].fit()

class EightPoints:
    def __init__(self,fsig,fsamp,Nhars=2):
        self.fsig  = fsig
        self.fsamp = fsamp
        self.Nhars = Nhars
        self.ats = -1*np.ones(8)
        self.ts = min(self.ats)
        self.Data = np.zeros(8,dtype=object)

    def setPoint(self,i,ch1,ch2,ch3,ch4,V1c,V2c,ts):
        self.Data[i] = FourChannels(self.fsig,self.fsamp,self.Nhars,ch1,ch2,ch3,ch4,V1c,V2c,i,ts)
        self.ats[i] =ts
        if min(self.ats)>0:
            self.ts = np.mean(self.ats)


    def calc(self):
        self.raw8 = np.zeros((8,4),dtype=complex)
        self.ctrl = np.zeros((8,2),dtype=complex)
        for i in range(8):
            phi = np.angle(self.Data[i].Data[0].Vc)
            cf = np.exp(-1j*phi)
            for j in range(4):
                self.raw8[i,j] =  self.Data[i].Data[j].Vc*cf
            self.ctrl[i,0] =self.Data[i].V1c
            self.ctrl[i,1] =self.Data[i].V2c
        self.ave4  = 0.5*(self.raw8[::2,:]+self.raw8[1::2,:])
        self.ctrla = 0.5*(self.ctrl[::2,:]+self.ctrl[1::2,:])
        self.Cir = np.zeros(4,dtype=object)
        for i in range(4):
            self.Cir[i] = R2FMath.FourCplxPts(self.ave4[:,i])

        self.circlepar = np.hstack(( np.sqrt(self.Cir[2].circlepar[0]**2+self.Cir[2].circlepar[1]**2),self.Cir[2].circlepar[2],\
                                    np.sqrt(self.Cir[3].circlepar[0]**2+self.Cir[3].circlepar[1]**2),self.Cir[3].circlepar[2]))

        self.pars3,self.vals3,self.errs3,self.Chi3,self.Cov3 = R2FMath.mycomplexfit(self.ctrla[:,0],\
                             self.ave4[:,2])
        self.pars4,self.vals4,self.errs4,self.Chi4,self.Cov4 = R2FMath.mycomplexfit(self.ctrla[:,0],\
                             self.ave4[:,3])

        MCs = 20000
        parMC3=np.random.multivariate_normal( self.pars3,self.Cov3,size =MCs)
        parMC4=np.random.multivariate_normal( self.pars4,self.Cov4,size =MCs)

        result3=-(parMC3[:,0]+1j*parMC3[:,1])/(parMC3[:,2]+1j*parMC3[:,3])                
        result4=-(parMC4[:,0]+1j*parMC4[:,1])/(parMC4[:,2]+1j*parMC4[:,3])        
        self.Vz3 =  np.mean(result3)        
        self.Vz4 =  np.mean(result4)
        self.Vz3e =  0.5* (np.percentile(np.real(result3),84.134) -np.percentile(np.real(result3),15.866)) + \
                0.5j*(np.percentile(np.imag(result3),84.134) -np.percentile(np.imag(result3),15.866))
        self.Vz4e =  0.5* (np.percentile(np.real(result4),84.134) -np.percentile(np.real(result4),15.866)) + \
                0.5j*(np.percentile(np.imag(result4),84.134) -np.percentile(np.imag(result4),15.866))
        
        ratio =10
        self.V1m = self.ave4[:,0]
        self.V2m = self.ave4[:,1]
        self.V3m = self.ave4[:,2]
        self.V4m = self.ave4[:,3]
        self.u =self.V1m/self.V2m+ratio
        self.v3 = self.V3m/self.V2m
        self.v4 = self.V4m/self.V2m
        self.fp3, self.fv3,self.fe3,self.C23,self.Cov3 = R2FMath.mycomplexfit(self.v3,self.u)
        self.fp4, self.fv4,self.fe4,self.C24,self.Cov4 = R2FMath.mycomplexfit(self.v4,self.u)
        self.gain3  = self.fp3[2]+1j*self.fp3[3]
        self.gain4  = self.fp4[2]+1j*self.fp4[3]
        self.alpha3 = np.real(-self.u+ self.gain3*self.v3)/ratio
        self.beta3  = np.imag(-self.u+ self.gain3*self.v3)/ratio
        self.alpha4 = np.real(-self.u+ self.gain4*self.v4)/ratio
        self.beta4  = np.imag(-self.u+ self.gain4*self.v4)/ratio
        self.alphamean3 = np.mean(self.alpha3)
        self.betamean3 = np.mean(self.beta3)
        self.alphamean4 = np.mean(self.alpha4)
        self.betamean4 = np.mean(self.beta4)
        self.setGoodFlag()

    def setGoodFlag(self):
        self.goodData=True
        if self.circlepar[0]>1:
            self.goodData=False
        if self.circlepar[2]>5:
            self.goodData=False


        
class AllData():
    def __init__(self):
        self.mydict={}

    def append(self,ND:EightPoints):
        f = ND.fsig
        if f not in self.mydict:
            self.mydict[f] = []
        self.mydict[f].append(ND)
    
    def countf(self,f):
        if f not in self.mydict:
            return 0
        else:
            return len(self.mydict[f])
    
    def getCircles(self,f,t0):
        t=[]
        ret=[]
        for a in self.mydict[f]:
            t.append(a.ts-t0)
            ret.append(a.circlepar)
        return np.array(t),np.array(ret)
        
        


    def geta3(self,f,t0):
        t=[]
        ret=[]
        for a in self.mydict[f]:
            t.append(a.ts-t0)
            ret.append(a.alphamean3)
        return np.array(t),np.array(ret)

    def getb3(self,f,t0):
        t=[]
        ret=[]
        for a in self.mydict[f]:
            t.append(a.ts-t0)
            ret.append(a.betamean3)
        return np.array(t),np.array(ret)

    def geta4(self,f,t0):
        t=[]
        ret=[]
        for a in self.mydict[f]:
            t.append(a.ts-t0)
            ret.append(a.alphamean4)
        return np.array(t),np.array(ret)
    
    def getb4(self,f,t0):
        t=[]
        ret=[]
        for a in self.mydict[f]:
            t.append(a.ts-t0)
            ret.append(a.betamean4)
        return np.array(t),np.array(ret)

        line = np.hstack(t1,np.real(V3),np.real(V3e),np.imag(V3),np.image(V3e),np.real(V4),np.real(V4e),np.imag(V4),np.image(V4e))
        if f not in self.mydict:
            self.myficyt[f]= 1



        




    

