import R2FMath
import numpy as np
import copy

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
    def __init__(self,fsig,fsamp,Nhars=2,g1=1,g2=1,ratio=10):
        self.Res={}
        self.Res['fsig']= fsig
        self.Res['ratio']= ratio
        self.Res['fsamp']= fsamp
        self.Res['Nhars']= Nhars
        self.Res['gain1']= g1
        self.Res['gain2']= g2      
        self.ats = -1*np.ones(8)
        self.Res['ts']= min(self.ats)      
        self.Data = np.zeros(8,dtype=object)

    def setPoint(self,i,ch1,ch2,ch3,ch4,V1c,V2c,ts):
        self.V2c=V2c
        self.Data[i] = FourChannels(self.Res['fsig'],self.Res['fsamp'],self.Res['Nhars'],\
                                    ch1,ch2,ch3,ch4,V1c,V2c,i,ts)
        self.ats[i] =ts
        if min(self.ats)>0:
            self.Res['ts'] = np.mean(self.ats)


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
        self.Ccirpar =[]
        la=['V1cplxcenter','V1cplxradius','V2cplxcenter','V2cplxradius',\
            'V3cplxcenter','V3cplxradius','V4cplxcenter','V4cplxradius']
        for i in range(4):
            self.Cir[i] = R2FMath.FourCplxPts(self.ave4[:,i])
            if i>0:
                self.Res[la[i*2]]= self.Cir[i].cplxcenter
                self.Res[la[i*2+1]]= self.Cir[i].cplxradius
        self.Res['V2setamp']= self.V2c
        self.V1fit= R2FMath.FourCplxPts(self.ctrla[:,0])        
        self.Res['V1setcplxcenter']= self.V1fit.cplxcenter
        self.Res['V1setcplxradius']= self.V1fit.cplxradius

        self.pars3,self.vals3,self.errs3,self.Chi3,self.Cov3 = R2FMath.mycomplexfit(self.ctrla[:,0],\
                             self.ave4[:,2])
        self.pars4,self.vals4,self.errs4,self.Chi4,self.Cov4 = R2FMath.mycomplexfit(self.ctrla[:,0],\
                             self.ave4[:,3])
        result3=-(self.pars3[0]+1j*self.pars3[1])/(self.pars3[2]+1j*self.pars3[3])                
        result4=-(self.pars4[0]+1j*self.pars4[1])/(self.pars4[2]+1j*self.pars4[3])        
        self.Res['Vz3']= result3
        self.Res['Vz4']= result4

        self.V1m = self.ave4[:,0]
        self.V2m = self.ave4[:,1]
        self.V3m = self.ave4[:,2]
        self.V4m = self.ave4[:,3]
        self.u =self.V1m/self.V2m+self.Res['ratio']
        self.v3 = self.V3m/self.V2m
        self.v4 = self.V4m/self.V2m
        self.fp3, self.fv3,self.fe3,self.C23,self.Cov3 = R2FMath.mycomplexfit(self.v3,self.u)
        self.fp4, self.fv4,self.fe4,self.C24,self.Cov4 = R2FMath.mycomplexfit(self.v4,self.u)
        self.gain3  = self.fp3[2]+1j*self.fp3[3]
        self.gain4  = self.fp4[2]+1j*self.fp4[3]
        self.alpha3 = np.real(-self.u+ self.gain3*self.v3)/self.Res['ratio']
        self.beta3  = np.imag(-self.u+ self.gain3*self.v3)/self.Res['ratio']
        self.alpha4 = np.real(-self.u+ self.gain4*self.v4)/self.Res['ratio']
        self.beta4  = np.imag(-self.u+ self.gain4*self.v4)/self.Res['ratio']
        self.Res['alpha3mean'] = np.mean(self.alpha3)
        self.Res['beta3mean']  = np.mean(self.beta3)
        self.Res['alpha4mean'] = np.mean(self.alpha4)
        self.Res['beta4mean']  = np.mean(self.beta4)
        self.setGoodFlag()

    def setGoodFlag(self):
        self.goodData=True
        if abs(self.Res['V3cplxcenter'])>1:
            self.goodData=False
        if abs(self.Res['V4cplxcenter'])>1:
            self.goodData=False
        
class AllData():
    def __init__(self):
        self.mydict={}

    def append(self,ND:EightPoints):
        f = ND.Res['fsig']
        if f not in self.mydict:
            self.mydict[f] = []
        self.mydict[f].append(ND)
    
    def deletekey(self,f):
        if f in self.mydict:
             del self.mydict[f]

    def count(self):
        return len(list(self.mydict))

    def countf(self,f):
        if f not in self.mydict:
            return 0
        else:
            return len(self.mydict[f])
    
    def getkeys(self,f,keys):
        retdict ={}
        L = len(self.mydict[f])
        if L==0:
            for k in keys:
                retdict[k]=np.array([])
                return retdict
        for k in keys:
            retdict[k]=np.empty(L,dtype=np.array(self.mydict[f][0]).dtype)
        for n,a in enumerate(self.mydict[f]):
            for k in keys:
                retdict[k][n] =a.Res[k]
        return retdict

    def getallkeys(self,f):
        k = list(self.mydict)
        return self.getkeys(f,k)

    def getAveVolts(self,f,t0=0):
        L = len(self.mydict[f])
        Nrows = np.shape(self.mydict[f][0].ave4)[0]  
        Ncols = 2+2*np.shape(self.mydict[f][0].ave4)[1]   +2*np.shape(self.mydict[f][0].ctrla)[1]
        #ret = np.empty(Nrows*L,Ncols)
        ret =[]
        for i in range(L):
            for  obj in self.mydict[f]:
                for j in range(Nrows):
                    line  = np.hstack((f,self.mydict[f][0].Res['ts']-t0))
                    for k in range(np.shape(obj.ave4)[1]):
                        line =np.hstack((line,obj.ave4[j,k].real,obj.ave4[j,k].imag))
                    for k in range(np.shape(obj.ctrla)[1]):
                        line =np.hstack((line,obj.ctrla[j,k].real,obj.ctrla[j,k].imag))
                ret.append(line)            
        ret = np.array(ret)
        return ret



    def getabf(self):
        f=list(self.mydict)
        ret=[]
        for ff in f:
            dummy=[]
            for a in self.mydict[ff]:
                dummy.append(np.hstack( (a.alphamean3,a.betamean3,a.alphamean4,a.betamean4)))
            dummy = np.array(dummy)
            da = dummy[:,2]-dummy[:,0]
            db = dummy[:,3]-dummy[:,1]
            if len(da)>2:
                sa =np.std(da,ddof=1)
            else:
                sa=0
            if len(db)>2:
                sb =np.std(da,ddof=1)
            else:
                sb=0
            ret.append(np.hstack((np.mean(da),np.mean(db),sa,sb)))

        return np.array(f),np.array(ret)

        




    

