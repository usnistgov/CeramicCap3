import R2FMath
import numpy as np
import copy
import scipy.optimize



def fit_sine_cplx(y, fsamp,fsig,fline=60, Nhars=1,use_hann=True):
    """
    Fits a sine wave with a DC offset to the data: 
    y = DC + cos_coeff*cos(wt) + sin_coeff*sin(wt)
    
    Parameters:
        y (array-like): The signal data.
        rf (float): Relative frequency (f0/fs).
        use_hann (bool): Whether to apply a Hanning window.
        
    Returns:
        tuple: (complex_amplitude, fitted_values, error_variance, residual_sum_of_squares)
    """
    y = np.asarray(y)
    n = len(y)
    
    rf = fsig/fsamp
    rlf = fline/fsamp
    # Pre-calculate angular frequency vector
    wt  =  2 * np.pi * np.arange(n) * rf
    wlf =  2 * np.pi * np.arange(n) * rlf

    
    # Build design matrix X directly using column_stack
    X = np.column_stack((np.ones(n), np.cos(wt), np.sin(wt),np.cos(wlf), np.sin(wlf)))
    cuhars=1
    while cuhars<Nhars:
        cuhars+=1
        X = np.column_stack((X, np.cos(cuhars*wt), np.sin(cuhars*wt)))

    if use_hann:
        # Generate window and apply weights
        w = np.hanning(n)
        # Broadcasting the window across X columns is faster and cleaner
        X_w = X * w[:, np.newaxis] 
        y_w = y * w
        
        # lstsq is numerically far more stable than solving normal equations (X^T * X)
        fit_pars, _, _, _ = np.linalg.lstsq(X_w, y_w, rcond=None)
    else:
        fit_pars, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        
    # Calculate fitted values and residuals
    fit_vals = X @ fit_pars
    residuals = y - fit_vals
    
    # Residual Sum of Squares (C2)
    rss = np.sum(residuals**2) 
    
    ndf = n - 3
    errv = np.sqrt(rss / ndf) if ndf > 0 else 0.0
    
    # Complex amplitude: cos_coeff - 1j * sin_coeff
    complex_amp = fit_pars[1] - 1j * fit_pars[2]
    
    return complex_amp, fit_vals, errv, rss


def get_f(y, fsamp, fsig_guess, fline_guess=60.0):
    n = len(y)
    
    def cost_function(params):
        fsig, fline = params
        return fit_sine_cplx(y, fsamp, fsig, fline)[3]
        
    initial_guess = [fsig_guess, fline_guess]
    
    # Define the tight search window using your n/(n±1) logic
    fsig_min = fsig_guess * n / (n + 1)
    fsig_max = fsig_guess * n / (n - 1)
    
    bounds = (
        (fsig_min, fsig_max),  # Strict cage for the main signal
        (59.5, 60.5)           # Strict cage for the line noise
    )
    
    result = scipy.optimize.minimize(cost_function, initial_guess, method='L-BFGS-B', bounds=bounds)
    
    if result.success:
        return result.x[0], result.x[1]
    else:
        print(f"Optimization failed: {result.message}")
        return fsig_guess, fline_guess
    


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
        self.fline = 60
        self.f0 = self.fsig/self.fsamp
        self.Nhars = Nhars

    def setf(self,fsig,fline):
        self.fsig,self.fline =  fsig,fline

    def findf(self):
        self.fsig,self.fline = get_f(self.data, self.fsamp, self.fsig, self.fline)
        return self.fsig,self.fline
        
    def fit(self):  #(y, fsamp,fsig,fline=60, use_hann=True):
        self.Vc, self.fv,errv,self.c2 = fit_sine_cplx(self.data,self.fsamp,self.fsig,self.fline,self.Nhars,True)


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
        fsig, fline =  self.Data[0].findf()
        for i in range(4):
            self.Data[i].setf(fsig, fline)
            self.Data[i].fit()

class NPoints:
    def __init__(self,fsig,fsamp,Nhars=2,g1=1,g2=1,ratio=10,N=8):
        self.N=N
        self.Res={}
        self.Res['fsig']= fsig
        self.Res['ratio']= ratio
        self.Res['fsamp']= fsamp
        self.Res['Nhars']= Nhars
        self.Res['gain1']= g1
        self.Res['gain2']= g2      
        self.ats = -1*np.ones(N)
        self.Res['ts']= min(self.ats)      
        self.Data = np.zeros(N,dtype=object)

    def setPoint(self,i,ch1,ch2,ch3,ch4,V1c,V2c,ts):
        self.V1c=V1c
        self.Data[i] = FourChannels(self.Res['fsig'],self.Res['fsamp'],self.Res['Nhars'],\
                                    ch1,ch2,ch3,ch4,V1c,V2c,i,ts)
        self.ats[i] =ts
        if min(self.ats)>0:
            self.Res['ts'] = np.mean(self.ats)

    def precalc(self):
        self.raw8 = np.zeros((self.N,4),dtype=complex)
        self.ctrl = np.zeros((self.N,2),dtype=complex)
        for i in range(self.N):
            phi = np.angle(self.Data[i].Data[0].Vc)
            cf = np.exp(-1j*phi)
            for j in range(4):
                self.raw8[i,j] =  self.Data[i].Data[j].Vc*cf
            self.ctrl[i,0] =self.Data[i].V1c
            self.ctrl[i,1] =self.Data[i].V2c
        self.ave4  = 0.5*(self.raw8[::2,:]+self.raw8[1::2,:])
        self.ctrla = 0.5*(self.ctrl[::2,:]+self.ctrl[1::2,:])

    def calc(self):
        self.precalc()
        self.Cir = np.zeros(self.N//2,dtype=object)
        self.Ccirpar =[]
        la=['V1cplxcenter','V1cplxradius','V1cplxangle','V2cplxcenter','V2cplxradius','V2cplxangle',\
            'V3cplxcenter','V3cplxradius','V3cplxangle','V4cplxcenter','V4cplxradius','V4cplxangle']
        for i in range(4):
            if i!=0:
                self.Cir[i] = R2FMath.FourPlusCplxPts(self.ave4[:,i])
                #R2FMath.FourCplxPts(self.ave4[:,i])
                self.Res[la[i*3]]= self.Cir[i].cplxcenter
                self.Res[la[i*3+1]]= self.Cir[i].cplxradius
                self.Res[la[i*3+2]]= self.Cir[i].cplxangle
            else:
                self.Res[la[i*3]]= np.mean(self.ave4[:,i])
                self.Res[la[i*3+1]]= 1e-20
                self.Res[la[i*3+2]]= 0

        self.Res['V1setamp']= self.V1c
        self.V2fit= R2FMath.FourPlusCplxPts(self.ctrla[:,1])        
        self.Res['V2setcplxcenter']= self.V2fit.cplxcenter
        self.Res['V2setcplxradius']= self.V2fit.cplxradius
        self.Res['V2setcplxangle']= self.V2fit.cplxangle


        self.pars3,self.vals3,self.errs3,self.Chi3,self.Cov3 = R2FMath.mycomplexfit(self.ctrla[:,1],\
                             self.ave4[:,2])
        self.pars4,self.vals4,self.errs4,self.Chi4,self.Cov4 = R2FMath.mycomplexfit(self.ctrla[:,1],\
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

    def append(self,ND:NPoints):
        f = ND.Res['fsig']
        if f not in self.mydict:
            self.mydict[f] = []
        self.mydict[f].append(ND)
    
    def deletekey(self,f):
        if f in self.mydict:
             del self.mydict[f]

    def count(self):
        return len(self.mydict)

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
            retdict[k]=np.empty(L,dtype=object)
        for n,a in enumerate(self.mydict[f]):
            for k in keys:
                retdict[k][n] =a.Res[k]
        return retdict

    def getallkeys(self,f):
        k = list(self.mydict)
        return self.getkeys(f,k)
    
    def getdictf(self, keys):
        retdict={}
        keys = list(keys)
        if 'fsig' not in keys:
            keys.append('fsig')
        allf = list(self.mydict)
        if len(allf)==0:
            for k in keys:
                retdict[k]=np.array([])
                return retdict
        else:
            for k in keys:
                retdict[k]=[]
            for f in allf:
                odict = self.getkeys(f,keys)
                for k in keys:
                    for item in odict[k]:
                        retdict[k].append(item)
            for k in keys:
                retdict[k] =np.array(retdict[k])
            return retdict


    def getAveVolts(self,f,t0=0):
        L = len(self.mydict[f])
        Nrows = np.shape(self.mydict[f][0].ave4)[0]  
        Ncols = 2+2*np.shape(self.mydict[f][0].ave4)[1]   +2*np.shape(self.mydict[f][0].ctrla)[1]
        #ret = np.empty(Nrows*L,Ncols)
        ret =[]
        
        obj = self.mydict[f][-1]
        Nrows = np.shape(obj.ave4)[0]
        for j in range(Nrows):
            line  = np.hstack((f,obj.Res['ts']-t0))
            for k in range(np.shape(obj.ave4)[1]):
                line =np.hstack((line,obj.ave4[j,k].real,obj.ave4[j,k].imag))
            for k in range(np.shape(obj.ctrla)[1]):
                line =np.hstack((line,obj.ctrla[j,k].real,obj.ctrla[j,k].imag))
            ret.append(line)            
        ret = np.array(ret)
        return ret


    def getRawPhasors(self,f,t0=0):
        L = len(self.mydict[f])
        Nrows = np.shape(self.mydict[f][0].raw8)[0]  
        Ncols = 2+2*np.shape(self.mydict[f][0].raw8)[1]   +2*np.shape(self.mydict[f][0].ctrla)[1]
        #ret = np.empty(Nrows*L,Ncols)
        ret =[]
        
        obj = self.mydict[f][-1]
        Nrows = np.shape(obj.raw8)[0]
        for j in range(Nrows):
            line  = np.hstack((f,obj.Res['ts']-t0))
            for k in range(np.shape(obj.raw8)[1]):
                line =np.hstack((line,obj.raw8[j,k].real,obj.raw8[j,k].imag))
            for k in range(np.shape(obj.ctrla)[1]):
                line =np.hstack((line,obj.ctrla[j//2,k].real,obj.ctrla[j//2,k].imag))
            ret.append(line)            
        ret = np.array(ret)
        return ret



    def getabf(self):
        f=list(self.mydict)
        ret=[]
        for ff in f:
            dummy=[]
            for a in self.mydict[ff]:
                dummy.append(np.hstack((a.Res['alpha3mean'], a.Res['beta3mean'], a.Res['alpha4mean'], a.Res['beta4mean'])))
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

        




    

