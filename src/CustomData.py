import numpy as np
import scipy.optimize



def fit_sine_cplx(y, fsamp, fsig, fline=60, Nhars=1, use_hann=True):
    """
    Fits a sine wave with a DC offset to the data:
    y = DC + cos_coeff*cos(wt) + sin_coeff*sin(wt) + [line harmonics] + [signal harmonics]

    Parameters:
        y (array-like): The signal data.
        fsamp (float): Sampling frequency.
        fsig (float): Signal frequency.
        fline (float): Line noise frequency.
        Nhars (int): Number of signal harmonics to fit (1 = fundamental only).
        use_hann (bool): Whether to apply a Hanning window consistently throughout.

    Returns:
        tuple: (complex_amplitude, fit_vals, errv, rss)
            complex_amplitude: cos_coeff - 1j * sin_coeff for the fundamental
            fit_vals: fitted values in the (possibly windowed) domain
            errv: RMS residual (sqrt of RSS / ndf)
            rss: residual sum of squares, in whichever domain was fitted
    """
    y = np.asarray(y, dtype=float)
    n = len(y)

    rf  = fsig  / fsamp
    rlf = fline / fsamp

    wt  = 2 * np.pi * np.arange(n) * rf
    wlf = 2 * np.pi * np.arange(n) * rlf

    cols = [np.ones(n), np.cos(wt), np.sin(wt), np.cos(wlf), np.sin(wlf), np.cos(2*wlf), np.sin(2*wlf)]
    for h in range(2, Nhars + 1):
        cols.extend([np.cos(h * wt), np.sin(h * wt)])
    X = np.column_stack(cols)

    n_params = 1 + 2 * Nhars + 4  # DC + signal harmonics + line fundamental + line 2nd harmonic
    ndf = n - n_params

    if use_hann:
        w = np.hanning(n)
        X_eff = X * w[:, np.newaxis]
        y_eff = y * w
    else:
        X_eff = X
        y_eff = y

    fit_pars, lstsq_res, _, _ = np.linalg.lstsq(X_eff, y_eff, rcond=None)

    fit_vals = X @ fit_pars
    rss = float(lstsq_res[0]) if len(lstsq_res) == 1 else float(np.sum((y_eff - X_eff @ fit_pars) ** 2))
    errv     = np.sqrt(rss / ndf) if ndf > 0 else 0.0

    complex_amp = fit_pars[1] - 1j * fit_pars[2]

    return complex_amp, fit_vals, errv, rss


def get_f(y, fsamp, fsig_guess, fline_guess=60.0, use_hann=True, Nhars=1):
    """
    Estimates the signal frequency by minimizing the residual sum of squares
    from fit_sine_cplx using Brent's bounded method (guaranteed convergence).

    Parameters:
        y (array-like): The signal data.
        fsamp (float): Sampling frequency.
        fsig_guess (float): Initial guess for the signal frequency.
        fline_guess (float): Line noise frequency (fixed during search).
        use_hann (bool): Passed through to fit_sine_cplx.
        Nhars (int): Number of signal harmonics — must match the final fit.

    Returns:
        tuple: (fsig, fline) best-fit frequencies.
    """
    y = np.asarray(y, dtype=float)
    n = len(y)

    fsig_min = fsig_guess * n / (n + 1)
    fsig_max = fsig_guess * n / (n - 1)

    res_sig = scipy.optimize.minimize_scalar(
        lambda fsig: fit_sine_cplx(y, fsamp, fsig, fline_guess, Nhars=Nhars, use_hann=use_hann)[3],
        bounds=(fsig_min, fsig_max),
        method='bounded',
    )
    best_fsig = res_sig.x

    res_line = scipy.optimize.minimize_scalar(
        lambda fline: fit_sine_cplx(y, fsamp, best_fsig, fline, Nhars=Nhars, use_hann=use_hann)[3],
        bounds=(fline_guess - 0.5, fline_guess + 0.5),
        method='bounded',
    )
    return best_fsig, res_line.x
        


class SampleData:
    def __init__(self,fsig,fsamp,data,Nhars=1):
        self.data  = np.array(data)
        self.fsig  = fsig
        self.fsamp = fsamp
        self.fline = 60
        self.Nhars = Nhars

    def setf(self,fsig,fline):
        self.fsig,self.fline =  fsig,fline

    def findf(self):
        self.fsig, self.fline = get_f(self.data, self.fsamp, self.fsig, self.fline, Nhars=self.Nhars)
        return self.fsig, self.fline

    def fit(self):  #(y, fsamp,fsig,fline=60, use_hann=True):
        self.Vc, self.fv,errv,self.c2 = fit_sine_cplx(self.data,self.fsamp,self.fsig,self.fline,self.Nhars,True)

    def strip_raw(self):
        self.data = None
        self.fv = None


class FourChannels:
    def __init__(self,fsig,fsamp,Nhars,ch1,ch2,ch3,ch4,V1c,V2c,i,ts):
        self.ts    = ts
        self.fsig  = fsig
        self.fsamp = fsamp
        self.i     = i
        self.Data  = []
        self.V1c   = V1c
        self.V2c   = V2c
        if ts<0:
            return
        self.Data.append(SampleData(fsig,fsamp,ch1,Nhars))
        self.Data.append(SampleData(fsig,fsamp,ch2,Nhars))
        self.Data.append(SampleData(fsig,fsamp,ch3,Nhars))
        self.Data.append(SampleData(fsig,fsamp,ch4,Nhars))
        fsig, fline =  self.Data[0].findf()
        for i in range(4):
            self.Data[i].setf(fsig, fline)
            self.Data[i].fit()

    def strip_raw(self):
        for sd in self.Data:
            sd.strip_raw()

class NPoints:
    def __init__(self,fsig,fsamp,Nhars=1,g1=1,g2=1,ratio=10,N=8):
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
        self.V1m = self.ave4[:, 0]
        self.V2m = self.ave4[:, 1]
        self.V3m = self.ave4[:, 2]
        self.V4m = self.ave4[:, 3]
        self.eta2 = self.V2m / self.V1m
        self.eta3 = self.V3m / self.V1m
        self.eta4 = self.V4m / self.V1m

        Xmat = np.column_stack([self.eta2, np.ones(len(self.eta2))])
        (m3, c3), _, _, _ = np.linalg.lstsq(Xmat, self.eta3, rcond=None)
        (m4, c4), _, _, _ = np.linalg.lstsq(Xmat, self.eta4, rcond=None)
        self.gamma3 = 1.0 / m3
        self.gamma4 = 1.0 / m4
        self.Res['Vz3'] = c3 / m3
        self.Res['Vz4'] = c4 / m4
        self.Res['gamma3'] = self.gamma3
        self.Res['gamma4'] = self.gamma4
        N2 = self.N // 2
        V4_raw = np.array([0.5*(self.Data[2*k].Data[3].Vc + self.Data[2*k+1].Data[3].Vc)
                           for k in range(N2)])
        V1rb_pts = self.ctrla[:, 0]
        Xmat = np.column_stack([V1rb_pts, np.ones(N2)])
        (m_bal, c_bal), _, _, _ = np.linalg.lstsq(Xmat, V4_raw, rcond=None)
        self.Res['V4fit_slope'] = m_bal
        self.Res['V4fit_intercept'] = c_bal
        self.Res['V1_balance'] = -c_bal / m_bal

        self.combined3 = self.gamma3 * self.eta3 - self.eta2
        self.combined4 = self.gamma4 * self.eta4 - self.eta2

        ratio = self.Res['ratio']
        self.alpha3 = np.real(self.combined3) / ratio - 1
        self.beta3  = np.imag(self.combined3) / ratio
        self.alpha4 = np.real(self.combined4) / ratio - 1
        self.beta4  = np.imag(self.combined4) / ratio
        self.Res['alpha3mean'] = np.mean(self.alpha3)
        self.Res['beta3mean']  = np.mean(self.beta3)
        self.Res['alpha4mean'] = np.mean(self.alpha4)
        self.Res['beta4mean']  = np.mean(self.beta4)
        self.Res['V1cReadback'] = self.V1c

        self.setGoodFlag()

    def setGoodFlag(self):
        self.goodData = not (
            np.any(~np.isfinite(np.abs(self.combined3))) or
            np.any(~np.isfinite(np.abs(self.combined4)))
        )

    def strip_raw(self):
        for fc in self.Data:
            if isinstance(fc, FourChannels):
                fc.strip_raw()


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
        if f not in self.mydict or len(self.mydict[f]) == 0:
            return {}
        keys = list(self.mydict[f][0].Res.keys())
        return self.getkeys(f, keys)
    
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


    def getRawPhasors(self,f,t0=0):
        ret = []
        for obj in self.mydict[f]:
            Nrows = np.shape(obj.raw8)[0]
            for j in range(Nrows):
                line = np.hstack((f, obj.Res['ts']-t0))
                for k in range(np.shape(obj.raw8)[1]):
                    line = np.hstack((line, obj.raw8[j,k].real, obj.raw8[j,k].imag))
                for k in range(np.shape(obj.ctrla)[1]):
                    line = np.hstack((line, obj.ctrla[j//2,k].real, obj.ctrla[j//2,k].imag))
                line = np.hstack((line, obj.Res['gain1'], obj.Res['gain2']))
                ret.append(line)
        return np.array(ret)





    

