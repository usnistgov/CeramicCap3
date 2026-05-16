import numpy as np
import sys
import os
sys.path.append('..')
# Now Python can see the 'src' folder!
import src.R2FMath as R2FMath

class oneCap:
    def __init__(self, bd, fn):
        data = np.loadtxt(os.path.join(bd,fn))
        self.ave = 0.5*(data[0:-1:2,:]+data[1::2,:])
        self.ana=[]
        self.meanV1 =[]
        self.elliV2 =[]
        self.elliV3 =[]
        self.elliV4 =[]
        self.elli2=[]
        self.elli3=[]
        self.elli4=[]
        self.allf = []
        self.aeta2 =[]
        self.aeta3 =[]
        self.aeta4 =[]
        for i in range(0, len(self.ave), 8):
            block = self.ave[i : i + 8,:]
            if len(block) < 8:
                print(f"Processing final partial block of size {len(block)}")
                break
            f = np.mean(block[:,0])
            self.allf.append(f)
            V1 = (block[:,2]+1j*block[:,3])
            V2 = (block[:,4]+1j*block[:,5])
            V3 = (block[:,6]+1j*block[:,7])
            V4 = (block[:,8]+1j*block[:,9])

            self.meanV1.append(np.mean(V1))
            self.elliV2.append(R2FMath.ComplexEllipse.fit_from_cmplx_points(V2))
            self.elliV3.append(R2FMath.ComplexEllipse.fit_from_cmplx_points(V3))
            self.elliV4.append(R2FMath.ComplexEllipse.fit_from_cmplx_points(V4))

            eta2 =  1000*V2/V1
            eta3 =  1000*V3/V1
            eta4 =  1000*V4/V1
            self.aeta2.append(eta2)
            self.aeta3.append(eta3)
            self.aeta4.append(eta4)

            elli2 = fitted_ellipse = R2FMath.ComplexEllipse.fit_from_cmplx_points(eta2)
            elli3 = fitted_ellipse = R2FMath.ComplexEllipse.fit_from_cmplx_points(eta3)
            elli4 = fitted_ellipse = R2FMath.ComplexEllipse.fit_from_cmplx_points(eta4)    
            self.elli2.append(elli2)
            self.elli3.append(elli3)
            self.elli4.append(elli4)
            #gain3 = (elli3.eta_cw/elli2.eta_cw+elli3.eta_ccw/elli2.eta_ccw)/2
            #gain4 = (elli4.eta_cw/elli2.eta_cw+elli4.eta_ccw/elli2.eta_ccw)/2
            
            

            mm, c = np.polyfit(eta2, eta3, 1)

            g_left = 1/mm#elli2.eta_cw/elli3.eta_cw


            #y2 = elli2.semi_major/elli4.semi_major
            #y1 = elli2.semi_minor/elli4.semi_minor
            #x2 = elli2.semi_major
            #x1 = elli2.semi_minor
            #g_right = y1 -(y2-y1)/(x2-x1)*x1

            mm, c = np.polyfit(eta2, eta4, 1)

            g_right = 1/mm#elli2.eta_cw/elli3.eta_cw

            #g_right = elli2.eta_cw/elli4.eta_cw


            
            ratio3raw = (g_left*elli3.eta_o-elli2.eta_o)/1000 # divide by 1000 to get rid off 000*V2/V1
            ratio4raw = (g_right*elli4.eta_o-elli2.eta_o)/1000 # divide by 1000 to get rid off 000*V2/V1
            
            ratio3 = ratio3raw/10-1  # divide by 1000 to get rid off 000*V2/V1
            ratio4 = ratio4raw/10-1  # divide by 1000 to get rid off 000*V2/V1



            dratio = ratio4-ratio3

            al_left     = ratio3.real        # This is alpha_32 - alpha_31 ~ alpha_32
            D_left      = -ratio3.imag       # This is D_32 - D_31 ~ D_32
            al_right    = ratio4.real        # This is alpha_42 - alpha_41 ~ alpha_42
            D_right     = -ratio4.imag       # This is D_42 - D_41 ~ D_42
            al_diff     = dratio.real
            D_diff      = -dratio.imag
            line = np.hstack((f, np.abs(g_left), np.angle(g_left),al_left,D_left ,np.abs(g_right) ,np.angle(g_right), al_right,D_right,
                              al_diff,D_diff,np.real(ratio4raw),np.imag(ratio4raw)))
             #                 0    1                2              3     4          5             6                    7    8
            self.ana.append(line)
        self.ana= np.array(self.ana)
        self.di={'f':0,\
                 'left gain(abs)':1,
                 'left gain(ang)':2,
                 'left alpha':3,
                 'left D':4,
                 'right gain(abs)':5,
                 'right gain(ang)':6,
                 'right alpha':7,
                 'right D':8,
                 'diff alpha':9,
                 'diff D ':10,
                 'rawratio4(re)':11,
                 'rawratio4(im)':12               
                 }
        self.ana_mean, self.ana_std =  self.average(self.ana)
        self.f = self.ana_mean[:,0]
        indices = self.f.argsort()
        self.f = self.f[indices]
        self.ana_mean= self.ana_mean[indices,:]


    def average(self,output):
        mydict = {}
        for line in output:
            f = line[0]
            if  f not in mydict:
                mydict[f] = np.array(line)
            else:
                mydict[f] = np.vstack((mydict[f] ,np.array(line)))
        means =[]
        stds =[]
        for f in list(mydict):
            if len(np.shape(mydict[f][1]))==1:
                means.append(np.median(mydict[f],axis=0))
                stds.append(np.std(mydict[f],axis=0,ddof=1))
            else:
                means.append(np.array(mydict[f]))
                stds.append(np.zeros_like(mydict[f]))
        means = np.array(means)
        stds = np.array(stds)

        return means,stds                      



class completeSet:
    """
    Assumes that the caps are sorted from smallest to largest. It can be 1,2,3,4,5 caps.
    Measurements may have different frequency sets; the intersection (rounded to nearest Hz)
    is used.

    C0: reference capacitance (farads) for the first cap in the chain.
    """
    def __init__(self, bds, fns, C0=100e-12, fmax=500000):
        self.C0 = C0
        self.myCaps = []
        for b, f in zip(bds, fns):
            self.myCaps.append(oneCap(b, f))
        self.di = self.myCaps[0].di
        self._all_rounded = [np.round(cap.ana_mean[:,0]).astype(int) for cap in self.myCaps]
        self.analyze(fmax)

    def analyze(self, fmax=500000):
        self.fmax = fmax
        common = self._all_rounded[0]
        for r in self._all_rounded[1:]:
            common = np.intersect1d(common, r)
        common = common[common <= fmax]

        def select_rows(cap, common_hz):
            cap_rounded = np.round(cap.ana_mean[:,0]).astype(int)
            idx = [np.where(cap_rounded == f)[0][0] for f in common_hz]
            return cap.ana_mean[idx, :]

        ana_means = [select_rows(cap, common) for cap in self.myCaps]
        self.ana_means = ana_means

        self.f = ana_means[0][:,0]
        self.w = 2*np.pi*self.f
        AbsCap=[]
        RelCap=[]
        D =[]
        R =[]
        D0 =[]
        R0 =[]
        ix = np.argmin((self.f-1000)**2)
        oldcplx = None

        for i, ana_mean in enumerate(ana_means):
            ratio4raw = (1 + ana_mean[:, 7] + 1j * ana_mean[:, 8]) * 10
            gamma = ratio4raw
            if i == 0:
                thiscplx = gamma * self.C0
            else:
                thiscplx = gamma * oldcplx
            thiscap = np.real(thiscplx)
            thisD   = np.imag(thiscplx) / thiscap
            AbsCap.append(thiscap)
            RelCap.append(thiscap - thiscap[ix])
            D.append(thisD)
            D0.append(thisD - thisD[ix])
            R.append(thisD / (self.w * thiscap))
            R0.append(thisD / (self.w * thiscap) - thisD[ix] / (self.w[ix] * thiscap[ix]))
            oldcplx = thiscplx
        self.AbsCap = np.array(AbsCap)
        self.RelCap = np.array(RelCap)
        self.D = np.array(D)
        self.D0 = np.array(D0)
        self.R = np.array(R)
        self.R0 = np.array(R0)


        






