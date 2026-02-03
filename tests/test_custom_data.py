#%%
import sys
import os
import pytest
import numpy as np
import matplotlib.pyplot as plt
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
import CustomData


def setAmp(co,NDpts=8,dV2=0.1,V1c= 9.9, V2c= -0.9):
        a = dV2*1.2
        b = dV2*0.8
        ang = (co//2)/NDpts*2*np.pi        
        theta = 30/180*np.pi
        V2cplxAmp = V2c+np.exp(1j*theta)*(a*np.cos(ang)+1j*b*np.sin(ang))
        V1cplxAmp = V1c
        return V1cplxAmp,V2cplxAmp



def simulator( V1cplxAmp, V2cplxAmp, fsig, fsamp, Nsamp,C31,C32,C41,C42):
    """
    Suppose signal is  V cos(2 pi fsig t) and t = i /fsamp, then the signal is
    V cos(2 pi fsig/fsmap i
    """
    VNoise = 1e-1
    gain3 = 1e6 # It goes from 1e3 V/A to 1e8 V/A
    gain4 = 1e6    
    t = np.arange(Nsamp)/fsamp
    wsig = 2 *np.pi*fsig
    Y31 = 1j*wsig*C31
    Y32 = 1j*wsig*C32
    Y41 = 1j*wsig*C41
    Y42 = 1j*wsig*C42
    V1 = V1cplxAmp* np.exp(1j*wsig*t) 
    V2 = V2cplxAmp* np.exp(1j*wsig*t)
    I31 = V1*Y31
    I32 = V2*Y32
    I3  = I31+I32
    I41 = V1*Y41
    I42 = V2*Y42
    I4  = I41+I42
    V3 = I3*gain3 + VNoise*np.random.randn(Nsamp)
    V4 = I4*gain4 + VNoise*np.random.randn(Nsamp)
    V1 +=  VNoise*np.random.randn(Nsamp)
    V2 +=  VNoise*np.random.randn(Nsamp)
    return np.real(V1),np.real(V2),np.real(V3),np.real(V4)
#%%

# %%
Npts = 16
fsig = 1000
fsamp = 800000
g1 = 1000
g2 =1000
rawN = CustomData.NPoints(fsig,fsamp,g1=g1,g2=g2,N=Npts)
Nsamp = int(fsamp/fsig) # This gives one period
# %%
allData = CustomData.AllData()

for co in range(16):
    V1rb,V2rb =  setAmp(co)
    ch1,ch2,ch3, ch4 = simulator(V1rb,V2rb, fsig, fsamp, Nsamp,1e-10,1e-9,1e-10,1e-9)
    rawN.setPoint(co,ch1,ch2,ch3,ch4,V2rb,V1rb,time.time())    
    #V1rb amd 2rb are switched here bc of the channel issue
print()
rawN.calc()
allData.append(rawN)
#%%
fig,ax = plt.subplots(2,2)
ax=ax.flatten()
for i in range(4):
    ax[i].plot(np.real(rawN.ave4[:,i]),np.imag(rawN.ave4[:,i]),'ro')
    if i>0:
        ax[i].scatter(**rawN.Cir[i].plot_mycircle(),marker='.', s=1,cmap='Greens')
ax[i].ticklabel_format(useOffset=False)


# %%
u=allData.getRawPhasors(fsig)
# %%
ax[0].plot(V1)
ax[1].plot(V2)
ax[2].plot(V3)
ax[3].plot(V4)

    

# %%
