import numpy as np
import scipy.linalg
import sys
import os

try:
    import R2FMath
except ImportError:
    _src = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src')
    sys.path.insert(0, _src)
    import R2FMath


def analyze_block(V1, V2, V3, V4):
    """
    Core phasor analysis for one ellipse block (8 switch-state-averaged points).
    V1..V4: arrays of 8 complex phasors.
    Returns al_left, D_left, al_right, D_right, g_left, g_right, ratio3raw, ratio4raw.
    """
    eta2 = 1000 * V2 / V1
    eta3 = 1000 * V3 / V1
    eta4 = 1000 * V4 / V1

    elli2 = R2FMath.ComplexEllipse.fit_from_cmplx_points(eta2)
    elli3 = R2FMath.ComplexEllipse.fit_from_cmplx_points(eta3)
    elli4 = R2FMath.ComplexEllipse.fit_from_cmplx_points(eta4)

    A = np.column_stack([eta2, np.ones(len(eta2))])
    g_left  = 1 / scipy.linalg.lstsq(A, eta3, lapack_driver='gelsy')[0][0]
    g_right = 1 / scipy.linalg.lstsq(A, eta4, lapack_driver='gelsy')[0][0]

    ratio3raw = (g_left  * elli3.eta_o - elli2.eta_o) / 1000
    ratio4raw = (g_right * elli4.eta_o - elli2.eta_o) / 1000
    ratio3 = ratio3raw / 10 - 1
    ratio4 = ratio4raw / 10 - 1

    return {
        'al_left':   ratio3.real,
        'D_left':   -ratio3.imag,
        'al_right':  ratio4.real,
        'D_right':  -ratio4.imag,
        'g_left':    g_left,
        'g_right':   g_right,
        'ratio3raw': ratio3raw,
        'ratio4raw': ratio4raw,
    }


class oneCap:
    def __init__(self, bd, fn):
        data = np.loadtxt(os.path.join(bd, fn))
        self.ave = 0.5 * (data[0:-1:2, :] + data[1::2, :])
        self.ana = []
        self.meanV1 = []
        self.elliV2 = []
        self.elliV3 = []
        self.elliV4 = []
        self.elli2 = []
        self.elli3 = []
        self.elli4 = []
        self.allf = []
        self.aeta2 = []
        self.aeta3 = []
        self.aeta4 = []
        for i in range(0, len(self.ave), 8):
            block = self.ave[i : i + 8, :]
            if len(block) < 8:
                print(f"Processing final partial block of size {len(block)}")
                break
            f = np.mean(block[:, 0])
            self.allf.append(f)
            V1 = block[:, 2] + 1j * block[:, 3]
            V2 = block[:, 4] + 1j * block[:, 5]
            V3 = block[:, 6] + 1j * block[:, 7]
            V4 = block[:, 8] + 1j * block[:, 9]

            self.meanV1.append(np.mean(V1))
            self.elliV2.append(R2FMath.ComplexEllipse.fit_from_cmplx_points(V2))
            self.elliV3.append(R2FMath.ComplexEllipse.fit_from_cmplx_points(V3))
            self.elliV4.append(R2FMath.ComplexEllipse.fit_from_cmplx_points(V4))

            eta2 = 1000 * V2 / V1
            eta3 = 1000 * V3 / V1
            eta4 = 1000 * V4 / V1
            self.aeta2.append(eta2)
            self.aeta3.append(eta3)
            self.aeta4.append(eta4)
            self.elli2.append(R2FMath.ComplexEllipse.fit_from_cmplx_points(eta2))
            self.elli3.append(R2FMath.ComplexEllipse.fit_from_cmplx_points(eta3))
            self.elli4.append(R2FMath.ComplexEllipse.fit_from_cmplx_points(eta4))

            res = analyze_block(V1, V2, V3, V4)
            g_left, g_right = res['g_left'], res['g_right']
            dratio = res['ratio4raw'] - res['ratio3raw']
            line = np.hstack((f,
                              np.abs(g_left),  np.angle(g_left),
                              res['al_left'],  res['D_left'],
                              np.abs(g_right), np.angle(g_right),
                              res['al_right'], res['D_right'],
                              dratio.real, -dratio.imag,
                              res['ratio4raw'].real, res['ratio4raw'].imag))
            self.ana.append(line)
        self.ana = np.array(self.ana)
        self.di = {
            'f':               0,
            'left gain(abs)':  1,
            'left gain(ang)':  2,
            'left alpha':      3,
            'left D':          4,
            'right gain(abs)': 5,
            'right gain(ang)': 6,
            'right alpha':     7,
            'right D':         8,
            'diff alpha':      9,
            'diff D ':        10,
            'rawratio4(re)':  11,
            'rawratio4(im)':  12,
        }
        self.ana_mean, self.ana_std = self.average(self.ana)
        self.f = self.ana_mean[:, 0]
        indices = self.f.argsort()
        self.f = self.f[indices]
        self.ana_mean = self.ana_mean[indices, :]

    def average(self, output):
        mydict = {}
        for line in output:
            f = line[0]
            if f not in mydict:
                mydict[f] = np.array(line)
            else:
                mydict[f] = np.vstack((mydict[f], np.array(line)))
        means = []
        stds = []
        for f in list(mydict):
            if len(np.shape(mydict[f][1])) == 1:
                means.append(np.median(mydict[f], axis=0))
                stds.append(np.std(mydict[f], axis=0, ddof=1))
            else:
                means.append(np.array(mydict[f]))
                stds.append(np.zeros_like(mydict[f]))
        return np.array(means), np.array(stds)


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
        self._all_rounded = [np.round(cap.ana_mean[:, 0]).astype(int) for cap in self.myCaps]
        self.analyze(fmax)

    def analyze(self, fmax=500000):
        self.fmax = fmax
        common = self._all_rounded[0]
        for r in self._all_rounded[1:]:
            common = np.intersect1d(common, r)
        common = common[common <= fmax]

        def select_rows(cap, common_hz):
            cap_rounded = np.round(cap.ana_mean[:, 0]).astype(int)
            idx = [np.where(cap_rounded == f)[0][0] for f in common_hz]
            return cap.ana_mean[idx, :]

        ana_means = [select_rows(cap, common) for cap in self.myCaps]
        self.ana_means = ana_means

        self.f = ana_means[0][:, 0]
        self.w = 2 * np.pi * self.f
        ix = np.argmin((self.f - 1000) ** 2)
        oldcplx = None

        AbsCap, RelCap, D, D0, R, R0 = [], [], [], [], [], []
        for i, ana_mean in enumerate(ana_means):
            ratio4raw = (1 + ana_mean[:, 7] + 1j * ana_mean[:, 8]) * 10
            gamma = ratio4raw
            thiscplx = gamma * self.C0 if i == 0 else gamma * oldcplx
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
