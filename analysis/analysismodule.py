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


def _parse_capval_pF(s):
    """Convert a capacitance label string to pF, e.g. '100pF'→100, '1nF'→1000."""
    s = s.strip()
    if s.endswith('pF'):
        return float(s[:-2])
    if s.endswith('nF'):
        return float(s[:-2]) * 1e3
    if s.endswith('uF') or s.endswith('μF'):
        return float(s[:-2]) * 1e6
    raise ValueError(f"Cannot parse capacitance label: {s!r}")


def analyze_block(V1, V2, V3, V4, fsig=None, C2_pF=None, correct_bias=False,
                  V1rb=None, V2rb=None):
    """
    Core phasor analysis for one ellipse block.
    V1..V4 : arrays of complex phasors in the phase-rotated frame (V1 is real).

    fsig        : signal frequency in Hz (required when C2_pF is given)
    C2_pF      : nominal C42 in pF.  When supplied, estimates I_bias by working
                  in the original (unrotated) frame using the SG readbacks V1rb/V2rb.
                  Model: V4_orig = a·V1rb + b  →  V4_Ibias = b − a·mean(V2rb).
                  V2rb is constant across the ellipse (fixed SG output), so the
                  cancellation is exact, unlike using the rotated DVM V2 channel.
    V1rb, V2rb  : SG readback phasors in the original frame (cols 10-11, 12-13 of
                  the data file).  Required for reliable I_bias estimation.
    correct_bias: subtract the estimated I_bias from V4 before computing the ratio.
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

    I_bias_pA = None
    V4_Ibias   = None
    if fsig is not None and C2_pF is not None and V1rb is not None and V2rb is not None:
        omega = 2 * np.pi * fsig
        # Work in the original (unrotated) frame.
        # V1_data is real (|V1_SG|); the rotation factor per point is V1rb/|V1rb|.
        rotation = V1rb / np.abs(V1rb)          # exp(+j·phi) for each ellipse point
        V4_orig  = V4 * rotation                 # unrotate V4

        # Fit V4_orig = a·V1rb + b.  V2rb ≈ const across the ellipse, so
        # V4_Ibias = b − a·mean(V2rb) cancels the V2 drive term exactly.
        A_bias = np.column_stack([V1rb, np.ones(len(V1rb))])
        params = np.linalg.lstsq(A_bias, V4_orig, rcond=None)[0]
        V4_Ibias = params[1] - params[0] * np.mean(V2rb)   # [V], original frame

        if omega * C2_pF < 2e6:   # ~3 kHz/100 pF — below this I_bias/I_C42 is extractable
            I_bias_pA = -V4_Ibias * (1j * omega * C2_pF) * g_right

        if correct_bias:
            # Subtract constant V4_Ibias in original frame, then re-rotate back
            V4 = (V4_orig - V4_Ibias) / rotation
            eta4  = 1000 * V4 / V1
            elli4 = R2FMath.ComplexEllipse.fit_from_cmplx_points(eta4)
            g_right = 1 / scipy.linalg.lstsq(A, eta4, lapack_driver='gelsy')[0][0]

    ratio3raw = (g_left  * elli3.eta_o - elli2.eta_o) / 1000
    ratio4raw = (g_right * elli4.eta_o - elli2.eta_o) / 1000
    ratio3 = ratio3raw / 10 - 1
    ratio4 = ratio4raw / 10 - 1

    return {
        'al_left':    ratio3.real,
        'D_left':    -ratio3.imag,
        'al_right':   ratio4.real,
        'D_right':   -ratio4.imag,
        'g_left':     g_left,
        'g_right':    g_right,
        'ratio3raw':  ratio3raw,
        'ratio4raw':  ratio4raw,
        'I_bias_pA':  I_bias_pA,   # None above ~3 kHz/100 pF — use V4_Ibias there
        'V4_Ibias':   V4_Ibias,    # [V] V4 offset due to I_bias, valid at all f
    }


class oneCap:
    """
    Load and analyse one capacitor position from one or more VOLT_ data files.

    bds, fns : str or list of str
        A single base-directory + filename pair, or parallel lists of
        base-directories and filenames.  When multiple files are given all
        blocks at the same frequency are pooled before averaging, so the
        resulting ana_err reflects run-to-run scatter rather than
        within-run scatter.
    """
    def __init__(self, bds, fns, C2_pF=None, correct_bias=False):
        if isinstance(bds, (str, os.PathLike)):
            bds = [bds]
        if isinstance(fns, (str, os.PathLike)):
            fns = [fns]

        self.ana   = []
        self.meanV1 = []
        self.elliV2 = []
        self.elliV3 = []
        self.elliV4 = []
        self.elli2  = []
        self.elli3  = []
        self.elli4  = []
        self.allf   = []
        self.aeta2  = []
        self.aeta3  = []
        self.aeta4  = []

        for bd, fn in zip(bds, fns):
            data = np.loadtxt(os.path.join(bd, fn))
            ave  = 0.5 * (data[0:-1:2, :] + data[1::2, :])
            for i in range(0, len(ave), 8):
                block = ave[i : i + 8, :]
                if len(block) < 8:
                    print(f"Skipping final partial block of size {len(block)} in {fn}")
                    break
                f   = np.mean(block[:, 0])
                V1  = block[:, 2] + 1j * block[:, 3]
                V2  = block[:, 4] + 1j * block[:, 5]
                V3  = block[:, 6] + 1j * block[:, 7]
                V4  = block[:, 8] + 1j * block[:, 9]
                # SG readbacks in original (unrotated) frame — used for I_bias fit
                V1rb = block[:, 10] + 1j * block[:, 11] if block.shape[1] > 10 else None
                V2rb = block[:, 12] + 1j * block[:, 13] if block.shape[1] > 12 else None

                self.allf.append(f)
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

                res = analyze_block(V1, V2, V3, V4,
                                    fsig=f, C2_pF=C2_pF, correct_bias=correct_bias,
                                    V1rb=V1rb, V2rb=V2rb)
                g_left, g_right = res['g_left'], res['g_right']
                dratio = res['ratio4raw'] - res['ratio3raw']
                def _cplx_or_nan(v):
                    return v if v is not None else complex(np.nan, np.nan)
                ibias  = _cplx_or_nan(res['I_bias_pA'])
                v4bias = _cplx_or_nan(res['V4_Ibias'])
                line = np.hstack((f,
                                  np.abs(g_left),  np.angle(g_left),
                                  res['al_left'],  res['D_left'],
                                  np.abs(g_right), np.angle(g_right),
                                  res['al_right'], res['D_right'],
                                  dratio.real, -dratio.imag,
                                  res['ratio4raw'].real, res['ratio4raw'].imag,
                                  ibias.real,  ibias.imag,
                                  v4bias.real, v4bias.imag))
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
            'I_bias_pA(re)':  13,   # NaN above ~3 kHz/100 pF — too noisy there
            'I_bias_pA(im)':  14,
            'V4_Ibias(re)':   15,   # [V] raw V4 offset; valid at all f
            'V4_Ibias(im)':   16,
        }
        self.ana_mean, self.ana_err = self.average(self.ana)
        self.f = self.ana_mean[:, 0]
        indices = self.f.argsort()
        self.f        = self.f[indices]
        self.ana_mean = self.ana_mean[indices, :]
        self.ana_err  = self.ana_err[indices, :]

    def average(self, output):
        mydict = {}
        for line in output:
            f = line[0]
            if f not in mydict:
                mydict[f] = np.array(line)
            else:
                mydict[f] = np.vstack((mydict[f], np.array(line)))
        means = []
        errs  = []
        for f in list(mydict):
            arr = mydict[f]
            if arr.ndim == 2:
                means.append(np.median(arr, axis=0))
                errs.append(np.std(arr, axis=0, ddof=1))
            else:
                # N=1: no spread estimable
                means.append(arr)
                errs.append(np.zeros_like(arr))
        return np.array(means), np.array(errs)


class completeSet:
    """
    Capacitor ladder from smallest to largest (1–5 caps).

    Preferred constructor: completeSet.from_runs() — see its docstring.

    Low-level constructor: pass parallel bds/fns lists, one element per cap.
    Each element may be a single string or a list of strings (multiple runs).
    """
    def __init__(self, bds, fns, C0=100, fmax=500000, C2_pF_list=None, correct_bias=False):
        self.C0 = C0
        self.myCaps = []
        for i, (b, f) in enumerate(zip(bds, fns)):
            c42 = C2_pF_list[i] if C2_pF_list is not None else None
            self.myCaps.append(oneCap(b, f, C2_pF=c42, correct_bias=correct_bias))
        self.di = self.myCaps[0].di
        self._all_rounded = [np.round(cap.ana_mean[:, 0]).astype(int) for cap in self.myCaps]
        self.analyze(fmax)

    @classmethod
    def from_runs(cls, base, entries, C0=100, ref_capval='100pF', fmax=500000,
                  correct_bias=False):
        """
        Convenience constructor — builds paths from a compact entry list.

        Parameters
        ----------
        base : str
            Root data directory, e.g. r'U:\\...\\CAPDATA'
        entries : list of [SN, capval, timestamp(s)]
            SN        — serial-number string, e.g. '1840J01469'
            capval    — cap label used in filenames, e.g. '1nF'
            timestamp — 'YYYYMMDD_HHMM' string, or a list of such strings
                        for multiple runs (pooled for better error bars)
        C0        : reference capacitance in pF (default 100 pF)
        ref_capval: label for C0 in filenames, e.g. '100pF'
        fmax      : maximum frequency passed to analyze()

        File path pattern:
            base / SN / YYMM / DD / VOLT_{capval}-{prev_capval}_{timestamp}.dat

        Examples
        --------
        Single run per cap:
            completeSet.from_runs(base, [
                ['1840J01469', '1nF',   '20260516_1414'],
                ['2519J00896', '10nF',  '20260516_2022'],
                ['2519J00896', '100nF', '20260515_1217'],
            ])

        Multiple runs per cap (pooled error bars):
            completeSet.from_runs(base, [
                ['1840J01469', '1nF',  ['20260516_1414', '20260501_1201']],
                ['2519J00896', '10nF', ['20260516_2022', '20260501_1401',
                                        '20260303_1533']],
            ])
        """
        bds = []
        fns = []
        C2_pF_list = []
        prev_capval = ref_capval
        for sn, capval, timestamps in entries:
            if isinstance(timestamps, str):
                timestamps = [timestamps]
            cap_bds, cap_fns = [], []
            for ts in timestamps:
                yymm = ts[2:6]   # YY+MM  e.g. '2605' from '20260516_1414'
                dd   = ts[6:8]   # DD     e.g. '16'
                cap_bds.append(os.path.join(base, sn, yymm, dd))
                cap_fns.append(f'VOLT_{capval}-{prev_capval}_{ts}.dat')
            bds.append(cap_bds)
            fns.append(cap_fns)
            C2_pF_list.append(_parse_capval_pF(prev_capval))
            prev_capval = capval
        return cls(bds, fns, C0=C0, fmax=fmax,
                   C2_pF_list=C2_pF_list, correct_bias=correct_bias)

    def analyze(self, fmax=500000):
        self.fmax = fmax
        common = self._all_rounded[0]
        for r in self._all_rounded[1:]:
            common = np.intersect1d(common, r)
        common = common[common <= fmax]

        def select_rows(cap, common_hz):
            cap_rounded = np.round(cap.ana_mean[:, 0]).astype(int)
            idx = [np.where(cap_rounded == f)[0][0] for f in common_hz]
            return cap.ana_mean[idx, :], cap.ana_err[idx, :]

        ana_data  = [select_rows(cap, common) for cap in self.myCaps]
        ana_means = [d[0] for d in ana_data]
        ana_errs  = [d[1] for d in ana_data]
        self.ana_means = ana_means

        self.f = ana_means[0][:, 0]
        self.w = 2 * np.pi * self.f
        ix = np.argmin((self.f - 1000) ** 2)
        oldcplx = None
        err_old_r = err_old_i = None

        AbsCap, RelCap, D, D0, R, R0 = [], [], [], [], [], []
        AbsCap_err, RelCap_err, D_err, D0_err, R_err, R0_err = [], [], [], [], [], []

        for i, (ana_mean, ana_err) in enumerate(zip(ana_means, ana_errs)):
            # col 7 = al_right = ratio4.real, col 8 = D_right = -ratio4.imag
            # so gamma.imag must use -D_right to get back +ratio4.imag
            ratio4raw = (1 + ana_mean[:, 7] - 1j * ana_mean[:, 8]) * 10
            gamma   = ratio4raw
            gamma_r = np.real(gamma)
            gamma_i = np.imag(gamma)
            err_gamma_r = 10 * ana_err[:, 7]
            err_gamma_i = 10 * ana_err[:, 8]

            if i == 0:
                thiscplx   = gamma * self.C0
                err_cplx_r = err_gamma_r * self.C0
                err_cplx_i = err_gamma_i * self.C0
            else:
                thiscplx = gamma * oldcplx
                old_r = np.real(oldcplx)
                old_i = np.imag(oldcplx)
                # Re(gamma * old) = gamma_r*old_r - gamma_i*old_i
                err_cplx_r = np.sqrt((old_r * err_gamma_r) ** 2 +
                                     (old_i * err_gamma_i) ** 2 +
                                     (gamma_r * err_old_r) ** 2 +
                                     (gamma_i * err_old_i) ** 2)
                # Im(gamma * old) = gamma_r*old_i + gamma_i*old_r
                err_cplx_i = np.sqrt((old_i * err_gamma_r) ** 2 +
                                     (old_r * err_gamma_i) ** 2 +
                                     (gamma_i * err_old_r) ** 2 +
                                     (gamma_r * err_old_i) ** 2)

            thiscap = np.real(thiscplx)
            thisD   = -np.imag(thiscplx) / thiscap   # Im(C) < 0 for positive D (C·(1-jD) convention)
            err_thiscap = err_cplx_r
            err_thisD   = np.sqrt(err_cplx_i ** 2 + (thisD * err_thiscap) ** 2) / thiscap
            err_R       = np.sqrt(err_thisD ** 2 +
                                  (thisD / thiscap * err_thiscap) ** 2) / (self.w * thiscap)

            AbsCap.append(thiscap)
            AbsCap_err.append(err_thiscap)
            RelCap.append(thiscap - thiscap[ix])
            RelCap_err.append(np.sqrt(err_thiscap ** 2 + err_thiscap[ix] ** 2))
            D.append(thisD)
            D_err.append(err_thisD)
            D0.append(thisD - thisD[ix])
            D0_err.append(np.sqrt(err_thisD ** 2 + err_thisD[ix] ** 2))
            R.append(thisD / (self.w * thiscap))
            R_err.append(err_R)
            R0.append(thisD / (self.w * thiscap) - thisD[ix] / (self.w[ix] * thiscap[ix]))
            R0_err.append(np.sqrt(err_R ** 2 + err_R[ix] ** 2))

            oldcplx   = thiscplx
            err_old_r = err_cplx_r
            err_old_i = err_cplx_i

        self.AbsCap     = np.array(AbsCap)
        self.RelCap     = np.array(RelCap)
        self.D          = np.array(D)
        self.D0         = np.array(D0)
        self.R          = np.array(R)
        self.R0         = np.array(R0)
        self.AbsCap_err = np.array(AbsCap_err)
        self.RelCap_err = np.array(RelCap_err)
        self.D_err      = np.array(D_err)
        self.D0_err     = np.array(D0_err)
        self.R_err      = np.array(R_err)
        self.R0_err     = np.array(R0_err)
