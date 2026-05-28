import numpy as np
import scipy.optimize



def fit_sine_cplx(y, fsamp, fsig, fline=60, Nhars=1, use_hann=True, chunk_periods=0):
    """
    Fits a sine wave with a DC offset to the data.

    When chunk_periods > 0, uses two-stage fitting:
      1. Global Hanning-windowed fit removes DC, line harmonics, and f/2 DDS spur.
      2. The signal-only residual is split into chunk_periods-long windows; each chunk
         gets its own [DC, cos, sin, cos/2, sin/2] fit. The array of per-chunk complex
         amplitudes is returned as a 5th element so that eta ratios can be computed
         per chunk before averaging (avoiding amplitude-averaging bias).

    Returns:
        4-tuple (complex_amp, fit_vals, errv, rss)                   when chunk_periods == 0
        5-tuple (complex_amp, fit_vals, errv, rss, chunk_amps)       when chunking succeeds
    """
    y = np.asarray(y, dtype=float)
    n = len(y)

    rf  = fsig  / fsamp
    rlf = fline / fsamp

    wt  = 2 * np.pi * np.arange(n) * rf
    wlf = 2 * np.pi * np.arange(n) * rlf

    cols = [np.ones(n), np.cos(wt), np.sin(wt),
            np.cos(wt / 2), np.sin(wt / 2),          # fsig/2 DDS spur — absorbed, discarded
            np.cos(wlf), np.sin(wlf), np.cos(2*wlf), np.sin(2*wlf)]
    for h in range(2, Nhars + 1):
        cols.extend([np.cos(h * wt), np.sin(h * wt)])

    # If any line harmonic falls within fline/2 of fsig, add it explicitly so it
    # doesn't leak into the signal estimate (e.g. 42×60=2520 Hz near 2512 Hz).
    n_extra_line = 0
    n_nearest = int(round(fsig / fline))
    for nh in range(max(3, n_nearest - 2), n_nearest + 3):
        fn = nh * fline
        if fsamp / n < abs(fn - fsig) < fline / 2:   # close but distinguishable
            wnh = 2 * np.pi * np.arange(n) * fn / fsamp
            cols.extend([np.cos(wnh), np.sin(wnh)])
            n_extra_line += 2

    X = np.column_stack(cols)

    n_params = 1 + 2 * Nhars + 4 + 2 + n_extra_line  # DC + signal harmonics + line harmonics + fsig/2 + nearby line
    ndf = n - n_params

    if use_hann:
        w = np.hanning(n)
        X_eff = X * w[:, np.newaxis]
        y_eff = y * w
    else:
        X_eff = X
        y_eff = y

    fit_pars, lstsq_res, _, _ = np.linalg.lstsq(X_eff, y_eff, rcond=None)

    # Check whether per-chunk fitting is feasible
    n_chunk_params = 5  # DC, cos(wt), sin(wt), cos(wt/2), sin(wt/2)
    chunk_size = int(round(chunk_periods * fsamp / fsig)) if chunk_periods > 0 else 0
    use_chunks = chunk_size > n_chunk_params + 1 and n // chunk_size >= 2

    if not use_chunks:
        fit_vals = X @ fit_pars
        rss = float(lstsq_res[0]) if len(lstsq_res) == 1 else float(np.sum((y_eff - X_eff @ fit_pars) ** 2))
        errv = np.sqrt(rss / ndf) if ndf > 0 else 0.0
        return fit_pars[1] - 1j * fit_pars[2], fit_vals, errv, rss

    # Stage 1: remove background (DC, line harmonics, f/2) leaving signal only
    bg_pars = fit_pars.copy()
    bg_pars[1] = 0.0   # zero cos(wt) fundamental
    bg_pars[2] = 0.0   # zero sin(wt) fundamental
    y_bg  = X @ bg_pars       # DC + line + f/2 + higher harmonics
    y_sig = y - y_bg           # signal-only residual, full record

    # Stage 2: fit per-chunk amplitudes
    n_chunks = n // chunk_size
    chunk_amps = np.zeros(n_chunks, dtype=complex)
    fit_sig    = np.zeros(n, dtype=float)
    rss_total  = 0.0
    ndf_total  = 0

    for k in range(n_chunks):
        i0  = k * chunk_size
        i1  = i0 + chunk_size
        idx = np.arange(i0, i1)
        wtc = 2 * np.pi * idx * rf
        Xc  = np.column_stack([np.ones(chunk_size),
                                np.cos(wtc), np.sin(wtc),
                                np.cos(wtc / 2), np.sin(wtc / 2)])
        yc = y_sig[i0:i1]
        if use_hann:
            wc = np.hanning(chunk_size)
            pars_c, _, _, _ = np.linalg.lstsq(Xc * wc[:, np.newaxis], yc * wc, rcond=None)
        else:
            pars_c, _, _, _ = np.linalg.lstsq(Xc, yc, rcond=None)
        chunk_amps[k] = pars_c[1] - 1j * pars_c[2]
        fit_sig[i0:i1] = Xc @ pars_c
        rss_total += float(np.sum((yc - Xc @ pars_c) ** 2))
        ndf_total += chunk_size - n_chunk_params

    # Tail samples not covered by whole chunks: use global fundamental
    tail = n_chunks * chunk_size
    if tail < n:
        tail_idx = np.arange(tail, n)
        fit_sig[tail:] = (fit_pars[1] * np.cos(2 * np.pi * tail_idx * rf) +
                          fit_pars[2] * np.sin(2 * np.pi * tail_idx * rf))

    complex_amp = np.mean(chunk_amps)
    fit_vals    = y_bg + fit_sig
    errv        = np.sqrt(rss_total / ndf_total) if ndf_total > 0 else 0.0

    return complex_amp, fit_vals, errv, rss_total, chunk_amps


def get_f(y, fsamp, fsig_guess, fline_guess=60.0, use_hann=True, Nhars=1):
    """
    Estimates the signal frequency by minimizing the residual sum of squares
    from fit_sine_cplx using Brent's bounded method (guaranteed convergence).
    chunk_periods is intentionally left at 0 here for speed.
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
    def __init__(self, fsig, fsamp, data, Nhars=1, chunk_periods=0):
        self.data          = np.array(data)
        self.fsig          = fsig
        self.fsamp         = fsamp
        self.fline         = 60
        self.Nhars         = Nhars
        self.chunk_periods = chunk_periods
        self.Vc_chunks     = None

    def setf(self, fsig, fline):
        self.fsig, self.fline = fsig, fline

    def findf(self):
        self.fsig, self.fline = get_f(self.data, self.fsamp, self.fsig, self.fline, Nhars=self.Nhars)
        return self.fsig, self.fline

    def fit(self):
        result = fit_sine_cplx(self.data, self.fsamp, self.fsig, self.fline,
                                self.Nhars, True, self.chunk_periods)
        self.Vc        = result[0]
        self.fv        = result[1]
        self.c2        = result[3]
        self.Vc_chunks = result[4] if len(result) > 4 else None

    def strip_raw(self):
        self.data = None
        self.fv   = None


class FourChannels:
    def __init__(self, fsig, fsamp, Nhars, ch1, ch2, ch3, ch4, V1c, V2c, i, ts,
                 chunk_periods=0):
        self.ts            = ts
        self.fsig          = fsig
        self.fsamp         = fsamp
        self.i             = i
        self.Data          = []
        self.V1c           = V1c
        self.V2c           = V2c
        self.phasor_chunks = None
        if ts < 0:
            return
        self.Data.append(SampleData(fsig, fsamp, ch1, Nhars, chunk_periods))
        self.Data.append(SampleData(fsig, fsamp, ch2, Nhars, chunk_periods))
        self.Data.append(SampleData(fsig, fsamp, ch3, Nhars, chunk_periods))
        self.Data.append(SampleData(fsig, fsamp, ch4, Nhars, chunk_periods))
        fsig, fline = self.Data[0].findf()
        # Clamp Nhars so no harmonic exceeds Nyquist (0.45 * fsamp as margin)
        nhars_used = max(1, min(Nhars, int(0.45 * fsamp / fsig)))
        for i in range(4):
            self.Data[i].Nhars = nhars_used
            self.Data[i].setf(fsig, fline)
            self.Data[i].fit()

        # Build per-chunk phasor matrix (n_chunks × 4): rotate each chunk by ch0's phase
        chunks0 = self.Data[0].Vc_chunks
        if chunks0 is not None and len(chunks0) > 1:
            n_ch = min(len(d.Vc_chunks) for d in self.Data if d.Vc_chunks is not None)
            self.phasor_chunks = np.zeros((n_ch, 4), dtype=complex)
            for k in range(n_ch):
                cf_k = np.exp(-1j * np.angle(chunks0[k]))
                for j in range(4):
                    self.phasor_chunks[k, j] = self.Data[j].Vc_chunks[k] * cf_k

    def strip_raw(self):
        for sd in self.Data:
            sd.strip_raw()


class NPoints:
    def __init__(self, fsig, fsamp, Nhars=1, g1=1, g2=1, ratio=10, N=8, chunk_periods=0):
        self.N             = N
        self.chunk_periods = chunk_periods
        self.Res = {}
        self.Res['fsig']  = fsig
        self.Res['ratio'] = ratio
        self.Res['fsamp'] = fsamp
        self.Res['Nhars'] = Nhars
        self.Res['gain1'] = g1
        self.Res['gain2'] = g2
        self.ats = -1 * np.ones(N)
        self.Res['ts']    = min(self.ats)
        self.Data = np.zeros(N, dtype=object)

    def setPoint(self, i, ch1, ch2, ch3, ch4, V1c, V2c, ts):
        self.V1c = V1c
        self.Data[i] = FourChannels(
            self.Res['fsig'], self.Res['fsamp'], self.Res['Nhars'],
            ch1, ch2, ch3, ch4, V1c, V2c, i, ts,
            chunk_periods=self.chunk_periods,
        )
        self.ats[i] = ts
        if min(self.ats) > 0:
            self.Res['ts'] = np.mean(self.ats)

    def precalc(self):
        self.raw8 = np.zeros((self.N, 4), dtype=complex)
        self.ctrl = np.zeros((self.N, 2), dtype=complex)
        for i in range(self.N):
            phi = np.angle(self.Data[i].Data[0].Vc)
            cf = np.exp(-1j * phi)
            for j in range(4):
                self.raw8[i, j] = self.Data[i].Data[j].Vc * cf
            self.ctrl[i, 0] = self.Data[i].V1c
            self.ctrl[i, 1] = self.Data[i].V2c
        self.ave4  = 0.5 * (self.raw8[::2, :] + self.raw8[1::2, :])
        self.ctrla = 0.5 * (self.ctrl[::2, :] + self.ctrl[1::2, :])

    def _build_fit_etas(self):
        """
        Returns (eta2, eta3, eta4) arrays for the ellipse lstsq.

        When per-chunk phasors are available, each ellipse point contributes n_chunks
        data rows instead of one: the per-chunk etas from the two paired measurements
        are averaged to cancel the relay-position offset, then concatenated across all
        ellipse points.  Falls back to the N/2 per-point etas when chunking was not used.
        """
        N2 = self.N // 2
        has_chunks = all(
            isinstance(self.Data[i], FourChannels) and self.Data[i].phasor_chunks is not None
            for i in range(self.N)
        )
        if not has_chunks:
            return self.eta2, self.eta3, self.eta4

        n_ch = min(len(self.Data[i].phasor_chunks) for i in range(self.N))
        if n_ch < 2:
            return self.eta2, self.eta3, self.eta4

        e2, e3, e4 = [], [], []
        for m in range(N2):
            pc_even = self.Data[2*m].phasor_chunks[:n_ch]    # (n_ch, 4)
            pc_odd  = self.Data[2*m+1].phasor_chunks[:n_ch]  # (n_ch, 4)
            v1_even = pc_even[:, 0]   # real and positive after per-chunk rotation
            v1_odd  = pc_odd[:, 0]
            # Compute etas per chunk per measurement, then average the paired estimates
            e2.append(0.5 * (pc_even[:, 1] / v1_even + pc_odd[:, 1] / v1_odd))
            e3.append(0.5 * (pc_even[:, 2] / v1_even + pc_odd[:, 2] / v1_odd))
            e4.append(0.5 * (pc_even[:, 3] / v1_even + pc_odd[:, 3] / v1_odd))
        return np.concatenate(e2), np.concatenate(e3), np.concatenate(e4)

    def calc(self):
        self.precalc()
        self.V1m = self.ave4[:, 0]
        self.V2m = self.ave4[:, 1]
        self.V3m = self.ave4[:, 2]
        self.V4m = self.ave4[:, 3]
        self.eta2 = self.V2m / self.V1m
        self.eta3 = self.V3m / self.V1m
        self.eta4 = self.V4m / self.V1m

        eta2_fit, eta3_fit, eta4_fit = self._build_fit_etas()

        Xmat = np.column_stack([eta2_fit, np.ones(len(eta2_fit))])
        (m3, c3), _, _, _ = np.linalg.lstsq(Xmat, eta3_fit, rcond=None)
        (m4, c4), _, _, _ = np.linalg.lstsq(Xmat, eta4_fit, rcond=None)
        self.gamma3 = 1.0 / m3
        self.gamma4 = 1.0 / m4
        self.Res['Vz3'] = c3 / m3
        self.Res['Vz4'] = c4 / m4
        self.Res['gamma3'] = self.gamma3
        self.Res['gamma4'] = self.gamma4
        N2 = self.N // 2
        V4_raw = np.array([0.5*(self.Data[2*k].Data[3].Vc + self.Data[2*k+1].Data[3].Vc)
                           for k in range(N2)])
        V1_meas = np.array([0.5*(self.Data[2*k].Data[0].Vc + self.Data[2*k+1].Data[0].Vc)
                            for k in range(N2)])
        eta4 = V4_raw / V1_meas
        V1rb_pts = self.ctrla[:, 0]
        V1rb_center = np.mean(V1rb_pts)
        Xmat = np.column_stack([V1rb_pts - V1rb_center, np.ones(N2)])
        (m_bal, c_bal), _, _, _ = np.linalg.lstsq(Xmat, eta4, rcond=None)
        step = c_bal / m_bal
        dV1_spread = np.std(V1rb_pts)
        step_limit = 5 * dV1_spread
        if np.abs(step) > step_limit:
            step = step * step_limit / np.abs(step)
        self.Res['eta4fit_slope'] = m_bal
        self.Res['eta4fit_intercept'] = c_bal
        self.Res['eta4_mean'] = np.mean(eta4)
        self.Res['V1rb_center'] = V1rb_center
        self.Res['V1_balance'] = V1rb_center - 0.5 * step

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

    def max_raw_amplitude(self, channel_idx):
        """Return max |raw sample| across all measurement points for the given channel index."""
        mx = 0.0
        for d in self.Data:
            if isinstance(d, FourChannels) and len(d.Data) > channel_idx:
                raw = d.Data[channel_idx].data
                if raw is not None:
                    mx = max(mx, float(np.max(np.abs(raw))))
        return mx

    def strip_raw(self):
        for fc in self.Data:
            if isinstance(fc, FourChannels):
                fc.strip_raw()


class AllData():
    def __init__(self):
        self.mydict = {}

    def append(self, ND: NPoints):
        f = ND.Res['fsig']
        if f not in self.mydict:
            self.mydict[f] = []
        self.mydict[f].append(ND)

    def deletekey(self, f):
        if f in self.mydict:
            del self.mydict[f]

    def count(self):
        return len(self.mydict)

    def countf(self, f):
        if f not in self.mydict:
            return 0
        else:
            return len(self.mydict[f])

    def getkeys(self, f, keys):
        retdict = {}
        L = len(self.mydict[f])
        if L == 0:
            for k in keys:
                retdict[k] = np.array([])
            return retdict
        for k in keys:
            retdict[k] = np.empty(L, dtype=object)
        for n, a in enumerate(self.mydict[f]):
            for k in keys:
                retdict[k][n] = a.Res[k]
        return retdict

    def getallkeys(self, f):
        if f not in self.mydict or len(self.mydict[f]) == 0:
            return {}
        keys = list(self.mydict[f][0].Res.keys())
        return self.getkeys(f, keys)

    def getdictf(self, keys):
        retdict = {}
        keys = list(keys)
        if 'fsig' not in keys:
            keys.append('fsig')
        allf = list(self.mydict)
        if len(allf) == 0:
            for k in keys:
                retdict[k] = np.array([])
            return retdict
        else:
            for k in keys:
                retdict[k] = []
            for f in allf:
                odict = self.getkeys(f, keys)
                for k in keys:
                    for item in odict[k]:
                        retdict[k].append(item)
            for k in keys:
                retdict[k] = np.array(retdict[k])
            return retdict

    def getRawPhasors(self, f, t0=0):
        ret = []
        for obj in self.mydict[f]:
            Nrows = np.shape(obj.raw8)[0]
            for j in range(Nrows):
                line = np.hstack((f, obj.Res['ts'] - t0))
                for k in range(np.shape(obj.raw8)[1]):
                    line = np.hstack((line, obj.raw8[j, k].real, obj.raw8[j, k].imag))
                for k in range(np.shape(obj.ctrla)[1]):
                    line = np.hstack((line, obj.ctrla[j//2, k].real, obj.ctrla[j//2, k].imag))
                line = np.hstack((line, obj.Res['gain1'], obj.Res['gain2']))
                ret.append(line)
        return np.array(ret)
