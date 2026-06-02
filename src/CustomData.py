import numpy as np
import scipy.optimize
import pandas as pd



def fit_sine_cplx(y, fsamp, fsig, fline=60, Nhars=1, use_hann=True, chunk_periods=0, cache=None):
    """
    Fits a sine wave with a DC offset to the data.

    When chunk_periods > 0, uses two-stage fitting:
      1. Global Hanning-windowed fit removes DC, line harmonics, and f/2 DDS spur.
      2. The signal-only residual is split into chunk_periods-long windows; each chunk
         gets its own [DC, cos, sin, cos/2, sin/2] fit. The array of per-chunk complex
         amplitudes is returned as a 5th element so that eta ratios can be computed
         per chunk before averaging (avoiding amplitude-averaging bias).

    When cache= is supplied (from build_fit_cache), the design-matrix construction and
    inversion are skipped; each call reduces to a matrix–vector multiply.

    Returns:
        4-tuple (complex_amp, fit_vals, errv, rss)                   when chunk_periods == 0
        5-tuple (complex_amp, fit_vals, errv, rss, chunk_amps)       when chunking succeeds
    """
    y = np.asarray(y, dtype=float)
    n = len(y)
    n_chunk_params = 5  # DC, cos(wt), sin(wt), cos(wt/2), sin(wt/2)

    if cache is not None:
        X, w, bg_pinv, rf = cache['X'], cache['w'], cache['bg_pinv'], cache['rf']
        fit_pars   = bg_pinv @ (y * w)
        use_chunks = 'chunk_pinvs' in cache
    else:
        rf  = fsig  / fsamp
        rlf = fline / fsamp
        wt  = 2 * np.pi * np.arange(n) * rf
        wlf = 2 * np.pi * np.arange(n) * rlf
        cols = [np.ones(n), np.cos(wt), np.sin(wt),
                np.cos(wt / 2), np.sin(wt / 2),
                np.cos(wlf), np.sin(wlf), np.cos(2*wlf), np.sin(2*wlf)]
        for h in range(2, Nhars + 1):
            cols.extend([np.cos(h * wt), np.sin(h * wt)])
        # If any line harmonic falls within fline/2 of fsig, add it explicitly so it
        # doesn't leak into the signal estimate (e.g. 42×60=2520 Hz near 2512 Hz).
        n_extra_line = 0
        n_nearest = int(round(fsig / fline))
        for nh in range(max(3, n_nearest - 2), n_nearest + 3):
            fn = nh * fline
            if fsamp / n < abs(fn - fsig) < fline / 2:
                wnh = 2 * np.pi * np.arange(n) * fn / fsamp
                cols.extend([np.cos(wnh), np.sin(wnh)])
                n_extra_line += 2
        X = np.column_stack(cols)
        w = np.hanning(n) if use_hann else np.ones(n)
        X_eff = X * w[:, np.newaxis]
        y_eff = y * w
        fit_pars, lstsq_res, _, _ = np.linalg.lstsq(X_eff, y_eff, rcond=None)
        chunk_size = int(round(chunk_periods * fsamp / fsig)) if chunk_periods > 0 else 0
        use_chunks = chunk_size > n_chunk_params + 1 and n // chunk_size >= 2

    fit_vals = X @ fit_pars
    rss  = float(np.sum((w * (y - fit_vals)) ** 2))
    ndf  = n - X.shape[1]
    errv = np.sqrt(rss / ndf) if ndf > 0 else 0.0

    if not use_chunks:
        return fit_pars[1] - 1j * fit_pars[2], fit_vals, errv, rss

    # Stage 1: remove background (DC, line harmonics, f/2) leaving signal only
    bg_pars    = fit_pars.copy()
    bg_pars[1] = 0.0
    bg_pars[2] = 0.0
    y_bg  = X @ bg_pars
    y_sig = y - y_bg

    # Stage 2: per-chunk amplitudes
    if cache is not None:
        n_chunks    = cache['n_chunks']
        chunk_size  = cache['chunk_size']
        chunk_pinvs = cache['chunk_pinvs']
    else:
        n_chunks    = n // chunk_size
        chunk_pinvs = None

    chunk_amps = np.zeros(n_chunks, dtype=complex)
    fit_sig    = np.zeros(n, dtype=float)
    rss_total  = 0.0
    ndf_total  = 0

    for k in range(n_chunks):
        i0 = k * chunk_size
        i1 = i0 + chunk_size
        yc = y_sig[i0:i1]
        if chunk_pinvs is not None:
            pinv_c, wc, Xc = chunk_pinvs[k]
            pars_c = pinv_c @ (yc * wc)
        else:
            idx = np.arange(i0, i1)
            wtc = 2 * np.pi * idx * rf
            Xc  = np.column_stack([np.ones(chunk_size),
                                    np.cos(wtc), np.sin(wtc),
                                    np.cos(wtc / 2), np.sin(wtc / 2)])
            wc = np.hanning(chunk_size) if use_hann else np.ones(chunk_size)
            pars_c, _, _, _ = np.linalg.lstsq(Xc * wc[:, np.newaxis], yc * wc, rcond=None)
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


def build_fit_cache(fsamp, fsig, fline, n, Nhars, use_hann=True, chunk_periods=0):
    """Precompute fit matrices for repeated calls at the same frequency and sample count.

    Builds the full design matrix X, the Hanning window w, and the pseudoinverse
    pinv(X * w[:, None]).  When chunk_periods > 0, also precomputes the per-chunk
    pseudoinverses so Stage-2 fits reduce to matrix–vector multiplies.

    Pass the returned dict as cache= to fit_sine_cplx.
    """
    rf  = fsig  / fsamp
    rlf = fline / fsamp
    wt  = 2 * np.pi * np.arange(n) * rf
    wlf = 2 * np.pi * np.arange(n) * rlf

    cols = [np.ones(n), np.cos(wt), np.sin(wt),
            np.cos(wt / 2), np.sin(wt / 2),
            np.cos(wlf), np.sin(wlf), np.cos(2*wlf), np.sin(2*wlf)]
    for h in range(2, Nhars + 1):
        cols.extend([np.cos(h * wt), np.sin(h * wt)])

    n_nearest = int(round(fsig / fline))
    for nh in range(max(3, n_nearest - 2), n_nearest + 3):
        fn = nh * fline
        if fsamp / n < abs(fn - fsig) < fline / 2:
            wnh = 2 * np.pi * np.arange(n) * fn / fsamp
            cols.extend([np.cos(wnh), np.sin(wnh)])

    X = np.column_stack(cols)
    w = np.hanning(n) if use_hann else np.ones(n)
    bg_pinv = np.linalg.pinv(X * w[:, np.newaxis])

    cache = {'X': X, 'w': w, 'bg_pinv': bg_pinv, 'rf': rf}

    n_chunk_params = 5
    chunk_size = int(round(chunk_periods * fsamp / fsig)) if chunk_periods > 0 else 0
    if chunk_size > n_chunk_params + 1 and n // chunk_size >= 2:
        n_chunks    = n // chunk_size
        chunk_pinvs = []
        for k in range(n_chunks):
            i0  = k * chunk_size
            i1  = i0 + chunk_size
            idx = np.arange(i0, i1)
            wtc = 2 * np.pi * idx * rf
            Xc  = np.column_stack([np.ones(chunk_size),
                                    np.cos(wtc), np.sin(wtc),
                                    np.cos(wtc / 2), np.sin(wtc / 2)])
            wc = np.hanning(chunk_size) if use_hann else np.ones(chunk_size)
            chunk_pinvs.append((np.linalg.pinv(Xc * wc[:, np.newaxis]), wc, Xc))
        cache['chunk_pinvs'] = chunk_pinvs
        cache['chunk_size']  = chunk_size
        cache['n_chunks']    = n_chunks

    return cache


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
        self.fsig, self.fline = get_f(self.data, self.fsamp, self.fsig, self.fline, Nhars=1)
        return self.fsig, self.fline

    def fit(self, cache=None):
        result = fit_sine_cplx(self.data, self.fsamp, self.fsig, self.fline,
                                self.Nhars, True, self.chunk_periods, cache=cache)
        self.Vc        = result[0]
        self.fv        = result[1]
        self.c2        = result[3]
        self.Vc_chunks = result[4] if len(result) > 4 else None

    def strip_raw(self):
        self.data = None
        self.fv   = None


class ThreeChannels:
    def __init__(self, fsig, fsamp, Nhars, ch1, ch2, ch3, V1c, V2c, i, ts,
                 chunk_periods=0, fit_cache=None):
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

        if fit_cache is not None and 'fsig' in fit_cache:
            # Reuse frequency and precomputed matrices from earlier point at same frequency
            fsig_r       = fit_cache['fsig']
            fline        = fit_cache['fline']
            nhars_used   = fit_cache['nhars_used']
            matrix_cache = fit_cache['matrix_cache']
        else:
            # First point: determine frequency and build matrices; use V2 (large signal)
            fsig_r, fline = self.Data[1].findf()
            nhars_used = max(1, min(Nhars, int(0.45 * fsamp / fsig_r)))
            matrix_cache = build_fit_cache(fsamp, fsig_r, fline, len(ch1), nhars_used,
                                           use_hann=True, chunk_periods=chunk_periods)
            if fit_cache is not None:
                fit_cache['fsig']         = fsig_r
                fit_cache['fline']        = fline
                fit_cache['nhars_used']   = nhars_used
                fit_cache['matrix_cache'] = matrix_cache

        for d in self.Data:
            d.Nhars = nhars_used
            d.setf(fsig_r, fline)
            d.fit(cache=matrix_cache)

        # Use V2 (Data[1]) as phase reference: rotate all Vc so V2 is real and positive
        v2_phase = np.exp(1j * np.angle(self.Data[1].Vc))
        for d in self.Data:
            d.Vc = d.Vc / v2_phase

        # Build per-chunk phasor matrix (n_chunks × 3): rotate each chunk by V2's phase
        chunks0 = self.Data[1].Vc_chunks
        if chunks0 is not None and len(chunks0) > 1:
            n_ch = min(len(d.Vc_chunks) for d in self.Data if d.Vc_chunks is not None)
            self.phasor_chunks = np.zeros((n_ch, 3), dtype=complex)
            for k in range(n_ch):
                cf_k = np.exp(-1j * np.angle(chunks0[k]))
                for j in range(3):
                    self.phasor_chunks[k, j] = self.Data[j].Vc_chunks[k] * cf_k

    def strip_raw(self):
        for sd in self.Data:
            sd.strip_raw()


class NPoints:
    def __init__(self, fsig, fsamp, Nhars=1, g2=1, ratio=10, N=8, chunk_periods=0, fit_cache=None):
        self.N             = N
        self.chunk_periods = chunk_periods
        self._fit_cache    = fit_cache if fit_cache is not None else {}
        self._n_filled     = 0
        self.Res = {}
        self.Res['fsig']  = fsig
        self.Res['ratio'] = ratio
        self.Res['fsamp'] = fsamp
        self.Res['Nhars'] = Nhars
        self.Res['gain2'] = g2
        self.ats = -1.0 * np.ones(N)
        self.Res['ts']    = -1.0
        self.Data = np.zeros(N, dtype=object)

    def setPoint(self, i, ch1, ch2, ch3, V1c, V2c, ts):
        self.V1c = V1c
        self.Data[i] = ThreeChannels(
            self.Res['fsig'], self.Res['fsamp'], self.Res['Nhars'],
            ch1, ch2, ch3, V1c, V2c, i, ts,
            chunk_periods=self.chunk_periods,
            fit_cache=self._fit_cache,
        )
        self.ats[i] = ts
        self._n_filled += 1
        if self.is_complete:
            self.Res['ts'] = float(np.mean(self.ats))

    @property
    def is_complete(self):
        return self._n_filled == self.N

    def precalc(self):
        self.raw3 = np.zeros((self.N, 3), dtype=complex)
        self.ctrl = np.zeros((self.N, 2), dtype=complex)
        for i in range(self.N):
            phi = np.angle(self.Data[i].Data[0].Vc)
            cf = np.exp(-1j * phi)
            for j in range(3):
                self.raw3[i, j] = self.Data[i].Data[j].Vc * cf
            self.ctrl[i, 0] = self.Data[i].V1c
            self.ctrl[i, 1] = self.Data[i].V2c
        self.ctrla = self.ctrl

    def _build_fit_etas(self):
        """
        Returns (eta2, eta3) arrays for the ellipse lstsq.

        When per-chunk phasors are available each ellipse point contributes n_chunks
        data rows, giving a denser lstsq fit.  Falls back to N per-point etas when
        chunking was not used.
        """
        has_chunks = all(
            isinstance(self.Data[i], ThreeChannels) and self.Data[i].phasor_chunks is not None
            for i in range(self.N)
        )
        if not has_chunks:
            return self.eta2, self.eta3

        n_ch = min(len(self.Data[i].phasor_chunks) for i in range(self.N))
        if n_ch < 2:
            return self.eta2, self.eta3

        e2, e3 = [], []
        for m in range(self.N):
            pc  = self.Data[m].phasor_chunks[:n_ch]
            v1  = pc[:, 0]
            e2.append(pc[:, 1] / v1)
            e3.append(pc[:, 2] / v1)
        return np.concatenate(e2), np.concatenate(e3)

    def calc(self):
        self.precalc()
        self.V1m = self.raw3[:, 0]
        self.V2m = self.raw3[:, 1]
        self.V3m = self.raw3[:, 2]
        self.eta2 = self.V2m / self.V1m
        self.eta3 = self.V3m / self.V1m

        eta2_fit, eta3_fit = self._build_fit_etas()

        Xmat = np.column_stack([eta2_fit, np.ones(len(eta2_fit))])
        (m4, c4), _, _, _ = np.linalg.lstsq(Xmat, eta3_fit, rcond=None)
        self.gamma3 = 1.0 / m4
        self.Res['Vz3'] = c4 / m4
        self.Res['gamma3'] = self.gamma3
        V4_raw  = np.array([self.Data[k].Data[2].Vc for k in range(self.N)])
        V1_meas = np.array([self.Data[k].Data[0].Vc for k in range(self.N)])
        eta3 = V4_raw / V1_meas
        V1rb_pts = self.ctrla[:, 0]
        V1rb_center = np.mean(V1rb_pts)
        Xmat = np.column_stack([V1rb_pts - V1rb_center, np.ones(self.N)])
        (m_bal, c_bal), _, _, _ = np.linalg.lstsq(Xmat, eta3, rcond=None)
        step = c_bal / m_bal
        dV1_spread = np.std(V1rb_pts)
        step_limit = 5 * dV1_spread
        if np.abs(step) > step_limit:
            step = step * step_limit / np.abs(step)
        self.Res['eta3fit_slope'] = m_bal
        self.Res['eta3fit_intercept'] = c_bal
        self.Res['eta3_mean'] = np.mean(eta3)
        self.Res['V1rb_center'] = V1rb_center
        self.Res['V1_balance'] = V1rb_center - 0.5 * step

        self.combined3 = self.gamma3 * self.eta3 - self.eta2

        ratio = self.Res['ratio']
        self.alpha3 = np.real(self.combined3) / ratio - 1
        self.beta3  = np.imag(self.combined3) / ratio
        self.Res['alpha3mean'] = np.mean(self.alpha3)
        self.Res['beta3mean']  = np.mean(self.beta3)
        self.Res['V1cReadback'] = self.V1c

        self.setGoodFlag()

    def setGoodFlag(self):
        self.goodData = not np.any(~np.isfinite(np.abs(self.combined3)))

    def max_raw_amplitude(self, channel_idx):
        """Return max |raw sample| across all measurement points for the given channel index."""
        mx = 0.0
        for d in self.Data:
            if isinstance(d, ThreeChannels) and len(d.Data) > channel_idx:
                raw = d.Data[channel_idx].data
                if raw is not None:
                    mx = max(mx, float(np.max(np.abs(raw))))
        return mx

    def strip_raw(self):
        for fc in self.Data:
            if isinstance(fc, ThreeChannels):
                fc.strip_raw()


class AllData():
    def __init__(self):
        # keyed by (round(fsig_hz), swpos); swpos 0=straight, 1=cross
        self.mydict = {}

    def _key(self, f, swpos):
        return (round(f), int(swpos))

    def append(self, nd: NPoints):
        key = self._key(nd.Res['fsig'], nd.Res.get('switch_pos', 0))
        if key not in self.mydict:
            self.mydict[key] = []
        self.mydict[key].append(nd)

    def deletekey(self, f):
        """Remove all entries at frequency f (both switch states)."""
        f_hz = round(f)
        for sw in (0, 1):
            self.mydict.pop((f_hz, sw), None)

    def freqs(self):
        """Sorted list of unique nominal frequencies (integer Hz) that have any data."""
        return sorted({k[0] for k in self.mydict})

    def entries(self, f):
        """All NPoints at frequency f: straight state first, then cross."""
        f_hz = round(f)
        result = []
        for sw in (0, 1):
            result.extend(self.mydict.get((f_hz, sw), []))
        return result

    def count(self):
        """Number of distinct frequencies with any data."""
        return len({k[0] for k in self.mydict})

    def countf(self, f):
        """Total NPoints at frequency f across both switch states."""
        return sum(len(self.mydict.get((round(f), sw), [])) for sw in (0, 1))

    def countf_sw(self, f, sw):
        return len(self.mydict.get(self._key(f, sw), []))

    def getkeys(self, f, keys):
        nds = self.entries(f)
        if not nds:
            return {k: np.array([]) for k in keys}
        retdict = {k: np.empty(len(nds), dtype=object) for k in keys}
        for n, nd in enumerate(nds):
            for k in keys:
                retdict[k][n] = nd.Res[k]
        return retdict

    def getallkeys(self, f):
        nds = self.entries(f)
        if not nds:
            return {}
        return self.getkeys(f, list(nds[0].Res.keys()))

    def getdictf(self, keys):
        keys = list(keys)
        if 'fsig' not in keys:
            keys.append('fsig')
        all_freqs = self.freqs()
        if not all_freqs:
            return {k: np.array([]) for k in keys}
        retdict = {k: [] for k in keys}
        for f in all_freqs:
            for k, v in self.getkeys(f, keys).items():
                retdict[k].extend(list(v))
        return {k: np.array(retdict[k]) for k in keys}

    def getRawPhasors(self, f, t0=0):
        COLS = ['frequency/Hz', 't/s',
                'reV1', 'imV1', 'reV2', 'imV2', 'reTZA', 'imTZA',
                'reV1set', 'imV1set', 'reV2set', 'imV2set',
                'gain2', 'fsamp', 'swpos']
        rows = []
        for obj in self.entries(f):
            for j in range(obj.raw3.shape[0]):
                row = [f, obj.Res['ts'] - t0]
                for k in range(obj.raw3.shape[1]):
                    row += [obj.raw3[j, k].real, obj.raw3[j, k].imag]
                for k in range(obj.ctrla.shape[1]):
                    row += [obj.ctrla[j, k].real, obj.ctrla[j, k].imag]
                row += [int(obj.Res['gain2']), int(obj.Res['fsamp']), int(obj.Res.get('switch_pos', 0))]
                rows.append(row)
        return pd.DataFrame(rows, columns=COLS)

    def getEllipseData(self, f):
        """One row per NPoints (one ellipse sweep): f, t, ratio, eta2, eta3, gain, gain2, fsamp, swpos."""
        rows = []
        for obj in self.entries(f):
            if not hasattr(obj, 'combined3'):
                continue
            rows.append({
                'f':     float(obj.Res['fsig']),
                't':     float(obj.Res['ts']),
                'ratio': complex(np.mean(obj.combined3)),
                'eta2':  complex(np.mean(obj.eta2)),
                'eta3':  complex(np.mean(obj.eta3)),
                'gain':  float(abs(obj.gamma3)),
                'gain2': int(obj.Res['gain2']),
                'fsamp': int(obj.Res['fsamp']),
                'swpos': int(obj.Res.get('switch_pos', 0)),
            })
        df = pd.DataFrame(rows)
        if not df.empty:
            for col in ('ratio', 'eta2', 'eta3'):
                df[col] = df[col].astype(complex)
        return df
