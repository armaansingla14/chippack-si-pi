"""Channel pulse response from S-parameters and eye / jitter analysis from it.

Uses the standard approach of IBIS-AMI style statistical analysis: the
frequency response (with TX/RX termination effects already in S21) is
converted to an impulse, convolved with a 1-UI pulse, and the resulting
cursors give a peak-distortion eye and a statistical eye at a target BER.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import erfcinv


def extrapolate_to_dc(f: np.ndarray, h: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Uniform grid from DC: magnitude held at the lowest point, phase linearly extrapolated."""
    df = f[1] - f[0]
    n0 = int(round(f[0] / df))
    if n0 == 0:
        return f, h
    fl = np.arange(n0) * df
    ph = np.unwrap(np.angle(h))
    slope = (ph[1] - ph[0]) / df
    phl = ph[0] + slope * (fl - f[0])
    hl = np.abs(h[0]) * np.exp(1j * phl)
    hl[0] = np.abs(h[0])
    return np.r_[fl, f], np.r_[hl, h]


def impulse_response(f: np.ndarray, h: np.ndarray, oversample: int = 4, window: bool = True):
    """Real impulse response (dt, h_t) from a one-sided response on a uniform grid."""
    f2, h2 = extrapolate_to_dc(f, h)
    if window:
        h2 = h2 * np.hanning(2 * h2.size)[h2.size:] ** 0.5
    n = oversample * h2.size
    dt = 1.0 / (2 * n * (f2[1] - f2[0]))
    ht = np.fft.irfft(h2, 2 * n) * (2 * n) * (f2[1] - f2[0])
    return dt, ht


def pulse_response(f: np.ndarray, h: np.ndarray, bit_rate: float, spu: int = 32, rise_ui: float = 0.0):
    """Response to a 1-UI rectangular pulse sampled at ``spu`` points per UI.

    Returns (t, p) with t in seconds. The pulse is formed in the frequency domain
    (sinc) so arbitrary rates work on any frequency grid.
    """
    ui = 1.0 / bit_rate
    f2, h2 = extrapolate_to_dc(f, h)
    df = f2[1] - f2[0]
    fs = spu / ui
    n_t = int(np.ceil(fs / df))
    fgrid = np.fft.rfftfreq(n_t, 1 / fs)
    hi = np.interp(fgrid, f2, h2.real, right=0) + 1j * np.interp(fgrid, f2, h2.imag, right=0)
    hi *= np.hanning(2 * fgrid.size)[fgrid.size:] ** 0.5  # tame truncation ringing
    pulse = ui * np.sinc(fgrid * ui) * np.exp(-1j * np.pi * fgrid * ui)
    if rise_ui > 0:
        sig = rise_ui * ui / 1.683
        pulse = pulse * np.exp(-0.5 * (2 * np.pi * fgrid * sig) ** 2)
    p = np.fft.irfft(hi * pulse, n_t) * fs
    # move 8 UI of 'negative time' (wrapped pre-cursor ringing) in front of the pulse
    p = np.roll(p, 8 * spu)
    return (np.arange(n_t) - 8 * spu) / fs, p


@dataclass
class EyeResult:
    peak_distortion_height: float  # worst-case (all ISI aligned) eye height at best phase, V per V swing
    stat_height: float  # eye height at target BER from the statistical eye
    stat_width_ui: float  # horizontal opening at the target BER where height > 0
    main_cursor: float
    isi_sum: float
    best_phase: int
    heights_by_phase: np.ndarray


def cursors_at_phase(p: np.ndarray, spu: int, peak: int, phase_offset: int, n_pre: int = 4, n_post: int = 60):
    """Cursors at UI spacing; indexing is circular because the FFT-based pulse is periodic."""
    idx = (peak + phase_offset + spu * np.arange(-n_pre, n_post + 1)) % p.size
    return p[idx], n_pre


def _isi_pdf(isi: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """PDF of sum_k c_k b_k with b_k = +/-1 equiprobable, on a uniform grid."""
    dx = grid[1] - grid[0]
    pdf = np.zeros(grid.size)
    pdf[np.argmin(np.abs(grid))] = 1.0
    for c in isi:
        sh = int(round(c / dx))
        if sh == 0:
            continue
        a = np.roll(pdf, sh)
        b = np.roll(pdf, -sh)
        pdf = 0.5 * (a + b)
    return pdf


def analyze_eye(p: np.ndarray, spu: int, noise_rms: float = 0.0, ber: float = 1e-12, n_bins: int = 4001,
                swing: float = 1.0) -> EyeResult:
    """NRZ eye from a pulse response (amplitude per unit swing)."""
    peak = int(np.argmax(p))
    q = np.sqrt(2) * erfcinv(2 * ber)
    heights = np.empty(spu)
    pd_heights = np.empty(spu)
    main0 = p[peak]
    span = 2 * np.sum(np.abs(p[::max(spu // 4, 1)])) + 10 * noise_rms + 1e-9
    grid = np.linspace(-span, span, n_bins)
    dx = grid[1] - grid[0]
    for k, off in enumerate(range(-spu // 2, spu - spu // 2)):
        c, i0 = cursors_at_phase(p, spu, peak, off)
        main = c[i0] * swing
        isi = np.delete(c, i0) * swing
        pd_heights[k] = 2 * (main - np.sum(np.abs(isi)))
        pdf = _isi_pdf(isi, grid)
        if noise_rms > 0:
            kk = np.arange(-int(7 * noise_rms / dx) - 1, int(7 * noise_rms / dx) + 2) * dx
            g = np.exp(-0.5 * (kk / noise_rms) ** 2)
            pdf = np.convolve(pdf, g / g.sum(), mode="same")
        cdf = np.cumsum(pdf)
        # lower edge of the '1' distribution at target BER (prob. of falling below)
        target = ber
        lo1 = np.interp(target, cdf, grid) + main if cdf[0] < target else grid[0] + main
        hi0 = -lo1  # symmetric for +/-1 signalling
        heights[k] = lo1 - hi0
    del q
    best = int(np.argmax(heights))
    open_ = heights > 0
    width = 0
    if open_[best]:
        width = 1
        for d in (1, -1):
            j = best
            while 0 <= j + d < spu and open_[j + d]:
                j += d
                width += 1
    c, i0 = cursors_at_phase(p, spu, peak, best - spu // 2)
    return EyeResult(float(pd_heights[best]), float(heights[best]), width / spu, float(c[i0]),
                     float(np.sum(np.abs(np.delete(c, i0)))), best, heights)


def ddj_from_pulse(p: np.ndarray, spu: int, n_hist: int = 8) -> float:
    """Data-dependent jitter (UI pk-pk) from zero crossings of all edge patterns.

    Builds the step response as cumulative sum of the pulse at UI spacing and
    evaluates crossing times of a 0->1 transition preceded by every ``n_hist``-bit
    history (exhaustive over 2**n_hist patterns).
    """
    peak = int(np.argmax(p))
    start = peak - spu  # region where the transition of interest happens
    seg = slice(start - spu // 2, start + 2 * spu)
    times = []
    for pat in range(2 ** n_hist):
        bits = [(pat >> i) & 1 for i in range(n_hist)]  # bits[i] = bit at i+1 UI before the edge bit 0
        # signal = sum_k a_k p(t - kT); current UI k=0 is '1', previous is '0', older from pattern
        a = np.r_[1.0, -1.0, 2.0 * np.array(bits) - 1.0]
        y = np.zeros(p.size)
        for k, ak in enumerate(a):
            sh = k * spu
            y[:p.size - sh] += ak * p[sh:]
        # add the far-future bits as alternating mean (0) -> ignored
        w = y[seg]
        s = np.nonzero((w[:-1] < 0) & (w[1:] >= 0))[0]
        if s.size:
            i = s[-1]
            times.append(i + w[i] / (w[i] - w[i + 1]))
    if not times:
        return float("nan")
    return float((max(times) - min(times)) / spu)
