"""Frequency grids, minimum-phase reconstruction, smoothing and filters."""

from __future__ import annotations

import numpy as np

from .constants import FS, N_FFT


def rfft_freqs(n: int = N_FFT, fs: float = FS) -> np.ndarray:
    """Bin frequencies of a real FFT of length ``n``."""
    return np.arange(n // 2 + 1) * fs / n


def log_grid(f_lo: float = 20.0, f_hi: float = 20_000.0, per_octave: int = 48) -> np.ndarray:
    """Logarithmically spaced display grid."""
    n = int(np.ceil(np.log2(f_hi / f_lo) * per_octave)) + 1
    return f_lo * 2.0 ** (np.arange(n) / per_octave)


def interp_db(f_src, db_src, f_dst, below=None, above=None):
    """Interpolate a dB curve in log-frequency.

    ``below``/``above`` give the slope (dB/octave) used to extend the curve outside the data,
    anchored at its first/last point; without them the end values are held.
    """
    f_src = np.asarray(f_src, float)
    db_src = np.asarray(db_src, float)
    f_dst = np.asarray(f_dst, float)
    out = np.interp(np.log2(np.maximum(f_dst, 1e-9)), np.log2(f_src), db_src)
    if below is not None:
        m = f_dst < f_src[0]
        out[m] = db_src[0] + below * np.log2(np.maximum(f_dst[m], 1e-9) / f_src[0])
    if above is not None:
        m = f_dst > f_src[-1]
        out[m] = db_src[-1] + above * np.log2(f_dst[m] / f_src[-1])
    return out


def minimum_phase(magnitude: np.ndarray) -> np.ndarray:
    """Minimum-phase spectrum with the given magnitude (real-FFT bins of an even-length grid).

    Homomorphic method: the real cepstrum of ``ln|H|`` is folded onto positive quefrencies,
    ``c_min[0] = c[0]``, ``c_min[n] = 2 c[n]`` for ``0 < n < N/2``, ``c_min[N/2] = c[N/2]``,
    and ``H_min = exp(FFT(c_min))``.
    """
    mag = np.maximum(np.asarray(magnitude, float), 1e-12 * np.max(magnitude))
    n = 2 * (len(mag) - 1)
    cep = np.fft.irfft(np.log(mag), n)
    fold = np.zeros(n)
    fold[0] = cep[0]
    fold[1 : n // 2] = 2.0 * cep[1 : n // 2]
    fold[n // 2] = cep[n // 2]
    return np.exp(np.fft.rfft(fold))


def minimum_phase_from_db(f_src, db_src, f_grid, below, above, oversample: int = 8):
    """Minimum-phase response on ``f_grid`` (an rfft grid) from a dB curve.

    The magnitude is interpolated on an ``oversample`` times finer grid first, so that the
    cepstrum does not alias, then decimated back to ``f_grid``.
    """
    n = 2 * (len(f_grid) - 1)
    fs = 2 * f_grid[-1]
    fine = rfft_freqs(n * oversample, fs)
    db = interp_db(f_src, db_src, fine, below=below, above=above)
    h = minimum_phase(10.0 ** (db / 20.0))
    return h[::oversample]


def smooth_power(freqs, power, fraction: float, out_freqs):
    """Fractional-octave smoothing of a power spectrum (rectangular window in log frequency).

    ``out[i]`` is the mean of ``power`` over the bins within ``out_freqs[i] * 2**(+-1/(2*fraction))``
    (``fraction = 12`` for 1/12 octave). ``power`` may be 2-D with frequency on the last axis.
    """
    freqs = np.asarray(freqs, float)
    power = np.asarray(power, float)
    out_freqs = np.asarray(out_freqs, float)
    half = 2.0 ** (1.0 / (2.0 * fraction))
    lo = np.searchsorted(freqs, out_freqs / half, side="left")
    hi = np.searchsorted(freqs, out_freqs * half, side="right")
    # At least one bin: the nearest one.
    nearest = np.clip(np.searchsorted(freqs, out_freqs), 0, len(freqs) - 1)
    empty = hi <= lo
    lo = np.where(empty, nearest, lo)
    hi = np.where(empty, nearest + 1, hi)
    csum = np.concatenate([np.zeros(power.shape[:-1] + (1,)), np.cumsum(power, axis=-1)], axis=-1)
    return (csum[..., hi] - csum[..., lo]) / (hi - lo)


def butter2(f, fc, kind):
    """Second-order Butterworth section evaluated at frequencies ``f``."""
    s = 1j * np.asarray(f, float) / fc
    den = s**2 + np.sqrt(2.0) * s + 1.0
    return (1.0 if kind == "lp" else s**2) / den


def lr4(f, fc, kind):
    """Fourth-order Linkwitz-Riley low-pass or high-pass (a squared Butterworth section)."""
    return butter2(f, fc, kind) ** 2


def raised_cosine_taper(f, f_start, f_stop):
    """1 below ``f_start``, 0 above ``f_stop``, raised cosine in between."""
    f = np.asarray(f, float)
    t = np.clip((f - f_start) / (f_stop - f_start), 0.0, 1.0)
    return 0.5 * (1.0 + np.cos(np.pi * t))


def db(x, floor=1e-20):
    """20 log10 |x|."""
    return 20.0 * np.log10(np.maximum(np.abs(x), floor))


def power_db(p, floor=1e-40):
    """10 log10 p."""
    return 10.0 * np.log10(np.maximum(p, floor))
