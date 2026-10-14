"""Error metrics shared by the figures and tables."""

from __future__ import annotations

import numpy as np

#: Evaluation band of every driver (Hz): its passband with the illustrative crossover, plus
#: about one octave of overlap.
BANDS = {"woofer": (50.0, 1000.0), "mid": (300.0, 6000.0), "tweeter": (2000.0, 20000.0),
         "system": (50.0, 20000.0)}

#: Directions where the rigorous response is more than this below the on-axis response at the
#: same frequency are left out of the statistics (deep nulls, the back of a dome at 15 kHz...).
LEVEL_FLOOR_DB = -20.0


def level_error(power, reference):
    """``10 log10(P / P_ref)``, the difference two plotted curves would show (dB)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return 10.0 * np.log10(power / reference)


def relevant(reference, freqs, angles, band, floor_db=LEVEL_FLOOR_DB):
    """Mask (angles x freqs) of the points that count: in band, and no more than ``floor_db``
    below the on-axis response."""
    on_axis = reference[list(angles).index(0)]
    in_band = (freqs >= band[0]) & (freqs <= band[1])
    return (reference >= on_axis[None, :] * 10.0 ** (floor_db / 10.0)) & in_band[None, :]


def summary(power, reference, freqs, angles, band, floor_db=LEVEL_FLOOR_DB):
    """RMS, 95th percentile and maximum of ``|level error|`` over the relevant points (dB)."""
    mask = relevant(reference, freqs, angles, band, floor_db)
    e = np.abs(level_error(power, reference))[mask]
    return {"rms": float(np.sqrt(np.mean(e**2))), "p95": float(np.percentile(e, 95)), "max": float(e.max())}


def per_angle_rms(power, reference, freqs, angles, band, floor_db=LEVEL_FLOOR_DB):
    """RMS level error over the relevant frequencies of each angle (dB)."""
    mask = relevant(reference, freqs, angles, band, floor_db)
    e = level_error(power, reference)
    out = np.full(len(angles), np.nan)
    for i in range(len(angles)):
        if mask[i].any():
            out[i] = np.sqrt(np.mean(e[i, mask[i]] ** 2))
    return out
