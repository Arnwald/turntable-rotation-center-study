"""An illustrative three-way crossover, to look at the summed loudspeaker.

Each driver is equalized to an acoustic target on axis at 1 m (rigorous setup):

    T_woofer  = LR4 low-pass at 450 Hz * 2nd-order Butterworth high-pass at 40 Hz
    T_mid     = LR4 high-pass at 450 Hz * LR4 low-pass at 3 kHz
    T_tweeter = LR4 high-pass at 3 kHz

The filter of driver ``i`` is the minimum-phase filter with magnitude ``L |T_i| / S_i``, where
``S_i`` is its 1/6-octave smoothed on-axis magnitude and ``L`` the target level (86 dB SPL at
1 m). Since ``T_i`` is itself minimum phase, the filter is ``L T_i`` times the minimum-phase
inverse of the smoothed driver: on axis the sum is close to the Linkwitz-Riley all-pass. The
same filters are used for every case and distance: only the measurements change.
"""

from __future__ import annotations

import numpy as np

from .constants import P_REF
from .signal import butter2, lr4, minimum_phase, smooth_power

CROSSOVER_HZ = (450.0, 3000.0)
TARGET_DB = 86.0


def targets(freqs, crossovers=CROSSOVER_HZ):
    lo, hi = crossovers
    f = np.maximum(freqs, 1e-3)
    return {
        "woofer": lr4(f, lo, "lp") * butter2(f, 40.0, "hp"),
        "mid": lr4(f, lo, "hp") * lr4(f, hi, "lp"),
        "tweeter": lr4(f, hi, "hp"),
    }


def design_filters(freqs, on_axis: dict[str, np.ndarray], crossovers=CROSSOVER_HZ, level_db=TARGET_DB):
    """Filters ``{driver: complex response}`` from the on-axis driver spectra at 1 m."""
    level = P_REF * 10.0 ** (level_db / 20.0)
    out = {}
    for name, target in targets(freqs, crossovers).items():
        f_eval = np.maximum(freqs, freqs[1])
        smooth = np.sqrt(smooth_power(freqs, np.abs(on_axis[name]) ** 2, 6, f_eval))
        smooth = np.maximum(smooth, 1e-6 * smooth.max())
        out[name] = minimum_phase(level * np.abs(target) / smooth)
        out[name][0] = 0.0
    return out


def system_sum(filters: dict, spectra: dict) -> np.ndarray:
    """Acoustic sum ``sum_i X_i p_i`` (works on arrays of any leading shape)."""
    return sum(filters[name] * spectra[name] for name in filters)
