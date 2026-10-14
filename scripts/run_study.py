"""Run every simulated measurement of the study and store the smoothed results.

Cases, for microphone distances of 1 m and 2 m, horizontal orbit -170 ... 180 degrees by 10:

* ``rigorous``: turntable axis through the reference point on the baffle;
* ``raw``: cabinet centered on the turntable, microphone moved sideways onto the reference axis;
* ``corrected``: ``raw`` after the rotation-center correction with the true offsets;
* variants of ``corrected`` with wrong or missing correction parameters, a 5-degree orbit, and
  the edge diffraction switched off.

Output: ``results/study.npz`` (1/24-octave smoothed power spectra on a 1/48-octave grid).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from rotcenter.correction import PlaneOffsets, Sample, correct_plane  # noqa: E402
from rotcenter.diffraction import direct_term  # noqa: E402
from rotcenter.geometry import CENTERED, RIGOROUS, Placement, correction_offsets, mic_positions  # noqa: E402
from rotcenter.signal import log_grid, smooth_power  # noqa: E402
from rotcenter.simulate import Loudspeaker  # noqa: E402
from rotcenter.system import design_filters, system_sum  # noqa: E402

ANGLES = np.arange(-170, 181, 10)
FINE_ANGLES = np.arange(-175, 181, 5)
DISTANCES = (1.0, 2.0)
DISPLAY = log_grid(20.0, 20_000.0, 48)
SMOOTHING = 24

#: Correction parameters entered by the user, as changes to the true ones (mm).
VARIANTS = {
    "lateral_omitted": dict(lateral_mm=lambda o: 0.0),
    "lateral_sign": dict(lateral_mm=lambda o: -o["lateral_mm"]),
    "woofer_offset_omitted": dict(driver_offset_mm=lambda o: 0.0),
    "depth_plus10": dict(depth_mm=lambda o: o["depth_mm"] + 10.0),
    "depth_minus10": dict(depth_mm=lambda o: o["depth_mm"] - 10.0),
    "lateral_plus10": dict(lateral_mm=lambda o: o["lateral_mm"] + 10.0),
    "lateral_minus10": dict(lateral_mm=lambda o: o["lateral_mm"] - 10.0),
    "distance_plus20": dict(mic_distance_mm=lambda o: o["mic_distance_mm"] + 20.0),
}


def corrected(ls, raw, angles, distance, change=None):
    """Spectra after the correction, ``{driver: (n_angles, n_freqs)}``."""
    out = {}
    for name, driver in ls.drivers.items():
        offsets = correction_offsets(CENTERED, driver.placement, distance)
        if change:
            for key, fn in change.items():
                if key == "driver_offset_mm" and name != "woofer":
                    continue
                offsets[key] = fn(offsets)
        irs = ls.impulse_responses(raw[name])
        fixed = correct_plane([Sample(float(a), ir) for a, ir in zip(angles, irs)], PlaneOffsets(**offsets), ls.fs)
        rows = [np.fft.rfft(c[: ls.n_fft], ls.n_fft) if c is not None else raw[name][i] for i, c in enumerate(fixed)]
        out[name] = np.stack(rows)
    return out


def direct_only_flat(ls, mics):
    """Direct sound only, every driver moved to the microphone height (no out-of-plane offset)."""
    out = {}
    for name, driver in ls.drivers.items():
        flat = type(driver)(driver.name, driver.label, Placement(driver.placement.x, float(mics[0][1])),
                            driver.radiator, driver.response)
        out[name] = np.stack([driver.response * direct_term(flat, m, ls.k) for m in mics])
    return out


class Store:
    def __init__(self, ls):
        self.ls = ls
        self.data = {"freqs": DISPLAY, "angles": ANGLES, "fine_angles": FINE_ANGLES}

    def power(self, spectra):
        return smooth_power(self.ls.freqs, np.abs(spectra) ** 2, SMOOTHING, DISPLAY)

    def add(self, key, spectra, reference=None):
        """Smoothed power of every driver (and of their difference with ``reference``)."""
        for name, s in spectra.items():
            self.data[f"{key}/{name}/power"] = self.power(s)
            if reference is not None:
                self.data[f"{key}/{name}/error"] = self.power(s - reference[name])

    def save(self, path):
        np.savez_compressed(path, **self.data)


def main():
    t0 = time.time()
    ls = Loudspeaker()
    print(f"loudspeaker ready ({time.time() - t0:.0f} s)")
    store = Store(ls)
    filters = None
    for distance in DISTANCES:
        tag = f"{distance:.0f}m"
        for diffraction in (True, False):
            dtag = "" if diffraction else "-direct"
            rigorous = ls.polar(mic_positions(RIGOROUS, ANGLES, distance), diffraction)
            raw = ls.polar(mic_positions(CENTERED, ANGLES, distance), diffraction)
            fixed = corrected(ls, raw, ANGLES, distance)
            if filters is None:
                # Crossover designed once, on the rigorous on-axis responses at 1 m.
                on_axis = {n: s[list(ANGLES).index(0)] for n, s in rigorous.items()}
                filters = design_filters(ls.freqs, on_axis)
                store.data["filters"] = np.stack([np.abs(filters[n]) for n in ls.drivers])
            cases = {"rigorous": rigorous, "raw": raw, "corrected": fixed}
            for name, spectra in cases.items():
                spectra = dict(spectra, system=system_sum(filters, spectra))
                ref = dict(rigorous, system=system_sum(filters, rigorous))
                store.add(f"{tag}{dtag}/{name}", spectra, ref)
            if diffraction:
                ref = dict(rigorous, system=system_sum(filters, rigorous))
                for variant, change in VARIANTS.items():
                    v = corrected(ls, raw, ANGLES, distance, change)
                    store.add(f"{tag}/{variant}", dict(v, system=system_sum(filters, v)), ref)
                # Time domain, a few angles (impulse responses of the tweeter and of the system).
                for name, spectra in cases.items():
                    for drv in ("tweeter", "mid", "woofer"):
                        store.data[f"{tag}/{name}/{drv}/ir"] = ls.impulse_responses(spectra[drv][[17, 26, 32, 35]])
                    sys_spec = system_sum(filters, spectra)[[17, 26, 32, 35]]
                    store.data[f"{tag}/{name}/system/ir"] = ls.impulse_responses(sys_spec)
            else:
                flat_r = direct_only_flat(ls, mic_positions(RIGOROUS, ANGLES, distance))
                flat_c = direct_only_flat(ls, mic_positions(CENTERED, ANGLES, distance))
                flat_fix = corrected(ls, flat_c, ANGLES, distance)
                store.add(f"{tag}-flat/rigorous", flat_r, flat_r)
                store.add(f"{tag}-flat/corrected", flat_fix, flat_r)
            print(f"{tag}{dtag} done ({time.time() - t0:.0f} s)")
        # Finer orbit (5 degrees), with diffraction: errors at the common 10-degree angles.
        rigorous5 = ls.polar(mic_positions(RIGOROUS, FINE_ANGLES, distance))
        raw5 = ls.polar(mic_positions(CENTERED, FINE_ANGLES, distance))
        fixed5 = corrected(ls, raw5, FINE_ANGLES, distance)
        keep = np.isin(FINE_ANGLES, ANGLES)
        r5 = {n: s[keep] for n, s in rigorous5.items()}
        f5 = {n: s[keep] for n, s in fixed5.items()}
        store.add(f"{tag}/step5", dict(f5, system=system_sum(filters, f5)),
                  dict(r5, system=system_sum(filters, r5)))
        print(f"{tag} 5-degree orbit done ({time.time() - t0:.0f} s)")
    out = ROOT / "results" / "study.npz"
    store.save(out)
    print(f"saved {out} ({out.stat().st_size / 1e6:.1f} MB, {time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
