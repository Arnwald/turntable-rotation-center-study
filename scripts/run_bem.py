"""The boundary-element check: the two-way variant, BEM against DED, rigorous against corrected.

Reads ``data/bem/two_way.npz`` (pressures of each driver per unit volume velocity), builds the
two-way from the datasheets and computes, at the BEM frequencies (1/48 octave, 20 Hz-16 kHz, no
smoothing) and positions, for 1 m and 2 m:

* ``rigorous``, ``raw`` and ``corrected`` data on the horizontal orbit, 10 and 5 degree steps
  (``*5``), each driver and the summed loudspeaker, for the BEM and for the DED model of the same
  loudspeaker; ``levels``: corrected with the powers interpolated instead of the complex responses;
* the correction variants of the three-way study (wrong parameters);
* the vertical orbit: ``v-rigorous`` for both models and, for the DED model only (the BEM was not
  solved there), ``v-raw`` and ``v-corrected`` with the loudspeaker on its side on the turntable.

Output: ``results/bem.npz`` (power spectra ``|p|^2``, smoothed to 1/24 octave).

Usage: ``python3 scripts/run_bem.py [DATA.npz [OUT.npz]]``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from rotcenter import twoway  # noqa: E402
from rotcenter.correction import PlaneOffsets, correct_spectra  # noqa: E402
from rotcenter.signal import smooth_power  # noqa: E402
from rotcenter.simulate import Loudspeaker  # noqa: E402
from run_study import VARIANTS  # noqa: E402

DISTANCES = (1.0, 2.0)


class Models:
    """The BEM and DED pressures of the two-way at the BEM frequencies."""

    def __init__(self, bem: twoway.BemData):
        self.bem = bem
        self.freqs = bem.freqs
        self.drivers = twoway.drivers(self.freqs)
        self.velocity = {n: twoway.volume_velocity(d, self.freqs) for n, d in self.drivers.items()}
        self.ded = Loudspeaker(freqs=self.freqs, drivers=self.drivers)

    def polar(self, model: str, setup: str, distance: float) -> dict:
        name = f"{setup}-{distance:.0f}m"
        if model == "bem":
            return self.bem.polar(name, self.velocity)
        mics = (self.bem.positions(name) if name in self.bem.sets
                else twoway.microphones(setup, distance, self.bem.angles))
        return self.ded.polar(mics)


def corrected(raw: dict, angles, distance, freqs, plane="h", change=None, levels=False) -> dict:
    """Each driver corrected with its own parameters (``change``: a variant's wrong ones;
    ``levels``: interpolating powers)."""
    out = {}
    for name, placement in twoway.PLACEMENTS.items():
        offsets = twoway.correction_offsets(placement, distance, plane)
        for key, fn in (change or {}).items():
            if key == "driver_offset_mm" and name != "woofer":
                continue
            offsets[key] = fn(offsets)
        out[name] = correct_spectra(raw[name], angles, PlaneOffsets(**offsets), freqs, levels=levels)
    return out


class Store:
    """Power spectra smoothed to 1/24 octave as in the three-way study: on the BEM grid (1/48
    octave), the mean of each frequency and its two neighbours."""

    def __init__(self, freqs):
        self.freqs = freqs
        self.data = {}

    def add(self, key, spectra, filters):
        spectra = dict(spectra, system=twoway.system_sum(filters, spectra))
        for name, s in spectra.items():
            self.data[f"{key}/{name}/power"] = smooth_power(self.freqs, np.abs(s) ** 2, 24, self.freqs)


def main():
    t0 = time.time()
    bem = twoway.BemData(Path(sys.argv[1])) if len(sys.argv) > 1 else twoway.BemData()
    freqs = bem.freqs
    models = Models(bem)
    print(f"models ready ({time.time() - t0:.0f} s)")
    store = Store(freqs)
    on_axis = models.polar("bem", "rigorous-h", 1.0)
    index0 = int(np.flatnonzero(bem.angles == 0.0)[0])
    filters = twoway.design_filters(freqs, {n: s[index0] for n, s in on_axis.items()})
    store.data.update(freqs=freqs, filters=np.stack([filters[n] for n in ("woofer", "tweeter")]),
                      velocity=np.stack([models.velocity[n] for n in ("woofer", "tweeter")]),
                      response=np.stack([models.drivers[n].response for n in ("woofer", "tweeter")]))
    for key, value in bem.solves.items():
        store.data[f"solve/{key}"] = value

    for model in ("bem", "ded"):
        for d in DISTANCES:
            tag = f"{model}/{d:.0f}m"
            rig = models.polar(model, "rigorous-h", d)
            raw = models.polar(model, "centered-h", d)
            for step, suffix in ((10, ""), (5, "5")):
                angles, r = twoway.orbit(rig, bem.angles, step)
                _, m = twoway.orbit(raw, bem.angles, step)
                store.data[f"angles{suffix}"] = angles
                fixed = corrected(m, angles, d, freqs)
                store.add(f"{tag}/rigorous{suffix}", r, filters)
                store.add(f"{tag}/raw{suffix}", m, filters)
                store.add(f"{tag}/corrected{suffix}", fixed, filters)
                store.add(f"{tag}/levels{suffix}", corrected(m, angles, d, freqs, levels=True), filters)
                if step == 10:
                    for variant, change in VARIANTS.items():
                        store.add(f"{tag}/{variant}", corrected(m, angles, d, freqs, change=change), filters)
            # Vertical orbit.
            angles, rv = twoway.orbit(models.polar(model, "rigorous-v", d), bem.angles, 10)
            store.add(f"{tag}/v-rigorous", rv, filters)
            if model == "ded":
                _, mv = twoway.orbit(models.polar(model, "centered-v", d), bem.angles, 10)
                store.add(f"{tag}/v-raw", mv, filters)
                store.add(f"{tag}/v-corrected", corrected(mv, angles, d, freqs, plane="v"), filters)
            print(f"{tag} done ({time.time() - t0:.0f} s)")
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "results" / "bem.npz"
    out.parent.mkdir(exist_ok=True)
    np.savez_compressed(out, **store.data)
    print(f"saved {out} ({out.stat().st_size / 1e6:.1f} MB, {time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
