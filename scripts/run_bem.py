"""The boundary-element check: the two-way variant, BEM against DED, rigorous against corrected.

Reads ``data/bem/two_way.npz`` (pressures of each driver per unit volume velocity), builds the
two-way from the datasheets and computes, at the BEM frequencies (1/48 octave, 20 Hz-16 kHz) and
positions, at every distance of the data, for the BEM and for the DED model of the same loudspeaker:

* on the horizontal orbit, ``rigorous``, ``raw`` and ``corrected`` data with 10, 5 and 2.5 degree
  steps (suffixes ``""``, ``"5"``, ``"2"``; the finer steps where the data have them), each driver
  and the summed loudspeaker; ``levels``: corrected with the powers interpolated instead of the
  complex responses;
* the correction variants of the three-way study (wrong parameters), 10 degree steps;
* the same on the vertical orbit (prefix ``v-``), measured with the loudspeaker on its side
  (``v-raw`` and ``v-corrected`` for the BEM where the data have the centered vertical orbit).

Output: ``results/bem.npz`` (power spectra ``|p|^2``, smoothed to 1/24 octave, single precision).

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

#: Angular steps (degrees) and the suffix of their cases.
STEPS = ((10.0, ""), (5.0, "5"), (2.5, "2"))


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
            self.data[f"{key}/{name}/power"] = smooth_power(self.freqs, np.abs(s) ** 2, 24,
                                                            self.freqs).astype(np.float32)


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

    base = 360.0 / len(bem.angles)
    steps = [(step, suffix) for step, suffix in STEPS if step >= base - 1e-9
             and abs(step / base - round(step / base)) < 1e-9]
    distances = sorted({distance for _, distance in bem.sets.values()})
    for step, suffix in steps:
        store.data[f"angles{suffix}"] = twoway.orbit({}, bem.angles, step)[0]
    for model in ("bem", "ded"):
        for d in distances:
            tag = f"{model}/{d:.0f}m"
            for plane, prefix in (("h", ""), ("v", "v-")):
                rig = models.polar(model, f"rigorous-{plane}", d)
                centered = f"centered-{plane}-{d:.0f}m"
                raw = models.polar(model, f"centered-{plane}", d) if model == "ded" or centered in bem.sets else None
                for step, suffix in steps:
                    angles, r = twoway.orbit(rig, bem.angles, step)
                    store.add(f"{tag}/{prefix}rigorous{suffix}", r, filters)
                    if raw is None:
                        continue
                    _, m = twoway.orbit(raw, bem.angles, step)
                    store.add(f"{tag}/{prefix}raw{suffix}", m, filters)
                    store.add(f"{tag}/{prefix}corrected{suffix}", corrected(m, angles, d, freqs, plane), filters)
                    store.add(f"{tag}/{prefix}levels{suffix}",
                              corrected(m, angles, d, freqs, plane, levels=True), filters)
                    if plane == "h" and step == 10.0:
                        for variant, change in VARIANTS.items():
                            store.add(f"{tag}/{variant}", corrected(m, angles, d, freqs, change=change), filters)
            print(f"{tag} done ({time.time() - t0:.0f} s)")
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "results" / "bem.npz"
    out.parent.mkdir(exist_ok=True)
    np.savez_compressed(out, **store.data)
    print(f"saved {out} ({out.stat().st_size / 1e6:.1f} MB, {time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
