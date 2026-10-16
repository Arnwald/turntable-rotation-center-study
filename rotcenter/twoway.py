"""The two-way variant, solved by the boundary-element method (BEM) and by the DED model.

Same cabinet as the three-way (440 x 290 x 290 mm, 24 L closed box) with two drivers on its
horizontal middle plane: the woofer (Dayton Audio SIG225-4) 85 mm left of the middle, the tweeter
(SICA LP 90.28/N92) 120 mm right of it, both at mid-height. The reference point is the tweeter
center, the microphone is at its height: every driver lies in the plane of the horizontal orbit.

Radiators, the same in both models:

* woofer: rigid recessed cone of radius ``sqrt(Sd / pi) = 83.8 mm``, 35 mm deep at ``r = 30 mm`` and
  40 mm at the center (the dust cap), its rim flush with the baffle;
* tweeter: the 28 mm dome of the three-way, 7.5 mm high, its rim flush with the baffle.

The BEM data (``data/bem/two_way.npz``) give, for each driver moving as a rigid body with a volume
velocity of 1 m^3/s, the complex pressure (Pa, time factor ``exp(+j omega t)``, the convention of
NumPy's FFT) at every microphone position, 20 Hz to 16 kHz. The drivers' volume velocities come from
the datasheets: ``U(f) = H(f) / q(f)``, where ``H`` is the on-axis half-space response of
:mod:`rotcenter.drivers` and ``q = j omega rho A(0, k) / (2 pi r_ref)`` the on-axis far-field pressure
of the same radiator in an infinite baffle per unit volume velocity. In an infinite baffle the BEM
pressure would then be ``H D(theta) (r_ref / r) exp(-j k r)`` in the far field: the direct term of
the DED model, so that both models share the same drivers and differ only by the radiation and the
diffraction.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np

from .constants import C, FS, P_REF, R_REF, RHO
from .directivity import Radiator
from .drivers import TWEETER_RADIATOR, Driver, tweeter_response, woofer_response
from .geometry import CENTERED, RIGOROUS, Placement, Turntable, mic_positions
from .signal import butter2, lr4, minimum_phase, rfft_freqs, smooth_power

DATA = Path(__file__).resolve().parent.parent / "data" / "bem" / "two_way.npz"

WOOFER = Placement(-0.085, 0.0)
TWEETER = Placement(0.120, 0.0)
PLACEMENTS = {"woofer": WOOFER, "tweeter": TWEETER}
#: The reference point: the tweeter center (the rotation axis of the rigorous setup).
REFERENCE = TWEETER

#: The cone's profile ``(r, z)``, m, ``z`` from the rim plane.
CONE_PROFILE = ((0.0838, 0.0), (0.030, -0.035), (0.0, -0.040))
WOOFER_RADIATOR = dict(radius=0.0838, profile=CONE_PROFILE)

#: Illustrative crossover: Linkwitz-Riley 4th order at 2 kHz, woofer high-passed at 40 Hz.
CROSSOVER_HZ = 2000.0
TARGET_DB = 86.0

#: Evaluation bands (Hz): one octave of overlap around the crossover.
BANDS = {"woofer": (50.0, 4000.0), "tweeter": (1000.0, 16000.0), "system": (50.0, 16000.0)}

#: Grid on which the drivers' minimum-phase responses are built before being sampled at any
#: frequency (0.18 Hz).
FINE_N = 1 << 18
#: Grid of the filters' minimum-phase equalizers: 192 kHz, so that the extension beyond the BEM
#: grid (16 kHz) can run two octaves before the magnitude levels off (0.18 Hz).
EQ_FS = 192_000.0
EQ_N = 1 << 20

#: The vertical orbit in the centered setup, the loudspeaker on its side: the axis runs through the
#: cabinet's center, along x; the tweeter (and the woofer) lie on the plane through that axis.
CENTERED_VERTICAL = Turntable("centered-vertical", 0.0, -0.145)


def sample_response(f_src, h, f_dst):
    """A smooth complex response sampled at ``f_dst``: log-magnitude and unwrapped phase,
    interpolated linearly (``f_src`` an rfft grid, DC left out)."""
    f_src, h = np.asarray(f_src)[1:], np.asarray(h)[1:]
    log_mag = np.interp(f_dst, f_src, np.log(np.maximum(np.abs(h), 1e-300)))
    phase = np.interp(f_dst, f_src, np.unwrap(np.angle(h)))
    return np.exp(log_mag + 1j * phase)


def minimum_phase_at(f_src, db_src, f_dst, below=0.0, above=0.0, reach=2.0):
    """Minimum-phase response with the dB magnitude ``db_src`` (given at ``f_src``) at ``f_dst``.

    Outside ``f_src`` the magnitude follows the slopes ``below`` and ``above`` (dB/octave) for
    ``reach`` octaves, then stays level."""
    fine = rfft_freqs(EQ_N, EQ_FS)
    octaves = np.log2(np.maximum(fine, 1e-9) / f_src[0])
    top = np.log2(f_src[-1] / f_src[0])
    db = np.interp(octaves, np.log2(np.asarray(f_src) / f_src[0]), db_src)
    db = np.where(octaves < 0, db_src[0] + below * np.maximum(octaves, -reach), db)
    db = np.where(octaves > top, db_src[-1] + above * np.minimum(octaves - top, reach), db)
    return sample_response(fine, minimum_phase(10.0 ** (db / 20.0)), f_dst)


def radiators():
    return {"woofer": Radiator(**WOOFER_RADIATOR), "tweeter": Radiator(**TWEETER_RADIATOR)}


def drivers(freqs) -> dict[str, Driver]:
    """The two drivers with their datasheet responses (half space, 1 m, 2.83 V) at ``freqs``."""
    fine = rfft_freqs(FINE_N, FS)
    rad = radiators()
    return {
        "woofer": Driver("woofer", "SIG225-4 (8 in cone)", WOOFER, rad["woofer"],
                         sample_response(fine, woofer_response(fine), freqs)),
        "tweeter": Driver("tweeter", "LP 90.28/N92 (1.1 in dome)", TWEETER, rad["tweeter"],
                          sample_response(fine, tweeter_response(fine), freqs)),
    }


def volume_velocity(driver: Driver, freqs) -> np.ndarray:
    """``U = H / q`` (m^3/s): the volume velocity that gives the datasheet response on axis in an
    infinite baffle (far field, same radiator)."""
    k = 2.0 * np.pi * np.asarray(freqs, float) / C
    q = 1j * 2.0 * np.pi * np.asarray(freqs, float) * RHO * driver.radiator.on_axis(k) / (2.0 * np.pi * R_REF)
    return driver.response / q


class BemData:
    """The boundary-element pressures and the microphone sets they were computed at."""

    def __init__(self, path: Path = DATA):
        d = np.load(path)
        self.freqs = d["freqs"]
        self.points = d["points"]
        self.sources = [str(s) for s in d["sources"]]
        self.transfer = {name: d["pressure"][i].astype(complex) for i, name in enumerate(self.sources)}
        self.angles = d["angles"]
        self.sets = {str(n): (int(f), float(dist)) for n, f, dist in zip(d["set_names"], d["set_first"],
                                                                          d["set_distance"])}
        self.solves = {k[len("solve_"):]: d[k] for k in d.files if k.startswith("solve_")}

    def positions(self, name: str) -> np.ndarray:
        first, _ = self.sets[name]
        return self.points[first : first + len(self.angles)]

    def polar(self, name: str, velocities: dict) -> dict[str, np.ndarray]:
        """Pressures of each driver on one set, ``{driver: (n_angles, n_freqs)}``."""
        first, _ = self.sets[name]
        span = slice(first, first + len(self.angles))
        return {n: velocities[n][None, :] * self.transfer[n][:, span].T for n in velocities}


def orbit(spectra: dict, angles, step: float):
    """The orbit every ``step`` degrees, from ``-180 + step`` to 180 (180 is the -180 point)."""
    angles = np.asarray(angles, float)
    wanted = -180.0 + step * np.arange(1, int(round(360.0 / step)) + 1)
    index = [int(np.flatnonzero(np.isclose(angles, a if a != 180 else -180.0))[0]) for a in wanted]
    return wanted, {n: s[index] for n, s in spectra.items()}


def targets(freqs, crossover=CROSSOVER_HZ):
    f = np.maximum(np.asarray(freqs, float), 1e-3)
    return {"woofer": lr4(f, crossover, "lp") * butter2(f, 40.0, "hp"), "tweeter": lr4(f, crossover, "hp")}


def _end_slopes(freqs, db, span_octaves=1.0 / 6.0):
    """Slopes (dB/octave) of a dB curve over its first and last ``span_octaves``."""
    octaves = np.log2(np.asarray(freqs, float) / freqs[0])
    lo = octaves <= span_octaves
    hi = octaves >= octaves[-1] - span_octaves
    return np.polyfit(octaves[lo], db[lo], 1)[0], np.polyfit(octaves[hi], db[hi], 1)[0]


def design_filters(freqs, on_axis: dict, crossover=CROSSOVER_HZ, level_db=TARGET_DB):
    """Filters ``L T_i E_i`` that bring each driver's on-axis response to its acoustic target, as for
    the three-way, on any frequency grid: ``T_i`` the target (minimum phase, exact), ``E_i`` the
    minimum-phase inverse of the driver's 1/6-octave smoothed on-axis magnitude, extended beyond
    the grid with its end slopes so that its phase is right up to the ends."""
    level = P_REF * 10.0 ** (level_db / 20.0)
    out = {}
    for name, target in targets(freqs, crossover).items():
        smooth = np.sqrt(smooth_power(freqs, np.abs(on_axis[name]) ** 2, 6, freqs))
        inverse_db = -20.0 * np.log10(np.maximum(smooth, 1e-12 * smooth.max()))
        offset = inverse_db.max()
        below, above = _end_slopes(freqs, inverse_db)
        inverse = minimum_phase_at(freqs, inverse_db - offset, freqs, below=below, above=above)
        out[name] = level * target * inverse * 10.0 ** (offset / 20.0)
    return out


def system_sum(filters: dict, spectra: dict) -> np.ndarray:
    return sum(filters[n] * spectra[n] for n in filters)


def microphones(setup: str, distance: float, angles) -> np.ndarray:
    """Microphone positions of the four orbits of the two-way (speaker frame, m)."""
    angles = np.asarray(angles, float)
    if setup == "rigorous-h":
        return mic_positions(RIGOROUS, angles, distance, height=0.0)
    if setup == "centered-h":
        return mic_positions(CENTERED, angles, distance, height=0.0)
    t = np.radians(angles)
    if setup == "rigorous-v":
        y, z = distance * np.sin(t), distance * np.cos(t)
    elif setup == "centered-v":
        a = CENTERED_VERTICAL.axis_z
        y, z = (distance - a) * np.sin(t), a + (distance - a) * np.cos(t)
    else:
        raise ValueError(setup)
    return np.stack([np.full_like(t, REFERENCE.x), y, z], axis=-1)


def correction_offsets(driver: Placement, distance: float, plane: str = "h") -> dict:
    """Parameters of the correction for the centered setup (mm): horizontal orbit (``"h"``), or
    vertical orbit with the loudspeaker on its side (``"v"``), whose axis is level with both
    drivers (no lateral offset, no driver offset)."""
    if plane == "h":
        return {"depth_mm": -CENTERED.axis_z * 1e3, "lateral_mm": (REFERENCE.x - CENTERED.axis_x) * 1e3,
                "driver_offset_mm": (driver.x - REFERENCE.x) * 1e3, "mic_distance_mm": distance * 1e3}
    return {"depth_mm": -CENTERED_VERTICAL.axis_z * 1e3, "lateral_mm": 0.0,
            "driver_offset_mm": (driver.y - REFERENCE.y) * 1e3, "mic_distance_mm": distance * 1e3}


def flat_drivers(drv: dict) -> dict:
    """The same drivers with a flat response: the DED transfer per unit on-axis response."""
    return {n: replace(d, response=np.ones_like(d.response)) for n, d in drv.items()}
