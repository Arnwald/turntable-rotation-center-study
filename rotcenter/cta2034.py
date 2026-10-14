"""CTA-2034-A spatial averages from the horizontal orbit only.

Only the horizontal orbit is simulated (a turntable measurement), so the curves below are the
horizontal-plane estimates of the standard's curves:

* Listening Window: RMS average of 0, +-10, +-20, +-30 degrees (the standard adds +-10 degrees
  vertical);
* Early Reflections: RMS average of the front-wall (0, +-10, +-20, +-30), side-wall (+-40 ... +-80)
  and rear-wall (+-90, 180) groups, each itself an RMS average (the standard adds the floor
  and ceiling bounces, taken on the vertical orbit);
* Sound Power: weighted RMS average of the 36 horizontal directions with the standard's weights
  (equivalent to assuming the vertical orbit identical to the horizontal one);
* Predicted In-Room Response: ``PIR^2 = 0.12 LW^2 + 0.44 ER^2 + 0.44 SP^2`` (squared pressures).

The sound power weights of the standard (Appendix C) are the solid angles of the bands between
the "spherical squares" of half-width ``theta - 5`` and ``theta + 5`` degrees centered on the
reference axis: ``w(theta) = S(theta + 5) - S(theta - 5)``, ``S(alpha) = 4 C - 2 pi`` with ``C`` the
corner angle of the square bounded by the great circles at +-alpha horizontally and vertically.
The weights are symmetric about 90 degrees and the 90-degree band is counted twice.

All inputs are power spectra ``|p|^2`` keyed by the horizontal angle in degrees.
"""

from __future__ import annotations

import math

import numpy as np

LISTENING_WINDOW = (0, 10, -10, 20, -20, 30, -30)
FRONT_WALL = (0, 10, -10, 20, -20, 30, -30)
SIDE_WALL = (40, -40, 50, -50, 60, -60, 70, -70, 80, -80)
REAR_WALL = (90, -90, 180)


def _square_solid_angle(alpha_deg: float) -> float:
    a = math.radians(alpha_deg)
    gamma = math.acos(math.cos(a) * math.cos(a))
    corner_a = math.atan(math.sin(a) / math.tan(a))
    corner = math.acos(-math.cos(corner_a) ** 2 + math.sin(corner_a) ** 2 * math.cos(gamma))
    return 4.0 * corner - 2.0 * math.pi


def sound_power_weights() -> dict[int, float]:
    """Weight of every angle of the 10-degree orbit, 0 ... 180 (symmetric about 90)."""
    edges = [i * 10 + 5 for i in range(9)] + [90]
    areas = [_square_solid_angle(e) for e in edges]
    w = [areas[0]] + [areas[i] - areas[i - 1] for i in range(1, 10)]
    w[9] *= 2.0
    out = {}
    for i, weight in enumerate(w):
        out[10 * i] = weight
        out[180 - 10 * i] = weight
    return out


SP_WEIGHTS = sound_power_weights()


def _fold(angle: float) -> int:
    a = int(round(angle)) % 360
    return a - 360 if a > 180 else a


def _rms(power: dict, angles) -> np.ndarray:
    return np.mean([power[_fold(a)] for a in angles], axis=0)


def horizontal_cta2034(power: dict) -> dict[str, np.ndarray]:
    """Horizontal-only LW, ER (and its groups), SP, PIR and DI as power spectra (DI in dB)."""
    power = {_fold(a): np.asarray(p, float) for a, p in power.items()}
    if -180 in power and 180 not in power:
        power[180] = power.pop(-180)
    lw = _rms(power, LISTENING_WINDOW)
    front = _rms(power, FRONT_WALL)
    side = _rms(power, SIDE_WALL)
    rear = _rms(power, REAR_WALL)
    er = (front + side + rear) / 3.0
    weights = {a: SP_WEIGHTS[abs(a)] for a in power}
    sp = sum(weights[a] * p for a, p in power.items()) / sum(weights.values())
    pir = 0.12 * lw + 0.44 * er + 0.44 * sp
    with np.errstate(divide="ignore", invalid="ignore"):
        di = 10.0 * np.log10(lw / sp)
    return {"on_axis": power[0], "lw": lw, "er_front": front, "er_side": side, "er_rear": rear,
            "er": er, "sp": sp, "pir": pir, "di": di}
