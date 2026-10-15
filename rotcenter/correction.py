"""Rotation-center correction of turntable measurements, one plane at a time.

A polar measurement should turn the loudspeaker about a reference point on its baffle, so that
the microphone stays on a circle of radius ``r0`` around that point. When the turntable axis is
elsewhere (behind the baffle, off to the side), the microphone travels on a circle around the
axis instead: off axis it is farther away, and seen from the driver it is at another angle.

The correction rebuilds, for every nominal angle, the response the microphone would have picked
up at the ideal position. Plane frame (m): rotation axis at the origin, ``y`` toward the front,
``x`` toward the + side (the side the microphone faces at +90 degrees).

* Reference point ``R = (lateral, depth)``, driver ``S = R + (driver_offset, 0)``, microphone
  distance ``r0`` (from the baffle at 0 degrees).
* Measurement at nominal angle ``theta``: microphone at the 0-degree position ``(R_x, R_y + r0)``
  rotated by ``theta`` about the axis, ``M_k = (x cos t + y sin t, -x sin t + y cos t)``.
* Ideal position: ``M* = R + r0 (sin theta, cos theta)``.
* Seen from the driver, every position has a direction ``phi = atan2(M_x - S_x, M_y - S_y)``
  (unwrapped next to the nominal angle) and a distance ``rho = |M - S|``.
* Point-source model around the driver, ``h(phi, rho, t) = g(phi, t - rho/c) / rho``: moving a
  measurement from ``rho_k`` to ``rho*`` is a delay ``(rho* - rho_k)/c`` (fractional, applied as a
  linear phase) and a gain ``rho_k / rho*``.
* The response in direction ``phi*`` is interpolated linearly in ``phi`` between the two
  measurements that bracket it, once both have been moved to ``rho*`` (so that they are in time
  and do not comb-filter). Outside the measured range the nearest measurement is used; a full
  circle wraps around. The 0-degree measurement is the reference and is never changed.

Assumptions: point source on the baffle plane, far-field directivity independent of distance,
everything in the plane of rotation (heights are ignored).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import C

#: Offsets below this (mm) count as zero.
NEGLIGIBLE_MM = 0.05


@dataclass(frozen=True)
class PlaneOffsets:
    depth_mm: float = 0.0
    lateral_mm: float = 0.0
    driver_offset_mm: float = 0.0
    mic_distance_mm: float = 1000.0

    def is_identity(self) -> bool:
        return abs(self.depth_mm) < NEGLIGIBLE_MM and abs(self.lateral_mm) < NEGLIGIBLE_MM


class PlaneGeometry:
    def __init__(self, o: PlaneOffsets):
        self.reference = (o.lateral_mm * 1e-3, o.depth_mm * 1e-3)
        self.source = (self.reference[0] + o.driver_offset_mm * 1e-3, self.reference[1])
        self.r0 = max(o.mic_distance_mm, 1.0) * 1e-3

    def seen_from_driver(self, mic, angle):
        dx, dy = mic[0] - self.source[0], mic[1] - self.source[1]
        phi = math.atan2(dx, dy)
        phi += 2.0 * math.pi * round((angle - phi) / (2.0 * math.pi))
        return phi, max(math.hypot(dx, dy), 1e-6)

    def measured(self, angle_deg):
        t = math.radians(angle_deg)
        x, y = self.reference[0], self.reference[1] + self.r0
        mic = (x * math.cos(t) + y * math.sin(t), -x * math.sin(t) + y * math.cos(t))
        return self.seen_from_driver(mic, t)

    def ideal(self, angle_deg):
        t = math.radians(angle_deg)
        mic = (self.reference[0] + self.r0 * math.sin(t), self.reference[1] + self.r0 * math.cos(t))
        return self.seen_from_driver(mic, t)


@dataclass
class Sample:
    angle_deg: float
    ir: np.ndarray
    #: A copy made by symmetry from the measurement at ``-angle_deg``.
    mirrored: bool = False


def _bracket(phis, phi):
    last = len(phis) - 1
    if phi <= phis[0]:
        return [(0, 1.0)]
    if phi >= phis[last]:
        return [(last, 1.0)]
    hi = min(int(np.searchsorted(phis, phi, side="right")), last)
    lo = hi - 1
    span = phis[hi] - phis[lo]
    if span < 1e-9:
        return [(lo, 1.0)]
    u = (phi - phis[lo]) / span
    return [(lo, 1.0 - u), (hi, u)]


def plan_correction(angles, offsets: PlaneOffsets, mirrored=None, present=None):
    """The correction as weights: ``{i: [(j, u, g, delay_s), ...]}``, meaning that the corrected
    response at ``angles[i]`` is ``sum u * g * h_j(t - delay_s)`` over measured responses ``h_j``:
    ``u`` the interpolation weights (they sum to 1), ``g = rho_j / rho*`` the gain and ``delay_s``
    the delay that move ``h_j`` to the intended distance.

    ``mirrored[j]``: ``h_j`` is a copy made by symmetry from the measurement at ``-angles[j]``;
    ``present[j]``: there is a measurement at ``angles[j]``. Angles that stand as they are (0
    degrees, no offset, no data) are left out.
    """
    n = len(angles)
    mirrored = [False] * n if mirrored is None else list(mirrored)
    present = [True] * n if present is None else list(present)
    if offsets.is_identity():
        return {}
    g = PlaneGeometry(offsets)

    supports = []
    for j, angle in enumerate(angles):
        if not present[j]:
            continue
        if mirrored[j]:
            phi, dist = g.measured(-angle)
            phi = -phi
        else:
            phi, dist = g.measured(angle)
        supports.append((phi, dist, j))
    if len(supports) < 2:
        return {}
    supports.sort(key=lambda s: s[0])

    gaps = [b[0] - a[0] for a, b in zip(supports, supports[1:])]
    span = supports[-1][0] - supports[0][0]
    if len(supports) >= 3 and span >= 2.0 * math.pi - 1.5 * max(gaps):
        wrapped = [(p - 2.0 * math.pi, d, j) for p, d, j in supports]
        wrapped += supports
        wrapped += [(p + 2.0 * math.pi, d, j) for p, d, j in supports]
        supports = sorted(wrapped, key=lambda s: s[0])
    phis = np.array([s[0] for s in supports])

    plans = {}
    for i, angle in enumerate(angles):
        if not present[i] or abs(angle) <= 1e-9:
            continue
        phi, dist = g.ideal(angle)
        plans[i] = [(supports[k][2], u, supports[k][1] / dist, (dist - supports[k][1]) / C)
                    for k, u in _bracket(phis, phi)]
    return plans


def correct_plane(samples: list[Sample], offsets: PlaneOffsets, fs: float):
    """Corrected IR for every sample (``None`` where it stands as it is: 0 degrees, no offset)."""
    out = [None] * len(samples)
    if offsets.is_identity() or fs <= 0:
        return out
    plans = plan_correction([s.angle_deg for s in samples], offsets,
                            [s.mirrored for s in samples], [len(s.ir) > 0 for s in samples])
    if not plans:
        return out

    max_len = max(len(s.ir) for s in samples)
    max_shift = max(abs(d) * fs for parts in plans.values() for _, _, _, d in parts)
    pre = int(math.ceil(max_shift)) + 1
    n = 1 << int(math.ceil(math.log2(max_len + 2 * pre + 64)))
    bins = np.arange(n, dtype=float)
    signed = np.where(bins <= n // 2, bins, bins - n)

    for index, parts in plans.items():
        acc = np.zeros(n, complex)
        out_len = len(samples[index].ir)
        for j, u, g, delay_s in parts:
            ir = samples[j].ir
            delay = delay_s * fs
            spectrum = np.fft.fft(ir, n)
            acc += spectrum * np.exp(-2j * np.pi * signed * (delay + pre) / n) * u * g
            out_len = max(out_len, int(math.ceil(len(ir) + max(delay, 0.0))))
        ir = np.fft.ifft(acc).real
        out[index] = ir[pre : pre + min(out_len, n - pre)]
    return out


def correct_spectra(spectra, angles, offsets: PlaneOffsets, freqs, mirrored=None, levels: bool = False):
    """The same correction on complex spectra (rows: ``angles``) at any frequencies ``freqs`` (Hz):
    each delay is the exact phase ``exp(-j 2 pi f delay)``. Rows that stand as they are are copied.

    ``levels``: interpolate the re-propagated powers instead of the complex responses,
    ``|H*|^2 = sum u_k |g_k H_k|^2`` (``g_k = rho_k / rho*``), keeping the phase of the complex
    interpolation. Where the two measurements that bracket a direction are not coherent (an
    interference pattern finer than the angular step), their complex mean cancels in part; the
    mean of their powers does not."""
    spectra = np.asarray(spectra, complex)
    freqs = np.asarray(freqs, float)
    out = spectra.copy()
    for i, parts in plan_correction(angles, offsets, mirrored).items():
        moved = [(u, g * np.exp(-2j * np.pi * freqs * delay) * spectra[j]) for j, u, g, delay in parts]
        out[i] = sum(u * h for u, h in moved)
        if levels and len(moved) > 1:
            out[i] = np.sqrt(sum(u * np.abs(h) ** 2 for u, h in moved)) * np.exp(1j * np.angle(out[i]))
    return out
