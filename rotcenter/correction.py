"""Rotation-center correction of turntable measurements, one plane at a time.

A polar measurement should turn the loudspeaker about a reference point on its baffle, so that
the microphone stays on a circle of radius ``r0`` around that point. When the turntable axis is
elsewhere (behind the baffle, off to the side), the microphone travels on a circle around the
axis instead: off axis it is farther away, and seen from the driver it is at another angle.

The correction rebuilds, for every nominal angle, the response the microphone would have picked
up at the ideal position.

* **Geometry, per driver, in 3D.** A plane of rotation is spanned by two unit vectors: ``front``,
  the direction of the microphone at 0 degrees, and ``side``, where positive angles take it. The
  turntable axis is perpendicular to the plane. A measurement at angle ``theta`` was taken at the
  0-degree position rotated by ``theta`` about that axis; its ideal position is on the circle of
  radius ``r*`` around the reference point. Seen from the driver's acoustic center ``S``, every
  position has a direction ``phi`` in the plane (``atan2`` of its ``side`` and ``front``
  components, unwrapped next to the nominal angle) and a distance ``rho = |M - S|`` in 3D: a driver
  above or below the microphone's height is farther away than its projection on the plane.
* **Point-source model** around the driver, ``h(phi, rho, t) = g(phi, t - rho/c) / rho``: moving a
  measurement from ``rho_k`` to ``rho*`` is a delay ``(rho* - rho_k)/c`` (fractional, applied as a
  linear phase) and a gain ``rho_k / rho*``.
* **Angular interpolation.** The response in direction ``phi*`` is interpolated linearly in
  ``phi`` between the two measurements that bracket it, once both have been moved to ``rho*`` (so
  that they are in time and do not comb-filter). The powers are interpolated,
  ``|H*|^2 = sum u_k |g_k H_k|^2``, with the phase of the complex interpolation
  ``sum u_k g_k H_k``: where the angular step is too coarse for the pattern (interference ripples
  finer than the step, behind the cabinet at high frequency), the two measurements are no longer
  in phase and their complex mean cancels in part, while the mean of their powers does not. Where
  they are in phase, both give the same response. Outside the measured range the nearest
  measurement is used; a full circle wraps around. A measurement taken where it should have been
  (0 degrees, usually) stands as it is.
* **At another distance.** The circle of the targets can be larger than the one measured on (the
  calculation distance; CTA-2034-A measures at 2 m): ``rho*`` and ``phi*`` are then taken from that
  circle, and the targets are brought back to the level and the time of flight of the measurement
  distance (``r*/r0`` on the level, ``r* - r0`` on the time of flight), so that the reference point
  keeps its measurements as they are and only the differences between drivers and directions
  change. A driver out of the plane (the woofer under the tweeter, in the horizontal plane) is then
  seen from above or below at another angle than it was measured at: the other plane's
  measurements give the change of its response between the two angles (:class:`Across`), a
  magnitude factor smoothed over a sixth of an octave and limited to +-6 dB. The driver's power in
  any direction is taken as its two planes' at the same angle off its axis, weighted by how far the
  direction turns toward each plane.

Assumptions: point source on the baffle plane, far-field directivity independent of distance (no
change of the diffraction with the microphone distance).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import C

#: Offsets below this (mm) count as zero.
NEGLIGIBLE_MM = 0.05
#: Half the width of the smoothing of the out-of-plane factor (octaves): it follows how the driver
#: beams, not the fine ripples of two measurements 10 degrees apart.
ACROSS_SMOOTHING_OCT = 1.0 / 12.0
#: The largest out-of-plane factor, up or down (6 dB).
ACROSS_MAX_FACTOR = 2.0


@dataclass(frozen=True)
class PlaneOffsets:
    """The correction's inputs in one plane (mm). Plane frame: rotation axis at the origin, ``z``
    toward the front, ``x`` toward the + side (where the microphone is at +90 degrees), ``y`` out of
    the plane.

    ``depth_mm``: rotation axis to front baffle, positive when the axis is behind the baffle.
    ``lateral_mm``: rotation axis to reference axis, positive toward the + side.
    ``driver_offset_mm``: reference axis to the driver's acoustic center, positive toward the + side.
    ``mic_distance_mm``: front baffle to microphone at 0 degrees.
    ``out_of_plane_mm``: the driver's acoustic center out of the microphone's plane, toward ``y``
    (above the microphone's height in the horizontal plane).
    ``calculation_distance_mm``: radius of the circle the results are rebuilt on (``None``: the
    microphone distance).
    """

    depth_mm: float = 0.0
    lateral_mm: float = 0.0
    driver_offset_mm: float = 0.0
    mic_distance_mm: float = 1000.0
    out_of_plane_mm: float = 0.0
    calculation_distance_mm: float | None = None


def _unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


class PlaneSetup:
    """One plane of rotation in 3D (m, any origin): where the driver is, where the microphone was at
    each angle (``axis``: a point of the turntable axis, perpendicular to the plane; ``mic0``: the
    microphone at 0 degrees) and where it should have been (the circle of ``target_radius`` around
    ``target_center``). ``reference_radius``: the distance the levels and times of flight are kept
    at, the measurement distance."""

    def __init__(self, side, front, source, axis, mic0, target_center, target_radius, reference_radius=None):
        self.side = _unit(side)
        self.front = _unit(front)
        self.source = np.asarray(source, float)
        self.axis = np.asarray(axis, float)
        self.mic0 = np.asarray(mic0, float)
        self.target_center = np.asarray(target_center, float)
        self.target_radius = float(target_radius)
        self.reference_radius = float(target_radius if reference_radius is None else reference_radius)

    @classmethod
    def from_offsets(cls, o: PlaneOffsets) -> PlaneSetup:
        """The plane described by the correction's inputs (see :class:`PlaneOffsets`)."""
        reference = np.array([o.lateral_mm, 0.0, o.depth_mm]) * 1e-3
        r0 = max(o.mic_distance_mm, 1.0) * 1e-3
        radius = r0 if o.calculation_distance_mm is None else max(o.calculation_distance_mm, 1.0) * 1e-3
        source = reference + np.array([o.driver_offset_mm, o.out_of_plane_mm, 0.0]) * 1e-3
        return cls((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), source, (0.0, 0.0, 0.0), reference + [0.0, 0.0, r0],
                   reference, radius, r0)

    def measured_position(self, angle_deg):
        """Where the microphone was at ``angle_deg``: its 0-degree position rotated about the axis."""
        t = math.radians(angle_deg)
        d = self.mic0 - self.axis
        a, b = d @ self.side, d @ self.front
        out_of_plane = d - a * self.side - b * self.front
        s, c = math.sin(t), math.cos(t)
        return self.axis + (a * c + b * s) * self.side + (-a * s + b * c) * self.front + out_of_plane

    def target_position(self, angle_deg):
        """Where the microphone should have been at ``angle_deg``."""
        t = math.radians(angle_deg)
        return self.target_center + self.target_radius * (math.sin(t) * self.side + math.cos(t) * self.front)

    def seen_from_driver(self, point, angle):
        """Direction (radians, in the plane, unwrapped next to ``angle``) and 3D distance of ``point``
        seen from the driver."""
        d = point - self.source
        phi = math.atan2(d @ self.side, d @ self.front)
        phi += 2.0 * math.pi * round((angle - phi) / (2.0 * math.pi))
        return phi, max(float(np.linalg.norm(d)), 1e-6)

    def measured(self, angle_deg):
        return self.seen_from_driver(self.measured_position(angle_deg), math.radians(angle_deg))

    def ideal(self, angle_deg):
        return self.seen_from_driver(self.target_position(angle_deg), math.radians(angle_deg))

    def elevation(self, point, direction):
        """How far out of the plane ``point`` is seen from the driver, toward ``direction`` (radians)."""
        d = point - self.source
        return math.atan2(d @ direction, math.hypot(d @ self.side, d @ self.front))

    def direction_at(self, phi, elevation, direction):
        """The direction (from the driver) at ``phi`` in the plane and ``elevation`` toward ``direction``."""
        return (math.cos(elevation) * (math.sin(phi) * self.side + math.cos(phi) * self.front)
                + math.sin(elevation) * np.asarray(direction, float))

    def off_axis(self, d, direction):
        """``d`` as its angle off the driver's axis (``front``), signed toward the plane's side and
        toward ``direction``, and how far it turns toward ``direction`` (squared sine of its angle
        around the axis, from the plane)."""
        length = max(float(np.linalg.norm(d)), 1e-12)
        front, side, across = d @ self.front / length, d @ self.side / length, d @ direction / length
        theta = math.acos(min(max(front, -1.0), 1.0))
        around = math.hypot(side, across)
        share = (across / around) ** 2 if around > 1e-12 else 0.0
        return (-theta if side < 0 else theta), (-theta if across < 0 else theta), share

    def reference_scaling(self):
        """Level factor and time of flight (m) that bring the targets back to the reference distance."""
        return self.target_radius / self.reference_radius, self.target_radius - self.reference_radius

    def is_ideal_at(self, angle_deg):
        """The measurement at ``angle_deg`` was taken where it should have been: in the target's
        direction, at the target's distance brought back to the reference one."""
        measured = self.measured_position(angle_deg) - self.source
        target = self.target_position(angle_deg) - self.source
        level, time = self.reference_scaling()
        tol = NEGLIGIBLE_MM * 1e-3
        return (np.linalg.norm(measured - target / level) < tol
                and abs(np.linalg.norm(target) - np.linalg.norm(measured) - time) < tol)

    def is_identity(self):
        return all(self.is_ideal_at(a) for a in (45.0, 90.0, 135.0, 180.0, -90.0, 0.0))


def setup_of(plane) -> PlaneSetup:
    return plane if isinstance(plane, PlaneSetup) else PlaneSetup.from_offsets(plane)


@dataclass
class Sample:
    angle_deg: float
    ir: np.ndarray
    #: A copy made by symmetry from the measurement at ``-angle_deg``.
    mirrored: bool = False
    #: The angle the turntable really stopped at, when known (where it sits among the supports).
    measured_deg: float | None = None


@dataclass
class Across:
    """The other plane's measurements of the same driver, for what the correction sees out of its
    plane: ``spectra`` (rows: ``angles``; or IRs, for :func:`correct_plane`), taken where ``setup``
    says, and ``direction``, the other plane's + side in this plane's frame. With fewer than two of
    them, the driver is taken as symmetric about its axis (this plane's measurements)."""

    angles: list
    spectra: object
    setup: object
    direction: tuple = (0.0, 1.0, 0.0)
    mirrored: list | None = None


@dataclass
class Plan:
    """The corrected response at ``angles[index]``: ``sum u * g * h_j(t - delay_s)`` over ``parts``
    ``(j, u, g, delay_s)``, and, when the target is seen at another angle across the plane than the
    measurements it is interpolated from, the off-axis description of both directions."""

    index: int
    parts: list
    across: tuple | None = None


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


def _wrap_around(rows):
    """``rows`` (their direction first) sorted by direction, with copies 2 pi below and above when
    they cover a full circle (one missing step allowed)."""
    rows = sorted(rows, key=lambda r: r[0])
    if len(rows) < 3:
        return rows
    gaps = [b[0] - a[0] for a, b in zip(rows, rows[1:])]
    if rows[-1][0] - rows[0][0] < 2.0 * math.pi - 1.5 * max(gaps):
        return rows
    def shifted(turns):
        return [(r[0] + turns * 2.0 * math.pi, *r[1:]) for r in rows]
    return shifted(-1) + rows + shifted(1)


def _supports(angles, setup, mirrored, present, taken, direction):
    """(phi, rho, elevation, index) of every measurement, sorted by direction, wrapped around for a
    full circle (one missing step allowed)."""
    supports = []
    for j, angle in enumerate(angles):
        if not present[j]:
            continue
        t = angle if taken[j] is None else taken[j]
        position, mirror = (-t, -1.0) if mirrored[j] else (t, 1.0)
        phi, dist = setup.measured(position)
        supports.append((phi * mirror, dist, setup.elevation(setup.measured_position(position), direction), j))
    return _wrap_around(supports)


def plan_correction(angles, measured, targets=None, mirrored=None, present=None, taken=None,
                    direction=(0.0, 0.0, 0.0)) -> list[Plan]:
    """The correction of one plane as weights (see :class:`Plan`).

    ``measured``: where the measurements were taken (a :class:`PlaneSetup` or :class:`PlaneOffsets`);
    ``targets``: where they should have been (``None``: the same setup). ``mirrored[j]``: ``angles[j]``
    is a copy made by symmetry from the measurement at ``-angles[j]``; ``present[j]``: there is a
    measurement at ``angles[j]``; ``taken[j]``: the angle it was really taken at, when known;
    ``direction``: across the plane, toward the other plane's + side (zero: not looked at). Angles
    that stand as they are are left out."""
    measured = setup_of(measured)
    same = targets is None
    targets = measured if same else setup_of(targets)
    n = len(angles)
    mirrored = [False] * n if mirrored is None else list(mirrored)
    present = [True] * n if present is None else list(present)
    taken = [None] * n if taken is None else list(taken)
    exact = [t is None or abs(t - a) < 1e-9 for t, a in zip(taken, angles)]
    if all(exact) and same and measured.is_identity():
        return []
    direction = np.asarray(direction, float)
    supports = _supports(angles, measured, mirrored, present, taken, direction)
    if len(supports) < 2:
        return []
    phis = np.array([s[0] for s in supports])
    level_scale, time_offset = targets.reference_scaling()

    plans = []
    for i, angle in enumerate(angles):
        if not present[i] or (exact[i] and same and measured.is_ideal_at(angle)):
            continue
        phi, dist = targets.ideal(angle)
        bracket = _bracket(phis, phi)
        parts = [(supports[k][3], u, supports[k][1] / dist * level_scale,
                  (dist - supports[k][1] - time_offset) / C) for k, u in bracket]
        elevation = sum(u * supports[k][2] for k, u in bracket)
        target = targets.target_position(angle)
        across = None
        if abs(targets.elevation(target, direction) - elevation) > 1e-4:
            interpolated = targets.direction_at(phi, elevation, direction)
            across = (targets.off_axis(interpolated, direction), targets.off_axis(target - targets.source, direction))
        plans.append(Plan(i, parts, across))
    return plans


def smooth_octaves(freqs, power, half_oct):
    """``power`` (last axis: ``freqs``, increasing) averaged over ``+-half_oct`` octaves around each
    frequency."""
    freqs = np.asarray(freqs, float)
    lo = np.searchsorted(freqs, freqs * 2.0 ** -half_oct * (1 - 1e-12), side="left")
    hi = np.searchsorted(freqs, freqs * 2.0 ** half_oct * (1 + 1e-12), side="right")
    hi = np.maximum(hi, np.arange(len(freqs)) + 1)
    csum = np.concatenate([np.zeros(power.shape[:-1] + (1,)), np.cumsum(power, axis=-1)], axis=-1)
    return (csum[..., hi] - csum[..., lo]) / (hi - lo)


class Pattern:
    """A plane's measurements of the driver as smoothed powers brought to 1 m, by direction, wrapped
    around for a full circle as the supports are."""

    def __init__(self, power, angles, setup, mirrored=None):
        mirrored = [False] * len(angles) if mirrored is None else mirrored
        rows = []
        for p, angle, mirror in zip(power, angles, mirrored):
            phi, dist = setup.measured(-angle) if mirror else setup.measured(angle)
            rows.append((-phi if mirror else phi, p * dist**2))
        rows = _wrap_around(rows)
        self.phis = np.array([r[0] for r in rows])
        self.power = np.stack([r[1] for r in rows])

    def at(self, phi):
        """The power between the two directions around ``phi`` (the nearest one outside them)."""
        return sum(u * self.power[k] for k, u in _bracket(self.phis, phi))


def directivity_factor(own: Pattern, across: Pattern, start, end):
    """Magnitude factor per frequency from direction ``start`` to direction ``end``, each
    ``(angle in this plane, angle in the other, share of the other)``."""
    def power(d):
        own_angle, across_angle, share = d
        return (1.0 - share) * own.at(own_angle) + share * across.at(across_angle)
    a, b = power(start), power(end)
    with np.errstate(divide="ignore", invalid="ignore"):
        factor = np.sqrt(b / a)
    factor = np.where((a > 0) & (b > 0), factor, 1.0)
    return np.clip(factor, 1.0 / ACROSS_MAX_FACTOR, ACROSS_MAX_FACTOR)


def _patterns(own_power, angles, measured, mirrored, present, across: Across | None, across_power):
    """The driver's patterns in this plane and in the other one (this plane's again when the other
    has fewer than two measurements)."""
    keep = [j for j in range(len(angles)) if present[j]]
    own = Pattern(own_power[keep], [angles[j] for j in keep], measured, [mirrored[j] for j in keep])
    if across is not None and across_power is not None and len(across.angles) >= 2:
        other = Pattern(across_power, list(across.angles), setup_of(across.setup), across.mirrored)
    else:
        other = own
    return own, other


def correct_plane(samples: list[Sample], measured, fs: float, targets=None, across: Across | None = None,
                  levels: bool = True):
    """Corrected IR for every sample (``None`` where it stands as it is: 0 degrees, no offset).

    ``measured``, ``targets``: see :func:`plan_correction`; ``across``: the other plane's IRs (see
    :class:`Across`); ``levels``: interpolate the powers (``False``: the complex responses)."""
    out = [None] * len(samples)
    if fs <= 0:
        return out
    measured = setup_of(measured)
    direction = np.zeros(3) if across is None else np.asarray(across.direction, float)
    plans = plan_correction([s.angle_deg for s in samples], measured, targets, [s.mirrored for s in samples],
                            [len(s.ir) > 0 for s in samples], [s.measured_deg for s in samples], direction)
    if not plans:
        return out

    across_len = 0 if across is None else max((len(ir) for ir in across.spectra), default=0)
    max_len = max(max(len(s.ir) for s in samples), across_len)
    max_shift = max(abs(d) * fs for p in plans for _, _, _, d in p.parts)
    pre = int(math.ceil(max_shift)) + 1
    n = 1 << int(math.ceil(math.log2(max_len + 2 * pre + 64)))
    bins = np.arange(n, dtype=float)
    signed = np.where(bins <= n // 2, bins, bins - n)
    folded = np.minimum(np.arange(n), n - np.arange(n))

    factors = None
    if across is not None and any(p.across is not None for p in plans):
        freqs = np.arange(n // 2 + 1) * fs / n

        def powers(rows):
            return smooth_octaves(freqs, np.abs(np.fft.rfft(rows, n)) ** 2, ACROSS_SMOOTHING_OCT)
        present = [len(s.ir) > 0 for s in samples]
        own_rows = np.stack([np.pad(s.ir, (0, max_len - len(s.ir))) if len(s.ir) else np.zeros(max_len)
                             for s in samples])
        other = (np.stack([np.pad(ir, (0, max_len - len(ir))) for ir in across.spectra])
                 if len(across.angles) >= 2 else None)
        own, other = _patterns(powers(own_rows), [s.angle_deg for s in samples], measured,
                               [s.mirrored for s in samples], present, across,
                               None if other is None else powers(other))
        factors = {p.index: directivity_factor(own, other, *p.across)[folded] for p in plans if p.across}

    for plan in plans:
        acc = np.zeros(n, complex)
        power = np.zeros(n)
        out_len = len(samples[plan.index].ir)
        for j, u, g, delay_s in plan.parts:
            ir = samples[j].ir
            delay = delay_s * fs
            moved = np.fft.fft(ir, n) * np.exp(-2j * np.pi * signed * (delay + pre) / n) * g
            acc += u * moved
            power += u * np.abs(moved) ** 2
            out_len = max(out_len, int(math.ceil(len(ir) + max(delay, 0.0))))
        if levels and len(plan.parts) > 1:
            acc = _with_power(acc, power)
        if factors is not None and plan.index in factors:
            acc *= factors[plan.index]
        ir = np.fft.ifft(acc).real
        out[plan.index] = ir[pre : pre + min(out_len, n - pre)]
    return out


def _with_power(acc, power):
    """The complex response ``acc`` given the magnitude ``sqrt(power)`` (its phase kept)."""
    norm = np.abs(acc)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(norm > 1e-300, acc * (np.sqrt(power) / norm), np.sqrt(power) + 0j)


def correct_spectra(spectra, angles, measured, freqs, mirrored=None, targets=None, across: Across | None = None,
                    levels: bool = True, taken=None):
    """The same correction on complex spectra (rows: ``angles``) at any frequencies ``freqs`` (Hz):
    each delay is the exact phase ``exp(-j 2 pi f delay)``. Rows that stand as they are are copied.

    ``across``: the other plane's spectra at the same frequencies (see :class:`Across`); ``levels``:
    interpolate the powers (``False``: the complex responses); ``taken``: the angles the measurements
    were really taken at, when known."""
    spectra = np.asarray(spectra, complex)
    freqs = np.asarray(freqs, float)
    measured = setup_of(measured)
    out = spectra.copy()
    n = len(angles)
    mirrored = [False] * n if mirrored is None else list(mirrored)
    direction = np.zeros(3) if across is None else np.asarray(across.direction, float)
    plans = plan_correction(angles, measured, targets, mirrored, None, taken, direction)
    factors = {}
    if across is not None and any(p.across is not None for p in plans):
        def powers(rows):
            return smooth_octaves(freqs, np.abs(np.asarray(rows, complex)) ** 2, ACROSS_SMOOTHING_OCT)
        own, other = _patterns(powers(spectra), list(angles), measured, mirrored, [True] * n, across,
                               powers(across.spectra) if len(across.angles) >= 2 else None)
        factors = {p.index: directivity_factor(own, other, *p.across) for p in plans if p.across}
    for plan in plans:
        moved = [(u, g * np.exp(-2j * np.pi * freqs * delay) * spectra[j]) for j, u, g, delay in plan.parts]
        acc = sum(u * h for u, h in moved)
        if levels and len(moved) > 1:
            acc = _with_power(acc, sum(u * np.abs(h) ** 2 for u, h in moved))
        out[plan.index] = acc * factors.get(plan.index, 1.0)
    return out
