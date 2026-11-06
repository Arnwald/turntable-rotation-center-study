"""The correction in 3D, at another distance, with powers interpolated, and across the plane."""

import numpy as np
import pytest

from rotcenter.constants import C
from rotcenter.correction import (Across, PlaneOffsets, PlaneSetup, Sample, correct_plane, correct_spectra,
                                  plan_correction)

FS = 48_000.0
LEN = 4096
N = 16_384
ANGLES = [float(a) for a in range(-170, 181, 10)]
NINETY = [float(a) for a in range(-90, 91, 10)]


def rippled_ir_at(phi, level_dist, time_dist, echo_m=0.0):
    """Band-limited IR of a directional point source seen ``phi`` (radians) off its axis: amplitude
    ``(0.55 + 0.45 cos phi) / level_dist``, a first-order low-pass whose corner falls off axis,
    arriving ``time_dist / c`` after a 2 ms pre-delay; with ``echo_m``, two baffle-edge echoes
    whose delays change with the direction, so that their ripples get finer with frequency."""
    f = np.fft.rfftfreq(N, 1.0 / FS)
    taper = np.where(f < 14_000, 1.0, np.where(f < 18_000, 0.5 * (1 + np.cos(np.pi * (f - 14_000) / 4_000)), 0.0))
    corner = 2_000.0 + 6_000.0 * (1.0 + np.cos(phi))
    spectrum = taper / (1 + 1j * f / corner) * (0.55 + 0.45 * np.cos(phi)) / level_dist
    spectrum = spectrum * np.exp(-2j * np.pi * f * (0.002 + time_dist / C))
    if echo_m > 0.0:
        echoes = ((0.8, echo_m * (1 - np.sin(phi)) / C), (0.6, 0.6 * echo_m * (1 + np.sin(phi)) / C))
        spectrum = spectrum * (1 + sum(level * np.exp(-2j * np.pi * f * delay) for level, delay in echoes))
    return np.fft.irfft(spectrum, N)[:LEN]


def source_ir(phi, dist, echo_m=0.0):
    return rippled_ir_at(phi, dist, dist, echo_m)


def error_db(a, b):
    n = max(len(a), len(b))
    a, b = np.pad(a, (0, n - len(a))), np.pad(b, (0, n - len(b)))
    return 10 * np.log10(np.sum((a - b) ** 2) / np.sum(b**2))


def samples(angles, irs):
    return [Sample(a, ir) for a, ir in zip(angles, irs)]


def measure(setup, angles, echo_m=0.0):
    return [source_ir(*setup.measured(a), echo_m) for a in angles]


def calculated_at(radius, source):
    """Measured at 1 m about the reference point (the origin), results calculated at ``radius``,
    the driver at ``source``."""
    return PlaneSetup((1, 0, 0), (0, 0, 1), source, (0, 0, 0), (0, 0, 1.0), (0, 0, 0), radius, 1.0)


def expected_at(setup, angle, theta):
    """What ``setup``'s target at ``angle`` must be, for a source whose response depends on the
    angle ``theta(direction)``: seen from the target circle, brought back to the reference distance."""
    d = setup.target_position(angle) - setup.source
    dist = np.linalg.norm(d)
    level, time = setup.reference_scaling()
    return rippled_ir_at(theta(d), dist / level, dist - time)


def test_at_another_distance_the_reference_point_stays_as_measured():
    setup = calculated_at(2.0, (0, 0, 0))
    assert setup.is_identity()
    out = correct_plane(samples(NINETY, measure(setup, NINETY)), setup, FS)
    assert all(x is None for x in out)


def test_at_another_distance_a_driver_off_the_reference_point_is_seen_from_there():
    """A woofer 200 mm to the side of the reference point, measured at 1 m, calculated at 2 m: every
    angle, 0 degrees too, rebuilt as heard from the 2 m circle, with the level and time of flight of
    the reference point at 1 m."""
    setup = calculated_at(2.0, (0.2, 0, 0))
    assert not setup.is_identity()
    out = correct_plane(samples(NINETY, measure(setup, NINETY)), setup, FS)
    phi = lambda d: np.arctan2(d[0], d[2])  # noqa: E731
    for i, a in enumerate(NINETY):
        assert error_db(out[i], expected_at(setup, a, phi)) < -25.0, a


def test_offsets_with_a_calculation_distance_describe_the_same_setup():
    o = PlaneOffsets(depth_mm=145.0, lateral_mm=120.0, driver_offset_mm=-205.0, mic_distance_mm=1000.0,
                     calculation_distance_mm=2000.0)
    g = PlaneSetup.from_offsets(o)
    assert g.reference_scaling() == pytest.approx((2.0, 1.0))
    assert g.target_position(0.0) == pytest.approx([0.12, 0.0, 2.145])
    assert g.measured_position(0.0) == pytest.approx([0.12, 0.0, 1.145])


def test_a_driver_under_the_plane_is_seen_from_another_angle_at_another_distance():
    """A woofer 250 mm under the reference point, measured in the horizontal plane at the
    reference point's height (seen from 14 degrees above at 1 m) and calculated at 2 m (7 degrees
    above): its vertical measurements give the change of its response with that angle."""
    source = (0.0, -0.25, 0.0)
    horizontal = calculated_at(2.0, source)
    vertical = PlaneSetup((0, 1, 0), (0, 0, 1), source, (0, 0, 0), (0, 0, 1.0), (0, 0, 0), 2.0, 1.0)
    off_axis = lambda d: np.arccos(np.clip(d[2] / np.linalg.norm(d), -1, 1))  # noqa: E731

    def measure_3d(setup):
        out = []
        for a in NINETY:
            d = setup.measured_position(a) - np.asarray(source)
            out.append(rippled_ir_at(off_axis(d), np.linalg.norm(d), np.linalg.norm(d)))
        return out

    h, v = measure_3d(horizontal), measure_3d(vertical)
    across = Across(NINETY, v, vertical, (0, 1, 0))
    with_ = correct_plane(samples(NINETY, h), horizontal, FS, across=across)
    without = correct_plane(samples(NINETY, h), horizontal, FS)
    for i, a in enumerate(NINETY):
        expected = expected_at(horizontal, a, off_axis)
        e_with, e_without = error_db(with_[i], expected), error_db(without[i], expected)
        assert e_with < -40.0, (a, e_with)
        if abs(a) < 90:
            assert e_with < e_without - 5.0, (a, e_with, e_without)


def test_a_driver_out_of_the_plane_of_the_microphone():
    """A woofer 250 mm below the microphone's height, the axis 180 mm behind the baffle: the
    distances are 3D."""
    o = PlaneOffsets(depth_mm=180.0, out_of_plane_mm=-250.0)
    g = PlaneSetup.from_offsets(o)
    theta = lambda d: np.arctan2(d[0], d[2])  # noqa: E731
    out = correct_plane(samples(ANGLES, measure(g, ANGLES)), o, FS)
    for i, a in enumerate(ANGLES):
        if a == 0:
            assert out[i] is None
            continue
        assert error_db(out[i], expected_at(g, a, theta)) < -25.0, a
    # Taken as in the plane, the woofer is rebuilt with the wrong distances.
    flat = correct_plane(samples(ANGLES, measure(g, ANGLES)), PlaneOffsets(depth_mm=180.0), FS)
    i = ANGLES.index(180.0)
    assert error_db(flat[i], expected_at(g, 180.0, theta)) > error_db(out[i], expected_at(g, 180.0, theta)) + 10


def band_levels(ir, bands, n=8192):
    spectrum = np.fft.rfft(ir, n)
    return np.array([10 * np.log10(np.sum(np.abs(spectrum[lo:hi]) ** 2)) for lo, hi in bands])


def test_interpolating_powers_avoids_notches_where_the_step_is_too_coarse():
    """Two edge echoes whose delays change by 0.15-0.25 ms every 10 degrees near the axis: above
    2-3 kHz, 10-degree steps undersample their ripples."""
    echo = 0.5
    o = PlaneOffsets(depth_mm=150.0)
    g = PlaneSetup.from_offsets(o)
    s = samples(ANGLES, measure(g, ANGLES, echo))
    levels = correct_plane(s, o, FS)
    complex_ = correct_plane(s, o, FS, levels=False)
    n = 8192
    to_bin = lambda f: int(f / FS * n)  # noqa: E731
    bands = [(to_bin(f), to_bin(f * 2 ** (1 / 6))) for f in 2000.0 * 2 ** (np.arange(16) / 6)]

    def errors(out):
        e = []
        for i, a in enumerate(ANGLES):
            if a != 0:
                e.extend(band_levels(out[i], bands) - band_levels(source_ir(*g.ideal(a), echo), bands))
        return np.array(e)

    lev, cpx = errors(levels), errors(complex_)
    assert cpx.mean() < -0.7
    assert abs(lev.mean()) < 0.5 * abs(cpx.mean())
    assert np.sqrt(np.mean(lev**2)) < 0.92 * np.sqrt(np.mean(cpx**2))


def test_a_measurement_taken_off_its_nominal_angle():
    """The table stopped at 31.5 degrees for the 30-degree measurement: rebuilt at 30, even with
    the axis on the reference point."""
    g = PlaneSetup.from_offsets(PlaneOffsets())
    angles, taken = [0.0, 30.0, 60.0], [0.0, 31.5, 60.0]
    irs = [source_ir(*g.measured(t)) for t in taken]
    out = correct_plane([Sample(a, ir, measured_deg=t) for a, t, ir in zip(angles, taken, irs)], g, FS)
    assert out[0] is None and out[2] is None
    ideal = source_ir(*g.ideal(30.0))
    assert error_db(out[1], ideal) < error_db(irs[1], ideal) - 6.0


def test_half_circle_with_symmetry_copies():
    """Measured from 0 to 180 degrees and completed by symmetry: the copies stand for the mirror
    image of their source."""
    o = PlaneOffsets(depth_mm=150.0)
    g = PlaneSetup.from_offsets(o)
    real = [float(a) for a in range(0, 181, 10)]
    irs = measure(g, real)
    s = samples(real, irs) + [Sample(-a, ir, mirrored=True) for a, ir in zip(real[1:18], irs[1:18])]
    out = correct_plane(s, o, FS)
    for i, sample in enumerate(s[1:], start=1):
        assert error_db(out[i], source_ir(*g.ideal(sample.angle_deg))) < -25.0, sample.angle_deg
    plus, minus = out[3], out[len(real) + 2]
    assert s[3].angle_deg == -s[len(real) + 2].angle_deg
    assert error_db(minus, plus) < -100.0


def test_spectra_and_impulse_responses_give_the_same_correction_across_the_plane():
    """The frequency-domain version, out-of-plane factor included, against the time-domain one."""
    source = (0.0, -0.25, 0.0)
    horizontal = PlaneSetup((1, 0, 0), (0, 0, 1), source, (0, 0, -0.15), (0, 0, 1.0), (0, 0, 0), 2.0, 1.0)
    vertical = PlaneSetup((0, 1, 0), (0, 0, 1), source, (0, 0, -0.15), (0, 0, 1.0), (0, 0, 0), 2.0, 1.0)
    h = [source_ir(*horizontal.measured(a)) for a in NINETY]
    v = [source_ir(*vertical.measured(a)) for a in NINETY]
    out = correct_plane(samples(NINETY, h), horizontal, FS, across=Across(NINETY, v, vertical))
    n = 8192
    freqs = np.fft.rfftfreq(n, 1.0 / FS)
    spectra = lambda irs: np.stack([np.fft.rfft(ir, n) for ir in irs])  # noqa: E731
    fixed = correct_spectra(spectra(h), NINETY, horizontal, freqs, across=Across(NINETY, spectra(v), vertical))
    band = (freqs > 200) & (freqs < 12_000)
    for i in range(len(NINETY)):
        ref = np.fft.rfft(out[i], n)
        level = 20 * np.log10(np.abs(fixed[i][band]) / np.abs(ref[band]))
        # The smoothing of the factor runs on two grids (8192 and 16384 points).
        assert np.max(np.abs(level)) < 0.05, NINETY[i]


def test_plans_carry_the_level_and_time_of_the_reference_distance():
    g = PlaneSetup.from_offsets(PlaneOffsets(depth_mm=145.0, mic_distance_mm=1000.0, calculation_distance_mm=4000.0))
    plans = {p.index: p for p in plan_correction(ANGLES, g)}
    zero = ANGLES.index(0.0)
    assert zero not in plans
    for i, a in enumerate(ANGLES):
        if a == 0:
            continue
        phi, dist = g.ideal(a)
        for j, u, gain, delay in plans[i].parts:
            _, rho = g.measured(ANGLES[j])
            assert gain == pytest.approx(rho / dist * 4.0)
            assert delay == pytest.approx((dist - rho - 3.0) / C)
