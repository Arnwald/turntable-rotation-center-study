import numpy as np
import pytest

from rotcenter.constants import C
from rotcenter.correction import PlaneGeometry, PlaneOffsets, Sample, correct_plane, correct_spectra
from rotcenter.geometry import CENTERED, MID, REFERENCE, RIGOROUS, TWEETER, WOOFER, correction_offsets, mic_positions

ANGLES = np.arange(-170, 181, 10)
FS = 48_000.0


def test_rigorous_microphone_stays_on_its_circle():
    for d in (1.0, 2.0):
        m = mic_positions(RIGOROUS, ANGLES, d)
        r = np.hypot(m[:, 0] - REFERENCE.x, m[:, 2])
        assert np.allclose(r, d)
        assert np.allclose(m[:, 1], REFERENCE.y)


def test_centered_microphone_turns_about_the_cabinet_center():
    d = 1.0
    m = mic_positions(CENTERED, ANGLES, d)
    r = np.hypot(m[:, 0] - CENTERED.axis_x, m[:, 2] - CENTERED.axis_z)
    assert np.allclose(r, np.hypot(REFERENCE.x, d + 0.145))
    zero = m[list(ANGLES).index(0)]
    assert np.allclose(zero, [REFERENCE.x, REFERENCE.y, d])
    # +90 degrees: the microphone is on the tweeter side (+x).
    assert m[list(ANGLES).index(90), 0] > 1.0


def test_correction_parameters_of_the_centered_setup():
    o = correction_offsets(CENTERED, WOOFER, 1.0)
    assert o["depth_mm"] == pytest.approx(145.0)
    assert o["lateral_mm"] == pytest.approx(120.0)
    assert o["driver_offset_mm"] == pytest.approx(-205.0)
    assert correction_offsets(CENTERED, TWEETER, 2.0)["driver_offset_mm"] == pytest.approx(0.0)
    assert correction_offsets(CENTERED, MID, 2.0)["mic_distance_mm"] == pytest.approx(2000.0)


@pytest.mark.parametrize("driver", [WOOFER, MID, TWEETER])
@pytest.mark.parametrize("distance", [1.0, 2.0])
def test_plane_model_matches_the_3d_kinematics(driver, distance):
    """The correction's 2D model of the measured positions (microphone shifted sideways onto the
    reference axis) must agree with the independent 3D kinematics of the simulation."""
    g = PlaneGeometry(PlaneOffsets(**correction_offsets(CENTERED, driver, distance)))
    mics = mic_positions(CENTERED, ANGLES, distance)
    for a, m in zip(ANGLES, mics):
        phi, rho = g.measured(a)
        dx, dz = m[0] - driver.x, m[2]
        assert rho == pytest.approx(np.hypot(dx, dz), abs=1e-12)
        assert np.sin(phi) == pytest.approx(dx / np.hypot(dx, dz), abs=1e-12)
        assert np.cos(phi) == pytest.approx(dz / np.hypot(dx, dz), abs=1e-12)
    rigorous = mic_positions(RIGOROUS, ANGLES, distance)
    for a, m in zip(ANGLES, rigorous):
        phi, rho = g.ideal(a)
        assert rho == pytest.approx(np.hypot(m[0] - driver.x, m[2]), abs=1e-12)


def pulse(delay_s, gain, n=4096, corner=None, fs=FS):
    """Band-limited pulse (raised-cosine spectrum below 16 kHz) with an optional low-pass."""
    f = np.fft.rfftfreq(n, 1.0 / fs)
    taper = np.where(f < 12_000, 1.0, np.where(f < 16_000, 0.5 * (1 + np.cos(np.pi * (f - 12_000) / 4_000)), 0.0))
    h = taper * np.exp(-2j * np.pi * f * delay_s) * gain
    if corner is not None:
        h = h / (1 + 1j * f / corner)
    return np.fft.irfft(h, n)


def measured_set(offsets, directivity):
    g = PlaneGeometry(offsets)
    irs = []
    for a in ANGLES:
        phi, rho = g.measured(a)
        gain, corner = directivity(phi)
        irs.append(pulse(0.002 + rho / C, gain / rho, corner=corner))
    return irs


def ideal(offsets, angle, directivity):
    phi, rho = PlaneGeometry(offsets).ideal(angle)
    gain, corner = directivity(phi)
    return pulse(0.002 + rho / C, gain / rho, corner=corner)


def error_db(a, b):
    n = max(len(a), len(b))
    a = np.pad(a, (0, n - len(a)))
    b = np.pad(b, (0, n - len(b)))
    return 10 * np.log10(np.sum((a - b) ** 2) / np.sum(b**2))


def test_identity_and_reference_are_left_alone():
    o = PlaneOffsets(driver_offset_mm=-205.0)
    assert all(x is None for x in correct_plane([Sample(a, np.ones(8)) for a in ANGLES], o, FS))
    o = PlaneOffsets(depth_mm=145.0, lateral_mm=120.0)
    out = correct_plane([Sample(a, ir) for a, ir in zip(ANGLES, measured_set(o, lambda p: (1.0, None)))], o, FS)
    assert out[list(ANGLES).index(0)] is None
    assert all(x is not None for i, x in enumerate(out) if ANGLES[i] != 0)


@pytest.mark.parametrize("driver_offset", [0.0, -205.0])
def test_omnidirectional_point_source_is_corrected_exactly(driver_offset):
    o = PlaneOffsets(depth_mm=145.0, lateral_mm=120.0, driver_offset_mm=driver_offset, mic_distance_mm=1000.0)
    omni = lambda phi: (1.0, None)  # noqa: E731
    irs = measured_set(o, omni)
    out = correct_plane([Sample(a, ir) for a, ir in zip(ANGLES, irs)], o, FS)
    for i, a in enumerate(ANGLES):
        if a == 0:
            continue
        assert error_db(out[i][:4096], ideal(o, a, omni)) < -60.0


def test_directional_point_source_is_corrected_closely():
    o = PlaneOffsets(depth_mm=145.0, lateral_mm=120.0, driver_offset_mm=0.0, mic_distance_mm=1000.0)
    beam = lambda phi: (0.55 + 0.45 * np.cos(phi), 2000.0 + 6000.0 * (1 + np.cos(phi)))  # noqa: E731
    irs = measured_set(o, beam)
    out = correct_plane([Sample(a, ir) for a, ir in zip(ANGLES, irs)], o, FS)
    for i, a in enumerate(ANGLES):
        if a == 0:
            continue
        before = error_db(irs[i], ideal(o, a, beam))
        after = error_db(out[i][:4096], ideal(o, a, beam))
        assert after < -25.0 and after < before - 10.0, (a, before, after)


def test_frequency_domain_correction_matches_the_time_domain():
    o = PlaneOffsets(depth_mm=145.0, lateral_mm=120.0, driver_offset_mm=-205.0, mic_distance_mm=1000.0)
    beam = lambda phi: (0.55 + 0.45 * np.cos(phi), 2000.0 + 6000.0 * (1 + np.cos(phi)))  # noqa: E731
    irs = measured_set(o, beam)
    out = correct_plane([Sample(a, ir) for a, ir in zip(ANGLES, irs)], o, FS)
    n = 8192
    freqs = np.fft.rfftfreq(n, 1.0 / FS)
    spectra = np.stack([np.fft.rfft(ir, n) for ir in irs])
    fixed = correct_spectra(spectra, ANGLES, o, freqs)
    for i in range(len(ANGLES)):
        ref = spectra[i] if out[i] is None else np.fft.rfft(out[i], n)
        # Up to the band-limited tails the time-domain version cuts at the ends of its window.
        assert np.max(np.abs(fixed[i] - ref)) < 2e-5 * np.max(np.abs(ref))


@pytest.mark.parametrize("driver_offset", [0.0, -205.0])
def test_frequency_domain_correction_of_a_point_source_is_exact(driver_offset):
    o = PlaneOffsets(depth_mm=145.0, lateral_mm=120.0, driver_offset_mm=driver_offset, mic_distance_mm=2000.0)
    g = PlaneGeometry(o)
    freqs = 20.0 * 2.0 ** (np.arange(400) / 40.0)
    source = lambda phi, rho: np.exp(-2j * np.pi * freqs * rho / C) / rho  # noqa: E731
    spectra = np.stack([source(*g.measured(a)) for a in ANGLES])
    fixed = correct_spectra(spectra, ANGLES, o, freqs)
    for i, a in enumerate(ANGLES):
        ideal_spectrum = source(*g.ideal(a)) if a != 0 else spectra[i]
        assert np.max(np.abs(fixed[i] - ideal_spectrum)) < 1e-12 * np.max(np.abs(ideal_spectrum))
