import numpy as np
import pytest

from rotcenter.constants import C
from rotcenter.cta2034 import SP_WEIGHTS, horizontal_cta2034
from rotcenter.diffraction import EdgeDiffraction, direct_term, edge_propagation
from rotcenter.directivity import Radiator, piston
from rotcenter.drivers import Driver, closed_box_alignment
from rotcenter.geometry import TWEETER, WOOFER
from rotcenter.signal import minimum_phase, rfft_freqs


def test_piston_on_axis_and_first_null():
    assert piston(np.array([0.0]))[0] == 1.0
    assert abs(piston(np.array([3.8317059702075125]))[0]) < 1e-12


def test_flat_dome_is_the_piston_and_domes_are_wider():
    k = 2 * np.pi * np.array([3000.0, 8000.0, 15000.0]) / C
    for theta in np.radians([20.0, 50.0, 80.0]):
        flat = Radiator(0.02, 0.0).directivity(k, theta)
        nearly_flat = Radiator(0.02, 1e-7).directivity(k, theta)
        assert np.allclose(flat, nearly_flat, atol=1e-5)
    dome, disc = Radiator(0.014, 0.0075), Radiator(0.014, 0.0)
    k = 2 * np.pi * 15000.0 / C
    assert abs(dome.directivity(np.array([k]), np.radians(60))[0]) > abs(disc.directivity(np.array([k]), np.radians(60))[0])
    assert np.allclose(dome.directivity(np.array([1.0, k]), 0.0), 1.0)


def test_closed_box_alignment():
    fc, qtc = closed_box_alignment()
    assert fc == pytest.approx(59.83, abs=0.01)
    assert qtc == pytest.approx(0.643, abs=0.001)


def test_minimum_phase_of_a_minimum_phase_filter():
    # Discrete minimum-phase system: zeros and poles inside the unit circle.
    n = 4096
    z = np.exp(1j * 2 * np.pi * np.arange(n // 2 + 1) / n)
    h = (1 - 0.5 / z) * (1 + 0.3 / z) / ((1 - 0.9 / z) * (1 - 0.6 / z + 0.25 / z**2))
    rebuilt = minimum_phase(np.abs(h))
    assert np.max(np.abs(rebuilt - h) / np.abs(h)) < 1e-6


@pytest.fixture(scope="module")
def small_drivers():
    freqs = rfft_freqs(2048)
    flat = np.ones_like(freqs, complex)
    woofer = Driver("woofer", "", WOOFER, Radiator(0.0838), flat)
    tweeter = Driver("tweeter", "", TWEETER, Radiator(0.014, 0.0075), flat)
    return freqs, {d.name: (d, EdgeDiffraction(d, freqs)) for d in (woofer, tweeter)}


def test_ded_baffle_step_on_axis(small_drivers):
    freqs, drivers = small_drivers
    k = 2 * np.pi * freqs / C
    for driver, edges in drivers.values():
        mic = np.array([driver.placement.x, driver.placement.y, 30.0])
        h = (direct_term(driver, mic, k) + edges.edge_term(edge_propagation(edges.points, mic, k))) * 30.0
        low = 20 * np.log10(abs(h[np.argmin(abs(freqs - 20))]))
        high = 20 * np.log10(abs(h[np.argmin(abs(freqs - 15000))]))
        assert low == pytest.approx(-6.0, abs=0.1)
        assert high == pytest.approx(0.0, abs=0.5)


def test_ded_edge_dipoles_vanish_in_the_baffle_plane(small_drivers):
    freqs, drivers = small_drivers
    k = 2 * np.pi * freqs / C
    _, edges = drivers["tweeter"]
    assert np.max(np.abs(edges.edge_term(edge_propagation(edges.points, np.array([3.0, 0.0, 0.0]), k)))) < 1e-9


def test_sound_power_weights():
    assert SP_WEIGHTS[0] == SP_WEIGHTS[180]
    assert SP_WEIGHTS[30] == SP_WEIGHTS[150]
    assert all(SP_WEIGHTS[a] < SP_WEIGHTS[a + 10] for a in range(0, 90, 10))


def test_cta_curves_of_an_omnidirectional_source():
    p = np.array([2.0, 3.0])
    curves = horizontal_cta2034({a: p for a in range(-170, 181, 10)})
    for name in ("lw", "er", "sp", "pir"):
        assert np.allclose(curves[name], p)
    assert np.allclose(curves["di"], 0.0)
