import json
from pathlib import Path

import numpy as np
import pytest

from rotcenter import twoway
from rotcenter.constants import C, P_REF, R_REF, RHO
from rotcenter.correction import PlaneGeometry, PlaneOffsets

ROOT = Path(__file__).resolve().parent.parent
FREQS = 20.0 * 800.0 ** (np.arange(464) / 463)


def test_microphones_are_the_points_of_the_bem_data():
    layout = json.loads((ROOT / "data" / "bem" / "points.json").read_text())
    points = np.array(layout["points"])
    for s in layout["sets"]:
        setup = s["name"].rsplit("-", 1)[0]
        mics = twoway.microphones(setup, s["distance"], s["angles_deg"])
        assert np.allclose(mics, points[s["first"] : s["first"] + s["count"]], atol=1e-12)


@pytest.mark.parametrize("distance", [1.0, 2.0])
def test_correction_parameters_match_the_kinematics(distance):
    angles = np.arange(-170, 181, 10)
    for plane, setup, rigorous in (("h", "centered-h", "rigorous-h"), ("v", "centered-v", "rigorous-v")):
        g = PlaneGeometry(PlaneOffsets(**twoway.correction_offsets(twoway.TWEETER, distance, plane)))
        measured = twoway.microphones(setup, distance, angles)
        ideal = twoway.microphones(rigorous, distance, angles)
        axis = 0 if plane == "h" else 1
        for a, m, i in zip(angles, measured, ideal):
            for (phi, rho), point in ((g.measured(a), m), (g.ideal(a), i)):
                d = point - np.array([twoway.TWEETER.x, twoway.TWEETER.y, 0.0])
                assert rho == pytest.approx(np.hypot(d[axis], d[2]), abs=1e-12)
                assert np.sin(phi) == pytest.approx(d[axis] / np.hypot(d[axis], d[2]), abs=1e-12)


def test_volume_velocity_gives_the_datasheet_response_in_an_infinite_baffle():
    drivers = twoway.drivers(FREQS)
    k = 2 * np.pi * FREQS / C
    for name, driver in drivers.items():
        u = twoway.volume_velocity(driver, FREQS)
        far_field = 1j * 2 * np.pi * FREQS * RHO * u * driver.radiator.on_axis(k) / (2 * np.pi * R_REF)
        assert np.allclose(far_field, driver.response)
    # The datasheets' passband levels: about 91 dB (woofer, 2.83 V), 90 dB (tweeter, 1 W).
    spl = {n: 20 * np.log10(np.abs(d.response) / P_REF) for n, d in drivers.items()}
    assert 88 < np.median(spl["woofer"][(FREQS > 200) & (FREQS < 1000)]) < 95
    assert 86 < np.median(spl["tweeter"][(FREQS > 3000) & (FREQS < 10000)]) < 95


def test_filters_equalize_a_flat_driver_to_its_target():
    flat = {"woofer": np.ones_like(FREQS, complex), "tweeter": np.ones_like(FREQS, complex)}
    filters = twoway.design_filters(FREQS, flat)
    level = P_REF * 10 ** (twoway.TARGET_DB / 20)
    for name, target in twoway.targets(FREQS).items():
        relevant = np.abs(target) > 1e-3
        assert np.allclose(filters[name][relevant], level * target[relevant], rtol=2e-3)
    total = twoway.system_sum(filters, flat)
    band = (FREQS > 200) & (FREQS < 15000)
    assert np.allclose(20 * np.log10(np.abs(total[band]) / level), 0.0, atol=0.02)


def test_filters_invert_smooth_drivers_and_sum_flat():
    s = 1j * FREQS
    woofer = (s / 60) ** 2 / ((s / 60) ** 2 + s / 60 / 0.7 + 1) / (1 + s / 3000) ** 2 * 0.02
    tweeter = (s / 900) ** 2 / ((s / 900) ** 2 + s / 900 / 0.8 + 1) / (1 + s / 20000) * 0.04
    drivers = {"woofer": woofer, "tweeter": tweeter}
    filters = twoway.design_filters(FREQS, drivers)
    level = P_REF * 10 ** (twoway.TARGET_DB / 20)
    for name, target in twoway.targets(FREQS).items():
        ratio = filters[name] * drivers[name] / (level * target)
        band = (FREQS > 150) & (FREQS < 15000)
        assert np.max(np.abs(20 * np.log10(np.abs(ratio[band])))) < 0.05
        # The phase is that of the target where the drivers overlap (the equalizer's own minimum
        # phase drifts by a few degrees, from the unknown response beyond the grid).
        overlap = (FREQS > 1000) & (FREQS < 3500)
        assert np.max(np.abs(np.angle(ratio[overlap], deg=True))) < 10.0
    total = twoway.system_sum(filters, drivers)
    band = (FREQS > 150) & (FREQS < 15000)
    assert np.max(np.abs(20 * np.log10(np.abs(total[band]) / level))) < 0.05


def test_orbit_takes_every_step_and_wraps_180():
    angles = np.arange(-180, 180, 5.0)
    spectra = {"x": angles[:, None] * np.ones((1, 3))}
    wanted, picked = twoway.orbit(spectra, angles, 10)
    assert list(wanted) == list(range(-170, 181, 10))
    assert np.allclose(picked["x"][:-1, 0], wanted[:-1])
    assert picked["x"][-1, 0] == -180.0
