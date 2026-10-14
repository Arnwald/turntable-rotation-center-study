"""Driver models: on-axis responses that follow the datasheets, and radiator geometry.

Each driver is described by its on-axis free-field response in an infinite baffle (half
space), at 1 m for 2.83 V, ``H(f)`` (complex, minimum phase), and by a radiator
(piston or dome) giving its directivity.

* The datasheet SPL curves (extracted from the vector graphics of the PDF datasheets) are
  taken as the half-space responses. Each is interpolated in log-frequency, extended outside
  the published range with fixed slopes, and given a minimum phase.
* The woofer is moved from the datasheet's (near infinite-baffle) loading to the 24 L closed
  box by the ratio of two second-order high-pass functions,
  ``H_cb(s) / H_ib(s)`` with ``H(s) = s^2 / (s^2 + s w0/Q + w0^2)``:
  ``(fs, Qts)`` for the datasheet, ``(fc, Qtc) = (fs, Qts) * sqrt(1 + Vas/Vb)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .constants import P_REF
from .directivity import Radiator
from .geometry import MID, TWEETER, WOOFER, Placement
from .signal import minimum_phase_from_db, raised_cosine_taper, rfft_freqs

DATA = Path(__file__).resolve().parent.parent / "data" / "datasheets"

#: Thiele-Small parameters of the woofer (Dayton Audio SIG225-4 datasheet).
SIG225 = dict(fs=40.0, qts=0.43, qes=0.45, qms=8.87, vas_l=29.7, re=3.9, le_mh=0.49, bl=9.0,
              mms_g=36.9, cms_mm_n=0.43, sd_cm2=220.8, spl=91.4)

#: Closed box net volume (L).
BOX_LITERS = 24.0

#: Radiators. Woofer: flat piston of area Sd. Tweeter: 28 mm silk dome (1.1 in voice coil),
#: assumed 7.5 mm high. Midrange: the DM-7600 sits in a horn-shaped faceplate; a 75 mm dome
#: alone is much wider than the datasheet's 30-degree curve, so an equivalent cap fitted to that
#: curve is used (120 mm across, 17.5 mm high: the mouth of the faceplate).
WOOFER_RADIATOR = dict(radius=float(np.sqrt(SIG225["sd_cm2"] * 1e-4 / np.pi)), height=0.0)
TWEETER_RADIATOR = dict(radius=0.014, height=0.0075)
MID_RADIATOR = dict(radius=0.060, height=0.0175)


@dataclass
class Driver:
    name: str
    label: str
    placement: Placement
    radiator: Radiator
    #: On-axis half-space response at 1 m for 2.83 V (Pa), on the rfft grid.
    response: np.ndarray


def load_datasheet(name: str):
    data = np.loadtxt(DATA / f"{name}.csv", delimiter=",", skiprows=2)
    return data[:, 0], data[:, 1]


def second_order_highpass(f, f0, q):
    s = 1j * np.asarray(f, float) / f0
    return s**2 / (s**2 + s / q + 1.0)


def closed_box_alignment(fs=SIG225["fs"], qts=SIG225["qts"], vas_l=SIG225["vas_l"], vb_l=BOX_LITERS):
    """System resonance and Q of the woofer in the closed box (lossless box)."""
    ratio = np.sqrt(1.0 + vas_l / vb_l)
    return fs * ratio, qts * ratio


def _response_from_datasheet(name, freqs, below, above):
    """Minimum-phase response (Pa at 1 m) from a datasheet SPL curve, with fixed-slope
    extensions (dB/octave) below and above the published range."""
    f, spl = load_datasheet(name)
    pa_db = spl + 20.0 * np.log10(P_REF)
    return minimum_phase_from_db(f, pa_db, freqs, below=below, above=above)


def woofer_response(freqs):
    """Datasheet response moved into the closed box: ``H_ds(f) * H_cb(f) / H_ib(f)``.

    Both high-pass functions are minimum phase, so the product stays minimum phase."""
    fc, qtc = closed_box_alignment()
    h = _response_from_datasheet("dayton_sig225_4_spl", freqs, below=12.0, above=-24.0)
    f = np.maximum(freqs, 1e-3)
    return h * second_order_highpass(f, fc, qtc) / second_order_highpass(f, SIG225["fs"], SIG225["qts"])


def mid_response(freqs):
    return _response_from_datasheet("hivi_dm7600_spl_0deg", freqs, below=12.0, above=-12.0)


def tweeter_response(freqs):
    return _response_from_datasheet("sica_lp9028_n92_spl", freqs, below=12.0, above=-12.0)


def build_drivers(freqs=None) -> dict[str, Driver]:
    """The three drivers on the default rfft grid, band-limited below 23.5 kHz."""
    freqs = rfft_freqs() if freqs is None else freqs
    taper = raised_cosine_taper(freqs, 20_500.0, 23_500.0)
    specs = [
        ("woofer", "SIG225-4 (8 in cone)", WOOFER, WOOFER_RADIATOR, woofer_response),
        ("mid", "DM-7600 (3 in dome)", MID, MID_RADIATOR, mid_response),
        ("tweeter", "LP 90.28/N92 (1.1 in dome)", TWEETER, TWEETER_RADIATOR, tweeter_response),
    ]
    out = {}
    for name, label, placement, radiator, response in specs:
        h = response(freqs) * taper
        h[0] = 0.0
        out[name] = Driver(name, label, placement, Radiator(**radiator), h)
    return out
