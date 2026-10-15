"""Cabinet edge diffraction with the Distributed Edge Dipole (DED) model.

The DED model (Urban, Heil, Pignon, Combet, Bauman, JAES 52(10), 2004) replaces the
diffraction by the edges of the baffle with a line of dipoles along those edges, driven by
the pressure the driver produces there. In the form used here (extended to a finite-size
source), the pressure of one driver at the microphone ``M`` is

    p(M, f) = H(f) * [ P_dir(M, f) + P_edge(M, f) ]

    P_dir  = D(theta, f) * (1 + cos theta) / 2 * (r_ref / r) * exp(-j k r)
    P_edge = sum_e  Q_e(f) * F_e(M) * (r_ref / r2_e) * exp(-j k r2_e)

* ``H(f)``: on-axis half-space response at ``r_ref = 1 m``;
* ``r, theta``: distance and angle from the driver axis to ``M``; ``D``: radiator directivity;
  ``(1 + cos theta)/2``: obliquity of the direct field (1 on axis, 1/2 in the baffle plane,
  0 behind);
* ``e``: elements of length ``dl_e`` along the four edges of the front baffle, at ``E_e``;
  ``r2_e = |M - E_e|``;
* ``F_e = -cos(theta_e)``, ``theta_e`` the angle between the baffle normal and ``M - E_e``:
  the edge dipole, zero in the baffle plane and changing sign behind it;
* ``Q_e = C * dl_e * G(rho_e, f) / A0(f)`` with the driving term

      G(rho, f) = sum_s w_s K_s exp(-j k r1_s) / r1_s,   K_s = (1 + cos theta_s) / 2,

  summed over a quadrature of the radiating surface (weights ``w_s`` summing to 1, ``r1_s``
  from the surface point to the edge point, ``theta_s`` from the baffle normal), ``rho_e``
  the distance from the driver center to the edge point and ``A0`` the radiator's on-axis
  factor (1 for a piston);
* ``C = (1/2) / sum_e dl_e G(rho_e, 0)``: normalization that makes ``P_dir + P_edge = 1/2``
  on axis at low frequency in the far field. The cabinet then radiates 6 dB below the
  half-space response at low frequency (full-space radiation) and reaches it at high
  frequency: the baffle step.

All paths keep their true length: the edge waves arrive ``(r1 + r2 - r) / c`` after the direct
sound, at the actual microphone position. This finite-distance geometry is exactly what a
point-source correction cannot reproduce.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

from .constants import C, R_REF
from .geometry import CABINET, Cabinet


def front_edges(cabinet: Cabinet = CABINET, step: float = 0.002):
    """Midpoints and lengths of the elements along the front baffle perimeter (z = 0)."""
    w, h = cabinet.width / 2, cabinet.height / 2
    corners = [(-w, -h), (w, -h), (w, h), (-w, h)]
    points, lengths = [], []
    for (x0, y0), (x1, y1) in zip(corners, corners[1:] + corners[:1]):
        length = np.hypot(x1 - x0, y1 - y0)
        n = int(np.ceil(length / step))
        t = (np.arange(n) + 0.5) / n
        points.append(np.stack([x0 + t * (x1 - x0), y0 + t * (y1 - y0), np.zeros(n)], axis=-1))
        lengths.append(np.full(n, length / n))
    return np.concatenate(points), np.concatenate(lengths)


class EdgeDiffraction:
    """Precomputed edge strengths ``Q_e(f)`` of one driver."""

    def __init__(self, driver, freqs, cabinet: Cabinet = CABINET, step: float = 0.002,
                 coarse_df: float = 50.0):
        self.points, self.lengths = front_edges(cabinet, step)
        sx, sy = driver.placement.x, driver.placement.y
        rho_e = np.hypot(self.points[:, 0] - sx, self.points[:, 1] - sy)
        radiator = driver.radiator
        a = radiator.radius

        # Quadrature of the radiating surface: Gauss-Legendre in r, midpoints in phi over half a
        # turn (the edge point is put on the +x axis, the integrand is even in phi). Enough nodes
        # for the phase spread k_max * 2a across the surface at the top frequency.
        spread = 2.0 * np.pi * freqs[-1] / C * 2.0 * a
        n_r = int(np.ceil(0.7 * spread + 16))
        n_phi = int(np.ceil(0.5 * spread + 16))
        x, wx = np.polynomial.legendre.leggauss(n_r)
        r = 0.5 * a * (x + 1.0)
        wr = 0.5 * a * wx
        z = radiator.surface(r)
        phi = (np.arange(n_phi) + 0.5) * np.pi / n_phi
        rr, pp = np.meshgrid(r, phi, indexing="ij")
        zz = np.broadcast_to(z[:, None], rr.shape)
        w = (wr[:, None] * rr) * (np.pi / n_phi) * 2.0 / (np.pi * a**2)
        rr, pp, zz, w = rr.ravel(), pp.ravel(), zz.ravel(), w.ravel()

        # G(rho, f) * exp(+j k rho) on a 1 mm grid of rho and a coarse frequency grid: once the
        # common delay rho/c is removed it varies slowly in both.
        rho_grid = np.arange(rho_e.min() - 0.002, rho_e.max() + 0.003, 0.001)
        f_coarse = np.arange(0.0, freqs[-1] + coarse_df, coarse_df)
        k_coarse = 2.0 * np.pi * f_coarse / C
        g_tab = np.empty((len(rho_grid), len(f_coarse)), complex)
        for j, rho in enumerate(rho_grid):
            r1 = np.sqrt((rho - rr * np.cos(pp)) ** 2 + (rr * np.sin(pp)) ** 2 + zz**2)
            amp = w * 0.5 * (1.0 - zz / r1) / r1
            g_tab[j] = np.exp(-1j * np.outer(k_coarse, r1 - rho)) @ amp

        # To the edge points (linear in rho), then to the fine frequency grid (cubic).
        idx = np.clip(np.searchsorted(rho_grid, rho_e) - 1, 0, len(rho_grid) - 2)
        t = (rho_e - rho_grid[idx]) / (rho_grid[idx + 1] - rho_grid[idx])
        g_edges = (1.0 - t)[:, None] * g_tab[idx] + t[:, None] * g_tab[idx + 1]
        self.normalization = 0.5 / np.sum(self.lengths * g_edges[:, 0].real)
        re = CubicSpline(f_coarse, g_edges.real, axis=1)(freqs)
        im = CubicSpline(f_coarse, g_edges.imag, axis=1)(freqs)
        k = 2.0 * np.pi * freqs / C
        on_axis = radiator.on_axis(k)
        q = (self.normalization * self.lengths)[:, None] * (re + 1j * im) * np.exp(-1j * np.outer(rho_e, k))
        self.strength = (q / on_axis[None, :]).astype(np.complex64)

    def edge_term(self, propagation: np.ndarray) -> np.ndarray:
        """``P_edge`` given the edge-to-microphone factors ``F_e r_ref / r2_e exp(-j k r2_e)``."""
        return np.einsum("ef,ef->f", self.strength, propagation)


def edge_propagation(points, mic, k) -> np.ndarray:
    """``F_e(M) * r_ref / r2_e * exp(-j k r2_e)`` for every edge element, shape ``(n_e, n_f)``."""
    d = mic[None, :] - points
    r2 = np.linalg.norm(d, axis=1)
    dipole = -d[:, 2] / r2
    return (dipole * R_REF / r2)[:, None] * np.exp(-1j * np.outer(r2, k))


def direct_term(driver, mic, k) -> np.ndarray:
    """``P_dir``: directivity, obliquity, spreading and delay of the direct sound."""
    d = mic - np.array([driver.placement.x, driver.placement.y, 0.0])
    r = np.linalg.norm(d)
    cos_t = d[2] / r
    theta = np.arccos(np.clip(cos_t, -1.0, 1.0))
    return driver.radiator.directivity(k, theta) * 0.5 * (1.0 + cos_t) * (R_REF / r) * np.exp(-1j * k * r)
