"""Far-field directivity of the radiators.

* Cone: rigid flat circular piston of radius ``a`` in an infinite baffle,
  ``D(theta) = 2 J1(k a sin theta) / (k a sin theta)``.
* Dome: rigid convex spherical cap (base radius ``a``, height ``h``) moving axially. The far
  field is the Rayleigh integral over the projected disc, each ring of radius ``r`` being
  advanced by the cap height ``z(r)`` above the rim plane:

      A(theta, k) = (2 / a^2) * int_0^a J0(k r sin theta) exp(j k z(r) cos theta) r dr,
      z(r) = sqrt(R^2 - r^2) - (R - h),   R = (a^2 + h^2) / (2 h),

  and ``D(theta) = A(theta, k) / A(0, k)``: unity on axis, complex off axis because the
  protruding center moves the apparent source. ``h = 0`` gives the flat piston.

Both are normalized to 1 on their own axis; the on-axis response is the driver model's job.
"""

from __future__ import annotations

import numpy as np
from scipy.special import j0, j1


def piston(x):
    """``2 J1(x) / x`` with its limit 1 at ``x = 0``."""
    x = np.asarray(x, float)
    out = np.ones_like(x)
    m = np.abs(x) > 1e-8
    out[m] = 2.0 * j1(x[m]) / x[m]
    return out


class Radiator:
    """A piston (``height = 0``) or a dome, with Gauss-Legendre nodes over its radius."""

    def __init__(self, radius: float, height: float = 0.0, nodes: int = 96):
        self.radius = float(radius)
        self.height = float(height)
        x, w = np.polynomial.legendre.leggauss(nodes)
        self.r = 0.5 * radius * (x + 1.0)
        self.w = 0.5 * radius * w
        if height > 0.0:
            sphere = (radius**2 + height**2) / (2.0 * height)
            self.z = np.sqrt(sphere**2 - self.r**2) - (sphere - height)
        else:
            self.z = np.zeros_like(self.r)

    @property
    def is_piston(self) -> bool:
        return self.height <= 0.0

    def aperture(self, k, theta):
        """Unnormalized far-field factor ``A(theta, k)``; ``k`` and ``theta`` broadcast."""
        k = np.asarray(k, float)[..., None]
        theta = np.asarray(theta, float)[..., None]
        if self.is_piston:
            return piston(k[..., 0] * self.radius * np.sin(theta[..., 0])).astype(complex)
        terms = j0(k * self.r * np.sin(theta)) * np.exp(1j * k * self.z * np.cos(theta))
        return (2.0 / self.radius**2) * np.sum(terms * self.r * self.w, axis=-1)

    def on_axis(self, k):
        """``A(0, k)``: 1 for the piston, the phase-spread loss of the dome otherwise."""
        return self.aperture(k, np.zeros_like(np.asarray(k, float)))

    def directivity(self, k, theta):
        """``D(theta) = A(theta, k) / A(0, k)``."""
        if self.is_piston:
            return self.aperture(k, theta)
        return self.aperture(k, theta) / self.on_axis(k)
