"""Loudspeaker layout and turntable kinematics.

Speaker frame (metres): origin at the center of the front baffle, x horizontal and positive to
the right as seen from the listening position, y up, z normal to the baffle and positive toward
the listener. The cabinet occupies ``-depth <= z <= 0``; every driver sits on the baffle plane
``z = 0``.

The loudspeaker is a horizontal three-way monitor: the woofer is on the left, the midrange and
the tweeter are stacked on the right (tweeter on top). The *reference point* is on the baffle,
on the vertical line through the midrange and tweeter centers, halfway between them; the
*reference axis* is the baffle normal through that point. The microphone is always on the
reference axis at 0 degrees.

Angles: a positive turntable angle brings the microphone, seen from the loudspeaker, toward
+x (the midrange/tweeter side).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Cabinet:
    width: float = 0.440
    height: float = 0.290
    depth: float = 0.290


@dataclass(frozen=True)
class Placement:
    """Center of a driver on the baffle (m)."""

    x: float
    y: float


CABINET = Cabinet()

#: Driver centers. Woofer frame 224.8 mm, midrange flange 156 mm, tweeter flange 90 mm: about
#: 20 mm of baffle around every flange.
WOOFER = Placement(-0.085, 0.0)
MID = Placement(0.120, -0.0485)
TWEETER = Placement(0.120, 0.0815)

#: Reference point on the baffle: midway between midrange and tweeter.
REFERENCE = Placement(0.120, (MID.y + TWEETER.y) / 2)


@dataclass(frozen=True)
class Turntable:
    """A vertical rotation axis, located by its (x, z) position in the speaker frame."""

    name: str
    axis_x: float
    axis_z: float


#: The textbook setup: the axis goes through the reference point on the baffle.
RIGOROUS = Turntable("rigorous", REFERENCE.x, 0.0)

#: The small-turntable setup: the cabinet is centered on the turntable.
CENTERED = Turntable("centered", 0.0, -CABINET.depth / 2)


def mic_positions(turntable: Turntable, angles_deg, distance: float, height: float | None = None):
    """Microphone positions in the speaker frame, shape ``(len(angles), 3)``.

    At 0 degrees the microphone is on the reference axis, ``distance`` in front of the baffle, at
    ``height`` (the reference point's height by default). Rotating the loudspeaker by ``theta``
    about the turntable axis moves the (fixed) microphone, seen from the loudspeaker, by
    ``-theta`` about that axis, i.e. toward +x for positive angles:

        x' = a_x + (x0 - a_x) cos(theta) + (z0 - a_z) sin(theta)
        z' = a_z - (x0 - a_x) sin(theta) + (z0 - a_z) cos(theta)
    """
    theta = np.radians(np.atleast_1d(np.asarray(angles_deg, dtype=float)))
    y = REFERENCE.y if height is None else height
    x0 = REFERENCE.x - turntable.axis_x
    z0 = distance - turntable.axis_z
    x = turntable.axis_x + x0 * np.cos(theta) + z0 * np.sin(theta)
    z = turntable.axis_z - x0 * np.sin(theta) + z0 * np.cos(theta)
    return np.stack([x, np.full_like(x, y), z], axis=-1)


def correction_offsets(turntable: Turntable, driver: Placement, distance: float) -> dict:
    """Parameters of the rotation-center correction for one driver in the horizontal plane (mm).

    ``depth``: rotation axis to front baffle, positive when the axis is behind the baffle.
    ``lateral``: rotation axis to reference axis, positive toward the + side.
    ``driver_offset``: reference axis to the driver's acoustic center, positive toward the + side.
    ``mic_distance``: front baffle to microphone at 0 degrees.
    """
    return {
        "depth_mm": -turntable.axis_z * 1e3,
        "lateral_mm": (REFERENCE.x - turntable.axis_x) * 1e3,
        "driver_offset_mm": (driver.x - REFERENCE.x) * 1e3,
        "mic_distance_mm": distance * 1e3,
    }
