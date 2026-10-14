"""2D and 3D schematics of the loudspeaker and of the two turntable setups.

Outputs (paper/figures): baffle_front.pdf, setups_top.pdf, setups_3d.pdf.
The 3D views use a small orthographic renderer: visible faces only (the cabinet is convex),
and ray casting to split the lines (rotation axis, microphone stand) into visible and hidden
parts.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from rotcenter import plotstyle as ps  # noqa: E402
from rotcenter.geometry import CABINET, CENTERED, MID, REFERENCE, RIGOROUS, TWEETER, WOOFER, mic_positions  # noqa: E402

FIG = ROOT / "paper" / "figures"
W, H, D = CABINET.width, CABINET.height, CABINET.depth
TURNTABLE_R = 0.20

FACE_FRONT = "#e9e7e1"
FACE_TOP = "#f5f4f0"
FACE_SIDE = "#d6d4cc"
EDGE = "#3d3c39"
CONE = "#2f2f2d"
FRAME = "#a9a79f"
DISC_TOP = "#eef2f7"
DISC_SIDE = "#c9d3df"
HIDDEN = "#b3b1aa"
AXIS_DASH = (0, (7, 2.2, 1.5, 2.2))

# Driver drawing radii (m): (outer frame/flange, cone or horn, dome or dust cap)
DRAW = {
    "woofer": (WOOFER, [(0.1124, FRAME), (0.089, CONE), (0.022, "#555451")]),
    "mid": (MID, [(0.078, FRAME), (0.060, "#77756f"), (0.0375, CONE)]),
    "tweeter": (TWEETER, [(0.045, FRAME), (0.018, "#77756f"), (0.014, CONE)]),
}


# ---------------------------------------------------------------------------------------------
# Front view of the baffle
# ---------------------------------------------------------------------------------------------

def dim_line(ax, p0, p1, text, offset=(0, 0), color=ps.INK_2, text_kw=None, fs=6.5):
    arrow = FancyArrowPatch(p0, p1, arrowstyle="<|-|>", mutation_scale=6, lw=0.6, color=color,
                            shrinkA=0, shrinkB=0)
    ax.add_patch(arrow)
    mid = ((p0[0] + p1[0]) / 2 + offset[0], (p0[1] + p1[1]) / 2 + offset[1])
    ax.text(*mid, text, ha="center", va="center", fontsize=fs, color=color,
            bbox=dict(fc="white", ec="none", pad=0.6), **(text_kw or {}))


def baffle_front():
    fig, ax = plt.subplots(figsize=(6.3, 4.3))
    ax.set_aspect("equal")
    ax.axis("off")
    mm = 1e3
    ax.add_patch(Rectangle((-W / 2 * mm, -H / 2 * mm), W * mm, H * mm, fc=FACE_FRONT, ec=EDGE, lw=1.0))
    for name, (pl, rings) in DRAW.items():
        for r, color in rings:
            ax.add_patch(Circle((pl.x * mm, pl.y * mm), r * mm, fc=color, ec=EDGE, lw=0.5))
    # Equivalent radiators used in the model.
    ax.add_patch(Circle((WOOFER.x * mm, WOOFER.y * mm), 83.8, fc="none", ec="white", lw=0.8, ls=(0, (3, 2))))
    ax.add_patch(Circle((MID.x * mm, MID.y * mm), 60.0, fc="none", ec="white", lw=0.8, ls=(0, (3, 2))))
    # Rotation axes (vertical lines in this view).
    ax.plot([REFERENCE.x * mm] * 2, [-H / 2 * mm - 35, H / 2 * mm + 30], color=ps.SLOTS[0], lw=1.1, ls=AXIS_DASH)
    ax.plot([0, 0], [-H / 2 * mm - 35, H / 2 * mm + 30], color=ps.SLOTS[1], lw=1.1, ls=AXIS_DASH)
    ax.text(REFERENCE.x * mm + 5, H / 2 * mm + 32, "rigorous axis\n(on the baffle)", color=ps.SLOTS[0],
            fontsize=6.5, ha="left", va="bottom")
    ax.text(-5, H / 2 * mm + 32, "centered axis\n(145 mm behind the baffle)", color="#c4501f",
            fontsize=6.5, ha="right", va="bottom")
    # Reference point and microphone height.
    ax.plot(REFERENCE.x * mm, REFERENCE.y * mm, marker="+", ms=9, mew=1.4, color="#e34948")
    ax.annotate("reference point\n(microphone axis)", (REFERENCE.x * mm, REFERENCE.y * mm),
                xytext=(REFERENCE.x * mm + 125, REFERENCE.y * mm + 5), fontsize=6.5, color=ps.INK,
                ha="left", va="center", arrowprops=dict(arrowstyle="-", lw=0.6, color=ps.INK_2))
    labels = {"woofer": "woofer\nSIG225-4\n(8 in cone)", "mid": "midrange\nDM-7600\n(3 in dome)",
              "tweeter": "tweeter\nLP 90.28/N92\n(1.1 in dome)"}
    spots = {"woofer": (-300, 90), "mid": (330, -95), "tweeter": (330, 110)}
    for name, (pl, rings) in DRAW.items():
        ax.annotate(labels[name], (pl.x * mm + (0 if name == "woofer" else 30), pl.y * mm),
                    xytext=spots[name], fontsize=6.5, ha="center", va="center", color=ps.INK,
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=ps.INK_2))
    # Dimensions (mm).
    y0 = -H / 2 * mm - 22
    dim_line(ax, (-220, y0 - 22), (220, y0 - 22), "440")
    dim_line(ax, (-220, y0), (WOOFER.x * mm, y0), "135")
    dim_line(ax, (WOOFER.x * mm, y0), (REFERENCE.x * mm, y0), "205")
    dim_line(ax, (REFERENCE.x * mm, y0), (220, y0), "100")
    x0 = -W / 2 * mm - 22
    dim_line(ax, (x0 - 22, -145), (x0 - 22, 145), "290", text_kw=dict(rotation=90))
    xr = W / 2 * mm + 22
    for (ya, yb, t) in [(-145, MID.y * mm, "96.5"), (MID.y * mm, REFERENCE.y * mm, "65"),
                        (REFERENCE.y * mm, TWEETER.y * mm, "65"), (TWEETER.y * mm, 145, "63.5")]:
        dim_line(ax, (xr, ya), (xr, yb), t, offset=(13, 0), text_kw=dict(rotation=90))
    for y in (MID.y * mm, REFERENCE.y * mm, TWEETER.y * mm):
        ax.plot([REFERENCE.x * mm, xr + 4], [y, y], color=ps.AXIS, lw=0.4)
    ax.set_xlim(-360, 420)
    ax.set_ylim(-235, 215)
    fig.savefig(FIG / "baffle_front.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Top views
# ---------------------------------------------------------------------------------------------

def top_panel(ax, turntable, title):
    mm = 1e3
    ax.set_aspect("equal")
    ax.axis("off")
    ax_x, ax_z = turntable.axis_x * mm, turntable.axis_z * mm
    ax.add_patch(Circle((ax_x, ax_z), TURNTABLE_R * mm, fc=DISC_TOP, ec="#8fa3bb", lw=0.8))
    ax.add_patch(Rectangle((-W / 2 * mm, -D * mm), W * mm, D * mm, fc=FACE_TOP, ec=EDGE, lw=1.0, alpha=0.92))
    # Drivers seen from above (on the baffle line).
    for x, half, depth in [(WOOFER.x, 84, 60), (MID.x, 60, 45), (TWEETER.x, 22, 30)]:
        ax.add_patch(Polygon([(x * mm - half, 0), (x * mm + half, 0), (x * mm + 0.35 * half, -depth),
                              (x * mm - 0.35 * half, -depth)], fc=CONE, ec="none", alpha=0.85))
    # Rotation axis.
    ax.add_patch(Circle((ax_x, ax_z), 16, fc="none", ec=ps.INK, lw=0.9, zorder=5))
    ax.plot([ax_x - 22, ax_x + 22], [ax_z, ax_z], color=ps.INK, lw=0.7, zorder=5)
    ax.plot([ax_x, ax_x], [ax_z - 22, ax_z + 22], color=ps.INK, lw=0.7, zorder=5)
    # Reference axis and microphone (distance not to scale).
    ref_x = REFERENCE.x * mm
    mic_z = 520
    ax.plot([ref_x, ref_x], [0, mic_z - 70], color=ps.INK_2, lw=0.7, ls=(0, (4, 2)))
    zig = np.array([[0, 0], [-8, 8], [8, 24], [0, 32]]) + [ref_x, mic_z - 70]
    ax.plot(zig[:, 0], zig[:, 1], color=ps.INK_2, lw=0.7)
    ax.plot([ref_x, ref_x], [mic_z - 38, mic_z - 14], color=ps.INK_2, lw=0.7, ls=(0, (4, 2)))
    ax.add_patch(Circle((ref_x, mic_z), 13, fc=ps.SLOTS[0], ec="white", lw=1.0))
    ax.text(ref_x + 22, mic_z, "microphone at 0°\n$d$ = 1 m or 2 m\nfrom the baffle", fontsize=6.3,
            va="center", color=ps.INK)
    ax.plot(ref_x, 0, marker="o", ms=3.0, color="#e34948", zorder=6)
    # Direction of positive angles: the microphone, seen from the speaker, turns toward +x.
    r = 330
    arc_t = np.radians(np.linspace(8, 32, 30))
    arc = np.stack([ax_x + r * np.sin(arc_t), ax_z + r * np.cos(arc_t)], axis=1)
    if turntable is CENTERED:
        arc = np.stack([ax_x + (r + 140) * np.sin(arc_t), ax_z + (r + 140) * np.cos(arc_t)], axis=1)
    ax.add_patch(FancyArrowPatch(arc[0], arc[-1], connectionstyle="arc3,rad=-0.12", arrowstyle="-|>",
                                 mutation_scale=7, lw=0.8, color=ps.SLOTS[0]))
    ax.text(arc[-1][0] + 8, arc[-1][1] - 10, "+θ (mic toward\nthe tweeter side)", fontsize=6.0, color=ps.SLOTS[0], va="top")
    ax.set_title(title, fontsize=7.5, loc="left")
    return ax_x, ax_z


def setups_top():
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.7))
    mm = 1e3
    top_panel(axes[0], RIGOROUS, "(a) rigorous: axis through the reference point")
    axes[0].annotate("axis on the baffle,\nbetween midrange\nand tweeter", (REFERENCE.x * mm + 14, -10),
                     xytext=(REFERENCE.x * mm + 120, -150), fontsize=6.0, color=ps.INK_2,
                     arrowprops=dict(arrowstyle="-", lw=0.5, color=ps.INK_2))
    ax_x, ax_z = top_panel(axes[1], CENTERED, "(b) small turntable: cabinet centered")
    ref_x = REFERENCE.x * mm
    dim_line(axes[1], (-W / 2 * mm - 25, ax_z), (-W / 2 * mm - 25, 0), "145", text_kw=dict(rotation=90), color="#c4501f")
    axes[1].plot([-W / 2 * mm - 30, ax_x], [ax_z, ax_z], color=ps.AXIS, lw=0.4)
    dim_line(axes[1], (ax_x, -D * mm - 55), (ref_x, -D * mm - 55), "120", color="#c4501f")
    axes[1].plot([ax_x, ax_x], [ax_z, -D * mm - 60], color=ps.AXIS, lw=0.4)
    axes[1].plot([ref_x, ref_x], [-D * mm, -D * mm - 60], color=ps.AXIS, lw=0.4)
    axes[1].text(-W / 2 * mm - 60, ax_z / 2, "depth", rotation=90, fontsize=6.0, color="#c4501f", va="center", ha="center")
    axes[1].text((ax_x + ref_x) / 2, -D * mm - 80, "lateral", fontsize=6.0, color="#c4501f", ha="center")
    dim_line(axes[1], (WOOFER.x * mm, 40), (ref_x, 40), "−205 (woofer offset)", color=ps.INK_2, fs=5.8)
    axes[1].plot([WOOFER.x * mm, WOOFER.x * mm], [0, 46], color=ps.AXIS, lw=0.4)
    for ax in axes:
        ax.set_xlim(-330, 470)
        ax.set_ylim(-470, 580)
    fig.savefig(FIG / "setups_top.pdf")
    plt.close(fig)


def positions_to_scale():
    """Microphone positions relative to the loudspeaker for both setups (d = 1 m, to scale)."""
    fig, ax = plt.subplots(figsize=(3.4, 3.6))
    ax.set_aspect("equal")
    angles = np.arange(-170, 181, 10)
    rig = mic_positions(RIGOROUS, angles, 1.0)
    cen = mic_positions(CENTERED, angles, 1.0)
    for a, b in zip(rig, cen):
        ax.plot([a[0], b[0]], [a[2], b[2]], color=ps.AXIS, lw=0.5, zorder=1)
    ax.plot(rig[:, 0], rig[:, 2], "o", ms=3.2, color=ps.SLOTS[0], label="rigorous (axis on baffle)", zorder=3)
    ax.plot(cen[:, 0], cen[:, 2], "s", ms=3.0, color=ps.SLOTS[1], label="centered cabinet", zorder=3)
    ax.add_patch(Rectangle((-W / 2, -D), W, D, fc=FACE_TOP, ec=EDGE, lw=0.8, zorder=2))
    ax.plot(REFERENCE.x, 0, "o", ms=2.5, color="#e34948", zorder=4)
    ax.plot(0, -D / 2, "+", ms=6, mew=1.0, color=ps.INK, zorder=4)
    for a in (0, 90, 180, -90):
        i = list(angles).index(a)
        txt = f"{a}°".replace("-", "\u2212")
        ax.annotate(txt, rig[i, [0, 2]], xytext=((rig[i, 0] - REFERENCE.x) * 1.32 + REFERENCE.x, rig[i, 2] * 1.22 + (0.05 if a == 0 else 0)),
                    fontsize=6.3, color=ps.INK, ha="center", va="center")
    ax.set_xlabel("x (m), + toward the tweeter side")
    ax.set_ylabel("z (m), toward the listener")
    ax.set_xlim(-1.55, 1.75)
    ax.set_ylim(-1.65, 1.45)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=6.3, handletextpad=0.3)
    ax.set_title("measurement positions seen from\nthe loudspeaker, $d$ = 1 m, 10° steps", fontsize=7.5, loc="left")
    fig.savefig(FIG / "positions_to_scale.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# 3D views
# ---------------------------------------------------------------------------------------------

class Ortho:
    """Orthographic camera; world frame X right, Y away from the listener, Z up (m)."""

    def __init__(self, view=(-0.50, -1.0, 0.50)):
        self.v = np.asarray(view, float) / np.linalg.norm(view)  # toward the camera
        right = np.cross([0.0, 0.0, 1.0], self.v)
        self.right = right / np.linalg.norm(right)
        self.up = np.cross(self.v, self.right)

    def p(self, pts):
        pts = np.atleast_2d(pts)
        return np.stack([pts @ self.right, pts @ self.up], axis=-1)


class Scene:
    def __init__(self, turntable, lift=0.34, disc_t=0.035):
        self.cam = Ortho()
        self.lift = lift
        self.disc_t = disc_t
        # World coordinates of the turntable axis.
        self.ax_xy = (turntable.axis_x, -turntable.axis_z)
        self.box_min = np.array([-W / 2, 0.0, lift])
        self.box_max = np.array([W / 2, D, lift + H])

    def world(self, x, y, z_speaker=0.0):
        """Speaker-frame point (x, y up, z toward listener) to world."""
        return np.array([x, -z_speaker, y + H / 2 + self.lift])

    def hidden(self, pt):
        """Is ``pt`` hidden from the camera by the cabinet or the turntable?"""
        for s in np.linspace(0.004, 3.0, 1500):
            q = pt + s * self.cam.v
            if np.all(q > self.box_min + 1e-6) and np.all(q < self.box_max - 1e-6):
                return True
            dx, dy = q[0] - self.ax_xy[0], q[1] - self.ax_xy[1]
            if dx * dx + dy * dy < TURNTABLE_R**2 and 0.0 < q[2] < self.disc_t:
                return True
        return False

    def line(self, ax, p0, p1, color, lw, ls, hidden_color=HIDDEN, n=240, zorder=6, draw_hidden=True):
        pts = np.linspace(p0, p1, n)
        flags = np.array([self.hidden(q) for q in pts])
        start = 0
        for i in range(1, n + 1):
            if i == n or flags[i] != flags[start]:
                seg = self.cam.p(pts[start:i + (0 if i == n else 1)])
                if flags[start]:
                    if draw_hidden:
                        ax.plot(seg[:, 0], seg[:, 1], color=hidden_color, lw=lw * 0.7, ls=ls, zorder=zorder - 3)
                else:
                    ax.plot(seg[:, 0], seg[:, 1], color=color, lw=lw, ls=ls, zorder=zorder)
                start = i

    def polygon(self, ax, pts, **kw):
        ax.add_patch(Polygon(self.cam.p(np.asarray(pts)), closed=True, **kw))

    def disc(self, ax):
        cx, cy = self.ax_xy
        t = np.linspace(0, 2 * np.pi, 200)
        ring = lambda z: np.stack([cx + TURNTABLE_R * np.cos(t), cy + TURNTABLE_R * np.sin(t), np.full_like(t, z)], 1)
        top, bottom = ring(self.disc_t), ring(0.0)
        # Side band facing the camera.
        facing = (np.cos(t) * self.cam.v[0] + np.sin(t) * self.cam.v[1]) > 0
        idx = np.where(facing)[0]
        # Make the facing arc contiguous (it may wrap around t = 0).
        if idx.size and idx[0] == 0 and idx[-1] == len(t) - 1:
            gap = np.where(np.diff(idx) > 1)[0]
            if gap.size:
                idx = np.concatenate([idx[gap[0] + 1:], idx[:gap[0] + 1]])
        band = np.concatenate([top[idx], bottom[idx][::-1]])
        self.polygon(ax, band, fc=DISC_SIDE, ec="#7f93ad", lw=0.6, zorder=1)
        self.polygon(ax, top, fc=DISC_TOP, ec="#7f93ad", lw=0.7, zorder=1)
        # Spokes on the top, to read it as a turntable.
        for a in np.radians([0, 60, 120]):
            p0 = [cx + 0.92 * TURNTABLE_R * np.cos(a), cy + 0.92 * TURNTABLE_R * np.sin(a), self.disc_t]
            p1 = [cx - 0.92 * TURNTABLE_R * np.cos(a), cy - 0.92 * TURNTABLE_R * np.sin(a), self.disc_t]
            q = self.cam.p(np.array([p0, p1]))
            ax.plot(q[:, 0], q[:, 1], color="#d3dce7", lw=0.6, zorder=1)

    def cabinet(self, ax):
        x0, y0, z0 = self.box_min
        x1, y1, z1 = self.box_max
        front = [[x0, y0, z0], [x1, y0, z0], [x1, y0, z1], [x0, y0, z1]]
        top = [[x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]]
        left = [[x0, y0, z0], [x0, y1, z0], [x0, y1, z1], [x0, y0, z1]]
        for face, color in ((left, FACE_SIDE), (top, FACE_TOP), (front, FACE_FRONT)):
            self.polygon(ax, face, fc=color, ec=EDGE, lw=0.8, zorder=3, joinstyle="round")
        t = np.linspace(0, 2 * np.pi, 120)
        for name, (pl, rings) in DRAW.items():
            for r, color in rings:
                pts = np.stack([pl.x + r * np.cos(t), np.zeros_like(t), pl.y + H / 2 + self.lift + r * np.sin(t)], 1)
                self.polygon(ax, pts, fc=color, ec=EDGE, lw=0.35, zorder=4)

    def microphone(self, ax, distance_draw):
        ref = self.world(REFERENCE.x, REFERENCE.y, distance_draw)
        tip = ref + np.array([0, 0.035, 0])
        tail = ref - np.array([0, 0.09, 0])
        q = self.cam.p(np.array([tail, tip]))
        ax.plot(q[:, 0], q[:, 1], color="#55606e", lw=4.2, solid_capstyle="round", zorder=9)
        qt = self.cam.p(tip)[0]
        ax.plot(qt[0], qt[1], "o", ms=4.5, color="#2b2f36", zorder=10)
        # Stand down to the floor.
        foot = np.array([ref[0], ref[1] - 0.04, 0.0])
        q = self.cam.p(np.array([ref - np.array([0, 0.04, 0]), foot]))
        ax.plot(q[:, 0], q[:, 1], color="#7c838c", lw=0.9, zorder=8)
        for a in (0, 120, 240):
            leg = foot + 0.07 * np.array([np.cos(np.radians(a + 30)), np.sin(np.radians(a + 30)), 0])
            q = self.cam.p(np.array([foot + [0, 0, 0.08], leg]))
            ax.plot(q[:, 0], q[:, 1], color="#7c838c", lw=0.7, zorder=8)
        return ref

    def label(self, ax, pt3, text, dxy, ha="left", color=ps.INK):
        p = self.cam.p(pt3)[0]
        ax.annotate(text, p, xytext=(p[0] + dxy[0], p[1] + dxy[1]), fontsize=6.3, ha=ha, va="center", color=color,
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=ps.INK_2, shrinkA=1, shrinkB=1), zorder=12)


def draw_3d(ax, turntable, title):
    sc = Scene(turntable)
    ax.set_aspect("equal")
    ax.axis("off")
    sc.disc(ax)
    sc.cabinet(ax)
    # Rotation axis: from the turntable center to above the cabinet. Hidden parts (inside the
    # cabinet or behind it) in a lighter tone.
    cx, cy = sc.ax_xy
    bottom = np.array([cx, cy, sc.disc_t])
    top = np.array([cx, cy, sc.box_max[2] + 0.16])
    sc.line(ax, bottom, top, color="#c62828", lw=1.4, ls=AXIS_DASH, hidden_color="#e9a3a3", zorder=7)
    q = sc.cam.p(bottom)[0]
    ax.plot(*q, "o", ms=3.2, color="#c62828", zorder=8)
    if turntable is CENTERED:
        for z in (sc.box_max[2],):
            q = sc.cam.p(np.array([cx, cy, z]))[0]
            ax.plot(*q, "o", ms=3.2, mfc="white", mec="#c62828", mew=1.0, zorder=8)
    # Reference axis and microphone (distance not to scale).
    mic = sc.microphone(ax, 0.55)
    ref_pt = sc.world(REFERENCE.x, REFERENCE.y, 0.0)
    sc.line(ax, ref_pt, mic + np.array([0, 0.09, 0]), color=ps.SLOTS[0], lw=0.9, ls=(0, (2, 1.6)), zorder=7,
            draw_hidden=False)
    q = sc.cam.p(ref_pt)[0]
    ax.plot(*q, "o", ms=3.0, color="#e34948", mec="white", mew=0.5, zorder=9)
    sc.label(ax, top - np.array([0, 0, 0.02]), "rotation axis", (0.07, 0.0), color="#c62828")
    sc.label(ax, np.array([cx + TURNTABLE_R * 0.75, cy - TURNTABLE_R * 0.6, sc.disc_t * 0.5]), "turntable", (0.20, -0.06))
    sc.label(ax, mic, "microphone at 0°", (0.12, 0.10))
    sc.label(ax, ref_pt, "reference point", (0.42, 0.17))
    ax.set_title(title, fontsize=7.5, loc="left")
    ax.set_xlim(-0.72, 0.98)
    ax.set_ylim(-0.38, 0.98)


def setups_3d():
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
    draw_3d(axes[0], RIGOROUS, "(a) rigorous: the axis runs along the baffle,\n     through midrange and tweeter")
    draw_3d(axes[1], CENTERED, "(b) small turntable: the axis runs through\n     the middle of the cabinet")
    fig.subplots_adjust(wspace=0.04, left=0.01, right=0.99, top=0.9, bottom=0.01)
    fig.savefig(FIG / "setups_3d.pdf")
    fig.savefig(FIG / "setups_3d.png", dpi=200)
    plt.close(fig)


def main():
    ps.apply()
    plt.rcParams["axes.grid"] = False
    FIG.mkdir(parents=True, exist_ok=True)
    baffle_front()
    setups_top()
    positions_to_scale()
    setups_3d()


if __name__ == "__main__":
    main()
