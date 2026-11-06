"""Figures of the two articles drawn from the study (a blog post and a magazine article).

Run ``scripts/run_study.py`` and ``scripts/run_bem.py`` first. Each figure is written to
``articles/figures/<name>.svg`` (web) and ``<name>.pdf`` and ``<name>.png`` (print, 300 dpi), with
larger type than the report's figures.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import make_schematics as sch  # noqa: E402
from rotcenter import plotstyle as ps  # noqa: E402
from rotcenter import twoway  # noqa: E402
from rotcenter.constants import C, P_REF  # noqa: E402
from rotcenter.correction import PlaneOffsets, PlaneSetup  # noqa: E402
from rotcenter.cta2034 import cta2034, horizontal_cta2034  # noqa: E402
from rotcenter.geometry import CABINET, CENTERED, RIGOROUS, Placement  # noqa: E402
from rotcenter.metrics import BANDS, level_error, per_angle_rms, relevant  # noqa: E402
from rotcenter.signal import smooth_power  # noqa: E402

OUT = ROOT / "articles" / "figures"
STUDY = np.load(ROOT / "results" / "study.npz")
BEM = np.load(ROOT / "results" / "bem.npz")
F3, A3 = STUDY["freqs"], STUDY["angles"]
F2 = BEM["freqs"]
ANGLES = {"": BEM["angles"], "5": BEM["angles5"], "2": BEM["angles2"]}
W, H, D = CABINET.width * 1e3, CABINET.height * 1e3, CABINET.depth * 1e3
SPL0 = 10.0 * np.log10(P_REF**2)

RIG, RAW, COR, ALT = ps.SLOTS[0], ps.SLOTS[1], ps.SLOTS[2], ps.SLOTS[6]
CASE = {"rigorous": ("axis through the reference point", RIG, "-", 1.7),
        "raw": ("cabinet centered, uncorrected", RAW, (0, (1.2, 1.4)), 1.6),
        "corrected": ("cabinet centered, corrected", COR, (0, (5, 2)), 1.5)}


def style():
    ps.apply()
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 9.5,
        "axes.titlesize": 10,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.8,
        "lines.linewidth": 1.5,
        "svg.fonttype": "none",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.svg")
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=300)
    plt.close(fig)


def scale_text(fig, k):
    for t in fig.findobj(mpl.text.Text):
        t.set_fontsize(t.get_fontsize() * k)


def spl(p):
    return ps.db10(p) - SPL0


def octave_smooth(freqs, power, fraction):
    return smooth_power(freqs, power, fraction, freqs)


def bem(case, tag="1m", drv="system", model="bem"):
    return BEM[f"{model}/{tag}/{case}/{drv}/power"].astype(float)


def study(key, drv="system"):
    return STUDY[f"{key}/{drv}/power"]


def suffix(case):
    return case[-1] if case[-1] in "52" and not case.endswith("m") else ""


def at_ten(p, angles):
    return {int(round(a)): p[i] for i, a in enumerate(angles) if np.isclose(a % 10, 0) or np.isclose(a % 10, 10)}


# ---------------------------------------------------------------------------------------------
# Setups
# ---------------------------------------------------------------------------------------------

TWO_WAY_DRAW = {
    "woofer": (twoway.WOOFER, [(0.1124, sch.FRAME), (0.089, sch.CONE), (0.022, "#555451")]),
    "tweeter": (twoway.TWEETER, [(0.045, sch.FRAME), (0.018, "#77756f"), (0.014, sch.CONE)]),
}


def setups_3d(kind):
    """The two setups in perspective, for the three-way (``"three"``) or the two-way (``"two"``)."""
    draw, reference = sch.DRAW, sch.REFERENCE
    through = "through midrange and tweeter"
    if kind == "two":
        sch.DRAW, sch.REFERENCE = TWO_WAY_DRAW, twoway.REFERENCE
        through = "through the tweeter"
    plt.rcParams["axes.grid"] = False
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
    sch.draw_3d(axes[0], RIGOROUS, f"(a) the axis runs along the baffle,\n     {through}")
    sch.draw_3d(axes[1], CENTERED, "(b) small turntable: the axis runs\n     through the middle of the cabinet")
    fig.subplots_adjust(wspace=0.04, left=0.01, right=0.99, top=0.9, bottom=0.01)
    scale_text(fig, 1.32)
    for ax in axes:
        ax.title.set_fontsize(10)
    sch.DRAW, sch.REFERENCE = draw, reference
    save(fig, f"setups-{kind}-way")
    plt.rcParams["axes.grid"] = True


def dim(ax, p0, p1, text, color, offset=(0, 0), rotation=0, fs=9):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="<|-|>", mutation_scale=8, lw=0.8, color=color,
                                 shrinkA=0, shrinkB=0))
    mid = ((p0[0] + p1[0]) / 2 + offset[0], (p0[1] + p1[1]) / 2 + offset[1])
    ax.text(*mid, text, ha="center", va="center", fontsize=fs, color=color, rotation=rotation,
            bbox=dict(fc="white", ec="none", pad=0.8))


def top_view():
    """Setup (b) seen from above, with what the correction needs."""
    plt.rcParams["axes.grid"] = False
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.set_aspect("equal")
    ax.axis("off")
    ax_x, ax_z = 0.0, -145.0
    ref_x = 120.0
    accent = "#c4501f"
    ax.add_patch(Circle((ax_x, ax_z), 200, fc=sch.DISC_TOP, ec="#8fa3bb", lw=0.9))
    ax.annotate("turntable", (ax_x - 120, ax_z - 165), xytext=(ax_x - 300, ax_z - 250), fontsize=9,
                color="#6d7f96", ha="center", va="center", arrowprops=dict(arrowstyle="-", lw=0.6, color="#8fa3bb"))
    ax.add_patch(Rectangle((-W / 2, -D), W, D, fc=sch.FACE_TOP, ec=sch.EDGE, lw=1.1, alpha=0.93))
    for x, half, depth in ((-85.0, 84, 60), (ref_x, 22, 30)):
        ax.add_patch(Polygon([(x - half, 0), (x + half, 0), (x + 0.35 * half, -depth), (x - 0.35 * half, -depth)],
                             fc=sch.CONE, ec="none", alpha=0.85))
    ax.text(-85, -95, "woofer", fontsize=9, ha="center", color=ps.INK)
    ax.text(ref_x, -65, "tweeter", fontsize=9, ha="center", color=ps.INK)
    # Rotation axis.
    ax.add_patch(Circle((ax_x, ax_z), 15, fc="white", ec=ps.INK, lw=1.0, zorder=5))
    ax.plot([ax_x - 22, ax_x + 22], [ax_z, ax_z], color=ps.INK, lw=0.8, zorder=5)
    ax.plot([ax_x, ax_x], [ax_z - 22, ax_z + 22], color=ps.INK, lw=0.8, zorder=5)
    ax.annotate("rotation axis", (ax_x + 14, ax_z - 12), xytext=(W / 2 + 40, ax_z - 90), fontsize=9,
                color=ps.INK, ha="left", va="center", arrowprops=dict(arrowstyle="-", lw=0.6, color=ps.INK_2))
    # Reference axis and microphone (distance not to scale).
    mic_z = 600
    ax.plot([ref_x, ref_x], [0, mic_z - 80], color=ps.INK_2, lw=0.8, ls=(0, (4, 2)))
    zig = np.array([[0, 0], [-9, 9], [9, 27], [0, 36]]) + [ref_x, mic_z - 80]
    ax.plot(zig[:, 0], zig[:, 1], color=ps.INK_2, lw=0.8)
    ax.plot([ref_x, ref_x], [mic_z - 44, mic_z - 15], color=ps.INK_2, lw=0.8, ls=(0, (4, 2)))
    ax.add_patch(Circle((ref_x, mic_z), 15, fc=ps.SLOTS[0], ec="white", lw=1.0, zorder=6))
    ax.text(ref_x + 28, mic_z, "microphone at 0°,\non the reference axis", fontsize=9, va="center", color=ps.INK)
    ax.plot(ref_x, 0, marker="o", ms=4.5, color="#e34948", zorder=6)
    ax.annotate("reference point", (ref_x + 4, 4), xytext=(ref_x + 120, 95), fontsize=9, color=ps.INK,
                arrowprops=dict(arrowstyle="-", lw=0.6, color=ps.INK_2))
    # The four inputs.
    x_dim = -W / 2 - 40
    ax.plot([x_dim - 6, ax_x - 18], [ax_z, ax_z], color=ps.AXIS, lw=0.5)
    ax.plot([x_dim - 6, -W / 2], [0, 0], color=ps.AXIS, lw=0.5)
    dim(ax, (x_dim, ax_z), (x_dim, 0), "depth\n145 mm", accent, offset=(-62, 0))
    y_dim = -D - 60
    ax.plot([ax_x, ax_x], [ax_z - 18, y_dim - 8], color=ps.AXIS, lw=0.5)
    ax.plot([ref_x, ref_x], [-D, y_dim - 8], color=ps.AXIS, lw=0.5)
    dim(ax, (ax_x, y_dim), (ref_x, y_dim), "lateral offset +120 mm", accent, offset=(60, -36))
    ax.plot([-85, -85], [0, 62], color=ps.AXIS, lw=0.5)
    dim(ax, (-85, 50), (ref_x, 50), "woofer offset −205 mm", ps.INK_2, offset=(0, 26))
    x_mic = -W / 2 + 40
    ax.plot([x_mic - 8, ref_x - 18], [mic_z, mic_z], color=ps.AXIS, lw=0.5)
    dim(ax, (x_mic, 0), (x_mic, mic_z), "", ps.INK_2)
    ax.text(x_mic - 14, 330, "microphone distance\n(1 m or 2 m)", fontsize=9, color=ps.INK_2, ha="right", va="center")
    # Positive angles.
    r = 520
    t = np.radians(np.linspace(14, 34, 30))
    arc = np.stack([ax_x + r * np.sin(t), ax_z + r * np.cos(t)], axis=1)
    ax.add_patch(FancyArrowPatch(arc[0], arc[-1], connectionstyle="arc3,rad=-0.12", arrowstyle="-|>",
                                 mutation_scale=10, lw=1.0, color=ps.SLOTS[0]))
    ax.text(arc[-1][0] + 12, arc[-1][1] - 8, "+θ: the microphone,\nseen from the speaker,\nmoves to the + side",
            fontsize=9, color=ps.SLOTS[0], va="top")
    ax.set_xlim(-430, 560)
    ax.set_ylim(-470, 660)
    save(fig, "top-view")
    plt.rcParams["axes.grid"] = True


def errors_geometry():
    """What the centered setup does to the measurements, seen from the drivers."""
    angles = np.arange(-180, 181, 1.0)
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))
    for dist, ls in ((1.0, "-"), (2.0, (0, (4, 2)))):
        for name, placement, color in (("tweeter", twoway.TWEETER, ps.SLOTS[0]), ("woofer", twoway.WOOFER, ps.SLOTS[1])):
            g = PlaneSetup.from_offsets(PlaneOffsets(**twoway.correction_offsets(placement, dist, "h")))
            meas = np.array([g.measured(a) for a in angles])
            ideal = np.array([g.ideal(a) for a in angles])
            axes[0].plot(angles, 20 * np.log10(ideal[:, 1] / meas[:, 1]), color=color, ls=ls, label=f"{name}, {dist:.0f} m")
            axes[1].plot(angles, (meas[:, 1] - ideal[:, 1]) / C * 1e3, color=color, ls=ls)
            axes[2].plot(angles, np.degrees(meas[:, 0] - ideal[:, 0]), color=color, ls=ls)
    for ax, title in zip(axes, ("level (dB)", "extra delay (ms)", "direction (°)")):
        ax.set_xlim(-180, 180)
        ax.set_xticks(np.arange(-180, 181, 90))
        ax.set_xlabel("nominal angle (°)")
        ax.set_title(title, loc="left")
    fig.tight_layout(w_pad=1.0, rect=(0, 0.12, 1, 1))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, handlelength=2.6)
    save(fig, "errors-geometry")


# ---------------------------------------------------------------------------------------------
# The two-way, BEM
# ---------------------------------------------------------------------------------------------

def angle_index(angles, a):
    a = -180.0 if a == 180 and not np.isclose(angles, 180).any() else a
    return int(np.flatnonzero(np.isclose(angles, a))[0])


def responses(tag="1m", picks=(60, 120, 150)):
    """The summed two-way at a few angles: rigorous, uncorrected, corrected (10-degree steps)."""
    fig, axes = plt.subplots(1, len(picks), figsize=(7.2, 2.9), sharey=True)
    for ax, a in zip(axes, picks):
        i = angle_index(ANGLES[""], a)
        for case in ("rigorous", "raw", "corrected"):
            label, color, ls, lw = CASE[case]
            p = octave_smooth(F2, bem(case, tag)[i], 48 / 4)
            ax.semilogx(F2, spl(p), color=color, ls=ls, lw=lw, label=label)
        ps.freq_axis(ax, 100, 16000, sparse=False)
        ax.set_xticks([100, 1000, 10000])
        ax.set_xticklabels(["100", "1k", "10k"])
        ax.set_title(f"{a}°", loc="left")
        ax.set_xlabel("frequency (Hz)")
    axes[0].set_ylabel(f"SPL at {tag[0]} m (dB)")
    axes[0].set_ylim(50, 92)
    fig.tight_layout(w_pad=0.6, rect=(0, 0.1, 1, 1))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, handlelength=2.6)
    save(fig, f"responses-{tag}")


def cta_both(case_h, case_v, tag="1m", model="bem"):
    sh, sv = suffix(case_h), suffix(case_v)
    return cta2034(at_ten(bem(case_h, tag, model=model), ANGLES[sh]), at_ten(bem(case_v, tag, model=model), ANGLES[sv]))


def spatial(tag="1m", upright=False):
    """PIR and sound power from both orbits, against the rigorous ones. ``upright``: the two-way
    read as an ordinary two-way standing on end, the orbit through both drivers its vertical one
    (measured lying on its side) and the other its horizontal one."""
    def both(h, v):
        return cta_both(v, h, tag) if upright else cta_both(h, v, tag)
    ref = both("rigorous", "v-rigorous")
    curves = (("raw", "v-raw", "uncorrected", RAW, (0, (1.2, 1.4)), 1.6),
              ("corrected", "v-corrected", "corrected, 10° steps", COR, "-", 1.5),
              ("corrected2", "v-corrected2", "corrected, 2.5° steps", ALT, (0, (5, 2)), 1.3))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
    for ax, key, title in zip(axes, ("pir", "sp"), ("Predicted In-Room Response", "sound power")):
        for h, v, label, color, ls, lw in curves:
            c = both(h, v)
            ax.semilogx(F2, ps.db10(c[key] / ref[key]), color=color, ls=ls, lw=lw, label=label)
        ax.axhline(0, color=ps.AXIS, lw=0.8)
        ps.freq_axis(ax, 100, 16000)
        ax.set_xticks([100, 1000, 10000])
        ax.set_xticklabels(["100", "1k", "10k"])
        ax.set_title(title, loc="left")
        ax.set_xlabel("frequency (Hz)")
    axes[0].set_ylabel("error (dB)")
    axes[0].set_ylim(-1.6, 0.8)
    fig.tight_layout(w_pad=0.6, rect=(0, 0.1, 1, 1))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, handlelength=2.6)
    save(fig, f"spatial-{'upright-' if upright else ''}{tag}")


def residual_angle(tags=("1m", "2m")):
    """RMS level error of the summed two-way at each angle, horizontal orbit."""
    curves = (("raw", "uncorrected, 10° steps", "#c9c7bf", "-", 1.3),
              ("complex", "corrected, complex responses interpolated", RAW, (0, (1.2, 1.4)), 1.6),
              ("corrected", "corrected, 10° steps", COR, "-", 1.6),
              ("corrected2", "corrected, 2.5° steps", ALT, (0, (5, 2)), 1.4))
    fig, axes = plt.subplots(1, len(tags), figsize=(7.2, 2.9), sharey=True)
    for ax, tag in zip(axes, tags):
        for case, label, color, ls, lw in curves:
            angles = ANGLES[suffix(case)]
            rig = bem("rigorous" + suffix(case), tag)
            e = per_angle_rms(bem(case, tag), rig, F2, angles, twoway.BANDS["system"])
            a = np.concatenate([[-180], angles])
            e = np.concatenate([e[-1:], e])
            ax.plot(a, e, color=color, ls=ls, lw=lw, label=label)
        ax.set_xlim(-180, 180)
        ax.set_xticks(np.arange(-180, 181, 60))
        ax.set_xlabel("angle (°), negative on the woofer's side")
        ax.set_title(f"microphone at {tag[0]} m", loc="left")
        ax.set_ylim(0, 4)
    axes[0].set_ylabel("RMS error, 50 Hz–16 kHz (dB)")
    fig.tight_layout(w_pad=0.6, rect=(0, 0.17, 1, 1))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, handlelength=2.6)
    save(fig, "residual-angle")


def step_convergence():
    """The residual against the angular step, at 1, 2 and 4 m."""
    steps = (10.0, 5.0, 2.5)
    sfx = ("", "5", "2")
    colors = dict(zip(("1m", "2m", "4m"), ps.SLOTS))
    lines = (("bem", "corrected", "BEM, powers interpolated", "-", "o"),
             ("bem", "complex", "BEM, complex responses interpolated", (0, (1.2, 1.4)), "s"),
             ("ded", "corrected", "edge-dipole model, powers", (0, (5, 2)), "^"))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for tag in ("1m", "2m", "4m"):
        for model, kind, _, ls, marker in lines:
            rms, pir = [], []
            for s in sfx:
                angles = ANGLES[s]
                p, r = bem(kind + s, tag, model=model), bem("rigorous" + s, tag, model=model)
                mask = relevant(r, F2, angles, twoway.BANDS["system"])
                rms.append(float(np.sqrt(np.mean(level_error(p, r)[mask] ** 2))))
                c = cta_both(kind + s, "v-" + kind + s, tag, model)
                ref = cta_both("rigorous" + s, "v-rigorous" + s, tag, model)
                band = (F2 >= 100) & (F2 <= 16000)
                pir.append(float(np.max(np.abs(ps.db10(c["pir"] / ref["pir"]))[band])))
            axes[0].plot(steps, rms, color=colors[tag], ls=ls, marker=marker, ms=4)
            axes[1].plot(steps, pir, color=colors[tag], ls=ls, marker=marker, ms=4)
    for ax, title in zip(axes, ("summed two-way, RMS error (dB)", "PIR from both orbits, largest error (dB)")):
        ax.set_xscale("log")
        ax.set_xticks(steps)
        ax.set_xticklabels([f"{s:g}°" for s in steps])
        ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        ax.set_xlim(11.5, 2.2)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("angular step")
        ax.set_title(title, loc="left")
    handles = [Line2D([], [], color=colors[t], lw=1.8) for t in ("1m", "2m", "4m")]
    handles += [Line2D([], [], color=ps.INK_2, ls=ls, marker=m, ms=4) for _, _, _, ls, m in lines]
    labels = ["1 m", "2 m", "4 m"] + [label for _, _, label, _, _ in lines]
    fig.tight_layout(w_pad=1.0, rect=(0, 0.2, 1, 1))
    fig.legend(handles, labels, loc="lower center", ncol=3, handlelength=2.6)
    save(fig, "step-convergence")


def lateral_offset():
    """The lateral offset entered wrong, with the powers or the complex responses interpolated:
    the two-way (BEM) at 1 m, horizontal orbit, 10-degree steps."""
    rows = (("raw", None, "no correction"),
            ("corrected", "complex", "right lateral offset"),
            ("lateral_omitted", "lateral_omitted-complex", "lateral offset left at 0"),
            ("lateral_sign", "lateral_sign-complex", "lateral offset, wrong sign"))
    ref = horizontal_cta2034(at_ten(bem("rigorous"), ANGLES[""]))
    band = (F2 >= 100) & (F2 <= 16000)

    def numbers(case):
        p, r = bem(case), bem("rigorous")
        mask = relevant(r, F2, ANGLES[""], twoway.BANDS["system"])
        rms = float(np.sqrt(np.mean(level_error(p, r)[mask] ** 2)))
        c = horizontal_cta2034(at_ten(p, ANGLES[""]))
        return rms, float(np.max(np.abs(ps.db10(c["pir"] / ref["pir"]))[band]))

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=True)
    y = np.arange(len(rows))[::-1]
    for k, (which, color, label) in enumerate(((0, COR, "powers interpolated"), (1, RAW, "complex responses interpolated"))):
        values = [numbers(r[which] if r[which] else r[0]) for r in rows]
        for ax, j in zip(axes, (0, 1)):
            v = [x[j] for x in values]
            ax.barh(y + (0.19 if k == 0 else -0.19), v, height=0.36, color=color, label=label, edgecolor="white", lw=0.5)
            for yy, val in zip(y, v):
                ax.text(val + 0.05, yy + (0.19 if k == 0 else -0.19), f"{val:.2f}", va="center", fontsize=8,
                        color=ps.INK_2)
    for ax, title in zip(axes, ("summed two-way, RMS error (dB)", "PIR, largest error (dB)")):
        ax.set_title(title, loc="left")
        ax.grid(axis="y", visible=False)
        ax.set_xlim(0, ax.get_xlim()[1] * 1.12)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([r[2] for r in rows])
    fig.tight_layout(w_pad=0.8, rect=(0, 0.12, 1, 1))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2)
    save(fig, "lateral-offset")


def distance():
    """Results at 2 m: from the 1 m data, as measured or rebuilt, and from 2 m data (2.5-degree steps)."""
    ref = cta_both("rigorous2", "v-rigorous2", "2m")
    curves = (("1m", "rigorous2", "v-rigorous2", 4.0, "rigorous at 1 m, used as measured", RIG, (0, (1.2, 1.4))),
              ("1m", "corrected@2m2", "v-corrected@2m2", 4.0, "centered at 1 m, rebuilt at 2 m", RAW, "-"),
              ("2m", "corrected2", "v-corrected2", 1.0, "centered at 2 m, corrected", COR, (0, (5, 2))))
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    for tag, h, v, scale, label, color, ls in curves:
        c = cta_both(h, v, tag)
        ax.semilogx(F2, ps.db10(c["pir"] / (scale * ref["pir"])), color=color, ls=ls, lw=1.6, label=label)
    ax.axhline(0, color=ps.AXIS, lw=0.8)
    ps.freq_axis(ax, 50, 16000)
    ax.set_ylim(-1.2, 1.2)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("error (dB)")
    ax.set_title("PIR at 2 m from both orbits, against the rigorous 2 m data (2.5° steps)", loc="left")
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "distance")


def bem_vs_ded(tag="1m", freqs=(3000.0, 6350.0, 10000.0)):
    """The tweeter's horizontal directivity, BEM against the edge-dipole model (rigorous, 5-degree steps)."""
    angles = ANGLES["5"]
    i0 = angle_index(angles, 0)
    fig, axes = plt.subplots(1, len(freqs), figsize=(7.2, 2.7), sharey=True)
    for ax, fr in zip(axes, freqs):
        j = int(np.argmin(np.abs(F2 - fr)))
        for model, label, color, ls in (("bem", "boundary elements", RIG, "-"), ("ded", "edge-dipole model", RAW, (0, (5, 2)))):
            p = bem("rigorous5", tag, "tweeter", model)
            a = np.concatenate([[-180], angles])
            v = ps.db10(p[:, j] / p[i0, j])
            ax.plot(a, np.concatenate([v[-1:], v]), color=color, ls=ls, lw=1.5, label=label)
        ax.set_xlim(-180, 180)
        ax.set_xticks(np.arange(-180, 181, 90))
        ax.set_ylim(-30, 5)
        ax.set_title(f"{F2[j] / 1000:.2g} kHz", loc="left")
        ax.set_xlabel("angle (°)")
    axes[0].set_ylabel("tweeter, level re 0° (dB)")
    fig.tight_layout(w_pad=0.6, rect=(0, 0.1, 1, 1))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, handlelength=2.6)
    save(fig, "bem-vs-ded")


# ---------------------------------------------------------------------------------------------
# The three-way, edge-dipole model
# ---------------------------------------------------------------------------------------------

def three_way_maps(tag="1m"):
    """The summed three-way: rigorous polar map, and the level error before and after the correction."""
    from rotcenter.metrics import BANDS as B3
    ref = study(f"{tag}/rigorous")
    on = ref[list(A3).index(0)]
    a = np.concatenate([[-180], A3])
    full = lambda rows: np.concatenate([rows[-1:], rows], axis=0)  # noqa: E731
    mid = np.sqrt(F3[1:] * F3[:-1])
    fe = np.concatenate([[F3[0] ** 2 / mid[0]], mid, [F3[-1] ** 2 / mid[-1]]])
    ae = np.concatenate([[a[0] - 5], (a[1:] + a[:-1]) / 2, [a[-1] + 5]])
    mask = ~relevant(ref, F3, A3, B3["system"])
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0), sharey=True)
    m0 = axes[0].pcolormesh(fe, ae, full(ps.db10(ref / on[None, :])), cmap=ps.SEQUENTIAL, vmin=-30, vmax=3,
                            rasterized=True)
    axes[0].set_title("level re 0° (dB)", loc="left")
    for ax, case, title in ((axes[1], "raw", "uncorrected: error (dB)"), (axes[2], "corrected", "corrected: error (dB)")):
        err = level_error(study(f"{tag}/{case}"), ref)
        m1 = ax.pcolormesh(fe, ae, full(err), cmap=ps.DIVERGING, vmin=-3, vmax=3, rasterized=True)
        ax.pcolormesh(fe, ae, full(np.where(mask, 1.0, np.nan)), cmap=mpl.colors.ListedColormap(["white"]),
                      alpha=0.85, rasterized=True)
        ax.set_title(title, loc="left")
    for ax in axes:
        ps.freq_axis(ax, 50, 20000, sparse=True)
        ax.set_xticks([100, 1000, 10000])
        ax.set_xticklabels(["100", "1k", "10k"])
        ax.set_xlabel("frequency (Hz)")
        ax.grid(False, which="both")
        ax.set_ylim(-180, 180)
        ax.set_yticks(np.arange(-180, 181, 60))
    axes[0].set_ylabel("angle (°)")
    fig.tight_layout(w_pad=0.4)
    fig.colorbar(m0, ax=axes[0], location="bottom", pad=0.25, aspect=18)
    fig.colorbar(m1, ax=axes[1:], location="bottom", pad=0.25, aspect=36, shrink=0.7)
    save(fig, f"three-way-maps-{tag}")


def three_way_spatial():
    """PIR and sound power of the three-way, horizontal orbit, against the rigorous ones."""
    def cta(key):
        p = study(key)
        return horizontal_cta2034({int(x): p[i] for i, x in enumerate(A3)})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
    for ax, tag in zip(axes, ("1m", "2m")):
        ref = cta(f"{tag}/rigorous")
        for case, label, color, ls, lw in (("raw", "uncorrected", RAW, (0, (1.2, 1.4)), 1.6),
                                           ("corrected", "corrected", COR, "-", 1.6)):
            c = cta(f"{tag}/{case}")
            ax.semilogx(F3, ps.db10(c["pir"] / ref["pir"]), color=color, ls=ls, lw=lw, label=f"PIR, {label}")
            ax.semilogx(F3, ps.db10(c["sp"] / ref["sp"]), color=color, ls=ls, lw=0.9, alpha=0.6,
                        label=f"sound power, {label}")
        ax.axhline(0, color=ps.AXIS, lw=0.8)
        ps.freq_axis(ax, 50, 20000)
        ax.set_xticks([100, 1000, 10000])
        ax.set_xticklabels(["100", "1k", "10k"])
        ax.set_xlabel("frequency (Hz)")
        ax.set_title(f"microphone at {tag[0]} m", loc="left")
    axes[0].set_ylabel("error (dB)")
    axes[0].set_ylim(-1.6, 0.6)
    fig.tight_layout(w_pad=0.6, rect=(0, 0.17, 1, 1))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, handlelength=2.6)
    save(fig, "three-way-spatial")


def main():
    style()
    setups_3d("three")
    setups_3d("two")
    top_view()
    errors_geometry()
    responses("1m")
    spatial("1m")
    spatial("1m", upright=True)
    residual_angle()
    step_convergence()
    lateral_offset()
    distance()
    bem_vs_ded()
    three_way_maps()
    three_way_spatial()


if __name__ == "__main__":
    main()
