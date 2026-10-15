"""Figures and LaTeX tables of the boundary-element check, from results/bem.npz.

Run ``scripts/run_bem.py`` first.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from rotcenter import plotstyle as ps  # noqa: E402
from rotcenter import twoway  # noqa: E402
from rotcenter.constants import C, P_REF, R_REF  # noqa: E402
from rotcenter.cta2034 import cta2034, horizontal_cta2034  # noqa: E402
from rotcenter.geometry import CABINET  # noqa: E402
from rotcenter.metrics import level_error, per_angle_rms, relevant, summary  # noqa: E402
from rotcenter.signal import smooth_power  # noqa: E402

FIG = ROOT / "paper" / "figures"
TAB = ROOT / "paper" / "tables"
DATA = np.load(sys.argv[1] if len(sys.argv) > 1 else ROOT / "results" / "bem.npz")
F = DATA["freqs"]
ANG = DATA["angles"]
ANG5 = DATA["angles5"]
DRIVERS = ("woofer", "tweeter")
ALL = DRIVERS + ("system",)
LABELS = {"woofer": "woofer", "tweeter": "tweeter", "system": "2-way sum"}
MODELS = {"bem": "BEM", "ded": "DED"}
CASES = ("rigorous", "raw", "corrected")
BANDS = twoway.BANDS
SPL0 = 10.0 * np.log10(P_REF**2)
#: Band of the CTA-2034-A deviations (Hz).
CTA_BAND = (F >= 100) & (F <= 16000)


def spl(power):
    return ps.db10(power) - SPL0


def get(model, tag, case, drv):
    return DATA[f"{model}/{tag}/{case}/{drv}/power"]


def angles_of(case):
    return ANG5 if case.endswith("5") else ANG


def full_circle(rows, angles):
    """Rows for ``-180 + step ... 180`` -> ``-180 ... 180`` (the 180-degree row repeated)."""
    return np.concatenate([[-180], angles]), np.concatenate([rows[-1:], rows], axis=0)


def freq_edges(f):
    mid = np.sqrt(f[1:] * f[:-1])
    return np.concatenate([[f[0] ** 2 / mid[0]], mid, [f[-1] ** 2 / mid[-1]]])


def angle_edges(a):
    step = a[1] - a[0]
    return np.concatenate([[a[0] - step / 2], (a[1:] + a[:-1]) / 2, [a[-1] + step / 2]])


def angle_axis(ax, label=True):
    ax.set_ylim(-180, 180)
    ax.set_yticks(np.arange(-180, 181, 60))
    if label:
        ax.set_ylabel("angle (°)")


def draw_map(ax, values, angles, vmin, vmax, cmap, mask=None, lo=100):
    a, values = full_circle(values, angles)
    mesh = ax.pcolormesh(freq_edges(F), angle_edges(a), values, cmap=cmap, vmin=vmin, vmax=vmax, shading="flat",
                         rasterized=True)
    if mask is not None:
        _, mask = full_circle(mask, angles)
        ax.pcolormesh(freq_edges(F), angle_edges(a), np.where(mask, 1.0, np.nan),
                      cmap=plt.matplotlib.colors.ListedColormap(["white"]), alpha=0.85, shading="flat",
                      rasterized=True)
    ps.freq_axis(ax, lo, 16000, sparse=True)
    angle_axis(ax)
    ax.grid(False, which="both")
    return mesh


def band_lines(ax, name):
    for edge in BANDS[name]:
        if 100 < edge < 16000:
            ax.axvline(edge, color=ps.INK_2, lw=0.5, ls=(0, (2, 2)))


def legend_cases(ax, cases=CASES, **kw):
    handles = [Line2D([], [], color=ps.CASE_COLORS[c], ls=ps.CASE_STYLES[c], lw=1.3) for c in cases]
    ax.legend(handles, [ps.CASE_LABELS[c] for c in cases], **kw)


# ---------------------------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------------------------

def fig_model():
    mm = 1e3
    W, H = CABINET.width * mm, CABINET.height * mm
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw=dict(width_ratios=[1.45, 1]))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((-W / 2, -H / 2), W, H, fc="#e9e7e1", ec="#3d3c39", lw=1.0))
    wx, tx = twoway.WOOFER.x * mm, twoway.TWEETER.x * mm
    ax.add_patch(Circle((wx, 0), 83.8, fc="#2f2f2d", ec="#3d3c39", lw=0.5))
    ax.add_patch(Circle((wx, 0), 30.0, fc="#555451", ec="#3d3c39", lw=0.5))
    ax.add_patch(Circle((tx, 0), 14.0, fc="#2f2f2d", ec="#3d3c39", lw=0.5))
    ax.plot([-W / 2 - 10, W / 2 + 10], [0, 0], color=ps.SLOTS[0], lw=0.8, ls=(0, (6, 2, 1.5, 2)))
    ax.text(W / 2 + 12, 0, "plane of the\nhorizontal orbit\n(symmetry plane)", fontsize=6.0, color=ps.SLOTS[0],
            va="center", ha="left")
    ax.plot(tx, 0, marker="+", ms=9, mew=1.4, color="#e34948")
    ax.annotate("reference point:\ntweeter center", (tx, 0), xytext=(tx + 10, 92), fontsize=6.2, color=ps.INK,
                ha="center", va="center", arrowprops=dict(arrowstyle="-", lw=0.6, color=ps.INK_2))
    ax.text(wx, -H / 2 + 22, "SIG225-4: rigid recessed cone", fontsize=6.2, color=ps.INK, ha="center", va="center")
    ax.text(tx, -H / 2 + 22, "LP 90.28/N92:\nrigid dome", fontsize=6.2, color=ps.INK, ha="center", va="center")
    y = -H / 2 - 14
    for x0, x1, t in ((-W / 2, wx, "135"), (wx, tx, "205"), (tx, W / 2, "100")):
        ax.annotate("", (x0, y), (x1, y), arrowprops=dict(arrowstyle="<|-|>", lw=0.6, color=ps.INK_2,
                                                           mutation_scale=6, shrinkA=0, shrinkB=0))
        ax.text((x0 + x1) / 2, y - 12, t, fontsize=6.2, color=ps.INK_2, ha="center", va="center")
    ax.text(0, y - 34, "front 440 × 290 mm, depth 290 mm; dimensions in mm", fontsize=6.2, color=ps.INK_2,
            ha="center", va="center")
    ax.set_xlim(-W / 2 - 12, W / 2 + 95)
    ax.set_ylim(-H / 2 - 56, H / 2 + 8)
    ax.set_title("front of the two-way", loc="left")

    rad = twoway.radiators()
    r = np.linspace(0, 1, 400)
    for name, color, x0, text in (("woofer", ps.SLOTS[0], 0.0, "cone, $S_d$ = 220.8 cm$^2$"),
                                  ("tweeter", ps.SLOTS[2], 125.0, "dome, 28 mm,\n7.5 mm high")):
        rr = r * rad[name].radius
        z = rad[name].surface(rr) * mm
        bx.plot(np.concatenate([x0 - rr[::-1] * mm, x0 + rr * mm]), np.concatenate([z[::-1], z]), color=color, lw=1.4)
        bx.text(x0, -48 if name == "woofer" else -6, text, fontsize=6.0, color=color, ha="center", va="top")
    bx.axhline(0, color=ps.MUTED, lw=0.8)
    bx.text(-100, 1.5, "baffle", fontsize=6.0, color=ps.INK_2, va="bottom")
    bx.set_aspect("equal")
    bx.set_xlim(-102, 150)
    bx.set_ylim(-62, 22)
    bx.set_xlabel("radius (mm)")
    bx.set_ylabel("z (mm)")
    bx.set_title("radiating surfaces (sections, to scale)", loc="left")
    fig.tight_layout(w_pad=1.0)
    fig.savefig(FIG / "bem_model.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# BEM against DED
# ---------------------------------------------------------------------------------------------

def half_space(name, distance):
    """``|H D(theta) r_ref / r|^2`` at the 0-degree position: the driver in an infinite baffle."""
    drv = twoway.drivers(F)[name]
    mic = np.array([twoway.REFERENCE.x, 0.0, distance])
    d = mic - np.array([drv.placement.x, drv.placement.y, 0.0])
    r = np.linalg.norm(d)
    k = 2 * np.pi * F / C
    p = drv.response * drv.radiator.directivity(k, np.arccos(d[2] / r)) * R_REF / r
    return smooth_power(F, np.abs(p) ** 2, 24, F)


def fig_onaxis():
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6), sharey=True)
    i0 = list(ANG5).index(0)
    for ax, tag in zip(axes, ("1m", "2m")):
        for name in DRIVERS:
            ref = half_space(name, float(tag[0]))
            for model, ls in (("bem", "-"), ("ded", "--")):
                p = get(model, tag, "rigorous5", name)[i0]
                ax.semilogx(F, ps.db10(p / ref), color=ps.DRIVER_COLORS[name], ls=ls, lw=1.1 if model == "bem" else 0.9,
                            label=f"{name}, {MODELS[model]}")
        ax.axhline(0, color=ps.AXIS, lw=0.6)
        ax.axhline(-6, color=ps.AXIS, lw=0.6, ls=":")
        ps.freq_axis(ax, 20, 16000)
        ax.set_ylim(-10, 5)
        ax.set_xlabel("frequency (Hz)")
        ax.set_title(f"cabinet effect at 0°, {tag[0]} m (re half space)", loc="left")
    axes[0].set_ylabel("dB")
    axes[1].legend(loc="lower right", fontsize=6.0, ncol=2)
    fig.tight_layout(w_pad=0.8)
    fig.savefig(FIG / "bem_onaxis.pdf")
    plt.close(fig)


def fig_maps(tag="1m"):
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.4), sharex=True, sharey=True)
    for row, name in enumerate(DRIVERS):
        ref = get("bem", tag, "rigorous5", name)
        on = ref[list(ANG5).index(0)]
        norm = ps.db10(ref / on[None, :])
        m0 = draw_map(axes[row, 0], norm, ANG5, -30, 3, ps.SEQUENTIAL, lo=50)
        a, full = full_circle(norm, ANG5)
        cs = axes[row, 0].contour(F, a, full, levels=[-20, -10, -6], colors=["white"], linewidths=0.5)
        axes[row, 0].clabel(cs, fmt="%d", fontsize=5.5, inline_spacing=1)
        ded = get("ded", tag, "rigorous5", name)
        diff = ps.db10(ded / ded[list(ANG5).index(0)][None, :]) - norm
        mask = ~relevant(ref, F, ANG5, BANDS[name])
        m1 = draw_map(axes[row, 1], diff, ANG5, -6, 6, ps.DIVERGING, mask, lo=50)
        band_lines(axes[row, 1], name)
        axes[row, 0].set_title(f"{name}, BEM, {tag[0]} m: level re 0°", loc="left")
        axes[row, 1].set_title(f"{name}: DED − BEM, levels re 0°", loc="left")
        axes[row, 1].set_ylabel("")
    for ax in axes[1]:
        ax.set_xlabel("frequency (Hz)")
    fig.tight_layout(h_pad=0.8, w_pad=0.4)
    cb0 = fig.colorbar(m0, ax=axes[:, 0], location="bottom", shrink=0.9, pad=0.12, aspect=25)
    cb0.set_label("dB re 0°")
    cb1 = fig.colorbar(m1, ax=axes[:, 1], location="bottom", shrink=0.9, pad=0.12, aspect=25)
    cb1.set_label("level difference (dB); faded: outside the band or below −20 dB re 0°")
    fig.savefig(FIG / f"bem_vs_ded_maps_{tag}.pdf")
    plt.close(fig)


POLAR_FREQS = (3000.0, 6350.0, 10000.0)


def fig_polars(tag="1m", name="tweeter"):
    """Level re 0° against angle at a few frequencies: BEM against DED (top), and what the
    correction makes of the BEM data (bottom)."""
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.2), sharex=True, sharey=True)
    i0 = list(ANG5).index(0)
    for col, fr in enumerate(POLAR_FREQS):
        j = int(np.argmin(np.abs(F - fr)))

        def curve(model, case, angles):
            p = get(model, tag, case, name)
            ref0 = get(model, tag, "rigorous5", name)[i0, j]
            return full_circle(ps.db10(p[:, j] / ref0), angles)

        ax = axes[0, col]
        for model, color, ls, lw in (("bem", ps.SLOTS[0], "-", 1.3), ("ded", ps.SLOTS[1], "--", 1.1)):
            a, v = curve(model, "rigorous5", ANG5)
            ax.plot(a, v, color=color, ls=ls, lw=lw, label=f"{MODELS[model]}, rigorous")
        ax.set_title(f"{name}, {F[j] / 1000:.2g} kHz, {tag[0]} m", loc="left")
        ax = axes[1, col]
        a, v = curve("bem", "rigorous5", ANG5)
        ax.plot(a, v, color=ps.CASE_COLORS["rigorous"], lw=1.3, label="BEM, rigorous")
        a, v = curve("bem", "corrected", ANG)
        ax.plot(a, v, color=ps.CASE_COLORS["corrected"], ls="--", lw=1.0, marker="o", ms=2.0,
                label="BEM, corrected, 10° steps")
        a, v = curve("bem", "corrected5", ANG5)
        ax.plot(a, v, color=ps.SLOTS[3], ls=":", lw=1.1, label="BEM, corrected, 5° steps")
    for ax in axes.ravel():
        ax.set_xlim(-180, 180)
        ax.set_xticks(np.arange(-180, 181, 60))
        ax.set_ylim(-30, 6)
    for ax in axes[1]:
        ax.set_xlabel("angle (°)")
    for ax in axes[:, 0]:
        ax.set_ylabel("level re 0° (dB)")
    fig.tight_layout(h_pad=0.6, w_pad=0.3)
    axes[0, 0].legend(loc="lower center", fontsize=6.0)
    axes[1, 0].legend(loc="lower center", fontsize=6.0)
    fig.savefig(FIG / f"bem_polars_{tag}.pdf")
    plt.close(fig)


def ded_vs_bem_numbers():
    """RMS of the DED - BEM level difference over the relevant points, front and rear halves."""
    out = {}
    for tag in ("1m", "2m"):
        for name in ALL:
            ref = get("bem", tag, "rigorous5", name)
            ded = get("ded", tag, "rigorous5", name)
            i0 = list(ANG5).index(0)
            err = level_error(ded / ded[i0][None, :], ref / ref[i0][None, :])
            mask = relevant(ref, F, ANG5, BANDS[name])
            front = np.abs(ANG5)[:, None] <= 90
            row = {}
            for half, sel in (("front", front), ("rear", ~front)):
                e = err[mask & sel]
                row[half] = (float(np.sqrt(np.mean(e**2))), float(np.percentile(np.abs(e), 95)))
            out[(tag, name)] = row
    return out


# ---------------------------------------------------------------------------------------------
# The correction, on the BEM data
# ---------------------------------------------------------------------------------------------

def fig_error_maps(tag):
    fig, axes = plt.subplots(3, 2, figsize=(6.2, 5.6), sharex=True, sharey=True)
    for row, name in enumerate(ALL):
        ref = get("bem", tag, "rigorous", name)
        mask = ~relevant(ref, F, ANG, BANDS[name])
        for col, case in enumerate(("raw", "corrected")):
            err = level_error(get("bem", tag, case, name), ref)
            mesh = draw_map(axes[row, col], err, ANG, -3, 3, ps.DIVERGING, mask)
            band_lines(axes[row, col], name)
            what = "uncorrected" if case == "raw" else "corrected"
            axes[row, col].set_title(f"{LABELS[name]}, {what} − rigorous", loc="left")
            if col:
                axes[row, col].set_ylabel("")
    for ax in axes[-1]:
        ax.set_xlabel("frequency (Hz)")
    fig.tight_layout(h_pad=0.6, w_pad=0.4)
    cb = fig.colorbar(mesh, ax=axes, location="bottom", shrink=0.7, pad=0.08, aspect=40)
    cb.set_label(f"level error at {tag[0]} m, BEM (dB); dotted lines: evaluation band; faded: below −20 dB re 0°")
    fig.savefig(FIG / f"bem_error_maps_{tag}.pdf")
    plt.close(fig)


def fig_residuals():
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 3.7), sharex=True, sharey="row")
    curves = (("bem", "corrected", "BEM, 10° steps", ps.SLOTS[0], "-"),
              ("bem", "levels", "BEM, 10° steps, levels interpolated", ps.SLOTS[3], "-."),
              ("bem", "corrected5", "BEM, 5° steps", ps.SLOTS[1], "--"),
              ("ded", "corrected", "DED, 10° steps", ps.SLOTS[2], ":"))
    for row, tag in enumerate(("1m", "2m")):
        for col, name in enumerate(ALL):
            ax = axes[row, col]
            for model, case, label, color, ls in curves:
                rig = "rigorous5" if case.endswith("5") else "rigorous"
                angles = angles_of(case)
                e = per_angle_rms(get(model, tag, case, name), get(model, tag, rig, name), F, angles, BANDS[name])
                a, e = full_circle(e, angles)
                ax.plot(a, e, color=color, ls=ls, lw=1.1, marker="o", ms=1.8 if case.endswith("5") else 2.2,
                        label=label)
            ax.set_xlim(-180, 180)
            ax.set_xticks(np.arange(-180, 181, 90))
            ax.set_title(f"{LABELS[name]}, {tag[0]} m", loc="left")
            if row == 1:
                ax.set_xlabel("angle (°)")
        axes[row, 0].set_ylabel("RMS level error (dB)")
    fig.tight_layout(h_pad=0.6, w_pad=0.3)
    axes[0, 0].legend(loc="upper center", fontsize=6.0, handlelength=2.2)
    fig.savefig(FIG / "bem_residuals.pdf")
    plt.close(fig)


def cta_h(model, tag, case, name="system"):
    p = get(model, tag, case, name)
    return horizontal_cta2034({int(a): p[i] for i, a in enumerate(ANG)})


def cta_hv(model, tag, h_case, v_case, name="system"):
    ph, pv = get(model, tag, h_case, name), get(model, tag, v_case, name)
    return cta2034({int(a): ph[i] for i, a in enumerate(ANG)}, {int(a): pv[i] for i, a in enumerate(ANG)})


def deviation(c, ref, key):
    if key in ("di", "erdi"):
        return float(np.max(np.abs(c[key] - ref[key])[CTA_BAND]))
    return float(np.max(np.abs(ps.db10(c[key] / ref[key]))[CTA_BAND]))


def fig_pir():
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.4), sharex=True, gridspec_kw=dict(height_ratios=[1.35, 1]))
    for col, tag in enumerate(("1m", "2m")):
        gain = 20 * np.log10(float(tag[0]))
        ref = cta_h("bem", tag, "rigorous")
        ax = axes[0, col]
        for case in CASES:
            c = cta_h("bem", tag, case)
            for curve, lw in (("pir", 1.2), ("sp", 0.9), ("lw", 0.9)):
                ax.semilogx(F, spl(c[curve]) + gain, color=ps.CASE_COLORS[case], ls=ps.CASE_STYLES[case], lw=lw)
        for curve, at, dy in (("lw", 300, 1.2), ("pir", 1200, 1.2), ("sp", 2500, -3.0)):
            v = spl(ref[curve]) + gain
            i = np.argmin(np.abs(F - at))
            ax.text(at, v[i] + dy, {"lw": "listening window", "pir": "PIR", "sp": "sound power"}[curve],
                    fontsize=6.2, color=ps.INK, ha="center")
        ps.freq_axis(ax, 50, 16000)
        ax.set_ylim(62, 92)
        ax.set_title(f"BEM, microphone at {tag[0]} m (levels referred to 1 m)", loc="left")
        ax = axes[1, col]
        for case, ls in (("raw", ":"), ("corrected", "--")):
            c = cta_h("bem", tag, case)
            color = ps.CASE_COLORS[case]
            ax.semilogx(F, ps.db10(c["pir"] / ref["pir"]), color=color, ls=ls, lw=1.2, label=f"PIR, {ps.CASE_LABELS[case]}")
            ax.semilogx(F, ps.db10(c["sp"] / ref["sp"]), color=color, ls=ls, lw=0.7, alpha=0.7,
                        label=f"sound power, {ps.CASE_LABELS[case]}")
        ax.axhline(0, color=ps.AXIS, lw=0.6)
        ps.freq_axis(ax, 50, 16000)
        ax.set_ylim(-1.6, 0.8)
        ax.set_xlabel("frequency (Hz)")
        ax.set_title("difference with the rigorous case (dB)", loc="left")
    axes[0, 0].set_ylabel("SPL (dB)")
    axes[1, 0].set_ylabel("dB")
    legend_cases(axes[0, 0], loc="lower left", fontsize=6.0)
    axes[1, 1].legend(loc="lower left", fontsize=5.8)
    fig.tight_layout(h_pad=0.6, w_pad=0.8)
    fig.savefig(FIG / "bem_pir.pdf")
    plt.close(fig)


def fig_spinorama():
    """The two-way's CTA-2034-A curves from both BEM orbits, 1 m."""
    c = cta_hv("bem", "1m", "rigorous", "v-rigorous")
    fig, ax = plt.subplots(figsize=(7.0, 3.0))
    for key, label, color, lw in (("on_axis", "on axis", ps.INK, 1.2), ("lw", "listening window", ps.SLOTS[0], 1.1),
                                  ("er", "early reflections", ps.SLOTS[1], 1.1), ("sp", "sound power", ps.SLOTS[2], 1.1),
                                  ("pir", "PIR", ps.SLOTS[3], 1.3)):
        ax.semilogx(F, spl(c[key]), color=color, lw=lw, label=label)
    ax2 = ax.twinx()
    ax2.semilogx(F, c["di"], color=ps.SLOTS[6], lw=0.9, ls="--", label="sound power DI")
    ax2.semilogx(F, c["erdi"], color=ps.SLOTS[4], lw=0.9, ls="--", label="early reflections DI")
    ax2.set_ylim(-5, 25)
    ax2.set_ylabel("DI (dB)")
    ax2.grid(False, which="both")
    ps.freq_axis(ax, 20, 16000)
    ax.set_ylim(55, 95)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("SPL at 1 m (dB)")
    ax.set_title("two-way, BEM, both orbits at 1 m: CTA-2034-A curves", loc="left")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="lower center", ncol=4, fontsize=6.0)
    fig.tight_layout()
    fig.savefig(FIG / "bem_spinorama.pdf")
    plt.close(fig)


VARIANT_LABELS = {
    "raw": "no correction",
    "corrected": "true parameters",
    "corrected5": "true parameters, 5° steps",
    "levels": "levels interpolated",
    "levels5": "levels interpolated, 5° steps",
    "lateral_omitted": "lateral offset omitted (0)",
    "lateral_sign": "lateral offset, wrong sign",
    "woofer_offset_omitted": "woofer offset omitted (0)",
    "depth_plus10": "depth +10 mm",
    "depth_minus10": "depth −10 mm",
    "lateral_plus10": "lateral +10 mm",
    "lateral_minus10": "lateral −10 mm",
    "distance_plus20": "mic distance +20 mm",
}


def error_summary(model, tag, case, name):
    rig = "rigorous5" if case.endswith("5") else "rigorous"
    return summary(get(model, tag, case, name), get(model, tag, rig, name), F, angles_of(case), BANDS[name])


def variant_numbers(model="bem"):
    out = {}
    for tag in ("1m", "2m"):
        ref = cta_h(model, tag, "rigorous")
        for v in VARIANT_LABELS:
            row = {name: error_summary(model, tag, v, name) for name in ALL}
            if v.endswith("5"):
                # The 5-degree correction, read at the 10-degree angles for the CTA-2034-A curves.
                keep = np.isin(ANG5, ANG)
                p = get(model, tag, v, "system")[keep]
                c = horizontal_cta2034({int(a): p[i] for i, a in enumerate(ANG5[keep])})
            else:
                c = cta_h(model, tag, v)
            row["pir"] = deviation(c, ref, "pir")
            out[(tag, v)] = row
    return out


def tables(numbers, ded_numbers, dvb):
    TAB.mkdir(parents=True, exist_ok=True)
    with open(TAB / "bem_errors.tex", "w") as fh:
        for tag in ("1m", "2m"):
            for model in ("bem", "ded"):
                for name in ALL:
                    cells = []
                    for case in ("raw", "corrected"):
                        s = error_summary(model, tag, case, name)
                        cells.append(f"{s['rms']:.2f} & {s['p95']:.2f} & {s['max']:.1f}")
                    for case in ("levels", "corrected5"):
                        s = error_summary(model, tag, case, name)
                        cells.append(f"{s['rms']:.2f} & {s['max']:.1f}")
                    fh.write(f"{tag[0]}\\,m & {MODELS[model]} & {LABELS[name]} & " + " & ".join(cells) + " \\\\\n")
                fh.write("\\addlinespace\n")
    with open(TAB / "bem_variants.tex", "w") as fh:
        for v, label in VARIANT_LABELS.items():
            cells = []
            for tag in ("1m", "2m"):
                r = numbers[(tag, v)]
                cells.append(f"{r['woofer']['rms']:.2f} & {r['tweeter']['rms']:.2f} & {r['system']['rms']:.2f} & "
                             f"{r['pir']:.2f}")
            fh.write(f"{label} & " + " & ".join(cells) + " \\\\\n")
    with open(TAB / "bem_cta.tex", "w") as fh:
        for tag in ("1m", "2m"):
            for model in ("bem", "ded"):
                ref = cta_h(model, tag, "rigorous")
                for case in ("raw", "corrected"):
                    c = cta_h(model, tag, case)
                    cells = [f"{deviation(c, ref, k):.2f}" for k in ("lw", "er", "sp", "pir", "di")]
                    what = "uncorrected" if case == "raw" else "corrected"
                    fh.write(f"{tag[0]}\\,m & {MODELS[model]} & {what} & " + " & ".join(cells) + " \\\\\n")
    # Both orbits: BEM (vertical orbit rigorous) and DED (vertical orbit measured on its side too).
    rows = (("bem", "raw", "v-rigorous", "uncorrected", "rigorous"),
            ("bem", "corrected", "v-rigorous", "corrected", "rigorous"),
            ("ded", "raw", "v-raw", "uncorrected", "uncorrected"),
            ("ded", "corrected", "v-corrected", "corrected", "corrected"))
    with open(TAB / "bem_cta2.tex", "w") as fh:
        for tag in ("1m", "2m"):
            for model, hc, vc, hl, vl in rows:
                ref = cta_hv(model, tag, "rigorous", "v-rigorous")
                c = cta_hv(model, tag, hc, vc)
                cells = [f"{deviation(c, ref, k):.2f}" for k in ("lw", "er", "sp", "pir", "di", "erdi")]
                fh.write(f"{tag[0]}\\,m & {MODELS[model]} & {hl} & {vl} & " + " & ".join(cells) + " \\\\\n")
            fh.write("\\addlinespace\n")
    with open(TAB / "bem_dedvsbem.tex", "w") as fh:
        for tag in ("1m", "2m"):
            for name in ALL:
                r = dvb[(tag, name)]
                fh.write(f"{tag[0]}\\,m & {LABELS[name]} & {r['front'][0]:.2f} & {r['front'][1]:.2f} & "
                         f"{r['rear'][0]:.2f} & {r['rear'][1]:.2f} \\\\\n")

    # Key numbers for the text.
    macros = {}
    for model, nums in (("bem", numbers), ("ded", ded_numbers)):
        for tag in ("1m", "2m"):
            t = "One" if tag == "1m" else "Two"
            for v, key in (("raw", "Raw"), ("corrected", "Cor"), ("corrected5", "Five"), ("levels", "Lev"),
                           ("levels5", "LevFive"), ("lateral_omitted", "LatZero"), ("lateral_sign", "LatSign"),
                           ("woofer_offset_omitted", "WfZero")):
                r = nums[(tag, v)]
                macros[f"{model}Pir{key}{t}"] = f"{r['pir']:.2f}"
                for name, short in (("system", "Sys"), ("tweeter", "Tw"), ("woofer", "Wf")):
                    macros[f"{model}{short}{key}{t}"] = f"{r[name]['rms']:.2f}"
                    macros[f"{model}{short}max{key}{t}"] = f"{r[name]['max']:.1f}"
    for tag in ("1m", "2m"):
        t = "One" if tag == "1m" else "Two"
        for name, short in (("system", "Sys"), ("tweeter", "Tw"), ("woofer", "Wf")):
            macros[f"dvb{short}Front{t}"] = f"{dvb[(tag, name)]['front'][0]:.2f}"
            macros[f"dvb{short}Rear{t}"] = f"{dvb[(tag, name)]['rear'][0]:.2f}"
        ref = cta_hv("bem", tag, "rigorous", "v-rigorous")
        macros[f"bemPirBothCor{t}"] = f"{deviation(cta_hv('bem', tag, 'corrected', 'v-rigorous'), ref, 'pir'):.2f}"
        macros[f"bemPirBothRaw{t}"] = f"{deviation(cta_hv('bem', tag, 'raw', 'v-rigorous'), ref, 'pir'):.2f}"
        ref = cta_hv("ded", tag, "rigorous", "v-rigorous")
        macros[f"dedPirBothCor{t}"] = f"{deviation(cta_hv('ded', tag, 'corrected', 'v-corrected'), ref, 'pir'):.2f}"
        macros[f"dedPirBothRaw{t}"] = f"{deviation(cta_hv('ded', tag, 'raw', 'v-raw'), ref, 'pir'):.2f}"
    solve = {k[len("solve/"):]: DATA[k] for k in DATA.files if k.startswith("solve/")}
    macros["bemFreqs"] = f"{len(F)}"
    if "unknowns" in solve:
        macros["bemUnknownsMax"] = f"{int(solve['unknowns'].max()):,}".replace(",", "\\,")
        macros["bemTrianglesMax"] = f"{int(solve['triangles'].max()):,}".replace(",", "\\,")
        macros["bemGigabytesMax"] = f"{solve['bytes'].max() / 1e9:.1f}"
    if "seconds" in solve:
        macros["bemHours"] = f"{solve['seconds'].sum() / 3600:.1f}"
        macros["bemSecondsTop"] = f"{solve['seconds'][np.argmax(F)]:.0f}"
    with open(TAB / "bem_numbers.tex", "w") as fh:
        for name, value in macros.items():
            fh.write(f"\\newcommand{{\\{name}}}{{{value}}}\n")


def main():
    ps.apply()
    FIG.mkdir(parents=True, exist_ok=True)
    fig_model()
    fig_onaxis()
    for tag in ("1m", "2m"):
        fig_maps(tag)
        fig_error_maps(tag)
    fig_polars("1m")
    fig_residuals()
    fig_pir()
    fig_spinorama()
    tables(variant_numbers("bem"), variant_numbers("ded"), ded_vs_bem_numbers())


if __name__ == "__main__":
    main()
