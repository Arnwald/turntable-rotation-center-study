"""Figures and LaTeX tables of the study, from results/study.npz.

Run ``scripts/run_study.py`` first.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from rotcenter import plotstyle as ps  # noqa: E402
from rotcenter.constants import C, P_REF  # noqa: E402
from rotcenter.correction import PlaneGeometry, PlaneOffsets  # noqa: E402
from rotcenter.cta2034 import SP_WEIGHTS, horizontal_cta2034  # noqa: E402
from rotcenter.directivity import Radiator  # noqa: E402
from rotcenter.drivers import (MID_RADIATOR, SIG225, build_drivers, closed_box_alignment,  # noqa: E402
                               load_datasheet)
from rotcenter.geometry import CENTERED, MID, TWEETER, WOOFER, correction_offsets  # noqa: E402
from rotcenter.metrics import BANDS, LEVEL_FLOOR_DB, level_error, per_angle_rms, relevant, summary  # noqa: E402
from rotcenter.signal import log_grid, rfft_freqs, smooth_power  # noqa: E402

FIG = ROOT / "paper" / "figures"
TAB = ROOT / "paper" / "tables"
DATA = np.load(ROOT / "results" / "study.npz")
F = DATA["freqs"]
ANG = DATA["angles"]
DRIVERS = ("woofer", "mid", "tweeter")
CASES = ("rigorous", "raw", "corrected")
SPL0 = 10.0 * np.log10(P_REF**2)


def spl(power):
    return ps.db10(power) - SPL0


def get(key, drv, what="power"):
    return DATA[f"{key}/{drv}/{what}"]


def full_circle(rows):
    """Rows for -170 ... 180 -> -180 ... 180 (the 180-degree row repeated at -180)."""
    return np.concatenate([rows[-1:], rows], axis=0)


ANG_FULL = np.concatenate([[-180], ANG])


def angle_axis(ax, label=True):
    ax.set_ylim(-180, 180)
    ax.set_yticks(np.arange(-180, 181, 60))
    if label:
        ax.set_ylabel("angle (°)")


def freq_edges(f):
    mid = np.sqrt(f[1:] * f[:-1])
    return np.concatenate([[f[0] ** 2 / mid[0]], mid, [f[-1] ** 2 / mid[-1]]])


def angle_edges(a):
    return np.concatenate([[a[0] - 5], (a[1:] + a[:-1]) / 2, [a[-1] + 5]])


def legend_cases(ax, cases=CASES, **kw):
    handles = [Line2D([], [], color=ps.CASE_COLORS[c], ls=ps.CASE_STYLES[c], lw=1.3) for c in cases]
    ax.legend(handles, [ps.CASE_LABELS[c] for c in cases], **kw)


# ---------------------------------------------------------------------------------------------
# Geometry of the measurement errors
# ---------------------------------------------------------------------------------------------

def fig_kinematics():
    angles = np.arange(-180, 181, 1.0)
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.4))
    rows = []
    for di, dist in enumerate((1.0, 2.0)):
        for name, pl, color in (("tweeter / midrange", TWEETER, ps.SLOTS[0]), ("woofer", WOOFER, ps.SLOTS[1])):
            g = PlaneGeometry(PlaneOffsets(**correction_offsets(CENTERED, pl, dist)))
            meas = np.array([g.measured(a) for a in angles])
            ideal = np.array([g.ideal(a) for a in angles])
            level = 20 * np.log10(ideal[:, 1] / meas[:, 1])
            delay = (meas[:, 1] - ideal[:, 1]) / C * 1e3
            direction = np.degrees(meas[:, 0] - ideal[:, 0])
            ls = "-" if dist == 1.0 else "--"
            lab = f"{name}, {dist:.0f} m"
            axes[0].plot(angles, level, ls=ls, color=color, label=lab)
            axes[1].plot(angles, delay, ls=ls, color=color)
            axes[2].plot(angles, direction, ls=ls, color=color)
            rows.append((name, dist, level.min(), level.max(), delay.min(), delay.max(),
                         direction.min(), direction.max()))
    titles = ("level error of the raw data (dB)", "arrival delay (ms)", "direction error (°)")
    for ax, t in zip(axes, titles):
        ax.set_xlim(-180, 180)
        ax.set_xticks(np.arange(-180, 181, 90))
        ax.set_xlabel("nominal angle (°)")
        ax.set_title(t, loc="left")
    fig.tight_layout(w_pad=1.2, rect=(0, 0.1, 1, 1))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=6.5, handlelength=2.4)
    fig.savefig(FIG / "kinematics.pdf")
    plt.close(fig)
    with open(TAB / "kinematics.tex", "w") as fh:
        for name, dist, lmin, lmax, tmin, tmax, amin, amax in rows:
            fh.write(f"{name} & {dist:.0f}\\,m & ${lmin:+.2f}$ to ${lmax:+.2f}$ & ${tmin:+.2f}$ to ${tmax:+.2f}$ "
                     f"& ${amin:+.1f}$ to ${amax:+.1f}$ \\\\\n")


# ---------------------------------------------------------------------------------------------
# Driver models
# ---------------------------------------------------------------------------------------------

def fig_driver_models():
    freqs = rfft_freqs()
    drivers = build_drivers(freqs)
    fd = log_grid(20, 20000, 48)
    datasheets = {"woofer": "dayton_sig225_4_spl", "mid": "hivi_dm7600_spl_0deg", "tweeter": "sica_lp9028_n92_spl"}
    titles = {"woofer": "woofer, Dayton Audio SIG225-4", "mid": "midrange, HiVi DM-7600",
              "tweeter": "tweeter, SICA LP 90.28/N92"}
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.5), sharey=True)
    for ax, name in zip(axes, DRIVERS):
        f, d = load_datasheet(datasheets[name])
        ax.semilogx(f, d, color=ps.MUTED, lw=2.2, alpha=0.55, label="datasheet, on axis")
        model = np.interp(np.log(fd), np.log(freqs[1:]), 20 * np.log10(np.abs(drivers[name].response[1:]) / P_REF))
        ax.semilogx(fd, model, color=ps.SLOTS[0], lw=1.1, label="model, on axis")
        if name == "mid":
            f30, d30 = load_datasheet("hivi_dm7600_spl_30deg")
            ax.semilogx(f30, d30, color=ps.SLOTS[1], lw=2.2, alpha=0.45, label="datasheet, 30°")
            rad = Radiator(**MID_RADIATOR)
            k = 2 * np.pi * freqs[1:] / C
            m30 = 20 * np.log10(np.abs(drivers[name].response[1:] * rad.directivity(k, np.radians(30))) / P_REF)
            ax.semilogx(fd, np.interp(np.log(fd), np.log(freqs[1:]), m30), color=ps.SLOTS[1], lw=1.1, ls="--",
                        label="model, 30°")
        if name == "woofer":
            fc, qtc = closed_box_alignment()
            ax.text(25, 101, f"24 L closed box:\n$f_c$ = {fc:.1f} Hz, $Q_{{tc}}$ = {qtc:.2f}", fontsize=6.3,
                    color=ps.INK_2, va="top")
        ps.freq_axis(ax, sparse=True)
        ax.set_ylim(55, 105)
        ax.set_title(titles[name], loc="left")
        ax.set_xlabel("frequency (Hz)")
        ax.legend(loc="lower right", fontsize=6.0)
    axes[0].set_ylabel("SPL at 1 m, 2.83 V, half space (dB)")
    fig.tight_layout(w_pad=0.6)
    fig.savefig(FIG / "driver_models.pdf")
    plt.close(fig)


def fig_inbox():
    from rotcenter.simulate import Loudspeaker
    ls = Loudspeaker()
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    for name in DRIVERS:
        drv = ls.drivers[name]
        own_axis = np.array([drv.placement.x, drv.placement.y, 1.0])
        inside = smooth_power(ls.freqs, np.abs(ls.pressures(own_axis)[name]) ** 2, 24, F)
        half = smooth_power(ls.freqs, np.abs(drv.response) ** 2, 24, F)
        color = ps.DRIVER_COLORS[name]
        axes[0].semilogx(F, spl(half), color=color, lw=0.8, alpha=0.45)
        axes[0].semilogx(F, spl(inside), color=color, lw=1.2, label=ps.DRIVER_LABELS[name])
        axes[1].semilogx(F, ps.db10(inside / half), color=color, lw=1.1, label=ps.DRIVER_LABELS[name])
    i0 = list(ANG).index(0)
    axes[1].semilogx(F, spl(get("1m/rigorous", "system")[i0]) - 86.0, color=ps.DRIVER_COLORS["system"], lw=1.1,
                     label="3-way sum on the reference axis (re 86 dB)")
    ps.freq_axis(axes[0])
    axes[0].set_ylim(55, 102)
    axes[0].set_ylabel("SPL at 1 m, 2.83 V (dB)")
    axes[0].set_title("on each driver's axis at 1 m", loc="left")
    axes[0].text(22, 99, "thick: in the cabinet\nthin: half space (datasheet)", fontsize=6.0, color=ps.INK_2, va="top")
    axes[0].legend(loc="lower right", fontsize=6.5)
    ps.freq_axis(axes[1])
    axes[1].set_ylim(-8, 4)
    axes[1].axhline(0, color=ps.AXIS, lw=0.6)
    axes[1].set_ylabel("dB")
    axes[1].set_title("cabinet effect (DED), same positions", loc="left")
    axes[1].legend(loc="lower right", fontsize=6.0)
    for ax in axes:
        ax.set_xlabel("frequency (Hz)")
    fig.tight_layout(w_pad=1.0)
    fig.savefig(FIG / "inbox_onaxis.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Polar maps
# ---------------------------------------------------------------------------------------------

def _map(ax, values, vmin, vmax, cmap, mask=None):
    mesh = ax.pcolormesh(freq_edges(F), angle_edges(ANG_FULL), values, cmap=cmap, vmin=vmin, vmax=vmax,
                         shading="flat", rasterized=True)
    if mask is not None:
        shade = np.where(mask, 1.0, np.nan)
        ax.pcolormesh(freq_edges(F), angle_edges(ANG_FULL), shade, cmap=plt.matplotlib.colors.ListedColormap(["white"]),
                      alpha=0.85, shading="flat", rasterized=True)
    ps.freq_axis(ax, 100, 20000, sparse=True)
    angle_axis(ax)
    ax.grid(False, which="both")
    return mesh


def fig_system_maps():
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.3), sharex=True, sharey=True)
    for row, tag in enumerate(("1m", "2m")):
        ref = get(f"{tag}/rigorous", "system")
        on = ref[list(ANG).index(0)]
        norm = full_circle(ps.db10(ref / on[None, :]))
        m0 = _map(axes[row, 0], norm, -30, 3, ps.SEQUENTIAL)
        cs = axes[row, 0].contour(F, ANG_FULL, norm, levels=[-20, -10, -6], colors=["white"], linewidths=0.5, linestyles="solid")
        axes[row, 0].clabel(cs, fmt="%d", fontsize=5.5, inline_spacing=1)
        mask = full_circle(~relevant(ref, F, ANG, BANDS["system"]))
        for col, case in ((1, "raw"), (2, "corrected")):
            err = full_circle(level_error(get(f"{tag}/{case}", "system"), ref))
            m1 = _map(axes[row, col], err, -3, 3, ps.DIVERGING, mask)
        axes[row, 0].set_title(f"rigorous, {tag[0]} m: level re on axis", loc="left")
        axes[row, 1].set_title(f"uncorrected − rigorous, {tag[0]} m", loc="left")
        axes[row, 2].set_title(f"corrected − rigorous, {tag[0]} m", loc="left")
        for ax in axes[row, 1:]:
            ax.set_ylabel("")
    for ax in axes[1]:
        ax.set_xlabel("frequency (Hz)")
    for ax in axes[0]:
        ax.set_xlabel("")
    fig.tight_layout(h_pad=0.8, w_pad=0.4)
    cb0 = fig.colorbar(m0, ax=axes[:, 0], location="bottom", shrink=0.9, pad=0.13, aspect=25)
    cb0.set_label("dB re on axis")
    cb1 = fig.colorbar(m1, ax=axes[:, 1:], location="bottom", shrink=0.6, pad=0.13, aspect=40)
    cb1.set_label("level error (dB); faded: more than 20 dB below the on-axis level")
    fig.savefig(FIG / "system_maps.pdf")
    plt.close(fig)


def fig_driver_error_maps(tag):
    fig, axes = plt.subplots(3, 2, figsize=(6.2, 5.6), sharex=True, sharey=True)
    for row, name in enumerate(DRIVERS):
        ref = get(f"{tag}/rigorous", name)
        mask = full_circle(~relevant(ref, F, ANG, BANDS[name]))
        for col, case in enumerate(("raw", "corrected")):
            err = full_circle(level_error(get(f"{tag}/{case}", name), ref))
            mesh = _map(axes[row, col], err, -3, 3, ps.DIVERGING, mask)
            lo, hi = BANDS[name]
            for edge in (lo, hi):
                if 100 < edge < 20000:
                    axes[row, col].axvline(edge, color=ps.INK_2, lw=0.5, ls=(0, (2, 2)))
            what = "uncorrected" if case == "raw" else "corrected"
            axes[row, col].set_title(f"{ps.DRIVER_LABELS[name]}, {what} − rigorous", loc="left")
            if col:
                axes[row, col].set_ylabel("")
    for ax in axes[-1]:
        ax.set_xlabel("frequency (Hz)")
    fig.tight_layout(h_pad=0.6, w_pad=0.4)
    cb = fig.colorbar(mesh, ax=axes, location="bottom", shrink=0.7, pad=0.08, aspect=40)
    cb.set_label(f"level error at {tag[0]} m (dB); dotted lines: evaluation band; faded: below −20 dB re on axis")
    fig.savefig(FIG / f"driver_error_maps_{tag}.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Responses at selected angles
# ---------------------------------------------------------------------------------------------

def fig_selected(tag, name="system"):
    sel = (30, 60, 90, 120, 150, 180, -60, -120)
    fig, axes = plt.subplots(2, 4, figsize=(7.0, 3.6), sharex=True, sharey=True)
    for ax, a in zip(axes.ravel(), sel):
        i = list(ANG).index(a)
        for case in CASES:
            p = smooth_power(F, get(f"{tag}/{case}", name)[i], 12, F)
            ax.semilogx(F, spl(p), color=ps.CASE_COLORS[case], ls=ps.CASE_STYLES[case], lw=1.1)
        ps.freq_axis(ax, 100, 20000, sparse=True)
        ax.set_title(("180°" if a == 180 else f"{a:+d}°").replace("-", "−"), loc="left")
    lo = 50 if name == "system" else 40
    for ax in axes[:, 0]:
        ax.set_ylabel("SPL (dB)")
    for ax in axes[1]:
        ax.set_xlabel("frequency (Hz)")
    top = spl(get(f"{tag}/rigorous", name)).max()
    axes[0, 0].set_ylim(top - 42, top + 4)
    fig.tight_layout(h_pad=0.6, w_pad=0.3)
    legend_cases(axes[0, 0], loc="lower left", fontsize=6.0, handlelength=2.4)
    fig.savefig(FIG / f"selected_{name}_{tag}.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# What is left after the correction: diffraction, angular step
# ---------------------------------------------------------------------------------------------

def fig_residuals():
    fig, axes = plt.subplots(2, 4, figsize=(7.0, 3.7), sharex=True, sharey="row")
    variants = (("corrected", "with diffraction, 10° steps", ps.SLOTS[0], "-"),
                ("step5", "with diffraction, 5° steps", ps.SLOTS[1], "--"),
                ("corrected@direct", "direct sound only, 10° steps", ps.SLOTS[2], ":"))
    for row, tag in enumerate(("1m", "2m")):
        for col, name in enumerate(DRIVERS + ("system",)):
            ax = axes[row, col]
            for key, label, color, ls in variants:
                if key.endswith("@direct"):
                    ref = get(f"{tag}-direct/rigorous", name)
                    pw = get(f"{tag}-direct/corrected", name)
                else:
                    ref = get(f"{tag}/rigorous", name)
                    pw = get(f"{tag}/{key}", name)
                e = per_angle_rms(pw, ref, F, ANG, BANDS[name])
                ax.plot(ANG_FULL, full_circle(e), color=color, ls=ls, lw=1.1, marker="o", ms=2.2, label=label)
            ax.set_xlim(-180, 180)
            ax.set_xticks(np.arange(-180, 181, 90))
            ax.set_title(f"{ps.DRIVER_LABELS[name]}, {tag[0]} m", loc="left")
            if row == 1:
                ax.set_xlabel("angle (°)")
        axes[row, 0].set_ylabel("RMS level error (dB)")
    axes[0, 0].set_ylim(0, 3.2)
    axes[1, 0].set_ylim(0, 3.2)
    fig.tight_layout(h_pad=0.6, w_pad=0.3)
    axes[0, 0].legend(loc="upper center", fontsize=6.0, handlelength=2.2)
    fig.savefig(FIG / "residuals.pdf")
    plt.close(fig)


def fig_impulses():
    fs = 48000.0
    names = {1: "90°", 2: "150°", 3: "180°"}
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 3.6), sharex="col")
    for row, name in enumerate(("tweeter", "mid")):
        irs = {c: DATA[f"1m/{c}/{name}/ir"] for c in CASES}
        for col, k in enumerate((1, 2, 3)):
            ax = axes[row, col]
            t = np.arange(irs["rigorous"].shape[1]) / fs * 1e3
            peak = np.argmax(np.abs(irs["rigorous"][k]))
            for case in CASES:
                ax.plot(t, irs[case][k] * 1e3, color=ps.CASE_COLORS[case], ls=ps.CASE_STYLES[case], lw=0.9)
            t0 = t[peak]
            ax.set_xlim(t0 - 0.35, t0 + 1.65)
            ax.set_title(f"{ps.DRIVER_LABELS[name]}, {names[k]}, 1 m", loc="left")
            if row == 1:
                ax.set_xlabel("time after emission (ms)")
        axes[row, 0].set_ylabel("pressure (mPa)")
    fig.tight_layout(h_pad=0.6, w_pad=0.6)
    legend_cases(axes[0, 0], loc="upper right", fontsize=5.8, handlelength=2.2)
    fig.savefig(FIG / "impulses.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# CTA-2034-A curves and Predicted In-Room Response
# ---------------------------------------------------------------------------------------------

def cta(key, name="system", angles=ANG):
    p = get(key, name)
    return horizontal_cta2034({int(a): p[i] for i, a in enumerate(angles)})


def fig_pir():
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.4), sharex=True, gridspec_kw=dict(height_ratios=[1.35, 1]))
    for col, tag in enumerate(("1m", "2m")):
        ref = cta(f"{tag}/rigorous")
        ax = axes[0, col]
        for case in CASES:
            c = cta(f"{tag}/{case}")
            for curve, lw, off in (("pir", 1.2, 0.0), ("sp", 0.9, 0.0), ("lw", 0.9, 0.0)):
                ax.semilogx(F, spl(c[curve]) + 20 * np.log10(float(tag[0])), color=ps.CASE_COLORS[case],
                            ls=ps.CASE_STYLES[case], lw=lw)
        for curve, y in (("lw", None), ("pir", None), ("sp", None)):
            v = spl(ref[curve]) + 20 * np.log10(float(tag[0]))
            i = np.argmin(np.abs(F - 650))
            ax.text(650, v[i] + 1.2, {"lw": "listening window", "pir": "PIR", "sp": "sound power"}[curve],
                    fontsize=6.2, color=ps.INK, ha="center")
        ps.freq_axis(ax, 50, 20000)
        ax.set_ylim(62, 92)
        ax.set_title(f"microphone at {tag[0]} m (levels referred to 1 m)", loc="left")
        ax = axes[1, col]
        for case, color, ls in (("raw", ps.CASE_COLORS["raw"], ":"), ("corrected", ps.CASE_COLORS["corrected"], "--")):
            c = cta(f"{tag}/{case}")
            ax.semilogx(F, ps.db10(c["pir"] / ref["pir"]), color=color, ls=ls, lw=1.2,
                        label=f"PIR, {ps.CASE_LABELS[case]}")
            ax.semilogx(F, ps.db10(c["sp"] / ref["sp"]), color=color, ls=ls, lw=0.7, alpha=0.7,
                        label=f"sound power, {ps.CASE_LABELS[case]}")
        ax.axhline(0, color=ps.AXIS, lw=0.6)
        ps.freq_axis(ax, 50, 20000)
        ax.set_ylim(-1.6, 0.6)
        ax.set_xlabel("frequency (Hz)")
        ax.set_title("difference with the rigorous case (dB)", loc="left")
    axes[0, 0].set_ylabel("SPL (dB)")
    axes[1, 0].set_ylabel("dB")
    legend_cases(axes[0, 0], loc="lower left", fontsize=6.0)
    axes[1, 1].legend(loc="lower left", fontsize=5.8, ncol=1)
    fig.tight_layout(h_pad=0.6, w_pad=0.8)
    fig.savefig(FIG / "pir.pdf")
    plt.close(fig)


VARIANT_LABELS = {
    "raw": "no correction",
    "corrected": "true parameters",
    "step5": "true parameters, 5° steps",
    "lateral_omitted": "lateral offset omitted (0)",
    "lateral_sign": "lateral offset, wrong sign",
    "woofer_offset_omitted": "woofer offset omitted (0)",
    "depth_plus10": "depth +10 mm",
    "depth_minus10": "depth −10 mm",
    "lateral_plus10": "lateral +10 mm",
    "lateral_minus10": "lateral −10 mm",
    "distance_plus20": "mic distance +20 mm",
}


def variant_numbers():
    out = {}
    band = (F >= 100) & (F <= 20000)
    for tag in ("1m", "2m"):
        ref_cta = cta(f"{tag}/rigorous")
        for v in VARIANT_LABELS:
            row = {}
            for name in DRIVERS + ("system",):
                row[name] = summary(get(f"{tag}/{v}", name), get(f"{tag}/rigorous", name), F, ANG, BANDS[name])
            c = cta(f"{tag}/{v}")
            row["pir"] = float(np.max(np.abs(ps.db10(c["pir"] / ref_cta["pir"]))[band]))
            row["sp"] = float(np.max(np.abs(ps.db10(c["sp"] / ref_cta["sp"]))[band]))
            row["di"] = float(np.max(np.abs(c["di"] - ref_cta["di"])[band]))
            out[(tag, v)] = row
    return out


def fig_sensitivity(numbers):
    order = ["raw", "corrected", "step5", "lateral_omitted", "lateral_sign", "woofer_offset_omitted",
             "depth_plus10", "depth_minus10", "lateral_plus10", "lateral_minus10", "distance_plus20"]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), sharey=True)
    y = np.arange(len(order))[::-1]
    for ax, metric, title in ((axes[0], "sys_rms", "RMS level error of the 3-way sum (dB)"),
                              (axes[1], "pir", "largest PIR error, 100 Hz–20 kHz (dB)")):
        for k, (tag, color) in enumerate((("1m", ps.SLOTS[0]), ("2m", ps.SLOTS[1]))):
            vals = [numbers[(tag, v)]["system"]["rms"] if metric == "sys_rms" else numbers[(tag, v)]["pir"] for v in order]
            ax.barh(y + (0.2 if k == 0 else -0.2), vals, height=0.38, color=color, label=f"{tag[0]} m",
                    edgecolor="white", linewidth=0.6)
            for yy, val in zip(y, vals):
                ax.text(val + 0.04, yy + (0.2 if k == 0 else -0.2), f"{val:.2f}", va="center", fontsize=5.6, color=ps.INK_2)
        ax.set_title(title, loc="left")
        ax.grid(axis="y", visible=False)
        ax.set_xlim(0, max(4.6, ax.get_xlim()[1]))
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([VARIANT_LABELS[v] for v in order], fontsize=6.5)
    axes[1].legend(loc="lower right", fontsize=6.5)
    fig.tight_layout(w_pad=0.6)
    fig.savefig(FIG / "sensitivity.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------------------------

def tables(numbers):
    TAB.mkdir(parents=True, exist_ok=True)
    # Error summary: raw / corrected, with and without diffraction.
    with open(TAB / "errors.tex", "w") as fh:
        for tag in ("1m", "2m"):
            for dtag, dlabel in (("", "DED"), ("-direct", "direct only")):
                for name in DRIVERS + ("system",):
                    ref = get(f"{tag}{dtag}/rigorous", name)
                    cells = []
                    for case in ("raw", "corrected"):
                        s = summary(get(f"{tag}{dtag}/{case}", name), ref, F, ANG, BANDS[name])
                        cells.append(f"{s['rms']:.2f} & {s['p95']:.2f} & {s['max']:.1f}")
                    fh.write(f"{tag[0]}\\,m & {dlabel} & {ps.DRIVER_LABELS[name]} & " + " & ".join(cells) + " \\\\\n")
                fh.write("\\addlinespace\n")
    # Variants.
    with open(TAB / "variants.tex", "w") as fh:
        for v, label in VARIANT_LABELS.items():
            cells = []
            for tag in ("1m", "2m"):
                r = numbers[(tag, v)]
                cells.append(f"{r['woofer']['rms']:.2f} & {r['mid']['rms']:.2f} & {r['tweeter']['rms']:.2f} & "
                             f"{r['system']['rms']:.2f} & {r['pir']:.2f}")
            fh.write(f"{label} & " + " & ".join(cells) + " \\\\\n")
    # CTA-2034-A curves: largest deviations.
    band = (F >= 100) & (F <= 20000)
    with open(TAB / "cta.tex", "w") as fh:
        for tag in ("1m", "2m"):
            for dtag, dlabel in (("", "DED"), ("-direct", "direct only")):
                ref = cta(f"{tag}{dtag}/rigorous")
                for case in ("raw", "corrected"):
                    c = cta(f"{tag}{dtag}/{case}")
                    cells = []
                    for k in ("lw", "er", "sp", "pir"):
                        cells.append(f"{np.max(np.abs(ps.db10(c[k] / ref[k]))[band]):.2f}")
                    cells.append(f"{np.max(np.abs(c['di'] - ref['di'])[band]):.2f}")
                    what = "uncorrected" if case == "raw" else "corrected"
                    fh.write(f"{tag[0]}\\,m & {dlabel} & {what} & " + " & ".join(cells) + " \\\\\n")
    # Sound power weights.
    total = sum(SP_WEIGHTS[abs(int(a))] for a in ANG)
    with open(TAB / "weights.tex", "w") as fh:
        cells = [f"{a}° & {SP_WEIGHTS[a] / total:.4f}" for a in range(0, 100, 10)]
        for i in range(0, 10, 5):
            fh.write(" & ".join(cells[i:i + 5]) + " \\\\\n")
    # Key numbers for the text.
    fc, qtc = closed_box_alignment()
    with open(TAB / "numbers.tex", "w") as fh:
        fh.write(f"\\newcommand{{\\fcbox}}{{{fc:.1f}}}\n\\newcommand{{\\qtcbox}}{{{qtc:.2f}}}\n")
        for tag in ("1m", "2m"):
            t = "One" if tag == "1m" else "Two"
            for v, key in (("raw", "Raw"), ("corrected", "Cor"), ("lateral_omitted", "LatZero"), ("lateral_sign", "LatSign"),
                           ("step5", "Five")):
                r = numbers[(tag, v)]
                fh.write(f"\\newcommand{{\\pir{key}{t}}}{{{r['pir']:.2f}}}\n")
                fh.write(f"\\newcommand{{\\sys{key}{t}}}{{{r['system']['rms']:.2f}}}\n")
                fh.write(f"\\newcommand{{\\sysmax{key}{t}}}{{{r['system']['max']:.1f}}}\n")
                fh.write(f"\\newcommand{{\\tw{key}{t}}}{{{r['tweeter']['rms']:.2f}}}\n")
                fh.write(f"\\newcommand{{\\twmax{key}{t}}}{{{r['tweeter']['max']:.1f}}}\n")
                fh.write(f"\\newcommand{{\\wf{key}{t}}}{{{r['woofer']['rms']:.2f}}}\n")
                fh.write(f"\\newcommand{{\\wfmax{key}{t}}}{{{r['woofer']['max']:.1f}}}\n")


def main():
    ps.apply()
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    fig_kinematics()
    fig_driver_models()
    fig_inbox()
    fig_system_maps()
    for tag in ("1m", "2m"):
        fig_driver_error_maps(tag)
        fig_selected(tag)
    fig_selected("1m", "tweeter")
    fig_residuals()
    fig_impulses()
    fig_pir()
    numbers = variant_numbers()
    fig_sensitivity(numbers)
    tables(numbers)


if __name__ == "__main__":
    main()
