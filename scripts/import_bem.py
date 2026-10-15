"""Convert the boundary-element solver's report into ``data/bem/two_way.npz``.

Usage: ``python3 scripts/import_bem.py REPORT.json [OUT.npz]``

The report (written by Huygens' ``huygens-bench turntable --points data/bem/points.json --out``)
holds the frequencies, the microphone positions, one record per solved frequency (mesh size, time)
and ``pressure[source][frequency][point] = [re, im]``: Pa for 1 m^3/s of volume velocity out of
each driver, time factor ``exp(+j omega t)``. The sets of microphone positions are those of
``data/bem/points.json``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent


def main():
    report_path = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "data" / "bem" / "two_way.npz"
    report = json.loads(report_path.read_text())
    layout = json.loads((ROOT / "data" / "bem" / "points.json").read_text())

    points = np.array(report["points"], float)
    if not np.allclose(points, np.array(layout["points"], float), atol=1e-12):
        raise SystemExit("the report was not computed on data/bem/points.json")
    freqs = np.array(report["frequencies"], float)
    pressure = np.array(report["pressure"], float)
    pressure = pressure[..., 0] + 1j * pressure[..., 1]
    if pressure.shape != (len(report["sources"]), len(freqs), len(points)):
        raise SystemExit(f"unexpected pressure array {pressure.shape}")
    if not np.all(np.isfinite(pressure)):
        raise SystemExit("non-finite pressures")

    angles = np.array(layout["sets"][0]["angles_deg"], float)
    for s in layout["sets"]:
        if not np.array_equal(np.array(s["angles_deg"], float), angles):
            raise SystemExit("every set must have the same angles")
    solves = sorted(report["solves"], key=lambda s: s["frequency"])
    if not np.allclose([s["frequency"] for s in solves], freqs):
        raise SystemExit("the solve records do not match the frequencies")

    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        freqs=freqs,
        points=points,
        sources=np.array(report["sources"]),
        pressure=pressure,
        angles=angles,
        set_names=np.array([s["name"] for s in layout["sets"]]),
        set_first=np.array([s["first"] for s in layout["sets"]]),
        set_distance=np.array([s["distance"] for s in layout["sets"]], float),
        solve_triangles=np.array([s["triangles"] for s in solves]),
        solve_unknowns=np.array([s["unknowns"] for s in solves]),
        solve_seconds=np.array([s["seconds"] for s in solves], float),
        solve_bytes=np.array([s["bytes"] for s in solves], float),
    )
    print(f"{out}: {len(freqs)} frequencies, {len(points)} points, {len(report['sources'])} sources, "
          f"{out.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
