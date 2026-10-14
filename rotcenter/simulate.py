"""Free-field pressure of each driver at given microphone positions."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from .constants import C, FS, N_FFT
from .diffraction import EdgeDiffraction, direct_term, edge_propagation, front_edges
from .drivers import build_drivers
from .signal import rfft_freqs

CACHE = Path(__file__).resolve().parent.parent / "results" / "cache"


class Loudspeaker:
    """The three drivers in the cabinet, with or without edge diffraction."""

    def __init__(self, n_fft: int = N_FFT, fs: float = FS, edge_step: float = 0.002):
        self.n_fft = n_fft
        self.fs = fs
        self.freqs = rfft_freqs(n_fft, fs)
        self.k = 2.0 * np.pi * self.freqs / C
        self.drivers = build_drivers(self.freqs)
        self.edge_points = front_edges(step=edge_step)[0]
        self.edges = {name: self._edges(d, edge_step) for name, d in self.drivers.items()}

    def _edges(self, driver, step):
        """Edge strengths of one driver, cached on disk (they take a while to compute)."""
        r = driver.radiator
        key = (f"v1-{driver.name}-{r.radius}-{r.height}-{driver.placement}-{step}-{self.n_fft}-{self.fs}-"
               f"{np.abs(driver.response).sum():.9e}")
        path = CACHE / f"edges-{driver.name}-{hashlib.sha1(key.encode()).hexdigest()[:10]}.npz"
        if path.exists():
            data = np.load(path)
            edges = EdgeDiffraction.__new__(EdgeDiffraction)
            edges.points, edges.lengths = data["points"], data["lengths"]
            edges.normalization = float(data["normalization"])
            edges.strength = data["strength"]
            return edges
        edges = EdgeDiffraction(driver, self.freqs, step=step)
        CACHE.mkdir(parents=True, exist_ok=True)
        np.savez(path, points=edges.points, lengths=edges.lengths,
                 normalization=edges.normalization, strength=edges.strength)
        return edges

    def pressures(self, mic, diffraction: bool = True) -> dict[str, np.ndarray]:
        """Complex pressure spectrum of every driver at ``mic`` (speaker frame, m)."""
        mic = np.asarray(mic, float)
        prop = edge_propagation(self.edge_points, mic, self.k) if diffraction else None
        out = {}
        for name, driver in self.drivers.items():
            p = direct_term(driver, mic, self.k)
            if diffraction:
                p = p + self.edges[name].edge_term(prop)
            out[name] = driver.response * p
        return out

    def polar(self, mics, diffraction: bool = True) -> dict[str, np.ndarray]:
        """Spectra at several positions: ``{driver: array (n_positions, n_freqs)}``."""
        rows = [self.pressures(m, diffraction) for m in np.asarray(mics, float)]
        return {name: np.stack([r[name] for r in rows]) for name in self.drivers}

    def impulse_responses(self, spectra: np.ndarray) -> np.ndarray:
        """Real impulse responses (last axis = time) from rfft spectra."""
        return np.fft.irfft(spectra, self.n_fft, axis=-1)
