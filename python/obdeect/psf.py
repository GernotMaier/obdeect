"""Point-spread function and effective area from traced photons."""

from __future__ import annotations

from typing import Any

import numpy as np

from obdeect.model import DETECTED, STATUS_NAMES


def containment(
    x: np.ndarray, y: np.ndarray, weights: np.ndarray | None = None, fraction: float = 0.8
) -> dict:
    """Centroid and the radius around it containing `fraction` of the weight."""
    if len(x) == 0:
        raise ValueError("no photons on the focal surface")
    w = np.ones_like(x) if weights is None else np.asarray(weights, dtype=float)
    total = w.sum()
    if not total > 0:  # Every detected photon has zero weight: use the geometry.
        w, total = np.ones_like(x), float(len(x))
    cx, cy = np.dot(w, x) / total, np.dot(w, y) / total
    r = np.hypot(x - cx, y - cy)
    order = np.argsort(r)
    cumulative = np.cumsum(w[order])
    index = min(int(np.searchsorted(cumulative, fraction * total)), len(r) - 1)
    return {
        "centroid_x_m": float(cx),
        "centroid_y_m": float(cy),
        "radius_m": float(r[order][index]),
    }


def summary(
    arrivals: dict[str, np.ndarray], sampling_area_m2: float, fraction: float = 0.8
) -> dict[str, Any]:
    """Status counts, effective area and focal-plane PSF of one traced source.

    The effective area is the detected fraction of the sampled photons times the
    sampled area, weighted by the traced optical weights.
    """
    status = arrivals["status"]
    detected = status == DETECTED
    counts = {name: int(np.count_nonzero(status == i)) for i, name in enumerate(STATUS_NAMES)}
    result: dict[str, Any] = {
        "photons": int(len(status)),
        "status": {k: v for k, v in counts.items() if v},
        "effective_area_m2": float(arrivals["weight"].sum() / len(status) * sampling_area_m2),
    }
    if detected.any():
        position = arrivals["position_m"][detected]
        psf = containment(position[:, 0], position[:, 1], arrivals["weight"][detected], fraction)
        result.update(psf)
        result[f"d{round(100 * fraction)}_m"] = 2 * psf.pop("radius_m")
    return result
