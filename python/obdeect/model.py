"""Optical models: build, store and trace."""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np

from obdeect import _core
from obdeect.compiler import compile_model, model_hash
from obdeect.formats import ModelError
from obdeect.repository import read_telescope

STATUS_NAMES: tuple[str, ...] = _core.STATUS_NAMES
DETECTED = STATUS_NAMES.index("detected")


class OpticalModel:
    """A compiled optical model, ready for tracing.

    The model is immutable; its SHA-256 identifies it in every result.
    """

    def __init__(self, document: dict[str, Any]):
        isolated_document = deepcopy(document)
        if isolated_document.get("sha256") != model_hash(isolated_document):
            raise ModelError("optical model content does not match its sha256")
        text = json.dumps(isolated_document, ensure_ascii=False, sort_keys=True)
        try:
            self._native = _core.Model(text)
        except ValueError as error:
            raise ModelError(str(error)) from error
        self._document = isolated_document

    @property
    def document(self) -> dict[str, Any]:
        """Return a copy of the compiled document."""
        return deepcopy(self._document)

    @classmethod
    def from_records(
        cls,
        records: Mapping[str, dict[str, Any]],
        read_file: Callable[[str], str],
        *,
        configuration_seed: int,
        zenith_deg: float | None = None,
        propagation_group_index: float = 1.0,
        provenance: dict[str, Any] | None = None,
    ) -> OpticalModel:
        """Compile simulation-models parameter records (see ``compiler.compile_model``)."""
        return cls(
            compile_model(
                records,
                read_file,
                configuration_seed=configuration_seed,
                zenith_deg=zenith_deg,
                propagation_group_index=propagation_group_index,
                provenance=provenance,
            )
        )

    @classmethod
    def from_simulation_models(
        cls,
        root: Path,
        model: str,
        version: str,
        *,
        configuration_seed: int,
        zenith_deg: float | None = None,
    ) -> OpticalModel:
        """Compile a telescope model from a simulation-models checkout."""
        telescope = read_telescope(Path(root), model, version)
        return cls.from_records(
            telescope.records,
            telescope.read_file,
            configuration_seed=configuration_seed,
            zenith_deg=zenith_deg,
            propagation_group_index=telescope.propagation_group_index,
            provenance=telescope.provenance,
        )

    @classmethod
    def load(cls, path: Path) -> OpticalModel:
        return cls(json.loads(Path(path).read_text(encoding="utf-8")))

    def save(self, path: Path) -> None:
        Path(path).write_text(
            json.dumps(self._document, ensure_ascii=False, indent=1) + "\n", encoding="utf-8"
        )

    @property
    def sha256(self) -> str:
        return self._document["sha256"]

    @property
    def has_structure(self) -> bool:
        """True if the model's own geometry shadows instead of ``telescope_transmission``."""
        return self._document["provenance"]["structure"]

    @property
    def sampling_area_m2(self) -> float:
        """Area of the disk star photons are sampled on (sim_telarray convention)."""
        return math.pi * self._document["sampling"]["disk_radius_m"] ** 2

    @property
    def focal_surface_z_m(self) -> float:
        """Height of the focal plane, or of the focal-surface vertex, in the telescope frame."""
        focal = self.document["optics"]["focal_surface"]
        return focal["asphere"]["vertex_z_m"] if "asphere" in focal else focal["plane_z_m"]

    def emitter_photons(
        self,
        emitted: float,
        *,
        seed: int,
        position_m: tuple[float, float, float],
        intensity: tuple[np.ndarray, np.ndarray],
        axis: tuple[float, float, float] = (0.0, 0.0, 1.0),
        pulse: tuple[np.ndarray, np.ndarray] | None = None,
        wavelength_nm: float = 400.0,
        disk: tuple[tuple[float, float, float], tuple[float, float, float], float] | None = None,
        first_id: int = 0,
    ) -> dict[str, np.ndarray]:
        """Photons of a point emitter (a calibration light source) that cross a disk.

        The emitter at `position_m` (telescope frame) sends `emitted` photons
        into 4 pi with the relative intensity per solid angle `intensity`
        (angles from `axis` [deg], values; linear in 1 - cos(angle), zero
        beyond the last angle) and emission times from the relative density
        `pulse` (times [ns], values; None: all at 0). Photons are made only
        on their way through `disk` (centre, unit normal, radius [m]; default
        the star-light disk). Their weights are in [0, 1], and their number is
        Poisson distributed, so keeping each photon with the probability of
        its weight gives the photons crossing the disk; the sum of weights
        estimates their mean number. Start positions are the emitter, times
        are emission times.
        """
        if disk is None:
            sampling = self.document["sampling"]
            disk = ((0.0, 0.0, sampling["disk_z_m"]), (0.0, 0.0, 1.0), sampling["disk_radius_m"])
        angle_deg, values = (np.asarray(column, dtype=float) for column in intensity)
        pulse_ns, pulse_values = (
            (np.empty(0), np.empty(0))
            if pulse is None
            else (np.asarray(column, dtype=float) for column in pulse)
        )
        centre, normal, radius = disk
        source = _core.Emitter(
            np.asarray(position_m, dtype=float),
            np.asarray(axis, dtype=float),
            1.0 - np.cos(np.radians(angle_deg)),
            values,
            pulse_ns,
            pulse_values,
            np.asarray(centre, dtype=float),
            np.asarray(normal, dtype=float),
            float(radius),
            wavelength_nm,
        )
        count = np.random.default_rng([seed, first_id]).poisson(
            emitted * source.candidates_per_photon()
        )
        return source.sample(first_id, int(count), seed)

    def star_photons(
        self,
        count: int,
        *,
        seed: int,
        towards_source: tuple[float, float, float] = (0.0, 0.0, 1.0),
        distance_m: float = math.inf,
        wavelength_nm: float = 400.0,
        first_id: int = 0,
    ) -> dict[str, np.ndarray]:
        """Photons from a star or a point source, as in sim_telarray's star light.

        `towards_source` is the unit vector to the source in the telescope frame
        (+z along the optical axis, towards the sky); `distance_m` is the
        source distance from the dish vertex (the centre of the sampled disk,
        as in sim_telarray). Times are zero for the wavefront at that point.
        """
        sampling = self._document["sampling"]
        direction = np.asarray(towards_source, dtype=float)
        return _core.sample_source(
            direction / np.linalg.norm(direction),
            distance_m,
            sampling["disk_radius_m"],
            sampling["disk_z_m"],
            sampling["start_z_m"],
            wavelength_nm,
            first_id,
            count,
            seed,
        )

    def trace(
        self,
        photons: Mapping[str, np.ndarray],
        *,
        seed: int,
        apply_responses: bool = True,
        threads: int = 0,
        direct: bool = False,
    ) -> dict[str, np.ndarray]:
        """Trace photons; returns per-photon arrays (see docs/USER_GUIDE.md).

        `photons` holds position_m (n, 3), direction (n, 3), wavelength_nm,
        time_ns, weight and photon_id. With apply_responses False only geometry
        acts (weights unchanged); `threads` 0 uses every core. `direct`
        photons start inside the telescope on their way to the focal surface
        (e.g. from a flasher): mirrors and telescope transmission do not act,
        the obscurers of that leg do.
        """

        def column(name: str, dtype: type) -> np.ndarray:
            return np.ascontiguousarray(photons[name], dtype=dtype)

        return self._native.trace(
            column("position_m", np.float64),
            column("direction", np.float64),
            column("wavelength_nm", np.float64),
            column("time_ns", np.float64),
            column("weight", np.float64),
            column("photon_id", np.uint64),
            seed=seed,
            apply_responses=apply_responses,
            threads=threads,
            direct=direct,
        )
