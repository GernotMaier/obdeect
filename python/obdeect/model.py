"""Optical models: build, store and trace."""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Iterator, Mapping, Sequence
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
ThreadPool = _core.ThreadPool


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
        focal = self._document["optics"]["focal_surface"]
        return focal["asphere"]["vertex_z_m"] if "asphere" in focal else focal["plane_z_m"]

    @property
    def propagation_group_index(self) -> float:
        """Group index of air the photon times use."""
        return self._document["optics"]["propagation_group_index"]

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
        source = self._emitter(position_m, intensity, axis, pulse, wavelength_nm, disk)
        return source.sample(first_id, _emitter_count(source, emitted, seed, first_id), seed)

    def trace_emitter(self, emitted: float, **options: Any) -> dict[str, np.ndarray]:
        """Make and trace the photons of a point emitter (see `trace_emitter_blocks`).

        Returns the photons of all blocks reaching the focal surface with a
        positive weight.
        """
        blocks = list(self.trace_emitter_blocks(emitted, **options))
        return {name: np.concatenate([block[name] for block in blocks]) for name in blocks[0]}

    def trace_emitter_blocks(
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
        direct: bool = False,
        threads: int = 0,
        block_size: int = 250_000,
    ) -> Iterator[dict[str, np.ndarray]]:
        """Make and trace the photons of a point emitter, `block_size` photons at a time.

        The photons are those of `emitter_photons`, traced as by `trace`;
        memory is bounded by one block. Yields, for each block (at least
        one), the photons reaching the focal surface with a positive weight:
        position_m, direction, time_ns and weight of the trace, and
        wavelength_nm. The photons do not depend on `block_size`.
        """
        if block_size < 1:
            raise ValueError("block_size must be positive")
        source = self._emitter(position_m, intensity, axis, pulse, wavelength_nm, disk)
        count = _emitter_count(source, emitted, seed, first_id)
        for start in range(0, count, block_size) or [0]:
            photons = source.sample(first_id + start, min(block_size, count - start), seed)
            traced = self.trace(photons, seed=seed, direct=direct, threads=threads)
            kept = (traced["status"] == DETECTED) & (traced["weight"] > 0.0)
            yield {name: traced[name][kept] for name in ("position_m", "direction", "time_ns")} | {
                "weight": traced["weight"][kept],
                "wavelength_nm": photons["wavelength_nm"][kept],
            }

    def _emitter(self, position_m, intensity, axis, pulse, wavelength_nm, disk) -> _core.Emitter:
        if disk is None:
            sampling = self._document["sampling"]
            disk = ((0.0, 0.0, sampling["disk_z_m"]), (0.0, 0.0, 1.0), sampling["disk_radius_m"])
        angle_deg, values = (np.asarray(column, dtype=float) for column in intensity)
        pulse_ns, pulse_values = (
            (np.empty(0), np.empty(0))
            if pulse is None
            else (np.asarray(column, dtype=float) for column in pulse)
        )
        centre, normal, radius = disk
        return _core.Emitter(
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
        outputs: Sequence[str] | None = None,
        detected_only: bool = False,
        classify_losses: bool = True,
        pool: ThreadPool | None = None,
    ) -> dict[str, np.ndarray]:
        """Trace photons; returns per-photon arrays (see docs/USER_GUIDE.md).

        `photons` holds position_m (n, 3), direction (n, 3), wavelength_nm,
        time_ns, weight and photon_id. With apply_responses False only geometry
        acts (weights unchanged); `threads` 0 uses every core. `direct`
        photons start inside the telescope on their way to the focal surface
        (e.g. from a flasher): mirrors and telescope transmission do not act,
        the obscurers of that leg do.

        `outputs` lists the columns to return (None: all of them); asking only
        for what is used saves memory traffic. With `detected_only` just the
        photons with status ``detected`` are returned, in input order, and the
        extra column ``index`` (their row in the input) becomes available.
        An unknown column raises ``ValueError``.

        With `classify_losses` False a photon that misses the primary mirror
        is reported as ``missed_primary`` even where it would be ``blocked``
        before missing it; detected photons are the same either way.

        With a `pool` (a `ThreadPool`) its threads run the work and `threads`
        is ignored; the results are the same.
        """
        arrays = _photon_arrays(photons)
        return self._native.trace(
            arrays["position_m"],
            arrays["direction"],
            arrays["wavelength_nm"],
            arrays["time_ns"],
            arrays["weight"],
            arrays["photon_id"],
            seed=seed,
            apply_responses=apply_responses,
            threads=threads,
            direct=direct,
            outputs=None if outputs is None else list(outputs),
            detected_only=detected_only,
            classify_losses=classify_losses,
            pool=pool,
        )

    def trace_bunches(
        self,
        bunches: Mapping[str, np.ndarray],
        *,
        seed: int,
        frame_rotation: np.ndarray,
        frame_offset_m: Sequence[float],
        length_unit_m: float,
        plane_z: float,
        acceptance: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
        wavelength_density: tuple[np.ndarray, np.ndarray] | None = None,
        roulette: bool = True,
        classify_losses: bool = True,
        pool: ThreadPool | None = None,
        threads: int = 0,
    ) -> dict[str, np.ndarray]:
        """Trace photon bunches of an air shower (see docs/USER_GUIDE.md).

        `bunches` holds x, y, cx, cy, time_ns, zem, size, wavelength_nm (all
        float32; wavelength 0 means "sample one from `wavelength_density`")
        and bunch_id (uint64, unique per bunch for one `seed`). `frame_rotation`
        (3 x 3, input to telescope frame), `frame_offset_m`, `length_unit_m`
        (metres per input length unit) and `plane_z` (input units) place the
        input plane in the telescope frame. `acceptance` is
        ``(wavelength_nm, zem_m, values)`` with values (n_wavelength, n_zem) in
        [0, 1]; `wavelength_density` is ``(nm, values)``.

        Each photon of a bunch is thinned by the acceptance and the model's
        response bound, traced, and, with `roulette`, kept with the probability
        weight / bound (survivors have weight 1; without it they carry that
        ratio as weight). Returns the surviving photons as columns position_m,
        direction, time_ns, weight, wavelength_nm, bunch_index (input row) and
        sub_photon, in input order. Invalid input raises ``ValueError`` before
        any tracing. The result does not depend on `threads` or `pool`.
        """
        dtypes = {
            "x": np.float32,
            "y": np.float32,
            "cx": np.float32,
            "cy": np.float32,
            "time_ns": np.float32,
            "zem": np.float32,
            "size": np.float32,
            "wavelength_nm": np.float32,
            "bunch_id": np.uint64,
        }
        arrays = {
            name: np.ascontiguousarray(bunches[name], dtype=dtype) for name, dtype in dtypes.items()
        }
        rotation = np.ascontiguousarray(frame_rotation, dtype=np.float64)
        if rotation.shape != (3, 3):
            raise ValueError("frame_rotation must be a 3 x 3 matrix")
        offset = np.ascontiguousarray(frame_offset_m, dtype=np.float64)
        if acceptance is not None:
            wavelengths, heights, values = (np.asarray(a, dtype=np.float64) for a in acceptance)
            if wavelengths.ndim != 1 or heights.ndim != 1:
                raise ValueError("acceptance axes must be one-dimensional")
            if values.shape != (wavelengths.size, heights.size):
                raise ValueError("acceptance values must have shape (n_wavelength, n_zem)")
            acceptance = (
                np.ascontiguousarray(wavelengths),
                np.ascontiguousarray(heights),
                np.ascontiguousarray(values).ravel(),
            )
        if wavelength_density is not None:
            nodes, weights = (np.ascontiguousarray(a, dtype=np.float64) for a in wavelength_density)
            if nodes.ndim != 1 or nodes.shape != weights.shape:
                raise ValueError("wavelength_density needs two equally long 1-D arrays")
            wavelength_density = (nodes, weights)
        return self._native.trace_bunches(
            arrays["x"],
            arrays["y"],
            arrays["cx"],
            arrays["cy"],
            arrays["time_ns"],
            arrays["zem"],
            arrays["size"],
            arrays["wavelength_nm"],
            arrays["bunch_id"],
            rotation.ravel(),
            offset,
            float(length_unit_m),
            float(plane_z),
            seed=seed,
            acceptance=acceptance,
            wavelength_density=wavelength_density,
            roulette=roulette,
            classify_losses=classify_losses,
            pool=pool if pool is not None else ThreadPool(threads),
        )


def _photon_arrays(photons: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """The photon columns as the contiguous arrays of the native tracer."""
    dtypes = {
        "position_m": np.float64,
        "direction": np.float64,
        "wavelength_nm": np.float64,
        "time_ns": np.float64,
        "weight": np.float64,
        "photon_id": np.uint64,
    }
    return {
        name: np.ascontiguousarray(photons[name], dtype=dtype) for name, dtype in dtypes.items()
    }


def trace_many(
    requests: Sequence[tuple[OpticalModel, Mapping[str, np.ndarray], Mapping[str, Any]]],
    *,
    seed: int,
    pool: ThreadPool,
) -> list[dict[str, np.ndarray]]:
    """Trace the photons of several models in one job on `pool`.

    Each request is ``(model, photons, options)``; `options` may hold
    `apply_responses`, `direct`, `outputs`, `detected_only` and
    `classify_losses` with the meaning and defaults of `OpticalModel.trace`.
    The result is the list of the dicts that one `OpticalModel.trace` call per
    request returns, with the same values; the blocks and models share the
    pool's threads, which keeps them busy when single blocks are small.
    """
    allowed = {"apply_responses", "direct", "outputs", "detected_only", "classify_losses"}
    items = []
    for model, photons, options in requests:
        unknown = set(options) - allowed
        if unknown:
            raise TypeError(f"unknown trace options: {', '.join(sorted(unknown))}")
        outputs = options.get("outputs")
        items.append((
            model._native,  # noqa: SLF001
            _photon_arrays(photons),
            options.get("apply_responses", True),
            options.get("direct", False),
            None if outputs is None else list(outputs),
            options.get("detected_only", False),
            options.get("classify_losses", True),
        ))
    return _core.trace_many(items, seed=seed, pool=pool)


def _emitter_count(source: _core.Emitter, emitted: float, seed: int, first_id: int) -> int:
    """Poisson number of candidate photons of an emitter (see ``OpticalModel.emitter_photons``)."""
    return int(
        np.random.default_rng([seed, first_id]).poisson(emitted * source.candidates_per_photon())
    )
