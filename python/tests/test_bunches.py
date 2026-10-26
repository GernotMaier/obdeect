"""Photon bunches: thinning, roulette, fractional sizes, frames and invariances."""

import math

import numpy as np
import pytest
from obdeect import OpticalModel, ThreadPool
from obdeect.compiler import model_hash
from scipy.stats import ks_2samp
from telescopes import record, segmented

IDENTITY = np.eye(3)
ORIGIN = (0.0, 0.0, 0.0)
SEED = 20261026
TOWARDS = (0.02, 0.01, math.sqrt(1 - 0.02**2 - 0.01**2))

TRACE_COLUMNS = ("position_m", "direction", "time_ns", "weight")


def compiled(**changes):
    records, files = segmented(**changes)
    return OpticalModel.from_records(records, files.__getitem__, configuration_seed=1)


def edited(model, change):
    """A copy of `model` whose document `change(document)` modified."""
    document = model.document
    change(document)
    document["sha256"] = model_hash(document)
    return OpticalModel(document)


def constant_reflectivity(value):
    def change(document):
        document["optics"]["primary"]["reflectivity"] = {
            "wavelength_min_nm": 300.0,
            "wavelength_step_nm": 100.0,
            "incidence_deg": [0.0],
            "values": [[value], [value]],
        }

    return change


def graded_responses(document):
    """Reflectivity falling with wavelength and a degradation gradient over each facet."""
    primary = document["optics"]["primary"]
    primary["reflectivity"] = {
        "wavelength_min_nm": 300.0,
        "wavelength_step_nm": 100.0,
        "incidence_deg": [0.0, 5.0],
        "values": [[0.9, 0.8], [0.8, 0.7], [0.7, 0.6]],
    }
    axis = [-0.6, 0.0, 0.6]
    primary["degradation"] = {
        "x_m": axis,
        "y_m": axis,
        "values": [[0.3, 0.3, 0.3], [0.6, 0.6, 0.6], [0.95, 0.95, 0.95]],
        "outside": "clamp",
    }


def star_bunches(model, count, *, first_id=0, size=1.0, wavelength_nm=400.0, towards=TOWARDS):
    """Star photons as bunches (identity frame, unit 1) and the plane height of their input."""
    photons = model.star_photons(count, seed=SEED, first_id=first_id, towards_source=towards)
    plane_z = float(photons["position_m"][0, 2])
    assert np.all(photons["position_m"][:, 2] == plane_z)
    bunches = {
        "x": photons["position_m"][:, 0].astype(np.float32),
        "y": photons["position_m"][:, 1].astype(np.float32),
        "cx": photons["direction"][:, 0].astype(np.float32),
        "cy": photons["direction"][:, 1].astype(np.float32),
        "time_ns": photons["time_ns"].astype(np.float32),
        "zem": np.full(count, 1000.0, dtype=np.float32),
        "size": np.full(count, size, dtype=np.float32),
        "wavelength_nm": np.full(count, wavelength_nm, dtype=np.float32),
        "bunch_id": np.arange(first_id, first_id + count, dtype=np.uint64),
    }
    return bunches, plane_z


def trace_bunches(model, bunches, plane_z, **options):
    options.setdefault("seed", SEED)
    return model.trace_bunches(
        bunches,
        frame_rotation=options.pop("frame_rotation", IDENTITY),
        frame_offset_m=options.pop("frame_offset_m", ORIGIN),
        length_unit_m=options.pop("length_unit_m", 1.0),
        plane_z=plane_z,
        **options,
    )


def as_photons(bunches, plane_z):
    """The photons of identity-frame bunches, as `trace` takes them (float32 values widened)."""
    n = len(bunches["x"])
    cx = bunches["cx"].astype(np.float64)
    cy = bunches["cy"].astype(np.float64)
    cz = -np.sqrt(np.maximum(0.0, 1.0 - cx * cx - cy * cy))
    position = np.stack(
        [
            bunches["x"].astype(np.float64),
            bunches["y"].astype(np.float64),
            np.full(n, plane_z),
        ],
        axis=1,
    )
    return {
        "position_m": position,
        "direction": np.stack([cx, cy, cz], axis=1),
        "wavelength_nm": bunches["wavelength_nm"].astype(np.float64),
        "time_ns": bunches["time_ns"].astype(np.float64),
        "weight": np.ones(n),
        "photon_id": bunches["bunch_id"],
    }


def constant_acceptance(value):
    return (np.array([100.0, 1000.0]), np.array([0.0, 1.0e6]), np.full((2, 2), value))


def test_without_responses_bunches_equal_traced_photons_bitwise():
    model = compiled(telescope_transmission=None)
    bunches, plane_z = star_bunches(model, 20000)
    out = trace_bunches(model, bunches, plane_z)
    ref = model.trace(as_photons(bunches, plane_z), seed=SEED, detected_only=True)
    assert len(out["weight"]) == len(ref["weight"]) > 1000
    for name in TRACE_COLUMNS:
        assert np.array_equal(out[name], ref[name]), name
    assert np.array_equal(out["wavelength_nm"], np.full(len(ref["weight"]), 400.0))
    assert np.all(out["sub_photon"] == 0)
    assert out["bunch_index"].dtype == np.uint64 and out["sub_photon"].dtype == np.uint32
    detected = model.trace(
        as_photons(bunches, plane_z), seed=SEED, detected_only=True, outputs=["index"]
    )
    assert np.array_equal(out["bunch_index"], detected["index"])


def test_expected_count_with_reflectivity_and_acceptance():
    model = edited(compiled(telescope_transmission=None), constant_reflectivity(0.8))
    n = 200_000
    bunches, plane_z = star_bunches(model, n)
    geometry = model.trace(
        as_photons(bunches, plane_z), seed=SEED, apply_responses=False, detected_only=True
    )
    expected = 0.3 * 0.8 * len(geometry["weight"])
    out = trace_bunches(model, bunches, plane_z, acceptance=constant_acceptance(0.3))
    count = len(out["weight"])
    print(f"expected count: {count} vs {expected:.1f}, limit {5 * math.sqrt(expected):.1f}")
    assert abs(count - expected) < 5 * math.sqrt(expected)
    assert np.all(out["weight"] == 1.0)
    # The thinning already contains the response bound, so the ratio is 1 and the
    # roulette keeps everything; without it the weight is that ratio.
    kept = trace_bunches(
        model, bunches, plane_z, acceptance=constant_acceptance(0.3), roulette=False
    )
    assert len(kept["weight"]) == count
    assert np.allclose(kept["weight"], 1.0, rtol=2e-12, atol=0.0)


def test_fractional_sizes():
    model = compiled(telescope_transmission=None)
    n = 100_000
    bunches, plane_z = star_bunches(model, n)
    geometry = model.trace(
        as_photons(bunches, plane_z), seed=SEED, apply_responses=False, detected_only=True
    )
    detected = len(geometry["weight"])
    expected = 2.5 * detected

    whole = {**bunches, "size": np.full(n, 2.5, dtype=np.float32)}
    count_whole = len(trace_bunches(model, whole, plane_z)["weight"])

    # The same photons as bunches of 2.0 and of 0.5 (other ids).
    split = {
        name: np.concatenate([values, values]) for name, values in bunches.items() if name != "size"
    }
    split["size"] = np.concatenate([np.full(n, 2.0), np.full(n, 0.5)]).astype(np.float32)
    split["bunch_id"] = np.arange(2 * n, dtype=np.uint64) + np.uint64(10_000_000)
    count_split = len(trace_bunches(model, split, plane_z)["weight"])

    limit = 5 * math.sqrt(expected)
    print(f"fractional sizes: {count_whole}, {count_split} vs {expected:.1f}, limit {limit:.1f}")
    assert abs(count_whole - expected) < limit
    assert abs(count_split - expected) < limit
    assert abs(count_whole - count_split) < limit


def test_mean_rows_per_bunch_equal_the_size():
    # Vertical light on the mirror at x = 1.25 m, away from the camera shadow: all bunches hit.
    model = compiled(telescope_transmission=None)
    rng = np.random.default_rng(5)
    n = 100_000
    for size in (2.5, 0.5, 3.2):
        bunches = {
            "x": (1.25 + rng.uniform(-0.3, 0.3, n)).astype(np.float32),
            "y": rng.uniform(-0.3, 0.3, n).astype(np.float32),
            "cx": np.zeros(n, dtype=np.float32),
            "cy": np.zeros(n, dtype=np.float32),
            "time_ns": np.zeros(n, dtype=np.float32),
            "zem": np.full(n, 8000.0, dtype=np.float32),
            "size": np.full(n, size, dtype=np.float32),
            "wavelength_nm": np.full(n, 400.0, dtype=np.float32),
            "bunch_id": np.arange(n, dtype=np.uint64),
        }
        out = trace_bunches(model, bunches, 30.0)
        mean = len(out["weight"]) / n
        limit = 5 * math.sqrt(size * n) / n
        print(f"mean rows per bunch of size {size}: {mean:.5f}, limit {limit:.5f}")
        assert abs(mean - np.float32(size)) < limit
        assert np.max(out["sub_photon"]) == math.ceil(size) - 1
        assert np.all(np.diff(out["bunch_index"].astype(np.int64)) >= 0)


def test_survivors_follow_the_weighted_trace_distribution():
    model = edited(compiled(telescope_transmission=None), graded_responses)
    n = 400_000
    bunches, plane_z = star_bunches(model, n, wavelength_nm=450.0)
    out = trace_bunches(model, bunches, plane_z)
    # Reference: other photons, traced with weights, then kept with probability = weight.
    reference_bunches, reference_z = star_bunches(
        model, n, first_id=50_000_000, wavelength_nm=450.0
    )
    traced = model.trace(as_photons(reference_bunches, reference_z), seed=SEED, detected_only=True)
    assert np.all(traced["weight"] <= 1.0)
    keep = np.random.default_rng(77).random(len(traced["weight"])) < traced["weight"]
    assert len(out["weight"]) > 20_000 and keep.sum() > 20_000
    for axis in (0, 1):
        statistic = ks_2samp(out["position_m"][:, axis], traced["position_m"][keep, axis])
        print(f"KS on axis {axis}: D={statistic.statistic:.5f}, p={statistic.pvalue:.4f}")
        assert statistic.pvalue > 1e-3
    # The mean weight also agrees: survivors / traced photons.
    expected = traced["weight"].sum()
    assert abs(len(out["weight"]) - expected) < 5 * math.sqrt(expected)


def test_threads_pool_and_splitting_do_not_change_the_result():
    model = edited(
        compiled(
            telescope_transmission=None,
            mirror_reflection_random_angle=record(
                "mirror_reflection_random_angle", [0.01, 0.2, 0.05], ["deg", "null", "deg"]
            ),
        ),
        graded_responses,
    )
    n = 30_000
    bunches, plane_z = star_bunches(model, n, size=2.5, wavelength_nm=0.0)
    density = (np.array([300.0, 600.0]), np.array([1.0, 1.0]))
    options = {"acceptance": constant_acceptance(0.9), "wavelength_density": density}
    one = trace_bunches(model, bunches, plane_z, threads=1, **options)
    assert len(one["weight"]) > 1000
    four = trace_bunches(model, bunches, plane_z, threads=4, **options)
    pool = ThreadPool(3)
    pooled = trace_bunches(model, bunches, plane_z, pool=pool, **options)
    for other in (four, pooled):
        assert set(other) == set(one)
        for name in one:
            assert np.array_equal(one[name], other[name]), name

    cut = 12345
    first = {name: values[:cut] for name, values in bunches.items()}
    second = {name: values[cut:] for name, values in bunches.items()}
    a = trace_bunches(model, first, plane_z, pool=pool, **options)
    b = trace_bunches(model, second, plane_z, pool=pool, **options)
    b["bunch_index"] = b["bunch_index"] + np.uint64(cut)
    for name in one:
        assert np.array_equal(one[name], np.concatenate([a[name], b[name]])), name

    # Sub-photons of one bunch draw their own scatter: same bunch, other position.
    index = one["bunch_index"].astype(np.int64)
    rows0 = np.flatnonzero(one["sub_photon"] == 0)
    rows1 = np.flatnonzero(one["sub_photon"] == 1)
    common, at0, at1 = np.intersect1d(index[rows0], index[rows1], return_indices=True)
    assert len(common) > 100
    assert np.all(one["position_m"][rows0[at0], 0] != one["position_m"][rows1[at1], 0])


def test_invalid_input_is_rejected_before_tracing():
    model = compiled(telescope_transmission=None)
    bunches, plane_z = star_bunches(model, 100)
    slanted = IDENTITY.copy()
    slanted[0, 0] = 1.0 + 1e-9
    with pytest.raises(ValueError, match="orthonormal"):
        trace_bunches(model, bunches, plane_z, frame_rotation=slanted)
    with pytest.raises(ValueError, match="length unit"):
        trace_bunches(model, bunches, plane_z, length_unit_m=0.0)
    with pytest.raises(ValueError, match="acceptance"):
        trace_bunches(model, bunches, plane_z, acceptance=constant_acceptance(1.01))
    with pytest.raises(ValueError, match="acceptance"):
        trace_bunches(model, bunches, plane_z, acceptance=constant_acceptance(-0.01))
    repeated = (np.array([100.0, 100.0]), np.array([0.0, 1.0e6]), np.full((2, 2), 0.5))
    with pytest.raises(ValueError, match="acceptance"):
        trace_bunches(model, bunches, plane_z, acceptance=repeated)
    unsampled = {**bunches, "wavelength_nm": np.zeros(100, dtype=np.float32)}
    with pytest.raises(ValueError, match="density"):
        trace_bunches(
            model, unsampled, plane_z, wavelength_density=(np.array([0.0, 500.0]), np.ones(2))
        )
    with pytest.raises(ValueError, match="wavelength"):
        trace_bunches(model, unsampled, plane_z)
    for name, bad in (("size", 0.0), ("size", np.nan), ("x", np.inf), ("cx", 1.5)):
        broken = {key: values.copy() for key, values in bunches.items()}
        broken[name][50] = bad
        with pytest.raises(ValueError):
            trace_bunches(model, broken, plane_z)
    with pytest.raises(ValueError, match="same length"):
        trace_bunches(model, {**bunches, "y": bunches["y"][:-1]}, plane_z)


def test_sampled_wavelengths_follow_the_density():
    model = compiled(telescope_transmission=None)
    nodes = np.linspace(300.0, 600.0, 301)
    values = 1.0 / nodes**2
    n = 300_000
    bunches, plane_z = star_bunches(model, n, wavelength_nm=0.0)
    out = trace_bunches(model, bunches, plane_z, wavelength_density=(nodes, values))

    # Exact moments of the piecewise linear density (Simpson is exact up to cubics).
    a, b = nodes[:-1], nodes[1:]
    mid = 0.5 * (a + b)
    density_mid = 0.5 * (values[:-1] + values[1:])

    def integral(power):
        return np.sum(
            (b - a)
            / 6
            * (a**power * values[:-1] + 4 * mid**power * density_mid + b**power * values[1:])
        )

    norm = integral(0)
    mean = integral(1) / norm
    sigma = math.sqrt(integral(2) / norm - mean**2)
    count = len(out["wavelength_nm"])
    sampled = out["wavelength_nm"].mean()
    print(f"wavelength mean {sampled:.4f} vs {mean:.4f}, limit {5 * sigma / math.sqrt(count):.4f}")
    assert count > 10_000
    assert abs(sampled - mean) < 5 * sigma / math.sqrt(count)
    assert out["wavelength_nm"].min() >= 300.0 and out["wavelength_nm"].max() <= 600.0


def test_frame_places_the_input_plane():
    # Rotating the input frame by 90 degrees about z and shifting it moves the photons.
    model = compiled(telescope_transmission=None)
    bunches, plane_z = star_bunches(model, 20000)
    base = trace_bunches(model, bunches, plane_z)
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    # Input (x, y) = rotation^T (telescope x, y): the same photons, now in a turned frame.
    turned = {
        **bunches,
        "x": bunches["y"].copy(),
        "y": -bunches["x"],
        "cx": bunches["cy"].copy(),
        "cy": -bunches["cx"],
    }
    out = trace_bunches(model, turned, plane_z, frame_rotation=rotation)
    assert np.array_equal(out["bunch_index"], base["bunch_index"])
    assert np.allclose(out["position_m"], base["position_m"], atol=1e-9)
    # Lengths in centimetres.
    centimetres = {
        **bunches,
        "x": bunches["x"] * np.float32(100.0),
        "y": bunches["y"] * np.float32(100.0),
    }
    scaled = trace_bunches(model, centimetres, plane_z * 100.0, length_unit_m=0.01)
    assert abs(len(scaled["weight"]) - len(base["weight"])) < 0.01 * len(base["weight"])
