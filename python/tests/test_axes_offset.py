"""Axes offset of a telescope: the optical origin displaced from the rotation axes."""

import math

import numpy as np
import pytest
from obdeect import OpticalModel
from obdeect.formats import ModelError
from telescopes import record, segmented

SEED = 20261029
THETA_DEG = 20.0
OFFSET_CM = 160.0
ORIGIN = (0.0, 0.0, 0.0)
COS, SIN = math.cos(math.radians(THETA_DEG)), math.sin(math.radians(THETA_DEG))
# Ground frame to telescope frame at zenith THETA_DEG and azimuth 0 (sim_telarray convention):
# p_telescope = ROTATION @ p_ground. The ground offset O = (1.6, 0, 0) m, the sim_telarray
# axes_offsets (160 cm, 0), has the telescope-frame image ROTATION @ O.
ROTATION = np.array([[COS, 0.0, -SIN], [0.0, 1.0, 0.0], [SIN, 0.0, COS]])
GROUND_OFFSET_M = np.array([OFFSET_CM / 100.0, 0.0, 0.0])
EXPECTED_OFFSET_M = ROTATION @ GROUND_OFFSET_M  # (O cos, 0, O sin)


def compiled(axes_offsets=None, zenith_deg=THETA_DEG):
    """Single-mirror telescope without transmission, with an optional axes_offsets record (cm)."""
    changes = {"telescope_transmission": None}
    if axes_offsets is not None:
        changes["axes_offsets"] = record("axes_offsets", axes_offsets, ["cm", "cm"])
    records, files = segmented(**changes)
    return OpticalModel.from_records(
        records, files.__getitem__, configuration_seed=1, zenith_deg=zenith_deg
    )


@pytest.fixture(scope="module")
def reference():
    return compiled()


@pytest.fixture(scope="module")
def offset_model():
    return compiled([OFFSET_CM, 0.0])


def ground_bunches(count, *, first_id=0):
    """Bunches with random positions in [-0.5, 0.5] m and small directions (ground frame)."""
    rng = np.random.default_rng(SEED)
    return {
        "x": rng.uniform(-0.5, 0.5, count).astype(np.float32),
        "y": rng.uniform(-0.5, 0.5, count).astype(np.float32),
        "cx": rng.uniform(-0.02, 0.02, count).astype(np.float32),
        "cy": rng.uniform(-0.02, 0.02, count).astype(np.float32),
        "time_ns": np.zeros(count, dtype=np.float32),
        "zem": np.full(count, 1000.0, dtype=np.float32),
        "size": np.ones(count, dtype=np.float32),
        "wavelength_nm": np.full(count, 400.0, dtype=np.float32),
        "bunch_id": np.arange(first_id, first_id + count, dtype=np.uint64),
    }


def trace_ground(model, bunches, height_m, rotation=ROTATION, offset_m=ORIGIN):
    return model.trace_bunches(
        bunches,
        seed=SEED,
        frame_rotation=rotation,
        frame_offset_m=offset_m,
        length_unit_m=1.0,
        plane_z=height_m,
    )


def telescope_photons(bunches, height_m, rotation, offset_m):
    """Photons of ground-frame bunches as `trace` takes them: R P + offset, R d."""
    n = len(bunches["x"])
    cx = bunches["cx"].astype(np.float64)
    cy = bunches["cy"].astype(np.float64)
    ground_direction = np.stack([cx, cy, -np.sqrt(1.0 - cx * cx - cy * cy)], axis=1)
    ground_position = np.stack(
        [
            bunches["x"].astype(np.float64),
            bunches["y"].astype(np.float64),
            np.full(n, height_m),
        ],
        axis=1,
    )
    return {
        "position_m": ground_position @ rotation.T + np.asarray(offset_m, dtype=np.float64),
        "direction": ground_direction @ rotation.T,
        "wavelength_nm": bunches["wavelength_nm"].astype(np.float64),
        "time_ns": bunches["time_ns"].astype(np.float64),
        "weight": np.ones(n),
        "photon_id": bunches["bunch_id"],
    }


def test_axes_offset_is_stored_in_the_telescope_frame(offset_model):
    assert offset_model.axes_offset_m == pytest.approx(tuple(EXPECTED_OFFSET_M))
    assert offset_model.document["optics"]["axes_offset_m"] == pytest.approx(
        list(EXPECTED_OFFSET_M)
    )


def test_axes_offset_is_zero_without_the_record(reference):
    assert reference.axes_offset_m == (0.0, 0.0, 0.0)
    assert "axes_offset_m" not in reference.document["optics"]


def test_zero_axes_offset_gives_the_same_model_as_no_record(reference):
    model = compiled([0.0, 0.0])
    assert model.axes_offset_m == (0.0, 0.0, 0.0)
    assert model.document["optics"] == reference.document["optics"]


def test_axes_offset_follows_the_zenith_angle():
    model = compiled([OFFSET_CM, 0.0], zenith_deg=0.0)
    assert model.axes_offset_m == pytest.approx((1.6, 0.0, 0.0))


def test_non_zero_perpendicular_component_is_not_supported():
    with pytest.raises(NotImplementedError, match="perpendicular"):
        compiled([OFFSET_CM, 10.0])


def test_axes_offset_needs_the_zenith_angle():
    with pytest.raises(ModelError, match="zenith"):
        compiled([OFFSET_CM, 0.0], zenith_deg=None)


def test_bunch_at_the_array_origin_is_traced_at_r_of_minus_offset(reference, offset_model):
    # The bunch starts at the ground origin of the input plane (x = y = 0, z = H); its photon
    # is at R (P - O) = R (0, 0, H) - R O in the telescope frame.
    height_m = 30.0
    bunches = ground_bunches(3000)
    bunches["x"][:] = 0.0
    bunches["y"][:] = 0.0
    out = trace_ground(offset_model, bunches, height_m)
    expected = telescope_photons(bunches, height_m, ROTATION, -EXPECTED_OFFSET_M)
    ref = reference.trace(expected, seed=SEED, detected_only=True, outputs=["position_m", "index"])
    assert len(out["weight"]) == len(ref["position_m"]) > 0
    assert np.array_equal(out["bunch_index"], ref["index"])
    np.testing.assert_allclose(out["position_m"], ref["position_m"], atol=1e-9)


def test_model_offset_equals_explicit_frame_offset(reference, offset_model):
    # The offset model equals a model without offset whose frame offset is -R O, the ground
    # offset O seen through the frame rotation.
    height_m = 30.0
    bunches = ground_bunches(5000)
    explicit = trace_ground(reference, bunches, height_m, offset_m=-ROTATION @ GROUND_OFFSET_M)
    model = trace_ground(offset_model, bunches, height_m)
    assert len(explicit["weight"]) > 0
    assert np.array_equal(model["bunch_index"], explicit["bunch_index"])
    for name in ("position_m", "direction", "time_ns", "weight", "wavelength_nm"):
        assert np.array_equal(model[name], explicit[name]), name


def test_input_plane_height_is_taken_at_the_offset_origin(offset_model):
    # The input origin of the plane is at height R_zz H - O_z in the telescope frame.
    top = offset_model.max_incoming_obscurer_z_m()
    bunches = ground_bunches(100)
    # Plane height in the telescope frame: R_zz H - O_z; it must reach the obscurer top.
    limit = (top + EXPECTED_OFFSET_M[2]) / ROTATION[2, 2]
    with pytest.raises(ValueError, match="input plane below incoming obscurer"):
        trace_ground(offset_model, bunches, limit - 0.05)
    assert "weight" in trace_ground(offset_model, bunches, limit + 0.05)


def test_star_tracing_is_unaffected_by_the_offset(reference, offset_model):
    # Star photons are sampled and traced in the telescope frame, where the offset does not act.
    stars = reference.star_photons(4000, seed=SEED, towards_source=(0.02, 0.01, 0.99))
    stars_offset = offset_model.star_photons(4000, seed=SEED, towards_source=(0.02, 0.01, 0.99))
    for name in stars:
        np.testing.assert_array_equal(stars[name], stars_offset[name], err_msg=name)
    traced = reference.trace(stars, seed=SEED)
    traced_offset = offset_model.trace(stars, seed=SEED)
    for name in traced:
        np.testing.assert_array_equal(traced[name], traced_offset[name], err_msg=name)
