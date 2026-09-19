#include "obdeect/ctao_models.hpp"
#include "obdeect/ctao_trace.hpp"
#include "obdeect/model_import.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {
void require(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}
}  // namespace

int main() {
  using namespace obdeect;
  // T-MODEL-001: each supported CTAO optical family has a stable identifier,
  // versioned provenance and a nonzero physical entrance pupil.
  for (const auto name : {"LST", "MST", "SST", "SCT"}) {
    const auto model = ctao_reference_model(name);
    require(model.has_value(), "known CTAO model name must resolve");
    require(!model->identifier.empty() && !model->simulation_models_version.empty(),
            "reference model must retain provenance");
    require(model->primary_outer_radius_m > 0.0 && model->focal_plane_radius_m > 0.0,
            "reference model must have physical apertures");
  }
  require(!ctao_reference_model("NOT-A-CTAO-TELESCOPE"), "unknown model must fail closed");
  auto invalid_model = lst_reference_model();
  invalid_model.focal_plane_radius_m = 0.0;
  require(!is_valid(invalid_model), "invalid reference aperture must fail closed");
  const auto imported = import_ctao_reference_model(
      "LST", {"LSTN-design", "6.3.0", "sha256:fixture-for-test"});
  require(imported.has_value() && imported->import.requires_segment_list,
          "import boundary must preserve explicit provenance and requirements");
  require(!import_ctao_reference_model("LST", {"", "6.3.0", "hash"}),
          "provenance-free import must fail closed");

  // T-MODEL-002: a parabolic LST focuses an on-axis parallel ray at its
  // declared focal plane through the common production trace entry point.
  const auto lst = lst_reference_model();
  const auto lst_path = trace_ctao_reference({{4.0, 0.0, 60.0}, {0.0, 0.0, -1.0}}, 1, lst);
  require(lst_path.status == PhotonStatus::detected, "LST on-axis ray must be detected");
  require(lst_path.point_count == 3, "single reflector path has entrance, M1, camera");
  require(std::hypot(lst_path.points_m[2].x, lst_path.points_m[2].y) < 1.0e-10,
          "LST paraboloid must focus on optical axis");
  const auto invalid_path = trace_ctao_reference({{NAN, 0.0, 60.0}, {0.0, 0.0, -1.0}}, 99, lst);
  require(invalid_path.status == PhotonStatus::invalid_input, "non-finite input must fail closed");

  // T-MODEL-003: the MST central sphere obeys its declared f=R/2 reference
  // relation and uses the same tracer as all other telescope types.
  const auto mst = mst_reference_model();
  const auto mst_path = trace_ctao_reference({{0.0, 0.0, 60.0}, {0.0, 0.0, -1.0}}, 2, mst);
  require(mst_path.status == PhotonStatus::detected, "MST central ray must be detected");
  require(std::abs(mst_path.points_m[2].z - mst.focal_plane_z_m) < 1.0e-12,
          "MST focal hit must be on declared camera plane");

  // T-MODEL-004: the two-mirror families are represented as such; a complete
  // segment/camera import is required before their production claim is made.
  const auto sst = sst_reference_model();
  const auto sct = sct_reference_model();
  require(sst.secondary.has_value() && sct.secondary.has_value(), "SC models require secondary surfaces");
  require(sst.secondary->outer_radius_m > 0.0 && sct.secondary->outer_radius_m > 0.0,
          "SC secondary apertures must be physical");

  std::cout << "CTAO reference-model tests passed\n";
}
