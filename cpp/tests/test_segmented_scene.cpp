#include "obdeect/trace.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

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
  const ModelProvenance provenance{"generic-fixture", "1.0.0", std::string(64, 'a')};
  const ImportedFacet facet{7, {1.0, 2.0, 3.0}, {0.0, 0.0, 1.0}, 1.2, 16.0,
                             FacetShape::hexagon_flat_y, {1.0, 0.0, 0.0}};
  const auto scene = compile_segmented_scene({provenance, {facet}});
  require(scene.has_value() && scene->primary_facets[0].shape == FacetShape::hexagon_flat_y,
          "T-IR-001: compiler preserves generic facet geometry");
  require(!compile_segmented_scene({provenance, {facet, facet}}),
          "T-IR-002: duplicate facet IDs fail closed");
  auto invalid = facet;
  invalid.unit_normal = {0.0, 0.0, 0.0};
  require(!compile_segmented_scene({provenance, {invalid}}),
          "T-IR-003: missing alignment normals fail closed");

  auto missing_orientation = facet;
  missing_orientation.unit_tangent_u = {};
  require(!compile_segmented_scene({provenance, {missing_orientation}}),
          "T-IR-006: polygon orientation cannot be inferred");

  // T-SEG-001: shape bounds use the documented diameter convention and do
  // not turn a hexagon into a circular aperture.
  require(contains_facet_point(facet, {1.0, 2.59, 3.0}), "point inside hexagonal flat boundary");
  require(!contains_facet_point(facet, {1.68, 2.59, 3.0}), "hexagonal corner outside aperture");
  const ImportedFacet square{8, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0, 16.0,
                              FacetShape::square, {1.0, 0.0, 0.0}};
  const ImportedFacet circle{9, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0, 16.0,
                              FacetShape::circle};
  require(contains_facet_point(square, {0.99, 0.99, 0.0}) &&
              !contains_facet_point(square, {1.01, 0.0, 0.0}),
          "square uses diameter as side length");
  require(contains_facet_point(circle, {0.0, 0.99, 0.0}) &&
              !contains_facet_point(circle, {1.01, 0.0, 0.0}),
          "circle uses diameter as diameter");

  // T-SEG-002: the closest finite panel wins, and reflection exposes the
  // supplied panel normal rather than an inferred dish normal.
  const ImportedFacet far{10, {0.0, 0.0, 2.0}, {0.0, 0.0, 1.0}, 2.0, 16.0,
                           FacetShape::circle};
  const ImportedFacet near{11, {0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}, 2.0, 16.0,
                            FacetShape::circle};
  const auto transport_scene = compile_segmented_scene({provenance, {near, far}});
  require(transport_scene.has_value(), "transport fixture compiles");
  const Ray on_axis{{0.0, 0.0, 5.0}, {0.0, 0.0, -1.0}};
  const auto hit = intersect_segmented_primary(on_axis, *transport_scene);
  require(hit.has_value() && hit->facet_id == far.id && std::abs(hit->distance_m - 3.0) < 1.e-12,
          "nearest finite panel is selected");

  const std::vector<Vec3> positions{{0.0, 0.0, 5.0}, {3.0, 0.0, 5.0}};
  const std::vector<Vec3> directions{{0.0, 0.0, -1.0}, {0.0, 0.0, -1.0}};
  const std::vector<double> wavelength{400.0, 400.0};
  const std::vector<double> time{7.0, 11.0};
  const std::vector<double> weight{2.0, 3.0};
  const std::vector<std::uint64_t> ids{41, 42};
  const PhotonBlockView input{positions, directions, wavelength, time, weight, ids};
  const auto result = trace(*transport_scene, input);
  require(result.photons.status[0] == PhotonStatus::no_detector &&
              result.photons.status[1] == PhotonStatus::missed_primary,
          "T-SEG-003: transport reports explicit no-detector and primary-miss terminals");
  require(result.photons.position_m[0].z == 2.0 && result.photons.direction[0].z > 0.999999999 &&
              std::abs(result.photons.optical_path_m[0] - 3.0) < 1.e-12 && result.photons.weight[0] == 0.0,
          "T-SEG-004: reflection state and unaccepted terminal weight are returned");
  require(result.summary.status_count[static_cast<std::size_t>(PhotonStatus::no_detector)] == 1 &&
              result.summary.status_count[static_cast<std::size_t>(PhotonStatus::missed_primary)] == 1,
          "T-SEG-005: terminal status accounting closes");

  // T-SEG-006: detector surfaces are explicit, finite and generic. The
  // nearest valid post-reflection surface wins; path, time and weight include
  // the complete physical path through the terminal detector intersection.
  const ImportedDetectorSurface farther_detector{21, {0.0, 0.0, 4.0}, {0.0, 0.0, 1.0}, 2.0,
                                                 FacetShape::circle};
  const ImportedDetectorSurface nearer_detector{22, {0.0, 0.0, 3.0}, {0.0, 0.0, 1.0}, 2.0,
                                                FacetShape::circle};
  const auto detected_scene = compile_segmented_scene({provenance, {near, far},
                                                        {farther_detector, nearer_detector}});
  require(detected_scene.has_value(), "scene with explicit detectors compiles");
  const auto detected = trace(*detected_scene, input);
  constexpr double speed_of_light_m_per_ns = 0.299792458;
  require(detected.photons.status[0] == PhotonStatus::detected &&
              detected.photons.surface_id[0] == nearer_detector.id &&
              detected.photons.position_m[0].z == 3.0 && detected.photons.direction[0].z > 0.999999999 &&
              std::abs(detected.photons.optical_path_m[0] - 4.0) < 1.e-12 &&
              std::abs(detected.photons.time_ns[0] - (7.0 + 4.0 / speed_of_light_m_per_ns)) < 1.e-12 &&
              detected.photons.weight[0] == 2.0,
          "detector hit preserves physical arrival state and input weight");
  require(detected.photons.status[1] == PhotonStatus::missed_primary &&
              detected.photons.surface_id[1] == PhotonResultBlock::kNoSurfaceId &&
              detected.summary.status_count[static_cast<std::size_t>(PhotonStatus::detected)] == 1 &&
              std::abs(detected.summary.detected_weight - 2.0) < 1.e-12,
          "detector terminal status and summary close");

  auto duplicate_detector = nearer_detector;
  duplicate_detector.id = near.id;
  require(!compile_segmented_scene({provenance, {near, far}, {duplicate_detector}}),
          "T-SEG-007: detector IDs must not collide with facet IDs");
  auto invalid_detector = nearer_detector;
  invalid_detector.diameter_m = 0.0;
  require(!compile_segmented_scene({provenance, {near, far}, {invalid_detector}}),
          "T-SEG-008: invalid detector geometry fails closed");
  std::cout << "segmented scene tests passed\n";
}
