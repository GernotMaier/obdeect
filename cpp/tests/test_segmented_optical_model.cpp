#include "obdeect/detector_assignment.hpp"
#include "obdeect/trace.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {
void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}
} // namespace

int main() {
  using namespace obdeect;
  const ModelProvenance provenance{"generic-fixture", "1.0.0", std::string(64, 'a')};
  const ImportedFacet facet{
      7, {1.0, 2.0, 3.0}, {0.0, 0.0, 1.0}, 1.2, 16.0, FacetShape::hexagon_flat_y, {1.0, 0.0, 0.0}};
  const auto optical_model = compile_segmented_optical_model({provenance, {facet}});
  require(optical_model.has_value() &&
              optical_model->primary_facets[0].shape == FacetShape::hexagon_flat_y,
          "T-IR-001: compiler preserves generic facet geometry");
  require(!compile_segmented_optical_model({provenance, {facet, facet}}),
          "T-IR-002: duplicate facet IDs fail closed");
  auto invalid_response_model = *optical_model;
  invalid_response_model.pixel_responses = std::make_shared<const PixelResponses>();
  require(!is_valid(invalid_response_model),
          "T-IR-007: compiled model validation checks attached pixel responses");
  auto invalid_assignment_model = *optical_model;
  invalid_assignment_model.detector_assignment =
      std::make_shared<const CompiledDetectorAssignmentGrid>();
  require(!is_valid(invalid_assignment_model),
          "compiled model rejects an invalid detector assignment grid");
  auto invalid = facet;
  invalid.unit_normal = {0.0, 0.0, 0.0};
  require(!compile_segmented_optical_model({provenance, {invalid}}),
          "T-IR-003: missing alignment normals fail closed");

  auto missing_orientation = facet;
  missing_orientation.unit_tangent_u = {};
  require(!compile_segmented_optical_model({provenance, {missing_orientation}}),
          "T-IR-006: polygon orientation cannot be inferred");

  const Ray invalid_position{{std::numeric_limits<double>::quiet_NaN(), 0, 10}, {0, 0, -1}};
  require(trace_segmented_path(invalid_position, 0, 400, *optical_model).status ==
              PhotonStatus::invalid_input,
          "scalar transport rejects nonfinite entrance coordinates");
  require(trace_segmented_path({{1, 2, 10}, {0, 0, -1}}, 0, 0, *optical_model).status ==
              PhotonStatus::invalid_input,
          "scalar transport rejects nonphysical wavelength even without response tables");

  // T-SEG-001: shape bounds use the documented diameter convention and do
  // not turn a hexagon into a circular aperture.
  require(contains_facet_point(facet, {1.0, 2.59, 3.0}), "point inside hexagonal flat boundary");
  require(!contains_facet_point(facet, {1.68, 2.59, 3.0}), "hexagonal corner outside aperture");
  const ImportedFacet square{8,    {0.0, 0.0, 0.0},    {0.0, 0.0, 1.0}, 2.0,
                             16.0, FacetShape::square, {1.0, 0.0, 0.0}};
  const ImportedFacet circle{9, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 2.0, 16.0, FacetShape::circle};
  require(contains_facet_point(square, {0.99, 0.99, 0.0}) &&
              !contains_facet_point(square, {1.01, 0.0, 0.0}),
          "square uses diameter as side length");
  require(contains_facet_point(circle, {0.0, 0.99, 0.0}) &&
              !contains_facet_point(circle, {1.01, 0.0, 0.0}),
          "circle uses diameter as diameter");
  const ImportedDetectorSurface square_detector{12,  {0.0, 0.0, 2.0},    {0.0, 0.0, 1.0},
                                                2.0, FacetShape::square, {1.0, 0.0, 0.0}};
  require(contains_detector_point(square_detector, {0.99, 0.99, 2.0}) &&
              !contains_detector_point(square_detector, {1.01, 0.0, 2.0}),
          "detector and facet use the same square aperture rule");
  auto curved = circle;
  curved.curvature_radius_m = 20.0;
  const auto curved_hit = intersect_segmented_facet({{0.5, 0.0, 5.0}, {0.0, 0.0, -1.0}}, curved);
  require(curved_hit && curved_hit->point_m.z > 0.0 && curved_hit->unit_normal.x > 0.0,
          "spherical panel uses its focal-length curvature and aperture");
  const auto curved_reflection =
      curved_hit ? reflect_specular({0.0, 0.0, -1.0}, curved_hit->unit_normal) : std::nullopt;
  require(curved_reflection.has_value(), "spherical panel reflection is defined");
  const double paraxial_distance = (10.0 - curved_hit->point_m.z) / curved_reflection->z;
  const Vec3 paraxial_hit = curved_hit->point_m + *curved_reflection * paraxial_distance;
  require(std::abs(paraxial_hit.x) < 1.e-3,
          "concave spherical panel sends parallel rays towards its paraxial focus");

  // T-SEG-002: the closest finite panel wins, and reflection exposes the
  // supplied panel normal rather than an inferred dish normal.
  const ImportedFacet far{10, {0.0, 0.0, 2.0}, {0.0, 0.0, 1.0}, 2.0, 16.0, FacetShape::circle};
  const ImportedFacet near{11, {0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}, 2.0, 16.0, FacetShape::circle};
  const auto transport_optical_model = compile_segmented_optical_model({provenance, {near, far}});
  require(transport_optical_model.has_value(), "transport fixture compiles");
  const Ray on_axis{{0.0, 0.0, 5.0}, {0.0, 0.0, -1.0}};
  const auto hit = intersect_segmented_primary(on_axis, *transport_optical_model);
  require(hit.has_value() && hit->facet_id == far.id && std::abs(hit->distance_m - 3.0) < 1.e-12,
          "nearest finite panel is selected");

  // Coincident edges must select the same physical identity for any input order.
  auto overlapping_facet = far;
  overlapping_facet.id = 12;
  const ImportedDetectorSurface detector_a{20, {0, 0, 4}, {0, 0, 1}, 2, FacetShape::circle};
  auto imaging = compile_segmented_optical_model({provenance, {near}, {detector_a}});
  require(imaging.has_value(), "imaging plane fixture compiles");
  const Ray outside_pixel{{3, 0, 5}, {0, 0, -1}};
  require(!intersect_detector_surfaces(outside_pixel, *imaging),
          "finite detector rejects photons outside its aperture");
  imaging->imaging_plane_z_m = 4;
  const auto image_hit = intersect_detector_surfaces(outside_pixel, *imaging);
  require(image_hit && image_hit->distance_m == 1 && image_hit->point_m.x == 3 &&
              image_hit->point_m.z == 4,
          "imaging plane records photons beyond the camera boundary at the declared plane");
  require(!intersect_detector_surfaces({{0, 0, 5}, {1, 0, 0}}, *imaging) &&
              !intersect_detector_surfaces({{0, 0, 5}, {0, 0, 1}}, *imaging),
          "imaging plane rejects parallel and backward intersections");
  auto detector_b = detector_a;
  detector_b.id = 21;
  const ImportedCylinderObscurer cylinder_a{30, {0, 0, 3}, {0, 0, 4}, 1};
  auto cylinder_b = cylinder_a;
  cylinder_b.id = 31;
  for (const bool reverse : {false, true}) {
    const auto tied = compile_segmented_optical_model(
        {provenance,
         reverse ? std::vector<ImportedFacet>{overlapping_facet, far}
                 : std::vector<ImportedFacet>{far, overlapping_facet},
         reverse ? std::vector<ImportedDetectorSurface>{detector_b, detector_a}
                 : std::vector<ImportedDetectorSurface>{detector_a, detector_b},
         reverse ? std::vector<ImportedCylinderObscurer>{cylinder_b, cylinder_a}
                 : std::vector<ImportedCylinderObscurer>{cylinder_a, cylinder_b}});
    require(tied.has_value(), "coincident-boundary optical model compiles");
    require(intersect_segmented_primary(on_axis, *tied)->facet_id == far.id,
            "coincident facets select lowest ID independent of ordering");
    require(intersect_detector_surfaces(on_axis, *tied)->surface_id == detector_a.id,
            "coincident detectors select lowest ID independent of ordering");
    require(intersect_cylinder_obscurers_unchecked(on_axis, *tied)->surface_id == cylinder_a.id,
            "coincident obscurers select lowest ID independent of ordering");
  }

  const std::vector<Vec3> positions{{0.0, 0.0, 5.0}, {3.0, 0.0, 5.0}};
  const std::vector<Vec3> directions{{0.0, 0.0, -1.0}, {0.0, 0.0, -1.0}};
  const std::vector<double> wavelength{400.0, 400.0};
  const std::vector<double> time{7.0, 11.0};
  const std::vector<double> weight{2.0, 3.0};
  const std::vector<std::uint64_t> ids{41, 42};
  const PhotonBlockView input{positions, directions, wavelength, time, weight, ids};
  const auto result = trace(*transport_optical_model, input);
  require(result.photons.status[0] == PhotonStatus::no_detector &&
              result.photons.status[1] == PhotonStatus::missed_primary,
          "T-SEG-003: transport reports explicit no-detector and primary-miss terminals");
  require(result.photons.position_m[0].z == 2.0 && result.photons.direction[0].z > 0.999999999 &&
              std::abs(result.photons.optical_path_m[0] - 3.0) < 1.e-12 &&
              result.photons.weight[0] == 0.0,
          "T-SEG-004: reflection state and unaccepted terminal weight are returned");
  require(result.summary.status_count[static_cast<std::size_t>(PhotonStatus::no_detector)] == 1 &&
              result.summary.status_count[static_cast<std::size_t>(PhotonStatus::missed_primary)] ==
                  1,
          "T-SEG-005: terminal status accounting closes");

  // T-SEG-006: detector surfaces are explicit, finite and generic. The
  // nearest valid post-reflection surface wins; path, time and weight include
  // the complete physical path through the terminal detector intersection.
  const ImportedDetectorSurface farther_detector{
      21, {0.0, 0.0, 4.0}, {0.0, 0.0, 1.0}, 2.0, FacetShape::circle};
  const ImportedDetectorSurface nearer_detector{
      22, {0.0, 0.0, 3.0}, {0.0, 0.0, 1.0}, 2.0, FacetShape::circle};
  const auto detected_optical_model = compile_segmented_optical_model(
      {provenance, {near, far}, {farther_detector, nearer_detector}});
  require(detected_optical_model.has_value(), "optical_model with explicit detectors compiles");
  const auto detected = trace(*detected_optical_model, input);
  auto air_model = *detected_optical_model;
  air_model.propagation_group_index = 1.00023;
  const auto air_result = trace(air_model, input);
  require(std::abs(air_result.photons.time_ns[0] - (7.0 + 4.0 * air_model.propagation_group_index /
                                                              kSpeedOfLightMPerNs)) < 1.e-12 &&
              air_result.photons.position_m[0].z == detected.photons.position_m[0].z &&
              air_result.photons.weight[0] == detected.photons.weight[0],
          "ambient group index changes flight time without changing geometry or throughput");
  air_model.propagation_group_index = 0.9;
  require(!is_valid(air_model), "invalid ambient group index rejects");
  constexpr double speed_of_light_m_per_ns = 0.299792458;
  require(detected.photons.status[0] == PhotonStatus::detected &&
              detected.photons.surface_id[0] == nearer_detector.id &&
              detected.photons.position_m[0].z == 3.0 &&
              detected.photons.direction[0].z > 0.999999999 &&
              std::abs(detected.photons.optical_path_m[0] - 4.0) < 1.e-12 &&
              std::abs(detected.photons.time_ns[0] - (7.0 + 4.0 / speed_of_light_m_per_ns)) <
                  1.e-12 &&
              detected.photons.weight[0] == 2.0,
          "detector hit preserves physical arrival state and input weight");
  require(detected.photons.status[1] == PhotonStatus::missed_primary &&
              detected.photons.surface_id[1] == PhotonResultBlock::kNoSurfaceId &&
              detected.summary.status_count[static_cast<std::size_t>(PhotonStatus::detected)] ==
                  1 &&
              std::abs(detected.summary.detected_weight - 2.0) < 1.e-12,
          "detector terminal status and summary close");

  auto duplicate_detector = nearer_detector;
  duplicate_detector.id = near.id;
  require(!compile_segmented_optical_model({provenance, {near, far}, {duplicate_detector}}),
          "T-SEG-007: detector IDs must not collide with facet IDs");
  auto invalid_detector = nearer_detector;
  invalid_detector.diameter_m = 0.0;
  require(!compile_segmented_optical_model({provenance, {near, far}, {invalid_detector}}),
          "T-SEG-008: invalid detector geometry fails closed");

  // T-SEG-009: model-provided opaque structure is part of the optical_model and wins
  // whenever it is closer than a mirror or detector intersection.
  const ImportedCylinderObscurer obscurer{30, {0.0, 0.0, 4.0}, {0.0, 0.0, 5.0}, 0.2};
  const auto obscured_optical_model =
      compile_segmented_optical_model({provenance, {near, far}, {nearer_detector}, {obscurer}});
  require(obscured_optical_model.has_value(), "optical_model with cylinder obscurer compiles");
  const auto obscured = trace(*obscured_optical_model, input);
  require(obscured.photons.status[0] == PhotonStatus::blocked_obscurer &&
              obscured.photons.surface_id[0] == obscurer.id &&
              obscured.photons.position_m[0].z == 4.0 &&
              std::abs(obscured.photons.optical_path_m[0] - 1.0) < 1.e-12 &&
              std::abs(obscured.photons.time_ns[0] - (7.0 + 1.0 / speed_of_light_m_per_ns)) <
                  1.e-12,
          "opaque cylinder records its incoming terminal state");

  auto invalid_response = SpectralResponse{{300.0, 200.0}, {0.7, 0.9}};
  require(!compile_segmented_optical_model({provenance, {facet}, {}, {}, invalid_response}),
          "unordered spectral response fails during optical model compilation");

  // A linear function of wavelength and incidence is reproduced exactly by the
  // bilinear table, including its four measured boundary knots.
  const SpectralResponse angular_response{{300, 500}, {0.4, 0.6, 0.6, 0.8}, {0, 60}};
  require(angular_response.is_valid() && std::abs(*angular_response.at(400, 30) - 0.6) < 1.e-12 &&
              angular_response.at(300, 0) == 0.4 && angular_response.at(500, 60) == 0.8 &&
              !angular_response.at(299, 30) && !angular_response.at(400, 61) &&
              !angular_response.at(400, std::numeric_limits<double>::quiet_NaN()),
          "mirror response interpolates the complete physical incidence/wavelength domain");
  auto invalid_angular = angular_response;
  invalid_angular.response.pop_back();
  require(!invalid_angular.is_valid(), "incomplete response grid fails validation");
  invalid_angular = angular_response;
  invalid_angular.incidence_angle_deg[1] = 91;
  require(!invalid_angular.is_valid(), "incidence response uses angles to the surface normal");
  auto angled_model = *detected_optical_model;
  angled_model.primary_reflectivity = angular_response;
  const auto angular_detected = trace(angled_model, input);
  require(std::abs(angular_detected.photons.weight[0] - 1.0) < 1.e-12,
          "normal incidence uses measured zero-degree response in transport");

  auto tilted = near;
  tilted.unit_normal = {0.5, 0, std::sqrt(0.75)};
  const auto tilted_model =
      compile_segmented_optical_model({provenance, {tilted}, {}, {}, angular_response});
  require(tilted_model.has_value(), "tilted response fixture compiles");
  const auto angled_path = trace_segmented_path({{0, 0, 5}, {0, 0, -1}}, 0, 400, *tilted_model);
  require(std::abs(angled_path.incidence_primary_deg - 30) < 1.e-12 &&
              std::abs(angled_path.surviving_throughput - 0.6) < 1.e-12,
          "transport uses incidence to the actual selected facet normal");

  const SpectralResponse reflectivity{{300.0, 500.0}, {0.7, 0.9}};
  require(reflectivity.is_valid() && std::abs(*reflectivity.at(400.0) - 0.8) < 1.e-12 &&
              !reflectivity.at(250.0),
          "spectral reflectivity interpolates only within declared range");
  auto lossy_model = *detected_optical_model;
  lossy_model.primary_reflectivity = reflectivity;
  const auto lossy = trace(lossy_model, input);
  require(std::abs(lossy.photons.weight[0] - 1.6) < 1.e-12 &&
              std::abs(lossy.summary.detected_weight - 1.6) < 1.e-12 &&
              std::abs(lossy.summary.response_loss_weight - 0.4) < 1.e-12 &&
              std::abs(lossy.summary.lost_weight - 3.4) < 1.e-12 &&
              std::abs(lossy.summary.input_weight - lossy.summary.detected_weight -
                       lossy.summary.lost_weight) < 1.e-12,
          "partial mirror loss and primary miss close independently of counts");
  lossy_model.detector_surfaces.clear();
  const auto lossy_miss = trace(lossy_model, input);
  require(
      std::abs(lossy_miss.summary.response_loss_weight - 0.4) < 1.e-12 &&
          std::abs(lossy_miss.summary.lost_weight - 5.0) < 1.e-12 &&
          std::abs(lossy_miss.summary
                       .terminal_loss_weight[static_cast<std::size_t>(PhotonStatus::no_detector)] -
                   1.6) < 1.e-12,
      "response loss is accounted even when the reflected photon misses detection");
  // Incoming housing masks must not clip rays on their way to the image.
  auto shadow_model = *detected_optical_model;
  const ImportedDetectorSurface housing{30,  {1, 2, 8},          {0, 0, 1},
                                        0.4, FacetShape::square, {1, 0, 0}};
  auto housing_planes = compile_detector_planes({housing});
  require(housing_planes.has_value(), "finite camera housing compiles");
  shadow_model.incoming_obscurer_planes =
      std::make_shared<const CompiledDetectorPlanes>(std::move(*housing_planes));
  const auto blocked = trace_segmented_path({{1, 2, 10}, {0, 0, -1}}, 0, 400, shadow_model);
  require(blocked.status == PhotonStatus::blocked_obscurer && blocked.point_count == 2 &&
              blocked.terminal_surface_id == 30 && std::abs(blocked.path_length_m - 2) < 1.e-12,
          "housing records its actual incoming intersection and loss");
  const auto below = trace_segmented_path({{1, 2, 7}, {0, 0, -1}}, 0, 400, shadow_model);
  const auto unmasked =
      trace_segmented_path({{1, 2, 7}, {0, 0, -1}}, 0, 400, *detected_optical_model);
  require(below.status == unmasked.status && below.path_length_m == unmasked.path_length_m,
          "incoming-only housing cannot intercept reflected rays");
  std::cout << "segmented optical_model tests passed\n";
}
