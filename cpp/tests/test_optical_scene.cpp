#include "obdeect/trace.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
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
  const ModelProvenance provenance{"geometry-oracle", "1", std::string(64, 'a')};
  const OpticalSurfaceRecord mirror{10, {{0.0, 0.0, 0.0}}, FacetShape::square, 2.0,
                                    SurfaceRole::mirror};
  const OpticalSurfaceRecord detector{20, {{0.0, 0.0, 4.0}}, FacetShape::square, 2.0,
                                      SurfaceRole::detector};
  const auto scene = compile_optical_scene({provenance, {detector, mirror}, 4});
  require(scene.has_value(), "T-GEO-001: unordered finite scene compiles");
  const std::vector<Vec3> positions{{0.0, 0.0, 3.0}, {1.1, 0.0, 3.0}};
  const std::vector<Vec3> directions{{0.0, 0.0, -1.0}, {0.0, 0.0, -1.0}};
  const std::vector<double> wavelength{400.0, 400.0};
  const std::vector<double> time{5.0, 5.0};
  const std::vector<double> weight{2.0, 3.0};
  const std::vector<std::uint64_t> ids{1, 2};
  const PhotonBlockView input{positions, directions, wavelength, time, weight, ids};
  const auto detected = trace(*scene, input);
  require(detected.photons.status[1] == PhotonStatus::escaped_scene &&
              detected.photons.surface_id[1] == PhotonResultBlock::kNoSurfaceId,
          "T-GEO-003: ray through finite-mask gap escapes");
  require(detected.photons.status[0] == PhotonStatus::detected &&
              detected.photons.surface_id[0] == detector.id &&
              std::abs(detected.photons.optical_path_m[0] - 7.0) < 1e-12 &&
              std::abs(detected.photons.time_ns[0] - (5.0 + 7.0 / kSpeedOfLightMPerNs)) < 1e-12 &&
              detected.photons.weight[0] == 2.0,
          "T-GEO-004: reflected path reaches detector with complete vacuum path and time");

  const double cosine = std::cos(std::numbers::pi / 8.0);
  const double sine = std::sin(std::numbers::pi / 8.0);
  const OpticalSurfaceRecord tilted_mirror{11, {{0.0, 0.0, 0.0}, {cosine, 0.0, -sine},
                                                {0.0, 1.0, 0.0}, {sine, 0.0, cosine}},
                                           FacetShape::circle, 2.0, SurfaceRole::mirror};
  const OpticalSurfaceRecord blocker{30, {{1.0, 0.0, 1.0}}, FacetShape::circle, 0.5,
                                     SurfaceRole::obscurer};
  const OpticalSurfaceRecord angled_detector{21, {{2.0, 0.0, 2.0}}, FacetShape::circle, 0.5,
                                             SurfaceRole::detector};
  const auto angled = compile_optical_scene({provenance, {angled_detector, tilted_mirror, blocker}, 4});
  require(angled.has_value(), "tilted scene compiles");
  const auto blocked = trace(*angled, input);
  require(blocked.photons.status[0] == PhotonStatus::blocked_obscurer &&
              blocked.photons.surface_id[0] == blocker.id &&
              std::abs(blocked.photons.position_m[0].x - 1.0) < 1e-12 &&
              std::abs(blocked.photons.optical_path_m[0] - (3.0 + std::sqrt(2.0))) < 1e-12,
          "T-GEO-002: nearest post-mirror blocker wins and is attributed");

  const OpticalSurfaceRecord front_blocker{31, {{0.0, 0.0, 2.0}}, FacetShape::circle, 2.0,
                                           SurfaceRole::obscurer};
  const auto front = compile_optical_scene({provenance, {mirror, detector, front_blocker}, 4});
  require(front.has_value(), "front blocker scene compiles");
  const auto pre_mirror = trace(*front, input);
  require(pre_mirror.photons.status[0] == PhotonStatus::blocked_obscurer &&
              pre_mirror.photons.surface_id[0] == front_blocker.id &&
              pre_mirror.photons.optical_path_m[0] == 1.0,
          "T-GEO-005: globally nearest pre-mirror obstruction wins");

  const OpticalSurfaceRecord upper_mirror{40, {{0.0, 0.0, 2.0}}, FacetShape::circle, 2.0,
                                          SurfaceRole::mirror};
  const auto cavity = compile_optical_scene({provenance, {upper_mirror, mirror}, 2});
  require(cavity.has_value(), "bounded mirror scene compiles");
  const std::vector<Vec3> cavity_position{{0.0, 0.0, 1.0}};
  const std::vector<Vec3> cavity_direction{{0.0, 0.0, -1.0}};
  const std::vector<double> cavity_scalar{400.0};
  const std::vector<double> cavity_time{0.0};
  const std::vector<double> cavity_weight{1.0};
  const std::vector<std::uint64_t> cavity_id{3};
  const PhotonBlockView cavity_input{cavity_position, cavity_direction, cavity_scalar, cavity_time,
                                     cavity_weight, cavity_id};
  const auto limited = trace(*cavity, cavity_input);
  require(limited.photons.status[0] == PhotonStatus::interaction_limit &&
              limited.photons.surface_id[0] == upper_mirror.id &&
              std::abs(limited.photons.optical_path_m[0] - 3.0) < 1e-12,
          "T-GEO-006: repeated reflections stop at the declared cap");

  require(!compile_optical_scene({provenance, {mirror, mirror}, 4}),
          "T-GEO-007: duplicate surface identities fail closed");
  require(!compile_optical_scene({provenance, {mirror}, kMaximumSceneInteractions + 1}),
          "T-GEO-009: interaction cap is bounded");
  auto invalid = mirror;
  invalid.frame.z_axis = {0.0, 0.0, -1.0};
  require(!compile_optical_scene({provenance, {invalid}, 4}),
          "T-GEO-008: non-rigid or left-handed transform fails closed");
  std::cout << "optical scene tests passed\n";
}
