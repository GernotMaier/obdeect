#include "obdeect/toy_mst.hpp"

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
  ToyMstConfig bare{};
  bare.include_structure = false;

  // T-TOY-001: the on-axis parallel ray reaches the paraxial focal screen.
  const auto central = trace_toy_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 1, bare);
  require(central.status == PhotonStatus::detected, "central ray must be detected without structure");
  require(central.point_count == 3, "detected path has source, mirror and screen vertices");
  require(std::abs(central.points_m[2].x) < 1e-12 && std::abs(central.points_m[2].y) < 1e-12,
          "central ray must focus on the optical axis");
  require(central.path_length_m > 20.0, "path length must be positive and accumulated");

  // T-OBS-001: the MST camera is a pre-M1 obstruction on the optical axis.
  ToyMstConfig structured{};
  const auto camera = trace_toy_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 2, structured);
  require(camera.status == PhotonStatus::blocked_camera, "camera must shadow the central incoming ray");
  require(camera.point_count == 2, "blocked path terminates at the obstruction");

  // T-OBS-004: camera and supports are closed solids.  This ray clears the
  // incoming camera aperture, reflects from M1, then enters the camera rear.
  const auto post_reflection_camera = trace_toy_mst({{0.6, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 8, structured);
  require(post_reflection_camera.status == PhotonStatus::blocked_camera,
          "camera must also obstruct the reflected segment");
  require(post_reflection_camera.point_count == 3,
          "post-reflection obstruction retains entrance, M1, and obstruction vertices");
  require(std::abs(post_reflection_camera.points_m[2].z -
                   (structured.focal_length_m - structured.camera_half_depth_m)) < 1e-12,
          "reflected ray must stop at the rear camera cap");

  // T-OBS-005: a closed cylinder catches axial hits on its end caps, which a
  // side-wall-only intersection would miss.
  const Ray axial_cylinder{{0.0, 0.0, 2.0}, {0.0, 0.0, -1.0}};
  const auto cap_hit = intersect_closed_finite_cylinder(axial_cylinder, {0.0, 0.0, 0.0},
                                                         {0.0, 0.0, 1.0}, 0.5);
  require(cap_hit && std::abs(*cap_hit - 1.0) < 1e-12, "closed cylinder must intersect its end cap");

  // T-OBS-002: an off-axis ray remains traceable through the MST structure.
  const auto outer = trace_toy_mst({{2.0, 1.0, 20.0}, {0.0, 0.0, -1.0}}, 3, structured);
  require(outer.status == PhotonStatus::detected || outer.status == PhotonStatus::blocked_mast ||
              outer.status == PhotonStatus::blocked_camera,
          "off-axis ray must either reach the screen or hit a physical obstruction");

  // T-SRC-001: artificial blue source emits unit, pupil-bounded directions.
  const auto rays = parallel_blue_cherenkov_rays(256, structured);
  require(rays.size() == 256, "source photon count");
  for (const auto& ray : rays) {
    require(std::abs(norm(ray.direction) - 1.0) < 1e-14, "source direction must be unit length");
    require(ray.position_m.x * ray.position_m.x + ray.position_m.y * ray.position_m.y <=
                structured.mirror_aperture_radius_m * structured.mirror_aperture_radius_m + 1e-12,
            "source samples must stay inside pupil");
  }
  require(parallel_blue_cherenkov_rays(0, structured).empty(), "zero source count must be empty");

  // T-ABI-001: invalid vectors/configuration fail closed rather than divide by zero.
  const auto zero_direction = trace_toy_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, 0.0}}, 4, bare);
  require(zero_direction.status == PhotonStatus::invalid_input, "zero direction must be rejected");
  ToyMstConfig invalid = bare;
  invalid.mirror_aperture_radius_m = invalid.mirror_radius_m + 1.0;
  const auto invalid_scene = trace_toy_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 5, invalid);
  require(invalid_scene.status == PhotonStatus::invalid_input, "impossible spherical cap must be rejected");
  invalid = bare;
  invalid.camera_half_depth_m = invalid.focal_length_m;
  const auto inverted_camera = trace_toy_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 6, invalid);
  require(inverted_camera.status == PhotonStatus::invalid_input,
          "camera support endpoint behind the primary must be rejected");

  // T-OBS-003: zero-radius camera/masts are disabled, not zero-area blockers.
  ToyMstConfig no_obstruction = structured;
  no_obstruction.camera_radius_m = 0.0;
  no_obstruction.mast_radius_m = 0.0;
  const auto unobscured = trace_toy_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 7, no_obstruction);
  require(unobscured.status == PhotonStatus::detected, "zero-radius obstructions must not block");

  std::cout << "obdeect core tests passed\n";
}
