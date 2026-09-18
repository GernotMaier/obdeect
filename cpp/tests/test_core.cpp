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

  // T-OBS-002: an off-axis ray remains traceable through the MST structure.
  const auto outer = trace_toy_mst({{2.0, 1.0, 20.0}, {0.0, 0.0, -1.0}}, 3, structured);
  require(outer.status == PhotonStatus::detected || outer.status == PhotonStatus::blocked_mast,
          "off-axis ray must either reach the screen or hit a physical mast");

  // T-SRC-001: artificial blue source emits unit, pupil-bounded directions.
  const auto rays = parallel_blue_cherenkov_rays(256, structured);
  require(rays.size() == 256, "source photon count");
  for (const auto& ray : rays) {
    require(std::abs(norm(ray.direction) - 1.0) < 1e-14, "source direction must be unit length");
    require(ray.position_m.x * ray.position_m.x + ray.position_m.y * ray.position_m.y <=
                structured.mirror_aperture_radius_m * structured.mirror_aperture_radius_m + 1e-12,
            "source samples must stay inside pupil");
  }

  std::cout << "obdeect core tests passed\n";
}
