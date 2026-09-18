#include "obdeect/sources.hpp"

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

  // T-SRC-011: a zenith star is a parallel 400-nm source over the pupil.
  const auto star = star_photons(128, 6.0, {});
  require(star.size() == 128, "star count");
  for (const auto& photon : star) {
    require(std::abs(norm(photon.ray.direction) - 1.0) < 1e-14, "star directions are unit");
    require(std::abs(photon.ray.direction.x) < 1e-14 && std::abs(photon.ray.direction.y) < 1e-14,
            "on-axis star is parallel to z");
    require(photon.wavelength_nm == 400.0 && photon.weight == 1.0, "star photon metadata");
  }

  // T-SRC-012: a finite illuminator has pupil-dependent direction and 1/r^2 weight.
  const auto illuminator = illuminator_photons(128, 6.0, {{0.0, 0.0, 30.0}, 400.0, 2.0});
  require(illuminator.size() == 128, "illuminator count");
  require(std::abs(illuminator.front().ray.direction.x - illuminator.back().ray.direction.x) > 1e-4,
          "finite illuminator directions differ across pupil");
  require(illuminator.front().weight > 0.0 && illuminator.back().weight > 0.0,
          "finite illuminator weights are positive");
  require(illuminator.front().time_ns == 0.0 && illuminator.back().time_ns == 0.0,
          "illuminator emission time is not preloaded with flight time");

  // T-SRC-013: laser directions stay inside configured divergence cone.
  constexpr double divergence_rad = 0.02;
  const auto laser = laser_photons(128, 6.0, {{0.0, 0.0, -1.0}, 50.0, divergence_rad, 400.0});
  require(laser.size() == 128, "laser count");
  for (const auto& photon : laser) {
    require(std::acos(std::clamp(-photon.ray.direction.z, -1.0, 1.0)) <= divergence_rad + 1e-12,
            "laser ray remains in divergence cone");
  }

  std::cout << "source tests passed\n";
}
