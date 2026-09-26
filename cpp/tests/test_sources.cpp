#include "obdeect/sources.hpp"
#include "obdeect/photon_input.hpp"
#include <algorithm>
#include <array>

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
  MemoryPhotonReader reader(star, {1, 2, 3, 4});
  std::array<OpticalPhoton, 7> batch{};
  std::size_t consumed = 0;
  for (;;) {
    const auto result = reader.read(batch);
    require(result.context.event_id == 2 && result.context.telescope_id == 4,
            "batch context is preserved");
    for (std::size_t i = 0; i < result.count; ++i)
      require(batch[i].photon_id == consumed + i, "chunking preserves photon identity and order");
    consumed += result.count;
    if (result.eof) break;
  }
  require(consumed == star.size() && reader.read(batch).count == 0, "EOF loses no photons");
  require(star.size() == 128, "star count");
  for (const auto& photon : star) {
    require(std::abs(norm(photon.ray.direction) - 1.0) < 1e-14, "star directions are unit");
    require(std::abs(photon.ray.direction.x) < 1e-14 && std::abs(photon.ray.direction.y) < 1e-14,
            "on-axis star is parallel to z");
    require(photon.wavelength_nm == 400.0 && photon.weight == 1.0, "star photon metadata");
  }

  // T-SRC-011b: an off-axis plane wave carries the phase across its launch
  // plane in emission time, rather than pretending every sample is emitted
  // simultaneously from a horizontal sheet.
  const auto tilted_star = star_photons(128, 6.0, {0.02, 0.0, 50.0, 400.0});
  require(tilted_star.front().time_ns != tilted_star.back().time_ns,
          "off-axis star has launch-plane phase times");

  // T-SRC-012: a finite illuminator has pupil-dependent direction and 1/r^2 weight.
  const auto illuminator = illuminator_photons(128, 6.0, {{0.0, 0.0, 30.0}, 400.0, 2.0});
  require(illuminator.size() == 128, "illuminator count");
  require(std::abs(illuminator.front().ray.direction.x - illuminator.back().ray.direction.x) > 1e-4,
          "finite illuminator directions differ across pupil");
  require(illuminator.front().weight > 0.0 && illuminator.back().weight > 0.0,
          "finite illuminator weights are positive");
  require(illuminator.front().time_ns == 0.0 && illuminator.back().time_ns == 0.0,
          "illuminator emission time is not preloaded with flight time");
  const double expected_weight_sum = 2.0 * 36.0 / 4.0;
  double illuminator_weight_sum = 0.0;
  for (const auto& photon : illuminator) illuminator_weight_sum += photon.weight;
  require(std::abs(illuminator_weight_sum - expected_weight_sum / (30.0 * 30.0)) < 2e-3,
          "illuminator weight includes projected pupil area and sample count");

  // T-SRC-013: laser directions stay inside configured divergence cone.
  constexpr double divergence_rad = 0.02;
  const auto laser = laser_photons(
      128, 6.0, {{0.0, 0.0, -1.0}, {1.5, -0.5, 50.0}, divergence_rad, 400.0});
  require(laser.size() == 128, "laser count");
  for (const auto& photon : laser) {
    require(std::acos(std::clamp(-photon.ray.direction.z, -1.0, 1.0)) <= divergence_rad + 1e-12,
            "laser ray remains in divergence cone");
  }
  require(std::abs(dot(laser.front().ray.position_m - laser.back().ray.position_m,
                       laser.front().ray.direction)) < 0.3,
          "laser launch samples use a plane normal to the beam");
  require(std::abs(laser.front().ray.position_m.x - 1.5) < 6.0,
          "laser launch plane retains configured transverse origin");

  std::cout << "source tests passed\n";
}
