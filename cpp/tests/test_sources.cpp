#include "obdeect/cli_parse.hpp"
#include "obdeect/photon_input.hpp"
#include "obdeect/sources.hpp"
#include <algorithm>
#include <array>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>

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
  std::vector<double> spectrum{400};
  for (const auto input : {"300,", ",300", "300,,400", ""}) {
    require(!parse_wavelengths_nm(input, spectrum), "empty spectral entries are rejected");
    require(spectrum == std::vector<double>{400}, "invalid spectra do not modify caller data");
  }

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
    if (result.eof)
      break;
  }
  require(consumed == star.size() && reader.read(batch).count == 0, "EOF loses no photons");
  require(star.size() == 128, "star count");
  for (const auto &photon : star) {
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
  // A plane perpendicular to the propagation direction is one wavefront.
  // Propagating every launch sample to it must give the same arrival time.
  const Vec3 direction = tilted_star.front().ray.direction;
  const Vec3 reference{50.0 * direction.x / direction.z, 50.0 * direction.y / direction.z, 50.0};
  const double reference_plane = dot(direction, reference);
  for (const auto &photon : tilted_star) {
    const double distance = reference_plane - dot(direction, photon.ray.position_m);
    const double arrival_ns = photon.time_ns + distance / kSpeedOfLightMPerNs;
    require(std::abs(arrival_ns) < 1e-12, "off-axis star samples share one wavefront");
  }

  // T-SRC-011c: field angle changes the direction, not the telescope that
  // receives the bundle.  Even a distant source must keep its pupil samples
  // on the primary aperture rather than shifting them by distance*tan(angle).
  const auto distant_star = star_photons(
      128, 6.0, {0.5 * std::numbers::pi / 180.0, -0.2 * std::numbers::pi / 180.0, 10000.0, 400.0});
  require(distant_star.size() == 128, "distant off-axis star count");
  for (const auto &photon : distant_star) {
    const double distance = -photon.ray.position_m.z / photon.ray.direction.z;
    const Vec3 primary_plane_hit = photon.ray.position_m + photon.ray.direction * distance;
    require(std::hypot(primary_plane_hit.x, primary_plane_hit.y) <= 6.1,
            "distant off-axis star remains on telescope pupil");
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
  const double expected_weight_sum = 2.0 * 36.0 / 4.0;
  double illuminator_weight_sum = 0.0;
  for (const auto &photon : illuminator)
    illuminator_weight_sum += photon.weight;
  require(std::abs(illuminator_weight_sum - expected_weight_sum / (30.0 * 30.0)) < 2e-3,
          "illuminator weight includes projected pupil area and sample count");

  // T-SRC-013: laser directions stay inside configured divergence cone.
  constexpr double divergence_rad = 0.02;
  const auto laser =
      laser_photons(128, 6.0, {{0.0, 0.0, -1.0}, {1.5, -0.5, 50.0}, divergence_rad, 400.0});
  require(laser.size() == 128, "laser count");
  for (const auto &photon : laser) {
    require(std::acos(std::clamp(-photon.ray.direction.z, -1.0, 1.0)) <= divergence_rad + 1e-12,
            "laser ray remains in divergence cone");
  }
  require(std::abs(dot(laser.front().ray.position_m - laser.back().ray.position_m,
                       laser.front().ray.direction)) < 0.3,
          "laser launch samples use a plane normal to the beam");
  require(std::abs(laser.front().ray.position_m.x - 1.5) < 1e-12,
          "laser launch plane retains configured transverse origin");

  // T-SRC-014: a pulse preserves any source-specific phase already present
  // and adds a deterministic emission-time distribution.
  const auto original_star_time = tilted_star.front().time_ns;
  auto pulsed_star = tilted_star;
  require(apply_top_hat_emission_times(pulsed_star, 4.0, 8.0), "valid top-hat pulse");
  require(std::abs(pulsed_star.front().time_ns - (original_star_time + 4.0 + 8.0 / 256.0)) < 1e-14,
          "pulse adds start and deterministic width without losing star phase");
  require(!apply_top_hat_emission_times(pulsed_star, 0.0, -1.0),
          "negative pulse width is rejected");

  std::cout << "source tests passed\n";
}
