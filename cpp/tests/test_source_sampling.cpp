#include "obdeect/source_sampling.hpp"

#include <array>
#include <cassert>

int main() {
  using namespace obdeect;
  const StarSource star{0.02, -0.01, 100, 400};
  const auto sampler = [&](std::uint64_t id) { return sample_star(id, 21, 2, star); };
  std::array<OpticalPhoton, 8> photons;
  assert(fill_source(photons, 0, sampler));
  std::array<OpticalPhoton, 3> block;
  assert(fill_source(block, 3, sampler));
  for (std::size_t i = 0; i < block.size(); ++i) {
    assert(block[i].photon_id == photons[i + 3].photon_id);
    assert(block[i].ray.position_m.x == photons[i + 3].ray.position_m.x);
    assert(block[i].time_ns == photons[i + 3].time_ns);
  }
  assert(sampler(4)->ray.position_m.x != sample_star(4, 22, 2, star)->ray.position_m.x);
  assert(!fill_source(block, std::numeric_limits<std::uint64_t>::max(), sampler));
  assert(!sample_star(0, 1, 0, star));
  assert(!sample_star(0, 1, 1, FiniteStarSource{{0, 0, 40}, 50, 400, 0}));

  const FiniteStarSource finite{{0, 0, 10000}, 50, 400, 2};
  const auto ray = sample_star(3, 21, 2, finite);
  assert(ray);
  assert(std::abs(ray->ray.position_m.z - finite.entrance_z_m) < 1e-10);
  const Vec3 target = sampled_disk(3, 21, 2);
  const double flight = -ray->ray.position_m.z / ray->ray.direction.z;
  assert(norm(ray->ray.position_m + ray->ray.direction * flight - target) < 1e-10);
  assert(std::abs(ray->time_ns + flight / kSpeedOfLightMPerNs -
                  (2 + norm(target - finite.position_m) / kSpeedOfLightMPerNs)) < 1e-8);
  const auto far = sample_star(3, 21, 2, FiniteStarSource{{0, 0, 1e9}, 50, 400, 0});
  assert(far && std::abs(far->ray.direction.x) < 2e-9 && far->ray.direction.z < -0.999999);

  const PointIlluminator illuminator{{0, 0, 10}, 400, 100};
  double captured = 0;
  for (std::size_t id = 0; id < 10000; ++id) {
    const auto photon = sample_illuminator(id, 9, 1, illuminator, 10000);
    assert(photon);
    captured += photon->weight;
  }
  const double expected = 100 * 0.5 * (1 - 10 / std::sqrt(101.0));
  assert(std::abs(captured - expected) < expected * 0.002);
  assert(!sample_illuminator(0, 9, 1, illuminator, 0));
  // Normalization remains a floating-point weight for the full identity range.
  const auto single = sample_illuminator(0, 9, 1, illuminator, 1);
  const auto large = sample_illuminator(0, 9, 1, illuminator, std::uint64_t{1} << 62);
  assert(single && large && std::isfinite(large->weight));
  assert(large->weight == single->weight / static_cast<double>(std::uint64_t{1} << 62));

  const LaserSource laser{{0, 0, -1}, {0, 0, 50}, 0.1, 400};
  double mean_cosine = 0;
  constexpr std::size_t count = 10000;
  for (std::size_t id = 0; id < count; ++id) {
    const auto photon = sample_laser(id, 4, 0.05, laser);
    assert(photon && std::hypot(photon->ray.position_m.x, photon->ray.position_m.y) <= 0.05);
    assert(-photon->ray.direction.z >= std::cos(0.1));
    mean_cosine -= photon->ray.direction.z / count;
  }
  assert(std::abs(mean_cosine - 0.5 * (1 + std::cos(0.1))) < 5e-5);
  const auto collimated = sample_laser(1, 4, 0.05, LaserSource{});
  assert(collimated && collimated->ray.direction.z == -1);
  assert(!sample_laser(1, 4, -1, laser));
}
