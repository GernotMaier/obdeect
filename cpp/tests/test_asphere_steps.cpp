// Newton steps of the aspheric intersection (built with OBDEECT_BENCH_COUNTERS).
// With a single proven root the first guess inside the ray's aperture range is
// tried first; a convex secondary-like asphere needs fewer than 3 steps per call
// on average (4.5 when the guess clamped to the start of the range goes first).
#include "check.hpp"

#include "obdeect/random.hpp"
#include "obdeect/surfaces.hpp"

#include <cstdint>
#include <vector>

using namespace obdeect;

namespace {

void secondary_like_steps() {
  Asphere a{4, 3.1084, 0.0, 0.9, 1.0, {0.0, -0.2335, -0.1, 0.05}};
  a.prepare();
  Asphere::counters() = {};
  std::size_t hits = 0;
  bool failed = false;
  for (std::uint64_t k = 0; k < 20000; ++k) {
    const auto u = uniform4(41, k, Stream::source);
    const Vec3 origin{2.0 * (u[0] - 0.5), 2.0 * (u[1] - 0.5), 0.2 * u[2]};
    const Vec3 target{0.9 * (u[3] - 0.5), 0.9 * (u[2] - 0.5), 3.1};
    const Ray ray{origin, normalized(target - origin)};
    hits += a.intersect(ray, kMinDistanceM, [](const Vec3 &) { return true; }, failed).has_value();
    CHECK(!failed);
  }
  const auto &c = Asphere::counters();
  CHECK(hits > 5000);
  CHECK(c.calls == 20000);
  const double steps = static_cast<double>(c.newton_steps) / static_cast<double>(c.calls);
  std::fprintf(stderr, "newton steps per call %.3f, hits %zu\n", steps, hits);
  CHECK(steps <= 3.0);
}

} // namespace

void all_tests() { secondary_like_steps(); }

int main() { return run_checks(all_tests); }
