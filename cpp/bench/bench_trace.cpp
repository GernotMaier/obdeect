// Single-thread cost of the tracing stages of an optical model.
//
//   bench_trace <model.json> <disk_radius_m> <disk_z_m> <start_z_m> [photons=2000000]
//
// Photons are sampled before the timing starts (star on the optical axis at
// 10 km, 400 nm, seed 3). Times are ns per photon, best of three. Not run by
// ctest.
#include "obdeect/model_json.hpp"
#include "obdeect/source.hpp"
#include "obdeect/trace.hpp"

#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>

using namespace obdeect;
using clk = std::chrono::steady_clock;

namespace {

// Best of three runs of body(), in ns per item.
template <class Body> double best_ns(std::size_t items, const Body &body) {
  double best = 1e30;
  for (int run = 0; run < 3; ++run) {
    const auto start = clk::now();
    body();
    best = std::min(best, std::chrono::duration<double, std::nano>(clk::now() - start).count());
  }
  return best / static_cast<double>(items);
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 5) {
    std::fprintf(stderr, "usage: %s model.json disk_radius_m disk_z_m start_z_m [photons]\n",
                 argv[0]);
    return 2;
  }
  std::ifstream file(argv[1]);
  if (!file) {
    std::fprintf(stderr, "cannot read %s\n", argv[1]);
    return 2;
  }
  std::stringstream text;
  text << file.rdbuf();
  const OpticalModel model = load_model(text.str());
  const PointSource source{{0, 0, 1},
                           1e4,
                           std::strtod(argv[2], nullptr),
                           std::strtod(argv[3], nullptr),
                           std::strtod(argv[4], nullptr),
                           400.0};
  const std::size_t n = argc > 5 ? std::strtoull(argv[5], nullptr, 10) : 2'000'000;
  std::vector<SourcePhoton> photons(n);
  for (std::size_t i = 0; i < n; ++i)
    photons[i] = sample_source(source, i, 3);

  double checksum = 0.0;
  std::array<std::size_t, kStatusCount> counts{};
  const double trace_ns = best_ns(n, [&] {
    counts.fill(0);
    for (std::size_t i = 0; i < n; ++i) {
      const auto r = trace_photon(model, photons[i].ray, 400.0, 0.0, 1.0, i, 3, false);
      ++counts[static_cast<std::size_t>(r.status)];
      checksum += r.position.x + r.time_ns;
    }
  });
  std::printf("trace_photon        %8.1f ns", trace_ns);
  for (std::size_t s = 0; s < kStatusCount; ++s)
    if (counts[s])
      std::printf("  %s %.4f", std::string(to_string(static_cast<Status>(s))).c_str(),
                  static_cast<double>(counts[s]) / static_cast<double>(n));
  std::printf("\n");

  bool failed = false;
  std::size_t hits = 0;
  std::vector<char> hit(n);
  const double primary_ns = best_ns(n, [&] {
    hits = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const auto h = detail::intersect_mirror(model.primary, photons[i].ray, failed);
      hit[i] = h.has_value();
      hits += hit[i];
      if (h)
        checksum += h->hit.t;
    }
  });
  std::printf("primary intersect   %8.1f ns  (hit %.4f)\n", primary_ns,
              static_cast<double>(hits) / static_cast<double>(n));

  if (!model.primary.facets.empty()) {
    std::size_t candidates = 0;
    for (std::size_t i = 0; i < n; ++i)
      model.primary.grid.candidates(photons[i].ray, model.primary.facets.size(),
                                    [&](std::size_t) { ++candidates; });
    std::printf("  candidate facets per ray %.2f\n",
                static_cast<double>(candidates) / static_cast<double>(n));
  }
  if (model.primary.asphere) {
    for (const char want : {char{1}, char{0}}) {
      std::size_t k = 0;
      for (std::size_t i = 0; i < n; ++i)
        k += hit[i] == want;
      const double ns = best_ns(k, [&] {
        for (std::size_t i = 0; i < n; ++i)
          if (hit[i] == want)
            checksum += detail::intersect_mirror(model.primary, photons[i].ray, failed).has_value();
      });
      std::printf("  asphere %-4s      %8.1f ns per ray (%zu rays)\n", want ? "hit" : "miss", ns,
                  k);
    }
  }

  std::size_t blocked = 0;
#ifdef OBDEECT_BENCH_COUNTERS
  Obscurers::counters() = {};
#endif
  const double obscurer_ns = best_ns(n, [&] {
    blocked = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (const auto b = model.obscurers.nearest(photons[i].ray, Leg::incoming, 1e30)) {
        ++blocked;
        checksum += b->t;
      }
  });
  std::printf("obscurers incoming  %8.1f ns  (blocked %.4f, %zu items)\n", obscurer_ns,
              static_cast<double>(blocked) / static_cast<double>(n),
              model.obscurers.items().size());
#ifdef OBDEECT_BENCH_COUNTERS
  std::printf("%s\n", Obscurers::counter_report(3 * n).c_str());
#endif
  std::printf("checksum %.17g\n", checksum);
  return 0;
}
