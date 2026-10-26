// Single-thread cost of the tracing stages of an optical model.
//
//   bench_trace <model.json> <disk_radius_m> <disk_z_m> <start_z_m> [photons=2000000]
//               [classify_losses=1]
//
// classify_losses 0 times the trace_photon lines with classify_losses false.
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
    std::fprintf(stderr,
                 "usage: %s model.json disk_radius_m disk_z_m start_z_m [photons] "
                 "[classify_losses]\n",
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
  const bool classify = argc > 6 ? std::strtol(argv[6], nullptr, 10) != 0 : true;
  std::vector<SourcePhoton> photons(n);
  for (std::size_t i = 0; i < n; ++i)
    photons[i] = sample_source(source, i, 3);

  double checksum = 0.0;
  std::array<std::size_t, kStatusCount> counts{};
  const double trace_ns = best_ns(n, [&] {
    counts.fill(0);
    for (std::size_t i = 0; i < n; ++i) {
      const auto r =
          trace_photon(model, photons[i].ray, 400.0, 0.0, 1.0, i, 3, false, false, true, classify);
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
#ifdef OBDEECT_BENCH_COUNTERS
  SegmentGrid::counters() = {};
#endif
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
#ifdef OBDEECT_BENCH_COUNTERS
  // Per primary hit of the timed runs (best_ns runs the body three times).
  if (hits && SegmentGrid::counters().calls)
    std::printf("  %s\n", SegmentGrid::counter_report(3 * hits).c_str());
#endif

  if (!model.primary.facets.empty()) {
    std::size_t candidates = 0;
    for (std::size_t i = 0; i < n; ++i)
      model.primary.grid.candidates(photons[i].ray, model.primary.facets.size(),
                                    [&](std::size_t) { ++candidates; });
    std::printf("  candidate facets per ray %.2f\n",
                static_cast<double>(candidates) / static_cast<double>(n));
  }
  // Photons that miss the primary (no hit above): the primary intersection alone.
  std::vector<std::size_t> misses;
  for (std::size_t i = 0; i < n; ++i)
    if (!hit[i])
      misses.push_back(i);
  if (!misses.empty()) {
    const double miss_ns = best_ns(misses.size(), [&] {
      for (const std::size_t i : misses)
        checksum += detail::intersect_mirror(model.primary, photons[i].ray, failed).has_value();
    });
    std::size_t miss_candidates = 0;
    if (!model.primary.facets.empty())
      for (const std::size_t i : misses)
        model.primary.grid.candidates(photons[i].ray, model.primary.facets.size(),
                                      [&](std::size_t) { ++miss_candidates; });
    std::printf("primary miss        %8.1f ns  (%zu rays, %.2f candidate facets per ray)\n",
                miss_ns, misses.size(),
                static_cast<double>(miss_candidates) / static_cast<double>(misses.size()));
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
  // Round-2 measures. Each one collects its inputs first and times only the calls.
  const auto focal_limit = [&](const Ray &r) -> std::optional<double> {
    if (model.focal_surface.asphere) {
      bool f = false;
      if (const auto h = model.focal_surface.asphere->intersect(
              r, kMinDistanceM, [](const Vec3 &) { return true; }, f))
        return h->t;
      return std::nullopt;
    }
    if (r.direction.z == 0.0)
      return std::nullopt;
    const double t = (model.focal_surface.plane_z - r.origin.z) / r.direction.z;
    if (t > kMinDistanceM)
      return t;
    return std::nullopt;
  };

  const double responses_ns = best_ns(n, [&] {
    for (std::size_t i = 0; i < n; ++i) {
      const auto r =
          trace_photon(model, photons[i].ray, 400.0, 0.0, 1.0, i, 3, true, false, true, classify);
      checksum += r.weight + r.time_ns;
    }
  });
  std::printf("trace_photon (responses) %8.1f ns\n", responses_ns);

  // Rays that hit the primary, with their hit records.
  std::vector<detail::MirrorHit> primary_hits;
  std::vector<std::size_t> primary_index;
  for (std::size_t i = 0; i < n; ++i)
    if (const auto h = detail::intersect_mirror(model.primary, photons[i].ray, failed)) {
      primary_hits.push_back(*h);
      primary_index.push_back(i);
    }
  if (!primary_hits.empty()) {
    const double reflect_ns = best_ns(primary_hits.size(), [&] {
      for (std::size_t k = 0; k < primary_hits.size(); ++k) {
        const std::size_t i = primary_index[k];
        const Vec3 d = detail::reflect_at(model.primary, primary_hits[k], photons[i].ray.direction,
                                          i, 3, Stream::primary_scatter);
        checksum += d.x;
      }
    });
    std::printf("reflect_at per hit  %8.1f ns  (%zu hits)\n", reflect_ns, primary_hits.size());
  }

  // Reflected rays after the primary, with the distance to their next surface.
  std::vector<Ray> after_primary_rays, after_secondary_rays;
#ifdef OBDEECT_BENCH_COUNTERS
  std::vector<Ray> reflected_primary_rays; // All of them, also those that miss the next surface.
#endif
  std::vector<double> after_primary_limits, after_secondary_limits;
  for (std::size_t k = 0; k < primary_hits.size(); ++k) {
    const std::size_t i = primary_index[k];
    const auto &mh = primary_hits[k];
    const Vec3 d1 = detail::reflect_at(model.primary, mh, photons[i].ray.direction, i, 3,
                                       Stream::primary_scatter);
    const Ray r1{mh.hit.point, d1};
#ifdef OBDEECT_BENCH_COUNTERS
    reflected_primary_rays.push_back(r1);
#endif
    std::optional<double> limit1;
    std::optional<detail::MirrorHit> m2;
    if (model.secondary) {
      m2 = detail::intersect_mirror(*model.secondary, r1, failed);
      if (m2)
        limit1 = m2->hit.t;
    } else {
      limit1 = focal_limit(r1);
    }
    if (!limit1)
      continue;
    after_primary_rays.push_back(r1);
    after_primary_limits.push_back(*limit1);
    if (!model.secondary)
      continue;
    const Vec3 d2 = detail::reflect_at(*model.secondary, *m2, d1, i, 3, Stream::secondary_scatter);
    const Ray r2{m2->hit.point, d2};
    if (const auto limit2 = focal_limit(r2)) {
      after_secondary_rays.push_back(r2);
      after_secondary_limits.push_back(*limit2);
    }
  }
  const auto time_obscurers = [&](const std::vector<Ray> &rays, const std::vector<double> &limits,
                                  Leg leg, const char *label) {
    if (rays.empty())
      return;
#ifdef OBDEECT_BENCH_COUNTERS
    Obscurers::counters() = {};
#endif
    const double ns = best_ns(rays.size(), [&] {
      for (std::size_t j = 0; j < rays.size(); ++j)
        if (const auto b = model.obscurers.nearest(rays[j], leg, limits[j]))
          checksum += b->t;
    });
    std::printf("obscurers %-16s %8.1f ns  (%zu rays)\n", label, ns, rays.size());
#ifdef OBDEECT_BENCH_COUNTERS
    std::printf("  %s\n", Obscurers::counter_report(3 * rays.size()).c_str());
#endif
  };
  time_obscurers(after_primary_rays, after_primary_limits, Leg::after_primary, "after primary");
  if (model.secondary)
    time_obscurers(after_secondary_rays, after_secondary_limits, Leg::after_secondary,
                   "after secondary");

  // Full trace cost of the photons that miss the primary (geometry only, as the
  // trace_photon line).
  std::vector<std::size_t> missed;
  for (std::size_t i = 0; i < n; ++i)
    if (trace_photon(model, photons[i].ray, 400.0, 0.0, 1.0, i, 3, false).status ==
        Status::missed_primary)
      missed.push_back(i);
  if (!missed.empty()) {
    const double miss_ns = best_ns(missed.size(), [&] {
      for (const std::size_t i : missed)
        checksum += trace_photon(model, photons[i].ray, 400.0, 0.0, 1.0, i, 3, false).path_length_m;
    });
    std::printf("miss cost           %8.1f ns  (%zu rays)\n", miss_ns, missed.size());
  }

#ifdef OBDEECT_BENCH_COUNTERS
  // Newton steps per call of each aspheric surface, for the rays that reach it:
  // photons (primary), all rays reflected by the primary (secondary, or the focal
  // surface without a secondary) and rays after the secondary (focal surface).
  const auto report_steps = [&](const char *label, const std::vector<Ray> &rays, auto intersect) {
    Asphere::counters() = {};
    for (const Ray &r : rays)
      checksum += intersect(r);
    const auto &c = Asphere::counters();
    if (c.calls)
      std::printf("  newton steps per call %-9s %.3f (%llu calls)\n", label,
                  static_cast<double>(c.newton_steps) / static_cast<double>(c.calls),
                  static_cast<unsigned long long>(c.calls));
  };
  const auto through_mirror = [&](const Mirror &m) {
    return [&](const Ray &r) {
      return detail::intersect_mirror(m, r, failed).has_value() ? 1.0 : 0.0;
    };
  };
  if (model.primary.asphere) {
    std::vector<Ray> rays;
    for (const auto &photon : photons)
      rays.push_back(photon.ray);
    report_steps("primary", rays, through_mirror(model.primary));
  }
  if (model.secondary && model.secondary->asphere)
    report_steps("secondary", reflected_primary_rays, through_mirror(*model.secondary));
  if (model.focal_surface.asphere)
    report_steps("focal", model.secondary ? after_secondary_rays : reflected_primary_rays,
                 [&](const Ray &r) {
                   return model.focal_surface.asphere
                                  ->intersect(
                                      r, kMinDistanceM, [](const Vec3 &) { return true; }, failed)
                                  .has_value()
                              ? 1.0
                              : 0.0;
                 });
#endif
  std::printf("checksum %.17g\n", checksum);
  return 0;
}
