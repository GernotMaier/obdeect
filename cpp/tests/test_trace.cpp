// End-to-end transport through analytic telescopes.
#include "check.hpp"

#include "obdeect/bunches.hpp"
#include "obdeect/source.hpp"
#include "obdeect/trace.hpp"

#include <array>
#include <atomic>
#include <cstring>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace obdeect;

namespace {

constexpr double kInfinity = std::numeric_limits<double>::infinity();

// One large spherical facet of focal length f; focal plane at z = f.
OpticalModel single_facet(double f, double half_width) {
  OpticalModel m;
  m.primary.facets.push_back({1,
                              {0, 0, 0},
                              {0, 0, 1},
                              {1, 0, 0},
                              {0, 1, 0},
                              {1, 0, 0},
                              {0, 1, 0},
                              half_width,
                              2 * f,
                              Aperture::circle});
  m.primary.grid = FacetGrid(m.primary.facets);
  m.focal_surface = {2, f, std::nullopt};
  return m;
}

// Paraboloid primary (f = 10, hole 1 m), flat secondary at z = 5 folding the
// focus to the vertex plane z = 0, where an asphere focal surface sits.
OpticalModel folded_paraboloid() {
  OpticalModel m;
  m.primary.asphere = Asphere{1, 0.0, 1.0, 6.0, 1.0, {0.0, 1.0 / 40.0}};
  Mirror secondary;
  secondary.asphere = Asphere{2, 5.0, 0.0, 3.0, 1.0, {0.0}};
  m.secondary = secondary;
  m.focal_surface = {3, 0.0, Asphere{3, 0.0, 0.0, 0.5, 1.0, {0.0}}};
  return m;
}

TraceResult shoot(const OpticalModel &m, Vec3 origin, Vec3 direction, std::uint64_t id = 0,
                  bool responses = true) {
  return trace_photon(m, {origin, normalized(direction)}, 400.0, 0.0, 1.0, id, 1, responses);
}

void spherical_mirror() {
  const auto m = single_facet(10.0, 0.5);
  // Paraxial on-axis rays focus at the axis (spherical aberration ~ h^2/(4R)).
  const auto r = shoot(m, {0.05, 0.0, 30.0}, {0, 0, -1});
  CHECK(r.status == Status::detected);
  CHECK(std::abs(r.position.x) < 1e-5 && std::abs(r.position.y) < 1e-12);
  CHECK_NEAR(r.position.z, 10.0, 1e-12);
  CHECK(r.primary_surface == 1 && r.terminal_surface == 2);
  CHECK_NEAR(r.primary_incidence_deg, std::asin(0.05 / 20.0) * 57.29577951308232, 1e-9);
  CHECK(shoot(m, {2.0, 0.0, 30.0}, {0, 0, -1}).status == Status::missed_primary);
  // Off-axis plane wave: the image centroid moves by f tan(theta).
  const double theta = 0.01;
  double x = 0.0;
  int n = 0;
  for (std::uint64_t i = 0; i < 2000; ++i) {
    const auto u = uniform4(2, i, Stream::source);
    const Vec3 start{0.4 * (u[0] - 0.5) + 30 * std::tan(theta), 0.4 * (u[1] - 0.5), 30.0};
    const auto hit = shoot(m, start, {-std::sin(theta), 0.0, -std::cos(theta)});
    if (hit.status == Status::detected) {
      x += hit.position.x;
      ++n;
    }
  }
  CHECK(n > 1500);
  CHECK_NEAR(x / n, -10.0 * std::tan(theta), 2e-5);
}

void folded_dual_mirror() {
  auto m = folded_paraboloid();
  // A paraboloid images an on-axis plane wave perfectly.
  for (const double r0 : {1.5, 3.0, 5.5}) {
    const auto r = shoot(m, {r0, 0.3, 30.0}, {0, 0, -1});
    CHECK(r.status == Status::detected);
    CHECK(std::hypot(r.position.x, r.position.y) < 1e-9);
    CHECK_NEAR(r.position.z, 0.0, 1e-12);
    CHECK_NEAR(r.secondary_point.z, 5.0, 1e-12);
    CHECK(r.terminal_surface == 3);
    // A point on the paraboloid lies sag + f from the focus, and the flat
    // fold preserves that distance: 30 - sag down, then sag + 10.
    CHECK_NEAR(r.path_length_m, 40.0, 1e-9);
  }
  CHECK(shoot(m, {0.5, 0.0, 30.0}, {0, 0, -1}).status ==
        Status::missed_primary); // Through the hole.
  // The back of the secondary blocks light from the sky when required.
  m.secondary_blocks_incoming = true;
  const auto blocked = shoot(m, {2.0, 0.0, 30.0}, {0, 0, -1});
  CHECK(blocked.status == Status::blocked && blocked.terminal_surface == 2);
  CHECK(shoot(m, {4.0, 0.0, 30.0}, {0, 0, -1}).status == Status::detected);
  CHECK(shoot(m, {0.2, 0.0, 30.0}, {0.05, 0, -1}).status == Status::blocked);
}

void obscurer_legs() {
  auto m = single_facet(10.0, 0.5);
  Obscurer rod;
  rod.kind = Obscurer::Kind::cylinder;
  rod.id = 9;
  rod.a = {-1, 0.2, 5};
  rod.b = {1, 0.2, 5};
  rod.radius_a = 0.05;
  rod.bound_a = rod.a;
  rod.bound_b = rod.b;
  rod.bound_radius = rod.radius_a;
  rod.legs = Leg::incoming;
  m.obscurers = Obscurers({rod}, {});
  // On the way in the ray passes the rod at y = 0.2; after reflection it
  // passes near the axis at z = 5, well clear of the rod.
  auto r = shoot(m, {0.0, 0.2, 30.0}, {0, 0, -1});
  CHECK(r.status == Status::blocked && r.terminal_surface == 9);
  CHECK_NEAR(r.position.z, 5.05, 1e-12);
  CHECK(r.weight == 0.0);
  // A ray at y = -0.4 reflects towards the focus and crosses y = -0.2 at
  // z = 5: a rod there absorbs it only on the after-primary leg.
  rod.a = {-1, -0.2, 5};
  rod.b = {1, -0.2, 5};
  rod.bound_a = rod.a;
  rod.bound_b = rod.b;
  rod.legs = Leg::after_primary;
  m.obscurers = Obscurers({rod}, {});
  r = shoot(m, {0.0, -0.4, 30.0}, {0, 0, -1});
  CHECK(r.status == Status::blocked); // Reflected ray passes y = -0.2 at z = 5.
  rod.legs = Leg::incoming;
  m.obscurers = Obscurers({rod}, {});
  CHECK(shoot(m, {0.0, -0.4, 30.0}, {0, 0, -1}).status == Status::detected);
}

// Photons from inside the telescope (a flasher) go straight to the focal
// surface: no mirror and no telescope transmission, only the obscurers of
// the leg to the focal surface.
void direct_photons() {
  auto m = single_facet(10.0, 0.5);
  m.telescope_transmission = TelescopeTransmission{0.5, 0.0, 1.0, 2.0, 1.0};
  m.primary.reflectivity = Response{300.0, 100.0, 3, {0.0}, {0.9, 0.8, 0.7}};
  const auto flash = [&m](double y) {
    return trace_photon(m, {{0.0, y, 0.5}, {0.0, 0.0, 1.0}}, 400.0, 1.0, 1.0, 0, 1, true, true);
  };
  auto r = flash(0.1);
  CHECK(r.status == Status::detected);
  CHECK_NEAR(r.position.z, 10.0, 1e-12);
  CHECK_NEAR(r.path_length_m, 9.5, 1e-12);
  CHECK_NEAR(r.weight, 1.0, 0.0);
  CHECK(r.primary_surface == kNoSurface);
  Obscurer rod;
  rod.kind = Obscurer::Kind::cylinder;
  rod.id = 9;
  rod.a = {-1, 0.1, 5};
  rod.b = {1, 0.1, 5};
  rod.radius_a = 0.05;
  rod.bound_a = rod.a;
  rod.bound_b = rod.b;
  rod.bound_radius = rod.radius_a;
  rod.legs = Leg::after_primary;
  m.obscurers = Obscurers({rod}, {});
  CHECK(flash(0.1).status == Status::blocked);
  CHECK(flash(0.3).status == Status::detected);
  rod.legs = Leg::incoming;
  m.obscurers = Obscurers({rod}, {});
  CHECK(flash(0.1).status == Status::detected);
  // Dual mirror: from between the mirrors down to the focal surface at z = 0.
  r = trace_photon(folded_paraboloid(), {{0.1, 0.0, 2.0}, {0.0, 0.0, -1.0}}, 400.0, 0.0, 1.0, 0, 1,
                   true, true);
  CHECK(r.status == Status::detected && r.terminal_surface == 3);
  CHECK_NEAR(r.position.z, 0.0, 1e-9);
}

// Without loss classification a photon that misses the primary is reported as
// missed_primary, also where it would be blocked on the incoming leg first.
// Every other photon, and every detected one, is the same with either setting.
void classify_losses_option() {
  auto m = single_facet(10.0, 0.5);
  Obscurer rod;
  rod.kind = Obscurer::Kind::cylinder;
  rod.id = 9;
  rod.a = {-1, 0.2, 5};
  rod.b = {1, 0.2, 5};
  rod.radius_a = 0.05;
  rod.bound_a = rod.a;
  rod.bound_b = rod.b;
  rod.bound_radius = rod.radius_a;
  rod.legs = Leg::incoming;
  m.obscurers = Obscurers({rod}, {});

  const auto same_double = [](double a, double b) {
    return a == b || (std::isnan(a) && std::isnan(b));
  };
  const auto same_vec = [&](const Vec3 &a, const Vec3 &b) {
    return same_double(a.x, b.x) && same_double(a.y, b.y) && same_double(a.z, b.z);
  };
  const auto same = [&](const TraceResult &a, const TraceResult &b) {
    return a.status == b.status && a.terminal_surface == b.terminal_surface &&
           a.primary_surface == b.primary_surface && a.secondary_surface == b.secondary_surface &&
           same_vec(a.position, b.position) && same_vec(a.direction, b.direction) &&
           same_vec(a.primary_point, b.primary_point) &&
           same_vec(a.secondary_point, b.secondary_point) &&
           same_double(a.path_length_m, b.path_length_m) && same_double(a.time_ns, b.time_ns) &&
           same_double(a.weight, b.weight) &&
           same_double(a.primary_incidence_deg, b.primary_incidence_deg) &&
           same_double(a.secondary_incidence_deg, b.secondary_incidence_deg) &&
           same_double(a.focal_incidence_deg, b.focal_incidence_deg);
  };

  // Photons across the dish and the rod (|x| < 1 is the rod, |x| > 0.5 misses the dish).
  std::vector<Ray> rays;
  for (int ix = 0; ix <= 48; ++ix)
    for (int iy = 0; iy <= 8; ++iy)
      for (const double tilt : {0.0, 0.03})
        rays.push_back(
            {{-1.2 + 0.05 * ix, -0.4 + 0.1 * iy, 30.0}, normalized(Vec3{tilt, 0.0, -1.0})});
  std::size_t detected = 0, differing = 0, all_differ_as_expected = 0;
  std::vector<Status> fast_status(rays.size());
  for (std::size_t i = 0; i < rays.size(); ++i) {
    const auto full = trace_photon(m, rays[i], 400.0, 0.0, 1.0, i, 1, false, false, true, true);
    const auto fast = trace_photon(m, rays[i], 400.0, 0.0, 1.0, i, 1, false, false, true, false);
    fast_status[i] = fast.status;
    if (full.status == fast.status) {
      CHECK(same(full, fast));
    } else {
      ++differing;
      if (full.status == Status::blocked && fast.status == Status::missed_primary)
        ++all_differ_as_expected;
    }
    detected += full.status == Status::detected;
  }
  CHECK(differing > 0 && detected > 0 && rays.size() == 2 * 49 * 9);
  CHECK(all_differ_as_expected == differing);

  // The block tracer passes the option on, with one and with a pool.
  const std::size_t n = rays.size();
  std::vector<double> position(3 * n), direction(3 * n), wavelength(n, 400.0), time(n, 0.0),
      weight(n, 1.0);
  std::vector<std::uint64_t> id(n);
  for (std::size_t i = 0; i < n; ++i) {
    id[i] = i;
    position[3 * i] = rays[i].origin.x;
    position[3 * i + 1] = rays[i].origin.y;
    position[3 * i + 2] = rays[i].origin.z;
    direction[3 * i] = rays[i].direction.x;
    direction[3 * i + 1] = rays[i].direction.y;
    direction[3 * i + 2] = rays[i].direction.z;
  }
  const PhotonBlock in{position.data(),
                       direction.data(),
                       wavelength.data(),
                       time.data(),
                       weight.data(),
                       id.data(),
                       n};
  std::vector<std::uint8_t> status(n);
  ResultBlock out{};
  out.status = status.data();
  const auto check_statuses = [&] {
    for (std::size_t i = 0; i < n; ++i)
      CHECK(status[i] == static_cast<std::uint8_t>(fast_status[i]));
  };
  CHECK(trace_block(m, in, out, 1, false, 2, false, false, false) == n);
  check_statuses();
  ThreadPool pool(2);
  CHECK(trace_block(m, in, out, 1, false, pool, false, false, false) == n);
  check_statuses();
}

void responses_and_timing() {
  auto m = single_facet(10.0, 0.5);
  m.primary.reflectivity = Response{300.0, 100.0, 3, {0.0}, {0.9, 0.8, 0.7}};
  m.telescope_transmission = TelescopeTransmission{0.5, 0.0, 1.0, 2.0, 1.0};
  m.propagation_group_index = 1.0003;
  const auto r = shoot(m, {0.1, 0.0, 30.0}, {0, 0, -1});
  CHECK_NEAR(r.weight, 0.5 * 0.8, 1e-12); // 400 nm on the 300/400/500 nm grid.
  CHECK_NEAR(r.time_ns, r.path_length_m * 1.0003 / kSpeedOfLightMPerNs, 1e-9);
  CHECK_NEAR(shoot(m, {0.1, 0.0, 30.0}, {0, 0, -1}, 0, false).weight, 1.0, 0.0);
  CHECK(trace_photon(m, {{0, 0, 30}, {0, 0, -1}}, 0.0, 0.0, 1.0, 0, 1).status ==
        Status::invalid_input);
  CHECK(trace_photon(m, {{0, 0, 30}, {0, 0, -1}}, kInfinity, 0.0, 1.0, 0, 1).status ==
        Status::invalid_input);
  CHECK(trace_photon(m, {{0, 0, 30}, {0, 0, -1}}, 400.0, 0.0, kInfinity, 0, 1).status ==
        Status::invalid_input);
  CHECK(trace_photon(m, {{0, 0, 30}, {0, 0, -1.1}}, 400.0, 0.0, 1.0, 0, 1).status ==
        Status::invalid_input);
}

// RMS deviation of reflected rays from a flat mirror per axis equals sigma.
void scatter_widths() {
  for (const bool slopes : {false, true}) {
    OpticalModel m;
    if (slopes) {
      m.primary.asphere = Asphere{1, 0.0, 0.0, 10.0, 1.0, {0.0}};
    } else {
      m.primary.facets.push_back({1,
                                  {0, 0, 0},
                                  {0, 0, 1},
                                  {1, 0, 0},
                                  {0, 1, 0},
                                  {1, 0, 0},
                                  {0, 1, 0},
                                  10.0,
                                  0.0,
                                  Aperture::circle});
      m.primary.grid = FacetGrid(m.primary.facets);
    }
    m.primary.scatter = Scatter{0.002, 0.0, 0.0, slopes};
    m.focal_surface = {2, 10.0, std::nullopt};
    double sxx = 0.0;
    const int n = 50000;
    for (int i = 0; i < n; ++i) {
      const auto r = shoot(m, {0.0, 0.0, 30.0}, {0, 0, -1}, static_cast<std::uint64_t>(i));
      sxx += (r.direction.x / r.direction.z) * (r.direction.x / r.direction.z);
    }
    CHECK_NEAR(std::sqrt(sxx / n), 0.002, 0.00004);
  }
}

// The sub-photon index selects other scatter draws; 0 is the default.
void sub_photon_index() {
  for (const bool dual : {false, true}) {
    auto m = dual ? folded_paraboloid() : single_facet(10.0, 0.5);
    m.primary.scatter = Scatter{0.002, 0.0, 0.0, dual};
    if (dual)
      m.secondary->scatter = Scatter{0.002, 0.0, 0.0, true};
    const Ray ray{{dual ? 3.0 : 0.05, 0.1, 30.0}, {0, 0, -1}};
    int differ = 0;
    for (std::uint64_t id = 0; id < 100; ++id) {
      const auto a = trace_photon(m, ray, 400.0, 0.0, 1.0, id, 7);
      const auto b = trace_photon(m, ray, 400.0, 0.0, 1.0, id, 7, true, false, true, true, 0);
      const auto c = trace_photon(m, ray, 400.0, 0.0, 1.0, id, 7, true, false, true, true, 1);
      CHECK(a.status == Status::detected && c.status == Status::detected);
      CHECK(std::memcmp(&a.position, &b.position, sizeof(Vec3)) == 0);
      CHECK(std::memcmp(&a.direction, &b.direction, sizeof(Vec3)) == 0);
      if (std::memcmp(&a.position, &c.position, sizeof(Vec3)) != 0)
        ++differ;
    }
    CHECK(differ == 100);
  }
}

void block_invariance() {
  auto m = single_facet(10.0, 0.5);
  m.primary.scatter = Scatter{0.001, 0.3, 0.004, false};
  const std::size_t n = 50000;
  std::vector<double> position(3 * n), direction(3 * n), wavelength(n, 400.0), time(n, 0.0),
      weight(n, 1.0);
  std::vector<std::uint64_t> id(n);
  const PointSource star{{0, 0, 1}, kInfinity, 0.6, 0.0, 30.0, 400.0};
  for (std::size_t i = 0; i < n; ++i) {
    id[i] = 1000 + i;
    const auto p = sample_source(star, id[i], 4);
    position[3 * i] = p.ray.origin.x;
    position[3 * i + 1] = p.ray.origin.y;
    position[3 * i + 2] = p.ray.origin.z;
    direction[3 * i] = p.ray.direction.x;
    direction[3 * i + 1] = p.ray.direction.y;
    direction[3 * i + 2] = p.ray.direction.z;
  }
  struct Out {
    std::vector<std::uint8_t> status;
    std::vector<std::uint32_t> terminal, m1, m2;
    std::vector<double> pos, dir, path, t, w, p1, p2, a1, a2, af;
    explicit Out(std::size_t k)
        : status(k), terminal(k), m1(k), m2(k), pos(3 * k), dir(3 * k), path(k), t(k), w(k),
          p1(3 * k), p2(3 * k), a1(k), a2(k), af(k) {}
    ResultBlock block() {
      return {status.data(), terminal.data(), m1.data(), m2.data(), pos.data(),
              dir.data(),    path.data(),     t.data(),  w.data(),  p1.data(),
              p2.data(),     a1.data(),       a2.data(), af.data()};
    }
  };
  const PhotonBlock in{position.data(),
                       direction.data(),
                       wavelength.data(),
                       time.data(),
                       weight.data(),
                       id.data(),
                       n};
  Out one(n);
  trace_block(m, in, one.block(), 4, true, 1);
  // Threads take chunks dynamically; every field must still agree bit for bit.
  for (const unsigned threads : {3U, 4U, 8U}) {
    Out many(n);
    trace_block(m, in, many.block(), 4, true, threads);
    const auto same_bits = [](const auto &x, const auto &y) {
      return x.size() == y.size() &&
             std::memcmp(x.data(), y.data(), x.size() * sizeof(x[0])) == 0; // NaN-safe.
    };
    CHECK(same_bits(one.status, many.status) && same_bits(one.terminal, many.terminal) &&
          same_bits(one.m1, many.m1) && same_bits(one.m2, many.m2) &&
          same_bits(one.pos, many.pos) && same_bits(one.dir, many.dir) &&
          same_bits(one.path, many.path) && same_bits(one.t, many.t) && same_bits(one.w, many.w) &&
          same_bits(one.p1, many.p1) && same_bits(one.p2, many.p2) && same_bits(one.a1, many.a1) &&
          same_bits(one.a2, many.a2) && same_bits(one.af, many.af));
    // The persistent pool gives the same bits as the threads.
    ThreadPool pool(threads);
    Out pooled(n);
    trace_block(m, in, pooled.block(), 4, true, pool);
    CHECK(same_bits(one.status, pooled.status) && same_bits(one.terminal, pooled.terminal) &&
          same_bits(one.m1, pooled.m1) && same_bits(one.m2, pooled.m2) &&
          same_bits(one.pos, pooled.pos) && same_bits(one.dir, pooled.dir) &&
          same_bits(one.path, pooled.path) && same_bits(one.t, pooled.t) &&
          same_bits(one.w, pooled.w) && same_bits(one.p1, pooled.p1) &&
          same_bits(one.p2, pooled.p2) && same_bits(one.a1, pooled.a1) &&
          same_bits(one.a2, pooled.a2) && same_bits(one.af, pooled.af));
  }
  // Detected-only output with a few columns: the detected rows of the full trace, in order.
  std::size_t expected = 0;
  for (std::size_t i = 0; i < n; ++i)
    expected += one.status[i] == static_cast<std::uint8_t>(Status::detected);
  CHECK(expected > 100 && expected < n);
  for (const unsigned threads : {1U, 3U, 8U}) {
    std::vector<std::uint8_t> status(n);
    std::vector<double> pos(3 * n);
    std::vector<std::uint64_t> index(n);
    ResultBlock partial;
    partial.status = status.data();
    partial.position = pos.data();
    partial.index = index.data();
    CHECK(trace_block(m, in, partial, 4, true, threads, false, true) == expected);
    std::size_t row = 0;
    bool same_rows = true;
    for (std::size_t i = 0; i < n; ++i)
      if (one.status[i] == static_cast<std::uint8_t>(Status::detected)) {
        same_rows = same_rows && index[row] == i && status[row] == one.status[i] &&
                    std::memcmp(&pos[3 * row], &one.pos[3 * i], 3 * sizeof(double)) == 0;
        ++row;
      }
    CHECK(same_rows);
    // The same on a pool: equal rows.
    ThreadPool pool(threads);
    std::vector<std::uint8_t> pooled_status(n);
    std::vector<double> pooled_pos(3 * n);
    std::vector<std::uint64_t> pooled_index(n);
    ResultBlock pooled;
    pooled.status = pooled_status.data();
    pooled.position = pooled_pos.data();
    pooled.index = pooled_index.data();
    CHECK(trace_block(m, in, pooled, 4, true, pool, false, true) == expected);
    CHECK(std::memcmp(pooled_status.data(), status.data(), expected) == 0 &&
          std::memcmp(pooled_pos.data(), pos.data(), 3 * expected * sizeof(double)) == 0 &&
          std::memcmp(pooled_index.data(), index.data(), expected * sizeof(std::uint64_t)) == 0);
  }
  // Without incidence-angle outputs the angles are not computed, yet an angle-dependent
  // reflectivity still sees them: the weights are the same.
  auto angular = m;
  angular.primary.reflectivity =
      Response{300.0, 100.0, 3, {0.0, 60.0}, {0.9, 0.5, 0.8, 0.4, 0.7, 0.3}};
  Out with_angles(n);
  trace_block(angular, in, with_angles.block(), 4, true, 2);
  std::vector<double> only_weight(n);
  ResultBlock weight_only;
  weight_only.weight = only_weight.data();
  trace_block(angular, in, weight_only, 4, true, 2);
  CHECK(only_weight == with_angles.w);
  CHECK(with_angles.w != one.w);
  // Reversed order gives the same result per photon id.
  std::vector<std::size_t> order(n);
  std::iota(order.rbegin(), order.rend(), std::size_t{0});
  std::vector<double> rp(3 * n), rd(3 * n);
  std::vector<std::uint64_t> rid(n);
  for (std::size_t k = 0; k < n; ++k) {
    for (int c = 0; c < 3; ++c) {
      rp[3 * k + c] = position[3 * order[k] + c];
      rd[3 * k + c] = direction[3 * order[k] + c];
    }
    rid[k] = id[order[k]];
  }
  Out reversed(n);
  trace_block(m,
              {rp.data(), rd.data(), wavelength.data(), time.data(), weight.data(), rid.data(), n},
              reversed.block(), 4, true, 3);
  bool same = true;
  for (std::size_t k = 0; k < n; ++k)
    for (int c = 0; c < 3; ++c)
      same = same && reversed.pos[3 * k + c] == one.pos[3 * order[k] + c];
  CHECK(same);
  // Several blocks and models in one job equal one trace_block call each, with
  // sizes that are no multiple of the chunk size, an empty block and detected_only.
  struct Case {
    const OpticalModel *model;
    std::size_t size;
    bool responses, direct, detected_only;
  };
  const std::array<Case, 5> cases{
      Case{&m, n, true, false, false}, Case{&angular, 10001, true, false, true},
      Case{&m, 0, true, false, false}, Case{&angular, 777, false, false, false},
      Case{&m, 4096, true, false, true}};
  std::vector<std::unique_ptr<Out>> separate, together;
  std::vector<BlockRequest> requests;
  std::vector<std::size_t> expected_rows;
  for (const Case &c : cases) {
    PhotonBlock part = in;
    part.size = c.size;
    separate.push_back(std::make_unique<Out>(c.size));
    together.push_back(std::make_unique<Out>(c.size));
    expected_rows.push_back(trace_block(*c.model, part, separate.back()->block(), 4, c.responses,
                                        2U, c.direct, c.detected_only));
    requests.push_back(
        {c.model, part, together.back()->block(), c.responses, c.direct, c.detected_only});
  }
  for (const unsigned threads : {1U, 2U, 8U}) {
    ThreadPool pool(threads);
    for (auto &t : together)
      t = std::make_unique<Out>(t->status.size());
    for (std::size_t b = 0; b < cases.size(); ++b)
      requests[b].out = together[b]->block();
    CHECK(trace_blocks(requests, 4, pool) == expected_rows);
    for (std::size_t b = 0; b < cases.size(); ++b) {
      const std::size_t rows = expected_rows[b];
      CHECK(std::memcmp(separate[b]->status.data(), together[b]->status.data(), rows) == 0);
      CHECK(std::memcmp(separate[b]->pos.data(), together[b]->pos.data(),
                        3 * rows * sizeof(double)) == 0);
      CHECK(std::memcmp(separate[b]->w.data(), together[b]->w.data(), rows * sizeof(double)) == 0);
    }
  }
  CHECK(expected_rows[1] > 0 && expected_rows[1] < 10001 && expected_rows[2] == 0);
}

// The pool runs every task exactly once, job after job, and shuts down cleanly.
void thread_pool() {
  constexpr std::size_t tasks = 1'000'000;
  std::vector<std::atomic<std::uint8_t>> seen(tasks);
  for (const unsigned threads : {1U, 2U, 8U}) {
    ThreadPool pool(threads);
    CHECK(pool.size() == threads);
    for (auto &s : seen)
      s.store(0, std::memory_order_relaxed);
    for (int job = 0; job < 100; ++job) // 100 consecutive jobs of 10^4 tasks: 10^6 in all.
      pool.run(tasks / 100, [&](std::size_t c) {
        seen[static_cast<std::size_t>(job) * (tasks / 100) + c].fetch_add(1);
      });
    bool once = true;
    for (const auto &s : seen)
      once = once && s.load() == 1;
    CHECK(once);
    pool.run(0, [](std::size_t) {}); // Empty job.
  }
  for (int i = 0; i < 50; ++i) { // Construct and destroy, with and without work.
    ThreadPool pool(1 + static_cast<unsigned>(i % 8));
    if (i % 2) {
      std::atomic<std::size_t> sum{0};
      pool.run(64, [&](std::size_t c) { sum.fetch_add(c); });
      CHECK(sum.load() == 64 * 63 / 2);
    }
  }
}

// An exception in a task reaches the caller; the pool stays usable.
void thread_pool_exceptions() {
  for (const unsigned threads : {1U, 4U}) {
    ThreadPool pool(threads);
    for (int round = 0; round < 3; ++round) {
      bool caught = false;
      try {
        pool.run(1000, [](std::size_t c) {
          if (c == 500)
            throw std::runtime_error("task failed");
        });
      } catch (const std::runtime_error &e) {
        caught = std::string(e.what()) == "task failed";
      }
      CHECK(caught);
      std::atomic<std::size_t> count{0}, sum{0};
      pool.run(1000, [&](std::size_t c) {
        count.fetch_add(1);
        sum.fetch_add(c);
      });
      CHECK(count.load() == 1000 && sum.load() == 1000 * 999 / 2);
    }
  }
}

template <class F> bool throws_invalid_argument(F &&f) {
  try {
    f();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}

// trace_bunches rejects a wavelength density reaching 0 nm and bad acceptance axes.
void bunch_option_validation() {
  const auto model = single_facet(10.0, 0.5);
  const float zero = 0.0F, one = 1.0F;
  const std::uint64_t id = 0;
  BunchBlock in;
  in.x = in.y = in.cx = in.cy = in.time_ns = in.zem = &zero;
  in.size = &one;
  in.wavelength_nm = &zero;
  in.bunch_id = &id;
  in.n = 1;
  BunchFrame frame;
  frame.plane_z = 20.0;
  ThreadPool pool(1);
  const auto allocate = [](std::size_t) { return BunchResult{}; };
  const auto run = [&](const BunchOptions &options) {
    trace_bunches(model, in, frame, options, 1, pool, allocate);
  };

  const LinearDensity from_zero({0.0, 500.0}, {1.0, 1.0}), from_below({-1.0, 500.0}, {1.0, 1.0});
  BunchOptions options;
  options.wavelength_density = &from_zero;
  CHECK(throws_invalid_argument([&] { run(options); }));
  options.wavelength_density = &from_below;
  CHECK(throws_invalid_argument([&] { run(options); }));

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::vector<double> ok{100.0, 1000.0};
  for (const auto &[x, y] : {std::pair{std::vector<double>{100.0, 100.0}, ok},
                             std::pair{std::vector<double>{1000.0, 100.0}, ok},
                             std::pair{ok, std::vector<double>{0.0, nan}},
                             std::pair{std::vector<double>{100.0, kInfinity}, ok}}) {
    const Map2D acceptance{x, y, std::vector<double>(4, 0.5), false};
    BunchOptions bad;
    bad.acceptance = &acceptance;
    CHECK(throws_invalid_argument([&] { run(bad); }));
  }

  // An unsupported telescope transmission is rejected with ModelError before any bunch is traced.
  auto unsupported = model;
  unsupported.telescope_transmission = TelescopeTransmission{1.0, -1.0, 0.1, 2.0, 1.0};
  bool model_error = false;
  try {
    trace_bunches(unsupported, in, frame, BunchOptions{}, 1, pool, allocate);
  } catch (const ModelError &) {
    model_error = true;
  }
  CHECK(model_error);
}

} // namespace

void all_tests() {
  spherical_mirror();
  folded_dual_mirror();
  obscurer_legs();
  classify_losses_option();
  direct_photons();
  responses_and_timing();
  scatter_widths();
  sub_photon_index();
  block_invariance();
  thread_pool();
  thread_pool_exceptions();
  bunch_option_validation();
}

int main() { return run_checks(all_tests); }
