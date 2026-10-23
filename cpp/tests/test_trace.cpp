// End-to-end transport through analytic telescopes.
#include "check.hpp"

#include "obdeect/source.hpp"
#include "obdeect/trace.hpp"

#include <numeric>
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

void block_invariance() {
  auto m = single_facet(10.0, 0.5);
  m.primary.scatter = Scatter{0.001, 0.3, 0.004, false};
  const std::size_t n = 20000;
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
  Out one(n), many(n);
  trace_block(m, in, one.block(), 4, true, 1);
  trace_block(m, in, many.block(), 4, true, 4);
  CHECK(one.pos == many.pos && one.dir == many.dir && one.status == many.status);
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
}

} // namespace

void all_tests() {
  spherical_mirror();
  folded_dual_mirror();
  obscurer_legs();
  direct_photons();
  responses_and_timing();
  scatter_widths();
  block_invariance();
}

int main() { return run_checks(all_tests); }
