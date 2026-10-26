// The response bound of an optical model against traced weights.
#include "check.hpp"

#include "obdeect/bounds.hpp"
#include "obdeect/trace.hpp"

#include <algorithm>
#include <vector>

using namespace obdeect;

namespace {

// One spherical facet of focal length f; focal plane at z = f.
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

// Paraboloid primary, flat secondary folding the focus to z = 0.
OpticalModel folded_paraboloid() {
  OpticalModel m;
  m.primary.asphere = Asphere{1, 0.0, 1.0, 6.0, 1.0, {0.0, 1.0 / 40.0}};
  Mirror secondary;
  secondary.asphere = Asphere{2, 5.0, 0.0, 3.0, 1.0, {0.0}};
  m.secondary = secondary;
  m.focal_surface = {3, 0.0, Asphere{3, 0.0, 0.0, 0.5, 1.0, {0.0}}};
  return m;
}

// Wavelengths 250, 400, 550, 700 nm; incidence 0, p, 2p, 4p degrees. The
// response is flat up to p degrees and falls beyond, so photons near normal
// incidence reach the bound.
Response angle_response(double p) {
  const auto row = [](double v) { return std::vector<double>{v, v, 0.6 * v, 0.3 * v}; };
  Response r{250.0, 150.0, 4, {0.0, p, 2 * p, 4 * p}, {}};
  for (const double v : {0.90, 0.95, 0.85, 0.70})
    r.values.insert(r.values.end(), row(v).begin(), row(v).end());
  return r;
}

// Map over [-half, half]^2: `hi` in the inner part, `lo` on the outer ring.
Map2D map(double lo, double hi, double half) {
  Map2D d;
  d.x_m = {-half, -0.5 * half, 0.5 * half, half};
  d.y_m = d.x_m;
  for (std::size_t i = 0; i < 4; ++i)
    for (std::size_t j = 0; j < 4; ++j)
      d.values.push_back(i == 0 || i == 3 || j == 0 || j == 3 ? lo : hi);
  d.zero_outside = true;
  return d;
}

// Photons over the dish (tilts up to 0.03 rad) at several wavelengths: no
// weight exceeds the bound, and some come close to it.
void check_model(const OpticalModel &m, double radius, std::uint64_t seed) {
  constexpr int n = 100000;
  for (const double lambda : {250.0, 300.0, 333.3, 400.0, 550.0, 700.0}) {
    const double bound = response_bound(m, lambda);
    double largest = 0.0;
    for (int i = 0; i < n; ++i) {
      const auto u = uniform4(seed, static_cast<std::uint64_t>(i), Stream::source);
      const double r = radius * std::sqrt(u[0]), phi = 6.283185307179586 * u[1];
      // Every fifth photon on axis, to reach the transmission maximum.
      const double tilt = i % 5 == 0 ? 0.0 : 0.03 * u[2];
      const double psi = 6.283185307179586 * u[3];
      const Vec3 direction = normalized({tilt * std::cos(psi), tilt * std::sin(psi), -1.0});
      const auto t = trace_photon(m, {{r * std::cos(phi), r * std::sin(phi), 30.0}, direction},
                                  lambda, 0.0, 1.0, static_cast<std::uint64_t>(i), seed);
      if (t.status != Status::detected)
        continue;
      CHECK(t.weight <= bound * 1.0);
      largest = std::max(largest, t.weight);
    }
    if (!(largest >= 0.9 * bound))
      std::fprintf(stderr, "bound loose at %g nm: largest weight %g, bound %g\n", lambda, largest,
                   bound);
    CHECK(largest >= 0.9 * bound);
  }
}

void bounds_single_mirror() {
  auto m = single_facet(10.0, 0.5);
  m.primary.reflectivity = angle_response(2.0);
  check_model(m, 0.45, 11);

  m.primary.degradation = map(0.8, 1.0, 0.6);
  check_model(m, 0.45, 12);

  m.telescope_transmission = TelescopeTransmission{0.9, 3.0, 0.02, 2.0, 1.5};
  check_model(m, 0.45, 13);
}

void bounds_dual_mirror() {
  auto m = folded_paraboloid();
  m.primary.reflectivity = angle_response(5.0);
  m.secondary->reflectivity = Response{250.0, 150.0, 4, {0.0}, {0.9, 0.95, 0.92, 0.7}};
  check_model(m, 5.5, 21);

  m.secondary->degradation = map(0.85, 1.0, 3.0);
  m.telescope_transmission = TelescopeTransmission{0.95, 1.0, 0.03, 2.0, 1.0};
  check_model(m, 5.5, 22);
}

void bounds_exact_and_errors() {
  // No responses at all: the bound is 1 (within the rounding margin).
  const auto plain = single_facet(10.0, 0.5);
  CHECK_NEAR(response_bound(plain, 400.0), 1.0, 1e-11);
  // Wavelength-only table: exact interpolated value, edge values outside.
  auto m = single_facet(10.0, 0.5);
  m.primary.reflectivity = Response{300.0, 100.0, 3, {0.0}, {0.9, 0.8, 0.7}};
  CHECK_NEAR(response_bound(m, 350.0), 0.85, 1e-11);
  CHECK_NEAR(response_bound(m, 100.0), 0.9, 1e-11);
  CHECK_NEAR(response_bound(m, 900.0), 0.7, 1e-11);
  // A transmission whose maximum is unknown is rejected.
  for (const auto &t : {TelescopeTransmission{1.0, -1.0, 0.1, 2.0, 1.0},
                        TelescopeTransmission{1.0, 1.0, 0.0, 2.0, 1.0},
                        TelescopeTransmission{1.0, 1.0, 0.1, 0.0, 1.0},
                        TelescopeTransmission{1.0, 1.0, 0.1, 2.0, -1.0}}) {
    m.telescope_transmission = t;
    bool threw = false;
    try {
      (void)response_bound(m, 400.0);
    } catch (const ModelError &) {
      threw = true;
    }
    CHECK(threw);
  }
}

} // namespace

void all_tests() {
  bounds_exact_and_errors();
  bounds_single_mirror();
  bounds_dual_mirror();
}

int main() { return run_checks(all_tests); }
