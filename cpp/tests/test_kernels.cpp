// Analytic tests of the random generator, tables, surfaces and obscurers.
#include "check.hpp"

#include "obdeect/model.hpp"
#include "obdeect/random.hpp"
#include "obdeect/source.hpp"

#include <optional>
#include <stdexcept>
#include <vector>

using namespace obdeect;

namespace {

void random_numbers() {
  // Random123 known-answer tests for Philox-4x32-10.
  CHECK((philox4x32({0, 0, 0, 0}, {0, 0}) ==
         std::array<std::uint32_t, 4>{0x6627e8d5U, 0xe169c58dU, 0xbc57ac4cU, 0x9b00dbd8U}));
  CHECK((philox4x32({0xffffffffU, 0xffffffffU, 0xffffffffU, 0xffffffffU},
                    {0xffffffffU, 0xffffffffU}) ==
         std::array<std::uint32_t, 4>{0x408f276dU, 0x41c83b0eU, 0xa20bc7c6U, 0x6d5451fdU}));
  CHECK((philox4x32({0x243f6a88U, 0x85a308d3U, 0x13198a2eU, 0x03707344U},
                    {0xa4093822U, 0x299f31d0U}) ==
         std::array<std::uint32_t, 4>{0xd16cfe09U, 0x94fdccebU, 0x5001e420U, 0x24126ea1U}));
  // Draws are reproducible, in (0, 1), and differ between streams and photons.
  const auto a = uniform4(7, 11, Stream::primary_scatter);
  CHECK(a == uniform4(7, 11, Stream::primary_scatter));
  CHECK(a != uniform4(7, 11, Stream::secondary_scatter));
  CHECK(a != uniform4(7, 12, Stream::primary_scatter));
  CHECK(a != uniform4(8, 11, Stream::primary_scatter));
  double mean = 0.0;
  for (std::uint64_t i = 0; i < 100000; ++i)
    for (const double u : uniform4(1, i, Stream::source)) {
      CHECK(u > 0.0 && u < 1.0);
      mean += u / 400000.0;
    }
  CHECK_NEAR(mean, 0.5, 0.002);
}

void tables() {
  Response r{400.0, 100.0, 3, {0.0}, {0.2, 0.4, 0.8}};
  CHECK_NEAR(r.at(450.0, 10.0), 0.3, 1e-15);
  CHECK_NEAR(r.at(100.0, 0.0), 0.2, 1e-15); // Clamped below.
  CHECK_NEAR(r.at(900.0, 0.0), 0.8, 1e-15); // Clamped above.
  Response angular{400.0, 100.0, 2, {0.0, 60.0}, {1.0, 0.5, 0.8, 0.2}};
  CHECK_NEAR(angular.at(450.0, 30.0), 0.5 * (0.75 + 0.5), 1e-15);
  CHECK_NEAR(angular.at(400.0, 90.0), 0.5, 1e-15);
  Map2D m{{0.0, 1.0}, {0.0, 2.0}, {0.0, 1.0, 0.5, 1.0}, false};
  CHECK_NEAR(m.at(0.5, 1.0), 0.625, 1e-15);
  CHECK_NEAR(m.at(-3.0, 0.0), 0.0, 1e-15);
  m.zero_outside = true;
  CHECK_NEAR(m.at(1.5, 1.0), 0.0, 1e-15);
}

void apertures() {
  CHECK(inside_aperture(Aperture::circle, 1.0, 0.7, 0.7));
  CHECK(!inside_aperture(Aperture::circle, 1.0, 0.8, 0.7));
  CHECK(inside_aperture(Aperture::square, 1.0, 0.99, -0.99));
  // A flat-x hexagon has flats at u = +-half and corners on the v axis.
  CHECK(inside_aperture(Aperture::hexagon_flat_x, 1.0, 0.99, 0.0));
  CHECK(!inside_aperture(Aperture::hexagon_flat_x, 1.0, 1.01, 0.0));
  CHECK(inside_aperture(Aperture::hexagon_flat_x, 1.0, 0.0, 1.15));
  CHECK(!inside_aperture(Aperture::hexagon_flat_x, 1.0, 0.0, 1.16));
  CHECK(inside_aperture(Aperture::hexagon_flat_y, 1.0, 1.15, 0.0));
  CHECK(!inside_aperture(Aperture::hexagon_flat_y, 1.0, 0.0, 1.01));
  CHECK_NEAR(circumradius(Aperture::hexagon_flat_x, 1.0), 2.0 / std::sqrt(3.0), 1e-15);
}

Facet facet(std::uint32_t id, Vec3 centre, double radius, double half_width, Aperture shape) {
  return {id,        centre,    {0, 0, 1},  {1, 0, 0}, {0, 1, 0},
          {1, 0, 0}, {0, 1, 0}, half_width, radius,    shape};
}

void facets() {
  const Facet f = facet(3, {0, 0, 0}, 20.0, 0.5, Aperture::hexagon_flat_y);
  // A paraxial ray parallel to the axis reflects through the focus R/2.
  const Ray ray{{0.01, 0.0, 30.0}, {0, 0, -1}};
  const auto hit = f.intersect(ray, kMinDistanceM);
  CHECK(hit && hit->id == 3);
  const Vec3 r = reflect(ray.direction, hit->normal);
  const double t = -hit->point.x / r.x;
  CHECK_NEAR(hit->point.z + t * r.z, 10.0, 1e-5);
  CHECK(!f.intersect({{0.0, 0.6, 30.0}, {0, 0, -1}}, kMinDistanceM)); // Outside the hexagon.
  CHECK(!f.intersect({{0.0, 0.0, -5.0}, {0, 0, 1}}, kMinDistanceM));  // Back face.
  const Facet flat = facet(4, {0, 0, 1}, 0.0, 0.5, Aperture::circle);
  const auto flat_hit = flat.intersect({{0.1, 0.1, 5.0}, normalized({0.01, 0, -1})}, kMinDistanceM);
  CHECK(flat_hit && std::abs(flat_hit->point.z - 1.0) < 1e-15);

  // The facet grid returns the same nearest facet as testing all facets.
  std::vector<Facet> layout;
  std::uint32_t id = 0;
  for (int i = -5; i <= 5; ++i)
    for (int j = -5; j <= 5; ++j)
      layout.push_back({id++,
                        {1.2 * i + 0.6 * (j % 2), 1.04 * j, 0.02 * (i * i + j * j)},
                        normalized({-0.02 * i, -0.02 * j, 1}),
                        {},
                        {},
                        {},
                        {},
                        0.6,
                        30.0,
                        Aperture::hexagon_flat_y});
  for (auto &l : layout) {
    l.u = normalized(cross({0, 1, 0}, l.normal));
    l.v = cross(l.normal, l.u);
  }
  const FacetGrid grid(layout);
  int compared = 0;
  for (std::uint64_t k = 0; k < 20000; ++k) {
    const auto u = uniform4(5, k, Stream::source);
    const Ray probe{{14 * (u[0] - 0.5), 14 * (u[1] - 0.5), 40.0},
                    normalized({0.2 * (u[2] - 0.5), 0.2 * (u[3] - 0.5), -1})};
    std::optional<SurfaceHit> all, gridded;
    for (const auto &l : layout)
      if (auto h = l.intersect(probe, kMinDistanceM); h && (!all || h->t < all->t))
        all = h;
    grid.candidates(probe, layout.size(), [&](std::size_t i) {
      if (auto h = layout[i].intersect(probe, kMinDistanceM); h && (!gridded || h->t < gridded->t))
        gridded = h;
    });
    CHECK(all.has_value() == gridded.has_value());
    if (all && gridded) {
      CHECK(all->id == gridded->id);
      ++compared;
    }
  }
  CHECK(compared > 10000);
}

// First sign change of the sag residual along the ray, by dense sampling.
double brute_force_root(const Asphere &a, const Ray &ray, double t0, double t1) {
  const auto g = [&](double t) {
    const Vec3 p = ray.at(t);
    return p.z - a.sag(std::hypot(p.x, p.y));
  };
  const int n = 200000;
  double previous = g(t0);
  for (int i = 1; i <= n; ++i) {
    const double t = t0 + (t1 - t0) * i / n, value = g(t);
    if ((previous < 0) != (value < 0)) {
      double lo = t - (t1 - t0) / n, hi = t;
      for (int k = 0; k < 200; ++k) {
        const double mid = 0.5 * (lo + hi);
        ((g(mid) < 0) == (g(lo) < 0) ? lo : hi) = mid;
      }
      return 0.5 * (lo + hi);
    }
    previous = value;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

void aspheres() {
  const auto accept = [](const Vec3 &) { return true; };
  bool failed = false;
  // Paraboloid z = r^2 / (4 f): compare with the closed-form root.
  const double f = 10.0;
  const Asphere parabola{1, 0.0, 0.0, 6.0, 1.0, {0.0, 1.0 / (4 * f)}};
  for (std::uint64_t k = 0; k < 2000; ++k) {
    const auto u = uniform4(9, k, Stream::source);
    const Ray ray{{10 * (u[0] - 0.5), 10 * (u[1] - 0.5), 20.0},
                  normalized({0.3 * (u[2] - 0.5), 0.3 * (u[3] - 0.5), -1})};
    const auto hit = parabola.intersect(ray, kMinDistanceM, accept, failed);
    CHECK(!failed);
    const double a =
        (ray.direction.x * ray.direction.x + ray.direction.y * ray.direction.y) / (4 * f);
    const double b = (ray.origin.x * ray.direction.x + ray.origin.y * ray.direction.y) / (2 * f) -
                     ray.direction.z;
    const double c =
        (ray.origin.x * ray.origin.x + ray.origin.y * ray.origin.y) / (4 * f) - ray.origin.z;
    const double t = (-b - std::sqrt(b * b - 4 * a * c)) / (2 * a);
    const Vec3 p = ray.at(t);
    if (std::hypot(p.x, p.y) <= 6.0) {
      CHECK(hit.has_value());
      if (hit)
        CHECK_NEAR(hit->t, t, 1e-9);
    }
  }
  // An SST-like primary (physical-centimetre polynomial, here in metres) and
  // an SCT-like high-order secondary: the solver finds the first crossing.
  std::vector<double> sst_cm{0.0,          0.000608051,  -7.36248e-10, 5.82125e-15,
                             -6.76447e-20, -3.89542e-24, -5.28036e-29, 2.99107e-34,
                             4.39153e-38,  6.17433e-43,  -2.73586e-47};
  std::vector<double> sst_m;
  for (std::size_t i = 0; i < sst_cm.size(); ++i)
    sst_m.push_back(sst_cm[i] * std::pow(0.01, 1.0 - 2.0 * static_cast<double>(i)));
  const Asphere sst{2, 0.0, 0.2, 2.12, 1.0, sst_m};
  std::vector<double> sct;
  for (const double c : {1.5, 0.416667, 0.145816, -0.712012, 4.17685, -23.1617, 118.844, -520.501,
                         1802.98, -4605.3, 8019.1, -8422.9, 4004.61})
    sct.push_back(c * 5.5863);
  const Asphere sct_like{3, 0.0, 0.0, 2.7, 5.5863, sct};
  for (const Asphere *surface : {&sst, &sct_like}) {
    for (std::uint64_t k = 0; k < 200; ++k) {
      const auto u = uniform4(13, k, Stream::source);
      // Include near-axial rays from 10 km, where precision is hardest.
      const double angle = k % 4 == 0 ? 1e-6 * u[2] : 0.05 * (u[2] - 0.5);
      const Vec3 d = normalized({angle, 0.03 * (u[3] - 0.5) * (k % 4 != 0), -1});
      const Vec3 target{4 * (u[0] - 0.5), 4 * (u[1] - 0.5), 0.0};
      const double distance = k % 2 ? 1e4 : 30.0;
      const Ray ray{target - d * distance, d};
      const auto hit = surface->intersect(ray, kMinDistanceM, accept, failed);
      CHECK(!failed);
      const double expected = brute_force_root(*surface, ray, distance - 40.0, distance + 40.0);
      const Vec3 p = ray.at(expected);
      if (!std::isnan(expected) && std::hypot(p.x, p.y) <= surface->outer_radius &&
          std::hypot(p.x, p.y) >= surface->inner_radius) {
        CHECK(hit.has_value());
        if (hit) {
          CHECK_NEAR(hit->t, expected, 1e-6);
          CHECK_NEAR(hit->point.z, surface->sag(std::hypot(hit->point.x, hit->point.y)), 1e-9);
        }
      }
    }
  }
}

// An SST-like asphere with a central hole and hexagonal segments separated by
// gaps (pitch 0.45 m, diameter 0.4 m).
struct SegmentedAsphere {
  Asphere asphere;
  std::vector<Segment> segments;
  std::vector<Segment::TangentFrame> frames;
  SegmentGrid grid;

  [[nodiscard]] bool accept(const Vec3 &p) const {
    for (std::size_t i = 0; i < segments.size(); ++i)
      if (segments[i].contains(frames[i], p))
        return true;
    return false;
  }

  // The lowest id of the segments containing p (a brute-force loop).
  [[nodiscard]] std::optional<std::uint32_t> lowest_id(const Vec3 &p) const {
    std::optional<std::uint32_t> id;
    for (std::size_t i = 0; i < segments.size(); ++i)
      if (segments[i].contains(frames[i], p) && (!id || segments[i].id < *id))
        id = segments[i].id;
    return id;
  }
};

SegmentedAsphere sst_like_segmented() {
  std::vector<double> sst_cm{0.0,          0.000608051,  -7.36248e-10, 5.82125e-15,
                             -6.76447e-20, -3.89542e-24, -5.28036e-29, 2.99107e-34,
                             4.39153e-38,  6.17433e-43,  -2.73586e-47};
  SegmentedAsphere m;
  m.asphere = Asphere{2, 0.0, 0.2, 2.12, 1.0, {}};
  for (std::size_t i = 0; i < sst_cm.size(); ++i)
    m.asphere.coefficients.push_back(sst_cm[i] *
                                     std::pow(0.01, 1.0 - 2.0 * static_cast<double>(i)));
  m.asphere.prepare();
  std::uint32_t id = 3;
  for (int row = -4; row <= 4; ++row)
    for (int column = -4; column <= 4; ++column) {
      const double x = 0.45 * (column + 0.5 * (row & 1)), y = 0.45 * 0.866 * row;
      if (std::hypot(x, y) < 0.3 || std::hypot(x, y) > 1.9)
        continue;
      Segment s;
      s.id = id++;
      s.shape = Segment::Shape::hexagon;
      s.centre_x = x;
      s.centre_y = y;
      s.half_width = 0.2;
      m.segments.push_back(s);
      m.frames.push_back(s.tangent_frame(m.asphere));
    }
  m.grid = SegmentGrid(m.segments, m.frames, m.asphere.outer_radius);
  return m;
}

void slope_bound() {
  const auto m = sst_like_segmented();
  CHECK(std::isfinite(m.asphere.max_slope) && m.asphere.max_slope > 0.0);
  std::vector<double> sct;
  for (const double c : {1.5, 0.416667, 0.145816, -0.712012, 4.17685, -23.1617, 118.844, -520.501,
                         1802.98, -4605.3, 8019.1, -8422.9, 4004.61})
    sct.push_back(c * 5.5863);
  const Asphere sct_like{3, 0.0, 0.0, 2.7, 5.5863, sct};
  for (const Asphere *surface : {&m.asphere, &sct_like}) {
    const double bound = surface->slope_bound();
    double largest = 0.0;
    for (int k = 0; k <= 1000; ++k) {
      const double slope = std::abs(surface->slope(surface->outer_radius * k / 1000.0));
      CHECK(slope <= bound);
      largest = std::max(largest, slope);
    }
    CHECK(bound >= largest);
    CHECK(bound <= 1e4 * std::max(largest, 1e-3));
  }
  CHECK(Asphere{}.max_slope == std::numeric_limits<double>::infinity());
}

// The single-root shortcut must give exactly the hits and misses of the
// search without it.
void single_root_shortcut() {
  const auto m = sst_like_segmented();
  Asphere reference = m.asphere;
  reference.max_slope = std::numeric_limits<double>::infinity();
  const auto accept = [&](const Vec3 &p) { return m.accept(p); };
  std::size_t hits = 0, misses = 0, shortcut_rays = 0;
  for (std::uint64_t k = 0; k < 200000; ++k) {
    const auto u = uniform4(21, k, Stream::source);
    // Mostly steep rays from far above, some tilted by up to 60 degrees.
    const double tilt = (k % 5 == 0 ? 1.05 : 0.15) * u[2];
    const double azimuth = 6.283185307179586 * u[3];
    const Vec3 d{std::sin(tilt) * std::cos(azimuth), std::sin(tilt) * std::sin(azimuth),
                 -std::cos(tilt)};
    const double radius = 2.3 * std::sqrt(u[0]), phi = 6.283185307179586 * u[1];
    const Vec3 target{radius * std::cos(phi), radius * std::sin(phi), 0.0};
    const Ray ray{target - d * (k % 2 ? 3.0 : 1e4), d};
    bool failed = false, failed_reference = false;
    const auto a = m.asphere.intersect(ray, kMinDistanceM, accept, failed);
    const auto b = reference.intersect(ray, kMinDistanceM, accept, failed_reference);
    CHECK(failed == failed_reference);
    CHECK(a.has_value() == b.has_value());
    if (a && b) {
      CHECK(a->t == b->t);
      CHECK(a->point.x == b->point.x && a->point.y == b->point.y && a->point.z == b->point.z);
      ++hits;
    } else {
      ++misses;
    }
    shortcut_rays += std::abs(d.z) > m.asphere.max_slope * std::hypot(d.x, d.y);
  }
  CHECK(hits > 20000 && misses > 20000 && shortcut_rays > 20000);
}

// The grid lookup finds the lowest-id segment that the brute-force loop finds,
// for every segment shape.
void segment_grid() {
  auto m = sst_like_segmented();
  // Mixed shapes with ids out of order: sectors (one wrapping past 2 pi), a
  // polygon, a square, a circle and a rotated hexagon.
  m.segments.clear();
  m.frames.clear();
  const auto add = [&](Segment s, std::uint32_t id) {
    s.id = id;
    m.segments.push_back(s);
    m.frames.push_back(s.tangent_frame(m.asphere));
  };
  Segment sector;
  sector.shape = Segment::Shape::annular_sector;
  sector.inner_radius = 0.3;
  sector.outer_radius = 1.9;
  sector.start = 5.5;
  sector.span = 1.5;
  sector.gap = 0.02;
  add(sector, 40);
  sector.start = 0.4;
  sector.span = 2.5;
  sector.inner_radius = 0.9;
  sector.gap_at_start = true;
  add(sector, 12);
  sector.start = 3.0;
  sector.span = 0.5;
  sector.inner_radius = 0.2;
  add(sector, 9);
  Segment polygon;
  polygon.shape = Segment::Shape::polygon;
  polygon.vertices = {{-1.0, -0.2}, {0.3, -1.5}, {0.9, -0.4}, {-0.2, 0.3}};
  add(polygon, 30);
  Segment other;
  other.shape = Segment::Shape::square;
  other.centre_x = 1.2;
  other.centre_y = 0.9;
  other.half_width = 0.35;
  other.rotation = 0.4;
  add(other, 5);
  other.shape = Segment::Shape::circle;
  other.centre_x = -1.3;
  other.centre_y = 1.0;
  add(other, 7);
  other.shape = Segment::Shape::hexagon;
  other.centre_x = 0.1;
  other.centre_y = 1.4;
  other.half_width = 0.3;
  add(other, 2);
  // Overlapping footprints: the lowest id must win.
  other.centre_x = 0.2;
  add(other, 1);
  m.grid = SegmentGrid(m.segments, m.frames, m.asphere.outer_radius);
  std::size_t found = 0, overlapping = 0;
  for (std::uint64_t k = 0; k < 300000; ++k) {
    const auto u = uniform4(31, k, Stream::source);
    const double outer = m.asphere.outer_radius;
    const double x = outer * (2.0 * u[0] - 1.0), y = outer * (2.0 * u[1] - 1.0);
    const Vec3 p{x, y, m.asphere.sag(std::hypot(x, y))};
    const auto expected = m.lowest_id(p);
    const auto actual = m.grid.find(m.segments, m.frames, p);
    CHECK(actual == expected);
    found += expected.has_value();
    overlapping += expected && *expected == 1;
  }
  CHECK(found > 30000 && overlapping > 100);
  CHECK(!m.grid.find(m.segments, m.frames, {5.0, 0.0, 0.0}).has_value());
  // The hexagonal mirror, on points near its segment edges as well.
  const auto hexagons = sst_like_segmented();
  for (std::uint64_t k = 0; k < 300000; ++k) {
    const auto u = uniform4(32, k, Stream::source);
    const double x = 2.12 * (2.0 * u[0] - 1.0), y = 2.12 * (2.0 * u[1] - 1.0);
    const Vec3 p{x, y, hexagons.asphere.sag(std::hypot(x, y))};
    CHECK(hexagons.grid.find(hexagons.segments, hexagons.frames, p) == hexagons.lowest_id(p));
  }
}

void segments() {
  const Asphere flat{1, 0.0, 0.0, 5.0, 1.0, {0.0}};
  Segment hexagon;
  hexagon.shape = Segment::Shape::hexagon;
  hexagon.centre_x = 1.0;
  hexagon.half_width = 0.5;
  CHECK(hexagon.contains(flat, {1.49, 0.0, 0.0}));
  CHECK(!hexagon.contains(flat, {1.51, 0.0, 0.0}));
  Segment sector;
  sector.shape = Segment::Shape::annular_sector;
  sector.inner_radius = 1.0;
  sector.outer_radius = 2.0;
  sector.span = 1.0;
  sector.gap = 0.05;
  CHECK(sector.contains(flat, {1.5, 0.5, 0.0}));
  CHECK(!sector.contains(flat, {1.5, -0.1, 0.0}));
  CHECK(!sector.contains(flat, {0.5, 0.1, 0.0}));
  Segment triangle;
  triangle.shape = Segment::Shape::polygon;
  triangle.vertices = {{0, 0}, {1, 0}, {0, 1}};
  CHECK(triangle.contains(flat, {0.2, 0.2, 0.0}));
  CHECK(!triangle.contains(flat, {0.8, 0.8, 0.0}));
  // On a curved mirror the precomputed tangent frame gives the same footprint.
  const Asphere curved{1, 0.0, 0.0, 5.0, 1.0, {0.0, 0.05, -0.001}};
  hexagon.centre_y = 0.6;
  hexagon.rotation = 0.3;
  const auto frame = hexagon.tangent_frame(curved);
  for (double x = 0.4; x < 1.6; x += 0.01)
    for (double y = 0.0; y < 1.2; y += 0.01) {
      const Vec3 p{x, y, curved.sag(std::hypot(x, y))};
      CHECK(hexagon.contains(frame, p) == hexagon.contains(curved, p));
    }
  CHECK(hexagon.contains(frame, {1.0, 0.6, curved.sag(std::hypot(1.0, 0.6))}));
}

void obscurers() {
  std::vector<Plane> planes;
  Obscurer cylinder;
  cylinder.kind = Obscurer::Kind::cylinder;
  cylinder.a = {0, 0, 0};
  cylinder.b = {0, 0, 10};
  cylinder.radius_a = 0.5;
  CHECK_NEAR(*cylinder.intersect({{3, 0, 5}, {-1, 0, 0}}, planes, kMinDistanceM), 2.5, 1e-12);
  CHECK_NEAR(*cylinder.intersect({{0.2, 0, 20}, {0, 0, -1}}, planes, kMinDistanceM), 10.0,
             1e-12); // End cap.
  CHECK(!cylinder.intersect({{3, 0, 11}, {-1, 0, 0}}, planes, kMinDistanceM));

  Obscurer cube;
  cube.kind = Obscurer::Kind::polyhedron;
  for (const Vec3 n : {Vec3{1, 0, 0}, Vec3{-1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, -1, 0}, Vec3{0, 0, 1},
                       Vec3{0, 0, -1}})
    planes.push_back({n, 1.0});
  cube.plane_end = 6;
  CHECK_NEAR(*cube.intersect({{5, 0.5, 0.5}, {-1, 0, 0}}, planes, kMinDistanceM), 4.0, 1e-12);
  CHECK(!cube.intersect({{5, 1.5, 0}, {-1, 0, 0}}, planes, kMinDistanceM));
  CHECK_NEAR(*cube.intersect({{0, 0, 0}, {1, 0, 0}}, planes, kMinDistanceM), kMinDistanceM,
             1e-15); // Starts inside.

  Obscurer plate;
  plate.kind = Obscurer::Kind::plate;
  plate.centre = {0, 0, 2};
  plate.normal = {0, 0, 1};
  plate.u = {1, 0, 0};
  plate.v = {0, 1, 0};
  plate.half_width = 1.0;
  plate.shape = Aperture::square;
  CHECK_NEAR(*plate.intersect({{0.9, -0.9, 5}, {0, 0, -1}}, planes, kMinDistanceM), 3.0, 1e-12);
  CHECK(!plate.intersect({{1.1, 0, 5}, {0, 0, -1}}, planes, kMinDistanceM));

  Obscurer cone;
  cone.kind = Obscurer::Kind::frustum;
  cone.a = {0, 0, 0};
  cone.b = {0, 0, 2};
  cone.radius_a = 1.0;
  cone.radius_b = 0.5;
  cone.solid = true;
  // Horizontal ray at z = 1 meets the wall at radius 0.75.
  CHECK_NEAR(*cone.intersect({{2, 0, 1}, {-1, 0, 0}}, planes, kMinDistanceM), 1.25, 1e-12);
  CHECK_NEAR(*cone.intersect({{0.2, 0, 5}, {0, 0, -1}}, planes, kMinDistanceM), 3.0,
             1e-12); // Top face.
  cone.solid = false;
  cone.thickness = 0.1;
  CHECK(!cone.intersect({{0.2, 0, 5}, {0, 0, -1}}, planes, kMinDistanceM)); // Open tube.
}

// The per-leg hierarchy finds the same nearest obscurer as testing them all.
void obscurer_hierarchy() {
  std::vector<Obscurer> rods;
  for (std::uint32_t i = 0; i < 300; ++i) {
    const auto u = uniform4(21, i, Stream::source), w = uniform4(22, i, Stream::source);
    Obscurer rod;
    rod.kind = Obscurer::Kind::cylinder;
    rod.id = 1000 - i;
    rod.legs = static_cast<std::uint8_t>(1 + i % 7);
    rod.a = {20 * u[0] - 10, 20 * u[1] - 10, 15 * u[2]};
    rod.b = rod.a + Vec3{2 * w[0] - 1, 2 * w[1] - 1, 2 * w[2] - 1};
    rod.radius_a = 0.05 + 0.2 * u[3];
    rod.bound_a = rod.a;
    rod.bound_b = rod.b;
    rod.bound_radius = rod.radius_a;
    rods.push_back(rod);
  }
  const Obscurers set(rods, {});
  int hits = 0;
  for (std::uint64_t k = 0; k < 20000; ++k) {
    const auto u = uniform4(23, k, Stream::source), w = uniform4(24, k, Stream::source);
    const Ray ray{{24 * u[0] - 12, 24 * u[1] - 12, 20.0},
                  normalized({0.6 * w[0] - 0.3, 0.6 * w[1] - 0.3, -1.0})};
    const double t_max = 5.0 + 20.0 * w[2];
    const Leg leg = k % 3 == 0   ? Leg::incoming
                    : k % 3 == 1 ? Leg::after_primary
                                 : Leg::after_secondary;
    std::optional<Obscurers::Hit> expected;
    for (const auto &rod : rods)
      if (rod.legs & leg)
        if (const auto t = rod.intersect(ray, {}, kMinDistanceM))
          if (*t < t_max &&
              (!expected || *t < expected->t || (*t == expected->t && rod.id < expected->id)))
            expected = Obscurers::Hit{*t, rod.id};
    const auto found = set.nearest(ray, leg, t_max);
    CHECK(found.has_value() == expected.has_value());
    if (found && expected) {
      CHECK(found->id == expected->id);
      CHECK(found->t == expected->t);
      ++hits;
    }
  }
  CHECK(hits > 500); // About 900 rays meet a rod.
  // A polyhedron's capsule encloses its vertices, hence the solid.
  const std::vector<Vec3> prism{
      {5.5, 3.0, 1.9}, {5.6, 0.1, 1.7}, {5.5, 0.1, 2.0}, {5.5, 3.0, 2.2}, {5.4, 3.0, 1.9}};
  Obscurer beam;
  bound_vertices(beam, prism);
  for (const auto &v : prism) {
    const Vec3 axis = normalized(beam.bound_b - beam.bound_a);
    const double along = dot(v - beam.bound_a, axis);
    CHECK(along >= -1e-12 && along <= norm(beam.bound_b - beam.bound_a) + 1e-12);
    CHECK(norm(v - beam.bound_a - axis * along) <= beam.bound_radius + 1e-12);
  }
  CHECK(beam.bound_radius < 0.4); // A sphere would need 1.5 m.
  CHECK(!Obscurers().nearest({{0, 0, 1}, {0, 0, -1}}, Leg::incoming, 10.0));
}

// Random slender beams (polyhedra with a vertex box) and rods: the grid for steep
// incoming rays, the hierarchy for every leg and a brute-force loop give the
// same nearest hit.
void obscurer_grid() {
  std::vector<Obscurer> items;
  std::vector<Plane> planes;
  for (std::uint32_t i = 0; i < 150; ++i) {
    const auto u = uniform4(41, i, Stream::source), w = uniform4(42, i, Stream::source);
    Obscurer o;
    o.id = 5000 - i;
    o.legs = static_cast<std::uint8_t>(1 + i % 7);
    if (i % 2) {
      o.kind = Obscurer::Kind::cylinder;
      o.a = {16 * u[0] - 8, 16 * u[1] - 8, 12 * u[2]};
      o.b = o.a + Vec3{4 * w[0] - 2, 4 * w[1] - 2, 4 * w[2] - 2};
      o.radius_a = 0.03 + 0.1 * u[3];
      o.bound_a = o.a;
      o.bound_b = o.b;
      o.bound_radius = o.radius_a;
    } else {
      // An oriented box: axes from two random directions.
      const Vec3 e1 = normalized({w[0] - 0.5, w[1] - 0.5, w[2] - 0.5});
      const Vec3 e2 = normalized(cross(e1, {0.3, 0.5 + u[3], 0.1}));
      const Vec3 e3 = cross(e1, e2);
      const Vec3 centre{16 * u[0] - 8, 16 * u[1] - 8, 12 * u[2]};
      const double h1 = 0.5 + 1.5 * u[2], h2 = 0.02 + 0.08 * w[3], h3 = 0.02 + 0.08 * u[3];
      std::vector<Vec3> corners;
      for (int a : {-1, 1})
        for (int b : {-1, 1})
          for (int c : {-1, 1})
            corners.push_back(centre + e1 * (a * h1) + e2 * (b * h2) + e3 * (c * h3));
      o.kind = Obscurer::Kind::polyhedron;
      o.plane_begin = static_cast<std::uint32_t>(planes.size());
      for (const auto &[axis, h] : {std::pair{e1, h1}, std::pair{e2, h2}, std::pair{e3, h3}})
        for (const double sign : {-1.0, 1.0})
          planes.push_back({axis * sign, dot(axis * sign, centre) + h});
      o.plane_end = static_cast<std::uint32_t>(planes.size());
      bound_vertices(o, corners);
    }
    items.push_back(o);
  }
  const Obscurers set(items, planes);
  std::size_t hits = 0, steep = 0;
  for (std::uint64_t k = 0; k < 60000; ++k) {
    const auto u = uniform4(43, k, Stream::source), w = uniform4(44, k, Stream::source);
    // Mostly steep rays, some shallow ones.
    const double tilt = k % 4 == 0 ? 1.5 : 0.8;
    const Ray ray{{20 * u[0] - 10, 20 * u[1] - 10, 20.0},
                  normalized({tilt * (w[0] - 0.5), tilt * (w[1] - 0.5), -1.0})};
    const double t_max = 8.0 + 20.0 * w[2];
    const Leg leg = k % 3 == 0   ? Leg::incoming
                    : k % 3 == 1 ? Leg::after_primary
                                 : Leg::after_secondary;
    std::optional<Obscurers::Hit> expected;
    for (const auto &item : items)
      if (item.legs & leg)
        if (const auto t = item.intersect(ray, planes, kMinDistanceM))
          if (*t < t_max &&
              (!expected || *t < expected->t || (*t == expected->t && item.id < expected->id)))
            expected = Obscurers::Hit{*t, item.id};
    for (const auto &found :
         {set.nearest(ray, leg, t_max), set.nearest_in_hierarchy(ray, leg, t_max)}) {
      CHECK(found.has_value() == expected.has_value());
      if (found && expected) {
        CHECK(found->id == expected->id);
        CHECK(found->t == expected->t);
      }
    }
    hits += expected.has_value();
    steep += leg == Leg::incoming && std::abs(ray.direction.z) >= 0.5;
  }
  CHECK(hits > 1500 && steep > 10000);
}

// The plane pre-test rejects only rays that the sphere test also rejects, and
// the tighter grid cells lose no facet.
void facet_prefilter() {
  std::size_t hits = 0;
  for (const Aperture shape :
       {Aperture::circle, Aperture::hexagon_flat_x, Aperture::square, Aperture::hexagon_flat_y})
    for (const double radius : {0.0, 40.0, 3.0}) { // Flat, gently and strongly curved.
      Facet facet;
      facet.id = 7;
      facet.centre = {0.3, -0.2, 1.0};
      facet.normal = normalized({0.2, 0.1, 1.0});
      facet.u = normalized(cross(facet.normal, {0.0, 1.0, 0.0}));
      facet.v = cross(facet.normal, facet.u);
      facet.shape = shape;
      facet.half_width = 0.5;
      facet.radius = radius;
      Facet checked = facet;
      checked.prepare();
      CHECK(checked.reach > 0.0 && (radius > 0.0 ? checked.sag > 0.0 : checked.sag == 0.0));
      for (std::uint64_t k = 0; k < 40000; ++k) {
        const auto u = uniform4(51 + k % 3, k, Stream::source);
        // Rays aimed near the facet from various directions, including shallow ones.
        const Vec3 target = facet.centre + facet.u * (1.6 * (u[0] - 0.5)) +
                            facet.v * (1.6 * (u[1] - 0.5)) + facet.normal * (0.1 * (u[2] - 0.5));
        const Vec3 d = normalized({1.4 * (u[3] - 0.5) + 0.1, 1.4 * (u[2] - 0.5), -1.0});
        const double distance = k % 2 ? 50.0 : 1e4;
        const Ray ray{target - d * distance, d};
        const auto a = checked.intersect(ray, kMinDistanceM),
                   b = facet.intersect(ray, kMinDistanceM);
        CHECK(a.has_value() == b.has_value());
        if (a && b) {
          CHECK(a->t == b->t && a->id == b->id);
          ++hits;
        }
      }
    }
  CHECK(hits > 20000);
  // Tighter cells: the grid still offers every facet that a ray hits.
  std::vector<Facet> facets;
  for (int row = -6; row <= 6; ++row)
    for (int column = -6; column <= 6; ++column) {
      Facet f;
      f.id = static_cast<std::uint32_t>(facets.size());
      f.centre = {0.9 * (column + 0.5 * (row & 1)), 0.78 * row,
                  0.001 * (row * row + column * column)};
      f.normal = normalized({-0.02 * f.centre.x, -0.02 * f.centre.y, 1.0});
      f.u = normalized(cross(f.normal, {0.0, 1.0, 0.0}));
      f.v = cross(f.normal, f.u);
      f.shape = Aperture::hexagon_flat_x;
      f.half_width = 0.45;
      f.radius = 30.0;
      f.prepare();
      facets.push_back(f);
    }
  const FacetGrid grid(facets);
  std::size_t found = 0;
  for (std::uint64_t k = 0; k < 40000; ++k) {
    const auto u = uniform4(61, k, Stream::source);
    const Vec3 d = normalized({0.6 * (u[2] - 0.5), 0.6 * (u[3] - 0.5), -1.0});
    const Ray ray{{12.0 * (u[0] - 0.5) - d.x * 20.0, 10.0 * (u[1] - 0.5) - d.y * 20.0, 20.0}, d};
    std::vector<char> offered(facets.size(), 0);
    grid.candidates(ray, facets.size(), [&](std::size_t i) { offered[i] = 1; });
    for (std::size_t i = 0; i < facets.size(); ++i)
      if (facets[i].intersect(ray, kMinDistanceM)) {
        CHECK(offered[i]);
        ++found;
      }
  }
  CHECK(found > 10000);
}

void sources() {
  // A plane wave fills the disk uniformly and starts on the start plane.
  PointSource star{
      normalized({0.01, 0.0, 1.0}), std::numeric_limits<double>::infinity(), 5.0, 0.0, 30.0, 400.0};
  double inside = 0.0;
  for (std::uint64_t i = 0; i < 100000; ++i) {
    const auto p = sample_source(star, i, 3);
    CHECK_NEAR(p.ray.origin.z, 30.0, 1e-12);
    const Vec3 target = p.ray.at(-p.ray.origin.z / p.ray.direction.z);
    CHECK(std::hypot(target.x, target.y) <= 5.0 + 1e-9);
    inside += std::hypot(target.x, target.y) < 2.5;
  }
  CHECK_NEAR(inside / 100000.0, 0.25, 0.005);
  // A nearby point source emits from its own position below the start plane;
  // its distance is measured from the disk centre (sim_telarray).
  PointSource lamp{normalized({0.6, 0, 0.8}), 20.0, 1.0, 1.5, 30.0, 400.0};
  const auto p = sample_source(lamp, 0, 3);
  CHECK_NEAR(p.ray.origin.x, 12.0, 1e-12);
  CHECK_NEAR(p.ray.origin.z, 17.5, 1e-12);
  CHECK_NEAR(p.time_ns, -20.0 / kSpeedOfLightMPerNs,
             1e-9); // Wavefront reaches the disk centre at t = 0.
  // A star's wavefront reaches the disk centre at t = 0.
  star.disk_z_m = 1.5;
  const auto s = sample_source(star, 7, 3);
  CHECK_NEAR(dot(s.ray.origin - Vec3{0, 0, 1.5}, s.ray.direction), s.time_ns * kSpeedOfLightMPerNs,
             1e-9);
}

void densities() {
  // Density 2x on [0, 1], 2 on [1, 3]: integral 5, CDF x^2 / 5 below 1.
  const LinearDensity d({0.0, 1.0, 3.0}, {0.0, 2.0, 2.0});
  CHECK_NEAR(d.integral(), 5.0, 1e-15);
  CHECK_NEAR(d.maximum(), 2.0, 0.0);
  CHECK_NEAR(d.at(0.5), 1.0, 1e-15);
  CHECK_NEAR(d.at(3.5), 0.0, 0.0);
  CHECK_NEAR(d.sample(0.1), std::sqrt(0.5), 1e-15);
  CHECK_NEAR(d.sample(0.6), 2.0, 1e-15);
  CHECK_NEAR(d.sample(1.0), 3.0, 1e-15);
  CHECK_NEAR(LinearDensity().sample(0.3), 0.0, 0.0); // Empty: point mass at zero.
  bool rejected = false;
  try {
    const LinearDensity bad({1.0, 0.0}, {1.0, 1.0});
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  CHECK(rejected);
}

void emitters() {
  // An isotropic emitter 10 m above the centre of a 2 m disk: the photons
  // crossing it per emitted photon are the disk's solid angle over 4 pi.
  Emitter lamp{{0.0, 0.0, 10.0},
               {0.0, 0.0, -1.0},
               LinearDensity({0.0, 2.0}, {1.0, 1.0}),
               LinearDensity({-1.0, 1.0}, {1.0, 1.0}),
               {0.0, 0.0, 0.0},
               {0.0, 0.0, 1.0},
               2.0,
               355.0};
  CHECK_NEAR(lamp.candidates_per_photon(), 4.0 / (4.0 * 100.0), 1e-15);
  const std::size_t n = 200000;
  double weights = 0.0, time = 0.0, time2 = 0.0;
  for (std::uint64_t i = 0; i < n; ++i) {
    const auto p = sample_emitter(lamp, i, 5);
    CHECK(p.weight > 0.0 && p.weight <= 1.0);
    CHECK_NEAR(norm(p.ray.origin - lamp.position), 0.0, 0.0);
    const Vec3 target = p.ray.at(-p.ray.origin.z / p.ray.direction.z);
    CHECK(std::hypot(target.x, target.y) <= 2.0 + 1e-9);
    weights += p.weight;
    time += p.time_ns;
    time2 += p.time_ns * p.time_ns;
  }
  const double solid_angle = 2.0 * 3.141592653589793 * (1.0 - 10.0 / std::sqrt(104.0));
  CHECK_NEAR(weights / n * lamp.candidates_per_photon(), solid_angle / (4.0 * 3.141592653589793),
             2e-5);
  CHECK_NEAR(time / n, 0.0, 0.003); // Flat pulse on [-1, 1]: mean 0, variance 1/3.
  CHECK_NEAR(time2 / n, 1.0 / 3.0, 0.003);
  // A 5 deg filled cone lights only the inner part of the disk.
  const double u5 = 1.0 - std::cos(5.0 / 57.29577951308232);
  lamp.intensity = LinearDensity({0.0, u5, u5 * 1.000001}, {1.0, 1.0, 0.0});
  for (std::uint64_t i = 0; i < 10000; ++i) {
    const auto p = sample_emitter(lamp, i, 5);
    const double angle = std::acos(-p.ray.direction.z) * 57.29577951308232;
    CHECK(angle < 5.0 + 1e-6 ? p.weight > 0.0 : p.weight == 0.0);
  }
}

void transmission() {
  const TelescopeTransmission t{0.9, 0.5, 0.1, 2.0, 1.0};
  CHECK_NEAR(t.at({0, 0, -1}), 0.9, 1e-15);
  CHECK_NEAR(t.at(normalized({0.1, 0, -1})),
             0.9 / (1 + 0.5 * std::pow(0.1 / std::sqrt(1.01) / 0.1, 2)), 1e-12);
}

} // namespace

void all_tests() {
  random_numbers();
  tables();
  apertures();
  facets();
  facet_prefilter();
  aspheres();
  slope_bound();
  single_root_shortcut();
  segment_grid();
  segments();
  obscurers();
  obscurer_hierarchy();
  obscurer_grid();
  sources();
  densities();
  emitters();
  transmission();
}

int main() { return run_checks(all_tests); }
