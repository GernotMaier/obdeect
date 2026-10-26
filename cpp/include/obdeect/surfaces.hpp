#pragma once

#include "obdeect/vec.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>
#ifdef OBDEECT_BENCH_COUNTERS
#include <string>
#endif

namespace obdeect {

// Planar aperture outlines. `half_width` is half the circle diameter or half
// the flat-to-flat width (sim_telarray mirror shape types 0, 2, 1 and 3).
enum class Aperture : std::uint8_t { circle, square, hexagon_flat_x, hexagon_flat_y };

[[nodiscard]] inline bool inside_aperture(Aperture shape, double half_width, double u, double v) {
  constexpr double s = 0.86602540378443864676; // sqrt(3)/2
  switch (shape) {
  case Aperture::circle:
    return u * u + v * v <= half_width * half_width;
  case Aperture::square:
    return std::abs(u) <= half_width && std::abs(v) <= half_width;
  case Aperture::hexagon_flat_x:
    return std::abs(u) <= half_width && std::abs(0.5 * u + s * v) <= half_width &&
           std::abs(0.5 * u - s * v) <= half_width;
  case Aperture::hexagon_flat_y:
    return std::abs(v) <= half_width && std::abs(0.5 * v + s * u) <= half_width &&
           std::abs(0.5 * v - s * u) <= half_width;
  }
  return false;
}

// Radius of the circle that encloses an aperture.
[[nodiscard]] inline double circumradius(Aperture shape, double half_width) {
  switch (shape) {
  case Aperture::circle:
    return half_width;
  case Aperture::square:
    return half_width * std::sqrt(2.0);
  case Aperture::hexagon_flat_x:
  case Aperture::hexagon_flat_y:
    return half_width * 2.0 / std::sqrt(3.0);
  }
  return half_width;
}

struct SurfaceHit {
  double t{};
  Vec3 point;
  Vec3 normal; // Unit normal on the reflecting side.
  std::uint32_t id{};
};

// A spherical (or flat, radius 0) mirror facet. `centre` is the vertex on the
// optical surface, `normal` points from the vertex towards the centre of
// curvature, and (u, v) span the aperture plane.
struct Facet {
  std::uint32_t id{};
  Vec3 centre, normal, u, v;
  Vec3 map_u, map_v; // Coordinate axes of the reflectivity map.
  double half_width{};
  double radius{}; // Radius of curvature; 0 for a flat facet.
  Aperture shape{Aperture::circle};
  // Radius of the circle around the vertex that encloses the aperture, and the
  // sag of the cap at that radius. reach == 0 (unset) disables the plane pre-test.
  double reach{}, sag{};

  void prepare() {
    reach = circumradius(shape, half_width);
    sag = radius > 0.0 ? radius - std::sqrt(std::max(0.0, radius * radius - reach * reach)) : 0.0;
  }

  [[nodiscard]] std::optional<SurfaceHit> intersect(const Ray &ray, double t_min) const {
    if (radius == 0.0) {
      const double denominator = dot(ray.direction, normal);
      if (denominator >= 0.0) // Parallel or reaching the back face.
        return std::nullopt;
      const double t = dot(centre - ray.origin, normal) / denominator;
      if (!(t > t_min))
        return std::nullopt;
      const Vec3 p = ray.at(t), offset = p - centre;
      if (!inside_aperture(shape, half_width, dot(offset, u), dot(offset, v)))
        return std::nullopt;
      return SurfaceHit{t, p, normal, id};
    }
    if (reach > 0.0) {
      // Reject at the facet plane. Every point of the cap inside the aperture
      // lies within `reach` of the vertex sideways and at most `sag` from the
      // plane along the normal; a ray that moves h along the normal moves at
      // most h sqrt(1 - dn^2) / |dn| sideways. So a ray farther from the vertex
      // than reach + sag sqrt(1 - dn^2) / |dn| where it crosses the plane
      // cannot meet the cap.
      const double dn = dot(ray.direction, normal);
      if (dn != 0.0) {
        const double tq = dot(centre - ray.origin, normal) / dn;
        const Vec3 q = ray.origin + ray.direction * tq - centre;
        const double along = dot(q, normal);
        const double limit =
            reach + sag * std::sqrt(std::max(0.0, 1.0 - dn * dn)) / std::abs(dn) + 1e-9;
        if (dot(q, q) - along * along > limit * limit)
          return std::nullopt;
      }
    }
    const Vec3 sphere = centre + normal * radius;
    const Vec3 oc = ray.origin - sphere;
    const double b = dot(oc, ray.direction);
    const double discriminant = b * b - (dot(oc, oc) - radius * radius);
    if (discriminant < 0.0)
      return std::nullopt;
    const double root = std::sqrt(discriminant);
    for (const double t : {-b - root, -b + root}) {
      if (!(t > t_min))
        continue;
      const Vec3 p = ray.at(t);
      const Vec3 inward = (sphere - p) * (1.0 / radius);
      // Only the concave cap around the vertex reflects, and only from its
      // front; back-face crossings are not mirror hits (as in sim_telarray).
      if (dot(inward, normal) <= 0.0 || dot(ray.direction, inward) >= 0.0)
        continue;
      const Vec3 offset = p - centre;
      if (inside_aperture(shape, half_width, dot(offset, u), dot(offset, v)))
        return SurfaceHit{t, p, inward, id};
    }
    return std::nullopt;
  }
};

// Uniform grid over the facets' footprints in the x-y plane. A ray is tested
// only against facets in the cells its path crosses between the lowest and
// highest facet heights; near-horizontal rays test all facets.
class FacetGrid {
public:
  // Cell size as a fraction of the smallest facet diameter.
  static constexpr double kCellScale = 0.25;

  FacetGrid() = default;
  explicit FacetGrid(const std::vector<Facet> &facets) {
    if (facets.empty())
      return;
    double reach = 0.0, cell = std::numeric_limits<double>::infinity();
    z_low_ = std::numeric_limits<double>::infinity();
    z_high_ = -z_low_;
    x0_ = y0_ = std::numeric_limits<double>::infinity();
    double x1 = -x0_, y1 = -y0_;
    for (const auto &f : facets) {
      const double r = circumradius(f.shape, f.half_width);
      // The facet's sag over its footprint bounds its height range.
      const double sag =
          f.radius > 0.0 ? f.radius - std::sqrt(std::max(0.0, f.radius * f.radius - r * r)) : 0.0;
      z_low_ = std::min(z_low_, f.centre.z - r - sag);
      z_high_ = std::max(z_high_, f.centre.z + r + sag);
      reach = std::max(reach, r);
      cell = std::min(cell, 2.0 * r);
      x0_ = std::min(x0_, f.centre.x - r);
      y0_ = std::min(y0_, f.centre.y - r);
      x1 = std::max(x1, f.centre.x + r);
      y1 = std::max(y1, f.centre.y + r);
    }
    cell_ = std::max(cell * kCellScale, 1e-3);
    nx_ = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil((x1 - x0_) / cell_)));
    ny_ = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil((y1 - y0_) / cell_)));
    std::vector<std::vector<std::uint32_t>> cells(nx_ * ny_);
    for (std::uint32_t i = 0; i < facets.size(); ++i) {
      const auto &f = facets[i];
      const double r = circumradius(f.shape, f.half_width);
      // The footprint lies within the circle of radius r (plus the sag of a
      // curved facet) around the centre: skip the cells of the bounding square
      // that this circle does not reach.
      const double sag =
          f.radius > 0.0 ? f.radius - std::sqrt(std::max(0.0, f.radius * f.radius - r * r)) : 0.0;
      const double reach = r + sag;
      for (std::size_t ix = index(f.centre.x - r, x0_, nx_); ix <= index(f.centre.x + r, x0_, nx_);
           ++ix)
        for (std::size_t iy = index(f.centre.y - r, y0_, ny_);
             iy <= index(f.centre.y + r, y0_, ny_); ++iy) {
          const double dx = std::max({x0_ + cell_ * static_cast<double>(ix) - f.centre.x,
                                      f.centre.x - x0_ - cell_ * static_cast<double>(ix + 1), 0.0});
          const double dy = std::max({y0_ + cell_ * static_cast<double>(iy) - f.centre.y,
                                      f.centre.y - y0_ - cell_ * static_cast<double>(iy + 1), 0.0});
          if (dx * dx + dy * dy <= reach * reach)
            cells[ix * ny_ + iy].push_back(i);
        }
    }
    offsets_.assign(1, 0);
    for (const auto &c : cells) {
      members_.insert(members_.end(), c.begin(), c.end());
      offsets_.push_back(members_.size());
    }
  }

  // Calls `visit(facet_index)` for every candidate facet (possibly twice).
  template <class Visit>
  void candidates(const Ray &ray, std::size_t facet_count, Visit visit) const {
    const double dz = ray.direction.z;
    if (std::abs(dz) < 0.05 || nx_ == 0) {
      for (std::size_t i = 0; i < facet_count; ++i)
        visit(i);
      return;
    }
    const Vec3 a = ray.at((z_low_ - ray.origin.z) / dz), b = ray.at((z_high_ - ray.origin.z) / dz);
    const std::size_t ix0 = index(std::min(a.x, b.x), x0_, nx_),
                      ix1 = index(std::max(a.x, b.x), x0_, nx_);
    const std::size_t iy0 = index(std::min(a.y, b.y), y0_, ny_),
                      iy1 = index(std::max(a.y, b.y), y0_, ny_);
    if (std::max(a.x, b.x) < x0_ || std::min(a.x, b.x) > x0_ + cell_ * nx_ ||
        std::max(a.y, b.y) < y0_ || std::min(a.y, b.y) > y0_ + cell_ * ny_)
      return;
    for (std::size_t ix = ix0; ix <= ix1; ++ix)
      for (std::size_t iy = iy0; iy <= iy1; ++iy)
        for (std::size_t k = offsets_[ix * ny_ + iy]; k < offsets_[ix * ny_ + iy + 1]; ++k)
          visit(members_[k]);
  }

private:
  [[nodiscard]] std::size_t index(double value, double origin, std::size_t count) const {
    const double cell = (value - origin) / cell_;
    if (!(cell > 0.0))
      return 0;
    return std::min(count - 1, static_cast<std::size_t>(cell));
  }

  double x0_{}, y0_{}, cell_{1.0}, z_low_{}, z_high_{};
  std::size_t nx_{}, ny_{};
  std::vector<std::size_t> offsets_;
  std::vector<std::uint32_t> members_;
};

// c + value * q; with OBDEECT_USE_FMA a fused multiply-add, which changes results
// at the rounding level (see CMakeLists.txt).
[[nodiscard]] inline double multiply_add(double c, double value, double q) {
#ifdef OBDEECT_USE_FMA
  return std::fma(value, q, c);
#else
  return c + value * q;
#endif
}

// Rotationally symmetric surface z(r) = vertex + sum_i c_i (r / scale)^(2i),
// within inner_radius <= r <= outer_radius, in the telescope frame.
struct Asphere {
  std::uint32_t id{};
  double vertex_z{}, inner_radius{}, outer_radius{}, scale{1.0};
  std::vector<double> coefficients; // c_0 ... c_n in metres.
  // Upper bound of |dz/dr| on [0, outer_radius]; infinity disables the
  // single-root shortcut in intersect().
  double max_slope{std::numeric_limits<double>::infinity()};
  // Vertex height and rim curvature of the first-guess paraboloid; NaN until
  // cache_guess() is called, then intersect() computes them per ray.
  double guess_z0{std::numeric_limits<double>::quiet_NaN()}, guess_k{guess_z0};
  // 1 / scale^2 and the squared aperture radii; NaN until prepare() is called,
  // then intersect() computes them per ray (as for the first-guess constants).
  double inv_scale2{guess_z0}, inner_radius2{guess_z0}, outer_radius2{guess_z0};

  // Computes max_slope, the first-guess constants and the squared quantities;
  // call after the shape is set.
  void prepare() {
    max_slope = slope_bound();
    inv_scale2 = 1.0 / (scale * scale);
    inner_radius2 = inner_radius * inner_radius;
    outer_radius2 = outer_radius * outer_radius;
    guess_z0 = sag(0.0);
    guess_k = (sag(outer_radius) - guess_z0) / (outer_radius * outer_radius);
  }

#ifdef OBDEECT_BENCH_COUNTERS
  // Benchmark-only counters (single thread); never compiled into the library.
  struct Counters {
    std::uint64_t calls{}, newton_steps{};
  };
  static Counters &counters() {
    static Counters c;
    return c;
  }
#endif

  // Rigorous bound sum_i 2 i |c_i| rho^(2i-1) / scale, rho = outer_radius / scale.
  // Every term grows with r, so the value at the outer radius bounds |slope(r)|
  // everywhere inside. A relative margin covers the rounding of the sum.
  [[nodiscard]] double slope_bound() const {
    const double rho = outer_radius / scale;
    double bound = 0.0;
    for (std::size_t i = coefficients.size(); i-- > 1;)
      bound = bound * rho * rho + 2.0 * static_cast<double>(i) * std::abs(coefficients[i]);
    return bound * rho / scale * (1.0 + 1e-12);
  }

  [[nodiscard]] double sag(double r) const {
    const double q = (r / scale) * (r / scale);
    double value = 0.0;
    for (auto c = coefficients.rbegin(); c != coefficients.rend(); ++c)
      value = value * q + *c;
    return vertex_z + value;
  }

  // dz/dr.
  [[nodiscard]] double slope(double r) const {
    const double rho = r / scale, q = rho * rho;
    double value = 0.0;
    for (std::size_t i = coefficients.size(); i-- > 1;)
      value = value * q + 2.0 * static_cast<double>(i) * coefficients[i];
    return value * rho / scale;
  }

  // The same sag in q = (x^2 + y^2) / scale^2 (Horner in q).
  [[nodiscard]] double sag_q(double q) const {
    double value = 0.0;
    for (auto c = coefficients.rbegin(); c != coefficients.rend(); ++c)
      value = multiply_add(*c, value, q);
    return vertex_z + value;
  }

  // (dz/dr) / r as a polynomial in q: (1 / scale^2) sum_{i>=1} 2 i c_i q^(i-1).
  // No special case at r = 0: the slope term is then k * 0 = 0.
  [[nodiscard]] double slope_over_r_q(double q) const {
    double value = 0.0;
    for (std::size_t i = coefficients.size(); i-- > 1;)
      value = multiply_add(2.0 * static_cast<double>(i) * coefficients[i], value, q);
    return value * (std::isnan(inv_scale2) ? 1.0 / (scale * scale) : inv_scale2);
  }

  [[nodiscard]] Vec3 normal_at(const Vec3 &p) const {
    const double inv = std::isnan(inv_scale2) ? 1.0 / (scale * scale) : inv_scale2;
    return normal_from_q(p, (p.x * p.x + p.y * p.y) * inv);
  }

  // The normal at p for q = (p.x^2 + p.y^2) / scale^2 already known.
  [[nodiscard]] Vec3 normal_from_q(const Vec3 &p, double q) const {
    const double k = slope_over_r_q(q);
    return normalized(Vec3{-k * p.x, -k * p.y, 1.0});
  }

  // Nearest intersection beyond t_min inside the annulus for which
  // accept(point) holds. Sets `failed` if the root search did not converge.
  template <class Accept>
  [[nodiscard]] std::optional<SurfaceHit> intersect(const Ray &ray, double t_min, Accept accept,
                                                    bool &failed) const {
    failed = false;
#ifdef OBDEECT_BENCH_COUNTERS
    ++counters().calls;
#endif
    // Work from a point near the vertex plane: distant sources otherwise
    // lose precision in every residual.
    const double shift = std::max(0.0, dot(Vec3{0.0, 0.0, vertex_z} - ray.origin, ray.direction));
    const Ray local{ray.at(shift), ray.direction};
    const double local_min = t_min - shift;
    // Everything in rr = p.x^2 + p.y^2 and q = rr / scale^2: no hypot, sqrt or
    // division per step.
    const double inv_s2 = std::isnan(inv_scale2) ? 1.0 / (scale * scale) : inv_scale2;
    const double inner2 = std::isnan(inner_radius2) ? inner_radius * inner_radius : inner_radius2;
    const double outer2 = std::isnan(outer_radius2) ? outer_radius * outer_radius : outer_radius2;
    // Residual and derivative at the point p = local.at(t) with rr = p.x^2 + p.y^2.
    const auto residual_at = [&](const Vec3 &p, double rr) { return p.z - sag_q(rr * inv_s2); };
    const auto derivative_at = [&](const Vec3 &p, double rr) {
      const Vec3 &d = local.direction;
      return d.z - slope_over_r_q(rr * inv_s2) * (p.x * d.x + p.y * d.y);
    };
    const auto residual = [&](double t) {
      const Vec3 p = local.at(t);
      return residual_at(p, p.x * p.x + p.y * p.y);
    };
    // A root t at the point p = local.at(t), rr = p.x^2 + p.y^2.
    struct Root {
      double t;
      Vec3 p;
      double rr;
    };
    const auto accept_root = [&](const Root &root) -> std::optional<SurfaceHit> {
      if (!(root.t > local_min))
        return std::nullopt;
      if (root.rr < inner2 || root.rr > outer2 || !accept(root.p))
        return std::nullopt;
      return SurfaceHit{root.t + shift, root.p, normal_from_q(root.p, root.rr * inv_s2), id};
    };
    const auto root_at = [&](double t) {
      const Vec3 p = local.at(t);
      return Root{t, p, p.x * p.x + p.y * p.y};
    };
    const auto newton = [&](double t, double lo, double hi) -> std::optional<Root> {
      for (int iteration = 0; iteration < 50; ++iteration) {
#ifdef OBDEECT_BENCH_COUNTERS
        ++counters().newton_steps;
#endif
        const Vec3 p = local.at(t);
        const double rr = p.x * p.x + p.y * p.y;
        const double g = residual_at(p, rr);
        if (std::abs(g) <= kTolerance)
          return Root{t, p, rr};
        const double gd = derivative_at(p, rr);
        if (gd == 0.0 || !std::isfinite(gd))
          return std::nullopt;
        const double next = t - g / gd;
        if (!(next >= lo && next <= hi) || next == t)
          return std::abs(g) <= 100 * kTolerance ? std::optional<Root>{Root{t, p, rr}}
                                                 : std::nullopt;
        t = next;
      }
      return std::nullopt;
    };

    // Range of t over which the ray lies inside the outer aperture cylinder.
    const Vec3 &o = local.origin, &d = local.direction;
    const double a2 = d.x * d.x + d.y * d.y, b2 = o.x * d.x + o.y * d.y;
    const double c2 = o.x * o.x + o.y * o.y - outer2;
    double lo = local_min, hi = std::numeric_limits<double>::infinity();
    if (a2 > 0.0) {
      const double disc = b2 * b2 - a2 * c2;
      if (disc < 0.0)
        return std::nullopt;
      const double root = std::sqrt(disc);
      lo = std::max(lo, (-b2 - root) / a2);
      hi = (-b2 + root) / a2;
      if (hi < lo)
        return std::nullopt;
    } else if (c2 > 0.0) {
      return std::nullopt;
    }

    // Along the ray the residual is g(t) = p_z(t) - sag(r(t)), with derivative
    // g'(t) = d_z - slope(r) (p_xy . d_xy) / r. As |p_xy . d_xy| / r <= |d_xy|,
    // g' has the sign of d_z whenever |d_z| > max_slope |d_xy|: g is strictly
    // monotone and the surface has at most one intersection, the only root.
    // A root that is not accepted then ends the search, and so does the same
    // sign of g at both ends of [lo, hi] when Newton found no root.
    const bool single_root = std::abs(d.z) > max_slope * std::sqrt(a2);

    // First guess: the paraboloid through the vertex and the outer rim.
    const double z0 = std::isnan(guess_z0) ? sag(0.0) : guess_z0;
    const double k = std::isnan(guess_k)
                         ? (sag(outer_radius) - sag(0.0)) / (outer_radius * outer_radius)
                         : guess_k;
    const double qa = -k * a2, qb = d.z - 2.0 * k * b2, qc = o.z - z0 - k * (o.x * o.x + o.y * o.y);
    std::array<double, 2> guesses{std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::quiet_NaN()};
    if (std::abs(qa) > 1e-14 * std::abs(qb)) {
      const double disc = qb * qb - 4.0 * qa * qc;
      if (disc >= 0.0) {
        const double q = -0.5 * (qb + std::copysign(std::sqrt(disc), qb));
        guesses = {q / qa, qc / q};
      }
    } else if (qb != 0.0) {
      guesses[0] = -qc / qb;
    }
    // Ascending, NaN last.
    if (!std::isnan(guesses[1]) && (std::isnan(guesses[0]) || guesses[1] < guesses[0]))
      std::swap(guesses[0], guesses[1]);
    // With the single root proven there is at most one intersection, so the
    // start point cannot change which root is found, only how fast Newton gets
    // there: try a guess inside [lo, hi] before one outside it (which would be
    // clamped to the end of the range). Otherwise keep the ascending order.
    if (single_root) {
      const auto inside = [&](double g) { return g >= lo && (!std::isfinite(hi) || g <= hi); };
      if (!std::isnan(guesses[1]) && !inside(guesses[0]) && inside(guesses[1]))
        std::swap(guesses[0], guesses[1]);
    }
    const double span = std::isfinite(hi) ? hi - lo : 0.0;
    for (const double guess : guesses) {
      if (std::isnan(guess))
        continue;
      const double start = std::isfinite(hi) ? std::clamp(guess, lo, hi) : std::max(guess, lo);
      const double upper = std::isfinite(hi) ? hi : start + 1e3;
      if (const auto t = newton(start, lo, upper)) {
        if (auto hit = accept_root(*t))
          return hit;
        if (single_root)
          return std::nullopt;
      }
    }
    if (single_root && std::isfinite(hi)) {
      const double g_lo = residual(lo), g_hi = residual(hi);
      if (std::isfinite(g_lo) && std::isfinite(g_hi) &&
          ((g_lo > 0.0 && g_hi > 0.0) || (g_lo < 0.0 && g_hi < 0.0)))
        return std::nullopt;
    }
    // Fallback: bracket sign changes along the aperture interval.
    if (!(span > 0.0))
      return std::nullopt;
    constexpr int samples = 64;
    double t0 = lo, g0 = residual(t0);
    for (int i = 1; i <= samples; ++i) {
      const double t1 = lo + span * i / samples, g1 = residual(t1);
      if (!std::isfinite(g0) || !std::isfinite(g1)) {
        failed = true;
        return std::nullopt;
      }
      if (g0 == 0.0 || (g0 < 0.0) != (g1 < 0.0)) {
        double a = t0, b = t1, ga = g0;
        for (int iteration = 0; iteration < 200 && b - a > 1e-13 * std::max(1.0, std::abs(b));
             ++iteration) {
          const double m = 0.5 * (a + b), gm = residual(m);
          if ((gm < 0.0) == (ga < 0.0)) {
            a = m;
            ga = gm;
          } else {
            b = m;
          }
        }
        if (auto hit = accept_root(root_at(0.5 * (a + b))))
          return hit;
      }
      t0 = t1;
      g0 = g1;
    }
    return std::nullopt;
  }

  static constexpr double kTolerance = 1e-10; // Metres along z.
};

// Footprint of one segment of an aspheric mirror (sim_telarray segmentation).
struct Segment {
  enum class Shape : std::uint8_t { hexagon, square, circle, annular_sector, polygon };
  std::uint32_t id{};
  Shape shape{};
  double centre_x{}, centre_y{}, half_width{}, rotation{};       // Hexagon, square, circle.
  double inner_radius{}, outer_radius{}, start{}, span{}, gap{}; // Annular sector.
  bool gap_at_start{};
  std::vector<std::array<double, 2>> vertices; // Convex polygon.

  // Tangent plane at the segment centre, in which hexagons, squares and circles
  // keep their physical flat-to-flat width. It depends only on the mirror and
  // the segment, so an optical model computes it once (Mirror::segment_frames).
  // `reach2` bounds the squared x-y distance from the centre of a point on the
  // mirror inside the segment: the tangent plane only shrinks in projection, and
  // hexagon, square and circle stay within 1.5 half widths of their centre.
  struct TangentFrame {
    double cx{1.0}, sx{}, cosine{1.0}, sine{}, centre_z{}, cos_rotation{1.0}, sin_rotation{};
    double reach2{std::numeric_limits<double>::infinity()};
  };

  [[nodiscard]] TangentFrame tangent_frame(const Asphere &mirror) const {
    const double rc = std::hypot(centre_x, centre_y);
    const double s = mirror.slope(rc), cosine = 1.0 / std::sqrt(1.0 + s * s);
    return {rc > 0.0 ? centre_x / rc : 1.0,
            rc > 0.0 ? centre_y / rc : 0.0,
            cosine,
            s * cosine,
            mirror.sag(rc),
            std::cos(rotation),
            std::sin(rotation),
            2.25 * half_width * half_width};
  }

  [[nodiscard]] bool contains(const Asphere &mirror, const Vec3 &p) const {
    return contains(tangent_frame(mirror), p);
  }

  [[nodiscard]] bool contains(const TangentFrame &f, const Vec3 &p) const {
    switch (shape) {
    case Shape::polygon: {
      double sign = 0.0;
      for (std::size_t i = 0; i < vertices.size(); ++i) {
        const auto &a = vertices[i], &b = vertices[(i + 1) % vertices.size()];
        const double side = (b[0] - a[0]) * (p.y - a[1]) - (b[1] - a[1]) * (p.x - a[0]);
        if (side == 0.0)
          continue;
        if (sign != 0.0 && (side < 0.0) != (sign < 0.0))
          return false;
        sign = side;
      }
      return true;
    }
    case Shape::annular_sector: {
      const double r = std::hypot(p.x, p.y);
      constexpr double two_pi = 6.283185307179586476925;
      double angle = std::atan2(p.y, p.x) - start;
      angle -= two_pi * std::floor(angle / two_pi);
      const double from_gap = gap_at_start ? angle * r : (span - angle) * r;
      return r >= inner_radius && r < outer_radius && angle < span && r >= gap && from_gap >= gap;
    }
    default: {
      const double dx = p.x - centre_x, dy = p.y - centre_y, dz = p.z - f.centre_z;
      if (dx * dx + dy * dy > f.reach2)
        return false;
      const double radial = (dx * f.cx + dy * f.sx) * f.cosine + dz * f.sine,
                   azimuthal = -dx * f.sx + dy * f.cx;
      const double x = radial * f.cx - azimuthal * f.sx, y = radial * f.sx + azimuthal * f.cx;
      const double u = x * f.cos_rotation + y * f.sin_rotation;
      const double v = -x * f.sin_rotation + y * f.cos_rotation;
      if (shape == Shape::circle)
        return u * u + v * v <= half_width * half_width;
      if (shape == Shape::square)
        return std::abs(u) <= half_width && std::abs(v) <= half_width;
      return inside_aperture(Aperture::hexagon_flat_x, half_width, u, v);
    }
    }
  }
};

// Uniform x-y grid over the square [-outer_radius, outer_radius]^2 of an
// aspheric mirror. Each cell lists, sorted by segment id, the segments whose
// conservative x-y bounding box overlaps it, so that a point is tested against
// the few segments near it instead of all of them. A cell whose result is the
// same for every point it can receive carries a verdict, which find() returns
// without testing any segment.
class SegmentGrid {
public:
  static constexpr std::size_t kCells = 128; // Cells per side.

#ifdef OBDEECT_BENCH_COUNTERS
  // Benchmark-only counters (single thread); never compiled into the library.
  struct Counters {
    std::uint64_t calls{}, contains_calls{}, decided{};
  };
  static Counters &counters() {
    static Counters c;
    return c;
  }
  static std::string counter_report(std::size_t hits) {
    const auto &c = counters();
    return "segment find calls " + std::to_string(static_cast<double>(c.calls) / hits) +
           ", segment contains calls " +
           std::to_string(static_cast<double>(c.contains_calls) / hits) +
           " (per primary hit, all runs), decided by a verdict " +
           std::to_string(c.calls ? static_cast<double>(c.decided) / c.calls : 0.0) +
           " of the find calls";
  }
#endif

  SegmentGrid() = default;
  SegmentGrid(const std::vector<Segment> &segments,
              const std::vector<Segment::TangentFrame> &frames, const Asphere &mirror)
      : outer_(mirror.outer_radius),
        inverse_cell_(static_cast<double>(kCells) / (2.0 * mirror.outer_radius)) {
    std::vector<std::uint32_t> by_id(segments.size());
    for (std::uint32_t i = 0; i < by_id.size(); ++i)
      by_id[i] = i;
    std::sort(by_id.begin(), by_id.end(), [&](std::uint32_t a, std::uint32_t b) {
      return segments[a].id < segments[b].id || (segments[a].id == segments[b].id && a < b);
    });
    std::array<std::vector<std::uint32_t>, kCells * kCells> cells;
    for (const std::uint32_t i : by_id) {
      const auto box = bounding_box(segments[i], frames[i]);
      // Hexagon, square and circle accept only points within sqrt(reach2) of
      // their centre: skip the corner cells of the box that are farther away.
      const bool round = segments[i].shape != Segment::Shape::polygon &&
                         segments[i].shape != Segment::Shape::annular_sector;
      const double reach = std::sqrt(frames[i].reach2) + 1e-9;
      for (std::size_t ix = index(box[0]); ix <= index(box[1]); ++ix)
        for (std::size_t iy = index(box[2]); iy <= index(box[3]); ++iy) {
          if (round) {
            const double dx = distance_to_cell(segments[i].centre_x, ix),
                         dy = distance_to_cell(segments[i].centre_y, iy);
            if (dx * dx + dy * dy > reach * reach)
              continue;
          }
          cells[ix * kCells + iy].push_back(i);
        }
    }
    offsets_.assign(1, 0);
    for (const auto &c : cells) {
      members_.insert(members_.end(), c.begin(), c.end());
      offsets_.push_back(static_cast<std::uint32_t>(members_.size()));
    }
    make_verdicts(segments, frames, mirror);
  }

  // Number of cells with a verdict (for tests and diagnostics).
  [[nodiscard]] std::size_t decided_cells() const {
    return static_cast<std::size_t>(
        std::count_if(verdicts_.begin(), verdicts_.end(),
                      [](const Verdict &v) { return v.kind != Kind::undecided; }));
  }

  // Id of the lowest-id segment that contains p, or nullopt if none does or p
  // lies outside the grid.
  [[nodiscard]] std::optional<std::uint32_t> find(const std::vector<Segment> &segments,
                                                  const std::vector<Segment::TangentFrame> &frames,
                                                  const Vec3 &p) const {
#ifdef OBDEECT_BENCH_COUNTERS
    ++counters().calls;
#endif
    if (members_.empty() || !(std::abs(p.x) <= outer_ && std::abs(p.y) <= outer_))
      return std::nullopt;
    const std::size_t cell = index(p.x) * kCells + index(p.y);
    const Verdict &verdict = verdicts_[cell];
    if (verdict.kind != Kind::undecided && p.z >= verdict.z_low && p.z <= verdict.z_high) {
#ifdef OBDEECT_BENCH_COUNTERS
      ++counters().decided;
#endif
      if (verdict.kind == Kind::none)
        return std::nullopt;
      return verdict.id;
    }
    for (std::uint32_t k = offsets_[cell]; k < offsets_[cell + 1]; ++k) {
      const std::uint32_t i = members_[k];
#ifdef OBDEECT_BENCH_COUNTERS
      ++counters().contains_calls;
#endif
      if (segments[i].contains(frames[i], p))
        return segments[i].id;
    }
    return std::nullopt;
  }

private:
  // Result of find() for every point in the cell's x-y box with z in
  // [z_low, z_high]: no segment, the id of the segment returned, or undecided.
  enum class Kind : std::uint32_t { undecided, none, sole };
  struct Verdict {
    float z_low{}, z_high{}; // Rounded inwards: never wider than the range proven.
    std::uint32_t id{};
    Kind kind{Kind::undecided};
  };
  enum class Cover { none, all, mixed }; // No / every / some point of a box is inside.

  // Margin (m) of all verdict tests, far above the rounding of Segment::contains.
  static constexpr double kMargin = 1e-9;

  // For each cell: the first member (in order) that contains the whole box
  // decides, provided all members before it contain none of it; a cell whose
  // members contain none of it has no result. Anything else stays undecided.
  void make_verdicts(const std::vector<Segment> &segments,
                     const std::vector<Segment::TangentFrame> &frames, const Asphere &mirror) {
    verdicts_.assign(kCells * kCells, Verdict{});
    if (!std::isfinite(mirror.max_slope))
      return;
    const double cell = 1.0 / inverse_cell_;
    for (std::size_t ix = 0; ix < kCells; ++ix)
      for (std::size_t iy = 0; iy < kCells; ++iy) {
        const std::size_t c = ix * kCells + iy;
        // x-y box of the cell, widened for the rounding of index().
        const double x0 = static_cast<double>(ix) * cell - outer_ - kMargin,
                     x1 = static_cast<double>(ix + 1) * cell - outer_ + kMargin,
                     y0 = static_cast<double>(iy) * cell - outer_ - kMargin,
                     y1 = static_cast<double>(iy + 1) * cell - outer_ + kMargin;
        // Mirror points in the box have |r - r_c| <= h and |dz/dr| <= max_slope.
        const double rc = std::hypot(0.5 * (x0 + x1), 0.5 * (y0 + y1)),
                     h = 0.5 * std::hypot(x1 - x0, y1 - y0), dz = mirror.max_slope * h + 1e-6,
                     zc = mirror.sag(rc);
        const double z0 = zc - dz, z1 = zc + dz;
        Verdict v{inward_float(z0, true), inward_float(z1, false), 0, Kind::none};
        for (std::uint32_t k = offsets_[c]; k < offsets_[c + 1]; ++k) {
          const std::uint32_t i = members_[k];
          const Cover cover = classify(segments[i], frames[i], x0, x1, y0, y1, z0, z1);
          if (cover == Cover::none)
            continue;
          v.kind = cover == Cover::all ? Kind::sole : Kind::undecided;
          v.id = segments[i].id;
          break;
        }
        verdicts_[c] = v;
      }
  }

  // The float nearest to `value` that is not below it (`up`) or not above it.
  static float inward_float(double value, bool up) {
    float f = static_cast<float>(value);
    if (up ? static_cast<double>(f) < value : static_cast<double>(f) > value)
      f = std::nextafter(f, up ? std::numeric_limits<float>::infinity()
                               : -std::numeric_limits<float>::infinity());
    return f;
  }

  // Whether Segment::contains accepts all, none or some of the box (Cover::mixed
  // when unsure). Hexagon, square and circle map the box affinely to (u, v), so
  // the images of its 8 corners bound every point of it.
  static Cover classify(const Segment &s, const Segment::TangentFrame &f, double x0, double x1,
                        double y0, double y1, double z0, double z1) {
    if (s.shape == Segment::Shape::annular_sector)
      return Cover::mixed;
    if (s.shape == Segment::Shape::polygon)
      return classify_polygon(s, x0, x1, y0, y1);
    // Nearest and farthest x-y distance from the centre to the box.
    const double dx = std::max({x0 - s.centre_x, s.centre_x - x1, 0.0}),
                 dy = std::max({y0 - s.centre_y, s.centre_y - y1, 0.0});
    const double reach = std::sqrt(f.reach2);
    if (dx * dx + dy * dy > (reach + kMargin) * (reach + kMargin))
      return Cover::none;
    const double fx = std::max(std::abs(x0 - s.centre_x), std::abs(x1 - s.centre_x)),
                 fy = std::max(std::abs(y0 - s.centre_y), std::abs(y1 - s.centre_y));
    const bool within_reach = fx * fx + fy * fy <= (reach - kMargin) * (reach - kMargin);
    // Same formulas as Segment::contains.
    std::array<std::array<double, 2>, 8> uv;
    std::size_t n = 0;
    for (const double x : {x0, x1})
      for (const double y : {y0, y1})
        for (const double z : {z0, z1}) {
          const double ddx = x - s.centre_x, ddy = y - s.centre_y, ddz = z - f.centre_z;
          const double radial = (ddx * f.cx + ddy * f.sx) * f.cosine + ddz * f.sine,
                       azimuthal = -ddx * f.sx + ddy * f.cx;
          const double px = radial * f.cx - azimuthal * f.sx, py = radial * f.sx + azimuthal * f.cx;
          uv[n++] = {px * f.cos_rotation + py * f.sin_rotation,
                     -px * f.sin_rotation + py * f.cos_rotation};
        }
    const double w = s.half_width;
    constexpr double sq = 0.86602540378443864676;
    // Linear forms whose absolute value must stay below w (hexagon, square),
    // or, for the circle, the 8 directions of the enclosing octagon.
    std::array<std::array<double, 2>, 8> forms;
    std::size_t m = 0;
    if (s.shape == Segment::Shape::hexagon)
      for (const auto &a :
           std::array<std::array<double, 2>, 3>{{{1.0, 0.0}, {0.5, sq}, {0.5, -sq}}})
        forms[m++] = a;
    else {
      forms[m++] = {1.0, 0.0};
      forms[m++] = {0.0, 1.0};
      if (s.shape == Segment::Shape::circle) {
        forms[m++] = {0.7071067811865475, 0.7071067811865475};
        forms[m++] = {0.7071067811865475, -0.7071067811865475};
      }
    }
    bool all = within_reach;
    for (std::size_t j = 0; j < m; ++j) {
      double lo = std::numeric_limits<double>::infinity(), hi = -lo;
      for (const auto &q : uv) {
        const double value = forms[j][0] * q[0] + forms[j][1] * q[1];
        lo = std::min(lo, value);
        hi = std::max(hi, value);
      }
      // The form range lies wholly outside [-w, w] on one side: nothing inside.
      if (lo > w + kMargin || hi < -w - kMargin)
        return Cover::none;
      all = all && lo >= -w + kMargin && hi <= w - kMargin;
    }
    if (s.shape == Segment::Shape::circle) {
      // The disk is convex: all corner images inside it means every point is.
      for (const auto &q : uv)
        all = all && q[0] * q[0] + q[1] * q[1] <= (w - kMargin) * (w - kMargin);
    }
    return all ? Cover::all : Cover::mixed;
  }

  // Convex polygon, x-y only: as Segment::contains, whose edge sides share one
  // sign inside. A non-convex or degenerate polygon is never decided.
  static Cover classify_polygon(const Segment &s, double x0, double x1, double y0, double y1) {
    const auto &vs = s.vertices;
    const std::size_t n = vs.size();
    if (n < 3)
      return Cover::mixed;
    const auto side = [&](std::size_t i, double x, double y) {
      const auto &a = vs[i], &b = vs[(i + 1) % n];
      return (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]);
    };
    double area = 0.0;
    for (std::size_t i = 0; i < n; ++i)
      area += vs[i][0] * vs[(i + 1) % n][1] - vs[(i + 1) % n][0] * vs[i][1];
    if (!(std::abs(area) > 0.0))
      return Cover::mixed;
    const double orientation = area > 0.0 ? 1.0 : -1.0;
    for (std::size_t i = 0; i < n; ++i) { // Convexity: every turn has the area's sign.
      const auto &a = vs[i], &b = vs[(i + 1) % n], &c = vs[(i + 2) % n];
      if (orientation * ((b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])) <= 0.0)
        return Cover::mixed;
    }
    bool all = true;
    for (std::size_t i = 0; i < n; ++i) {
      const auto &a = vs[i], &b = vs[(i + 1) % n];
      const double margin = kMargin * std::hypot(b[0] - a[0], b[1] - a[1]);
      double lo = std::numeric_limits<double>::infinity(), hi = -lo;
      for (const double x : {x0, x1})
        for (const double y : {y0, y1}) {
          const double value = orientation * side(i, x, y);
          lo = std::min(lo, value);
          hi = std::max(hi, value);
        }
      if (hi < -margin) // Wholly outside this edge: some other edge is positive.
        return Cover::none;
      all = all && lo > margin;
    }
    return all ? Cover::all : Cover::mixed;
  }

  [[nodiscard]] std::size_t index(double value) const {
    const double cell = (value + outer_) * inverse_cell_;
    if (!(cell > 0.0))
      return 0;
    return std::min(kCells - 1, static_cast<std::size_t>(cell));
  }

  // Distance along one axis from `value` to the cell with index `i` (0 inside it).
  [[nodiscard]] double distance_to_cell(double value, std::size_t i) const {
    const double low = static_cast<double>(i) / inverse_cell_ - outer_,
                 high = static_cast<double>(i + 1) / inverse_cell_ - outer_;
    return std::max({low - value, value - high, 0.0});
  }

  // {x_min, x_max, y_min, y_max} of the points Segment::contains can accept,
  // widened by a margin far above the rounding of the tests.
  static std::array<double, 4> bounding_box(const Segment &s, const Segment::TangentFrame &frame) {
    constexpr double margin = 1e-9;
    double x0 = std::numeric_limits<double>::infinity(), x1 = -x0, y0 = x0, y1 = -x0;
    const auto include = [&](double x, double y) {
      x0 = std::min(x0, x);
      x1 = std::max(x1, x);
      y0 = std::min(y0, y);
      y1 = std::max(y1, y);
    };
    switch (s.shape) {
    case Segment::Shape::polygon:
      for (const auto &v : s.vertices)
        include(v[0], v[1]);
      break;
    case Segment::Shape::annular_sector: {
      constexpr double two_pi = 6.283185307179586476925;
      const double end = s.start + s.span;
      for (const double radius : {s.inner_radius, s.outer_radius})
        for (const double angle : {s.start, end})
          include(radius * std::cos(angle), radius * std::sin(angle));
      // The outer arc reaches +-outer_radius on an axis inside the angular span.
      for (int axis = 0; axis < 4; ++axis) {
        double from_start = 0.25 * axis * two_pi - s.start;
        from_start -= two_pi * std::floor(from_start / two_pi);
        if (from_start < s.span)
          include(s.outer_radius * ((axis == 0) - (axis == 2)),
                  s.outer_radius * ((axis == 1) - (axis == 3)));
      }
      break;
    }
    default: {
      // Segment::contains rejects points farther than sqrt(reach2) from the centre.
      const double reach = std::sqrt(frame.reach2);
      include(s.centre_x - reach, s.centre_y - reach);
      include(s.centre_x + reach, s.centre_y + reach);
    }
    }
    return {x0 - margin, x1 + margin, y0 - margin, y1 + margin};
  }

  double outer_{}, inverse_cell_{1.0};
  std::vector<std::uint32_t> offsets_, members_;
  std::vector<Verdict> verdicts_;
};

} // namespace obdeect
