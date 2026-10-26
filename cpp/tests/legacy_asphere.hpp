#pragma once
// Reference copy of the Asphere intersection before the r^2 formulation (hypot,
// sag(r), slope(r)) and with the first-guess paraboloid through the vertex and
// the outer rim. Test-only: it checks that later reformulations of
// Asphere::intersect keep the surface they find.
#include "obdeect/surfaces.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>

namespace legacy {

using namespace obdeect;

template <class Accept>
std::optional<SurfaceHit> intersect(const Asphere &a, const Ray &ray, double t_min, Accept accept,
                                    bool &failed) {
  failed = false;
  const double tolerance = 1e-10;
  const double shift = std::max(0.0, dot(Vec3{0.0, 0.0, a.vertex_z} - ray.origin, ray.direction));
  const Ray local{ray.at(shift), ray.direction};
  const double local_min = t_min - shift;
  const auto residual_at = [&](const Vec3 &p, double r) { return p.z - a.sag(r); };
  const auto derivative_at = [&](const Vec3 &p, double r) {
    const Vec3 &d = local.direction;
    return r > 0.0 ? d.z - a.slope(r) * (p.x * d.x + p.y * d.y) / r : d.z;
  };
  const auto residual = [&](double t) {
    const Vec3 p = local.at(t);
    return residual_at(p, std::hypot(p.x, p.y));
  };
  struct Root {
    double t;
    Vec3 p;
    double r;
  };
  const auto normal_at = [&](const Vec3 &p, double r) {
    if (r == 0.0)
      return Vec3{0.0, 0.0, 1.0};
    const double s = a.slope(r);
    return normalized(Vec3{-s * p.x / r, -s * p.y / r, 1.0});
  };
  const auto accept_root = [&](const Root &root) -> std::optional<SurfaceHit> {
    if (!(root.t > local_min))
      return std::nullopt;
    if (root.r < a.inner_radius || root.r > a.outer_radius || !accept(root.p))
      return std::nullopt;
    return SurfaceHit{root.t + shift, root.p, normal_at(root.p, root.r), a.id};
  };
  const auto root_at = [&](double t) {
    const Vec3 p = local.at(t);
    return Root{t, p, std::hypot(p.x, p.y)};
  };
  const auto newton = [&](double t, double lo, double hi) -> std::optional<Root> {
    for (int iteration = 0; iteration < 50; ++iteration) {
      const Vec3 p = local.at(t);
      const double r = std::hypot(p.x, p.y);
      const double g = residual_at(p, r);
      if (std::abs(g) <= tolerance)
        return Root{t, p, r};
      const double gd = derivative_at(p, r);
      if (gd == 0.0 || !std::isfinite(gd))
        return std::nullopt;
      const double next = t - g / gd;
      if (!(next >= lo && next <= hi) || next == t)
        return std::abs(g) <= 100 * tolerance ? std::optional<Root>{Root{t, p, r}} : std::nullopt;
      t = next;
    }
    return std::nullopt;
  };
  const Vec3 &o = local.origin, &d = local.direction;
  const double a2 = d.x * d.x + d.y * d.y, b2 = o.x * d.x + o.y * d.y;
  const double c2 = o.x * o.x + o.y * o.y - a.outer_radius * a.outer_radius;
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
  const double max_slope = a.slope_bound();
  const bool single_root = std::abs(d.z) > max_slope * std::sqrt(a2);
  const double z0 = a.sag(0.0);
  const double k = (a.sag(a.outer_radius) - z0) / (a.outer_radius * a.outer_radius);
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
  if (!std::isnan(guesses[1]) && (std::isnan(guesses[0]) || guesses[1] < guesses[0]))
    std::swap(guesses[0], guesses[1]);
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
      double lo_b = t0, hi_b = t1, ga = g0;
      for (int iteration = 0;
           iteration < 200 && hi_b - lo_b > 1e-13 * std::max(1.0, std::abs(hi_b)); ++iteration) {
        const double m = 0.5 * (lo_b + hi_b), gm = residual(m);
        if ((gm < 0.0) == (ga < 0.0)) {
          lo_b = m;
          ga = gm;
        } else {
          hi_b = m;
        }
      }
      if (auto hit = accept_root(root_at(0.5 * (lo_b + hi_b))))
        return hit;
    }
    t0 = t1;
    g0 = g1;
  }
  return std::nullopt;
}

} // namespace legacy
