#pragma once

#include "obdeect/model.hpp"
#include "obdeect/random.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string_view>
#include <thread>
#include <vector>

namespace obdeect {

enum class Status : std::uint8_t {
  detected,             // Reached the focal surface.
  missed_primary,       // Never met the primary mirror.
  missed_secondary,     // Left the primary without meeting the secondary.
  missed_focal_surface, // Never met the focal surface.
  blocked,              // Absorbed by an obscurer (terminal_surface names it).
  invalid_input,        // Non-finite input, wavelength <= 0 or non-unit direction.
  numerical_failure,    // An intersection did not converge.
};
inline constexpr std::size_t kStatusCount = 7;

[[nodiscard]] constexpr std::string_view to_string(Status status) {
  constexpr std::string_view names[kStatusCount] = {
      "detected", "missed_primary", "missed_secondary", "missed_focal_surface",
      "blocked",  "invalid_input",  "numerical_failure"};
  return names[static_cast<std::size_t>(status)];
}

inline constexpr std::uint32_t kNoSurface = 0xFFFFFFFFU;

// Result of one photon. Mirror hits not reached keep NaN coordinates and
// kNoSurface; `weight` is the input weight times the applied responses.
struct TraceResult {
  Status status{Status::invalid_input};
  std::uint32_t terminal_surface{kNoSurface};
  Vec3 position, direction;
  double path_length_m{}, time_ns{}, weight{};
  std::uint32_t primary_surface{kNoSurface}, secondary_surface{kNoSurface};
  Vec3 primary_point{nan(), nan(), nan()}, secondary_point{nan(), nan(), nan()};
  double primary_incidence_deg{nan()}, secondary_incidence_deg{nan()}, focal_incidence_deg{nan()};

  static constexpr double nan() { return std::numeric_limits<double>::quiet_NaN(); }
};

namespace detail {

constexpr double kDegree = 57.29577951308232087680;

inline double incidence_deg(const Vec3 &direction, const Vec3 &normal) {
  return std::acos(std::min(1.0, std::abs(dot(direction, normal)))) * kDegree;
}

struct MirrorHit {
  SurfaceHit hit;
  std::uint32_t facet_index{};
};

// Nearest hit on a mirror (lowest id on exact ties). `failed` reports a
// non-converged aspheric intersection.
inline std::optional<MirrorHit> intersect_mirror(const Mirror &mirror, const Ray &ray,
                                                 bool &failed) {
  failed = false;
  if (mirror.asphere) {
    const auto &a = *mirror.asphere;
    const auto inside = [&](std::size_t i, const Vec3 &p) {
      return mirror.segments[i].contains(mirror.segment_frames[i], p);
    };
    const auto accept = [&](const Vec3 &p) {
      if (mirror.segments.empty())
        return true;
      for (std::size_t i = 0; i < mirror.segments.size(); ++i)
        if (inside(i, p))
          return true;
      return false;
    };
    auto hit = a.intersect(ray, kMinDistanceM, accept, failed);
    if (!hit)
      return std::nullopt;
    for (std::size_t i = 0; i < mirror.segments.size(); ++i) // Lowest segment id at the point.
      if ((hit->id == a.id || mirror.segments[i].id < hit->id) && inside(i, hit->point))
        hit->id = mirror.segments[i].id;
    return MirrorHit{*hit, 0};
  }
  std::optional<MirrorHit> nearest;
  mirror.grid.candidates(ray, mirror.facets.size(), [&](std::size_t i) {
    const auto &f = mirror.facets[i];
    if (const auto hit = f.intersect(ray, kMinDistanceM))
      if (!nearest || hit->t < nearest->hit.t ||
          (hit->t == nearest->hit.t && hit->id < nearest->hit.id))
        nearest = MirrorHit{*hit, static_cast<std::uint32_t>(i)};
  });
  return nearest;
}

// Rotation of v about the unit axis k by `angle` (Rodrigues).
inline Vec3 rotate(const Vec3 &v, const Vec3 &k, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  return v * c + cross(k, v) * s + k * (dot(k, v) * (1.0 - c));
}

// Reflected direction, with the mirror's random reflection error.
inline Vec3 reflect_at(const Mirror &mirror, const MirrorHit &mh, const Vec3 &direction,
                       std::uint64_t photon_id, std::uint64_t seed, Stream stream) {
  const Vec3 &n = mh.hit.normal;
  if (!mirror.scatter)
    return reflect(direction, n);
  const auto &s = *mirror.scatter;
  const auto u = uniform4(seed, photon_id, stream);
  const double sigma = u[0] < s.fraction2 ? s.sigma2_rad : s.sigma1_rad;
  if (sigma == 0.0)
    return reflect(direction, n);
  const double radius = sigma * std::sqrt(-2.0 * std::log(u[1]));
  const double phi = 6.283185307179586476925 * u[2];
  const double a1 = radius * std::cos(phi), a2 = radius * std::sin(phi);
  if (s.surface_slopes) {
    // Graph-surface normal (-dz/dx, -dz/dy, 1): perturb each slope by half
    // the angle, so the reflected ray deviates by about the full angle.
    const Vec3 tilted{n.x / n.z + 0.5 * a1, n.y / n.z + 0.5 * a2, 1.0};
    return reflect(direction, normalized(tilted));
  }
  // Rotate the reflected ray about two axes in the mirror's tangent plane.
  const Vec3 e1 = mirror.facets.empty() ? normalized(std::abs(n.z) < 0.9 ? cross(Vec3{0, 0, 1}, n)
                                                                         : cross(Vec3{1, 0, 0}, n))
                                        : mirror.facets[mh.facet_index].u;
  const Vec3 e2 = cross(n, e1);
  return normalized(rotate(rotate(reflect(direction, n), e2, a1), e1, a2));
}

inline double mirror_response(const Mirror &mirror, const MirrorHit &mh, double wavelength_nm,
                              double incidence) {
  double r = mirror.reflectivity ? mirror.reflectivity->at(wavelength_nm, incidence) : 1.0;
  if (mirror.degradation) {
    const Vec3 &p = mh.hit.point;
    if (mirror.facets.empty()) {
      r *= mirror.degradation->at(p.x, p.y);
    } else {
      const auto &f = mirror.facets[mh.facet_index];
      r *= mirror.degradation->at(dot(p - f.centre, f.map_u), dot(p - f.centre, f.map_v));
    }
  }
  return r;
}

} // namespace detail

// Trace one photon through the optical model. With `apply_responses` false
// the weight stays at its input value (geometry-only transport). A `direct`
// photon starts inside the telescope on the leg to the focal surface (after
// the last mirror), e.g. from a flasher: no mirror and no telescope
// transmission act, only the obscurers of that leg.
[[nodiscard]] inline TraceResult trace_photon(const OpticalModel &model, Ray ray,
                                              double wavelength_nm, double time_ns, double weight,
                                              std::uint64_t photon_id, std::uint64_t seed,
                                              bool apply_responses = true, bool direct = false) {
  using namespace detail;
  TraceResult out;
  out.position = ray.origin;
  out.direction = ray.direction;
  const double length = norm(ray.direction);
  if (!finite(ray.origin) || !std::isfinite(length) || std::abs(length - 1.0) > 1e-6 ||
      !(wavelength_nm > 0.0) || !std::isfinite(wavelength_nm) || !std::isfinite(time_ns) ||
      !(weight >= 0.0) || !std::isfinite(weight))
    return out;
  ray.direction = ray.direction * (1.0 / length);
  out.weight = weight;
  double path = 0.0;
  const auto finish = [&](Status status, const Vec3 &position, std::uint32_t surface) {
    out.status = status;
    out.position = position;
    out.direction = ray.direction;
    out.terminal_surface = surface;
    out.path_length_m = path;
    out.time_ns = time_ns + path * model.propagation_group_index / kSpeedOfLightMPerNs;
    if (status != Status::detected)
      out.weight = 0.0;
    return out;
  };
  // A blocked photon stops at the obscurer; it is absorbed there.
  const auto blocked_on = [&](Leg leg, double t_limit) -> std::optional<TraceResult> {
    if (const auto b = model.obscurers.nearest(ray, leg, t_limit)) {
      path += b->t;
      return finish(Status::blocked, ray.at(b->t), b->id);
    }
    return std::nullopt;
  };

  bool failed = false;
  if (!direct) {
    if (apply_responses && model.telescope_transmission)
      out.weight *= model.telescope_transmission->at(ray.direction);

    // Source to primary mirror.
    const auto m1 = intersect_mirror(model.primary, ray, failed);
    if (failed)
      return finish(Status::numerical_failure, ray.origin, kNoSurface);
    double t_m1 = m1 ? m1->hit.t : std::numeric_limits<double>::infinity();
    if (model.secondary_blocks_incoming) {
      const auto back = intersect_mirror(*model.secondary, ray, failed);
      if (failed)
        return finish(Status::numerical_failure, ray.origin, kNoSurface);
      if (back && back->hit.t < t_m1) {
        if (auto b = blocked_on(Leg::incoming, back->hit.t))
          return *b;
        path += back->hit.t;
        return finish(Status::blocked, back->hit.point, back->hit.id);
      }
    }
    if (auto b = blocked_on(Leg::incoming, t_m1))
      return *b;
    if (!m1)
      return finish(Status::missed_primary, ray.origin, kNoSurface);

    path += m1->hit.t;
    out.primary_point = m1->hit.point;
    out.primary_surface = m1->hit.id;
    out.primary_incidence_deg = incidence_deg(ray.direction, m1->hit.normal);
    if (apply_responses)
      out.weight *= mirror_response(model.primary, *m1, wavelength_nm, out.primary_incidence_deg);
    ray = {m1->hit.point,
           reflect_at(model.primary, *m1, ray.direction, photon_id, seed, Stream::primary_scatter)};

    if (model.secondary) {
      const auto m2 = intersect_mirror(*model.secondary, ray, failed);
      if (failed)
        return finish(Status::numerical_failure, ray.origin, kNoSurface);
      if (auto b = blocked_on(Leg::after_primary,
                              m2 ? m2->hit.t : std::numeric_limits<double>::infinity()))
        return *b;
      if (!m2)
        return finish(Status::missed_secondary, ray.origin, kNoSurface);
      path += m2->hit.t;
      out.secondary_point = m2->hit.point;
      out.secondary_surface = m2->hit.id;
      out.secondary_incidence_deg = incidence_deg(ray.direction, m2->hit.normal);
      if (apply_responses)
        out.weight *=
            mirror_response(*model.secondary, *m2, wavelength_nm, out.secondary_incidence_deg);
      ray = {m2->hit.point, reflect_at(*model.secondary, *m2, ray.direction, photon_id, seed,
                                       Stream::secondary_scatter)};
    }
  }
  const Leg leg = model.secondary ? Leg::after_secondary : Leg::after_primary;

  // To the focal surface.
  std::optional<SurfaceHit> focal;
  const auto &fs = model.focal_surface;
  if (fs.asphere) {
    focal = fs.asphere->intersect(ray, kMinDistanceM, [](const Vec3 &) { return true; }, failed);
    if (failed)
      return finish(Status::numerical_failure, ray.origin, kNoSurface);
  } else if (ray.direction.z != 0.0) {
    const double t = (fs.plane_z - ray.origin.z) / ray.direction.z;
    if (t > kMinDistanceM)
      focal = SurfaceHit{t, ray.at(t), Vec3{0.0, 0.0, 1.0}, fs.id};
  }
  if (auto b = blocked_on(leg, focal ? focal->t : std::numeric_limits<double>::infinity()))
    return *b;
  if (!focal)
    return finish(Status::missed_focal_surface, ray.origin, kNoSurface);
  path += focal->t;
  out.focal_incidence_deg = incidence_deg(ray.direction, focal->normal);
  return finish(Status::detected, focal->point, fs.id);
}

// Structure-of-arrays photon block. Vectors are (n, 3) row-major.
struct PhotonBlock {
  const double *position{}, *direction{}, *wavelength_nm{}, *time_ns{}, *weight{};
  const std::uint64_t *photon_id{};
  std::size_t size{};
};

struct ResultBlock {
  std::uint8_t *status{};
  std::uint32_t *terminal_surface{}, *primary_surface{}, *secondary_surface{};
  double *position{}, *direction{}, *path_length_m{}, *time_ns{}, *weight{};
  double *primary_point{}, *secondary_point{};
  double *primary_incidence_deg{}, *secondary_incidence_deg{}, *focal_incidence_deg{};
};

// Run body(begin, end) over [0, n) in contiguous chunks on up to `threads`
// threads (at least 4096 items each).
template <class Body> void parallel_for(std::size_t n, unsigned threads, const Body &body) {
  constexpr std::size_t kMinimumPerThread = 4096;
  const std::size_t workers = std::max<std::size_t>(
      1, std::min<std::size_t>(threads, (n + kMinimumPerThread - 1) / kMinimumPerThread));
  if (workers == 1) {
    body(std::size_t{0}, n);
    return;
  }
  std::vector<std::thread> pool;
  const std::size_t chunk = (n + workers - 1) / workers;
  for (std::size_t w = 0; w < workers; ++w)
    pool.emplace_back(body, w * chunk, std::min(n, (w + 1) * chunk));
  for (auto &t : pool)
    t.join();
}

// Trace a block on `threads` threads. Every photon writes only its own row,
// and its random numbers depend only on (seed, photon_id), so results do not
// depend on the thread count or the order of photons.
inline void trace_block(const OpticalModel &model, const PhotonBlock &in, const ResultBlock &out,
                        std::uint64_t seed, bool apply_responses, unsigned threads,
                        bool direct = false) {
  parallel_for(in.size, threads, [&](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      const Ray ray{{in.position[3 * i], in.position[3 * i + 1], in.position[3 * i + 2]},
                    {in.direction[3 * i], in.direction[3 * i + 1], in.direction[3 * i + 2]}};
      const auto r = trace_photon(model, ray, in.wavelength_nm[i], in.time_ns[i], in.weight[i],
                                  in.photon_id[i], seed, apply_responses, direct);
      const auto put = [i](double *target, const Vec3 &v) {
        target[3 * i] = v.x;
        target[3 * i + 1] = v.y;
        target[3 * i + 2] = v.z;
      };
      out.status[i] = static_cast<std::uint8_t>(r.status);
      out.terminal_surface[i] = r.terminal_surface;
      out.primary_surface[i] = r.primary_surface;
      out.secondary_surface[i] = r.secondary_surface;
      put(out.position, r.position);
      put(out.direction, r.direction);
      put(out.primary_point, r.primary_point);
      put(out.secondary_point, r.secondary_point);
      out.path_length_m[i] = r.path_length_m;
      out.time_ns[i] = r.time_ns;
      out.weight[i] = r.weight;
      out.primary_incidence_deg[i] = r.primary_incidence_deg;
      out.secondary_incidence_deg[i] = r.secondary_incidence_deg;
      out.focal_incidence_deg[i] = r.focal_incidence_deg;
    }
  });
}

} // namespace obdeect
