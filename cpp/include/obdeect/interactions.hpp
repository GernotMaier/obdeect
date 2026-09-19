#pragma once

#include "obdeect/math.hpp"

#include <optional>
#include <algorithm>

namespace obdeect {

[[nodiscard]] inline std::optional<Vec3> reflect_specular(const Vec3& incident, const Vec3& normal) {
  const auto unit_incident = normalised_checked(incident);
  const auto unit_normal = normalised_checked(normal);
  if (!unit_incident || !unit_normal) return std::nullopt;
  return normalised_checked(*unit_incident - *unit_normal * (2.0 * dot(*unit_incident, *unit_normal)));
}

struct DielectricInterfaceResult {
  Vec3 reflected;
  std::optional<Vec3> transmitted;
  double reflectance_s{};
  double reflectance_p{};
};

// Lossless, isotropic dielectric boundary. The normal points into the
// incident medium; indices are phase indices at the photon wavelength.
// T_s,p = 1 - R_s,p is an energy fraction, not a field amplitude.
[[nodiscard]] inline std::optional<DielectricInterfaceResult> dielectric_interface(
    Vec3 incident, Vec3 normal, double n_incident, double n_transmitted) {
  const auto d = normalised_checked(incident);
  const auto n = normalised_checked(normal);
  if (!d || !n || !std::isfinite(n_incident) || !std::isfinite(n_transmitted) ||
      n_incident <= 0 || n_transmitted <= 0 || dot(*d, *n) > 0) return std::nullopt;
  const Vec3 reflected = *d - *n * (2 * dot(*d, *n));
  if (n_incident == n_transmitted)
    return DielectricInterfaceResult{reflected, *d, 0, 0};
  const double ci = std::clamp(-dot(*d, *n), 0.0, 1.0);
  const double eta = n_incident / n_transmitted;
  const double st2 = eta * eta * (1 - ci * ci);
  if (!std::isfinite(st2)) return std::nullopt;
  if (st2 >= 1) return DielectricInterfaceResult{reflected, std::nullopt, 1, 1};
  const double ct = std::sqrt(1 - st2);
  const auto transmitted = normalised_checked(*d * eta + *n * (eta * ci - ct));
  if (!transmitted) return std::nullopt;
  const double rs = (n_incident * ci - n_transmitted * ct) /
                    (n_incident * ci + n_transmitted * ct);
  const double rp = (n_transmitted * ci - n_incident * ct) /
                    (n_transmitted * ci + n_incident * ct);
  return DielectricInterfaceResult{reflected, transmitted, rs * rs, rp * rp};
}

[[nodiscard]] inline std::optional<double> bulk_transmission(double distance_m,
                                                            double absorption_per_m) {
  if (!std::isfinite(distance_m) || !std::isfinite(absorption_per_m) ||
      distance_m < 0 || absorption_per_m < 0) return std::nullopt;
  return std::exp(-absorption_per_m * distance_m);
}

}  // namespace obdeect
