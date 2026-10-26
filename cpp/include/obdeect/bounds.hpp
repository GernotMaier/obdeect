#pragma once

#include "obdeect/model.hpp"

#include <algorithm>
#include <cstddef>

namespace obdeect {

namespace detail {

// Upper bound of Response::at(wavelength, any angle): the angle interpolation
// is a convex combination of node values, so the largest node of the two
// interpolated wavelength rows bounds it. Same wavelength arithmetic as
// Response::at.
inline double response_bound(const Response &r, double wavelength_nm) {
  const double position = std::clamp((wavelength_nm - r.wavelength_min_nm) / r.wavelength_step_nm,
                                     0.0, static_cast<double>(r.wavelength_count - 1));
  const auto w = std::min(static_cast<std::size_t>(position), r.wavelength_count - 1);
  const double fw = position - static_cast<double>(w);
  const std::size_t w1 = std::min(w + 1, r.wavelength_count - 1);
  const std::size_t na = r.incidence_deg.size();
  double bound = r.values[w * na] * (1 - fw) + r.values[w1 * na] * fw;
  for (std::size_t j = 1; j < na; ++j)
    bound = std::max(bound, r.values[w * na + j] * (1 - fw) + r.values[w1 * na + j] * fw);
  return bound;
}

inline double response_bound(const Mirror &mirror, double wavelength_nm) {
  double bound = 1.0;
  if (mirror.reflectivity)
    bound *= response_bound(*mirror.reflectivity, wavelength_nm);
  // Bilinear interpolation never exceeds the largest node; outside the map the
  // value is either the edge value or 0.
  if (mirror.degradation)
    bound *=
        *std::max_element(mirror.degradation->values.begin(), mirror.degradation->values.end());
  return bound;
}

} // namespace detail

// Throws ModelError for a telescope transmission whose maximum is not known.
inline void validate_bound_support(const OpticalModel &model) {
  if (const auto &t = model.telescope_transmission)
    if (!(t->amplitude >= 0.0 && t->power > 0.0 && t->angular_scale_rad > 0.0 &&
          t->outer_power >= 0.0))
      throw ModelError("response_bound: unsupported telescope transmission");
}

// An upper bound of the weight factor that trace_photon with apply_responses
// can apply to a photon of this wavelength: reflectivities, degradation maps
// and the telescope transmission. Throws ModelError for a transmission whose
// maximum is not known.
[[nodiscard]] inline double response_bound(const OpticalModel &model, double wavelength_nm) {
  validate_bound_support(model);
  double bound = detail::response_bound(model.primary, wavelength_nm);
  if (model.secondary)
    bound *= detail::response_bound(*model.secondary, wavelength_nm);
  if (const auto &t = model.telescope_transmission)
    bound *= t->on_axis;
  return bound * (1.0 + 1e-12);
}

} // namespace obdeect
