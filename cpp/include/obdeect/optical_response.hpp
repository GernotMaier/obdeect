#pragma once

#include "obdeect/interpolation.hpp"
#include "obdeect/tables.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <vector>

namespace obdeect {

struct WavelengthSampling {
  double width_nm{}, offset_nm{};
  std::uint32_t first_bin{}, last_bin{};
  std::optional<double> projection_angle_deg{};

  [[nodiscard]] bool is_valid() const {
    return std::isfinite(width_nm) && width_nm > 0 && std::isfinite(offset_nm) && offset_nm >= 0 &&
           offset_nm < width_nm && first_bin <= last_bin &&
           (!projection_angle_deg || (std::isfinite(*projection_angle_deg) &&
                                      *projection_angle_deg >= 0 && *projection_angle_deg <= 90));
  }
};

struct SpectralResponse {
  std::vector<double> wavelength_nm;
  // Wavelength-major response values; an empty angle axis denotes a 1D table.
  std::vector<double> response;
  std::vector<double> incidence_angle_deg{};
  TableInterpolation interpolation{};
  std::optional<WavelengthSampling> wavelength_sampling{};
  std::vector<double> spectral_envelope{};
  TableInterpolation envelope_interpolation{};
  bool relative_to_envelope{};

  [[nodiscard]] bool is_valid() const {
    if (relative_to_envelope &&
        (incidence_angle_deg.empty() || spectral_envelope.size() != wavelength_nm.size()))
      return false;
    if (wavelength_sampling &&
        (!wavelength_sampling->is_valid() ||
         (!incidence_angle_deg.empty() &&
          (spectral_envelope.size() != wavelength_nm.size() ||
           !envelope_interpolation.is_valid(wavelength_nm.empty() ? 0 : wavelength_nm.size() - 1) ||
           std::any_of(spectral_envelope.begin(), spectral_envelope.end(), [](double value) {
             return !std::isfinite(value) || value < 0 || value > 1;
           })))))
      return false;
    const std::size_t angles = incidence_angle_deg.empty() ? 1 : incidence_angle_deg.size();
    if (!interpolation.is_valid(wavelength_nm.empty() ? 0 : wavelength_nm.size() - 1,
                                !incidence_angle_deg.empty()))
      return false;
    if (interpolation.y_log && (incidence_angle_deg.empty() || incidence_angle_deg.front() <= 0))
      return false;
    if (interpolation.value_log &&
        std::any_of(response.begin(), response.end(), [](double v) { return v <= 0; }))
      return false;
    if (wavelength_nm.size() < 2 || (!incidence_angle_deg.empty() && angles < 2) ||
        response.size() / angles != wavelength_nm.size() || response.size() % angles != 0)
      return false;
    for (std::size_t i = 0; i < wavelength_nm.size(); ++i)
      if (!std::isfinite(wavelength_nm[i]) || wavelength_nm[i] <= 0 ||
          (i && wavelength_nm[i] <= wavelength_nm[i - 1]))
        return false;
    for (std::size_t i = 0; i < incidence_angle_deg.size(); ++i)
      if (!std::isfinite(incidence_angle_deg[i]) || incidence_angle_deg[i] < 0 ||
          incidence_angle_deg[i] > 90 ||
          (i && incidence_angle_deg[i] <= incidence_angle_deg[i - 1]))
        return false;
    return std::all_of(response.begin(), response.end(), [](double value) {
      return std::isfinite(value) && value >= 0 && value <= 1;
    });
  }

  // The caller has validated the immutable table once, before transport.
  [[nodiscard]] std::optional<double> at_continuous_unchecked(double wavelength,
                                                              double angle_deg = 0) const {
    if (!std::isfinite(wavelength) || wavelength <= 0 || !std::isfinite(angle_deg) ||
        angle_deg < 0 || angle_deg > 90)
      return std::nullopt;
    if (incidence_angle_deg.empty())
      return interpolate_curve_unchecked(wavelength_nm, response, wavelength, interpolation);
    if (wavelength < wavelength_nm.front() || wavelength > wavelength_nm.back() ||
        angle_deg < incidence_angle_deg.front() || angle_deg > incidence_angle_deg.back()) {
      if (interpolation.boundary == TableBoundary::reject)
        return std::nullopt;
      if (interpolation.boundary == TableBoundary::zero)
        return 0;
      wavelength = std::clamp(wavelength, wavelength_nm.front(), wavelength_nm.back());
      angle_deg = std::clamp(angle_deg, incidence_angle_deg.front(), incidence_angle_deg.back());
    }
    const auto bracket = [](const std::vector<double> &axis, double value) {
      const auto upper = std::upper_bound(axis.begin(), axis.end(), value);
      return std::min(static_cast<std::size_t>(upper - axis.begin()), axis.size() - 1);
    };
    const auto w = bracket(wavelength_nm, wavelength);
    const double wf =
        table_fraction(wavelength, wavelength_nm[w - 1], wavelength_nm[w], interpolation.x_log);
    const std::size_t angles = incidence_angle_deg.empty() ? 1 : incidence_angle_deg.size();
    const auto spectral = [&](std::size_t a) {
      const auto value = [&](std::size_t i) {
        return interpolation.value_log ? std::log(response[i]) : response[i];
      };
      return value((w - 1) * angles + a) * (1 - wf) + value(w * angles + a) * wf;
    };
    if (incidence_angle_deg.empty())
      return spectral(0);
    if (angle_deg < incidence_angle_deg.front() || angle_deg > incidence_angle_deg.back())
      return std::nullopt;
    const auto a = bracket(incidence_angle_deg, angle_deg);
    const double af = table_fraction(angle_deg, incidence_angle_deg[a - 1], incidence_angle_deg[a],
                                     interpolation.y_log);
    const double value = spectral(a - 1) * (1 - af) + spectral(a) * af;
    return interpolation.value_log ? std::exp(value) : value;
  }

  [[nodiscard]] std::optional<double> at_unchecked(double wavelength, double angle_deg = 0) const {
    if (relative_to_envelope) {
      const auto physical = at_continuous_unchecked(wavelength, angle_deg);
      const auto envelope = interpolate_curve_unchecked(wavelength_nm, spectral_envelope,
                                                        wavelength, envelope_interpolation);
      if (!physical || !envelope)
        return std::nullopt;
      return *envelope > 0 ? *physical / *envelope : 0;
    }
    if (!wavelength_sampling)
      return at_continuous_unchecked(wavelength, angle_deg);
    if (!std::isfinite(wavelength) || wavelength <= 0 || !std::isfinite(angle_deg) ||
        angle_deg < 0 || angle_deg > 90)
      return std::nullopt;
    const auto &sampling = *wavelength_sampling;
    const double bin = std::floor(wavelength / sampling.width_nm);
    if (bin < sampling.first_bin || bin > sampling.last_bin)
      return 0;
    const double sampled = bin * sampling.width_nm + sampling.offset_nm;
    if (incidence_angle_deg.empty())
      return at_continuous_unchecked(sampled, angle_deg);
    const auto physical = at_continuous_unchecked(wavelength, angle_deg);
    const auto nominal = interpolate_curve_unchecked(wavelength_nm, spectral_envelope, wavelength,
                                                     envelope_interpolation);
    const auto upstream = sampling.projection_angle_deg
                              ? at_continuous_unchecked(sampled, *sampling.projection_angle_deg)
                              : interpolate_curve_unchecked(wavelength_nm, spectral_envelope,
                                                            sampled, envelope_interpolation);
    if (!physical || !nominal || !upstream)
      return std::nullopt;
    return *nominal > 0 ? *physical * *upstream / *nominal : 0;
  }

  [[nodiscard]] std::optional<double> at(double wavelength, double angle_deg = 0) const {
    return is_valid() ? at_unchecked(wavelength, angle_deg) : std::nullopt;
  }
};

// Measured angular concentrator efficiency, distinct from simulated guide
// walls. Owning arrays live in the immutable compiled optical model.
struct CameraIncidenceResponse {
  std::vector<double> incidence_angle_deg;
  std::vector<double> response;
  TableInterpolation interpolation{};

  [[nodiscard]] Table1DView view() const { return {incidence_angle_deg, response}; }
  [[nodiscard]] bool is_valid() const {
    return interpolation.is_valid(incidence_angle_deg.empty() ? 0
                                                              : incidence_angle_deg.size() - 1) &&
           (!interpolation.x_log ||
            (!incidence_angle_deg.empty() && incidence_angle_deg.front() > 0)) &&
           (!interpolation.value_log ||
            std::all_of(response.begin(), response.end(), [](double v) { return v > 0; })) &&
           view().is_valid() && incidence_angle_deg.front() >= 0 &&
           incidence_angle_deg.back() <= 90 &&
           std::all_of(response.begin(), response.end(), [](double value) { return value <= 1; });
  }
};

// Complete measured entrance response. Apply this product exactly once at the
// optical detector boundary; do not also simulate the same filter/guide loss.
// The file loader requires explicit measured_complete_response semantics.
struct CameraResponse {
  double camera_transmission{1};
  std::optional<SpectralResponse> camera_filter;
  std::optional<CameraIncidenceResponse> lightguide_efficiency;

  [[nodiscard]] bool is_valid() const {
    return std::isfinite(camera_transmission) && camera_transmission >= 0 &&
           camera_transmission <= 1 && (!camera_filter || camera_filter->is_valid()) &&
           (!lightguide_efficiency || lightguide_efficiency->is_valid());
  }

  // Model validity is established once during compilation. Querying the
  // immutable response is allocation-free and rejects all extrapolation.
  [[nodiscard]] std::optional<double> at_unchecked(double wavelength_nm,
                                                   double incidence_angle_deg) const {
    if (!std::isfinite(wavelength_nm) || wavelength_nm <= 0 ||
        !std::isfinite(incidence_angle_deg) || incidence_angle_deg < 0 || incidence_angle_deg > 90)
      return std::nullopt;
    double transmission = camera_transmission;
    if (camera_filter) {
      const auto filter = camera_filter->at_unchecked(wavelength_nm, incidence_angle_deg);
      if (!filter)
        return std::nullopt;
      transmission *= *filter;
    }
    if (lightguide_efficiency) {
      const auto guide = interpolate_curve_unchecked(
          lightguide_efficiency->incidence_angle_deg, lightguide_efficiency->response,
          incidence_angle_deg, lightguide_efficiency->interpolation);
      if (!guide)
        return std::nullopt;
      transmission *= *guide;
    }
    return transmission;
  }

  [[nodiscard]] std::optional<double> at(double wavelength_nm, double incidence_angle_deg) const {
    return is_valid() ? at_unchecked(wavelength_nm, incidence_angle_deg) : std::nullopt;
  }
};

} // namespace obdeect
