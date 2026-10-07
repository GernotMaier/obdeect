#pragma once

#include "obdeect/tables.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

namespace obdeect {

struct SpectralResponse {
  std::vector<double> wavelength_nm;
  // Wavelength-major response values; an empty angle axis denotes a 1D table.
  std::vector<double> response;
  std::vector<double> incidence_angle_deg{};

  [[nodiscard]] bool is_valid() const {
    const std::size_t angles = incidence_angle_deg.empty() ? 1 : incidence_angle_deg.size();
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
  [[nodiscard]] std::optional<double> at_unchecked(double wavelength, double angle_deg = 0) const {
    if (!std::isfinite(wavelength) || wavelength < wavelength_nm.front() ||
        wavelength > wavelength_nm.back() || !std::isfinite(angle_deg) || angle_deg < 0 ||
        angle_deg > 90)
      return std::nullopt;
    const auto bracket = [](const std::vector<double> &axis, double value) {
      const auto upper = std::upper_bound(axis.begin(), axis.end(), value);
      return std::min(static_cast<std::size_t>(upper - axis.begin()), axis.size() - 1);
    };
    const auto w = bracket(wavelength_nm, wavelength);
    const double wf =
        (wavelength - wavelength_nm[w - 1]) / (wavelength_nm[w] - wavelength_nm[w - 1]);
    const std::size_t angles = incidence_angle_deg.empty() ? 1 : incidence_angle_deg.size();
    const auto spectral = [&](std::size_t a) {
      return response[(w - 1) * angles + a] * (1 - wf) + response[w * angles + a] * wf;
    };
    if (incidence_angle_deg.empty())
      return spectral(0);
    if (angle_deg < incidence_angle_deg.front() || angle_deg > incidence_angle_deg.back())
      return std::nullopt;
    const auto a = bracket(incidence_angle_deg, angle_deg);
    const double af = (angle_deg - incidence_angle_deg[a - 1]) /
                      (incidence_angle_deg[a] - incidence_angle_deg[a - 1]);
    return spectral(a - 1) * (1 - af) + spectral(a) * af;
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

  [[nodiscard]] Table1DView view() const { return {incidence_angle_deg, response}; }
  [[nodiscard]] bool is_valid() const {
    return view().is_valid() && incidence_angle_deg.front() >= 0 &&
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
      const auto guide = lightguide_efficiency->view().interpolate_unchecked(incidence_angle_deg);
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
