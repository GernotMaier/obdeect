#pragma once

#include "obdeect/photon_input.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace obdeect {

// Direct-beam optical depth tabulated against wavelength and altitude above
// sea level. The first altitude is the observation level, where depth is zero.
// It describes removal from the direct beam (absorption + out-scattering),
// not photons scattered into the telescope.
class AtmosphereTransmissionTable {
 public:
  AtmosphereTransmissionTable(std::vector<double> wavelength_nm,
                              std::vector<double> altitude_m,
                              std::vector<double> optical_depth)
      : wavelength_nm_(std::move(wavelength_nm)), altitude_m_(std::move(altitude_m)),
        optical_depth_(std::move(optical_depth)) {
    if (wavelength_nm_.size() < 2 || altitude_m_.size() < 2 ||
        optical_depth_.size() != wavelength_nm_.size() * altitude_m_.size())
      throw std::invalid_argument("invalid atmospheric optical-depth grid shape");
    for (std::size_t i = 0; i < wavelength_nm_.size(); ++i)
      if (!std::isfinite(wavelength_nm_[i]) || wavelength_nm_[i] <= 0 ||
          (i && wavelength_nm_[i] <= wavelength_nm_[i - 1]))
        throw std::invalid_argument("invalid atmospheric wavelength grid");
    for (std::size_t i = 0; i < altitude_m_.size(); ++i)
      if (!std::isfinite(altitude_m_[i]) || altitude_m_[i] <= 0 ||
          (i && altitude_m_[i] <= altitude_m_[i - 1]))
        throw std::invalid_argument("invalid atmospheric altitude grid");
    for (std::size_t w = 0; w < wavelength_nm_.size(); ++w) {
      if (optical_depth_[w * altitude_m_.size()] != 0.0)
        throw std::invalid_argument("optical depth at observation level must be zero");
      for (std::size_t h = 0; h < altitude_m_.size(); ++h) {
        const double value = optical_depth_[w * altitude_m_.size() + h];
        if (!std::isfinite(value) || value < 0 ||
            (h && value + 1e-10 < optical_depth_[w * altitude_m_.size() + h - 1]))
          throw std::invalid_argument("invalid atmospheric optical depth");
      }
    }
  }

  // Read the MODTRAN-derived table format used by sim_telarray. Data values
  // are optical depth, not transmission; 99999 is an invalid-grid sentinel.
  [[nodiscard]] static AtmosphereTransmissionTable from_simtel_file(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::invalid_argument("cannot open atmospheric table: " + path);
    std::vector<double> altitude_m, wavelength_nm, optical_depth;
    std::string line;
    while (std::getline(input, line)) {
      if (line.rfind("# H2=", 0) == 0) {
        const auto h1 = line.find("H1=");
        if (h1 == std::string::npos || !altitude_m.empty())
          throw std::invalid_argument("invalid atmospheric table altitude header");
        const auto h2_text = line.substr(5, h1 - 5);
        altitude_m.push_back(1000.0 * std::stod(h2_text));
        std::istringstream heights(line.substr(h1 + 3));
        double height_km = 0;
        while (heights >> height_km) altitude_m.push_back(1000.0 * height_km);
      } else if (!line.empty() && line.front() != '#') {
        if (altitude_m.size() < 2)
          throw std::invalid_argument("atmospheric table data precedes altitude header");
        std::istringstream row(line);
        double wavelength = 0;
        if (!(row >> wavelength)) throw std::invalid_argument("invalid atmospheric wavelength row");
        wavelength_nm.push_back(wavelength);
        optical_depth.push_back(0.0);
        double value = 0;
        for (std::size_t h = 1; h < altitude_m.size(); ++h) {
          if (!(row >> value)) throw std::invalid_argument("incomplete atmospheric wavelength row");
          optical_depth.push_back(value);
        }
        if (row >> value) throw std::invalid_argument("extra atmospheric depth column");
      }
    }
    if (input.bad()) throw std::invalid_argument("error reading atmospheric table");
    return {std::move(wavelength_nm), std::move(altitude_m), std::move(optical_depth)};
  }

  [[nodiscard]] double observation_altitude_m() const noexcept { return altitude_m_.front(); }

  [[nodiscard]] double depth(double wavelength_nm, double altitude_m) const {
    if (!std::isfinite(wavelength_nm) || !std::isfinite(altitude_m) ||
        wavelength_nm < wavelength_nm_.front() || wavelength_nm > wavelength_nm_.back() ||
        altitude_m < altitude_m_.front() || altitude_m > altitude_m_.back())
      throw std::out_of_range("atmospheric wavelength or altitude outside table");
    const auto wi = lower_interval(wavelength_nm_, wavelength_nm);
    const auto hi = lower_interval(altitude_m_, altitude_m);
    const double w = (wavelength_nm - wavelength_nm_[wi]) /
                     (wavelength_nm_[wi + 1] - wavelength_nm_[wi]);
    // Match the reference interpolation coordinate for altitude. In the first
    // interval the lower endpoint remains the observation altitude.
    const double h = (std::log(altitude_m) - std::log(altitude_m_[hi])) /
                     (std::log(altitude_m_[hi + 1]) - std::log(altitude_m_[hi]));
    const auto corner = [&](std::size_t x, std::size_t y) {
      const double value = optical_depth_[x * altitude_m_.size() + y];
      if (value >= 99999.0) throw std::out_of_range("invalid atmospheric table cell");
      return value;
    };
    const double lower = std::lerp(corner(wi, hi), corner(wi + 1, hi), w);
    const double upper = std::lerp(corner(wi, hi + 1), corner(wi + 1, hi + 1), w);
    return std::lerp(lower, upper, h);
  }

  [[nodiscard]] double direct_survival(const OpticalPhoton& photon,
                                       const PhotonBatchContext& context) const {
    if (!valid_photon(photon) || !valid_context(context) || photon.wavelength_nm <= 0)
      throw std::invalid_argument("direct extinction needs a resolved physical photon");
    if (photon.ray.direction.z >= -0.2)
      throw std::invalid_argument("plane-parallel extinction is invalid for near-horizontal/upward rays");
    const double arrival_altitude = observation_altitude_m() + context.telescope_position_m.z +
                                    photon.ray.position_m.z;
    const bool has_height = std::isfinite(photon.emission_height_m);
    const bool has_distance = std::isfinite(photon.emission_distance_m);
    if (has_height == has_distance)
      throw std::invalid_argument("exactly one emission height or distance is required");
    const double emission_altitude = has_height
        ? photon.emission_height_m
        : arrival_altitude - photon.ray.direction.z * photon.emission_distance_m;
    if (emission_altitude < arrival_altitude)
      throw std::invalid_argument("emission below arrival requires a path-integral model");
    const double delta = depth(photon.wavelength_nm, emission_altitude) -
                         depth(photon.wavelength_nm, arrival_altitude);
    if (delta < -1e-10) throw std::runtime_error("negative atmospheric optical-depth difference");
    return std::exp(-std::max(0.0, delta) / -photon.ray.direction.z);
  }

 private:
  [[nodiscard]] static std::size_t lower_interval(const std::vector<double>& grid, double value) {
    const auto upper = std::upper_bound(grid.begin(), grid.end(), value);
    return upper == grid.begin() ? 0 :
           std::min(static_cast<std::size_t>(upper - grid.begin() - 1), grid.size() - 2);
  }

  std::vector<double> wavelength_nm_, altitude_m_, optical_depth_;
};

}  // namespace obdeect
