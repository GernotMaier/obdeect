#pragma once

#include "obdeect/math.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

namespace obdeect {

enum class PhotonStatus : std::uint8_t {
  detected,
  blocked_camera,
  blocked_mast,
  missed_primary,
  missed_screen,
  invalid_input,
};

inline std::string_view to_string(PhotonStatus status) {
  switch (status) {
    case PhotonStatus::detected: return "detected";
    case PhotonStatus::blocked_camera: return "blocked_camera";
    case PhotonStatus::blocked_mast: return "blocked_mast";
    case PhotonStatus::missed_primary: return "missed_primary";
    case PhotonStatus::missed_screen: return "missed_screen";
    case PhotonStatus::invalid_input: return "invalid_input";
  }
  return "unknown";
}

// Non-owning SoA boundary suitable for a future NumPy/nanobind view. Every
// field must have the same length and direction vectors are validated once.
struct PhotonBlockView {
  std::span<const Vec3> position_m;
  std::span<const Vec3> direction;
  std::span<const double> wavelength_nm;
  std::span<const double> time_ns;
  std::span<const double> weight;
  std::span<const std::uint64_t> photon_id;

  [[nodiscard]] bool is_consistent() const {
    const std::size_t size = position_m.size();
    return direction.size() == size && wavelength_nm.size() == size && time_ns.size() == size &&
           weight.size() == size && photon_id.size() == size;
  }
};

struct PhotonResultBlock {
  std::vector<Vec3> position_m;
  std::vector<Vec3> direction;
  std::vector<double> optical_path_m;
  std::vector<double> time_ns;
  std::vector<double> weight;
  std::vector<PhotonStatus> status;

  explicit PhotonResultBlock(std::size_t size = 0)
      : position_m(size), direction(size), optical_path_m(size), time_ns(size), weight(size), status(size) {}
};

struct PathRecord {
  std::uint64_t photon_id{};
  double wavelength_nm{400.0};
  PhotonStatus status{PhotonStatus::missed_primary};
  std::array<Vec3, 3> points_m{};
  std::uint8_t point_count{};
  double path_length_m{};
};

}  // namespace obdeect
