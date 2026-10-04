#pragma once

#include "obdeect/math.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string_view>
#include <vector>

namespace obdeect {

enum class PhotonStatus : std::uint8_t {
  detected,
  blocked_camera,
  blocked_mast,
  missed_primary,
  missed_secondary,
  missed_screen,
  no_detector,
  invalid_input,
  blocked_obscurer,
  escaped_optical_model,
  interaction_limit,
  intersection_failure,
  absorbed_material,
  escaped_material,
  count,
};

constexpr std::size_t kPhotonStatusCount = static_cast<std::size_t>(PhotonStatus::count);

inline std::string_view to_string(PhotonStatus status) {
  switch (status) {
  case PhotonStatus::detected:
    return "detected";
  case PhotonStatus::blocked_camera:
    return "blocked_camera";
  case PhotonStatus::blocked_mast:
    return "blocked_mast";
  case PhotonStatus::missed_primary:
    return "missed_primary";
  case PhotonStatus::missed_secondary:
    return "missed_secondary";
  case PhotonStatus::missed_screen:
    return "missed_screen";
  case PhotonStatus::no_detector:
    return "no_detector";
  case PhotonStatus::invalid_input:
    return "invalid_input";
  case PhotonStatus::blocked_obscurer:
    return "blocked_obscurer";
  case PhotonStatus::escaped_optical_model:
    return "escaped_optical_model";
  case PhotonStatus::interaction_limit:
    return "interaction_limit";
  case PhotonStatus::intersection_failure:
    return "intersection_failure";
  case PhotonStatus::absorbed_material:
    return "absorbed_material";
  case PhotonStatus::escaped_material:
    return "escaped_material";
  case PhotonStatus::count:
    break;
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
  std::vector<double> geometric_path_m;
  std::vector<double> time_ns;
  std::vector<double> weight;
  std::vector<PhotonStatus> status;
  // Terminal optical surface, or kNoSurfaceId if none was reached.
  std::vector<std::uint32_t> surface_id;

  static constexpr std::uint32_t kNoSurfaceId = std::numeric_limits<std::uint32_t>::max();

  explicit PhotonResultBlock(std::size_t size = 0)
      : position_m(size), direction(size), optical_path_m(size), geometric_path_m(size),
        time_ns(size), weight(size), status(size), surface_id(size, kNoSurfaceId) {}
};

enum class OpticalInteractionKind { mirror, detector, obscurer, refractive_interface };

template <std::size_t PointCapacity> struct BasicPathRecord {
  std::uint64_t photon_id{};
  double wavelength_nm{400.0};
  PhotonStatus status{PhotonStatus::missed_primary};
  // entrance, M1, M2 (when present), focal plane.  Single-reflector paths
  // simply use the first three entries.
  std::array<Vec3, PointCapacity> points_m{};
  std::uint8_t point_count{};
  double path_length_m{};
  Vec3 final_direction{};
  double incidence_primary_deg{};
  double incidence_secondary_deg{};
  double incidence_focal_deg{};
  // Optical response applied before the terminal geometry loss. At a detector
  // this is the arriving fraction; otherwise it is lost at the terminal.
  double surviving_throughput{1.0};
  std::array<OpticalInteractionKind, PointCapacity - 1> interaction_kinds{};
  std::array<Vec3, PointCapacity - 1> interaction_normals{};
  std::array<Vec3, PointCapacity - 1> interaction_incoming_directions{};
  std::array<Vec3, PointCapacity - 1> interaction_outgoing_directions{};
  std::array<double, PointCapacity - 1> interaction_throughput = [] {
    std::array<double, PointCapacity - 1> values{};
    values.fill(1.0);
    return values;
  }();
  bool material_transport{};
  double optical_path_m{};
  double group_delay_ns{};
  std::array<double, PointCapacity - 1> interaction_optical_path_m{};
  std::array<double, PointCapacity - 1> interaction_group_delay_ns{};
  std::uint32_t terminal_surface_id{PhotonResultBlock::kNoSurfaceId};
  std::array<std::uint32_t, PointCapacity - 1> interaction_surface_ids = [] {
    std::array<std::uint32_t, PointCapacity - 1> values{};
    values.fill(PhotonResultBlock::kNoSurfaceId);
    return values;
  }();
};

using PathRecord = BasicPathRecord<4>;

} // namespace obdeect
