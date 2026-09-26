#pragma once

#include "obdeect/photon_input.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace obdeect {

struct EventioInputLimits {
  std::size_t max_block_bytes{512ULL * 1024 * 1024};
  std::size_t max_bunches_per_telescope{8ULL * 1024 * 1024};
};

struct EventioRunInfo {
  double wavelength_lower_nm{};
  double wavelength_upper_nm{};
  bool ceffic{};
  bool atmext{};
  bool refraction{};
  bool has_event_header{};
  bool reuse_weight_known{};
  double array_time_offset_ns{};
  double observation_altitude_m{};
  bool curved{};
};

// One output row represents one CORSIKA photon bunch. Its weight is the
// fractional photon count, and wavelength zero remains unresolved. This
// adapter deliberately rejects CEFFIC/photoelectron bunches; they cannot be
// passed through the ordinary atmospheric and optical efficiency stages.
// The external EventIO C library is needed only by this optional target.
class EventioPhotonReader final : public PhotonReader {
 public:
  explicit EventioPhotonReader(const std::string& path, EventioInputLimits limits = {});
  ~EventioPhotonReader() override;
  EventioPhotonReader(const EventioPhotonReader&) = delete;
  EventioPhotonReader& operator=(const EventioPhotonReader&) = delete;
  EventioPhotonReader(EventioPhotonReader&&) noexcept;
  EventioPhotonReader& operator=(EventioPhotonReader&&) noexcept;

  PhotonReadResult read(std::span<OpticalPhoton> destination) override;
  // Metadata for the most recently returned batch (stable across lookahead).
  [[nodiscard]] const EventioRunInfo& run_info() const noexcept;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace obdeect
