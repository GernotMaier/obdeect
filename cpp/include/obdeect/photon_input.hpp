#pragma once

#include "obdeect/sources.hpp"
#include <algorithm>
#include <span>
#include <stdexcept>

namespace obdeect {

// One batch belongs to one event, array reuse and telescope, in telescope-local
// SI coordinates. Adapters resolve direction/frame and spectral conventions
// before exposing OpticalPhoton. Zero/negative wavelength is not silently blue.
struct PhotonBatchContext {
  std::uint64_t run_id{};
  std::uint64_t event_id{};
  std::uint64_t array_id{};
  std::uint64_t telescope_id{};
};

struct PhotonReadResult {
  PhotonBatchContext context;
  std::size_t count{};
  bool eof{};
};

// Virtual dispatch occurs once per batch, never per ray. Caller owns/reuses the
// destination storage; errors throw (never masquerade as EOF). A non-EOF read
// must return count > 0. Event boundaries must not be mixed within a batch.
class PhotonReader {
 public:
  virtual ~PhotonReader() = default;
  virtual PhotonReadResult read(std::span<OpticalPhoton> destination) = 0;
};

// First adapter: deterministic test/calibration sources. Source storage must
// outlive this reader; no photon copy beyond the requested batch is allocated.
class MemoryPhotonReader final : public PhotonReader {
 public:
  MemoryPhotonReader(std::span<const OpticalPhoton> source, PhotonBatchContext context)
      : source_(source), context_(context) {}

  PhotonReadResult read(std::span<OpticalPhoton> destination) override {
    if (destination.empty()) throw std::invalid_argument("photon batch capacity must be positive");
    const auto count = std::min(destination.size(), source_.size() - offset_);
    std::copy_n(source_.begin() + offset_, count, destination.begin());
    offset_ += count;
    return {context_, count, offset_ == source_.size()};
  }

 private:
  std::span<const OpticalPhoton> source_;
  PhotonBatchContext context_;
  std::size_t offset_{};
};

}  // namespace obdeect
