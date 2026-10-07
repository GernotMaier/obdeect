#pragma once

#include <cstdint>

namespace obdeect {

// Counter-based dimensions: changing order, batch size or total count does not
// change a ray with the same identity and seed. No mutable random engine exists.
[[nodiscard]] inline double source_uniform(std::uint64_t id, std::uint64_t seed,
                                           std::uint64_t dimension) {
  std::uint64_t value = id ^ (seed + 0x9e3779b97f4a7c15ULL * (dimension + 1));
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  value ^= value >> 31;
  return static_cast<double>(value >> 11) * 0x1.0p-53;
}

} // namespace obdeect
