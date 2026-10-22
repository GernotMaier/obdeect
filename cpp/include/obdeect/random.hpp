#pragma once

#include <array>
#include <cstdint>

namespace obdeect {

// Philox-4x32-10 (Salmon et al., SC'11): a counter-based generator. A draw
// depends only on (seed, photon id, stream, block), never on the order in
// which photons are processed, so results are identical for any block size
// and thread count.
[[nodiscard]] constexpr std::array<std::uint32_t, 4> philox4x32(std::array<std::uint32_t, 4> c,
                                                                std::array<std::uint32_t, 2> k) {
  constexpr std::uint32_t m0 = 0xD2511F53U, m1 = 0xCD9E8D57U;
  constexpr std::uint32_t w0 = 0x9E3779B9U, w1 = 0xBB67AE85U;
  for (int round = 0; round < 10; ++round) {
    const std::uint64_t p0 = std::uint64_t{m0} * c[0];
    const std::uint64_t p1 = std::uint64_t{m1} * c[2];
    c = {static_cast<std::uint32_t>(p1 >> 32) ^ c[1] ^ k[0], static_cast<std::uint32_t>(p1),
         static_cast<std::uint32_t>(p0 >> 32) ^ c[3] ^ k[1], static_cast<std::uint32_t>(p0)};
    k = {k[0] + w0, k[1] + w1};
  }
  return c;
}

// Independent random streams of one photon.
enum class Stream : std::uint32_t { source = 0, primary_scatter = 1, secondary_scatter = 2 };

// Four uniform deviates in the open interval (0, 1).
[[nodiscard]] constexpr std::array<double, 4> uniform4(std::uint64_t seed, std::uint64_t photon_id,
                                                       Stream stream, std::uint32_t block = 0) {
  const auto r = philox4x32(
      {static_cast<std::uint32_t>(photon_id), static_cast<std::uint32_t>(photon_id >> 32),
       static_cast<std::uint32_t>(stream), block},
      {static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)});
  constexpr double scale = 1.0 / 4294967296.0;
  return {(r[0] + 0.5) * scale, (r[1] + 0.5) * scale, (r[2] + 0.5) * scale, (r[3] + 0.5) * scale};
}

} // namespace obdeect
