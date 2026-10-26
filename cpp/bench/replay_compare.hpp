// Limit-based comparison of two replay record sets (tier G, S or T).
//
// Used by replay.cpp (`compare --tier` and `compare --limits`) and by
// cpp/tests/test_replay_compare.cpp. Standard library only, no files.
//
// A photon is a flip when its status or any surface id differs. Flips are
// counted per 10^6 photons and their values are not compared. For every other
// photon, each field's |difference| is collected; the field passes when its
// maximum is within the limit (weight: relative to the first file's weight).
// Quantiles are nearest-rank over the compared photons.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <vector>

namespace replay {

inline constexpr std::size_t kFieldCount = 9;
// Field order of Record::value.
inline constexpr std::array<const char *, kFieldCount> kFieldNames = {
    "position.x",  "position.y",    "position.z", "direction.x", "direction.y",
    "direction.z", "path_length_m", "time_ns",    "weight"};

// One photon of a replay file: status, the three surface ids and the nine values.
struct Record {
  std::uint8_t status{};
  std::uint32_t terminal{}, primary{}, secondary{};
  std::array<double, kFieldCount> value{};
};

// Limits in the units of the table in OBDEECT-PERFORMANCE-PLAN-3-CLASS-B.md.
// A limit of zero requires bit-identical values.
struct Limits {
  double flips_per_million;
  double position_m;    // Each position component, absolute.
  double direction;     // Each direction component, absolute.
  double path_length_m; // Absolute.
  double time_ns;       // Absolute.
  double weight_relative;

  [[nodiscard]] double field(std::size_t k) const {
    if (k < 3)
      return position_m;
    if (k < 6)
      return direction;
    if (k == 6)
      return path_length_m;
    return k == 7 ? time_ns : weight_relative;
  }
};

inline constexpr Limits kTierG{5.0, 1e-8, 1e-9, 1e-8, 1e-7, 1e-9};
inline constexpr Limits kTierS{5.0, 1e-10, 1e-13, 1e-10, 1e-9, 1e-13};
inline constexpr Limits kTierT{0.0, 0.0, 0.0, 0.0, 0.0, 1e-7};

// Tier by its letter ('G', 'S' or 'T').
[[nodiscard]] inline std::optional<Limits> tier_limits(char name) {
  switch (name) {
  case 'G':
    return kTierG;
  case 'S':
    return kTierS;
  case 'T':
    return kTierT;
  default:
    return std::nullopt;
  }
}

struct FieldStats {
  std::size_t bitwise_different{};
  double median{}, p99{}, p999{}, max{};
  double limit{};
  bool within{true};
};

struct Comparison {
  std::size_t photons{};
  std::size_t flips{};
  // Differences of status, terminal, primary and secondary surface ids.
  std::array<std::size_t, 4> id_differences{};
  double flips_per_million{};
  bool flips_within{true};
  std::array<FieldStats, kFieldCount> fields{};
  bool within{true};
};

[[nodiscard]] inline bool same_bits(double a, double b) {
  return std::memcmp(&a, &b, sizeof a) == 0;
}

// Nearest-rank quantile of sorted values (q in (0, 1]); NaN for no values.
[[nodiscard]] inline double nearest_rank(const std::vector<double> &sorted, double q) {
  if (sorted.empty())
    return std::numeric_limits<double>::quiet_NaN();
  const auto rank = static_cast<std::size_t>(std::ceil(q * static_cast<double>(sorted.size())));
  return sorted[std::clamp<std::size_t>(rank, 1, sorted.size()) - 1];
}

// Compares two equally long record sets. Both sets are assumed to be the same
// photons in the same order.
[[nodiscard]] inline Comparison compare_limits(const std::vector<Record> &a,
                                               const std::vector<Record> &b, const Limits &limits) {
  Comparison out;
  out.photons = std::min(a.size(), b.size());
  std::array<std::vector<double>, kFieldCount> diffs;
  for (auto &d : diffs)
    d.reserve(out.photons);
  for (std::size_t i = 0; i < out.photons; ++i) {
    const Record &x = a[i], &y = b[i];
    const bool discrete = x.status != y.status || x.terminal != y.terminal ||
                          x.primary != y.primary || x.secondary != y.secondary;
    out.id_differences[0] += x.status != y.status;
    out.id_differences[1] += x.terminal != y.terminal;
    out.id_differences[2] += x.primary != y.primary;
    out.id_differences[3] += x.secondary != y.secondary;
    if (discrete) {
      ++out.flips;
      continue;
    }
    for (std::size_t k = 0; k < kFieldCount; ++k) {
      const double u = x.value[k], v = y.value[k];
      double d = 0.0; // Bit-identical values (including NaN of missed photons) differ by 0.
      if (!same_bits(u, v)) {
        ++out.fields[k].bitwise_different;
        d = std::abs(u - v);
        if (k == 8 && d != 0.0)
          d = u != 0.0 ? d / std::abs(u) : std::numeric_limits<double>::infinity();
        // A NaN difference fails any limit: it is mapped to infinity.
        if (std::isnan(d))
          d = std::numeric_limits<double>::infinity();
      }
      diffs[k].push_back(d);
    }
  }
  out.flips_per_million =
      out.photons ? 1e6 * static_cast<double>(out.flips) / static_cast<double>(out.photons) : 0.0;
  out.flips_within = out.flips_per_million <= limits.flips_per_million;
  out.within = out.flips_within;
  for (std::size_t k = 0; k < kFieldCount; ++k) {
    auto &d = diffs[k];
    std::sort(d.begin(), d.end());
    FieldStats &s = out.fields[k];
    s.limit = limits.field(k);
    s.median = nearest_rank(d, 0.5);
    s.p99 = nearest_rank(d, 0.99);
    s.p999 = nearest_rank(d, 0.999);
    s.max = d.empty() ? 0.0 : d.back();
    // A zero limit needs bit-identical values, not only zero differences.
    s.within = s.max <= s.limit && (s.limit > 0.0 || s.bitwise_different == 0);
    out.within = out.within && s.within;
  }
  return out;
}

} // namespace replay
