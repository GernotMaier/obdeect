// The limit comparison of the replay tool (cpp/bench/replay_compare.hpp).
#include "check.hpp"

#include "replay_compare.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace {

using replay::Limits;
using replay::Record;

// Typical magnitudes of the nine fields, so that a limit is tested at a realistic scale.
constexpr std::array<double, replay::kFieldCount> kBase = {0.9,  -0.4, 1.2,   0.3, 0.5,
                                                           -0.8, 40.0, 100.0, 0.7};

std::vector<Record> sample(std::size_t n) {
  std::vector<Record> records(n);
  for (std::size_t i = 0; i < n; ++i) {
    records[i].status = static_cast<std::uint8_t>(i % 3);
    records[i].terminal = static_cast<std::uint32_t>(i % 5);
    records[i].primary = 1;
    records[i].secondary = 2;
    records[i].value = kBase;
    records[i].value[0] += 1e-3 * static_cast<double>(i % 7);
  }
  return records;
}

// Copy of `records` with field k of photon i set to `value`.
std::vector<Record> with_value(std::vector<Record> records, std::size_t i, std::size_t k,
                               double value) {
  records[i].value[k] = value;
  return records;
}

bool passes(const std::vector<Record> &a, const std::vector<Record> &b, const Limits &limits) {
  return replay::compare_limits(a, b, limits).within;
}

void identical_passes_every_tier() {
  const auto records = sample(1000);
  for (const char name : {'G', 'S', 'T'}) {
    const auto limits = replay::tier_limits(name);
    CHECK(limits.has_value());
    if (limits)
      CHECK(passes(records, records, *limits));
  }
  CHECK(!replay::tier_limits('X').has_value());
}

void nan_positions_of_missed_photons() {
  auto a = sample(10), b = sample(10);
  for (auto *set : {&a, &b})
    for (std::size_t k = 0; k < 3; ++k)
      (*set)[4].value[k] = std::numeric_limits<double>::quiet_NaN();
  // Bit-identical NaNs are equal; a NaN against a number fails every tier.
  CHECK(passes(a, b, replay::kTierT));
  b[4].value[0] = 0.5;
  CHECK(!passes(a, b, replay::kTierG));
  CHECK(!passes(a, b, replay::kTierT));
}

// A change of 0.5 limits passes, a change of 1.5 limits fails, in every field of tiers G and S.
// Weight limits are relative, so the change is 0.5 or 1.5 limits times the weight.
void each_field_just_inside_and_outside_its_limit() {
  for (const Limits *tier : {&replay::kTierG, &replay::kTierS}) {
    for (std::size_t k = 0; k < replay::kFieldCount; ++k) {
      const double scale = k == 8 ? kBase[8] : 1.0;
      const double limit = tier->field(k) * scale;
      const auto a = sample(20);
      CHECK(passes(a, with_value(a, 7, k, kBase[k] + 0.5 * limit), *tier));
      CHECK(!passes(a, with_value(a, 7, k, kBase[k] + 1.5 * limit), *tier));
    }
  }
}

// Tier T: geometry must be bit-identical, only the weight may change (by 1e-7 relative).
void tier_t_geometry_must_be_bit_identical() {
  const auto a = sample(20);
  for (std::size_t k = 0; k < 8; ++k) {
    const double up = std::nextafter(kBase[k], std::numeric_limits<double>::infinity());
    CHECK(!passes(a, with_value(a, 7, k, up), replay::kTierT));
  }
  CHECK(passes(a, with_value(a, 7, 8, kBase[8] * (1.0 + 0.5e-7)), replay::kTierT));
  CHECK(!passes(a, with_value(a, 7, 8, kBase[8] * (1.0 + 1.5e-7)), replay::kTierT));
}

void flips_per_million_limit() {
  // 10^6 photons: 5 flips is 5 per 10^6 (within tier G), 6 flips is not.
  const std::size_t n = 1'000'000;
  const auto a = sample(n);
  auto b = a;
  for (std::size_t i = 0; i < 5; ++i)
    b[i].status = static_cast<std::uint8_t>(a[i].status + 1);
  const auto five = replay::compare_limits(a, b, replay::kTierG);
  CHECK(five.flips == 5);
  CHECK(five.flips_within);
  CHECK(five.within);
  b[5].status = static_cast<std::uint8_t>(a[5].status + 1);
  const auto six = replay::compare_limits(a, b, replay::kTierG);
  CHECK(six.flips == 6);
  CHECK(!six.flips_within);
  CHECK(!six.within);
  // Tier T allows no flips at all.
  CHECK(!passes(a, b, replay::kTierT));
}

void flipped_photons_are_not_compared_by_value() {
  const auto a = sample(10);
  auto b = a;
  b[3].status = static_cast<std::uint8_t>(a[3].status + 1);
  b[3].value[0] += 1.0; // Would fail the position limit if the photon were compared.
  const auto c = replay::compare_limits(a, b, replay::kTierG);
  CHECK(c.flips == 1);
  CHECK(c.fields[0].max == 0.0);
  CHECK(c.fields[0].within);
  CHECK(c.within == c.flips_within);
}

void surface_id_change_is_a_flip() {
  const auto a = sample(10);
  auto b = a;
  b[2].secondary = 9;
  const auto c = replay::compare_limits(a, b, replay::kTierS);
  CHECK(c.flips == 1);
  CHECK(c.id_differences[3] == 1);
  CHECK(c.id_differences[2] == 0);
}

void quantiles_are_nearest_rank() {
  const auto a = sample(100);
  auto b = a;
  // Time differences of 1, 2, ..., 100 times 1e-12 ns (time is about 100 ns; ulp 1.4e-14).
  for (std::size_t i = 0; i < 100; ++i)
    b[i].value[7] += static_cast<double>(i + 1) * 1e-12;
  const auto c = replay::compare_limits(a, b, replay::kTierG);
  CHECK_NEAR(c.fields[7].median, 50e-12, 2e-14);
  CHECK_NEAR(c.fields[7].p99, 99e-12, 2e-14);
  CHECK_NEAR(c.fields[7].p999, 100e-12, 2e-14);
  CHECK_NEAR(c.fields[7].max, 100e-12, 2e-14);
}

} // namespace

void all_tests() {
  identical_passes_every_tier();
  nan_positions_of_missed_photons();
  each_field_just_inside_and_outside_its_limit();
  tier_t_geometry_must_be_bit_identical();
  flips_per_million_limit();
  flipped_photons_are_not_compared_by_value();
  surface_id_change_is_a_flip();
  quantiles_are_nearest_rank();
}

int main() { return run_checks(all_tests); }
