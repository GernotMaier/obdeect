// Bit-level regression check of the tracer.
//
//   replay write <model.json> <out.bin>
//   replay compare [--class-b] <a.bin> <b.bin>
//
// `write` traces six fixed sets of 10^6 photons (star at 10 km and at
// infinity, 0, 1 and 3 degrees off axis; responses on; seed 7; photon ids
// increasing across the sets) and writes one record per photon, field by
// field. The sampling disk is read from the model's "sampling" entry.
// `compare` exits 0 only if both files are identical, or, with --class-b, if
// they agree within the class-B limits (status and surface ids may differ for
// at most 5 photons per 10^6; position 1e-9 m; direction 1e-12; time 1e-9 ns;
// weight 1e-12 relative). Only API of the core on `main` is used, so the same
// source builds against the reference headers. Not run by ctest.
#include "obdeect/json.hpp"
#include "obdeect/model_json.hpp"
#include "obdeect/source.hpp"
#include "obdeect/trace.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace obdeect;

namespace {

constexpr std::size_t kPhotonsPerSet = 1'000'000;
constexpr std::uint64_t kSeed = 7;
constexpr double kDistancesM[] = {1e4, std::numeric_limits<double>::infinity()};
constexpr double kOffAxisDeg[] = {0.0, 1.0, 3.0};

// Record layout: status u8, three u32 surface ids, 6 + 3 doubles.
constexpr std::size_t kRecordBytes = 1 + 3 * 4 + 9 * 8;

struct Record {
  std::uint8_t status;
  std::uint32_t terminal, primary, secondary;
  std::array<double, 9> value; // position(3), direction(3), path, time, weight.
};

void put(std::vector<char> &buffer, const void *data, std::size_t size) {
  const auto *bytes = static_cast<const char *>(data);
  buffer.insert(buffer.end(), bytes, bytes + size);
}

std::string read_text(const char *path) {
  std::ifstream file(path);
  if (!file) {
    std::fprintf(stderr, "cannot read %s\n", path);
    std::exit(2);
  }
  std::stringstream text;
  text << file.rdbuf();
  return text.str();
}

int write_command(const char *model_path, const char *out_path) {
  const std::string text = read_text(model_path);
  const OpticalModel model = load_model(text);
  const json::Value root = json::Parser::parse(text);
  const auto number = [&](const char *key) {
    return std::get<double>(root.find("sampling")->find(key)->data);
  };
  const double radius = number("disk_radius_m"), disk_z = number("disk_z_m"),
               start_z = number("start_z_m");
  std::ofstream out(out_path, std::ios::binary);
  std::uint64_t next_id = 0;
  std::vector<char> buffer;
  for (const double distance : kDistancesM)
    for (const double degrees : kOffAxisDeg) {
      const double theta = degrees * 3.14159265358979323846 / 180.0;
      const PointSource source{
          {std::sin(theta), 0.0, std::cos(theta)}, distance, radius, disk_z, start_z, 400.0};
      buffer.clear();
      for (std::size_t k = 0; k < kPhotonsPerSet; ++k, ++next_id) {
        const auto photon = sample_source(source, next_id, kSeed);
        const auto r =
            trace_photon(model, photon.ray, 400.0, photon.time_ns, 1.0, next_id, kSeed, true);
        const auto status = static_cast<std::uint8_t>(r.status);
        put(buffer, &status, 1);
        put(buffer, &r.terminal_surface, 4);
        put(buffer, &r.primary_surface, 4);
        put(buffer, &r.secondary_surface, 4);
        const double values[9] = {r.position.x,    r.position.y,  r.position.z,
                                  r.direction.x,   r.direction.y, r.direction.z,
                                  r.path_length_m, r.time_ns,     r.weight};
        put(buffer, values, sizeof values);
      }
      out.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    }
  return out ? 0 : 2;
}

std::vector<char> read_binary(const char *path) {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    std::fprintf(stderr, "cannot read %s\n", path);
    std::exit(2);
  }
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

Record record_at(const std::vector<char> &data, std::size_t i) {
  const char *p = data.data() + i * kRecordBytes;
  Record r;
  r.status = static_cast<std::uint8_t>(p[0]);
  std::memcpy(&r.terminal, p + 1, 4);
  std::memcpy(&r.primary, p + 5, 4);
  std::memcpy(&r.secondary, p + 9, 4);
  std::memcpy(r.value.data(), p + 13, sizeof r.value);
  return r;
}

bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

int compare_command(bool class_b, const char *path_a, const char *path_b) {
  const auto a = read_binary(path_a), b = read_binary(path_b);
  if (a.size() != b.size() || a.size() % kRecordBytes != 0) {
    std::printf("DIFFERENT SIZE: %zu and %zu bytes\n", a.size(), b.size());
    return 1;
  }
  const std::size_t n = a.size() / kRecordBytes;
  static const char *names[9] = {"position.x",    "position.y",  "position.z",
                                 "direction.x",   "direction.y", "direction.z",
                                 "path_length_m", "time_ns",     "weight"};
  std::size_t status_diff = 0, terminal_diff = 0, primary_diff = 0, secondary_diff = 0;
  std::size_t id_status_diff = 0, printed = 0;
  std::array<std::size_t, 9> bit_diff{};
  std::array<double, 9> max_diff{};
  bool within = true;
  for (std::size_t i = 0; i < n; ++i) {
    const Record x = record_at(a, i), y = record_at(b, i);
    const bool discrete = x.status != y.status || x.terminal != y.terminal ||
                          x.primary != y.primary || x.secondary != y.secondary;
    status_diff += x.status != y.status;
    terminal_diff += x.terminal != y.terminal;
    primary_diff += x.primary != y.primary;
    secondary_diff += x.secondary != y.secondary;
    if (discrete) {
      ++id_status_diff;
      if (printed++ < 20)
        std::printf("photon %zu: status %u vs %u, terminal %u vs %u, primary %u vs %u, secondary "
                    "%u vs %u\n",
                    i, x.status, y.status, x.terminal, y.terminal, x.primary, y.primary,
                    x.secondary, y.secondary);
      continue; // Values of a photon that took another path are not comparable.
    }
    for (std::size_t k = 0; k < 9; ++k) {
      const double u = x.value[k], v = y.value[k];
      if (same_bits(u, v))
        continue;
      ++bit_diff[k];
      const double d = std::abs(u - v);
      if (std::isnan(d)) {
        within = false;
        max_diff[k] = std::numeric_limits<double>::infinity();
        continue;
      }
      max_diff[k] = std::max(max_diff[k], d);
      const double limit = k < 3 ? 1e-9 : k < 6 ? 1e-12 : k == 6 ? 1e-9 : k == 7 ? 1e-9 : 1e-12;
      const double scale = k == 8 ? std::abs(u) : 1.0;
      if (d > limit * scale)
        within = false;
    }
  }
  std::printf("photons %zu\n", n);
  std::printf("status differences %zu, terminal_surface %zu, primary_surface %zu, "
              "secondary_surface %zu (photons with any: %zu)\n",
              status_diff, terminal_diff, primary_diff, secondary_diff, id_status_diff);
  std::size_t any_value = 0;
  for (std::size_t k = 0; k < 9; ++k) {
    std::printf("%-14s bitwise differences %9zu  max |diff| %.3e\n", names[k], bit_diff[k],
                max_diff[k]);
    any_value += bit_diff[k];
  }
  const std::size_t allowed = 5 * (n / 1'000'000);
  if (class_b) {
    const bool ok = within && id_status_diff <= allowed;
    std::printf("CLASS B: %s (limit %zu discrete differences)\n", ok ? "within limits" : "FAILED",
                allowed);
    return ok ? 0 : 1;
  }
  const bool identical = id_status_diff == 0 && any_value == 0;
  std::printf("%s\n", identical ? "IDENTICAL" : "DIFFERENT");
  return identical ? 0 : 1;
}

} // namespace

int main(int argc, char **argv) {
  if (argc == 4 && std::strcmp(argv[1], "write") == 0)
    return write_command(argv[2], argv[3]);
  if (argc == 4 && std::strcmp(argv[1], "compare") == 0)
    return compare_command(false, argv[2], argv[3]);
  if (argc == 5 && std::strcmp(argv[1], "compare") == 0 && std::strcmp(argv[2], "--class-b") == 0)
    return compare_command(true, argv[3], argv[4]);
  std::fprintf(stderr, "usage: %s write <model.json> <out.bin> | compare [--class-b] <a> <b>\n",
               argv[0]);
  return 2;
}
