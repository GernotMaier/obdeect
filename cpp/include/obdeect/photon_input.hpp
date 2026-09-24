#pragma once

#include "obdeect/sources.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace obdeect {

// One batch belongs to one event, array reuse and telescope, in telescope-local
// SI coordinates. Adapters resolve direction/frame and spectral conventions
// before exposing OpticalPhoton. Zero/negative wavelength is not silently blue.
struct PhotonBatchContext {
  std::uint64_t run_id{};
  std::uint64_t event_id{};
  std::uint64_t array_id{};
  std::uint64_t telescope_id{};
  Vec3 telescope_position_m{};
  // Multiplicative reuse factor supplied by the upstream array simulation.
  // It is metadata, not a photon count, and must never be rounded.
  double array_reuse_weight{1.0};
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

// A deliberately small, format-neutral interchange reader.  It makes CSV
// fixtures and calibration exports usable without putting an EventIO decoder
// in the core.  Required columns are:
// run_id,event_id,array_id,telescope_id,photon_id,x_m,y_m,z_m,dx,dy,dz,
// wavelength_nm,time_ns,weight.  Optional columns preserve raw-bunch and
// telescope metadata: bunch_id,emission_height_m,emission_distance_m,
// telescope_x_m,telescope_y_m,telescope_z_m,array_reuse_weight.
//
// One returned batch never crosses context boundaries.  Zero wavelength is
// preserved as the explicit unspecified-spectrum sentinel.  Malformed input
// throws invalid_argument; it is never reported as EOF.
class CsvPhotonReader final : public PhotonReader {
 public:
  explicit CsvPhotonReader(const std::string& path) : input_(path) {
    if (!input_) throw std::invalid_argument("cannot open photon CSV: " + path);
    initialise_header();
  }

  PhotonReadResult read(std::span<OpticalPhoton> destination) override {
    if (destination.empty()) throw std::invalid_argument("photon batch capacity must be positive");
    const auto first = next_row();
    if (!first) return {{}, 0, true};

    const PhotonBatchContext context = first->context;
    destination[0] = first->photon;
    std::size_t count = 1;
    while (count < destination.size()) {
      const auto row = next_row();
      if (!row) return {context, count, true};
      if (!same_context(context, row->context)) {
        pending_ = *row;
        return {context, count, false};
      }
      destination[count++] = row->photon;
    }

    // Read one row ahead so eof is accurate even when the final batch exactly
    // fills caller storage.  The row is retained for the next call.
    const auto row = next_row();
    if (!row) return {context, count, true};
    pending_ = *row;
    return {context, count, false};
  }

 private:
  struct Row {
    PhotonBatchContext context;
    OpticalPhoton photon;
  };

  static constexpr std::string_view kRequiredColumns[] = {
      "run_id", "event_id", "array_id", "telescope_id", "photon_id", "x_m", "y_m", "z_m",
      "dx", "dy", "dz", "wavelength_nm", "time_ns", "weight"};

  void initialise_header() {
    std::string header;
    if (!std::getline(input_, header)) throw std::invalid_argument("photon CSV is missing a header");
    const auto columns = split(header);
    for (std::size_t index = 0; index < columns.size(); ++index) {
      if (columns[index].empty() || !header_index_.emplace(columns[index], index).second)
        throw std::invalid_argument("photon CSV has an empty or duplicate column name");
    }
    for (const auto name : kRequiredColumns)
      if (!header_index_.contains(std::string{name}))
        throw std::invalid_argument("photon CSV is missing required column: " + std::string{name});
  }

  [[nodiscard]] std::optional<Row> next_row() {
    if (pending_) return std::exchange(pending_, std::nullopt);
    std::string line;
    while (std::getline(input_, line)) {
      ++line_number_;
      if (line.empty() || line.front() == '#') continue;
      return parse_row(split(line));
    }
    if (input_.bad()) throw std::invalid_argument("error while reading photon CSV");
    return std::nullopt;
  }

  [[nodiscard]] Row parse_row(const std::vector<std::string>& fields) const {
    if (fields.size() != header_index_.size()) throw invalid("wrong number of fields");
    const auto integer = [this, &fields](std::string_view name) {
      const auto& text = fields.at(header_index_.at(std::string{name}));
      std::size_t consumed{};
      try {
        // std::stoull accepts a leading minus sign and converts it modulo the
        // unsigned range. IDs must fail closed instead of silently wrapping.
        if (text.empty() || text.front() == '-') throw invalid("invalid integer in " + std::string{name});
        const auto value = std::stoull(text, &consumed);
        if (consumed != text.size()) throw invalid("invalid integer in " + std::string{name});
        return static_cast<std::uint64_t>(value);
      } catch (const std::exception&) { throw invalid("invalid integer in " + std::string{name}); }
    };
    const auto number = [this, &fields](std::string_view name, double fallback) {
      const auto found = header_index_.find(std::string{name});
      if (found == header_index_.end()) return fallback;
      const auto& text = fields.at(found->second);
      std::size_t consumed{};
      try {
        const double value = std::stod(text, &consumed);
        if (consumed != text.size() || !std::isfinite(value)) throw invalid("invalid number in " + std::string{name});
        return value;
      } catch (const std::exception&) { throw invalid("invalid number in " + std::string{name}); }
    };

    PhotonBatchContext context{integer("run_id"), integer("event_id"), integer("array_id"), integer("telescope_id"),
                               {number("telescope_x_m", 0.0), number("telescope_y_m", 0.0),
                                number("telescope_z_m", 0.0)},
                               number("array_reuse_weight", 1.0)};
    OpticalPhoton photon{{{number("x_m", 0.0), number("y_m", 0.0), number("z_m", 0.0)},
                          {number("dx", 0.0), number("dy", 0.0), number("dz", 0.0)}},
                         integer("photon_id"), number("wavelength_nm", 0.0), number("time_ns", 0.0),
                         number("weight", 0.0),
                         header_index_.contains("bunch_id") ? integer("bunch_id") : 0,
                         number("emission_height_m", std::numeric_limits<double>::quiet_NaN()),
                         number("emission_distance_m", std::numeric_limits<double>::quiet_NaN())};
    if (!normalised_checked(photon.ray.direction) || photon.wavelength_nm < 0.0 || photon.weight < 0.0 ||
        context.array_reuse_weight < 0.0)
      throw invalid("invalid photon or batch metadata");
    return {context, photon};
  }

  [[nodiscard]] std::invalid_argument invalid(const std::string& message) const {
    return std::invalid_argument("photon CSV line " + std::to_string(line_number_) + ": " + message);
  }

  static std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> fields;
    const std::size_t record_end = !line.empty() && line.back() == '\r' ? line.size() - 1 : line.size();
    std::size_t start{};
    for (std::size_t comma = line.find(',', start);
         comma != std::string::npos && comma < record_end; comma = line.find(',', start)) {
      fields.push_back(line.substr(start, comma - start));
      start = comma + 1;
    }
    fields.push_back(line.substr(start, record_end - start));
    return fields;
  }

  static bool same_context(const PhotonBatchContext& left, const PhotonBatchContext& right) {
    return left.run_id == right.run_id && left.event_id == right.event_id && left.array_id == right.array_id &&
           left.telescope_id == right.telescope_id && left.telescope_position_m.x == right.telescope_position_m.x &&
           left.telescope_position_m.y == right.telescope_position_m.y && left.telescope_position_m.z == right.telescope_position_m.z &&
           left.array_reuse_weight == right.array_reuse_weight;
  }

  std::ifstream input_;
  std::unordered_map<std::string, std::size_t> header_index_;
  std::optional<Row> pending_;
  std::size_t line_number_{1};
};

}  // namespace obdeect
