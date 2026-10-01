#pragma once

#include "obdeect/axisymmetric_optics.hpp"
#include "obdeect/json.hpp"
#include "obdeect/photon_buffer.hpp"
#include "obdeect/segmented_optical_model.hpp"

#include <cmath>
#include <algorithm>
#include <array>
#include <charconv>
#include <cstdint>
#include <fstream>
#include <limits>
#include <numbers>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace obdeect {

struct AxisymmetricOpticalModel {
  ModelProvenance provenance;
  AxisymmetricMirror primary;
  AxisymmetricMirror secondary;
  AxisymmetricMirror detector;
  std::optional<SpectralResponse> primary_reflectivity;
  std::optional<SpectralResponse> secondary_reflectivity;
};

namespace detail {

inline void append_json_string(std::string& output, std::string_view value) {
  static constexpr char hex[] = "0123456789abcdef";
  output.push_back('"');
  for (const unsigned char character : value) {
    switch (character) {
      case '"': output += "\\\""; break;
      case '\\': output += "\\\\"; break;
      case '\b': output += "\\b"; break;
      case '\f': output += "\\f"; break;
      case '\n': output += "\\n"; break;
      case '\r': output += "\\r"; break;
      case '\t': output += "\\t"; break;
      default:
        if (character < 0x20) {
          output += "\\u00";
          output.push_back(hex[character >> 4]);
          output.push_back(hex[character & 0x0f]);
        } else {
          output.push_back(static_cast<char>(character));
        }
        break;
    }
  }
  output.push_back('"');
}

inline void append_canonical_json(const json::Value& value, std::string& output,
                                  bool omit_root_hash = false, bool root = true) {
  switch (value.kind) {
    case json::Value::Kind::null: output += "null"; break;
    case json::Value::Kind::boolean: output += value.boolean ? "true" : "false"; break;
    case json::Value::Kind::number: {
      if (!value.number_text.empty()) {
        output += value.number_text;
      } else {
        std::array<char, 64> buffer{};
        const auto result = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value.number);
        if (result.ec != std::errc{}) return;
        output.append(buffer.data(), result.ptr);
      }
      break;
    }
    case json::Value::Kind::string: append_json_string(output, value.string); break;
    case json::Value::Kind::array:
      output.push_back('[');
      for (std::size_t index = 0; index < value.array.size(); ++index) {
        if (index != 0) output.push_back(',');
        append_canonical_json(value.array[index], output, false, false);
      }
      output.push_back(']');
      break;
    case json::Value::Kind::object: {
      std::vector<const std::pair<std::string, json::Value>*> entries;
      entries.reserve(value.object.size());
      for (const auto& entry : value.object) {
        if (root && omit_root_hash && entry.first == "optical_model_sha256") continue;
        entries.push_back(&entry);
      }
      std::sort(entries.begin(), entries.end(), [](const auto* left, const auto* right) {
        return left->first < right->first;
      });
      output.push_back('{');
      for (std::size_t index = 0; index < entries.size(); ++index) {
        if (index != 0) output.push_back(',');
        append_json_string(output, entries[index]->first);
        output.push_back(':');
        append_canonical_json(entries[index]->second, output, false, false);
      }
      output.push_back('}');
      break;
    }
  }
}

[[nodiscard]] inline std::string canonical_json_without_hash(const json::Value& root) {
  std::string result;
  append_canonical_json(root, result, true);
  return result;
}

class Sha256 {
 public:
  Sha256() : state_{0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
                    0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U} {}

  void update(std::string_view input) {
    for (const unsigned char byte : input) {
      buffer_[buffer_size_++] = byte;
      if (buffer_size_ == buffer_.size()) {
        transform();
        bit_count_ += 512;
        buffer_size_ = 0;
      }
    }
  }

  [[nodiscard]] std::string finish() {
    const std::uint64_t original_bits = bit_count_ + buffer_size_ * 8;
    buffer_[buffer_size_++] = 0x80;
    if (buffer_size_ > 56) {
      while (buffer_size_ < 64) buffer_[buffer_size_++] = 0;
      transform();
      buffer_size_ = 0;
    }
    while (buffer_size_ < 56) buffer_[buffer_size_++] = 0;
    for (int shift = 56; shift >= 0; shift -= 8) buffer_[buffer_size_++] = static_cast<unsigned char>(original_bits >> shift);
    transform();

    static constexpr char hex[] = "0123456789abcdef";
    std::string result;
    result.reserve(64);
    for (const std::uint32_t word : state_) {
      for (int shift = 28; shift >= 0; shift -= 4) result.push_back(hex[(word >> shift) & 0x0f]);
    }
    return result;
  }

 private:
  static constexpr std::array<std::uint32_t, 64> constants_ = {
      0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U,
      0x923f82a4U, 0xab1c5ed5U, 0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
      0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U, 0xe49b69c1U, 0xefbe4786U,
      0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
      0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U,
      0x06ca6351U, 0x14292967U, 0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
      0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U, 0xa2bfe8a1U, 0xa81a664bU,
      0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
      0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU,
      0x5b9cca4fU, 0x682e6ff3U, 0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
      0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};

  static constexpr std::uint32_t rotate_right(std::uint32_t value, int shift) {
    return (value >> shift) | (value << (32 - shift));
  }

  void transform() {
    std::array<std::uint32_t, 64> schedule{};
    for (std::size_t index = 0; index < 16; ++index) {
      schedule[index] = (static_cast<std::uint32_t>(buffer_[4 * index]) << 24) |
                        (static_cast<std::uint32_t>(buffer_[4 * index + 1]) << 16) |
                        (static_cast<std::uint32_t>(buffer_[4 * index + 2]) << 8) |
                        static_cast<std::uint32_t>(buffer_[4 * index + 3]);
    }
    for (std::size_t index = 16; index < schedule.size(); ++index) {
      const auto first = rotate_right(schedule[index - 15], 7) ^ rotate_right(schedule[index - 15], 18) ^
                         (schedule[index - 15] >> 3);
      const auto second = rotate_right(schedule[index - 2], 17) ^ rotate_right(schedule[index - 2], 19) ^
                          (schedule[index - 2] >> 10);
      schedule[index] = schedule[index - 16] + first + schedule[index - 7] + second;
    }
    auto working = state_;
    for (std::size_t index = 0; index < schedule.size(); ++index) {
      const auto sigma_one = rotate_right(working[4], 6) ^ rotate_right(working[4], 11) ^ rotate_right(working[4], 25);
      const auto choice = (working[4] & working[5]) ^ (~working[4] & working[6]);
      const auto temp_one = working[7] + sigma_one + choice + constants_[index] + schedule[index];
      const auto sigma_zero = rotate_right(working[0], 2) ^ rotate_right(working[0], 13) ^ rotate_right(working[0], 22);
      const auto majority = (working[0] & working[1]) ^ (working[0] & working[2]) ^ (working[1] & working[2]);
      const auto temp_two = sigma_zero + majority;
      working[7] = working[6];
      working[6] = working[5];
      working[5] = working[4];
      working[4] = working[3] + temp_one;
      working[3] = working[2];
      working[2] = working[1];
      working[1] = working[0];
      working[0] = temp_one + temp_two;
    }
    for (std::size_t index = 0; index < state_.size(); ++index) state_[index] += working[index];
  }

  std::array<std::uint32_t, 8> state_{};
  std::array<unsigned char, 64> buffer_{};
  std::size_t buffer_size_{};
  std::uint64_t bit_count_{};
};

[[nodiscard]] inline std::string sha256(std::string_view input) {
  Sha256 digest;
  digest.update(input);
  return digest.finish();
}

[[nodiscard]] inline const json::Value* field(const json::Value& value, std::string_view name) {
  return value.find(name);
}

[[nodiscard]] inline const std::string* string_field(const json::Value& value, std::string_view name) {
  const auto* result = field(value, name);
  return result && result->kind == json::Value::Kind::string ? &result->string : nullptr;
}

[[nodiscard]] inline std::optional<double> number_field(const json::Value& value, std::string_view name) {
  const auto* result = field(value, name);
  if (!result || result->kind != json::Value::Kind::number || !std::isfinite(result->number)) return std::nullopt;
  return result->number;
}

[[nodiscard]] inline std::optional<std::uint32_t> uint_field(const json::Value& value, std::string_view name) {
  const auto number = number_field(value, name);
  if (!number || *number < 0.0 || *number > std::numeric_limits<std::uint32_t>::max() ||
      std::floor(*number) != *number) return std::nullopt;
  return static_cast<std::uint32_t>(*number);
}

[[nodiscard]] inline std::optional<Vec3> vec3_field(const json::Value& value, std::string_view name) {
  const auto* result = field(value, name);
  if (!result || result->kind != json::Value::Kind::array || result->array.size() != 3) return std::nullopt;
  Vec3 vector{};
  double* components[] = {&vector.x, &vector.y, &vector.z};
  for (std::size_t index = 0; index < 3; ++index) {
    if (result->array[index].kind != json::Value::Kind::number || !std::isfinite(result->array[index].number)) return std::nullopt;
    *components[index] = result->array[index].number;
  }
  return vector;
}

[[nodiscard]] inline std::optional<FacetShape> facet_shape(const json::Value& value) {
  const auto* shape = string_field(value, "shape");
  if (!shape) return std::nullopt;
  if (*shape == "circle") return FacetShape::circle;
  if (*shape == "hexagon_flat_y" || *shape == "hexagon") return FacetShape::hexagon_flat_y;
  if (*shape == "square") return FacetShape::square;
  if (*shape == "hexagon_flat_x") return FacetShape::hexagon_flat_x;
  return std::nullopt;
}

[[nodiscard]] inline std::optional<ModelProvenance> provenance(const json::Value& root) {
  const auto* source = field(root, "provenance");
  if (!source || source->kind != json::Value::Kind::object) return std::nullopt;
  const auto* model = string_field(*source, "model");
  const auto* version = string_field(*source, "model_version");
  const auto* hash = string_field(root, "optical_model_sha256");
  if (!model || !version || !hash) return std::nullopt;
  ModelProvenance result{*model, *version, *hash};
  return has_valid_provenance(result) ? std::optional{std::move(result)} : std::nullopt;
}

[[nodiscard]] inline std::optional<SpectralResponse> response_field(const json::Value& value, std::string_view name) {
  const auto* response = field(value, name);
  if (!response || response->kind != json::Value::Kind::array) return std::nullopt;
  SpectralResponse result;
  for (const auto& entry : response->array) {
    if (entry.kind != json::Value::Kind::object) return std::nullopt;
    const auto wavelength = number_field(entry, "wavelength_nm");
    const auto response_value = number_field(entry, "response");
    if (!wavelength || !response_value) return std::nullopt;
    result.wavelength_nm.push_back(*wavelength);
    result.response.push_back(*response_value);
  }
  return result.is_valid() ? std::optional{std::move(result)} : std::nullopt;
}

[[nodiscard]] inline std::optional<AxisymmetricMirror> axisymmetric_surface(const json::Value& value) {
  const auto vertex = number_field(value, "vertex_z_m");
  const auto inner = number_field(value, "inner_radius_m");
  const auto outer = number_field(value, "outer_radius_m");
  const auto scale = number_field(value, "radial_scale_m");
  const auto* coefficients = field(value, "coefficient_m");
  if (!vertex || !inner || !outer || !scale || !coefficients || coefficients->kind != json::Value::Kind::array ||
      coefficients->array.size() != 13) return std::nullopt;
  AxisymmetricMirror result{};
  result.vertex_z_m = *vertex;
  result.inner_radius_m = *inner;
  result.outer_radius_m = *outer;
  result.surface.radial_scale_m = *scale;
  for (std::size_t index = 0; index < result.surface.coefficient_m.size(); ++index) {
    const auto& coefficient = coefficients->array[index];
    if (coefficient.kind != json::Value::Kind::number || !std::isfinite(coefficient.number)) return std::nullopt;
    result.surface.coefficient_m[index] = coefficient.number;
  }
  return is_valid(result) ? std::optional{result} : std::nullopt;
}

[[nodiscard]] inline std::optional<json::Value> read_json(const std::string& path) {
  std::ifstream input(path);
  if (!input) return std::nullopt;
  std::stringstream buffer;
  buffer << input.rdbuf();
  const auto root = json::Parser{buffer.str()}.parse();
  if (!root || root->kind != json::Value::Kind::object || !string_field(*root, "format") ||
      *string_field(*root, "format") != "obdeect.compiled-optical-model.v1") return std::nullopt;
  const auto* hash = string_field(*root, "optical_model_sha256");
  if (!hash || hash->size() != 64 ||
      std::any_of(hash->begin(), hash->end(), [](const char character) {
        return !((character >= '0' && character <= '9') ||
                 (character >= 'a' && character <= 'f'));
      }) ||
      sha256(canonical_json_without_hash(*root)) != *hash) {
    return std::nullopt;
  }
  return root;
}

}  // namespace detail

// Read the single versioned JSON artifact emitted by obdeect-compile-optical-model.
// Its trace_model section is directly consumable; the remainder is retained for
// provenance, validation and inspection.
inline std::optional<CompiledSegmentedOpticalModel> read_segmented_optical_model(const std::string& path) {
  const auto root = detail::read_json(path);
  if (!root) return std::nullopt;
  const auto model_provenance = detail::provenance(*root);
  const auto* trace = detail::field(*root, "trace_model");
  if (!model_provenance || !trace || trace->kind != json::Value::Kind::object ||
      !detail::string_field(*trace, "kind") || *detail::string_field(*trace, "kind") != "segmented") return std::nullopt;
  const auto* facet_values = detail::field(*trace, "primary_facets");
  const auto* detector_values = detail::field(*trace, "detector_surfaces");
  const auto* obscurer_values = detail::field(*trace, "cylinder_obscurers");
  if (!facet_values || !detector_values || !obscurer_values || facet_values->kind != json::Value::Kind::array ||
      detector_values->kind != json::Value::Kind::array || obscurer_values->kind != json::Value::Kind::array) return std::nullopt;
  std::vector<ImportedFacet> facets;
  for (const auto& value : facet_values->array) {
    if (value.kind != json::Value::Kind::object) return std::nullopt;
    const auto id = detail::uint_field(value, "id");
    const auto centre = detail::vec3_field(value, "centre_m");
    const auto normal = detail::vec3_field(value, "normal");
    const auto tangent = detail::vec3_field(value, "tangent");
    const auto diameter = detail::number_field(value, "diameter_m");
    const auto focal_length = detail::number_field(value, "focal_length_m");
    const auto shape = detail::facet_shape(value);
    if (!id || !centre || !normal || !tangent || !diameter || !focal_length || !shape) return std::nullopt;
    facets.push_back({*id, *centre, *normal, *diameter, *focal_length, *shape, *tangent, 2.0 * *focal_length});
  }
  std::vector<ImportedDetectorSurface> detectors;
  for (const auto& value : detector_values->array) {
    if (value.kind != json::Value::Kind::object) return std::nullopt;
    const auto id = detail::uint_field(value, "id");
    const auto centre = detail::vec3_field(value, "centre_m");
    const auto normal = detail::vec3_field(value, "normal");
    const auto tangent = detail::vec3_field(value, "tangent");
    const auto diameter = detail::number_field(value, "diameter_m");
    const auto shape = detail::facet_shape(value);
    if (!id || !centre || !normal || !tangent || !diameter || !shape) return std::nullopt;
    detectors.push_back({*id, *centre, *normal, *diameter, *shape, *tangent});
  }
  std::vector<ImportedCylinderObscurer> obscurers;
  for (const auto& value : obscurer_values->array) {
    if (value.kind != json::Value::Kind::object) return std::nullopt;
    const auto id = detail::uint_field(value, "id");
    const auto first = detail::vec3_field(value, "first_endpoint_m");
    const auto second = detail::vec3_field(value, "second_endpoint_m");
    const auto diameter = detail::number_field(value, "diameter_m");
    if (!id || !first || !second || !diameter) return std::nullopt;
    obscurers.push_back({*id, *first, *second, *diameter});
  }
  const auto* reflectivity_value = detail::field(*trace, "primary_reflectivity");
  const auto reflectivity = detail::response_field(*trace, "primary_reflectivity");
  if ((reflectivity_value && reflectivity_value->kind != json::Value::Kind::array) ||
      (reflectivity_value && !reflectivity_value->array.empty() && !reflectivity) ||
      facets.empty() || detectors.empty()) return std::nullopt;
  return compile_segmented_optical_model({*model_provenance, std::move(facets), std::move(detectors),
                                          std::move(obscurers), std::move(reflectivity)});
}

[[nodiscard]] inline std::optional<AxisymmetricOpticalModel> read_axisymmetric_optical_model(const std::string& path) {
  const auto root = detail::read_json(path);
  if (!root) return std::nullopt;
  const auto model_provenance = detail::provenance(*root);
  const auto* trace = detail::field(*root, "trace_model");
  if (!model_provenance || !trace || trace->kind != json::Value::Kind::object ||
      !detail::string_field(*trace, "kind") || *detail::string_field(*trace, "kind") != "axisymmetric") return std::nullopt;
  const auto* primary = detail::field(*trace, "primary");
  const auto* secondary = detail::field(*trace, "secondary");
  const auto* detector = detail::field(*trace, "detector");
  if (!primary || !secondary || !detector) return std::nullopt;
  const auto primary_surface = detail::axisymmetric_surface(*primary);
  const auto secondary_surface = detail::axisymmetric_surface(*secondary);
  const auto detector_surface = detail::axisymmetric_surface(*detector);
  const auto* primary_response = detail::field(*trace, "primary_reflectivity");
  const auto* secondary_response = detail::field(*trace, "secondary_reflectivity");
  const auto primary_reflectivity = detail::response_field(*trace, "primary_reflectivity");
  const auto secondary_reflectivity = detail::response_field(*trace, "secondary_reflectivity");
  if ((primary_response && primary_response->kind != json::Value::Kind::array) ||
      (secondary_response && secondary_response->kind != json::Value::Kind::array) ||
      (primary_response && !primary_response->array.empty() && !primary_reflectivity) ||
      (secondary_response && !secondary_response->array.empty() && !secondary_reflectivity) ||
      !primary_surface || !secondary_surface || !detector_surface)
    return std::nullopt;
  return AxisymmetricOpticalModel{*model_provenance, *primary_surface, *secondary_surface, *detector_surface,
                                  std::move(primary_reflectivity), std::move(secondary_reflectivity)};
}

[[nodiscard]] inline PathRecord trace_axisymmetric_optical_model(const Ray& input, std::uint64_t photon_id,
                                                          const AxisymmetricOpticalModel& optical_model) {
  PathRecord record{};
  record.photon_id = photon_id;
  record.points_m[0] = input.position_m;
  record.point_count = 1;
  const auto direction = normalised_checked(input.direction);
  if (!direction) { record.status = PhotonStatus::invalid_input; return record; }
  Ray ray{input.position_m, *direction};
  const auto primary = intersect_axisymmetric_mirror(ray, optical_model.primary);
  if (!primary) { record.status = PhotonStatus::missed_primary; record.final_direction = ray.direction; return record; }
  record.points_m[1] = primary->point_m;
  record.point_count = 2;
  record.path_length_m = primary->distance_m;
  record.incidence_primary_deg = std::acos(std::clamp(std::abs(dot(ray.direction, primary->unit_normal)), 0.0, 1.0)) * 180.0 / std::numbers::pi;
  const auto after_primary = reflect(ray, *primary);
  if (!after_primary) { record.status = PhotonStatus::invalid_input; return record; }
  ray = *after_primary;
  const auto secondary = intersect_axisymmetric_mirror(ray, optical_model.secondary);
  if (!secondary) { record.status = PhotonStatus::missed_screen; record.final_direction = ray.direction; return record; }
  record.points_m[2] = secondary->point_m;
  record.point_count = 3;
  record.path_length_m += secondary->distance_m;
  record.incidence_secondary_deg = std::acos(std::clamp(std::abs(dot(ray.direction, secondary->unit_normal)), 0.0, 1.0)) * 180.0 / std::numbers::pi;
  const auto after_secondary = reflect(ray, *secondary);
  if (!after_secondary) { record.status = PhotonStatus::invalid_input; return record; }
  ray = *after_secondary;
  const auto detector = intersect_axisymmetric_mirror(ray, optical_model.detector);
  if (!detector) { record.status = PhotonStatus::missed_screen; record.final_direction = ray.direction; return record; }
  record.points_m[3] = detector->point_m;
  record.point_count = 4;
  record.path_length_m += detector->distance_m;
  record.incidence_focal_deg = std::acos(std::clamp(std::abs(dot(ray.direction, detector->unit_normal)), 0.0, 1.0)) * 180.0 / std::numbers::pi;
  record.final_direction = ray.direction;
  record.status = PhotonStatus::detected;
  return record;
}

}  // namespace obdeect
