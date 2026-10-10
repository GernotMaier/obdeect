#pragma once

#include "obdeect/pixel_response.hpp"

#include "obdeect/axisymmetric_optics.hpp"
#include "obdeect/axisymmetric_segments.hpp"
#include "obdeect/detector_assignment.hpp"
#include "obdeect/detector_planes.hpp"
#include "obdeect/json.hpp"
#include "obdeect/optical_model.hpp"
#include "obdeect/photon_buffer.hpp"
#include "obdeect/segmented_optical_model.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace obdeect {

struct AxisymmetricOpticalModel {
  ModelProvenance provenance;
  AxisymmetricMirror primary;
  AxisymmetricMirror secondary;
  AxisymmetricMirror detector;
  std::optional<SpectralResponse> primary_reflectivity;
  std::optional<SpectralResponse> secondary_reflectivity;
  std::vector<AxisymmetricSegment> primary_segments{};
  std::vector<AxisymmetricSegment> secondary_segments{};
  std::uint32_t primary_surface_id{PhotonResultBlock::kNoSurfaceId};
  std::uint32_t secondary_surface_id{PhotonResultBlock::kNoSurfaceId};
  std::uint32_t detector_surface_id{PhotonResultBlock::kNoSurfaceId};
  bool block_incoming_secondary{};
  std::optional<TelescopeTransmission> telescope_transmission{};
  double propagation_group_index{1.0};
  std::optional<CompiledDetectorPlanes> detector_planes{};
  std::optional<CompiledDetectorAssignmentGrid> detector_assignment{};
  bool camera_degradation_in_detector_frame{};
  std::optional<MirrorScatter> primary_scatter{}, secondary_scatter{};
  std::optional<CameraResponse> camera_response{};
  std::optional<CompiledDetectorPlanes> primary_to_secondary_planes{};
  std::vector<ImportedCylinderObscurer> primary_to_secondary_cylinders{};
  std::optional<CompiledDetectorPlanes> incoming_obscurer_planes{};
  std::vector<OpaqueSurface> opaque_obscurers{};
  std::shared_ptr<const PixelResponses> pixel_responses{};
  std::optional<SpatialResponse> primary_degradation{}, secondary_degradation{},
      camera_degradation{};
};

namespace detail {

template <class Output> inline void append_json_string(Output &output, std::string_view value) {
  static constexpr char hex[] = "0123456789abcdef";
  output.push_back('"');
  for (const unsigned char character : value) {
    switch (character) {
    case '"':
      output += "\\\"";
      break;
    case '\\':
      output += "\\\\";
      break;
    case '\b':
      output += "\\b";
      break;
    case '\f':
      output += "\\f";
      break;
    case '\n':
      output += "\\n";
      break;
    case '\r':
      output += "\\r";
      break;
    case '\t':
      output += "\\t";
      break;
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

template <class Output>
inline void append_canonical_json(const json::Value &value, Output &output,
                                  bool omit_root_hash = false, bool root = true) {
  switch (value.kind()) {
  case json::Value::Kind::null:
    output += "null";
    break;
  case json::Value::Kind::boolean:
    output += value.boolean() ? "true" : "false";
    break;
  case json::Value::Kind::number: {
    if (!value.number_text().empty()) {
      output += value.number_text();
    } else {
      std::array<char, 64> buffer{};
      const auto result =
          std::to_chars(buffer.data(), buffer.data() + buffer.size(), value.number());
      if (result.ec != std::errc{})
        return;
      output.append(buffer.data(), result.ptr);
    }
    break;
  }
  case json::Value::Kind::string:
    append_json_string(output, value.string());
    break;
  case json::Value::Kind::array:
    output.push_back('[');
    for (std::size_t index = 0; index < value.array().size(); ++index) {
      if (index != 0)
        output.push_back(',');
      append_canonical_json(value.array()[index], output, false, false);
    }
    output.push_back(']');
    break;
  case json::Value::Kind::object: {
    std::vector<const std::pair<std::string, json::Value> *> entries;
    entries.reserve(value.object().size());
    for (const auto &entry : value.object()) {
      if (root && omit_root_hash && entry.first == "optical_model_sha256")
        continue;
      entries.push_back(&entry);
    }
    std::sort(entries.begin(), entries.end(),
              [](const auto *left, const auto *right) { return left->first < right->first; });
    output.push_back('{');
    for (std::size_t index = 0; index < entries.size(); ++index) {
      if (index != 0)
        output.push_back(',');
      append_json_string(output, entries[index]->first);
      output.push_back(':');
      append_canonical_json(entries[index]->second, output, false, false);
    }
    output.push_back('}');
    break;
  }
  }
}

[[nodiscard]] inline std::string canonical_json_without_hash(const json::Value &root) {
  std::string result;
  append_canonical_json(root, result, true);
  return result;
}

class Sha256 {
public:
  Sha256()
      : state_{0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
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
      while (buffer_size_ < 64)
        buffer_[buffer_size_++] = 0;
      transform();
      buffer_size_ = 0;
    }
    while (buffer_size_ < 56)
      buffer_[buffer_size_++] = 0;
    for (int shift = 56; shift >= 0; shift -= 8)
      buffer_[buffer_size_++] = static_cast<unsigned char>(original_bits >> shift);
    transform();

    static constexpr char hex[] = "0123456789abcdef";
    std::string result;
    result.reserve(64);
    for (const std::uint32_t word : state_) {
      for (int shift = 28; shift >= 0; shift -= 4)
        result.push_back(hex[(word >> shift) & 0x0f]);
    }
    return result;
  }

private:
  static constexpr std::array<std::uint32_t, 64> constants_ = {
      0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U, 0x923f82a4U,
      0xab1c5ed5U, 0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U, 0x72be5d74U, 0x80deb1feU,
      0x9bdc06a7U, 0xc19bf174U, 0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU,
      0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU, 0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
      0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U, 0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU,
      0x53380d13U, 0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U, 0xa2bfe8a1U, 0xa81a664bU,
      0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U, 0x19a4c116U,
      0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
      0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U, 0x90befffaU, 0xa4506cebU, 0xbef9a3f7U,
      0xc67178f2U};

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
      const auto first = rotate_right(schedule[index - 15], 7) ^
                         rotate_right(schedule[index - 15], 18) ^ (schedule[index - 15] >> 3);
      const auto second = rotate_right(schedule[index - 2], 17) ^
                          rotate_right(schedule[index - 2], 19) ^ (schedule[index - 2] >> 10);
      schedule[index] = schedule[index - 16] + first + schedule[index - 7] + second;
    }
    auto working = state_;
    for (std::size_t index = 0; index < schedule.size(); ++index) {
      const auto sigma_one =
          rotate_right(working[4], 6) ^ rotate_right(working[4], 11) ^ rotate_right(working[4], 25);
      const auto choice = (working[4] & working[5]) ^ (~working[4] & working[6]);
      const auto temp_one = working[7] + sigma_one + choice + constants_[index] + schedule[index];
      const auto sigma_zero =
          rotate_right(working[0], 2) ^ rotate_right(working[0], 13) ^ rotate_right(working[0], 22);
      const auto majority =
          (working[0] & working[1]) ^ (working[0] & working[2]) ^ (working[1] & working[2]);
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
    for (std::size_t index = 0; index < state_.size(); ++index)
      state_[index] += working[index];
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

// Feed canonical tokens directly into SHA-256 instead of allocating a second
// serialization of the complete optical-model artifact.
class CanonicalHashOutput {
public:
  void operator+=(std::string_view token) { digest_.update(token); }
  void push_back(char character) { digest_.update(std::string_view{&character, 1}); }
  void append(const char *first, const char *last) {
    digest_.update(std::string_view{first, static_cast<std::size_t>(last - first)});
  }
  std::string finish() { return digest_.finish(); }

private:
  Sha256 digest_;
};

[[nodiscard]] inline std::string canonical_json_hash(const json::Value &root) {
  CanonicalHashOutput output;
  append_canonical_json(root, output, true);
  return output.finish();
}

[[nodiscard]] inline const json::Value *field(const json::Value &value, std::string_view name) {
  return value.find(name);
}

[[nodiscard]] inline bool fields_supported(const json::Value &value,
                                           std::initializer_list<std::string_view> allowed) {
  return value.kind() == json::Value::Kind::object &&
         std::all_of(value.object().begin(), value.object().end(), [&](const auto &field) {
           return std::find(allowed.begin(), allowed.end(), field.first) != allowed.end();
         });
}

[[nodiscard]] inline const std::string *string_field(const json::Value &value,
                                                     std::string_view name) {
  const auto *result = field(value, name);
  return result && result->kind() == json::Value::Kind::string ? &result->string() : nullptr;
}

[[nodiscard]] inline std::optional<double> number_field(const json::Value &value,
                                                        std::string_view name) {
  const auto *result = field(value, name);
  if (!result || result->kind() != json::Value::Kind::number || !std::isfinite(result->number()))
    return std::nullopt;
  return result->number();
}

[[nodiscard]] inline std::optional<std::uint32_t> uint_field(const json::Value &value,
                                                             std::string_view name) {
  const auto number = number_field(value, name);
  if (!number || *number < 0.0 || *number > std::numeric_limits<std::uint32_t>::max() ||
      std::floor(*number) != *number)
    return std::nullopt;
  return static_cast<std::uint32_t>(*number);
}

[[nodiscard]] inline std::optional<Vec3> vec3_field(const json::Value &value,
                                                    std::string_view name) {
  const auto *result = field(value, name);
  if (!result || result->kind() != json::Value::Kind::array || result->array().size() != 3)
    return std::nullopt;
  Vec3 vector{};
  double *components[] = {&vector.x, &vector.y, &vector.z};
  for (std::size_t index = 0; index < 3; ++index) {
    if (result->array()[index].kind() != json::Value::Kind::number ||
        !std::isfinite(result->array()[index].number()))
      return std::nullopt;
    *components[index] = result->array()[index].number();
  }
  return vector;
}

[[nodiscard]] inline std::optional<FacetShape> facet_shape(const json::Value &value) {
  const auto *shape = string_field(value, "shape");
  if (!shape)
    return std::nullopt;
  if (*shape == "circle")
    return FacetShape::circle;
  if (*shape == "hexagon_flat_y" || *shape == "hexagon")
    return FacetShape::hexagon_flat_y;
  if (*shape == "square")
    return FacetShape::square;
  if (*shape == "hexagon_flat_x")
    return FacetShape::hexagon_flat_x;
  return std::nullopt;
}

[[nodiscard]] inline std::optional<ModelProvenance> provenance(const json::Value &root) {
  const auto *source = field(root, "provenance");
  if (!source || source->kind() != json::Value::Kind::object)
    return std::nullopt;
  const auto *model = string_field(*source, "model");
  const auto *version = string_field(*source, "model_version");
  const auto *hash = string_field(root, "optical_model_sha256");
  if (!model || !version || !hash)
    return std::nullopt;
  ModelProvenance result{*model, *version, *hash};
  return has_valid_provenance(result) ? std::optional{std::move(result)} : std::nullopt;
}

[[nodiscard]] inline std::optional<TableInterpolation>
interpolation_field(const json::Value &value) {
  if (!fields_supported(value,
                        {"boundary", "scheme", "x_log", "y_log", "value_log", "coefficients"}))
    return std::nullopt;
  TableInterpolation result;
  const auto *boundary = string_field(value, "boundary");
  const auto scheme = uint_field(value, "scheme");
  if (!boundary || !scheme || *scheme > 4)
    return std::nullopt;
  if (*boundary == "clamp")
    result.boundary = TableBoundary::clamp;
  else if (*boundary == "zero")
    result.boundary = TableBoundary::zero;
  else if (*boundary != "reject")
    return std::nullopt;
  result.scheme = *scheme;
  for (const auto &[name, destination] :
       {std::pair{"x_log", &result.x_log}, std::pair{"y_log", &result.y_log},
        std::pair{"value_log", &result.value_log}}) {
    const auto *flag = value.find(name);
    if (!flag || flag->kind() != json::Value::Kind::boolean)
      return std::nullopt;
    *destination = flag->boolean();
  }
  if (const auto *coefficients = value.find("coefficients")) {
    if (coefficients->kind() != json::Value::Kind::array)
      return std::nullopt;
    for (const auto &row : coefficients->array()) {
      if (row.kind() != json::Value::Kind::array || row.array().size() != 4)
        return std::nullopt;
      std::array<double, 4> values;
      for (std::size_t i = 0; i < 4; ++i) {
        if (row.array()[i].kind() != json::Value::Kind::number ||
            !std::isfinite(row.array()[i].number()))
          return std::nullopt;
        values[i] = row.array()[i].number();
      }
      result.coefficients.push_back(values);
    }
  }
  return result;
}

[[nodiscard]] inline std::optional<SpectralResponse> response_field(const json::Value &value,
                                                                    std::string_view name) {
  const auto *response = field(value, name);
  if (!response || response->kind() != json::Value::Kind::array)
    return std::nullopt;
  SpectralResponse result;
  const bool angular =
      !response->array().empty() && response->array().front().find("incidence_angle_deg");
  std::map<double, std::map<double, double>> grid;
  for (const auto &entry : response->array()) {
    if (!fields_supported(entry, {"wavelength_nm", "response", "incidence_angle_deg",
                                  "interpolation", "wavelength_sampling"}))
      return std::nullopt;
    if (const auto *sampling = entry.find("wavelength_sampling")) {
      if (&entry != &response->array().front() ||
          !fields_supported(*sampling, {"width_nm", "offset_nm", "first_bin", "last_bin",
                                        "projection_angle_deg"}))
        return std::nullopt;
      const auto width = number_field(*sampling, "width_nm"),
                 offset = number_field(*sampling, "offset_nm");
      const auto first = uint_field(*sampling, "first_bin"),
                 last = uint_field(*sampling, "last_bin");
      if (!width || !offset || !first || !last)
        return std::nullopt;
      result.wavelength_sampling = WavelengthSampling{*width, *offset, *first, *last};
      if (sampling->find("projection_angle_deg")) {
        const auto angle = number_field(*sampling, "projection_angle_deg");
        if (!angle)
          return std::nullopt;
        result.wavelength_sampling->projection_angle_deg = *angle;
      }
    }
    if (const auto *options = entry.find("interpolation")) {
      if (&entry != &response->array().front())
        return std::nullopt;
      const auto parsed = interpolation_field(*options);
      if (!parsed)
        return std::nullopt;
      result.interpolation = *parsed;
    }
    const auto wavelength = number_field(entry, "wavelength_nm");
    const auto value = number_field(entry, "response");
    const auto angle = number_field(entry, "incidence_angle_deg");
    if (!wavelength || !value || angular != (entry.find("incidence_angle_deg") != nullptr) ||
        (angular && !angle))
      return std::nullopt;
    if (angular) {
      if (!grid[*wavelength].emplace(*angle, *value).second)
        return std::nullopt;
    } else {
      result.wavelength_nm.push_back(*wavelength);
      result.response.push_back(*value);
    }
  }
  for (const auto &[wavelength, values] : grid) {
    if (result.incidence_angle_deg.empty())
      for (const auto &[angle, value] : values)
        result.incidence_angle_deg.push_back(angle);
    if (values.size() != result.incidence_angle_deg.size())
      return std::nullopt;
    result.wavelength_nm.push_back(wavelength);
    std::size_t index = 0;
    for (const auto &[angle, value] : values) {
      if (angle != result.incidence_angle_deg[index++])
        return std::nullopt;
      result.response.push_back(value);
    }
    result.spectral_envelope.push_back(
        std::max_element(values.begin(), values.end(), [](const auto &a, const auto &b) {
          return a.second < b.second;
        })->second);
  }
  result.envelope_interpolation = result.interpolation;
  result.envelope_interpolation.y_log = false;
  return result.is_valid() ? std::optional{std::move(result)} : std::nullopt;
}

[[nodiscard]] inline std::optional<std::vector<OpaqueSurface>>
opaque_fields(const json::Value &trace) {
  std::vector<OpaqueSurface> result;
  const auto *values = trace.find("opaque_obscurers");
  if (!values)
    return result;
  if (values->kind() != json::Value::Kind::array)
    return std::nullopt;
  for (const auto &value : values->array()) {
    const auto id = uint_field(value, "id");
    const auto *shape = string_field(value, "shape");
    if (!id || !shape)
      return std::nullopt;
    OpaqueSurface surface;
    surface.id = *id;
    if (*shape == "quadrilateral") {
      if (!fields_supported(value, {"id", "shape", "vertices_m"}))
        return std::nullopt;
      const auto *vertices = value.find("vertices_m");
      if (!vertices || vertices->kind() != json::Value::Kind::array ||
          vertices->array().size() != 4)
        return std::nullopt;
      for (std::size_t i = 0; i < 4; ++i) {
        const auto &vertex = vertices->array()[i];
        if (vertex.kind() != json::Value::Kind::array || vertex.array().size() != 3)
          return std::nullopt;
        for (const auto &v : vertex.array())
          if (v.kind() != json::Value::Kind::number)
            return std::nullopt;
        surface.vertices[i] = {vertex.array()[0].number(), vertex.array()[1].number(),
                               vertex.array()[2].number()};
      }
    } else if (*shape == "hollow_frustum" || *shape == "solid_frustum") {
      if (!fields_supported(value, {"id", "shape", "first_endpoint_m", "second_endpoint_m",
                                    "first_radius_m", "second_radius_m", "thickness_m"}))
        return std::nullopt;
      const auto first = vec3_field(value, "first_endpoint_m"),
                 second = vec3_field(value, "second_endpoint_m");
      const auto r1 = number_field(value, "first_radius_m"),
                 r2 = number_field(value, "second_radius_m"),
                 thickness = number_field(value, "thickness_m");
      if (!first || !second || !r1 || !r2 || !thickness)
        return std::nullopt;
      surface.shape = *shape == "hollow_frustum" ? OpaqueSurface::Shape::hollow_frustum
                                                 : OpaqueSurface::Shape::solid_frustum;
      surface.first = *first;
      surface.second = *second;
      surface.first_radius_m = *r1;
      surface.second_radius_m = *r2;
      surface.thickness_m = *thickness;
    } else
      return std::nullopt;
    if (!is_valid(surface))
      return std::nullopt;
    result.push_back(surface);
  }
  return result;
}

[[nodiscard]] inline std::optional<SpatialResponse> spatial_field(const json::Value &parent,
                                                                  std::string_view name) {
  const auto *value = parent.find(name);
  if (!value || !fields_supported(*value, {"x_m", "y_m", "response", "clip", "x_basis", "y_basis",
                                           "interpolation"}))
    return std::nullopt;
  SpatialResponse result;
  for (const auto &[name, destination] :
       {std::pair{"x_m", &result.x_m}, std::pair{"y_m", &result.y_m},
        std::pair{"response", &result.response}}) {
    const auto *values = value->find(name);
    if (!values || values->kind() != json::Value::Kind::array)
      return std::nullopt;
    for (const auto &v : values->array()) {
      if (v.kind() != json::Value::Kind::number)
        return std::nullopt;
      destination->push_back(v.number());
    }
  }
  for (const auto &[name, destination] :
       {std::pair{"x_basis", &result.x_basis}, std::pair{"y_basis", &result.y_basis}})
    if (value->find(name)) {
      const auto vector = vec3_field(*value, name);
      if (!vector)
        return std::nullopt;
      *destination = *vector;
    }
  if (const auto *options = value->find("interpolation")) {
    const auto parsed = interpolation_field(*options);
    if (!parsed)
      return std::nullopt;
    result.interpolation = *parsed;
  }
  if (const auto *clip = value->find("clip")) {
    if (clip->kind() != json::Value::Kind::boolean)
      return std::nullopt;
    result.clip = clip->boolean();
  }
  return result.is_valid() ? std::optional{std::move(result)} : std::nullopt;
}

[[nodiscard]] inline std::optional<MirrorScatter> scatter_field(const json::Value &parent,
                                                                std::string_view name) {
  const auto *value = parent.find(name);
  if (!value || !fields_supported(*value, {"sigma1_rad", "sigma2_rad", "fraction2", "method"}))
    return std::nullopt;
  const auto first = number_field(*value, "sigma1_rad"),
             second = number_field(*value, "sigma2_rad"),
             fraction = number_field(*value, "fraction2");
  const auto *method = string_field(*value, "method");
  if (!first || !second || !fraction || !method ||
      (*method != "outgoing_angles" && *method != "surface_slopes"))
    return std::nullopt;
  MirrorScatter result{*first, *fraction, *second,
                       *method == "outgoing_angles" ? MirrorScatterMethod::outgoing_angles
                                                    : MirrorScatterMethod::surface_slopes};
  return result.is_valid() ? std::optional{result} : std::nullopt;
}

[[nodiscard]] inline std::optional<TelescopeTransmission>
telescope_transmission_field(const json::Value &parent) {
  const auto *value = parent.find("telescope_transmission");
  if (!value || !fields_supported(
                    *value, {"on_axis", "amplitude", "angular_scale_rad", "power", "outer_power"}))
    return {};
  const auto t0 = number_field(*value, "on_axis"), a = number_field(*value, "amplitude"),
             scale = number_field(*value, "angular_scale_rad"), p = number_field(*value, "power"),
             q = number_field(*value, "outer_power");
  if (!t0 || !a || !scale || !p || !q)
    return {};
  TelescopeTransmission result{*t0, *a, *scale, *p, *q};
  return result.is_valid() ? std::optional{result} : std::nullopt;
}

[[nodiscard]] inline std::optional<CameraResponse> camera_response_field(const json::Value &parent,
                                                                         std::string_view name) {
  const auto *value = parent.find(name);
  if (!value || !fields_supported(*value, {"semantics", "camera_transmission", "camera_filter",
                                           "lightguide_efficiency"}))
    return std::nullopt;
  const auto *semantics = string_field(*value, "semantics");
  const auto transmission = number_field(*value, "camera_transmission");
  if (!semantics || *semantics != "measured_complete_response" || !transmission)
    return std::nullopt;
  CameraResponse result;
  result.camera_transmission = *transmission;
  if (value->find("camera_filter")) {
    result.camera_filter = response_field(*value, "camera_filter");
    if (!result.camera_filter)
      return std::nullopt;
  }
  if (const auto *knots = value->find("lightguide_efficiency")) {
    if (knots->kind() != json::Value::Kind::array)
      return std::nullopt;
    CameraIncidenceResponse guide;
    for (const auto &knot : knots->array()) {
      if (const auto *options = knot.find("interpolation")) {
        if (&knot != &knots->array().front())
          return std::nullopt;
        const auto parsed = interpolation_field(*options);
        if (!parsed)
          return std::nullopt;
        guide.interpolation = *parsed;
      }
      const auto angle = number_field(knot, "incidence_angle_deg"),
                 response = number_field(knot, "response");
      if (!angle || !response ||
          !fields_supported(knot, {"incidence_angle_deg", "response", "interpolation"}))
        return std::nullopt;
      guide.incidence_angle_deg.push_back(*angle);
      guide.response.push_back(*response);
    }
    result.lightguide_efficiency = std::move(guide);
  }
  return result.is_valid() ? std::optional{std::move(result)} : std::nullopt;
}

[[nodiscard]] inline std::optional<std::vector<ImportedDetectorSurface>>
detector_fields(const json::Value &trace, std::string_view name = "detector_surfaces") {
  const auto *entries = trace.find(name);
  if (!entries)
    return std::vector<ImportedDetectorSurface>{};
  if (entries->kind() != json::Value::Kind::array)
    return std::nullopt;
  std::vector<ImportedDetectorSurface> result;
  for (const auto &entry : entries->array()) {
    if (!fields_supported(entry, {"id", "centre_m", "normal", "tangent", "diameter_m", "shape",
                                  "enabled", "source_pixel_id", "source_type_id",
                                  "assignment_radius_m", "response_y_sign"}))
      return std::nullopt;
    for (const auto name : {"source_pixel_id", "source_type_id"})
      if (entry.find(name) && !uint_field(entry, name))
        return std::nullopt;
    if (const auto *enabled = entry.find("enabled"))
      if (enabled->kind() != json::Value::Kind::boolean || !enabled->boolean())
        return std::nullopt;
    const auto id = uint_field(entry, "id");
    const auto centre = vec3_field(entry, "centre_m");
    const auto normal = vec3_field(entry, "normal");
    const auto tangent = vec3_field(entry, "tangent");
    const auto diameter = number_field(entry, "diameter_m");
    const auto shape = facet_shape(entry);
    if (!id || *id == PhotonResultBlock::kNoSurfaceId || !centre || !normal || !tangent ||
        !diameter || !shape)
      return std::nullopt;
    result.push_back({*id, *centre, *normal, *diameter, *shape, *tangent});
    if (entry.find("assignment_radius_m")) {
      const auto radius = number_field(entry, "assignment_radius_m");
      if (!radius || *radius <= 0)
        return std::nullopt;
      result.back().assignment_radius_m = *radius;
    }
    if (entry.find("response_y_sign")) {
      const auto sign = number_field(entry, "response_y_sign");
      if (!sign || (*sign != 1 && *sign != -1))
        return std::nullopt;
      result.back().response_y_sign = *sign;
    }
  }
  return result;
}

[[nodiscard]] inline std::optional<bool>
detector_map_frame_field(const json::Value &trace,
                         std::span<const ImportedDetectorSurface> surfaces, bool has_map) {
  const auto *flag = trace.find("camera_degradation_in_detector_frame");
  if (!flag)
    return false;
  if (flag->kind() != json::Value::Kind::boolean)
    return std::nullopt;
  if (flag->boolean() && (!has_map || surfaces.empty()))
    return std::nullopt;
  if (flag->boolean())
    for (const auto &surface : surfaces) {
      const auto u = normalised_checked(surface.unit_tangent_u),
                 n = normalised_checked(surface.unit_normal);
      if (!u || !n || std::abs(dot(*u, *n)) > kEpsilon)
        return std::nullopt;
    }
  return flag->boolean();
}

[[nodiscard]] inline std::optional<CompiledDetectorAssignmentGrid>
detector_assignment_field(const json::Value &trace,
                          std::span<const ImportedDetectorSurface> surfaces) {
  const auto *value = trace.find("detector_assignment");
  if (!value || !fields_supported(*value, {"nx", "ny", "x_low_m", "x_high_m", "y_low_m", "y_high_m",
                                           "x_basis", "y_basis", "reference_plane_z_m"}))
    return std::nullopt;
  const auto nx = uint_field(*value, "nx"), ny = uint_field(*value, "ny");
  const auto xl = number_field(*value, "x_low_m"), xh = number_field(*value, "x_high_m"),
             yl = number_field(*value, "y_low_m"), yh = number_field(*value, "y_high_m");
  const auto u = vec3_field(*value, "x_basis"), v = vec3_field(*value, "y_basis");
  if (!nx || !ny || !xl || !xh || !yl || !yh || !u || !v)
    return std::nullopt;
  DetectorAssignmentGrid grid{*nx, *ny, *xl, *xh, *yl, *yh, *u, *v};
  if (value->find("reference_plane_z_m")) {
    grid.reference_plane_z_m = number_field(*value, "reference_plane_z_m");
    if (!grid.reference_plane_z_m)
      return std::nullopt;
  }
  return CompiledDetectorAssignmentGrid::compile(grid, surfaces);
}

[[nodiscard]] inline std::optional<AxisymmetricMirror>
axisymmetric_surface(const json::Value &value) {
  if (!fields_supported(value, {"vertex_z_m", "inner_radius_m", "outer_radius_m", "radial_scale_m",
                                "coefficient_m"}))
    return std::nullopt;
  const auto vertex = number_field(value, "vertex_z_m");
  const auto inner = number_field(value, "inner_radius_m");
  const auto outer = number_field(value, "outer_radius_m");
  const auto scale = number_field(value, "radial_scale_m");
  const auto *coefficients = field(value, "coefficient_m");
  if (!vertex || !inner || !outer || !scale || !coefficients ||
      coefficients->kind() != json::Value::Kind::array || coefficients->array().size() != 13)
    return std::nullopt;
  AxisymmetricMirror result{};
  result.vertex_z_m = *vertex;
  result.inner_radius_m = *inner;
  result.outer_radius_m = *outer;
  result.surface.radial_scale_m = *scale;
  for (std::size_t index = 0; index < result.surface.coefficient_m.size(); ++index) {
    const auto &coefficient = coefficients->array()[index];
    if (coefficient.kind() != json::Value::Kind::number || !std::isfinite(coefficient.number()))
      return std::nullopt;
    result.surface.coefficient_m[index] = coefficient.number();
  }
  return is_valid(result) ? std::optional{result} : std::nullopt;
}

[[nodiscard]] inline std::optional<PixelResponses> pixel_response_field(const json::Value &trace) {
  const auto *value = trace.find("pixel_responses");
  if (!value || !fields_supported(*value, {"tables", "bindings"}))
    return std::nullopt;
  const auto *tables = value->find("tables"), *bindings = value->find("bindings");
  if (!tables || !bindings || tables->kind() != json::Value::Kind::array ||
      bindings->kind() != json::Value::Kind::array)
    return std::nullopt;
  PixelResponses result;
  for (const auto &entry : tables->array()) {
    if (!fields_supported(entry, {"id", "method", "tangent_bin_width", "wavelength_bin_width_nm",
                                  "wavelength_bin_origin_nm", "angular_efficiency",
                                  "spectral_correction", "transparency", "wall_reflectivity"}))
      return std::nullopt;
    const auto id = uint_field(entry, "id");
    const auto *method = string_field(entry, "method");
    if (!id || !method || (*method != "measured" && *method != "single_reflection"))
      return std::nullopt;
    PixelResponseTable table;
    table.id = *id;
    table.method = *method == "measured" ? PixelResponseTable::Method::measured
                                         : PixelResponseTable::Method::single_reflection;
    for (const auto &[name, destination] :
         {std::pair{"tangent_bin_width", &table.tangent_bin_width},
          std::pair{"wavelength_bin_width_nm", &table.wavelength_bin_width_nm},
          std::pair{"wavelength_bin_origin_nm", &table.wavelength_bin_origin_nm},
          std::pair{"transparency", &table.transparency},
          std::pair{"wall_reflectivity", &table.wall_reflectivity}})
      if (entry.find(name)) {
        const auto number = number_field(entry, name);
        if (!number)
          return std::nullopt;
        *destination = *number;
      }
    for (const auto &[name, destination] :
         {std::pair{"angular_efficiency", &table.angular_efficiency},
          std::pair{"spectral_correction", &table.spectral_correction}})
      if (const auto *values = entry.find(name)) {
        if (values->kind() != json::Value::Kind::array)
          return std::nullopt;
        for (const auto &v : values->array()) {
          if (v.kind() != json::Value::Kind::number)
            return std::nullopt;
          destination->push_back(v.number());
        }
      }
    result.tables.push_back(std::move(table));
  }
  for (const auto &entry : bindings->array()) {
    if (!fields_supported(entry, {"id", "table_id", "cathode"}))
      return std::nullopt;
    const auto id = uint_field(entry, "id"), table = uint_field(entry, "table_id");
    if (!id || !table)
      return std::nullopt;
    PixelResponseBinding binding{*id, *table};
    if (const auto *plane = entry.find("cathode")) {
      if (!fields_supported(*plane, {"shape", "centre_m", "normal", "tangent", "diameter_m"}))
        return std::nullopt;
      const auto centre = vec3_field(*plane, "centre_m"), normal = vec3_field(*plane, "normal"),
                 tangent = vec3_field(*plane, "tangent");
      const auto diameter = number_field(*plane, "diameter_m");
      const auto shape = facet_shape(*plane);
      if (!centre || !normal || !tangent || !diameter || !shape)
        return std::nullopt;
      binding.cathode = ImportedDetectorSurface{*id, *centre, *normal, *diameter, *shape, *tangent};
    }
    result.bindings.push_back(std::move(binding));
  }
  std::sort(result.tables.begin(), result.tables.end(),
            [](const auto &a, const auto &b) { return a.id < b.id; });
  std::sort(result.bindings.begin(), result.bindings.end(),
            [](const auto &a, const auto &b) { return a.detector_id < b.detector_id; });
  return result.is_valid() ? std::optional{std::move(result)} : std::nullopt;
}

[[nodiscard]] inline std::optional<std::vector<AxisymmetricSegment>>
segment_fields(const json::Value &trace, std::string_view name) {
  const auto *entries = trace.find(name);
  if (!entries)
    return std::vector<AxisymmetricSegment>{};
  if (entries->kind() != json::Value::Kind::array)
    return std::nullopt;
  std::vector<AxisymmetricSegment> result;
  for (const auto &entry : entries->array()) {
    const auto id = uint_field(entry, "id");
    const auto *shape = string_field(entry, "shape");
    if (!id || *id == PhotonResultBlock::kNoSurfaceId || !shape)
      return std::nullopt;
    AxisymmetricSegment segment{};
    segment.id = *id;
    constexpr double radians_per_degree = std::numbers::pi / 180;
    if (*shape == "hexagon" || *shape == "square" || *shape == "circle") {
      if (!fields_supported(entry, {"id", "shape", "centre_xy_m", "diameter_m", "rotation_deg"}))
        return std::nullopt;
      const auto *centre = entry.find("centre_xy_m");
      const auto diameter = number_field(entry, "diameter_m");
      const auto rotation = number_field(entry, "rotation_deg");
      if (!centre || centre->kind() != json::Value::Kind::array || centre->array().size() != 2 ||
          centre->array()[0].kind() != json::Value::Kind::number ||
          centre->array()[1].kind() != json::Value::Kind::number ||
          !std::isfinite(centre->array()[0].number()) ||
          !std::isfinite(centre->array()[1].number()) || !diameter || *diameter <= 0 || !rotation)
        return std::nullopt;
      segment.shape = *shape == "hexagon"  ? AxisymmetricSegmentShape::hexagon
                      : *shape == "square" ? AxisymmetricSegmentShape::square
                                           : AxisymmetricSegmentShape::circle;
      segment.centre_m = {centre->array()[0].number(), centre->array()[1].number(), 0};
      segment.diameter_m = *diameter;
      segment.rotation_rad = *rotation * radians_per_degree;
    } else if (*shape == "polygon") {
      if (!fields_supported(entry, {"id", "shape", "vertices_xy_m"}))
        return std::nullopt;
      const auto *vertices = entry.find("vertices_xy_m");
      if (!vertices || vertices->kind() != json::Value::Kind::array || vertices->array().size() < 3)
        return std::nullopt;
      segment.shape = AxisymmetricSegmentShape::polygon;
      for (const auto &point : vertices->array()) {
        if (point.kind() != json::Value::Kind::array || point.array().size() != 2)
          return std::nullopt;
        for (const auto &v : point.array())
          if (v.kind() != json::Value::Kind::number || !std::isfinite(v.number()))
            return std::nullopt;
        segment.vertices.push_back({point.array()[0].number(), point.array()[1].number(), 0});
      }
      double sign = 0;
      for (std::size_t i = 0; i < segment.vertices.size(); ++i) {
        const auto a = segment.vertices[(i + 1) % segment.vertices.size()] - segment.vertices[i];
        const auto b = segment.vertices[(i + 2) % segment.vertices.size()] -
                       segment.vertices[(i + 1) % segment.vertices.size()];
        const double turn = cross(a, b).z;
        if (std::abs(turn) <= kEpsilon || (sign && sign * turn < 0))
          return std::nullopt;
        sign = turn;
      }
    } else if (*shape == "annular_sector") {
      if (!fields_supported(entry, {"id", "shape", "inner_radius_m", "outer_radius_m", "start_deg",
                                    "span_deg", "gap_m", "gap_at_start"}))
        return std::nullopt;
      const auto inner = number_field(entry, "inner_radius_m"),
                 outer = number_field(entry, "outer_radius_m");
      const auto start = number_field(entry, "start_deg"), span = number_field(entry, "span_deg");
      const auto gap = number_field(entry, "gap_m");
      if (!inner || !outer || !start || !span || !gap || *inner < 0 || *outer <= *inner ||
          *span <= 0 || *span > 360 || *gap < 0)
        return std::nullopt;
      segment.shape = AxisymmetricSegmentShape::annular_sector;
      segment.inner_radius_m = *inner;
      segment.outer_radius_m = *outer;
      segment.start_rad = *start * radians_per_degree;
      segment.span_rad = *span * radians_per_degree;
      segment.gap_m = *gap;
      if (const auto *edge = entry.find("gap_at_start")) {
        if (edge->kind() != json::Value::Kind::boolean)
          return std::nullopt;
        segment.gap_at_start = edge->boolean();
      }
    } else
      return std::nullopt;
    result.push_back(segment);
  }
  return result;
}

[[nodiscard]] inline std::optional<json::Value> read_json(const std::string &path) {
  std::optional<json::Value> root;
  {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input)
      return std::nullopt;
    const auto size = input.tellg();
    if (size < 0)
      return std::nullopt;
    std::string document(static_cast<std::size_t>(size), '\0');
    input.seekg(0);
    if (!input.read(document.data(), static_cast<std::streamsize>(document.size())))
      return std::nullopt;
    root = json::Parser{document}.parse();
  }
  if (!root || root->kind() != json::Value::Kind::object || !string_field(*root, "format") ||
      *string_field(*root, "format") != "obdeect.compiled-optical-model.v1")
    return std::nullopt;
  const auto *hash = string_field(*root, "optical_model_sha256");
  if (!hash || hash->size() != 64 ||
      std::any_of(hash->begin(), hash->end(),
                  [](const char character) {
                    return !((character >= '0' && character <= '9') ||
                             (character >= 'a' && character <= 'f'));
                  }) ||
      canonical_json_hash(*root) != *hash) {
    return std::nullopt;
  }
  return root;
}

// Read the single versioned JSON artifact emitted by obdeect-compile-optical-model.
// Its trace_model section is directly consumable; the remainder is retained for
// provenance, validation and inspection.
inline std::optional<CompiledSegmentedOpticalModel>
segmented_optical_model_from_json(const json::Value &root) {
  const auto model_provenance = detail::provenance(root);
  const auto *trace = detail::field(root, "trace_model");
  if (!model_provenance || !trace || trace->kind() != json::Value::Kind::object ||
      !detail::string_field(*trace, "kind") ||
      *detail::string_field(*trace, "kind") != "segmented" ||
      !detail::fields_supported(
          *trace, {"kind", "primary_facets", "detector_surfaces", "cylinder_obscurers",
                   "primary_reflectivity", "primary_scatter", "camera_response",
                   "incoming_obscurer_planes", "opaque_obscurers", "primary_degradation",
                   "camera_degradation", "primary_degradation_in_facet_frame", "pixel_responses",
                   "telescope_transmission", "propagation_group_index", "detector_assignment",
                   "camera_degradation_in_detector_frame"}))
    return std::nullopt;
  const auto *facet_values = detail::field(*trace, "primary_facets");
  const auto *detector_values = detail::field(*trace, "detector_surfaces");
  const auto *obscurer_values = detail::field(*trace, "cylinder_obscurers");
  if (!facet_values || !detector_values || !obscurer_values ||
      facet_values->kind() != json::Value::Kind::array ||
      detector_values->kind() != json::Value::Kind::array ||
      obscurer_values->kind() != json::Value::Kind::array)
    return std::nullopt;
  std::vector<ImportedFacet> facets;
  for (const auto &value : facet_values->array()) {
    if (!detail::fields_supported(value, {"id", "centre_m", "normal", "tangent", "diameter_m",
                                          "focal_length_m", "shape", "response_basis_u",
                                          "response_basis_v"}))
      return std::nullopt;
    const auto id = detail::uint_field(value, "id");
    const auto centre = detail::vec3_field(value, "centre_m");
    const auto normal = detail::vec3_field(value, "normal");
    const auto tangent = detail::vec3_field(value, "tangent");
    const auto diameter = detail::number_field(value, "diameter_m");
    const auto focal_length = detail::number_field(value, "focal_length_m");
    const auto shape = detail::facet_shape(value);
    if (!id || *id == PhotonResultBlock::kNoSurfaceId || !centre || !normal || !tangent ||
        !diameter || !focal_length || !shape)
      return std::nullopt;
    facets.push_back(
        {*id, *centre, *normal, *diameter, *focal_length, *shape, *tangent, 2.0 * *focal_length});
    if (value.find("response_basis_u") || value.find("response_basis_v")) {
      const auto u = detail::vec3_field(value, "response_basis_u"),
                 v = detail::vec3_field(value, "response_basis_v");
      if (!u || !v || std::abs(norm(*u) - 1) > kEpsilon || std::abs(norm(*v) - 1) > kEpsilon ||
          std::abs(dot(*u, *v)) > kEpsilon || std::abs(dot(*u, *normal)) > kEpsilon ||
          std::abs(dot(*v, *normal)) > kEpsilon)
        return std::nullopt;
      facets.back().response_basis_u = *u;
      facets.back().response_basis_v = *v;
    }
  }
  const auto detectors = detail::detector_fields(*trace);
  if (!detectors)
    return std::nullopt;
  std::vector<ImportedCylinderObscurer> obscurers;
  for (const auto &value : obscurer_values->array()) {
    if (!detail::fields_supported(value,
                                  {"id", "first_endpoint_m", "second_endpoint_m", "diameter_m"}))
      return std::nullopt;
    const auto id = detail::uint_field(value, "id");
    const auto first = detail::vec3_field(value, "first_endpoint_m");
    const auto second = detail::vec3_field(value, "second_endpoint_m");
    const auto diameter = detail::number_field(value, "diameter_m");
    if (!id || *id == PhotonResultBlock::kNoSurfaceId || !first || !second || !diameter)
      return std::nullopt;
    obscurers.push_back({*id, *first, *second, *diameter});
  }
  const auto *reflectivity_value = detail::field(*trace, "primary_reflectivity");
  const auto reflectivity = detail::response_field(*trace, "primary_reflectivity");
  if ((reflectivity_value && reflectivity_value->kind() != json::Value::Kind::array) ||
      (reflectivity_value && !reflectivity_value->array().empty() && !reflectivity) ||
      facets.empty() || detectors->empty())
    return std::nullopt;
  auto compiled = compile_segmented_optical_model({*model_provenance, std::move(facets), *detectors,
                                                   std::move(obscurers), std::move(reflectivity)});
  if (!compiled)
    return std::nullopt;
  auto planes = compile_detector_planes(*detectors);
  if (!planes)
    return std::nullopt;
  compiled->detector_planes = std::make_shared<const CompiledDetectorPlanes>(std::move(*planes));
  const auto incoming = detail::detector_fields(*trace, "incoming_obscurer_planes");
  if (!incoming)
    return std::nullopt;
  if (!incoming->empty()) {
    auto shadows = compile_detector_planes(*incoming);
    if (!shadows)
      return std::nullopt;
    for (const auto &surface : *incoming) {
      for (const auto &facet : compiled->primary_facets)
        if (facet.id == surface.id)
          return std::nullopt;
      for (const auto &detector : *detectors)
        if (detector.id == surface.id)
          return std::nullopt;
      for (const auto &cylinder : compiled->cylinder_obscurers)
        if (cylinder.id == surface.id)
          return std::nullopt;
    }
    compiled->incoming_obscurer_planes =
        std::make_shared<const CompiledDetectorPlanes>(std::move(*shadows));
  }
  auto opaque = detail::opaque_fields(*trace);
  if (!opaque)
    return std::nullopt;
  compiled->opaque_obscurers = std::move(*opaque);
  if (!is_valid(*compiled))
    return std::nullopt;
  if (compiled->incoming_obscurer_planes)
    for (const auto &a : compiled->incoming_obscurer_planes->surfaces())
      for (const auto &b : compiled->opaque_obscurers)
        if (a.id == b.id)
          return std::nullopt;
  for (const auto &[name, destination] :
       {std::pair{"primary_degradation", &compiled->primary_degradation},
        std::pair{"camera_degradation", &compiled->camera_degradation}}) {
    if (trace->find(name)) {
      *destination = detail::spatial_field(*trace, name);
      if (!*destination)
        return std::nullopt;
    }
  }
  if (const auto *local = trace->find("primary_degradation_in_facet_frame")) {
    if (local->kind() != json::Value::Kind::boolean)
      return std::nullopt;
    compiled->primary_degradation_in_facet_frame = local->boolean();
    if (local->boolean())
      for (const auto &facet : compiled->primary_facets)
        if (norm(facet.response_basis_u) <= kEpsilon || norm(facet.response_basis_v) <= kEpsilon)
          return std::nullopt;
  }
  if (trace->find("pixel_responses")) {
    auto responses = detail::pixel_response_field(*trace);
    if (!responses)
      return std::nullopt;
    if (responses->bindings.size() != detectors->size())
      return std::nullopt;
    for (const auto &binding : responses->bindings)
      if (std::none_of(detectors->begin(), detectors->end(),
                       [&](const auto &d) { return d.id == binding.detector_id; }))
        return std::nullopt;
    if (!responses->is_valid())
      return std::nullopt;
    compiled->pixel_responses = std::make_shared<const PixelResponses>(std::move(*responses));
  }
  compiled->primary_scatter = detail::scatter_field(*trace, "primary_scatter");
  compiled->camera_response = detail::camera_response_field(*trace, "camera_response");
  const auto camera_frame = detail::detector_map_frame_field(
      *trace, *detectors, compiled->camera_degradation.has_value());
  if (!camera_frame)
    return std::nullopt;
  compiled->camera_degradation_in_detector_frame = *camera_frame;
  if (trace->find("detector_assignment")) {
    auto grid = detail::detector_assignment_field(*trace, *detectors);
    if (!grid || !grid->description().reference_plane_z_m)
      return std::nullopt;
    compiled->detector_assignment =
        std::make_shared<const CompiledDetectorAssignmentGrid>(std::move(*grid));
  }
  compiled->telescope_transmission = detail::telescope_transmission_field(*trace);
  if (trace->find("propagation_group_index")) {
    const auto index = detail::number_field(*trace, "propagation_group_index");
    if (!index || *index < 1)
      return std::nullopt;
    compiled->propagation_group_index = *index;
  }
  if (trace->find("telescope_transmission") && !compiled->telescope_transmission)
    return std::nullopt;
  if ((trace->find("primary_scatter") && !compiled->primary_scatter) ||
      (compiled->primary_scatter &&
       compiled->primary_scatter->method != MirrorScatterMethod::outgoing_angles) ||
      (trace->find("camera_response") && !compiled->camera_response))
    return std::nullopt;
  if (!is_valid(*compiled))
    return std::nullopt;
  return compiled;
}

[[nodiscard]] inline std::optional<AxisymmetricOpticalModel>
axisymmetric_optical_model_from_json(const json::Value &root) {
  const auto model_provenance = detail::provenance(root);
  const auto *trace = detail::field(root, "trace_model");
  if (!model_provenance || !trace || trace->kind() != json::Value::Kind::object ||
      !detail::string_field(*trace, "kind") ||
      *detail::string_field(*trace, "kind") != "axisymmetric" ||
      !detail::fields_supported(*trace, {"kind",
                                         "primary",
                                         "secondary",
                                         "detector",
                                         "primary_segments",
                                         "secondary_segments",
                                         "primary_surface_id",
                                         "secondary_surface_id",
                                         "detector_surface_id",
                                         "block_incoming_secondary",
                                         "detector_surfaces",
                                         "primary_reflectivity",
                                         "secondary_reflectivity",
                                         "primary_scatter",
                                         "secondary_scatter",
                                         "camera_response",
                                         "primary_to_secondary_planes",
                                         "primary_to_secondary_cylinders",
                                         "incoming_obscurer_planes",
                                         "opaque_obscurers",
                                         "primary_degradation",
                                         "secondary_degradation",
                                         "camera_degradation",
                                         "pixel_responses",
                                         "telescope_transmission",
                                         "propagation_group_index",
                                         "detector_assignment",
                                         "camera_degradation_in_detector_frame"}))
    return std::nullopt;
  const auto *primary = detail::field(*trace, "primary");
  const auto *secondary = detail::field(*trace, "secondary");
  const auto *detector = detail::field(*trace, "detector");
  if (!primary || !secondary || !detector)
    return std::nullopt;
  const auto primary_surface = detail::axisymmetric_surface(*primary);
  const auto secondary_surface = detail::axisymmetric_surface(*secondary);
  const auto detector_surface = detail::axisymmetric_surface(*detector);
  const auto *primary_response = detail::field(*trace, "primary_reflectivity");
  const auto *secondary_response = detail::field(*trace, "secondary_reflectivity");
  const auto primary_reflectivity = detail::response_field(*trace, "primary_reflectivity");
  const auto secondary_reflectivity = detail::response_field(*trace, "secondary_reflectivity");
  if ((primary_response && primary_response->kind() != json::Value::Kind::array) ||
      (secondary_response && secondary_response->kind() != json::Value::Kind::array) ||
      (primary_response && !primary_response->array().empty() && !primary_reflectivity) ||
      (secondary_response && !secondary_response->array().empty() && !secondary_reflectivity) ||
      !primary_surface || !secondary_surface || !detector_surface)
    return std::nullopt;
  auto primary_segments = detail::segment_fields(*trace, "primary_segments");
  auto secondary_segments = detail::segment_fields(*trace, "secondary_segments");
  if (!primary_segments || !secondary_segments)
    return std::nullopt;
  std::unordered_set<std::uint32_t> segment_ids;
  for (const auto *segments : {&*primary_segments, &*secondary_segments})
    for (const auto &segment : *segments)
      if (!segment_ids.insert(segment.id).second)
        return std::nullopt;
  AxisymmetricOpticalModel model{*model_provenance,
                                 *primary_surface,
                                 *secondary_surface,
                                 *detector_surface,
                                 std::move(primary_reflectivity),
                                 std::move(secondary_reflectivity)};
  model.primary_segments = std::move(*primary_segments);
  model.secondary_segments = std::move(*secondary_segments);
  for (const auto &name : {"primary_surface_id", "secondary_surface_id", "detector_surface_id"})
    if (trace->find(name)) {
      const auto id = detail::uint_field(*trace, name);
      if (!id || *id == PhotonResultBlock::kNoSurfaceId)
        return std::nullopt;
    }
  model.primary_surface_id =
      detail::uint_field(*trace, "primary_surface_id").value_or(PhotonResultBlock::kNoSurfaceId);
  model.secondary_surface_id =
      detail::uint_field(*trace, "secondary_surface_id").value_or(PhotonResultBlock::kNoSurfaceId);
  model.detector_surface_id =
      detail::uint_field(*trace, "detector_surface_id").value_or(PhotonResultBlock::kNoSurfaceId);
  for (const auto id :
       {model.primary_surface_id, model.secondary_surface_id, model.detector_surface_id})
    if (id != PhotonResultBlock::kNoSurfaceId && !segment_ids.insert(id).second)
      return std::nullopt;
  if (const auto *shadow = trace->find("block_incoming_secondary")) {
    if (shadow->kind() != json::Value::Kind::boolean)
      return std::nullopt;
    model.block_incoming_secondary = shadow->boolean();
  }
  const auto detectors = detail::detector_fields(*trace);
  if (!detectors)
    return std::nullopt;
  if (!detectors->empty()) {
    for (const auto &surface : *detectors)
      if (!segment_ids.insert(surface.id).second)
        return std::nullopt;
    model.detector_planes = compile_detector_planes(*detectors);
    if (!model.detector_planes)
      return std::nullopt;
  }
  if (trace->find("detector_assignment")) {
    model.detector_assignment = detail::detector_assignment_field(*trace, *detectors);
    if (!model.detector_assignment || model.detector_assignment->description().reference_plane_z_m)
      return std::nullopt;
  }
  model.primary_scatter = detail::scatter_field(*trace, "primary_scatter");
  const auto obscurer_planes = detail::detector_fields(*trace, "primary_to_secondary_planes");
  if (!obscurer_planes)
    return std::nullopt;
  if (!obscurer_planes->empty()) {
    for (const auto &plane : *obscurer_planes)
      if (!segment_ids.insert(plane.id).second)
        return std::nullopt;
    model.primary_to_secondary_planes = compile_detector_planes(*obscurer_planes);
    if (!model.primary_to_secondary_planes)
      return std::nullopt;
  }
  if (const auto *entries = trace->find("primary_to_secondary_cylinders")) {
    if (entries->kind() != json::Value::Kind::array)
      return std::nullopt;
    for (const auto &entry : entries->array()) {
      if (!detail::fields_supported(entry,
                                    {"id", "first_endpoint_m", "second_endpoint_m", "diameter_m"}))
        return std::nullopt;
      const auto id = detail::uint_field(entry, "id");
      const auto first = detail::vec3_field(entry, "first_endpoint_m");
      const auto second = detail::vec3_field(entry, "second_endpoint_m");
      const auto diameter = detail::number_field(entry, "diameter_m");
      if (!id || *id == PhotonResultBlock::kNoSurfaceId || !first || !second || !diameter ||
          !segment_ids.insert(*id).second)
        return std::nullopt;
      const ImportedCylinderObscurer cylinder{*id, *first, *second, *diameter};
      if (!is_valid(cylinder))
        return std::nullopt;
      model.primary_to_secondary_cylinders.push_back(cylinder);
    }
  }
  for (const auto &[name, destination] :
       {std::pair{"primary_degradation", &model.primary_degradation},
        std::pair{"secondary_degradation", &model.secondary_degradation},
        std::pair{"camera_degradation", &model.camera_degradation}}) {
    if (trace->find(name)) {
      *destination = detail::spatial_field(*trace, name);
      if (!*destination)
        return std::nullopt;
    }
  }
  if (trace->find("pixel_responses")) {
    auto responses = detail::pixel_response_field(*trace);
    if (!responses || responses->bindings.size() != detectors->size())
      return std::nullopt;
    for (const auto &binding : responses->bindings)
      if (std::none_of(detectors->begin(), detectors->end(),
                       [&](const auto &d) { return d.id == binding.detector_id; }))
        return std::nullopt;
    if (!responses->is_valid())
      return std::nullopt;
    model.pixel_responses = std::make_shared<const PixelResponses>(std::move(*responses));
  }
  auto opaque = detail::opaque_fields(*trace);
  if (!opaque)
    return std::nullopt;
  for (const auto &surface : *opaque)
    if (!segment_ids.insert(surface.id).second)
      return std::nullopt;
  model.opaque_obscurers = std::move(*opaque);
  const auto incoming_planes = detail::detector_fields(*trace, "incoming_obscurer_planes");
  if (!incoming_planes)
    return std::nullopt;
  if (!incoming_planes->empty()) {
    for (const auto &plane : *incoming_planes)
      if (!segment_ids.insert(plane.id).second)
        return std::nullopt;
    model.incoming_obscurer_planes = compile_detector_planes(*incoming_planes);
    if (!model.incoming_obscurer_planes)
      return std::nullopt;
  }
  model.secondary_scatter = detail::scatter_field(*trace, "secondary_scatter");
  model.camera_response = detail::camera_response_field(*trace, "camera_response");
  const auto camera_frame =
      detail::detector_map_frame_field(*trace, *detectors, model.camera_degradation.has_value());
  if (!camera_frame)
    return std::nullopt;
  model.camera_degradation_in_detector_frame = *camera_frame;
  model.telescope_transmission = detail::telescope_transmission_field(*trace);
  if (trace->find("propagation_group_index")) {
    const auto index = detail::number_field(*trace, "propagation_group_index");
    if (!index || *index < 1)
      return std::nullopt;
    model.propagation_group_index = *index;
  }
  if (trace->find("telescope_transmission") && !model.telescope_transmission)
    return std::nullopt;
  if ((trace->find("primary_scatter") && !model.primary_scatter) ||
      (trace->find("secondary_scatter") && !model.secondary_scatter) ||
      (trace->find("camera_response") && !model.camera_response))
    return std::nullopt;
  return model;
}

// Explicit nonsequential media/interfaces are generic data. Existing CTAO
// compatibility response tables cannot establish these physical interfaces.
[[nodiscard]] inline std::optional<CompiledOpticalModel>
nonsequential_optical_model_from_json(const json::Value &root) {
  const auto provenance = detail::provenance(root);
  const auto *trace = root.find("trace_model");
  using detail::fields_supported;
  if (!provenance || !trace || !detail::string_field(*trace, "kind") ||
      *detail::string_field(*trace, "kind") != "nonsequential" ||
      !fields_supported(
          *trace, {"kind", "materials", "surfaces", "max_interactions", "entrance_medium_id"}))
    return std::nullopt;
  const auto *materials = trace->find("materials");
  const auto *surfaces = trace->find("surfaces");
  const auto limit = detail::uint_field(*trace, "max_interactions");
  const auto entrance_medium = detail::uint_field(*trace, "entrance_medium_id");
  if (!materials || materials->kind() != json::Value::Kind::array || !surfaces ||
      surfaces->kind() != json::Value::Kind::array || !limit || !entrance_medium)
    return std::nullopt;
  ImportedOpticalModel model{*provenance, {}, *limit, {}, *entrance_medium};
  const auto curve = [&](const json::Value &material,
                         std::string_view name) -> std::optional<MaterialCurve> {
    const auto *values = material.find(name);
    if (!values || values->kind() != json::Value::Kind::array)
      return std::nullopt;
    MaterialCurve result;
    for (const auto &knot : values->array()) {
      const auto wavelength = detail::number_field(knot, "wavelength_nm");
      const auto value = detail::number_field(knot, "value");
      if (!wavelength || !value || !fields_supported(knot, {"wavelength_nm", "value"}))
        return std::nullopt;
      result.wavelength_nm.push_back(*wavelength);
      result.value.push_back(*value);
    }
    return result;
  };
  for (const auto &material : materials->array()) {
    const auto id = detail::uint_field(material, "id");
    auto phase = curve(material, "phase_index");
    auto group = curve(material, "group_index");
    auto absorption = curve(material, "absorption_per_m");
    if (!id || !phase || !group || !absorption ||
        !fields_supported(material, {"id", "phase_index", "group_index", "absorption_per_m"}))
      return std::nullopt;
    model.materials.push_back({*id, std::move(*phase), std::move(*group), std::move(*absorption)});
  }
  for (const auto &surface : surfaces->array()) {
    const auto id = detail::uint_field(surface, "id");
    const auto diameter = detail::number_field(surface, "diameter_m");
    const auto shape = detail::facet_shape(surface);
    const auto *role = detail::string_field(surface, "role");
    const auto *frame = surface.find("frame");
    if (!id || !diameter || !shape || !role || !frame ||
        !fields_supported(surface, {"id", "diameter_m", "shape", "role", "frame", "sag",
                                    "front_medium_id", "back_medium_id", "transmission",
                                    "transmission_semantics", "material_id", "reflectivity"}) ||
        !fields_supported(*frame, {"origin_m", "x_axis", "y_axis", "z_axis"}))
      return std::nullopt;
    const auto origin = detail::vec3_field(*frame, "origin_m");
    const auto x = detail::vec3_field(*frame, "x_axis");
    const auto y = detail::vec3_field(*frame, "y_axis");
    const auto z = detail::vec3_field(*frame, "z_axis");
    if (!origin || !x || !y || !z)
      return std::nullopt;
    OpticalSurfaceRecord record{*id, {*origin, *x, *y, *z}, *shape, *diameter};
    if (*role == "mirror")
      record.role = SurfaceRole::mirror;
    else if (*role == "detector")
      record.role = SurfaceRole::detector;
    else if (*role == "obscurer")
      record.role = SurfaceRole::obscurer;
    else if (*role == "refractive_interface")
      record.role = SurfaceRole::refractive_interface;
    else
      return std::nullopt;
    if (const auto *sag = surface.find("sag")) {
      const auto parsed = detail::axisymmetric_surface(*sag);
      if (!parsed || parsed->vertex_z_m != 0 || parsed->outer_radius_m != *diameter / 2 ||
          !fields_supported(*sag, {"vertex_z_m", "inner_radius_m", "outer_radius_m",
                                   "radial_scale_m", "coefficient_m"}))
        return std::nullopt;
      record.sag = parsed->surface;
      record.inner_radius_m = parsed->inner_radius_m;
    }
    for (const auto &[name, destination] : {std::pair{"front_medium_id", &record.front_medium_id},
                                            std::pair{"back_medium_id", &record.back_medium_id}}) {
      if (surface.find(name)) {
        const auto medium = detail::uint_field(surface, name);
        if (!medium)
          return std::nullopt;
        *destination = *medium;
      } else if (record.role == SurfaceRole::refractive_interface)
        return std::nullopt;
    }
    if (surface.find("transmission")) {
      record.transmission = detail::response_field(surface, "transmission");
      const auto *semantics = detail::string_field(surface, "transmission_semantics");
      if (!record.transmission || !semantics)
        return std::nullopt;
      if (*semantics == "additional_coating")
        record.transmission_semantics = InterfaceTransmissionSemantics::additional_coating;
      else if (*semantics == "complete_interface")
        record.transmission_semantics = InterfaceTransmissionSemantics::complete_interface;
      else
        return std::nullopt;
    } else if (surface.find("transmission_semantics"))
      return std::nullopt;
    for (const auto name : {"transmission", "reflectivity"})
      if (const auto *response = surface.find(name)) {
        if (response->kind() != json::Value::Kind::array)
          return std::nullopt;
        for (const auto &knot : response->array())
          if (!fields_supported(
                  knot, {"wavelength_nm", "response", "incidence_angle_deg", "interpolation"}))
            return std::nullopt;
      }
    if (surface.find("reflectivity")) {
      record.reflectivity = detail::response_field(surface, "reflectivity");
      if (!record.reflectivity)
        return std::nullopt;
    }
    if (surface.find("material_id")) {
      const auto material = detail::uint_field(surface, "material_id");
      if (!material)
        return std::nullopt;
      record.material_id = *material;
    }
    model.surfaces.push_back(std::move(record));
  }
  return compile_optical_model(model);
}

} // namespace detail

struct LoadedOpticalModel {
  std::optional<CompiledSegmentedOpticalModel> segmented;
  std::optional<AxisymmetricOpticalModel> axisymmetric;
  std::optional<CompiledOpticalModel> nonsequential;
  bool production_ready{};
  std::optional<double> imaging_plane_z_m{};
  std::uint64_t detector_configuration_seed{};
};

// Verify the artifact once, compile its declared trace model and release the
// JSON tree before tracing. Readiness is covered by the same content hash.
[[nodiscard]] inline std::optional<LoadedOpticalModel> read_optical_model(const std::string &path) {
  const auto root = detail::read_json(path);
  if (!root)
    return std::nullopt;
  const auto *trace = root->find("trace_model");
  const auto *kind = trace ? detail::string_field(*trace, "kind") : nullptr;
  if (!kind)
    return std::nullopt;
  LoadedOpticalModel loaded;
  if (const auto *seeds = root->find("random_seeds")) {
    if (seeds->kind() != json::Value::Kind::object ||
        !detail::fields_supported(*seeds, {"detector_configuration_seed"}))
      return std::nullopt;
    const auto seed = detail::uint_field(*seeds, "detector_configuration_seed");
    if (!seed)
      return std::nullopt;
    loaded.detector_configuration_seed = *seed;
  }
  loaded.imaging_plane_z_m = detail::number_field(*root, "detector_vertex_z_m");
  if (*kind == "segmented")
    loaded.segmented = detail::segmented_optical_model_from_json(*root);
  else if (*kind == "axisymmetric")
    loaded.axisymmetric = detail::axisymmetric_optical_model_from_json(*root);
  else if (*kind == "nonsequential") {
    if (auto model = detail::nonsequential_optical_model_from_json(*root))
      loaded.nonsequential.emplace(std::move(*model));
  }
  if (!loaded.segmented && !loaded.axisymmetric && !loaded.nonsequential)
    return std::nullopt;
  const auto *report = root->find("report");
  const auto *readiness = report ? report->find("production_trace_ready") : nullptr;
  const auto *blockers = report ? report->find("trace_blockers") : nullptr;
  loaded.production_ready =
      readiness && readiness->kind() == json::Value::Kind::boolean && readiness->boolean() &&
      blockers && blockers->kind() == json::Value::Kind::array && blockers->array().empty();
  return loaded;
}

[[nodiscard]] inline std::optional<CompiledSegmentedOpticalModel>
read_segmented_optical_model(const std::string &path) {
  const auto root = detail::read_json(path);
  return root ? detail::segmented_optical_model_from_json(*root) : std::nullopt;
}

[[nodiscard]] inline std::optional<AxisymmetricOpticalModel>
read_axisymmetric_optical_model(const std::string &path) {
  const auto root = detail::read_json(path);
  return root ? detail::axisymmetric_optical_model_from_json(*root) : std::nullopt;
}

[[nodiscard]] inline std::optional<CompiledOpticalModel>
read_nonsequential_optical_model(const std::string &path) {
  const auto root = detail::read_json(path);
  return root ? detail::nonsequential_optical_model_from_json(*root) : std::nullopt;
}

[[nodiscard]] inline PathRecord
trace_axisymmetric_optical_model(const Ray &input, std::uint64_t photon_id,
                                 const AxisymmetricOpticalModel &optical_model,
                                 double wavelength_nm = 400.0, std::uint64_t ray_tracing_seed = 0) {
  PathRecord record{};
  record.photon_id = photon_id;
  record.wavelength_nm = wavelength_nm;
  record.points_m[0] = input.position_m;
  record.point_count = 1;
  const auto direction = normalised_checked(input.direction);
  if (!direction || !std::isfinite(input.position_m.x) || !std::isfinite(input.position_m.y) ||
      !std::isfinite(input.position_m.z) || !std::isfinite(wavelength_nm) || wavelength_nm <= 0) {
    record.status = PhotonStatus::invalid_input;
    return record;
  }
  if (optical_model.telescope_transmission)
    record.surviving_throughput = *optical_model.telescope_transmission->at_unchecked(*direction);
  Ray ray{input.position_m, *direction};
  bool primary_failure = false;
  const auto primary =
      intersect_segmented_asphere(ray, optical_model.primary, optical_model.primary_segments,
                                  optical_model.primary_surface_id, &primary_failure);
  if (primary_failure) {
    record.status = PhotonStatus::intersection_failure;
    record.final_direction = ray.direction;
    return record;
  }
  {
    DetectorSurfaceHit shadow{};
    bool has_shadow = false;
    if (optical_model.incoming_obscurer_planes) {
      const auto incoming_shadow = optical_model.incoming_obscurer_planes->intersect(ray);
      if (incoming_shadow) {
        shadow = *incoming_shadow;
        has_shadow = true;
      }
    }
    const auto opaque = intersect_opaque_surfaces(ray, optical_model.opaque_obscurers);
    if (opaque && (!has_shadow || opaque->distance_m < shadow.distance_m)) {
      shadow = DetectorSurfaceHit{opaque->surface_id, opaque->distance_m, opaque->point_m,
                                  opaque->unit_normal};
      has_shadow = true;
    }
    if (has_shadow && (!primary || shadow.distance_m < primary->distance_m)) {
      record.points_m[1] = shadow.point_m;
      record.point_count = 2;
      record.path_length_m = shadow.distance_m;
      record.status = PhotonStatus::blocked_obscurer;
      record.final_direction = ray.direction;
      record.terminal_surface_id = shadow.surface_id;
      record.interaction_surface_ids[0] = shadow.surface_id;
      record.interaction_kinds[0] = OpticalInteractionKind::obscurer;
      record.interaction_normals[0] = shadow.unit_normal;
      record.interaction_incoming_directions[0] = ray.direction;
      record.interaction_outgoing_directions[0] = ray.direction;
      return record;
    }
  }
  if (optical_model.block_incoming_secondary) {
    bool shadow_failure = false;
    const auto shadow =
        intersect_segmented_asphere(ray, optical_model.secondary, optical_model.secondary_segments,
                                    optical_model.secondary_surface_id, &shadow_failure);
    if (shadow_failure) {
      record.status = PhotonStatus::intersection_failure;
      record.final_direction = ray.direction;
      return record;
    }
    if (shadow && (!primary || shadow->distance_m < primary->distance_m)) {
      record.points_m[1] = shadow->point_m;
      record.point_count = 2;
      record.path_length_m = shadow->distance_m;
      record.status = PhotonStatus::blocked_obscurer;
      record.final_direction = ray.direction;
      record.terminal_surface_id =
          *axisymmetric_segment_id(optical_model.secondary_segments, optical_model.secondary,
                                   shadow->point_m, optical_model.secondary_surface_id);
      record.interaction_surface_ids[0] = record.terminal_surface_id;
      record.interaction_kinds[0] = OpticalInteractionKind::obscurer;
      record.interaction_normals[0] = shadow->unit_normal;
      record.interaction_incoming_directions[0] = ray.direction;
      record.interaction_outgoing_directions[0] = ray.direction;
      return record;
    }
  }
  if (!primary) {
    record.status = PhotonStatus::missed_primary;
    record.final_direction = ray.direction;
    return record;
  }
  record.points_m[1] = primary->point_m;
  record.terminal_surface_id =
      *axisymmetric_segment_id(optical_model.primary_segments, optical_model.primary,
                               primary->point_m, optical_model.primary_surface_id);
  record.interaction_surface_ids[0] = record.terminal_surface_id;
  record.interaction_normals[0] = primary->unit_normal;
  record.interaction_incoming_directions[0] = ray.direction;
  record.point_count = 2;
  record.path_length_m = primary->distance_m;
  record.incidence_primary_deg =
      std::acos(std::clamp(std::abs(dot(ray.direction, primary->unit_normal)), 0.0, 1.0)) * 180.0 /
      std::numbers::pi;
  const auto after_primary = optical_model.primary_scatter ? [&]() -> std::optional<Ray> {
    const auto direction = reflect_with_scatter(ray.direction, primary->unit_normal, {1, 0, 0},
                                                *optical_model.primary_scatter, photon_id,
                                                record.terminal_surface_id, ray_tracing_seed);
    return direction ? std::optional{Ray{primary->point_m, *direction}} : std::nullopt;
  }()
      : reflect(ray, *primary);
  if (!after_primary) {
    record.status = PhotonStatus::invalid_input;
    return record;
  }
  ray = *after_primary;
  record.final_direction = ray.direction;
  record.interaction_outgoing_directions[0] = ray.direction;
  const auto primary_response = optical_model.primary_reflectivity
                                    ? optical_model.primary_reflectivity->at_unchecked(
                                          wavelength_nm, record.incidence_primary_deg)
                                    : std::optional<double>{1.0};
  if (!primary_response) {
    record.status = PhotonStatus::invalid_input;
    return record;
  }
  record.surviving_throughput *= *primary_response;
  if (optical_model.primary_degradation) {
    const auto response = optical_model.primary_degradation->at_point_unchecked(primary->point_m);
    if (!response) {
      record.status = PhotonStatus::invalid_input;
      return record;
    }
    record.surviving_throughput *= *response;
  }
  record.interaction_throughput[0] = record.surviving_throughput;
  bool secondary_failure = false;
  const auto secondary =
      intersect_segmented_asphere(ray, optical_model.secondary, optical_model.secondary_segments,
                                  optical_model.secondary_surface_id, &secondary_failure);
  if (secondary_failure) {
    record.status = PhotonStatus::intersection_failure;
    record.final_direction = ray.direction;
    return record;
  }
  std::optional<DetectorSurfaceHit> obscurer{};
  if (optical_model.primary_to_secondary_planes)
    obscurer = optical_model.primary_to_secondary_planes->intersect(ray);
  for (const auto &cylinder : optical_model.primary_to_secondary_cylinders) {
    const auto distance = intersect_closed_finite_cylinder(
        ray, cylinder.first_endpoint_m, cylinder.second_endpoint_m, cylinder.diameter_m * 0.5);
    if (!distance)
      continue;
    if (obscurer) {
      const auto &current = *obscurer;
      if (*distance > current.distance_m ||
          (*distance == current.distance_m && cylinder.id >= current.surface_id))
        continue;
    }
    const Vec3 point = ray.position_m + ray.direction * *distance;
    const Vec3 axis = *normalised_checked(cylinder.second_endpoint_m - cylinder.first_endpoint_m);
    const double projection = dot(point - cylinder.first_endpoint_m, axis);
    const double length = norm(cylinder.second_endpoint_m - cylinder.first_endpoint_m);
    const Vec3 normal =
        projection <= kEpsilon ? axis * -1
        : projection >= length - kEpsilon
            ? axis
            : *normalised_checked(point - cylinder.first_endpoint_m - axis * projection);
    obscurer = DetectorSurfaceHit{cylinder.id, *distance, point, normal};
  }
  const auto opaque = intersect_opaque_surfaces(ray, optical_model.opaque_obscurers);
  if (opaque && (!obscurer || opaque->distance_m < obscurer->distance_m))
    obscurer = DetectorSurfaceHit{opaque->surface_id, opaque->distance_m, opaque->point_m,
                                  opaque->unit_normal};
  if (obscurer) {
    const auto hit = *obscurer;
    if (!secondary || hit.distance_m < secondary->distance_m) {
      record.points_m[2] = hit.point_m;
      record.point_count = 3;
      record.path_length_m += hit.distance_m;
      record.status = PhotonStatus::blocked_obscurer;
      record.terminal_surface_id = hit.surface_id;
      record.interaction_surface_ids[1] = hit.surface_id;
      record.interaction_kinds[1] = OpticalInteractionKind::obscurer;
      record.interaction_normals[1] = hit.unit_normal;
      record.interaction_incoming_directions[1] = ray.direction;
      record.interaction_outgoing_directions[1] = ray.direction;
      record.interaction_throughput[1] = record.surviving_throughput;
      return record;
    }
  }
  if (!secondary) {
    record.status = PhotonStatus::missed_secondary;
    record.final_direction = ray.direction;
    return record;
  }
  record.points_m[2] = secondary->point_m;
  record.terminal_surface_id =
      *axisymmetric_segment_id(optical_model.secondary_segments, optical_model.secondary,
                               secondary->point_m, optical_model.secondary_surface_id);
  record.interaction_surface_ids[1] = record.terminal_surface_id;
  record.interaction_normals[1] = secondary->unit_normal;
  record.interaction_incoming_directions[1] = ray.direction;
  record.interaction_throughput[1] = record.surviving_throughput;
  record.point_count = 3;
  record.path_length_m += secondary->distance_m;
  record.incidence_secondary_deg =
      std::acos(std::clamp(std::abs(dot(ray.direction, secondary->unit_normal)), 0.0, 1.0)) *
      180.0 / std::numbers::pi;
  const auto after_secondary = optical_model.secondary_scatter ? [&]() -> std::optional<Ray> {
    const auto direction = reflect_with_scatter(ray.direction, secondary->unit_normal, {1, 0, 0},
                                                *optical_model.secondary_scatter, photon_id,
                                                record.terminal_surface_id, ray_tracing_seed);
    return direction ? std::optional{Ray{secondary->point_m, *direction}} : std::nullopt;
  }()
      : reflect(ray, *secondary);
  if (!after_secondary) {
    record.status = PhotonStatus::invalid_input;
    return record;
  }
  ray = *after_secondary;
  record.final_direction = ray.direction;
  record.interaction_outgoing_directions[1] = ray.direction;
  const auto secondary_response = optical_model.secondary_reflectivity
                                      ? optical_model.secondary_reflectivity->at_unchecked(
                                            wavelength_nm, record.incidence_secondary_deg)
                                      : std::optional<double>{1.0};
  if (!secondary_response) {
    record.status = PhotonStatus::invalid_input;
    return record;
  }
  record.surviving_throughput *= *secondary_response;
  if (optical_model.secondary_degradation) {
    const auto response =
        optical_model.secondary_degradation->at_point_unchecked(secondary->point_m);
    if (!response) {
      record.status = PhotonStatus::invalid_input;
      return record;
    }
    record.surviving_throughput *= *response;
  }
  record.interaction_throughput[1] = record.surviving_throughput;
  bool detector_failure = false;
  std::optional<AxisymmetricHit> detector;
  std::optional<Vec3> detector_local_position;
  std::uint32_t detector_surface_id = optical_model.detector_surface_id;
  if (optical_model.detector_planes) {
    std::optional<DetectorSurfaceHit> hit;
    if (optical_model.detector_assignment) {
      const auto focal = intersect_segmented_asphere(
          ray, optical_model.detector, {}, optical_model.detector_surface_id, &detector_failure);
      hit =
          focal ? optical_model.detector_assignment->intersect(ray, focal->point_m) : std::nullopt;
    } else
      hit = optical_model.detector_planes->intersect(ray);
    if (hit) {
      detector_local_position = hit->local_position_m;
      detector = AxisymmetricHit{hit->distance_m, hit->point_m, hit->unit_normal};
      detector_surface_id = hit->surface_id;
    }
  } else
    detector = intersect_segmented_asphere(ray, optical_model.detector, {},
                                           optical_model.detector_surface_id, &detector_failure);
  if (detector_failure) {
    record.status = PhotonStatus::intersection_failure;
    record.final_direction = ray.direction;
    return record;
  }
  const auto final_obscurer = intersect_opaque_surfaces(ray, optical_model.opaque_obscurers);
  if (final_obscurer && (!detector || final_obscurer->distance_m < detector->distance_m)) {
    const auto &hit = *final_obscurer;
    record.points_m[3] = hit.point_m;
    record.point_count = 4;
    record.path_length_m += hit.distance_m;
    record.status = PhotonStatus::blocked_obscurer;
    record.terminal_surface_id = hit.surface_id;
    record.interaction_surface_ids[2] = hit.surface_id;
    record.interaction_kinds[2] = OpticalInteractionKind::obscurer;
    record.interaction_normals[2] = hit.unit_normal;
    record.interaction_incoming_directions[2] = ray.direction;
    record.interaction_outgoing_directions[2] = ray.direction;
    record.interaction_throughput[2] = record.surviving_throughput;
    return record;
  }
  if (!detector) {
    record.status = PhotonStatus::missed_screen;
    record.final_direction = ray.direction;
    return record;
  }
  record.points_m[3] = detector->point_m;
  record.terminal_surface_id = detector_surface_id;
  record.interaction_surface_ids[2] = record.terminal_surface_id;
  record.interaction_kinds[2] = OpticalInteractionKind::detector;
  record.interaction_normals[2] = detector->unit_normal;
  record.interaction_incoming_directions[2] = ray.direction;
  record.interaction_outgoing_directions[2] = ray.direction;
  record.interaction_throughput[2] = record.surviving_throughput;
  record.point_count = 4;
  record.path_length_m += detector->distance_m;
  record.incidence_focal_deg =
      std::acos(std::clamp(std::abs(dot(ray.direction, detector->unit_normal)), 0.0, 1.0)) * 180.0 /
      std::numbers::pi;
  record.final_direction = ray.direction;
  if (optical_model.camera_degradation) {
    const auto point = optical_model.camera_degradation_in_detector_frame
                           ? detector_local_position
                           : std::optional{detector->point_m};
    const auto response =
        !point ? std::nullopt : optical_model.camera_degradation->at_point_unchecked(*point);
    if (!response) {
      record.status = PhotonStatus::invalid_input;
      return record;
    }
    record.surviving_throughput *= *response;
  }
  record.interaction_throughput[2] = record.surviving_throughput;
  if (optical_model.pixel_responses) {
    const auto response = optical_model.pixel_responses->at_unchecked(
        detector_surface_id, wavelength_nm, record.incidence_focal_deg,
        {detector->point_m, ray.direction});
    if (!response) {
      record.status = PhotonStatus::invalid_input;
      return record;
    }
    record.surviving_throughput *= *response;
    record.interaction_throughput[2] = record.surviving_throughput;
  }
  record.status = PhotonStatus::detected;
  if (optical_model.camera_response) {
    const auto response =
        optical_model.camera_response->at_unchecked(wavelength_nm, record.incidence_focal_deg);
    if (!response)
      record.status = PhotonStatus::invalid_input;
    else
      record.surviving_throughput *= *response;
    record.interaction_throughput[2] = record.surviving_throughput;
  }
  return record;
}

} // namespace obdeect
