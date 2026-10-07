#include "obdeect/optical_model_file.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>

std::string with_valid_hash(std::string document) {
  const auto root = obdeect::json::Parser{document}.parse();
  assert(root);
  const auto hash = obdeect::detail::sha256(obdeect::detail::canonical_json_without_hash(*root));
  const std::string marker = "\"optical_model_sha256\":\"" + std::string(64, 'a') + "\"";
  const auto position = document.find(marker);
  assert(position != std::string::npos);
  document.replace(position, marker.size(), "\"optical_model_sha256\":\"" + hash + "\"");
  return document;
}

template <class Reader>
void assert_unknown_fields_rejected(const char *path, Reader read,
                                    std::initializer_list<std::string_view> components) {
  std::ifstream input(path);
  const std::string original((std::istreambuf_iterator<char>(input)),
                             std::istreambuf_iterator<char>());
  for (const auto component : components) {
    auto root = obdeect::json::Parser{original}.parse();
    assert(root);
    auto *target = const_cast<obdeect::json::Value *>(root->find("trace_model"));
    if (!component.empty())
      target = const_cast<obdeect::json::Value *>(target->find(component));
    assert(target);
    if (target->kind == obdeect::json::Value::Kind::array)
      target = &target->array.front();
    target->object.emplace_back("unsupported_loss", obdeect::json::Value{});
    auto *hash = const_cast<obdeect::json::Value *>(root->find("optical_model_sha256"));
    hash->string = obdeect::detail::sha256(obdeect::detail::canonical_json_without_hash(*root));
    std::string document;
    obdeect::detail::append_canonical_json(*root, document);
    {
      std::ofstream output(path);
      output << document;
    }
    assert(!read(path));
  }
  std::ofstream output(path);
  output << original;
}

int main() {
  assert(obdeect::detail::sha256("") ==
         "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
  std::size_t response_case = 0;
  for (
      const auto document :
      {R"({"response":[{"wavelength_nm":300,"incidence_angle_deg":0,"response":0.4},{"wavelength_nm":300,"incidence_angle_deg":60,"response":0.6},{"wavelength_nm":500,"incidence_angle_deg":0,"response":0.6},{"wavelength_nm":500,"incidence_angle_deg":60,"response":0.8}]})",
       R"({"response":[{"wavelength_nm":300,"incidence_angle_deg":0,"response":0.4},{"wavelength_nm":300,"incidence_angle_deg":60,"response":0.6},{"wavelength_nm":500,"incidence_angle_deg":0,"response":0.6}]})",
       R"({"response":[{"wavelength_nm":300,"incidence_angle_deg":0,"response":0.4},{"wavelength_nm":500,"response":0.8}]})"}) {
    const auto root = obdeect::json::Parser{document}.parse();
    assert(root);
    const auto response = obdeect::detail::response_field(*root, "response");
    const bool complete = response_case++ == 0;
    assert(response.has_value() == complete);
    if (response)
      assert(std::abs(*response->at(400, 30) - 0.6) < 1.e-12);
  }
  const char *segmented_path = "obdeect-optical-model-test.json";
  {
    std::ofstream output(segmented_path);
    output << with_valid_hash(
        R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","provenance":{"model":"LSTN-design","model_version":"7.0.0"},"trace_model":{"kind":"segmented","primary_facets":[{"id":0,"shape":"circle","centre_m":[0,0,0],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":1,"focal_length_m":10}],"detector_surfaces":[{"id":1,"shape":"circle","centre_m":[0,0,10],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":1}],"cylinder_obscurers":[{"id":2,"first_endpoint_m":[0,0,2],"second_endpoint_m":[0,0,3],"diameter_m":0.1}],"primary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.9}]}})");
  }
  const auto segmented = obdeect::read_segmented_optical_model(segmented_path);
  assert_unknown_fields_rejected(segmented_path, obdeect::read_segmented_optical_model,
                                 {"", "primary_facets", "detector_surfaces", "cylinder_obscurers"});
  if (!segmented || segmented->primary_facets.size() != 1 ||
      segmented->detector_surfaces.size() != 1 || segmented->cylinder_obscurers.size() != 1 ||
      !segmented->primary_reflectivity ||
      std::abs(*segmented->primary_reflectivity->at(400) - 0.85) > 1.e-12) {
    std::remove(segmented_path);
    return 1;
  }
  {
    std::ifstream input(segmented_path);
    std::string tampered((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    const auto position = tampered.find("\"diameter_m\":1");
    assert(position != std::string::npos);
    tampered.replace(position, std::string("\"diameter_m\":1").size(), "\"diameter_m\":1.5");
    std::ofstream output(segmented_path);
    output << tampered;
  }
  if (obdeect::read_segmented_optical_model(segmented_path))
    return 1;
  std::remove(segmented_path);

  const char *optional_response_path = "obdeect-optical-model-optional-response-test.json";
  {
    std::ofstream output(optional_response_path);
    output << with_valid_hash(
        R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","provenance":{"model":"LSTN-design","model_version":"7.0.0"},"trace_model":{"kind":"segmented","primary_facets":[{"id":0,"shape":"circle","centre_m":[0,0,0],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":1,"focal_length_m":10}],"detector_surfaces":[{"id":1,"shape":"circle","centre_m":[0,0,10],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":1}],"cylinder_obscurers":[]}})");
  }
  const auto optional_response = obdeect::read_segmented_optical_model(optional_response_path);
  std::remove(optional_response_path);
  if (!optional_response || optional_response->primary_reflectivity)
    return 1;

  const char *python_artifact_path = "obdeect-python-hash-artifact-test.json";
  {
    std::ofstream output(python_artifact_path);
    output
        << R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"edb26c34fa09b514e7e751d40cc948dc19928256091a4916a9acea5dddb011e1","provenance":{"model":"python-model","model_version":"1.0"},"trace_model":{"cylinder_obscurers":[],"detector_surfaces":[{"centre_m":[0.0,0.0,10.0],"diameter_m":1.0,"id":1,"normal":[0.0,0.0,1.0],"shape":"circle","tangent":[1.0,0.0,0.0]}],"kind":"segmented","primary_facets":[{"centre_m":[0.0,0.0,0.0],"diameter_m":1.0,"focal_length_m":10.0,"id":0,"normal":[0.0,0.0,1.0],"shape":"circle","tangent":[1.0,0.0,0.0]}],"primary_reflectivity":[{"response":0.8,"wavelength_nm":300.0},{"response":0.9,"wavelength_nm":500.0}]}})";
  }
  const auto python_artifact = obdeect::read_segmented_optical_model(python_artifact_path);
  std::remove(python_artifact_path);
  if (!python_artifact)
    return 1;

  const char *axisymmetric_path = "obdeect-axisymmetric-optical-model-test.json";
  {
    std::ofstream output(axisymmetric_path);
    output << with_valid_hash(
        R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","provenance":{"model":"SSTS-design","model_version":"7.0.0"},"trace_model":{"kind":"axisymmetric","primary":{"vertex_z_m":0,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"secondary":{"vertex_z_m":2,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"detector":{"vertex_z_m":1,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"primary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.9}],"secondary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.9}]}})");
  }
  const auto axisymmetric = obdeect::read_axisymmetric_optical_model(axisymmetric_path);
  assert_unknown_fields_rejected(axisymmetric_path, obdeect::read_axisymmetric_optical_model,
                                 {"", "primary", "secondary", "detector"});
  if (!axisymmetric) {
    std::remove(axisymmetric_path);
    return 1;
  }
  assert(obdeect::trace_axisymmetric_optical_model(
             {{std::numeric_limits<double>::quiet_NaN(), 0, 10}, {0, 0, -1}}, 0, *axisymmetric)
             .status == obdeect::PhotonStatus::invalid_input);
  auto unit_response = *axisymmetric;
  unit_response.primary_reflectivity.reset();
  unit_response.secondary_reflectivity.reset();
  assert(obdeect::trace_axisymmetric_optical_model({{0, 0, 10}, {0, 0, -1}}, 0, unit_response, 0)
             .status == obdeect::PhotonStatus::invalid_input);
  auto masked = *axisymmetric;
  masked.primary_surface_id = 10;
  masked.secondary_surface_id = 11;
  masked.detector_surface_id = 12;
  obdeect::AxisymmetricSegment sector{};
  sector.id = 21;
  sector.shape = obdeect::AxisymmetricSegmentShape::annular_sector;
  sector.inner_radius_m = 0.5;
  sector.outer_radius_m = 2;
  sector.span_rad = std::numbers::pi / 2;
  masked.primary_segments = {sector};
  auto overlapping_sector = sector;
  overlapping_sector.id = 20;
  const std::array sectors{sector, overlapping_sector};
  const std::array reversed_sectors{overlapping_sector, sector};
  assert(obdeect::axisymmetric_segment_id(sectors, masked.primary, {1, 0, 0}, 10) == 20);
  assert(obdeect::axisymmetric_segment_id(reversed_sectors, masked.primary, {1, 0, 0}, 10) == 20);
  auto mask_hit = obdeect::trace_axisymmetric_optical_model({{1, 0, 10}, {0, 0, -1}}, 2, masked);
  assert(mask_hit.status == obdeect::PhotonStatus::detected &&
         mask_hit.interaction_surface_ids[0] == 21 && mask_hit.interaction_surface_ids[1] == 11 &&
         mask_hit.terminal_surface_id == 12);
  auto mask_miss = obdeect::trace_axisymmetric_optical_model({{-1, 0, 10}, {0, 0, -1}}, 2, masked);
  assert(mask_miss.status == obdeect::PhotonStatus::missed_primary && mask_miss.point_count == 1);
  masked.primary_segments[0].gap_m = 0.1;
  assert(!obdeect::contains_axisymmetric_segment(
      masked.primary_segments[0], masked.primary,
      {std::cos(std::numbers::pi / 2 - 0.05), std::sin(std::numbers::pi / 2 - 0.05), 0}));
  assert(obdeect::contains_axisymmetric_segment(
      masked.primary_segments[0], masked.primary,
      {std::cos(std::numbers::pi / 2 - 0.15), std::sin(std::numbers::pi / 2 - 0.15), 0}));
  auto reversed_gap = masked.primary_segments[0];
  reversed_gap.gap_at_start = true;
  assert(!obdeect::contains_axisymmetric_segment(reversed_gap, masked.primary,
                                                 {std::cos(0.05), std::sin(0.05), 0}));
  assert(obdeect::contains_axisymmetric_segment(reversed_gap, masked.primary,
                                                {std::cos(0.15), std::sin(0.15), 0}));
  sector.id = 22;
  sector.start_rad = std::numbers::pi / 2;
  masked.secondary_segments = {sector};
  auto secondary_miss =
      obdeect::trace_axisymmetric_optical_model({{1, 0, 10}, {0, 0, -1}}, 2, masked);
  assert(secondary_miss.status == obdeect::PhotonStatus::missed_secondary &&
         secondary_miss.terminal_surface_id == 21);
  masked.secondary_segments.clear();
  masked.block_incoming_secondary = true;
  auto shadow = obdeect::trace_axisymmetric_optical_model({{1, 0, 10}, {0, 0, -1}}, 2, masked);
  assert(shadow.status == obdeect::PhotonStatus::blocked_obscurer && shadow.point_count == 2 &&
         shadow.terminal_surface_id == 11 && std::abs(shadow.path_length_m - 8) < 1.e-12);
  obdeect::AxisymmetricSegment hexagon{};
  hexagon.id = 30;
  hexagon.shape = obdeect::AxisymmetricSegmentShape::hexagon;
  hexagon.diameter_m = 2;
  assert(obdeect::contains_axisymmetric_segment(hexagon, masked.primary, {0, 1.1, 0}));
  assert(!obdeect::contains_axisymmetric_segment(hexagon, masked.primary, {1.01, 0, 0}));
  hexagon.rotation_rad = std::numbers::pi / 2;
  assert(!obdeect::contains_axisymmetric_segment(hexagon, masked.primary, {0, 1.1, 0}));
  const auto malformed_masks =
      obdeect::json::Parser{R"({"primary_segments":[{"id":2,"shape":"unknown"}]})"}.parse();
  assert(malformed_masks && !obdeect::detail::segment_fields(*malformed_masks, "primary_segments"));
  auto numerical = *axisymmetric;
  numerical.primary.surface.coefficient_m[12] = std::numeric_limits<double>::max();
  auto solver_failure =
      obdeect::trace_axisymmetric_optical_model({{1.5, 0, 10}, {0, 0, -1}}, 3, numerical);
  assert(solver_failure.status == obdeect::PhotonStatus::intersection_failure);
  const auto record =
      obdeect::trace_axisymmetric_optical_model({{0, 0, 10}, {0, 0, -1}}, 1, *axisymmetric);
  // Explicit IDs must distinguish continuous surfaces and finite segments.
  for (
      const std::string extra :
      {R"(,"primary_surface_id":10,"secondary_surface_id":10)",
       R"(,"detector_surface_id":4294967295)",
       R"(,"detector_surface_id":21,"primary_segments":[{"id":21,"shape":"hexagon","centre_xy_m":[0,0],"diameter_m":1,"rotation_deg":0}])"}) {
    std::ifstream input(axisymmetric_path);
    std::string document((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    auto original = obdeect::json::Parser{document}.parse();
    assert(original);
    const auto hash = *obdeect::detail::string_field(*original, "optical_model_sha256");
    document.replace(document.find(hash), hash.size(), std::string(64, 'a'));
    document.insert(document.size() - 2, extra);
    {
      std::ofstream output(axisymmetric_path);
      output << with_valid_hash(document);
    }
    assert(!obdeect::read_axisymmetric_optical_model(axisymmetric_path));
    // Restore the original for the next independent invalid-ID case.
    std::string restored;
    obdeect::detail::append_canonical_json(*original, restored);
    std::ofstream restore(axisymmetric_path);
    restore << restored;
  }
  std::remove(axisymmetric_path);
  if (record.status != obdeect::PhotonStatus::detected || record.point_count != 4)
    return 1;
}
