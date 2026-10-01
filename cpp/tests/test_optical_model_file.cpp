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

int main() {
  assert(obdeect::detail::sha256("") ==
         "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
  const char* segmented_path = "obdeect-optical-model-test.json";
  {
    std::ofstream output(segmented_path);
    output << with_valid_hash(R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","provenance":{"model":"LSTN-design","model_version":"7.0.0"},"trace_model":{"kind":"segmented","primary_facets":[{"id":0,"shape":"circle","centre_m":[0,0,0],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":1,"focal_length_m":10}],"detector_surfaces":[{"id":1,"shape":"circle","centre_m":[0,0,10],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":1}],"cylinder_obscurers":[{"id":2,"first_endpoint_m":[0,0,2],"second_endpoint_m":[0,0,3],"diameter_m":0.1}],"primary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.9}]}})");
  }
  const auto segmented = obdeect::read_segmented_optical_model(segmented_path);
  if (!segmented || segmented->primary_facets.size() != 1 || segmented->detector_surfaces.size() != 1 ||
      segmented->cylinder_obscurers.size() != 1 || !segmented->primary_reflectivity ||
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
  if (obdeect::read_segmented_optical_model(segmented_path)) return 1;
  std::remove(segmented_path);

  const char* optional_response_path = "obdeect-optical-model-optional-response-test.json";
  {
    std::ofstream output(optional_response_path);
    output << with_valid_hash(R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","provenance":{"model":"LSTN-design","model_version":"7.0.0"},"trace_model":{"kind":"segmented","primary_facets":[{"id":0,"shape":"circle","centre_m":[0,0,0],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":1,"focal_length_m":10}],"detector_surfaces":[{"id":1,"shape":"circle","centre_m":[0,0,10],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":1}],"cylinder_obscurers":[]}})");
  }
  const auto optional_response = obdeect::read_segmented_optical_model(optional_response_path);
  std::remove(optional_response_path);
  if (!optional_response || optional_response->primary_reflectivity) return 1;

  const char* python_artifact_path = "obdeect-python-hash-artifact-test.json";
  {
    std::ofstream output(python_artifact_path);
    output << R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"edb26c34fa09b514e7e751d40cc948dc19928256091a4916a9acea5dddb011e1","provenance":{"model":"python-model","model_version":"1.0"},"trace_model":{"cylinder_obscurers":[],"detector_surfaces":[{"centre_m":[0.0,0.0,10.0],"diameter_m":1.0,"id":1,"normal":[0.0,0.0,1.0],"shape":"circle","tangent":[1.0,0.0,0.0]}],"kind":"segmented","primary_facets":[{"centre_m":[0.0,0.0,0.0],"diameter_m":1.0,"focal_length_m":10.0,"id":0,"normal":[0.0,0.0,1.0],"shape":"circle","tangent":[1.0,0.0,0.0]}],"primary_reflectivity":[{"response":0.8,"wavelength_nm":300.0},{"response":0.9,"wavelength_nm":500.0}]}})";
  }
  const auto python_artifact = obdeect::read_segmented_optical_model(python_artifact_path);
  std::remove(python_artifact_path);
  if (!python_artifact) return 1;

  const char* axisymmetric_path = "obdeect-axisymmetric-optical-model-test.json";
  {
    std::ofstream output(axisymmetric_path);
    output << with_valid_hash(R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","provenance":{"model":"SSTS-design","model_version":"7.0.0"},"trace_model":{"kind":"axisymmetric","primary":{"vertex_z_m":0,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"secondary":{"vertex_z_m":2,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"detector":{"vertex_z_m":1,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"primary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.9}],"secondary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.9}]}})");
  }
  const auto axisymmetric = obdeect::read_axisymmetric_optical_model(axisymmetric_path);
  if (!axisymmetric) {
    std::remove(axisymmetric_path);
    return 1;
  }
  const auto record = obdeect::trace_axisymmetric_optical_model({{0, 0, 10}, {0, 0, -1}}, 1, *axisymmetric);
  std::remove(axisymmetric_path);
  if (record.status != obdeect::PhotonStatus::detected || record.point_count != 4) return 1;
}
