#include "obdeect/optical_model_file.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>

int main() {
  const char* segmented_path = "obdeect-optical-model-test.json";
  {
    std::ofstream output(segmented_path);
    output << R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","provenance":{"model":"LSTN-design","model_version":"7.0.0"},"trace_model":{"kind":"segmented","primary_facets":[{"id":0,"shape":"circle","centre_m":[0,0,0],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":1,"focal_length_m":10}],"detector_surfaces":[{"id":1,"shape":"circle","centre_m":[0,0,10],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":1}],"cylinder_obscurers":[{"id":2,"first_endpoint_m":[0,0,2],"second_endpoint_m":[0,0,3],"diameter_m":0.1}],"primary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.9}]}})";
  }
  const auto segmented = obdeect::read_segmented_optical_model(segmented_path);
  if (!segmented || segmented->primary_facets.size() != 1 || segmented->detector_surfaces.size() != 1 ||
      segmented->cylinder_obscurers.size() != 1 || !segmented->primary_reflectivity ||
      std::abs(*segmented->primary_reflectivity->at(400) - 0.85) > 1.e-12) {
    std::remove(segmented_path);
    return 1;
  }
  std::remove(segmented_path);

  const char* axisymmetric_path = "obdeect-axisymmetric-optical-model-test.json";
  {
    std::ofstream output(axisymmetric_path);
    output << R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","provenance":{"model":"SSTS-design","model_version":"7.0.0"},"trace_model":{"kind":"axisymmetric","primary":{"vertex_z_m":0,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"secondary":{"vertex_z_m":2,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"detector":{"vertex_z_m":1,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"primary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.9}],"secondary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.9}]}})";
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
