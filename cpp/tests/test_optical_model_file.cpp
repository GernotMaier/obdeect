#include "obdeect/optical_model_file.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>

int main() {
  const char* path = "obdeect-native-optical_model-test.csv";
  {
    std::ofstream output(path);
    output << "obdeect-optical-model-v1\n"
              "provenance,LSTN-design,7.0.0,aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n"
              "surface_id,role,shape,cx_m,cy_m,cz_m,nx,ny,nz,tx,ty,tz,diameter_m,focal_length_m\n"
              "0,mirror,circle,0,0,0,0,0,1,1,0,0,1,10\n"
              "1,detector,circle,0,0,10,0,0,1,1,0,0,1,0\n"
              "obscurer_cylinder,2,0,0,2,0,0,3,0.1\n"
              "primary_reflectivity,300,0.8\n"
              "primary_reflectivity,500,0.9\n";
  }
  const auto optical_model = obdeect::read_native_optical_model(path);
  if (!optical_model || optical_model->primary_facets.size() != 1 || optical_model->detector_surfaces.size() != 1 ||
      optical_model->cylinder_obscurers.size() != 1 || !optical_model->primary_reflectivity ||
      std::abs(*optical_model->primary_reflectivity->at(400) - 0.85) > 1.e-12) {
    std::remove(path);
    return 1;
  }
  std::remove(path);

  const char* axisymmetric_path = "obdeect-axisymmetric-optical-model-test.csv";
  {
    std::ofstream output(axisymmetric_path);
    output << "obdeect-axisymmetric-optical-model-v1\n"
              "provenance,SSTS-design,7.0.0,aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n"
              "role,vertex_z_m,inner_radius_m,outer_radius_m,radial_scale_m,c0_m,c1_m,c2_m,c3_m,c4_m,c5_m,c6_m,c7_m,c8_m,c9_m,c10_m,c11_m,c12_m\n"
              "primary,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0\n"
              "secondary,2,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0\n"
              "detector,1,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0\n";
  }
  const auto axisymmetric_optical_model = obdeect::read_native_axisymmetric_optical_model(axisymmetric_path);
  if (!axisymmetric_optical_model) {
    std::remove(axisymmetric_path);
    return 1;
  }
  const auto record = obdeect::trace_axisymmetric_optical_model({{0, 0, 10}, {0, 0, -1}}, 1, *axisymmetric_optical_model);
  std::remove(axisymmetric_path);
  if (record.status != obdeect::PhotonStatus::detected || record.point_count != 4) return 1;
}
