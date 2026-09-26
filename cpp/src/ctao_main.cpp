#include "obdeect/cli_parse.hpp"
#include "obdeect/ctao_models.hpp"
#include "obdeect/ctao_trace.hpp"
#include "obdeect/sources.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>

namespace {
void usage() {
  std::cout << "Usage: obdeect_ctao --telescope LST|MST|SST|SCT [--photons N] [--output paths.csv]\n"
            << "                     [--field-x-deg D] [--field-y-deg D] [--wavelength-nm D]\n"
            << "Trace monochromatic parallel reference photons through an axisymmetric CTAO optical prescription.\n";
}
}  // namespace

int main(int argc, char** argv) {
  using obdeect::parse_finite_double;
  using obdeect::parse_positive_size;
  std::size_t photon_count = 10000;
  std::string telescope_name;
  std::string output_path{"ctao_paths.csv"};
  double field_x_deg = 0.0;
  double field_y_deg = 0.0;
  double wavelength_nm = 400.0;
  for (int index = 1; index < argc; ++index) {
    const std::string argument{argv[index]};
    if (argument == "--telescope" && index + 1 < argc) { telescope_name = argv[++index]; continue; }
    if (argument == "--photons" && index + 1 < argc && parse_positive_size(argv[++index], photon_count)) continue;
    if (argument == "--output" && index + 1 < argc) { output_path = argv[++index]; continue; }
    if (argument == "--field-x-deg" && index + 1 < argc && parse_finite_double(argv[++index], field_x_deg)) continue;
    if (argument == "--field-y-deg" && index + 1 < argc && parse_finite_double(argv[++index], field_y_deg)) continue;
    if (argument == "--wavelength-nm" && index + 1 < argc && parse_finite_double(argv[++index], wavelength_nm) &&
        wavelength_nm > 0.0) {
      continue;
    }
    usage();
    return 2;
  }
  const auto model = obdeect::ctao_reference_model(telescope_name);
  if (!model) {
    std::cerr << "--telescope must be one of LST, MST, SST, SCT\n";
    return 2;
  }
  const double radians_per_degree = std::numbers::pi / 180.0;
  const auto photons = obdeect::star_photons(
      photon_count, model->primary_outer_radius_m,
      {field_x_deg * radians_per_degree, field_y_deg * radians_per_degree, 60.0, wavelength_nm});
  std::ofstream output{output_path};
  if (!output) {
    std::cerr << "Cannot write " << output_path << '\n';
    return 1;
  }
  output << std::setprecision(17)
         << "photon_id,wavelength_nm,emission_time_ns,source_weight,throughput,status,point_count,path_length_m,"
            "x0_m,y0_m,z0_m,x1_m,y1_m,z1_m,"
            "x2_m,y2_m,z2_m,x3_m,y3_m,z3_m\n";
  std::size_t detected = 0;
  for (const auto& photon : photons) {
    auto path = obdeect::trace_ctao_reference(photon.ray, photon.photon_id, *model);
    path.wavelength_nm = photon.wavelength_nm;
    detected += path.status == obdeect::PhotonStatus::detected;
    const double throughput = path.status == obdeect::PhotonStatus::detected ? 1.0 : 0.0;
    output << path.photon_id << ',' << path.wavelength_nm << ',' << photon.time_ns << ',' << photon.weight << ','
           << throughput << ',' << obdeect::to_string(path.status) << ',' << static_cast<int>(path.point_count)
           << ',' << path.path_length_m;
    for (const auto& point : path.points_m) output << ',' << point.x << ',' << point.y << ',' << point.z;
    output << '\n';
  }
  std::cout << model->identifier << " reference optical trace: " << photon_count << ' ' << wavelength_nm
            << "-nm photons\n"
            << "  detected: " << detected << " (" << 100.0 * detected / photon_count << "%)\n"
            << "  paths: " << output_path << '\n';
  return 0;
}
