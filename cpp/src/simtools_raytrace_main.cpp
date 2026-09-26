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
  std::cout << "Usage: obdeect-simtools-raytrace --telescope LST|MST|SST|SCT [options]\n"
            << "  --source star|illuminator|laser  (default: star)\n"
            << "  --photons N --output FILE --field-x-deg D --field-y-deg D\n"
            << "  --distance-m D --wavelength-nm D --divergence-deg D\n"
            << "  --source-x-m D --source-y-m D --source-z-m D\n"
            << "  --direction-x D --direction-y D --direction-z D  (laser axis)\n"
            << "Reference optical CLI: imported production scenes are not used yet.\n";
}

bool value(int& index, int argc, char** argv, const char* option, std::string& result) {
  if (std::string{argv[index]} != option || index + 1 >= argc) return false;
  result = argv[++index];
  return true;
}
}  // namespace

int main(int argc, char** argv) {
  using obdeect::parse_finite_double;
  using obdeect::parse_positive_size;
  std::string telescope;
  std::string source{"star"};
  std::string output_path{"obdeect_paths.csv"};
  std::size_t photons_count = 10000;
  double field_x_deg = 0.0, field_y_deg = 0.0, distance_m = 10000.0;
  double wavelength_nm = 400.0, divergence_deg = 0.0;
  double source_x_m = 0.0, source_y_m = 0.0, source_z_m = 50.0;
  double direction_x = 0.0, direction_y = 0.0, direction_z = -1.0;
  for (int index = 1; index < argc; ++index) {
    std::string value_string;
    if (value(index, argc, argv, "--telescope", telescope) ||
        value(index, argc, argv, "--source", source) ||
        value(index, argc, argv, "--output", output_path)) {
      continue;
    }
    if (std::string{argv[index]} == "--photons" && index + 1 < argc &&
        parse_positive_size(argv[++index], photons_count)) continue;
    auto number = [&](const char* option, double& target) {
      return std::string{argv[index]} == option && index + 1 < argc &&
             parse_finite_double(argv[++index], target);
    };
    if (number("--field-x-deg", field_x_deg) || number("--field-y-deg", field_y_deg) ||
        number("--distance-m", distance_m) || number("--wavelength-nm", wavelength_nm) ||
        number("--divergence-deg", divergence_deg) || number("--source-x-m", source_x_m) ||
        number("--source-y-m", source_y_m) || number("--source-z-m", source_z_m) ||
        number("--direction-x", direction_x) || number("--direction-y", direction_y) ||
        number("--direction-z", direction_z)) continue;
    usage();
    return 2;
  }
  const auto model = obdeect::ctao_reference_model(telescope);
  if (!model || (source != "star" && source != "illuminator" && source != "laser") ||
      wavelength_nm <= 0.0 || distance_m <= 0.0 || divergence_deg < 0.0 || divergence_deg >= 90.0) {
    usage();
    return 2;
  }
  const double radians_per_degree = std::numbers::pi / 180.0;
  std::vector<obdeect::OpticalPhoton> input;
  if (source == "star") {
    input = obdeect::star_photons(photons_count, model->primary_outer_radius_m,
                                  {field_x_deg * radians_per_degree, field_y_deg * radians_per_degree,
                                   distance_m, wavelength_nm});
  } else if (source == "illuminator") {
    input = obdeect::illuminator_photons(
        photons_count, model->primary_outer_radius_m,
        {{source_x_m, source_y_m, source_z_m}, wavelength_nm, 1.0});
  } else {
    input = obdeect::laser_photons(
        photons_count, model->primary_outer_radius_m,
        {{direction_x, direction_y, direction_z}, {source_x_m, source_y_m, source_z_m},
         divergence_deg * radians_per_degree, wavelength_nm});
  }
  if (input.size() != photons_count) {
    std::cerr << "source configuration produced no valid photons\n";
    return 1;
  }
  std::ofstream output{output_path};
  if (!output) {
    std::cerr << "cannot write " << output_path << '\n';
    return 1;
  }
  output << std::setprecision(17)
         << "contract_version,photon_id,source_kind,wavelength_nm,emission_time_ns,source_weight,throughput,status,point_count,path_length_m,"
            "incidence_primary_deg,incidence_secondary_deg,incidence_focal_deg,";
  for (int point = 0; point < 4; ++point)
    output << "x" << point << "_m,y" << point << "_m,z" << point << "_m" << (point == 3 ? '\n' : ',');
  std::size_t detected = 0;
  for (const auto& photon : input) {
    const auto path = obdeect::trace_ctao_reference(photon.ray, photon.photon_id, *model);
    detected += path.status == obdeect::PhotonStatus::detected;
    const double throughput = path.status == obdeect::PhotonStatus::detected ? 1.0 : 0.0;
    output << "obdeect-arrival-v1," << path.photon_id << ',' << source << ',' << photon.wavelength_nm << ',' << photon.time_ns << ',' << photon.weight << ','
           << throughput << ',' << obdeect::to_string(path.status) << ',' << static_cast<int>(path.point_count) << ','
           << path.path_length_m << ',' << path.incidence_primary_deg << ',' << path.incidence_secondary_deg << ','
           << path.incidence_focal_deg;
    for (const auto& point : path.points_m) output << ',' << point.x << ',' << point.y << ',' << point.z;
    output << '\n';
  }
  std::cout << "reference " << model->identifier << ": " << photons_count << " " << source
            << " photons, detected " << detected << "\n";
  return 0;
}
