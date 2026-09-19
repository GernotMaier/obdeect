#include "obdeect/sources.hpp"
#include "obdeect/toy_mst.hpp"

#include <charconv>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>

namespace {

bool parse_size(const char* value, std::size_t& output) {
  const std::string text{value};
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), output);
  return error == std::errc{} && end == text.data() + text.size() && output > 0;
}

bool parse_double(const char* value, double& output) {
  const std::string text{value};
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), output);
  return error == std::errc{} && end == text.data() + text.size() && std::isfinite(output);
}

bool parse_source(const std::string& value, obdeect::ArtificialSourceKind& output) {
  if (value == "star") {
    output = obdeect::ArtificialSourceKind::star;
    return true;
  }
  if (value == "illuminator") {
    output = obdeect::ArtificialSourceKind::illuminator;
    return true;
  }
  if (value == "laser") {
    output = obdeect::ArtificialSourceKind::laser;
    return true;
  }
  return false;
}

void usage() {
  std::cout << "Usage: obdeect_toy [--photons N] [--output paths.csv] [--no-structure]\n"
            << "                    [--source star|illuminator|laser] [--field-x-deg D]\n"
            << "                    [--field-y-deg D] [--distance-m D] [--divergence-deg D]\n"
            << "Trace 400-nm artificial photons through a simple MST-inspired spherical\n"
            << "mirror, camera shadow and four mast supports.\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::size_t photon_count = 10000;
  std::string output_path{"toy_mst_paths.csv"};
  obdeect::ToyMstConfig config{};
  obdeect::ArtificialSourceKind source_kind = obdeect::ArtificialSourceKind::star;
  double field_x_deg = 0.0;
  double field_y_deg = 0.0;
  double distance_m = 50.0;
  double divergence_deg = 0.0;
  for (int index = 1; index < argc; ++index) {
    const std::string argument{argv[index]};
    if (argument == "--photons" && index + 1 < argc && parse_size(argv[++index], photon_count)) continue;
    if (argument == "--output" && index + 1 < argc) {
      output_path = argv[++index];
      continue;
    }
    if (argument == "--no-structure") {
      config.include_structure = false;
      continue;
    }
    if (argument == "--source" && index + 1 < argc && parse_source(argv[++index], source_kind)) continue;
    if (argument == "--field-x-deg" && index + 1 < argc && parse_double(argv[++index], field_x_deg)) continue;
    if (argument == "--field-y-deg" && index + 1 < argc && parse_double(argv[++index], field_y_deg)) continue;
    if (argument == "--distance-m" && index + 1 < argc && parse_double(argv[++index], distance_m) &&
        distance_m > 0.0) {
      continue;
    }
    if (argument == "--divergence-deg" && index + 1 < argc && parse_double(argv[++index], divergence_deg) &&
        divergence_deg >= 0.0) {
      continue;
    }
    usage();
    return 2;
  }

  std::ofstream output{output_path};
  if (!output) {
    std::cerr << "Cannot write " << output_path << '\n';
    return 1;
  }
  output << std::setprecision(17);
  output << "photon_id,wavelength_nm,emission_time_ns,source_weight,throughput,status,point_count,path_length_m,"
            "x0_m,y0_m,z0_m,x1_m,y1_m,z1_m,x2_m,y2_m,z2_m,x3_m,y3_m,z3_m\n";

  std::size_t detected = 0;
  std::size_t camera_blocked = 0;
  std::size_t mast_blocked = 0;
  const double radians_per_degree = std::numbers::pi / 180.0;
  std::vector<obdeect::OpticalPhoton> photons;
  if (source_kind == obdeect::ArtificialSourceKind::star) {
    photons = obdeect::star_photons(photon_count, config.mirror_aperture_radius_m,
                                    {field_x_deg * radians_per_degree, field_y_deg * radians_per_degree,
                                     distance_m, 400.0});
  } else if (source_kind == obdeect::ArtificialSourceKind::illuminator) {
    photons = obdeect::illuminator_photons(photon_count, config.mirror_aperture_radius_m,
                                            {{0.0, 0.0, distance_m}, 400.0, 1.0});
  } else {
    photons = obdeect::laser_photons(photon_count, config.mirror_aperture_radius_m,
                                      {{0.0, 0.0, -1.0}, distance_m,
                                       divergence_deg * radians_per_degree, 400.0});
  }
  if (photons.size() != photon_count) {
    std::cerr << "Cannot generate the requested source photons\n";
    return 1;
  }
  for (const auto& photon : photons) {
    auto record = obdeect::trace_toy_mst(photon.ray, photon.photon_id, config);
    record.wavelength_nm = photon.wavelength_nm;
    detected += record.status == obdeect::PhotonStatus::detected;
    camera_blocked += record.status == obdeect::PhotonStatus::blocked_camera;
    mast_blocked += record.status == obdeect::PhotonStatus::blocked_mast;
    const double throughput = record.status == obdeect::PhotonStatus::detected ? 1.0 : 0.0;
    output << record.photon_id << ',' << record.wavelength_nm << ',' << photon.time_ns << ',' << photon.weight << ','
           << throughput << ',' << obdeect::to_string(record.status) << ',' << static_cast<int>(record.point_count)
           << ',' << record.path_length_m;
    for (const auto& point : record.points_m) output << ',' << point.x << ',' << point.y << ',' << point.z;
    output << '\n';
  }

  std::cout << "MST baseline: " << photon_count << " artificial 400-nm " << obdeect::to_string(source_kind)
            << " photons\n"
            << "  detected: " << detected << " (" << 100.0 * detected / photon_count << "%)\n"
            << "  camera shadow: " << camera_blocked << "\n"
            << "  mast shadow: " << mast_blocked << "\n"
            << "  paths: " << output_path << '\n';
}
