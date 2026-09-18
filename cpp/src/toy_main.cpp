#include "obdeect/toy_mst.hpp"

#include <charconv>
#include <fstream>
#include <iostream>
#include <string>

namespace {

bool parse_size(const char* value, std::size_t& output) {
  const std::string text{value};
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), output);
  return error == std::errc{} && end == text.data() + text.size() && output > 0;
}

void usage() {
  std::cout << "Usage: obdeect_toy [--photons N] [--output paths.csv] [--no-structure]\n"
            << "Trace 400-nm parallel Cherenkov photons through a simple MST-inspired\n"
            << "spherical mirror, camera shadow and four mast supports.\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::size_t photon_count = 10000;
  std::string output_path{"toy_mst_paths.csv"};
  obdeect::ToyMstConfig config{};
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
    usage();
    return 2;
  }

  std::ofstream output{output_path};
  if (!output) {
    std::cerr << "Cannot write " << output_path << '\n';
    return 1;
  }
  output << "photon_id,wavelength_nm,status,point_count,path_length_m,x0_m,y0_m,z0_m,x1_m,y1_m,z1_m,x2_m,y2_m,z2_m\n";

  std::size_t detected = 0;
  std::size_t camera_blocked = 0;
  std::size_t mast_blocked = 0;
  const auto rays = obdeect::parallel_blue_cherenkov_rays(photon_count, config);
  for (std::size_t index = 0; index < rays.size(); ++index) {
    const auto record = obdeect::trace_toy_mst(rays[index], index, config);
    detected += record.status == obdeect::PhotonStatus::detected;
    camera_blocked += record.status == obdeect::PhotonStatus::blocked_camera;
    mast_blocked += record.status == obdeect::PhotonStatus::blocked_mast;
    output << record.photon_id << ',' << record.wavelength_nm << ',' << obdeect::to_string(record.status) << ','
           << static_cast<int>(record.point_count) << ',' << record.path_length_m;
    for (const auto& point : record.points_m) output << ',' << point.x << ',' << point.y << ',' << point.z;
    output << '\n';
  }

  std::cout << "MST baseline: " << photon_count << " artificial 400-nm photons\n"
            << "  detected: " << detected << " (" << 100.0 * detected / photon_count << "%)\n"
            << "  camera shadow: " << camera_blocked << "\n"
            << "  mast shadow: " << mast_blocked << "\n"
            << "  paths: " << output_path << '\n';
}
