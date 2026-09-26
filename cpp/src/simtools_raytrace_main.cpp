#include "obdeect/cli_parse.hpp"
#include "obdeect/ctao_models.hpp"
#include "obdeect/ctao_trace.hpp"
#include "obdeect/interactions.hpp"
#include "obdeect/scene_file.hpp"
#include "obdeect/segmented_scene.hpp"
#include "obdeect/sources.hpp"

#include <fstream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>

namespace {
void usage() {
  std::cout << "Usage: obdeect-simtools-raytrace [--telescope LST|MST|SST|SCT] [options]\n"
            << "  --scene-file FILE  model-derived obdeect-scene-v1 surface table\n"
            << "  --source star|illuminator|laser  (default: star)\n"
            << "  --photons N --output FILE --field-x-deg D --field-y-deg D\n"
            << "  --distance-m D --wavelength-nm D --divergence-deg D\n"
            << "  --source-x-m D --source-y-m D --source-z-m D\n"
            << "  --direction-x D --direction-y D --direction-z D  (laser axis)\n"
            << "  Without --scene-file, use the built-in CTAO reference prescription.\n";
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
  std::string scene_file;
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
        value(index, argc, argv, "--scene-file", scene_file) ||
        value(index, argc, argv, "--scene", scene_file) ||
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
  const auto model = scene_file.empty() ? obdeect::ctao_reference_model(telescope) : std::nullopt;
  const auto imported_scene = scene_file.empty() ? std::optional<obdeect::CompiledSegmentedScene>{}
                                                 : obdeect::read_native_scene(scene_file);
  if ((!model && !imported_scene) || (source != "star" && source != "illuminator" && source != "laser") ||
      wavelength_nm <= 0.0 || distance_m <= 0.0 || divergence_deg < 0.0 || divergence_deg >= 90.0) {
    usage();
    return 2;
  }
  const double radians_per_degree = std::numbers::pi / 180.0;
  const double pupil_radius = imported_scene
                                  ? [&] {
                                      double radius = 0.0;
                                      for (const auto& facet : imported_scene->primary_facets)
                                        radius = std::max(radius, std::hypot(facet.centre_m.x, facet.centre_m.y) +
                                                                  facet.diameter_m * 0.5);
                                      return radius;
                                    }()
                                  : model->primary_outer_radius_m;
  if (!std::isfinite(pupil_radius) || pupil_radius <= 0.0) return 1;
  std::vector<obdeect::OpticalPhoton> input;
  if (source == "star") {
    input = obdeect::star_photons(photons_count, pupil_radius,
                                  {field_x_deg * radians_per_degree, field_y_deg * radians_per_degree,
                                   distance_m, wavelength_nm});
  } else if (source == "illuminator") {
    input = obdeect::illuminator_photons(
        photons_count, pupil_radius,
        {{source_x_m, source_y_m, source_z_m}, wavelength_nm, 1.0});
  } else {
    input = obdeect::laser_photons(
        photons_count, pupil_radius,
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
    obdeect::PathRecord path{};
    path.photon_id = photon.photon_id;
    path.wavelength_nm = photon.wavelength_nm;
    path.points_m[0] = photon.ray.position_m;
    path.point_count = 1;
    if (imported_scene) {
      const auto direction = obdeect::normalised_checked(photon.ray.direction);
      if (!direction) {
        path.status = obdeect::PhotonStatus::invalid_input;
      } else {
        const obdeect::Ray ray{photon.ray.position_m, *direction};
        const auto primary_hit = obdeect::intersect_segmented_primary_unchecked(ray, *imported_scene);
        if (!primary_hit) {
          path.status = obdeect::PhotonStatus::missed_primary;
          path.final_direction = *direction;
        } else {
          path.points_m[1] = primary_hit->point_m;
          path.point_count = 2;
          path.path_length_m = primary_hit->distance_m;
          path.incidence_primary_deg = std::acos(std::clamp(std::abs(obdeect::dot(*direction, primary_hit->unit_normal)), 0.0, 1.0)) *
                                       180.0 / std::numbers::pi;
          const auto reflected = obdeect::reflect_specular(*direction, primary_hit->unit_normal);
          if (!reflected) {
            path.status = obdeect::PhotonStatus::invalid_input;
          } else {
            const auto detector_hit = obdeect::intersect_detector_surfaces_unchecked(
                obdeect::Ray{primary_hit->point_m, *reflected}, *imported_scene);
            if (!detector_hit) {
              path.status = obdeect::PhotonStatus::missed_screen;
              path.final_direction = *reflected;
            } else {
              path.points_m[2] = detector_hit->point_m;
              path.point_count = 3;
              path.path_length_m += detector_hit->distance_m;
              path.final_direction = *reflected;
              path.incidence_focal_deg = std::acos(std::clamp(std::abs(obdeect::dot(*reflected, detector_hit->unit_normal)), 0.0, 1.0)) *
                                         180.0 / std::numbers::pi;
              path.status = obdeect::PhotonStatus::detected;
            }
          }
        }
      }
    } else {
      path = obdeect::trace_ctao_reference(photon.ray, photon.photon_id, *model);
      path.wavelength_nm = photon.wavelength_nm;
    }
    detected += path.status == obdeect::PhotonStatus::detected;
    const double throughput = path.status == obdeect::PhotonStatus::detected ? 1.0 : 0.0;
    output << "obdeect-arrival-v1," << path.photon_id << ',' << source << ',' << photon.wavelength_nm << ',' << photon.time_ns << ',' << photon.weight << ','
           << throughput << ',' << obdeect::to_string(path.status) << ',' << static_cast<int>(path.point_count) << ','
           << path.path_length_m << ',' << path.incidence_primary_deg << ',' << path.incidence_secondary_deg << ','
           << path.incidence_focal_deg;
    for (const auto& point : path.points_m) output << ',' << point.x << ',' << point.y << ',' << point.z;
    output << '\n';
  }
  if (imported_scene) {
    std::cout << "native scene " << scene_file << ": " << photons_count << " " << source
              << " photons, detected " << detected << "\n";
  } else {
    std::cout << "reference " << model->identifier << ": " << photons_count << " " << source
              << " photons, detected " << detected << "\n";
  }
  return 0;
}
