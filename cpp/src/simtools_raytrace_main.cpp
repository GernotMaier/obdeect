#include "obdeect/cli_parse.hpp"
#include "obdeect/ctao_models.hpp"
#include "obdeect/ctao_trace.hpp"
#include "obdeect/interactions.hpp"
#include "obdeect/optical_model_file.hpp"
#include "obdeect/segmented_optical_model.hpp"
#include "obdeect/sources.hpp"

#include <fstream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

namespace {
void usage() {
  std::cout << "Usage: obdeect-simtools-raytrace [--telescope LST|MST|SST|SCT] [options]\n"
            << "  --optical-model FILE  compiled obdeect optical-model JSON\n"
            << "  --source star|illuminator|laser  (default: star)\n"
            << "  --photons N --output FILE --field-x-deg D --field-y-deg D\n"
            << "  --distance-m D --wavelength-nm N[,N...] --divergence-deg D --panel-id N\n"
            << "  --source-x-m D --source-y-m D --source-z-m D\n"
            << "  --screen-x-m D --screen-y-m D --screen-z-m D --screen-radius-m D\n"
            << "  --direction-x D --direction-y D --direction-z D  (laser axis)\n"
            << "  --emission-time-ns D --pulse-width-ns D  (deterministic top-hat pulse)\n"
            << "  --optical-model is required for model-derived CTAO simulations; it accepts general panel geometry.\n"
            << "  Without it, this runs an analytic diagnostic prescription only.\n";
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
  std::string optical_model_path;
  std::string source{"star"};
  std::string output_path{"obdeect_paths.csv"};
  std::size_t photons_count = 10000;
  double field_x_deg = 0.0, field_y_deg = 0.0, distance_m = 10000.0;
  std::vector<double> wavelengths_nm{400.0};
  double divergence_deg = 0.0;
  double source_x_m = 0.0, source_y_m = 0.0, source_z_m = 50.0;
  double direction_x = 0.0, direction_y = 0.0, direction_z = -1.0;
  double emission_time_ns = 0.0, pulse_width_ns = 0.0;
  double screen_x_m = 0.0, screen_y_m = 0.0, screen_z_m = 0.0, screen_radius_m = 0.0;
  bool screen_x_set = false, screen_y_set = false, screen_z_set = false, screen_radius_set = false;
  std::optional<std::uint32_t> panel_id;
  for (int index = 1; index < argc; ++index) {
    std::string value_string;
    if (value(index, argc, argv, "--telescope", telescope) ||
        value(index, argc, argv, "--optical-model", optical_model_path) ||
        value(index, argc, argv, "--source", source) ||
        value(index, argc, argv, "--output", output_path)) {
      continue;
    }
    if (std::string{argv[index]} == "-h" || std::string{argv[index]} == "--help") {
      usage();
      return 0;
    }
    if (std::string{argv[index]} == "--photons" && index + 1 < argc &&
        parse_positive_size(argv[++index], photons_count)) continue;
    if (std::string{argv[index]} == "--wavelength-nm" && index + 1 < argc &&
        obdeect::parse_wavelengths_nm(argv[++index], wavelengths_nm)) continue;
    if (std::string{argv[index]} == "--panel-id" && index + 1 < argc) {
      std::uint32_t parsed{};
      if (!obdeect::parse_uint32(argv[++index], parsed)) {
        usage();
        return 2;
      }
      panel_id = parsed;
      continue;
    }
    auto number = [&](const char* option, double& target, bool* provided = nullptr) {
      if (std::string{argv[index]} != option || index + 1 >= argc ||
          !parse_finite_double(argv[index + 1], target)) {
        return false;
      }
      ++index;
      if (provided) *provided = true;
      return true;
    };
    if (number("--field-x-deg", field_x_deg) || number("--field-y-deg", field_y_deg) ||
        number("--distance-m", distance_m) ||
        number("--divergence-deg", divergence_deg) || number("--source-x-m", source_x_m) ||
        number("--source-y-m", source_y_m) || number("--source-z-m", source_z_m) ||
        number("--direction-x", direction_x) || number("--direction-y", direction_y) ||
        number("--direction-z", direction_z) || number("--emission-time-ns", emission_time_ns) ||
        number("--pulse-width-ns", pulse_width_ns) ||
        number("--screen-x-m", screen_x_m, &screen_x_set) ||
        number("--screen-y-m", screen_y_m, &screen_y_set) ||
        number("--screen-z-m", screen_z_m, &screen_z_set) ||
        number("--screen-radius-m", screen_radius_m, &screen_radius_set)) continue;
    usage();
    return 2;
  }
  const auto model = optical_model_path.empty() ? obdeect::ctao_reference_model(telescope) : std::nullopt;
  auto imported_optical_model = optical_model_path.empty() ? std::optional<obdeect::CompiledSegmentedOpticalModel>{}
                                           : obdeect::read_segmented_optical_model(optical_model_path);
  const auto axisymmetric_optical_model = optical_model_path.empty() ? std::optional<obdeect::AxisymmetricOpticalModel>{}
                                                      : obdeect::read_axisymmetric_optical_model(optical_model_path);
  if ((!model && !imported_optical_model && !axisymmetric_optical_model) ||
      (source != "star" && source != "illuminator" && source != "laser") ||
      distance_m <= 0.0 || divergence_deg < 0.0 || divergence_deg >= 90.0 || pulse_width_ns < 0.0 ||
      (panel_id && !imported_optical_model) ||
      ((screen_x_set || screen_y_set || screen_z_set || screen_radius_set) &&
       !(screen_x_set && screen_y_set && screen_z_set && screen_radius_set)) ||
      (screen_radius_set && screen_radius_m <= 0.0)) {
    usage();
    return 2;
  }
  if (panel_id) {
    auto& facets = imported_optical_model->primary_facets;
    facets.erase(std::remove_if(facets.begin(), facets.end(), [&](const auto& facet) {
      return facet.id != *panel_id;
    }), facets.end());
    if (facets.empty()) {
      std::cerr << "panel " << *panel_id << " is not present in " << optical_model_path << '\n';
      return 2;
    }
  }
  if (screen_x_set) {
    if (!imported_optical_model) {
      std::cerr << "a custom screen requires --optical-model\n";
      return 2;
    }
    std::uint32_t maximum_id = 0;
    for (const auto& facet : imported_optical_model->primary_facets) maximum_id = std::max(maximum_id, facet.id);
    for (const auto& detector : imported_optical_model->detector_surfaces)
      maximum_id = std::max(maximum_id, detector.id);
    if (maximum_id == std::numeric_limits<std::uint32_t>::max()) {
      std::cerr << "cannot allocate a custom screen ID\n";
      return 2;
    }
    imported_optical_model->detector_surfaces = {{maximum_id + 1,
                                          {screen_x_m, screen_y_m, screen_z_m},
                                          {0.0, 0.0, 1.0},
                                          2.0 * screen_radius_m,
                                          obdeect::FacetShape::circle,
                                          {1.0, 0.0, 0.0}}};
  }
  const double radians_per_degree = std::numbers::pi / 180.0;
  const double pupil_radius = imported_optical_model
                                  ? [&] {
                                      double radius = 0.0;
                                      for (const auto& facet : imported_optical_model->primary_facets)
                                        radius = std::max(radius, std::hypot(facet.centre_m.x, facet.centre_m.y) +
                                                                  facet.diameter_m * 0.5);
                                      return radius;
                                    }()
                                  : axisymmetric_optical_model ? axisymmetric_optical_model->primary.outer_radius_m
                                                       : model->primary_outer_radius_m;
  if (!std::isfinite(pupil_radius) || pupil_radius <= 0.0) return 1;
  std::vector<obdeect::OpticalPhoton> input;
  if (source == "star") {
    input = obdeect::star_photons(photons_count, pupil_radius,
                                  {field_x_deg * radians_per_degree, field_y_deg * radians_per_degree,
                                   distance_m, wavelengths_nm.front()});
  } else if (source == "illuminator") {
    input = obdeect::illuminator_photons(
        photons_count, pupil_radius,
        {{source_x_m, source_y_m, source_z_m}, wavelengths_nm.front(), 1.0});
  } else {
    input = obdeect::laser_photons(
        photons_count, pupil_radius,
        {{direction_x, direction_y, direction_z}, {source_x_m, source_y_m, source_z_m},
         divergence_deg * radians_per_degree, wavelengths_nm.front()});
  }
  if (input.size() != photons_count) {
    std::cerr << "source configuration produced no valid photons\n";
    return 1;
  }
  for (std::size_t index = 0; index < input.size(); ++index)
    input[index].wavelength_nm = wavelengths_nm[index % wavelengths_nm.size()];
  if (!obdeect::apply_top_hat_emission_times(input, emission_time_ns, pulse_width_ns)) {
    std::cerr << "invalid emission-time or pulse-width value\n";
    return 2;
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
    if (imported_optical_model) {
      const auto direction = obdeect::normalised_checked(photon.ray.direction);
      if (!direction) {
        path.status = obdeect::PhotonStatus::invalid_input;
      } else {
        const obdeect::Ray ray{photon.ray.position_m, *direction};
        const auto incoming_obscurer = obdeect::intersect_cylinder_obscurers_unchecked(ray, *imported_optical_model);
        const auto primary_hit = obdeect::intersect_segmented_primary_unchecked(ray, *imported_optical_model);
        if (incoming_obscurer && (!primary_hit || incoming_obscurer->distance_m < primary_hit->distance_m)) {
          path.status = obdeect::PhotonStatus::blocked_obscurer;
          path.final_direction = *direction;
        } else if (!primary_hit) {
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
                obdeect::Ray{primary_hit->point_m, *reflected}, *imported_optical_model);
            const obdeect::Ray reflected_ray{primary_hit->point_m, *reflected};
            const auto outgoing_obscurer =
                obdeect::intersect_cylinder_obscurers_unchecked(reflected_ray, *imported_optical_model);
            if (outgoing_obscurer && (!detector_hit || outgoing_obscurer->distance_m < detector_hit->distance_m)) {
              path.status = obdeect::PhotonStatus::blocked_obscurer;
              path.path_length_m += outgoing_obscurer->distance_m;
              path.final_direction = *reflected;
            } else if (!detector_hit) {
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
    } else if (axisymmetric_optical_model) {
      path = obdeect::trace_axisymmetric_optical_model(photon.ray, photon.photon_id, *axisymmetric_optical_model);
      path.wavelength_nm = photon.wavelength_nm;
    } else {
      path = obdeect::trace_ctao_reference(photon.ray, photon.photon_id, *model);
      path.wavelength_nm = photon.wavelength_nm;
    }
    double throughput = path.status == obdeect::PhotonStatus::detected ? 1.0 : 0.0;
    if (path.status == obdeect::PhotonStatus::detected && imported_optical_model && imported_optical_model->primary_reflectivity) {
      const auto response = imported_optical_model->primary_reflectivity->at(photon.wavelength_nm);
      if (!response) {
        path.status = obdeect::PhotonStatus::invalid_input;
        throughput = 0.0;
      } else {
        throughput = *response;
      }
    }
    if (path.status == obdeect::PhotonStatus::detected && axisymmetric_optical_model) {
      for (const auto* response : {axisymmetric_optical_model->primary_reflectivity ? &*axisymmetric_optical_model->primary_reflectivity : nullptr,
                                   axisymmetric_optical_model->secondary_reflectivity ? &*axisymmetric_optical_model->secondary_reflectivity : nullptr}) {
        if (!response) continue;
        const auto value = response->at(photon.wavelength_nm);
        if (!value) {
          path.status = obdeect::PhotonStatus::invalid_input;
          throughput = 0.0;
          break;
        }
        throughput *= *value;
      }
    }
    detected += path.status == obdeect::PhotonStatus::detected;
    output << "obdeect-arrival-v1," << path.photon_id << ',' << source << ',' << photon.wavelength_nm << ',' << photon.time_ns << ',' << photon.weight << ','
           << throughput << ',' << obdeect::to_string(path.status) << ',' << static_cast<int>(path.point_count) << ','
           << path.path_length_m << ',' << path.incidence_primary_deg << ',' << path.incidence_secondary_deg << ','
           << path.incidence_focal_deg;
    for (const auto& point : path.points_m) output << ',' << point.x << ',' << point.y << ',' << point.z;
    output << '\n';
  }
  if (imported_optical_model || axisymmetric_optical_model) {
    std::cout << "optical model " << optical_model_path << ": " << photons_count << " " << source
              << " photons, detected " << detected;
    if (panel_id) std::cout << ", panel " << *panel_id;
    std::cout << "\n";
  } else {
    std::cout << "reference " << model->identifier << ": " << photons_count << " " << source
              << " photons, detected " << detected << "\n";
  }
  return 0;
}
