#include "obdeect/cli_parse.hpp"
#include "obdeect/ctao_models.hpp"
#include "obdeect/ctao_trace.hpp"
#include "obdeect/interaction_csv.hpp"
#include "obdeect/interactions.hpp"
#include "obdeect/optical_model_file.hpp"
#include "obdeect/photon_input.hpp"
#include "obdeect/segmented_optical_model.hpp"
#include "obdeect/segmented_path.hpp"
#include "obdeect/source_sampling.hpp"
#include "obdeect/sources.hpp"
#include <charconv>

#include <memory>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

namespace {
void usage() {
  std::cout
      << "Usage: obdeect-simtools-raytrace [--telescope LST|MST|SST|SCT] [options]\n"
      << "  --optical-model FILE  compiled obdeect optical-model JSON\n"
      << "  --require-production-ready  reject optical models lacking a completed production "
         "report\n"
      << "  --photon-input FILE --input-block-size N  replay telescope-local CSV photons\n"
      << "  --interactions-output FILE --interaction-record-limit N --interaction-byte-limit N\n"
      << "  --source star|illuminator|laser  (default: star)\n"
      << "  --photons N --output FILE --field-x-deg D --field-y-deg D\n"
      << "  --distance-m D --wavelength-nm N[,N...] --divergence-deg D --panel-id N\n"
      << "  --source-x-m D --source-y-m D --source-z-m D\n"
      << "  --screen-x-m D --screen-y-m D --screen-z-m D --screen-radius-m D\n"
      << "  --direction-x D --direction-y D --direction-z D  (laser axis)\n"
      << "  --emission-time-ns D --pulse-width-ns D  (deterministic top-hat pulse)\n"
      << "  --optical-model is required for model-derived CTAO simulations; it accepts "
         "general panel geometry.\n"
      << "  Without it, this runs an analytic diagnostic prescription only.\n";
}

bool aliases(const std::string &first, const std::string &second) {
  if (first.empty() || second.empty())
    return false;
  std::error_code error;
  if (std::filesystem::equivalent(first, second, error))
    return true;
  error.clear();
  const auto canonical_first = std::filesystem::weakly_canonical(first, error);
  if (error)
    return false;
  const auto canonical_second = std::filesystem::weakly_canonical(second, error);
  return !error && canonical_first == canonical_second;
}

bool value(int &index, int argc, char **argv, const char *option, std::string &result) {
  if (std::string{argv[index]} != option || index + 1 >= argc)
    return false;
  result = argv[++index];
  return true;
}
} // namespace

int main(int argc, char **argv) {
  using obdeect::parse_finite_double;
  using obdeect::parse_positive_size;
  std::optional<std::uint64_t> source_normalization_photons;
  double emitted_weight = 1;
  bool emitted_weight_set = false;
  std::optional<std::uint64_t> source_seed;
  std::uint64_t first_photon_id = 0;
  std::string star_mode = "plane-wave";
  double beam_radius_m = 0, entrance_z_m = 50;
  bool beam_radius_set = false, first_id_set = false;
  bool require_production_ready = false;
  bool production_ready = false;
  std::string interactions_output_path;
  std::size_t interaction_record_limit = 10000, interaction_byte_limit = 4 * 1024 * 1024;
  std::string photon_input_path;
  std::size_t input_block_size = 4096;
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
        value(index, argc, argv, "--star-mode", star_mode) ||
        value(index, argc, argv, "--source", source) ||
        value(index, argc, argv, "--interactions-output", interactions_output_path) ||
        value(index, argc, argv, "--photon-input", photon_input_path) ||
        value(index, argc, argv, "--output", output_path)) {
      continue;
    }
    if (std::string{argv[index]} == "--require-production-ready") {
      require_production_ready = true;
      continue;
    }
    if (std::string{argv[index]} == "-h" || std::string{argv[index]} == "--help") {
      usage();
      return 0;
    }
    if ((std::string{argv[index]} == "--source-seed" ||
         std::string{argv[index]} == "--first-photon-id" ||
         std::string{argv[index]} == "--source-normalization-photons") &&
        index + 1 < argc) {
      const bool normalization_option =
          std::string{argv[index]} == "--source-normalization-photons";
      const bool seed_option = std::string{argv[index]} == "--source-seed";
      const std::string text = argv[++index];
      std::uint64_t parsed{};
      const auto result = std::from_chars(text.data(), text.data() + text.size(), parsed);
      if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
        usage();
        return 2;
      }
      if (normalization_option) {
        if (parsed == 0) {
          usage();
          return 2;
        }
        source_normalization_photons = parsed;
      } else if (seed_option)
        source_seed = parsed;
      else {
        first_photon_id = parsed;
        first_id_set = true;
      }
      continue;
    }
    if (std::string{argv[index]} == "--photons" && index + 1 < argc &&
        parse_positive_size(argv[++index], photons_count))
      continue;
    if (std::string{argv[index]} == "--input-block-size" && index + 1 < argc &&
        parse_positive_size(argv[++index], input_block_size))
      continue;
    if (std::string{argv[index]} == "--interaction-record-limit" && index + 1 < argc &&
        parse_positive_size(argv[++index], interaction_record_limit))
      continue;
    if (std::string{argv[index]} == "--interaction-byte-limit" && index + 1 < argc &&
        parse_positive_size(argv[++index], interaction_byte_limit))
      continue;
    if (std::string{argv[index]} == "--wavelength-nm" && index + 1 < argc &&
        obdeect::parse_wavelengths_nm(argv[++index], wavelengths_nm))
      continue;
    if (std::string{argv[index]} == "--panel-id" && index + 1 < argc) {
      std::uint32_t parsed{};
      if (!obdeect::parse_uint32(argv[++index], parsed)) {
        usage();
        return 2;
      }
      panel_id = parsed;
      continue;
    }
    auto number = [&](const char *option, double &target, bool *provided = nullptr) {
      if (std::string{argv[index]} != option || index + 1 >= argc ||
          !parse_finite_double(argv[index + 1], target)) {
        return false;
      }
      ++index;
      if (provided)
        *provided = true;
      return true;
    };
    if (number("--emitted-weight", emitted_weight, &emitted_weight_set) ||
        number("--beam-radius-m", beam_radius_m, &beam_radius_set) ||
        number("--entrance-z-m", entrance_z_m) || number("--field-x-deg", field_x_deg) ||
        number("--field-y-deg", field_y_deg) || number("--distance-m", distance_m) ||
        number("--divergence-deg", divergence_deg) || number("--source-x-m", source_x_m) ||
        number("--source-y-m", source_y_m) || number("--source-z-m", source_z_m) ||
        number("--direction-x", direction_x) || number("--direction-y", direction_y) ||
        number("--direction-z", direction_z) || number("--emission-time-ns", emission_time_ns) ||
        number("--pulse-width-ns", pulse_width_ns) ||
        number("--screen-x-m", screen_x_m, &screen_x_set) ||
        number("--screen-y-m", screen_y_m, &screen_y_set) ||
        number("--screen-z-m", screen_z_m, &screen_z_set) ||
        number("--screen-radius-m", screen_radius_m, &screen_radius_set))
      continue;
    usage();
    return 2;
  }
  const auto model =
      optical_model_path.empty() ? obdeect::ctao_reference_model(telescope) : std::nullopt;
  auto imported_optical_model = optical_model_path.empty()
                                    ? std::optional<obdeect::CompiledSegmentedOpticalModel>{}
                                    : obdeect::read_segmented_optical_model(optical_model_path);
  const auto axisymmetric_optical_model =
      optical_model_path.empty() ? std::optional<obdeect::AxisymmetricOpticalModel>{}
                                 : obdeect::read_axisymmetric_optical_model(optical_model_path);
  if ((!model && !imported_optical_model && !axisymmetric_optical_model) ||
      (source != "star" && source != "illuminator" && source != "laser") || distance_m <= 0.0 ||
      divergence_deg < 0.0 || divergence_deg >= 90.0 || pulse_width_ns < 0.0 ||
      (panel_id && !imported_optical_model) ||
      ((screen_x_set || screen_y_set || screen_z_set || screen_radius_set) &&
       !(screen_x_set && screen_y_set && screen_z_set && screen_radius_set)) ||
      (screen_radius_set && screen_radius_m <= 0.0)) {
    usage();
    return 2;
  }
  const bool sampled_source = source_seed.has_value() || first_id_set || star_mode == "finite" ||
                              beam_radius_set || source_normalization_photons.has_value() ||
                              emitted_weight_set;
  if ((star_mode != "plane-wave" && star_mode != "finite") ||
      (star_mode == "finite" && source != "star") ||
      (beam_radius_set && (source != "laser" || beam_radius_m <= 0)) ||
      (source != "illuminator" && (source_normalization_photons || emitted_weight_set)) ||
      emitted_weight < 0 || (!photon_input_path.empty() && sampled_source) || entrance_z_m < 0 ||
      (star_mode == "finite" && distance_m <= entrance_z_m) ||
      (sampled_source &&
       photons_count - 1 > std::numeric_limits<std::uint64_t>::max() - first_photon_id)) {
    std::cerr << "invalid keyed source, finite star, beam radius, or photon ID range\n";
    return 2;
  }
  if (aliases(output_path, photon_input_path) || aliases(output_path, optical_model_path) ||
      aliases(interactions_output_path, photon_input_path) ||
      aliases(interactions_output_path, optical_model_path) ||
      aliases(output_path, interactions_output_path)) {
    std::cerr << "output files must not alias each other, the photon input, or the optical model\n";
    return 2;
  }
  if (!optical_model_path.empty()) {
    const auto root = obdeect::detail::read_json(optical_model_path);
    const auto report = root ? root->find("report") : nullptr;
    const auto readiness = report ? report->find("production_trace_ready") : nullptr;
    const auto blockers = report ? report->find("trace_blockers") : nullptr;
    production_ready =
        readiness && readiness->kind == obdeect::json::Value::Kind::boolean && readiness->boolean &&
        blockers && blockers->kind == obdeect::json::Value::Kind::array && blockers->array.empty();
  }
  // A reviewed artifact's qualification does not cover diagnostic geometry overrides.
  production_ready = production_ready && !panel_id && !screen_x_set;
  if (require_production_ready && !production_ready) {
    std::cerr << "optical model has no completed production readiness report or retains production "
                 "blockers\n";
    return 2;
  }
  if (!interactions_output_path.empty() &&
      ((!imported_optical_model && !axisymmetric_optical_model) ||
       interactions_output_path == output_path)) {
    std::cerr
        << "interaction diagnostics require a compiled optical model and a separate output file\n";
    return 2;
  }
  if (panel_id) {
    auto &facets = imported_optical_model->primary_facets;
    facets.erase(std::remove_if(facets.begin(), facets.end(),
                                [&](const auto &facet) { return facet.id != *panel_id; }),
                 facets.end());
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
    for (const auto &facet : imported_optical_model->primary_facets)
      maximum_id = std::max(maximum_id, facet.id);
    for (const auto &detector : imported_optical_model->detector_surfaces)
      maximum_id = std::max(maximum_id, detector.id);
    for (const auto &obscurer : imported_optical_model->cylinder_obscurers)
      maximum_id = std::max(maximum_id, obscurer.id);
    if (maximum_id >= obdeect::PhotonResultBlock::kNoSurfaceId - 1) {
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
  if (!std::isfinite(pupil_radius) || pupil_radius <= 0.0)
    return 1;
  std::vector<obdeect::OpticalPhoton> input;
  if (sampled_source) {
    input.resize(std::min(input_block_size, photons_count));
  } else if (!photon_input_path.empty()) {
    source = "replay";
    input.resize(input_block_size);
  } else if (source == "star") {
    input =
        obdeect::star_photons(photons_count, pupil_radius,
                              {field_x_deg * radians_per_degree, field_y_deg * radians_per_degree,
                               distance_m, wavelengths_nm.front()});
  } else if (source == "illuminator") {
    input = obdeect::illuminator_photons(
        photons_count, pupil_radius,
        {{source_x_m, source_y_m, source_z_m}, wavelengths_nm.front(), 1.0});
  } else {
    input = obdeect::laser_photons(photons_count, pupil_radius,
                                   {{direction_x, direction_y, direction_z},
                                    {source_x_m, source_y_m, source_z_m},
                                    divergence_deg * radians_per_degree,
                                    wavelengths_nm.front()});
  }
  if (!sampled_source && photon_input_path.empty() && input.size() != photons_count) {
    std::cerr << "source configuration produced no valid photons\n";
    return 1;
  }
  if (!sampled_source && photon_input_path.empty())
    for (std::size_t index = 0; index < input.size(); ++index)
      input[index].wavelength_nm = wavelengths_nm[index % wavelengths_nm.size()];
  if (!sampled_source && photon_input_path.empty() &&
      !obdeect::apply_top_hat_emission_times(input, emission_time_ns, pulse_width_ns)) {
    std::cerr << "invalid emission-time or pulse-width value\n";
    return 2;
  }
  std::unique_ptr<obdeect::InteractionCsvWriter> interaction_writer;
  try {
    if (!interactions_output_path.empty())
      interaction_writer = std::make_unique<obdeect::InteractionCsvWriter>(
          interactions_output_path, interaction_record_limit, interaction_byte_limit);
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
  std::ofstream output{output_path};
  if (!output) {
    std::cerr << "cannot write " << output_path << '\n';
    return 1;
  }
  output << std::setprecision(17)
         << "contract_version,photon_id,source_kind,wavelength_nm,emission_time_ns,source_weight,"
            "throughput,status,point_count,path_length_m,"
            "incidence_primary_deg,incidence_secondary_deg,incidence_focal_deg,";
  for (int point = 0; point < 4; ++point)
    output << "x" << point << "_m,y" << point << "_m,z" << point << "_m" << ',';
  output << "run_id,event_id,array_id,telescope_id,bunch_id,arrival_time_ns,terminal_surface_id,"
            "final_dx,final_dy,final_dz,interaction_surface_ids,response_loss_fraction,terminal_"
            "loss_fraction\n";
  std::size_t detected = 0;
  std::size_t traced_count = 0;
  std::unique_ptr<obdeect::CsvPhotonReader> reader;
  try {
    if (!photon_input_path.empty())
      reader = std::make_unique<obdeect::CsvPhotonReader>(photon_input_path);
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
  bool complete = false;
  while (!complete) {
    obdeect::PhotonBatchContext context{};
    if (sampled_source) {
      input.resize(std::min(input_block_size, photons_count - traced_count));
      const bool generated = obdeect::fill_source(
          input, first_photon_id + traced_count,
          [&](std::uint64_t id) -> std::optional<obdeect::OpticalPhoton> {
            const auto seed = source_seed.value_or(0);
            std::optional<obdeect::OpticalPhoton> photon;
            if (source == "star" && star_mode == "finite") {
              photon = obdeect::sample_star(
                  id, seed, pupil_radius,
                  obdeect::FiniteStarSource{
                      {-distance_m * std::tan(field_x_deg * radians_per_degree),
                       -distance_m * std::tan(field_y_deg * radians_per_degree), distance_m},
                      entrance_z_m,
                      wavelengths_nm.front(),
                      emission_time_ns});
            } else if (source == "star") {
              photon =
                  obdeect::sample_star(id, seed, pupil_radius,
                                       obdeect::StarSource{field_x_deg * radians_per_degree,
                                                           field_y_deg * radians_per_degree,
                                                           distance_m, wavelengths_nm.front()});
              if (photon)
                photon->time_ns += emission_time_ns;
            } else if (source == "illuminator") {
              photon = obdeect::sample_illuminator(
                  id, seed, pupil_radius,
                  obdeect::PointIlluminator{
                      {source_x_m, source_y_m, source_z_m}, wavelengths_nm.front(), emitted_weight},
                  source_normalization_photons.value_or(photons_count));
              if (photon)
                photon->time_ns += emission_time_ns;
            } else {
              photon = obdeect::sample_laser(
                  id, seed, beam_radius_set ? beam_radius_m : pupil_radius,
                  obdeect::LaserSource{{direction_x, direction_y, direction_z},
                                       {source_x_m, source_y_m, source_z_m},
                                       divergence_deg * radians_per_degree,
                                       wavelengths_nm.front()});
              if (photon)
                photon->time_ns += emission_time_ns;
            }
            if (photon) {
              const auto wavelength_index = static_cast<std::size_t>(
                  obdeect::source_uniform(id, seed, 4) * wavelengths_nm.size());
              photon->wavelength_nm = wavelengths_nm[wavelength_index];
              photon->time_ns += pulse_width_ns * (obdeect::source_uniform(id, seed, 5) - 0.5);
            }
            return photon;
          });
      if (!generated) {
        std::cerr << "invalid keyed source configuration\n";
        return 2;
      }
      complete = traced_count + input.size() == photons_count;
    } else if (reader) {
      input.resize(input_block_size);
      try {
        const auto batch = reader->read(input);
        context = batch.context;
        input.resize(batch.count);
        complete = batch.eof;
      } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 2;
      }
    } else
      complete = true;
    traced_count += input.size();
    for (const auto &photon : input) {
      obdeect::PathRecord path{};
      path.photon_id = photon.photon_id;
      path.wavelength_nm = photon.wavelength_nm;
      path.points_m[0] = photon.ray.position_m;
      path.point_count = 1;
      if (imported_optical_model) {
        path = obdeect::trace_segmented_path(photon.ray, photon.photon_id, photon.wavelength_nm,
                                             *imported_optical_model);
      } else if (axisymmetric_optical_model) {
        path = obdeect::trace_axisymmetric_optical_model(photon.ray, photon.photon_id,
                                                         *axisymmetric_optical_model);
        path.wavelength_nm = photon.wavelength_nm;
      } else {
        path = obdeect::trace_ctao_reference(photon.ray, photon.photon_id, *model);
        path.wavelength_nm = photon.wavelength_nm;
      }
      double throughput = path.status == obdeect::PhotonStatus::detected ? 1.0 : 0.0;
      if (path.status == obdeect::PhotonStatus::detected && imported_optical_model)
        throughput = path.surviving_throughput;
      if (path.status == obdeect::PhotonStatus::detected && axisymmetric_optical_model)
        throughput = path.surviving_throughput;
      if (interaction_writer) {
        try {
          interaction_writer->write(context, photon, path);
        } catch (const std::exception &error) {
          std::cerr << error.what() << '\n';
          return 1;
        }
      }
      detected += path.status == obdeect::PhotonStatus::detected;
      output << "obdeect-arrival-v1," << path.photon_id << ',' << source << ','
             << photon.wavelength_nm << ',' << photon.time_ns << ',' << photon.weight << ','
             << throughput << ',' << obdeect::to_string(path.status) << ','
             << static_cast<int>(path.point_count) << ',' << path.path_length_m << ','
             << path.incidence_primary_deg << ',' << path.incidence_secondary_deg << ','
             << path.incidence_focal_deg;
      for (const auto &point : path.points_m)
        output << ',' << point.x << ',' << point.y << ',' << point.z;
      output << ',' << context.run_id << ',' << context.event_id << ',' << context.array_id << ','
             << context.telescope_id << ',' << photon.bunch_id << ','
             << photon.time_ns + path.path_length_m / obdeect::kSpeedOfLightMPerNs << ','
             << path.terminal_surface_id << ',' << path.final_direction.x << ','
             << path.final_direction.y << ',' << path.final_direction.z << ',';
      for (std::size_t interaction = 0; interaction + 1 < path.point_count; ++interaction) {
        if (interaction != 0)
          output << ';';
        output << path.interaction_surface_ids[interaction];
      }
      output << ',' << 1.0 - path.surviving_throughput << ','
             << (path.status == obdeect::PhotonStatus::detected ? 0.0 : path.surviving_throughput)
             << '\n';
    }
  }
  photons_count = traced_count;
  if (interaction_writer)
    std::cout << "interaction diagnostics: " << interaction_writer->records_written() << " records"
              << (interaction_writer->truncated() ? " (storage cap reached)" : "") << '\n';
  if (imported_optical_model || axisymmetric_optical_model) {
    std::cout << "optical model " << optical_model_path << ": " << photons_count << " " << source
              << " photons, detected " << detected
              << (production_ready ? ", production readiness declared"
                                   : ", nominal optical transport; production unvalidated");
    if (panel_id)
      std::cout << ", panel " << *panel_id;
    std::cout << "\n";
  } else {
    std::cout << "reference " << model->identifier << ": " << photons_count << " " << source
              << " photons, detected " << detected << "\n";
  }
  return 0;
}
