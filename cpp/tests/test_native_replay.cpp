#include "obdeect/optical_model_file.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << message << '\n';
    std::exit(1);
  }
}
std::string quoted(const std::string &value) {
  std::string result = "\"";
  for (const auto character : value) {
    if (character == '\\' || character == '"' || character == '$' || character == '`')
      result += '\\';
    result += character;
  }
  return result + "\"";
}
std::string contents(const std::filesystem::path &path) {
  std::ifstream input(path);
  return {(std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>()};
}
std::vector<std::string> fields(const std::string &line) {
  std::vector<std::string> result;
  std::size_t start = 0;
  while (true) {
    const auto end = line.find(',', start);
    result.push_back(line.substr(start, end - start));
    if (end == std::string::npos)
      return result;
    start = end + 1;
  }
}
} // namespace

int main(int argc, char **argv) {
  require(argc == 2, "native executable argument required");
  const auto directory = std::filesystem::current_path() / "native-replay-test";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);
  const auto model = directory / "model.json", photons = directory / "photons.csv";
  auto root =
      obdeect::json::Parser{
          R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"","provenance":{"model":"generic-replay-fixture","model_version":"1"},"report":{"production_trace_ready":false,"trace_blockers":["nominal fixture"]},"trace_model":{"kind":"segmented","primary_facets":[{"id":7,"shape":"circle","centre_m":[0,0,0],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":2,"focal_length_m":10}],"detector_surfaces":[{"id":8,"shape":"circle","centre_m":[0,0,3],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":2}],"cylinder_obscurers":[{"id":9,"first_endpoint_m":[-0.1,0,0.9],"second_endpoint_m":[-0.1,0,1.1],"diameter_m":0.05}],"primary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.8}]}})"}
          .parse();
  require(root.has_value(), "fixture parses");
  auto *hash = const_cast<obdeect::json::Value *>(root->find("optical_model_sha256"));
  hash->string() = obdeect::detail::sha256(obdeect::detail::canonical_json_without_hash(*root));
  {
    std::string json;
    obdeect::detail::append_canonical_json(*root, json);
    std::ofstream output(model);
    output << json;
  }
  {
    std::ofstream output(photons);
    output.precision(17);
    output << "run_id,event_id,array_id,telescope_id,photon_id,x_m,y_m,z_m,dx,dy,dz,wavelength_nm,"
              "time_ns,weight,bunch_id\n"
           << "11,12,13,14,41,0,0,5,0,0,-1,400,7,2,99\n"
           << "11,12,13,14,42,0.5,0,5," << -0.1 / std::sqrt(1.01) << ",0," << -1 / std::sqrt(1.01)
           << ",400,9,3,100\n"
           << "11,12,13,14,43,3,0,5,0,0,-1,400,11,4,101\n";
  }
  const auto run = [&](int block_size, const std::filesystem::path &result, bool gate = false,
                       const std::string &extra = "") {
    const std::string command = quoted(argv[1]) + " --optical-model " + quoted(model.string()) +
                                " --photon-input " + quoted(photons.string()) +
                                " --input-block-size " + std::to_string(block_size) + " --output " +
                                quoted(result.string()) +
                                (gate ? " --require-production-ready" : "") + extra;
    return std::system(command.c_str());
  };
  const auto original_photons = contents(photons), original_model = contents(model);
  require(run(1, photons) != 0 && contents(photons) == original_photons,
          "arrival output cannot overwrite frozen photon input");
  const auto photon_link = directory / "photons-hardlink.csv";
  std::filesystem::create_hard_link(photons, photon_link);
  require(run(1, photon_link) != 0 && contents(photons) == original_photons,
          "hard-linked arrival output cannot overwrite input");
  const auto model_link = directory / "model-symlink.json";
  std::filesystem::create_symlink(model.filename(), model_link);
  require(run(1, model_link) != 0 && contents(model) == original_model,
          "symlinked arrival output cannot overwrite the optical model");
  require(run(1, directory / "protected.csv", false,
              " --interactions-output " + quoted(photons.string())) != 0 &&
              contents(photons) == original_photons,
          "interaction output cannot overwrite input");
  const auto same_output = directory / "same-output.csv";
  require(run(1, same_output, false,
              " --interactions-output " + quoted((directory / "." / "same-output.csv").string())) !=
                  0 &&
              !std::filesystem::exists(same_output),
          "lexical output aliases are rejected before either file is opened");
  const auto first = directory / "first.csv", second = directory / "second.csv";
  require(run(1, first) == 0 && run(2, second) == 0,
          "native replay succeeds across block boundaries");
  require(contents(first) == contents(second), "block size cannot alter frozen photon results");
  std::ifstream input(first);
  std::string line;
  std::getline(input, line);
  const auto header = fields(line);
  std::unordered_map<std::string, std::size_t> column;
  for (std::size_t i = 0; i < header.size(); ++i)
    column.emplace(header[i], i);
  std::size_t row_index = 0;
  while (std::getline(input, line)) {
    const auto row = fields(line);
    const auto get = [&](const char *name) -> const std::string & {
      return row.at(column.at(name));
    };
    require(get("source_kind") == "replay" && get("run_id") == "11" && get("event_id") == "12" &&
                get("telescope_id") == "14",
            "replay identity preserved");
    require(get("photon_id") == std::to_string(41 + row_index), "photon ID survives replay");
    require(get("sampling_area_m2").empty(), "replay cannot invent a sampling-area normalization");
    require(get("source_weight") == std::to_string(2 + row_index), "source weight survives replay");
    if (row_index == 0) {
      require(get("status") == "detected" &&
                  std::abs(std::stod(get("throughput")) - 0.8) < 1.e-12 &&
                  std::abs(std::stod(get("response_loss_fraction")) - 0.2) < 1.e-12 &&
                  get("terminal_loss_fraction") == "0" && get("terminal_surface_id") == "8" &&
                  get("interaction_surface_ids") == "7;8",
              "detector arrival retains response and surfaces");
    } else if (row_index == 1) {
      require(get("status") == "blocked_obscurer" &&
                  std::abs(std::stod(get("response_loss_fraction")) - 0.2) < 1.e-12 &&
                  std::abs(std::stod(get("terminal_loss_fraction")) - 0.8) < 1.e-12 &&
                  get("terminal_surface_id") == "9" && get("point_count") == "3" &&
                  get("interaction_surface_ids") == "7;9",
              "outgoing obscurer retains its endpoint and component");
      const double x0 = std::stod(get("x0_m")), z0 = std::stod(get("z0_m"));
      const double x1 = std::stod(get("x1_m")), z1 = std::stod(get("z1_m"));
      const double x2 = std::stod(get("x2_m")), z2 = std::stod(get("z2_m"));
      require(std::abs(std::stod(get("path_length_m")) - std::hypot(x1 - x0, z1 - z0) -
                       std::hypot(x2 - x1, z2 - z1)) < 1.e-12,
              "recorded loss path ends at actual interaction");
    } else
      require(get("status") == "missed_primary", "off-aperture input retains physical miss");
    ++row_index;
  }
  require(row_index == 3, "all input photons written once");
  require(run(1, directory / "missing-plane.csv", false, " --focal-surface-image") != 0,
          "segmented imaging requires an explicit focal plane");
  root->object().emplace_back("detector_vertex_z_m",
                              obdeect::json::Value{obdeect::json::Value::Number{3.0, "3"}});
  hash = const_cast<obdeect::json::Value *>(root->find("optical_model_sha256"));
  hash->string() = obdeect::detail::sha256(obdeect::detail::canonical_json_without_hash(*root));
  {
    std::string json;
    obdeect::detail::append_canonical_json(*root, json);
    std::ofstream output(model);
    output << json;
  }
  const auto segmented_image_model = contents(model);
  const auto segmented_image_csv = directory / "segmented-image.csv";
  require(run(1, segmented_image_csv, false, " --focal-surface-image") == 0 &&
              contents(model) == segmented_image_model,
          "segmented imaging preserves the optical-model artifact");
  std::ifstream segmented_image_input(segmented_image_csv);
  std::getline(segmented_image_input, line);
  std::getline(segmented_image_input, line);
  const auto segmented_image_row = fields(line);
  require(segmented_image_row[column.at("status")] == "detected" &&
              segmented_image_row[column.at("z2_m")] == "3" &&
              segmented_image_row[column.at("detector_boundary")] == "continuous_focal_surface" &&
              std::abs(std::stod(segmented_image_row[column.at("throughput")]) - 1.0) < 1.e-12,
          "segmented imaging preserves mirror weights and the declared focal-plane position");
  const auto star_area_csv = directory / "star-area.csv";
  const std::string star_area_command =
      quoted(argv[1]) + " --optical-model " + quoted(model.string()) +
      " --source star --photons 1 --sampling-radius-m 2 --output " + quoted(star_area_csv.string());
  require(std::system(star_area_command.c_str()) == 0, "star sampling-area fixture runs");
  std::ifstream star_area_input(star_area_csv);
  std::getline(star_area_input, line);
  std::getline(star_area_input, line);
  require(std::abs(std::stod(fields(line)[column.at("sampling_area_m2")]) - 4 * std::acos(-1.0)) <
              1.e-12,
          "native star output records the exact sampled pupil area");
  const auto keyed_star_area_csv = directory / "keyed-star-area.csv";
  const std::string keyed_star_area_command =
      quoted(argv[1]) + " --optical-model " + quoted(model.string()) +
      " --source star --photons 1 --sampling-radius-m 2 --source-seed 7 --output " +
      quoted(keyed_star_area_csv.string());
  require(std::system(keyed_star_area_command.c_str()) == 0,
          "keyed star sampling-area fixture runs");
  std::ifstream keyed_star_area_input(keyed_star_area_csv);
  std::getline(keyed_star_area_input, line);
  std::getline(keyed_star_area_input, line);
  require(std::abs(std::stod(fields(line)[column.at("sampling_area_m2")]) - 4 * std::acos(-1.0)) <
              1.e-12,
          "keyed star output records the sampled pupil area");
  const auto default_star_area_csv = directory / "default-star-area.csv";
  const std::string default_star_area_command =
      quoted(argv[1]) + " --optical-model " + quoted(model.string()) +
      " --source star --photons 1 --output " + quoted(default_star_area_csv.string());
  require(std::system(default_star_area_command.c_str()) == 0,
          "default star sampling-area fixture runs");
  std::ifstream default_star_area_input(default_star_area_csv);
  std::getline(default_star_area_input, line);
  std::getline(default_star_area_input, line);
  require(std::abs(std::stod(fields(line)[column.at("sampling_area_m2")]) -
                   std::acos(-1.0) * 1.2 * 1.2) < 1.e-12,
          "default star sampling area follows the 1.2-times mirror convention");
  const auto interactions = directory / "interactions.csv";
  require(run(2, directory / "with-diagnostics.csv", false,
              " --interactions-output " + quoted(interactions.string()) +
                  " --interaction-record-limit 1 --interaction-byte-limit 4096") == 0,
          "optional interaction records run through native command");
  require(contents(first) == contents(directory / "with-diagnostics.csv"),
          "diagnostic caps cannot change optical results");
  std::ifstream diagnostics(interactions);
  std::getline(diagnostics, line);
  const auto diagnostic_header = fields(line);
  std::unordered_map<std::string, std::size_t> diagnostic_column;
  for (std::size_t i = 0; i < diagnostic_header.size(); ++i)
    diagnostic_column[diagnostic_header[i]] = i;
  require(static_cast<bool>(std::getline(diagnostics, line)), "one sampled record is stored");
  const auto diagnostic_row = fields(line);
  require(diagnostic_row[diagnostic_column.at("surface_id")] == "7" &&
              diagnostic_row[diagnostic_column.at("kind")] == "mirror" &&
              std::abs(std::stod(diagnostic_row[diagnostic_column.at("outgoing_weight")]) - 1.6) <
                  1.e-12 &&
              std::abs(std::stod(diagnostic_row[diagnostic_column.at("in_dz")]) + 1) < 1.e-12 &&
              std::abs(std::stod(diagnostic_row[diagnostic_column.at("out_dz")]) - 1) < 1.e-12,
          "recorded mirror interaction retains physical normal/directions and attenuation");
  require(!std::getline(diagnostics, line), "interaction count cap is enforced");
  require(run(2, directory / "byte-capped.csv", false,
              " --interactions-output " + quoted(interactions.string()) +
                  " --interaction-record-limit 100 --interaction-byte-limit 450") == 0 &&
              std::filesystem::file_size(interactions) <= 450,
          "interaction byte cap includes header and every record");
  const auto generate = [&](const std::filesystem::path &result, int count, int block,
                            std::uint64_t first_id, const std::string &options) {
    const std::string command = quoted(argv[1]) + " --optical-model " + quoted(model.string()) +
                                " --output " + quoted(result.string()) +
                                " --source-seed 12345 --photons " + std::to_string(count) +
                                " --input-block-size " + std::to_string(block) +
                                " --first-photon-id " + std::to_string(first_id) + options;
    return std::system(command.c_str());
  };
  const auto generated_full = directory / "generated-full.csv";
  const auto generated_blocked = directory / "generated-blocked.csv";
  const auto generated_subset = directory / "generated-subset.csv";
  const auto resolved_source = directory / "resolved-source.csv";
  const auto replayed_source = directory / "replayed-source.csv";
  for (const auto &options :
       {std::string(" --source star --wavelength-nm 350,450 --pulse-width-ns 3"),
        std::string(" --source star --star-mode finite --distance-m 60 --entrance-z-m 50"),
        std::string(" --source laser --beam-radius-m 0.05 --divergence-deg 1"),
        std::string(
            " --source illuminator --emitted-weight 200 --source-normalization-photons 17")}) {
    require(generate(generated_full, 17, 17, 0,
                     options + " --photon-output " + quoted(resolved_source.string())) == 0 &&
                generate(generated_blocked, 17, 2, 0, options) == 0,
            "explicit source mode generates reusable blocks");
    require(contents(generated_full) == contents(generated_blocked),
            "source geometry/time/spectrum/weights independent of generation blocks");
    const std::string source = options.find("--source laser") != std::string::npos ? "laser"
                               : options.find("--source illuminator") != std::string::npos
                                   ? "illuminator"
                                   : "star";
    const auto replay_command = quoted(argv[1]) + " --optical-model " + quoted(model.string()) +
                                " --photon-input " + quoted(resolved_source.string()) +
                                " --source " + source + " --input-block-size 3 --output " +
                                quoted(replayed_source.string());
    require(std::system(replay_command.c_str()) == 0, "resolved photon replay succeeds");
    std::ifstream generated_rows(generated_full), replayed_rows(replayed_source);
    std::string generated_line, replayed_line;
    require(std::getline(generated_rows, generated_line) &&
                std::getline(replayed_rows, replayed_line) && generated_line == replayed_line,
            "generated and replayed arrival headers agree");
    for (int row = 0; row < 17; ++row) {
      require(std::getline(generated_rows, generated_line) &&
                  std::getline(replayed_rows, replayed_line),
              "every resolved photon is replayed");
      auto generated_fields = fields(generated_line), replayed_fields = fields(replayed_line);
      require(generated_fields.at(2) == source && replayed_fields.at(2) == "replay",
              "source provenance distinguishes generation from resolved replay");
      generated_fields.at(2) = "replay";
      const auto sampling_area = column.at("sampling_area_m2");
      if (source == "star") {
        require(!generated_fields.at(sampling_area).empty() &&
                    replayed_fields.at(sampling_area).empty(),
                "generated stars carry launch-area metadata; replay inputs do not invent it");
      }
      generated_fields.at(sampling_area).clear();
      replayed_fields.at(sampling_area).clear();
      require(generated_fields == replayed_fields,
              "resolved replay exactly preserves all physical results and identities");
    }
    require(!std::getline(replayed_rows, replayed_line), "replay has no extra photons");
    require(generate(generated_subset, 7, 2, 5, options) == 0,
            "source partition with explicit photon identities succeeds");
    std::ifstream full(generated_full), subset(generated_subset);
    std::getline(full, line);
    std::getline(subset, line);
    for (int skipped = 0; skipped < 5; ++skipped)
      std::getline(full, line);
    std::string expected;
    for (int row = 0; row < 7; ++row) {
      require(
          static_cast<bool>(std::getline(full, expected)) &&
              static_cast<bool>(std::getline(subset, line)) && expected == line,
          "fixed source identity preserves complete optical result across count/partition changes");
    }
  }
  require(run(2, directory / "source-protected.csv", false,
              " --photon-output " + quoted(photons.string())) != 0 &&
              contents(photons) == original_photons,
          "resolved photon output cannot overwrite its input");
  require(run(2, directory / "gated.csv", true) != 0, "production gate rejects nominal model");
  const auto screen = directory / "screen.csv";
  require(run(2, screen, false,
              " --screen-x-m 0 --screen-y-m 0 --screen-z-m 3 --screen-radius-m 1") == 0,
          "custom detector screen succeeds");
  require(contents(screen).find("7;10") != std::string::npos,
          "custom screen ID follows obscurers as well as optical surfaces");
  auto *report = const_cast<obdeect::json::Value *>(root->find("report"));
  const_cast<obdeect::json::Value *>(report->find("production_trace_ready"))->boolean() = true;
  const_cast<obdeect::json::Value *>(report->find("trace_blockers"))->array().clear();
  hash->string() = obdeect::detail::sha256(obdeect::detail::canonical_json_without_hash(*root));
  {
    std::string json;
    obdeect::detail::append_canonical_json(*root, json);
    std::ofstream output(model);
    output << json;
  }
  require(run(2, directory / "qualified.csv", true) == 0,
          "production gate accepts fixture with a completed report");
  require(run(2, directory / "modified.csv", true, " --panel-id 7") != 0,
          "qualification does not carry over to a diagnostic panel subset");
  require(run(2, directory / "modified-screen.csv", true,
              " --screen-x-m 0 --screen-y-m 0 --screen-z-m 3 --screen-radius-m 1") != 0,
          "qualification does not carry over to a custom detector screen");
  // Native and sidecar paths use phase optical distance and group delay for a
  // declared finite window; vacuum distance is retained separately.
  const auto window_root =
      obdeect::json::Parser{
          R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"","provenance":{"model":"synthetic-finite-window","model_version":"1"},"trace_model":{"kind":"nonsequential","max_interactions":8,"entrance_medium_id":4294967295,"materials":[{"id":1,"phase_index":[{"wavelength_nm":300,"value":1.5},{"wavelength_nm":500,"value":1.5}],"group_index":[{"wavelength_nm":300,"value":1.8},{"wavelength_nm":500,"value":1.8}],"absorption_per_m":[{"wavelength_nm":300,"value":0.2},{"wavelength_nm":500,"value":0.2}]}],"surfaces":[{"id":1,"role":"refractive_interface","shape":"circle","diameter_m":4,"front_medium_id":4294967295,"back_medium_id":1,"frame":{"origin_m":[0,0,2],"x_axis":[1,0,0],"y_axis":[0,1,0],"z_axis":[0,0,1]}},{"id":2,"role":"refractive_interface","shape":"circle","diameter_m":4,"front_medium_id":1,"back_medium_id":4294967295,"frame":{"origin_m":[0,0,1],"x_axis":[1,0,0],"y_axis":[0,1,0],"z_axis":[0,0,1]}},{"id":3,"role":"detector","shape":"circle","diameter_m":4,"frame":{"origin_m":[0,0,0],"x_axis":[1,0,0],"y_axis":[0,1,0],"z_axis":[0,0,1]}}]}})"}
          .parse();
  require(window_root.has_value(), "synthetic finite window schema parses");
  const_cast<obdeect::json::Value *>(window_root->find("optical_model_sha256"))->string() =
      obdeect::detail::sha256(obdeect::detail::canonical_json_without_hash(*window_root));
  {
    std::string json;
    obdeect::detail::append_canonical_json(*window_root, json);
    std::ofstream output(model);
    output << json;
  }
  require(obdeect::read_nonsequential_optical_model(model.string()).has_value(),
          "explicit finite-window optical model loads");
  const auto window_csv = directory / "window.csv",
             window_interactions = directory / "window-interactions.csv";
  require(run(2, window_csv, false,
              " --interactions-output " + quoted(window_interactions.string())) == 0,
          "native finite-window replay succeeds");
  std::ifstream window_input(window_csv);
  std::getline(window_input, line);
  const auto window_header = fields(line);
  std::getline(window_input, line);
  const auto window_row = fields(line);
  const auto window_value = [&](const char *name) {
    const auto found = std::find(window_header.begin(), window_header.end(), name);
    require(found != window_header.end(), "window CSV field exists");
    return window_row.at(found - window_header.begin());
  };
  require(window_value("status") == "detected" &&
              window_value("interaction_surface_ids") == "1;2;3" &&
              std::abs(std::stod(window_value("path_length_m")) - 5) < 1.e-12 &&
              std::abs(std::stod(window_value("optical_path_m")) - 5.5) < 1.e-12 &&
              std::abs(std::stod(window_value("arrival_time_ns")) -
                       (7 + 5.8 / obdeect::kSpeedOfLightMPerNs)) < 1.e-12,
          "native window preserves geometric, phase and group transport separately");
  require(contents(window_interactions).find("refractive_interface") != std::string::npos,
          "actual finite refractive interactions are recorded");
  require(run(1, directory / "unsupported-image.csv", false, " --focal-surface-image") != 0,
          "focal-surface imaging rejects non-axisymmetric optical models");
  auto image_root =
      obdeect::json::Parser{
          R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"","provenance":{"model":"synthetic-imaging-boundary","model_version":"1"},"report":{"production_trace_ready":true,"trace_blockers":[]},"trace_model":{"kind":"axisymmetric","primary_surface_id":10,"secondary_surface_id":11,"detector_surface_id":12,"block_incoming_secondary":false,"primary":{"vertex_z_m":0,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"secondary":{"vertex_z_m":2,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"detector":{"vertex_z_m":1,"inner_radius_m":0,"outer_radius_m":2,"radial_scale_m":1,"coefficient_m":[0,0,0,0,0,0,0,0,0,0,0,0,0]},"detector_surfaces":[{"id":13,"shape":"circle","centre_m":[1,0,1],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":0.1}],"primary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.8}],"secondary_reflectivity":[{"wavelength_nm":300,"response":0.9},{"wavelength_nm":500,"response":0.9}]}})"}
          .parse();
  require(image_root.has_value(), "imaging boundary fixture parses");
  const_cast<obdeect::json::Value *>(image_root->find("optical_model_sha256"))->string() =
      obdeect::detail::sha256(obdeect::detail::canonical_json_without_hash(*image_root));
  {
    std::string json;
    obdeect::detail::append_canonical_json(*image_root, json);
    std::ofstream output(model);
    output << json;
  }
  const auto image_model = contents(model);
  const auto image_csv = directory / "image.csv";
  require(run(1, image_csv, false, " --focal-surface-image") == 0 && contents(model) == image_model,
          "continuous focal-surface diagnostic preserves the compiled optical model");
  std::ifstream image_input(image_csv);
  std::getline(image_input, line);
  const auto image_header = fields(line);
  std::getline(image_input, line);
  const auto image_row = fields(line);
  const auto image_value = [&](const char *name) {
    const auto found = std::find(image_header.begin(), image_header.end(), name);
    require(found != image_header.end(), "imaging CSV field exists");
    return image_row.at(found - image_header.begin());
  };
  require(image_value("status") == "detected" &&
              image_value("detector_boundary") == "continuous_focal_surface" &&
              std::abs(std::stod(image_value("throughput")) - 1.0) < 1.e-12,
          "imaging diagnostic ignores finite pixel acceptance and preserves mirror response");
  require(run(1, directory / "gated-image.csv", true, " --focal-surface-image") != 0,
          "production gate rejects an overridden imaging boundary");
  std::filesystem::remove_all(directory);
}
