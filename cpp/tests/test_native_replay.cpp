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
  std::filesystem::create_directories(directory);
  const auto model = directory / "model.json", photons = directory / "photons.csv";
  auto root =
      obdeect::json::Parser{
          R"({"format":"obdeect.compiled-optical-model.v1","optical_model_sha256":"","provenance":{"model":"generic-replay-fixture","model_version":"1"},"report":{"production_trace_ready":false,"trace_blockers":["nominal fixture"]},"trace_model":{"kind":"segmented","primary_facets":[{"id":7,"shape":"circle","centre_m":[0,0,0],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":2,"focal_length_m":10}],"detector_surfaces":[{"id":8,"shape":"circle","centre_m":[0,0,3],"normal":[0,0,1],"tangent":[1,0,0],"diameter_m":2}],"cylinder_obscurers":[{"id":9,"first_endpoint_m":[-0.1,0,0.9],"second_endpoint_m":[-0.1,0,1.1],"diameter_m":0.05}],"primary_reflectivity":[{"wavelength_nm":300,"response":0.8},{"wavelength_nm":500,"response":0.8}]}})"}
          .parse();
  require(root.has_value(), "fixture parses");
  auto *hash = const_cast<obdeect::json::Value *>(root->find("optical_model_sha256"));
  hash->string = obdeect::detail::sha256(obdeect::detail::canonical_json_without_hash(*root));
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
  for (const auto &options :
       {std::string(" --source star --wavelength-nm 350,450 --pulse-width-ns 3"),
        std::string(" --source star --star-mode finite --distance-m 60 --entrance-z-m 50"),
        std::string(" --source laser --beam-radius-m 0.05 --divergence-deg 1"),
        std::string(
            " --source illuminator --emitted-weight 200 --source-normalization-photons 17")}) {
    require(generate(generated_full, 17, 17, 0, options) == 0 &&
                generate(generated_blocked, 17, 2, 0, options) == 0,
            "explicit source mode generates reusable blocks");
    require(contents(generated_full) == contents(generated_blocked),
            "source geometry/time/spectrum/weights independent of generation blocks");
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
  require(run(2, directory / "gated.csv", true) != 0, "production gate rejects nominal model");
  const auto screen = directory / "screen.csv";
  require(run(2, screen, false,
              " --screen-x-m 0 --screen-y-m 0 --screen-z-m 3 --screen-radius-m 1") == 0,
          "custom detector screen succeeds");
  require(contents(screen).find("7;10") != std::string::npos,
          "custom screen ID follows obscurers as well as optical surfaces");
  auto *report = const_cast<obdeect::json::Value *>(root->find("report"));
  const_cast<obdeect::json::Value *>(report->find("production_trace_ready"))->boolean = true;
  const_cast<obdeect::json::Value *>(report->find("trace_blockers"))->array.clear();
  hash->string = obdeect::detail::sha256(obdeect::detail::canonical_json_without_hash(*root));
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
  std::filesystem::remove_all(directory);
}
