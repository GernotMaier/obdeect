// Boundary-only driver: link with an explicitly instrumented sim_telarray build.
#include "obdeect/photon_input.hpp"
#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <tuple>

extern "C" {
struct ReferenceMeasurement {
  int status, mirror, loss_line;
  double position[3], direction[3], time_ns, relative_efficiency;
  double primary[3], secondary[3], primary_cosine, secondary_cosine;
  double propagation_group_index, primary_envelope, secondary_envelope;
  double upstream_optical_efficiency;
  int pixel_status, pixel_id;
  double pixel_time_ns, camera_absolute_efficiency, pixel_x_cm, pixel_y_cm;
};
void obdeect_reference_measure(void *, unsigned, const double *, const double *, double, double,
                               double, ReferenceMeasurement *);

// Called only after the simulator has completed its own configuration setup.
int obdeect_reference_replay(void *optics, unsigned count) {
  try {
    const char *input = std::getenv("OBDEECT_REPLAY_INPUT");
    const char *output = std::getenv("OBDEECT_REPLAY_OUTPUT");
    const char *index_text = std::getenv("OBDEECT_REPLAY_TELESCOPE_INDEX");
    if (!input || !output || !index_text)
      throw std::invalid_argument("replay input, output and telescope index must be explicit");
    std::string index_string(index_text);
    std::size_t consumed{};
    const auto index = std::stoul(index_string, &consumed);
    if (consumed != index_string.size() || index_string.starts_with('-') || index >= count)
      throw std::invalid_argument("invalid zero-based replay telescope index");
    const auto source_path = std::filesystem::weakly_canonical(input);
    const auto output_path = std::filesystem::weakly_canonical(output);
    if (source_path == output_path || (std::filesystem::exists(output_path) &&
                                       std::filesystem::equivalent(source_path, output_path)))
      throw std::invalid_argument("replay output aliases photon input");
    obdeect::CsvPhotonReader reader(input);
    std::ofstream rows(output);
    if (!rows)
      throw std::invalid_argument("cannot open reference replay output");
    rows << "run_id,event_id,array_id,telescope_id,bunch_id,photon_id,status,"
            "camera_x_m,camera_y_m,camera_z_m,camera_dx,camera_dy,camera_dz,"
            "travel_time_ns,arrival_time_ns,relative_efficiency,wavelength_nm,source_weight,"
            "mirror_index,loss_source_line,primary_x_m,primary_y_m,primary_z_m,"
            "secondary_x_m,secondary_y_m,secondary_z_m,primary_cosine,secondary_cosine,"
            "propagation_group_index,primary_envelope,secondary_envelope,"
            "upstream_optical_efficiency,pixel_status,pixel_id,pixel_time_ns,"
            "camera_absolute_efficiency,pixel_x_cm,pixel_y_cm\n";
    rows << std::setprecision(17);
    std::array<obdeect::OpticalPhoton, 1024> photons;
    using Identity = std::tuple<std::uint64_t, std::uint64_t, std::uint64_t, std::uint64_t,
                                std::uint64_t, std::uint64_t>;
    std::set<Identity> identities;
    std::optional<std::uint64_t> telescope;
    for (;;) {
      const auto batch = reader.read(photons);
      const auto &c = batch.context;
      if (batch.count && telescope && *telescope != c.telescope_id)
        throw std::invalid_argument("one replay must contain one input telescope identity");
      if (batch.count)
        telescope = c.telescope_id;
      for (std::size_t i = 0; i < batch.count; ++i) {
        const auto &photon = photons[i];
        if (!identities
                 .emplace(c.run_id, c.event_id, c.array_id, c.telescope_id, photon.bunch_id,
                          photon.photon_id)
                 .second)
          throw std::invalid_argument("duplicate replay photon identity");
        if (photon.wavelength_nm <= 0)
          throw std::invalid_argument("reference replay requires resolved positive wavelength");
        const double p[] = {100 * photon.ray.position_m.x, 100 * photon.ray.position_m.y,
                            100 * photon.ray.position_m.z};
        const double d[] = {photon.ray.direction.x, photon.ray.direction.y, photon.ray.direction.z};
        const double distance =
            std::isfinite(photon.emission_distance_m) ? 100 * photon.emission_distance_m : 1e30;
        ReferenceMeasurement result{};
        obdeect_reference_measure(optics, static_cast<unsigned>(index), p, d, photon.wavelength_nm,
                                  distance, photon.weight, &result);
        if (result.status < 0)
          throw std::runtime_error("reference optical trace failed");
        rows << c.run_id << ',' << c.event_id << ',' << c.array_id << ',' << c.telescope_id << ','
             << photon.bunch_id << ',' << photon.photon_id << ','
             << (result.status == 0 ? "hit_camera" : "lost_unclassified") << ',';
        // The legacy ABI does not expose valid terminal fields for lost photons.
        for (double value : result.position) {
          if (result.status == 0)
            rows << value * 0.01;
          rows << ',';
        }
        for (double value : result.direction)
          if (result.status == 0)
            rows << value << ',';
          else
            rows << ',';
        if (result.status == 0)
          rows << result.time_ns << ',' << photon.time_ns + result.time_ns << ','
               << result.relative_efficiency;
        else
          rows << ",,";
        rows << ',' << photon.wavelength_nm << ',' << photon.weight << ',';
        if (result.status == 0)
          rows << result.mirror;
        rows << ',' << result.loss_line << ',';
        for (const auto *point : {result.primary, result.secondary})
          for (int axis = 0; axis < 3; ++axis) {
            if (std::isfinite(point[axis]))
              rows << point[axis] * 0.01;
            rows << ',';
          }
        if (std::isfinite(result.primary[0]))
          rows << result.primary_cosine;
        rows << ',';
        if (std::isfinite(result.secondary[0]))
          rows << result.secondary_cosine;
        rows << ',' << result.propagation_group_index << ',' << result.primary_envelope << ','
             << result.secondary_envelope << ',' << result.upstream_optical_efficiency << ','
             << result.pixel_status << ',' << result.pixel_id << ',' << result.pixel_time_ns << ','
             << result.camera_absolute_efficiency << ',' << result.pixel_x_cm << ','
             << result.pixel_y_cm << '\n';
      }
      if (batch.eof)
        break;
    }
    rows.flush();
    if (!rows)
      throw std::runtime_error("reference output write failed");
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "reference replay: " << error.what() << '\n';
    return 2;
  }
}
}
