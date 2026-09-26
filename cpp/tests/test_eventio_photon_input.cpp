#include "obdeect/eventio_photon_input.hpp"

extern "C" {
#include "initial.h"
#include "io_basic.h"
#include "mc_tel.h"
}

#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

void flush(IO_BUFFER* buffer) {
  // Ending a top-level EventIO item writes the block automatically.
  require(std::fflush(buffer->output_file) == 0, "cannot flush EventIO test file");
}

struct Fixture {
  std::filesystem::path path = std::filesystem::temp_directory_path() / "obdeect-corsika-eventio-test.dat";
  Fixture(bool ceffic = false, bool compact = false) {
    std::FILE* file = std::fopen(path.string().c_str(), "wb");
    require(file != nullptr, "cannot open test file");
    IO_BUFFER* buffer = allocate_io_buffer(32768);
    require(buffer != nullptr, "cannot allocate test buffer");
    buffer->output_file = file;
    std::array<real, 273> run{};
    std::memcpy(run.data(), "RUNH", 4);
    run[1] = 7;
    run[4] = 1;
    run[5] = 200000;
    require(write_tel_block(buffer, IO_TYPE_MC_RUNH, 7, run.data(), 273) == 0, "run header");
    flush(buffer);

    double tx[2] = {100.0, 200.0}, ty[2] = {300.0, 400.0};
    double tz[2] = {0.0, 50.0}, tr[2] = {1000.0, 1200.0};
    require(write_tel_pos(buffer, 2, tx, ty, tz, tr) == 0, "telescope positions");
    flush(buffer);

    std::array<real, 273> event{};
    std::memcpy(event.data(), "EVTH", 4);
    event[1] = 8;
    event[76] = ceffic ? 4 : 8;
    event[95] = 300;
    event[96] = 600;
    require(write_tel_block(buffer, IO_TYPE_MC_EVTH, 8, event.data(), 273) == 0, "event header");
    flush(buffer);
    double ax[1] = {10}, ay[1] = {20}, aw[1] = {1.5};
    require(write_tel_offset_w(buffer, 1, 42, ax, ay, aw) == 0, "array offsets");
    flush(buffer);

    IO_ITEM_HEADER array{};
    require(begin_write_tel_array(buffer, &array, 0) == 0, "array begin");
    struct bunch full[2]{};
    full[0] = {2.5f, 100.0f, 200.0f, 0.0f, 0.0f, 12.0f, 1200000.0f,
               ceffic ? -1.0f : 0.0f};
    full[1] = {1.25f, 200.0f, 300.0f, 0.3f, 0.4f, 13.0f, 1300000.0f,
               ceffic ? -1.0f : 400.0f};
    if (compact) {
      struct compact_bunch packed[2]{};
      packed[0] = {250, 1000, 2000, 0, 0, 120, 6079,
                   static_cast<short>(ceffic ? -1 : 0)};
      packed[1] = {125, 2000, 3000, 9000, 12000, 130, 6114,
                   static_cast<short>(ceffic ? -1 : 400)};
      require(write_tel_compact_photons(buffer, 0, 0, 3.75, packed, 2, 0, nullptr) == 0,
              "compact photons");
    } else {
      require(write_tel_photons(buffer, 0, 0, 3.75, full, 2, 0, nullptr) == 0,
              "full photons");
    }
    require(end_write_tel_array(buffer, &array) == 0, "array end");
    flush(buffer);

    // A second array reuse uses the separate top-level telescope blocks.
    require(write_tel_array_head(buffer, &array, 0) == 0, "split array begin");
    flush(buffer);
    struct bunch3d spatial[1]{};
    spatial[0] = {0.75f, 100.0f, 0.0f, 50.0f, 0.0f, 0.0f, -1.0f,
                  20.0f, 100000.0f, ceffic ? -1.0f : 350.0f};
    require(write_tel_photons3d(buffer, 0, 1, 0.75, spatial, 1, 0, nullptr) == 0,
            "3D photons");
    flush(buffer);
    require(write_tel_array_end(buffer, &array, 0) == 0, "split array end");
    flush(buffer);
    std::array<real, 273> event_end{};
    std::memcpy(event_end.data(), "EVTE", 4);
    event_end[1] = 8;
    require(write_tel_block(buffer, IO_TYPE_MC_EVTE, 8, event_end.data(), 273) == 0,
            "event end");
    flush(buffer);
    std::array<real, 3> run_end{};
    std::memcpy(run_end.data(), "RUNE", 4);
    run_end[1] = 7;
    require(write_tel_block(buffer, IO_TYPE_MC_RUNE, 7, run_end.data(), 3) == 0,
            "run end");
    flush(buffer);
    std::fclose(file);
    free_io_buffer(buffer);
  }
  ~Fixture() { std::filesystem::remove(path); }
};

}  // namespace

int main() {
  try {
    Fixture input;
    obdeect::EventioPhotonReader reader(input.path.string());
    std::array<obdeect::OpticalPhoton, 1> batch{};
    auto first = reader.read(batch);
    require(first.count == 1 && !first.eof, "first chunk");
    require(first.context.run_id == 7 && first.context.event_id == 8 &&
                first.context.array_id == 0 && first.context.telescope_id == 0,
            "run and telescope IDs");
    require(first.context.array_reuse_weight == 1.5 &&
                first.context.telescope_position_m.x == 1.0,
            "array and telescope metadata");
    require(batch[0].ray.position_m.x == 1.0 && batch[0].ray.position_m.y == 2.0 &&
                batch[0].ray.direction.z == -1.0 && batch[0].weight == 2.5 &&
                batch[0].wavelength_nm == 0.0 && batch[0].emission_height_m == 12000.0,
            "full bunch conversion and unknown wavelength");
    require(reader.run_info().wavelength_lower_nm == 300 &&
                reader.run_info().wavelength_upper_nm == 600 &&
                reader.run_info().atmext && reader.run_info().reuse_weight_known &&
                reader.run_info().observation_altitude_m == 2000,
            "run options and wavelength band");
    auto second = reader.read(batch);
    require(second.count == 1 && !second.eof && batch[0].weight == 1.25 &&
                std::abs(batch[0].ray.direction.z + std::sqrt(0.75)) < 1e-6,
            "second full bunch and downward direction");
    auto third = reader.read(batch);
    require(third.count == 1 && third.eof && third.context.telescope_id == 1 &&
                batch[0].ray.position_m.z == 0.5 && batch[0].emission_distance_m == 1000.0 &&
                batch[0].wavelength_nm == 350.0,
            "split 3D bunch and EOF");
    require(reader.read(batch).count == 0, "stable EOF");

    Fixture compact_input(false, true);
    obdeect::EventioPhotonReader compact_reader(compact_input.path.string());
    const auto compact = compact_reader.read(batch);
    require(compact.count == 1 && !compact.eof && batch[0].weight == 2.5 &&
                batch[0].ray.position_m.x == 1.0 && batch[0].time_ns == 12.0 &&
                std::abs(batch[0].emission_height_m - 12000.0) < 20.0,
            "compact decoding and quantisation");

    Fixture preselected(true);
    obdeect::EventioPhotonReader forbidden(preselected.path.string());
    bool rejected = false;
    try { (void) forbidden.read(batch); }
    catch (const std::runtime_error&) { rejected = true; }
    require(rejected, "preselected CEFFIC input must fail closed");
    std::cout << "CORSIKA EventIO reader integration tests passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
