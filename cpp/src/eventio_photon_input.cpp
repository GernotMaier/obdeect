#include "obdeect/eventio_photon_input.hpp"

extern "C" {
#include "initial.h"
#include "io_basic.h"
#include "mc_tel.h"
}
#include "fileopen.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace obdeect {
namespace {

constexpr int kMaxTelescopes = 4096;
constexpr int kMaxArrays = 4096;

[[nodiscard]] std::uint64_t integer_id(double value, const char* name) {
  if (!std::isfinite(value) || value < 0 || value > 9007199254740991.0 ||
      std::floor(value) != value)
    throw std::runtime_error(std::string("invalid EventIO ") + name);
  return static_cast<std::uint64_t>(value);
}

[[nodiscard]] Vec3 downward_direction(double cx, double cy) {
  const double transverse = cx * cx + cy * cy;
  if (!std::isfinite(transverse) || transverse > 1.0 + 2e-4)
    throw std::runtime_error("invalid EventIO 2D direction cosines");
  // Compact direction cosines are quantised independently. The small excess
  // at the horizon is a format effect; all other invalid vectors fail.
  const double scale = transverse > 1.0 ? 1.0 / std::sqrt(transverse) : 1.0;
  const Vec3 direction{cx * scale, cy * scale,
                       -std::sqrt(std::max(0.0, 1.0 - transverse * scale * scale))};
  return direction;
}

}  // namespace

struct EventioPhotonReader::Impl {
  explicit Impl(const std::string& input_path, EventioInputLimits requested_limits)
      : path(input_path), limits(requested_limits), telescope_x(kMaxTelescopes),
        telescope_y(kMaxTelescopes), telescope_z(kMaxTelescopes), telescope_r(kMaxTelescopes),
        array_x(kMaxArrays), array_y(kMaxArrays), array_weight(kMaxArrays) {
    if (limits.max_block_bytes < 32 ||
        limits.max_block_bytes > static_cast<std::size_t>(std::numeric_limits<long>::max()) ||
        limits.max_bunches_per_telescope == 0 ||
        limits.max_bunches_per_telescope > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::invalid_argument("invalid EventIO input limits");
    buffer = allocate_io_buffer(32768);
    if (!buffer) throw std::runtime_error("cannot allocate EventIO input buffer");
    buffer->max_length = static_cast<long>(limits.max_block_bytes);
    file = fileopen(path.c_str(), READ_BINARY);
    if (!file) {
      free_io_buffer(buffer);
      buffer = nullptr;
      throw std::runtime_error("cannot open CORSIKA EventIO file: " + path);
    }
    buffer->input_file = file;
  }

  ~Impl() {
    if (file) fileclose(file);
    if (buffer) free_io_buffer(buffer);
  }

  [[noreturn]] void fail(const std::string& reason) const {
    throw std::runtime_error("CORSIKA EventIO " + path + ": " + reason);
  }

  void read_header(int type) {
    std::array<real, 273> values{};
    if (read_tel_block(buffer, type, values.data(), static_cast<int>(values.size())) != 0)
      fail("cannot decode run/event header");
    if (type == IO_TYPE_MC_RUNH) {
      if (have_run || have_event) fail("run header before prior run ended");
      run_id = integer_id(values[1], "run ID");
      have_run = true;
      have_event = false;
      telescope_count = 0;
      array_count = 0;
      info = {};
      const auto levels = integer_id(values[4], "observation-level count");
      if (levels == 0 || levels > 10) fail("invalid observation-level count");
      info.observation_altitude_m = 0.01 * values[4 + levels];
      if (!std::isfinite(info.observation_altitude_m) || info.observation_altitude_m <= 0)
        fail("invalid observation altitude");
    } else {
      if (!have_run) fail("event header before run header");
      if (have_event) fail("event header before prior event ended");
      event_id = integer_id(values[1], "event ID");
      have_event = true;
      array_count = 0;
      info.reuse_weight_known = false;
      if (!std::isfinite(options_value) || options_value < 0 || options_value > 65535 ||
          std::floor(options_value) != options_value)
        fail("invalid IACT options in event header");
      const auto options = static_cast<unsigned>(options_value) & 0x3ffU;
      info.ceffic = (options & 0x04U) != 0;
      info.atmext = (options & 0x08U) != 0;
      info.refraction = (options & 0x10U) != 0;
      info.curved = (options & 0x40U) != 0;
      info.wavelength_lower_nm = values[95];
      info.wavelength_upper_nm = values[96];
      if (!std::isfinite(info.wavelength_lower_nm) ||
          !std::isfinite(info.wavelength_upper_nm) ||
          info.wavelength_lower_nm <= 0 ||
          info.wavelength_lower_nm >= info.wavelength_upper_nm)
        fail("invalid Cherenkov wavelength band in event header");
      info.has_event_header = true;
    }
  }

  void read_footer(int type) {
    std::array<real, 273> values{};
    if (read_tel_block(buffer, type, values.data(), static_cast<int>(values.size())) != 0)
      fail("cannot decode run/event footer");
    if (type == IO_TYPE_MC_EVTE) {
      if (!have_event || integer_id(values[1], "event end ID") != event_id || split_array)
        fail("mismatched event end");
      have_event = false;
    } else {
      if (!have_run || have_event || integer_id(values[1], "run end ID") != run_id)
        fail("mismatched run end");
      have_run = false;
    }
  }

  void read_positions() {
    if (read_tel_pos(buffer, kMaxTelescopes, &telescope_count, telescope_x.data(),
                     telescope_y.data(), telescope_z.data(), telescope_r.data()) != 0)
      fail("cannot decode telescope positions or telescope limit exceeded");
    for (int i = 0; i < telescope_count; ++i)
      if (!std::isfinite(telescope_x[i]) || !std::isfinite(telescope_y[i]) ||
          !std::isfinite(telescope_z[i]) || !std::isfinite(telescope_r[i]) ||
          telescope_r[i] < 0)
        fail("invalid telescope position or radius");
  }

  void read_offsets(unsigned version) {
    if (read_tel_offset_w(buffer, kMaxArrays, &array_count, &info.array_time_offset_ns,
                          array_x.data(), array_y.data(), array_weight.data()) != 0)
      fail("cannot decode array offsets or array limit exceeded");
    info.reuse_weight_known = version == 1;
    if (!std::isfinite(info.array_time_offset_ns)) fail("invalid array time offset");
    for (int i = 0; i < array_count; ++i)
      if (!std::isfinite(array_x[i]) || !std::isfinite(array_y[i]) ||
          (info.reuse_weight_known && (!std::isfinite(array_weight[i]) || array_weight[i] < 0)))
        fail("invalid array offset or weight");
  }

  void load_photons(int type, int expected_array) {
    if (!have_event || telescope_count == 0 || array_count == 0)
      fail("photon block lacks event, telescope or array metadata");
    int count = 0, array = -1, telescope = -1;
    double total = 0;
    const int probe = type == IO_TYPE_MC_PHOTONS
                          ? read_tel_photons(buffer, 0, &array, &telescope, &total, nullptr, &count)
                          : read_tel_photons3d(buffer, 0, &array, &telescope, &total, nullptr, &count);
    if (probe != -10) fail("cannot inspect photon bunch count");
    if (count < 0 || static_cast<std::size_t>(count) > limits.max_bunches_per_telescope)
      fail("photon bunch count exceeds limit");
    if (array < 0 || array >= array_count || telescope < 0 || telescope >= telescope_count ||
        array != expected_array)
      fail("invalid photon array/telescope ID");
    if (!std::isfinite(total) || total < 0) fail("invalid telescope photon total");
    bunches.clear();
    bunches3d.clear();
    const auto allocation = static_cast<std::size_t>(std::max(1, count));
    if (type == IO_TYPE_MC_PHOTONS) {
      bunches.resize(allocation);
      if (read_tel_photons(buffer, count, &array, &telescope, &total, bunches.data(), &count) != 0)
        fail("cannot decode 2D photon bunches");
    } else {
      bunches3d.resize(allocation);
      if (read_tel_photons3d(buffer, count, &array, &telescope, &total, bunches3d.data(), &count) != 0)
        fail("cannot decode 3D photon bunches");
    }
    context = {run_id, event_id, static_cast<std::uint64_t>(array),
               static_cast<std::uint64_t>(telescope),
               {0.01 * telescope_x[telescope], 0.01 * telescope_y[telescope],
                0.01 * telescope_z[telescope]},
               info.reuse_weight_known ? array_weight[array] : 1.0};
    if (!valid_context(context)) fail("invalid photon batch context");
    current_type = type;
    bunch_count = static_cast<std::size_t>(count);
    bunch_index = 0;
    prepared = bunch_count > 0;
  }

  [[nodiscard]] OpticalPhoton convert(std::size_t index) {
    double x, y, z, cx, cy, cz, time, emission, photons, wavelength;
    const bool is_3d = current_type == IO_TYPE_MC_PHOTONS3D;
    if (is_3d) {
      const auto& bunch = bunches3d[index];
      x = bunch.x; y = bunch.y; z = bunch.z;
      cx = bunch.cx; cy = bunch.cy; cz = bunch.cz;
      time = bunch.ctime; emission = bunch.dist;
      photons = bunch.photons; wavelength = bunch.lambda;
    } else {
      const auto& bunch = bunches[index];
      x = bunch.x; y = bunch.y; z = 0.0;
      cx = bunch.cx; cy = bunch.cy; cz = 0.0;
      time = bunch.ctime; emission = bunch.zem;
      photons = bunch.photons; wavelength = bunch.lambda;
    }
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) ||
        !std::isfinite(time) || !std::isfinite(emission) || emission < 0 ||
        !std::isfinite(photons) || photons < 0 || !std::isfinite(wavelength))
      fail("invalid values in photon bunch");
    if (wavelength < 0 || info.ceffic)
      fail("CEFFIC/photoelectron bunches require a separate input mode");
    if (wavelength > 0 && (wavelength < info.wavelength_lower_nm ||
                           wavelength > info.wavelength_upper_nm))
      fail("out-of-band wavelength or unsupported source marker");
    const Vec3 direction = [&] {
      if (!is_3d) return downward_direction(cx, cy);
      const double norm_squared = cx * cx + cy * cy + cz * cz;
      if (!std::isfinite(norm_squared) || std::abs(norm_squared - 1.0) > 2e-4)
        fail("invalid 3D photon direction");
      const double inverse = 1.0 / std::sqrt(norm_squared);
      return Vec3{cx * inverse, cy * inverse, cz * inverse};
    }();
    OpticalPhoton result{{{0.01 * x, 0.01 * y, 0.01 * z}, direction},
                         next_photon_id++, wavelength, time, photons,
                         static_cast<std::uint64_t>(index),
                         is_3d ? std::numeric_limits<double>::quiet_NaN() : 0.01 * emission,
                         is_3d ? 0.01 * emission : std::numeric_limits<double>::quiet_NaN()};
    if (!valid_photon(result)) fail("invalid converted photon bunch");
    return result;
  }

  void prepare_next() {
    prepared = false;
    for (;;) {
      if (in_array) {
        const int type = next_subitem_type(buffer);
        if (type == -1) fail("invalid array subitem");
        if (type == -2) {
          if (end_read_tel_array(buffer, &array_header) != 0) fail("invalid array block end");
          in_array = false;
          continue;
        }
        if (type == IO_TYPE_MC_PHOTONS || type == IO_TYPE_MC_PHOTONS3D) {
          load_photons(type, active_array);
          if (prepared) return;
        } else if (skip_subitem(buffer) != 0) {
          fail("cannot skip array subitem");
        }
        continue;
      }
      IO_ITEM_HEADER header{};
      const int found = find_io_block(buffer, &header);
      if (found == -2) {
        if (split_array || in_array || have_event || have_run)
          fail("truncated CORSIKA EventIO run or event at end of file");
        eof = true;
        return;
      }
      if (found != 0) fail("cannot find next EventIO block");
      if (header.length > limits.max_block_bytes) fail("EventIO block exceeds byte limit");
      const int type = static_cast<int>(header.type);
      const bool relevant = type == IO_TYPE_MC_RUNH || type == IO_TYPE_MC_EVTH ||
                            type == IO_TYPE_MC_TELPOS || type == IO_TYPE_MC_TELOFF ||
                            type == IO_TYPE_MC_TELARRAY || type == IO_TYPE_MC_TELARRAY_HEAD ||
                            type == IO_TYPE_MC_TELARRAY_END || type == IO_TYPE_MC_PHOTONS ||
                            type == IO_TYPE_MC_PHOTONS3D || type == IO_TYPE_MC_EVTE ||
                            type == IO_TYPE_MC_RUNE;
      if (!relevant) {
        if (skip_io_block(buffer, &header) != 0) fail("cannot skip EventIO block");
        continue;
      }
      if (read_io_block(buffer, &header) != 0) fail("cannot read EventIO block");
      switch (type) {
        case IO_TYPE_MC_RUNH: read_header(type); break;
        case IO_TYPE_MC_EVTH: read_header(type); break;
        case IO_TYPE_MC_TELPOS: read_positions(); break;
        case IO_TYPE_MC_TELOFF: read_offsets(header.version); break;
        case IO_TYPE_MC_TELARRAY:
          if (split_array) fail("nested array inside split array");
          if (begin_read_tel_array(buffer, &array_header, &active_array) != 0)
            fail("cannot begin telescope array");
          if (active_array < 0 || active_array >= array_count) fail("invalid array block ID");
          in_array = true;
          break;
        case IO_TYPE_MC_TELARRAY_HEAD:
          if (split_array) fail("nested split array headers");
          split_array = true;
          active_array = static_cast<int>(header.ident);
          if (active_array < 0 || active_array >= array_count) fail("invalid split array ID");
          break;
        case IO_TYPE_MC_TELARRAY_END:
          if (!split_array || header.ident != active_array) fail("mismatched split array end");
          split_array = false;
          break;
        case IO_TYPE_MC_PHOTONS:
        case IO_TYPE_MC_PHOTONS3D:
          // A top-level 999/999 item contains ground particles, not photons.
          if (split_array) {
            load_photons(type, active_array);
            if (prepared) return;
          }
          break;
        case IO_TYPE_MC_EVTE: read_footer(type); break;
        case IO_TYPE_MC_RUNE: read_footer(type); break;
        default: break;
      }
    }
  }

  std::string path;
  EventioInputLimits limits;
  std::FILE* file{};
  IO_BUFFER* buffer{};
  EventioRunInfo info{};
  EventioRunInfo last_batch_info{};
  std::uint64_t run_id{}, event_id{}, next_photon_id{};
  bool have_run{}, have_event{}, eof{}, prepared{}, in_array{}, split_array{};
  int telescope_count{}, array_count{}, active_array{-1}, current_type{};
  IO_ITEM_HEADER array_header{};
  PhotonBatchContext context{};
  std::vector<double> telescope_x, telescope_y, telescope_z, telescope_r;
  std::vector<double> array_x, array_y, array_weight;
  std::vector<struct bunch> bunches;
  std::vector<struct bunch3d> bunches3d;
  std::size_t bunch_count{}, bunch_index{};
};

EventioPhotonReader::EventioPhotonReader(const std::string& path, EventioInputLimits limits)
    : impl_(std::make_unique<Impl>(path, limits)) {}
EventioPhotonReader::~EventioPhotonReader() = default;
EventioPhotonReader::EventioPhotonReader(EventioPhotonReader&&) noexcept = default;
EventioPhotonReader& EventioPhotonReader::operator=(EventioPhotonReader&&) noexcept = default;

PhotonReadResult EventioPhotonReader::read(std::span<OpticalPhoton> destination) {
  if (destination.empty()) throw std::invalid_argument("photon batch capacity must be positive");
  for (;;) {
    if (!impl_->prepared && !impl_->eof) impl_->prepare_next();
    if (impl_->eof) return {{}, 0, true};
    const PhotonBatchContext context = impl_->context;
    const EventioRunInfo batch_info = impl_->info;
    std::size_t written = 0;
    while (written < destination.size() && impl_->bunch_index < impl_->bunch_count) {
      const auto index = impl_->bunch_index++;
      const auto wavelength = impl_->current_type == IO_TYPE_MC_PHOTONS
                                  ? impl_->bunches[index].lambda : impl_->bunches3d[index].lambda;
      if (wavelength >= 9000.0) continue;  // IACTEXT emitter marker, not light.
      destination[written++] = impl_->convert(index);
    }
    if (impl_->bunch_index == impl_->bunch_count) {
      impl_->prepared = false;
      impl_->prepare_next();
    }
    // Marker-only blocks can produce no optical bunches.
    if (written > 0 || impl_->eof) {
      impl_->last_batch_info = batch_info;
      return {context, written, impl_->eof};
    }
  }
}

const EventioRunInfo& EventioPhotonReader::run_info() const noexcept {
  return impl_->last_batch_info;
}

}  // namespace obdeect
