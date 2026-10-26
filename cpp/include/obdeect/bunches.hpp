#pragma once

#include "obdeect/bounds.hpp"
#include "obdeect/source.hpp"
#include "obdeect/trace.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace obdeect {

// Photon bunches of an air shower (SoA, n entries). A bunch is `size` photons
// that share position, direction, time, emission height and wavelength. Lengths
// are in input units (BunchFrame::length_unit_m); the photons move towards
// -z in the input frame: cz = -sqrt(max(0, 1 - cx^2 - cy^2)).
struct BunchBlock {
  const float *x{}, *y{}, *cx{}, *cy{}, *time_ns{}, *zem{}, *size{}, *wavelength_nm{};
  const std::uint64_t *bunch_id{}; // Unique per bunch for one seed.
  std::size_t n{};
};

// Input frame to telescope frame: p = R (u x, u y, u plane_z) + offset and
// d = R (cx, cy, cz). `rotation` is row-major and must be orthonormal.
struct BunchFrame {
  double rotation[9]{1, 0, 0, 0, 1, 0, 0, 0, 1};
  double offset_m[3]{};
  double length_unit_m{1.0}; // Metres per input length unit.
  double plane_z{};          // Input units.
};

struct BunchOptions {
  // Acceptance outside the telescope as a function of wavelength [nm] (x) and
  // emission height [m] (y), values in [0, 1]; null means 1.
  const Map2D *acceptance{};
  const LinearDensity *wavelength_density{}; // Over nm, for bunches of wavelength 0.
  bool roulette{true};                       // Survivors get weight 1 (else the weight ratio).
  bool classify_losses{true};
};

// Output columns, with room for the rows the allocator was asked for. Vectors
// are (rows, 3) row-major.
struct BunchResult {
  double *position_m{}, *direction{}, *time_ns{}, *weight{}, *wavelength_nm{};
  std::uint64_t *bunch_index{}; // Input row of the bunch.
  std::uint32_t *sub_photon{};  // Photon number within the bunch.
};

namespace detail {

inline void validate_frame(const BunchFrame &f) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      double s = 0.0;
      for (int k = 0; k < 3; ++k)
        s += f.rotation[3 * i + k] * f.rotation[3 * j + k];
      if (!(std::abs(s - (i == j ? 1.0 : 0.0)) < 1e-12))
        throw std::invalid_argument("bunch frame: rotation is not orthonormal");
    }
  if (!(f.length_unit_m > 0.0) || !std::isfinite(f.length_unit_m) || !std::isfinite(f.plane_z) ||
      !std::isfinite(f.offset_m[0]) || !std::isfinite(f.offset_m[1]) ||
      !std::isfinite(f.offset_m[2]))
    throw std::invalid_argument("bunch frame: invalid length unit, plane or offset");
}

inline void validate_acceptance(const Map2D &a) {
  if (a.x_m.empty() || a.y_m.empty() || a.values.size() != a.x_m.size() * a.y_m.size())
    throw std::invalid_argument("bunch acceptance: inconsistent map");
  for (const double v : a.values)
    if (!(v >= 0.0 && v <= 1.0))
      throw std::invalid_argument("bunch acceptance: values must lie in [0, 1]");
  for (const auto *axis : {&a.x_m, &a.y_m})
    for (std::size_t i = 0; i < axis->size(); ++i)
      if (!std::isfinite((*axis)[i]) || (i > 0 && !((*axis)[i] > (*axis)[i - 1])))
        throw std::invalid_argument("bunch acceptance: axes must be finite and increase");
}

// Largest bunch size: the sub-photon index is 32-bit.
inline constexpr double kMaxBunchSize = 4.0e9;

// Photons in a bunch of this size.
[[nodiscard]] inline std::size_t sub_photons(float size) {
  return static_cast<std::size_t>(std::ceil(static_cast<double>(size)));
}

// The work of one trace_bunches call.
class BunchTask {
public:
  BunchTask(const OpticalModel &model, const BunchBlock &in, const BunchFrame &frame,
            const BunchOptions &options, std::uint64_t seed)
      : model_(model), in_(in), frame_(frame), options_(options), seed_(seed) {}

  // Rows that chunk c can write at most (ceil of the sizes); sets `bad` for an
  // invalid bunch.
  std::size_t rows_max(std::size_t begin, std::size_t end, bool &bad) const {
    std::size_t rows = 0;
    for (std::size_t i = begin; i < end; ++i) {
      const double size = in_.size[i];
      const double lambda = in_.wavelength_nm[i];
      if (!(size > 0.0 && size <= kMaxBunchSize) || !std::isfinite(lambda) || lambda < 0.0 ||
          (lambda == 0.0 && !options_.wavelength_density) || !std::isfinite(in_.x[i]) ||
          !std::isfinite(in_.y[i]) || !std::isfinite(in_.cx[i]) || !std::isfinite(in_.cy[i]) ||
          !(static_cast<double>(in_.cx[i]) * in_.cx[i] +
                static_cast<double>(in_.cy[i]) * in_.cy[i] <=
            1.0 + 1e-6) ||
          !std::isfinite(in_.time_ns[i]) || !std::isfinite(in_.zem[i]))
        bad = true;
      else
        rows += sub_photons(in_.size[i]);
    }
    return rows;
  }

  // Trace the bunches [begin, end), writing rows from `row` on; returns the
  // number of rows written.
  std::size_t run(std::size_t begin, std::size_t end, const BunchResult &out,
                  std::size_t row) const {
    const std::size_t first = row;
    const double u = frame_.length_unit_m;
    const double *r = frame_.rotation;
    // Position and direction in the telescope frame, computed for the first
    // photon of a bunch that survives the thinning.
    const auto transform = [&](std::size_t i) {
      const double x = u * in_.x[i], y = u * in_.y[i], z = u * frame_.plane_z;
      const double cx = in_.cx[i], cy = in_.cy[i];
      const double cz = -std::sqrt(std::max(0.0, 1.0 - cx * cx - cy * cy));
      return Ray{{r[0] * x + r[1] * y + r[2] * z + frame_.offset_m[0],
                  r[3] * x + r[4] * y + r[5] * z + frame_.offset_m[1],
                  r[6] * x + r[7] * y + r[8] * z + frame_.offset_m[2]},
                 {r[0] * cx + r[1] * cy + r[2] * cz, r[3] * cx + r[4] * cy + r[5] * cz,
                  r[6] * cx + r[7] * cy + r[8] * cz}};
    };
    // Bound and acceptance of the last wavelength (and height): showers have
    // many bunches of one wavelength.
    double cached_lambda = -1.0, cached_zem = 0.0, cached_bound = 0.0, cached_acceptance = 1.0;
    const auto response = [&](double lambda, double zem_m) {
      if (lambda != cached_lambda) {
        cached_lambda = lambda;
        cached_bound = response_bound(model_, lambda);
        cached_zem = std::numeric_limits<double>::quiet_NaN();
      }
      if (zem_m != cached_zem) {
        cached_zem = zem_m;
        cached_acceptance = options_.acceptance ? options_.acceptance->at(lambda, zem_m) : 1.0;
      }
    };
    for (std::size_t i = begin; i < end; ++i) {
      const double zem_m = u * in_.zem[i];
      const double size = in_.size[i], time = in_.time_ns[i];
      const std::uint64_t id = in_.bunch_id[i];
      const std::size_t n_sub = sub_photons(in_.size[i]);
      const bool fixed = in_.wavelength_nm[i] > 0.0;
      double lambda = in_.wavelength_nm[i];
      if (fixed)
        response(lambda, zem_m);
      bool placed = false;
      Ray ray;
      for (std::size_t k = 0; k < n_sub; ++k) {
        const auto sub = static_cast<std::uint32_t>(k);
        const double frac = static_cast<double>(k + 1) <= size ? 1.0 : size - std::floor(size);
        if (!fixed) {
          lambda =
              options_.wavelength_density->sample(uniform4(seed_, id, Stream::wavelength, sub)[0]);
          response(lambda, zem_m);
        }
        const double bound = cached_bound;
        const double keep = frac * cached_acceptance * bound;
        const auto draw = uniform4(seed_, id, Stream::thinning, sub);
        if (!(draw[0] < keep))
          continue;
        if (!placed) {
          ray = transform(i);
          placed = true;
        }
        const auto t = trace_photon(model_, ray, lambda, time, 1.0, id, seed_, true, false, false,
                                    options_.classify_losses, sub);
        if (t.status != Status::detected)
          continue;
        const double ratio = std::min(1.0, t.weight / bound);
        if (options_.roulette && !(draw[1] < ratio))
          continue;
        out.position_m[3 * row] = t.position.x;
        out.position_m[3 * row + 1] = t.position.y;
        out.position_m[3 * row + 2] = t.position.z;
        out.direction[3 * row] = t.direction.x;
        out.direction[3 * row + 1] = t.direction.y;
        out.direction[3 * row + 2] = t.direction.z;
        out.time_ns[row] = t.time_ns;
        out.weight[row] = options_.roulette ? 1.0 : ratio;
        out.wavelength_nm[row] = lambda;
        out.bunch_index[row] = i;
        out.sub_photon[row] = sub;
        ++row;
      }
    }
    return row - first;
  }

private:
  const OpticalModel &model_;
  BunchBlock in_;
  BunchFrame frame_;
  BunchOptions options_;
  std::uint64_t seed_;
};

} // namespace detail

// Trace photon bunches through the optical model. Each photon of a bunch is
// kept with probability size-fraction x acceptance x response bound, traced
// with its own random numbers (seed, bunch_id, sub-photon), and, with
// `roulette`, kept again with probability weight / bound, so that the expected
// weight is the weight a full trace would give. The results do not depend on
// the thread count, the chunking or the order of the bunches; bunch ids must be
// unique per seed.
//
// Invalid input (a non-orthonormal frame, acceptance outside [0, 1], a bunch of
// wavelength 0 without a density, bad sizes) throws std::invalid_argument
// before any tracing; so does a transmission without a known maximum
// (ModelError).
//
// `allocate(rows)` is called once, after validation, with the largest number
// of rows the call can write, and returns the output columns for that many
// rows. Returns the rows written; they come first, in input order, bunch by
// bunch and sub-photon by sub-photon.
template <class Allocate>
std::size_t trace_bunches(const OpticalModel &model, const BunchBlock &in, const BunchFrame &frame,
                          const BunchOptions &options, std::uint64_t seed, ThreadPool &pool,
                          Allocate &&allocate) {
  detail::validate_frame(frame);
  if (options.acceptance)
    detail::validate_acceptance(*options.acceptance);
  if (options.wavelength_density && !(options.wavelength_density->lower() > 0.0))
    throw std::invalid_argument("bunch wavelength density: nodes must start above 0 nm");
  validate_bound_support(model); // Throws ModelError for an unsupported transmission.
  const detail::BunchTask task(model, in, frame, options, seed);

  // First pass: the most rows each chunk can write (and the input check).
  const std::size_t chunks = (in.n + kChunkSize - 1) / kChunkSize;
  std::vector<std::size_t> start(chunks + 1, 0), found(chunks, 0);
  std::vector<char> bad(chunks, 0);
  pool.run(chunks, [&](std::size_t c) {
    bool b = false;
    start[c + 1] = task.rows_max(c * kChunkSize, std::min(in.n, (c + 1) * kChunkSize), b);
    bad[c] = b;
  });
  if (std::find(bad.begin(), bad.end(), char{1}) != bad.end())
    throw std::invalid_argument("bunches: invalid size, position, direction, time or wavelength "
                                "(wavelength 0 needs a wavelength density)");
  for (std::size_t c = 0; c < chunks; ++c)
    start[c + 1] += start[c];

  const BunchResult out = allocate(start.back());

  // Second pass: chunk c writes from start[c] on.
  pool.run(chunks, [&](std::size_t c) {
    found[c] = task.run(c * kChunkSize, std::min(in.n, (c + 1) * kChunkSize), out, start[c]);
  });

  std::vector<detail::Column> columns;
  detail::add_column(columns, out.position_m, 3);
  detail::add_column(columns, out.direction, 3);
  detail::add_column(columns, out.time_ns, 1);
  detail::add_column(columns, out.weight, 1);
  detail::add_column(columns, out.wavelength_nm, 1);
  detail::add_column(columns, out.bunch_index, 1);
  detail::add_column(columns, out.sub_photon, 1);
  start.pop_back();
  return detail::compact_chunks(columns, start, found);
}

} // namespace obdeect
