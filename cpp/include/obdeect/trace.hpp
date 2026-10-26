#pragma once

#include "obdeect/model.hpp"
#include "obdeect/random.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <limits>
#include <mutex>
#include <optional>
#include <string_view>
#include <thread>
#include <vector>

namespace obdeect {

enum class Status : std::uint8_t {
  detected,             // Reached the focal surface.
  missed_primary,       // Never met the primary mirror.
  missed_secondary,     // Left the primary without meeting the secondary.
  missed_focal_surface, // Never met the focal surface.
  blocked,              // Absorbed by an obscurer (terminal_surface names it).
  invalid_input,        // Non-finite input, wavelength <= 0 or non-unit direction.
  numerical_failure,    // An intersection did not converge.
};
inline constexpr std::size_t kStatusCount = 7;

[[nodiscard]] constexpr std::string_view to_string(Status status) {
  constexpr std::string_view names[kStatusCount] = {
      "detected", "missed_primary", "missed_secondary", "missed_focal_surface",
      "blocked",  "invalid_input",  "numerical_failure"};
  return names[static_cast<std::size_t>(status)];
}

inline constexpr std::uint32_t kNoSurface = 0xFFFFFFFFU;

// Result of one photon. Mirror hits not reached keep NaN coordinates and
// kNoSurface; `weight` is the input weight times the applied responses.
struct TraceResult {
  Status status{Status::invalid_input};
  std::uint32_t terminal_surface{kNoSurface};
  Vec3 position, direction;
  double path_length_m{}, time_ns{}, weight{};
  std::uint32_t primary_surface{kNoSurface}, secondary_surface{kNoSurface};
  Vec3 primary_point{nan(), nan(), nan()}, secondary_point{nan(), nan(), nan()};
  double primary_incidence_deg{nan()}, secondary_incidence_deg{nan()}, focal_incidence_deg{nan()};

  static constexpr double nan() { return std::numeric_limits<double>::quiet_NaN(); }
};

namespace detail {

constexpr double kDegree = 57.29577951308232087680;

inline double incidence_deg(const Vec3 &direction, const Vec3 &normal) {
  return std::acos(std::min(1.0, std::abs(dot(direction, normal)))) * kDegree;
}

struct MirrorHit {
  SurfaceHit hit;
  std::uint32_t facet_index{};
};

// Nearest hit on a mirror (lowest id on exact ties). `failed` reports a
// non-converged aspheric intersection.
inline std::optional<MirrorHit> intersect_mirror(const Mirror &mirror, const Ray &ray,
                                                 bool &failed) {
  failed = false;
  if (mirror.asphere) {
    const auto &a = *mirror.asphere;
    // The lowest segment id at the accepted point: accept() is the last call
    // that returned true when intersect() returns, and it saw the hit point.
    std::uint32_t segment_id = a.id;
    const auto accept = [&](const Vec3 &p) {
      if (mirror.segments.empty())
        return true;
      const auto id = mirror.segment_grid.find(mirror.segments, mirror.segment_frames, p);
      if (id)
        segment_id = *id;
      return id.has_value();
    };
    auto hit = a.intersect(ray, kMinDistanceM, accept, failed);
    if (!hit)
      return std::nullopt;
    hit->id = segment_id;
    return MirrorHit{*hit, 0};
  }
  std::optional<MirrorHit> nearest;
  mirror.grid.candidates(ray, mirror.facets.size(), [&](std::size_t i) {
    const auto &f = mirror.facets[i];
    if (const auto hit = f.intersect(ray, kMinDistanceM))
      if (!nearest || hit->t < nearest->hit.t ||
          (hit->t == nearest->hit.t && hit->id < nearest->hit.id))
        nearest = MirrorHit{*hit, static_cast<std::uint32_t>(i)};
  });
  return nearest;
}

// Rotation of v about the unit axis k by `angle` (Rodrigues).
inline Vec3 rotate(const Vec3 &v, const Vec3 &k, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  return v * c + cross(k, v) * s + k * (dot(k, v) * (1.0 - c));
}

// Reflected direction, with the mirror's random reflection error.
inline Vec3 reflect_at(const Mirror &mirror, const MirrorHit &mh, const Vec3 &direction,
                       std::uint64_t photon_id, std::uint64_t seed, Stream stream,
                       std::uint32_t sub_photon = 0) {
  const Vec3 &n = mh.hit.normal;
  if (!mirror.scatter)
    return reflect(direction, n);
  const auto &s = *mirror.scatter;
  const auto u = uniform4(seed, photon_id, stream, sub_photon);
  const double sigma = u[0] < s.fraction2 ? s.sigma2_rad : s.sigma1_rad;
  if (sigma == 0.0)
    return reflect(direction, n);
  const double radius = sigma * std::sqrt(-2.0 * std::log(u[1]));
  const double phi = 6.283185307179586476925 * u[2];
  const double a1 = radius * std::cos(phi), a2 = radius * std::sin(phi);
  if (s.surface_slopes) {
    // Graph-surface normal (-dz/dx, -dz/dy, 1): perturb each slope by half
    // the angle, so the reflected ray deviates by about the full angle.
    const Vec3 tilted{n.x / n.z + 0.5 * a1, n.y / n.z + 0.5 * a2, 1.0};
    return reflect(direction, normalized(tilted));
  }
  // Rotate the reflected ray about two axes in the mirror's tangent plane.
  const Vec3 e1 = mirror.facets.empty() ? normalized(std::abs(n.z) < 0.9 ? cross(Vec3{0, 0, 1}, n)
                                                                         : cross(Vec3{1, 0, 0}, n))
                                        : mirror.facets[mh.facet_index].u;
  const Vec3 e2 = cross(n, e1);
  return normalized(rotate(rotate(reflect(direction, n), e2, a1), e1, a2));
}

// Does the mirror's reflectivity depend on the incidence angle?
inline bool angle_matters(const Mirror &mirror) {
  return mirror.reflectivity && mirror.reflectivity->incidence_deg.size() > 1;
}

inline double mirror_response(const Mirror &mirror, const MirrorHit &mh, double wavelength_nm,
                              double incidence) {
  double r = mirror.reflectivity ? mirror.reflectivity->at(wavelength_nm, incidence) : 1.0;
  if (mirror.degradation) {
    const Vec3 &p = mh.hit.point;
    if (mirror.facets.empty()) {
      r *= mirror.degradation->at(p.x, p.y);
    } else {
      const auto &f = mirror.facets[mh.facet_index];
      r *= mirror.degradation->at(dot(p - f.centre, f.map_u), dot(p - f.centre, f.map_v));
    }
  }
  return r;
}

} // namespace detail

// Trace one photon through the optical model. With `apply_responses` false
// the weight stays at its input value (geometry-only transport). A `direct`
// photon starts inside the telescope on the leg to the focal surface (after
// the last mirror), e.g. from a flasher: no mirror and no telescope
// transmission act, only the obscurers of that leg. With `incidence_angles`
// false the incidence angles stay NaN (saving an arccosine each) except
// where a response table depends on them. With `classify_losses` false a
// photon that misses the primary is reported as missed_primary even where it
// would be blocked on the incoming leg (before missing the primary); detected
// photons are the same either way. The default keeps the full classification.
// `sub_photon` numbers the photons that share one `photon_id` (the random
// counter's block); 0 gives the random numbers of a single photon.
[[nodiscard]] inline TraceResult
trace_photon(const OpticalModel &model, Ray ray, double wavelength_nm, double time_ns,
             double weight, std::uint64_t photon_id, std::uint64_t seed,
             bool apply_responses = true, bool direct = false, bool incidence_angles = true,
             bool classify_losses = true, std::uint32_t sub_photon = 0) {
  using namespace detail;
  TraceResult out;
  out.position = ray.origin;
  out.direction = ray.direction;
  const double length = norm(ray.direction);
  if (!finite(ray.origin) || !std::isfinite(length) || std::abs(length - 1.0) > 1e-6 ||
      !(wavelength_nm > 0.0) || !std::isfinite(wavelength_nm) || !std::isfinite(time_ns) ||
      !(weight >= 0.0) || !std::isfinite(weight))
    return out;
  ray.direction = ray.direction * (1.0 / length);
  out.weight = weight;
  double path = 0.0;
  const auto finish = [&](Status status, const Vec3 &position, std::uint32_t surface) {
    out.status = status;
    out.position = position;
    out.direction = ray.direction;
    out.terminal_surface = surface;
    out.path_length_m = path;
    out.time_ns = time_ns + path * model.propagation_group_index / kSpeedOfLightMPerNs;
    if (status != Status::detected)
      out.weight = 0.0;
    return out;
  };
  // A blocked photon stops at the obscurer; it is absorbed there.
  const auto blocked_on = [&](Leg leg, double t_limit) -> std::optional<TraceResult> {
    if (const auto b = model.obscurers.nearest(ray, leg, t_limit)) {
      path += b->t;
      return finish(Status::blocked, ray.at(b->t), b->id);
    }
    return std::nullopt;
  };

  bool failed = false;
  if (!direct) {
    if (apply_responses && model.telescope_transmission)
      out.weight *= model.telescope_transmission->at(ray.direction);

    // Source to primary mirror.
    const auto m1 = intersect_mirror(model.primary, ray, failed);
    if (failed)
      return finish(Status::numerical_failure, ray.origin, kNoSurface);
    double t_m1 = m1 ? m1->hit.t : std::numeric_limits<double>::infinity();
    if (!classify_losses && !m1)
      return finish(Status::missed_primary, ray.origin, kNoSurface);
    if (model.secondary_blocks_incoming) {
      const auto back = intersect_mirror(*model.secondary, ray, failed);
      if (failed)
        return finish(Status::numerical_failure, ray.origin, kNoSurface);
      if (back && back->hit.t < t_m1) {
        if (auto b = blocked_on(Leg::incoming, back->hit.t))
          return *b;
        path += back->hit.t;
        return finish(Status::blocked, back->hit.point, back->hit.id);
      }
    }
    if (auto b = blocked_on(Leg::incoming, t_m1))
      return *b;
    if (!m1)
      return finish(Status::missed_primary, ray.origin, kNoSurface);

    path += m1->hit.t;
    out.primary_point = m1->hit.point;
    out.primary_surface = m1->hit.id;
    if (incidence_angles || angle_matters(model.primary))
      out.primary_incidence_deg = incidence_deg(ray.direction, m1->hit.normal);
    if (apply_responses)
      out.weight *=
          mirror_response(model.primary, *m1, wavelength_nm,
                          std::isnan(out.primary_incidence_deg) ? 0.0 : out.primary_incidence_deg);
    ray = {m1->hit.point, reflect_at(model.primary, *m1, ray.direction, photon_id, seed,
                                     Stream::primary_scatter, sub_photon)};

    if (model.secondary) {
      const auto m2 = intersect_mirror(*model.secondary, ray, failed);
      if (failed)
        return finish(Status::numerical_failure, ray.origin, kNoSurface);
      if (auto b = blocked_on(Leg::after_primary,
                              m2 ? m2->hit.t : std::numeric_limits<double>::infinity()))
        return *b;
      if (!m2)
        return finish(Status::missed_secondary, ray.origin, kNoSurface);
      path += m2->hit.t;
      out.secondary_point = m2->hit.point;
      out.secondary_surface = m2->hit.id;
      if (incidence_angles || angle_matters(*model.secondary))
        out.secondary_incidence_deg = incidence_deg(ray.direction, m2->hit.normal);
      if (apply_responses)
        out.weight *= mirror_response(
            *model.secondary, *m2, wavelength_nm,
            std::isnan(out.secondary_incidence_deg) ? 0.0 : out.secondary_incidence_deg);
      ray = {m2->hit.point, reflect_at(*model.secondary, *m2, ray.direction, photon_id, seed,
                                       Stream::secondary_scatter, sub_photon)};
    }
  }
  const Leg leg = model.secondary ? Leg::after_secondary : Leg::after_primary;

  // To the focal surface.
  std::optional<SurfaceHit> focal;
  const auto &fs = model.focal_surface;
  if (fs.asphere) {
    focal = fs.asphere->intersect(ray, kMinDistanceM, [](const Vec3 &) { return true; }, failed);
    if (failed)
      return finish(Status::numerical_failure, ray.origin, kNoSurface);
  } else if (ray.direction.z != 0.0) {
    const double t = (fs.plane_z - ray.origin.z) / ray.direction.z;
    if (t > kMinDistanceM)
      focal = SurfaceHit{t, ray.at(t), Vec3{0.0, 0.0, 1.0}, fs.id};
  }
  if (auto b = blocked_on(leg, focal ? focal->t : std::numeric_limits<double>::infinity()))
    return *b;
  if (!focal)
    return finish(Status::missed_focal_surface, ray.origin, kNoSurface);
  path += focal->t;
  if (incidence_angles)
    out.focal_incidence_deg = incidence_deg(ray.direction, focal->normal);
  return finish(Status::detected, focal->point, fs.id);
}

// Structure-of-arrays photon block. Vectors are (n, 3) row-major.
struct PhotonBlock {
  const double *position{}, *direction{}, *wavelength_nm{}, *time_ns{}, *weight{};
  const std::uint64_t *photon_id{};
  std::size_t size{};
};

// Output columns. A null pointer means "not requested": nothing is written.
// Vectors are (n, 3) row-major. `index` receives the input row of each output
// row (useful with detected_only).
struct ResultBlock {
  std::uint8_t *status{};
  std::uint32_t *terminal_surface{}, *primary_surface{}, *secondary_surface{};
  double *position{}, *direction{}, *path_length_m{}, *time_ns{}, *weight{};
  double *primary_point{}, *secondary_point{};
  double *primary_incidence_deg{}, *secondary_incidence_deg{}, *focal_incidence_deg{};
  std::uint64_t *index{};
};

// Chunk size of parallel_for: photons cost very different amounts, so threads
// take fixed chunks from a shared counter instead of one static share each.
inline constexpr std::size_t kChunkSize = 4096;

// Run body(begin, end) over [0, n) in chunks of kChunkSize items (the last one
// shorter) on up to `threads` threads. Each worker takes the next free chunk
// until none is left; the body must write only the items of its range.
template <class Body> void parallel_for(std::size_t n, unsigned threads, const Body &body) {
  const std::size_t chunks = (n + kChunkSize - 1) / kChunkSize;
  const std::size_t workers = std::max<std::size_t>(1, std::min<std::size_t>(threads, chunks));
  if (workers == 1) {
    body(std::size_t{0}, n);
    return;
  }
  std::atomic<std::size_t> next{0};
  const auto work = [&] {
    for (std::size_t c = next.fetch_add(1, std::memory_order_relaxed); c < chunks;
         c = next.fetch_add(1, std::memory_order_relaxed))
      body(c * kChunkSize, std::min(n, (c + 1) * kChunkSize));
  };
  std::vector<std::thread> pool;
  for (std::size_t w = 1; w < workers; ++w)
    pool.emplace_back(work);
  work(); // The calling thread works too.
  for (auto &t : pool)
    t.join();
}

// A fixed set of worker threads that runs one job at a time; the calling
// thread works too. Owned by the caller (no global pool). A task must never
// call `run` on the same pool from inside a running job (it would deadlock).
class ThreadPool {
public:
  explicit ThreadPool(unsigned threads) {
    try {
      for (unsigned i = 1; i < std::max(1U, threads); ++i)
        workers_.emplace_back([this] { loop(); });
    } catch (...) { // Joinable threads must not be destroyed.
      {
        std::lock_guard<std::mutex> lock(mutex_);
        stop_ = true;
      }
      wake_.notify_all();
      for (auto &w : workers_)
        w.join();
      throw;
    }
  }
  ~ThreadPool() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_ = true;
    }
    wake_.notify_all();
    for (auto &w : workers_)
      w.join();
  }
  ThreadPool(const ThreadPool &) = delete;
  ThreadPool &operator=(const ThreadPool &) = delete;

  [[nodiscard]] unsigned size() const { return static_cast<unsigned>(workers_.size()) + 1; }

  // Calls task(c) for every c in [0, count) and returns when all calls are done.
  // If `task` throws, the remaining calls are skipped and the first exception
  // is rethrown here once all workers have left the job.
  template <class Task> void run(std::size_t count, const Task &task) {
    if (count == 0)
      return;
    std::lock_guard<std::mutex> one_job(run_mutex_);
    Job job;
    job.count = count;
    job.call = [&task](std::size_t c) { task(c); };
    {
      std::lock_guard<std::mutex> lock(mutex_);
      job_ = &job;
      ++generation_;
    }
    wake_.notify_all();
    work(job);
    std::unique_lock<std::mutex> lock(mutex_);
    job_ = nullptr; // Workers that have not joined yet will not touch the job.
    done_.wait(lock, [&] { return job.active == 0; });
    if (job.error)
      std::rethrow_exception(job.error);
  }

private:
  struct Job {
    std::size_t count{};
    std::function<void(std::size_t)> call;
    std::atomic<std::size_t> next{0};
    std::size_t active{0}; // Workers inside work(job); guarded by mutex_.
    std::exception_ptr error;
    std::mutex error_mutex;
  };
  static void work(Job &job) {
    for (std::size_t c = job.next.fetch_add(1); c < job.count; c = job.next.fetch_add(1))
      try {
        job.call(c);
      } catch (...) {
        {
          std::lock_guard<std::mutex> lock(job.error_mutex);
          if (!job.error)
            job.error = std::current_exception();
        }
        job.next.store(job.count); // Skip the remaining chunks.
      }
  }
  void loop() {
    std::uint64_t seen = 0;
    for (;;) {
      Job *job = nullptr;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        wake_.wait(lock, [&] { return stop_ || (job_ != nullptr && generation_ != seen); });
        if (stop_)
          return;
        seen = generation_;
        job = job_;
        ++job->active;
      }
      work(*job);
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (--job->active == 0)
          done_.notify_all();
      }
    }
  }
  std::vector<std::thread> workers_;
  std::mutex mutex_, run_mutex_;
  std::condition_variable wake_, done_;
  Job *job_{nullptr};
  std::uint64_t generation_{0};
  bool stop_{false};
};

// The same on a persistent pool: the chunks of kChunkSize items are the pool's
// tasks.
template <class Body> void parallel_for(std::size_t n, ThreadPool &pool, const Body &body) {
  const std::size_t chunks = (n + kChunkSize - 1) / kChunkSize;
  pool.run(chunks, [&](std::size_t c) { body(c * kChunkSize, std::min(n, (c + 1) * kChunkSize)); });
}

namespace detail {

// An output column for compact_chunks: its data and the bytes of one row.
struct Column {
  void *data;
  std::size_t bytes;
};

// Add a column of `per_row` elements per row; a null pointer is not requested.
template <class T> void add_column(std::vector<Column> &columns, T *pointer, std::size_t per_row) {
  if (pointer)
    columns.push_back({pointer, per_row * sizeof(T)});
}

// Chunk c wrote found[c] rows from row chunk_start[c] on (starts increasing);
// move them down in order so that all rows are contiguous, and return their
// number. Moving in order never overwrites rows not yet moved.
inline std::size_t compact_chunks(const std::vector<Column> &columns,
                                  const std::vector<std::size_t> &chunk_start,
                                  const std::vector<std::size_t> &found) {
  std::size_t total = 0;
  for (std::size_t c = 0; c < found.size(); ++c) {
    if (total != chunk_start[c])
      for (const Column &column : columns) {
        auto *bytes = static_cast<char *>(column.data);
        std::memmove(bytes + total * column.bytes, bytes + chunk_start[c] * column.bytes,
                     found[c] * column.bytes);
      }
    total += found[c];
  }
  return total;
}

// One block to trace: the per-chunk work and the final compaction, shared by
// every way of running the chunks (threads, pool, several blocks in one job).
class BlockTask {
public:
  BlockTask(const OpticalModel &model, const PhotonBlock &in, const ResultBlock &out,
            std::uint64_t seed, bool apply_responses, bool direct, bool detected_only,
            bool classify_losses)
      : model_(model), in_(in), out_(out), seed_(seed), apply_responses_(apply_responses),
        direct_(direct), detected_only_(detected_only), classify_losses_(classify_losses),
        wants_angles_(out.primary_incidence_deg || out.secondary_incidence_deg ||
                      out.focal_incidence_deg),
        found_(detected_only ? (in.size + kChunkSize - 1) / kChunkSize : 0, 0) {}

  [[nodiscard]] std::size_t size() const { return in_.size; }

  // Trace the photons [begin, end); the range must not cross a kChunkSize
  // boundary unless it starts at a chunk boundary (as parallel_for gives).
  void run(std::size_t begin, std::size_t end) {
    if (!detected_only_) {
      for (std::size_t i = begin; i < end; ++i)
        write_row(i, i, trace_one(i));
      return;
    }
    // Each chunk writes its detected rows, in input order, from the chunk's
    // first row on and records their number.
    for (std::size_t first = begin; first < end; first += kChunkSize) {
      std::size_t row = first;
      for (std::size_t i = first; i < std::min(end, first + kChunkSize); ++i) {
        const auto r = trace_one(i);
        if (r.status == Status::detected)
          write_row(row++, i, r);
      }
      found_[first / kChunkSize] = row - first;
    }
  }

  // Number of rows written; with detected_only the chunks are moved down.
  std::size_t finish() const {
    if (!detected_only_)
      return in_.size;
    std::vector<Column> columns;
    add_column(columns, out_.status, 1);
    add_column(columns, out_.terminal_surface, 1);
    add_column(columns, out_.primary_surface, 1);
    add_column(columns, out_.secondary_surface, 1);
    add_column(columns, out_.position, 3);
    add_column(columns, out_.direction, 3);
    add_column(columns, out_.path_length_m, 1);
    add_column(columns, out_.time_ns, 1);
    add_column(columns, out_.weight, 1);
    add_column(columns, out_.primary_point, 3);
    add_column(columns, out_.secondary_point, 3);
    add_column(columns, out_.primary_incidence_deg, 1);
    add_column(columns, out_.secondary_incidence_deg, 1);
    add_column(columns, out_.focal_incidence_deg, 1);
    add_column(columns, out_.index, 1);
    std::vector<std::size_t> chunk_start(found_.size());
    for (std::size_t c = 0; c < found_.size(); ++c)
      chunk_start[c] = c * kChunkSize;
    return compact_chunks(columns, chunk_start, found_);
  }

private:
  [[nodiscard]] TraceResult trace_one(std::size_t i) const {
    const Ray ray{{in_.position[3 * i], in_.position[3 * i + 1], in_.position[3 * i + 2]},
                  {in_.direction[3 * i], in_.direction[3 * i + 1], in_.direction[3 * i + 2]}};
    return trace_photon(model_, ray, in_.wavelength_nm[i], in_.time_ns[i], in_.weight[i],
                        in_.photon_id[i], seed_, apply_responses_, direct_, wants_angles_,
                        classify_losses_);
  }

  void write_row(std::size_t row, std::size_t input, const TraceResult &r) const {
    const auto put = [row](double *target, const Vec3 &v) {
      if (!target)
        return;
      target[3 * row] = v.x;
      target[3 * row + 1] = v.y;
      target[3 * row + 2] = v.z;
    };
    const auto set = [row](auto *target, auto value) {
      if (target)
        target[row] = value;
    };
    set(out_.status, static_cast<std::uint8_t>(r.status));
    set(out_.terminal_surface, r.terminal_surface);
    set(out_.primary_surface, r.primary_surface);
    set(out_.secondary_surface, r.secondary_surface);
    put(out_.position, r.position);
    put(out_.direction, r.direction);
    put(out_.primary_point, r.primary_point);
    put(out_.secondary_point, r.secondary_point);
    set(out_.path_length_m, r.path_length_m);
    set(out_.time_ns, r.time_ns);
    set(out_.weight, r.weight);
    set(out_.primary_incidence_deg, r.primary_incidence_deg);
    set(out_.secondary_incidence_deg, r.secondary_incidence_deg);
    set(out_.focal_incidence_deg, r.focal_incidence_deg);
    set(out_.index, static_cast<std::uint64_t>(input));
  }

  const OpticalModel &model_;
  PhotonBlock in_;
  ResultBlock out_;
  std::uint64_t seed_;
  bool apply_responses_, direct_, detected_only_, classify_losses_, wants_angles_;
  std::vector<std::size_t> found_; // Detected rows per chunk (detected_only).
};

} // namespace detail

// Trace a block on `threads` threads. Every photon writes only its own row,
// and its random numbers depend only on (seed, photon_id), so results do not
// depend on the thread count or the order of photons. Returns the number of
// rows written: in.size, or with `detected_only` the number of detected
// photons, whose rows come first, in input order (the rest of the output is
// left unwritten).
inline std::size_t trace_block(const OpticalModel &model, const PhotonBlock &in,
                               const ResultBlock &out, std::uint64_t seed, bool apply_responses,
                               unsigned threads, bool direct = false, bool detected_only = false,
                               bool classify_losses = true) {
  detail::BlockTask task(model, in, out, seed, apply_responses, direct, detected_only,
                         classify_losses);
  parallel_for(in.size, threads, [&](std::size_t begin, std::size_t end) { task.run(begin, end); });
  return task.finish();
}

// The same on a persistent pool; the results are equal.
inline std::size_t trace_block(const OpticalModel &model, const PhotonBlock &in,
                               const ResultBlock &out, std::uint64_t seed, bool apply_responses,
                               ThreadPool &pool, bool direct = false, bool detected_only = false,
                               bool classify_losses = true) {
  detail::BlockTask task(model, in, out, seed, apply_responses, direct, detected_only,
                         classify_losses);
  parallel_for(in.size, pool, [&](std::size_t begin, std::size_t end) { task.run(begin, end); });
  return task.finish();
}

// One block of trace_blocks.
struct BlockRequest {
  const OpticalModel *model{};
  PhotonBlock in;
  ResultBlock out;
  bool apply_responses{true}, direct{false}, detected_only{false}, classify_losses{true};
};

// Trace several blocks, each with its own model, in one pool job: all chunks of
// all blocks form one flat list. Returns the rows written per block; the
// results equal one trace_block call per block.
inline std::vector<std::size_t> trace_blocks(const std::vector<BlockRequest> &requests,
                                             std::uint64_t seed, ThreadPool &pool) {
  std::vector<detail::BlockTask> tasks;
  tasks.reserve(requests.size());
  std::vector<std::size_t> first_chunk(requests.size() + 1, 0); // Chunks before each block.
  for (std::size_t b = 0; b < requests.size(); ++b) {
    const BlockRequest &r = requests[b];
    tasks.emplace_back(*r.model, r.in, r.out, seed, r.apply_responses, r.direct, r.detected_only,
                       r.classify_losses);
    first_chunk[b + 1] = first_chunk[b] + (r.in.size + kChunkSize - 1) / kChunkSize;
  }
  pool.run(first_chunk.back(), [&](std::size_t c) {
    const std::size_t b = static_cast<std::size_t>(
        std::upper_bound(first_chunk.begin(), first_chunk.end(), c) - first_chunk.begin() - 1);
    const std::size_t begin = (c - first_chunk[b]) * kChunkSize;
    tasks[b].run(begin, std::min(tasks[b].size(), begin + kChunkSize));
  });
  std::vector<std::size_t> rows;
  rows.reserve(tasks.size());
  for (const auto &task : tasks)
    rows.push_back(task.finish());
  return rows;
}

} // namespace obdeect
