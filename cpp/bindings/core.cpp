#include "obdeect/bunches.hpp"
#include "obdeect/model_json.hpp"
#include "obdeect/source.hpp"
#include "obdeect/trace.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

template <class T> using In1 = nb::ndarray<const T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using In3 = nb::ndarray<const double, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>;

// A NumPy array that owns its buffer, left uninitialised: the tracer writes
// every element.
template <class T> nb::object owned_array(std::size_t rows, std::size_t columns = 0) {
  auto *data = new T[rows * (columns ? columns : 1)];
  nb::capsule owner(data, [](void *p) noexcept { delete[] static_cast<T *>(p); });
  if (columns)
    return nb::ndarray<nb::numpy, T>(data, {rows, columns}, owner).cast();
  return nb::ndarray<nb::numpy, T>(data, {rows}, owner).cast();
}

template <class T> T *data_of(const nb::object &array) {
  return static_cast<T *>(nb::cast<nb::ndarray<nb::numpy, T>>(array).data());
}

unsigned thread_count(unsigned requested) {
  return requested ? requested : std::max(1U, std::thread::hardware_concurrency());
}

class Model {
public:
  explicit Model(const std::string &json) {
    try {
      model_ = obdeect::load_model(json);
    } catch (const obdeect::ModelError &error) {
      throw nb::value_error(error.what());
    }
  }

  // Output columns of trace(): name, element type and values per row.
  static constexpr const char *kColumns[] = {"status",
                                             "terminal_surface",
                                             "position_m",
                                             "direction",
                                             "path_length_m",
                                             "time_ns",
                                             "weight",
                                             "primary_point_m",
                                             "primary_surface",
                                             "primary_incidence_deg",
                                             "secondary_point_m",
                                             "secondary_surface",
                                             "secondary_incidence_deg",
                                             "focal_incidence_deg"};

  // The output columns of one trace: the Python arrays and the block that
  // points into them. Only the requested columns are allocated.
  struct Columns {
    nb::dict result;
    std::vector<std::string> requested;
    obdeect::ResultBlock block;
  };

  static Columns allocate(std::size_t n, const std::optional<std::vector<std::string>> &outputs,
                          bool detected_only) {
    std::vector<std::string> valid(std::begin(kColumns), std::end(kColumns));
    valid.emplace_back("index");
    Columns columns;
    if (outputs) {
      columns.requested = *outputs;
    } else {
      columns.requested.assign(std::begin(kColumns), std::end(kColumns));
    }
    const auto &requested = columns.requested;
    for (const auto &name : requested) {
      if (std::find(valid.begin(), valid.end(), name) == valid.end()) {
        std::string message = "unknown output column '" + name + "'; valid columns are:";
        for (const auto &v : valid)
          message += " " + v;
        throw nb::value_error(message.c_str());
      }
      if (name == "index" && !detected_only)
        throw nb::value_error("the output column 'index' needs detected_only=True");
    }
    const auto wants = [&](const char *name) {
      return std::find(requested.begin(), requested.end(), name) != requested.end();
    };
    nb::dict &result = columns.result;
    // The requested columns are allocated; the others stay null (not computed).
    const auto make = [&](const char *name, auto *null, std::size_t width = 0) {
      using T = std::remove_pointer_t<decltype(null)>;
      if (!wants(name))
        return static_cast<T *>(nullptr);
      nb::object array = owned_array<T>(n, width);
      result[name] = array;
      return data_of<T>(array);
    };
    constexpr std::uint8_t *u8 = nullptr;
    constexpr std::uint32_t *u32 = nullptr;
    constexpr std::uint64_t *u64 = nullptr;
    constexpr double *f64 = nullptr;
    const auto status = make("status", u8);
    const auto surface = make("terminal_surface", u32);
    const auto pos = make("position_m", f64, 3);
    const auto dir = make("direction", f64, 3);
    const auto path = make("path_length_m", f64);
    const auto time = make("time_ns", f64);
    const auto w = make("weight", f64);
    const auto m1 = make("primary_point_m", f64, 3);
    const auto m1_id = make("primary_surface", u32);
    const auto m1_angle = make("primary_incidence_deg", f64);
    const auto m2 = make("secondary_point_m", f64, 3);
    const auto m2_id = make("secondary_surface", u32);
    const auto m2_angle = make("secondary_incidence_deg", f64);
    const auto focal_angle = make("focal_incidence_deg", f64);
    const auto index = make("index", u64);
    columns.block = {status, surface, m1_id, m2_id,    pos,      dir,         path, time,
                     w,      m1,      m2,    m1_angle, m2_angle, focal_angle, index};
    return columns;
  }

  // With detected_only the arrays become views of the first `count` rows.
  static nb::dict finish(Columns &columns, std::size_t count, bool detected_only) {
    if (detected_only)
      for (const auto &name : columns.requested) {
        const nb::object array = columns.result[name.c_str()];
        columns.result[name.c_str()] = array[nb::slice(
            std::ptrdiff_t{0}, static_cast<std::ptrdiff_t>(count), std::ptrdiff_t{1})];
      }
    return columns.result;
  }

  static obdeect::PhotonBlock photons(const In3 &position, const In3 &direction,
                                      const In1<double> &wavelength_nm, const In1<double> &time_ns,
                                      const In1<double> &weight,
                                      const In1<std::uint64_t> &photon_id) {
    const std::size_t n = position.shape(0);
    if (direction.shape(0) != n || wavelength_nm.shape(0) != n || time_ns.shape(0) != n ||
        weight.shape(0) != n || photon_id.shape(0) != n)
      throw nb::value_error("photon arrays must have the same length");
    return {position.data(),
            direction.data(),
            wavelength_nm.data(),
            time_ns.data(),
            weight.data(),
            photon_id.data(),
            n};
  }

  // `pool`, if given, runs the chunks and `threads` is ignored.
  nb::dict trace(In3 position, In3 direction, In1<double> wavelength_nm, In1<double> time_ns,
                 In1<double> weight, In1<std::uint64_t> photon_id, std::uint64_t seed,
                 bool apply_responses, unsigned threads, bool direct,
                 std::optional<std::vector<std::string>> outputs, bool detected_only,
                 bool classify_losses, obdeect::ThreadPool *pool) const {
    const obdeect::PhotonBlock in =
        photons(position, direction, wavelength_nm, time_ns, weight, photon_id);
    Columns columns = allocate(in.size, outputs, detected_only);
    std::size_t count = 0;
    {
      nb::gil_scoped_release release;
      count = pool ? obdeect::trace_block(model_, in, columns.block, seed, apply_responses, *pool,
                                          direct, detected_only, classify_losses)
                   : obdeect::trace_block(model_, in, columns.block, seed, apply_responses,
                                          thread_count(threads), direct, detected_only,
                                          classify_losses);
    }
    return finish(columns, count, detected_only);
  }

  // Photon bunches (see obdeect::trace_bunches). Acceptance and density are
  // given as arrays; the output arrays are views of the rows written.
  nb::dict
  trace_bunches(In1<float> x, In1<float> y, In1<float> cx, In1<float> cy, In1<float> time_ns,
                In1<float> zem, In1<float> size, In1<float> wavelength_nm,
                In1<std::uint64_t> bunch_id, In1<double> rotation, In1<double> offset_m,
                double length_unit_m, double plane_z, std::uint64_t seed,
                const std::optional<std::tuple<In1<double>, In1<double>, In1<double>>> &acceptance,
                const std::optional<std::tuple<In1<double>, In1<double>>> &density, bool roulette,
                bool classify_losses, obdeect::ThreadPool &pool) const {
    const std::size_t n = x.shape(0);
    for (const std::size_t length :
         {y.shape(0), cx.shape(0), cy.shape(0), time_ns.shape(0), zem.shape(0), size.shape(0),
          wavelength_nm.shape(0), bunch_id.shape(0)})
      if (length != n)
        throw nb::value_error("bunch arrays must have the same length");
    if (rotation.shape(0) != 9 || offset_m.shape(0) != 3)
      throw nb::value_error("rotation needs nine elements and offset_m three");
    obdeect::BunchFrame frame;
    std::copy_n(rotation.data(), 9, frame.rotation);
    std::copy_n(offset_m.data(), 3, frame.offset_m);
    frame.length_unit_m = length_unit_m;
    frame.plane_z = plane_z;

    obdeect::Map2D acceptance_map;
    obdeect::LinearDensity wavelength_density;
    obdeect::BunchOptions options;
    options.roulette = roulette;
    options.classify_losses = classify_losses;
    try {
      if (acceptance) {
        const auto &[wavelengths, heights, values] = *acceptance;
        acceptance_map.x_m.assign(wavelengths.data(), wavelengths.data() + wavelengths.shape(0));
        acceptance_map.y_m.assign(heights.data(), heights.data() + heights.shape(0));
        acceptance_map.values.assign(values.data(), values.data() + values.shape(0));
        if (!std::is_sorted(acceptance_map.x_m.begin(), acceptance_map.x_m.end()) ||
            !std::is_sorted(acceptance_map.y_m.begin(), acceptance_map.y_m.end()))
          throw std::invalid_argument("bunch acceptance: axes must increase");
        options.acceptance = &acceptance_map;
      }
      if (density) {
        const auto &[nodes, weights] = *density;
        wavelength_density =
            obdeect::LinearDensity({nodes.data(), nodes.data() + nodes.shape(0)},
                                   {weights.data(), weights.data() + weights.shape(0)});
        options.wavelength_density = &wavelength_density;
      }
    } catch (const std::invalid_argument &error) {
      throw nb::value_error(error.what());
    }

    const obdeect::BunchBlock block{x.data(),        y.data(),   cx.data(),   cy.data(),
                                    time_ns.data(),  zem.data(), size.data(), wavelength_nm.data(),
                                    bunch_id.data(), n};
    static constexpr const char *names[] = {"position_m",    "direction",   "time_ns",   "weight",
                                            "wavelength_nm", "bunch_index", "sub_photon"};
    nb::dict result;
    std::size_t count = 0;
    try {
      nb::gil_scoped_release release;
      count =
          obdeect::trace_bunches(model_, block, frame, options, seed, pool, [&](std::size_t rows) {
            nb::gil_scoped_acquire acquire; // The arrays are Python objects.
            const auto position = owned_array<double>(rows, 3);
            const auto direction = owned_array<double>(rows, 3);
            const auto time = owned_array<double>(rows), weight = owned_array<double>(rows);
            const auto wavelength = owned_array<double>(rows);
            const auto index = owned_array<std::uint64_t>(rows);
            const auto sub = owned_array<std::uint32_t>(rows);
            result["position_m"] = position;
            result["direction"] = direction;
            result["time_ns"] = time;
            result["weight"] = weight;
            result["wavelength_nm"] = wavelength;
            result["bunch_index"] = index;
            result["sub_photon"] = sub;
            return obdeect::BunchResult{data_of<double>(position),   data_of<double>(direction),
                                        data_of<double>(time),       data_of<double>(weight),
                                        data_of<double>(wavelength), data_of<std::uint64_t>(index),
                                        data_of<std::uint32_t>(sub)};
          });
    } catch (const obdeect::ModelError &error) {
      throw nb::value_error(error.what());
    }
    for (const char *name : names) {
      const nb::object array = result[name];
      result[name] = array[nb::slice(std::ptrdiff_t{0}, static_cast<std::ptrdiff_t>(count),
                                     std::ptrdiff_t{1})];
    }
    return result;
  }

  [[nodiscard]] const obdeect::OpticalModel &optical_model() const { return model_; }

private:
  obdeect::OpticalModel model_;
};

// Photon arrays of `count` rows filled by photon(i) -> (ray, time, weight).
template <class Photon>
nb::dict photon_arrays(std::size_t count, std::uint64_t first_id, double wavelength_nm,
                       const Photon &photon) {
  const auto position = owned_array<double>(count, 3), direction = owned_array<double>(count, 3);
  const auto time = owned_array<double>(count), wavelength = owned_array<double>(count);
  const auto weight = owned_array<double>(count);
  const auto id = owned_array<std::uint64_t>(count);
  double *p = data_of<double>(position), *d = data_of<double>(direction),
         *t = data_of<double>(time);
  double *l = data_of<double>(wavelength), *w = data_of<double>(weight);
  std::uint64_t *k = data_of<std::uint64_t>(id);
  {
    nb::gil_scoped_release release;
    obdeect::parallel_for(count, thread_count(0), [&](std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; ++i) {
        const auto [ray, time_ns, photon_weight] = photon(first_id + i);
        p[3 * i] = ray.origin.x;
        p[3 * i + 1] = ray.origin.y;
        p[3 * i + 2] = ray.origin.z;
        d[3 * i] = ray.direction.x;
        d[3 * i + 1] = ray.direction.y;
        d[3 * i + 2] = ray.direction.z;
        t[i] = time_ns;
        l[i] = wavelength_nm;
        w[i] = photon_weight;
        k[i] = first_id + i;
      }
    });
  }
  nb::dict result;
  result["position_m"] = position;
  result["direction"] = direction;
  result["time_ns"] = time;
  result["wavelength_nm"] = wavelength;
  result["weight"] = weight;
  result["photon_id"] = id;
  return result;
}

obdeect::Vec3 vec3(const In1<double> &v, const char *name) {
  if (v.shape(0) != 3)
    throw nb::value_error((std::string(name) + " must have three components").c_str());
  const obdeect::Vec3 result{v(0), v(1), v(2)};
  if (!obdeect::finite(result))
    throw nb::value_error((std::string(name) + " must be finite").c_str());
  return result;
}

nb::dict sample_source(In1<double> towards_source, double distance_m, double disk_radius_m,
                       double disk_z_m, double start_z_m, double wavelength_nm,
                       std::uint64_t first_id, std::size_t count, std::uint64_t seed) {
  if (towards_source.shape(0) != 3)
    throw nb::value_error("towards_source must have three components");
  const obdeect::Vec3 towards{towards_source(0), towards_source(1), towards_source(2)};
  if (std::abs(obdeect::norm(towards) - 1.0) > 1e-9 || towards.z <= 0.0)
    throw nb::value_error("towards_source must be a unit vector pointing upwards");
  if (!(distance_m > 0.0) || !(disk_radius_m > 0.0) || !(wavelength_nm > 0.0) ||
      !std::isfinite(disk_z_m) || !std::isfinite(start_z_m))
    throw nb::value_error("distance, disk radius and wavelength must be positive; heights finite");
  const obdeect::PointSource source{towards,  distance_m, disk_radius_m,
                                    disk_z_m, start_z_m,  wavelength_nm};
  return photon_arrays(count, first_id, wavelength_nm, [&](std::uint64_t id) {
    const auto photon = obdeect::sample_source(source, id, seed);
    return std::tuple{photon.ray, photon.time_ns, 1.0};
  });
}

// A point emitter (see obdeect::Emitter).
class Emitter {
public:
  Emitter(In1<double> position_m, In1<double> axis, In1<double> intensity_u, In1<double> intensity,
          In1<double> pulse_time_ns, In1<double> pulse, In1<double> disk_centre_m,
          In1<double> disk_normal, double disk_radius_m, double wavelength_nm)
      : emitter_{vec3(position_m, "position_m"),
                 vec3(axis, "axis"),
                 density(intensity_u, intensity),
                 density(pulse_time_ns, pulse),
                 vec3(disk_centre_m, "disk_centre_m"),
                 vec3(disk_normal, "disk_normal"),
                 disk_radius_m,
                 wavelength_nm} {
    if (intensity_u.shape(0) == 0)
      throw nb::value_error("the intensity table must not be empty");
    if (std::abs(obdeect::norm(emitter_.axis) - 1.0) > 1e-9 ||
        std::abs(obdeect::norm(emitter_.disk_normal) - 1.0) > 1e-9)
      throw nb::value_error("axis and disk_normal must be unit vectors");
    if (!(disk_radius_m > 0.0) || !(wavelength_nm > 0.0) || !(emitter_.disk_distance_m() > 0.0))
      throw nb::value_error(
          "disk radius and wavelength must be positive and the emitter off the disk plane");
  }

  [[nodiscard]] double candidates_per_photon() const { return emitter_.candidates_per_photon(); }

  [[nodiscard]] nb::dict sample(std::uint64_t first_id, std::size_t count,
                                std::uint64_t seed) const {
    return photon_arrays(count, first_id, emitter_.wavelength_nm, [&](std::uint64_t id) {
      const auto photon = obdeect::sample_emitter(emitter_, id, seed);
      return std::tuple{photon.ray, photon.time_ns, photon.weight};
    });
  }

private:
  static obdeect::LinearDensity density(const In1<double> &x, const In1<double> &y) {
    if (x.shape(0) == 0 && y.shape(0) == 0)
      return {};
    try {
      return {{x.data(), x.data() + x.shape(0)}, {y.data(), y.data() + y.shape(0)}};
    } catch (const std::invalid_argument &error) {
      throw nb::value_error(error.what());
    }
  }

  obdeect::Emitter emitter_;
};

// One request of trace_many: (Model, photons, apply_responses, direct, outputs,
// detected_only, classify_losses), where photons holds the contiguous arrays that
// Model.trace takes.
using ManyItem =
    std::tuple<Model *, nb::dict, bool, bool, std::optional<std::vector<std::string>>, bool, bool>;

nb::list trace_many(const std::vector<ManyItem> &items, std::uint64_t seed,
                    obdeect::ThreadPool &pool) {
  struct Prepared {
    In3 position, direction;
    In1<double> wavelength_nm, time_ns, weight;
    In1<std::uint64_t> photon_id;
    Model::Columns columns;
    bool detected_only;
  };
  std::vector<Prepared> prepared;
  std::vector<obdeect::BlockRequest> requests;
  prepared.reserve(items.size());
  for (const auto &[model, photon_dict, apply_responses, direct, outputs, detected_only,
                    classify_losses] : items) {
    const auto get = [&photon_dict](const char *name) { return photon_dict[name]; };
    Prepared p{nb::cast<In3>(get("position_m")),
               nb::cast<In3>(get("direction")),
               nb::cast<In1<double>>(get("wavelength_nm")),
               nb::cast<In1<double>>(get("time_ns")),
               nb::cast<In1<double>>(get("weight")),
               nb::cast<In1<std::uint64_t>>(get("photon_id")),
               {},
               detected_only};
    const obdeect::PhotonBlock in =
        Model::photons(p.position, p.direction, p.wavelength_nm, p.time_ns, p.weight, p.photon_id);
    p.columns = Model::allocate(in.size, outputs, detected_only);
    requests.push_back({&model->optical_model(), in, p.columns.block, apply_responses, direct,
                        detected_only, classify_losses});
    prepared.push_back(std::move(p));
  }
  std::vector<std::size_t> rows;
  {
    nb::gil_scoped_release release;
    rows = obdeect::trace_blocks(requests, seed, pool);
  }
  nb::list results;
  for (std::size_t i = 0; i < prepared.size(); ++i)
    results.append(Model::finish(prepared[i].columns, rows[i], prepared[i].detected_only));
  return results;
}

} // namespace

NB_MODULE(_core, module) {
  module.doc() = "obdeect native optical transport";
  nb::class_<obdeect::ThreadPool>(module, "ThreadPool")
      .def(
          "__init__",
          [](obdeect::ThreadPool *self, unsigned threads) {
            new (self) obdeect::ThreadPool(thread_count(threads));
          },
          "threads"_a = 0U)
      .def_prop_ro("size", &obdeect::ThreadPool::size);
  nb::class_<Model>(module, "Model")
      .def(nb::init<const std::string &>(), "json"_a)
      .def("trace", &Model::trace, "position_m"_a, "direction"_a, "wavelength_nm"_a, "time_ns"_a,
           "weight"_a, "photon_id"_a, nb::kw_only(), "seed"_a, "apply_responses"_a = true,
           "threads"_a = 0U, "direct"_a = false, "outputs"_a = nb::none(),
           "detected_only"_a = false, "classify_losses"_a = true, "pool"_a = nb::none())
      .def("trace_bunches", &Model::trace_bunches, "x"_a, "y"_a, "cx"_a, "cy"_a, "time_ns"_a,
           "zem"_a, "size"_a, "wavelength_nm"_a, "bunch_id"_a, "rotation"_a, "offset_m"_a,
           "length_unit_m"_a, "plane_z"_a, nb::kw_only(), "seed"_a, "acceptance"_a = nb::none(),
           "wavelength_density"_a = nb::none(), "roulette"_a = true, "classify_losses"_a = true,
           "pool"_a);
  module.def("trace_many", &trace_many, "items"_a, nb::kw_only(), "seed"_a, "pool"_a);
  nb::class_<Emitter>(module, "Emitter")
      .def(nb::init<In1<double>, In1<double>, In1<double>, In1<double>, In1<double>, In1<double>,
                    In1<double>, In1<double>, double, double>(),
           "position_m"_a, "axis"_a, "intensity_u"_a, "intensity"_a, "pulse_time_ns"_a, "pulse"_a,
           "disk_centre_m"_a, "disk_normal"_a, "disk_radius_m"_a, "wavelength_nm"_a)
      .def("candidates_per_photon", &Emitter::candidates_per_photon)
      .def("sample", &Emitter::sample, "first_id"_a, "count"_a, "seed"_a);
  module.def("sample_source", &sample_source, "towards_source"_a, "distance_m"_a, "disk_radius_m"_a,
             "disk_z_m"_a, "start_z_m"_a, "wavelength_nm"_a, "first_id"_a, "count"_a, "seed"_a);
  nb::list names;
  for (std::size_t i = 0; i < obdeect::kStatusCount; ++i)
    names.append(nb::str(std::string(obdeect::to_string(static_cast<obdeect::Status>(i))).c_str()));
  module.attr("STATUS_NAMES") = nb::tuple(names);
  module.attr("MODEL_FORMAT") = std::string(obdeect::kModelFormat);
}
