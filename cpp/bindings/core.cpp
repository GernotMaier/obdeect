#include "obdeect/model_json.hpp"
#include "obdeect/source.hpp"
#include "obdeect/trace.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <memory>
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

  nb::dict trace(In3 position, In3 direction, In1<double> wavelength_nm, In1<double> time_ns,
                 In1<double> weight, In1<std::uint64_t> photon_id, std::uint64_t seed,
                 bool apply_responses, unsigned threads, bool direct) const {
    const std::size_t n = position.shape(0);
    if (direction.shape(0) != n || wavelength_nm.shape(0) != n || time_ns.shape(0) != n ||
        weight.shape(0) != n || photon_id.shape(0) != n)
      throw nb::value_error("photon arrays must have the same length");
    nb::dict result;
    const auto add = [&](const char *name, nb::object array) {
      result[name] = array;
      return array;
    };
    const auto status = add("status", owned_array<std::uint8_t>(n));
    const auto surface = add("terminal_surface", owned_array<std::uint32_t>(n));
    const auto pos = add("position_m", owned_array<double>(n, 3));
    const auto dir = add("direction", owned_array<double>(n, 3));
    const auto path = add("path_length_m", owned_array<double>(n));
    const auto time = add("time_ns", owned_array<double>(n));
    const auto w = add("weight", owned_array<double>(n));
    const auto m1 = add("primary_point_m", owned_array<double>(n, 3));
    const auto m1_id = add("primary_surface", owned_array<std::uint32_t>(n));
    const auto m1_angle = add("primary_incidence_deg", owned_array<double>(n));
    const auto m2 = add("secondary_point_m", owned_array<double>(n, 3));
    const auto m2_id = add("secondary_surface", owned_array<std::uint32_t>(n));
    const auto m2_angle = add("secondary_incidence_deg", owned_array<double>(n));
    const auto focal_angle = add("focal_incidence_deg", owned_array<double>(n));
    const obdeect::PhotonBlock in{position.data(),
                                  direction.data(),
                                  wavelength_nm.data(),
                                  time_ns.data(),
                                  weight.data(),
                                  photon_id.data(),
                                  n};
    const obdeect::ResultBlock out{data_of<std::uint8_t>(status),
                                   data_of<std::uint32_t>(surface),
                                   data_of<std::uint32_t>(m1_id),
                                   data_of<std::uint32_t>(m2_id),
                                   data_of<double>(pos),
                                   data_of<double>(dir),
                                   data_of<double>(path),
                                   data_of<double>(time),
                                   data_of<double>(w),
                                   data_of<double>(m1),
                                   data_of<double>(m2),
                                   data_of<double>(m1_angle),
                                   data_of<double>(m2_angle),
                                   data_of<double>(focal_angle)};
    {
      nb::gil_scoped_release release;
      obdeect::trace_block(model_, in, out, seed, apply_responses, thread_count(threads), direct);
    }
    return result;
  }

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

} // namespace

NB_MODULE(_core, module) {
  module.doc() = "obdeect native optical transport";
  nb::class_<Model>(module, "Model")
      .def(nb::init<const std::string &>(), "json"_a)
      .def("trace", &Model::trace, "position_m"_a, "direction"_a, "wavelength_nm"_a, "time_ns"_a,
           "weight"_a, "photon_id"_a, nb::kw_only(), "seed"_a, "apply_responses"_a = true,
           "threads"_a = 0U, "direct"_a = false);
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
