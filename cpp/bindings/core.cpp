#include "obdeect/model_json.hpp"
#include "obdeect/source.hpp"
#include "obdeect/trace.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <memory>
#include <string>
#include <thread>
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
                 bool apply_responses, unsigned threads) const {
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
      obdeect::trace_block(model_, in, out, seed, apply_responses, thread_count(threads));
    }
    return result;
  }

private:
  obdeect::OpticalModel model_;
};

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
  nb::dict result;
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
        const auto photon = obdeect::sample_source(source, first_id + i, seed);
        p[3 * i] = photon.ray.origin.x;
        p[3 * i + 1] = photon.ray.origin.y;
        p[3 * i + 2] = photon.ray.origin.z;
        d[3 * i] = photon.ray.direction.x;
        d[3 * i + 1] = photon.ray.direction.y;
        d[3 * i + 2] = photon.ray.direction.z;
        t[i] = photon.time_ns;
        l[i] = wavelength_nm;
        w[i] = 1.0;
        k[i] = first_id + i;
      }
    });
  }
  result["position_m"] = position;
  result["direction"] = direction;
  result["time_ns"] = time;
  result["wavelength_nm"] = wavelength;
  result["weight"] = weight;
  result["photon_id"] = id;
  return result;
}

} // namespace

NB_MODULE(_core, module) {
  module.doc() = "obdeect native optical transport";
  nb::class_<Model>(module, "Model")
      .def(nb::init<const std::string &>(), "json"_a)
      .def("trace", &Model::trace, "position_m"_a, "direction"_a, "wavelength_nm"_a, "time_ns"_a,
           "weight"_a, "photon_id"_a, nb::kw_only(), "seed"_a, "apply_responses"_a = true,
           "threads"_a = 0U);
  module.def("sample_source", &sample_source, "towards_source"_a, "distance_m"_a, "disk_radius_m"_a,
             "disk_z_m"_a, "start_z_m"_a, "wavelength_nm"_a, "first_id"_a, "count"_a, "seed"_a);
  nb::list names;
  for (std::size_t i = 0; i < obdeect::kStatusCount; ++i)
    names.append(nb::str(std::string(obdeect::to_string(static_cast<obdeect::Status>(i))).c_str()));
  module.attr("STATUS_NAMES") = nb::tuple(names);
  module.attr("MODEL_FORMAT") = std::string(obdeect::kModelFormat);
}
