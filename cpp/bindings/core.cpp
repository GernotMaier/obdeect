#include "obdeect/material_path.hpp"
#include "obdeect/optical_model_file.hpp"
#include "obdeect/segmented_path.hpp"
#include "obdeect/sources.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <memory>

namespace nb = nanobind;
using namespace nb::literals;
using Coordinates =
    nb::ndarray<nb::numpy, const double, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>;
using Scalars = nb::ndarray<nb::numpy, const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using Ids = nb::ndarray<nb::numpy, const std::uint64_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

namespace {
struct Output {
  std::vector<double> position, direction, path, optical_path, time, weight, throughput,
      response_loss, terminal_loss;
  std::vector<std::uint64_t> ids;
  std::vector<std::uint32_t> surfaces;
  std::vector<std::uint8_t> status;
  explicit Output(std::size_t n)
      : position(3 * n), direction(3 * n), path(n), optical_path(n), time(n), weight(n),
        throughput(n), response_loss(n), terminal_loss(n), ids(n), surfaces(n), status(n) {}
};

template <class T>
nb::object array(std::vector<T> &values, const nb::capsule &owner, std::size_t n,
                 bool coordinates = false) {
  if (coordinates)
    return nb::ndarray<nb::numpy, T>(values.data(), {n, 3}, owner).cast();
  return nb::ndarray<nb::numpy, T>(values.data(), {n}, owner).cast();
}

class OpticalModel {
public:
  explicit OpticalModel(const std::string &path) {
    if (auto loaded = obdeect::read_optical_model(path)) {
      segmented_ = std::move(loaded->segmented);
      axisymmetric_ = std::move(loaded->axisymmetric);
      if (loaded->nonsequential)
        nonsequential_.emplace(std::move(*loaded->nonsequential));
    }
    if (!segmented_ && !axisymmetric_ && !nonsequential_)
      throw nb::value_error("Cannot load optical model: invalid schema, content hash or geometry");
  }

  std::string hash() const {
    return segmented_      ? segmented_->provenance.content_hash
           : axisymmetric_ ? axisymmetric_->provenance.content_hash
                           : nonsequential_->provenance().content_hash;
  }

  nb::dict trace(Coordinates positions, Coordinates directions, Scalars wavelengths, Scalars times,
                 Scalars weights, Ids ids) const {
    const auto count = positions.shape(0);
    if (directions.shape(0) != count || wavelengths.shape(0) != count || times.shape(0) != count ||
        weights.shape(0) != count || ids.shape(0) != count)
      throw nb::value_error("Photon arrays must have the same length");
    std::vector<std::uint64_t> identities;
    if (count)
      identities.assign(ids.data(), ids.data() + count);
    std::sort(identities.begin(), identities.end());
    if (std::adjacent_find(identities.begin(), identities.end()) != identities.end())
      throw nb::value_error("Photon IDs must be unique within the batch");
    for (std::size_t i = 0; i < count; ++i) {
      const obdeect::Vec3 direction{directions(i, 0), directions(i, 1), directions(i, 2)};
      if (!std::isfinite(positions(i, 0)) || !std::isfinite(positions(i, 1)) ||
          !std::isfinite(positions(i, 2)) || !std::isfinite(obdeect::norm(direction)) ||
          std::abs(obdeect::norm(direction) - 1) > 1e-12 ||
          !obdeect::is_valid_wavelength(wavelengths(i)) || !std::isfinite(times(i)) ||
          !std::isfinite(weights(i)) || weights(i) < 0)
        throw nb::value_error(
            "Photon positions, unit directions, wavelengths, times and weights must be valid");
    }
    auto output = std::make_unique<Output>(count);
    {
      nb::gil_scoped_release release;
      for (std::size_t i = 0; i < count; ++i) {
        const obdeect::Ray ray{{positions(i, 0), positions(i, 1), positions(i, 2)},
                               {directions(i, 0), directions(i, 1), directions(i, 2)}};
        const auto store_result = [&](const auto &result) {
          const auto point = result.points_m[result.point_count - 1];
          output->position[3 * i] = point.x;
          output->position[3 * i + 1] = point.y;
          output->position[3 * i + 2] = point.z;
          output->direction[3 * i] = result.final_direction.x;
          output->direction[3 * i + 1] = result.final_direction.y;
          output->direction[3 * i + 2] = result.final_direction.z;
          output->path[i] = result.path_length_m;
          output->optical_path[i] =
              result.material_transport ? result.optical_path_m : result.path_length_m;
          output->time[i] = times(i) + (result.material_transport
                                            ? result.group_delay_ns
                                            : result.path_length_m / obdeect::kSpeedOfLightMPerNs);
          output->throughput[i] =
              result.status == obdeect::PhotonStatus::detected ? result.surviving_throughput : 0;
          output->weight[i] = weights(i) * output->throughput[i];
          output->response_loss[i] = weights(i) * (1 - result.surviving_throughput);
          output->terminal_loss[i] = result.status == obdeect::PhotonStatus::detected
                                         ? 0
                                         : weights(i) * result.surviving_throughput;
          output->ids[i] = ids(i);
          output->surfaces[i] = result.terminal_surface_id;
          output->status[i] = static_cast<std::uint8_t>(result.status);
        };
        if (nonsequential_)
          store_result(
              obdeect::trace_material_path<false>(ray, ids(i), wavelengths(i), *nonsequential_));
        else {
          const auto result =
              segmented_ ? obdeect::trace_segmented_path(ray, ids(i), wavelengths(i), *segmented_)
                         : obdeect::trace_axisymmetric_optical_model(ray, ids(i), *axisymmetric_,
                                                                     wavelengths(i));
          store_result(result);
        }
      }
    }
    Output *storage = output.get();
    const nb::capsule owner(storage,
                            [](void *pointer) noexcept { delete static_cast<Output *>(pointer); });
    output.release();
    nb::dict result;
    result["position_m"] = array(storage->position, owner, count, true);
    result["direction"] = array(storage->direction, owner, count, true);
    result["path_length_m"] = array(storage->path, owner, count);
    result["optical_path_m"] = array(storage->optical_path, owner, count);
    result["arrival_time_ns"] = array(storage->time, owner, count);
    result["optical_weight"] = array(storage->weight, owner, count);
    result["throughput"] = array(storage->throughput, owner, count);
    result["response_loss_weight"] = array(storage->response_loss, owner, count);
    result["terminal_loss_weight"] = array(storage->terminal_loss, owner, count);
    result["photon_id"] = array(storage->ids, owner, count);
    result["terminal_surface_id"] = array(storage->surfaces, owner, count);
    result["status"] = array(storage->status, owner, count);
    return result;
  }

private:
  std::optional<obdeect::CompiledSegmentedOpticalModel> segmented_;
  std::optional<obdeect::AxisymmetricOpticalModel> axisymmetric_;
  std::optional<obdeect::CompiledOpticalModel> nonsequential_;
};
} // namespace

NB_MODULE(_core, module) {
  module.doc() = "Bulk nominal optical transport through an immutable compiled optical model";
  nb::class_<OpticalModel>(module, "OpticalModel")
      .def(nb::init<const std::string &>(), "path"_a)
      .def_prop_ro("optical_model_sha256", &OpticalModel::hash)
      .def("trace", &OpticalModel::trace, "position_m"_a.noconvert(), "direction"_a.noconvert(),
           "wavelength_nm"_a.noconvert(), "emission_time_ns"_a.noconvert(),
           "source_weight"_a.noconvert(), "photon_id"_a.noconvert());
  nb::list statuses;
  for (std::size_t i = 0; i < obdeect::kPhotonStatusCount; ++i)
    statuses.append(nb::str(obdeect::to_string(static_cast<obdeect::PhotonStatus>(i)).data()));
  module.attr("STATUS_NAMES") = nb::tuple(statuses);
}
