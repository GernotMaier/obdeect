#include "obdeect/segmented_path.hpp"

#include <cassert>

int main() {
  using namespace obdeect;
  {
    ImportedSegmentedOpticalModel input{
        {"local-response-fixture", "1", std::string(64, 'a')},
        {{7, {0.5, 0.5, 0}, {0, 0, 1}, 2, 10, FacetShape::circle, {1, 0, 0}}},
        {{8, {0.5, 0.5, 3}, {0, 0, 1}, 2, FacetShape::circle, {1, 0, 0}, 0, -1}},
        {},
        SpectralResponse{{300, 500}, {0.8, 0.8}}};
    auto model = compile_segmented_optical_model(input);
    assert(model);
    model->camera_degradation =
        SpatialResponse{{-0.1, 0.1}, {-0.1, 0.1}, {0.25, 0.75, 0.25, 0.75}, true};
    const Ray ray{{0.5, 0.55, 5}, {0, 0, -1}};
    assert(trace_segmented_path(ray, 1, 400, *model).surviving_throughput == 0);
    model->camera_degradation_in_detector_frame = true;
    assert(is_valid(*model));
    model->pixel_responses = std::make_shared<const PixelResponses>();
    assert(!is_valid(*model));
    model->pixel_responses.reset();
    assert(std::abs(trace_segmented_path(ray, 1, 400, *model).surviving_throughput - 0.3) < 1e-15);
    model->detector_surfaces[0].response_y_sign = 1;
    assert(std::abs(trace_segmented_path(ray, 1, 400, *model).surviving_throughput - 0.5) < 1e-15);
  }
  {
    SpectralResponse sampled{{300, 500}, {0.2, 0.8}};
    sampled.wavelength_sampling = WavelengthSampling{1, 0, 200, 999};
    assert(sampled.is_valid());
    assert(std::abs(*sampled.at(400.9) - 0.5) < 1e-15);
    sampled.wavelength_sampling->offset_nm = 0.5;
    assert(std::abs(*sampled.at(400.9) - 0.5015) < 1e-15);
    assert(*sampled.at(199.9) == 0 && *sampled.at(1000) == 0);
    sampled.wavelength_sampling->offset_nm = 1;
    assert(!sampled.is_valid());
    SpectralResponse angular{{300, 500}, {0.2, 0.4, 0.4, 0.8}, {0, 60}};
    angular.wavelength_sampling = WavelengthSampling{1, 0, 200, 999};
    angular.spectral_envelope = {0.4, 0.8};
    assert(angular.is_valid());
    assert(std::abs(*angular.at(400.9, 30) - 0.45) < 1e-15);
    angular.wavelength_sampling->offset_nm = 0.5;
    assert(std::abs(*angular.at(400.9, 30) - 0.45075) < 1e-15);
    angular.wavelength_sampling->offset_nm = 0;
    angular.wavelength_sampling->projection_angle_deg = 0;
    assert(std::abs(*angular.at(400.9, 30) - 0.225) < 1e-15);
    angular.relative_to_envelope = true;
    assert(std::abs(*angular.at(400.9, 30) - 0.75) < 1e-15);
    angular.wavelength_sampling.reset();
    assert(angular.is_valid());
    assert(std::abs(*angular.at(400.9, 30) - 0.75) < 1e-15);
    angular.envelope_interpolation.scheme = 2;
    assert(!angular.is_valid());
    angular.envelope_interpolation.scheme = 1;
    angular.spectral_envelope[0] = 1.1;
    assert(!angular.is_valid());
  }
  const Vec3 incoming{0, 0, -1}, normal{0, 0, 1}, tangent{1, 0, 0};
  MirrorScatter scatter{0.001, 0, 0.004, MirrorScatterMethod::outgoing_angles};
  assert(scatter.is_valid());
  double first_sum{}, second_sum{}, squared_sum{};
  for (std::uint64_t id = 0; id < 20000; ++id) {
    const auto reflected = reflect_with_scatter(incoming, normal, tangent, scatter, id, 7, 123);
    const auto repeat = reflect_with_scatter(incoming, normal, tangent, scatter, id, 7, 123);
    assert(reflected && repeat && norm(*reflected) > 1 - 1e-14 && norm(*reflected) < 1 + 1e-14);
    assert(reflected->x == repeat->x && reflected->y == repeat->y && reflected->z == repeat->z);
    first_sum += reflected->x;
    second_sum += reflected->y;
    squared_sum += reflected->x * reflected->x + reflected->y * reflected->y;
  }
  assert(reflect_with_scatter(incoming, normal, tangent, scatter, 1, 7, 123)->x !=
         reflect_with_scatter(incoming, normal, tangent, scatter, 1, 7, 124)->x);
  assert(std::abs(first_sum / 20000) < 3e-5 && std::abs(second_sum / 20000) < 3e-5);
  assert(std::abs(squared_sum / 20000 - 2e-6) < 8e-8);
  scatter.sigma1_rad = 0;
  assert(reflect_with_scatter(incoming, normal, tangent, scatter, 1, 7, 123)->z == 1);
  scatter.fraction2 = 1;
  assert(reflect_with_scatter(incoming, normal, tangent, scatter, 1, 7, 123)->z != 1);
  scatter.method = MirrorScatterMethod::surface_slopes;
  assert(reflect_with_scatter(incoming, normal, tangent, scatter, 1, 7, 123));
  scatter.fraction2 = 2;
  assert(!scatter.is_valid());

  ImportedSegmentedOpticalModel input{{"synthetic-response", "1", std::string(64, 'a')},
                                      {{7, {}, normal, 2, 10, FacetShape::circle, tangent, 20}},
                                      {{8, {0, 0, 3}, normal, 2, FacetShape::circle, tangent}},
                                      {},
                                      SpectralResponse{{300, 500}, {0.8, 0.8}}};
  input.camera_response = CameraResponse{0.5, SpectralResponse{{300, 500}, {0.8, 0.8}},
                                         CameraIncidenceResponse{{0, 90}, {0.9, 0.9}}};
  auto model = compile_segmented_optical_model(input);
  assert(model);
  const auto path = trace_segmented_path({{0, 0, 5}, incoming}, 1, 400, *model);
  assert(path.status == PhotonStatus::detected);
  assert(std::abs(path.surviving_throughput - 0.8 * 0.5 * 0.8 * 0.9) < 1e-15);
  assert(path.interaction_throughput[1] == path.surviving_throughput);
  input.camera_response->camera_transmission = 0;
  model = compile_segmented_optical_model(input);
  assert(model);
  const auto zero = trace_segmented_path({{0, 0, 5}, incoming}, 1, 400, *model);
  assert(zero.status == PhotonStatus::detected && zero.surviving_throughput == 0);
  input.primary_scatter = MirrorScatter{0.001, 0, 0, MirrorScatterMethod::outgoing_angles};
  model = compile_segmented_optical_model(input);
  const auto scattered = trace_segmented_path({{0, 0, 5}, incoming}, 1, 400, *model, 123);
  assert(scattered.status == PhotonStatus::detected && scattered.final_direction.x != 0);
  input.primary_scatter->sigma1_rad = -1;
  assert(!compile_segmented_optical_model(input));
}
