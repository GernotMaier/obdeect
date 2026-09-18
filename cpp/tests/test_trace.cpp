#include "obdeect/trace.hpp"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

void require(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}

}  // namespace

int main() {
  using namespace obdeect;
  ToyMstConfig configuration{};
  configuration.include_structure = false;
  const auto scene = compile_scene(configuration);
  require(scene.has_value(), "valid toy scene compiles");

  const std::vector<Vec3> positions{{0.0, 0.0, 20.0}, {10.0, 0.0, 20.0}};
  const std::vector<Vec3> directions{{0.0, 0.0, -1.0}, {0.0, 0.0, -1.0}};
  const std::vector<double> wavelength{400.0, 400.0};
  const std::vector<double> time{3.0, 4.0};
  const std::vector<double> weight{2.0, 3.0};
  const std::vector<std::uint64_t> ids{10, 11};
  const PhotonBlockView input{positions, directions, wavelength, time, weight, ids};
  const auto result = trace(*scene, input);
  require(result.photons.status[0] == PhotonStatus::detected, "on-axis ray detected");
  require(result.photons.status[1] == PhotonStatus::missed_primary, "outside aperture ray missed");
  require(result.photons.weight[0] == 2.0 && result.photons.weight[1] == 0.0,
          "terminal optical weights retained exactly once");
  require(result.photons.direction[0].z > 0.0, "detected ray exposes reflected direction");
  require(result.photons.time_ns[0] > time[0], "trace adds geometric flight time once");
  require(result.summary.status_count[static_cast<std::size_t>(PhotonStatus::detected)] == 1,
          "summary detected count");
  require(result.summary.status_count[static_cast<std::size_t>(PhotonStatus::missed_primary)] == 1,
          "summary missed count");

  const PhotonBlockView malformed{positions, directions, wavelength, time, weight, {ids.data(), 1}};
  const auto invalid = trace(*scene, malformed);
  require(invalid.photons.status[0] == PhotonStatus::invalid_input,
          "inconsistent public block fails closed");
  require(invalid.summary.status_count[static_cast<std::size_t>(PhotonStatus::invalid_input)] == positions.size(),
          "invalid block summary closes status count");

  std::cout << "trace tests passed\n";
}
