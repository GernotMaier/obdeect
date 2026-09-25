#include "obdeect/abi.hpp"
#include "obdeect/frames.hpp"
#include "obdeect/tables.hpp"

#include <cmath>
#include <cstdlib>
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

  // T-TAB-001: linear interpolation and out-of-range rejection.
  const std::vector<double> wavelength{300.0, 400.0, 500.0};
  const std::vector<double> reflectivity{0.8, 0.9, 0.7};
  const Table1DView table{wavelength, reflectivity};
  require(table.is_valid(), "monotonic response table valid");
  const auto interpolated = table.interpolate(450.0);
  require(interpolated.has_value() && std::abs(*interpolated - 0.8) < 1e-14,
          "linear wavelength interpolation");
  require(!table.interpolate(250.0), "response table must not extrapolate silently");

  // T-ABI-002: every SoA field is required and unit direction is enforced.
  const std::vector<Vec3> position{{0.0, 0.0, 20.0}};
  const std::vector<Vec3> direction{{0.0, 0.0, -1.0}};
  const std::vector<double> photon_wavelength{400.0};
  const std::vector<double> time{0.0};
  const std::vector<double> weight{1.0};
  const std::vector<std::uint64_t> id{1};
  const PhotonBlockView valid{position, direction, photon_wavelength, time, weight, id};
  require(validate_photon_block(valid), "valid SoA block accepted");
  const PhotonBlockView inconsistent{position, direction, photon_wavelength, time, weight, {}};
  require(!validate_photon_block(inconsistent), "inconsistent SoA block rejected");
  const std::vector<Vec3> nonunit_direction{{0.0, 0.0, -2.0}};
  const PhotonBlockView nonunit{position, nonunit_direction, photon_wavelength, time, weight, id};
  require(!validate_photon_block(nonunit), "non-unit direction rejected");

  // T-FRAME-001: rotation and translation preserve a physical ray through a
  // parent/local round trip. Reflections and scaled axes are rejected.
  const RigidFrame frame{{2.0, -3.0, 5.0}, {0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
                         {0.0, 0.0, 1.0}};
  require(frame.is_valid(), "right-handed rigid frame accepted");
  const Ray local{{1.0, 2.0, 3.0}, {0.0, 0.0, -1.0}};
  const Ray parent = frame.ray_to_parent(local);
  const Ray round_trip = frame.ray_from_parent(parent);
  require(parent.position_m.x == 0.0 && parent.position_m.y == -2.0 &&
              parent.position_m.z == 8.0 && round_trip.position_m.x == local.position_m.x &&
              round_trip.position_m.y == local.position_m.y &&
              round_trip.position_m.z == local.position_m.z &&
              round_trip.direction.z == local.direction.z,
          "ray round trip preserves position and direction");
  RigidFrame reflected = frame;
  reflected.z_axis = {0.0, 0.0, -1.0};
  require(!reflected.is_valid(), "left-handed frame rejected");

  std::cout << "contract tests passed\n";
}
