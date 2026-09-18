#include "obdeect/abi.hpp"
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

  std::cout << "contract tests passed\n";
}
