#include "obdeect/eventio_transport.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
}  // namespace

int main() {
  using namespace obdeect;
  try {
    const auto path = std::filesystem::temp_directory_path() / "obdeect-atmosphere-test.dat";
    {
      std::ofstream output(path);
      output << "# H2= 2.0, H1= 3.0 4.0\n300 1 2\n600 1 2\n";
    }
    const auto table = AtmosphereTransmissionTable::from_simtel_file(path.string());
    std::filesystem::remove(path);
    require(table.observation_altitude_m() == 2000 && table.depth(450, 3000) == 1,
            "table parser and knot interpolation");
    PhotonBatchContext context{7, 8, 0, 0, {0, 0, 0}, 1};
    OpticalPhoton bunch{{{0, 0, 0}, {0, 0, -1}}, 10, 0, 12, 2.5, 0, 3000,
                         std::numeric_limits<double>::quiet_NaN()};
    EventioRunInfo info{300, 600, false, true, false, true, false, 0};
    info.observation_altitude_m = 2000;
    auto children = resolve_eventio_spectrum(bunch, info, 16, 1234);
    auto repeated = resolve_eventio_spectrum(bunch, info, 16, 1234);
    double total = 0;
    for (std::size_t i = 0; i < children.size(); ++i) {
      total += children[i].weight;
      require(children[i].wavelength_nm >= 300 && children[i].wavelength_nm <= 600 &&
                  children[i].wavelength_nm == repeated[i].wavelength_nm &&
                  children[i].photon_id == repeated[i].photon_id,
              "stable in-band reciprocal-wavelength stratification");
    }
    require(std::abs(total - 2.5) < 1e-12, "spectral splitting conserves weight");
    const double removed = attenuate_eventio_direct_beam(children, context, info, table);
    require(std::abs(removed - 2.5 * (1 - std::exp(-1))) < 1e-12,
            "direct extinction uses optical depth once");

    OpticalPhoton spatial{{{0, 0, 0}, {0, 0, -1}}, 11, 450, 12, 1, 0,
                          std::numeric_limits<double>::quiet_NaN(), 1000};
    require(std::abs(table.direct_survival(spatial, context) - std::exp(-1)) < 1e-12,
            "3D emission distance produces expected altitude");
    info.ceffic = true;
    bool rejected = false;
    try { (void) attenuate_eventio_direct_beam(children, context, info, table); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "preselected CEFFIC light cannot be attenuated twice");
    info.ceffic = false;
    info.refraction = true;
    rejected = false;
    try { (void) attenuate_eventio_direct_beam(children, context, info, table); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "refracted rays require a path-integral atmosphere model");
    std::cout << "atmosphere transport tests passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
