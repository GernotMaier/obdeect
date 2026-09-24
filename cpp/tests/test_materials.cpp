#include "obdeect/interactions.hpp"
#include "obdeect/materials.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>

void require(bool ok, const char* message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(1); }
}

int main() {
  using namespace obdeect;
  const auto normal = dielectric_interface({0, 0, -1}, {0, 0, 1}, 1, 1.5);
  require(normal && normal->transmitted, "normal incidence transmits");
  require(std::abs(normal->reflectance_s - 0.04) < 1e-12 &&
          std::abs(normal->reflectance_p - 0.04) < 1e-12, "Fresnel normal incidence");
  const double angle = std::atan(1.5);
  const auto brewster = dielectric_interface(
      {std::sin(angle), 0, -std::cos(angle)}, {0, 0, 1}, 1, 1.5);
  require(brewster && brewster->transmitted && brewster->reflectance_p < 1e-24,
          "Brewster p reflection vanishes");
  require(std::abs(brewster->transmitted->x * 1.5 - std::sin(angle)) < 1e-12,
          "Snell law tangential component");
  const auto tir = dielectric_interface({std::sqrt(0.75), 0, -0.5}, {0, 0, 1}, 1.5, 1);
  require(tir && !tir->transmitted && tir->reflectance_s == 1, "total internal reflection");
  require(!dielectric_interface({0, 0, -1}, {0, 0, -1}, 1, 1.5),
          "incorrect normal orientation rejected");
  require(std::abs(*bulk_transmission(2, 0.5) - std::exp(-1.0)) < 1e-12,
          "Beer Lambert absorption");
  require(!bulk_transmission(-1, 0.5), "negative length rejected");

  const std::vector<double> wavelength_nm{300.0, 500.0};
  const std::vector<double> index{1.5, 1.6};
  const std::vector<double> absorption{0.5, 0.0};
  const std::vector<double> coating_transmission{0.8, 0.9};
  const SlabMaterial glass{{wavelength_nm, index}, {wavelength_nm, absorption},
                           Table1DView{wavelength_nm, coating_transmission}};
  require(glass.is_valid(), "wavelength-dependent slab material is valid");

  const auto window = transmit_through_slab({0, 0, -1}, {0, 0, 1}, {1, 2, 3}, 1.0, 2.0, 400.0, glass);
  require(window.status == SlabTransportStatus::transmitted, "window entry and exit transmit");
  require(std::abs(window.position_m.x - 1.0) < 1e-12 && std::abs(window.position_m.y - 2.0) < 1e-12 &&
              std::abs(window.position_m.z - 1.0) < 1e-12,
          "normal-incidence window exits at the far face");
  require(std::abs(window.direction.x) < 1e-12 && std::abs(window.direction.y) < 1e-12 &&
              std::abs(window.direction.z + 1.0) < 1e-12,
          "parallel slab restores direction");
  require(std::abs(window.geometric_path_m - 2.0) < 1e-12 && std::abs(window.optical_path_m - 3.1) < 1e-12,
          "window records internal geometric and optical paths");
  const double fresnel_t = 1.0 - std::pow((1.55 - 1.0) / (1.55 + 1.0), 2);
  const double expected_t = fresnel_t * fresnel_t * 0.85 * std::exp(-0.5);
  require(std::abs(window.transmission - expected_t) < 1e-12,
          "Fresnel, wavelength coating, and Beer-Lambert losses multiply");

  const auto oblique = transmit_through_slab({0.5, 0, -std::sqrt(0.75)}, {0, 0, 1}, {0, 0, 0},
                                             1.0, 0.01, 500.0, glass);
  require(oblique.status == SlabTransportStatus::transmitted && oblique.position_m.x > 0.0 &&
              std::abs(oblique.direction.x - 0.5) < 1e-12,
          "oblique window laterally displaces but restores the outgoing direction");

  const std::vector<double> opaque{0.0, 0.0};
  const SlabMaterial black_glass{{wavelength_nm, index}, {wavelength_nm, absorption},
                                 Table1DView{wavelength_nm, opaque}};
  require(transmit_through_slab({0, 0, -1}, {0, 0, 1}, {0, 0, 0}, 1.0, 0.01, 400.0, black_glass).status ==
              SlabTransportStatus::absorbed,
          "zero transmission is an attributed material loss");
  require(transmit_through_slab({std::sqrt(0.75), 0, -0.5}, {0, 0, 1}, {0, 0, 0}, 1.8, 0.01, 300.0, glass).status ==
              SlabTransportStatus::total_internal_reflection,
          "total internal reflection has an explicit material status");
  require(transmit_through_slab({0, 0, -1}, {0, 0, 1}, {0, 0, 0}, 1.0, 0.01, 700.0, glass).status ==
              SlabTransportStatus::wavelength_out_of_range,
          "unmapped wavelength is rejected explicitly");
}
