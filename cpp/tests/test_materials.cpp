#include "obdeect/interactions.hpp"
#include <cstdlib>
#include <iostream>

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
}
