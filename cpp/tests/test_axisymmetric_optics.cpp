#include "obdeect/axisymmetric_optics.hpp"
#include "obdeect/ctao_optical_specs.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

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

  // T-ASP-001: an ideal paraboloid sends an on-axis parallel ray to f.
  constexpr double focal_length_m = 28.0;
  const AxisymmetricMirror lst_primary{0.0, 0.0, 12.0, paraboloid(focal_length_m)};
  const Ray input{{4.0, 0.0, 60.0}, {0.0, 0.0, -1.0}};
  const auto primary_hit = intersect_axisymmetric_mirror(input, lst_primary);
  require(primary_hit.has_value(), "paraboloid must intersect incoming ray");
  const auto reflected = reflect(input, *primary_hit);
  require(reflected.has_value(), "paraboloid reflection must be defined");
  const auto focal_t = intersect_plane_z(*reflected, focal_length_m);
  require(focal_t.has_value(), "reflected paraboloid ray must cross focal plane");
  const Vec3 focal_hit = reflected->position_m + reflected->direction * *focal_t;
  require(std::hypot(focal_hit.x, focal_hit.y) < 1e-10, "ideal paraboloid must focus on axis");

  // T-ASP-002: aperture holes and invalid profiles fail closed.
  const AxisymmetricMirror annular{0.0, 1.0, 12.0, paraboloid(focal_length_m)};
  require(!intersect_axisymmetric_mirror({{0.0, 0.0, 60.0}, {0.0, 0.0, -1.0}}, annular),
          "central hole must reject ray");
  const AxisymmetricMirror invalid{0.0, 2.0, 1.0, paraboloid(focal_length_m)};
  require(!intersect_axisymmetric_mirror(input, invalid), "invalid aperture must reject ray");

  // T-SC-001: conversion is exact at the reference radius and supplied SC
  // profiles remain finite over their primary apertures.
  const std::array<double, 13> centimetre_profile{0.0, 0.25};
  const auto metre_profile = centimetre_even_polynomial_to_metres(centimetre_profile);
  require(std::abs(metre_profile.sag(1.0) - 25.0) < 1e-12, "cm polynomial conversion at 1 m");
  const auto normalised_profile = centimetre_normalised_radius_polynomial_to_metres({0.0, 25.0}, 2.0);
  require(std::abs(normalised_profile.sag(2.0) - 0.25) < 1e-12,
          "normalised-radius conversion at reference radius");
  const auto ssts = ssts_design_reference();
  const auto scts = scts_design_reference();
  require(std::isfinite(ssts.primary.sag(0.5 * ssts.primary_diameter_m)), "SSTS primary profile finite");
  require(std::isfinite(scts.secondary.sag(0.5 * scts.secondary_diameter_m)), "SCTS secondary profile finite");
  require(kMstNectarCam.family == TelescopeOpticalFamily::mst_modified_davies_cotton,
          "MST reference must retain Davies-Cotton family");

  std::cout << "axisymmetric optics tests passed\n";
}
