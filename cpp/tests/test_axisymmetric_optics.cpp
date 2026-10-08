#include "obdeect/axisymmetric_optics.hpp"
#include "obdeect/ctao_optical_specs.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>

namespace {

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}

} // namespace

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

  const AxisymmetricMirror horizontal_surface{0, 0, 4, paraboloid(1)};
  const auto horizontal = intersect_axisymmetric_mirror({{0, 0, 1}, {1, 0, 0}}, horizontal_surface);
  require(horizontal && std::abs(horizontal->distance_m - 2) < 1.e-10,
          "horizontal ray intersects finite paraboloid");
  const auto multiple = intersect_axisymmetric_mirror({{-3, 0, 1}, {1, 0, 0}}, horizontal_surface);
  require(multiple && std::abs(multiple->distance_m - 1) < 1.e-10,
          "nearest of two horizontal paraboloid roots wins");
  const auto grazing = intersect_axisymmetric_mirror({{-3, 2, 1}, {1, 0, 0}}, horizontal_surface);
  require(grazing && std::abs(grazing->distance_m - 3) < 1.e-9,
          "double root at grazing incidence is retained");
  const auto parallel =
      intersect_axisymmetric_mirror({{0, 0, 1}, {1, 0, -1.e-14}}, horizontal_surface);
  require(parallel && std::abs(parallel->distance_m - 2) < 1.e-10,
          "near horizontal ray does not depend on vertex-plane seed");
  const AxisymmetricMirror holed{0, 1, 4, paraboloid(1)};
  require(!intersect_axisymmetric_mirror({{-3, 0, 0.01}, {1, 0, 0}}, holed),
          "roots in central hole are physical misses");
  auto extreme = horizontal_surface;
  extreme.surface.coefficient_m[12] = std::numeric_limits<double>::max();
  bool numerical_failure = false;
  const auto failed = intersect_axisymmetric_mirror(
      {{1.5, 0, 10}, {0, 0, -1}}, extreme, kEpsilon, [](const Vec3 &) { return true; },
      &numerical_failure);
  require(!failed && numerical_failure,
          "overflow in a supported aperture is explicit numerical failure");
  // T-SC-001: conversion is exact at the reference radius and supplied SC
  // profiles remain finite over their primary apertures.
  const std::array<double, 13> centimetre_profile{0.0, 0.25};
  const auto metre_profile = centimetre_even_polynomial_to_metres(centimetre_profile);
  require(std::abs(metre_profile.sag(1.0) - 25.0) < 1e-12, "cm polynomial conversion at 1 m");
  const auto normalised_profile =
      centimetre_normalised_radius_polynomial_to_metres({0.0, 25.0}, 2.0);
  require(std::abs(normalised_profile.sag(2.0) - 0.25) < 1e-12,
          "normalised-radius conversion at reference radius");
  const auto ssts = ssts_design_reference();
  const auto scts = scts_design_reference();
  require(std::abs(scts.secondary.sag(0.0) - 8.37945) < 1e-12,
          "SCT secondary vertex uses sag scale R, not centimetres or aperture radius");
  const auto scaled = reference_radius_polynomial({1.5, 0.25}, 4.0);
  require(std::abs(scaled.sag(4.0) - 7.0) < 1e-12 &&
              std::abs(scaled.radial_slope(4.0) - 0.5) < 1e-12,
          "reference radius scales both height and radial slope");
  const AxisymmetricMirror offset_plane{0, 0, 2, reference_radius_polynomial({1.5}, 4.0)};
  const auto offset_hit = intersect_axisymmetric_mirror({{0, 0, 1}, {0, 0, 1}}, offset_plane);
  require(offset_hit && std::abs(offset_hit->point_m.z - 6.0) < 1e-12,
          "nonzero constant sag must enter intersection seed");
  require(std::isfinite(ssts.primary.sag(0.5 * ssts.primary_diameter_m)),
          "SSTS primary profile finite");
  require(std::isfinite(scts.secondary.sag(0.5 * scts.secondary_diameter_m)),
          "SCTS secondary profile finite");
  require(kMstNectarCam.family == TelescopeOpticalFamily::mst_modified_davies_cotton,
          "MST reference must retain Davies-Cotton family");

  // T-ASP-003: distant oblique rays retain local intersection precision.
  // The analytic paraboloid contains (1, 0, 0.25), regardless of launch distance.
  for (const double angle_deg : {0.0, 0.5, 1.0, 2.0, 3.0}) {
    const double angle = angle_deg * std::numbers::pi / 180.0;
    const Vec3 direction{std::sin(angle), 0, -std::cos(angle)};
    for (const double launch_distance : {10.0, 1.e4, 1.e7}) {
      const Vec3 target{1, 0, 0.25};
      const Ray distant{target - direction * launch_distance, direction};
      bool failure = false;
      const auto hit = intersect_axisymmetric_mirror(
          distant, horizontal_surface, kEpsilon, [](const Vec3 &) { return true; }, &failure);
      require(hit && !failure, "distant off-axis paraboloid intersection must not fail");
      // Input coordinates at 10,000 km have nanometre-scale representation error.
      require(norm(hit->point_m - target) < 1.e-8,
              "distant launch must preserve the analytic local intersection");
      require(std::abs(hit->distance_m - launch_distance) < 1.e-8,
              "distant intersection must preserve full travel distance");
      require(std::abs(hit->point_m.z - horizontal_surface.surface.sag(
                                            std::hypot(hit->point_m.x, hit->point_m.y))) < 1.e-10,
              "distant intersection retains the original surface residual tolerance");
    }
  }
  // Versioned SSTS reference coefficients reproduce the 7.0.0 distant-ray
  // failure: expansion alone left the local sag residual above 1e-10 m.
  const AxisymmetricMirror sst_primary{0, 0, ssts.primary_diameter_m * 0.5, ssts.primary};
  bool failed_distant_asphere = false;
  const auto distant_asphere = intersect_axisymmetric_mirror(
      {{-524076.55125975102, 0.12236330347260739, 1.e7},
       {0.052335956242943842, 0, -0.99862953475457394}},
      sst_primary, kEpsilon, [](const Vec3 &) { return true; }, &failed_distant_asphere);
  require(distant_asphere && !failed_distant_asphere,
          "high-order distant asphere root must satisfy the unchanged residual tolerance");

  std::cout << "axisymmetric optics tests passed\n";
}
