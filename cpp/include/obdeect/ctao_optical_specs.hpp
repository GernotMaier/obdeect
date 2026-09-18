#pragma once

#include "obdeect/axisymmetric_optics.hpp"

namespace obdeect {

enum class TelescopeOpticalFamily : std::uint8_t { lst_parabolic, mst_modified_davies_cotton, sst_sc, sct_sc };

struct SingleReflectorReference {
  TelescopeOpticalFamily family;
  double focal_length_m;
  double facet_focal_length_m;
  double dish_sphere_radius_m;
  bool parabolic_dish;
};

// Values pinned from simulation-models 6.3.0.  These records describe the
// dish/facet convention only; actual LST/MST ray tracing must additionally load
// the model's mirror list and camera/structure data.
inline constexpr SingleReflectorReference kLstNorthDesign{
    TelescopeOpticalFamily::lst_parabolic, 28.0, 0.0, 12.5, true};
inline constexpr SingleReflectorReference kMstNectarCam{
    TelescopeOpticalFamily::mst_modified_davies_cotton, 16.0, 16.083, 9.15, false};

struct SchwarzschildCouderReference {
  TelescopeOpticalFamily family;
  double focal_length_m;
  double primary_diameter_m;
  double secondary_diameter_m;
  EvenPolynomialSurface primary;
  EvenPolynomialSurface secondary;
};

// Coefficients are converted from the simulation-models 6.3.0 SSTS/SCTS
// records. They are surfaces only; segment layout, reference radii, holes and
// mirror vertex transforms remain model-import responsibilities.
[[nodiscard]] inline SchwarzschildCouderReference ssts_design_reference() {
  return {TelescopeOpticalFamily::sst_sc,
          2.15,
          4.241,
          1.80,
          centimetre_even_polynomial_to_metres(
              {0.0, 0.000608051, -7.36248e-10, 5.82125e-15, -6.76447e-20, -3.89542e-24,
               -5.28036e-29, 2.99107e-34, 4.39153e-38, 6.17433e-43, -2.73586e-47, 0.0, 0.0}),
          centimetre_even_polynomial_to_metres(
              {310.84, 0.00229358, 2.8273e-08, -2.7689e-12, 8.80066e-17, 3.37315e-21,
               -1.02973e-25, -6.72882e-30, -3.06437e-34, 3.1718e-38, -3.71217e-43, 0.0, 0.0})};
}

[[nodiscard]] inline SchwarzschildCouderReference scts_design_reference() {
  return {TelescopeOpticalFamily::sct_sc,
          5.5863,
          9.6638,
          5.4166,
          centimetre_even_polynomial_to_metres(
              {0.0, 0.111112, -0.00698726, -0.00206487, -0.00689219, 0.0301911, -0.119762,
               0.319791, -0.602077, 0.777846, -0.661167, 0.333439, -0.0764291}),
          centimetre_even_polynomial_to_metres(
              {1.5, 0.416667, 0.145816, -0.712012, 4.17685, -23.1617, 118.844, -520.501,
               1802.98, -4605.3, 8019.1, -8422.9, 4004.61})};
}

}  // namespace obdeect
