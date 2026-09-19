#include "obdeect/facets.hpp"

#include <array>
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
  constexpr std::array<double, 2> wavelength_nm{300.0, 500.0};
  constexpr std::array<double, 2> reflectivity{0.80, 0.90};
  const Table1DView coating{wavelength_nm, reflectivity};
  const CircularFacet near{4, {0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}, 1.0, coating};
  const CircularFacet far{9, {0.0, 0.0, 2.0}, {0.0, 0.0, 1.0}, 1.0, coating};

  // T-FACET-001: finite aperture, wavelength response and specular direction.
  const Ray on_axis{{0.0, 0.0, 5.0}, {0.0, 0.0, -1.0}};
  const auto hit = intersect_facet(on_axis, near, 400.0);
  require(hit.has_value() && hit->facet_id == 4, "on-axis ray must hit finite facet");
  require(std::abs(hit->reflectivity - 0.85) < 1.0e-12, "coating must interpolate in wavelength");
  const auto reflected = reflect(on_axis, *hit);
  require(reflected.has_value() && reflected->direction.z > 0.999999999,
          "facet must reflect an on-axis ray upward");

  // T-FACET-002: nearest physical facet wins; misses and out-of-table
  // wavelengths fail closed rather than applying an extrapolated coating.
  const auto nearest = intersect_faceted_mirror(on_axis, {{far, near}}, 400.0);
  require(nearest.has_value() && nearest->facet_id == 9, "nearest facet must win");
  require(!intersect_facet({{2.0, 0.0, 5.0}, {0.0, 0.0, -1.0}}, near, 400.0),
          "ray outside finite facet aperture must miss");
  require(!intersect_facet(on_axis, near, 700.0), "out-of-table wavelength must fail closed");
  const std::array<double, 2> unphysical{1.1, 0.9};
  auto invalid = near;
  invalid.reflectivity = {wavelength_nm, unphysical};
  require(!is_valid(invalid), "mirror coating cannot amplify photon weight");

  std::cout << "facet tests passed\n";
}
