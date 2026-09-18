#include "obdeect/toy_mst.hpp"

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

  // T-LIT-001. Garczarczyk, "MST MC parameters" (2017), p. 19 gives
  // central facet R = 32.14 m and f = 16.07 m = R/2. This is a one-facet
  // analytic validation only; it is not the 86-facet modified Davies-Cotton
  // telescope described on p. 2 of that document.
  ToyMstConfig central_facet{32.14, 0.6, 16.07, 0.3, 0.0, 0.0, 0.0, false};
  require(std::abs(central_facet.mirror_radius_m / 2.0 - central_facet.focal_length_m) < 1e-12,
          "published central-facet paraxial focus must equal R/2");

  const auto ray = trace_toy_mst({{0.0, 0.0, 50.0}, {0.0, 0.0, -1.0}}, 0, central_facet);
  require(ray.status == PhotonStatus::detected, "published central-facet ray must reach focal plane");
  require(std::abs(ray.points_m[2].x) < 1e-12 && std::abs(ray.points_m[2].y) < 1e-12,
          "published central-facet on-axis ray must land at focal-plane origin");
  std::cout << "literature central-facet validation passed\n";
}
