#include "obdeect/obscurers.hpp"
#include "obdeect/pixel_response.hpp"
#include "obdeect/spatial_response.hpp"

#include "obdeect/artificial_mst.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}

} // namespace

int main() {
  {
    obdeect::TelescopeTransmission response{0.96, 0.2, 0.1, 2, 1};
    require(response.is_valid(), "telescope transmission validates");
    require(*response.at_unchecked({0, 0, -1}) == 0.96, "on-axis structural transmission");
    require(std::abs(*response.at_unchecked({0.1, 0, -std::sqrt(0.99)}) - 0.8) < 1.e-12,
            "angular structural transmission uses sine of incidence");
    response.outer_power = 3;
    require(std::abs(*response.at_unchecked({0.1, 0, -std::sqrt(0.99)}) - 0.96 / 1.728) < 1.e-12,
            "outer transmission exponent");
    response.amplitude = -1;
    require(!response.is_valid(), "invalid transmission coefficient rejects");
  }
  {
    using namespace obdeect;
    const std::array<double, 2> axis{1, 10}, values{0.1, 1};
    TableInterpolation options{TableBoundary::reject};
    require(!interpolate_curve_unchecked(axis, values, 0, options), "table rejects outside range");
    options.boundary = TableBoundary::clamp;
    require(*interpolate_curve_unchecked(axis, values, 0, options) == 0.1,
            "table clamps outside range");
    options.boundary = TableBoundary::zero;
    require(*interpolate_curve_unchecked(axis, values, 0, options) == 0,
            "table clips outside range");
    options.x_log = options.value_log = true;
    require(std::abs(*interpolate_curve_unchecked(axis, values, std::sqrt(10.), options) -
                     std::sqrt(0.1)) < 1.e-12,
            "table logarithmic interpolation");
    options.x_log = options.value_log = false;
    options.scheme = 0;
    require(*interpolate_curve_unchecked(axis, values, 5, options) == 0.1,
            "nearest neighbour lower interval");
    options.scheme = 2;
    options.coefficients = {{0.1, 0, 0.9, 0}};
    require(options.is_valid(1), "compiled polynomial interval validates");
    require(std::abs(*interpolate_curve_unchecked(axis, values, 5.5, options) - 0.325) < 1.e-12,
            "table evaluates compiled polynomial");
    CameraIncidenceResponse invalid;
    require(!invalid.is_valid(), "empty angular response rejects safely");
    invalid.interpolation.x_log = true;
    require(!invalid.is_valid(), "empty logarithmic angular response rejects safely");
  }
  {
    using namespace obdeect;
    PixelResponseTable table;
    table.id = 2;
    table.tangent_bin_width = 0.1;
    table.angular_efficiency = {0.5, 0.8};
    table.wavelength_bin_width_nm = 1;
    table.wavelength_bin_origin_nm = 0.5;
    table.spectral_correction = {0.4, 0.6};
    PixelResponses response{{table}, {{7, 2, {}}}};
    require(response.is_valid(), "measured pixel response validates");
    require(std::abs(*response.at_unchecked(7, 1.1, 0, {}) - 0.3) < 1.e-12,
            "measured spectral integer bin");
    require(*response.at_unchecked(7, 1.1, 45, {}) == 0, "angular table boundary");
    table.method = PixelResponseTable::Method::single_reflection;
    table.angular_efficiency.clear();
    table.spectral_correction.clear();
    table.transparency = 0.9;
    table.wall_reflectivity = 0.7;
    ImportedDetectorSurface cathode;
    cathode.id = 8;
    cathode.diameter_m = 1;
    cathode.unit_normal = {0, 0, 1};
    cathode.unit_tangent_u = {1, 0, 0};
    response = {{table}, {{7, 2, cathode}}};
    require(response.is_valid(), "legacy pixel response validates");
    require(*response.at_unchecked(7, 400, 0, {{0, 0, 0}, {0, 0, -1}}) == 0.9,
            "zero-depth cathode accepts direct light");
    require(std::abs(*response.at_unchecked(7, 400, 0, {{1, 0, 0}, {0, 0, -1}}) - 0.63) < 1.e-12,
            "single reflection applies wall loss");
  }
  {
    obdeect::SpatialResponse map{{-1, 1}, {-1, 1}, {0.2, 0.4, 0.6, 0.8}};
    require(bool(map.is_valid()), "manual optical primitive");
    require(bool(std::abs(*map.at_unchecked(0, 0) - 0.5) < 1.e-12), "manual optical primitive");
    require(bool(std::abs(*map.at_unchecked(2, 0) - 0.7) < 1.e-12), "manual optical primitive");
    map.clip = true;
    require(bool(*map.at_unchecked(2, 0) == 0), "manual optical primitive");
    map.clip = false;
    map.interpolation.boundary = obdeect::TableBoundary::reject;
    require(!map.at_unchecked(2, 0), "spatial map rejects outside range");
    map.response.pop_back();
    require(bool(!map.is_valid()), "manual optical primitive");
  }

  {
    using namespace obdeect;
    OpaqueSurface baffle;
    baffle.shape = OpaqueSurface::Shape::hollow_frustum;
    baffle.first = {0, 0, 1};
    baffle.second = {0, 0, 3};
    baffle.first_radius_m = baffle.second_radius_m = 1;
    baffle.thickness_m = 0.1;
    require(bool(is_valid(baffle)), "manual optical primitive");
    require(bool(!intersect_opaque_surface({{0, 0, 4}, {0, 0, -1}}, baffle)),
            "manual optical primitive");
    const auto end = intersect_opaque_surface({{1.05, 0, 4}, {0, 0, -1}}, baffle);
    require(bool(end && std::abs(end->distance_m - 1) < 1.e-12), "manual optical primitive");
    const auto side = intersect_opaque_surface({{2, 0, 2}, {-1, 0, 0}}, baffle);
    require(bool(side && std::abs(side->distance_m - 0.9) < 1.e-12), "manual optical primitive");
    baffle.second_radius_m = 2;
    baffle.thickness_m = 0;
    const auto cone = intersect_opaque_surface({{3, 0, 2}, {-1, 0, 0}}, baffle);
    require(bool(cone && std::abs(cone->distance_m - 1.5) < 1.e-12), "manual optical primitive");
    baffle.shape = OpaqueSurface::Shape::solid_frustum;
    const auto closed_end = intersect_opaque_surface({{0, 0, 4}, {0, 0, -1}}, baffle);
    require(closed_end && std::abs(closed_end->distance_m - 1) < 1.e-12,
            "solid structural rod has a closed end");
    baffle.thickness_m = 0.1;
    require(!is_valid(baffle), "solid rod rejects hollow wall thickness");
    OpaqueSurface plate;
    plate.vertices = {Vec3{-1, -1, 2}, Vec3{1, -1, 2}, Vec3{1, 1, 2}, Vec3{-1, 1, 2}};
    require(bool(is_valid(plate)), "manual optical primitive");
    require(bool(intersect_opaque_surface({{0, 0, 4}, {0, 0, -1}}, plate)),
            "manual optical primitive");
    require(bool(!intersect_opaque_surface({{2, 0, 4}, {0, 0, -1}}, plate)),
            "manual optical primitive");
    plate.vertices[3].z = 3;
    require(bool(!is_valid(plate)), "manual optical primitive");
  }

  using namespace obdeect;
  ArtificialMstConfig bare{};
  bare.include_structure = false;

  // T-REFERENCE-001: the on-axis parallel ray reaches the paraxial focal screen.
  const auto central = trace_artificial_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 1, bare);
  require(central.status == PhotonStatus::detected,
          "central ray must be detected without structure");
  require(central.point_count == 3, "detected path has source, mirror and screen vertices");
  require(std::abs(central.points_m[2].x) < 1e-12 && std::abs(central.points_m[2].y) < 1e-12,
          "central ray must focus on the optical axis");
  require(central.path_length_m > 20.0, "path length must be positive and accumulated");

  // T-OBS-001: the MST camera is a pre-M1 obstruction on the optical axis.
  ArtificialMstConfig structured{};
  const auto camera = trace_artificial_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 2, structured);
  require(camera.status == PhotonStatus::blocked_camera,
          "camera must shadow the central incoming ray");
  require(camera.point_count == 2, "blocked path terminates at the obstruction");

  // T-OBS-004: camera and supports are closed solids.  This ray clears the
  // incoming camera aperture, reflects from M1, then enters the camera rear.
  const auto post_reflection_camera =
      trace_artificial_mst({{0.6, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 8, structured);
  require(post_reflection_camera.status == PhotonStatus::blocked_camera,
          "camera must also obstruct the reflected segment");
  require(post_reflection_camera.point_count == 3,
          "post-reflection obstruction retains entrance, M1, and obstruction vertices");
  require(std::abs(post_reflection_camera.points_m[2].z -
                   (structured.focal_length_m - structured.camera_half_depth_m)) < 1e-12,
          "reflected ray must stop at the rear camera cap");

  // T-OBS-005: a closed cylinder catches axial hits on its end caps, which a
  // side-wall-only intersection would miss.
  const Ray axial_cylinder{{0.0, 0.0, 2.0}, {0.0, 0.0, -1.0}};
  const auto cap_hit =
      intersect_closed_finite_cylinder(axial_cylinder, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 0.5);
  require(cap_hit && std::abs(*cap_hit - 1.0) < 1e-12,
          "closed cylinder must intersect its end cap");

  // T-OBS-002: an off-axis ray remains traceable through the MST structure.
  const auto outer = trace_artificial_mst({{2.0, 1.0, 20.0}, {0.0, 0.0, -1.0}}, 3, structured);
  require(outer.status == PhotonStatus::detected || outer.status == PhotonStatus::blocked_mast ||
              outer.status == PhotonStatus::blocked_camera,
          "off-axis ray must either reach the screen or hit a physical obstruction");

  // T-SRC-001: artificial blue source emits unit, pupil-bounded directions.
  const auto rays = parallel_blue_cherenkov_rays(256, structured);
  require(rays.size() == 256, "source photon count");
  for (const auto &ray : rays) {
    require(std::abs(norm(ray.direction) - 1.0) < 1e-14, "source direction must be unit length");
    require(ray.position_m.x * ray.position_m.x + ray.position_m.y * ray.position_m.y <=
                structured.mirror_aperture_radius_m * structured.mirror_aperture_radius_m + 1e-12,
            "source samples must stay inside pupil");
  }
  require(parallel_blue_cherenkov_rays(0, structured).empty(), "zero source count must be empty");

  // T-ABI-001: invalid vectors/configuration fail closed rather than divide by zero.
  const auto zero_direction = trace_artificial_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, 0.0}}, 4, bare);
  require(zero_direction.status == PhotonStatus::invalid_input, "zero direction must be rejected");
  ArtificialMstConfig invalid = bare;
  invalid.mirror_aperture_radius_m = invalid.mirror_radius_m + 1.0;
  const auto invalid_optical_model =
      trace_artificial_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 5, invalid);
  require(invalid_optical_model.status == PhotonStatus::invalid_input,
          "impossible spherical cap must be rejected");
  invalid = bare;
  invalid.camera_half_depth_m = invalid.focal_length_m;
  const auto inverted_camera =
      trace_artificial_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 6, invalid);
  require(inverted_camera.status == PhotonStatus::invalid_input,
          "camera support endpoint behind the primary must be rejected");

  // T-OBS-003: zero-radius camera/masts are disabled, not zero-area blockers.
  ArtificialMstConfig no_obstruction = structured;
  no_obstruction.camera_radius_m = 0.0;
  no_obstruction.mast_radius_m = 0.0;
  const auto unobscured =
      trace_artificial_mst({{0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}}, 7, no_obstruction);
  require(unobscured.status == PhotonStatus::detected, "zero-radius obstructions must not block");

  std::cout << "obdeect core tests passed\n";
}
