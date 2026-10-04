#include "obdeect/trace.hpp"

#include <array>
#include <cassert>

int main() {
  using namespace obdeect;
  const ModelProvenance provenance{"finite-window-fixture", "1", std::string(64, 'a')};
  const OpticalMaterialRecord glass{
      1, {{300, 500}, {1.5, 1.5}}, {{300, 500}, {1.8, 1.8}}, {{300, 500}, {0.2, 0.2}}};
  OpticalSurfaceRecord entry{
      1, {{0, 0, 2}}, FacetShape::circle, 4, SurfaceRole::refractive_interface};
  entry.back_medium_id = 1;
  auto exit = entry;
  exit.id = 2;
  exit.frame.origin_m.z = 1;
  exit.front_medium_id = 1;
  exit.back_medium_id = kVacuumMediumId;
  const OpticalSurfaceRecord detector{3, {}, FacetShape::circle, 4, SurfaceRole::detector};
  const auto window = compile_optical_model({provenance, {entry, exit, detector}, 8, {glass}});
  assert(window);
  const auto normal = trace_material_path({{0, 0, 3}, {0, 0, -1}}, 1, 400, *window);
  assert(normal.status == PhotonStatus::detected && normal.point_count == 4);
  assert(std::abs(normal.path_length_m - 3) < 1e-12);
  assert(std::abs(normal.optical_path_m - 3.5) < 1e-12);
  assert(std::abs(normal.group_delay_ns - 3.8 / kSpeedOfLightMPerNs) < 1e-12);
  assert(std::abs(normal.surviving_throughput - 0.96 * 0.96 * std::exp(-0.2)) < 1e-12);
  assert(normal.interaction_surface_ids[0] == 1 && normal.interaction_surface_ids[1] == 2);

  // Snell's law and lateral displacement of a finite parallel window.
  const auto oblique =
      trace_material_path({{0, 0, 3}, {0.5, 0, -std::sqrt(0.75)}}, 2, 400, *window);
  assert(oblique.status == PhotonStatus::detected);
  const double internal_sine = 0.5 / 1.5;
  const double internal_cosine = std::sqrt(1 - internal_sine * internal_sine);
  assert(std::abs(oblique.points_m[2].x - oblique.points_m[1].x - internal_sine / internal_cosine) <
         1e-12);
  assert(std::abs(oblique.final_direction.x - 0.5) < 1e-12);
  assert(std::abs(oblique.final_direction.z + std::sqrt(0.75)) < 1e-12);

  // Curved finite interfaces use their actual polynomial normals. Index one
  // makes them optically invisible without changing their recorded geometry.
  auto vacuum_material = glass;
  vacuum_material.phase_index.value = {1, 1};
  vacuum_material.group_index.value = {1, 1};
  vacuum_material.absorption_per_m.value = {0, 0};
  entry.sag = paraboloid(4);
  exit.sag = paraboloid(8);
  const auto curved =
      compile_optical_model({provenance, {entry, exit, detector}, 8, {vacuum_material}});
  assert(curved);
  const auto curved_path = trace_material_path({{0.5, 0, 3}, {0, 0, -1}}, 3, 400, *curved);
  assert(curved_path.status == PhotonStatus::detected && curved_path.surviving_throughput == 1);
  assert(std::abs(curved_path.points_m[1].z - (2 + 0.25 / 16)) < 1e-12);
  assert(std::abs(curved_path.points_m[2].z - (1 + 0.25 / 32)) < 1e-12);
  assert(std::abs(curved_path.optical_path_m - 3) < 1e-12);
  const auto refracting_curved =
      compile_optical_model({provenance, {entry, exit, detector}, 8, {glass}});
  assert(refracting_curved);
  const auto refracted = trace_material_path({{0.5, 0, 3}, {0, 0, -1}}, 4, 400, *refracting_curved);
  assert(refracted.status == PhotonStatus::detected);
  const auto &n = refracted.interaction_normals[0];
  const auto &d = refracted.interaction_outgoing_directions[0];
  assert(std::abs(norm(cross(Vec3{0, 0, -1}, n)) - 1.5 * norm(cross(d, n))) < 1e-12);

  // Complete-interface data replace the same analytic Fresnel loss once.
  entry.sag.reset();
  exit.sag.reset();
  entry.transmission = SpectralResponse{{300, 500}, {0.8, 0.8}};
  entry.transmission_semantics = InterfaceTransmissionSemantics::complete_interface;
  const auto measured = compile_optical_model({provenance, {entry, exit, detector}, 8, {glass}});
  assert(measured);
  assert(std::abs(
             trace_material_path({{0, 0, 3}, {0, 0, -1}}, 5, 400, *measured).surviving_throughput -
             0.8 * 0.96 * std::exp(-0.2)) < 1e-12);
  entry.transmission_semantics.reset();
  assert(!compile_optical_model({provenance, {entry, exit}, 8, {glass}}));
  entry.transmission.reset();

  auto wrong_medium = exit;
  wrong_medium.front_medium_id = 2;
  assert(!compile_optical_model({provenance, {entry, wrong_medium}, 8, {glass}}));
  assert(!compile_optical_model({provenance, {entry, exit}, 8, {glass, glass}}));
  auto invalid_glass = glass;
  invalid_glass.group_index.value[0] = 0;
  assert(!compile_optical_model({provenance, {entry, exit}, 8, {invalid_glass}}));
  assert(trace_material_path({{0, 0, 3}, {0, 0, -1}}, 6, 600, *window).status ==
         PhotonStatus::invalid_input);

  // An unclosed finite window cannot silently become vacuum after its last
  // known interface. TIR preserves the medium and all power in a clear glass.
  const auto open_window = compile_optical_model({provenance, {entry}, 8, {glass}});
  assert(open_window);
  assert(trace_material_path({{0, 0, 3}, {0, 0, -1}}, 7, 400, *open_window).status ==
         PhotonStatus::escaped_material);
  auto clear = glass;
  clear.absorption_per_m.value = {0, 0};
  auto tir_exit = exit;
  tir_exit.diameter_m = 20;
  auto upper_detector = detector;
  upper_detector.frame.origin_m.z = 3;
  upper_detector.diameter_m = 20;
  const auto tir_model =
      compile_optical_model({provenance, {tir_exit, upper_detector}, 8, {clear}, 1});
  assert(tir_model);
  const auto tir = trace_material_path({{0, 0, 2}, {std::sqrt(0.75), 0, -0.5}}, 8, 400, *tir_model);
  assert(tir.status == PhotonStatus::detected && tir.final_direction.z == 0.5 &&
         tir.surviving_throughput == 1);

  const std::array<Vec3, 1> position{{{0, 0, 3}}}, direction{{{0, 0, -1}}};
  const std::array<double, 1> wavelength{400}, time{5}, weight{2};
  const std::array<std::uint64_t, 1> ids{9};
  const auto bulk = trace(*window, {position, direction, wavelength, time, weight, ids});
  assert(bulk.photons.status[0] == PhotonStatus::detected);
  assert(std::abs(bulk.photons.optical_path_m[0] - 3.5) < 1e-12);
  assert(std::abs(bulk.photons.time_ns[0] - (5 + 3.8 / kSpeedOfLightMPerNs)) < 1e-12);
  assert(std::abs(bulk.summary.input_weight - bulk.summary.detected_weight -
                  bulk.summary.lost_weight) < 1e-12);

  auto strongly_absorbing = clear;
  strongly_absorbing.absorption_per_m.value = {1e6, 1e6};
  const auto dark_detector =
      compile_optical_model({provenance, {upper_detector}, 8, {strongly_absorbing}, 1});
  assert(dark_detector);
  const auto dark = trace_material_path({{0, 0, 2}, {0, 0, 1}}, 10, 400, *dark_detector);
  assert(dark.status == PhotonStatus::detected && dark.surviving_throughput == 0);

  auto coated_mirror = entry;
  coated_mirror.role = SurfaceRole::mirror;
  coated_mirror.back_medium_id = kVacuumMediumId;
  coated_mirror.reflectivity = SpectralResponse{{300, 500}, {0.7, 0.7}};
  const auto reflective = compile_optical_model({provenance, {coated_mirror, upper_detector}, 8});
  assert(reflective);
  const auto reflected = trace_material_path({{0, 0, 2.5}, {0, 0, -1}}, 11, 400, *reflective);
  assert(reflected.status == PhotonStatus::detected && reflected.surviving_throughput == 0.7);
  assert(reflected.final_direction.z == 1);
}
