#include "obdeect/detector_assignment.hpp"
#include "obdeect/detector_planes.hpp"
#include "obdeect/optical_response.hpp"

#include <cstdlib>
#include <iostream>
#include <new>

namespace {
bool count_allocations = false;
std::size_t allocations = 0;
void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}
} // namespace

void *operator new(std::size_t size) {
  if (count_allocations)
    ++allocations;
  if (void *pointer = std::malloc(size ? size : 1))
    return pointer;
  throw std::bad_alloc();
}
void operator delete(void *pointer) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::size_t) noexcept { std::free(pointer); }

int main() {
  {
    using namespace obdeect;
    std::vector<ImportedDetectorSurface> surfaces{
        {17, {0, 0, 1}, {0, 0, 1}, 2, FacetShape::square, {1, 0, 0}, 1},
        {18, {0, 0, 2}, {0, 0, 1}, 2, FacetShape::square, {1, 0, 0}, 1}};
    DetectorAssignmentGrid description{4, 4, -2, 2, -2, 2, {1, 0, 0}, {0, 1, 0}, 0};
    const auto grid = CompiledDetectorAssignmentGrid::compile(description, surfaces);
    require(grid.has_value(), "explicit assignment grid compiles");
    const Ray ray{{0, 0, 3}, {0, 0, -1}};
    const auto nearest = compile_detector_planes(surfaces)->intersect(ray);
    require(nearest->surface_id == 18, "physical intersection selects nearest plane");
    surfaces.clear();
    allocations = 0;
    count_allocations = true;
    const auto assigned = grid->intersect(ray, {0, 0, 0});
    const auto rejected = grid->intersect(ray, {-1.9, -1.9, 0});
    count_allocations = false;
    require(assigned && assigned->surface_id == 17 && assigned->distance_m == 2,
            "assignment keeps declared candidate order and immutable surface data");
    require(!rejected, "projected candidate cell can reject a physical intersection");
    require(allocations == 0, "assignment path does not allocate");

    const std::vector<ImportedDetectorSurface> boundary_surface{
        {19, {1, 1, 1}, {0, 0, 1}, 2, FacetShape::square, {1, 0, 0}, 1}};
    const DetectorAssignmentGrid boundary_description{2, 2, 0, 2, 0, 2, {1, 0, 0}, {0, 1, 0}, 0};
    const auto boundary_grid =
        CompiledDetectorAssignmentGrid::compile(boundary_description, boundary_surface);
    require(boundary_grid.has_value(), "grid containing the aperture boundary compiles");
    require(boundary_grid->intersect({{2, 2, 0}, {0, 0, 1}}, {2, 2, 0}).has_value(),
            "upper grid corner clamps into the final cell and accepts an aperture edge hit");
  }
  using namespace obdeect;
  const ImportedDetectorSurface left{11, {-0.6, 0, 3}, {0, 0, 1}, 1, FacetShape::square, {1, 0, 0}};
  const ImportedDetectorSurface right{12, {0.6, 0, 3}, {0, 0, 1}, 1, FacetShape::square, {1, 0, 0}};
  const auto pixels = compile_detector_planes({left, right});
  require(pixels.has_value(), "compile finite pixel entrances");
  require(!pixels->intersect({{0, 0, 0}, {0, 0, 1}}), "physical pixel gap remains a miss");
  const auto hit = pixels->intersect({{-0.6, 0, 0}, {0, 0, 2}});
  require(hit && hit->surface_id == 11 && hit->distance_m == 3,
          "finite entrance returns actual pixel ID and physical path length");
  require(!pixels->intersect({{-0.6, 0, 0}, {1, 0, 0}}), "parallel rays do not hit");
  require(!pixels->intersect({{-0.6, 0, 4}, {0, 0, 1}}), "surfaces behind launch do not hit");

  const double q = std::sqrt(0.5);
  const ImportedDetectorSurface tilted{21, {2, -1, 5},         {q, 0, q},
                                       1,  FacetShape::square, {q, 0, -q}};
  const auto inclined = compile_detector_planes({tilted});
  const Ray axial{{2, -1, 0}, {0, 0, 1}};
  const auto inclined_hit = inclined->intersect(axial);
  require(inclined_hit && inclined_hit->distance_m == 5 && inclined_hit->unit_normal.x == q,
          "translated tilted plane intersection uses supplied frame");
  const Vec3 edge = tilted.centre_m + tilted.unit_tangent_u * 0.49;
  require(inclined->intersect({edge - tilted.unit_normal * 2, tilted.unit_normal}).has_value(),
          "tilted local aperture includes a point just inside its square edge");
  const Vec3 outside = tilted.centre_m + tilted.unit_tangent_u * 0.51;
  require(!inclined->intersect({outside - tilted.unit_normal * 2, tilted.unit_normal}),
          "tilted local aperture excludes a point beyond its square edge");
  auto nearer = tilted;
  nearer.id = 22;
  nearer.centre_m.z = 4;
  auto coincident = tilted;
  coincident.id = 20;
  const auto overlap = compile_detector_planes({tilted, coincident});
  require(overlap->intersect(axial)->surface_id == 20,
          "equal-distance tie uses stable smallest ID");
  const auto depth = compile_detector_planes({tilted, nearer});
  require(depth->intersect(axial)->surface_id == 22, "nearest physical entrance wins");

  std::vector<ImportedDetectorSurface> grid;
  for (std::uint32_t id = 0; id < 256; ++id)
    grid.push_back({id,
                    {static_cast<double>(id % 16), static_cast<double>(id / 16), 3},
                    {0, 0, 1},
                    0.8,
                    static_cast<FacetShape>(id % 4),
                    {1, 0, 0}});
  const auto indexed = compile_detector_planes(grid);
  std::reverse(grid.begin(), grid.end());
  const auto reordered = compile_detector_planes(grid);
  count_allocations = true;
  for (std::uint32_t id = 0; id < 256; ++id) {
    for (const double offset : {0.0, 0.35, 0.45, 0.65}) {
      const Ray ray{
          {static_cast<double>(id % 16) + offset, static_cast<double>(id / 16) + offset, 0},
          {0, 0, 1}};
      const auto accelerated = indexed->intersect(ray);
      const auto shuffled = reordered->intersect(ray);
      std::optional<DetectorSurfaceHit> expected;
      for (const auto &surface : grid) {
        const auto candidate = intersect_detector_surface_unchecked(ray, surface);
        if (candidate && (!expected || candidate->distance_m < expected->distance_m ||
                          (candidate->distance_m == expected->distance_m &&
                           candidate->surface_id < expected->surface_id)))
          expected = candidate;
      }
      require(accelerated.has_value() == expected.has_value() &&
                  shuffled.has_value() == expected.has_value(),
              "BVH agrees with exhaustive apertures");
      if (expected)
        require(accelerated->surface_id == expected->surface_id &&
                    shuffled->surface_id == expected->surface_id &&
                    accelerated->distance_m == expected->distance_m,
                "aperture shape and ordering preserve exact indexed hit");
    }
  }
  count_allocations = false;
  require(allocations == 0, "per-photon entrance traversal does not allocate");
  require(!compile_detector_planes({}), "empty entrance set fails closed");
  require(!compile_detector_planes({left, left}), "duplicate component IDs fail closed");
  auto invalid = left;
  invalid.unit_tangent_u = {0, 0, 1};
  require(!compile_detector_planes({invalid}), "nonorthogonal entrance frame fails closed");
  invalid = left;
  invalid.id = std::numeric_limits<std::uint32_t>::max();
  require(!compile_detector_planes({invalid}), "reserved component ID fails closed");
  const CameraResponse response{0.8, SpectralResponse{{300, 500}, {0.4, 0.6, 0.8, 1.0}, {0, 60}},
                                CameraIncidenceResponse{{0, 60}, {0.9, 0.3}}};
  require(response.is_valid(), "complete measured camera response validates");
  const auto throughput = response.at(400, 30);
  require(throughput && std::abs(*throughput - 0.8 * 0.7 * 0.6) < 1e-15,
          "bilinear filter and linear guide multiply scalar transmission once");
  require(response.at(300, 0) == 0.8 * 0.4 * 0.9 && response.at(500, 60) == 0.8 * 1.0 * 0.3,
          "response includes exact wavelength and incidence table endpoints");
  require(!response.at(299, 30) && !response.at(400, 61) && !response.at(0, 30),
          "response never extrapolates spectral or angular data");
  const CameraResponse spectral_only{0.5, SpectralResponse{{300, 500}, {0.2, 0.8}}, {}};
  require(spectral_only.at(400, 45) == 0.25, "one-dimensional filter needs no angle grid");
  const CameraResponse zero{0, {}, {}};
  require(zero.at(400, 45) == 0, "zero optical response remains an explicit zero");
  auto invalid_response = response;
  invalid_response.lightguide_efficiency->response[0] = 1.1;
  require(!invalid_response.is_valid() && !invalid_response.at(400, 30),
          "guide efficiencies above unity fail closed");
  invalid_response = response;
  invalid_response.lightguide_efficiency->incidence_angle_deg[1] = 91;
  require(!invalid_response.is_valid(), "guide angle above 90 degrees fails closed");
  invalid_response = response;
  invalid_response.camera_transmission = -0.1;
  require(!invalid_response.is_valid(), "negative scalar response fails closed");
  count_allocations = true;
  for (std::size_t index = 0; index < 1024; ++index)
    require(response.at_unchecked(400, 30).has_value(), "valid unchecked response query");
  count_allocations = false;
  require(allocations == 0, "per-photon camera response never allocates");
}
