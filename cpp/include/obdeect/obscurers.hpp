#pragma once

#include "obdeect/surfaces.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace obdeect {

// Path segments on which an obscurer absorbs light.
enum Leg : std::uint8_t {
  incoming = 1,        // Source to primary mirror.
  after_primary = 2,   // Primary to secondary mirror or focal surface.
  after_secondary = 4, // Secondary mirror to focal surface.
};

struct Plane {
  Vec3 normal; // Outward unit normal; the solid is dot(normal, x) <= offset.
  double offset{};
};

// An opaque body. Only the fields of its `kind` are used.
struct Obscurer {
  enum class Kind : std::uint8_t { cylinder, polyhedron, plate, frustum };
  std::uint32_t id{};
  std::uint8_t legs{};
  Kind kind{};
  // Cylinder and frustum axis from `a` to `b`.
  Vec3 a, b;
  double radius_a{}, radius_b{}, thickness{}; // Frustum radii; cylinder uses radius_a.
  bool solid{};                               // Frustum: solid cone instead of a wall.
  // Plate.
  Vec3 centre, normal, u, v;
  double half_width{};
  Aperture shape{};
  // Polyhedron faces, indices into the model's plane list.
  std::uint32_t plane_begin{}, plane_end{};
  // Bounding capsule: every point lies within `bound_radius` of the segment
  // from `bound_a` to `bound_b`.
  Vec3 bound_a, bound_b;
  double bound_radius{};
  // Optional axis-aligned box that also encloses the body (a polyhedron's
  // vertices), tighter than the capsule for slender beams.
  Vec3 box_low, box_high;
  bool has_box{};

  [[nodiscard]] std::optional<double> intersect(const Ray &ray, const std::vector<Plane> &planes,
                                                double t_min) const {
    switch (kind) {
    case Kind::cylinder:
      return intersect_cylinder(ray, t_min);
    case Kind::polyhedron:
      return intersect_polyhedron(ray, planes, t_min);
    case Kind::plate: {
      const double denominator = dot(ray.direction, normal);
      if (denominator == 0.0)
        return std::nullopt;
      const double t = dot(centre - ray.origin, normal) / denominator;
      if (!(t > t_min))
        return std::nullopt;
      const Vec3 offset = ray.at(t) - centre;
      return inside_aperture(shape, half_width, dot(offset, u), dot(offset, v))
                 ? std::optional<double>{t}
                 : std::nullopt;
    }
    case Kind::frustum:
      return intersect_frustum(ray, t_min);
    }
    return std::nullopt;
  }

private:
  // Closed cylinder: side wall plus both end caps.
  [[nodiscard]] std::optional<double> intersect_cylinder(const Ray &ray, double t_min) const {
    const Vec3 axis = b - a, offset = ray.origin - a;
    const double axis2 = dot(axis, axis), da = dot(ray.direction, axis), oa = dot(offset, axis);
    const double qa = 1.0 - da * da / axis2;
    const double qb = dot(ray.direction, offset) - da * oa / axis2;
    const double qc = dot(offset, offset) - oa * oa / axis2 - radius_a * radius_a;
    double best = std::numeric_limits<double>::infinity();
    if (qa > 0.0) {
      const double disc = qb * qb - qa * qc;
      if (disc >= 0.0) {
        const double root = std::sqrt(disc);
        for (const double t : {(-qb - root) / qa, (-qb + root) / qa}) {
          const double along = oa + t * da;
          if (t > t_min && along >= 0.0 && along <= axis2)
            best = std::min(best, t);
        }
      }
    }
    if (da != 0.0)
      for (const double cap : {0.0, axis2}) {
        const double t = (cap - oa) / da;
        const Vec3 radial = offset + ray.direction * t - axis * (cap / axis2);
        if (t > t_min && dot(radial, radial) <= radius_a * radius_a)
          best = std::min(best, t);
      }
    return std::isfinite(best) ? std::optional<double>{best} : std::nullopt;
  }

  // Convex polyhedron by slab clipping; a ray starting inside is absorbed at once.
  [[nodiscard]] std::optional<double>
  intersect_polyhedron(const Ray &ray, const std::vector<Plane> &planes, double t_min) const {
    double enter = -std::numeric_limits<double>::infinity(),
           leave = std::numeric_limits<double>::infinity();
    for (std::uint32_t i = plane_begin; i < plane_end; ++i) {
      const double denominator = dot(planes[i].normal, ray.direction);
      const double distance = planes[i].offset - dot(planes[i].normal, ray.origin);
      if (denominator == 0.0) {
        if (distance < 0.0)
          return std::nullopt;
      } else if (denominator > 0.0) {
        leave = std::min(leave, distance / denominator);
      } else {
        enter = std::max(enter, distance / denominator);
      }
      if (enter > leave)
        return std::nullopt;
    }
    if (!(leave > t_min))
      return std::nullopt;
    return std::max(enter, t_min);
  }

  // Cone frustum between radii radius_a (at a) and radius_b (at b): an open
  // wall of radial `thickness` with annular end faces, or a solid cone.
  [[nodiscard]] std::optional<double> intersect_frustum(const Ray &ray, double t_min) const {
    const Vec3 delta = b - a;
    const double length = norm(delta);
    const Vec3 axis = delta * (1.0 / length), offset = ray.origin - a;
    const double pa = dot(offset, axis), da = dot(ray.direction, axis);
    const Vec3 rp = offset - axis * pa, rd = ray.direction - axis * da;
    const double slope = (radius_b - radius_a) / length;
    double best = std::numeric_limits<double>::infinity();
    const auto consider = [&](double t) {
      if (t > t_min)
        best = std::min(best, t);
    };
    for (const double extra : {0.0, thickness}) {
      if (solid && extra > 0.0)
        break;
      const double r0 = radius_a + extra;
      const double qa = dot(rd, rd) - slope * slope * da * da;
      const double qb = dot(rd, rp) - (r0 + slope * pa) * slope * da;
      const double qc = dot(rp, rp) - (r0 + slope * pa) * (r0 + slope * pa);
      const auto wall = [&](double t) {
        const double along = pa + t * da;
        if (along >= 0.0 && along <= length)
          consider(t);
      };
      if (std::abs(qa) > 1e-15) {
        const double disc = qb * qb - qa * qc;
        if (disc >= 0.0) {
          wall((-qb - std::sqrt(disc)) / qa);
          wall((-qb + std::sqrt(disc)) / qa);
        }
      } else if (qb != 0.0) {
        wall(-qc / (2.0 * qb));
      }
      if (thickness == 0.0)
        break;
    }
    if ((solid || thickness > 0.0) && da != 0.0)
      for (const double along : {0.0, length}) {
        const double t = (along - pa) / da;
        const Vec3 radial = rp + rd * t;
        const double r = radius_a + slope * along, r2 = dot(radial, radial);
        if ((solid || r2 >= r * r) && r2 <= (r + thickness) * (r + thickness))
          consider(t);
      }
    return std::isfinite(best) ? std::optional<double>{best} : std::nullopt;
  }
};

// Bounding capsule of a polyhedron: along its longest vertex-to-vertex
// distance, wide enough for every vertex.
inline void bound_vertices(Obscurer &o, const std::vector<Vec3> &vertices) {
  Vec3 first = vertices.at(0), second = vertices[0];
  for (const auto &p : vertices)
    for (const auto &q : vertices)
      if (norm(q - p) > norm(second - first)) {
        first = p;
        second = q;
      }
  const Vec3 axis = norm(second - first) > 0.0 ? normalized(second - first) : Vec3{0, 0, 1};
  double low = 0.0, high = 0.0;
  o.bound_radius = 0.0;
  for (const auto &vertex : vertices) {
    const double along = dot(vertex - first, axis);
    low = std::min(low, along);
    high = std::max(high, along);
    o.bound_radius = std::max(o.bound_radius, norm(vertex - first - axis * along));
  }
  o.bound_a = first + axis * low;
  o.bound_b = first + axis * high;
  o.box_low = o.box_high = vertices.at(0);
  for (const auto &vertex : vertices) {
    o.box_low = {std::min(o.box_low.x, vertex.x), std::min(o.box_low.y, vertex.y),
                 std::min(o.box_low.z, vertex.z)};
    o.box_high = {std::max(o.box_high.x, vertex.x), std::max(o.box_high.y, vertex.y),
                  std::max(o.box_high.z, vertex.z)};
  }
  o.has_box = true;
}

// The obscurers of an optical model, with a bounding-box hierarchy per leg so
// that a ray tests only the obscurers near its path. Each bounding capsule is
// cut into short pieces, so that long, slanted masts get tight boxes.
// Immutable once built.
class Obscurers {
public:
  struct Hit {
    double t;
    std::uint32_t id;
  };

  Obscurers() = default;
  Obscurers(std::vector<Obscurer> items, std::vector<Plane> planes)
      : items_(std::move(items)), planes_(std::move(planes)) {
    for (std::size_t leg = 0; leg < trees_.size(); ++leg) {
      auto &tree = trees_[leg];
      for (std::uint32_t i = 0; i < items_.size(); ++i)
        if (items_[i].legs & (1U << leg))
          add_pieces(tree, i, kPieceLengthRadii);
      if (!tree.pieces.empty())
        build(tree, 0, static_cast<std::uint32_t>(tree.pieces.size()), 0);
    }
    build_grid();
  }

#ifdef OBDEECT_BENCH_COUNTERS
  // Benchmark-only counters (single thread); never compiled into the library.
  struct Counters {
    std::uint64_t node_visits{}, exact_tests{};
  };
  static Counters &counters() {
    static Counters c;
    return c;
  }
  static std::string counter_report(std::size_t photons) {
    const auto &c = counters();
    return "node visits " + std::to_string(static_cast<double>(c.node_visits) / photons) +
           ", exact tests " + std::to_string(static_cast<double>(c.exact_tests) / photons) +
           " (per nearest() call, all runs)";
  }
#endif

  [[nodiscard]] const std::vector<Obscurer> &items() const { return items_; }
  [[nodiscard]] const std::vector<Plane> &planes() const { return planes_; }

  // Nearest obscurer on `leg` with kMinDistanceM < t < t_max; the lowest id
  // wins exact ties, so the result does not depend on the hierarchy.
  [[nodiscard]] std::optional<Hit> nearest(const Ray &ray, Leg leg, double t_max) const {
    if (leg == Leg::incoming && grid_.nx > 0 && std::abs(ray.direction.z) >= kGridMinDz)
      return nearest_in_grid(ray, t_max);
    return nearest_in_hierarchy(ray, leg, t_max);
  }

  // The same through the bounding-box hierarchy alone (the reference for the grid).
  [[nodiscard]] std::optional<Hit> nearest_in_hierarchy(const Ray &ray, Leg leg,
                                                        double t_max) const {
    const auto &tree = trees_[leg == Leg::incoming ? 0 : leg == Leg::after_primary ? 1 : 2];
    std::optional<Hit> best;
    if (tree.nodes.empty())
      return best;
    const Vec3 inverse{1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z};
    std::array<std::uint32_t, kMaxDepth + 2> stack{};
    std::size_t top = 0;
    stack[top++] = 0;
    // Items already tested by this call (pieces of one obscurer often share
    // a leaf or neighbouring leaves), used as a ring.
    std::array<std::uint32_t, kMailbox> tested;
    tested.fill(kNone);
    std::size_t next_slot = 0;
    while (top > 0) {
      const Node &node = tree.nodes[stack[--top]];
#ifdef OBDEECT_BENCH_COUNTERS
      ++counters().node_visits;
#endif
      if (!enters(node.low, node.high, ray, inverse, best ? best->t : t_max))
        continue;
      if (node.count == 0) {
        // Push the far child first so the near one is popped first. The left
        // child follows its parent; the right one is at `first`.
        const auto left = static_cast<std::uint32_t>(&node - tree.nodes.data()) + 1;
        const bool left_near = component(ray.direction, node.axis) >= 0.0;
        stack[top++] = left_near ? node.first : left;
        stack[top++] = left_near ? left : node.first;
        continue;
      }
      for (std::uint32_t k = node.first; k < node.first + node.count; ++k) {
        const Piece &piece = tree.pieces[k];
        if (seen(tested, piece.item) ||
            !enters(piece.low, piece.high, ray, inverse, best ? best->t : t_max))
          continue;
        tested[next_slot] = piece.item;
        next_slot = (next_slot + 1) % kMailbox;
#ifdef OBDEECT_BENCH_COUNTERS
        ++counters().exact_tests;
#endif
        const Obscurer &o = items_[piece.item];
        if (const auto t = o.intersect(ray, planes_, kMinDistanceM))
          if (*t < t_max && (!best || *t < best->t || (*t == best->t && o.id < best->id)))
            best = Hit{*t, o.id};
      }
    }
    return best;
  }

private:
  static constexpr std::uint32_t kNone = 0xFFFFFFFFU;
  static constexpr std::size_t kMailbox = 8;      // Items remembered per nearest() call.
  static constexpr std::size_t kMaxDepth = 64;    // Tree depth bound (stack size).
  static constexpr std::size_t kSahDepth = 32;    // Deeper levels split at the median.
  static constexpr std::uint32_t kLeafPieces = 4; // Largest leaf.
  static constexpr double kGridCellM = 0.125;     // Cell size of the incoming-leg grid.
  static constexpr double kGridMinDz = 0.5;       // Rays steeper than this use the grid.
  static constexpr double kGridPieceRadii = 4.0;  // Piece length for the grid, in radii.
  static constexpr std::size_t kGridMaxCells = 1U << 20;
  static constexpr std::size_t kGridMinItems = 16; // Fewer items: the hierarchy is faster.
  static constexpr std::size_t kBins = 12;         // Bins of the SAH split.
  static constexpr double kPieceLengthRadii = 2.0; // Longest piece, in capsule radii.
  struct Piece {
    Vec3 low, high;
    std::uint32_t item{};
  };
  // Leaf: `count` pieces from pieces[first]. Inner node (count 0): left child
  // at the next index, right child at `first`.
  struct Node {
    Vec3 low, high;
    std::uint32_t first{}, count{};
    std::uint8_t axis{}; // Split axis of an inner node.
  };
  struct Tree {
    std::vector<Node> nodes;
    std::vector<Piece> pieces;
  };

  // Pieces of the capsule no longer than a few capsule radii.
  void add_pieces(Tree &tree, std::uint32_t item, double piece_radii) {
    const Obscurer &o = items_[item];
    const double r = o.bound_radius, length = norm(o.bound_b - o.bound_a);
    const auto count = static_cast<std::uint32_t>(
        std::clamp(std::ceil(length / (piece_radii * std::max(r, 1e-6))), 1.0, 64.0));
    for (std::uint32_t i = 0; i < count; ++i) {
      const Vec3 p = o.bound_a + (o.bound_b - o.bound_a) * (static_cast<double>(i) / count);
      const Vec3 q = o.bound_a + (o.bound_b - o.bound_a) * (static_cast<double>(i + 1) / count);
      Piece piece{Vec3{std::min(p.x, q.x) - r, std::min(p.y, q.y) - r, std::min(p.z, q.z) - r},
                  Vec3{std::max(p.x, q.x) + r, std::max(p.y, q.y) + r, std::max(p.z, q.z) + r},
                  item};
      if (o.has_box) { // Both boxes enclose the body, so does their intersection.
        constexpr double margin = 1e-9;
        const Vec3 low{std::max(piece.low.x, o.box_low.x - margin),
                       std::max(piece.low.y, o.box_low.y - margin),
                       std::max(piece.low.z, o.box_low.z - margin)};
        const Vec3 high{std::min(piece.high.x, o.box_high.x + margin),
                        std::min(piece.high.y, o.box_high.y + margin),
                        std::min(piece.high.z, o.box_high.z + margin)};
        if (low.x <= high.x && low.y <= high.y && low.z <= high.z) {
          piece.low = low;
          piece.high = high;
        } else {
          piece.low = o.box_low - Vec3{margin, margin, margin};
          piece.high = o.box_high + Vec3{margin, margin, margin};
        }
      }
      tree.pieces.push_back(piece);
    }
  }

  // Uniform x-y grid over the incoming-leg obscurers: each cell lists the items
  // with a capsule piece whose x-y box overlaps it.
  struct Grid {
    double x0{}, y0{}, cell{kGridCellM}, z_low{}, z_high{};
    std::size_t nx{}, ny{};
    std::vector<std::uint32_t> offsets;
    std::vector<Piece> members; // Pieces by cell, each with its box.
  };

  void build_grid() {
    Tree pieces;
    std::size_t incoming = 0;
    for (std::uint32_t i = 0; i < items_.size(); ++i)
      if (items_[i].legs & Leg::incoming) {
        add_pieces(pieces, i, kGridPieceRadii);
        ++incoming;
      }
    if (incoming < kGridMinItems)
      return;
    Vec3 lo = pieces.pieces[0].low, hi = pieces.pieces[0].high;
    for (const Piece &p : pieces.pieces) {
      lo = {std::min(lo.x, p.low.x), std::min(lo.y, p.low.y), std::min(lo.z, p.low.z)};
      hi = {std::max(hi.x, p.high.x), std::max(hi.y, p.high.y), std::max(hi.z, p.high.z)};
    }
    Grid g;
    g.x0 = lo.x;
    g.y0 = lo.y;
    g.z_low = lo.z;
    g.z_high = hi.z;
    // Enlarge the cells if the structure is so wide that the grid would be huge.
    while (std::ceil((hi.x - lo.x) / g.cell) * std::ceil((hi.y - lo.y) / g.cell) >
           static_cast<double>(kGridMaxCells))
      g.cell *= 2.0;
    g.nx = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil((hi.x - lo.x) / g.cell)));
    g.ny = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil((hi.y - lo.y) / g.cell)));
    std::vector<std::vector<Piece>> cells(g.nx * g.ny);
    for (const Piece &p : pieces.pieces)
      for (std::size_t ix = cell_index(p.low.x, g.x0, g.cell, g.nx);
           ix <= cell_index(p.high.x, g.x0, g.cell, g.nx); ++ix)
        for (std::size_t iy = cell_index(p.low.y, g.y0, g.cell, g.ny);
             iy <= cell_index(p.high.y, g.y0, g.cell, g.ny); ++iy) {
          cells[ix * g.ny + iy].push_back(p);
        }
    g.offsets.assign(1, 0);
    for (const auto &c : cells) {
      g.members.insert(g.members.end(), c.begin(), c.end());
      g.offsets.push_back(static_cast<std::uint32_t>(g.members.size()));
    }
    grid_ = std::move(g);
  }

  [[nodiscard]] static std::size_t cell_index(double value, double origin, double cell,
                                              std::size_t count) {
    const double c = (value - origin) / cell;
    if (!(c > 0.0))
      return 0;
    return std::min(count - 1, static_cast<std::size_t>(c));
  }

  // Nearest hit among the items in the cells crossed by the ray between the
  // lowest and highest height of the structure; nothing else can be hit.
  [[nodiscard]] std::optional<Hit> nearest_in_grid(const Ray &ray, double t_max) const {
    const Grid &g = grid_;
    const Vec3 inverse{1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z};
    // Only the cells to visit depend on these heights, so they carry a margin.
    constexpr double margin = 1e-9;
    const double ta = (g.z_low - margin - ray.origin.z) * inverse.z,
                 tb = (g.z_high + margin - ray.origin.z) * inverse.z;
    const double t0 = std::max(std::min(ta, tb), kMinDistanceM),
                 t1 = std::min(std::max(ta, tb), t_max);
    std::optional<Hit> best;
    if (!(t0 <= t1))
      return best;
    const Vec3 a = ray.at(t0), b = ray.at(t1);
    const double x_max = g.x0 + g.cell * static_cast<double>(g.nx),
                 y_max = g.y0 + g.cell * static_cast<double>(g.ny);
    if (std::max(a.x, b.x) < g.x0 || std::min(a.x, b.x) > x_max || std::max(a.y, b.y) < g.y0 ||
        std::min(a.y, b.y) > y_max)
      return best;
    std::array<std::uint32_t, kMailbox> tested;
    tested.fill(kNone);
    std::size_t next_slot = 0;
    for (std::size_t ix = cell_index(std::min(a.x, b.x), g.x0, g.cell, g.nx);
         ix <= cell_index(std::max(a.x, b.x), g.x0, g.cell, g.nx); ++ix)
      for (std::size_t iy = cell_index(std::min(a.y, b.y), g.y0, g.cell, g.ny);
           iy <= cell_index(std::max(a.y, b.y), g.y0, g.cell, g.ny); ++iy)
        for (std::uint32_t k = g.offsets[ix * g.ny + iy]; k < g.offsets[ix * g.ny + iy + 1]; ++k) {
          const Piece &piece = g.members[k];
          const std::uint32_t item = piece.item;
          if (seen(tested, item) ||
              !enters(piece.low, piece.high, ray, inverse, best ? best->t : t_max))
            continue;
          tested[next_slot] = item;
          next_slot = (next_slot + 1) % kMailbox;
#ifdef OBDEECT_BENCH_COUNTERS
          ++counters().exact_tests;
#endif
          const Obscurer &o = items_[item];
          if (const auto t = o.intersect(ray, planes_, kMinDistanceM))
            if (*t < t_max && (!best || *t < best->t || (*t == best->t && o.id < best->id)))
              best = Hit{*t, o.id};
        }
    return best;
  }

  // Is `item` in the mailbox? A plain loop: the ring holds only a few entries.
  [[nodiscard]] static bool seen(const std::array<std::uint32_t, kMailbox> &ring,
                                 std::uint32_t item) {
    bool found = false;
    for (const std::uint32_t entry : ring)
      found |= entry == item;
    return found;
  }

  [[nodiscard]] static double component(const Vec3 &v, std::size_t axis) {
    return axis == 0 ? v.x : axis == 1 ? v.y : v.z;
  }

  [[nodiscard]] static double area(const Vec3 &low, const Vec3 &high) {
    const Vec3 e = high - low;
    return e.x * e.y + e.y * e.z + e.z * e.x;
  }

  // Bounding-box hierarchy of pieces [begin, end): leaves of up to kLeafPieces,
  // binned surface-area-heuristic splits along the axis of widest centroid
  // spread, the median when all centroids share a bin or below kSahDepth.
  // Ties are broken by item index, so the build is deterministic.
  void build(Tree &tree, std::uint32_t begin, std::uint32_t end, std::size_t depth) {
    if (depth > kMaxDepth)
      throw std::length_error("obscurer hierarchy deeper than the traversal stack");
    const auto index = static_cast<std::uint32_t>(tree.nodes.size());
    tree.nodes.push_back({});
    Vec3 lo = tree.pieces[begin].low, hi = tree.pieces[begin].high;
    Vec3 clo = (lo + hi) * 0.5, chi = clo;
    for (std::uint32_t k = begin; k < end; ++k) {
      const Piece &piece = tree.pieces[k];
      const Vec3 c = (piece.low + piece.high) * 0.5;
      lo = {std::min(lo.x, piece.low.x), std::min(lo.y, piece.low.y), std::min(lo.z, piece.low.z)};
      hi = {std::max(hi.x, piece.high.x), std::max(hi.y, piece.high.y),
            std::max(hi.z, piece.high.z)};
      clo = {std::min(clo.x, c.x), std::min(clo.y, c.y), std::min(clo.z, c.z)};
      chi = {std::max(chi.x, c.x), std::max(chi.y, c.y), std::max(chi.z, c.z)};
    }
    tree.nodes[index].low = lo;
    tree.nodes[index].high = hi;
    if (end - begin <= kLeafPieces) {
      tree.nodes[index].first = begin;
      tree.nodes[index].count = end - begin;
      return;
    }
    const Vec3 spread = chi - clo;
    const std::size_t axis = spread.x >= spread.y && spread.x >= spread.z ? 0
                             : spread.y >= spread.z                       ? 1
                                                                          : 2;
    tree.nodes[index].axis = static_cast<std::uint8_t>(axis);
    const auto key = [axis](const Piece &piece) {
      return 0.5 * (component(piece.low, axis) + component(piece.high, axis));
    };
    const double c0 = component(clo, axis), width = component(chi, axis) - c0;
    std::uint32_t middle = begin + (end - begin) / 2;
    bool split = false;
    if (depth < kSahDepth && width > 0.0) {
      const auto bin_of = [&](const Piece &piece) {
        return std::min(kBins - 1, static_cast<std::size_t>((key(piece) - c0) / width * kBins));
      };
      std::array<std::uint32_t, kBins> count{};
      std::array<Vec3, kBins> bin_low, bin_high;
      bin_low.fill({1e300, 1e300, 1e300});
      bin_high.fill({-1e300, -1e300, -1e300});
      for (std::uint32_t k = begin; k < end; ++k) {
        const Piece &piece = tree.pieces[k];
        const std::size_t b = bin_of(piece);
        ++count[b];
        bin_low[b] = {std::min(bin_low[b].x, piece.low.x), std::min(bin_low[b].y, piece.low.y),
                      std::min(bin_low[b].z, piece.low.z)};
        bin_high[b] = {std::max(bin_high[b].x, piece.high.x), std::max(bin_high[b].y, piece.high.y),
                       std::max(bin_high[b].z, piece.high.z)};
      }
      // Cost of splitting after bin b: left bins 0..b, right bins b+1...
      std::array<double, kBins> left_cost{};
      Vec3 l = {1e300, 1e300, 1e300}, h = {-1e300, -1e300, -1e300};
      std::uint32_t n = 0;
      for (std::size_t b = 0; b + 1 < kBins; ++b) {
        if (count[b]) {
          l = {std::min(l.x, bin_low[b].x), std::min(l.y, bin_low[b].y),
               std::min(l.z, bin_low[b].z)};
          h = {std::max(h.x, bin_high[b].x), std::max(h.y, bin_high[b].y),
               std::max(h.z, bin_high[b].z)};
          n += count[b];
        }
        left_cost[b] = n ? area(l, h) * n : 0.0;
      }
      double best_cost = std::numeric_limits<double>::infinity();
      std::size_t best_bin = 0;
      std::uint32_t left_count = 0;
      l = {1e300, 1e300, 1e300};
      h = {-1e300, -1e300, -1e300};
      n = 0;
      std::array<std::uint32_t, kBins> left_n{};
      for (std::size_t b = 0, m = 0; b + 1 < kBins; ++b) {
        m += count[b];
        left_n[b] = static_cast<std::uint32_t>(m);
      }
      for (std::size_t b = kBins - 1; b-- > 0;) {
        const std::size_t r = b + 1;
        if (count[r]) {
          l = {std::min(l.x, bin_low[r].x), std::min(l.y, bin_low[r].y),
               std::min(l.z, bin_low[r].z)};
          h = {std::max(h.x, bin_high[r].x), std::max(h.y, bin_high[r].y),
               std::max(h.z, bin_high[r].z)};
          n += count[r];
        }
        if (left_n[b] == 0 || n == 0)
          continue;
        const double cost = left_cost[b] + area(l, h) * n;
        if (cost < best_cost) {
          best_cost = cost;
          best_bin = b;
          left_count = left_n[b];
        }
      }
      if (left_count > 0 && left_count < end - begin) {
        std::stable_partition(tree.pieces.begin() + begin, tree.pieces.begin() + end,
                              [&](const Piece &piece) { return bin_of(piece) <= best_bin; });
        middle = begin + left_count;
        split = true;
      }
    }
    if (!split)
      std::nth_element(tree.pieces.begin() + begin, tree.pieces.begin() + middle,
                       tree.pieces.begin() + end, [&](const Piece &a, const Piece &b) {
                         return key(a) < key(b) || (key(a) == key(b) && a.item < b.item);
                       });
    build(tree, begin, middle, depth + 1);
    tree.nodes[index].first = static_cast<std::uint32_t>(tree.nodes.size());
    build(tree, middle, end, depth + 1);
  }

  // Does the ray enter the box before t_max? (Slab test. A zero direction
  // component gives infinite slab distances; a ray lying exactly in a face
  // plane may be rejected, as it can at most touch the capsule inside.)
  [[nodiscard]] static bool enters(const Vec3 &low, const Vec3 &high, const Ray &ray,
                                   const Vec3 &inverse, double t_max) {
    const double tx0 = (low.x - ray.origin.x) * inverse.x,
                 tx1 = (high.x - ray.origin.x) * inverse.x;
    const double ty0 = (low.y - ray.origin.y) * inverse.y,
                 ty1 = (high.y - ray.origin.y) * inverse.y;
    const double tz0 = (low.z - ray.origin.z) * inverse.z,
                 tz1 = (high.z - ray.origin.z) * inverse.z;
    double enter = kMinDistanceM, leave = t_max;
    enter = std::max(enter, std::min(tx0, tx1));
    leave = std::min(leave, std::max(tx0, tx1));
    enter = std::max(enter, std::min(ty0, ty1));
    leave = std::min(leave, std::max(ty0, ty1));
    enter = std::max(enter, std::min(tz0, tz1));
    leave = std::min(leave, std::max(tz0, tz1));
    return enter <= leave;
  }

  std::vector<Obscurer> items_;
  std::vector<Plane> planes_;
  std::array<Tree, 3> trees_;
  Grid grid_;
};

} // namespace obdeect
