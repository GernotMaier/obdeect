#pragma once

#include "obdeect/surfaces.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
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
          add_pieces(tree, i);
      if (!tree.pieces.empty())
        build(tree, 0, static_cast<std::uint32_t>(tree.pieces.size()));
    }
  }

  [[nodiscard]] const std::vector<Obscurer> &items() const { return items_; }
  [[nodiscard]] const std::vector<Plane> &planes() const { return planes_; }

  // Nearest obscurer on `leg` with kMinDistanceM < t < t_max; the lowest id
  // wins exact ties, so the result does not depend on the hierarchy.
  [[nodiscard]] std::optional<Hit> nearest(const Ray &ray, Leg leg, double t_max) const {
    const auto &tree = trees_[leg == Leg::incoming ? 0 : leg == Leg::after_primary ? 1 : 2];
    std::optional<Hit> best;
    if (tree.nodes.empty())
      return best;
    const Vec3 inverse{1.0 / ray.direction.x, 1.0 / ray.direction.y, 1.0 / ray.direction.z};
    std::array<std::uint32_t, 64> stack{};
    std::size_t top = 0;
    stack[top++] = 0;
    std::uint32_t last = kNone; // Pieces of one obscurer often share a leaf.
    while (top > 0) {
      const Node &node = tree.nodes[stack[--top]];
      if (!enters(node.low, node.high, ray, inverse, best ? best->t : t_max))
        continue;
      if (node.count == 0) {
        stack[top++] = node.first; // Right child; the left one follows its parent.
        stack[top++] = static_cast<std::uint32_t>(&node - tree.nodes.data()) + 1;
        continue;
      }
      for (std::uint32_t k = node.first; k < node.first + node.count; ++k) {
        const Piece &piece = tree.pieces[k];
        if (piece.item == last ||
            !enters(piece.low, piece.high, ray, inverse, best ? best->t : t_max))
          continue;
        last = piece.item;
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
  struct Piece {
    Vec3 low, high;
    std::uint32_t item{};
  };
  // Leaf: `count` pieces from pieces[first]. Inner node (count 0): left child
  // at the next index, right child at `first`.
  struct Node {
    Vec3 low, high;
    std::uint32_t first{}, count{};
  };
  struct Tree {
    std::vector<Node> nodes;
    std::vector<Piece> pieces;
  };

  // Pieces of the capsule no longer than a few capsule radii.
  void add_pieces(Tree &tree, std::uint32_t item) {
    const Obscurer &o = items_[item];
    const double r = o.bound_radius, length = norm(o.bound_b - o.bound_a);
    const auto count = static_cast<std::uint32_t>(
        std::clamp(std::ceil(length / (4.0 * std::max(r, 1e-6))), 1.0, 64.0));
    for (std::uint32_t i = 0; i < count; ++i) {
      const Vec3 p = o.bound_a + (o.bound_b - o.bound_a) * (static_cast<double>(i) / count);
      const Vec3 q = o.bound_a + (o.bound_b - o.bound_a) * (static_cast<double>(i + 1) / count);
      tree.pieces.push_back(
          {Vec3{std::min(p.x, q.x) - r, std::min(p.y, q.y) - r, std::min(p.z, q.z) - r},
           Vec3{std::max(p.x, q.x) + r, std::max(p.y, q.y) + r, std::max(p.z, q.z) + r}, item});
    }
  }

  void build(Tree &tree, std::uint32_t begin, std::uint32_t end) {
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
    if (end - begin <= 4) {
      tree.nodes[index].first = begin;
      tree.nodes[index].count = end - begin;
      return;
    }
    // Median split along the widest spread of the piece centres.
    const Vec3 spread = chi - clo;
    const int axis = spread.x >= spread.y && spread.x >= spread.z ? 0
                     : spread.y >= spread.z                       ? 1
                                                                  : 2;
    const auto key = [axis](const Piece &piece) {
      const Vec3 c = (piece.low + piece.high) * 0.5;
      return axis == 0 ? c.x : axis == 1 ? c.y : c.z;
    };
    const std::uint32_t middle = begin + (end - begin) / 2;
    std::nth_element(tree.pieces.begin() + begin, tree.pieces.begin() + middle,
                     tree.pieces.begin() + end, [&](const Piece &a, const Piece &b) {
                       return key(a) < key(b) || (key(a) == key(b) && a.item < b.item);
                     });
    build(tree, begin, middle);
    tree.nodes[index].first = static_cast<std::uint32_t>(tree.nodes.size());
    build(tree, middle, end);
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
};

} // namespace obdeect
