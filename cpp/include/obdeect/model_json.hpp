#pragma once

#include "obdeect/json.hpp"
#include "obdeect/model.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace obdeect {

inline constexpr std::string_view kModelFormat = "obdeect.optical-model.v2";

namespace detail {

// A JSON value together with its path, for precise error messages.
class Node {
public:
  Node(const json::Value &value, std::string path) : value_(&value), path_(std::move(path)) {}

  [[noreturn]] void fail(const std::string &what) const { throw ModelError(path_ + ": " + what); }

  [[nodiscard]] bool is_null() const { return value_->kind() == json::Value::Kind::null; }

  [[nodiscard]] const json::Value::Object &
  object(std::initializer_list<std::string_view> keys) const {
    if (value_->kind() != json::Value::Kind::object)
      fail("expected an object");
    const auto &members = std::get<json::Value::Object>(value_->data);
    for (const auto &[name, value] : members)
      if (std::find(keys.begin(), keys.end(), name) == keys.end())
        fail("unknown key '" + name + "'");
    return members;
  }
  [[nodiscard]] Node operator[](std::string_view key) const {
    const json::Value *member = value_->find(key);
    if (member == nullptr)
      fail("missing key '" + std::string(key) + "'");
    return {*member, path_ + "." + std::string(key)};
  }
  [[nodiscard]] bool has(std::string_view key) const {
    const json::Value *member = value_->find(key);
    return member != nullptr && member->kind() != json::Value::Kind::null;
  }

  [[nodiscard]] std::vector<Node> array() const {
    if (value_->kind() != json::Value::Kind::array)
      fail("expected an array");
    std::vector<Node> items;
    const auto &values = std::get<json::Value::Array>(value_->data);
    for (std::size_t i = 0; i < values.size(); ++i)
      items.emplace_back(values[i], path_ + "[" + std::to_string(i) + "]");
    return items;
  }
  [[nodiscard]] double number() const {
    if (value_->kind() != json::Value::Kind::number)
      fail("expected a number");
    return std::get<double>(value_->data);
  }
  [[nodiscard]] double positive() const {
    const double value = number();
    if (!(value > 0.0))
      fail("expected a positive number");
    return value;
  }
  [[nodiscard]] double non_negative() const {
    const double value = number();
    if (!(value >= 0.0))
      fail("expected a non-negative number");
    return value;
  }
  [[nodiscard]] double fraction() const {
    const double value = number();
    if (!(value >= 0.0 && value <= 1.0))
      fail("expected a value in [0, 1]");
    return value;
  }
  [[nodiscard]] std::uint32_t id() const {
    const double value = number();
    if (!(value >= 0.0 && value < 4294967295.0) || std::floor(value) != value)
      fail("expected an unsigned 32-bit identifier");
    return static_cast<std::uint32_t>(value);
  }
  [[nodiscard]] bool boolean() const {
    if (value_->kind() != json::Value::Kind::boolean)
      fail("expected true or false");
    return std::get<bool>(value_->data);
  }
  [[nodiscard]] const std::string &string() const {
    if (value_->kind() != json::Value::Kind::string)
      fail("expected a string");
    return std::get<std::string>(value_->data);
  }
  [[nodiscard]] std::vector<double> numbers(std::size_t minimum = 1) const {
    std::vector<double> values;
    for (const auto &item : array())
      values.push_back(item.number());
    if (values.size() < minimum)
      fail("expected at least " + std::to_string(minimum) + " numbers");
    return values;
  }
  [[nodiscard]] std::vector<double> increasing() const {
    auto values = numbers();
    for (std::size_t i = 1; i < values.size(); ++i)
      if (!(values[i] > values[i - 1]))
        fail("expected strictly increasing values");
    return values;
  }
  [[nodiscard]] Vec3 vec3() const {
    const auto values = numbers(3);
    if (values.size() != 3)
      fail("expected three numbers");
    return {values[0], values[1], values[2]};
  }
  [[nodiscard]] Vec3 unit() const {
    const Vec3 value = vec3();
    if (std::abs(norm(value) - 1.0) > 1e-9)
      fail("expected a unit vector");
    return normalized(value);
  }
  // A 2-D table given as rows; checks the shape and the value range.
  [[nodiscard]] std::vector<double> matrix(std::size_t rows, std::size_t columns,
                                           bool fractions) const {
    const auto items = array();
    if (items.size() != rows)
      fail("expected " + std::to_string(rows) + " rows");
    std::vector<double> values;
    values.reserve(rows * columns);
    for (const auto &row : items) {
      const auto entries = row.array();
      if (entries.size() != columns)
        row.fail("expected " + std::to_string(columns) + " values");
      for (const auto &entry : entries)
        values.push_back(fractions ? entry.fraction() : entry.non_negative());
    }
    return values;
  }

private:
  const json::Value *value_;
  std::string path_;
};

inline Aperture aperture(const Node &node) {
  const auto &name = node.string();
  if (name == "circle")
    return Aperture::circle;
  if (name == "square")
    return Aperture::square;
  if (name == "hexagon_flat_x")
    return Aperture::hexagon_flat_x;
  if (name == "hexagon_flat_y")
    return Aperture::hexagon_flat_y;
  node.fail("unknown aperture shape '" + name + "'");
}

// Unit tangent orthogonal to `normal`.
inline Vec3 tangent(const Node &node, const Vec3 &normal) {
  const Vec3 t = node.unit();
  if (std::abs(dot(t, normal)) > 1e-9)
    node.fail("tangent is not perpendicular to the normal");
  return normalized(t - normal * dot(t, normal));
}

inline Response response(const Node &node) {
  (void)node.object({"wavelength_min_nm", "wavelength_step_nm", "incidence_deg", "values"});
  Response r;
  r.wavelength_min_nm = node["wavelength_min_nm"].positive();
  r.wavelength_step_nm = node["wavelength_step_nm"].positive();
  r.incidence_deg = node["incidence_deg"].increasing();
  const auto rows = node["values"].array();
  r.wavelength_count = rows.size();
  if (r.wavelength_count < 1)
    node["values"].fail("expected at least one row");
  r.values = node["values"].matrix(r.wavelength_count, r.incidence_deg.size(), true);
  return r;
}

inline Map2D map2d(const Node &node) {
  (void)node.object({"x_m", "y_m", "values", "outside"});
  Map2D m;
  m.x_m = node["x_m"].increasing();
  m.y_m = node["y_m"].increasing();
  m.values = node["values"].matrix(m.x_m.size(), m.y_m.size(), true);
  const auto outside = node["outside"].string();
  if (outside != "clamp" && outside != "zero")
    node["outside"].fail("expected 'clamp' or 'zero'");
  m.zero_outside = outside == "zero";
  return m;
}

inline Scatter scatter(const Node &node) {
  (void)node.object({"sigma1_rad", "fraction2", "sigma2_rad", "method"});
  Scatter s{node["sigma1_rad"].non_negative(), node["fraction2"].fraction(),
            node["sigma2_rad"].non_negative()};
  const auto method = node["method"].string();
  if (method != "outgoing_angles" && method != "surface_slopes")
    node["method"].fail("expected 'outgoing_angles' or 'surface_slopes'");
  s.surface_slopes = method == "surface_slopes";
  return s;
}

inline Asphere asphere(const Node &node, std::uint32_t id) {
  (void)node.object(
      {"vertex_z_m", "inner_radius_m", "outer_radius_m", "radial_scale_m", "coefficients_m"});
  Asphere a;
  a.id = id;
  a.vertex_z = node["vertex_z_m"].number();
  a.inner_radius = node["inner_radius_m"].non_negative();
  a.outer_radius = node["outer_radius_m"].positive();
  a.scale = node["radial_scale_m"].positive();
  a.coefficients = node["coefficients_m"].numbers();
  if (a.inner_radius >= a.outer_radius)
    node.fail("inner radius must be smaller than the outer radius");
  a.prepare();
  return a;
}

inline Segment segment(const Node &node) {
  constexpr double degree = 3.14159265358979323846 / 180.0;
  Segment s;
  const auto shape = node["shape"].string();
  if (shape == "annular_sector") {
    (void)node.object({"id", "shape", "inner_radius_m", "outer_radius_m", "start_deg", "span_deg",
                       "gap_m", "gap_at_start"});
    s.shape = Segment::Shape::annular_sector;
    s.inner_radius = node["inner_radius_m"].non_negative();
    s.outer_radius = node["outer_radius_m"].positive();
    s.start = node["start_deg"].number() * degree;
    s.span = node["span_deg"].positive() * degree;
    s.gap = node["gap_m"].non_negative();
    s.gap_at_start = node["gap_at_start"].boolean();
  } else if (shape == "polygon") {
    (void)node.object({"id", "shape", "vertices_xy_m"});
    s.shape = Segment::Shape::polygon;
    for (const auto &vertex : node["vertices_xy_m"].array()) {
      const auto xy = vertex.numbers(2);
      if (xy.size() != 2)
        vertex.fail("expected two numbers");
      s.vertices.push_back({xy[0], xy[1]});
    }
    if (s.vertices.size() < 3)
      node.fail("a polygon needs at least three vertices");
  } else {
    (void)node.object({"id", "shape", "centre_xy_m", "diameter_m", "rotation_deg"});
    if (shape == "hexagon")
      s.shape = Segment::Shape::hexagon;
    else if (shape == "square")
      s.shape = Segment::Shape::square;
    else if (shape == "circle")
      s.shape = Segment::Shape::circle;
    else
      node["shape"].fail("unknown segment shape '" + shape + "'");
    const auto xy = node["centre_xy_m"].numbers(2);
    s.centre_x = xy[0];
    s.centre_y = xy[1];
    s.half_width = 0.5 * node["diameter_m"].positive();
    s.rotation = node["rotation_deg"].number() * degree;
  }
  s.id = node["id"].id();
  return s;
}

inline Mirror mirror(const Node &node, std::vector<std::uint32_t> &ids) {
  (void)node.object(
      {"facets", "asphere", "asphere_id", "segments", "reflectivity", "scatter", "degradation"});
  Mirror m;
  if (node.has("facets") == node.has("asphere"))
    node.fail("expected exactly one of 'facets' and 'asphere'");
  if (node.has("facets")) {
    for (const auto &f : node["facets"].array()) {
      (void)f.object({"id", "centre_m", "normal", "tangent", "map_tangent", "shape", "diameter_m",
                      "radius_of_curvature_m"});
      Facet facet;
      facet.id = f["id"].id();
      facet.centre = f["centre_m"].vec3();
      facet.normal = f["normal"].unit();
      facet.u = tangent(f["tangent"], facet.normal);
      facet.v = cross(facet.normal, facet.u);
      facet.map_u = f.has("map_tangent") ? tangent(f["map_tangent"], facet.normal) : facet.u;
      facet.map_v = cross(facet.normal, facet.map_u);
      facet.shape = aperture(f["shape"]);
      facet.half_width = 0.5 * f["diameter_m"].positive();
      facet.radius = f["radius_of_curvature_m"].non_negative();
      facet.prepare();
      ids.push_back(facet.id);
      m.facets.push_back(facet);
    }
    if (m.facets.empty())
      node["facets"].fail("expected at least one facet");
    m.grid = FacetGrid(m.facets);
  } else {
    m.asphere = asphere(node["asphere"], node["asphere_id"].id());
    ids.push_back(m.asphere->id);
    if (node.has("segments"))
      for (const auto &s : node["segments"].array()) {
        m.segments.push_back(segment(s));
        m.segment_frames.push_back(m.segments.back().tangent_frame(*m.asphere));
        ids.push_back(m.segments.back().id);
      }
    m.segment_grid = SegmentGrid(m.segments, m.segment_frames, *m.asphere);
  }
  if (node.has("reflectivity"))
    m.reflectivity = response(node["reflectivity"]);
  if (node.has("scatter"))
    m.scatter = scatter(node["scatter"]);
  if (node.has("degradation"))
    m.degradation = map2d(node["degradation"]);
  return m;
}

inline std::uint8_t legs(const Node &node) {
  std::uint8_t mask = 0;
  for (const auto &item : node.array()) {
    const auto &name = item.string();
    if (name == "incoming")
      mask |= Leg::incoming;
    else if (name == "after_primary")
      mask |= Leg::after_primary;
    else if (name == "after_secondary")
      mask |= Leg::after_secondary;
    else
      item.fail("unknown leg '" + name + "'");
  }
  if (mask == 0)
    node.fail("expected at least one leg");
  return mask;
}

inline Obscurer obscurer(const Node &node, std::vector<Plane> &planes) {
  Obscurer o;
  const auto shape = node["shape"].string();
  if (shape == "cylinder") {
    (void)node.object({"id", "name", "group", "legs", "shape", "start_m", "end_m", "diameter_m"});
    o.kind = Obscurer::Kind::cylinder;
    o.a = node["start_m"].vec3();
    o.b = node["end_m"].vec3();
    o.radius_a = 0.5 * node["diameter_m"].positive();
    if (norm(o.b - o.a) == 0.0)
      node.fail("cylinder has zero length");
    o.bound_a = o.a;
    o.bound_b = o.b;
    o.bound_radius = o.radius_a;
    // The box of a closed cylinder: the end points widened by r sqrt(1 - n_i^2).
    const Vec3 n = normalized(o.b - o.a);
    const Vec3 e{o.radius_a * std::sqrt(std::max(0.0, 1.0 - n.x * n.x)),
                 o.radius_a * std::sqrt(std::max(0.0, 1.0 - n.y * n.y)),
                 o.radius_a * std::sqrt(std::max(0.0, 1.0 - n.z * n.z))};
    o.box_low = {std::min(o.a.x, o.b.x) - e.x, std::min(o.a.y, o.b.y) - e.y,
                 std::min(o.a.z, o.b.z) - e.z};
    o.box_high = {std::max(o.a.x, o.b.x) + e.x, std::max(o.a.y, o.b.y) + e.y,
                  std::max(o.a.z, o.b.z) + e.z};
    o.has_box = true;
  } else if (shape == "polyhedron") {
    (void)node.object({"id", "name", "group", "legs", "shape", "vertices_m", "planes"});
    o.kind = Obscurer::Kind::polyhedron;
    std::vector<Vec3> vertices;
    for (const auto &vertex : node["vertices_m"].array())
      vertices.push_back(vertex.vec3());
    if (vertices.size() < 4)
      node["vertices_m"].fail("expected at least four vertices");
    o.plane_begin = static_cast<std::uint32_t>(planes.size());
    for (const auto &plane : node["planes"].array()) {
      const auto values = plane.numbers(4);
      if (values.size() != 4)
        plane.fail("expected [nx, ny, nz, offset]");
      const Vec3 n{values[0], values[1], values[2]};
      if (std::abs(norm(n) - 1.0) > 1e-9)
        plane.fail("expected a unit normal");
      planes.push_back({normalized(n), values[3]});
    }
    o.plane_end = static_cast<std::uint32_t>(planes.size());
    if (o.plane_end - o.plane_begin < 4)
      node["planes"].fail("expected at least four planes");
    for (const auto &vertex : vertices)
      for (std::uint32_t i = o.plane_begin; i < o.plane_end; ++i)
        if (dot(planes[i].normal, vertex) > planes[i].offset + 1e-6)
          node.fail("vertex lies outside a face plane");
    bound_vertices(o, vertices);
  } else if (shape == "plate") {
    (void)node.object({"id", "name", "group", "legs", "shape", "centre_m", "normal", "tangent",
                       "aperture", "diameter_m"});
    o.kind = Obscurer::Kind::plate;
    o.centre = node["centre_m"].vec3();
    o.normal = node["normal"].unit();
    o.u = tangent(node["tangent"], o.normal);
    o.v = cross(o.normal, o.u);
    o.shape = aperture(node["aperture"]);
    o.half_width = 0.5 * node["diameter_m"].positive();
    o.bound_a = o.bound_b = o.centre;
    o.bound_radius = circumradius(o.shape, o.half_width);
  } else if (shape == "frustum") {
    (void)node.object({"id", "name", "group", "legs", "shape", "start_m", "end_m", "start_radius_m",
                       "end_radius_m", "thickness_m", "solid"});
    o.kind = Obscurer::Kind::frustum;
    o.a = node["start_m"].vec3();
    o.b = node["end_m"].vec3();
    o.radius_a = node["start_radius_m"].positive();
    o.radius_b = node["end_radius_m"].positive();
    o.thickness = node["thickness_m"].non_negative();
    o.solid = node["solid"].boolean();
    if (norm(o.b - o.a) == 0.0)
      node.fail("frustum has zero length");
    if (o.solid && o.thickness != 0.0)
      node.fail("a solid frustum has no wall thickness");
    o.bound_a = o.a;
    o.bound_b = o.b;
    o.bound_radius = std::max(o.radius_a, o.radius_b) + o.thickness;
  } else {
    node["shape"].fail("unknown obscurer shape '" + shape + "'");
  }
  o.id = node["id"].id();
  o.legs = legs(node["legs"]);
  (void)node["name"].string();
  (void)node["group"].string();
  return o;
}

} // namespace detail

// Parse and validate an optical model (format obdeect.optical-model.v2).
[[nodiscard]] inline OpticalModel load_model(std::string_view text) {
  using detail::Node;
  json::Value root_value;
  try {
    root_value = json::Parser::parse(text);
  } catch (const std::runtime_error &error) {
    throw ModelError(error.what());
  }
  const Node root{root_value, "model"};
  (void)root.object({"format", "sha256", "provenance", "notes", "sampling", "optics"});
  if (root["format"].string() != kModelFormat)
    root["format"].fail("expected '" + std::string(kModelFormat) + "'");
  const Node optics = root["optics"];
  (void)optics.object({"propagation_group_index", "telescope_transmission", "primary", "secondary",
                       "secondary_blocks_incoming", "focal_surface", "obscurers"});
  OpticalModel model;
  std::vector<std::uint32_t> ids;
  model.propagation_group_index = optics["propagation_group_index"].number();
  if (!(model.propagation_group_index >= 1.0))
    optics["propagation_group_index"].fail("expected a value >= 1");
  if (optics.has("telescope_transmission")) {
    const Node t = optics["telescope_transmission"];
    (void)t.object({"on_axis", "amplitude", "angular_scale_rad", "power", "outer_power"});
    model.telescope_transmission = TelescopeTransmission{
        t["on_axis"].fraction(), t["amplitude"].non_negative(), t["angular_scale_rad"].positive(),
        t["power"].positive(), t["outer_power"].positive()};
  }
  model.primary = detail::mirror(optics["primary"], ids);
  if (optics.has("secondary")) {
    model.secondary = detail::mirror(optics["secondary"], ids);
    if (!model.secondary->asphere)
      optics["secondary"].fail("a secondary mirror must be an asphere");
  }
  model.secondary_blocks_incoming = optics["secondary_blocks_incoming"].boolean();
  if (model.secondary_blocks_incoming && !model.secondary)
    optics["secondary_blocks_incoming"].fail("there is no secondary mirror");
  const Node focal = optics["focal_surface"];
  (void)focal.object({"id", "plane_z_m", "asphere"});
  model.focal_surface.id = focal["id"].id();
  ids.push_back(model.focal_surface.id);
  if (focal.has("plane_z_m") == focal.has("asphere"))
    focal.fail("expected exactly one of 'plane_z_m' and 'asphere'");
  if (focal.has("plane_z_m"))
    model.focal_surface.plane_z = focal["plane_z_m"].number();
  else
    model.focal_surface.asphere = detail::asphere(focal["asphere"], model.focal_surface.id);
  std::vector<Obscurer> obscurers;
  std::vector<Plane> planes;
  for (const auto &item : optics["obscurers"].array()) {
    obscurers.push_back(detail::obscurer(item, planes));
    ids.push_back(obscurers.back().id);
  }
  model.obscurers = Obscurers(std::move(obscurers), std::move(planes));
  std::unordered_set<std::uint32_t> unique;
  for (const auto id : ids)
    if (!unique.insert(id).second)
      throw ModelError("model.optics: surface id " + std::to_string(id) + " is used twice");
  return model;
}

} // namespace obdeect
