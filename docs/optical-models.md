# Optical models

## Frames and units

All lengths are metres, angles in degrees unless named `_rad`. The telescope
frame is sim_telarray's: +z along the optical axis towards the sky, the dish
vertex at z = -`mirror_offset`, and x, y as in the mirror lists and the camera
photon lists of sim_telarray (with the telescope pointing above the horizon, x
points down). `star_photons` samples a disk of radius 1.2 times the dish radius
around the optical axis, perpendicular to it and centred on the dish vertex, as
sim_telarray does for star light; `sampling_area_m2` is its area. The distance
of a point source is measured from the disk centre.

## File format

An optical model (`obdeect.optical-model.v2`) is one JSON document:

- `format`, `sha256` (of the document without it), `provenance` (parameter
  versions, file hashes, seed, zenith angle, structure flag), `notes`;
- `sampling`: the star-light disk;
- `optics`: `primary`, optional `secondary` (each either `facets` or an
  `asphere` with `segments`, a `reflectivity` table, `scatter`, optional
  `degradation` map), `focal_surface` (plane or asphere), `obscurers`
  (cylinder, polyhedron, plate, frustum; each acting on the legs `incoming`,
  `after_primary`, `after_secondary`), `telescope_transmission`,
  `propagation_group_index`.

The C++ loader accepts only known keys, checks every value and names the JSON
path of the first error. Surface and obscurer ids are unique.

For mirrors, obscurers and the telescope-structure table, see
[Optical elements](optical-elements.md).
