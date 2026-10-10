"""Build a shared-photon diagnostic reference without modifying the simulator checkout.

Compiler/linker argv come from an explicit JSON recipe. Placeholders are literal
argv entries or fragments: {source}, {object}, {driver}, {executable}, {include}.
Record all compiler flags and link objects from the matching simulator build.
The result records genuine hits/losses; it is not a production-normalized table.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

_ANCHOR = "   if ( (options.save_pe_amp & 0x01) )"
_BRIDGE = r"""
/* Diagnostic boundary instrumentation; combined simulator: GPL-3.0-or-later. */
struct ReferenceMeasurement {
   int status, mirror, loss_line;
   double position[3], direction[3], time_ns, relative_efficiency;
   double primary[3], secondary[3], primary_cosine, secondary_cosine;
   double propagation_group_index, primary_envelope, secondary_envelope;
   double upstream_optical_efficiency;
   int pixel_status, pixel_id;
   double pixel_time_ns, camera_absolute_efficiency, pixel_x_cm, pixel_y_cm;
};
extern int obdeect_reference_loss_line;
extern double r1_pos_x,r1_pos_y,r1_pos_z,r2_pos_x,r2_pos_y,r2_pos_z,dc_prm,dc_sec;
void obdeect_reference_measure(void *storage, unsigned index,
   const double *position, const double *direction, double wavelength,
   double distance, double weight, struct ReferenceMeasurement *measurement);
void obdeect_reference_measure(void *storage, unsigned index,
   const double *position, const double *direction, double wavelength,
   double distance, double weight, struct ReferenceMeasurement *measurement)
{
   struct telescope_optics *optics = ((struct telescope_optics *)storage)+index;
   struct sim_optical_input input;
   struct sim_optical_result result;
   memset(&input,0,sizeof(input));
   memcpy(input.position,position,sizeof(input.position));
   memcpy(input.direction,direction,sizeof(input.direction));
   transform_off(input.position,input.direction,&optics->tel_trans,-1);
   input.wavelength=wavelength; input.distance=distance; input.relative_efficiency=weight;
   r1_pos_x=r1_pos_y=r1_pos_z=r2_pos_x=r2_pos_y=r2_pos_z=NAN;
   dc_prm=dc_sec=NAN; obdeect_reference_loss_line=0;
   measurement->status=sim_optical_trace(optics,&input,&result);
   measurement->loss_line=obdeect_reference_loss_line;
   memcpy(measurement->position,result.position,sizeof(result.position));
   memcpy(measurement->direction,result.direction,sizeof(result.direction));
   measurement->time_ns=result.travel_time;
   measurement->relative_efficiency=result.relative_efficiency;
   measurement->propagation_group_index=29.9792458/airlightspeed;
   measurement->primary_envelope=optics->with_mirror_ref_2d ?
      rpolate_1d(optics->mirror_ref_2d,wavelength,-1) : 1.;
   measurement->secondary_envelope=optics->with_mirror2_ref_2d ?
      rpolate_1d(optics->mirror2_ref_2d,wavelength,-1) : 1.;
   measurement->upstream_optical_efficiency=(wavelength>=0. && wavelength<MAX_LAMBDA) ?
      obdeect_reference_electronics[index].optics_efficiency[(int)wavelength] : 0.;
   measurement->mirror=result.mirror;
   measurement->primary[0]=r1_pos_x; measurement->primary[1]=r1_pos_y;
   measurement->primary[2]=r1_pos_z; measurement->secondary[0]=r2_pos_x;
   measurement->secondary[1]=r2_pos_y; measurement->secondary[2]=r2_pos_z;
   measurement->primary_cosine=dc_prm; measurement->secondary_cosine=dc_sec;
   measurement->pixel_status=-1; measurement->pixel_id=-1;
   measurement->pixel_time_ns=NAN; measurement->camera_absolute_efficiency=0.;
   measurement->pixel_x_cm=measurement->pixel_y_cm=NAN;
   if (measurement->status==0 && wavelength>0. && wavelength<MAX_LAMBDA) {
      struct pm_camera *camera=obdeect_reference_camera+index;
      double p[3], d[3], x=result.position[0], y=result.position[1];
      double sx=result.direction[0]/result.direction[2], sy=result.direction[1]/result.direction[2];
      double time=result.travel_time, efficiency=result.relative_efficiency;
      int iwl=(int)wavelength, pixel;
      memcpy(p,result.position,sizeof(p)); memcpy(d,result.direction,sizeof(d));
      pixel=camera_hit(camera,optics,p,d,&x,&y,&sx,&sy,&time);
      measurement->pixel_id=pixel;
      if(pixel>=0 && pixel<camera->pixels) {
         if(optics->camera_degraded_map!=NULL) {
            double degradation=rpolate_2d(optics->camera_degraded_map,x,y,1);
            if(degradation<1.) efficiency*=degradation;
         }
         if(camera->with_filter_2d) {
            double maximum=rpolate_1d(camera->filter_trans_2d,wavelength,-1);
            double filter=rpolate_2d(camera->filter_trans_2d,wavelength,atan(hypot(sx,sy)),1);
            efficiency*=maximum!=0. ? filter/maximum : 0.;
         }
         measurement->pixel_status=cathode_hit(camera,pixel,x,y,sx,sy,iwl,&efficiency);
         measurement->pixel_time_ns=time;
         measurement->pixel_x_cm=x; measurement->pixel_y_cm=y;
         if(measurement->pixel_status>0)
            measurement->camera_absolute_efficiency=efficiency*
               measurement->upstream_optical_efficiency*camera->filter_trans[iwl]*
               optics->camera_transmission*optics->camera_degraded_efficiency;
      }
   }
}
"""


def instrument(root: Path, destination: Path) -> tuple[Path, Path]:
    """Instrument exact setup/loss boundaries, failing when source anchors changed."""
    main = root / "common/sim_telarray.c"
    imaging = root / "common/sim_imaging.c"
    text = main.read_text()
    if text.count(_ANCHOR) != 1 or '#include "sim_optical_backend.h"' not in text:
        raise ValueError("unsupported sim_telarray setup/optical backend boundary")
    declaration = (
        "extern int obdeect_reference_replay(void *, unsigned);\n"
        "static struct camera_electronics *obdeect_reference_electronics;\n"
        "static struct pm_camera *obdeect_reference_camera;\n"
    )
    text = text.replace(
        '#include "sim_optical_backend.h"', '#include "sim_optical_backend.h"\n' + declaration
    )
    text = text.replace(
        _ANCHOR,
        "   obdeect_reference_electronics=array.electronics;\n"
        "   obdeect_reference_camera=array.camera;\n"
        "   exit(obdeect_reference_replay(array.optics,array.max_tel));\n" + _ANCHOR,
    )
    main_output = destination / "simtel_replay_main.c"
    main_output.write_text(text + _BRIDGE)
    text = imaging.read_text()
    declarations = 'int obdeect_reference_loss_line=0;\n#line 1 "common/sim_imaging.c"\n'
    names = ("segmented", "paraboloid", "fresnel", "secondary")
    spans = []
    for name in names:
        symbol = f"trace_photon_{'with_' if name == 'secondary' else 'in_'}{name}"
        match = re.search(r"^int " + symbol + r"\s*\(", text, re.MULTILINE)
        if not match:
            raise ValueError(f"unsupported reference tracer: {symbol}")
        end = text.find("\n}\n", match.end())
        if end < 0:
            raise ValueError(f"unterminated reference tracer: {symbol}")
        spans.append((match.start(), end))
    for start, end in reversed(spans):
        block = text[start:end]
        if not re.search(r"return\s+-1\s*;", block):
            raise ValueError("reference tracer has no recognized loss boundary")
        block = re.sub(
            r"return\s+-1\s*;", "{ obdeect_reference_loss_line=__LINE__; return -1; }", block
        )
        text = text[:start] + block + text[end:]
    imaging_output = destination / "simtel_replay_imaging.c"
    imaging_output.write_text(declarations + text)
    return main_output, imaging_output


def build(root: Path, destination: Path, recipe: dict) -> dict:
    """Use an explicit, reproducible toolchain; retain hashes and build logs."""
    expected = {"compile_c", "compile_cxx", "link", "inputs"}
    if set(recipe) != expected:
        raise ValueError(f"build recipe requires exactly {sorted(expected)}")
    for key in ("compile_c", "compile_cxx", "link", "inputs"):
        if (
            not isinstance(recipe[key], list)
            or not recipe[key]
            or not all(isinstance(item, str) for item in recipe[key])
        ):
            raise ValueError(f"{key} must be a nonempty list of strings")
    root = root.resolve()
    destination = destination.resolve()
    if destination == root or root in destination.parents:
        raise ValueError("build outside the simulator source checkout")
    destination.mkdir(parents=True, exist_ok=True)
    sources = instrument(root, destination)
    repo = Path(__file__).resolve().parents[1]
    executable = destination / "simtel-replay"
    values = {
        "include": str(repo / "cpp/include"),
        "driver": str(destination / "driver.o"),
        "executable": str(executable),
        "main_object": str(destination / "main.o"),
        "imaging_object": str(destination / "imaging.o"),
    }
    inputs = [
        root / "common/sim_telarray.c",
        root / "common/sim_imaging.c",
        repo / "tools/simtel_replay.cpp",
        *map(Path, recipe["inputs"]),
    ]
    before = {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    commands = []
    for source, name in zip(sources, ("main", "imaging"), strict=True):
        arguments = {**values, "source": str(source), "object": str(destination / f"{name}.o")}
        commands.append([arg.format_map(arguments) for arg in recipe["compile_c"]])
    arguments = {
        **values,
        "source": str(repo / "tools/simtel_replay.cpp"),
        "object": values["driver"],
    }
    commands.append([arg.format_map(arguments) for arg in recipe["compile_cxx"]])
    commands.append([arg.format_map(values) for arg in recipe["link"]])
    with (destination / "build.log").open("w") as log:
        for command in commands:
            log.write(json.dumps(command) + "\n")
            log.flush()
            subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT, check=True)
    if any(
        hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest
        for path, digest in before.items()
    ):
        raise ValueError("reference build inputs changed during compilation")
    record = {
        "schema_version": 1,
        "production_normalized": False,
        "input_sha256": before,
        "commands": commands,
        "cwd": str(root),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
    }
    (destination / "build-record.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--simtel-root", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--build-recipe", required=True, type=Path)
    args = parser.parse_args()
    build(args.simtel_root, args.output_directory, json.loads(args.build_recipe.read_text()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
