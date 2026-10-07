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
   measurement->mirror=result.mirror;
   measurement->primary[0]=r1_pos_x; measurement->primary[1]=r1_pos_y;
   measurement->primary[2]=r1_pos_z; measurement->secondary[0]=r2_pos_x;
   measurement->secondary[1]=r2_pos_y; measurement->secondary[2]=r2_pos_z;
   measurement->primary_cosine=dc_prm; measurement->secondary_cosine=dc_sec;
}
"""


def instrument(root: Path, destination: Path) -> tuple[Path, Path]:
    """Instrument exact setup/loss boundaries, failing when source anchors changed."""
    main = root / "common/sim_telarray.c"
    imaging = root / "common/sim_imaging.c"
    text = main.read_text()
    if text.count(_ANCHOR) != 1 or '#include "sim_optical_backend.h"' not in text:
        raise ValueError("unsupported sim_telarray setup/optical backend boundary")
    declaration = "extern int obdeect_reference_replay(void *, unsigned);\n"
    text = text.replace(
        '#include "sim_optical_backend.h"', '#include "sim_optical_backend.h"\n' + declaration
    )
    text = text.replace(
        _ANCHOR, "   exit(obdeect_reference_replay(array.optics,array.max_tel));\n" + _ANCHOR
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
