"""Exact (50-digit) reference values of the aspheric sag and slope / r polynomials.

usage: python make_polynomial_reference.py [--points N] model.json [model.json ...] > ref.txt

For every aspheric surface found in the optical-model files, evaluates at N
values of q = (r / scale)^2 (default 2000, r uniform up to the outer radius):

    sag(q)        = vertex_z + sum_i c_i q^i
    slope/r (q)   = (1 / scale^2) sum_{i>=1} 2 i c_i q^(i-1)

with mpmath at 50 digits, from the double-precision inputs taken as exact.
Doubles are written as C99 hexadecimal floats (exact); each reference as hi + lo
doubles (hexadecimal) and then as a 50-digit decimal.
Columns: q sag_hi sag_lo slope_hi slope_lo sag slope. The coefficients of real
telescope models must not be committed: write the output outside the repository.
"""

import argparse
import json
import sys

import mpmath


def surfaces(node, path=""):
    """Yield (path, asphere dict) for every object that has coefficients_m."""
    if isinstance(node, dict):
        if "coefficients_m" in node:
            yield path, node
        for key, value in node.items():
            yield from surfaces(value, f"{path}/{key}")
    elif isinstance(node, list):
        for value in node[:1]:  # Segment lists repeat the same asphere description.
            yield from surfaces(value, f"{path}[]")


def split(value):
    """The value as hi + lo doubles (hexadecimal): about 106 bits."""
    hi = float(value)
    return f"{hi.hex()} {float(value - mpmath.mpf(hi)).hex()}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("models", nargs="+")
    parser.add_argument("--points", type=int, default=2000)
    args = parser.parse_args()

    mpmath.mp.dps = 50
    print(f"# mpmath {mpmath.__version__}, python {sys.version.split()[0]}, dps 50")
    for name in args.models:
        print(f"# input {name}")
    for name in args.models:
        with open(name, encoding="utf-8") as stream:
            model = json.load(stream)
        for path, a in surfaces(model.get("optics", model)):
            c = [float(x) for x in a["coefficients_m"]]
            vertex, outer, scale = (
                float(a["vertex_z_m"]),
                float(a["outer_radius_m"]),
                float(a["radial_scale_m"]),
            )
            print(f"surface {name.rsplit('/', 1)[-1]}{path} {len(c)} {args.points}")
            print("coefficients " + " ".join(x.hex() for x in c))
            print(f"vertex {vertex.hex()} scale {scale.hex()}")
            for i in range(1, args.points + 1):
                r = outer * i / args.points
                q = (r / scale) * (r / scale)
                mq = mpmath.mpf(q)
                sag = mpmath.mpf(vertex) + mpmath.polyval([mpmath.mpf(x) for x in reversed(c)], mq)
                # d/dq of sum c_i q^i, times 1 / scale^2.
                deriv = [2 * k * mpmath.mpf(c[k]) for k in range(len(c) - 1, 0, -1)]
                slope = (
                    mpmath.polyval(deriv, mq) / (mpmath.mpf(scale) ** 2) if deriv else mpmath.mpf(0)
                )
                columns = [q.hex(), split(sag), split(slope), mpmath.nstr(sag, 50)]
                columns.append(mpmath.nstr(slope, 50))
                print(" ".join(columns))


if __name__ == "__main__":
    main()
