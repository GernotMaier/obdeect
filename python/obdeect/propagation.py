"""Resolve the uniform ambient refractive index from an explicit site profile."""

import bisect
import math

from obdeect.simtel_tables import interval_coefficients


def ambient_group_index(contents: str, observation_level_m: float) -> float:
    """Use the reference natural spline and 40,000-point logarithmic lookup grid."""
    lines = [
        line.split()
        for line in contents.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines or "altitude" not in lines[0] or "refractive_index" not in lines[0]:
        raise ValueError("atmospheric profile requires altitude and refractive_index columns")
    altitude_column = lines[0].index("altitude")
    index_column = lines[0].index("refractive_index")
    samples = [(float(row[altitude_column]) * 1000, float(row[index_column])) for row in lines[1:]]
    if len(samples) < 2 or any(
        not math.isfinite(h) or not math.isfinite(n) or n <= 0 for h, n in samples
    ):
        raise ValueError("atmospheric profile requires finite positive refractivity")
    heights = [height for height, _ in samples]
    if any(a >= b for a, b in zip(heights, heights[1:])) or not math.isfinite(observation_level_m):
        raise ValueError("invalid atmospheric altitude grid or observation level")
    if observation_level_m <= heights[0]:
        return 1 + samples[0][1]
    if observation_level_m >= heights[-1]:
        return 1
    if observation_level_m > heights[-2]:
        raise ValueError("observation level lies in the atmosphere's extrapolated upper layer")
    coefficients = interval_coefficients(heights, [math.log(n) for _, n in samples], 3)

    def spline(height):
        lower = min(len(heights) - 2, max(0, bisect.bisect_right(heights, height) - 1))
        fraction = (height - heights[lower]) / (heights[lower + 1] - heights[lower])
        return sum(
            coefficient * fraction**power for power, coefficient in enumerate(coefficients[lower])
        )

    step = (heights[-1] - heights[0]) / 39999
    coordinate = (observation_level_m - heights[0]) / step
    lower = int(coordinate)
    fraction = coordinate - lower
    first = spline(heights[0] + lower * step)
    second = spline(heights[0] + (lower + 1) * step)
    return 1 + math.exp(first + fraction * (second - first))
