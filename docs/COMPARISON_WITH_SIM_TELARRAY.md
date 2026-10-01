# How to compare obdeect with sim_telarray

For the current implementation status, see [STATUS.md](STATUS.md).

## What this document is for

This is the practical recipe for comparing obdeect and sim_telarray.

In one sentence:

> Give both programs the same photons, through equivalent optical models, and
> compare the same measurable results.

It is not a list of telescope models or a claim that the two programs already
agree. It describes the information that must be fixed before such a comparison
is meaningful.

## 1. Start with the same photons

Every input photon needs an unambiguous identity:

- run, event, array, telescope, bunch, and photon IDs;
- starting position and direction in the telescope coordinate system;
- wavelength, emission time, and weight.

The source description must also say what was simulated: for example, a star,
an air shower, or a calibration laser; its position and coordinate frame; its
energy or wavelength distribution; its time distribution; the random seed; and
how the photon weights were assigned.

The two programs must receive the same input sample and use the same coordinate
and unit conventions. Input files may use nanometres and nanoseconds at the CSV
boundary, but the comparison adapter must convert them to metres and seconds
before comparing results. A wavelength of zero is only an unresolved-input
marker; it must be replaced by a real spectrum before wavelength-dependent
transport is compared.

## 2. Describe the optical model used

The optical model must identify:

- the source model and the model-data files used;
- hashes of those input files;
- telescope and camera coordinate frames;
- mirror, obscurer, detector, and pixel IDs;
- mirror response or material data;
- detector geometry; and
- any required geometry or response that is still missing.

If required geometry or response data is missing, the optical model is not ready
for a production comparison. The validation must stop and report what is
missing; it must not silently replace it with an approximation.

## 3. Report the same result for every photon

For each input photon, both programs should report:

| Quantity | Meaning |
| --- | --- |
| Photon ID | Which input photon this row belongs to. |
| Status | What happened: detected, missed a surface, or was lost at a component. |
| Last position and direction | Where the ray stopped and which way it was travelling. |
| Surface or pixel ID | The last optical component involved. |
| Optical path and arrival time | How far the ray travelled and when it arrived. |
| Wavelength and weight | The photon properties at that point. |
| Interaction history | The ordered list of mirror, obscurer, and detector interactions. |

At minimum, a detected photon needs its focal-plane position, optical path
length, arrival time, wavelength, weight, and incidence angles. A lost photon
needs a useful loss category, such as missed primary mirror, blocked by an
obscurer, missed secondary, or failed to reach the detector.

Detector arrival, concentrator loss, and sensor conversion are different
outcomes and must not be folded into one generic “lost” category.

## 4. How we will compare the outputs

The comparison should proceed in this order:

1. Check that both programs saw the same photon IDs.
2. Compare the terminal status of every photon.
3. For photons detected by both programs, compare focal position, path length,
   arrival time, wavelength, and incidence angles.
4. Record the first photon and first quantity that disagree.
5. Report the declared tolerances and the largest residuals.

Do not compare histograms first. A matching histogram can hide different
photon-by-photon behaviour.

## Current implementation status

`CsvPhotonReader`, `MemoryPhotonReader`, and the optional
`EventioPhotonReader` use the shared `OpticalPhoton` batch type. EventIO tests
cover full, compact, and 3-D photon files. `RigidFrame` validates right-handed
coordinate transforms and their round trips.

The native tracer still writes a short, fixed path record. It does not yet
provide the complete interaction history, group arrival time, or pixel-level
result required above. Until those fields are available, the project can run
diagnostic comparisons, but it must not claim a complete obdeect/sim_telarray
reference comparison.
