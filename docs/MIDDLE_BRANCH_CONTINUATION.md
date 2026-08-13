# Middle-branch continuation

The lower middle branches are locally closer to graphs of `b = b(sigma)` than
to graphs of `sigma = sigma(b)`. They therefore require the opposite corrector
geometry from the side branches.

## Corrector geometries

- `horizontal_side` (default): predict in `(sigma,b)`, hold predicted `b`
  fixed, sweep in `sigma`, and select a threshold crossing.
- `vertical_middle`: predict in `(sigma,b)`, hold predicted `sigma` fixed,
  sweep in `b`, and select a threshold crossing.

`bisect_interval.m` supports both modes. The middle mode is selected by

```matlab
modelpar.corrector_mode = 'vertical_middle';
modelpar.db_local = 0.01;
modelpar.local_halfwidth_steps_b = 6;
modelpar.expand_factors_b = [1 2 4];
```

The middle-right branches select a downward crossing in increasing `b` and
track the crossing near the previous accepted value. The middle-left branches
select an upward crossing and use the strongest valid crossing, matching the
historical middle-left metadata.

All feature evaluations in one corrector sweep use the same `ic_det`, so the
sampled response is not warm-started along the local `b` grid.

## Branch tasks

Both main scripts use the existing `TASK` convention:

| TASK | Branch | Corrector |
|---:|---|---|
| 0 | right-up | horizontal in `sigma` |
| 1 | right-down | horizontal in `sigma` |
| 2 | left-up | horizontal in `sigma` |
| 3 | left-down | horizontal in `sigma` |
| 4 | middle-right | vertical in `b` |
| 5 | middle-left | vertical in `b` |

From the `Brusselator Continuation` directory, examples for macOS are:

```sh
TASK=4 /Applications/MATLAB_R2026a.app/bin/matlab -batch "Main_spiral"
TASK=5 /Applications/MATLAB_R2026a.app/bin/matlab -batch "Main_spiral"
TASK=4 /Applications/MATLAB_R2026a.app/bin/matlab -batch "Main_target"
TASK=5 /Applications/MATLAB_R2026a.app/bin/matlab -batch "Main_target"
```

Each new run writes to a timestamped directory under
`outputs/continuation/`. It does not overwrite the root-level historical
`.mat` files. Runs started with an older copy of the driver may still finish
in a root-level `full_continuation_*` directory; move that whole directory
only after MATLAB has returned to the prompt.

## Defaults restored from saved metadata

The genuine historical middle-left files contain metadata from an earlier
vertical corrector that is absent from the current source. The restored middle
defaults follow that metadata:

- Spiral: arclength step `0.003`, `db_local = 0.01`, and six grid steps on
  each side. Middle-right uses a tracked downward crossing; middle-left uses
  the historical strongest upward crossing.
- Target: arclength step `0.002`, the same side-dependent vertical crossing
  logic, and the thesis threshold `F_thr = 8`.

The approximate direction seed is used only to form the first predictor and is
not included in the returned `p_hist_mid`. Every reported point is therefore a
corrector output.

## Historical-data warning

Do not treat all four root-level middle files as equivalent raw continuation
output:

- `middle_left.mat` contains evaluated continuation history and uses the
  left spiral threshold `0.5`.
- `target_middle_left.mat` contains evaluated continuation history but uses
  target threshold `10`, not the thesis threshold `8`.
- `middle_right.mat` is explicitly marked as a synthetic, manually assembled
  plotting curve.
- `target_middle_right.mat` is explicitly marked as a synthetic, linearly
  interpolated plotting curve.

The root-level target side branches also mix thresholds `8` and `10`. New
publication runs should regenerate every target branch at one declared
threshold before the curves are combined.

## Validation

Run the geometry tests without solving the PDE:

```matlab
addpath(genpath(pwd)); test_corrector_geometry
```

For each completed middle run:

1. Plot the accepted curve against the manually swept reference points.
2. Re-evaluate several accepted points with a finer independent vertical
   sweep, preferably `db <= 0.002`.
3. Inspect spacetime patterns below, at, and above the selected crossing.
4. Confirm that the temporal-coherence gate is not determining the crossing.
5. Only after these checks, copy the validated result into a publication-data
   directory and update figure scripts to use it explicitly.
