# Feature-Based Continuation in the One-Dimensional Brusselator

This repository contains the MATLAB code and numerical data used for the
feature-based continuation study of pattern transitions in a one-dimensional
Brusselator. The computations distinguish wave/stripe-like,
spiral/source-defect-like, and target-like regimes in the parameter plane
`(sigma,b)`.

The continued side branches are threshold level sets of features extracted
from late-time PDE simulations. The lower middle transition estimates shown
in the manuscript were obtained separately from vertical parameter sweeps and
direct inspection of spacetime plots. The historical-data notes in
[`docs/MIDDLE_BRANCH_CONTINUATION.md`](docs/MIDDLE_BRANCH_CONTINUATION.md)
explain this distinction and identify synthetic plotting files.

## Requirements

- MATLAB (the current scripts were tested with MATLAB R2026a)
- A machine with sufficient memory and runtime for repeated PDE simulations
- Optional: a Slurm cluster for the provided CCV job scripts

No absolute local paths are required by the main continuation scripts.

## Main files

- `Main_spiral.m`: spiral/source-defect transition runs
- `Main_target.m`: target transition runs
- `continuation.m`: predictor-corrector continuation driver
- `bisect_interval.m`: horizontal and vertical local sweep correctors
- `solve_brusselator_1d.m`: one-dimensional Brusselator simulation
- `feature_evaluation_left.m` and `feature_evaluation_right.m`: branch-adapted
  spiral features
- `feature_evaluation_target.m`: target feature
- `jobs/`: Slurm submission scripts
- `scripts/figures/`: manuscript figure-generation scripts
- `scripts/sweeps/`: horizontal and vertical validation sweeps
- `outputs/ccv_*_20260810/`: completed middle-branch CCV runs

The root-level `.mat` files contain the seeds and historical continuation data
used by the plotting scripts. See the warning in
[`docs/MIDDLE_BRANCH_CONTINUATION.md`](docs/MIDDLE_BRANCH_CONTINUATION.md)
before interpreting the middle-branch files.

## Running a continuation

Start MATLAB in the repository root. The environment variable `TASK` selects
the branch:

| `TASK` | Branch |
|---:|---|
| 0 | right-up |
| 1 | right-down |
| 2 | left-up |
| 3 | left-down |
| 4 | middle-right |
| 5 | middle-left |

For example, from a shell:

```sh
TASK=0 matlab -batch "Main_spiral"
TASK=2 matlab -batch "Main_target"
```

New runs are written to timestamped folders under `outputs/continuation/`.
Those local run products are ignored by Git.

On Slurm, submit the middle branches with:

```sh
mkdir -p logs
sbatch jobs/ccv_middle_spiral.sh
sbatch jobs/ccv_middle_target.sh
```

## Validation test

The corrector-geometry test does not solve the PDE:

```matlab
addpath(genpath(pwd));
test_corrector_geometry
```

## Reproducibility notes

The feature threshold, branch direction, step size, and corrector geometry are
recorded in the two main scripts and in the saved run metadata. Because mixed
and irregular patterns occur in the lower parameter region, a smooth scalar
threshold crossing should not automatically be interpreted as the intended
physical pattern transition. Representative spacetime plots and independent
sweeps should be checked as described in the middle-branch documentation.

## License

No software license has been selected yet. The files are publicly available
for inspection and reproducibility, but reuse rights are not granted beyond
those provided by applicable law.
