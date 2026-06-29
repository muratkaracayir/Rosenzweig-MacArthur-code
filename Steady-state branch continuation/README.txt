# Steady-state branch continuation

This folder contains a minimal MATLAB workflow for constructing and continuing
stationary, spatially nonhomogeneous steady states of the one-dimensional
Rosenzweig--MacArthur reaction--diffusion system

```text
u_t = d1 u_xx + u(1 - u/k) - m u v/(u + 1),
v_t = d2 v_xx - theta v + m u v/(u + 1),
```

on the interval `0 <= x <= ell*pi`, with homogeneous Neumann boundary
conditions.  The public workflow is specialized to Test case 3 / Example 3.12
from the manuscript and uses the theta-fixed parameter route.

The main tasks are:

1. construct one corrected nonconstant steady state near a steady-state
   bifurcation point;
2. continue a steady-state branch by pseudo-arclength continuation (PALC);
3. plot representative steady-state profiles in two panels.

The workflow is intentionally small.  It is not a full reproduction archive for
all manuscript figures.

## Folder layout

```text
Steady-state branch continuation/
  README.md
  config/
    case_ex312.m
  run/
    start_steady_branch.m
    continue_steady_branch.m
  continuation/
    palc_steady_branch.m
  solver/
    newton_solver_steady.m
    residual_and_jacobian_steady.m
  figures/
    plot_ss_profiles.m
  data/
    SS_Ex3.12_mode5_thetaFixed_lambda0.37726.mat
```

The `data/` folder contains one precomputed corrected steady state used by the
default continuation driver.  New `SS_*.mat` and `PALC_*.mat` files are created
there when the workflow is run.

## Requirements

The code is written for standard MATLAB syntax and avoids newer plotting
features such as `tiledlayout`, `sgtitle`, and `yline`.  It was organized with
older MATLAB compatibility in mind.

No external continuation package is required.  The steady residual and Jacobian
are assembled in `solver/residual_and_jacobian_steady.m`, and the PALC corrector
is implemented inside `continuation/palc_steady_branch.m`.

## Quick start

Start MATLAB in the top-level folder

```text
Steady-state branch continuation/
```

or navigate there using `cd`.

### 1. Construct the default starting steady state

Run:

```matlab
start_steady_branch
```

This script constructs a small-amplitude initial guess near the selected
steady-state bifurcation point, applies Newton correction, and saves the
corrected nonconstant steady state to `data/`.

The saved file has the form

```text
SS_Ex3.12_mode5_thetaFixed_lambda0.37726.mat
```

up to MATLAB's numeric formatting.  The file contains one top-level structure
named `SS`.

### 2. Continue the branch by PALC

Run:

```matlab
continue_steady_branch
```

The default driver loads

```matlab
steadyFile = 'SS_Ex3.12_mode5_thetaFixed_lambda0.37726.mat';
```

from the `data/` folder and passes it to

```matlab
PALC = palc_steady_branch(UV0, cfg0, steady, stepIdx0, opts);
```

The loaded state does not have to be the original branch anchor.  It only has
to be a corrected steady state on the desired branch, with compatible `UV`,
`cfg`, and `steady` metadata.

### 3. Plot steady-state profiles

Run:

```matlab
plot_ss_profiles
```

from the `figures/` folder, or add that folder to the MATLAB path and call the
function.  The plotter scans `data/` for `SS_*.mat` files and produces a
two-panel figure with `u(x)` on the left and `v(x)` on the right.

## Configuration file: `case_ex312.m`

The main configuration file is

```matlab
config/case_ex312.m
```

It defines the Test case 3 / Example 3.12 parameter set:

```matlab
cfg.d1  = 1.0;
cfg.d2  = 1.0;
cfg.k   = 3.0;
cfg.ell = 35.0;
cfg.Nx  = 128;
cfg.L   = cfg.ell * pi;
```

The workflow uses the theta-fixed route

```matlab
theta0 = 0.003;

cfg.route.type = 'thetaFixed';
cfg.route.fun  = @(lambdaValue) deal( ...
    theta0 * (1 + 1./lambdaValue), ...
    theta0);
```

Thus

```matlab
m(lambda)     = theta0 * (1 + 1/lambda),
theta(lambda) = theta0.
```

The field

```matlab
cfg.route.type = 'thetaFixed';
```

is kept as metadata for filenames and saved structures.  It is not a general
route selector in this public workflow.

The coexistence equilibrium is stored as

```matlab
cfg.u_star = cfg.lambda;
cfg.v_star = cfg.lambda * (1 - cfg.lambda/cfg.k) / cfg.theta;
```

The route metadata also includes the endpoint values used by the seed and by
PALC fold-suppression logic:

```matlab
cfg.route.lambdaL
cfg.route.lambdaR
cfg.route.anchorDesc
```

For the default file, these describe the mode-5 branch used in the public
example.

## The seed side and the anchor lambda value

The seed settings are stored in

```matlab
cfg.seed
```

The key fields are:

```matlab
cfg.seed.eps0
cfg.seed.lambdaBif
cfg.seed.side
cfg.seed.lambda0
```

The default bifurcation point is

```matlab
cfg.seed.lambdaBif = cfg.route.lambdaL;
```

The nearby Newton parameter value `cfg.seed.lambda0` is chosen from a grid with
spacing

```matlab
h = 1e-5;
cfg.seed.lambdaGrid = h;
```

and depends on `cfg.seed.side`:

```matlab
if cfg.seed.side >= 0
    cfg.seed.lambda0 = (floor(cfg.seed.lambdaBif/h) + 1) * h;
else
    cfg.seed.lambda0 = (ceil(cfg.seed.lambdaBif/h) - 1) * h;
end
```

Therefore:

```text
cfg.seed.side >= 0  gives  cfg.seed.lambda0 > cfg.seed.lambdaBif,
cfg.seed.side <  0  gives  cfg.seed.lambda0 < cfg.seed.lambdaBif.
```

The same `side` also enters the initial profile perturbation.  In
`start_steady_branch.m`, the initial guess is built as

```matlab
U0 = cfg0.u_star + side * eps0 * phi * q(1);
V0 = cfg0.v_star + side * eps0 * phi * q(2);
```

where

```matlab
phi(x) = cos(n*pi*x/L) = cos(n*x/ell)
```

is the Neumann mode and `q` is a null vector of the mode-reduced steady
bifurcation matrix at `cfg.seed.lambdaBif`.

If the Newton correction fails or starts on the undesired local branch, change

```matlab
cfg.seed.side
```

from `-1` to `+1`, or vice versa.

## Starting a branch: `start_steady_branch.m`

The driver

```matlab
run/start_steady_branch.m
```

performs the branch-starting step.  Its main user settings are near the top:

```matlab
caseFun   = @case_ex312;
modeIndex = 5;
lambdaForConfig = 0.377251393220232;
doSaveSS = true;
```

The value `lambdaForConfig` only instantiates the configuration.  The actual
bifurcation point and nearby Newton parameter are read from

```matlab
cfg.seed.lambdaBif
cfg.seed.lambda0
```

The script contains local helpers for:

```text
make_steady_seed
compute_seed_amp_eq
kern_vec_2x2
apply_route
snap_lambda
```

so it does not require a standalone seed-construction file or a standalone
`apply_route.m`.

The Newton correction is performed by

```matlab
[SS, UV, resNorm] = newton_solver_steady(UV0, cfg0, steady, doSaveSS);
```

A corrected state is saved only if Newton satisfies the residual tolerance and
the result is not essentially the homogeneous equilibrium.

The saved structure has fields including:

```matlab
SS.stepIdx
SS.lambda
SS.UV
SS.cfg
SS.steady
SS.resNorm
SS.amp
SS.newton
```

For a Newton-produced starting state,

```matlab
SS.stepIdx = 0;
```

This is useful for later plotting and sheet reconstruction.

## Newton solver and steady residual

The Newton solver is

```matlab
solver/newton_solver_steady.m
```

Its fixed internal settings are:

```matlab
maxIter = 100;
resTol  = 1e-10;
stepTol = 1e-10;
damping = 1.0;
ampMin  = 1e-8;
```

It calls

```matlab
residual_and_jacobian_steady(UV, cfg)
```

to assemble the steady residual and sparse Jacobian.

The residual routine uses the endpoint grid

```matlab
x = linspace(0, cfg.L, cfg.Nx)
```

and imposes homogeneous Neumann boundary conditions through a DCT-I spectral
Laplacian.  It returns the stacked residual

```matlab
R = [Ru; Rv]
```

with

```matlab
Ru = d1*u_xx + u*(1-u/k) - m*u*v/(u+1),
Rv = d2*v_xx - theta*v + m*u*v/(u+1).
```

When requested, it also returns the sparse Jacobian with respect to the stacked
unknown vector `[u; v]`.

## Continuing a branch: `continue_steady_branch.m`

The public continuation driver is

```matlab
run/continue_steady_branch.m
```

It adds the `config/`, `continuation/`, and `solver/` folders to the MATLAB
path, loads one saved `SS` file from `data/`, extracts

```matlab
UV0    = SS.UV;
cfg0   = SS.cfg;
steady = SS.steady;
```

and calls the PALC engine.

The default input file is

```matlab
steadyFile = 'SS_Ex3.12_mode5_thetaFixed_lambda0.37726.mat';
```

The default starting global step index is

```matlab
stepIdx0 = 0;
```

If continuing from a later saved PALC point, change `steadyFile` and set
`stepIdx0` consistently with the loaded point.

The default continuation settings in the driver are:

```matlab
opts = struct();
opts.nSteps = 100;
opts.ds     = -1e-2;

opts.saveEvery     = 10;
opts.verbose       = true;
opts.savePALCAtEnd = true;
```

The sign and magnitude of `opts.ds` are important:

```text
abs(opts.ds)  controls the initial pseudo-arclength step size,
sign(opts.ds) controls the initial lambda direction.
```

In `palc_steady_branch.m`, the initial tangent is oriented so that:

```text
opts.ds > 0  initially increases lambda,
opts.ds < 0  initially decreases lambda.
```

After folds are encountered, lambda need not remain monotone; PALC follows the
branch in pseudo-arclength.

The option

```matlab
opts.nSteps = 100;
```

sets the maximum number of accepted PALC steps attempted by the driver.  The
PALC routine may stop earlier if a stopping criterion is triggered.

## PALC engine: `palc_steady_branch.m`

The PALC engine is

```matlab
continuation/palc_steady_branch.m
```

with interface

```matlab
PALC = palc_steady_branch(UV0, cfg0, steady, stepIdx0, opts)
```

It does not call `newton_solver_steady.m`.  Instead, it performs its own
augmented Newton correction internally and calls

```matlab
residual_and_jacobian_steady(UV, cfg)
```

directly.

The most important PALC options are:

```matlab
opts.nSteps
opts.ds
opts.dsMin
opts.dsMax
opts.maxCorr
opts.rTol
opts.sTol
opts.saveEvery
opts.savePALCAtEnd
```

Their roles are:

```text
opts.nSteps        maximum number of continuation steps;
opts.ds            signed initial pseudo-arclength step;
opts.dsMin         minimum allowed absolute step size;
opts.dsMax         maximum allowed absolute step size;
opts.maxCorr       maximum augmented Newton iterations per correction;
opts.rTol          residual tolerance for the PALC corrector;
opts.sTol          arclength constraint tolerance;
opts.saveEvery     save one SS file every this many accepted steps;
opts.savePALCAtEnd save the final PALC structure at termination.
```

Additional default diagnostics and stopping controls include:

```matlab
opts.specModesK
opts.specJumpThreshL1
opts.stopOnAbruptSpectrum
opts.stopOnDupLambdaMode
opts.dupTolLambda
opts.stopOnCollapsedEq
opts.collapseAmpThresh
opts.collapsePersistSteps
opts.stopOnDsUnderflow
opts.tailCheckEnabled
opts.tailFracStart
opts.tailEnergyThresh
```

By default, the code also uses endpoint values from

```matlab
cfg.route.lambdaL
cfg.route.lambdaR
```

to suppress artificial fold reports extremely close to known branch endpoints.
The relevant tolerances are

```matlab
opts.foldSuppLamTol
opts.foldSuppAmpTol
```

The PALC output structure stores branch-level arrays such as:

```matlab
PALC.lambdas
PALC.resNorms
PALC.ampEq
PALC.ampEqNorm
PALC.diag
PALC.folds
PALC.stopReason
PALC.stopMeta
```

At the end, the routine also adds a compact table

```matlab
PALC.tbl
```

with columns:

```text
step, lambda, ampEq, ampEqNorm, modeDominant, resNorm
```

when the PALC structure is saved.

## Files generated by PALC

PALC saves individual steady states in `data/` with filenames of the form

```text
SS_thetaFixed_step00010_mode5_plus_lambda0.42.mat
SS_thetaFixed_step00020_mode5_minus_lambda0.55.mat
```

The exact mode, sign label, step number, and lambda value depend on the run.

Each PALC-saved `SS` structure includes

```matlab
SS.stepIdx
SS.lambda
SS.UV
SS.cfg
SS.steady
SS.ampEq
SS.ampEqNorm
```

The field

```matlab
SS.stepIdx
```

records the global continuation step and is used by the profile plotter in
advanced mode to reconstruct continuation order.

At termination, PALC saves a branch-level file of the form

```text
PALC_thetaFixed_step00100_end_lambda0.52.mat
```

when

```matlab
opts.savePALCAtEnd = true;
```

## Profile plotter: `plot_ss_profiles.m`

The profile plotter is

```matlab
figures/plot_ss_profiles.m
```

It scans `data/` for files named

```text
SS_*.mat
```

and expects each valid file to contain a top-level structure named `SS` with
at least

```matlab
SS.lambda
SS.UV
SS.cfg
```

The main user settings are:

```matlab
targetLambda = 0.375;
advancedMode = true;
lambdaTol = 1e-3;
sheetTurnTol = 1e-10;
```

The plotter produces a two-panel figure:

```text
left panel:  u(x),
right panel: v(x).
```

### Simple mode

Set

```matlab
advancedMode = false;
```

The script selects the single file whose `SS.lambda` is closest to
`targetLambda`.  It plots this profile only if

```matlab
abs(SS.lambda - targetLambda) <= lambdaTol.
```

Otherwise it stops with a clear tolerance error.

### Advanced mode

Set

```matlab
advancedMode = true;
```

The script loads all valid `SS_*.mat` files, orders them by `SS.stepIdx` when
all step indices are finite and unique, splits the ordered lambda sequence into
sheets by detecting turning points in lambda, and selects the profile closest
to `targetLambda` on each sheet.

A sheet is plotted only if its selected profile satisfies

```matlab
abs(SS.lambda - targetLambda) <= lambdaTol.
```

If `SS.stepIdx` is missing or not unique for all files, the script falls back
to file timestamps and warns that sheet detection may be unreliable.

For the most reliable advanced-mode plots, use a clean `data/` folder containing
only steady states from the branch being visualized.

## Changing the example mode or branch

To start from another mode, for example mode 7, change the top of
`start_steady_branch.m`:

```matlab
modeIndex = 7;
```

Then update the relevant configuration metadata in `case_ex312.m`:

```matlab
cfg.route.lambdaL
cfg.route.lambdaR
cfg.route.anchorDesc
cfg.seed.lambdaBif
cfg.seed.side
```

The PALC routine itself is mode-agnostic.  It can start from any corrected
steady state on the desired branch, not only from the initial anchor produced
by `start_steady_branch.m`.

To continue from a different saved state, change the hardcoded filename in
`continue_steady_branch.m`:

```matlab
steadyFile = '...';
```

and update

```matlab
stepIdx0
```

when continuing from a saved PALC point whose step index is not zero.

## Recommended contents of `data/` for a public repository

For a minimal public repository, include only one default corrected starting
profile:

```text
data/
  SS_Ex3.12_mode5_thetaFixed_lambda0.37726.mat
```

Do not include generated PALC output files unless the repository is intended to
serve as a results archive.  The files

```text
SS_thetaFixed_step*.mat
PALC_thetaFixed_*.mat
```

are generated by running the continuation workflow.

## Common adjustments

### Change the initial PALC direction

In `continue_steady_branch.m`, change the sign of

```matlab
opts.ds
```

Use a positive value to initially increase lambda and a negative value to
initially decrease lambda.

### Change the continuation length

In `continue_steady_branch.m`, change

```matlab
opts.nSteps
```

For a quick test, use a small value such as `10`.  For a longer branch, increase
it.

### Save every accepted PALC point

Set

```matlab
opts.saveEvery = 1;
```

This produces many `SS_thetaFixed_step*.mat` files and is useful before running
`plot_ss_profiles.m` in advanced mode.

### Avoid too many saved files

Use a larger value such as

```matlab
opts.saveEvery = 10;
```

or set it to `0` to save only critical/final states according to the PALC save
logic.

### Make profile selection stricter or looser

In `plot_ss_profiles.m`, adjust

```matlab
lambdaTol
```

The default is strict:

```matlab
lambdaTol = 1e-3;
```

Increase it if no saved profiles are sufficiently close to `targetLambda`.

## Troubleshooting

### Newton converges to the homogeneous equilibrium

Try changing

```matlab
cfg.seed.side
```

or increasing/decreasing

```matlab
cfg.seed.eps0
```

in `case_ex312.m`.

### Newton fails from the initial seed

Try reducing

```matlab
cfg.seed.eps0
```

or switching the sign of

```matlab
cfg.seed.side
```

### PALC stops early because of duplicate lambda values

The default duplicate-lambda behavior is conservative.  To disable it, set

```matlab
opts.stopOnDupLambdaMode = 'none';
```

in `continue_steady_branch.m`.

### Advanced profile plotting gives a sheet warning

Make sure the `SS_*.mat` files were generated by the current workflow and
contain finite, unique `SS.stepIdx` values.  Also consider cleaning the `data/`
folder so it contains only files from the branch being plotted.
