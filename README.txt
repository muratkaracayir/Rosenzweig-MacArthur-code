# Hopf branch continuation

This folder contains a minimal public MATLAB workflow for computing spatially
nonhomogeneous periodic orbits of the one-dimensional diffusive
Rosenzweig--MacArthur predator--prey system by Hopf-branch continuation.

The code is intended as a compact companion to the manuscript, not as a full
private workbase. It contains the core continuation, Newton-correction, selected
single-orbit Floquet analysis, and quick-look plotting routines. It does not
include manuscript table generators, full figure-production scripts, direct PDE
time-stepping verification, or branchwise Floquet scans.

## Folder layout

```text
Hopf branch continuation/
  README.md
  config/
    case_ex262.m
  run/
    compute_anchor_orbit.m
    run_continue_branch.m
    run_floquet_orbit.m
  continuation/
    prepare_data.m
    continue_branch.m
  solver/
    newton_solver.m
    residual_and_jacobian.m
  floquet/
    make_floqCfg.m
    compute_floquet_IFRK4.m
    floquet_buildM_IFRK4.m
    floquet_processM.m
    floquet_orbitEval_trig.m
  utils/
    merge_UPO_files.m
  figures/
    plot_upo_colormaps.m
    plot_upo_branch.m
  data/
```

The `data/` folder is used for generated `.mat` files. It may be empty when the
repository is first cloned.

## Model and test case

The workflow uses the nondimensional diffusive Rosenzweig--MacArthur system

```text
u_t = d1 u_xx + u(1 - u/k) - m u v/(u + 1),
v_t = d2 v_xx - theta v + m u v/(u + 1),
```

on a one-dimensional interval with homogeneous Neumann boundary conditions.
The public configuration `config/case_ex262.m` corresponds to Example 2.6.2 in
the manuscript. The continuation parameter is

```text
lambda = theta/(m - theta),
```

and this public workflow uses only the theta-fixed route

```text
theta(lambda) = theta0,
m(lambda)     = theta0*(1 + 1/lambda).
```

Other parameter routes used during private experimentation are intentionally not
included.

## Requirements

The scripts are written for MATLAB and avoid newer plotting features such as
`tiledlayout`, `sgtitle`, and `yline`. They are intended to remain compatible
with older MATLAB versions used in the computations.

No external MATLAB packages are required.

## Basic workflow

### 1. Compute an anchor periodic orbit

Edit the user settings in

```text
run/compute_anchor_orbit.m
```

if needed, then run the script. The default settings are

```matlab
lambdaHopf = 0.5;
modeIndex  = 2;
doSave     = true;
```

The script constructs a small-amplitude Hopf initial guess using
`continuation/prepare_data.m`, then Newton-corrects it with
`solver/newton_solver.m`.

With `doSave = true`, the corrected orbit is saved in `data/` as a file of the
form

```text
UPO_Ex2.6.2_mode2_thetaFixed_lambda0.5001.mat
```

The exact lambda value in the filename is the corrected orbit parameter value.
A saved UPO file contains the variables

```text
UV, T, lambda, resNorm, hopf, ampUV
```

where `UV` is the packed space-time orbit, `T` is its period, and `ampUV` is the
amplitude diagnostic used by the branch plotter.

### 2. Continue the branch

Edit the user settings in

```text
run/run_continue_branch.m
```

The main choices are

```matlab
hopfLambda     = 0.5;
modeIndex      = 2;
startUPOLambda = 0.5001;

nSteps      = 50;
dLambda     = 0.01;
keepInitial = true;
```

The driver locates the saved starting UPO in `data/` and calls
`continuation/continue_branch.m`. The continuation result is saved in `data/` as
a branch file of the form

```text
Branch_Ex2.6.2_mode2_thetaFixed_lambda0.5001_50Steps_dL0.01.mat
```

This branch file contains one top-level variable:

```text
UPOs
```

with entries

```text
UPOs(1), UPOs(2), ..., UPOs(end).
```

Each entry stores one corrected periodic orbit and its diagnostics.

### 3. Optionally merge individual UPO files

If individual `UPO_*.mat` files have been generated and you want a compact master
branch file, run

```matlab
utils/merge_UPO_files
```

or call it with explicit options, for example

```matlab
merge_UPO_files( ...
    'filePattern', 'UPO_Ex2.6.2_mode2_thetaFixed_lambda*.mat', ...
    'outName', 'UPOMaster_Ex2.6.2_mode2_thetaFixed.mat');
```

The saved master file contains only one top-level variable:

```text
UPOs
```

The entries are sorted by lambda, and duplicate lambda values are removed up to
the specified tolerance.

### 4. Plot quick-look figures

Two lightweight plotting scripts are included in `figures/`.

#### Single-orbit colormaps

```text
figures/plot_upo_colormaps.m
```

This script loads one hardcoded UPO file from `data/` and plots a two-panel
colormap figure:

```text
left:  u(x,t)
right: v(x,t)
```

The default input file is

```matlab
upoFile = 'UPO_Ex2.6.2_mode2_thetaFixed_lambda2.mat';
```

Edit this line to inspect a different saved orbit.

#### Branch diagram

```text
figures/plot_upo_branch.m
```

This script loads one hardcoded branch/master file from `data/` and plots

```text
left:  amplitude vs lambda
right: period T vs lambda
```

The default input file is

```matlab
masterFile = 'UPOMaster_Ex2.6.2_mode2_thetaFixed.mat';
```

The script only requires the loaded file to contain a top-level `UPOs` variable.
Thus it can also be used with a `Branch_*.mat` file by changing `masterFile`.

Both figure scripts set `doSave = false` by default. Set `doSave = true` to save
`.fig` and `.png` copies.

### 5. Run Floquet analysis for one orbit

Edit the selected UPO filename and Floquet settings in

```text
run/run_floquet_orbit.m
```

For example:

```matlab
upoFile = 'UPO_Ex2.6.2_mode2_thetaFixed_lambda0.5001.mat';

doSave = true;

floqOverrides = struct( ...
    'mSub', 8, ...
    'nStore', 16, ...
    'verbose', false);
```

The driver loads the selected orbit, calls

```matlab
compute_floquet_IFRK4
```

and prints a compact stability summary, including the leading nontrivial
multiplier, the number of unstable directions, and a Floquet-stable or
Floquet-unstable verdict.

With `doSave = true`, the output is saved in `data/` as a file of the form

```text
Floquet_Ex2.6.2_mode2_thetaFixed_lambda0.5001.mat
```

The saved file contains

```text
floquet, cfgFloq, sourceFile, upoFile, lambda
```

The `floquet` struct contains method metadata, monodromy-building diagnostics,
stored multipliers, the identified trivial multiplier, stability counts, and
quality-control flags.

## Data conventions

A single saved UPO file contains

```text
UV, T, lambda, resNorm, hopf, ampUV
```

A branch or master file contains

```text
UPOs
```

where each `UPOs(k)` entry contains at least

```text
UV, T, lambda, hopf, ampUV
```

and usually also

```text
resNorm
```

This convention is intentionally simple: plotting and postprocessing scripts can
work directly from the `UPOs` array without requiring additional top-level
branch arrays.

## Notes

- This repository focuses on the theta-fixed Hopf branch from Test case 2.
- The quick-look plots are included for inspection, not for reproducing the
  manuscript figures exactly.
- Floquet analysis is provided for one selected orbit at a time. A branchwise
  Floquet scan is intentionally omitted from this minimal public workflow.
- The full monodromy matrix used in the Floquet computation can be large. Higher
  spatial or temporal resolutions may require substantial memory.
