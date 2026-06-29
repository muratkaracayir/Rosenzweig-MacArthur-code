# Rosenzweig--MacArthur reaction--diffusion workflows

This repository contains a minimal public MATLAB companion codebase for
selected numerical workflows associated with a study of the one-dimensional
Rosenzweig--MacArthur reaction--diffusion predator--prey system

```text
u_t = d1 u_xx + u(1 - u/k) - m u v/(u + 1),
v_t = d2 v_xx - theta v + m u v/(u + 1),
```

on an interval with homogeneous Neumann boundary conditions.  The code is
organized as four independent workflow folders.  Each folder has its own
`README.md` and `README.txt` with detailed instructions.

The repository is intentionally compact.  It contains the core scripts needed
to run representative computations, inspect the numerical output, and reproduce
basic diagnostic plots.  It is not a full private workbase and does not include
every parameter sweep, table generator, manuscript-figure script, or auxiliary
one-off utility used during development.

## Repository layout

```text
<repository root>/
  README.md
  README.txt

  Constant Equilibria/
    README.md
    README.txt
    config/
    run/
    solver/
    figures/
    data/

  Homogeneous limit cycle/
    README.md
    README.txt
    config/
    ode_solver/
    pde_solver/
    postprocess/
    run/
    figures/
    data/

  Hopf branch continuation/
    README.md
    README.txt
    config/
    run/
    continuation/
    solver/
    floquet/
    utils/
    figures/
    data/

  Steady-state branch continuation/
    README.md
    README.txt
    config/
    run/
    continuation/
    solver/
    figures/
    data/
```

The folder names are part of the public workflow organization.  Because several
folder names contain spaces, use quotes when changing directory from the MATLAB
command line, for example

```matlab
cd('Homogeneous limit cycle')
```

## Workflow summary

### 1. Constant Equilibria

This folder contains a Strang-splitting time-stepper for checking convergence
toward spatially homogeneous equilibria in a representative parameter regime.
The public entry point is

```matlab
run/run_constant_equilibrium
```

and the quick-look plotter is

```matlab
figures/plot_constant_history
```

The solver uses an endpoint grid, homogeneous Neumann boundary conditions, and a
DCT-I diffusion step.

### 2. Homogeneous limit cycle

This folder computes a spatially homogeneous ODE reference limit cycle and then
runs the PDE time-stepper to compare the late-time spatial mean with that ODE
reference orbit.  The main steps are

```matlab
ode_solver/compute_reference_orbit
run/run_homogeneous_periodic
```

followed by optional plots from

```matlab
figures/plot_ode_reference
figures/plot_pdemean_vs_ode
```

The default public run corresponds to Test case 2 / Example 2.6.2 with
`lambda = 5`, `Nx = 128`, and `dt = 0.001`.

### 3. Hopf branch continuation

This folder contains a compact Hopf-branch continuation workflow for spatially
nonhomogeneous periodic orbits.  It includes construction of an initial
small-amplitude Hopf orbit, Newton correction, branch continuation, selected
single-orbit Floquet analysis, and quick-look plots.

The main scripts are

```matlab
run/compute_anchor_orbit
run/run_continue_branch
run/run_floquet_orbit
```

with plotting scripts

```matlab
figures/plot_upo_colormaps
figures/plot_upo_branch
```

The public workflow uses the theta-fixed route

```text
theta(lambda) = theta0,
m(lambda)     = theta0*(1 + 1/lambda).
```

Other private parameter routes are intentionally not included.

### 4. Steady-state branch continuation

This folder contains the steady-state analogue of the continuation workflow. It
constructs a corrected nonconstant steady state near a steady-state bifurcation
point, continues the branch by pseudo-arclength continuation (PALC), and plots
representative steady profiles.

The main scripts are

```matlab
run/start_steady_branch
run/continue_steady_branch
```

and the PALC engine is

```matlab
continuation/palc_steady_branch
```

The profile plotter is

```matlab
figures/plot_ss_profiles
```

The default public example corresponds to Test case 3 / Example 3.12 and uses
the theta-fixed route.

## Requirements

The code is written for MATLAB.  The scripts avoid newer plotting commands such
as `tiledlayout`, `sgtitle`, and `yline`, and were organized with older MATLAB
compatibility in mind.

No external continuation package is required.  Newton correction, PALC, residual
assembly, Jacobian assembly, and the DCT-I diffusion steps are implemented in
the repository.

One part of the homogeneous limit-cycle postprocessing uses MATLAB spline
utilities such as `csape` and `fnval`.  If those functions are unavailable in
your MATLAB installation, the corresponding script will issue an error.  The
other workflows do not rely on those spline utilities.

## Quick start

Each workflow is meant to be run from its own folder.  For example,

```matlab
cd('Constant Equilibria')
run/run_constant_equilibrium
```

or

```matlab
cd('Steady-state branch continuation')
run/start_steady_branch
run/continue_steady_branch
```

The driver scripts add their local subfolders to the MATLAB path at runtime, so
you should not need to edit MATLAB's permanent path.

For details, open the workflow-specific README file in the folder you want to
run.

## Data folders and generated files

Each workflow has its own `data/` folder.  These folders may contain one or more
small precomputed files so that the default plotting or continuation scripts can
run immediately.  Additional `.mat` files are generated when the workflows are
run.

Generated output files can usually be deleted and regenerated from the scripts.
Large branch computations may create many files, especially in the continuation
workflows.  The public repository keeps only representative data needed for the
default examples.

## Reproducibility scope

This codebase is designed to make the main numerical workflows readable and
runnable.  It is not intended to be a complete archival snapshot of the full
research workbase.  In particular, it does not include every exploratory script,
all intermediate data, every manuscript figure generator, or every table
construction utility.

The workflow READMEs describe which scripts are entry points, which files are
expected as input, and which output files are generated.

## Citation

If you use this code, please cite the associated manuscript.  If a formal
software citation file is added to the repository, please also follow the
instructions in that file.

## Suggested reading order for new users

A practical way to explore the repository is:

1. Run `Constant Equilibria/` to see the basic time-stepper and data layout.
2. Run `Homogeneous limit cycle/` to compare PDE mean dynamics with an ODE
   reference orbit.
3. Run `Hopf branch continuation/` to inspect periodic-orbit continuation.
4. Run `Steady-state branch continuation/` to inspect steady-state PALC
   continuation and profile plotting.

The four folders are independent enough to be inspected separately, but this
order follows the increasing complexity of the numerical workflows.
