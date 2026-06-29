# Homogeneous limit cycles

This workflow contains the core MATLAB codes for the homogeneous limit-cycle
computations in **Test case 2** of the manuscript.  The test case corresponds
to Example 2.6.2 in the reference setting used in the code.  The default
public run uses

```text
lambda = 5,
Nx     = 128,
dt     = 0.001.
```

For this parameter set, the spatially homogeneous Hopf value is

```text
lambda_0^H = (k - 1)/2 = 8,
```

so the default value `lambda = 5` lies in the range where the homogeneous
ODE dynamics have a stable limit cycle.

## Folder structure

```text
Homogeneous limit cycles/
  README.md

  config/
    case_ex262.m

  ode_solver/
    compute_reference_orbit.m

  pde_solver/
    run_strang.m
    createIC.m
    strang_split.m
    diffusion_step_DCT.m
    reaction_step_MPRK2.m
    dct1_endpoints.m

  postprocess/
    build_ref_orbit_interpolant.m
    eval_ref_orbit.m
    postprocess_homogeneous_periodic_run.m

  run/
    run_homogeneous_periodic.m

  figures/
    plot_ode_reference.m
    plot_pdemean_vs_ode.m

  data/
    RefOrbit_*.mat
    PDEOrbit_*.mat
```

The `data/` folder is used for generated `.mat` files.  It may be empty when
the repository is first downloaded.

## Main workflow

Run the scripts from MATLAB in the following order.

### 1. Compute the ODE reference orbit

```matlab
cd('Homogeneous limit cycles/ode_solver')
ref = compute_reference_orbit();
```

This computes one period of the spatially homogeneous ODE limit cycle and,
by default, saves the output structure `ref` in the neighboring `data/`
folder as

```text
data/RefOrbit_Ex2.6.2_lambda_5.mat
```

This file is required by the PDE driver.

### 2. Run the PDE time-stepper and postprocess the result

```matlab
cd('../run')
run_homogeneous_periodic
```

The driver loads the ODE reference orbit, runs the Strang-splitting PDE
solver, computes postprocessing diagnostics, and saves the processed output
in `data/`.  With the default public settings, the expected output filename
is

```text
data/PDEOrbit_Ex2.6.2_lambda_5_Nx128_dt0.001.mat
```

If a file with the same name already exists, the driver may create a
versioned filename such as

```text
data/PDEOrbit_Ex2.6.2_lambda_5_Nx128_dt0.001_v2.mat
```

The saved PDE file contains, among other variables,

```text
cfg
run
results
postOpts
post
refFileName
refFile
```

The top-level driver owns file saving.  The solver and postprocessor
functions themselves do not save files.

## Quick-look plots

After the two workflow steps above, the following lightweight plotters can
be run from the `figures/` folder.

### ODE reference cycle

```matlab
cd('../figures')
plot_ode_reference
```

This loads

```text
data/RefOrbit_Ex2.6.2_lambda_5.mat
```

and plots the ODE reference cycle in the `(u,v)` phase plane.

### PDE mean tail versus ODE reference cycle

```matlab
plot_pdemean_vs_ode
```

This loads

```text
data/RefOrbit_Ex2.6.2_lambda_5.mat
data/PDEOrbit_Ex2.6.2_lambda_5_Nx128_dt0.001.mat
```

and plots the late-time PDE spatial-mean trajectory together with the ODE
reference cycle.  The ODE orbit is plotted with a dashed line, and the PDE
mean tail is plotted with a solid line.

The default simulation index is `simIndex = 2`.  To plot a different stored
initial condition, call for example

```matlab
plot_pdemean_vs_ode(1)
```

## File roles

- `config/case_ex262.m` defines the Test case 2 parameters for a chosen
  value of `lambda`.
- `ode_solver/compute_reference_orbit.m` computes and saves the ODE
  reference limit cycle.
- `run/run_homogeneous_periodic.m` is the public entry point for the PDE
  time-stepping workflow.
- `pde_solver/run_strang.m` advances the PDE system and records late-time
  diagnostics needed by the homogeneous limit-cycle workflow.
- `postprocess/postprocess_homogeneous_periodic_run.m` computes period and
  reference-orbit comparison diagnostics from the PDE output.
- `figures/plot_ode_reference.m` and `figures/plot_pdemean_vs_ode.m` provide
  simple one-panel visual checks.

## Notes

- The public entry scripts use relative paths based on the location of the
  running script, so they do not require editing MATLAB's permanent path.
  In particular, `run/run_homogeneous_periodic.m` adds `config/`,
  `ode_solver/`, `pde_solver/`, and `postprocess/` automatically.
- The PDE solver uses an endpoint grid with homogeneous Neumann boundary
  conditions and a DCT-I diffusion step.
- The reference-orbit interpolant uses MATLAB spline utilities such as
  `csape` and `fnval`.
- These scripts are intended as a minimal public code companion.  They are
  not the full private workbase used to create every manuscript figure and
  table.
