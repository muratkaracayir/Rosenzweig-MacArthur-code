# Constant-equilibrium time stepping

This workflow folder contains a minimal MATLAB implementation for time
stepping the one-dimensional diffusive Rosenzweig--MacArthur predator--prey
system and checking convergence toward a spatially homogeneous equilibrium in
a representative parameter regime.

This folder is one component of the larger code companion repository. The
planned top-level repository contains separate workflow folders for constant
equilibria, homogeneous limit cycles, nonhomogeneous limit cycles, and steady
states. The present README describes only the Constant equilibria/ workflow.

The code here is deliberately compact. It contains the core routines needed to
run the constant-equilibrium example, save the output, and make one quick-look
plot. It does not include the full private workbase, parameter sweeps,
convergence studies, or publication-figure scripts.

## Folder structure

Constant equilibria/
  README.md

  config/
    case_ex261.m

  run/
    run_constant_equilibrium.m

  solver/
    run_strang.m
    createIC.m
    strang_split.m
    diffusion_step_DCT.m
    reaction_step_MPRK2.m
    dct1_endpoints.m

  figures/
    plot_constant_history_quicklook.m

  data/
    saved output files

## Quick start

Open MATLAB and change directory to the Constant equilibria/ folder, namely
the folder containing this README file. Then run

run/run_constant_equilibrium

The driver adds the required config/ and solver/ folders to the MATLAB path
automatically.

With the default settings, the run writes a .mat file to data/ with name

Constant_Ex2.6.1_lambda_10_Nx128_dt0.0001_k3.mat

After running the default simulation, create a one-panel diagnostic plot by
running

figures/plot_constant_history_quicklook

The plotting script visualizes the time histories of the spatial envelopes

u_min(t), u_max(t), v_min(t), v_max(t)

for one representative initial condition. It also saves the resulting MATLAB
figure file in the figures/ folder.

## Main files

- config/case_ex261.m defines the parameter set and homogeneous equilibrium
  data used by the example.
- run/run_constant_equilibrium.m is the user-facing driver for this workflow.
- solver/run_strang.m runs the time-stepping experiment for the selected
  initial conditions.
- solver/strang_split.m performs one Strang-splitting time step.
- solver/diffusion_step_DCT.m applies the diffusion substep using a DCT-I
  representation compatible with homogeneous Neumann boundary conditions.
- solver/reaction_step_MPRK2.m applies the reaction substep.
- solver/dct1_endpoints.m implements the endpoint-grid DCT-I transform used
  by the diffusion step.
- figures/plot_constant_history_quicklook.m produces a simple diagnostic
  figure from the default saved output.

## Notes on output

The default driver uses

run.saveResults = true;
run.onlyICMode  = '';

Thus the full selected initial-condition suite is run and the result is saved
in data/. The option run.onlyICMode may be changed later for quick
interactive tests with a single initial-condition family, but the default
public setting is chosen so that a new user immediately obtains a saved output
file for plotting.

## Notes on scope

This folder is a curated core workflow, not a full reproducibility archive for
all computations in the manuscript. It is intended to make the main
constant-equilibrium time-stepping code readable, runnable, and easy to inspect.
