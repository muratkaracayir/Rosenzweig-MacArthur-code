%% run_homogeneous_periodic.m
% Top-level driver for the homogeneous limit-cycle PDE workflow.
%
% This script assumes the following workflow-folder layout:
%
%   Homogeneous limit cycles/
%     config/
%       case_ex262.m
%     ode_solver/
%       compute_reference_orbit.m
%     pde_solver/
%       run_strang.m and helper routines
%     postprocess/
%       build_ref_orbit_interpolant.m
%       postprocess_homogeneous_periodic_run.m
%     run/
%       run_homogeneous_periodic.m
%     data/
%       RefOrbit_*.mat and saved PDE output
%
% Before running this script, first run compute_reference_orbit.m to create
% the ODE reference-orbit file in data/. This driver then loads that
% reference orbit, runs the PDE time stepper, compares the PDE output with
% the reference orbit, and saves the processed result.
%
% For this example, the homogeneous Hopf value is lambda_0^H = 8. To observe
% convergence to the homogeneous limit cycle, choose 0 < lambda < 8.

clear; clc;

%% ------------------------------------------------------------
% 0) Paths relative to this run script
%% ------------------------------------------------------------
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
rootDir  = fileparts(thisDir);

dataDir = fullfile(rootDir, 'data');
if exist(dataDir, 'dir') ~= 7
    mkdir(dataDir);
end

addpath(fullfile(rootDir, 'config'));
addpath(fullfile(rootDir, 'ode_solver'));
addpath(fullfile(rootDir, 'pde_solver'));
addpath(fullfile(rootDir, 'postprocess'));

%% ------------------------------------------------------------
% 1) Case and parameter
%% ------------------------------------------------------------
lambda = 5;                  % Choose 0 < lambda < 8 for the limit-cycle run.
cfg = case_ex262(lambda);

%% ------------------------------------------------------------
% 2) Load ODE reference orbit
%% ------------------------------------------------------------
refFileName = 'RefOrbit_Ex2.6.2_lambda_5.mat';
refFile = fullfile(dataDir, refFileName);

S = load(refFile);

if ~isfield(S, 'ref')
    error('The reference-orbit MAT file must contain a variable named ref.');
end

ref = S.ref;

% Check that the loaded reference orbit matches the selected lambda.
if isfield(ref, 'lambda')
    if abs(ref.lambda - lambda) > 1e-12 * max(1, abs(lambda))
        error(['Reference-orbit file lambda (%.12g) does not match ', ...
               'driver lambda (%.12g).'], ref.lambda, lambda);
    end
end

% Build the helper structure used by the PDE stopping and postprocessing
% diagnostics.
refOrbit = build_ref_orbit_interpolant(ref);

%% ------------------------------------------------------------
% 3) PDE time-stepping parameters
%% ------------------------------------------------------------
run.Nx          = 128;
run.dt          = 1e-3;
run.Tfinal      = 1000;
run.outputEvery = 100;
run.verbose     = true;

% Initial-condition control:
%   ''                         -> run the full base suite
%   'equilibriumPerturbation'   -> run only this IC family
%   'smoothRandom'              -> run only this IC family
%   'largeSmooth'               -> run only this IC family
run.exploreIC    = false;
run.onlyICMode   = '';
run.onlyICSample = [];
run.saveResults  = true;

run.seed = 1;

%% ------------------------------------------------------------
% 4) Stopping settings
%% ------------------------------------------------------------
run.useStopping = true;

nStopPerPeriod = 100;
run.stopEvery  = max(1, round(ref.T / (nStopPerPeriod * run.dt)));

run.homTol          = 1e-100;
run.orbitTol        = 1e-80;
run.minStopTime     = 50;
run.nStopPass       = 3;
run.tailPeriods     = 3;
run.storeDenseTail  = true;
run.homScaleFloor   = 1;
run.orbitScaleFloor = 1;

%% ------------------------------------------------------------
% 5) Reference-orbit settings
%% ------------------------------------------------------------
run.refOrbit = refOrbit;

% false -> use the dense stored one-period sample in refOrbit
% true  -> use spline evaluation on a uniform phase-search grid
run.useRefInterpolant = true;
run.refPhaseSearchN   = 1000000;

%% ------------------------------------------------------------
% 6) Postprocessing options
%% ------------------------------------------------------------
postOpts.preferDenseTail = true;

% Period search near the ODE reference period.
postOpts.periodSearchRel = [0.8, 1.2];
postOpts.coarseSearchN   = 200;

% Reference-phase search.
postOpts.refSearchN = 200;

% Interpolation and scaling.
postOpts.interpMethod  = 'pchip';
postOpts.scaleFloor    = 1;
postOpts.minOverlapRel = 0.75;

% Convergence-time options based on stabilization of orbit-distance dips.
postOpts.convNDips     = 5;
postOpts.convRelTol    = 0.01;
postOpts.convUseInterp = true;

%% ------------------------------------------------------------
% 7) Run PDE simulation suite and postprocess
%% ------------------------------------------------------------
results = run_strang(cfg, run);
post = postprocess_homogeneous_periodic_run(results, postOpts);

%% ------------------------------------------------------------
% 8) Remove large fields already retained in post
%% ------------------------------------------------------------
if isfield(results.sim, 'denseTail')
    results.sim = rmfield(results.sim, 'denseTail');
end

if isfield(results.sim, 'orbitDistHistory')
    results.sim = rmfield(results.sim, 'orbitDistHistory');
end

%% ------------------------------------------------------------
% 9) Save output for full-suite runs
%% ------------------------------------------------------------
if isfield(run, 'saveResults') && run.saveResults && isempty(run.onlyICMode)
    saveName = sprintf('PDEOrbit_%s_lambda_%g_Nx%d_dt%g.mat', ...
        cfg.caseName, cfg.lambda, run.Nx, run.dt);

    savePath = fullfile(dataDir, saveName);

    if exist(savePath, 'file') == 2
        [~, baseName, ext] = fileparts(saveName);
        version = 2;

        while true
            candidateName = sprintf('%s_v%d%s', baseName, version, ext);
            candidatePath = fullfile(dataDir, candidateName);

            if exist(candidatePath, 'file') ~= 2
                savePath = candidatePath;
                break;
            end

            version = version + 1;
        end
    end

    save(savePath, 'cfg', 'run', 'results', 'postOpts', 'post', ...
        'refFileName', 'refFile');

    fprintf('\nSaved results to:\n%s\n', savePath);
end

%% ------------------------------------------------------------
% 10) Compact command-window summary
%% ------------------------------------------------------------
fprintf('\n');
fprintf('========================================\n');
fprintf('Completed run_homogeneous_periodic\n');
fprintf('Case          : %s\n', cfg.caseName);
fprintf('lambda        : %.12g\n', cfg.lambda);
fprintf('Nx            : %d\n', run.Nx);
fprintf('dt            : %.12g\n', run.dt);
fprintf('Tfinal        : %.12g\n', run.Tfinal);
fprintf('nIC           : %d\n', numel(results.sim));
fprintf('Ref file      : %s\n', refFileName);
fprintf('Use spline    : %d\n', run.useRefInterpolant);
fprintf('========================================\n\n');

for k = 1:numel(results.sim)
    sim = results.sim(k);

    fprintf(['IC %2d: %-28s  homErr = %9.3e   orbitDist = %9.3e', ...
             '   alarmTriggered = %d   stopTime = %9.3f\n'], ...
        k, sim.icInfo.label, sim.finalHomErr, sim.finalOrbitDist, ...
        sim.alarmTriggered, sim.stopTime);
end

clearvars -except results post
