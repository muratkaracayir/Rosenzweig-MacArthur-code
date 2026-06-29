%% run_constant_equilibrium.m
%RUN_CONSTANT_EQUILIBRIUM  Time stepping toward a constant equilibrium.
%
%   This script is the public entry point for the constant-equilibrium
%   time-stepping workflow in the repository.
%
%   Intended repository layout:
%
%       Constant equilibria/
%         config/
%           case_ex261.m
%         run/
%           run_constant_equilibrium.m
%         solver/
%           run_strang.m
%           ...
%
%   The script adds the required config/ and solver/ folders to the MATLAB
%   path relative to its own location, then runs a Strang-splitting time
%   integration for the selected parameter value.
%
%   The numerical settings below are intentionally modest and are meant as a
%   representative public example. Larger final times, smaller time steps, or
%   more initial data may be used for production computations.
%
%   MATLAB compatibility: written in a style compatible with MATLAB R2016a.

clear; clc;

%% ------------------------------------------------------------
% 0) Add required repository folders to the MATLAB path
%% ------------------------------------------------------------
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
rootDir  = fileparts(thisDir);

addpath(fullfile(rootDir, 'config'));
addpath(fullfile(rootDir, 'solver'));

%% ------------------------------------------------------------
% 1) Select test case and bifurcation parameter
%% ------------------------------------------------------------
% For case_ex261, all model parameters are positive. The value lambda is
% assumed positive. If 0 < lambda < k, the relevant nonnegative homogeneous
% equilibrium is the coexistence equilibrium (lambda, v_lambda). If
% lambda >= k, the relevant nonnegative homogeneous equilibrium is the
% prey-only state (k,0).
lambda = 12;

% Choose the active case here.
cfg = case_ex261(lambda);

%% ------------------------------------------------------------
% 2) Numerical / implementation parameters
%% ------------------------------------------------------------
run.Nx         = 128;
run.dt         = 1e-4;
run.Tfinal     = 200;

run.saveEvery   = 1000;
run.outputEvery = 1000;
run.probeStep   = 1;

run.exploreIC     = false;   % false -> base 3 ICs, true -> extended 11 ICs
run.onlyICMode    = '';
run.onlyICSample  = [];
run.seed          = 1;
run.verbose       = true;
run.saveResults   = true;

run.useStopping   = true;
run.stopEvery     = run.outputEvery;
run.eqTol         = 1e-8;
run.stepTol       = 1e-6;
run.minStopTime   = 1;
run.stopEqMode    = 'relative';   % 'absolute' or 'relative'
run.eqScaleFloor  = 1;            % used only when stopEqMode = 'relative'

%% ------------------------------------------------------------
% 3) Run simulation suite
%% ------------------------------------------------------------
results = run_strang(cfg, run);

%% ------------------------------------------------------------
% 4) Basic quick-look summary in the command window
%% ------------------------------------------------------------
fprintf('\n');
fprintf('========================================\n');
fprintf('Completed run_constant_equilibrium\n');
fprintf('Case      : %s\n', cfg.caseName);
fprintf('lambda    : %.12g\n', cfg.lambda);
fprintf('Nx        : %d\n', run.Nx);
fprintf('dt        : %.12g\n', run.dt);
fprintf('Tfinal    : %.12g\n', run.Tfinal);
fprintf('nIC       : %d\n', numel(results.sim));
fprintf('========================================\n\n');

% Determine the target constant equilibrium used in the final-distance
% diagnostic. At lambda = k, both formulas give the prey-only state.
if cfg.lambda < cfg.k
    targetU = cfg.u_star;
    targetV = cfg.v_star;
else
    targetU = cfg.k;
    targetV = 0;
end

for k = 1:numel(results.sim)
    Ufinal = results.sim(k).Ufinal;

    eqDistU = max(abs(Ufinal(1,:) - targetU));
    eqDistV = max(abs(Ufinal(2,:) - targetV));

    fprintf(['IC %2d: %-28s  eqDistU = %9.3e   eqDistV = %9.3e', ...
             '   stoppedEarly = %d   stopTime = %9.3f\n'], ...
        k, results.sim(k).icInfo.label, eqDistU, eqDistV, ...
        results.sim(k).stoppedEarly, results.sim(k).stopTime);
end

clearvars -except results
