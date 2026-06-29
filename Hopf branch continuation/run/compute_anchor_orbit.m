%% compute_anchor_orbit.m
% Compute one Newton-corrected periodic orbit from a Hopf initialization.
%
% This is the first entry point for the Hopf-branch continuation workflow.
% It prepares a small-amplitude mode-n Hopf guess and corrects it by Newton's
% method. If doSave is true, newton_solver saves the corrected orbit in data/.

%% Path setup
thisFile = mfilename('fullpath');
runDir   = fileparts(thisFile);
rootDir  = fileparts(runDir);

addpath(fullfile(rootDir, 'config'));
addpath(fullfile(rootDir, 'continuation'));
addpath(fullfile(rootDir, 'solver'));

%% User settings
lambdaHopf = 0.5;   % mode-2 Hopf point for Example 2.6.2
modeIndex  = 2;
doSave     = true;

%% Prepare Hopf initial guess
cfg = case_ex262(lambdaHopf);
[init, hopf] = prepare_data(cfg, modeIndex);

%% Newton correction
UV_init = [init.U0(:); init.V0(:)];
[UV, T, resNorm, ampUV] = newton_solver(UV_init, init.T_H, hopf, doSave);

%% Diagnostics
fprintf('\n=== Anchor orbit computation complete ===\n');
fprintf('Case:              %s\n', hopf.cfg.caseName);
fprintf('Mode:              %d\n', hopf.mode);
if isfield(hopf.cfg, 'route') && isfield(hopf.cfg.route, 'type')
    routeTag = hopf.cfg.route.type;
elseif isfield(hopf.cfg, 'route') && isfield(hopf.cfg.route, 'name')
    routeTag = hopf.cfg.route.name;
else
    routeTag = 'unknown';
end
fprintf('Route:             %s\n', routeTag);
fprintf('Corrected lambda:  %.12g\n', hopf.cfg.lambda);
fprintf('Corrected period:  %.12g\n', T);
fprintf('Residual norm:     %.3e\n', resNorm);
fprintf('Amplitude:         %.6e\n', ampUV);
fprintf('=========================================\n\n');
