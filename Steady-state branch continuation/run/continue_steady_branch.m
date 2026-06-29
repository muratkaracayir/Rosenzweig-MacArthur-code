function PALC = continue_steady_branch()
%CONTINUE_STEADY_BRANCH  Continue a steady-state branch by PALC.
%
%   PALC = continue_steady_branch()
%
% This public driver loads one corrected steady state from the workflow data
% folder and passes it to palc_steady_branch.m.  The loaded point does not
% have to be the branch anchor; it only has to be a corrected steady state
% on the desired branch.
%
% Intended location:
%   Steady-state branch continuation/run/continue_steady_branch.m
%
% Required folder layout:
%   config/
%   continuation/
%   solver/
%   data/

%% Path setup
thisFile = mfilename('fullpath');
runDir   = fileparts(thisFile);
rootDir  = fileparts(runDir);

dataDir = fullfile(rootDir, 'data');

addpath(fullfile(rootDir, 'config'));
addpath(fullfile(rootDir, 'continuation'));
addpath(fullfile(rootDir, 'solver'));

%% User settings
% The default file is the nontrivial steady state produced by
% start_steady_branch.m in the basic mode-5 workflow.
steadyFile = 'SS_Ex3.12_mode5_thetaFixed_lambda0.37725.mat';

% If continuing from a saved PALC point with a known global step number,
% change this value accordingly.  For a fresh start from a corrected steady
% state, use stepIdx0 = 0.
stepIdx0 = 0;

opts = struct();
opts.nSteps = 100;
opts.ds     = -1e-2;

% Optional commonly adjusted settings.  The remaining defaults are set in
% palc_steady_branch.m.
opts.saveEvery     = 10;
opts.verbose       = true;
opts.savePALCAtEnd = true;

%% Load initial corrected steady state
steadyPath = fullfile(dataDir, steadyFile);
if ~exist(steadyPath, 'file')
    error('continue_steady_branch:MissingInputFile', ...
        'Could not find the input steady-state file:\n  %s', steadyPath);
end

[UV0, cfg0, steady] = load_initial_state(steadyPath);

fprintf('Loaded initial steady state:\n');
fprintf('  file   = %s\n', steadyFile);
fprintf('  lambda = %.16g\n', cfg0.lambda);
if isfield(steady, 'mode') && ~isempty(steady.mode)
    fprintf('  mode   = %d\n', steady.mode);
end
fprintf('  nSteps = %d\n', opts.nSteps);
fprintf('  ds     = %.3e\n', opts.ds);

%% Continue branch
PALC = palc_steady_branch(UV0, cfg0, steady, stepIdx0, opts);
end

% =========================================================================
function [UV0, cfg0, steady] = load_initial_state(filePath)
%LOAD_INITIAL_STATE  Extract UV0, cfg0, and steady from a saved SS file.

S = load(filePath);

if isfield(S, 'SS')
    SS = S.SS;
else
    error('continue_steady_branch:MissingSS', ...
        'The input file must contain a top-level variable named SS.');
end

if ~isstruct(SS)
    error('continue_steady_branch:InvalidSS', ...
        'The variable SS must be a structure.');
end

requiredFields = {'UV', 'cfg', 'steady'};
for j = 1:numel(requiredFields)
    fieldName = requiredFields{j};
    if ~isfield(SS, fieldName)
        error('continue_steady_branch:MissingSSField', ...
            'SS.%s must be present in the input file.', fieldName);
    end
end

UV0    = SS.UV;
cfg0   = SS.cfg;
steady = SS.steady;

validate_loaded_state(UV0, cfg0, steady);
end

% =========================================================================
function validate_loaded_state(UV0, cfg0, steady)
%VALIDATE_LOADED_STATE  Basic structural checks before calling PALC.

if ~isstruct(cfg0)
    error('continue_steady_branch:InvalidCfg', ...
        'SS.cfg must be a structure.');
end

requiredCfgFields = {'Nx', 'lambda', 'k', 'd1', 'd2', 'theta'};
for j = 1:numel(requiredCfgFields)
    fieldName = requiredCfgFields{j};
    if ~isfield(cfg0, fieldName)
        error('continue_steady_branch:MissingCfgField', ...
            'SS.cfg.%s must be present.', fieldName);
    end
end

if ~isnumeric(cfg0.Nx) || ~isscalar(cfg0.Nx) || ...
        cfg0.Nx ~= round(cfg0.Nx) || cfg0.Nx < 2
    error('continue_steady_branch:InvalidNx', ...
        'SS.cfg.Nx must be an integer greater than or equal to 2.');
end

if ~isnumeric(UV0) || ndims(UV0) ~= 2 || ...
        size(UV0, 1) ~= cfg0.Nx || size(UV0, 2) ~= 2
    error('continue_steady_branch:InvalidUV', ...
        'SS.UV must be a numeric cfg.Nx-by-2 array.');
end

if any(~isfinite(UV0(:)))
    error('continue_steady_branch:NonfiniteUV', ...
        'SS.UV must contain only finite values.');
end

if ~isstruct(steady)
    error('continue_steady_branch:InvalidSteady', ...
        'SS.steady must be a structure.');
end
end
