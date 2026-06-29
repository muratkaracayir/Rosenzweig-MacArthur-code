%% run_continue_branch.m
% Continue a Hopf branch from a saved anchor periodic orbit.
%
% This is the public-facing driver for branch continuation. Edit the user
% settings below, then run this script from MATLAB. The continuation engine is
% continuation/continue_branch.m.

%% Path setup
thisFile = mfilename('fullpath');
runDir   = fileparts(thisFile);
rootDir  = fileparts(runDir);

dataDir = fullfile(rootDir, 'data');

addpath(fullfile(rootDir, 'config'));
addpath(fullfile(rootDir, 'continuation'));
addpath(fullfile(rootDir, 'solver'));

%% User settings
caseFun      = @case_ex262;
hopfLambda   = 0.5;      % Hopf point used to identify the test case
modeIndex    = 2;        % spatial mode
anchorLambda = 0.5001;   % lambda value of the saved anchor UPO

nSteps      = 50;        % number of continuation steps
dLambda     = 0.01;      % lambda increment per continuation step
keepInitial = true;      % include the anchor orbit in the saved branch

%% Build the expected anchor filename
cfg = caseFun(hopfLambda);
caseName = cfg.caseName;

if isfield(cfg, 'route') && isfield(cfg.route, 'type')
    routeTag = cfg.route.type;
elseif isfield(cfg, 'route') && isfield(cfg.route, 'name')
    routeTag = cfg.route.name;
else
    routeTag = 'thetaFixed';
end

expectedName = sprintf('UPO_%s_mode%d_%s_lambda%g.mat', ...
    caseName, modeIndex, routeTag, anchorLambda);
expectedFile = fullfile(dataDir, expectedName);

%% Resolve the anchor file
anchorFile = '';

if exist(expectedFile, 'file') == 2
    anchorFile = expectedFile;
else
    searchPattern = sprintf('UPO_%s_mode%d_%s_lambda*.mat', ...
        caseName, modeIndex, routeTag);
    candidates = dir(fullfile(dataDir, searchPattern));

    lambdaTol = 1e-10;
    bestMismatch = Inf;
    bestFile = '';

    for j = 1:numel(candidates)
        candidateFile = fullfile(dataDir, candidates(j).name);
        S = load(candidateFile, 'lambda');

        if isfield(S, 'lambda')
            mismatch = abs(S.lambda - anchorLambda);
            if mismatch < bestMismatch
                bestMismatch = mismatch;
                bestFile = candidateFile;
            end

            if mismatch <= lambdaTol
                anchorFile = candidateFile;
                break;
            end
        end
    end

    if isempty(anchorFile)
        fprintf('\nExpected anchor file was not found:\n  %s\n', expectedFile);
        if ~isempty(bestFile)
            fprintf('Closest matching file was:\n  %s\n', bestFile);
            fprintf('with lambda mismatch %.3e.\n', bestMismatch);
        end
        error('run_continue_branch:AnchorNotFound', ...
            'No saved UPO matching the requested case, mode, route, and anchor lambda was found.');
    end
end

%% Run continuation
fprintf('\n=== Hopf branch continuation ===\n');
fprintf('Case:          %s\n', caseName);
fprintf('Mode:          %d\n', modeIndex);
fprintf('Route:         %s\n', routeTag);
fprintf('Anchor file:   %s\n', anchorFile);
fprintf('Anchor lambda: %.12g\n', anchorLambda);
fprintf('Steps:         %d\n', nSteps);
fprintf('dLambda:       %.12g\n', dLambda);
fprintf('Keep initial:  %d\n', keepInitial);
fprintf('================================\n\n');

continue_branch(anchorFile, nSteps, dLambda, keepInitial);
