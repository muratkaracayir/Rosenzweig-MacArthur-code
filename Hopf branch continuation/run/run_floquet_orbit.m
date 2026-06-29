%% run_floquet_orbit.m
% Run Floquet analysis for one saved periodic orbit.
%
% This public driver loads one UPO file from the workflow data/ folder,
% computes the Floquet multipliers with the IFRK4 explicit-monodromy
% workflow, prints a compact stability summary, and optionally saves the
% Floquet output in data/.
%
% Place this file in:
%   Hopf branch continuation/run/run_floquet_orbit.m

%% Path setup
thisFile = mfilename('fullpath');
runDir   = fileparts(thisFile);
rootDir  = fileparts(runDir);

addpath(fullfile(rootDir, 'floquet'));

dataDir = fullfile(rootDir, 'data');
if exist(dataDir, 'dir') ~= 7
    error('run_floquet_orbit:MissingDataDir', ...
        'The data folder was not found: %s', dataDir);
end

%% User-editable choices
% Saved periodic orbit to analyze. The file is searched for in data/ unless
% an absolute or relative path to an existing file is supplied.
upoFile = 'UPO_Ex2.6.2_mode2_thetaFixed_2.mat';

% Save the Floquet result as Floquet_*.mat in data/.
doSave = true;

% Optional Floquet settings. Increase mSub for a more resolved monodromy
% computation, or increase nStore to keep more leading multipliers.
floqOverrides = struct( ...
    'mSub', 8, ...
    'nStore', 16, ...
    'verbose', false);

%% Locate and load the saved orbit
if exist(upoFile, 'file') == 2
    upoPath = upoFile;
else
    upoPath = fullfile(dataDir, upoFile);
end

if exist(upoPath, 'file') ~= 2
    error('run_floquet_orbit:MissingUPOFile', ...
        'Could not find the saved orbit file: %s', upoPath);
end

S = load(upoPath);
requiredVars = {'UV', 'T', 'hopf'};
for j = 1:numel(requiredVars)
    name = requiredVars{j};
    if ~isfield(S, name)
        error('run_floquet_orbit:MissingVariable', ...
            'The saved orbit file must contain variable "%s".', name);
    end
end

UV = S.UV;
T  = S.T;
hopf = S.hopf;

if ~isfield(hopf, 'cfg')
    error('run_floquet_orbit:MissingHopfCfg', ...
        'The loaded hopf struct must contain hopf.cfg.');
end
cfgOrbit = hopf.cfg;

%% Run Floquet analysis
cfgFloq = make_floqCfg(floqOverrides);

[~, upoBase, ~] = fileparts(upoPath);
runInfo = struct();
runInfo.sourceFile = upoPath;
runInfo.sourceOrbit = upoBase;

floquet = compute_floquet_IFRK4(UV, T, cfgOrbit, cfgFloq, runInfo);

%% Compact stability summary
muStored = floquet.eig.mu(:);
posTrivial = floquet.eig.posTrivial;
nontrivialMask = true(size(muStored));
if ~isempty(posTrivial) && posTrivial >= 1 && posTrivial <= numel(muStored)
    nontrivialMask(posTrivial) = false;
end

muNontrivial = muStored(nontrivialMask);
if isempty(muNontrivial)
    leadingMu = NaN;
else
    [~, idxLead] = max(abs(muNontrivial));
    leadingMu = muNontrivial(idxLead);
end

if floquet.summary.nUnstable == 0
    verdict = 'Floquet-stable';
else
    verdict = 'Floquet-unstable';
end

fprintf('\nFloquet analysis completed.\n');
fprintf('Source orbit: %s\n', upoBase);
fprintf('Leading nontrivial multiplier: %.16g%+.16gi  |mu| = %.16g\n', ...
    real(leadingMu), imag(leadingMu), abs(leadingMu));
fprintf('Number of unstable directions: %d\n', floquet.summary.nUnstable);
fprintf('Verdict: %s\n', verdict);

if ~floquet.qc.ok
    fprintf('QC warning(s):\n');
    for j = 1:numel(floquet.qc.reason)
        fprintf('  - %s\n', floquet.qc.reason{j});
    end
end

%% Optional save
if doSave
    if length(upoBase) >= 4 && strcmp(upoBase(1:4), 'UPO_')
        floquetBase = ['Floquet_' upoBase(5:end)];
    else
        floquetBase = ['Floquet_' upoBase];
    end
    floquetFile = [floquetBase '.mat'];
    floquetPath = fullfile(dataDir, floquetFile);

    sourceFile = upoPath; %#ok<NASGU>
    save(floquetPath, 'floquet', 'cfgFloq', 'sourceFile', 'upoFile');
    fprintf('Saved Floquet result: %s\n', floquetPath);
end
