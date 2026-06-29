function plot_ss_profiles()
%PLOT_SS_PROFILES  Two-panel steady-state profile plotter.
%
%   plot_ss_profiles()
%
% The function scans the workflow data folder for SS_*.mat files.  Each
% file is expected to contain a structure named SS with at least the fields
% lambda, UV, and cfg.
%
% Simple mode plots the single steady state whose lambda value is closest
% to targetLambda.
%
% Advanced mode reconstructs the continuation order using SS.stepIdx when
% available, splits the ordered states into sheets by detecting turning
% points in lambda, and plots the closest state to targetLambda on each
% sheet, provided that it lies within lambdaTol.
%
% Intended location:
%   Steady-state branch continuation/figures/plot_ss_profiles.m

clear; close all; clc;

%% User settings
targetLambda = 0.375;
advancedMode = true;

lambdaTol = 1e-3;
sheetTurnTol = 1e-10;

figPos     = [100, 100, 1100, 430];
fontSize   = 11;
lineWidth  = 1.6;

%% Locate folders
figuresDir = fileparts(mfilename('fullpath'));
rootDir    = fileparts(figuresDir);
dataDir    = fullfile(rootDir, 'data');

if ~exist(dataDir, 'dir')
    error('plot_ss_profiles:MissingDataFolder', ...
        'The data folder does not exist: %s', dataDir);
end

%% Load steady-state files
records = load_ss_records(dataDir);

if isempty(records)
    error('plot_ss_profiles:NoSteadyStates', ...
        'No valid SS_*.mat files were found in %s.', dataDir);
end

fprintf('Loaded %d steady-state files from:\n  %s\n', ...
    numel(records), dataDir);
fprintf('Target lambda: %.12g\n', targetLambda);

%% Select profiles
if advancedMode
    records = order_records_for_sheets(records);
    sheetIdx = split_into_sheets([records.lambda], sheetTurnTol);
    selected = select_one_profile_per_sheet(records, sheetIdx, ...
        targetLambda, lambdaTol);
else
    selected = select_single_profile(records, targetLambda, lambdaTol);
end

if isempty(selected)
    error('plot_ss_profiles:NoProfilesInTolerance', ...
        ['No selected profile satisfies |lambda - targetLambda| <= lambdaTol. ', ...
         'Try increasing lambdaTol or choosing a different targetLambda.']);
end

print_selection_summary(selected, targetLambda, lambdaTol);

%% Plot selected profiles
fig = figure('Color', 'w', 'Position', figPos);
movegui(fig, 'center');

axU = subplot(1, 2, 1); hold(axU, 'on'); grid(axU, 'on'); box(axU, 'on');
axV = subplot(1, 2, 2); hold(axV, 'on'); grid(axV, 'on'); box(axV, 'on');

C = lines(max(numel(selected), 1));
lineStyles = {'-', '--', '-.', ':'};
legendText = cell(numel(selected), 1);
legendHandles = zeros(numel(selected), 1);

for j = 1:numel(selected)
    rec = selected(j);
    [x, U, V] = unpack_profile(rec.SS);

    style = lineStyles{1 + mod(j - 1, numel(lineStyles))};

    hU = plot(axU, x, U, ...
        'LineWidth', lineWidth, ...
        'LineStyle', style, ...
        'Color', C(j, :));

    plot(axV, x, V, ...
        'LineWidth', lineWidth, ...
        'LineStyle', style, ...
        'Color', C(j, :));

    legendHandles(j) = hU;
    legendText{j} = make_legend_entry(rec, targetLambda);
end

plot_equilibrium_lines(axU, axV, selected(1).SS.cfg);

xlabel(axU, 'x', 'Interpreter', 'tex');
ylabel(axU, 'u(x)', 'Interpreter', 'tex');
title(axU, 'Prey profile', ...
    'Interpreter', 'tex', ...
    'FontWeight', 'normal');

xlabel(axV, 'x', 'Interpreter', 'tex');
ylabel(axV, 'v(x)', 'Interpreter', 'tex');
title(axV, 'Predator profile', ...
    'Interpreter', 'tex', ...
    'FontWeight', 'normal');

set([axU axV], 'FontSize', fontSize);

lg = legend(axU, legendHandles, legendText, ...
    'Location', 'best', ...
    'Interpreter', 'tex');
set(lg, 'Box', 'off', 'FontSize', fontSize - 1);

end

% ---------------------------------------------------------------------- %
function records = load_ss_records(dataDir)
%LOAD_SS_RECORDS  Load valid SS structures from SS_*.mat files.

files = dir(fullfile(dataDir, 'SS_*.mat'));
records = empty_record_array();

for k = 1:numel(files)
    filePath = fullfile(dataDir, files(k).name);
    S = load(filePath);

    if ~isfield(S, 'SS') || ~isstruct(S.SS)
        warning('plot_ss_profiles:SkippingFile', ...
            'Skipping %s because it does not contain a structure SS.', ...
            files(k).name);
        continue;
    end

    SS = S.SS;
    if ~is_valid_ss(SS)
        warning('plot_ss_profiles:SkippingFile', ...
            'Skipping %s because SS is missing required fields.', ...
            files(k).name);
        continue;
    end

    rec = struct();
    rec.fileName = files(k).name;
    rec.filePath = filePath;
    rec.fileDate = files(k).datenum;
    rec.SS       = SS;
    rec.lambda   = double(SS.lambda);
    rec.stepIdx  = get_step_index(SS);
    rec.sheetNo  = NaN;

    records(end + 1, 1) = rec; %#ok<AGROW>
end
end

% ---------------------------------------------------------------------- %
function tf = is_valid_ss(SS)
%IS_VALID_SS  Basic structural check for a saved steady-state record.

tf = isfield(SS, 'lambda') && isnumeric(SS.lambda) && ...
        isscalar(SS.lambda) && isfinite(SS.lambda) && ...
     isfield(SS, 'UV') && isnumeric(SS.UV) && ndims(SS.UV) == 2 && ...
        size(SS.UV, 2) == 2 && ...
     isfield(SS, 'cfg') && isstruct(SS.cfg);
end

% ---------------------------------------------------------------------- %
function stepIdx = get_step_index(SS)
%GET_STEP_INDEX  Return SS.stepIdx when available.

if isfield(SS, 'stepIdx') && isnumeric(SS.stepIdx) && ...
        isscalar(SS.stepIdx) && isfinite(SS.stepIdx)
    stepIdx = double(SS.stepIdx);
else
    stepIdx = NaN;
end
end

% ---------------------------------------------------------------------- %
function records = order_records_for_sheets(records)
%ORDER_RECORDS_FOR_SHEETS  Order records before sheet detection.

stepIdx = [records.stepIdx];

if all(isfinite(stepIdx)) && numel(unique(stepIdx)) == numel(stepIdx)
    [~, ord] = sort(stepIdx);
    records = records(ord);
    fprintf('Advanced mode: records ordered by SS.stepIdx.\n');
else
    warning('plot_ss_profiles:MissingStepIdx', ...
        ['Advanced mode is using file timestamps because not all SS files ', ...
         'contain unique finite SS.stepIdx values. Sheet detection may be ', ...
         'unreliable for files not saved in continuation order.']);
    [~, ord] = sort([records.fileDate]);
    records = records(ord);
end
end

% ---------------------------------------------------------------------- %
function sheetIdx = split_into_sheets(lambda, tol)
%SPLIT_INTO_SHEETS  Split an ordered lambda sequence at turning points.

n = numel(lambda);
if n == 0
    sheetIdx = cell(0, 1);
    return;
end

if n == 1
    sheetIdx = {1};
    return;
end

dl = diff(lambda(:));
signDL = zeros(size(dl));
signDL(dl >  tol) =  1;
signDL(dl < -tol) = -1;

breaks = [];
lastSign = 0;

for j = 1:numel(signDL)
    thisSign = signDL(j);

    if thisSign == 0
        continue;
    end

    if lastSign ~= 0 && thisSign ~= lastSign
        breaks(end + 1) = j + 1; %#ok<AGROW>
    end

    lastSign = thisSign;
end

starts = [1, breaks];
stops  = [breaks - 1, n];

sheetIdx = cell(numel(starts), 1);
for j = 1:numel(starts)
    sheetIdx{j} = starts(j):stops(j);
end

fprintf('Advanced mode: detected %d sheet(s).\n', numel(sheetIdx));
end

% ---------------------------------------------------------------------- %
function selected = select_single_profile(records, targetLambda, lambdaTol)
%SELECT_SINGLE_PROFILE  Select the globally closest profile within tolerance.

[~, idx] = min(abs([records.lambda] - targetLambda));
rec = records(idx);
rec.sheetNo = 1;

if abs(rec.lambda - targetLambda) <= lambdaTol
    selected = rec;
else
    selected = empty_record_array();
end
end

% ---------------------------------------------------------------------- %
function selected = select_one_profile_per_sheet(records, sheetIdx, ...
    targetLambda, lambdaTol)
%SELECT_ONE_PROFILE_PER_SHEET  Select closest in-tolerance profile per sheet.

selected = empty_record_array();

for s = 1:numel(sheetIdx)
    idxSheet = sheetIdx{s};
    lambdaSheet = [records(idxSheet).lambda];
    [~, loc] = min(abs(lambdaSheet - targetLambda));

    rec = records(idxSheet(loc));
    rec.sheetNo = s;

    if abs(rec.lambda - targetLambda) <= lambdaTol
        selected(end + 1, 1) = rec; %#ok<AGROW>
    else
        fprintf(['  sheet %d skipped: closest lambda = %.12g, ', ...
            '|dlam| = %.3e > lambdaTol = %.3e\n'], ...
            s, rec.lambda, abs(rec.lambda - targetLambda), lambdaTol);
    end
end
end

% ---------------------------------------------------------------------- %
function [x, U, V] = unpack_profile(SS)
%UNPACK_PROFILE  Extract x, u, and v from a saved SS structure.

UV = SS.UV;
U  = UV(:, 1);
V  = UV(:, 2);
Nx = size(UV, 1);

cfg = SS.cfg;

if isfield(cfg, 'x') && isnumeric(cfg.x) && numel(cfg.x) == Nx
    x = cfg.x(:);
elseif isfield(cfg, 'L') && isnumeric(cfg.L) && isscalar(cfg.L)
    x = linspace(0, cfg.L, Nx).';
else
    x = (1:Nx).';
end
end

% ---------------------------------------------------------------------- %
function plot_equilibrium_lines(axU, axV, cfg)
%PLOT_EQUILIBRIUM_LINES  Add homogeneous equilibrium reference levels.

if isfield(cfg, 'u_star') && isnumeric(cfg.u_star) && isscalar(cfg.u_star)
    xl = get(axU, 'XLim');
    plot(axU, xl, [cfg.u_star, cfg.u_star], 'k:', 'LineWidth', 0.8);
    set(axU, 'XLim', xl);
end

if isfield(cfg, 'v_star') && isnumeric(cfg.v_star) && isscalar(cfg.v_star)
    xl = get(axV, 'XLim');
    plot(axV, xl, [cfg.v_star, cfg.v_star], 'k:', 'LineWidth', 0.8);
    set(axV, 'XLim', xl);
end
end

% ---------------------------------------------------------------------- %
function txt = make_legend_entry(rec, targetLambda)
%MAKE_LEGEND_ENTRY  Legend text for a selected profile.

if isfinite(rec.sheetNo) && rec.sheetNo > 0
    prefix = sprintf('sheet %d', rec.sheetNo);
else
    prefix = 'profile';
end

delta = abs(rec.lambda - targetLambda);
txt = sprintf('%s: \\lambda=%.6g, |\\Delta\\lambda|=%.2e', ...
    prefix, rec.lambda, delta);
end

% ---------------------------------------------------------------------- %
function print_selection_summary(selected, targetLambda, lambdaTol)
%PRINT_SELECTION_SUMMARY  Print selected profiles and lambda distances.

fprintf('\nSelected profile(s):\n');
for j = 1:numel(selected)
    rec = selected(j);
    dlam = abs(rec.lambda - targetLambda);

    if isfinite(rec.stepIdx)
        stepText = sprintf('%g', rec.stepIdx);
    else
        stepText = 'missing';
    end

    fprintf('  sheet %d: lambda = %.12g, |dlam| = %.3e, stepIdx = %s, file = %s\n', ...
        rec.sheetNo, rec.lambda, dlam, stepText, rec.fileName);

    if dlam > lambdaTol
        error('plot_ss_profiles:InternalToleranceError', ...
            'Selected profile violates lambdaTol: |dlam| = %.3e > %.3e.', ...
            dlam, lambdaTol);
    end
end
end

% ---------------------------------------------------------------------- %
function records = empty_record_array()
%EMPTY_RECORD_ARRAY  Empty structure array used for loaded records.

records = struct( ...
    'fileName', {}, ...
    'filePath', {}, ...
    'fileDate', {}, ...
    'SS', {}, ...
    'lambda', {}, ...
    'stepIdx', {}, ...
    'sheetNo', {});
end
