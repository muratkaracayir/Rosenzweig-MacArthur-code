function plot_upo_branch()
%% plot_upo_branch.m
% Quick-look branch diagram for one theta-fixed UPO master file.
%
% Required variable in the selected .mat file:
%     UPOs
%
% Each UPOs(k) entry should contain at least:
%     lambda, T, ampUV

clear; close all; clc;

%% User settings
figuresDir = fileparts(mfilename('fullpath'));
rootDir    = fileparts(figuresDir);
dataDir    = fullfile(rootDir, 'data');

masterFile = 'UPOMaster_Ex2.6.2_mode=2_thetaFixed.mat';

doSave = false;
outFileBase = 'upo_branch';

lineWidth = 1.4;
markerSize = 4;
fontSize = 11;
figPos = [100, 100, 1100, 440];

%% Load branch data
S = load(fullfile(dataDir, masterFile), 'UPOs');

if ~isfield(S, 'UPOs')
    error('plot_upo_branch:MissingUPOs', ...
        'The selected master file must contain a top-level variable named UPOs.');
end

UPOs = S.UPOs;
numUPOs = numel(UPOs);

if numUPOs == 0
    error('plot_upo_branch:EmptyBranch', ...
        'The selected master file contains an empty UPOs array.');
end

lambda = nan(numUPOs, 1);
Tper   = nan(numUPOs, 1);
ampUV  = nan(numUPOs, 1);

for j = 1:numUPOs
    require_field(UPOs(j), 'lambda', j);
    require_field(UPOs(j), 'T', j);
    require_field(UPOs(j), 'ampUV', j);

    lambda(j) = UPOs(j).lambda;
    Tper(j)   = UPOs(j).T;
    ampUV(j)  = UPOs(j).ampUV;
end

modeLabel = infer_mode_label(UPOs(1));

%% Plot
fig = figure('Color', 'w', 'Position', figPos);
movegui(fig, 'center');

ax1 = subplot(1,2,1);
plot(ax1, lambda, ampUV, 'o-', ...
    'LineWidth', lineWidth, ...
    'MarkerSize', markerSize);
grid(ax1, 'on');
box(ax1, 'on');
xlabel(ax1, '\lambda', 'Interpreter', 'tex');
ylabel(ax1, 'Amplitude', 'Interpreter', 'tex');
title(ax1, [modeLabel ': amplitude'], ...
    'Interpreter', 'tex', ...
    'FontWeight', 'normal');
set(ax1, 'FontSize', fontSize);

ax2 = subplot(1,2,2);
plot(ax2, lambda, Tper, 'o-', ...
    'LineWidth', lineWidth, ...
    'MarkerSize', markerSize);
grid(ax2, 'on');
box(ax2, 'on');
xlabel(ax2, '\lambda', 'Interpreter', 'tex');
ylabel(ax2, 'Period T', 'Interpreter', 'tex');
title(ax2, [modeLabel ': period'], ...
    'Interpreter', 'tex', ...
    'FontWeight', 'normal');
set(ax2, 'FontSize', fontSize);

set(fig, 'Renderer', 'painters');

%% Optional save
if doSave
    savefig(fig, fullfile(figuresDir, [outFileBase, '.fig']));
    print(fig, fullfile(figuresDir, [outFileBase, '.png']), '-dpng', '-r300');
end

%% Local helper functions
function require_field(S, fieldName, idx)
    if ~isfield(S, fieldName)
        error('plot_upo_branch:MissingField', ...
            'UPOs(%d) must contain the field %s.', idx, fieldName);
    end
end

function modeLabel = infer_mode_label(upo)
    modeLabel = 'Theta-fixed branch';

    if isfield(upo, 'hopf') && isfield(upo.hopf, 'mode')
        modeLabel = sprintf('Theta-fixed mode %d branch', upo.hopf.mode);
    elseif isfield(upo, 'mode')
        modeLabel = sprintf('Theta-fixed mode %d branch', upo.mode);
    end
end
end