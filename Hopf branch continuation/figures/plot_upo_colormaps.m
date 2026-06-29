%% plot_upo_colormaps.m
% Quick-look two-panel colormap plot for one stored UPO from Test case 2.
%
% Required variables in the selected .mat file:
%     UV, T, hopf
%
% Optional variable:
%     lambda

clear; close all; clc;

%% User settings
figuresDir = fileparts(mfilename('fullpath'));
rootDir    = fileparts(figuresDir);
dataDir    = fullfile(rootDir, 'data');

upoFile = 'UPO_Ex2.6.2_mode2_thetaFixed_lambda2.mat';

doSave = false;
outFileBase = 'upo_colormaps';

figPos = [100, 100, 1100, 460];
fontSizeAxes = 11;

%% Load orbit
S = load(fullfile(dataDir, upoFile), 'UV', 'T', 'hopf', 'lambda');

UV   = S.UV;
T    = S.T;
hopf = S.hopf;

if isfield(S, 'lambda')
    lambda = S.lambda;
else
    lambda = hopf.cfg.lambda;
end

Nx = hopf.cfg.periodic_orbit.Nx;
Nt = hopf.cfg.periodic_orbit.Nt;
L  = hopf.cfg.L;

Ntot = Nx * Nt;
U = reshape(UV(1:Ntot), Nx, Nt);
V = reshape(UV(Ntot+1:end), Nx, Nt);

x = linspace(0, L, Nx);
t = linspace(0, T, Nt + 1);
t(end) = [];

modeIndex = hopf.mode;

%% Plot
fig = figure('Color', 'w', 'Position', figPos);
movegui(fig, 'center');

subplot(1,2,1);
imagesc(x, t, U.');
axis xy;
colormap(gca, parula);
colorbar;
xlabel('x', 'Interpreter', 'tex');
ylabel('t', 'Interpreter', 'tex');
title(sprintf('(a) Test case 2, mode %d, \\lambda = %.12g: u(x,t)', ...
    modeIndex, lambda), 'Interpreter', 'tex', 'FontWeight', 'normal');
set(gca, 'FontSize', fontSizeAxes);

subplot(1,2,2);
imagesc(x, t, V.');
axis xy;
colormap(gca, parula);
colorbar;
xlabel('x', 'Interpreter', 'tex');
ylabel('t', 'Interpreter', 'tex');
title(sprintf('(b) Test case 2, mode %d, \\lambda = %.12g: v(x,t)', ...
    modeIndex, lambda), 'Interpreter', 'tex', 'FontWeight', 'normal');
set(gca, 'FontSize', fontSizeAxes);

%% Optional save
if doSave
    savefig(fig, fullfile(dataDir, [outFileBase, '.fig']));
    print(fig, fullfile(dataDir, [outFileBase, '.png']), '-dpng', '-r300');
end
