function hFig = plot_ode_reference()
%PLOT_ODE_REFERENCE  Quick-look plot of the ODE reference cycle.
%
%   hFig = PLOT_ODE_REFERENCE() loads the ODE reference orbit for Test case 2
%   at lambda = 5 and plots the cycle in the (u,v) phase plane.
%
%   The script expects the file
%
%       data/RefOrbit_Ex2.6.2_lambda_5.mat
%
%   relative to the workflow folder. This file is produced by running
%
%       ode_solver/compute_reference_orbit.m
%
%   The saved .mat file must contain a top-level structure named ref with
%   fields uCycle and vCycle. If tCycle is available, it is used to place a
%   few phase markers along the orbit.
%
%   This file is intentionally lightweight: it is a public quick-look
%   plotter, not the manuscript figure-generation script.
%
%   Output:
%       hFig   Figure handle.

    %% Locate workflow folders
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    rootDir  = fileparts(thisDir);
    dataDir  = fullfile(rootDir, 'data');

    %% User-facing settings
    refFileName  = 'RefOrbit_Ex2.6.2_lambda_5.mat';
    saveFigure   = true;
    figureName   = 'ode_reference_lambda_5.fig';
    phaseMarks   = [0, 0.25, 0.50, 0.75];

    %% Load reference orbit
    refPath = fullfile(dataDir, refFileName);
    S = load(refPath);
    ref = S.ref;

    u = ref.uCycle(:);
    v = ref.vCycle(:);

    if isfield(ref, 'tCycle')
        t = ref.tCycle(:);
        s = (t - t(1)) / (t(end) - t(1));
    elseif isfield(ref, 'phase')
        s = ref.phase(:);
        if max(s) > 1 + 1e-12 || min(s) < -1e-12
            s = (s - s(1)) / (s(end) - s(1));
        end
    else
        s = linspace(0, 1, numel(u)).';
    end

    if isfield(ref, 'lambda')
        lambda = ref.lambda;
    else
        lambda = 5;
    end

    if isfield(ref, 'T')
        periodText = sprintf(', T = %.6g', ref.T);
    else
        periodText = '';
    end

    %% Plot
    hFig = figure('Name', 'ODE reference cycle', 'Color', 'w');
    ax = axes('Parent', hFig);
    hold(ax, 'on');
    box(ax, 'on');

    plot(ax, u, v, '-', 'LineWidth', 1.5);

    % Mark the starting point.
    plot(ax, u(1), v(1), 'o', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.0);

    % Mark a few normalized phases along the orbit.
    for q = 1:numel(phaseMarks)
        [~, idx] = min(abs(s - phaseMarks(q)));
        plot(ax, u(idx), v(idx), 'o', ...
            'MarkerSize', 5, ...
            'MarkerFaceColor', 0.80 * [1 1 1], ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.8);
        text(u(idx), v(idx), sprintf('  s = %.2g', phaseMarks(q)), ...
            'Parent', ax, ...
            'FontSize', 9, ...
            'VerticalAlignment', 'middle');
    end

    xlabel(ax, 'u');
    ylabel(ax, 'v');
    title(ax, sprintf('ODE reference cycle, Test case 2, lambda = %.3g%s', ...
        lambda, periodText), 'Interpreter', 'none');

    axis(ax, 'tight');
    grid(ax, 'on');

    %% Save figure
    if saveFigure
        savefig(hFig, fullfile(thisDir, figureName));
    end
end
