function plot_constant_history()
%PLOT_CONSTANT_HISTORY  Quick plot of one saved time-stepping run.
%
%   This script is a lightweight visualization companion for the public
%   constant-equilibrium time-stepping example. It loads one saved .mat file
%   from the repository data/ folder and plots the spatial minimum and
%   maximum histories of u and v for one representative initial condition.
%
%   Intended repository layout:
%
%       Constant equilibria/
%         data/
%           Constant_Ex2.6.1_lambda_10_Nx128_dt0.0001_k3.mat
%         figures/
%           plot_constant_history_quicklook.m
%
%   Before running this script, run
%
%       run/run_constant_equilibrium.m
%
%   with run.saveResults = true and run.onlyICMode = ''.
%
%   MATLAB compatibility: written in a style compatible with MATLAB R2016a.

    % ------------------------------------------------------------
    % User settings
    % ------------------------------------------------------------
    dataFileName = 'Constant_Ex2.6.1_lambda_10_Nx128_dt0.01_k3.mat';

    % In the default base suite, simIndex = 2 is the smooth-random initial
    % condition. Change this index to inspect another stored run.
    simIndex = 2;

    saveFigure = true;

    % ------------------------------------------------------------
    % Locate folders relative to this file
    % ------------------------------------------------------------
    thisFile   = mfilename('fullpath');
    figuresDir = fileparts(thisFile);
    projectDir = fileparts(figuresDir);
    dataDir    = fullfile(projectDir, 'data');

    % ------------------------------------------------------------
    % Load saved run
    % ------------------------------------------------------------
    S = load(fullfile(dataDir, dataFileName), 'cfg', 'results', 'summary');

    cfg     = S.cfg;
    results = S.results;
    H       = results.sim(simIndex).history;

    % Prefer the saved summary target. This keeps the plot consistent with
    % the equilibrium convention used by run_strang.
    if isfield(S, 'summary') && isfield(S.summary, 'targetU') && isfield(S.summary, 'targetV')
        targetU = S.summary.targetU;
        targetV = S.summary.targetV;
    else
        if cfg.lambda < cfg.k
            targetU = cfg.u_star;
            targetV = cfg.v_star;
        else
            targetU = cfg.k;
            targetV = 0;
        end
    end

    % ------------------------------------------------------------
    % One-panel quick-look figure
    % ------------------------------------------------------------
    fig = figure('Color', 'w', ...
                 'Name', 'Constant-equilibrium quick look', ...
                 'NumberTitle', 'off');

    [ax, ~, ~] = plotyy(H.t, H.umin, H.t, H.vmin); %#ok<PLOTYY>
    hold(ax(1), 'on');
    hold(ax(2), 'on');

    % Remove the two dummy lines created by plotyy.
    delete(get(ax(1), 'Children'));
    delete(get(ax(2), 'Children'));

    uColor = [0.00, 0.25, 0.75];
    vColor = [0.80, 0.30, 0.00];

    % u-envelope on left axis.
    plot(ax(1), H.t, H.umin, '-',  'LineWidth', 1.2, 'Color', uColor);
    plot(ax(1), H.t, H.umax, '--', 'LineWidth', 1.2, 'Color', uColor);

    % v-envelope on right axis.
    plot(ax(2), H.t, H.vmin, '-',  'LineWidth', 1.2, 'Color', vColor);
    plot(ax(2), H.t, H.vmax, '--', 'LineWidth', 1.2, 'Color', vColor);

    % Equilibrium levels.
    xRange = [H.t(1), H.t(end)];
    plot(ax(1), xRange, [targetU, targetU], ':', ...
         'LineWidth', 1.5, 'Color', uColor);
    plot(ax(2), xRange, [targetV, targetV], ':', ...
         'LineWidth', 1.5, 'Color', vColor);

    % Labels and title.
    xlabel(ax(1), 't');
    ylabel(ax(1), 'u envelope');
    ylabel(ax(2), 'v envelope');

    title(ax(1), sprintf('%s, \\lambda = %.6g, IC: %s', ...
          cfg.caseName, cfg.lambda, results.sim(simIndex).icInfo.label));

    set(ax(1), 'YColor', uColor, 'Box', 'on');
    set(ax(2), 'YColor', vColor, 'Box', 'off');
    grid(ax(1), 'on');

    xlim(ax(1), xRange);
    xlim(ax(2), xRange);

    legend(ax(1), {'u_{min}', 'u_{max}', 'u_*'}, ...
           'Location', 'best');

    % Add a small text annotation for the v-equilibrium, since plotyy cannot
    % combine the left- and right-axis objects into one simple legend.
    text('Parent', ax(2), ...
         'Units', 'normalized', ...
         'Position', [0.98, 0.08, 0], ...
         'String', 'v_{min}, v_{max}, v_* on right axis', ...
         'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'bottom', ...
         'Color', vColor, ...
         'FontSize', 9);

    % ------------------------------------------------------------
    % Optional save
    % ------------------------------------------------------------
    if saveFigure
        outFig = fullfile(figuresDir, 'constant_history_quicklook.fig');
        savefig(fig, outFig);
        fprintf('Saved figure to:\n  %s\n', outFig);
    end
end
