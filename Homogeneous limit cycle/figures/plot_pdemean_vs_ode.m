function hFig = plot_pdemean_vs_ode(simIndex)
%PLOT_PDEMEAN_VS_ODE  Quick-look plot of PDE mean tail and ODE reference.
%
%   hFig = PLOT_PDEMEAN_VS_ODE() loads the Test case 2 data at lambda = 5
%   and plots the late-time PDE spatial-mean trajectory together with the
%   ODE reference cycle in the (u,v) phase plane.
%
%   hFig = PLOT_PDEMEAN_VS_ODE(simIndex) uses the requested stored PDE
%   simulation. The default is simIndex = 2.
%
%   The script expects the files
%
%       data/RefOrbit_Ex2.6.2_lambda_5.mat
%       data/PDEOrbit_Ex2.6.2_lambda_5_Nx128_dt0.01.mat
%
%   relative to the workflow folder. The ODE file is produced by
%
%       ode_solver/compute_reference_orbit.m
%
%   and the PDE file is produced by
%
%       run/run_homogeneous_periodic.m
%
%   The ODE reference orbit is plotted with a dashed line, and the PDE
%   spatial-mean tail is plotted with a solid line. This file is
%   intentionally lightweight: it is a public quick-look plotter, not the
%   manuscript figure-generation script.
%
%   This function does not save any files.

    if nargin < 1 || isempty(simIndex)
        simIndex = 2;
    end

    %% ------------------------------------------------------------
    % Locate workflow folders
    %% ------------------------------------------------------------
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    rootDir  = fileparts(thisDir);
    dataDir  = fullfile(rootDir, 'data');

    %% ------------------------------------------------------------
    % Default files
    %% ------------------------------------------------------------
    refFileName = 'RefOrbit_Ex2.6.2_lambda_5.mat';
    pdeFileName = 'PDEOrbit_Ex2.6.2_lambda_5_Nx128_dt0.001.mat';

    refPath = fullfile(dataDir, refFileName);
    pdePath = fullfile(dataDir, pdeFileName);

    if exist(refPath, 'file') ~= 2
        error('Reference-orbit file not found: %s', refPath);
    end
    if exist(pdePath, 'file') ~= 2
        error('PDE-output file not found: %s', pdePath);
    end

    %% ------------------------------------------------------------
    % Load files
    %% ------------------------------------------------------------
    Sref = load(refPath);
    Spde = load(pdePath);

    if ~isfield(Sref, 'ref')
        error('The reference-orbit file must contain a top-level variable named ref.');
    end
    if ~isfield(Spde, 'results')
        error('The PDE-output file must contain a top-level variable named results.');
    end

    ref     = Sref.ref;
    results = Spde.results;

    if ~isfield(results, 'sim') || isempty(results.sim)
        error('results must contain a nonempty results.sim array.');
    end
    if ~isscalar(simIndex) || simIndex ~= round(simIndex) || ...
            simIndex < 1 || simIndex > numel(results.sim)
        error('simIndex must be an integer between 1 and %d.', numel(results.sim));
    end

    %% ------------------------------------------------------------
    % ODE reference orbit data
    %% ------------------------------------------------------------
    if isfield(ref, 'uCycle') && isfield(ref, 'vCycle')
        uRef = ref.uCycle(:);
        vRef = ref.vCycle(:);
    elseif isfield(ref, 'u') && isfield(ref, 'v')
        uRef = ref.u(:);
        vRef = ref.v(:);
    else
        error('ref must contain either uCycle/vCycle or u/v fields.');
    end

    if numel(uRef) ~= numel(vRef) || isempty(uRef)
        error('Reference-orbit arrays must be nonempty and have matching lengths.');
    end

    %% ------------------------------------------------------------
    % PDE mean-tail data
    %% ------------------------------------------------------------
    [tailData, tailLabel] = select_tail_data(Spde, simIndex);

    uPDE = tailData.ubar(:);
    vPDE = tailData.vbar(:);

    if numel(uPDE) ~= numel(vPDE) || isempty(uPDE)
        error('PDE mean-tail arrays must be nonempty and have matching lengths.');
    end

    goodRef = isfinite(uRef) & isfinite(vRef);
    goodPDE = isfinite(uPDE) & isfinite(vPDE);

    uRef = uRef(goodRef);
    vRef = vRef(goodRef);
    uPDE = uPDE(goodPDE);
    vPDE = vPDE(goodPDE);

    if isempty(uRef) || isempty(uPDE)
        error('No finite data remain for plotting.');
    end

    %% ------------------------------------------------------------
    % Plot
    %% ------------------------------------------------------------
    hFig = figure;
    ax = axes('Parent', hFig);
    hold(ax, 'on');

    plot(ax, uRef, vRef, '--', 'LineWidth', 1.5);
    plot(ax, uPDE, vPDE, '-',  'LineWidth', 1.5);

    plot(ax, uPDE(1),   vPDE(1),   's', 'MarkerSize', 6, ...
        'HandleVisibility', 'off');
    plot(ax, uPDE(end), vPDE(end), '^', 'MarkerSize', 6, ...
        'HandleVisibility', 'off');

    hold(ax, 'off');
    box(ax, 'on');
    grid(ax, 'on');

    xlabel(ax, 'u');
    ylabel(ax, 'v');
    legend(ax, {'ODE reference orbit', tailLabel}, 'Location', 'best');

    lambdaValue = get_lambda_value(Spde, ref);
    title(ax, sprintf('Test case 2: PDE mean tail vs ODE reference (\\lambda = %g)', ...
        lambdaValue));

end

% =====================================================================
% Local helpers
% =====================================================================
function [tailData, tailLabel] = select_tail_data(Spde, simIndex)

    % Prefer the postprocessed copy, because the public driver may remove
    % denseTail from results.sim before saving the MAT file.
    if isfield(Spde, 'post') && isfield(Spde.post, 'sim') && ...
            numel(Spde.post.sim) >= simIndex

        P = Spde.post.sim(simIndex);

        if isfield(P, 'tailData')
            if isfield(P.tailData, 'dense') && is_usable_tail(P.tailData.dense)
                tailData  = P.tailData.dense;
                tailLabel = 'PDE mean tail (dense)';
                return;
            end

            if isfield(P.tailData, 'tail') && is_usable_tail(P.tailData.tail)
                tailData  = P.tailData.tail;
                tailLabel = 'PDE mean tail';
                return;
            end
        end
    end

    % Fallback to the raw results structure.
    results = Spde.results;
    sim = results.sim(simIndex);

    if isfield(sim, 'denseTail') && is_usable_tail(sim.denseTail)
        tailData  = sim.denseTail;
        tailLabel = 'PDE mean tail (dense)';
        return;
    end

    if isfield(sim, 'tail') && is_usable_tail(sim.tail)
        tailData  = sim.tail;
        tailLabel = 'PDE mean tail';
        return;
    end

    error('No usable PDE mean-tail data found for simIndex = %d.', simIndex);

end

function tf = is_usable_tail(D)

    tf = isstruct(D) && ...
        isfield(D, 'ubar') && ~isempty(D.ubar) && ...
        isfield(D, 'vbar') && ~isempty(D.vbar);

end

function lambdaValue = get_lambda_value(Spde, ref)

    lambdaValue = 5;

    if isfield(Spde, 'results') && isfield(Spde.results, 'cfg') && ...
            isfield(Spde.results.cfg, 'lambda')
        lambdaValue = Spde.results.cfg.lambda;
    elseif isfield(Spde, 'cfg') && isfield(Spde.cfg, 'lambda')
        lambdaValue = Spde.cfg.lambda;
    elseif isfield(ref, 'lambda')
        lambdaValue = ref.lambda;
    end

end
