function ref = compute_reference_orbit()
% COMPUTE_REFERENCE_ORBIT  Compute an ODE reference limit cycle.
%
%   ref = COMPUTE_REFERENCE_ORBIT() integrates the spatially
%   homogeneous Rosenzweig--MacArthur ODE system with MATLAB's ODE45,
%   detects successive local maxima of u, and reconstructs one reference
%   period on a uniform time grid.
%
%   The parameter set is supplied by case_ex262(lambda).  For this example,
%   the spatially homogeneous Hopf value is lambda_0^H = (k - 1)/2 = 8.  To
%   observe convergence to the homogeneous limit cycle, choose
%   0 < lambda < 8.
%
%   The output ref is a structure containing the case data, period estimate,
%   sampled orbit, detected event states, and convergence diagnostics.
%
%   This file is intended to live in the ode_solver/ folder of the
%   Homogeneous limit cycles workflow.  It adds the neighboring config/
%   folder to the MATLAB path automatically.
%
%   See also ODE45, ODESET, CASE_EX262.

    %% Path setup
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    rootDir  = fileparts(thisDir);

    addpath(fullfile(rootDir, 'config'));

    %% User settings
    lambda = 5;

    % Initial condition for transient integration.  Leave empty to use the
    % automatic positive default below.
    y0User = [8; 3];

    % ODE45 tolerances for transient orbit finding.
    relTolTransient = 1e-12;
    absTolTransient = 1e-10;

    % ODE45 tolerances for final one-cycle reconstruction.
    relTolCycle = 1e-12;
    absTolCycle = 1e-12;

    % Maximum ODE45 step size.  This helps make event detection robust.
    maxStep = 0.001;

    % Transient integration controls.
    tChunk           = 100;
    maxTransientTime = 20e3;
    maxEvents        = 5000;

    % Convergence criteria for the detected sequence of maxima.
    stateTol     = 1e-10;
    periodTol    = 1e-10;
    closureTol   = 1e-8;
    scaleFloor   = 1.0;
    nConsecutive = 3;

    % Final one-cycle sampling.
    nCycleSample = 100000;

    % Progress output.  Since ODE45 is adaptive, this counter refers to
    % returned solver points rather than uniformly spaced physical times.
    outputEvery = 40000;
    verbose     = true;

    %% Case configuration
    cfg = case_ex262(lambda);

    lambdaHopf = (cfg.k - 1) / 2;
    if ~(cfg.lambda > 0 && cfg.lambda < lambdaHopf)
        error(['This reference-orbit computation is intended for ', ...
               '0 < lambda < (k - 1)/2 = %.12g. Current lambda = %.12g.'], ...
               lambdaHopf, cfg.lambda);
    end

    %% Initial condition
    if isempty(y0User)
        y0 = [0.8 * cfg.k; 0.2 * cfg.k / cfg.theta];
    else
        y0 = y0User(:);
    end

    if numel(y0) ~= 2 || any(~isfinite(y0)) || any(y0 <= 0)
        error('Initial condition y0 must be a finite positive 2-vector.');
    end

    %% Solver options
    optsTransient = odeset( ...
        'RelTol',  relTolTransient, ...
        'AbsTol',  absTolTransient, ...
        'MaxStep', maxStep, ...
        'Events',  @event_max_u);

    optsCycle = odeset( ...
        'RelTol',  relTolCycle, ...
        'AbsTol',  absTolCycle, ...
        'MaxStep', maxStep);

    %% Event storage and diagnostics
    tEvents       = zeros(1, maxEvents);
    yEvents       = zeros(2, maxEvents);
    periodVals    = nan(1, maxEvents);
    stateErrVals  = nan(1, maxEvents);
    periodErrVals = nan(1, maxEvents);

    nEvents = 0;
    nPass   = 0;
    nOutput = 0;

    convergedTransient = false;

    t0     = 0;
    yStart = y0;

    if verbose
        fprintf('\n');
        fprintf('========================================\n');
        fprintf('compute_reference_orbit_ode45\n');
        fprintf('========================================\n');
        fprintf('Case          : %s\n', cfg.caseName);
        fprintf('lambda        : %.12g\n', cfg.lambda);
        fprintf('Initial state : [%.12g, %.12g]^T\n', yStart(1), yStart(2));
        fprintf('stateTol      : %.3e\n', stateTol);
        fprintf('periodTol     : %.3e\n', periodTol);
        fprintf('closureTol    : %.3e\n', closureTol);
        fprintf('nConsecutive  : %d\n', nConsecutive);
        fprintf('========================================\n\n');
    end

    %% Transient integration loop
    while ~convergedTransient
        if t0 >= maxTransientTime
            break;
        end
        if nEvents >= maxEvents
            break;
        end

        tspan = [t0, min(t0 + tChunk, maxTransientTime)];

        [tSeg, ySeg, tEventSeg, yEventSeg] = ...
            ode45(@rhs_rm, tspan, yStart, optsTransient);

        if verbose && ~isempty(outputEvery) && isfinite(outputEvery) && outputEvery > 0
            for iOut = 2:numel(tSeg)
                nOutput = nOutput + 1;

                if mod(nOutput, outputEvery) == 0
                    fprintf('t = %12.6f   u = %.12g   v = %.12g\n', ...
                        tSeg(iOut), ySeg(iOut, 1), ySeg(iOut, 2));
                end
            end
        end

        t0     = tSeg(end);
        yStart = ySeg(end, :).';

        for jEvent = 1:numel(tEventSeg)
            te = tEventSeg(jEvent);
            ye = yEventSeg(jEvent, :).';

            % Skip duplicate events at chunk boundaries.
            if nEvents >= 1
                if abs(te - tEvents(nEvents)) <= 1e-12 * max(1, abs(te))
                    continue;
                end
            end

            nEvents = nEvents + 1;
            tEvents(nEvents)    = te;
            yEvents(:, nEvents) = ye;

            if verbose
                fprintf('Event %4d at t = %.12g : u = %.12g, v = %.12g\n', ...
                    nEvents, te, ye(1), ye(2));
            end

            % Three events give two successive period estimates.
            if nEvents >= 3
                Tcurr = tEvents(nEvents)     - tEvents(nEvents - 1);
                Tprev = tEvents(nEvents - 1) - tEvents(nEvents - 2);

                periodVals(nEvents) = Tcurr;

                stateErr = max( ...
                    abs(yEvents(:, nEvents) - yEvents(:, nEvents - 1)) ./ ...
                    max(abs(yEvents(:, nEvents)), scaleFloor));

                periodErr = abs(Tcurr - Tprev) / Tcurr;

                stateErrVals(nEvents)  = stateErr;
                periodErrVals(nEvents) = periodErr;

                if verbose
                    fprintf(['          Tcurr = %.12g, stateErr = %.3e, ', ...
                             'periodErr = %.3e\n'], ...
                        Tcurr, stateErr, periodErr);
                end

                if stateErr <= stateTol && periodErr <= periodTol
                    nPass = nPass + 1;
                else
                    nPass = 0;
                end

                if nPass >= nConsecutive
                    convergedTransient = true;
                    break;
                end
            end

            if nEvents >= maxEvents
                break;
            end
        end
    end

    %% Transient convergence check
    if ~convergedTransient
        error(['Transient orbit finding did not converge before hitting ', ...
               'a hard limit. Increase maxTransientTime or maxEvents, ', ...
               'or adjust the initial condition.']);
    end

    %% Accepted phase-zero state and period estimate
    y0Cycle = yEvents(:, nEvents);

    validPeriods = periodVals(isfinite(periodVals));
    nAvg = min(3, numel(validPeriods));
    T = mean(validPeriods(end - nAvg + 1:end));

    if verbose
        fprintf('\nAccepted phase-zero state:\n');
        fprintf('  y0Cycle = [%.12g, %.12g]^T\n', y0Cycle(1), y0Cycle(2));
        fprintf('Accepted period estimate:\n');
        fprintf('  T = %.12g\n\n', T);
    end

    %% Final one-cycle reconstruction on a uniform time grid
    tCycle = linspace(0, T, nCycleSample + 1);
    [tCycleOut, yCycleOut] = ode45(@rhs_rm, tCycle, y0Cycle, optsCycle);

    uCycle = yCycleOut(:, 1).';
    vCycle = yCycleOut(:, 2).';

    %% Closure check
    yEnd = yCycleOut(end, :).';

    closureErr = max(abs(yEnd - y0Cycle) ./ max(abs(y0Cycle), 1.0));

    if verbose
        fprintf('Final one-cycle closure error:\n');
        fprintf('  closureErr = %.3e\n\n', closureErr);
    end

    if closureErr > closureTol
        warning('compute_reference_orbit_ode45:closureDefect', ...
            ['Closure error %.3e exceeds closureTol %.3e. ', ...
             'The returned orbit may still be useful, but you may wish ', ...
             'to tighten tolerances or increase the transient criteria.'], ...
             closureErr, closureTol);
    end

    %% Trim storage and build output structure
    tEvents       = tEvents(1:nEvents);
    yEvents       = yEvents(:, 1:nEvents);
    periodVals    = periodVals(1:nEvents);
    stateErrVals  = stateErrVals(1:nEvents);
    periodErrVals = periodErrVals(1:nEvents);

    ref = struct();

    ref.cfg    = cfg;
    ref.lambda = cfg.lambda;

    ref.y0Transient = y0;
    ref.y0Cycle     = y0Cycle;
    ref.T           = T;

    ref.tCycle = tCycleOut(:).';
    ref.uCycle = uCycle;
    ref.vCycle = vCycle;
    ref.phase  = ref.tCycle / ref.T;

    ref.tEvents = tEvents;
    ref.yEvents = yEvents;

    ref.periodVals    = periodVals;
    ref.stateErrVals  = stateErrVals;
    ref.periodErrVals = periodErrVals;

    ref.closureErr = closureErr;

    ref.convergedTransient = convergedTransient;
    ref.convergedClosure   = (closureErr <= closureTol);
    ref.converged          = convergedTransient && (closureErr <= closureTol);

    ref.settings = struct( ...
        'relTolTransient',  relTolTransient, ...
        'absTolTransient',  absTolTransient, ...
        'relTolCycle',      relTolCycle, ...
        'absTolCycle',      absTolCycle, ...
        'maxStep',          maxStep, ...
        'tChunk',           tChunk, ...
        'maxTransientTime', maxTransientTime, ...
        'maxEvents',        maxEvents, ...
        'stateTol',         stateTol, ...
        'periodTol',        periodTol, ...
        'closureTol',       closureTol, ...
        'scaleFloor',       scaleFloor, ...
        'nConsecutive',     nConsecutive, ...
        'nCycleSample',     nCycleSample);

    if verbose
        fprintf('Done.\n');
        fprintf('convergedTransient = %d\n', ref.convergedTransient);
        fprintf('convergedClosure   = %d\n', ref.convergedClosure);
        fprintf('converged          = %d\n\n', ref.converged);
    end

    %% Nested functions
    function dydt = rhs_rm(~, y)
    %RHS_RM  Spatially homogeneous Rosenzweig--MacArthur ODE system.

        u = y(1);
        v = y(2);

        du = u * (1 - u / cfg.k) - cfg.m * u * v / (1 + u);
        dv = -cfg.theta * v + cfg.m * u * v / (1 + u);

        dydt = [du; dv];
    end

    function [value, isterminal, direction] = event_max_u(~, y)
    %EVENT_MAX_U  Detect local maxima of u by downward du/dt crossing.

        dydt = rhs_rm([], y);

        value      = dydt(1);
        isterminal = 0;
        direction  = -1;
    end
end
