function results = run_strang(cfg, run)
%RUN_STRANG  PDE time-stepper for the homogeneous limit-cycle workflow.
%
%   results = run_strang(cfg, run)
%
%   This routine evolves the one-dimensional diffusive Rosenzweig--MacArthur
%   system with Strang splitting. It is specialized to the homogeneous
%   limit-cycle workflow: in addition to time stepping, it records tail
%   histories and compares the spatial mean of the PDE solution with a
%   precomputed ODE reference orbit supplied in run.refOrbit.
%
%   Required fields of run:
%     Nx          number of endpoint grid points
%     dt          time step
%     Tfinal      final integration time
%     refOrbit    reference-orbit structure, normally produced by
%                 build_ref_orbit_interpolant.m
%
%   Common optional fields:
%     outputEvery         coarse-history output stride
%     exploreIC           use the larger initial-condition suite
%     onlyICMode          restrict to one IC family; empty runs all families
%     onlyICSample        select one sample inside onlyICMode
%     seed                random seed used by createIC.m
%     verbose             print progress information
%
%     useStopping         enable convergence/alarm checks
%     stopEvery           stride for reference-orbit distance checks
%     homTol              spatial-homogeneity tolerance
%     orbitTol            distance-to-reference-orbit tolerance
%     minStopTime         earliest allowed criterion alarm time
%     nStopPass           number of consecutive successful checks required
%
%     useRefInterpolant   if true, evaluate run.refOrbit with eval_ref_orbit;
%                         otherwise use its stored discrete samples s,u,v
%     tailPeriods         number of reference periods retained after alarm
%     storeTailSnapshots  store full-state snapshots in the retained tail
%     storeDenseTail      store dense spatial means after alarm
%
%   This function does not save files. Saving is handled by the public
%   driver run_homogeneous_periodic.m.
    %% ------------------------------------------------------------
    % Defaults and checks
    %% ------------------------------------------------------------
    if ~isfield(run, 'outputEvery'),         run.outputEvery = 100; end
    if ~isfield(run, 'exploreIC'),           run.exploreIC   = false; end
    if ~isfield(run, 'onlyICMode'),          run.onlyICMode  = ''; end
    if ~isfield(run, 'onlyICSample'),        run.onlyICSample = []; end
    if ~isfield(run, 'seed'),                run.seed        = 1; end
    if ~isfield(run, 'verbose'),             run.verbose     = true; end

    if ~isfield(run, 'useStopping'),         run.useStopping = true; end
    if ~isfield(run, 'stopEvery'),           run.stopEvery   = run.outputEvery; end
    if ~isfield(run, 'homTol'),              run.homTol      = 1e-6; end
    if ~isfield(run, 'orbitTol'),            run.orbitTol    = 1e-6; end
    if ~isfield(run, 'minStopTime'),         run.minStopTime = 0; end
    if ~isfield(run, 'nStopPass'),           run.nStopPass   = 3; end
    if ~isfield(run, 'homScaleFloor'),       run.homScaleFloor   = 1; end
    if ~isfield(run, 'orbitScaleFloor'),     run.orbitScaleFloor = 1; end

    if ~isfield(run, 'useRefInterpolant'),   run.useRefInterpolant = false; end
    if ~isfield(run, 'refPhaseSearchN'),     run.refPhaseSearchN   = 4000; end

    if ~isfield(run, 'tailPeriods'),         run.tailPeriods = 4; end
    if ~isfield(run, 'storeTailSnapshots'),  run.storeTailSnapshots = false; end
    if ~isfield(run, 'snapshotEvery'),       run.snapshotEvery = run.outputEvery; end

    if ~isfield(run, 'storeDenseTail'),      run.storeDenseTail = false; end

    if ~isfield(run, 'refOrbit')
        error('run.refOrbit is required.');
    end
    ref = run.refOrbit;

    if ~isfield(ref, 'T') || ~isscalar(ref.T) || ~isfinite(ref.T) || ref.T <= 0
        error('run.refOrbit.T must be a positive finite scalar.');
    end

    if ~isfield(cfg, 'L')
        cfg.L = cfg.ell * pi;
    end

    if ~isscalar(run.Nx) || run.Nx < 2 || run.Nx ~= round(run.Nx)
        error('run.Nx must be an integer >= 2.');
    end
    if ~isscalar(run.dt) || run.dt <= 0 || ~isfinite(run.dt)
        error('run.dt must be a positive finite scalar.');
    end
    if ~isscalar(run.Tfinal) || run.Tfinal <= 0 || ~isfinite(run.Tfinal)
        error('run.Tfinal must be a positive finite scalar.');
    end

    %% ------------------------------------------------------------
    % Grid and time
    %% ------------------------------------------------------------
    x  = linspace(0, cfg.L, run.Nx);
    Nt = round(run.Tfinal / run.dt);

    tailDuration = run.tailPeriods * ref.T;
    tForceAlarm  = max(0, run.Tfinal - tailDuration);

    %% ------------------------------------------------------------
    % Prepare reference phase-search arrays
    %% ------------------------------------------------------------
    if run.useRefInterpolant
        tSearch = linspace(0, ref.T, run.refPhaseSearchN + 1);
        tSearch = tSearch(1:end-1);

        [uSearch, vSearch] = eval_ref_orbit(ref, tSearch);
        sSearch = tSearch / ref.T;
    else
        if ~all(isfield(ref, {'s','u','v'}))
            error('Discrete refOrbit must contain fields s, u, and v.');
        end

        sSearch = ref.s(:).';
        uSearch = ref.u(:).';
        vSearch = ref.v(:).';

        if numel(sSearch) >= 2 && abs(sSearch(end) - 1) < 1e-12
            sSearch = sSearch(1:end-1);
            uSearch = uSearch(1:end-1);
            vSearch = vSearch(1:end-1);
        end
    end

    %% ------------------------------------------------------------
    % IC selection
    %% ------------------------------------------------------------
    if run.exploreIC
        nICfull = 11;
    else
        nICfull = 3;
    end

    allICIndices = 1:nICfull;

    if isempty(run.onlyICMode)
        if ~isempty(run.onlyICSample)
            error('run.onlyICSample can only be used together with run.onlyICMode.');
        end
        icIndices = allICIndices;
    else
        switch lower(run.onlyICMode)
            case 'equilibriumperturbation'
                icIndices = 1;

            case 'smoothrandom'
                if run.exploreIC
                    icIndices = 2:6;
                else
                    icIndices = 2;
                end

            case 'largesmooth'
                if run.exploreIC
                    icIndices = 7:11;
                else
                    icIndices = 3;
                end

            otherwise
                error('Unknown run.onlyICMode: %s', run.onlyICMode);
        end

        if ~isempty(run.onlyICSample)
            s = run.onlyICSample;
            if ~isscalar(s) || s ~= round(s) || s < 1 || s > numel(icIndices)
                error('run.onlyICSample must be an integer between 1 and %d.', numel(icIndices));
            end
            icIndices = icIndices(s);
        end
    end

    nIC = numel(icIndices);

    %% ------------------------------------------------------------
    % Output struct
    %% ------------------------------------------------------------
    results = struct();
    results.cfg = cfg;
    results.run = run;
    results.x   = x;

    simTemplate = struct( ...
        'icInfo',           [], ...
        'icIndex',          [], ...
        'U0',               [], ...
        'Ufinal',           [], ...
        'alarmTriggered',   false, ...
        'alarmReason',      '', ...
        'tAlarm',           [], ...
        'tForceAlarm',      [], ...
        'tailStartTime',    [], ...
        'stopTime',         [], ...
        'orbitDistHistory', [], ...
        'tail',             [], ...
        'denseTail',        [], ...
        'tailSnapshots',    [], ...
        'finalHomErr',      [], ...
        'finalOrbitDist',   [], ...
        'finalRefPhase',    [], ...
        'finalRefState',    [] );

    results.sim = repmat(simTemplate, 1, nIC);

    %% ------------------------------------------------------------
    % Main IC loop
    %% ------------------------------------------------------------
    for k = 1:nIC
        icIndex = icIndices(k);
        [U0, icInfo] = createIC(x, cfg, run, icIndex);

        U = U0;
        tNow = 0;
        passCount = 0;

        alarmTriggered = false;
        alarmReason = '';
        tAlarm = NaN;
        stopTime = run.Tfinal;
        targetStopTime = run.Tfinal;

        % Full temporary coarse history (trimmed to tail at the end)
        nHistMax = floor(Nt / run.outputEvery) + 5;
        hist = struct();
        hist.t         = zeros(1, nHistMax);
        hist.ubar      = zeros(1, nHistMax);
        hist.vbar      = zeros(1, nHistMax);
        hist.umin      = zeros(1, nHistMax);
        hist.umax      = zeros(1, nHistMax);
        hist.vmin      = zeros(1, nHistMax);
        hist.vmax      = zeros(1, nHistMax);
        hist.spreadU   = zeros(1, nHistMax);
        hist.spreadV   = zeros(1, nHistMax);
        hist.homErr    = zeros(1, nHistMax);
        hist.orbitDist = zeros(1, nHistMax);
        hist.refPhase  = zeros(1, nHistMax);

        iHist = 0;
        M = compute_metrics(U);
        append_history(0, M);
        
        % Orbit-distance history on the stopEvery schedule only
        nOrbitHistMax = floor(Nt / run.stopEvery) + 5;
        orbitHist = struct();
        orbitHist.t         = zeros(1, nOrbitHistMax);
        orbitHist.orbitDist = zeros(1, nOrbitHistMax);

        iOrbitHist = 0;
        append_orbit_history(0, M);

        % Optional full-state snapshots (trimmed to tail at the end)
        if run.storeTailSnapshots
            nSnapMax = floor(Nt / run.snapshotEvery) + 5;
            snapT = zeros(1, nSnapMax);
            snapU = cell(1, nSnapMax);
            iSnap = 1;
            snapT(iSnap) = 0;
            snapU{iSnap} = U;
        else
            snapT = [];
            snapU = {};
            iSnap = 0;
        end

        % Optional dense tail after alarm
        if run.storeDenseTail
            nDenseMax = Nt + 5;
            denseT    = zeros(1, nDenseMax);
            denseUbar = zeros(1, nDenseMax);
            denseVbar = zeros(1, nDenseMax);
            iDense = 0;
        else
            denseT = [];
            denseUbar = [];
            denseVbar = [];
            iDense = 0;
        end

        lastMetrics = M;
        lastMetricsTime = 0;

        % Forced alarm already at t = 0 if the requested tail spans the whole run
        if ~alarmTriggered && (tForceAlarm <= 0)
            alarmTriggered = true;
            alarmReason = 'forced';
            tAlarm = 0;
            targetStopTime = run.Tfinal;

            if run.storeDenseTail
                append_dense_tail(0, M.ubar, M.vbar);
            end

            if run.verbose
                fprintf('\n');
                fprintf('========================================\n');
                fprintf('IC %d / %d : %s\n', k, nIC, icInfo.label);
                fprintf('========================================\n');
                fprintf(['FORCED ALARM at t = %10.4f   homErr = %9.3e   ', ...
                         'orbitDist = %9.3e\n'], ...
                    tAlarm, M.homErr, M.orbitDist);
            end
        elseif run.verbose
            fprintf('\n');
            fprintf('========================================\n');
            fprintf('IC %d / %d : %s\n', k, nIC, icInfo.label);
            fprintf('========================================\n');
        end

        for n = 1:Nt
            U = strang_split(U, cfg, run.Nx, run.dt);
            tNow = n * run.dt;

            needOutput = (mod(n, run.outputEvery) == 0) || (n == Nt);
            needStop   = run.useStopping && (mod(n, run.stopEvery) == 0);
            needSnap   = run.storeTailSnapshots && ...
                         ((mod(n, run.snapshotEvery) == 0) || (n == Nt));

            % ------------------------------------------------------------
            % Online diagnostics
            % ------------------------------------------------------------
            M = [];

            % Expensive reference-orbit metric only on the stopping schedule
            if needStop
                M = compute_metrics(U);
                lastMetrics = M;
                lastMetricsTime = tNow;

                append_orbit_history(tNow, M);

                if ~alarmTriggered
                    if (tNow >= run.minStopTime) && ...
                       (M.homErr <= run.homTol) && ...
                       (M.orbitDist <= run.orbitTol)
                        passCount = passCount + 1;
                    else
                        passCount = 0;
                    end

                    if passCount >= run.nStopPass
                        alarmTriggered = true;
                        alarmReason = 'criterion';
                        tAlarm = tNow;
                        targetStopTime = min(run.Tfinal, tAlarm + tailDuration);

                        % Save immediate dense-tail start
                        if run.storeDenseTail
                            append_dense_tail(tNow, M.ubar, M.vbar);
                        end

                        % Save immediate coarse-history point if this is not already an output checkpoint
                        if ~needOutput
                            append_history(tNow, M);
                        end

                        if run.verbose
                            fprintf(['ALARM at t = %10.4f   homErr = %9.3e   ', ...
                                     'orbitDist = %9.3e\n'], ...
                                tNow, M.homErr, M.orbitDist);
                        end
                    end
                end
            end

            % Cheap output diagnostics on the output schedule.
            % orbitDist/refPhase are copied from the latest expensive stopping check.
            if needOutput
                if isempty(M)
                    M = compute_basic_metrics(U);
                    M.orbitDist = lastMetrics.orbitDist;
                    M.refPhase  = lastMetrics.refPhase;
                    M.refState  = lastMetrics.refState;
                end

                append_history(tNow, M);

                if run.verbose
                    fprintf(['t = %10.4f   homErr = %9.3e   orbitDist = %9.3e', ...
                             '   refPhase = %7.4f\n'], ...
                        tNow, M.homErr, M.orbitDist, M.refPhase);
                end
            end

            % Forced alarm once the run enters the final tail window
            if ~alarmTriggered && (tNow >= tForceAlarm)
                if abs(lastMetricsTime - tNow) > 1e-14 * max(1, tNow)
                    M = compute_metrics(U);
                    lastMetrics = M;
                    lastMetricsTime = tNow;
                else
                    M = lastMetrics;
                end

                alarmTriggered = true;
                alarmReason = 'forced';
                tAlarm = tNow;
                targetStopTime = run.Tfinal;

                if run.storeDenseTail
                    append_dense_tail(tNow, M.ubar, M.vbar);
                end

                if ~needOutput
                    append_history(tNow, M);
                end

                if run.verbose
                    fprintf(['FORCED ALARM at t = %10.4f   homErr = %9.3e   ', ...
                             'orbitDist = %9.3e\n'], ...
                        tNow, M.homErr, M.orbitDist);
                end
            end

            % Dense tail storage at every step after alarm
            if run.storeDenseTail && alarmTriggered
                if ~(abs(tNow - tAlarm) <= 1e-14 * max(1, tAlarm))
                    append_dense_tail(tNow, mean(U(1,:)), mean(U(2,:)));
                end
            end

            if needSnap
                iSnap = iSnap + 1;
                snapT(iSnap) = tNow;
                snapU{iSnap} = U;
            end

            % If alarm has triggered, stop after the requested extra tail
            if alarmTriggered && (tNow >= targetStopTime)
                stopTime = tNow;
                break;
            end
        end

        % Final metrics only if stale
        if abs(lastMetricsTime - tNow) > 1e-14 * max(1, tNow)
            lastMetrics = compute_metrics(U);
        end

        % --------------------------------------------------------
        % Tail extraction
        % --------------------------------------------------------
        if alarmTriggered
            tailStartTime = tAlarm;
        else
            tailStartTime = max(0, tNow - tailDuration);
        end

        histNames = fieldnames(hist);
        for j = 1:numel(histNames)
            hist.(histNames{j}) = hist.(histNames{j})(1:iHist);
        end
        
        orbitHist.t         = orbitHist.t(1:iOrbitHist);
        orbitHist.orbitDist = orbitHist.orbitDist(1:iOrbitHist);

        maskTail = hist.t >= tailStartTime - 1e-14 * max(1, tailStartTime);

        tail = struct();
        for j = 1:numel(histNames)
            tail.(histNames{j}) = hist.(histNames{j})(maskTail);
        end

        if run.storeDenseTail && alarmTriggered
            denseTail = struct();
            denseTail.t    = denseT(1:iDense);
            denseTail.ubar = denseUbar(1:iDense);
            denseTail.vbar = denseVbar(1:iDense);
        else
            denseTail = struct();
            denseTail.t    = [];
            denseTail.ubar = [];
            denseTail.vbar = [];
        end

        if run.storeTailSnapshots
            snapT = snapT(1:iSnap);
            snapU = snapU(1:iSnap);

            maskSnap = snapT >= tailStartTime - 1e-14 * max(1, tailStartTime);

            tailSnapshots = struct();
            tailSnapshots.t = snapT(maskSnap);
            tailSnapshots.U = snapU(maskSnap);
        else
            tailSnapshots = [];
        end

        % --------------------------------------------------------
        % Store this simulation
        % --------------------------------------------------------
        sim = simTemplate;
        sim.icInfo         = icInfo;
        sim.icIndex        = icIndex;
        sim.U0             = U0;
        sim.Ufinal         = U;
        sim.alarmTriggered = alarmTriggered;
        sim.alarmReason    = alarmReason;
        sim.tAlarm         = tAlarm;
        sim.tForceAlarm    = tForceAlarm;
        sim.tailStartTime  = tailStartTime;
        sim.stopTime       = stopTime;
        sim.orbitDistHistory = orbitHist;
        sim.tail           = tail;
        sim.denseTail      = denseTail;
        sim.tailSnapshots  = tailSnapshots;
        sim.finalHomErr    = lastMetrics.homErr;
        sim.finalOrbitDist = lastMetrics.orbitDist;
        sim.finalRefPhase  = lastMetrics.refPhase;
        sim.finalRefState  = lastMetrics.refState;

        results.sim(k) = sim;
    end

    % ============================================================
    % Local helpers
    % ============================================================
    function M = compute_basic_metrics(Uloc)
        u = Uloc(1,:);
        v = Uloc(2,:);

        M.ubar = mean(u);
        M.vbar = mean(v);

        M.umin = min(u);
        M.umax = max(u);
        M.vmin = min(v);
        M.vmax = max(v);

        M.spreadU = M.umax - M.umin;
        M.spreadV = M.vmax - M.vmin;

        M.homErr = max([ ...
            M.spreadU / max(abs(M.ubar), run.homScaleFloor), ...
            M.spreadV / max(abs(M.vbar), run.homScaleFloor) ]);
    end

    function M = compute_metrics(Uloc)
        M = compute_basic_metrics(Uloc);

        du = abs(M.ubar - uSearch) ./ max(abs(uSearch), run.orbitScaleFloor);
        dv = abs(M.vbar - vSearch) ./ max(abs(vSearch), run.orbitScaleFloor);

        distVals = max(du, dv);
        [M.orbitDist, idxMin] = min(distVals);

        M.refPhase = sSearch(idxMin);
        M.refState = [uSearch(idxMin); vSearch(idxMin)];
    end

    function append_history(tVal, M)
        iHist = iHist + 1;
        hist.t(iHist)         = tVal;
        hist.ubar(iHist)      = M.ubar;
        hist.vbar(iHist)      = M.vbar;
        hist.umin(iHist)      = M.umin;
        hist.umax(iHist)      = M.umax;
        hist.vmin(iHist)      = M.vmin;
        hist.vmax(iHist)      = M.vmax;
        hist.spreadU(iHist)   = M.spreadU;
        hist.spreadV(iHist)   = M.spreadV;
        hist.homErr(iHist)    = M.homErr;
        hist.orbitDist(iHist) = M.orbitDist;
        hist.refPhase(iHist)  = M.refPhase;
    end

    function append_orbit_history(tVal, M)
        iOrbitHist = iOrbitHist + 1;
        orbitHist.t(iOrbitHist)         = tVal;
        orbitHist.orbitDist(iOrbitHist) = M.orbitDist;
    end

    function append_dense_tail(tVal, ubarVal, vbarVal)
        iDense = iDense + 1;
        denseT(iDense)    = tVal;
        denseUbar(iDense) = ubarVal;
        denseVbar(iDense) = vbarVal;
    end
end
