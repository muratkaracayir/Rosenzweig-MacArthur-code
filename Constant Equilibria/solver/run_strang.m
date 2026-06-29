function results = run_strang(cfg, run)
%RUN_STRANG  Run Strang-splitting simulations for constant-equilibrium tests.
%
%   results = run_strang(cfg, run)
%
%   This routine is the main solver-level driver used by the public
%   constant-equilibrium workflow. It evolves the one-dimensional diffusive
%   Rosenzweig--MacArthur system from one or more prescribed initial
%   conditions and records basic convergence diagnostics toward the relevant
%   spatially homogeneous equilibrium.
%
%   Required cfg fields
%     cfg.d1, cfg.d2      diffusion coefficients
%     cfg.k              carrying capacity
%     cfg.theta          predator death rate
%     cfg.m              interaction strength
%     cfg.lambda         bifurcation parameter
%     cfg.u_star         coexistence prey value, usually lambda
%     cfg.v_star         coexistence predator value
%     cfg.caseName       short case identifier
%     cfg.ell or cfg.L   domain length parameter, with L = ell*pi
%
%   Required run fields
%     run.Nx             number of endpoint-including spatial grid points
%     run.dt             time step
%     run.Tfinal         final integration time
%
%   Optional run fields
%     run.saveEvery              snapshot stride, in time steps (default 100)
%     run.outputEvery            diagnostic stride, in time steps (default 1000)
%     run.probeStep              spatial subsampling stride for snapshots (default 1)
%     run.exploreIC              use extended initial-condition suite (default false)
%     run.seed                   random seed for reproducible ICs (default 1)
%     run.eqType                 'coexistence', 'preyonly', or auto if empty
%     run.verbose                print progress diagnostics (default true)
%     run.onlyICMode             restrict to one IC family (default '')
%     run.onlyICSample           restrict to one sample in that family (default [])
%     run.useStopping            enable early stopping (default false)
%     run.stopEvery              early-stopping check stride (default outputEvery)
%     run.eqTol                  equilibrium-distance tolerance (default 1e-6)
%     run.stepTol                step-difference tolerance (default 1e-8)
%     run.minStopTime            minimum time before stopping is allowed (default 0)
%     run.saveResults            save full suite to ../data (default false)
%     run.stopEqMode             'absolute' or 'relative' (default 'absolute')
%     run.eqScaleFloor           scale floor for relative stopping (default 1)
%     run.storeHistorySnapshots  store history/snapshot arrays (default true)
%
%   Equilibrium convention
%     If run.eqType is empty, the target is chosen automatically: for
%     lambda < k the target is the coexistence equilibrium (lambda,v_lambda),
%     while for lambda >= k the target is the prey-only equilibrium (k,0).
%
%   Dependencies in this folder
%     createIC, strang_split.

    % --- Required-field checks ---
    require_fields(run, {'Nx','dt','Tfinal'}, 'run');
    require_fields(cfg, {'d1','d2','k','theta','m','lambda', ...
                         'u_star','v_star','caseName'}, 'cfg');

    % --- Defaults ---
    if ~isfield(run, 'saveEvery'),   run.saveEvery   = 100;   end
    if ~isfield(run, 'outputEvery'), run.outputEvery = 1000;  end
    if ~isfield(run, 'probeStep'),   run.probeStep   = 1;     end
    if ~isfield(run, 'exploreIC'),   run.exploreIC   = false; end
    if ~isfield(run, 'seed'),        run.seed        = 1;     end
    if ~isfield(run, 'eqType'),      run.eqType      = '';    end
    if ~isfield(run, 'verbose'),     run.verbose     = true;  end
    if ~isfield(run, 'onlyICMode'),   run.onlyICMode   = '';  end
    if ~isfield(run, 'onlyICSample'), run.onlyICSample = [];  end
    if ~isfield(run, 'useStopping'), run.useStopping = false; end
    if ~isfield(run, 'stopEvery'),   run.stopEvery   = run.outputEvery; end
    if ~isfield(run, 'eqTol'),       run.eqTol       = 1e-6;  end
    if ~isfield(run, 'stepTol'),     run.stepTol     = 1e-8;  end
    if ~isfield(run, 'minStopTime'), run.minStopTime = 0;     end
    if ~isfield(run, 'saveResults'), run.saveResults = false; end
    if ~isfield(run, 'stopEqMode'),   run.stopEqMode   = 'absolute'; end
    if ~isfield(run, 'eqScaleFloor'), run.eqScaleFloor = 1; end
    if ~isfield(run, 'storeHistorySnapshots')
        run.storeHistorySnapshots = true;
    end

    validate_run_options(run);

    % --- Domain length ---
    if ~isfield(cfg, 'L')
        if ~isfield(cfg, 'ell')
            error('run_strang: cfg must contain either cfg.L or cfg.ell.');
        end
        cfg.L = cfg.ell * pi;
    end

    % --- Target constant equilibrium for diagnostics and stopping ---
    run.eqType = resolve_eq_type(cfg, run.eqType);
    [ueq, veq] = target_equilibrium(cfg, run.eqType);
    Ueq = [ueq * ones(1, run.Nx); veq * ones(1, run.Nx)];

    % --- Grid and time data ---
    x  = linspace(0, cfg.L, run.Nx);
    Nt = round(run.Tfinal / run.dt);

    idxProbe = 1:run.probeStep:run.Nx;
    nProbe   = numel(idxProbe);

    % --- Scaling for optional relative equilibrium stopping ---
    switch lower(run.stopEqMode)
        case 'absolute'
            UeqScale = [];

        case 'relative'
            UeqScale = max(abs(Ueq), run.eqScaleFloor * ones(size(Ueq)));

        otherwise
            error('run_strang: unknown run.stopEqMode: %s.', run.stopEqMode);
    end

    % --- Decide which ICs to run ---
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
                error('run_strang: unknown run.onlyICMode: %s.', run.onlyICMode);
        end

        if ~isempty(run.onlyICSample)
            s = run.onlyICSample;
            if ~isscalar(s) || s ~= round(s) || s < 1 || s > numel(icIndices)
                error('run.onlyICSample must be an integer between 1 and %d.', ...
                      numel(icIndices));
            end
            icIndices = icIndices(s);
        end
    end

    nIC = numel(icIndices);

    % --- Global metadata ---
    results.cfg    = cfg;
    results.run    = run;
    results.x      = x;
    results.xProbe = x(idxProbe);

    simTemplate = struct( ...
        'icType',       '', ...
        'icInfo',       [], ...
        'icIndex',      [], ...
        'U0',           [], ...
        'Ufinal',       [], ...
        'history',      [], ...
        'snapshots',    [], ...
        'stoppedEarly', false, ...
        'stopTime',     [] );

    results.sim = repmat(simTemplate, 1, nIC);

    % --- Loop over initial conditions ---
    for k = 1:nIC
        sim = simTemplate;
        icIndex = icIndices(k);

        [U0, icInfo] = createIC(x, cfg, run, icIndex);
        U = U0;

        prevStopCheckU = U;
        stoppedEarly   = false;
        stopTime       = run.Tfinal;

        % History / snapshot allocation
        if run.storeHistorySnapshots
            nHist = floor(Nt / run.outputEvery) + 1;
            history.t         = zeros(1, nHist);
            history.umin      = zeros(1, nHist);
            history.umax      = zeros(1, nHist);
            history.vmin      = zeros(1, nHist);
            history.vmax      = zeros(1, nHist);
            history.eqDistU   = zeros(1, nHist);
            history.eqDistV   = zeros(1, nHist);
            history.stepDistU = zeros(1, nHist);
            history.stepDistV = zeros(1, nHist);

            nSnaps = floor(Nt / run.saveEvery) + 1;
            snapshots.t = zeros(1, nSnaps);
            snapshots.u = zeros(nSnaps, nProbe);
            snapshots.v = zeros(nSnaps, nProbe);

            % Record initial state
            histCount = 1;
            snapCount = 1;
            prevHistU = U;

            history.t(histCount)         = 0;
            history.umin(histCount)      = min(U(1,:));
            history.umax(histCount)      = max(U(1,:));
            history.vmin(histCount)      = min(U(2,:));
            history.vmax(histCount)      = max(U(2,:));
            history.eqDistU(histCount)   = max(abs(U(1,:) - ueq));
            history.eqDistV(histCount)   = max(abs(U(2,:) - veq));
            history.stepDistU(histCount) = 0;
            history.stepDistV(histCount) = 0;
            histCount = histCount + 1;

            snapshots.t(snapCount)   = 0;
            snapshots.u(snapCount,:) = U(1, idxProbe);
            snapshots.v(snapCount,:) = U(2, idxProbe);
            snapCount = snapCount + 1;
        else
            history   = [];
            snapshots = [];
            histCount = [];
            snapCount = [];
            prevHistU = [];
        end

        % Time stepping
        for n = 1:Nt
            U = strang_split(U, cfg, run.Nx, run.dt);

            if mod(n, run.outputEvery) == 0
                tnow = n * run.dt;

                eqDistU_now = max(abs(U(1,:) - ueq));
                eqDistV_now = max(abs(U(2,:) - veq));

                if run.storeHistorySnapshots
                    history.t(histCount)         = tnow;
                    history.umin(histCount)      = min(U(1,:));
                    history.umax(histCount)      = max(U(1,:));
                    history.vmin(histCount)      = min(U(2,:));
                    history.vmax(histCount)      = max(U(2,:));
                    history.eqDistU(histCount)   = eqDistU_now;
                    history.eqDistV(histCount)   = eqDistV_now;
                    history.stepDistU(histCount) = max(abs(U(1,:) - prevHistU(1,:)));
                    history.stepDistV(histCount) = max(abs(U(2,:) - prevHistU(2,:)));

                    prevHistU = U;
                    histCount = histCount + 1;
                end

                if run.verbose
                    fprintf(['IC %2d/%2d  %-28s  t = %9.3f', ...
                             '  eqDistU = %9.3e  eqDistV = %9.3e\n'], ...
                        k, nIC, icInfo.label, tnow, eqDistU_now, eqDistV_now);
                end
            end

            if run.storeHistorySnapshots && mod(n, run.saveEvery) == 0
                tnow = n * run.dt;
                snapshots.t(snapCount)   = tnow;
                snapshots.u(snapCount,:) = U(1, idxProbe);
                snapshots.v(snapCount,:) = U(2, idxProbe);
                snapCount = snapCount + 1;
            end

            % Optional stopping criterion for convergence to the target equilibrium
            if run.useStopping && mod(n, run.stopEvery) == 0
                tnow = n * run.dt;

                switch lower(run.stopEqMode)
                    case 'absolute'
                        eqErr = max(max(abs(U - Ueq)));

                    case 'relative'
                        eqErr = max(max(abs(U - Ueq) ./ UeqScale));
                end

                stepErr = max(max(abs(U - prevStopCheckU)));

                if tnow >= run.minStopTime && ...
                   eqErr <= run.eqTol && stepErr <= run.stepTol
                    stoppedEarly = true;
                    stopTime = tnow;

                    if run.verbose
                        fprintf(['Stopping early at t = %.3f', ...
                                 '  eqErr = %.3e  stepErr = %.3e\n'], ...
                            tnow, eqErr, stepErr);
                    end

                    break;
                end

                prevStopCheckU = U;
            end
        end

        % Trim unused entries
        if run.storeHistorySnapshots
            history.t         = history.t(1:histCount-1);
            history.umin      = history.umin(1:histCount-1);
            history.umax      = history.umax(1:histCount-1);
            history.vmin      = history.vmin(1:histCount-1);
            history.vmax      = history.vmax(1:histCount-1);
            history.eqDistU   = history.eqDistU(1:histCount-1);
            history.eqDistV   = history.eqDistV(1:histCount-1);
            history.stepDistU = history.stepDistU(1:histCount-1);
            history.stepDistV = history.stepDistV(1:histCount-1);

            snapshots.t = snapshots.t(1:snapCount-1);
            snapshots.u = snapshots.u(1:snapCount-1, :);
            snapshots.v = snapshots.v(1:snapCount-1, :);
        else
            history   = [];
            snapshots = [];
        end

        % Store result for this IC
        sim.icType       = icInfo.mode;
        sim.icInfo       = icInfo;
        sim.icIndex      = icIndex;
        sim.U0           = U0;
        sim.Ufinal       = U;
        sim.history      = history;
        sim.snapshots    = snapshots;
        sim.stoppedEarly = stoppedEarly;
        sim.stopTime     = stopTime;

        results.sim(k) = sim;
    end

    results = save_results(cfg, run, results);
end

function results = save_results(cfg, run, results)
%SAVE_RESULTS  Save full-suite constant-equilibrium results if requested.
%
%   Saving occurs only when run.saveResults is true and run.onlyICMode is
%   empty. Files are written to ../data relative to this file, i.e. to the
%   project-root data folder when run_strang.m is located in solver/.

    if ~run.saveResults
        return;
    end

    if ~isempty(run.onlyICMode)
        return;
    end

    saveName = build_save_name(cfg, run, results);
    savePath = versioned_save_path(saveName);

    summary.caseName    = cfg.caseName;
    summary.fileName    = saveName;
    summary.lambda      = cfg.lambda;
    summary.Nx          = run.Nx;
    summary.dt          = run.dt;
    summary.Tfinal      = run.Tfinal;
    summary.nIC         = numel(results.sim);
    summary.useStopping = run.useStopping;

    [ueq, veq] = target_equilibrium(cfg, run.eqType);

    summary.targetU = ueq;
    summary.targetV = veq;

    summary.stoppedEarly = [results.sim.stoppedEarly];
    summary.eqDistU = zeros(1, numel(results.sim));
    summary.eqDistV = zeros(1, numel(results.sim));

    for i = 1:numel(results.sim)
        Ufinal = results.sim(i).Ufinal;
        summary.eqDistU(i) = max(abs(Ufinal(1,:) - ueq));
        summary.eqDistV(i) = max(abs(Ufinal(2,:) - veq));
    end

    save(savePath, 'cfg', 'run', 'results', 'summary');

    results.savePath = savePath;

    fprintf('Saved results to:\n%s\n', savePath);
end

function saveName = build_save_name(cfg, run, results)
%BUILD_SAVE_NAME  Construct the primary output file name.

    k = numel(results.sim);

    saveName = sprintf('Constant_%s_lambda_%g_Nx%d_dt%g_k%d.mat', ...
        cfg.caseName, cfg.lambda, run.Nx, run.dt, k);
end

function savePath = versioned_save_path(saveName)
%VERSIONED_SAVE_PATH  Return a non-overwriting path in ../data.
%
%   If saveName already exists, append _v2, _v3, ... before .mat.

    thisFileDir = fileparts(mfilename('fullpath'));
    dataDir = fullfile(thisFileDir, '..', 'data');

    if exist(dataDir, 'dir') ~= 7
        mkdir(dataDir);
    end

    savePath = fullfile(dataDir, saveName);

    if exist(savePath, 'file') ~= 2
        return;
    end

    [baseName, ext] = split_file_name(saveName);

    version = 2;
    while true
        candidate = fullfile(dataDir, sprintf('%s_v%d%s', baseName, version, ext));
        if exist(candidate, 'file') ~= 2
            savePath = candidate;
            return;
        end
        version = version + 1;
    end
end

function [baseName, ext] = split_file_name(fileName)
%SPLIT_FILE_NAME  Split file name into base name and extension.

    dotPos = find(fileName == '.', 1, 'last');

    if isempty(dotPos)
        baseName = fileName;
        ext = '';
    else
        baseName = fileName(1:dotPos-1);
        ext = fileName(dotPos:end);
    end
end

function eqType = resolve_eq_type(cfg, eqType)
%RESOLVE_EQ_TYPE  Choose or validate the target equilibrium type.

    if isempty(eqType)
        if cfg.lambda < cfg.k
            eqType = 'coexistence';
        else
            eqType = 'preyonly';
        end
        return;
    end

    switch lower(eqType)
        case 'coexistence'
            if cfg.lambda >= cfg.k
                error(['run_strang: coexistence target requested, but ', ...
                       'lambda >= k. Use run.eqType = ''preyonly''.']);
            end
            eqType = 'coexistence';

        case 'preyonly'
            eqType = 'preyonly';

        otherwise
            error('run_strang: unknown run.eqType: %s.', eqType);
    end
end

function [ueq, veq] = target_equilibrium(cfg, eqType)
%TARGET_EQUILIBRIUM  Return the homogeneous target equilibrium.

    switch lower(eqType)
        case 'coexistence'
            ueq = cfg.u_star;
            veq = cfg.v_star;

        case 'preyonly'
            ueq = cfg.k;
            veq = 0;

        otherwise
            error('run_strang: unknown equilibrium type: %s.', eqType);
    end
end

function require_fields(s, names, structName)
%REQUIRE_FIELDS  Error if any required field is missing.

    for j = 1:numel(names)
        if ~isfield(s, names{j})
            error('run_strang: missing required field %s.%s.', ...
                  structName, names{j});
        end
    end
end

function validate_run_options(run)
%VALIDATE_RUN_OPTIONS  Basic public-facing checks for run options.

    if ~isscalar(run.Nx) || run.Nx ~= round(run.Nx) || run.Nx < 2
        error('run_strang: run.Nx must be an integer at least 2.');
    end

    if ~isscalar(run.dt) || ~isreal(run.dt) || ~isfinite(run.dt) || run.dt <= 0
        error('run_strang: run.dt must be a positive real scalar.');
    end

    if ~isscalar(run.Tfinal) || ~isreal(run.Tfinal) || ...
       ~isfinite(run.Tfinal) || run.Tfinal <= 0
        error('run_strang: run.Tfinal must be a positive real scalar.');
    end

    integerFields = {'saveEvery','outputEvery','probeStep','stopEvery'};
    for j = 1:numel(integerFields)
        value = run.(integerFields{j});
        if ~isscalar(value) || value ~= round(value) || value < 1
            error('run_strang: run.%s must be a positive integer.', ...
                  integerFields{j});
        end
    end

    if run.probeStep > run.Nx
        error('run_strang: run.probeStep must not exceed run.Nx.');
    end
end
