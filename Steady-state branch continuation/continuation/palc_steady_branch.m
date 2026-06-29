function PALC = palc_steady_branch(UV0, cfg0, steady, stepIdx0, opts)
%PALC_STEADY_BRANCH  Pseudo-arclength continuation for steady branches.
%
%   PALC = palc_steady_branch(UV0, cfg0, steady, stepIdx0, opts)
%
% Continues a branch of stationary solutions for the theta-fixed
% Rosenzweig--MacArthur workflow.  The routine performs its own PALC
% predictor-corrector steps and calls residual_and_jacobian_steady.m
% from the neighboring solver folder.
%
% Inputs:
%   UV0      : initial steady state, stored as an Nx-by-2 array.
%   cfg0     : configuration structure at the initial lambda value.
%   steady   : metadata structure inherited from the branch anchor.
%   stepIdx0 : initial global step index.  If omitted, it is set to 0.
%   opts     : optional continuation settings.
%
% The sign of opts.ds controls the initial lambda direction.  Positive ds
% initially increases lambda, while negative ds initially decreases lambda.
% The full PALC structure is saved at the end, and the last accepted steady
% state is always saved.

    % -------------------- nargin handling --------------------
    if nargin < 4 || isempty(stepIdx0)
        stepIdx0 = 0;
    end
    if nargin < 5 || isempty(opts)
        opts = struct();
    end

    add_solver_path();

    % -------------------- defaults --------------------
    opts.nSteps   = get_def(opts,'nSteps', 1143);

    opts.ds       = get_def(opts,'ds', -5e-2);
    opts.dsMin    = get_def(opts,'dsMin', 1e-6);
    opts.dsMax    = get_def(opts,'dsMax', 5e-2);  % max |ds|

    opts.maxCorr  = get_def(opts,'maxCorr', 15);
    opts.rTol     = get_def(opts,'rTol', 1e-10);
    opts.sTol     = get_def(opts,'sTol', 1e-12);
    opts.bt       = get_def(opts,'bt', true);
    opts.btShrink = get_def(opts,'btShrink', 0.5);
    opts.btMin    = get_def(opts,'btMin', 1e-6);

    opts.fdLamRel = get_def(opts,'fdLamRel', 1e-7);
    opts.verbose  = get_def(opts,'verbose', true);

    % Saving SS control
    opts.saveEvery = get_def(opts,'saveEvery', 10); % 1 save each; 0 only critical; else periodic

    % Always save PALC at end
    opts.savePALCAtEnd = get_def(opts,'savePALCAtEnd', true);

    % Diagnostics
    opts.specModesK           = get_def(opts,'specModesK', 32);
    opts.specJumpThreshL1     = get_def(opts,'specJumpThreshL1', 0.40);
    opts.stopOnAbruptSpectrum = get_def(opts,'stopOnAbruptSpectrum', true);

    % Duplicate lambda stopping
    opts.stopOnDupLambdaMode = get_def(opts,'stopOnDupLambdaMode','hard'); % 'hard'|'none'
    opts.dupTolLambda        = get_def(opts,'dupTolLambda', 1e-10);

    % Collapse-to-equilibrium stopping
    opts.stopOnCollapsedEq    = get_def(opts,'stopOnCollapsedEq', true);
    opts.collapseAmpThresh    = get_def(opts,'collapseAmpThresh', 1e-8);
    opts.collapsePersistSteps = get_def(opts,'collapsePersistSteps', 3);

    % ds-underflow stopping
    opts.stopOnDsUnderflow    = get_def(opts,'stopOnDsUnderflow', true);

    % Tail energy sanity check (warn)
    opts.tailCheckEnabled  = get_def(opts,'tailCheckEnabled', true);
    opts.tailFracStart     = get_def(opts,'tailFracStart', 0.25);
    opts.tailEnergyThresh  = get_def(opts,'tailEnergyThresh', 1e-8);

    % Fold suppression near known branch endpoints.
    % By default, the endpoint values are read from cfg.route.lambdaL/R
    % when those metadata fields are present.
    opts.foldSuppLamTol = get_def(opts,'foldSuppLamTol', 2e-5);
    opts.foldSuppAmpTol = get_def(opts,'foldSuppAmpTol', 5e-4);
    opts = complete_fold_suppression_points(opts, cfg0);
    % --------------------------------------------------

    validate_inputs(UV0, cfg0, steady, stepIdx0, opts);

    cfg = cfg0;
    Nx = cfg.Nx;
    if size(UV0,1) ~= Nx || size(UV0,2) ~= 2
        error('palc_steady_branch:InvalidUV0', ...
            'UV0 must be an Nx-by-2 array with Nx = %d.', Nx);
    end

    if ~isfield(cfg,'L') || isempty(cfg.L)
        cfg.L = cfg.ell * pi;
    end
    if ~isfield(cfg,'x') || isempty(cfg.x)
        cfg.x = linspace(0, cfg.L, Nx)'; % includes endpoints
    end

    pack   = @(UVin) [UVin(:,1); UVin(:,2)];
    unpack = @(w) [w(1:Nx), w(Nx+1:end)];

    cfg = apply_thetaFixed(cfg, cfg.lambda);

    % Initial residual
    [R0, ~] = residual_and_jacobian_steady(UV0, cfg);
    res0 = norm(R0,2);
    if opts.verbose
        fprintf('PALC init: lambda=%.12g, ||F||=%.3e\n', cfg.lambda, res0);
    end

    % Initial tangent
    [tU, tlam] = compute_tangent(UV0, cfg, opts);

    ds = opts.ds;
    if ds == 0, error('opts.ds cannot be 0.'); end
    ds = sign(ds) * min(opts.dsMax, abs(ds));

    % Enforce INITIAL direction in lambda:
    % lam_pred = lam + ds*tlam, want sign(lam_pred-lam) = sign(ds)
    if sign(ds*tlam) ~= sign(ds)
        tU   = -tU;
        tlam = -tlam;
    end

    if opts.verbose
        fprintf('  Tangent init: tlam=%.3e, ds=%.3e\n', tlam, ds);
        if ds > 0
            fprintf('  Direction rule: ds>0 => lambda increases initially (ds*tlam>0 enforced).\n');
        else
            fprintf('  Direction rule: ds<0 => lambda decreases initially (ds*tlam<0 enforced).\n');
        end
    end

    % Allocate outputs (arrays)
    PALC = struct();
    PALC.caseName = get_case_name(cfg);
    PALC.route    = get_route_type(cfg);
    PALC.opts     = opts;
    PALC.stepIdx0 = stepIdx0;

    PALC.lambdas    = nan(opts.nSteps+1,1);
    PALC.resNorms   = nan(opts.nSteps+1,1);
    PALC.ampEq      = nan(opts.nSteps+1,1);
    PALC.ampEqNorm  = nan(opts.nSteps+1,1);

    % diag as struct-of-arrays
    PALC.diag = struct();
    PALC.diag.tlam            = nan(opts.nSteps+1,1);
    PALC.diag.dlam_ds_emp     = nan(opts.nSteps+1,1);
    PALC.diag.modeDominant    = nan(opts.nSteps+1,1);
    PALC.diag.specDistL1      = nan(opts.nSteps+1,1);
    PALC.diag.tailEnergyRatio = nan(opts.nSteps+1,1);
    PALC.diag.tailWarn        = false(opts.nSteps+1,1);

    % variable-length per step
    PALC.diag.modeAmp = cell(opts.nSteps+1,1);
    PALC.diag.specP   = cell(opts.nSteps+1,1);

    % fold log (store all folds)
    PALC.folds = struct('step',{},'lambdaEst',{},'lastLambda',{},'currLambda',{},'lastTlam',{},'currTlam',{},'alpha',{});

    % warnings
    PALC.warn = struct();
    PALC.warn.tailEnergy = false;
    PALC.warn.tailEnergySteps = [];
    PALC.warn.tailEnergyRatioMax = NaN;

    % stopping
    PALC.stopReason = 'none';
    PALC.stopMeta   = struct();

    % Step 0 store (seed)
    PALC.lambdas(1)  = cfg.lambda;
    PALC.resNorms(1) = res0;
    [a0, a0n] = compute_amp_eq_numbers(UV0, cfg);
    PALC.ampEq(1)     = a0;
    PALC.ampEqNorm(1) = a0n;

    [modeAmp0, modeDom0, specP0] = compute_mode_diagnostics(UV0, cfg, opts);
    PALC.diag.tlam(1)         = tlam;
    PALC.diag.modeDominant(1) = modeDom0;
    PALC.diag.modeAmp{1}      = modeAmp0;
    PALC.diag.specP{1}        = specP0;
    PALC.diag.specDistL1(1)   = NaN;

    if opts.tailCheckEnabled
        tailRatio0 = compute_tail_energy_ratio(UV0, cfg, opts);
        PALC.diag.tailEnergyRatio(1) = tailRatio0;
        if isfinite(tailRatio0) && tailRatio0 > opts.tailEnergyThresh
            PALC.diag.tailWarn(1) = true;
            PALC.warn.tailEnergy = true;
            PALC.warn.tailEnergySteps = [PALC.warn.tailEnergySteps; stepIdx0];
            PALC.warn.tailEnergyRatioMax = tailRatio0;
            warn_tail_energy(stepIdx0, cfg.lambda, tailRatio0, opts);
        end
    end

    % State
    UV  = UV0;
    lam = cfg.lambda;

    last_lam  = lam;
    last_tlam = tlam;

    prevSpecP    = specP0;
    prevModeDom  = modeDom0;

    seenPlus  = lam;
    seenMinus = lam;
    collapseCount = 0;

    step = 0;
    while step < opts.nSteps
        step = step + 1;
        gStep = stepIdx0 + step;  % global step id used in filenames, printing, tables

        if opts.verbose
            fprintf('\nPALC step %d/%d (ds=%.3e): lambda=%.12g\n', gStep, stepIdx0 + opts.nSteps, ds, lam);
        end

        % predictor
        w = pack(UV);
        w_pred   = w + ds*tU;
        lam_pred = lam + ds*tlam;

        cfg_pred = cfg;
        cfg_pred = apply_thetaFixed(cfg_pred, lam_pred);

        % corrector retries
        maxRetry = 50;
        retry = 0;
        ok = false;

        while ~ok
            retry = retry + 1;
            if retry > maxRetry
                % Ensure last state saved
                save_SS_data_thetaFixed(UV, cfg, steady, +1, prevModeDom, gStep-1);

                PALC.stopReason = 'corrector_failed';
                PALC.stopMeta.step = gStep;
                PALC.stopMeta.lambda = lam;
                PALC.stopMeta.note = 'maxRetry exceeded';
                PALC = trim_PALC_arrays(PALC, step); % last accepted index = step
                save_PALC_at_end(PALC, cfg, lam, gStep-1);
                return;
            end

            [w_new, lam_new, ok, resNorm, tU_new, tlam_new] = ...
                corrector_newton_aug(w_pred, lam_pred, tU, tlam, cfg_pred, opts, unpack);

            if ok, break; end

            ds_old = ds;
            ds = sign(ds) * max(opts.dsMin, 0.5*abs(ds));
            if opts.verbose
                fprintf('  Corrector failed. Shrinking |ds|: %.3e -> %.3e and retrying.\n', abs(ds_old), abs(ds));
            end

            if abs(ds) <= opts.dsMin + eps
                % Ensure last state saved
                % (last accepted is current UV/cfg/lam, which corresponds to step gStep-1)
                % Determine a sid from "no move": fallback to ds*tlam
                sid_last = +1; if ds*tlam < 0, sid_last = -1; end
                save_SS_data_thetaFixed(UV, cfg, steady, sid_last, prevModeDom, gStep-1);

                if opts.stopOnDsUnderflow
                    PALC.stopReason = 'ds_underflow';
                    PALC.stopMeta.step = gStep;
                    PALC.stopMeta.lambda = lam;
                    PALC.stopMeta.ds = ds;
                    if opts.verbose
                        fprintf('  >>> STOP: ds_underflow at step %d (|ds|=%.3e)\n', gStep, abs(ds));
                    end
                    PALC = trim_PALC_arrays(PALC, step);
                    save_PALC_at_end(PALC, cfg, lam, gStep-1);
                    return;
                else
                    PALC.stopReason = 'corrector_failed';
                    PALC.stopMeta.step = gStep;
                    PALC.stopMeta.lambda = lam;
                    PALC.stopMeta.ds = ds;
                    PALC = trim_PALC_arrays(PALC, step);
                    save_PALC_at_end(PALC, cfg, lam, gStep-1);
                    return;
                end
            end

            % rebuild predictor with shrunken ds
            w = pack(UV);
            w_pred   = w + ds*tU;
            lam_pred = lam + ds*tlam;
            cfg_pred = cfg;
            cfg_pred = apply_thetaFixed(cfg_pred, lam_pred);
        end

        % accept
        UV  = unpack(w_new);
        lam = lam_new;
        dlam_step = lam - last_lam;

        cfg = cfg_pred;
        cfg = apply_thetaFixed(cfg, lam);

        PALC.lambdas(step+1)  = lam;
        PALC.resNorms(step+1) = resNorm;

        [aEq, aEqN] = compute_amp_eq_numbers(UV, cfg);
        PALC.ampEq(step+1)     = aEq;
        PALC.ampEqNorm(step+1) = aEqN;

        dlam_ds_emp = (lam - last_lam)/ds;

        [modeAmp, modeDom, specP] = compute_mode_diagnostics(UV, cfg, opts);
        specDistL1 = NaN;
        if ~isempty(prevSpecP) && ~isempty(specP)
            specDistL1 = sum(abs(specP - prevSpecP));
        end

        PALC.diag.tlam(step+1)         = tlam;
        PALC.diag.dlam_ds_emp(step+1)  = dlam_ds_emp;
        PALC.diag.modeDominant(step+1) = modeDom;
        PALC.diag.modeAmp{step+1}      = modeAmp;
        PALC.diag.specP{step+1}        = specP;
        PALC.diag.specDistL1(step+1)   = specDistL1;

        if opts.verbose
            fprintf('  Accepted: lambda=%.12g, ||F||=%.3e, ampEqNorm=%.3e\n', lam, resNorm, aEqN);
            fprintf('  dlam/ds: empirical=%.3e, tangent(tlam)=%.3e\n', dlam_ds_emp, tlam);
            fprintf('  Dominant mode: %d\n', modeDom);
        end

        % Tail check
        tailEvent = false;
        if opts.tailCheckEnabled
            tailRatio = compute_tail_energy_ratio(UV, cfg, opts);
            PALC.diag.tailEnergyRatio(step+1) = tailRatio;
            if isfinite(tailRatio) && tailRatio > opts.tailEnergyThresh
                PALC.diag.tailWarn(step+1) = true;
                PALC.warn.tailEnergy = true;
                PALC.warn.tailEnergySteps = [PALC.warn.tailEnergySteps; gStep];
                if ~isfinite(PALC.warn.tailEnergyRatioMax)
                    PALC.warn.tailEnergyRatioMax = tailRatio;
                else
                    PALC.warn.tailEnergyRatioMax = max(PALC.warn.tailEnergyRatioMax, tailRatio);
                end
                warn_tail_energy(gStep, lam, tailRatio, opts);
                tailEvent = true;
            end
        end

        % Fold detection
        foldEvent = false;
        if (last_tlam ~= 0) && (tlam ~= 0) && (sign(tlam) ~= sign(last_tlam))
            if ~suppress_fold_here(lam, aEqN, opts)
                alpha = last_tlam/(last_tlam - tlam);
                lam_est = last_lam + alpha*(lam - last_lam);

                foldEvent = true;
                fprintf('  >>> FOLD DETECTED at step %d, est lambda %.12g\n', gStep, lam_est);

                kFold = numel(PALC.folds) + 1;
                PALC.folds(kFold).step = gStep;
                PALC.folds(kFold).lambdaEst = lam_est;
                PALC.folds(kFold).lastLambda = last_lam;
                PALC.folds(kFold).currLambda = lam;
                PALC.folds(kFold).lastTlam = last_tlam;
                PALC.folds(kFold).currTlam = tlam;
                PALC.folds(kFold).alpha = alpha;
            else
                if opts.verbose
                    fprintf('  (Fold suppressed near bifurcation; ampEqNorm=%.3e)\n', aEqN);
                end
            end
        end

        % Mode switch
        modeSwitchEvent = (~isnan(prevModeDom)) && (~isnan(modeDom)) && (modeDom ~= prevModeDom);
        if modeSwitchEvent
            fprintf('  >>> Dominant mode switch: %d -> %d at step %d (lambda=%.12g)\n', prevModeDom, modeDom, gStep, lam);
        end

        % Abrupt spectrum jump stop
        if opts.stopOnAbruptSpectrum && isfinite(specDistL1) && (~modeSwitchEvent)
            if specDistL1 > opts.specJumpThreshL1
                % Ensure last accepted SS saved
                sid_last = +1; if dlam_step < 0, sid_last = -1; end
                save_SS_data_thetaFixed(UV, cfg, steady, sid_last, modeDom, gStep);

                PALC.stopReason = 'abrupt_spectrum_jump';
                PALC.stopMeta.step = gStep;
                PALC.stopMeta.lambda = lam;
                PALC.stopMeta.distL1 = specDistL1;
                fprintf('  >>> STOP: abrupt_spectrum_jump (L1=%.3f > %.3f)\n', specDistL1, opts.specJumpThreshL1);
                PALC = trim_PALC_arrays(PALC, step+1);
                save_PALC_at_end(PALC, cfg, lam, gStep);
                return;
            end
        end

        % Sheet label based on actual lambda direction
        if dlam_step < 0
            sid = -1;
        elseif dlam_step > 0
            sid = +1;
        else
            sid = +1;
            if ds*tlam < 0, sid = -1; end
        end

        % Duplicate lambda stop (as before)
        if ~strcmpi(opts.stopOnDupLambdaMode,'none')
            if sid > 0
                [isDup, minDist] = is_duplicate_lambda(lam, seenPlus, opts.dupTolLambda);
                if isDup
                    % Ensure last accepted SS saved
                    save_SS_data_thetaFixed(UV, cfg, steady, sid, modeDom, gStep);

                    PALC.stopReason = 'duplicate_lambda_same_sheet';
                    PALC.stopMeta.step = gStep;
                    PALC.stopMeta.lambda = lam;
                    PALC.stopMeta.sheet = 'plus';
                    PALC.stopMeta.minDist = minDist;
                    fprintf('  >>> STOP: duplicate_lambda_same_sheet (+), |dlam|~%.3e\n', minDist);
                    PALC = trim_PALC_arrays(PALC, step+1);
                    save_PALC_at_end(PALC, cfg, lam, gStep);
                    return;
                else
                    seenPlus = [seenPlus; lam]; %#ok<AGROW>
                end
            else
                [isDup, minDist] = is_duplicate_lambda(lam, seenMinus, opts.dupTolLambda);
                if isDup
                    % Ensure last accepted SS saved
                    save_SS_data_thetaFixed(UV, cfg, steady, sid, modeDom, gStep);

                    PALC.stopReason = 'duplicate_lambda_same_sheet';
                    PALC.stopMeta.step = gStep;
                    PALC.stopMeta.lambda = lam;
                    PALC.stopMeta.sheet = 'minus';
                    PALC.stopMeta.minDist = minDist;
                    fprintf('  >>> STOP: duplicate_lambda_same_sheet (-), |dlam|~%.3e\n', minDist);
                    PALC = trim_PALC_arrays(PALC, step+1);
                    save_PALC_at_end(PALC, cfg, lam, gStep);
                    return;
                else
                    seenMinus = [seenMinus; lam]; %#ok<AGROW>
                end
            end
        end

        % Collapse stop
        if opts.stopOnCollapsedEq
            if isfinite(aEqN) && (aEqN < opts.collapseAmpThresh)
                collapseCount = collapseCount + 1;
            else
                collapseCount = 0;
            end
            if collapseCount >= opts.collapsePersistSteps
                % Ensure last accepted SS saved
                save_SS_data_thetaFixed(UV, cfg, steady, sid, modeDom, gStep);

                PALC.stopReason = 'collapsed_to_equilibrium';
                PALC.stopMeta.step = gStep;
                PALC.stopMeta.lambda = lam;
                PALC.stopMeta.ampEqNorm = aEqN;
                fprintf('  >>> STOP: collapsed_to_equilibrium (ampEqNorm=%.3e)\n', aEqN);
                PALC = trim_PALC_arrays(PALC, step+1);
                save_PALC_at_end(PALC, cfg, lam, gStep);
                return;
            end
        end

        % SS saving: periodic or genuinely critical
        % Tail-energy warnings should NOT override saveEvery.
        savedThisStep = false;
        criticalThisStep = foldEvent || modeSwitchEvent;

        if (opts.saveEvery > 0) && (mod(step, opts.saveEvery) == 0)
            save_SS_data_thetaFixed(UV, cfg, steady, sid, modeDom, gStep);
            savedThisStep = true;
        elseif criticalThisStep
            save_SS_data_thetaFixed(UV, cfg, steady, sid, modeDom, gStep);
            savedThisStep = true;
        end

        % Update history
        last_lam  = lam;
        last_tlam = tlam;
        prevSpecP = specP;
        prevModeDom = modeDom;

        % Update tangent
        tU   = tU_new;
        tlam = tlam_new;

        % Keep ds bounded, preserve sign
        ds = sign(ds) * min(opts.dsMax, abs(ds));
        if resNorm < 1e-12 && abs(ds) < opts.dsMax
            ds = sign(ds) * min(opts.dsMax, 1.2*abs(ds));
        end

        % If this is the last planned step, force-save last accepted state
        if step == opts.nSteps && ~savedThisStep
            save_SS_data_thetaFixed(UV, cfg, steady, sid, modeDom, gStep);
        end
    end

    PALC.stopReason = 'nSteps_reached';
    PALC.stopMeta.step = stepIdx0 + opts.nSteps;
    PALC.stopMeta.lambda = lam;

    PALC = trim_PALC_arrays(PALC, opts.nSteps+1);
    save_PALC_at_end(PALC, cfg, lam, stepIdx0 + opts.nSteps);
end

% ======= Local path setup =======
function add_solver_path()
    thisFile = mfilename('fullpath');
    continuationDir = fileparts(thisFile);
    rootDir = fileparts(continuationDir);
    solverDir = fullfile(rootDir, 'solver');

    if exist(solverDir, 'dir')
        addpath(solverDir);
    end
end

% ======= Input checks =======
function validate_inputs(UV0, cfg0, steady, stepIdx0, opts)
    if ~isnumeric(UV0) || ndims(UV0) ~= 2 || size(UV0, 2) ~= 2
        error('palc_steady_branch:InvalidUV0', ...
            'UV0 must be a numeric Nx-by-2 array.');
    end

    if ~isstruct(cfg0)
        error('palc_steady_branch:InvalidConfig', ...
            'cfg0 must be a structure.');
    end

    required = {'Nx','d1','d2','k','lambda'};
    for j = 1:numel(required)
        name = required{j};
        if ~isfield(cfg0, name) || isempty(cfg0.(name))
            error('palc_steady_branch:MissingConfigField', ...
                'cfg0.%s must be supplied.', name);
        end
    end

    if ~isfield(cfg0, 'L') && ~isfield(cfg0, 'ell')
        error('palc_steady_branch:MissingDomainLength', ...
            'cfg0.L or cfg0.ell must be supplied.');
    end

    if ~isstruct(steady)
        error('palc_steady_branch:InvalidSteadyMetadata', ...
            'steady must be a metadata structure.');
    end

    if ~isnumeric(stepIdx0) || ~isscalar(stepIdx0) || ~isfinite(stepIdx0)
        error('palc_steady_branch:InvalidStepIndex', ...
            'stepIdx0 must be a finite numeric scalar.');
    end

    if ~isstruct(opts)
        error('palc_steady_branch:InvalidOptions', ...
            'opts must be a structure.');
    end

    if isfield(opts, 'stopOnDupLambdaMode') && ...
            ~any(strcmpi(opts.stopOnDupLambdaMode, {'hard','none'}))
        error('palc_steady_branch:InvalidDuplicateMode', ...
            'opts.stopOnDupLambdaMode must be ''hard'' or ''none''.');
    end
end

% ======= Fold-suppression reference points =======
function opts = complete_fold_suppression_points(opts, cfg)
    if ~isfield(opts, 'bif') || isempty(opts.bif)
        opts.bif = struct();
    end

    if ~isstruct(opts.bif)
        error('palc_steady_branch:InvalidBifOptions', ...
            'opts.bif must be a structure when supplied.');
    end

    if ~isempty(fieldnames(opts.bif))
        return;
    end

    if isfield(cfg, 'route') && isstruct(cfg.route)
        if isfield(cfg.route, 'lambdaL') && isfinite_numeric_scalar(cfg.route.lambdaL)
            opts.bif.lambdaL = cfg.route.lambdaL;
        end
        if isfield(cfg.route, 'lambdaR') && isfinite_numeric_scalar(cfg.route.lambdaR)
            opts.bif.lambdaR = cfg.route.lambdaR;
        end
    end
end

function tf = isfinite_numeric_scalar(x)
    tf = isnumeric(x) && isscalar(x) && isfinite(x);
end

function routeType = get_route_type(cfg)
    if isfield(cfg, 'route') && isfield(cfg.route, 'type') && ...
            ~isempty(cfg.route.type)
        routeType = cfg.route.type;
    else
        routeType = 'thetaFixed';
    end
end

function caseName = get_case_name(cfg)
    if isfield(cfg, 'caseName') && ~isempty(cfg.caseName)
        caseName = cfg.caseName;
    else
        caseName = 'steady_state_branch';
    end
end

% ======= Save full PALC at end =======
function save_PALC_at_end(PALC, cfg, lambdaLast, stepLastGlobal)
    if isfield(PALC, 'opts') && isfield(PALC.opts, 'savePALCAtEnd') && ...
            ~PALC.opts.savePALCAtEnd
        if isfield(PALC, 'opts') && isfield(PALC.opts, 'verbose') && ...
                PALC.opts.verbose
            fprintf('PALC end save skipped because opts.savePALCAtEnd=false.\n');
        end
        return;
    end

    n = numel(PALC.lambdas);
    stepId = (PALC.stepIdx0 : PALC.stepIdx0 + n - 1).';

    PALC.tbl = table( ...
        stepId, ...
        PALC.lambdas(:), ...
        PALC.ampEq(:), ...
        PALC.ampEqNorm(:), ...
        PALC.diag.modeDominant(:), ...
        PALC.resNorms(:), ...
        'VariableNames', {'step','lambda','ampEq','ampEqNorm','modeDominant','resNorm'} );

    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(thisDir);
    dataFolder  = fullfile(projectRoot,'data');
    if ~exist(dataFolder,'dir'), mkdir(dataFolder); end

    fname = sprintf('PALC_thetaFixed_step%05d_end_lambda%.6g.mat', stepLastGlobal, lambdaLast);
    save(fullfile(dataFolder, fname), 'PALC');
    fprintf('PALC saved to %s\n', fullfile(dataFolder, fname));
end

% ======= Fold suppression near known bif points =======
function tf = suppress_fold_here(lam, ampEqNorm, opts)
    tf = false;

    if ampEqNorm > opts.foldSuppAmpTol
        return;
    end

    if ~isfield(opts,'bif') || isempty(opts.bif) || ~isstruct(opts.bif)
        return;
    end

    vals = [];
    fn = fieldnames(opts.bif);

    for k = 1:numel(fn)
        v = opts.bif.(fn{k});

        % collect any finite numeric entries, scalar or vector
        if isnumeric(v) && ~isempty(v)
            v = v(:);
            v = v(isfinite(v));
            vals = [vals; v]; %#ok<AGROW>
        end
    end

    if isempty(vals)
        return;
    end

    if min(abs(lam - vals)) < opts.foldSuppLamTol
        tf = true;
    end
end

% ======= Tangent computation =======
function [tU, tlam] = compute_tangent(UV, cfg, opts)
    Nx = cfg.Nx;
    w = [UV(:,1); UV(:,2)];
    [~, J] = residual_and_jacobian_steady(UV, cfg);
    J = sparse(J);
    Fl = finite_diff_Flambda(w, cfg, opts, Nx);

    tlam = 1.0;
    tU = J \ (-Fl);

    nrm = sqrt(norm(tU,2)^2 + tlam^2);
    tU = tU / nrm;
    tlam = tlam / nrm;
end

% ======= Corrector =======
function [w, lam, ok, resNorm, tU_new, tlam_new] = corrector_newton_aug(w, lam, tU, tlam, cfg, opts, unpack)
    Nx = cfg.Nx;
    ok = false;
    tU_new = tU; tlam_new = tlam;
    w_pred = w; lam_pred = lam;

    for it = 1:opts.maxCorr
        UV = unpack(w);
        cfg = apply_thetaFixed(cfg, lam);

        [R, J] = residual_and_jacobian_steady(UV, cfg);
        J = sparse(J);
        resNorm = norm(R,2);

        sRes = dot(w - w_pred, tU) + (lam - lam_pred)*tlam;

        if opts.verbose
            fprintf('   Corr %2d: ||F||=%.3e, |s|=%.3e\n', it, resNorm, abs(sRes));
        end

        if (resNorm < opts.rTol) && (abs(sRes) < opts.sTol)
            ok = true;
            break;
        end

        Fl = finite_diff_Flambda(w, cfg, opts, Nx);

        A = [J, Fl; tU.', tlam];
        rhs = -[R; sRes];

        step = A \ rhs;
        dw = step(1:2*Nx);
        dl = step(end);

        alpha = 1.0;
        merit0 = sqrt(resNorm^2 + sRes^2);

        if opts.bt
            while alpha >= opts.btMin
                w_try = w + alpha*dw;
                lam_try = lam + alpha*dl;

                cfg_try = apply_thetaFixed(cfg, lam_try);
                UV_try = unpack(w_try);

                R_try = residual_and_jacobian_steady(UV_try, cfg_try);
                s_try = dot(w_try - w_pred, tU) + (lam_try - lam_pred)*tlam;
                merit_try = sqrt(norm(R_try,2)^2 + s_try^2);

                if isfinite(merit_try) && (merit_try <= (1 - 1e-4*alpha)*merit0)
                    break;
                end
                alpha = alpha * opts.btShrink;
            end
        end

        w   = w + alpha*dw;
        lam = lam + alpha*dl;
    end

    cfg = apply_thetaFixed(cfg, lam);
    [tU_new, tlam_new] = compute_tangent(unpack(w), cfg, opts);

    if dot(tU_new, tU) + tlam_new*tlam < 0
        tU_new = -tU_new;
        tlam_new = -tlam_new;
    end

    cfg = apply_thetaFixed(cfg, lam);
    UV = unpack(w);
    Rfin = residual_and_jacobian_steady(UV, cfg);
    resNorm = norm(Rfin,2);
end

% ======= Finite difference F_lambda =======
function Fl = finite_diff_Flambda(w, cfg, opts, Nx)
    lam = cfg.lambda;
    h = opts.fdLamRel * max(1, abs(lam));
    lam_p = lam + h;
    lam_m = lam - h;

    if lam_m <= 1e-12
        lam_m = lam; lam_p = lam + h;
        oneSided = true;
    else
        oneSided = false;
    end

    unpack = @(wvec) [wvec(1:Nx), wvec(Nx+1:end)];
    UV = unpack(w);

    if ~oneSided
        cfgp = apply_thetaFixed(cfg, lam_p);
        cfgm = apply_thetaFixed(cfg, lam_m);
        Rp = residual_and_jacobian_steady(UV, cfgp);
        Rm = residual_and_jacobian_steady(UV, cfgm);
        Fl = (Rp - Rm) / (lam_p - lam_m);
    else
        cfgp = apply_thetaFixed(cfg, lam_p);
        Rp = residual_and_jacobian_steady(UV, cfgp);
        R0 = residual_and_jacobian_steady(UV, cfg);
        Fl = (Rp - R0) / (lam_p - lam);
    end
end

% ======= theta-fixed route update =======
function cfg = apply_thetaFixed(cfg, lambda)
    if ~isnumeric(lambda) || ~isscalar(lambda) || ~isfinite(lambda) || lambda <= 0
        error('palc_steady_branch:InvalidLambda', ...
            'lambda must be a positive finite scalar.');
    end

    if isfield(cfg, 'route') && isfield(cfg.route, 'type') && ...
            ~strcmp(cfg.route.type, 'thetaFixed')
        error('palc_steady_branch:UnsupportedRoute', ...
            'This public PALC workflow supports only cfg.route.type = ''thetaFixed''.');
    end

    cfg.lambda = lambda;

    if isfield(cfg, 'route') && isfield(cfg.route, 'fun') && ...
            isa(cfg.route.fun, 'function_handle')
        [mNow, thetaNow] = cfg.route.fun(lambda);
    else
        if isfield(cfg, 'route') && isfield(cfg.route, 'theta0') && ...
                ~isempty(cfg.route.theta0)
            thetaNow = cfg.route.theta0;
        elseif isfield(cfg, 'theta') && ~isempty(cfg.theta)
            thetaNow = cfg.theta;
        else
            error('palc_steady_branch:MissingTheta', ...
                'cfg.theta or cfg.route.theta0 must be supplied.');
        end
        mNow = thetaNow * (1 + 1./lambda);
    end

    if ~isnumeric(mNow) || ~isscalar(mNow) || ~isfinite(mNow) || mNow <= 0
        error('palc_steady_branch:InvalidM', ...
            'The theta-fixed route produced an invalid m value.');
    end
    if ~isnumeric(thetaNow) || ~isscalar(thetaNow) || ...
            ~isfinite(thetaNow) || thetaNow <= 0
        error('palc_steady_branch:InvalidTheta', ...
            'The theta-fixed route produced an invalid theta value.');
    end

    cfg.m      = mNow;
    cfg.theta  = thetaNow;
    cfg.u_star = lambda;

    k = cfg.k;
    cfg.v_star = ((k - lambda) * (1 + lambda)) / (k * cfg.m);
end

% ======= Amp numbers only =======
function [ampEq, ampEqNorm] = compute_amp_eq_numbers(UV, cfg)
    Nx = cfg.Nx;
    EQ = [cfg.u_star*ones(Nx,1), cfg.v_star*ones(Nx,1)];
    dUV = UV - EQ;

    ampEq = max( sqrt(dUV(:,1).^2 + dUV(:,2).^2) );

    eqMag = sqrt(cfg.u_star^2 + cfg.v_star^2);
    if eqMag > 0
        ampEqNorm = ampEq / eqMag;
    else
        ampEqNorm = NaN;
    end
end

% ======= Mode diagnostics =======
function [modeAmp, modeDom, specP] = compute_mode_diagnostics(UV, cfg, opts)
    Nx = cfg.Nx; L = cfg.L; x = cfg.x(:);
    K = min(opts.specModesK, Nx-1);
    if K < 1
        modeAmp = []; modeDom = NaN; specP = [];
        return;
    end

    du = UV(:,1) - cfg.u_star;
    dv = UV(:,2) - cfg.v_star;

    w = ones(Nx,1); w(1) = 0.5; w(end) = 0.5;
    dx = L/(Nx-1);

    modeAmp = zeros(K,1);
    for n = 1:K
        phi = cos(n*pi*x/L);
        au = dx * sum(w .* du .* phi);
        av = dx * sum(w .* dv .* phi);
        modeAmp(n) = sqrt(au^2 + av^2);
    end
    [~, modeDom] = max(modeAmp);

    E = sum(modeAmp);
    if E <= 0
        specP = modeAmp*0;
    else
        specP = modeAmp / E;
    end
end

% ======= Tail energy ratio =======
function tailRatio = compute_tail_energy_ratio(UV, cfg, opts)
    Nx = cfg.Nx;
    du = UV(:,1) - cfg.u_star;
    dv = UV(:,2) - cfg.v_star;

    cu = dct1_fast(du);
    cv = dct1_fast(dv);

    idx0 = 2:Nx;
    totalE = sum(cu(idx0).^2) + sum(cv(idx0).^2);
    if totalE <= 0 || ~isfinite(totalE)
        tailRatio = NaN;
        return;
    end

    tailStart = max(2, ceil(opts.tailFracStart * Nx));
    tailE = sum(cu(tailStart:Nx).^2) + sum(cv(tailStart:Nx).^2);
    tailRatio = tailE / totalE;
end

function c = dct1_fast(x)
    x = x(:);
    N = length(x);
    if N == 1, c = x; return; end
    y = [x; x(end-1:-1:2)];
    Y = real(fft(y));
    c = Y(1:N);
end

% ======= SS saver: step first in filename =======
function save_SS_data_thetaFixed(UV, cfg, steady, sid, modeDom, stepGlobal)
    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(thisDir);
    dataFolder  = fullfile(projectRoot,'data');
    if ~exist(dataFolder,'dir'), mkdir(dataFolder); end

    pm = 'plus'; if sid < 0, pm = 'minus'; end

    SS = struct();
    SS.stepIdx = stepGlobal;
    SS.lambda = cfg.lambda;
    SS.UV     = UV;
    SS.cfg    = cfg;
    SS.steady = steady;
    [aEq,aEqN] = compute_amp_eq_numbers(UV,cfg);
    SS.ampEq = aEq;
    SS.ampEqNorm = aEqN;

    fname = sprintf('SS_thetaFixed_step%05d_mode%d_%s_lambda%.6g.mat', stepGlobal, modeDom, pm, cfg.lambda);
    save(fullfile(dataFolder, fname), 'SS');
    fprintf('SS saved to %s\n', fullfile(dataFolder, fname));
end

% ======= Trimming =======
function PALC = trim_PALC_arrays(PALC, idx)
    if idx < 1, idx = 1; end
    PALC.lambdas   = PALC.lambdas(1:idx);
    PALC.resNorms  = PALC.resNorms(1:idx);
    PALC.ampEq     = PALC.ampEq(1:idx);
    PALC.ampEqNorm = PALC.ampEqNorm(1:idx);

    PALC.diag.tlam            = PALC.diag.tlam(1:idx);
    PALC.diag.dlam_ds_emp     = PALC.diag.dlam_ds_emp(1:idx);
    PALC.diag.modeDominant    = PALC.diag.modeDominant(1:idx);
    PALC.diag.specDistL1      = PALC.diag.specDistL1(1:idx);
    PALC.diag.tailEnergyRatio = PALC.diag.tailEnergyRatio(1:idx);
    PALC.diag.tailWarn        = PALC.diag.tailWarn(1:idx);

    PALC.diag.modeAmp = PALC.diag.modeAmp(1:idx);
    PALC.diag.specP   = PALC.diag.specP(1:idx);
end

% ======= Duplicate check =======
function [isDup, minDist] = is_duplicate_lambda(lam, seen, tol)
    if isempty(seen)
        isDup = false; minDist = Inf; return;
    end
    d = abs(seen - lam);
    minDist = min(d);
    isDup = (minDist < tol);
end

% ======= Tail warning =======
function warn_tail_energy(step, lam, tailRatio, opts)
    warning('PALC:TailEnergy', ...
        'Tail-energy warning at step %d (lambda=%.12g): tailRatio=%.3e > %.3e', ...
        step, lam, tailRatio, opts.tailEnergyThresh);
end

% ======= get_def =======
function v = get_def(S, name, def)
    if isfield(S,name) && ~isempty(S.(name))
        v = S.(name);
    else
        v = def;
    end
end