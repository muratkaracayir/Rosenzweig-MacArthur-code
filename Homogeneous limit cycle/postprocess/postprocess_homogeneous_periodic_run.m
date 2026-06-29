function post = postprocess_homogeneous_periodic_run(results, arg2, arg3)
%POSTPROCESS_HOMOGENEOUS_PERIODIC_RUN Postprocess homogeneous-periodic PDE runs.
%
%   post = postprocess_homogeneous_periodic_run(results)
%   post = postprocess_homogeneous_periodic_run(results, opts)
%   post = postprocess_homogeneous_periodic_run(results, simIndex)
%   post = postprocess_homogeneous_periodic_run(results, simIndex, opts)
%
%   This workflow-specific helper analyzes the output of the homogeneous
%   limit-cycle PDE time-stepping run. It compares the saved late-time PDE
%   trajectory with itself and with the ODE reference orbit stored in
%   results.run.refOrbit.
%
%   In suite mode, the function postprocesses all entries in results.sim and
%   returns post.sim(k). In single-run mode, it returns the postprocessing
%   structure for the requested simulation only.
%
%   The diagnostics include:
%     - an estimated period from saved late-time spatial means;
%     - a PDE-to-itself periodicity error, when dense-tail data are present;
%     - a PDE-to-reference-orbit error, when dense-tail data are present;
%     - a convergence-time diagnostic from period-window minima of the
%       orbit-distance history, when that history is present.
%
%   Optional fields in opts:
%     opts.preferDenseTail = true
%     opts.periodSearchRel = [0.8, 1.2]
%     opts.coarseSearchN   = 400
%     opts.refSearchN      = 400
%     opts.interpMethod    = 'pchip'
%     opts.scaleFloor      = 1
%     opts.minOverlapRel   = 0.75
%     opts.convNDips       = 5
%     opts.convRelTol      = 0.05
%     opts.convUseInterp   = true
%
%   This file belongs in the postprocess/ folder of the homogeneous
%   limit-cycle workflow. It does not advance either the ODE or PDE system.

    % ------------------------------------------------------------
    % Input parsing
    % ------------------------------------------------------------
    if nargin < 2
        simMode = 'all';
        opts = struct();
        simIndex = [];
    elseif nargin < 3
        if isstruct(arg2) || isempty(arg2)
            simMode = 'all';
            opts = arg2;
            simIndex = [];
        else
            simMode = 'one';
            simIndex = arg2;
            opts = struct();
        end
    else
        simMode = 'one';
        simIndex = arg2;
        opts = arg3;
    end

    if isempty(opts)
        opts = struct();
    end

    % ------------------------------------------------------------
    % Defaults
    % ------------------------------------------------------------
    if ~isfield(opts, 'preferDenseTail'),    opts.preferDenseTail = true; end
    if ~isfield(opts, 'periodSearchRel'),    opts.periodSearchRel = [0.8, 1.2]; end
    if ~isfield(opts, 'coarseSearchN'),      opts.coarseSearchN   = 400; end
    if ~isfield(opts, 'refSearchN'),         opts.refSearchN      = 400; end
    if ~isfield(opts, 'interpMethod'),       opts.interpMethod    = 'pchip'; end
    if ~isfield(opts, 'scaleFloor'),         opts.scaleFloor      = 1; end
    if ~isfield(opts, 'minOverlapRel'),      opts.minOverlapRel   = 0.75; end

    if ~isfield(opts, 'convNDips'),          opts.convNDips       = 5; end
    if ~isfield(opts, 'convRelTol'),         opts.convRelTol      = 0.05; end
    if ~isfield(opts, 'convUseInterp'),      opts.convUseInterp   = true; end

    % ------------------------------------------------------------
    % Basic validation
    % ------------------------------------------------------------
    if ~isstruct(results) || ~isfield(results, 'sim')
        error('results must be a struct containing results.sim.');
    end

    if ~isfield(results, 'run') || ~isfield(results.run, 'refOrbit')
        error('results.run.refOrbit is required.');
    end

    requiredRefFields = {'T', 's', 'u', 'v'};
    if ~all(isfield(results.run.refOrbit, requiredRefFields))
        error('results.run.refOrbit must contain fields T, s, u, and v.');
    end

    if strcmp(simMode, 'all')
        post = struct();
        post.meta = struct();
        post.meta.case   = get_field(results, 'cfg', 'caseName', '');
        post.meta.lambda = get_field(results, 'cfg', 'lambda', NaN);
        post.meta.dt     = get_field(results, 'run', 'dt', NaN);
        post.meta.Nx     = get_field(results, 'run', 'Nx', NaN);
        post.meta.nSim   = numel(results.sim);

        if isempty(results.sim)
            post.sim = struct([]);
        else
            firstSim = process_one(1);
            post.sim = repmat(firstSim, 1, numel(results.sim));
            post.sim(1) = firstSim;

            for k = 2:numel(results.sim)
                post.sim(k) = process_one(k);
            end
        end
    else
        if ~isscalar(simIndex) || simIndex ~= round(simIndex) || ...
           simIndex < 1 || simIndex > numel(results.sim)
            error('simIndex is out of range.');
        end
        post = process_one(simIndex);
    end

    % ============================================================
    % Internal single-run processor
    % ============================================================
    function out = process_one(k)
        sim  = results.sim(k);
        ref  = results.run.refOrbit;
        Tref = ref.T;

        % --------------------------------------------------------
        % Select late-time mean data for period estimation
        % --------------------------------------------------------
        useDenseForPeriod = false;
        if opts.preferDenseTail && isfield(sim, 'denseTail')
            D = sim.denseTail;
            if isfield(D, 't') && isfield(D, 'ubar') && isfield(D, 'vbar') && ...
               ~isempty(D.t)
                useDenseForPeriod = true;
            end
        end

        if useDenseForPeriod
            t = sim.denseTail.t(:);
            u = sim.denseTail.ubar(:);
            v = sim.denseTail.vbar(:);
            tailSource = 'dense';
        else
            if ~isfield(sim, 'tail') || ~all(isfield(sim.tail, {'t','ubar','vbar'}))
                error('No usable tail mean data found in sim.tail or sim.denseTail.');
            end
            t = sim.tail.t(:);
            u = sim.tail.ubar(:);
            v = sim.tail.vbar(:);
            tailSource = 'tail';
        end

        if numel(t) < 5
            error('Not enough tail points for period estimation.');
        end

        [t, ia] = unique(t, 'stable');
        u = u(ia);
        v = v(ia);

        if numel(t) < 5
            error('Tail data became too short after removing duplicate times.');
        end

        tailSpan = t(end) - t(1);
        if tailSpan <= 0
            error('Tail time span must be positive.');
        end

        % --------------------------------------------------------
        % Period estimate
        % --------------------------------------------------------
        tauMin = opts.periodSearchRel(1) * Tref;
        tauMax = opts.periodSearchRel(2) * Tref;

        minOverlapTime = opts.minOverlapRel * Tref;
        tauMaxEff = min(tauMax, tailSpan - minOverlapTime);

        if tauMaxEff <= tauMin
            error(['Tail is too short for the requested search interval and overlap. ', ...
                   'Increase the saved tail length or relax opts.minOverlapRel / opts.periodSearchRel.']);
        end

        Su = max(max(u) - min(u), opts.scaleFloor);
        Sv = max(max(v) - min(v), opts.scaleFloor);

        tauGrid = linspace(tauMin, tauMaxEff, opts.coarseSearchN);
        JGrid = zeros(size(tauGrid));

        for j = 1:numel(tauGrid)
            JGrid(j) = shift_objective(tauGrid(j));
        end

        [Jcoarse, idxBest] = min(JGrid);

        if numel(tauGrid) >= 3
            iL = max(idxBest - 1, 1);
            iR = min(idxBest + 1, numel(tauGrid));
            left  = tauGrid(iL);
            right = tauGrid(iR);

            if right > left
                [tauEst, JEst] = fminbnd(@shift_objective, left, right);
            else
                tauEst = tauGrid(idxBest);
                JEst   = JGrid(idxBest);
            end
        else
            tauEst = tauGrid(idxBest);
            JEst   = JGrid(idxBest);
        end

        % --------------------------------------------------------
        % Dense-tail diagnostics: PDE-to-itself and PDE-to-reference
        % --------------------------------------------------------
        haveDense = false;
        if isfield(sim, 'denseTail')
            D = sim.denseTail;
            if isfield(D, 't') && isfield(D, 'ubar') && isfield(D, 'vbar') && ...
               ~isempty(D.t)
                haveDense = true;
                td = D.t(:);
                ud = D.ubar(:);
                vd = D.vbar(:);

                [td, ia] = unique(td, 'stable');
                ud = ud(ia);
                vd = vd(ia);
            end
        end

        if haveDense
            Sud = max(max(ud) - min(ud), opts.scaleFloor);
            Svd = max(max(vd) - min(vd), opts.scaleFloor);

            % PDE-to-itself
            [tS, uS, vS, uShiftS, vShiftS] = build_shifted_overlap_dense(tauEst);
            dSelf = max([abs(uShiftS - uS) / Sud, abs(vShiftS - vS) / Svd], [], 2);

            selfOk   = true;
            selfD    = max(dSelf);
            selfDrms = sqrt(mean(dSelf.^2));

            % PDE-to-reference
            phiGrid = linspace(0, Tref, opts.refSearchN + 1);
            phiGrid = phiGrid(1:end-1);
            JR = zeros(size(phiGrid));

            for j = 1:numel(phiGrid)
                JR(j) = ref_objective(phiGrid(j));
            end

            [~, idxPhi] = min(JR);

            if numel(phiGrid) >= 3
                dphi = Tref / opts.refSearchN;
                left  = max(0, phiGrid(idxPhi) - dphi);
                right = min(Tref, phiGrid(idxPhi) + dphi);

                if right > left
                    [phiEst, JRef] = fminbnd(@ref_objective, left, right);
                else
                    phiEst = phiGrid(idxPhi);
                    JRef   = JR(idxPhi);
                end
            else
                phiEst = phiGrid(idxPhi);
                JRef   = JR(idxPhi);
            end

            [uR, vR, uRefFit, vRefFit] = build_ref_fit(phiEst);
            dRef = max([abs(uRefFit - uR) / Sud, abs(vRefFit - vR) / Svd], [], 2);

            refOk   = true;
            refD    = max(dRef);
            refDrms = sqrt(mean(dRef.^2));
        else
            selfOk   = false;
            selfD    = [];
            selfDrms = [];
            tS = []; uS = []; vS = []; uShiftS = []; vShiftS = []; dSelf = [];

            refOk   = false;
            refD    = [];
            refDrms = [];
            phiEst  = [];
            JRef    = [];
            phiGrid = [];
            JR      = [];
        end

        % --------------------------------------------------------
        % Convergence diagnostic from orbit-distance history
        % --------------------------------------------------------
        if isfield(sim, 'orbitDistHistory') && ...
           isstruct(sim.orbitDistHistory) && ...
           isfield(sim.orbitDistHistory, 't') && ...
           isfield(sim.orbitDistHistory, 'orbitDist') && ...
           ~isempty(sim.orbitDistHistory.t)

            orbitHist = sim.orbitDistHistory;
            th = orbitHist.t(:);
            dh = orbitHist.orbitDist(:);

            [th, ia] = unique(th, 'stable');
            dh = dh(ia);

            orbitHist.t = th;
            orbitHist.orbitDist = dh;

            conv = compute_convergence_from_history(th, dh, tauEst);
        else
            orbitHist = struct();
            orbitHist.t = [];
            orbitHist.orbitDist = [];

            conv = empty_conv_struct();
        end

        % --------------------------------------------------------
        % Package output for this simulation
        % --------------------------------------------------------
        out = struct();

        out.meta = struct();
        out.meta.case   = get_field(results, 'cfg', 'caseName', '');
        out.meta.lambda = get_field(results, 'cfg', 'lambda', NaN);
        out.meta.dt     = get_field(results, 'run', 'dt', NaN);
        out.meta.Nx     = get_field(results, 'run', 'Nx', NaN);
        out.meta.ic     = get_ic_label(sim);
        out.meta.alarm  = get_field(sim, [], 'alarmTriggered', []);
        out.meta.reason = get_field(sim, [], 'alarmReason', '');
        out.meta.tAlarm = get_field(sim, [], 'tAlarm', NaN);
        out.meta.tTail  = get_field(sim, [], 'tailStartTime', NaN);
        out.meta.tStop  = get_field(sim, [], 'stopTime', NaN);
        out.meta.source = tailSource;
        out.meta.Tref   = Tref;

        out.tailData = struct();
        out.tailData.tail  = get_field(sim, [], 'tail', []);
        out.tailData.dense = get_field(sim, [], 'denseTail', []);

        out.period = struct();
        out.period.Tref    = Tref;
        out.period.T       = tauEst;
        out.period.relErr  = abs(tauEst - Tref) / Tref;
        out.period.J       = JEst;
        out.period.J0      = Jcoarse;
        out.period.win     = [tauMin, tauMaxEff];
        out.period.overlap = minOverlapTime;

        out.self = struct();
        out.self.ok   = selfOk;
        out.self.D    = selfD;
        out.self.Drms = selfDrms;

        out.ref = struct();
        out.ref.ok   = refOk;
        out.ref.phi  = phiEst;
        out.ref.D    = refD;
        out.ref.Drms = refDrms;
        out.ref.J    = JRef;

        out.plotData = struct();

        out.plotData.period = struct();
        out.plotData.period.tau = tauGrid;
        out.plotData.period.J   = JGrid;

        out.plotData.self = struct();
        out.plotData.self.t  = tS;
        out.plotData.self.u  = uS;
        out.plotData.self.v  = vS;
        out.plotData.self.us = uShiftS;
        out.plotData.self.vs = vShiftS;
        out.plotData.self.d  = dSelf;

        out.plotData.ref = struct();
        out.plotData.ref.phi = phiGrid;
        out.plotData.ref.J   = JR;

        out.plotData.orbitDistHistory = orbitHist;
        out.plotData.conv = conv;

        % --------------------------------------------------------
        % Nested helpers
        % --------------------------------------------------------
        function J = shift_objective(tau)
            mask = (t <= t(end) - tau);

            t0 = t(mask);
            u0 = u(mask);
            v0 = v(mask);

            tq = t0 + tau;

            uq = interp1(t, u, tq, opts.interpMethod);
            vq = interp1(t, v, tq, opts.interpMethod);

            du = (uq - u0) / Su;
            dv = (vq - v0) / Sv;

            J = sqrt(mean(du.^2 + dv.^2));
        end

        function [t0, u0, v0, uq, vq] = build_shifted_overlap_dense(tau)
            mask = (td <= td(end) - tau);

            t0 = td(mask);
            u0 = ud(mask);
            v0 = vd(mask);

            tq = t0 + tau;

            uq = interp1(td, ud, tq, opts.interpMethod);
            vq = interp1(td, vd, tq, opts.interpMethod);
        end

        function J = ref_objective(phi)
            trel = td - td(1);
            [ur, vr] = eval_ref_time(trel + phi);

            du = (ur - ud) / Sud;
            dv = (vr - vd) / Svd;

            J = sqrt(mean(du.^2 + dv.^2));
        end

        function [u0, v0, ur, vr] = build_ref_fit(phi)
            u0 = ud;
            v0 = vd;

            trel = td - td(1);
            [ur, vr] = eval_ref_time(trel + phi);
        end

        function [ur, vr] = eval_ref_time(tq)
            s = mod(tq(:), Tref) / Tref;

            sRef = ref.s(:);
            uRef = ref.u(:);
            vRef = ref.v(:);

            if abs(sRef(end) - 1) > 1e-12
                sRef = [sRef; 1];
                uRef = [uRef; uRef(1)];
                vRef = [vRef; vRef(1)];
            end

            ur = interp1(sRef, uRef, s, opts.interpMethod);
            vr = interp1(sRef, vRef, s, opts.interpMethod);
        end

        function convLoc = compute_convergence_from_history(tHist, dHist, Tper)
            convLoc = empty_conv_struct();

            if numel(tHist) < 2 || ~isfinite(Tper) || Tper <= 0
                return;
            end

            t0 = tHist(1);
            tEnd = tHist(end);

            nWin = floor((tEnd - t0) / Tper);
            if nWin < 1
                return;
            end

            pp = pchip(tHist, dHist);

            tMin = NaN(nWin, 1);
            dMin = NaN(nWin, 1);
            tMinInterp = NaN(nWin, 1);
            dMinInterp = NaN(nWin, 1);

            for q = 1:nWin
                w0 = t0 + (q - 1) * Tper;
                w1 = t0 + q * Tper;

                if q < nWin
                    mask = (tHist >= w0) & (tHist < w1);
                else
                    mask = (tHist >= w0) & (tHist <= w1);
                end

                if ~any(mask)
                    continue;
                end

                tw = tHist(mask);
                dw = dHist(mask);

                % Sampled period-window minimum
                [dMin(q), idxLoc] = min(dw);
                tMin(q) = tw(idxLoc);

                % Interpolated period-window minimum
                if numel(tw) == 1 || w1 <= w0
                    tMinInterp(q) = tMin(q);
                    dMinInterp(q) = dMin(q);
                else
                    fwin = @(tt) ppval(pp, tt);
                    tStar = fminbnd(fwin, w0, w1);
                    dStar = fwin(tStar);

                    % Cheap safeguard: do not let the interpolated minimum
                    % exceed the sampled minimum in the same window.
                    dStar = min(dStar, dMin(q));

                    % Guard against tiny negative interpolation artefacts
                    dStar = max(dStar, 0);

                    tMinInterp(q) = tStar;
                    dMinInterp(q) = dStar;
                end
            end

            valid = ~isnan(tMin) & ~isnan(dMin) & ~isnan(tMinInterp) & ~isnan(dMinInterp);

            tMin = tMin(valid);
            dMin = dMin(valid);
            tMinInterp = tMinInterp(valid);
            dMinInterp = dMinInterp(valid);

            convLoc.tMin = tMin;
            convLoc.dMin = dMin;
            convLoc.tMinInterp = tMinInterp;
            convLoc.dMinInterp = dMinInterp;

            if opts.convUseInterp
                tUse = tMinInterp;
                dUse = dMinInterp;
            else
                tUse = tMin;
                dUse = dMin;
            end

            nBlock = opts.convNDips;
            convLoc.nDips = nBlock;

            % Need two nonoverlapping blocks of length nBlock
            if numel(dUse) < 2 * nBlock
                return;
            end

            nPairs = numel(dUse) - 2 * nBlock + 1;
            tBlock    = NaN(nPairs, 1);
            varBlock  = NaN(nPairs, 1);
            spanBlock = NaN(nPairs, 1);
            tolBlock  = NaN(nPairs, 1);
            passRelTolBlock = NaN(nPairs, 1);

            for p = 1:nPairs
                prevD = dUse(p : p + nBlock - 1);
                currD = dUse(p + nBlock : p + 2*nBlock - 1);

                prevT = tUse(p : p + nBlock - 1);
                currT = tUse(p + nBlock : p + 2*nBlock - 1);

                avgPrev = mean(prevD);
                avgCurr = mean(currD);

                tBlock(p)    = currT(end);
                varBlock(p)  = abs(avgCurr - avgPrev);
                spanBlock(p) = currT(end) - prevT(1);
                tolBlock(p)  = opts.convRelTol * avgCurr;
                passRelTolBlock(p) = required_passed_reltol(varBlock(p), avgCurr);

                if ~convLoc.ok && (varBlock(p) <= tolBlock(p))
                    convLoc.ok    = true;
                    convLoc.tConv = currT(end);
                    convLoc.d     = avgCurr;
                    convLoc.idx   = p;
                end
            end

            finiteRelTol = passRelTolBlock(isfinite(passRelTolBlock));
            if isempty(finiteRelTol)
                convLoc.passedRelTol = inf;
            else
                convLoc.passedRelTol = min(finiteRelTol);
            end
            convLoc.tBlock    = tBlock;
            convLoc.varBlock  = varBlock;
            convLoc.spanBlock = spanBlock;
            convLoc.tolBlock  = tolBlock;
        end
        
        function r = required_passed_reltol(varVal, avgVal)
            if avgVal > 0
                r = varVal / avgVal;
            elseif varVal == 0
                r = 0;
            else
                r = inf;
            end
        end

        function convLoc = empty_conv_struct()
            convLoc = struct();
            convLoc.ok          = false;
            convLoc.tConv       = [];
            convLoc.d           = [];
            convLoc.idx         = [];
            convLoc.nDips       = opts.convNDips;
            convLoc.relTol      = opts.convRelTol;
            convLoc.passedRelTol = inf;
            convLoc.tMin        = [];
            convLoc.dMin        = [];
            convLoc.tMinInterp  = [];
            convLoc.dMinInterp  = [];
            convLoc.tBlock      = [];
            convLoc.varBlock    = [];
            convLoc.spanBlock   = [];
            convLoc.tolBlock    = [];
        end
    end
end

% ==============================================================
% Small local utilities
% ==============================================================
function val = get_field(S, parent, name, defaultVal)
    if isempty(parent)
        if isfield(S, name)
            val = S.(name);
        else
            val = defaultVal;
        end
    else
        if isfield(S, parent) && isfield(S.(parent), name)
            val = S.(parent).(name);
        else
            val = defaultVal;
        end
    end
end

function label = get_ic_label(sim)
    label = '';
    if isfield(sim, 'icInfo') && isstruct(sim.icInfo) && isfield(sim.icInfo, 'label')
        label = sim.icInfo.label;
    end
end