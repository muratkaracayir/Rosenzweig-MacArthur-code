function [eigOut, summary, qc, info] = floquet_processM(M, T, cfgFloq, UV, cfgOrbit)
%FLOQUET_PROCESSM  Postprocess an explicit Floquet monodromy matrix.
%
%   [eigOut, summary, qc, info] = floquet_processM(M, T, cfgFloq, UV, cfgOrbit)
%
%   This routine computes the full eigendecomposition of the explicit
%   monodromy matrix M, identifies the trivial Floquet multiplier by
%   alignment with the time-shift direction q, stores the leading
%   multipliers by modulus, and returns basic stability and quality-control
%   diagnostics.
%
%   Inputs
%   ------
%   M        : square monodromy matrix.
%   T        : orbit period. Used here only for validation/provenance.
%   cfgFloq  : Floquet configuration produced by make_floqCfg.
%   UV       : packed orbit [U(:); V(:)] on an Nx-by-Nt periodic time grid.
%   cfgOrbit : orbit/model configuration, usually hopf.cfg.
%
%   Outputs
%   -------
%   eigOut  : stored multipliers, optional eigenvectors, and residuals.
%   summary : stability and neutral-direction diagnostics.
%   qc      : quality-control flags and thresholds.
%   info    : lightweight provenance information.
%
%   MATLAB R2016a compatible.

    validate_inputs(M, T, cfgFloq, UV, cfgOrbit);

    n = size(M, 1);

    % Neutral time-shift direction at the first time slice.
    q = compute_neutral_direction(UV, cfgOrbit);

    % Full eigendecomposition.
    tEig = tic;
    [V, D] = eig(M);
    eigTime = toc(tEig);
    mu = diag(D);

    % Eigenpair residuals for the full spectrum.
    resAll = zeros(n, 1);
    for j = 1:n
        x = V(:, j);
        nx = norm(x);
        if nx == 0
            resAll(j) = Inf;
        else
            resAll(j) = norm(M*x - mu(j)*x) / nx;
        end
    end

    % Normalize eigenvectors before q-alignment and optional storage.
    for j = 1:n
        nj = norm(V(:, j));
        if nj > 0
            V(:, j) = V(:, j) / nj;
        end
    end

    % Identify the trivial multiplier by maximum alignment with q.
    align = abs(V' * q);
    [~, idxTrivial] = max(align);
    muTrivial = mu(idxTrivial);
    muTrivialDist = abs(muTrivial - 1);

    % Also record the multiplier closest to one, for ambiguity diagnostics.
    [~, idxNearOne] = min(abs(mu - 1));
    muNearOne = mu(idxNearOne);

    % Neutral-direction diagnostics.
    Mq = M * q;
    eta_q = norm(Mq - q) / norm(q);
    mu_q = (q' * Mq) / (q' * q);
    r_q = norm(Mq - mu_q*q) / norm(q);
    r_triv = norm(Mq - muTrivial*q) / norm(q);

    % Store the leading multipliers by modulus, always including trivial.
    muAbs = abs(mu);
    [~, permAbsDesc] = sort(muAbs, 'descend');

    nStore = min(cfgFloq.nStore, n);
    idxStored = permAbsDesc(1:nStore);
    if ~any(idxStored == idxTrivial)
        idxStored = [idxStored(:); idxTrivial]; %#ok<AGROW>
    end

    muStored = mu(idxStored);
    resStored = resAll(idxStored);

    if cfgFloq.storeEigenvectors
        XStored = V(:, idxStored);
    else
        XStored = [];
    end

    posTrivialStored = find(idxStored == idxTrivial, 1, 'first');

    % Stability counts from the full spectrum, excluding the trivial pair.
    idxNontriv = true(n, 1);
    idxNontriv(idxTrivial) = false;

    muNon = mu(idxNontriv);
    muNonAbs = abs(muNon);

    if isempty(muNon)
        nUnstable = 0;
        gap = NaN;
        muGap = NaN;
    else
        nUnstable = sum(muNonAbs > (1 + cfgFloq.muTol));
        [gap, idxGap] = min(abs(muNonAbs - 1));
        muGap = muNon(idxGap);
    end

    maxResAll = max(resAll);

    % Quality-control flags.
    reasons = {};

    if muTrivialDist > cfgFloq.maxMuTrivialDist
        reasons{end+1} = 'trivialNotNear1'; %#ok<AGROW>
    end

    if ~isempty(muNon) && gap < cfgFloq.near1Tol
        reasons{end+1} = 'nearUnitCircle'; %#ok<AGROW>
    end

    if maxResAll > cfgFloq.eigResTol
        reasons{end+1} = 'eigResidualLarge'; %#ok<AGROW>
    end

    if abs(muTrivial - muNearOne) > cfgFloq.trivialAgreeTol
        reasons{end+1} = 'trivialAmbiguous'; %#ok<AGROW>
    end

    if r_triv > cfgFloq.rTrivTol
        reasons{end+1} = 'neutralDefectLarge'; %#ok<AGROW>
    end

    % Output structures.
    qc = struct();
    qc.reason = reasons;
    qc.ok = isempty(reasons);
    qc.needsRefine = ~isempty(reasons);
    qc.thresholds = struct( ...
        'muTol', cfgFloq.muTol, ...
        'near1Tol', cfgFloq.near1Tol, ...
        'maxMuTrivialDist', cfgFloq.maxMuTrivialDist, ...
        'eigResTol', cfgFloq.eigResTol, ...
        'trivialAgreeTol', cfgFloq.trivialAgreeTol, ...
        'rTrivTol', cfgFloq.rTrivTol);

    eigOut = struct();
    eigOut.mu = muStored;
    eigOut.abs = abs(muStored);
    eigOut.res = resStored;
    eigOut.X = XStored;
    eigOut.posTrivial = posTrivialStored;
    eigOut.muTrivial = muTrivial;
    eigOut.maxEigRes = maxResAll;

    summary = struct();
    summary.muTrivialDist = muTrivialDist;
    summary.nUnstable = nUnstable;
    summary.gapUnitCircle = gap;
    summary.muGap = muGap;
    summary.eta_q = eta_q;
    summary.mu_q = mu_q;
    summary.r_q = r_q;
    summary.r_triv = r_triv;

    info = struct();
    info.n = n;
    info.period = T;
    info.eigTime_sec = eigTime;
    info.q = struct('method', 'RHS-at-t0', 'tIndex', 1);

end

function validate_inputs(M, T, cfgFloq, UV, cfgOrbit)
%VALIDATE_INPUTS  Lightweight validation for Floquet postprocessing.

    if ~ismatrix(M) || size(M, 1) ~= size(M, 2)
        error('floquet_processM:BadM', 'M must be a square matrix.');
    end

    if ~isscalar(T) || ~isfinite(T) || T <= 0
        error('floquet_processM:BadPeriod', 'T must be a finite positive scalar.');
    end

    if ~isstruct(cfgFloq)
        error('floquet_processM:BadCfgFloq', 'cfgFloq must be a struct.');
    end

    requiredFloqFields = { ...
        'nStore', 'storeEigenvectors', 'muTol', 'near1Tol', ...
        'maxMuTrivialDist', 'eigResTol', 'trivialAgreeTol', 'rTrivTol'};
    for j = 1:numel(requiredFloqFields)
        field = requiredFloqFields{j};
        if ~isfield(cfgFloq, field)
            error('floquet_processM:MissingFloqField', ...
                'cfgFloq.%s is required.', field);
        end
    end

    if cfgFloq.nStore <= 0 || fix(cfgFloq.nStore) ~= cfgFloq.nStore
        error('floquet_processM:BadNStore', ...
            'cfgFloq.nStore must be a positive integer.');
    end

    if ~isvector(UV) || isempty(UV)
        error('floquet_processM:BadUV', 'UV must be a nonempty vector.');
    end

    if ~isstruct(cfgOrbit)
        error('floquet_processM:BadCfgOrbit', 'cfgOrbit must be a struct.');
    end

end

function q = compute_neutral_direction(UV, cfgOrbit)
%COMPUTE_NEUTRAL_DIRECTION  Compute q = dU/dt at the first time slice.
%
%   The neutral direction is computed from the semi-discrete PDE right-hand
%   side at t=0, using the same DCT-I Neumann Laplacian convention as the
%   periodic-orbit residual.

    require_orbit_fields(cfgOrbit);

    Nx = cfgOrbit.periodic_orbit.Nx;
    Nt = cfgOrbit.periodic_orbit.Nt;
    Ntot = Nx * Nt;

    if numel(UV) ~= 2*Ntot
        error('floquet_processM:BadUVLength', ...
            'UV must have length 2*Nx*Nt.');
    end

    if isfield(cfgOrbit, 'L') && ~isempty(cfgOrbit.L)
        Lx = cfgOrbit.L;
    elseif isfield(cfgOrbit, 'ell') && ~isempty(cfgOrbit.ell)
        Lx = cfgOrbit.ell * pi;
    else
        error('floquet_processM:MissingDomainLength', ...
            'cfgOrbit must define either L or ell.');
    end

    U = reshape(UV(1:Ntot), Nx, Nt);
    V = reshape(UV(Ntot+1:end), Nx, Nt);

    u0 = U(:, 1);
    v0 = V(:, 1);

    uxx0 = laplacian_vec_dct1(u0, Lx);
    vxx0 = laplacian_vec_dct1(v0, Lx);

    m = cfgOrbit.m;
    K = cfgOrbit.k;
    theta = cfgOrbit.theta;

    pred = (m .* u0 .* v0) ./ (u0 + 1);

    f = u0 .* (1 - u0./K) - pred;
    g = -theta .* v0 + pred;

    ut0 = cfgOrbit.d1 .* uxx0 + f;
    vt0 = cfgOrbit.d2 .* vxx0 + g;

    q = [ut0; vt0];
    nq = norm(q);
    if nq == 0
        error('floquet_processM:ZeroNeutralDirection', ...
            'The computed neutral direction has zero norm.');
    end
    q = q / nq;

end

function require_orbit_fields(cfgOrbit)
%REQUIRE_ORBIT_FIELDS  Check fields needed to compute the neutral direction.

    if ~isfield(cfgOrbit, 'periodic_orbit') || ...
            ~isfield(cfgOrbit.periodic_orbit, 'Nx') || ...
            ~isfield(cfgOrbit.periodic_orbit, 'Nt')
        error('floquet_processM:MissingGridFields', ...
            'cfgOrbit.periodic_orbit.Nx and Nt are required.');
    end

    fields = {'d1', 'd2', 'm', 'theta', 'k'};
    for j = 1:numel(fields)
        field = fields{j};
        if ~isfield(cfgOrbit, field) || isempty(cfgOrbit.(field))
            error('floquet_processM:MissingModelField', ...
                'cfgOrbit.%s is required.', field);
        end
    end

end

function uxx = laplacian_vec_dct1(u, Lx)
%LAPLACIAN_VEC_DCT1  DCT-I Neumann second derivative on endpoint grid.

    u = u(:);
    Nx = length(u);
    if Nx == 1
        uxx = 0*u;
        return;
    end

    c = dct1_fast(u);
    k = (0:Nx-1)';
    lambda = -(pi*k/Lx).^2;
    uxx = dct1_fast(c .* lambda) * (Nx-1)/2;

end

function a = dct1_fast(u)
%DCT1_FAST  Direct DCT-I helper matching the periodic-orbit residual.

    u = u(:);
    Nx = length(u);
    if Nx == 1
        a = u;
        return;
    end

    k = (0:Nx-1)';
    n = 1:Nx-2;
    C = cos(pi * k * n / (Nx-1));
    a = 0.5*u(1) + 0.5*u(end) * (-1).^k + C * u(2:end-1);
    a = a * 2/(Nx-1);

end
