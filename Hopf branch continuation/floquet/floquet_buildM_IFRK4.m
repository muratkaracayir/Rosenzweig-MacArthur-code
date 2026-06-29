function out = floquet_buildM_IFRK4(UV, T, cfg, opts)
%FLOQUET_BUILDM_IFRK4  Build the monodromy matrix with an IFRK4 scheme.
%
%   out = floquet_buildM_IFRK4(UV, T, cfg, opts)
%
%   This routine builds the explicit monodromy matrix for the variational
%   equation along one stored periodic orbit. The diffusion part is treated
%   by an integrating factor based on the same DCT-I Neumann Laplacian used
%   in the periodic-orbit residual, while the reaction Jacobian is advanced
%   with classical fourth-order Runge--Kutta stages.
%
%   Inputs
%   ------
%   UV   : packed orbit vector of length 2*Nx*Nt,
%          [u(:,1); ...; u(:,Nt); v(:,1); ...; v(:,Nt)].
%          The stored time grid excludes the endpoint t = T.
%   T    : orbit period.
%   cfg  : configuration struct containing d1, d2, k, m, theta, and either
%          cfg.L or cfg.ell. Grid sizes are read from cfg.periodic_orbit
%          when available, otherwise from cfg.Nx and cfg.Nt.
%   opts : optional struct with fields
%          mSub    : substeps per stored time interval, default 16;
%          verbose : print progress messages, default true.
%
%   Output
%   ------
%   out.M    : explicit monodromy matrix of size 2*Nx by 2*Nx.
%   out.diag : diagnostics, including time-step information and the neutral
%              direction residual ||Mq-q||/||q||.
%
%   MATLAB R2016a compatible.

    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    opts = local_complete_opts(opts);

    [Nx, Nt] = local_get_NxNt(cfg);
    local_validate_inputs(UV, T, cfg, opts, Nx, Nt);

    n = 2 * Nx;

    dtStore = T / Nt;
    mSub    = opts.mSub;
    dt      = dtStore / mSub;
    nSteps  = Nt * mSub;

    evalOrbit = floquet_orbitEval_trig(UV, T, cfg);

    Lx   = local_get_Lx(cfg);
    diff = local_make_diffusion_cache(Nx, Lx);

    M = zeros(n, n);
    tStart = tic;

    for j = 1:n
        z = zeros(n, 1);
        z(j) = 1;
        t = 0;

        for step = 1:nSteps
            z = local_ifrk4_step(z, t, dt, evalOrbit, cfg, Nx, diff);
            t = t + dt;
            if t >= T
                t = t - T;
            end
        end

        M(:, j) = z;

        if opts.verbose && (mod(j, 32) == 0 || j == n)
            fprintf('  built %d/%d columns\n', j, n);
        end
    end

    elapsed = toc(tStart);

    q = local_compute_q_from_UV(UV, T, cfg);
    Mq = M * q;
    eta_q = norm(Mq - q) / norm(q);
    mu_q  = (q' * Mq) / (q' * q);

    out = struct();
    out.M = M;
    out.diag = struct( ...
        'mSub', mSub, ...
        'dt', dt, ...
        'nSteps', nSteps, ...
        'elapsed_sec', elapsed, ...
        'eta_q', eta_q, ...
        'mu_q', mu_q);

    if opts.verbose
        fprintf('--- floquet_buildM_IFRK4 done ---\n');
        fprintf('eta_q = ||Mq-q||/||q|| = %.3e\n', eta_q);
        fprintf('Rayleigh(mu along q) = %.16g%+.3gi\n', real(mu_q), imag(mu_q));
    end
end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function opts = local_complete_opts(opts)
    if ~isfield(opts, 'mSub') || isempty(opts.mSub)
        opts.mSub = 16;
    end
    if ~isfield(opts, 'verbose') || isempty(opts.verbose)
        opts.verbose = true;
    end
end

function local_validate_inputs(UV, T, cfg, opts, Nx, Nt)
    if numel(UV) ~= 2 * Nx * Nt
        error('floquet_buildM_IFRK4:BadUVSize', ...
            'UV must have length 2*Nx*Nt.');
    end
    if ~isscalar(T) || T <= 0
        error('floquet_buildM_IFRK4:BadPeriod', ...
            'T must be a positive scalar.');
    end
    if opts.mSub < 1 || opts.mSub ~= round(opts.mSub)
        error('floquet_buildM_IFRK4:BadSubstepCount', ...
            'opts.mSub must be a positive integer.');
    end

    required = {'d1', 'd2', 'k', 'm', 'theta'};
    for i = 1:numel(required)
        if ~isfield(cfg, required{i})
            error('floquet_buildM_IFRK4:MissingConfigField', ...
                'cfg.%s is required.', required{i});
        end
    end

    if ~(isfield(cfg, 'L') && ~isempty(cfg.L)) && ~isfield(cfg, 'ell')
        error('floquet_buildM_IFRK4:MissingDomainLength', ...
            'The configuration must define cfg.L or cfg.ell.');
    end
end

function zNext = local_ifrk4_step(z, t0, dt, evalOrbit, cfg, Nx, diff)
% One integrating-factor RK4 step for the variational equation.

    w0  = z;
    k1w = local_Jreact_times_w(t0, w0, evalOrbit, cfg, Nx);
    k1  = k1w;

    z2  = z + (dt / 2) * k1;
    w2  = local_apply_diff_exp_block(z2, diff, cfg, +dt / 2, Nx);
    k2w = local_Jreact_times_w(t0 + 0.5 * dt, w2, evalOrbit, cfg, Nx);
    k2  = local_apply_diff_exp_block(k2w, diff, cfg, -dt / 2, Nx);

    z3  = z + (dt / 2) * k2;
    w3  = local_apply_diff_exp_block(z3, diff, cfg, +dt / 2, Nx);
    k3w = local_Jreact_times_w(t0 + 0.5 * dt, w3, evalOrbit, cfg, Nx);
    k3  = local_apply_diff_exp_block(k3w, diff, cfg, -dt / 2, Nx);

    z4  = z + dt * k3;
    w4  = local_apply_diff_exp_block(z4, diff, cfg, +dt, Nx);
    k4w = local_Jreact_times_w(t0 + dt, w4, evalOrbit, cfg, Nx);
    k4  = local_apply_diff_exp_block(k4w, diff, cfg, -dt, Nx);

    zNew = z + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    zNext = local_apply_diff_exp_block(zNew, diff, cfg, +dt, Nx);
end

function Jw = local_Jreact_times_w(t, w, evalOrbit, cfg, Nx)
% Apply the pointwise reaction Jacobian along the orbit to perturbation w.

    du = w(1:Nx);
    dv = w(Nx+1:end);

    [u, v] = evalOrbit(t);

    m = cfg.m;
    k = cfg.k;
    theta = cfg.theta;

    dfdu = 1 - 2 * u / k - (m * v) ./ (u + 1) + (m * u .* v) ./ (u + 1).^2;
    dfdv = -m * u ./ (u + 1);
    dgdu = (m * v) ./ (u + 1) - (m * u .* v) ./ (u + 1).^2;
    dgdv = -theta + (m * u) ./ (u + 1);

    du_t = dfdu .* du + dfdv .* dv;
    dv_t = dgdu .* du + dgdv .* dv;

    Jw = [du_t; dv_t];
end

function w = local_apply_diff_exp_block(w, diff, cfg, tau, Nx)
% Apply the block diffusion exponential to the prey and predator components.

    du = w(1:Nx);
    dv = w(Nx+1:end);

    du = local_apply_diff_exp_scalar(du, diff, cfg.d1, tau);
    dv = local_apply_diff_exp_scalar(dv, diff, cfg.d2, tau);

    w = [du; dv];
end

function y = local_apply_diff_exp_scalar(y, diff, d, tau)
% Apply exp(d*Dxx*tau) using the cached eigendecomposition of Dxx.

    c = diff.iV * y;
    c = exp(d * diff.lam * tau) .* c;
    y = diff.V * c;
end

function diff = local_make_diffusion_cache(Nx, Lx)
% Build the Neumann DCT-I Laplacian and cache its eigendecomposition.

    Dxx = zeros(Nx, Nx);
    for j = 1:Nx
        e = zeros(Nx, 1);
        e(j) = 1;
        Dxx(:, j) = local_laplacian_vec_dct1(e, Lx);
    end

    [V, Lam] = eig(Dxx);
    diff.V   = V;
    diff.iV  = inv(V);
    diff.lam = diag(Lam);
end

function uxx = local_laplacian_vec_dct1(u, Lx)
% Apply the same DCT-I Neumann Laplacian used in the orbit residual.

    u = u(:);
    Nx = length(u);
    if Nx == 1
        uxx = 0 * u;
        return
    end

    c = local_dct1_fast(u);
    k = (0:Nx-1)';
    lambda = -(pi * k / Lx).^2;
    c2 = c .* lambda;
    uxx = local_dct1_fast(c2) * (Nx - 1) / 2;
end

function a = local_dct1_fast(u)
% DCT-I transform with the scaling convention used by the residual code.

    u = u(:);
    Nx = length(u);
    if Nx == 1
        a = u;
        return
    end

    k = (0:Nx-1)';
    n = 1:Nx-2;
    C = cos(pi * k * n / (Nx - 1));
    a = 0.5 * u(1) + 0.5 * u(end) * (-1).^k + C * u(2:end-1);
    a = a * 2 / (Nx - 1);
end

function [Nx, Nt] = local_get_NxNt(cfg)
    if isfield(cfg, 'periodic_orbit') && isstruct(cfg.periodic_orbit) ...
            && isfield(cfg.periodic_orbit, 'Nx') && isfield(cfg.periodic_orbit, 'Nt')
        Nx = cfg.periodic_orbit.Nx;
        Nt = cfg.periodic_orbit.Nt;
    elseif isfield(cfg, 'Nx') && isfield(cfg, 'Nt')
        Nx = cfg.Nx;
        Nt = cfg.Nt;
    else
        error('floquet_buildM_IFRK4:MissingGridSize', ...
            'The configuration must define Nx and Nt.');
    end
end

function Lx = local_get_Lx(cfg)
    if isfield(cfg, 'L') && ~isempty(cfg.L)
        Lx = cfg.L;
    else
        Lx = cfg.ell * pi;
    end
end

function q = local_compute_q_from_UV(UV, T, cfg)
% Compute q = dU/dt at t = 0 using Fourier differentiation in time.

    [Nx, Nt] = local_get_NxNt(cfg);
    Ntot = Nx * Nt;

    U = reshape(UV(1:Ntot), Nx, Nt);
    V = reshape(UV(Ntot+1:end), Nx, Nt);

    k = 0:(Nt - 1);
    half = floor(Nt / 2);
    kSigned = k;
    kSigned(k > half) = k(k > half) - Nt;
    omegaK = (2 * pi / T) * kSigned;

    Ut = real(ifft(bsxfun(@times, 1i * omegaK, fft(U, [], 2)), [], 2));
    Vt = real(ifft(bsxfun(@times, 1i * omegaK, fft(V, [], 2)), [], 2));

    q = [Ut(:, 1); Vt(:, 1)];
    q = q / norm(q);
end
