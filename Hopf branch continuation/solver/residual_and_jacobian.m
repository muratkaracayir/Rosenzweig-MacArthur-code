function [R, J] = residual_and_jacobian(UV, hopf)
% RESIDUAL_AND_JACOBIAN  Residual and Jacobian for periodic-orbit Newton solves.
%
%   R = RESIDUAL_AND_JACOBIAN(UV, hopf) returns the periodic-orbit
%   residual with one phase condition appended.
%
%   [R, J] = RESIDUAL_AND_JACOBIAN(UV, hopf) also returns the Jacobian with
%   respect to the unknowns [UV; T].  The orbit is represented on a
%   cosine-spectral grid in space and a Fourier grid in time.

    require_hopf_fields(hopf);

    ops = setup_operators(hopf);

    [u, v] = reshape_solution(UV, ops.Nx, ops.Nt);
    ut  = time_derivative(u, ops);
    vt  = time_derivative(v, ops);
    uxx = laplacian(u, ops.Lx);
    vxx = laplacian(v, ops.Lx);

    [f, g, dfdu, dfdv, dgdu, dgdv] = reaction_terms(u, v, hopf);

    d1 = hopf.cfg.d1;
    d2 = hopf.cfg.d2;

    Ru = ut - d1 * uxx - f;
    Rv = vt - d2 * vxx - g;
    R = [Ru(:); Rv(:)];

    if nargout == 1
        R = add_phase_condition(R, [], UV, hopf, ops);
        return;
    end

    J = assemble_jacobian(dfdu, dfdv, dgdu, dgdv, hopf, ops);

    JT = period_column(u, v, hopf, ops);
    J = [J, JT];

    [R, J] = add_phase_condition(R, J, UV, hopf, ops);
end

% -------------------------------------------------------------------------
function require_hopf_fields(hopf)
    if ~isfield(hopf, 'cfg')
        error('residual_and_jacobian:MissingConfig', ...
              'The input hopf must contain hopf.cfg.');
    end

    if ~isfield(hopf, 'T') || isempty(hopf.T)
        error('residual_and_jacobian:MissingPeriod', ...
              'The input hopf must contain the current period hopf.T.');
    end

    cfg = hopf.cfg;
    if ~isfield(cfg, 'periodic_orbit') || ...
       ~isfield(cfg.periodic_orbit, 'Nx') || ...
       ~isfield(cfg.periodic_orbit, 'Nt')
        error('residual_and_jacobian:MissingGrid', ...
              'The configuration must define cfg.periodic_orbit.Nx and Nt.');
    end
end

% -------------------------------------------------------------------------
function ops = setup_operators(hopf)
    cfg = hopf.cfg;

    Nx = cfg.periodic_orbit.Nx;
    Nt = cfg.periodic_orbit.Nt;
    Lx = cfg.L;
    T  = hopf.T;

    x = linspace(0, Lx, Nx)';
    t = linspace(0, T, Nt + 1);
    t(end) = [];

    mu_x = -(((0:Nx-1)' * pi / Lx).^2);

    k = 0:(Nt - 1);
    half = floor(Nt / 2);
    k_signed = k;
    k_signed(k > half) = k(k > half) - Nt;
    omega_k = (2*pi / T) * k_signed;

    ops.Nx = Nx;
    ops.Nt = Nt;
    ops.Lx = Lx;
    ops.T  = T;
    ops.x  = x;
    ops.t  = t(:);
    ops.mu_x = mu_x;
    ops.omega_k = omega_k;

    Dxx = zeros(Nx, Nx);
    for j = 1:Nx
        e = zeros(Nx, 1);
        e(j) = 1;
        Dxx(:, j) = laplacian(e, Lx);
    end

    F = fft(eye(Nt));
    Finv = ifft(eye(Nt));
    Dt = real(Finv * (1i * diag(omega_k)) * F);

    ops.Dxx_sp = sparse(Dxx);
    ops.Dt_sp  = sparse(Dt);
end

% -------------------------------------------------------------------------
function uxx = laplacian(u, Lx)
    [Nx, Nt] = size(u);
    uxx = zeros(Nx, Nt);

    k = (0:Nx-1)';
    lambda = -(pi * k / Lx).^2;

    for j = 1:Nt
        c = dct1_fast(u(:, j));
        c2 = c .* lambda;
        uxx(:, j) = dct1_fast(c2) * (Nx - 1) / 2;
    end
end

% -------------------------------------------------------------------------
function a = dct1_fast(u)
    u = u(:);
    Nx = length(u);

    if Nx == 1
        a = u;
        return;
    end

    k = (0:Nx-1)';
    n = 1:(Nx-2);
    C = cos(pi * k * n / (Nx - 1));

    a = 0.5 * u(1) + 0.5 * u(end) * (-1).^k + C * u(2:end-1);
    a = a * 2 / (Nx - 1);
end

% -------------------------------------------------------------------------
function ut = time_derivative(u, ops)
    ut = real(ifft(bsxfun(@times, 1i * ops.omega_k, fft(u, [], 2)), [], 2));
end

% -------------------------------------------------------------------------
function [f, g, dfdu, dfdv, dgdu, dgdv] = reaction_terms(u, v, hopf)
    cfg = hopf.cfg;
    m = cfg.m;
    k = cfg.k;
    theta = cfg.theta;

    denom = u + 1;

    f = u .* (1 - u / k) - m * u .* v ./ denom;
    g = -theta * v + m * u .* v ./ denom;

    dfdu = 1 - 2 * u / k - (m * v) ./ denom + ...
           (m * u .* v) ./ denom.^2;
    dfdv = -m * u ./ denom;

    dgdu = (m * v) ./ denom - (m * u .* v) ./ denom.^2;
    dgdv = -theta + (m * u) ./ denom;
end

% -------------------------------------------------------------------------
function J = assemble_jacobian(dfdu, dfdv, dgdu, dgdv, hopf, ops)
    Nx = ops.Nx;
    Nt = ops.Nt;
    Ntot = Nx * Nt;

    Ix = speye(Nx);
    It = speye(Nt);

    Jtime = kron(ops.Dt_sp, Ix);
    Ju_space = -hopf.cfg.d1 * kron(It, ops.Dxx_sp);
    Jv_space = -hopf.cfg.d2 * kron(It, ops.Dxx_sp);

    Juu = Jtime + Ju_space - spdiags(dfdu(:), 0, Ntot, Ntot);
    Juv =                    - spdiags(dfdv(:), 0, Ntot, Ntot);
    Jvu =                    - spdiags(dgdu(:), 0, Ntot, Ntot);
    Jvv = Jtime + Jv_space - spdiags(dgdv(:), 0, Ntot, Ntot);

    J = [Juu, Juv;
         Jvu, Jvv];
end

% -------------------------------------------------------------------------
function JT = period_column(u, v, hopf, ops)
    T = hopf.T;
    domega_dT = -ops.omega_k / T;

    dUt_dT = real(ifft(bsxfun(@times, 1i * domega_dT, fft(u, [], 2)), [], 2));
    dVt_dT = real(ifft(bsxfun(@times, 1i * domega_dT, fft(v, [], 2)), [], 2));

    JT = [dUt_dT(:); dVt_dT(:)];
end

% -------------------------------------------------------------------------
function [R, J] = add_phase_condition(R, J, UV, hopf, ops)
    Nx = ops.Nx;
    Nt = ops.Nt;
    Ntot = Nx * Nt;

    [u_ref, v_ref] = phase_reference(hopf, ops);
    u_ref_t = time_derivative(u_ref, ops);
    v_ref_t = time_derivative(v_ref, ops);

    [u, v] = reshape_solution(UV, Nx, Nt);

    Phi = sum((u(:) - u_ref(:)) .* u_ref_t(:)) + ...
          sum((v(:) - v_ref(:)) .* v_ref_t(:));
    R = [R; Phi];

    if isempty(J)
        J = [];
        return;
    end

    phase_row = [u_ref_t(:).', v_ref_t(:).'];
    Phi_T = 0;

    J = [J; phase_row, Phi_T];
end

% -------------------------------------------------------------------------
function [u_ref, v_ref] = phase_reference(hopf, ops)
    Nx = ops.Nx;
    Nt = ops.Nt;
    Ntot = Nx * Nt;

    if isfield(hopf, 'phaseRefUV') && ~isempty(hopf.phaseRefUV)
        UVref = hopf.phaseRefUV;

        if numel(UVref) ~= 2 * Ntot
            error('residual_and_jacobian:PhaseReferenceSize', ...
                  'phaseRefUV length mismatch: expected %d, got %d.', ...
                  2*Ntot, numel(UVref));
        end

        u_ref = reshape(UVref(1:Ntot), Nx, Nt);
        v_ref = reshape(UVref(Ntot+1:end), Nx, Nt);
        return;
    end

    if ~isfield(hopf, 'psi') || ~isfield(hopf, 'omega') || ~isfield(hopf, 'mode')
        error('residual_and_jacobian:MissingPhaseReference', ...
              ['The input hopf must contain either hopf.phaseRefUV, or ', ...
               'the Hopf reference fields hopf.psi, hopf.omega, and hopf.mode.']);
    end

    psi = hopf.psi;
    omega = hopf.omega;
    t = ops.t(:);

    phi_x = cos(hopf.mode * pi * ops.x / ops.Lx);

    time_u = real(psi(1) * exp(1i * omega * t'));
    time_v = real(psi(2) * exp(1i * omega * t'));

    u_ref = phi_x * time_u;
    v_ref = phi_x * time_v;
end

% -------------------------------------------------------------------------
function [u, v] = reshape_solution(UV, Nx, Nt)
    Ntot = Nx * Nt;

    if numel(UV) ~= 2 * Ntot
        error('residual_and_jacobian:SolutionSize', ...
              'UV length mismatch: expected %d, got %d.', 2*Ntot, numel(UV));
    end

    u = reshape(UV(1:Ntot), Nx, Nt);
    v = reshape(UV(Ntot+1:end), Nx, Nt);
end
