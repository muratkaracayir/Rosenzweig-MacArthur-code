function [init, hopf] = prepare_data(cfg, n)
%PREPARE_DATA  Build a small-amplitude Hopf initial guess.
%
%   [init, hopf] = PREPARE_DATA(cfg, n) computes the mode-n Hopf eigenpair,
%   constructs a small-amplitude periodic-orbit initial guess, updates the
%   parameter value used for Newton correction, and saves the result in data/.
%
%   Inputs:
%       cfg : configuration struct, e.g. from case_ex262.m
%       n   : spatial mode index
%
%   Outputs:
%       init : initial periodic-orbit guess
%       hopf : Hopf/equilibrium data used by the Newton solver

    fprintf('\n=== Preparing initial guess for periodic solver ===\n');

    % Step 1: Compute Hopf eigenpair.
    require_route_fun(cfg);
    hopf = compute_hopf_eigenpair(cfg, n);
    fprintf('Hopf eigenpair computed for mode n = %d.\n', n);

    % Step 2: Construct small-amplitude initial guess.
    init = construct_initial_guess(hopf);
    fprintf('Initial periodic orbit constructed.\n');

    % Step 3: Move slightly away from the Hopf point for Newton correction.
    hopf.cfg = update_cfg_for_periodic_orbit(hopf.cfg);

    % Step 4: Save results.
    thisDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(thisDir);
    dataDir = fullfile(rootDir, 'data');
    if ~exist(dataDir, 'dir')
        mkdir(dataDir);
    end

    fname = fullfile(dataDir, sprintf('initial_guess_%s_mode%d.mat', ...
        cfg.caseName, n));
    save(fname, 'init', 'hopf');

    fprintf('Data saved to %s\n', fname);
    fprintf('Updated parameters: lambda = %.6f, m = %.6f\n', ...
        hopf.cfg.lambda, hopf.cfg.m);
    fprintf('===================================================\n\n');
end

function hopf = compute_hopf_eigenpair(cfg, n)
%COMPUTE_HOPF_EIGENPAIR  Compute the mode-n eigenpair at a Hopf point.

    % Parameters.
    L     = cfg.L;
    d1    = cfg.d1;
    d2    = cfg.d2;
    k     = cfg.k;
    theta = cfg.theta;
    m     = cfg.m;

    % Equilibrium.
    u_star = cfg.u_star;
    v_star = cfg.v_star;

    % Neumann mode cos(n*pi*x/L), with Laplacian eigenvalue (n*pi/L)^2.
    mu_n = (n * pi / L)^2;

    % Reaction Jacobian entries at (u_star, v_star).
    a11 = 1 - 2*u_star/k - (m * v_star) / (u_star + 1)^2;
    a12 = - (m * u_star) / (u_star + 1);
    a21 =   (m * v_star) / (u_star + 1)^2;
    a22 = -theta + (m * u_star) / (u_star + 1);

    % Mode-dependent linearization, including diffusion.
    L_n = [a11 - d1 * mu_n,  a12; ...
           a21,              a22 - d2 * mu_n];

    % Eigenvalues and eigenvectors.
    opts.disp = 0;
    if isfield(cfg, 'eigopts')
        fn = fieldnames(cfg.eigopts);
        for kf = 1:numel(fn)
            opts.(fn{kf}) = cfg.eigopts.(fn{kf});
        end
    end

    [W, D] = eig(L_n); %#ok<ASGLU>
    eigvals = diag(D);

    % Pick the eigenvalue with positive imaginary part.
    [~, idx] = max(imag(eigvals));
    omega_H = imag(eigvals(idx));
    psi = W(:, idx);
    psi = psi / norm(psi);

    hopf = struct();
    hopf.mode    = n;
    hopf.mu_n    = mu_n;
    hopf.L_n     = L_n;
    hopf.eigvals = eigvals;
    hopf.eigvecs = W;
    hopf.lambda  = cfg.lambda;
    hopf.omega   = omega_H;
    hopf.psi     = psi;
    hopf.T       = 2*pi / abs(omega_H);
    hopf.u_star  = u_star;
    hopf.v_star  = v_star;
    hopf.cfg     = cfg;
end

function init = construct_initial_guess(hopf)
%CONSTRUCT_INITIAL_GUESS  Construct a small-amplitude periodic-orbit guess.

    cfg = hopf.cfg;
    omega = hopf.omega;
    psi = hopf.psi;
    u_star = hopf.u_star;
    v_star = hopf.v_star;
    eps_amp = cfg.periodic_orbit.amplitude;

    % Space-time discretization.
    Nx = cfg.periodic_orbit.Nx;
    L  = cfg.L;
    x  = linspace(0, L, Nx).';
    Nt = cfg.periodic_orbit.Nt;
    T_H = hopf.T;
    t  = linspace(0, T_H, Nt + 1);
    t(end) = [];

    % Spatial eigenmode.
    n = hopf.mode;
    phi_n = cos(n*pi*x/L);

    % Complex eigenvector components.
    psi_u    = real(psi(1));
    psi_v    = real(psi(2));
    psi_u_im = imag(psi(1));
    psi_v_im = imag(psi(2));

    % Space-time fields.
    U0 = zeros(Nx, Nt);
    V0 = zeros(Nx, Nt);

    for kt = 1:Nt
        phase = omega * t(kt);
        U0(:, kt) = u_star + eps_amp * ...
            (psi_u * cos(phase) - psi_u_im * sin(phase)) .* phi_n;
        V0(:, kt) = v_star + eps_amp * ...
            (psi_v * cos(phase) - psi_v_im * sin(phase)) .* phi_n;
    end

    init.U0 = U0;
    init.V0 = V0;
    init.x = x;
    init.t = t;
    init.T_H = T_H;
end

function cfg = update_cfg_for_periodic_orbit(cfg)
%UPDATE_CFG_FOR_PERIODIC_ORBIT  Shift lambda and update route-dependent data.

    require_route_fun(cfg);

    delta_lambda = cfg.periodic_orbit.delta_lambda;
    cfg.lambda = cfg.lambda + delta_lambda;

    [m, theta] = cfg.route.fun(cfg.lambda);
    cfg.theta = theta;
    cfg.m = m;
    cfg.u_star = cfg.lambda;
    cfg.v_star = (cfg.k - cfg.lambda) * (1 + cfg.lambda) / (cfg.k * cfg.m);
end

function require_route_fun(cfg)
%REQUIRE_ROUTE_FUN  Check that the configuration supplies the parameter route.

    if ~isfield(cfg, 'route') || ...
            ~isfield(cfg.route, 'fun') || isempty(cfg.route.fun)
        error(['The configuration must define cfg.route.fun. ', ...
               'For the public workflow, this is set in case_ex262.m.']);
    end
end
