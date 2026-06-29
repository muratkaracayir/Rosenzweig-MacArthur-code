function Uf = diffusion_step_DCT(U, cfg, Nx, dt)
%DIFFUSION_STEP_DCT  Exact cosine-spectral diffusion step with Neumann BCs.
%
%   Uf = diffusion_step_DCT(U, cfg, Nx, dt)
%
%   Input
%     U   2 x Nx array of nodal values [u; v]
%     cfg structure with fields cfg.d1 and cfg.d2, and either cfg.ell or cfg.L
%     Nx  number of endpoint-including spatial grid points
%     dt  diffusion time to advance, usually Delta t/2 inside Strang splitting
%
%   Output
%     Uf  2 x Nx array after the exact modal diffusion step
%
%   Grid and modes
%     The grid is the endpoint-including uniform grid
%
%         x_j = j*L/(Nx-1),   j = 0,...,Nx-1,
%
%     with L = ell*pi. The Neumann cosine modes are
%
%         cos(q*x/ell) = cos(pi*j*q/(Nx-1)),   q = 0,...,Nx-1,
%
%     and the Laplacian eigenvalues are -(q/ell)^2.
%
%   Dependency in this folder
%     dct1_endpoints.

    % --- Basic checks ---
    if ~isnumeric(U) || ~isreal(U) || ndims(U) ~= 2 || size(U,1) ~= 2
        error('diffusion_step_DCT: U must be a real 2 x Nx array.');
    end

    if ~isscalar(Nx) || Nx < 2 || Nx ~= round(Nx)
        error('diffusion_step_DCT: Nx must be an integer at least 2.');
    end

    if size(U,2) ~= Nx
        error('diffusion_step_DCT: size(U,2) must equal Nx.');
    end

    if ~isscalar(dt) || ~isreal(dt) || ~isfinite(dt) || dt < 0
        error('diffusion_step_DCT: dt must be a nonnegative real scalar.');
    end

    if ~isstruct(cfg) || ~all(isfield(cfg, {'d1','d2'}))
        error('diffusion_step_DCT: cfg must contain fields d1 and d2.');
    end

    if cfg.d1 < 0 || cfg.d2 < 0
        error('diffusion_step_DCT: cfg.d1 and cfg.d2 must be nonnegative.');
    end

    if isfield(cfg, 'ell')
        ell = cfg.ell;
    elseif isfield(cfg, 'L')
        ell = cfg.L / pi;
    else
        error('diffusion_step_DCT: cfg must contain either ell or L.');
    end

    if ~isscalar(ell) || ~isreal(ell) || ~isfinite(ell) || ell <= 0
        error('diffusion_step_DCT: ell must be a positive real scalar.');
    end

    % --- Cosine mode indices and exact diffusion multipliers ---
    q = 0:(Nx-1);
    mu = (q / ell).^2;

    Eu = exp(-cfg.d1 * mu * dt);
    Ev = exp(-cfg.d2 * mu * dt);

    % --- Transform to cosine coefficients ---
    Uhat = dct1_endpoints(U, 'forward');

    % --- Exact modal diffusion step ---
    Uhat(1,:) = Eu .* Uhat(1,:);
    Uhat(2,:) = Ev .* Uhat(2,:);

    % --- Transform back to nodal values ---
    Uf = dct1_endpoints(Uhat, 'inverse');

    % --- Positivity check: abort if any nodal value became nonpositive ---
    bad = find(Uf <= 0, 1, 'first');
    if ~isempty(bad)
        [comp, j] = ind2sub(size(Uf), bad);

        if comp == 1
            species = 'u';
        else
            species = 'v';
        end

        warning('diffusion_step_DCT:nonpositiveValue', ...
            ['Diffusion half-step produced a nonpositive nodal value ', ...
             'in component %s at node j = %d.'], species, j);

        error('diffusion_step_DCT:nonpositiveValue', ...
            'Aborting because the diffusion half-step violated positivity.');
    end
end
