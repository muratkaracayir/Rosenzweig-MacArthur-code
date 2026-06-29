function Uf = diffusion_step_DCT(U, cfg, Nx, dt)
%DIFFUSION_STEP_DCT  Exact cosine-spectral diffusion step with Neumann BCs.
%
%   Uf = diffusion_step_DCT(U, cfg, Nx, dt)
%
%   Advances the diffusion part of the semidiscrete system over time dt on
%   an endpoint-including uniform grid. The transform is a DCT-I-style
%   cosine transform implemented by dct1_endpoints.m.
%
%   Inputs:
%     U    2 x Nx array of nodal values [u; v]
%     cfg  structure with fields d1 and d2, and either ell or L
%     Nx   number of endpoint grid points
%     dt   diffusion time, usually one half of the Strang time step
%
%   The cosine modes have Neumann Laplacian eigenvalues -(q/ell)^2,
%   q = 0,...,Nx-1, where L = ell*pi.
    % --- Basic checks ---
    if ~isnumeric(U) || ~isreal(U) || ndims(U) ~= 2 || size(U,1) ~= 2
        error('diffusion_step_DCT: U must be a real 2 x Nx array.');
    end

    if ~isscalar(Nx) || Nx < 1 || Nx ~= round(Nx)
        error('diffusion_step_DCT: Nx must be a positive integer.');
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