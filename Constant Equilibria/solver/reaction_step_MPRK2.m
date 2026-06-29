function Unew = reaction_step_MPRK2(U, cfg, dt)
%REACTION_STEP_MPRK2  Two-stage MPRK reaction step for the RM kinetics.
%
%   Unew = reaction_step_MPRK2(U, cfg, dt)
%
%   Applies the local reaction update used in the Strang-splitting scheme for
%   the Rosenzweig--MacArthur reaction terms
%
%       f(u,v) = u(1-u/k) - m*u*v/(1+u),
%       g(u,v) = -theta*v + m*u*v/(1+u).
%
%   Input
%     U    2 x Nx array after the first diffusion half-step
%     cfg  structure with fields k, m, and theta
%     dt   reaction-step size
%
%   Output
%     Unew 2 x Nx array after the reaction update
%
%   Notes
%     The update uses the explicit sequential formulas corresponding to the
%     manuscript's production-destruction discretization. The formulas assume
%     strictly positive input nodal values. This function checks for obvious
%     failures such as nonpositive inputs, nonpositive denominators, or
%     nonfinite outputs.

    % --- Basic input checks ---
    if ~isnumeric(U) || ndims(U) ~= 2 || size(U,1) ~= 2
        error('reaction_step_MPRK2: U must be a 2 x Nx numeric array.');
    end

    if ~isstruct(cfg) || ~all(isfield(cfg, {'k','m','theta'}))
        error('reaction_step_MPRK2: cfg must contain fields k, m, and theta.');
    end

    if ~isscalar(dt) || ~isreal(dt) || ~isfinite(dt) || dt <= 0
        error('reaction_step_MPRK2: dt must be a positive real scalar.');
    end

    u = U(1,:);
    v = U(2,:);

    if any(~isfinite(u)) || any(~isfinite(v))
        error('reaction_step_MPRK2: input state contains NaN or Inf values.');
    end

    if any(u <= 0) || any(v <= 0)
        error(['reaction_step_MPRK2: the MPRK reaction update assumes ', ...
               'strictly positive nodal values in U.']);
    end

    k = cfg.k;
    m = cfg.m;
    theta = cfg.theta;

    if ~isscalar(k) || ~isscalar(m) || ~isscalar(theta) || ...
       ~isfinite(k) || ~isfinite(m) || ~isfinite(theta) || ...
       k <= 0 || m <= 0 || theta <= 0
        error(['reaction_step_MPRK2: cfg.k, cfg.m, and cfg.theta must be ', ...
               'positive finite scalars.']);
    end

    % Local interaction term H(u,v) = m*u*v/(1+u)
    H0 = m .* u .* v ./ (1 + u);

    % ---------- First stage ----------
    % u^{(1)} = u / [1 - dt + dt*(u/k + H0/u)]
    den_u1 = 1 - dt + dt .* (u ./ k + H0 ./ u);
    if any(den_u1 <= 0)
        error(['reaction_step_MPRK2: nonpositive denominator encountered ', ...
               'in the stage-1 prey update. Try a smaller dt.']);
    end
    u1 = u ./ den_u1;

    % v^{(1)} = [v + dt*H0*(u1/u)] / (1 + theta*dt)
    den_v1 = 1 + theta * dt;
    if den_v1 <= 0
        error('reaction_step_MPRK2: nonpositive stage-1 predator denominator.');
    end
    v1 = (v + dt .* H0 .* (u1 ./ u)) ./ den_v1;

    if any(u1 <= 0) || any(v1 <= 0)
        error(['reaction_step_MPRK2: stage-1 values became nonpositive. ', ...
               'Try a smaller dt.']);
    end

    H1 = m .* u1 .* v1 ./ (1 + u1);

    % ---------- Second stage ----------
    % u^{DR} = u / [1 - dt/(2*u1)*(u + u1 - u^2/k - u1^2/k - H0 - H1)]
    bracket_u2 = u + u1 - (u.^2) ./ k - (u1.^2) ./ k - H0 - H1;
    den_u2 = 1 - (dt ./ (2 .* u1)) .* bracket_u2;
    if any(den_u2 <= 0)
        error(['reaction_step_MPRK2: nonpositive denominator encountered ', ...
               'in the stage-2 prey update. Try a smaller dt.']);
    end
    udr = u ./ den_u2;

    % v^{DR} = [v + dt/2*(H0+H1)*(udr/u1)] /
    %          [1 + theta*dt/2 * (v+v1)/v1]
    den_v2 = 1 + (theta * dt / 2) .* ((v + v1) ./ v1);
    if any(den_v2 <= 0)
        error(['reaction_step_MPRK2: nonpositive denominator encountered ', ...
               'in the stage-2 predator update. Try a smaller dt.']);
    end
    vdr = (v + (dt / 2) .* (H0 + H1) .* (udr ./ u1)) ./ den_v2;

    if any(~isfinite(udr)) || any(~isfinite(vdr))
        error('reaction_step_MPRK2: output contains NaN or Inf values.');
    end

    if any(udr <= 0) || any(vdr <= 0)
        error(['reaction_step_MPRK2: final reaction update became ', ...
               'nonpositive. Try a smaller dt.']);
    end

    Unew = [udr; vdr];
end
