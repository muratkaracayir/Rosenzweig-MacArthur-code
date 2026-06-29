function Unew = reaction_step_MPRK2(U, cfg, dt)
%REACTION_STEP_MPRK2  Positivity-oriented MPRK2 reaction update.
%
%   Unew = reaction_step_MPRK2(U, cfg, dt)
%
%   Advances the local Rosenzweig--MacArthur reaction subsystem over one
%   time step dt. The update is applied independently at each spatial grid
%   point and is used as the reaction part of the Strang-splitting scheme.
%
%   Inputs:
%     U    2 x Nx array with strictly positive nodal values [u; v]
%     cfg  structure with fields k, m, and theta
%     dt   reaction-step size
%
%   Output:
%     Unew  2 x Nx reaction-updated state
%
%   The formulas assume strictly positive input values. The routine checks
%   for nonpositive inputs and nonpositive update denominators, but it does
%   not independently establish a time-step restriction.
    % Basic input checks
    if ~isnumeric(U) || ndims(U) ~= 2 || size(U,1) ~= 2
        error('U must be a 2 x Nx numeric array.');
    end
    if ~isstruct(cfg) || ~all(isfield(cfg, {'k','m','theta'}))
        error('cfg must be a struct with fields k, m, and theta.');
    end
    if ~isscalar(dt) || ~isreal(dt) || ~isfinite(dt) || dt <= 0
        error('dt must be a positive real scalar.');
    end

    u = U(1,:);
    v = U(2,:);

    if any(~isfinite(u)) || any(~isfinite(v))
        error('Input state contains NaN or Inf values.');
    end
    if any(u <= 0) || any(v <= 0)
        error(['The MPRK reaction update assumes strictly positive nodal ', ...
               'values in U.']);
    end

    k = cfg.k;
    m = cfg.m;
    theta = cfg.theta;

    if ~isscalar(k) || ~isscalar(m) || ~isscalar(theta) || ...
       ~isfinite(k) || ~isfinite(m) || ~isfinite(theta) || k <= 0
        error('cfg.k, cfg.m, cfg.theta must be finite scalars, with k > 0.');
    end

    % Local interaction term H(u,v) = m*u*v/(1+u)
    H0 = m .* u .* v ./ (1 + u);

    % ---------- First stage ----------
    % u^{n,D,(1)} = u^{n,D} / [1 - dt + dt*(u^{n,D}/k + H0/u^{n,D})]
    den_u1 = 1 - dt + dt .* (u ./ k + H0 ./ u);
    if any(den_u1 <= 0)
        error(['Nonpositive denominator encountered in stage-1 prey update. ', ...
               'Try a smaller dt.']);
    end
    u1 = u ./ den_u1;

    % v^{n,D,(1)} = [v^{n,D} + dt*H0*(u1/u)] / (1 + theta*dt)
    den_v1 = 1 + theta * dt;
    if den_v1 <= 0
        error('Nonpositive denominator encountered in stage-1 predator update.');
    end
    v1 = (v + dt .* H0 .* (u1 ./ u)) ./ den_v1;

    if any(u1 <= 0) || any(v1 <= 0)
        error(['Stage-1 values became nonpositive. This indicates the time ', ...
               'step is outside the admissible regime for the manuscript formulas.']);
    end

    H1 = m .* u1 .* v1 ./ (1 + u1);

    % ---------- Second stage ----------
    % u^{n,DR} = u^{n,D} / [1 - dt/(2*u1)*(u + u1 - u^2/k - u1^2/k - H0 - H1)]
    bracket_u2 = u + u1 - (u.^2) ./ k - (u1.^2) ./ k - H0 - H1;
    den_u2 = 1 - (dt ./ (2 .* u1)) .* bracket_u2;
    if any(den_u2 <= 0)
        error(['Nonpositive denominator encountered in stage-2 prey update. ', ...
               'Try a smaller dt.']);
    end
    udr = u ./ den_u2;

    % v^{n,DR} = [v^{n,D} + dt/2*(H0+H1)*(udr/u1)] /
    %            [1 + theta*dt/2 * (v+v1)/v1]
    den_v2 = 1 + (theta * dt / 2) .* ((v + v1) ./ v1);
    if any(den_v2 <= 0)
        error(['Nonpositive denominator encountered in stage-2 predator update. ', ...
               'Try a smaller dt.']);
    end
    vdr = (v + (dt / 2) .* (H0 + H1) .* (udr ./ u1)) ./ den_v2;

    if any(~isfinite(udr)) || any(~isfinite(vdr))
        error('Output contains NaN or Inf values.');
    end
    if any(udr <= 0) || any(vdr <= 0)
        error(['Final reaction update became nonpositive. This indicates the ', ...
               'time step is outside the admissible regime for the manuscript formulas.']);
    end

    Unew = [udr; vdr];
end
