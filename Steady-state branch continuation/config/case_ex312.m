function cfg = case_ex312(lambda)
%CASE_EX312  Configuration for Test case 3, Example 3.12.
%
% This public configuration uses the theta-fixed route from the manuscript
% computations.  The field cfg.route.type is kept as a descriptive tag for
% filenames and metadata; it is not an editable route selector.

if nargin < 1 || isempty(lambda)
    error('case_ex312:MissingLambda', ...
        'A lambda value must be supplied.');
end

%% Physical and kinetic parameters
cfg.d1  = 1.0;
cfg.d2  = 1.0;
cfg.k   = 3.0;
cfg.ell = 35.0;

%% Primary parameter
cfg.lambda = lambda;

%% Spatial domain and discretization
cfg.Nx = 128;
cfg.L  = cfg.ell * pi;
cfg.x  = linspace(0, cfg.L, cfg.Nx)';
cfg.dx = cfg.x(2) - cfg.x(1);

%% Theta-fixed route lambda -> (m, theta)
theta0 = 0.003;

cfg.route = struct();
cfg.route.type       = 'thetaFixed';
cfg.route.theta0     = theta0;
cfg.route.lambdaL    = 0.377251393220232;
cfg.route.lambdaR    = 0.676167479810115;
cfg.route.anchorDesc = 'mode5-left-steady';
cfg.route.fun        = @(lambdaValue) deal( ...
    theta0 * (1 + 1./lambdaValue), ...
    theta0);

% Lambda snapping is retained for consistency with the continuation code.
cfg.route.snap = struct();
cfg.route.snap.enable = true;
cfg.route.snap.h      = 1e-8;
cfg.route.snap.mode   = 'round';

[cfg.m, cfg.theta] = cfg.route.fun(cfg.lambda);

%% Coexistence equilibrium
cfg.u_star = cfg.lambda;
cfg.v_star = cfg.lambda * (1 - cfg.lambda/cfg.k) / cfg.theta;

%% Seed options
cfg.seed = struct();
cfg.seed.eps0      = 8e-2;
cfg.seed.lambdaBif = cfg.route.lambdaL;
cfg.seed.side      = -1;   % switch to +1 if start_steady_branch fails

% Nearby grid point used for the initial Newton correction.
h = 1e-5;
cfg.seed.lambdaGrid = h;
if cfg.seed.side >= 0
    cfg.seed.lambda0 = (floor(cfg.seed.lambdaBif/h) + 1) * h;
else
    cfg.seed.lambda0 = (ceil(cfg.seed.lambdaBif/h) - 1) * h;
end

%% Continuation options
cfg.cont = struct();
cfg.cont.usePALC  = false;
cfg.cont.ds       = 1e-2;
cfg.cont.dsMin    = 1e-4;
cfg.cont.dsMax    = 5e-2;
cfg.cont.maxSteps = 200;
cfg.cont.param    = 'lambda';

%% Observable options
cfg.obs = struct();
cfg.obs.ampType = 'L2';

%% Metadata
cfg.caseName = 'Ex3.12';
cfg.workflow = 'steady-state branch continuation';

end
