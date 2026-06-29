function cfg = case_ex262(lambda)
%CASE_EX262  Configuration for Example 2.6.2.
%
%   cfg = case_ex262(lambda)
%
%   Minimal public configuration for the Hopf-branch continuation workflow.
%   The parameter route is the theta-fixed route used in the manuscript.

% --- Physical and kinetic parameters ---
cfg.d1    = 1.0;
cfg.d2    = 3.0;
cfg.k     = 17.0;
cfg.theta = 4.0;
cfg.ell   = 4*sqrt(85)/5;   % domain is (0, ell*pi)

% --- Primary bifurcation parameter ---
cfg.lambda = lambda;

% --- Theta-fixed parameter route ---
theta0 = cfg.theta;
cfg.route.type   = 'thetaFixed';
cfg.route.theta0 = theta0;
cfg.route.fun    = @(lambdaValue) deal( ...
    theta0 * (1 + 1./lambdaValue), ...
    theta0 );

[cfg.m, cfg.theta] = cfg.route.fun(cfg.lambda);

% --- Coexistence equilibrium ---
cfg.u_star = cfg.lambda;
cfg.v_star = cfg.u_star * (1 - cfg.u_star/cfg.k) / cfg.theta;

% --- Spatial domain and discretization ---
cfg.Nx = 128;
cfg.L  = cfg.ell * pi;
cfg.x  = linspace(0, cfg.L, cfg.Nx)';
cfg.dx = cfg.x(2) - cfg.x(1);

% --- Periodic-orbit discretization and initialization ---
cfg.periodic_orbit.Nx            = cfg.Nx;
cfg.periodic_orbit.Nt            = 48;
cfg.periodic_orbit.amplitude     = 1e-2;
cfg.periodic_orbit.amplitudeType = 2;

grid = 1e-4;
cfg.periodic_orbit.delta_lambda = ceil(cfg.lambda/grid)*grid - cfg.lambda;
if cfg.periodic_orbit.delta_lambda == 0
    cfg.periodic_orbit.delta_lambda = grid;
end

% --- Eigenvalue solver options ---
cfg.eigopts.disp  = 0;
cfg.eigopts.maxit = 200;
cfg.eigopts.tol   = 1e-10;

% --- Case name for saved data ---
cfg.caseName = 'Ex2.6.2';

end
