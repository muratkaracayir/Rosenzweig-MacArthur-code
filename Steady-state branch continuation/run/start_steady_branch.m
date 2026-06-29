function start_steady_branch()
%START_STEADY_BRANCH  Start one nonconstant steady branch.
%
% This driver constructs a small-amplitude steady-state seed near a selected
% steady bifurcation point, applies Newton correction at the nearby seed
% parameter value, and saves the resulting steady state to the workflow data
% folder when requested.
%
% Intended location in the public workflow:
%   Steady-state branch continuation/run/start_steady_branch.m

%% Path setup
thisFile = mfilename('fullpath');
runDir   = fileparts(thisFile);
rootDir  = fileparts(runDir);

addpath(fullfile(rootDir, 'config'));
addpath(fullfile(rootDir, 'solver'));

%% User settings
caseFun   = @case_ex312;
modeIndex = 5;

% This value is used only to instantiate the configuration. The branch seed
% itself uses cfg.seed.lambdaBif and cfg.seed.lambda0.
lambdaForConfig = 0.377251393220232;

doSaveSS          = true;

%% Build configuration
cfg = caseFun(lambdaForConfig);

% Use the bifurcation value stored in the configuration as the anchor point.
cfgBif = apply_route(cfg, cfg.seed.lambdaBif);

%% Prepare data folder
dataDir = fullfile(rootDir, 'data');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

%% Build initial seed near the bifurcation point
[steady, UV0, cfg0] = make_steady_seed(cfgBif, modeIndex);

fprintf('Seed constructed:\n');
fprintf('  mode      = %d\n', steady.mode);
fprintf('  lambdaBif = %.16g\n', steady.lambdaL);
fprintf('  lambda0   = %.16g\n', steady.lambda0);
fprintf('  eps0      = %.3g\n', steady.eps0);
fprintf('  side      = %+d\n', steady.side);
fprintf('  ampEq     = %.3e\n', steady.seedAmpEq);

%% Newton correction
[SS, UV, resNorm] = newton_solver_steady(UV0, cfg0, steady, doSaveSS); %#ok<ASGLU>

fprintf('Newton done:\n');
fprintf('  lambda    = %.16g\n', SS.lambda);
fprintf('  resNorm   = %.3e\n', resNorm);
fprintf('  ampEqNorm = %.3e\n', SS.amp.ampEqNorm);
end

% =========================================================================
function [steady, UV0, cfg0] = make_steady_seed(cfg, modeIndex)
%MAKE_STEADY_SEED  Build a small-amplitude steady seed near lambdaBif.
%
% The Neumann eigenfunction on x in [0,L], L = ell*pi, is
%     phi(x) = cos(n*pi*x/L) = cos(n*x/ell).

cfg0 = apply_route(cfg, cfg.seed.lambda0);

n = modeIndex;
phi = cos(n * pi * cfg0.x / cfg0.L);

cfgBif = apply_route(cfg, cfg.seed.lambdaBif);
q = kern_vec_2x2(cfgBif, n);
q = q / norm(q);
if q(1) < 0
    q = -q;
end

eps0 = cfg.seed.eps0;
side = cfg.seed.side;

U0 = cfg0.u_star + side * eps0 * phi * q(1);
V0 = cfg0.v_star + side * eps0 * phi * q(2);
UV0 = [U0, V0];

if any(UV0(:) <= 0)
    warning('start_steady_branch:NonpositiveSeed', ...
        'Initial steady seed has nonpositive entries; consider a smaller cfg.seed.eps0.');
end

seedAmpEq = compute_seed_amp_eq(UV0, cfg0);

steady = struct();
steady.mode      = n;
steady.lambdaL   = cfg.seed.lambdaBif;
steady.lambdaR   = cfg.route.lambdaR;
steady.lambda0   = cfg.seed.lambda0;
steady.side      = side;
steady.eps0      = eps0;
steady.qKern     = q;
steady.seedAmpEq = seedAmpEq;
steady.note      = cfg.route.anchorDesc;
end

% =========================================================================
function ampEq = compute_seed_amp_eq(UV, cfg)
%COMPUTE_SEED_AMP_EQ  Equilibrium-referenced seed amplitude.

EQ = [cfg.u_star * ones(cfg.Nx, 1), cfg.v_star * ones(cfg.Nx, 1)];
dUV = UV - EQ;
ampEq = max(sqrt(dUV(:,1).^2 + dUV(:,2).^2));
end

% =========================================================================
function q = kern_vec_2x2(cfgBif, n)
%KERN_VEC_2X2  Null vector for the mode-reduced steady bifurcation matrix.

d1    = cfgBif.d1;
d2    = cfgBif.d2;
k     = cfgBif.k;
theta = cfgBif.theta;
ell   = cfgBif.ell;
lam   = cfgBif.lambda;

p = (n / ell)^2;

A = lam * (k - 1 - 2 * lam) / (k * (1 + lam));
C = (k - lam) / (k * (1 + lam));

L = [A - d1 * p, -theta; ...
     C,          -d2 * p];

q1 = [-L(1,2); L(1,1)];
q2 = [-L(2,2); L(2,1)];

if norm(q1) >= norm(q2)
    q = q1;
else
    q = q2;
end

if norm(q) == 0
    error('start_steady_branch:ZeroKernelVector', ...
        'Kernel vector construction failed.');
end
end

function cfg = apply_route(cfg, lambdaNew)
%APPLY_ROUTE  Update configuration fields at a new lambda value.
%
%   cfg = apply_route(cfg, lambdaNew)
%
% This local helper is used by the steady-state branch workflow.  The public
% version supports the theta-fixed route used in Test case 3.  The route type
% is kept as descriptive metadata in cfg.route.type, while the actual map
% lambda -> (m, theta) is supplied by cfg.route.fun.

if nargin < 2
    error('apply_route:NotEnoughInputs', ...
        'Both cfg and lambdaNew must be supplied.');
end

if ~isnumeric(lambdaNew) || ~isscalar(lambdaNew) || ~isfinite(lambdaNew)
    error('apply_route:InvalidLambda', ...
        'lambdaNew must be a finite numeric scalar.');
end

if lambdaNew <= 0
    error('apply_route:InvalidLambda', ...
        'lambdaNew must be positive.');
end

if ~isfield(cfg, 'k') || ~isnumeric(cfg.k) || ~isscalar(cfg.k)
    error('apply_route:MissingK', ...
        'cfg.k must be defined as a numeric scalar.');
end

if ~isfield(cfg, 'route') || ~isstruct(cfg.route)
    error('apply_route:MissingRoute', ...
        'cfg.route must be defined in the configuration.');
end

if ~isfield(cfg.route, 'type')
    error('apply_route:MissingRouteType', ...
        'cfg.route.type must be defined.');
end

if ~strcmp(cfg.route.type, 'thetaFixed')
    error('apply_route:UnsupportedRoute', ...
        'This public workflow supports only cfg.route.type = ''thetaFixed''.');
end

if ~isfield(cfg.route, 'fun') || ~isa(cfg.route.fun, 'function_handle')
    error('apply_route:MissingRouteFunction', ...
        'cfg.route.fun must be defined as a function handle.');
end

lambdaNew = snap_lambda(lambdaNew, cfg);

[mNow, thetaNow] = cfg.route.fun(lambdaNew);

if ~isnumeric(mNow) || ~isscalar(mNow) || ~isfinite(mNow) || mNow <= 0
    error('apply_route:InvalidM', ...
        'The route map returned an invalid m value at lambda = %.16g.', ...
        lambdaNew);
end

if ~isnumeric(thetaNow) || ~isscalar(thetaNow) || ...
        ~isfinite(thetaNow) || thetaNow <= 0
    error('apply_route:InvalidTheta', ...
        'The route map returned an invalid theta value at lambda = %.16g.', ...
        lambdaNew);
end

cfg.lambda = lambdaNew;
cfg.m      = mNow;
cfg.theta  = thetaNow;

cfg.u_star = lambdaNew;
cfg.v_star = (cfg.k - lambdaNew) * (1 + lambdaNew) / (cfg.k * cfg.m);
end

function lambdaNew = snap_lambda(lambdaNew, cfg)
%SNAP_LAMBDA  Optionally snap lambda to the route grid.

if ~isfield(cfg, 'route') || ~isfield(cfg.route, 'snap')
    return;
end

snap = cfg.route.snap;

if ~isfield(snap, 'enable') || ~snap.enable
    return;
end

if ~isfield(snap, 'h') || ~isnumeric(snap.h) || ...
        ~isscalar(snap.h) || snap.h <= 0
    error('snap_lambda:InvalidStep', ...
        'cfg.route.snap.h must be a positive scalar when snapping is enabled.');
end

if isfield(snap, 'mode')
    snapMode = snap.mode;
else
    snapMode = 'round';
end

switch snapMode
    case 'round'
        lambdaNew = round(lambdaNew / snap.h) * snap.h;
    case 'floor'
        lambdaNew = floor(lambdaNew / snap.h) * snap.h;
    case 'ceil'
        lambdaNew = ceil(lambdaNew / snap.h) * snap.h;
    otherwise
        error('snap_lambda:UnknownMode', ...
            'Unknown cfg.route.snap.mode: %s', snapMode);
end
end
