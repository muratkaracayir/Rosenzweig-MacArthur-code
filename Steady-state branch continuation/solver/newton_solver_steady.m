function [SS, UV, resNorm] = newton_solver_steady(UV, cfg, steady, doSave)
%NEWTON_SOLVER_STEADY  Newton correction for one steady state profile.
%
%   [SS, UV, resNorm] = newton_solver_steady(UV, cfg, steady, doSave)
%
% Inputs
%   UV     : cfg.Nx-by-2 initial guess.  The two columns are u and v.
%   cfg    : configuration structure at the target lambda value.
%   steady : metadata structure for the steady branch seed.
%   doSave : optional logical flag.  If true, a converged nontrivial steady
%            state is saved to the workflow data folder.
%
% Outputs
%   SS      : steady-state solution structure.
%   UV      : corrected cfg.Nx-by-2 profile.
%   resNorm : final Euclidean norm of the residual.
%
% This solver expects residual_and_jacobian_steady.m to be on the MATLAB
% path.  The public workflow uses the theta-fixed route; route information
% is used here only for file naming.

if nargin < 4 || isempty(doSave)
    doSave = false;
end

validate_inputs(UV, cfg, steady, doSave);

% Newton settings.  These are deliberately local to keep the public workflow
% minimal and reproducible.
maxIter = 100;
resTol  = 1e-10;
stepTol = 1e-10;
damping = 1.0;
ampMin  = 1e-8;
verbose = true;

Nx = cfg.Nx;
pack   = @(UVnow) [UVnow(:, 1); UVnow(:, 2)];
unpack = @(w) [w(1:Nx), w(Nx+1:end)];

resNorm   = NaN;
deltaNorm = NaN;
converged = false;
stagnated = false;

for iter = 1:maxIter
    [R, J] = residual_and_jacobian_steady(UV, cfg);
    resNorm = norm(R, 2);

    if verbose
        fprintf('Iter %2d: ||R|| = %.3e\n', iter - 1, resNorm);
    end

    if resNorm < resTol
        converged = true;
        break;
    end

    delta = sparse(J) \ (-R);
    deltaNorm = norm(delta, 2);

    w  = pack(UV);
    w  = w + damping * delta;
    UV = unpack(w);

    Rnew = residual_and_jacobian_steady(UV, cfg);
    resNorm = norm(Rnew, 2);

    if verbose
        fprintf('          ||R_new|| = %.3e, ||delta|| = %.3e\n', ...
            resNorm, deltaNorm);
    end

    if resNorm < resTol
        converged = true;
        break;
    end

    if deltaNorm < stepTol
        stagnated = true;
        break;
    end
end

SS = make_solution_struct(UV, cfg, steady, resNorm, converged, ...
    stagnated, deltaNorm, maxIter, resTol, stepTol, damping);

if doSave
    if converged && SS.amp.ampEqNorm > ampMin
        save_SS_data(SS);
    elseif converged
        fprintf(['Converged to a near-equilibrium steady state; ', ...
            'not saving SS (ampEqNorm = %.3e).\n'], SS.amp.ampEqNorm);
    else
        fprintf(['Newton did not satisfy the residual tolerance; ', ...
            'not saving SS (resNorm = %.3e).\n'], resNorm);
    end
end
end

% ---------------------------------------------------------------------- %
function validate_inputs(UV, cfg, steady, doSave)
%VALIDATE_INPUTS  Check the inputs needed by the Newton solver.

if nargin < 3
    error('newton_solver_steady:NotEnoughInputs', ...
        'UV, cfg, and steady must be supplied.');
end

if ~isstruct(cfg)
    error('newton_solver_steady:InvalidCfg', ...
        'cfg must be a structure.');
end

requiredFields = {'Nx', 'lambda', 'caseName', 'u_star', 'v_star'};
for j = 1:numel(requiredFields)
    fieldName = requiredFields{j};
    if ~isfield(cfg, fieldName)
        error('newton_solver_steady:MissingField', ...
            'cfg.%s must be defined.', fieldName);
    end
end

if ~isnumeric(cfg.Nx) || ~isscalar(cfg.Nx) || ~isfinite(cfg.Nx) || ...
        cfg.Nx ~= round(cfg.Nx) || cfg.Nx < 2
    error('newton_solver_steady:InvalidNx', ...
        'cfg.Nx must be an integer greater than or equal to 2.');
end

numericFields = {'lambda', 'u_star', 'v_star'};
for j = 1:numel(numericFields)
    fieldName = numericFields{j};
    value = cfg.(fieldName);
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
        error('newton_solver_steady:InvalidField', ...
            'cfg.%s must be a finite numeric scalar.', fieldName);
    end
end

if ~ischar(cfg.caseName)
    error('newton_solver_steady:InvalidCaseName', ...
        'cfg.caseName must be a character vector.');
end

if ~isnumeric(UV) || ndims(UV) ~= 2 || ...
        size(UV, 1) ~= cfg.Nx || size(UV, 2) ~= 2
    error('newton_solver_steady:InvalidUV', ...
        'UV must be a numeric cfg.Nx-by-2 array.');
end

if any(~isfinite(UV(:)))
    error('newton_solver_steady:NonfiniteUV', ...
        'UV must contain only finite values.');
end

if ~isstruct(steady)
    error('newton_solver_steady:InvalidSteady', ...
        'steady must be a structure.');
end

if ~islogical(doSave) && ~(isnumeric(doSave) && isscalar(doSave))
    error('newton_solver_steady:InvalidDoSave', ...
        'doSave must be a logical scalar.');
end
end

% ---------------------------------------------------------------------- %
function SS = make_solution_struct(UV, cfg, steady, resNorm, converged, ...
    stagnated, deltaNorm, maxIter, resTol, stepTol, damping)
%MAKE_SOLUTION_STRUCT  Assemble the output steady-state structure.

SS = struct();
SS.lambda  = cfg.lambda;
SS.UV      = UV;
SS.cfg     = cfg;
SS.steady  = steady;
SS.resNorm = resNorm;
SS.amp     = compute_amp_eq(UV, cfg);

SS.newton = struct();
SS.newton.converged = converged;
SS.newton.stagnated = stagnated;
SS.newton.deltaNorm = deltaNorm;
SS.newton.maxIter   = maxIter;
SS.newton.resTol    = resTol;
SS.newton.stepTol   = stepTol;
SS.newton.damping   = damping;
end

% ---------------------------------------------------------------------- %
function amp = compute_amp_eq(UV, cfg)
%COMPUTE_AMP_EQ  Equilibrium-referenced amplitude diagnostics.
%
%   ampEq     = max_x ||(u,v) - (u*,v*)||_2
%   ampEqNorm = ampEq / ||(u*,v*)||_2

Nx = cfg.Nx;

ueq = cfg.u_star;
veq = cfg.v_star;

EQ  = [ueq * ones(Nx, 1), veq * ones(Nx, 1)];
dUV = UV - EQ;

ampEq = max(sqrt(dUV(:, 1).^2 + dUV(:, 2).^2));

eqMag = sqrt(ueq^2 + veq^2);
if eqMag > 0
    ampEqNorm = ampEq / eqMag;
else
    ampEqNorm = NaN;
end

amp = struct();
amp.ampEq     = ampEq;
amp.ampEqNorm = ampEqNorm;
amp.ampDesc   = ['ampEq = max_x ||(u,v)-(u*,v*)||_2; ', ...
    'ampEqNorm = ampEq / ||(u*,v*)||_2'];
end

% ---------------------------------------------------------------------- %
function save_SS_data(SS)
%SAVE_SS_DATA  Save one converged steady-state solution to data/.

cfg = SS.cfg;

rootDir = fileparts(fileparts(mfilename('fullpath')));
dataDir = fullfile(rootDir, 'data');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

lambda   = cfg.lambda;
routeTag = get_route_tag(cfg);
modeTag  = get_mode_tag(SS);

fileName = sprintf('SS_%s_%s%s_lambda%.7g.mat', ...
    cfg.caseName, modeTag, routeTag, lambda);
filePath = fullfile(dataDir, fileName);

save(filePath, 'SS');
fprintf('SS saved to:\n  %s\n', filePath);
end

% ---------------------------------------------------------------------- %
function modeTag = get_mode_tag(SS)
%GET_MODE_TAG  Return a filename tag for the spatial mode, if available.

modeTag = '';

if isfield(SS, 'steady') && isstruct(SS.steady) && ...
        isfield(SS.steady, 'mode') && ~isempty(SS.steady.mode)
    modeTag = sprintf('mode%d_', SS.steady.mode);
end
end

% ---------------------------------------------------------------------- %
function routeTag = get_route_tag(cfg)
%GET_ROUTE_TAG  Route tag for the public theta-fixed workflow.

if ~isfield(cfg, 'route') || ~isstruct(cfg.route) || ...
        ~isfield(cfg.route, 'type') || isempty(cfg.route.type)
    routeTag = 'routeUnknown';
    return;
end

if ~strcmp(cfg.route.type, 'thetaFixed')
    error('newton_solver_steady:UnsupportedRoute', ...
        'This public workflow supports only cfg.route.type = ''thetaFixed''.');
end

routeTag = cfg.route.type;
end
