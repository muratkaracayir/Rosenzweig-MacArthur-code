function [R, J] = residual_and_jacobian_steady(UV, cfg)
%RESIDUAL_AND_JACOBIAN_STEADY  Residual and Jacobian for the steady PDE.
%
%   R = residual_and_jacobian_steady(UV, cfg)
%   [R, J] = residual_and_jacobian_steady(UV, cfg)
%
% The endpoint grid is x = linspace(0, cfg.L, cfg.Nx).  Neumann boundary
% conditions are imposed through a DCT-I spectral Laplacian on this grid.
%
% Input
%   UV  : cfg.Nx-by-2 array whose columns are the steady profiles u and v.
%   cfg : configuration structure with fields
%           Nx, L, d1, d2, k, m, theta.
%
% Output
%   R   : stacked residual [Ru; Rv], where
%           Ru = d1*u_xx + u*(1-u/k) - m*u*v/(u+1),
%           Rv = d2*v_xx - theta*v + m*u*v/(u+1).
%   J   : sparse Jacobian of R with respect to [u; v].

validate_inputs(UV, cfg);

Nx = cfg.Nx;
Lx = cfg.L;

u = UV(:, 1);
v = UV(:, 2);

% Diffusion terms.
uxx = laplacian_dct1(u, Lx);
vxx = laplacian_dct1(v, Lx);

% Reaction terms and pointwise derivatives.
[f, g, fu, fv, gu, gv] = reaction_and_derivatives(u, v, cfg);

% Stacked residual.
Ru = cfg.d1 * uxx + f;
Rv = cfg.d2 * vxx + g;
R  = [Ru; Rv];

if nargout < 2
    return;
end

% Sparse Jacobian.  The cached Dxx matrix is built using the same
% DCT-I Laplacian routine used above for the residual.
Dxx = get_cached_Dxx(cfg);

J11 = cfg.d1 * Dxx + spdiags(fu, 0, Nx, Nx);
J12 =                 spdiags(fv, 0, Nx, Nx);
J21 =                 spdiags(gu, 0, Nx, Nx);
J22 = cfg.d2 * Dxx + spdiags(gv, 0, Nx, Nx);

J = [J11, J12;
     J21, J22];
end

% ---------------------------------------------------------------------- %
function validate_inputs(UV, cfg)
%VALIDATE_INPUTS  Check the fields required by the residual/Jacobian.

if nargin < 2
    error('residual_and_jacobian_steady:NotEnoughInputs', ...
        'Both UV and cfg must be supplied.');
end

if ~isstruct(cfg)
    error('residual_and_jacobian_steady:InvalidCfg', ...
        'cfg must be a structure.');
end

requiredFields = {'Nx', 'L', 'd1', 'd2', 'k', 'm', 'theta'};
for j = 1:numel(requiredFields)
    fieldName = requiredFields{j};
    if ~isfield(cfg, fieldName)
        error('residual_and_jacobian_steady:MissingField', ...
            'cfg.%s must be defined.', fieldName);
    end
    value = cfg.(fieldName);
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
        error('residual_and_jacobian_steady:InvalidField', ...
            'cfg.%s must be a finite numeric scalar.', fieldName);
    end
end

if cfg.Nx ~= round(cfg.Nx) || cfg.Nx < 2
    error('residual_and_jacobian_steady:InvalidNx', ...
        'cfg.Nx must be an integer greater than or equal to 2.');
end

positiveFields = {'L', 'd1', 'd2', 'k', 'm', 'theta'};
for j = 1:numel(positiveFields)
    fieldName = positiveFields{j};
    if cfg.(fieldName) <= 0
        error('residual_and_jacobian_steady:InvalidField', ...
            'cfg.%s must be positive.', fieldName);
    end
end

if ~isnumeric(UV) || ndims(UV) ~= 2
    error('residual_and_jacobian_steady:InvalidUV', ...
        'UV must be a numeric cfg.Nx-by-2 array.');
end

if size(UV, 1) ~= cfg.Nx || size(UV, 2) ~= 2
    error('residual_and_jacobian_steady:InvalidUVSize', ...
        'UV must have size cfg.Nx-by-2, with cfg.Nx = %d.', cfg.Nx);
end

if any(~isfinite(UV(:)))
    error('residual_and_jacobian_steady:NonfiniteUV', ...
        'UV must contain only finite values.');
end

if any(abs(UV(:, 1) + 1) < 100 * eps)
    error('residual_and_jacobian_steady:SingularReaction', ...
        'The reaction term is singular when u is close to -1.');
end
end

% ---------------------------------------------------------------------- %
function uxx = laplacian_dct1(u, Lx)
%LAPLACIAN_DCT1  Neumann spectral Laplacian on an endpoint grid.

u = u(:);
Nx = length(u);

coeffs = dct1_endpoint(u);
modeIndex = (0:Nx-1).';
lambda = -(pi * modeIndex / Lx).^2;

uxx = dct1_endpoint(lambda .* coeffs) * (Nx - 1) / 2;
end

% ---------------------------------------------------------------------- %
function coeffs = dct1_endpoint(u)
%DCT1_ENDPOINT  DCT-I transform with endpoint-grid normalization.
%
% With this convention, applying the transform twice gives the inverse up
% to the factor (Nx - 1)/2 used in laplacian_dct1.

u = u(:);
Nx = length(u);

if Nx == 1
    coeffs = u;
    return;
end

modeIndex = (0:Nx-1).';
interiorIndex = 1:Nx-2;

cosMatrix = cos(pi * modeIndex * interiorIndex / (Nx - 1));
coeffs = 0.5 * u(1) + 0.5 * u(end) * (-1).^modeIndex;
coeffs = coeffs + cosMatrix * u(2:end-1);
coeffs = coeffs * 2 / (Nx - 1);
end

% ---------------------------------------------------------------------- %
function [f, g, fu, fv, gu, gv] = reaction_and_derivatives(u, v, cfg)
%REACTION_AND_DERIVATIVES  Holling type-II kinetics and derivatives.

m = cfg.m;
k = cfg.k;
theta = cfg.theta;

den = u + 1;

f = u .* (1 - u / k) - m * u .* v ./ den;
g = -theta * v + m * u .* v ./ den;

fu = 1 - 2 * u / k - (m * v) ./ den + (m * u .* v) ./ den.^2;
fv = -m * u ./ den;

gu = (m * v) ./ den - (m * u .* v) ./ den.^2;
gv = -theta + (m * u) ./ den;
end

% ---------------------------------------------------------------------- %
function Dxx = get_cached_Dxx(cfg)
%GET_CACHED_DXX  Build and cache the DCT-I Laplacian matrix.

persistent cacheNx cacheL cacheDxx

Nx = cfg.Nx;
Lx = cfg.L;

cacheIsStale = isempty(cacheDxx) || isempty(cacheNx) || isempty(cacheL) || ...
    cacheNx ~= Nx || cacheL ~= Lx;

if cacheIsStale
    DxxDense = zeros(Nx, Nx);
    for j = 1:Nx
        basisVector = zeros(Nx, 1);
        basisVector(j) = 1;
        DxxDense(:, j) = laplacian_dct1(basisVector, Lx);
    end

    cacheDxx = sparse(DxxDense);
    cacheNx  = Nx;
    cacheL   = Lx;
end

Dxx = cacheDxx;
end
