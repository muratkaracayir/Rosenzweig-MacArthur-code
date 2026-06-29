function cfg = case_ex261(lambda)
%CASE_EX261  Yi--Wei--Shi (2009) Example 2.6.1 case definition.
%
%   cfg = case_ex261(lambda)
%
%   Return the model parameters and derived homogeneous coexistence
%   equilibrium for the Example 2.6.1 parameter set.
%
%   The input lambda is the main bifurcation parameter used in the
%   theta-fixed parametrization
%
%       m = theta * (1 + 1/lambda).
%
%   The parameters d1, d2, k, theta, ell, and lambda are assumed to be
%   positive. For the coexistence equilibrium to be positive, lambda must
%   additionally satisfy 0 < lambda < k. When lambda >= k, the relevant
%   nonnegative homogeneous equilibrium is the prey-only state (k,0).
%
%   This function returns only the model/case parameters and derived
%   homogeneous equilibrium quantities. Numerical discretization choices
%   (Nx, dt, Tfinal, IC mode, output controls, etc.) should be supplied
%   externally by the driver script.

    if nargin < 1
        error('You must provide the bifurcation parameter lambda.');
    end

    % --- Physical / model parameters ---
    % All physical/model parameters below are positive.
    cfg.d1    = 1.0;
    cfg.d2    = 3.0;
    cfg.k     = 17.0;
    cfg.theta = 4.0;
    cfg.ell   = 2*sqrt(119)/7;
    cfg.L     = cfg.ell * pi;

    % --- Main bifurcation parameter ---
    % The user supplies lambda. It is assumed positive; in this example,
    % 0 < lambda < k gives a positive coexistence equilibrium, while
    % lambda >= k corresponds to the prey-only equilibrium (k,0).
    cfg.lambda = lambda;
    cfg.m      = cfg.theta * (1 + 1/cfg.lambda);

    % --- Homogeneous coexistence equilibrium ---
    cfg.u_star = cfg.lambda;
    cfg.v_star = cfg.u_star * (1 - cfg.u_star / cfg.k) / cfg.theta;

    % --- Fixed case identifier ---
    cfg.caseName = 'Ex2.6.1';
end
