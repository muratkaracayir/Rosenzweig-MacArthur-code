function cfg = case_ex262(lambda)
%CASE_EX262  Yi--Wei--Shi (2009) Example 2.6.2 case definition.
%
%   cfg = case_ex262(lambda)
%
%   Return the model parameters for the Example 2.6.2 parameter set.
%   The input lambda is the main bifurcation parameter used in the
%   theta-fixed parametrization
%
%       m = theta * (1 + 1/lambda).
%
%   The input lambda is assumed positive. In this example, k = 17, so the
%   spatially homogeneous Hopf value is lambda_0^H = (k - 1)/2 = 8. To
%   observe convergence to the homogeneous limit cycle, choose
%   0 < lambda < 8.
%
%   This function returns only model/case parameters and reference
%   quantities associated with lambda. Numerical choices such as Nx, dt,
%   Tfinal, initial data, and output controls are supplied by the driver.

    if nargin < 1
        error('You must provide the bifurcation parameter lambda.');
    end

    % --- Physical / model parameters ---
    cfg.d1    = 1.0;
    cfg.d2    = 3.0;
    cfg.k     = 17.0;
    cfg.theta = 4.0;
    cfg.ell   = 4*sqrt(85)/5;
    cfg.L     = cfg.ell * pi;

    % --- Main bifurcation parameter ---
    cfg.lambda = lambda;
    cfg.m      = cfg.theta * (1 + 1/cfg.lambda);

    % --- Reference values associated with lambda ---
    cfg.u_star = cfg.lambda;
    cfg.v_star = cfg.u_star * (1 - cfg.u_star / cfg.k) / cfg.theta;

    % --- Fixed case identifier ---
    cfg.caseName = 'Ex2.6.2';
end
