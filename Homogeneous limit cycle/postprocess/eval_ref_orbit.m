function [uRef, vRef, Uref] = eval_ref_orbit(refInt, tQuery)
%EVAL_REF_ORBIT  Evaluate the periodic ODE reference orbit.
%
%   [uRef, vRef, Uref] = eval_ref_orbit(refInt, tQuery)
%
%   evaluates the periodic interpolant returned by
%   build_ref_orbit_interpolant at the physical times tQuery.
%
%   Input
%   -----
%   refInt is the structure returned by build_ref_orbit_interpolant and must
%   contain the fields
%
%       refInt.T
%       refInt.ppU
%       refInt.ppV
%
%   tQuery may be a scalar, vector, or array of physical query times.
%
%   Output
%   ------
%   uRef contains the prey component of the reference orbit at tQuery.
%   vRef contains the predator component of the reference orbit at tQuery.
%   Uref is the stacked 2-by-numel(tQuery) array
%
%       Uref = [uRef(:).'; vRef(:).'];
%
%   Periodicity is enforced by wrapping tQuery modulo refInt.T before the
%   spline is evaluated.
%
%   This routine uses fnval.

    requiredFields = {'T', 'ppU', 'ppV'};
    if ~isstruct(refInt) || ~all(isfield(refInt, requiredFields))
        error('Input refInt must contain the fields T, ppU, and ppV.');
    end

    if exist('fnval', 'file') ~= 2
        error('eval_ref_orbit requires fnval, which is not available.');
    end

    if ~isnumeric(tQuery) || ~isreal(tQuery)
        error('tQuery must be a real numeric array.');
    end

    T = refInt.T;
    if ~isscalar(T) || ~isfinite(T) || T <= 0
        error('refInt.T must be a positive finite scalar.');
    end

    % Wrap physical time to one period and convert to normalized phase.
    sQuery = mod(tQuery, T) / T;

    uRef = fnval(refInt.ppU, sQuery);
    vRef = fnval(refInt.ppV, sQuery);

    Uref = [uRef(:).'; vRef(:).'];
end
