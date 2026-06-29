function refInt = build_ref_orbit_interpolant(ref)
%BUILD_REF_ORBIT_INTERPOLANT  Build a periodic spline for an ODE orbit.
%
%   refInt = build_ref_orbit_interpolant(ref)
%
%   converts one saved reference cycle of the spatially homogeneous ODE
%   problem into a periodic interpolant.  The resulting structure is used by
%   the homogeneous limit-cycle postprocessing routines to evaluate the ODE
%   reference orbit at arbitrary PDE output times.
%
%   Input
%   -----
%   ref is the structure produced by compute_reference_orbit and must contain
%
%       ref.T        period of the reference orbit
%       ref.tCycle   time samples over one cycle, including the endpoint T
%       ref.uCycle   prey component over one cycle
%       ref.vCycle   predator component over one cycle
%
%   Output
%   ------
%   refInt is a structure with fields
%
%       refInt.T     period
%       refInt.s     normalized phase samples in [0,1]
%       refInt.u     prey samples on refInt.s
%       refInt.v     predator samples on refInt.s
%       refInt.ppU   periodic cubic spline for u(s)
%       refInt.ppV   periodic cubic spline for v(s)
%
%   The optional fields ref.lambda and ref.closureErr are copied when they
%   are present.
%
%   Notes
%   -----
%   The input cycle is assumed to include the duplicated endpoint at t = T.
%   This endpoint is removed before the normalized phase grid is closed again
%   at s = 1 for periodic spline construction.
%
%   This routine uses csape(...,'periodic').

    requiredFields = {'T', 'tCycle', 'uCycle', 'vCycle'};
    if ~isstruct(ref) || ~all(isfield(ref, requiredFields))
        error(['Input ref must be a structure containing the fields ', ...
               'T, tCycle, uCycle, and vCycle.']);
    end

    if exist('csape', 'file') ~= 2
        error(['build_ref_orbit_interpolant requires ', ...
               'csape(...,''periodic''), which is not available.']);
    end

    T = ref.T;
    t = ref.tCycle(:).';
    u = ref.uCycle(:).';
    v = ref.vCycle(:).';

    if ~isscalar(T) || ~isfinite(T) || T <= 0
        error('ref.T must be a positive finite scalar.');
    end

    if numel(t) < 3 || numel(u) ~= numel(t) || numel(v) ~= numel(t)
        error('ref.tCycle, ref.uCycle, and ref.vCycle must have matching lengths >= 3.');
    end

    % Remove the duplicated endpoint at t = T.
    t = t(1:end-1);
    u = u(1:end-1);
    v = v(1:end-1);

    % Convert physical time to normalized phase.
    s = t / T;

    % Close the phase interval explicitly for periodic spline construction.
    sFull = [s, 1];
    uFull = [u, u(1)];
    vFull = [v, v(1)];

    refInt = struct();
    refInt.T   = T;
    refInt.s   = sFull;
    refInt.u   = uFull;
    refInt.v   = vFull;
    refInt.ppU = csape(sFull, uFull, 'periodic');
    refInt.ppV = csape(sFull, vFull, 'periodic');

    if isfield(ref, 'lambda')
        refInt.lambda = ref.lambda;
    end

    if isfield(ref, 'closureErr')
        refInt.closureErr = ref.closureErr;
    end
end
