function floquet = compute_floquet_IFRK4(UV, T, cfgOrbit, cfgFloq, runInfo)
%COMPUTE_FLOQUET_IFRK4  Compute Floquet multipliers for one stored orbit.
%
%   floquet = compute_floquet_IFRK4(UV, T, cfgOrbit, cfgFloq, runInfo)
%
%   This routine is the high-level single-orbit Floquet workflow. It builds
%   the explicit monodromy matrix with floquet_buildM_IFRK4, postprocesses
%   the resulting multipliers with floquet_processM, and returns a compact
%   output structure containing method information, orbit metadata, selected
%   multipliers/eigenvectors, stability counts, and quality-control flags.
%
%   Inputs
%   ------
%   UV       : packed periodic orbit vector of length 2*Nx*Nt.
%   T        : orbit period.
%   cfgOrbit : orbit/model configuration, usually hopf.cfg.
%   cfgFloq  : Floquet configuration produced by make_floqCfg.
%   runInfo  : optional metadata struct, for example branchID, mode,
%              routeType, and s.
%
%   Output
%   ------
%   floquet : struct with fields info, build, orbit, eig, summary, and qc.
%
%   MATLAB R2016a compatible.

    if nargin < 5 || isempty(runInfo)
        runInfo = struct();
    end

    validate_inputs(UV, T, cfgOrbit, cfgFloq, runInfo);

    optsBuild = struct();
    optsBuild.mSub = cfgFloq.mSub;
    optsBuild.verbose = cfgFloq.verbose;

    outM = floquet_buildM_IFRK4(UV, T, cfgOrbit, optsBuild);

    [eigOut, summary, qc, infoSpec] = ...
        floquet_processM(outM.M, T, cfgFloq, UV, cfgOrbit);

    floquet = struct();
    floquet.info = make_info_block(cfgFloq, runInfo);
    floquet.build = make_build_block(cfgFloq, outM, infoSpec);
    floquet.orbit = make_orbit_block(UV, T, cfgOrbit);
    floquet.eig = eigOut;
    floquet.summary = summary;
    floquet.qc = qc;

    if cfgFloq.verbose
        print_summary(floquet, cfgFloq);
    end
end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function validate_inputs(UV, T, cfgOrbit, cfgFloq, runInfo)
    if ~isvector(UV) || isempty(UV)
        error('compute_floquet_IFRK4:BadUV', ...
            'UV must be a nonempty vector.');
    end
    if ~isscalar(T) || ~isfinite(T) || T <= 0
        error('compute_floquet_IFRK4:BadPeriod', ...
            'T must be a finite positive scalar.');
    end
    if ~isstruct(cfgOrbit)
        error('compute_floquet_IFRK4:BadCfgOrbit', ...
            'cfgOrbit must be a struct.');
    end
    if ~isstruct(cfgFloq)
        error('compute_floquet_IFRK4:BadCfgFloq', ...
            'cfgFloq must be a struct.');
    end
    if ~isstruct(runInfo)
        error('compute_floquet_IFRK4:BadRunInfo', ...
            'runInfo must be a struct when supplied.');
    end

    requiredFloq = {'methodTagBase', 'mSub', 'verbose'};
    require_fields(cfgFloq, requiredFloq, 'cfgFloq');

    requiredOrbit = {'lambda', 'm', 'theta', 'k', 'u_star', 'v_star'};
    require_fields(cfgOrbit, requiredOrbit, 'cfgOrbit');

    if ~isfield(cfgOrbit, 'periodic_orbit') || ...
            ~isfield(cfgOrbit.periodic_orbit, 'Nx') || ...
            ~isfield(cfgOrbit.periodic_orbit, 'Nt')
        error('compute_floquet_IFRK4:MissingGridFields', ...
            'cfgOrbit.periodic_orbit.Nx and Nt are required.');
    end

    Nx = cfgOrbit.periodic_orbit.Nx;
    Nt = cfgOrbit.periodic_orbit.Nt;
    if numel(UV) ~= 2 * Nx * Nt
        error('compute_floquet_IFRK4:BadUVLength', ...
            'UV must have length 2*Nx*Nt.');
    end

    if ~isfield(cfgOrbit, 'L') && ~isfield(cfgOrbit, 'ell')
        error('compute_floquet_IFRK4:MissingDomainLength', ...
            'cfgOrbit must define cfgOrbit.L or cfgOrbit.ell.');
    end
end

function require_fields(S, fields, structName)
    for j = 1:numel(fields)
        field = fields{j};
        if ~isfield(S, field)
            error('compute_floquet_IFRK4:MissingField', ...
                '%s.%s is required.', structName, field);
        end
    end
end

function info = make_info_block(cfgFloq, runInfo)
    info = struct();
    info.methodTag = sprintf('%s-mSub%d', cfgFloq.methodTagBase, cfgFloq.mSub);
    info.timestamp = datestr(now, 30);
    info.branchID  = get_field(runInfo, 'branchID', '');
    info.mode      = get_field(runInfo, 'mode', NaN);
    info.routeType = get_field(runInfo, 'routeType', '');
    info.s         = get_field(runInfo, 's', NaN);
end

function build = make_build_block(cfgFloq, outM, infoSpec)
    build = struct();
    build.mSub = cfgFloq.mSub;
    build.dt = outM.diag.dt;
    build.nSteps = outM.diag.nSteps;
    build.elapsed_sec = outM.diag.elapsed_sec;
    build.q = infoSpec.q;
    build.n = infoSpec.n;
    build.eigTime_sec = infoSpec.eigTime_sec;
end

function orbit = make_orbit_block(UV, T, cfgOrbit)
    Nx = cfgOrbit.periodic_orbit.Nx;
    Nt = cfgOrbit.periodic_orbit.Nt;
    Lx = get_domain_length(cfgOrbit);

    orbit = struct();
    orbit.lambda = cfgOrbit.lambda;
    orbit.T = T;

    params = struct();
    params.m = cfgOrbit.m;
    params.theta = cfgOrbit.theta;
    params.k = cfgOrbit.k;
    params.ustar = cfgOrbit.u_star;
    params.vstar = cfgOrbit.v_star;
    params.s = get_route_s(cfgOrbit);
    orbit.params = params;

    orbit.Nx = Nx;
    orbit.Nt = Nt;
    orbit.domain = struct('Lx', Lx, 'type', '(0,Lx) with Neumann BC');

    orbit.grid = struct();
    orbit.grid.type = 'DCT-I (Neumann cosine) grid';
    orbit.grid.x = linspace(0, Lx, Nx).';
    orbit.grid.t = (0:Nt-1).' * (T / Nt);

    orbit.packing = '[u(:,1..Nt); v(:,1..Nt)] on x-grid, no t=T slice';
    orbit.UV = UV(:);
end

function Lx = get_domain_length(cfgOrbit)
    if isfield(cfgOrbit, 'L') && ~isempty(cfgOrbit.L)
        Lx = cfgOrbit.L;
    else
        Lx = cfgOrbit.ell * pi;
    end
end

function s = get_route_s(cfgOrbit)
    if isfield(cfgOrbit, 'route') && isfield(cfgOrbit.route, 's')
        s = cfgOrbit.route.s;
    else
        s = NaN;
    end
end

function v = get_field(S, name, default)
    if isfield(S, name) && ~isempty(S.(name))
        v = S.(name);
    else
        v = default;
    end
end

function print_summary(floquet, cfgFloq)
    fprintf('--- compute_floquet_IFRK4 done ---\n');
    fprintf('methodTag: %s\n', floquet.info.methodTag);
    fprintf('muTrivial (q-aligned) = %.16g%+.16gi\n', ...
        real(floquet.eig.muTrivial), imag(floquet.eig.muTrivial));
    fprintf('nUnstable = %d (muTol=%.1e), gapUnitCircle=%.3e\n', ...
        floquet.summary.nUnstable, cfgFloq.muTol, ...
        floquet.summary.gapUnitCircle);
    fprintf('QC ok = %d\n', floquet.qc.ok);

    if ~floquet.qc.ok
        fprintf('QC reasons:\n');
        for j = 1:numel(floquet.qc.reason)
            fprintf('  - %s\n', floquet.qc.reason{j});
        end
    end
end
