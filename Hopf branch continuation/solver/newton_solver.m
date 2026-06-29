function [UV, T, resNorm, ampUV] = newton_solver(UV, T, hopf, doSave)
%NEWTON_SOLVER  Newton correction for a time-periodic orbit.
%
%   [UV, T, resNorm, ampUV] = NEWTON_SOLVER(UV, T, hopf, doSave)
%   corrects an initial periodic-orbit guess using the residual and Jacobian
%   assembled by residual_and_jacobian.m.
%
%   UV is the flattened vector [U(:); V(:)], where U and V are Nx-by-Nt.
%   T is the orbit period. The structure hopf supplies the configuration,
%   phase data, and current period used by the residual routine.
%
%   If doSave is true, the corrected orbit is saved to the workflow data/
%   folder after convergence.

    if nargin < 4
        doSave = false;
    end

    opts = newton_options(hopf);
    hopf.T = T;

    converged = false;
    deltaNorm = Inf;

    for iter = 1:opts.maxIter
        [R, J] = residual_and_jacobian(UV, hopf);

        delta = sparse(J) \ (-R);
        dUV = delta(1:end-1);
        dT  = delta(end);

        UV = UV + opts.alpha * dUV;
        T  = T  + opts.alpha * dT;
        hopf.T = T;

        R = residual_and_jacobian(UV, hopf);

        resNorm   = norm(R);
        deltaNorm = norm(delta);

        if opts.verbose
            resPDE = norm(R(1:end-1));
            resPHI = abs(R(end));
            fprintf('Iter %2d: ||R|| = %.3e, ||delta|| = %.3e\n', ...
                iter, resNorm, deltaNorm);
            fprintf('Iter %2d: ||R_PDE|| = %.3e, |Phi| = %.3e\n', ...
                iter, resPDE, resPHI);
        end

        if deltaNorm < opts.stepTol || resNorm < opts.resTol
            converged = true;
            break;
        end
    end

    R = residual_and_jacobian(UV, hopf);
    resNorm = norm(R);

    ampUV = spatial_amplitude_UV( ...
        UV, ...
        hopf.cfg.periodic_orbit.Nx, ...
        hopf.cfg.periodic_orbit.Nt, ...
        opts.amplitudeType, ...
        hopf.cfg.u_star, ...
        hopf.cfg.v_star);

    if ~converged && resNorm < opts.resTol
        converged = true;
    end

    if converged && doSave
        save_UPO_data(UV, hopf, resNorm, ampUV);
    elseif ~converged && opts.verbose
        warning('newton_solver:NoConvergence', ...
            'Newton iteration stopped with ||R|| = %.3e and ||delta|| = %.3e.', ...
            resNorm, deltaNorm);
    end
end

function opts = newton_options(hopf)
%NEWTON_OPTIONS  Numerical options for Newton correction.

    opts.maxIter = 10;
    opts.stepTol = 1e-10;
    opts.resTol  = 1e-10;
    opts.alpha   = 1.0;
    opts.verbose = true;

    if isfield(hopf, 'cfg') && isfield(hopf.cfg, 'newton')
        nopts = hopf.cfg.newton;
        if isfield(nopts, 'maxIter'), opts.maxIter = nopts.maxIter; end
        if isfield(nopts, 'stepTol'), opts.stepTol = nopts.stepTol; end
        if isfield(nopts, 'resTol'),  opts.resTol  = nopts.resTol;  end
        if isfield(nopts, 'alpha'),   opts.alpha   = nopts.alpha;   end
        if isfield(nopts, 'verbose'), opts.verbose = nopts.verbose; end
    end

    opts.amplitudeType = 2;
    if isfield(hopf, 'cfg') && ...
            isfield(hopf.cfg, 'periodic_orbit') && ...
            isfield(hopf.cfg.periodic_orbit, 'amplitudeType')
        opts.amplitudeType = hopf.cfg.periodic_orbit.amplitudeType;
    end
end

function amp = spatial_amplitude_UV(UV, Nx, Nt, amplitudeType, u_star, v_star)
%SPATIAL_AMPLITUDE_UV  Spatial nonhomogeneity amplitude for both species.
%
%   amplitudeType = 1: absolute RMS spatial variation over time.
%   amplitudeType = 2: RMS spatial variation normalized by equilibrium scale.
%   amplitudeType = 3: RMS spatial variation normalized by spatial-mean scale.

    if nargin < 4 || isempty(amplitudeType)
        amplitudeType = 2;
    end

    eps0 = 1e-14;
    Ntot = Nx * Nt;

    if numel(UV) ~= 2*Ntot
        error('newton_solver:UVLengthMismatch', ...
            'UV length mismatch: expected %d, got %d.', ...
            2*Ntot, numel(UV));
    end

    U = reshape(UV(1:Ntot), Nx, Nt);
    V = reshape(UV(Ntot+1:end), Nx, Nt);

    Ubar = mean(U, 1);
    Vbar = mean(V, 1);

    Uvar = U - repmat(Ubar, Nx, 1);
    Vvar = V - repmat(Vbar, Nx, 1);

    energy = sum(Uvar.^2, 1) + sum(Vvar.^2, 1);
    numerator = sqrt(mean(energy));

    ampType = parse_amplitude_type(amplitudeType);

    switch ampType
        case 1
            amp = numerator;

        case 2
            if nargin < 6 || isempty(u_star) || isempty(v_star)
                u_star = mean(Ubar);
                v_star = mean(Vbar);
            end
            denom = sqrt(Nx) * sqrt(u_star^2 + v_star^2) + eps0;
            amp = numerator / denom;

        case 3
            ubarRMS = sqrt(mean(Ubar.^2));
            vbarRMS = sqrt(mean(Vbar.^2));
            denom = sqrt(Nx) * sqrt(ubarRMS^2 + vbarRMS^2) + eps0;
            amp = numerator / denom;
    end
end

function ampType = parse_amplitude_type(amplitudeType)
%PARSE_AMPLITUDE_TYPE  Convert amplitude-type labels to numeric codes.

    if ischar(amplitudeType)
        key = lower(strtrim(amplitudeType));
        if any(strcmp(key, {'1', 'rms', 'abs', 'absolute'}))
            ampType = 1;
        elseif any(strcmp(key, {'2', 'rel_eq', 'eq', 'equilibrium'}))
            ampType = 2;
        elseif any(strcmp(key, {'3', 'rel_mean', 'mean', 'time_mean'}))
            ampType = 3;
        else
            error('newton_solver:UnknownAmplitudeType', ...
                'Unknown amplitudeType string: %s.', amplitudeType);
        end
    else
        ampType = amplitudeType;
    end

    if ~any(ampType == [1, 2, 3])
        error('newton_solver:UnknownAmplitudeType', ...
            'Unknown amplitudeType: %s.', num2str(amplitudeType));
    end
end

function save_UPO_data(UV, hopf, resNorm, ampUV)
%SAVE_UPO_DATA  Save a corrected periodic orbit for continuation.

    rootDir = workflow_root();
    dataDir = fullfile(rootDir, 'data');
    if ~exist(dataDir, 'dir')
        mkdir(dataDir);
    end

    lambda = hopf.cfg.lambda;
    routeTag = route_tag(hopf.cfg);
    fname = sprintf('UPO_%s_mode%d_%s_lambda%g.mat', ...
        hopf.cfg.caseName, hopf.mode, routeTag, lambda);

    T = hopf.T;
    save(fullfile(dataDir, fname), ...
        'UV', 'T', 'lambda', 'resNorm', 'hopf', 'ampUV');

    fprintf('UPO saved to %s\n', fullfile(dataDir, fname));
end

function rootDir = workflow_root()
%WORKFLOW_ROOT  Return the top-level workflow folder.

    thisFile = mfilename('fullpath');
    thisDir = fileparts(thisFile);
    rootDir = fileparts(thisDir);
end

function routeTag = route_tag(cfg)
%ROUTE_TAG  Descriptive route tag used in saved filenames.

    if isfield(cfg, 'route') && isfield(cfg.route, 'type')
        routeTag = cfg.route.type;
    elseif isfield(cfg, 'route') && isfield(cfg.route, 'name')
        routeTag = cfg.route.name;
    else
        routeTag = 'thetaFixed';
    end
end
