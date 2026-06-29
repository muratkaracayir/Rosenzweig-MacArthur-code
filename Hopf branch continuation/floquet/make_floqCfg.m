function cfg = make_floqCfg(overrides)
%MAKE_FLOQCFG  Create options for single-orbit Floquet analysis.
%
%   cfg = make_floqCfg()
%   cfg = make_floqCfg(overrides)
%
%   The returned struct configures the explicit-monodromy IFRK4 Floquet
%   workflow used by compute_floquet_IFRK4.  The optional input overrides is
%   a struct whose fields replace the defaults below.  This is a shallow
%   merge: nested structs, if supplied, are replaced as whole fields.
%
%   Example:
%       cfgFloq = make_floqCfg(struct('mSub', 8, ...
%                                    'nStore', 20, ...
%                                    'muTol', 1e-5));
%
%   MATLAB R2016a compatible.

    if nargin < 1 || isempty(overrides)
        overrides = struct();
    end
    if ~isstruct(overrides)
        error('make_floqCfg:InvalidInput', ...
            'Input overrides must be a struct.');
    end

    cfg = default_floquet_config();
    cfg = apply_overrides(cfg, overrides);
    validate_cfg(cfg);

    % Normalize logical flags after accepting numeric 0/1 overrides.
    cfg.verbose = logical(cfg.verbose);
    cfg.storeEigenvectors = logical(cfg.storeEigenvectors);
end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function cfg = default_floquet_config()
%DEFAULT_FLOQUET_CONFIG  Defaults for explicit-monodromy IFRK4 analysis.

    cfg = struct();

    % Method provenance.
    cfg.methodTagBase = 'IFRK4-explicitM-eig';
    cfg.mSub = 8;             % IFRK4 substeps per stored orbit time step.
    cfg.verbose = true;

    % Eigenvector storage policy.
    cfg.storeEigenvectors = true;
    cfg.nStore = 16;          % Leading multipliers by modulus are stored.

    % Stability and quality-control tolerances.
    cfg.muTol = 1e-6;                 % unstable if |mu| > 1 + muTol
    cfg.near1Tol = 1e-3;              % nontrivial multiplier near unit circle
    cfg.maxMuTrivialDist = 1e-3;      % distance of trivial multiplier from 1
    cfg.rTrivTol = 1e-3;              % neutral-direction residual threshold
    cfg.eigResTol = 1e-6;             % full eig residual threshold
    cfg.trivialAgreeTol = 1e-8;       % q-aligned vs closest-to-one agreement
end

function cfg = apply_overrides(cfg, overrides)
%APPLY_OVERRIDES  Shallow-merge user-supplied options.

    names = fieldnames(overrides);
    for j = 1:numel(names)
        name = names{j};
        cfg.(name) = overrides.(name);
    end
end

function validate_cfg(cfg)
%VALIDATE_CFG  Basic validation for the Floquet option struct.

    require_nonempty_char(cfg, 'methodTagBase');

    require_positive_integer(cfg, 'mSub');
    require_logical_scalar(cfg, 'verbose');

    require_logical_scalar(cfg, 'storeEigenvectors');
    require_positive_integer(cfg, 'nStore');

    require_positive_scalar(cfg, 'muTol');
    require_positive_scalar(cfg, 'near1Tol');
    require_positive_scalar(cfg, 'maxMuTrivialDist');
    require_positive_scalar(cfg, 'rTrivTol');
    require_positive_scalar(cfg, 'eigResTol');
    require_positive_scalar(cfg, 'trivialAgreeTol');
end

function require_positive_scalar(cfg, field)
    if ~isfield(cfg, field) || ~isscalar(cfg.(field)) || ...
            ~isfinite(cfg.(field)) || cfg.(field) <= 0
        error('make_floqCfg:InvalidValue', ...
            'Field "%s" must be a finite positive scalar.', field);
    end
end

function require_positive_integer(cfg, field)
    if ~isfield(cfg, field)
        error('make_floqCfg:MissingField', ...
            'Missing field "%s".', field);
    end

    value = cfg.(field);
    if ~isscalar(value) || ~isfinite(value) || value <= 0 || ...
            fix(value) ~= value
        error('make_floqCfg:InvalidValue', ...
            'Field "%s" must be a finite positive integer.', field);
    end
end

function require_nonempty_char(cfg, field)
    if ~isfield(cfg, field) || ~ischar(cfg.(field)) || isempty(cfg.(field))
        error('make_floqCfg:InvalidValue', ...
            'Field "%s" must be a nonempty char array.', field);
    end
end

function require_logical_scalar(cfg, field)
    if ~isfield(cfg, field) || ~isscalar(cfg.(field))
        error('make_floqCfg:InvalidValue', ...
            'Field "%s" must be a logical scalar.', field);
    end

    value = cfg.(field);
    if islogical(value)
        return;
    end
    if isnumeric(value) && isfinite(value) && (value == 0 || value == 1)
        return;
    end

    error('make_floqCfg:InvalidValue', ...
        'Field "%s" must be a logical scalar.', field);
end
