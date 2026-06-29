function continue_branch(fileName, nSteps, dLambda, keepInitial)
%CONTINUE_BRANCH  Continue a periodic-orbit branch from a saved anchor orbit.
%
%   continue_branch(fileName, nSteps, dLambda)
%   continue_branch(fileName, nSteps, dLambda, keepInitial)
%
%   fileName    Name of a .mat file in the workflow data/ folder, or a full
%               path to such a file. The file must contain UV, T, hopf, and
%               resNorm, as produced by compute_anchor_orbit/newton_solver.
%   nSteps      Number of continuation steps.
%   dLambda     Parameter increment used at each step.
%   keepInitial If true, include the loaded anchor orbit as the first entry
%               in the saved branch structure. The default is false.
%
%   The output is saved as a Branch_*.mat file in the workflow data/ folder.

    if nargin < 4 || isempty(keepInitial)
        keepInitial = false;
    end
    check_inputs(fileName, nSteps, dLambda);

    rootDir = workflow_root();
    dataDir = fullfile(rootDir, 'data');
    inputFile = resolve_input_file(fileName, dataDir);

    data = load_anchor_orbit(inputFile);

    UV      = data.UV;
    T       = data.T;
    hopf    = data.hopf;
    resNorm = data.resNorm;
    ampUV   = data.ampUV;

    lambdaStart = hopf.cfg.lambda;
    lambda      = lambdaStart;

    UPOs = initialize_storage(nSteps, keepInitial);
    idx0 = 0;

    if keepInitial
        UPOs = record_step(UPOs, 1, UV, T, lambda, resNorm, ampUV, hopf);
        idx0 = 1;
    end

    fprintf('Starting continuation from lambda = %.12g\n', lambdaStart);

    for step = 1:nSteps
        lambdaNext = lambda + dLambda;
        fprintf('\nStep %d/%d: lambda %.12g -> %.12g\n', ...
            step, nSteps, lambda, lambdaNext);

        lambda = lambdaNext;
        hopf = update_model_params(hopf, lambda);

        hopf.phaseRefUV = UV;
        saveLastOrbit = (step == nSteps);
        [UV, T, resNorm, ampUV] = newton_solver(UV, T, hopf, saveLastOrbit);
        hopf.T = T;

        if isfield(hopf, 'phaseRefUV')
            hopf = rmfield(hopf, 'phaseRefUV');
        end

        UPOs = record_step(UPOs, step + idx0, UV, T, lambda, resNorm, ampUV, hopf);
    end

    outputFile = save_branch_output(UPOs, rootDir, lambdaStart, nSteps, dLambda);

    fprintf('\nContinuation completed.\n');
    fprintf('Branch saved to %s\n', outputFile);
end

function check_inputs(fileName, nSteps, dLambda)
%CHECK_INPUTS  Basic input validation for the continuation driver.

    if ~ischar(fileName) || isempty(fileName)
        error('continue_branch:InvalidFileName', ...
            'fileName must be a nonempty character vector.');
    end

    if ~isnumeric(nSteps) || ~isscalar(nSteps) || nSteps < 1 || nSteps ~= round(nSteps)
        error('continue_branch:InvalidNSteps', ...
            'nSteps must be a positive integer.');
    end

    if ~isnumeric(dLambda) || ~isscalar(dLambda) || dLambda == 0
        error('continue_branch:InvalidDLambda', ...
            'dLambda must be a nonzero numeric scalar.');
    end
end

function data = load_anchor_orbit(inputFile)
%LOAD_ANCHOR_ORBIT  Load and validate the saved orbit used to start continuation.

    raw = load(inputFile);

    require_field(raw, 'UV', inputFile);
    require_field(raw, 'T', inputFile);
    require_field(raw, 'hopf', inputFile);
    require_field(raw, 'resNorm', inputFile);

    data.UV      = raw.UV;
    data.T       = raw.T;
    data.hopf    = raw.hopf;
    data.resNorm = raw.resNorm;

    if isfield(raw, 'ampUV')
        data.ampUV = raw.ampUV;
    else
        data.ampUV = NaN;
    end

    if ~isfield(data.hopf, 'cfg')
        error('continue_branch:MissingConfig', ...
            'The loaded hopf structure in %s must contain hopf.cfg.', inputFile);
    end

    require_route_fun(data.hopf.cfg);

    if ~isfield(data.hopf.cfg, 'lambda')
        error('continue_branch:MissingLambda', ...
            'The loaded configuration in %s must contain cfg.lambda.', inputFile);
    end
end

function require_field(s, fieldName, inputFile)
%REQUIRE_FIELD  Require a variable in a loaded MAT-file structure.

    if ~isfield(s, fieldName)
        error('continue_branch:MissingVariable', ...
            'The file %s must contain the variable %s.', inputFile, fieldName);
    end
end

function require_route_fun(cfg)
%REQUIRE_ROUTE_FUN  Require the parameter-route function in the configuration.

    if ~isfield(cfg, 'route') || ...
            ~isfield(cfg.route, 'fun') || isempty(cfg.route.fun)
        error('continue_branch:MissingRouteFun', ...
            'The configuration must define cfg.route.fun.');
    end
end

function hopf = update_model_params(hopf, lambda)
%UPDATE_MODEL_PARAMS  Update lambda-dependent parameters along the route.

    cfg = hopf.cfg;
    require_route_fun(cfg);

    [m, theta] = cfg.route.fun(lambda);

    cfg.lambda = lambda;
    cfg.m      = m;
    cfg.theta  = theta;
    cfg.u_star = lambda;
    cfg.v_star = lambda * (1 - lambda/cfg.k) / theta;

    hopf.cfg = cfg;
end

function UPOs = initialize_storage(nSteps, keepInitial)
%INITIALIZE_STORAGE  Allocate the branch output structure.

    nSlots = nSteps;
    if keepInitial
        nSlots = nSlots + 1;
    end

    template = struct( ...
        'UV',      [], ...
        'lambda',  [], ...
        'T',       [], ...
        'resNorm', [], ...
        'ampUV',   [], ...
        'hopf',    []);

    UPOs = repmat(template, nSlots, 1);
end

function UPOs = record_step(UPOs, k, UV, T, lambda, resNorm, ampUV, hopf)
%RECORD_STEP  Store one corrected periodic orbit in the branch structure.

    UPOs(k).UV      = UV;
    UPOs(k).T       = T;
    UPOs(k).lambda  = lambda;
    UPOs(k).resNorm = resNorm;
    UPOs(k).ampUV   = ampUV;
    UPOs(k).hopf    = hopf;
end

function outputFile = save_branch_output(UPOs, rootDir, lambdaStart, nSteps, dLambda)
%SAVE_BRANCH_OUTPUT  Save the computed branch to the workflow data/ folder.

    dataDir = fullfile(rootDir, 'data');
    if ~exist(dataDir, 'dir')
        mkdir(dataDir);
    end

    cfg = UPOs(1).hopf.cfg;
    routeTag = route_tag(cfg);

    outputName = sprintf('Branch_%s_mode%d_%s_lambda%g_%dSteps_dL%g.mat', ...
        cfg.caseName, UPOs(1).hopf.mode, routeTag, lambdaStart, nSteps, dLambda);
    outputFile = fullfile(dataDir, outputName);

    if exist(outputFile, 'file') == 2
        warning('continue_branch:OverwriteBranchFile', ...
            'File %s already exists. Overwriting.', outputFile);
    end

    save(outputFile, 'UPOs', '-v7.3');
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

function inputFile = resolve_input_file(fileName, dataDir)
%RESOLVE_INPUT_FILE  Resolve a MAT-file path for a saved periodic orbit.

    if exist(fileName, 'file') == 2
        inputFile = fileName;
        return;
    end

    inputFile = fullfile(dataDir, fileName);
    if exist(inputFile, 'file') == 2
        return;
    end

    [pathPart, baseName, extPart] = fileparts(fileName);
    if isempty(extPart)
        candidate = fullfile(dataDir, [fileName '.mat']);
        if exist(candidate, 'file') == 2
            inputFile = candidate;
            return;
        end

        if ~isempty(pathPart)
            candidate = fullfile(pathPart, [baseName '.mat']);
            if exist(candidate, 'file') == 2
                inputFile = candidate;
                return;
            end
        end
    end

    error('continue_branch:FileNotFound', ...
        'Could not find the input file %s.', fileName);
end

function rootDir = workflow_root()
%WORKFLOW_ROOT  Return the top-level workflow folder.

    thisFile = mfilename('fullpath');
    thisDir = fileparts(thisFile);
    rootDir = fileparts(thisDir);
end
