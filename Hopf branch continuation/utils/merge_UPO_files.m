function merge_UPO_files(varargin)
%MERGE_UPO_FILES  Merge saved UPO files into one branch master file.
%
%   merge_UPO_files()
%   merge_UPO_files('filePattern', 'UPO_Ex2.6.2_mode2_thetaFixed_lambda*.mat')
%   merge_UPO_files('outName', 'UPOMaster_Ex2.6.2_mode2_thetaFixed.mat')
%
% The function scans the workflow data folder, loads selected UPO_*.mat
% files, sorts the orbits by lambda, removes duplicate lambda values up to
% a tolerance, and saves a UPOMaster_*.mat file.
%
% Expected variables in each UPO file:
%   UV, T, lambda, hopf
%
% Optional variables:
%   resNorm, ampUV

opts = default_options();
opts = parse_options(opts, varargin{:});

rootDir = workflow_root();
dataDir = fullfile(rootDir, 'data');

files = dir(fullfile(dataDir, opts.filePattern));
if isempty(files)
    fprintf('No files matching %s found in %s\n', opts.filePattern, dataDir);
    return;
end

fprintf('Merging %d UPO files from %s\n', numel(files), dataDir);

UPOs = empty_upo_array();
sourceFiles = cell(0, 1);

for j = 1:numel(files)
    filePath = fullfile(dataDir, files(j).name);
    S = load(filePath);

    if ~has_required_upo_fields(S)
        warning('merge_UPO_files:SkippingFile', ...
            'Skipping %s because it does not contain UV, T, lambda, and hopf.', ...
            files(j).name);
        continue;
    end

    entry = make_upo_entry(S);
    UPOs(end + 1, 1) = entry; %#ok<AGROW>
    sourceFiles{end + 1, 1} = files(j).name; %#ok<AGROW>
end

if isempty(UPOs)
    warning('merge_UPO_files:NoEntries', 'No valid UPO entries were merged.');
    return;
end

[UPOs, sourceFiles] = sort_upos_by_lambda(UPOs, sourceFiles);
[UPOs, sourceFiles] = remove_duplicate_lambdas(UPOs, sourceFiles, opts.tol);

lambdas  = column_vector([UPOs.lambda]);
periods  = column_vector([UPOs.T]);
resNorms = collect_optional_numeric(UPOs, 'resNorm');
ampUVs   = collect_optional_numeric(UPOs, 'ampUV');
metadata = make_metadata(opts, dataDir, sourceFiles);

outPath = fullfile(dataDir, opts.outName);
save(outPath, 'UPOs', 'lambdas', 'periods', 'resNorms', ...
    'ampUVs', 'sourceFiles', 'metadata', '-v7.3');

fprintf('Saved %d merged UPO entries to %s\n', numel(UPOs), outPath);
end

function opts = default_options()
opts.filePattern = 'UPO_Ex2.6.2_mode2_thetaFixed_lambda*.mat';
opts.outName = 'UPOMaster_Ex2.6.2_mode2_thetaFixed.mat';
opts.tol = 1e-10;
end

function opts = parse_options(opts, varargin)
if mod(numel(varargin), 2) ~= 0
    error('merge_UPO_files:InvalidOptions', ...
        'Options must be supplied as name/value pairs.');
end

for j = 1:2:numel(varargin)
    name = lower(varargin{j});
    value = varargin{j + 1};

    switch name
        case 'filepattern'
            validate_char_option(name, value);
            opts.filePattern = value;
        case 'outname'
            validate_char_option(name, value);
            opts.outName = value;
        case 'tol'
            if ~isnumeric(value) || ~isscalar(value) || value < 0
                error('merge_UPO_files:InvalidTolerance', ...
                    'The duplicate tolerance must be a nonnegative scalar.');
            end
            opts.tol = value;
        otherwise
            error('merge_UPO_files:UnknownOption', ...
                'Unknown option: %s', varargin{j});
    end
end
end

function validate_char_option(name, value)
if ~ischar(value) || isempty(value)
    error('merge_UPO_files:InvalidOption', ...
        'The %s option must be a nonempty character vector.', name);
end
end

function rootDir = workflow_root()
thisFile = mfilename('fullpath');
postprocessDir = fileparts(thisFile);
rootDir = fileparts(postprocessDir);
end

function tf = has_required_upo_fields(S)
tf = isfield(S, 'UV') && isfield(S, 'T') && ...
    isfield(S, 'lambda') && isfield(S, 'hopf');
end

function UPOs = empty_upo_array()
UPOs = repmat(struct( ...
    'UV', [], ...
    'T', [], ...
    'lambda', [], ...
    'resNorm', [], ...
    'ampUV', [], ...
    'hopf', []), 0, 1);
end

function entry = make_upo_entry(S)
entry = struct();
entry.UV = S.UV;
entry.T = S.T;
entry.lambda = S.lambda;
entry.hopf = S.hopf;

entry.resNorm = [];
if isfield(S, 'resNorm')
    entry.resNorm = S.resNorm;
end

entry.ampUV = [];
if isfield(S, 'ampUV')
    entry.ampUV = S.ampUV;
end
end

function [UPOs, sourceFiles] = sort_upos_by_lambda(UPOs, sourceFiles)
[~, idx] = sort([UPOs.lambda]);
UPOs = UPOs(idx);
sourceFiles = sourceFiles(idx);
end

function [UPOs, sourceFiles] = remove_duplicate_lambdas(UPOs, sourceFiles, tol)
if numel(UPOs) <= 1
    return;
end

keep = true(numel(UPOs), 1);
lambdas = [UPOs.lambda];

for j = 2:numel(UPOs)
    if abs(lambdas(j) - lambdas(j - 1)) <= tol
        prevScore = residual_score(UPOs(j - 1));
        currScore = residual_score(UPOs(j));

        if currScore <= prevScore
            keep(j - 1) = false;
        else
            keep(j) = false;
        end
    end
end

UPOs = UPOs(keep);
sourceFiles = sourceFiles(keep);
end

function score = residual_score(entry)
if isfield(entry, 'resNorm') && isnumeric(entry.resNorm) && ~isempty(entry.resNorm)
    score = min(entry.resNorm(:));
else
    score = inf;
end
end

function values = collect_optional_numeric(UPOs, fieldName)
values = nan(numel(UPOs), 1);
for j = 1:numel(UPOs)
    if isfield(UPOs(j), fieldName) && isnumeric(UPOs(j).(fieldName)) && ...
            ~isempty(UPOs(j).(fieldName))
        fieldValue = UPOs(j).(fieldName);
        values(j) = fieldValue(1);
    end
end
end

function v = column_vector(v)
v = v(:);
end

function metadata = make_metadata(opts, dataDir, sourceFiles)
metadata = struct();
metadata.created = datestr(now, 30);
metadata.dataDir = dataDir;
metadata.filePattern = opts.filePattern;
metadata.outName = opts.outName;
metadata.duplicateTol = opts.tol;
metadata.numEntries = numel(sourceFiles);
metadata.sourceFiles = sourceFiles;
end
