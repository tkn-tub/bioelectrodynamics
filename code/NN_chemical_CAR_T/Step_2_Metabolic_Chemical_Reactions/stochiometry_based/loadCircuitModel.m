function [model, modelFile] = loadCircuitModel(fileName, modelVariable)
%LOADCIRCUITMODEL Find a MAT/SBML model with which(), then load without dialogs.
if nargin < 2, modelVariable = ''; end
fileName = char(fileName);
modelFile = which(fileName);
if isempty(modelFile) && exist(fileName, 'file') == 2
    modelFile = fileName;
end
if isempty(modelFile)
    error('Circuit:ModelNotFound', ...
        ['Model "%s" was not found. Call addingPathParentFolderByName(''code'') ' ...
         'from within your code tree, or provide the full model path.'], fileName);
end
% Convert a relative fallback to an absolute path.
if isempty(regexp(modelFile, '^(?:[A-Za-z]:[\\/]|/|\\\\)', 'once'))
    modelFile = fullfile(pwd, modelFile);
end
[~,~,ext] = fileparts(modelFile);
switch lower(ext)
    case '.mat'
        contents = load(modelFile);
        if ~isempty(modelVariable)
            if ~isfield(contents, modelVariable)
                error('Circuit:MATVariable', 'MAT file has no variable "%s".', modelVariable);
            end
            model = contents.(modelVariable);
        else
            names = fieldnames(contents);
            candidates = {};
            for k = 1:numel(names)
                v = contents.(names{k});
                if isstruct(v) && isscalar(v) && all(isfield(v, {'S','rxns','mets'}))
                    candidates{end+1} = names{k}; %#ok<AGROW>
                end
            end
            if numel(candidates) ~= 1
                error('Circuit:MATVariable', ...
                    'Found %d model structures. Set MODEL_VARIABLE explicitly.', numel(candidates));
            end
            model = contents.(candidates{1});
        end
    case {'.xml','.sbml'}
        if exist('readCbModel', 'file') ~= 2
            error('Circuit:COBRA', ...
                'SBML loading requires the COBRA Toolbox readCbModel function and SBML support.');
        end
        model = readCbModel(modelFile, 'fileType', 'SBML');
        % Retain EC annotations that an SBML importer may omit.
        model = addSBMLECAnnotations(model, modelFile);
    otherwise
        error('Circuit:ModelFormat', 'Use a .mat, .xml, or .sbml model.');
end
if ~isstruct(model) || ~isscalar(model) || ~all(isfield(model, {'S','rxns','mets'}))
    error('Circuit:ModelFormat', 'The loaded object is not a COBRA model structure.');
end
fprintf('Loaded model file: %s\n', modelFile);
end

function model = addSBMLECAnnotations(model, modelFile)
% Exact reaction-ID matching; no reaction-name inference or online lookup.
n = numel(model.rxns);
extra = repmat({''}, n, 1);
xml = fileread(modelFile);
blocks = regexp(xml, '(?s)<(?:\w+:)?reaction\b[^>]*>.*?</(?:\w+:)?reaction\s*>', 'match');
ids = model.rxns;
if ~iscell(ids), ids = cellstr(ids); end
lookup = containers.Map(ids(:), num2cell((1:n)'));
for b = 1:numel(blocks)
    opening = regexp(blocks{b}, '^<[^>]+>', 'match', 'once');
    tok = regexp(opening, '\bid\s*=\s*["'']([^"'']+)["'']', 'tokens', 'once');
    if isempty(tok), continue; end
    id = tok{1};
    candidates = {id, regexprep(id, '^R_', '')};
    idx = [];
    for k = 1:numel(candidates)
        if isKey(lookup, candidates{k}), idx = lookup(candidates{k}); break; end
    end
    if isempty(idx), continue; end
    ec = regexp(blocks{b}, ...
        '(?i)ec-code(?:/|%2F|:)([1-7]\.\d+\.\d+\.\d+)(?![\d.])', 'tokens');
    notes = regexp(blocks{b}, ...
        '(?i)\bEC(?:\s*(?:NUMBER[S]?|CODE[S]?))?\s*[:=]\s*([^<\r\n]+)', 'tokens');
    values = cellfun(@(x) x{1}, ec, 'UniformOutput', false);
    for k = 1:numel(notes)
        values = [values, regexp(notes{k}{1}, ...
            '(?<![\d.])[1-7]\.\d+\.\d+\.\d+(?![\d.])', 'match')]; %#ok<AGROW>
    end
    extra{idx} = strjoin(unique(values), '; ');
end
model.circuitSBMLEC = extra;
end
