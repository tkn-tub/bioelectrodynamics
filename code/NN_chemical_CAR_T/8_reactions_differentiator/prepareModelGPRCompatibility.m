function [model,hasGPR,gprSource] = prepareModelGPRCompatibility(model)
%PREPAREMODELGPRCOMPATIBILITY Normalize GPR fields for the motif mappers.
%
%   [model,hasGPR,gprSource] = prepareModelGPRCompatibility(model)
%
% Some SBML models loaded by readCbModel contain model.rules and model.genes
% but do not contain the optional human-readable model.grRules field. The
% differentiator, counter and NN mapping scripts use model.grRules both for
% GPR filtering and for reporting. This function creates a consistent
% nRxns-by-1 cell array in model.grRules and returns the corresponding
% logical hasGPR vector.
%
% Priority:
%   1) normalize an existing model.grRules field;
%   2) convert exact Boolean model.rules expressions using model.genes;
%   3) reconstruct gene sets from model.rxnGeneMat (Boolean logic unknown);
%   4) create empty rules and warn when the model has no gene associations.

if ~isstruct(model) || ~isfield(model,'rxns')
    error('GPRCompatibility:InvalidModel', ...
        'The input must be a COBRA model containing model.rxns.');
end

nRxns = numel(model.rxns);
normalizedGrRules = repmat({''},nRxns,1);
gprSource = "none";

%% Case 1: normalize an existing grRules field
if isfield(model,'grRules') && ~isempty(model.grRules)
    sourceRules = localTextColumn(model.grRules,'model.grRules');
    nCopy = min(nRxns,numel(sourceRules));
    normalizedGrRules(1:nCopy) = sourceRules(1:nCopy);
    gprSource = "model.grRules";

    if numel(sourceRules) ~= nRxns
        warning('GPRCompatibility:GrRulesLengthMismatch', ...
            ['model.grRules contains %d entries for %d reactions. ', ...
             'Missing entries were padded with empty rules.'], ...
            numel(sourceRules),nRxns);
    end

%% Case 2: reconstruct readable rules from rules + genes
elseif isfield(model,'rules') && ~isempty(model.rules) && ...
        isfield(model,'genes') && ~isempty(model.genes)

    machineRules = localTextColumn(model.rules,'model.rules');
    geneIds = localTextColumn(model.genes,'model.genes');

    if numel(machineRules) ~= nRxns
        warning('GPRCompatibility:RulesLengthMismatch', ...
            ['model.rules contains %d entries for %d reactions. ', ...
             'Missing reaction rules were treated as empty.'], ...
            numel(machineRules),nRxns);
    end

    nRulesToConvert = min(nRxns,numel(machineRules));
    for r = 1:nRulesToConvert
        expression = strtrim(machineRules{r});
        if isempty(expression)
            continue;
        end

        tokens = regexp(expression,'x\(([0-9]+)\)','tokens');
        if ~isempty(tokens)
            geneIndices = unique(cellfun( ...
                @(token) str2double(token{1}),tokens));

            invalidGeneIndices = geneIndices( ...
                ~isfinite(geneIndices) | geneIndices < 1 | ...
                geneIndices > numel(geneIds) | ...
                geneIndices ~= floor(geneIndices));
            if ~isempty(invalidGeneIndices)
                error('GPRCompatibility:InvalidGeneIndex', ...
                    ['Reaction %d contains gene index/indices outside ', ...
                     'model.genes: %s.'], ...
                    r,strjoin(string(invalidGeneIndices),', '));
            end

            for q = 1:numel(geneIndices)
                geneIdx = geneIndices(q);
                expression = strrep(expression, ...
                    sprintf('x(%d)',geneIdx),geneIds{geneIdx});
            end
        end

        expression = strrep(expression,'&',' and ');
        expression = strrep(expression,'|',' or ');
        expression = regexprep(expression,'\s+',' ');
        normalizedGrRules{r} = strtrim(expression);
    end
    gprSource = "model.rules + model.genes";

%% Case 3: retain gene membership from rxnGeneMat
elseif isfield(model,'rxnGeneMat') && ~isempty(model.rxnGeneMat) && ...
        isfield(model,'genes') && ~isempty(model.genes)

    geneIds = localTextColumn(model.genes,'model.genes');
    rxnGeneMat = model.rxnGeneMat;

    if size(rxnGeneMat,1) ~= nRxns && size(rxnGeneMat,2) == nRxns
        rxnGeneMat = rxnGeneMat.';
        warning('GPRCompatibility:TransposedRxnGeneMat', ...
            'model.rxnGeneMat was transposed to reactions-by-genes.');
    end

    if size(rxnGeneMat,1) ~= nRxns || ...
            size(rxnGeneMat,2) ~= numel(geneIds)
        error('GPRCompatibility:RxnGeneMatSizeMismatch', ...
            ['model.rxnGeneMat is %d-by-%d, but the model has %d ', ...
             'reactions and %d genes.'], ...
            size(rxnGeneMat,1),size(rxnGeneMat,2), ...
            nRxns,numel(geneIds));
    end

    for r = 1:nRxns
        geneIndices = find(rxnGeneMat(r,:) ~= 0);
        if ~isempty(geneIndices)
            % rxnGeneMat records membership but not AND/OR structure. "or"
            % is used only to create readable text; hasGPR and gene-overlap
            % tests depend on membership and therefore remain valid.
            normalizedGrRules{r} = strjoin(geneIds(geneIndices),' or ');
        end
    end
    gprSource = "model.rxnGeneMat + model.genes (gene-set fallback)";
    warning('GPRCompatibility:GeneLogicUnavailable', ...
        ['GPRs were reconstructed from rxnGeneMat. Gene membership is ', ...
         'preserved, but original AND/OR enzyme-complex logic is unavailable.']);

%% Case 4: no gene associations are available
else
    warning('GPRCompatibility:NoGeneAssociations', ...
        ['The model contains no usable grRules, rules+genes, or ', ...
         'rxnGeneMat+genes data. model.grRules was initialized as empty. ', ...
         'With REQUIRE_GPR=true, the motif search will return no results.']);
end

model.grRules = normalizedGrRules;

%% Build the logical enzyme-reaction mask used by Blocks 1-4
hasGPR = false(nRxns,1);
for r = 1:nRxns
    rule = strtrim(string(model.grRules{r}));
    hasGPR(r) = strlength(rule) > 0 && ...
        ~strcmpi(rule,"s0001") && ...
        ~strcmpi(rule,"spontaneous");
end

fprintf('GPR source: %s\n',gprSource);
fprintf('Enzyme-associated reactions: %d out of %d\n', ...
    sum(hasGPR),nRxns);

if ~any(hasGPR)
    warning('GPRCompatibility:NoEnzymeReactions', ...
        ['No enzyme-associated reactions were identified. Keep ', ...
         'REQUIRE_GPR=true for a strict search, or explicitly set it to ', ...
         'false for a topology-only exploratory search.']);
end
end

function values = localTextColumn(rawValues,fieldLabel)
% Convert char matrices, string arrays or cell arrays into column cellstr.
if ischar(rawValues)
    values = cellstr(rawValues);
elseif isstring(rawValues)
    values = cellstr(rawValues(:));
elseif iscell(rawValues)
    values = cell(numel(rawValues),1);
    for k = 1:numel(rawValues)
        value = string(rawValues{k});
        if ~isscalar(value)
            error('GPRCompatibility:NonScalarTextEntry', ...
                '%s entry %d is not scalar text.',fieldLabel,k);
        end
        values{k} = char(value);
    end
else
    error('GPRCompatibility:UnsupportedTextType', ...
        '%s must be a char array, string array or cell array.',fieldLabel);
end
values = values(:);
end
