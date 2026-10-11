function [hasEC,codes,fieldsUsed] = getReactionECAnnotations(model)
%GETREACTIONECANNOTATIONS Resolve complete EC codes without inferring them.
% At least one complete four-number EC code is required per reaction.
% Merge annotations from eccodes, rxnECNumbers and ecNumbers. Partial codes,
% missing values and placeholders fail. This check is independent of GPRs.
if ~isstruct(model) || ~isscalar(model) || ...
        ~all(isfield(model,{'S','rxns'}))
    error('EC:InvalidModel','A scalar model with S and rxns is required.');
end
nRxns = size(model.S,2);
if numel(model.rxns) ~= nRxns
    error('EC:ModelSize','The reaction IDs disagree with model.S.');
end
hasEC = false(nRxns,1);
codes = strings(nRxns,1);
fieldsUsed = strings(nRxns,1);
fields = {'eccodes','rxnECNumbers','ecNumbers'};
pattern = ['(?<![\w.\-])(?:[Ee][Cc]\s*:?\s*)?', ...
    '([1-7]\.[0-9]+\.[0-9]+\.[0-9]+)(?![\w.\-])'];
for f = 1:numel(fields)
    field = fields{f};
    if ~isfield(model,field) || isempty(model.(field)), continue; end
    values = model.(field);
    if ischar(values)
        if size(values,1) ~= nRxns
            error('EC:AnnotationSize','model.%s needs one entry per reaction.',field);
        end
        values = cellstr(values);
    elseif numel(values) ~= nRxns
        error('EC:AnnotationSize','model.%s needs one entry per reaction.',field);
    end
    for r = 1:nRxns
        if iscell(values), value = values{r}; else, value = values(r); end
        tokens = regexp(char(localText(value)),pattern,'tokens');
        if isempty(tokens), continue; end
        matches = cellfun(@(t) t{1},tokens,'UniformOutput',false);
        merged = unique([split(codes(r),"; ");string(matches(:))],'stable');
        merged = merged(strlength(merged) > 0);
        codes(r) = strjoin(merged,"; ");
        hasEC(r) = true;
        if strlength(fieldsUsed(r)) == 0
            fieldsUsed(r) = string(field);
        else
            fieldsUsed(r) = fieldsUsed(r) + "; " + string(field);
        end
    end
end
end

function out = localText(value)
out = "";
if iscell(value)
    parts = strings(numel(value),1);
    for k = 1:numel(value), parts(k) = localText(value{k}); end
    out = strjoin(parts,"; ");
elseif ischar(value) || isstring(value)
    parts = string(value);
    parts = parts(~ismissing(parts));
    out = strjoin(parts(:),"; ");
end
end
