function [isMembrane,labels] = classifyMinimalMembrane(model)
% Classify each model reaction once. Internal organelle membranes are retained.
nm = size(model.S,1); nr = size(model.S,2);
labels = strings(nm,1);
if isfield(model,'metComps') && numel(model.metComps) == nm && ...
        (isfield(model,'compNames') || isfield(model,'comps'))
    if isfield(model,'compNames'), names = model.compNames; else, names = model.comps; end
    if ischar(names), names = cellstr(names); end
    names = string(names(:)); idx = double(model.metComps(:));
    if any(~isfinite(idx) | idx < 1 | idx > numel(names) | idx ~= floor(idx))
        error('MinimalMembrane:Compartments','Invalid model.metComps.');
    end
    labels = lower(strtrim(names(idx)));
else
    ids = string(model.mets(:));
    for m = 1:nm
        token = regexp(char(ids(m)),'\[([^\]]+)\]$','tokens','once');
        if isempty(token), token = regexp(char(ids(m)),'_([A-Za-z][A-Za-z0-9]*)$','tokens','once'); end
        if isempty(token)
            error('MinimalMembrane:UnknownCompartment', ...
                'Cannot determine the compartment of %s. Supply metComps and compNames.',ids(m));
        end
        code = string(token{1}); label = lower(code);
        if isfield(model,'comps') && isfield(model,'compNames')
            codes = model.comps; names = model.compNames;
            if ischar(codes), codes = cellstr(codes); end
            if ischar(names), names = cellstr(names); end
            k = find(strcmpi(string(codes(:)),code),1);
            if ~isempty(k), label = lower(string(names(k))); end
        end
        labels(m) = strtrim(label);
    end
end
cyto = ismember(labels,["c","cytoplasm","cytosol","cytoplasmic"]);
extra = ismember(labels,["e","extracellular","extracellular space","extracellular region"]);
envelope = ismember(labels,["ce","pm","cell envelope","plasma membrane","cell membrane"]);
if ~any(cyto) || ~any(extra)
    error('MinimalMembrane:Labels', ...
        'Could not identify cytoplasm and extracellular space. Inspect compartment names before filtering.');
end
participates = model.S ~= 0;
isMembrane = full(any(participates(envelope,:),1) | ...
    (any(participates(cyto,:),1) & any(participates(extra,:),1))).';
if numel(isMembrane) ~= nr, error('MinimalMembrane:Size','Classification size mismatch.'); end
end
