function [D, meta] = prepareCircuitNetwork(model, cfg, p)
%PREPARECIRCUITNETWORK Build oriented incidence and annotation conflict tables.
required = {'S','rxns','mets','lb','ub'};
if ~isstruct(model) || ~isscalar(model) || ~all(isfield(model, required))
    error('Circuit:ModelFields', 'Model must have S, rxns, mets, lb, and ub.');
end
[nm,nr] = size(model.S);
if numel(model.mets) ~= nm || numel(model.rxns) ~= nr || ...
        numel(model.lb) ~= nr || numel(model.ub) ~= nr
    error('Circuit:ModelDimensions', 'Model field sizes do not agree with S.');
end
if any(~isfinite(nonzeros(model.S))) || any(isnan(model.lb(:))) || ...
        any(isnan(model.ub(:))) || any(model.lb(:) > model.ub(:))
    error('Circuit:ModelNumbers', 'S or reaction bounds contain invalid values.');
end
meta.metIDs = textColumn(model.mets,nm);
meta.rxnIDs = textColumn(model.rxns,nr);
if numel(unique(meta.metIDs)) ~= nm || numel(unique(meta.rxnIDs)) ~= nr
    error('Circuit:DuplicateIDs', 'Model metabolite and reaction IDs must be unique.');
end
meta.metNames = getText(model, {'metNames'}, nm, meta.metIDs);
meta.rxnNames = getText(model, {'rxnNames'}, nr, meta.rxnIDs);
meta.formulas = getText(model, {'metFormulas'}, nm, repmat({''},nm,1));
meta.gpr = getText(model, {'grRules'}, nr, repmat({''},nr,1));
meta.compartments = metaboliteCompartments(model, meta.metIDs);
meta.metKeys = meta.metIDs;
if strcmpi(cfg.METABOLITE_IDENTITY,'chemical')
    for k = 1:nm
        id = meta.metIDs{k}; comp = meta.compartments{k};
        base = regexprep(id, '\[[^\]]+\]$', '');
        if strcmp(base,id) && ~isempty(comp)
            for sep = {'_','.'}
                suffix = [sep{1},comp];
                if length(id) > length(suffix) && ...
                        strcmp(id(end-length(suffix)+1:end),suffix)
                    base = id(1:end-length(suffix)); break;
                end
            end
        end
        meta.metKeys{k} = base;
    end
end
[meta.keyNames,~,D.metGroup] = unique(meta.metKeys);
ng = numel(meta.keyNames);

ecFields = {'rxnECNumbers','eccodes','ecCodes','ecNumbers','rxnEC', ...
    'rxnECNumber','rxnECcodes','rxnEnzymeClassification','circuitSBMLEC'};
meta.ecTokens = cell(nr,1);
meta.geneTokens = cell(nr,1);
meta.ec = repmat({''},nr,1);
for j = 1:nr
    meta.ecTokens{j} = {};
end
for f = 1:numel(ecFields)
    if ~isfield(model,ecFields{f}), continue; end
    values = textColumn(model.(ecFields{f}),nr);
    for j = 1:nr
        ec = regexp(values{j}, '(?<![\d.])[1-7]\.\d+\.\d+\.\d+(?![\d.])', 'match');
        meta.ecTokens{j} = unique([meta.ecTokens{j},ec]);
    end
end
hasGeneMatrix = isfield(model,'rxnGeneMat') && isfield(model,'genes') && ...
    isequal(size(model.rxnGeneMat),[nr,numel(model.genes)]);
if hasGeneMatrix, geneIDs = textColumn(model.genes,numel(model.genes)); end
for j = 1:nr
    meta.ec{j} = strjoin(meta.ecTokens{j},'; ');
    genes = regexp(meta.gpr{j}, '[A-Za-z0-9_.:\-]+', 'match');
    genes = genes(~ismember(lowerCell(genes), ...
        {'and','or','not','true','false','spontaneous','unknown','na','nan'}));
    if hasGeneMatrix
        genes = [genes, geneIDs(find(model.rxnGeneMat(j,:)))']; %#ok<FNDSB>
    end
    % If only indexed COBRA rules are available, resolve x(i) via genes.
    if isempty(genes) && isfield(model,'rules') && isfield(model,'genes')
        rules = textColumn(model.rules,nr);
        ids = regexp(rules{j}, 'x\((\d+)\)', 'tokens');
        allGenes = textColumn(model.genes,numel(model.genes));
        for k = 1:numel(ids)
            index = str2double(ids{k}{1});
            if index >= 1 && index <= numel(allGenes)
                genes{end+1} = allGenes{index}; %#ok<AGROW>
            end
        end
    end
    meta.geneTokens{j} = unique(genes);
end
meta.enzymeTokens = cell(nr,1);
hasEC = ~cellfun(@isempty,meta.ecTokens);
hasGenes = ~cellfun(@isempty,meta.geneTokens);
curated = ~isempty(cfg.ENZYME_ID_FIELD);
customValues = repmat({''},nr,1);
if curated
    if ~isfield(model,cfg.ENZYME_ID_FIELD)
        error('Circuit:EnzymeField','Missing enzyme ID field: %s',cfg.ENZYME_ID_FIELD);
    end
    customValues = textColumn(model.(cfg.ENZYME_ID_FIELD),nr);
end
for j = 1:nr
    if curated
        enzymes = regexp(customValues{j}, '[^;,|\s]+', 'match');
        resources = cellfun(@(x) ['enzyme:',x],enzymes,'UniformOutput',false);
    else
        resources = cellfun(@(x) ['EC:',x],meta.ecTokens{j},'UniformOutput',false);
    end
    geneResources = cellfun(@(x) ['gene:',x],meta.geneTokens{j},'UniformOutput',false);
    meta.enzymeTokens{j} = unique([resources,geneResources]);
end
if curated
    % All selected reactions must have comparable, curated identities.
    evidence = ~cellfun(@isempty,strtrimCell(customValues));
elseif strcmpi(cfg.ENZYME_EVIDENCE,'EC')
    evidence = hasEC;
else
    evidence = hasEC | hasGenes;
end
notes = getText(model, {'rxnNotes'}, nr, repmat({''},nr,1));
spontaneous = false(nr,1);
for j = 1:nr
    spontaneous(j) = ~isempty(regexpi([meta.rxnNames{j},' ',notes{j},' ',meta.gpr{j}], ...
        'spontaneous|non[- ]?enzymatic|uncataly[sz]ed', 'once'));
end
S = sparse(model.S);
[row,col,value] = find(S);
keep = abs(value) > cfg.TOLERANCE;
S = sparse(row(keep),col(keep),value(keep),nm,nr);
twoSided = full(sum(S < 0,1)' > 0 & sum(S > 0,1)' > 0);
compOK = true(nr,1);
if ~isempty(cfg.COMPARTMENT)
    outside = ~strcmp(meta.compartments,cfg.COMPARTMENT);
    compOK = full(sum(spones(S(outside,:)),1)' == 0);
end
eligible = evidence & ~spontaneous & twoSided & compOK;
forward = find(eligible & model.ub(:) > cfg.TOLERANCE);
reverse = [];
if ~cfg.FORWARD_ONLY
    reverse = find(eligible & model.lb(:) < -cfg.TOLERANCE);
end
D.origRxn = [forward;reverse];
D.direction = [ones(numel(forward),1);-ones(numel(reverse),1)];
no = numel(D.origRxn);
OS = S(:,D.origRxn) * spdiags(D.direction,0,no,no);
D.incidence = spones(OS);
% GroupCount counts distinct model species of the same chemical per reaction.
groupMatrix = sparse(D.metGroup,(1:nm)',ones(nm,1),ng,nm);
D.groupCount = groupMatrix * D.incidence;
[ri,rj,sv] = find(OS);
if isempty(ri), validRole = false(0,1);
else, validRole = full(D.groupCount(sub2ind([ng,no],D.metGroup(ri),rj))) == 1;
end
D.sub = sparse(ri(sv < 0 & validRole),rj(sv < 0 & validRole), ...
    ones(nnz(sv < 0 & validRole),1),nm,no);
D.prod = sparse(ri(sv > 0 & validRole),rj(sv > 0 & validRole), ...
    ones(nnz(sv > 0 & validRole),1),nm,no);
D.conflict = resourceConflicts(meta.enzymeTokens,nr);
D.nm = nm; D.nr = nr; D.no = no; D.ng = ng; D.p = p;
D.exemptAux = matchingKeys(cfg.AUXILIARY_REUSE_ALLOWED,meta,D.metGroup);
D.allowedMet = true(nm,1);
if ~isempty(cfg.COMPARTMENT)
    D.allowedMet = D.allowedMet & strcmp(meta.compartments,cfg.COMPARTMENT);
end
excluded = matchingKeys(cfg.EXCLUDED_ROLE_METABOLITES,meta,D.metGroup);
D.allowedMet = D.allowedMet & ~ismember(D.metGroup,excluded);
D.activation = sparse(nm,no);
if strcmpi(cfg.NEURON_INPUT_MODE,'activation')
    pairs = cfg.ACTIVATION_PAIRS;
    if ~iscell(pairs) || isempty(pairs) || size(pairs,2) < 2
        error('Circuit:ActivationEvidence', ...
            'Activation mode requires ACTIVATION_PAIRS = {reactionID, metaboliteID; ...}.');
    end
    for k = 1:size(pairs,1)
        j = find(strcmp(meta.rxnIDs,char(pairs{k,1})),1);
        i = find(strcmp(meta.metIDs,char(pairs{k,2})),1);
        if isempty(j) || isempty(i)
            error('Circuit:ActivationEvidence','Unknown reaction/metabolite in activation pair %d.',k);
        end
        D.activation(i,D.origRxn == j) = 1;
    end
end
meta.originalReactions = nr;
meta.originalMetabolites = nm;
meta.annotatedECReactions = nnz(hasEC);
meta.annotatedGPRReactions = nnz(hasGenes);
meta.eligibleReactions = nnz(eligible);
meta.orientedReactions = no;
meta.noEvidence = nnz(~evidence);
meta.spontaneousReactions = nnz(spontaneous);
meta.boundaryReactions = nnz(~twoSided);
meta.curatedEnzymes = curated;
end

function C = resourceConflicts(tokens,n)
rr = []; names = {};
for j = 1:n
    t = tokens{j};
    rr = [rr;repmat(j,numel(t),1)]; %#ok<AGROW>
    names = [names;t(:)]; %#ok<AGROW>
end
if isempty(names), C = speye(n); return; end
[~,~,cc] = unique(names);
B = sparse(rr,cc,ones(numel(rr),1),n,max(cc));
C = spones(B*B' + speye(n)); % Same underlying reaction always conflicts.
end

function groups = matchingKeys(values,meta,metGroup)
if ischar(values) || isstring(values), values = cellstr(values); end
groups = [];
for k = 1:numel(values)
    mask = strcmp(meta.metIDs,char(values{k})) | strcmp(meta.metKeys,char(values{k}));
    if ~any(mask)
        error('Circuit:UnknownMetabolite', 'Unknown metabolite ID/base ID: %s',char(values{k}));
    end
    groups = [groups;metGroup(mask)]; %#ok<AGROW>
end
groups = unique(groups);
end

function comps = metaboliteCompartments(model,ids)
n = numel(ids); comps = repmat({''},n,1);
known = {};
if isfield(model,'comps'), known = textColumn(model.comps,numel(model.comps)); end
if isfield(model,'metComps') && numel(model.metComps) == n && ~isempty(known)
    for i = 1:n
        c = model.metComps(i);
        if isnumeric(c) && isfinite(c) && c == floor(c) && c >= 1 && c <= numel(known)
            comps{i} = known{c};
        end
    end
end
for i = 1:n
    if ~isempty(comps{i}), continue; end
    tok = regexp(ids{i}, '\[([^\]]+)\]$', 'tokens','once');
    if ~isempty(tok), comps{i} = tok{1}; continue; end
    for k = 1:numel(known)
        suffix = ['_',known{k}];
        if length(ids{i}) > length(suffix) && strcmp(ids{i}(end-length(suffix)+1:end),suffix)
            comps{i} = known{k}; break;
        end
    end
end
end

function out = getText(s,fields,n,fallback)
out = fallback;
for k = 1:numel(fields)
    if isfield(s,fields{k}) && ~isempty(s.(fields{k}))
        out = textColumn(s.(fields{k}),n); return;
    end
end
end

function out = textColumn(values,n)
if isempty(values), out = repmat({''},n,1); return; end
if ischar(values), values = cellstr(values); end
if isstring(values), values = cellstr(values); end
if isnumeric(values), values = num2cell(values); end
if ~iscell(values) || numel(values) ~= n
    error('Circuit:AnnotationSize','An annotation field is not a %d-element text array.',n);
end
out = cell(n,1);
for k = 1:n
    v = values{k};
    if isempty(v), out{k} = '';
    elseif iscell(v), out{k} = strjoin(textColumn(v,numel(v))','; ');
    elseif isnumeric(v), out{k} = strtrim(num2str(v));
    else, out{k} = strtrim(char(v));
    end
end
end
function out = lowerCell(x)
out = cellfun(@lower,x,'UniformOutput',false);
end
function out = strtrimCell(x)
out = cellfun(@strtrim,x,'UniformOutput',false);
end
