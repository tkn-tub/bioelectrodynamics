function [ranking, paths, rankInfo] = rankCircuitInputsByCAMP(model,motifs,cfg,p)
%RANKCIRCUITINPUTSBYCAMP Rank returned motifs by directed conversion distance.
% Distance counts annotated, bound-permitted substrate -> product reactions
% from intracellular 3',5'-cAMP to the exact model species assigned to m_in.
% This is a graph projection, not evidence of regulation, flux or kinetics.
% Other substrates need not be produced by this path. Selected signal
% chemicals and enzymes are protected against bridge reactions. Enzymes may
% repeat within a bridge; BridgeEnzymeConflict reports that limitation.
% No motifs are reordered or removed. Paths follow the table's sorted rows.
if nargin < 3 || isempty(cfg), cfg = struct(); end
cfg = completeConfiguration(cfg);
if nargin < 4 || isempty(p), p = circuitPattern(cfg.NEURON_INPUT_MODE); end
inputRole = find(strcmp(p.metNames,'m_in'),1);
if isempty(inputRole), error('Circuit:CAMPPattern','Pattern has no m_in role.'); end
[D,meta] = prepareCircuitNetwork(model,cfg,p);
nm = D.nm; no = D.no; n = numel(motifs);

% Base IDs always identify chemical aliases across compartments, even when
% the motif search was explicitly configured to use model-ID identity.
baseIDs = chemicalBaseIDs(meta);
[~,~,chemicalGroup] = unique(baseIDs);
sourceSynonyms = {'camp','35camp','cyclicamp','35cyclicamp', ...
    'adenosine35cyclicmonophosphate','adenosine35cyclicphosphate', ...
    'adenosinecyclic35monophosphate','35cyclicadenosinemonophosphate', ...
    'cyclicadenosinemonophosphate'};
idNorm = cellfun(@normalizeText,baseIDs,'UniformOutput',false);
nameNorm = cellfun(@normalizeText,meta.metNames,'UniformOutput',false);
autoSource = ismember(idNorm,sourceSynonyms) | ismember(nameNorm,sourceSynonyms);
% A conventional ID must never override an explicitly 2',3'-cyclic name.
is23 = cellfun(@isTwoThreeCAMP,baseIDs) | cellfun(@isTwoThreeCAMP,meta.metNames);
autoSource(is23) = false;
if isempty(cfg.CAMP_SOURCE_METABOLITES)
    sourceCandidates = find(autoSource); detection = 'automatic_exact_synonyms';
else
    sourceCandidates = resolveIDs(cfg.CAMP_SOURCE_METABOLITES,meta,baseIDs,false);
    if any(is23(sourceCandidates))
        error('Circuit:CAMPSourceIsomer','An explicit cAMP source is annotated as 2'',3''-cAMP.');
    end
    detection = 'explicit_model_ids';
end
sources = sourceCandidates;
if ~isempty(cfg.CAMP_SOURCE_COMPARTMENT)
    sources = sources(strcmp(meta.compartments(sources),cfg.CAMP_SOURCE_COMPARTMENT));
end
sources = sort(unique(sources));

currency = currencySynonyms();
excluded = ismember(idNorm,currency) | ismember(nameNorm,currency);
explicitExcluded = resolveIDs(cfg.CAMP_EXCLUDED_PATH_METABOLITES,meta,baseIDs,true);
excluded(explicitExcluded) = true;
% A chemical exclusion protects every compartment copy of that chemical.
excluded = ismember(chemicalGroup,unique(chemicalGroup(excluded)));

% Build the bridge graph from actual stoichiometry. D.sub/D.prod deliberately
% suppress chemical aliases for signal matching, which would lose transports.
[ri,rj,sv] = find(sparse(model.S));
keep = abs(sv) > cfg.TOLERANCE;
S = sparse(ri(keep),rj(keep),sv(keep),D.nm,D.nr);
OS = S(:,D.origRxn) * spdiags(D.direction,0,no,no);
sub = OS < 0; prod = OS > 0; incidence = spones(OS);
substrateOf = cell(nm,1); productsOf = cell(no,1);
[si,sj] = find(sub);
for k = 1:numel(si), substrateOf{si(k)}(end+1) = sj(k); end %#ok<AGROW>
[pi,pj] = find(prod);
for k = 1:numel(pi), productsOf{pj(k)}(end+1) = pi(k); end %#ok<AGROW>
transport = false(no,1); unknownComp = false(no,1);
for j = 1:no
    participating = find(incidence(:,j));
    comps = meta.compartments(participating);
    unknownComp(j) = any(cellfun(@isempty,comps));
    transport(j) = numel(unique(comps(~cellfun(@isempty,comps)))) > 1;
end
baseAllowed = true(no,1);
if ~cfg.CAMP_ALLOW_TRANSPORT, baseAllowed(transport) = false; end

rankInfo.status = 'complete';
rankInfo.sourceDetection = detection;
rankInfo.sourceCandidateMetaboliteIDs = meta.metIDs(sourceCandidates);
rankInfo.sourceIndices = sources;
rankInfo.sourceMetaboliteIndices = sources;
rankInfo.sourceMetaboliteIDs = meta.metIDs(sources);
rankInfo.metadata = meta;
rankInfo.sourceCompartment = cfg.CAMP_SOURCE_COMPARTMENT;
rankInfo.excludedPathMetaboliteIDs = meta.metIDs(excluded);
rankInfo.maxDepth = cfg.CAMP_MAX_DEPTH;
rankInfo.allowTransport = cfg.CAMP_ALLOW_TRANSPORT;
rankInfo.unknownCompartmentOrientations = nnz(unknownComp & baseAllowed);
rankInfo.motifsRanked = n;
rankInfo.scope = 'Only returned motifs; no global optimum is established.';
rankInfo.interpretation = ['Directed enzymatic conversion distance; other substrates ', ...
    'are assumed available. No receptor, regulatory, flux or kinetic inference.'];
rankInfo.bridgeEnzymePolicy = ['Bridge enzymes cannot conflict with selected circuit ', ...
    'enzymes. Conflicts among bridge steps are reported, not rejected.'];
rankInfo.depthInterpretation = ['A bounded failed search does not establish ', ...
    'disconnection outside the configured depth.'];

Rank = zeros(n,1); MotifIndex = (1:n)'; Distance = Inf(n,1);
InputMetabolite = cell(n,1); InputName = cell(n,1); Status = cell(n,1);
SourceMetabolite = repmat({''},n,1); PathReactionIDs = repmat({''},n,1);
PathDirections = repmat({''},n,1); PathMetabolites = repmat({''},n,1);
BridgeEnzymeConflict = false(n,1);
template = struct('metIndex',[],'reactionIndex',[],'direction',[], ...
    'distance',Inf,'status','','enzymeConflict',false,'unknownCompartments',false);
paths = repmat(template,n,1);
if n == 0, rankInfo.status = 'no_motifs';
elseif isempty(sources), rankInfo.status = 'source_missing'; end

for k = 1:n
    motif = motifs(k);
    validateMotif(motif,nm,D.nr,numel(p.metNames));
    target = motif.metIndex(inputRole);
    InputMetabolite{k} = meta.metIDs{target}; InputName{k} = meta.metNames{target};
    if ismember(target,sources)
        paths(k).metIndex = target;
        paths(k).distance = 0; paths(k).status = 'exact_input';
    elseif isempty(sources)
        paths(k).status = 'source_missing';
    elseif excluded(target)
        paths(k).status = 'excluded_input';
    else
        % Any species of another selected signal chemical must be absent
        % from every bridge reaction, including an auxiliary participation.
        signalGroups = unique(chemicalGroup(motif.metIndex));
        inputGroup = chemicalGroup(target);
        otherSignal = ismember(chemicalGroup,setdiff(signalGroups,inputGroup));
        forbiddenSignal = full(sum(incidence(otherSignal,:),1)' > 0);
        selectedRxns = unique(motif.reactionIndex(:));
        conflicting = full(sum(D.conflict(:,selectedRxns),2) > 0);
        allowed = baseAllowed & ~forbiddenSignal & ~conflicting(D.origRxn);
        targetFamily = chemicalGroup == inputGroup;
        touchesInput = full(sum(incidence(targetFamily,:),1)' > 0);
        % Input species cannot be used as an intermediate. A reaction that
        % touches this signal is usable only for arrival at the exact target.
        allowed = allowed & (~touchesInput | full(prod(target,:))');
        [metPath,oriPath] = shortestPath(sources,target,excluded,allowed, ...
            touchesInput,substrateOf,productsOf,cfg.CAMP_MAX_DEPTH,nm,no);
        if isempty(metPath)
            if isfinite(cfg.CAMP_MAX_DEPTH)
                paths(k).status = 'not_reached_within_depth';
            else
                paths(k).status = 'unreachable';
            end
        else
            paths(k).metIndex = metPath;
            paths(k).reactionIndex = D.origRxn(oriPath);
            paths(k).direction = D.direction(oriPath);
            paths(k).distance = numel(oriPath); paths(k).status = 'reachable';
            pathRxns = paths(k).reactionIndex;
            paths(k).enzymeConflict = nnz(triu(D.conflict(pathRxns,pathRxns),1)) > 0;
            paths(k).unknownCompartments = any(unknownComp(oriPath));
        end
    end
    Distance(k) = paths(k).distance; Status{k} = paths(k).status;
    BridgeEnzymeConflict(k) = paths(k).enzymeConflict;
    if ~isempty(paths(k).metIndex)
        SourceMetabolite{k} = meta.metIDs{paths(k).metIndex(1)};
        PathMetabolites{k} = strjoin(meta.metIDs(paths(k).metIndex),'; ');
        PathReactionIDs{k} = strjoin(meta.rxnIDs(paths(k).reactionIndex),'; ');
        dirs = arrayfun(@(d) sprintf('%+d',d),paths(k).direction,'UniformOutput',false);
        PathDirections{k} = strjoin(dirs,'; ');
    end
end

% Stable deterministic tie handling, while retaining the original motif ID.
[~,order] = sortrows([Distance,MotifIndex],[1 2]);
if n > 0
    sortedDistance = Distance(order); lastRank = 1;
    for k = 1:n
        if k == 1 || sortedDistance(k) ~= sortedDistance(k-1), lastRank = k; end
        Rank(order(k)) = lastRank;
    end
end
ranking = table(Rank,MotifIndex,InputMetabolite,InputName,Distance,Status, ...
    SourceMetabolite,PathReactionIDs,PathDirections,PathMetabolites,BridgeEnzymeConflict);
ranking = ranking(order,:); paths = paths(order);
end

function [metPath,oriPath] = shortestPath(sources,target,excluded,allowed, ...
    touchesInput,substrateOf,productsOf,maxDepth,nm,no)
% Multi-source breadth-first search. Each orientation expands once at its
% earliest accessible substrate. This is a substrate/product projection.
distance = Inf(nm,1); previousMet = zeros(nm,1); previousOri = zeros(nm,1);
usedOrientation = false(no,1); queue = zeros(nm,1);
queue(1:numel(sources)) = sources; head = 1; tail = numel(sources);
distance(sources) = 0; metPath = []; oriPath = [];
while head <= tail && ~isfinite(distance(target))
    current = queue(head); head = head+1;
    if distance(current) >= maxDepth, continue; end
    % Explicitly excluded compounds remain valid exact sources, but cannot
    % serve as currency intermediates. A source is allowed to start its path.
    if excluded(current) && ~ismember(current,sources), continue; end
    candidates = substrateOf{current};
    for j = candidates
        if ~allowed(j) || usedOrientation(j), continue; end
        usedOrientation(j) = true;
        products = productsOf{j};
        if touchesInput(j), products = target; end
        for next = products
            if excluded(next) || isfinite(distance(next)), continue; end
            distance(next) = distance(current)+1;
            previousMet(next) = current; previousOri(next) = j;
            tail = tail+1; queue(tail) = next;
        end
    end
end
if ~isfinite(distance(target)), return; end
metPath = target; node = target;
while previousMet(node) ~= 0
    oriPath = [previousOri(node);oriPath]; %#ok<AGROW>
    node = previousMet(node); metPath = [node;metPath]; %#ok<AGROW>
end
end

function cfg = completeConfiguration(cfg)
defaults = circuitSearchDefaults();
defaults.CAMP_SOURCE_METABOLITES = {};
defaults.CAMP_SOURCE_COMPARTMENT = 'c';
defaults.CAMP_ALLOW_TRANSPORT = false;
defaults.CAMP_EXCLUDED_PATH_METABOLITES = {};
defaults.CAMP_MAX_DEPTH = 6;
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(cfg,names{k}), cfg.(names{k}) = defaults.(names{k}); end
end
d = cfg.CAMP_MAX_DEPTH;
if ~isnumeric(d) || ~isscalar(d) || isnan(d) || d < 0 || (isfinite(d) && d ~= floor(d))
    error('Circuit:CAMPDepth','CAMP_MAX_DEPTH must be a nonnegative integer or Inf.');
end
flag = cfg.CAMP_ALLOW_TRANSPORT;
if ~(islogical(flag) || isnumeric(flag)) || ~isscalar(flag) || ~ismember(flag,[0 1])
    error('Circuit:CAMPTransport','CAMP_ALLOW_TRANSPORT must be true or false.');
end
comp = cfg.CAMP_SOURCE_COMPARTMENT;
if isstring(comp) && isscalar(comp), comp = char(comp); end
if ~ischar(comp) || size(comp,1) > 1
    error('Circuit:CAMPCompartment','CAMP_SOURCE_COMPARTMENT must be a compartment ID or ''''.');
end
cfg.CAMP_SOURCE_COMPARTMENT = comp;
cfg.CAMP_SOURCE_METABOLITES = textList(cfg.CAMP_SOURCE_METABOLITES);
cfg.CAMP_EXCLUDED_PATH_METABOLITES = textList(cfg.CAMP_EXCLUDED_PATH_METABOLITES);
end

function ids = chemicalBaseIDs(meta)
ids = meta.metIDs;
for k = 1:numel(ids)
    id = ids{k}; comp = meta.compartments{k};
    base = regexprep(id,'\[[^\]]+\]$','');
    if strcmp(base,id) && ~isempty(comp)
        suffixes = {['_',comp],['.',comp]};
        for j = 1:numel(suffixes)
            suffix = suffixes{j};
            if length(id) > length(suffix) && strcmp(id(end-length(suffix)+1:end),suffix)
                base = id(1:end-length(suffix)); break;
            end
        end
    end
    ids{k} = base;
end
end

function out = normalizeText(value)
value = lower(strtrim(char(value)));
value = regexprep(value,'^m_','');
value = regexprep(value,'\[[^\]]+\]$','');
value = regexprep(value,'\(\d*[+-]\)$','');
out = regexprep(value,'[^a-z0-9]','');
end

function yes = isTwoThreeCAMP(value)
v = normalizeText(value);
yes = ~isempty(strfind(v,'23')) && ...
    (~isempty(strfind(v,'camp')) || ~isempty(strfind(v,'cyclic'))); %#ok<STREMP>
end

function ids = resolveIDs(values,meta,baseIDs,allowBase)
ids = [];
for k = 1:numel(values)
    matches = strcmp(meta.metIDs,values{k});
    if allowBase, matches = matches | strcmp(baseIDs,values{k}); end
    if ~any(matches)
        error('Circuit:CAMPUnknownMetabolite','Unknown cAMP ranking metabolite ID: %s',values{k});
    end
    ids = [ids;find(matches)]; %#ok<AGROW>
end
ids = unique(ids);
end

function values = textList(values)
if isempty(values), values = {}; return; end
if ischar(values) || isstring(values), values = cellstr(values); end
if ~iscell(values) || any(~cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)),values))
    error('Circuit:CAMPList','cAMP metabolite selections must be text IDs in a cell array.');
end
values = cellfun(@char,values(:),'UniformOutput',false);
end

function validateMotif(motif,nm,nr,nRoles)
if ~isstruct(motif) || ~all(isfield(motif,{'metIndex','reactionIndex'})) || ...
        numel(motif.metIndex) ~= nRoles || isempty(motif.reactionIndex)
    error('Circuit:CAMPMotif','Each motif must contain its selected role species and reactions.');
end
mi = motif.metIndex(:); ri = motif.reactionIndex(:);
if any(~isfinite(mi) | mi < 1 | mi > nm | mi ~= floor(mi)) || ...
        any(~isfinite(ri) | ri < 1 | ri > nr | ri ~= floor(ri))
    error('Circuit:CAMPMotif','A motif contains invalid model indices.');
end
end

function synonyms = currencySynonyms()
% Exact normalized chemical names/IDs only. AMP and adenosine are retained
% deliberately so cAMP phosphodiesterase and downstream routes remain usable.
synonyms = {'h2o','water','h','h1','proton','protons','hydrogenion', ...
    'atp','adenosinetriphosphate','adenosine5triphosphate', ...
    'adp','adenosinediphosphate','adenosine5diphosphate', ...
    'pi','phosphate','orthophosphate','inorganicphosphate', ...
    'ppi','pyrophosphate','diphosphate','inorganicpyrophosphate', ...
    'gtp','guanosinetriphosphate','gdp','guanosinediphosphate', ...
    'utp','uridinetriphosphate','udp','uridinediphosphate', ...
    'ctp','cytidinetriphosphate','cdp','cytidinediphosphate', ...
    'itp','inosinetriphosphate','idp','inosinediphosphate', ...
    'nad','nadplus','nicotinamideadeninedinucleotide', ...
    'nadh','reducednicotinamideadeninedinucleotide', ...
    'nadp','nadpplus','nicotinamideadeninedinucleotidephosphate', ...
    'nadph','reducednicotinamideadeninedinucleotidephosphate', ...
    'fad','flavinadeninedinucleotide','fadh2', ...
    'fmn','flavinmononucleotide','fmnh2','coa','coenzymea', ...
    'o2','oxygen','dioxygen','co2','carbondioxide','hco3','bicarbonate', ...
    'na1','sodium','sodiumion','k','potassium','potassiumion', ...
    'cl','chloride','chlorideion'};
end
