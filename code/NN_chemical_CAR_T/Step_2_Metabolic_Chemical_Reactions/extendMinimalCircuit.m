function [motifs,diagnostics,paths,info,evidence] = extendMinimalCircuit(model,bad,paths,up,stage,o)
% Controlled reaction fragments plus explicit NONCONSUMING regulation links.
% Design mode assigns regulator links; evidence mode requires supplied links.
if ~ismember(stage,[3 4]), error('MinimalExtend:Stage','Use stage 3 or 4.'); end
if ~isfield(o,'mapCsLeak'), o.mapCsLeak = false; end
if ~isscalar(o.mapCsLeak) || ~ismember(o.mapCsLeak,[0 1])
    error('MinimalExtend:BufferMode','mapCsLeak must be true or false.');
end
if stage == 4 && ismember('r_CsLeak',up.Properties.VariableNames) ~= logical(o.mapCsLeak)
    error('MinimalExtend:BufferMode','Set MAP_COUNTER_CS_LEAK to the mode used by the loaded Block 3 checkpoint.');
end
for key = {'maxHits','poolCap','maxChecks','maxJoinChecks','optionsCap','nodeCap','nodesPerRole'}
    value = o.(key{1});
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || isnan(value) || ...
            value <= 0 || value ~= floor(value)
        error('MinimalExtend:Budget','%s must be a positive integer or inf.',key{1});
    end
end
if ~ismember(string(o.regulationMode),["design","evidence"])
    error('MinimalExtend:Mode','REGULATION_MODE must be design or evidence.');
end
evidence = loadMinimalRegulation(model,o.regulationFile);
graph = minimalDirectedGraph(model,bad,o.requireGPR);
fprintf('Block %d: %d eligible EC-annotated reactions; regulation mode %s.\n', ...
    stage,sum(graph.eligible),string(o.regulationMode));
paths = requireECAnnotatedAdapters(model,paths,getReactionECAnnotations(model),o.requireAdapterEC);
prior = minimalCircuitSchema(stage-1,o.mapCsLeak);
if ~all(ismember([prior.mets prior.rxns prior.dirs],up.Properties.VariableNames))
    error('MinimalExtend:OldSchema','Upstream motifs do not use the new minimal topology.');
end
indices = table2array(up(:,prior.mets));
if ~isnumeric(indices) || ~isreal(indices) || ...
        any(~isfinite(indices) | indices < 1 | indices > size(model.S,1) | indices ~= floor(indices),'all')
    error('MinimalExtend:ModelIndex','Upstream metabolite indices do not match this model.');
end
if stage == 3 || stage == 4
    % Physical aliases follow the updated diagram: m_d=Y1, m_d2=Y2.
    if ismember('m_d',up.Properties.VariableNames) && any(up.m_d ~= up.m_Y1)
        error('MinimalExtend:OutputAlias','m_d must identify Y1 in the current topology.');
    end
    if ismember('m_d1',up.Properties.VariableNames) && any(up.m_d1 ~= up.m_Y1)
        error('MinimalExtend:OutputAlias','m_d1 must identify Y1 in the current topology.');
    end
    if ismember('m_d2',up.Properties.VariableNames) && any(up.m_d2 ~= up.m_Y2)
        error('MinimalExtend:OutputAlias','m_d2 must identify the physical Block 2 output Y2.');
    end
    up.m_d2 = up.m_Y2;
    up.m_d = up.m_Y1;
    up.m_d1 = up.m_Y1;
end
original = (1:height(up)).';
if ~isempty(up)
    [up,index] = requireECAnnotatedRows(up,graph.eligible,prior.rxns,paths);
    original = original(index);
    keep = true(height(up),1);
    feasible = [graph.rxn,graph.dir];
    for k = 1:numel(prior.rxns)
        keep = keep & ismember([up.(prior.rxns{k}),up.(prior.dirs{k})],feasible,'rows');
    end
    up = up(keep,:); original = original(keep);
end
% Exclude signals/reactions present in EVERY upstream row before fragment
% caps. Otherwise early fragments can all reuse a fixed input or counter C.
graph.excludedMet = false(size(model.S,1),1);
graph.excludedEdge = false(numel(graph.rxn),1);
if ~isempty(up)
    universalM = up{1,prior.mets}; universalR = up{1,prior.rxns};
    for u = 2:height(up)
        universalM = intersect(universalM,up{u,prior.mets});
        universalR = intersect(universalR,up{u,prior.rxns});
        if isempty(universalM) && isempty(universalR), break; end
    end
    graph.excludedMet(universalM) = true;
    graph.excludedEdge = ismember(graph.rxn,universalR);
end
oldNames = up.Properties.VariableNames;
oldNames = oldNames(~ismember(oldNames,{'RegulatoryEvidenceComplete','KineticsValidated'}));
if stage == 3
    counterRoles = {'r_Cp','r_Cg','r_Cr'};
    if o.mapCsLeak, counterRoles{end+1} = 'r_CsLeak'; end
    newNames = [{'block2_idx','counter_idx','m_c','m_Cs','m_PC'},counterRoles,strrep(counterRoles,'r_','dir_')];
else
    newNames = {'block3_idx','nn_idx','m_Ni','m_Na','m_Nb','m_out', ...
        'm_PNi','m_PNa','m_PNb','m_PNo', ...
        'r_Ni','r_NiDecay','r_Na','r_NaDecay','r_Nb','r_NbDecay','r_No','r_NoDecay', ...
        'dir_Ni','dir_NiDecay','dir_Na','dir_NaDecay','dir_Nb','dir_NbDecay','dir_No','dir_NoDecay'};
end
names = [oldNames,newNames];
rows = zeros(min(1000,o.maxHits),numel(names));
supported = false(size(rows,1),1);
upValues = table2array(up(:,oldNames));
nHits = 0; tested = 0; rejected = 0; poolStats = zeros(1,4);
info = struct('stage',stage,'regulationMode',string(o.regulationMode), ...
    'counterCsLeakMapped',logical(o.mapCsLeak), ...
    'upstreamRows',height(up),'acceptedRows',0,'joinChecks',0, ...
    'joinLimitReached',false,'hitLimitReached',false,'poolTruncated',false, ...
    'regulationIsInferredFromS',false,'kineticsValidated',false);
clock = tic; pulse = tic;
if isempty(up)
    pool = zeros(0,1);
elseif stage == 3
    [pool,poolStats] = localCounterPool(graph,o,up,evidence);
    joinPool(pool,1:height(up));
elseif string(o.regulationMode) == "design"
    [pool,poolStats] = localNNPool(graph,o,0,evidence);
    joinPool(pool,1:height(up));
else
    % Search once for each distinct counter output, which fixes C->r_Ni.
    inputs = unique(up.m_c,'stable');
    pool = zeros(0,16);
    for a = 1:numel(inputs)
        [fragment,stat] = localNNPool(graph,o,inputs(a),evidence);
        poolStats = poolStats + stat;
        joinPool(fragment,find(up.m_c == inputs(a)));
        if nHits >= o.maxHits || tested >= o.maxJoinChecks, break; end
    end
end
motifs = array2table(rows(1:nHits,:),'VariableNames',names);
motifs.RegulatoryEvidenceComplete = supported(1:nHits);
motifs.KineticsValidated = false(nHits,1);
diagnostics = table(poolStats(1),poolStats(2),logical(poolStats(3)),logical(poolStats(4)), ...
    tested,rejected,nHits,'VariableNames',{'PoolChecks','PoolCandidates', ...
    'PoolLimitReached','PoolCheckLimitReached','JoinChecks','RejectedJoins','AcceptedJoins'});
info.acceptedRows = nHits; info.joinChecks = tested;
info.hitLimitReached = nHits >= o.maxHits; info.joinLimitReached = tested >= o.maxJoinChecks;
info.poolTruncated = any(poolStats(3:4) ~= 0); info.elapsedSeconds = toc(clock);
fprintf('Block %d: %d structural joins; %d have supplied regulatory evidence.\n', ...
    stage,nHits,sum(supported(1:nHits)));
if info.poolTruncated || info.joinLimitReached || info.hitLimitReached
    fprintf('Search stopped at a configured limit. Inspect SearchInfo and diagnostics.\n');
end
if string(o.regulationMode) == "design"
    fprintf('Regulatory assignments are design requirements; unannotated links remain pending.\n');
end

    function joinPool(pool,members)
        for candidateIndex = 1:size(pool,1)
            candidate = pool(candidateIndex,:);
            if stage == 3
                edges = candidate(4:end);
                extra = [candidate(1:3),graph.rxn(edges).',graph.dir(edges).'];
            else
                % Four production/removal pairs, each [signal,feed,prod,decay].
                nodes = reshape(candidate,4,4).';
                edges = reshape(nodes(:,3:4).',1,[]);
                extra = [nodes(:,1).',nodes(:,2).',graph.rxn(edges).',graph.dir(edges).'];
            end
            for u = members(:).'
                if tested >= o.maxJoinChecks || nHits >= o.maxHits, return; end
                tested = tested+1;
                if o.showProgress && toc(pulse) >= 5
                    fprintf('Block %d joins: %d checks | %d hits | %.1f min\n',stage,tested,nHits,toc(clock)/60);
                    pulse = tic;
                end
                values = [upValues(u,:),original(u),candidateIndex,extra];
                row = cell2struct(num2cell(values),names,2);
                if o.distinctAdapter
                    path = find([paths.pathIndex] == row.adapterPathIndex,1);
                    if isempty(path) || any(ismember(graph.rxn(edges),paths(path).reactionIndices))
                        rejected = rejected+1; continue;
                    end
                end
                % Distinct signals and reactions always; optional global cross-talk.
                schema = minimalCircuitSchema(stage,o.mapCsLeak);
                allM = cellfun(@(f) row.(f),schema.mets);
                allR = cellfun(@(f) row.(f),schema.rxns);
                if numel(unique(allM)) ~= numel(allM) || numel(unique(allR)) ~= numel(allR)
                    rejected = rejected+1; continue;
                end
                if o.crossTalk && ~checkMinimalCircuit(model,row,stage)
                    rejected = rejected+1; continue;
                end
                hasSupport = minimalRegulatoryEdges(row,stage,evidence);
                if string(o.regulationMode) == "evidence" && ~hasSupport
                    rejected = rejected+1; continue;
                end
                nHits = nHits+1;
                if nHits > size(rows,1)
                    capacity = min(o.maxHits,max(nHits,2*size(rows,1)));
                    rows(capacity,numel(names)) = 0; supported(capacity,1) = false;
                end
                rows(nHits,:) = values; supported(nHits) = hasSupport;
                if o.showProgress && toc(pulse) >= 5
                    fprintf('Block %d joins: %d checks | %d hits | %.1f min\n',stage,tested,nHits,toc(clock)/60);
                    pulse = tic;
                end
            end
            if nHits >= o.maxHits || tested >= o.maxJoinChecks, return; end
        end
    end
end

function [pool,stat] = localCounterPool(g,o,up,evidence)
% Cp: feed->C; Cg: C+Cs->waste; Cr: Cs->waste (regulated).
% Cs supply is external. An optional CsLeak implements basal removal.
numberEdges = 3+logical(o.mapCsLeak); width = 3+numberEdges;
pool = zeros(min(256,o.poolCap),width); hits = 0; checks = 0; cappedOptions = false;
hitCap = false; checkCap = false; pulse = tic;
regulatoryPairs = [];
if string(o.regulationMode) == "evidence"
    anchors = unique(up.m_d2);
    keys = evidence{:,{'MetIndex','RxnIndex','Direction','Effect'}};
    regulatoryPairs = false(size(anchors,1),numel(g.rxn));
    for a = 1:size(anchors,1)
        regulatoryPairs(a,:) = ismember([repmat(anchors(a),numel(g.rxn),1), ...
            g.rxn,g.dir,-ones(numel(g.rxn),1)],keys,'rows').';
    end
end
for gate = 1:numel(g.rxn)
    if g.excludedEdge(gate), continue; end
    substrates = g.substrates{gate};
    signals = substrates(~g.bad(substrates) & ~g.excludedMet(substrates));
    if o.allowCurrencyCs, csOptions = substrates(~g.excludedMet(substrates)); else, csOptions = signals; end
    if isempty(g.products{gate}), continue; end
    for c = signals
        for cs = csOptions
            if c == cs, continue; end
            if checks >= o.maxChecks, checkCap = true; finish(); return; end
            checks = checks+1;
            cp = localProducers(g,c,cs); cr = localRemovers(g,cs,c);
            cappedOptions = cappedOptions || any([numel(cp),numel(cr)] > o.optionsCap);
            cp = cp(1:min(end,o.optionsCap)); cr = cr(1:min(end,o.optionsCap));
            if o.mapCsLeak, basalOptions = cr; else, basalOptions = 0; end
            pairHits = 0;
            for producer = cp
                for remover = cr
                    if string(o.regulationMode) == "evidence" && ...
                            ~any(regulatoryPairs(:,producer) & regulatoryPairs(:,remover)), continue; end
                    for basal = basalOptions
                            if checks >= o.maxChecks, checkCap = true; finish(); return; end
                            checks = checks+1;
                            if o.showProgress && toc(pulse) >= 5
                                fprintf('Counter fragments: %d checks | %d candidates\n',checks,hits); pulse = tic;
                            end
                            edges = [producer gate remover];
                            intendedC = logical([0 1 0;0 1 1]);
                            intendedP = logical([1 0 0;0 0 0]);
                            if o.mapCsLeak
                                edges(end+1) = basal;
                                intendedC(:,end+1) = [false;true]; intendedP(:,end+1) = false;
                            end
                            if numel(unique(g.rxn(edges))) ~= numberEdges, continue; end
                            net = full(g.S([c cs],edges));
                            if ~isequal(net < 0,intendedC) || ~isequal(net > 0,intendedP), continue; end
                            hits = hits+1;
                            pairHits = pairHits+1;
                            if hits > size(pool,1), pool(min(o.poolCap,2*size(pool,1)),width) = 0; end
                            pool(hits,:) = [c cs,g.substrates{producer}(1),edges];
                            if hits >= o.poolCap, hitCap = true; finish(); return; end
                            if pairHits >= o.optionsCap, cappedOptions = true; break; end
                    end
                    if pairHits >= o.optionsCap, break; end
                end
                if pairHits >= o.optionsCap, break; end
            end
        end
    end
    if o.showProgress && toc(pulse) >= 5
        fprintf('Counter fragments: %.1f%% reactions | %d checks | %d candidates\n', ...
            100*gate/numel(g.rxn),checks,hits); pulse = tic;
    end
end
finish();
    function finish()
        pool = pool(1:hits,:);
        stat = [checks hits hitCap || cappedOptions checkCap];
    end
end

function producers = localProducers(g,m,other)
producers = g.producers{m};
keep = false(size(producers));
for k = 1:numel(producers)
    e = producers(k);
    keep(k) = ~g.excludedEdge(e) && ~isempty(g.substrates{e}) && g.S(other,e) == 0;
end
producers = producers(keep);
end

function removers = localRemovers(g,m,other)
removers = g.consumers{m};
keep = false(size(removers));
for k = 1:numel(removers)
    e = removers(k);
    keep(k) = ~g.excludedEdge(e) && ~isempty(g.products{e}) && g.S(other,e) == 0;
end
removers = removers(keep);
end

function [pool,stat] = localNNPool(g,o,c,evidence)
% Hill topology: C activates Ni production; Ni activates Na and Nb;
% Na and Nb activate output production. All four signals have removal.
units = zeros(min(256,o.nodeCap),4); nUnits = 0; checks = 0;
unitCap = false; hitCap = false; checkCap = false; optionsCap = false; pulse = tic;
for m = 1:numel(g.bad)
    if g.bad(m) || g.excludedMet(m), continue; end
    producers = g.producers{m}; removers = g.consumers{m};
    producers = producers(~g.excludedEdge(producers)); removers = removers(~g.excludedEdge(removers));
    producers = producers(cellfun(@(e) ~isempty(e),g.substrates(producers)));
    removers = removers(cellfun(@(e) ~isempty(e),g.products(removers)));
    optionsCap = optionsCap || numel(producers) > o.optionsCap || numel(removers) > o.optionsCap;
    producers = producers(1:min(end,o.optionsCap)); removers = removers(1:min(end,o.optionsCap));
    speciesUnits = 0;
    for producer = producers
        if string(o.regulationMode) == "evidence" && ...
                ~ismember([g.rxn(producer),g.dir(producer)], ...
                evidence{evidence.Effect == 1,{'RxnIndex','Direction'}},'rows'), continue; end
        for remover = removers
            if checks >= o.maxChecks, checkCap = true; break; end
            checks = checks+1;
            if g.rxn(producer) == g.rxn(remover), continue; end
            nUnits = nUnits+1;
            speciesUnits = speciesUnits+1;
            if nUnits > size(units,1), units(min(o.nodeCap,2*size(units,1)),4) = 0; end
            units(nUnits,:) = [m,g.substrates{producer}(1),producer,remover];
            if nUnits >= o.nodeCap, unitCap = true; break; end
            if speciesUnits >= o.optionsCap, optionsCap = true; break; end
        end
        if unitCap || checkCap || speciesUnits >= o.optionsCap, break; end
    end
    if unitCap || checkCap, break; end
end
units = units(1:nUnits,:);
pool = zeros(min(256,o.poolCap),16); hits = 0; selected = zeros(4,4);
% Role-specific filtering occurs before the NN combination limit.
walk(1);
pool = pool(1:hits,:);
stat = [checks hits hitCap || unitCap || optionsCap checkCap];

    function walk(role)
        if hits >= o.poolCap || checks >= o.maxChecks, return; end
        candidates = 1:nUnits;
        if string(o.regulationMode) == "evidence"
            if role == 1, regulators = c;
            elseif role <= 3, regulators = selected(1,1);
            else, regulators = selected(2:3,1).'; end
            allowed = true(numel(candidates),1);
            keys = evidence{:,{'MetIndex','RxnIndex','Direction','Effect'}};
            for regulator = regulators
                prod = units(:,3);
                allowed = allowed & ismember([repmat(regulator,nUnits,1),g.rxn(prod),g.dir(prod),ones(nUnits,1)],keys,'rows');
            end
            candidates = candidates(allowed);
        end
        if numel(candidates) > o.nodesPerRole, optionsCap = true; end
        candidates = localDiverseUnits(units,candidates,o.nodesPerRole);
        for u = candidates
            if checks >= o.maxChecks, checkCap = true; return; end
            checks = checks+1;
            selected(role,:) = units(u,:);
            met = selected(1:role,1).';
            edges = reshape(selected(1:role,3:4).',1,[]);
            if numel(unique(met)) ~= role || numel(unique(g.rxn(edges))) ~= 2*role, continue; end
            net = full(g.S(met,edges));
            intendedP = false(role,2*role); intendedC = intendedP;
            for k = 1:role, intendedP(k,2*k-1) = true; intendedC(k,2*k) = true; end
            if ~isequal(net > 0,intendedP) || ~isequal(net < 0,intendedC), continue; end
            if role < 4
                walk(role+1);
            else
                hits = hits+1;
                if hits > size(pool,1), pool(min(o.poolCap,2*size(pool,1)),16) = 0; end
                pool(hits,:) = reshape(selected.',1,[]);
                if hits >= o.poolCap, hitCap = true; return; end
            end
            if o.showProgress && toc(pulse) >= 5
                fprintf('NN fragments: %d checks | %d candidates\n',checks,hits); pulse = tic;
            end
            if hits >= o.poolCap || checkCap, return; end
        end
    end
end

function selected = localDiverseUnits(units,candidates,cap)
% Round-robin over signal metabolites: early rows must not fill all four
% neuron roles with dozens of variants of the same first metabolite.
selected = zeros(1,0);
metabolites = unique(units(candidates,1),'stable');
groups = cell(numel(metabolites),1);
for k = 1:numel(metabolites), groups{k} = candidates(units(candidates,1) == metabolites(k)); end
depth = 1;
while numel(selected) < cap
    added = false;
    for k = 1:numel(groups)
        if numel(groups{k}) >= depth
            selected(end+1) = groups{k}(depth); added = true;
            if numel(selected) >= cap, return; end
        end
    end
    if ~added, return; end
    depth = depth+1;
end
end
