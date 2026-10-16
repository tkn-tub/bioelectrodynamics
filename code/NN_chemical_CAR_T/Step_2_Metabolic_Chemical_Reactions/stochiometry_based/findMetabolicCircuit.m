function [motifs, info, p] = findMetabolicCircuit(model, options)
%FINDMETABOLICCIRCUIT Global constraint search for the two supplied diagrams.
% No optimization solver or old block files are needed. Direct connections
% are required; each reaction label maps to exactly one original reaction.
if nargin < 2, options = struct(); end
cfg = mergeOptions(options);
p = circuitPattern(cfg.NEURON_INPUT_MODE);
[D,meta] = prepareCircuitNetwork(model,cfg,p);
info.metadata = meta;
info.initialReactionCandidates = zeros(18,1);
info.initialMetaboliteCandidates = zeros(13,1);
info.nodes = 0; info.elapsedSeconds = 0; info.rootCandidates = 0;
info.parallelUsed = false; info.tasksCompleted = 0; info.tasksTotal = 0;
info.status = 'exhausted'; info.exhaustive = true; info.motifsReturned = 0;
info.oneMotifPerInput = cfg.ONE_MOTIF_PER_INPUT;
info.inputIdentity = cfg.METABOLITE_IDENTITY;
info.distinctInputsReturned = 0;
if cfg.ONE_MOTIF_PER_INPUT
    info.exhaustiveScope = 'One representative per feasible input, not all circuits for each input.';
else
    info.exhaustiveScope = 'All circuit assignments under the configured restrictions.';
end
motifs = emptyMotifs();
s.M = repmat({find(D.allowedMet)},13,1);
nin = full(sum(D.sub,1))'; nout = full(sum(D.prod,1))';
s.R = cell(18,1);
for j = 1:18
    s.R{j} = find(nin >= nnz(p.T(:,j) < 0) & nout >= nnz(p.T(:,j) > 0));
end
s = applyAnchors(s,D,meta,cfg);
info.initialReactionCandidates = cellfun(@numel,s.R);
info.initialMetaboliteCandidates = cellfun(@numel,s.M);
if cfg.VERBOSE
    fprintf('Model: %d reactions; %d metabolites.\n',D.nr,D.nm);
    fprintf('Complete EC annotations: %d; GPR annotations: %d.\n', ...
        meta.annotatedECReactions,meta.annotatedGPRReactions);
    fprintf('Eligible enzymatic conversions: %d; allowed orientations: %d.\n', ...
        meta.eligibleReactions,D.no);
    fprintf('Initial reaction candidate counts:\n');
    for j = 1:18, fprintf('  %-9s : %d\n',p.rxnNames{j},numel(s.R{j})); end
end
startClock = tic;
deadline = now + cfg.MAX_SECONDS/86400;
[s,ok,timedOut] = propagate(s,D,cfg,deadline);
if ~ok
    if timedOut, info.status = 'time_limit'; info.exhaustive = false; end
    info.elapsedSeconds = toc(startClock);
    if cfg.VERBOSE, fprintf('Initial constraint propagation: %s.\n',info.status); end
    return;
end
[kind,varIndex,values] = chooseVariable(s,D);
if isempty(values)
    if isCompleteMatch(s,D,cfg), motifs = makeMotif(s,D); end
    info.nodes = 1; info.motifsReturned = numel(motifs);
    info.distinctInputsReturned = numel(motifs);
    info.elapsedSeconds = toc(startClock); return;
end
if cfg.ONE_MOTIF_PER_INPUT
    % Partition by input identity, keeping all compartment aliases together.
    % Disjoint input groups go to disjoint workers: no duplicate input can
    % be found by another partition, and no post-search deduplication is needed.
    kind = 'M'; varIndex = find(strcmp(p.metNames,'m_in'),1);
    values = unique(D.metGroup(s.M{varIndex}),'stable');
    if cfg.VERBOSE
        fprintf('Distinct input groups to test: %d; stop after the first circuit per input.\n',numel(values));
    end
end
info.rootCandidates = numel(values);
pool = [];
if cfg.USE_PARALLEL
    if ~license('test','Distrib_Computing_Toolbox') || exist('parpool','file') ~= 2
        fprintf('Parallel toolbox unavailable: using serial search.\n');
    else
        pool = gcp('nocreate');
        if isempty(pool)
            if isempty(cfg.PARALLEL_PROFILE)
                if isempty(cfg.NUM_WORKERS), pool = parpool;
                else, pool = parpool(cfg.NUM_WORKERS); end
            else
                if isempty(cfg.NUM_WORKERS), pool = parpool(cfg.PARALLEL_PROFILE);
                else, pool = parpool(cfg.PARALLEL_PROFILE,cfg.NUM_WORKERS); end
            end
        elseif ~isempty(cfg.NUM_WORKERS) && pool.NumWorkers ~= cfg.NUM_WORKERS
            fprintf('Using existing pool with %d workers (requested %d).\n', ...
                pool.NumWorkers,cfg.NUM_WORKERS);
        end
    end
end
if isempty(pool)
    info.tasksTotal = 1;
    serialCfg = cfg; serialCfg.USE_PARALLEL = false;
    [motifs,stats] = runPartition(D,s,kind,varIndex,values,serialCfg,deadline);
    info.nodes = stats.nodes;
    info.tasksCompleted = 1;
    if stats.timedOut, info.status = 'time_limit'; info.exhaustive = false;
    elseif stats.hitLimit, info.status = 'motif_limit'; info.exhaustive = false;
    end
else
    info.parallelUsed = true;
    nTasks = min(numel(values),max(1,4*pool.NumWorkers));
    info.tasksTotal = nTasks;
    futures = parallel.FevalFuture.empty(0,1);
    for task = 1:nTasks
        % Disjoint root domains: every possible root assignment occurs once.
        subset = values(task:nTasks:end);
        futures(task,1) = parfeval(pool,@runPartition,2, ...
            D,s,kind,varIndex,subset,cfg,deadline);
    end
    cleanup = onCleanup(@() cancel(futures));
    processed = false(nTasks,1);
    anyLimit = false; anyTimeout = false;
    nextProgress = 30;
    while ~all(processed)
        for task = 1:nTasks
            if processed(task) || ~strcmp(futures(task).State,'finished'), continue; end
            if ~isempty(futures(task).Error)
                throw(futures(task).Error);
            end
            [found,stats] = fetchOutputs(futures(task));
            processed(task) = true;
            motifs = [motifs;found]; %#ok<AGROW>
            info.nodes = info.nodes + stats.nodes;
            info.tasksCompleted = info.tasksCompleted + 1;
            anyLimit = anyLimit || stats.hitLimit;
            anyTimeout = anyTimeout || stats.timedOut;
            if numel(motifs) >= cfg.MAX_MOTIFS, break; end
        end
        if numel(motifs) >= cfg.MAX_MOTIFS
            info.status = 'motif_limit'; info.exhaustive = false; break;
        end
        if now >= deadline
            % Workers stop at their next deadline check and return partial
            % results. Collect those results instead of discarding them.
            info.status = 'time_limit'; info.exhaustive = false;
        end
        if cfg.VERBOSE && toc(startClock) >= nextProgress
            fprintf('Search: %.0f s; %d/%d tasks complete; %d motifs received.\n', ...
                toc(startClock),nnz(processed),nTasks,numel(motifs));
            nextProgress = toc(startClock) + 30;
        end
        if ~all(processed), pause(0.1); end
    end
    if all(processed) && strcmp(info.status,'exhausted')
        if anyTimeout, info.status = 'time_limit'; info.exhaustive = false;
        elseif anyLimit, info.status = 'motif_limit'; info.exhaustive = false;
        end
    end
    clear cleanup; % Cancel unfinished futures, while preserving the user's pool.
end
if isfinite(cfg.MAX_MOTIFS) && numel(motifs) > cfg.MAX_MOTIFS
    motifs = motifs(1:cfg.MAX_MOTIFS);
    info.status = 'motif_limit'; info.exhaustive = false;
end
info.elapsedSeconds = toc(startClock);
info.motifsReturned = numel(motifs);
if ~isempty(motifs)
    inputRole = find(strcmp(p.metNames,'m_in'),1);
    inputs = arrayfun(@(m) m.metIndex(inputRole),motifs);
    info.distinctInputsReturned = numel(unique(D.metGroup(inputs)));
    if cfg.ONE_MOTIF_PER_INPUT
        assert(info.distinctInputsReturned == numel(motifs), ...
            'Circuit:DuplicateInputs','Distinct-input search returned a repeated input.');
    end
end
if cfg.VERBOSE
    fprintf('Search finished: %s; %d motifs; %.1f s; %d visited states.\n', ...
        info.status,numel(motifs),info.elapsedSeconds,info.nodes);
end
end

function cfg = mergeOptions(options)
cfg = circuitSearchDefaults();
names = fieldnames(options);
for k = 1:numel(names)
    if ~isfield(cfg,names{k}), error('Circuit:Option','Unknown option: %s',names{k}); end
    cfg.(names{k}) = options.(names{k});
end
choices = {'ENZYME_EVIDENCE',{'EC','EC_OR_GPR'}; ...
    'METABOLITE_SCOPE',{'all','roles'}; ...
    'METABOLITE_IDENTITY',{'chemical','model_id'}; ...
    'NEURON_INPUT_MODE',{'substrate','activation'}};
for k = 1:size(choices,1)
    if ~any(strcmpi(cfg.(choices{k,1}),choices{k,2}))
        error('Circuit:Option','Invalid value for %s.',choices{k,1});
    end
end
flag = cfg.ONE_MOTIF_PER_INPUT;
if ~(islogical(flag) || isnumeric(flag)) || ~isscalar(flag) || ~ismember(flag,[0 1])
    error('Circuit:Option','ONE_MOTIF_PER_INPUT must be true or false.');
end
if ~isnumeric(cfg.MAX_MOTIFS) || ~isscalar(cfg.MAX_MOTIFS) || ...
        isnan(cfg.MAX_MOTIFS) || cfg.MAX_MOTIFS < 1 || cfg.MAX_MOTIFS ~= floor(cfg.MAX_MOTIFS)
    error('Circuit:Option','MAX_MOTIFS must be a positive integer or Inf.');
end
if ~isnumeric(cfg.MAX_SECONDS) || ~isscalar(cfg.MAX_SECONDS) || ...
        isnan(cfg.MAX_SECONDS) || cfg.MAX_SECONDS <= 0
    error('Circuit:Option','MAX_SECONDS must be positive or Inf.');
end
if ~isscalar(cfg.TOLERANCE) || ~isfinite(cfg.TOLERANCE) || cfg.TOLERANCE <= 0
    error('Circuit:Option','TOLERANCE must be finite and positive.');
end
if ~isempty(cfg.NUM_WORKERS) && (~isnumeric(cfg.NUM_WORKERS) || ...
        ~isscalar(cfg.NUM_WORKERS) || ~isfinite(cfg.NUM_WORKERS) || ...
        cfg.NUM_WORKERS < 1 || cfg.NUM_WORKERS ~= floor(cfg.NUM_WORKERS))
    error('Circuit:Option','NUM_WORKERS must be empty or a positive integer.');
end
end

function s = applyAnchors(s,D,meta,cfg)
names = fieldnames(cfg.FIXED_METABOLITES);
for k = 1:numel(names)
    i = find(strcmp(D.p.metNames,names{k}),1);
    if isempty(i) && isfield(D.p,'legacyMetNames')
        i = find(strcmp(D.p.legacyMetNames,names{k}),1);
    end
    if isempty(i), error('Circuit:Anchor','Unknown metabolite role: %s',names{k}); end
    id = char(cfg.FIXED_METABOLITES.(names{k}));
    m = find(strcmp(meta.metIDs,id),1);
    if isempty(m), error('Circuit:Anchor','Unknown model metabolite: %s',id); end
    s.M{i} = intersect(s.M{i},m);
end
names = fieldnames(cfg.FIXED_REACTIONS);
for k = 1:numel(names)
    j = find(strcmp(D.p.rxnNames,names{k}),1);
    if isempty(j), error('Circuit:Anchor','Unknown reaction role: %s',names{k}); end
    id = char(cfg.FIXED_REACTIONS.(names{k}));
    r = find(strcmp(meta.rxnIDs,id),1);
    if isempty(r), error('Circuit:Anchor','Unknown model reaction: %s',id); end
    s.R{j} = intersect(s.R{j},find(D.origRxn == r));
end
end

function [s,ok,timedOut] = propagate(s,D,cfg,deadline)
% Positive/negative arcs get two-way support filtering. Absence constraints
% get forward checking when either end is assigned. All constraints are
% checked again at a complete assignment. Nothing is greedily committed.
ok = false; timedOut = false;
while true
    if now >= deadline, timedOut = true; return; end
    sizes = [cellfun(@numel,s.M);cellfun(@numel,s.R)];
    if any(sizes == 0), return; end
    before = sum(sizes);
    assignedM = find(cellfun(@numel,s.M) == 1);
    assignedR = find(cellfun(@numel,s.R) == 1);
    usedM = cellfun(@(x) x(1),s.M(assignedM));
    usedGroups = D.metGroup(usedM);
    if numel(unique(usedGroups)) ~= numel(usedGroups), return; end
    for i = 1:13
        if numel(s.M{i}) > 1
            s.M{i} = s.M{i}(~ismember(D.metGroup(s.M{i}),usedGroups));
        end
    end
    usedO = cellfun(@(x) x(1),s.R(assignedR));
    usedR = D.origRxn(usedO);
    if numel(unique(usedR)) ~= numel(usedR), return; end
    if nnz(D.conflict(usedR,usedR)) > numel(usedR), return; end
    for j = 1:18
        if numel(s.R{j}) > 1 && ~isempty(usedR)
            blocked = full(any(D.conflict(usedR,D.origRxn(s.R{j})),1))';
            s.R{j} = s.R{j}(~blocked);
        end
    end
    if any(cellfun(@isempty,s.M)) || any(cellfun(@isempty,s.R)), return; end
    [ii,jj] = find(D.p.T);
    for edge = 1:numel(ii)
        i = ii(edge); j = jj(edge);
        if D.p.T(i,j) < 0, Q = D.sub; else, Q = D.prod; end
        block = Q(s.M{i},s.R{j});
        keepM = full(any(block,2)); keepR = full(any(block,1))';
        s.M{i} = s.M{i}(keepM); s.R{j} = s.R{j}(keepR);
        if isempty(s.M{i}) || isempty(s.R{j}), return; end
    end
    % A role cannot secretly recur as an auxiliary species or compartment alias.
    for k = 1:numel(assignedM)
        i = assignedM(k); m = s.M{i}; g = D.metGroup(m);
        for j = 1:18
            if D.p.T(i,j) == 0
                keep = full(D.groupCount(g,s.R{j}) == 0)';
            elseif D.p.T(i,j) < 0
                keep = full(D.sub(m,s.R{j}) > 0)';
            else
                keep = full(D.prod(m,s.R{j}) > 0)';
            end
            s.R{j} = s.R{j}(keep);
            if isempty(s.R{j}), return; end
        end
    end
    for k = 1:numel(assignedR)
        j = assignedR(k); r = s.R{j};
        for i = 1:13
            if D.p.T(i,j) == 0
                keep = full(D.groupCount(D.metGroup(s.M{i}),r) == 0);
            elseif D.p.T(i,j) < 0
                keep = full(D.sub(s.M{i},r) > 0);
            else
                keep = full(D.prod(s.M{i},r) > 0);
            end
            s.M{i} = s.M{i}(keep);
            if isempty(s.M{i}), return; end
        end
    end
    if strcmpi(cfg.NEURON_INPUT_MODE,'activation')
        block = D.activation(s.M{10},s.R{16});
        keepM = full(any(block,2)); keepR = full(any(block,1))';
        s.M{10} = s.M{10}(keepM); s.R{16} = s.R{16}(keepR);
        if isempty(s.M{10}) || isempty(s.R{16}), return; end
    end
    if strcmpi(cfg.METABOLITE_SCOPE,'all') && numel(usedO) > 1
        % A compound shared by selected reactions must be a named role, or an
        % explicit auxiliary exemption. If only one role can cover it, bind it.
        repeated = find(full(sum(D.groupCount(:,usedO) > 0,2)) > 1);
        repeated = setdiff(repeated,D.exemptAux);
        for key = repeated(:)'
            roles = find(cellfun(@(x) any(D.metGroup(x) == key),s.M));
            if isempty(roles), return; end
            if numel(roles) == 1
                i = roles(1);
                s.M{i} = s.M{i}(D.metGroup(s.M{i}) == key);
            end
        end
    end
    afterSizes = [cellfun(@numel,s.M);cellfun(@numel,s.R)];
    if any(afterSizes == 0), return; end
    if sum(afterSizes) == before, ok = true; return; end
end
end

function [kind,index,values] = chooseVariable(s,D)
% Minimum remaining values; prefer high-degree variables on a tie.
sizes = [cellfun(@numel,s.M);cellfun(@numel,s.R)];
scores = sizes; scores(sizes <= 1) = Inf;
minimum = min(scores);
if isinf(minimum), kind = ''; index = []; values = []; return; end
ties = find(scores == minimum);
degree = [sum(D.p.T ~= 0,2);sum(D.p.T ~= 0,1)'];
[~,best] = max(degree(ties)); chosen = ties(best);
if chosen <= 13, kind = 'M'; index = chosen; values = s.M{index};
else, kind = 'R'; index = chosen-13; values = s.R{index}; end
end

function [motifs,stats] = runPartition(D,initial,kind,index,values,cfg,deadline)
motifs = emptyMotifs();
stats.nodes = 0; stats.timedOut = false; stats.hitLimit = false;
stats.rootsCompleted = 0;
progressClock = tic; nextProgress = 30;
inputFound = false;
for v = values(:)'
    inputFound = false;
    if stopped(), break; end
    child = initial;
    if cfg.ONE_MOTIF_PER_INPUT
        % v is an input group, not a model metabolite index. Try its allowed
        % species until the first valid circuit, then abandon the whole group.
        child.M{index} = initial.M{index}(D.metGroup(initial.M{index}) == v);
    else
        child.(kind){index} = v;
    end
    dfs(child);
    if ~stats.timedOut && ~stats.hitLimit
        stats.rootsCompleted = stats.rootsCompleted+1;
    end
end
    function yes = stopped()
        yes = false;
        if now >= deadline, stats.timedOut = true; yes = true;
        elseif numel(motifs) >= cfg.MAX_MOTIFS, stats.hitLimit = true; yes = true;
        elseif cfg.ONE_MOTIF_PER_INPUT && inputFound, yes = true;
        end
    end
    function dfs(state)
        if stopped(), return; end
        stats.nodes = stats.nodes+1;
        [state,consistent,timeout] = propagate(state,D,cfg,deadline);
        if timeout, stats.timedOut = true; return; end
        if ~consistent, return; end
        if cfg.VERBOSE && ~cfg.USE_PARALLEL && toc(progressClock) >= nextProgress
            fprintf('Search: %.0f s; %d visited states; %d motifs.\n', ...
                toc(progressClock),stats.nodes,numel(motifs));
            nextProgress = toc(progressClock)+30;
        end
        [nextKind,nextIndex,nextValues] = chooseVariable(state,D);
        if isempty(nextValues)
            if isCompleteMatch(state,D,cfg)
                motifs(end+1,1) = makeMotif(state,D); %#ok<AGROW>
                if cfg.ONE_MOTIF_PER_INPUT, inputFound = true; end
            end
            return;
        end
        for nextValue = nextValues(:)'
            if stopped(), return; end
            next = state; next.(nextKind){nextIndex} = nextValue;
            dfs(next);
        end
    end
end

function yes = isCompleteMatch(s,D,cfg)
yes = false;
if any(cellfun(@numel,s.M) ~= 1) || any(cellfun(@numel,s.R) ~= 1), return; end
mi = cellfun(@(x) x(1),s.M); oi = cellfun(@(x) x(1),s.R);
ri = D.origRxn(oi); groups = D.metGroup(mi);
if numel(unique(groups)) ~= 13 || numel(unique(ri)) ~= 18, return; end
if nnz(D.conflict(ri,ri)) ~= 18, return; end
actual = full(D.prod(mi,oi)-D.sub(mi,oi));
if ~isequal(actual,D.p.T), return; end
if ~isequal(full(D.groupCount(groups,oi)),abs(D.p.T)), return; end
if strcmpi(cfg.METABOLITE_SCOPE,'all')
    repeated = find(full(sum(D.groupCount(:,oi) > 0,2)) > 1);
    if ~isempty(setdiff(repeated,[groups;D.exemptAux])), return; end
end
if strcmpi(cfg.NEURON_INPUT_MODE,'activation') && D.activation(mi(10),oi(16)) == 0
    return;
end
yes = true;
end

function motif = makeMotif(s,D)
motif.metIndex = cellfun(@(x) x(1),s.M);
motif.orientedIndex = cellfun(@(x) x(1),s.R);
motif.reactionIndex = D.origRxn(motif.orientedIndex);
motif.direction = D.direction(motif.orientedIndex);
end
function motifs = emptyMotifs()
motifs = repmat(struct('metIndex',[],'orientedIndex',[], ...
    'reactionIndex',[],'direction',[]),0,1);
end
