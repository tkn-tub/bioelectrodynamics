function [motifs,diagnostics,paths,info] = mapFiveReactionDifferentiators(model,bad,paths,up,o)
% Stage 1: X->M1->Y1; X->F1; F1+Y1->products; F1->products.
% Stage 2: Y1->M2->Y2; Y1->Y2; M2+Y2->products; Y2->products.
% Each stage contains five distinct reactions, with independent rate tuning.
% Only selected reactions are tested for unintended internal consumption.
% External producers/consumers are not counted or prohibited.
% +1/-1 directions refer to model.S, subject to model.ub/model.lb.
validateattributes(model.S,{'numeric'},{'2d','real','finite'},mfilename,'model.S');
if ~ismember(o.stage,[1 2]), error('PolarityDiff:Stage','stage must be 1 or 2.'); end
[nm,nr] = size(model.S);
if numel(model.mets) ~= nm || numel(model.rxns) ~= nr
    error('MinimalDiff:ModelSize','S, mets and rxns have inconsistent sizes.');
end
if numel(bad) ~= nm, error('MinimalDiff:CurrencyMask','Incorrect isBadMet length.'); end
bad = logical(bad(:));
if ~isfield(o,'maxJoinChecks'), o.maxJoinChecks = 2e6; end
for name = {'maxHits','maxInputs','perInputCap','maxChecks','maxJoinChecks'}
    value = o.(name{1});
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
            isnan(value) || value <= 0 || value ~= floor(value)
        error('MinimalDiff:Budget','%s must be a positive integer or inf.',name{1});
    end
end
[ec,~] = getReactionECAnnotations(model);
eligible = ec;
if o.requireGPR
    [~,gpr] = prepareModelGPRCompatibility(model);
    eligible = eligible & gpr;
end
fprintf('Block %d eligible EC%s reactions: %d / %d\n', ...
    o.stage,localGPRText(o.requireGPR),sum(eligible),nr);
paths = requireECAnnotatedAdapters(model,paths,ec,o.requireAdapterEC);
if isempty(paths), error('MinimalDiff:NoAdapters','No eligible Block 0 paths are available.'); end
if ~isfield(paths,'pathIndex')
    for p = 1:numel(paths), paths(p).pathIndex = p; end
end
pathIds = double([paths.pathIndex]).';
validateattributes(pathIds,{'numeric'},{'positive','integer','finite'});
if numel(unique(pathIds)) ~= numel(pathIds)
    error('MinimalDiff:AdapterIdentity','Block 0 pathIndex values must be unique.');
end

% Build only feasible, annotated directed reactions. Never modify model.S.
lb = zeros(nr,1); ub = inf(nr,1);
if isfield(model,'lb'), lb = model.lb(:); end
if isfield(model,'ub'), ub = model.ub(:); end
if numel(lb) ~= nr || numel(ub) ~= nr || any(isnan(lb) | isnan(ub) | lb > ub)
    error('MinimalDiff:Bounds','Invalid model reaction bounds.');
end
forward = find(eligible & ub > 0);
reverse = find(eligible & lb < 0);
graph.rxn = [forward;reverse];
graph.dir = [ones(numel(forward),1);-ones(numel(reverse),1)];
graph.S = model.S(:,graph.rxn)*spdiags(graph.dir,0,numel(graph.dir),numel(graph.dir));
graph.bad = bad;
graph.consumers = cell(nm,1);
for m = 1:nm
    edges = find(graph.S(m,:) < 0); graph.consumers{m} = edges(:).';
end
graph.products = cell(numel(graph.rxn),1);
for e = 1:numel(graph.rxn), graph.products{e} = find(graph.S(:,e) > 0).'; end

stage1Names = {'adapterPathIndex','m_in','m_M1','m_Y1','m_F1_delay', ...
    'r_S1','r_S2','r_F','r_I','r_R','dir_S1','dir_S2','dir_F','dir_I','dir_R'};
original = [];
upValues = [];
if o.stage == 1
    inputs = unique(double([paths.anchorMetIdx]).','stable');
    names = stage1Names;
elseif o.stage == 2
    if ~istable(up) || ~all(ismember(stage1Names,up.Properties.VariableNames))
        error('MinimalDiff:OldBlock1', ...
            'Run the five-reaction Block 1 before Block 2. Three-reaction rows cannot be extended.');
    end
    if isempty(up), error('MinimalDiff:EmptyBlock1','Block 1 has no motifs to extend.'); end
    localValidateUpstream(up,model,stage1Names);
    [up,original] = requireECAnnotatedRows(up,eligible,{'r_S1','r_S2','r_F','r_I','r_R'},paths);
    % Recheck direction feasibility, including on a fresh checkpoint load.
    keep = true(height(up),1);
    for k = 1:5
        role = {'r_S1','r_S2','r_F','r_I','r_R'}; direction = strrep(role,'r_','dir_');
        r = up.(role{k}); d = up.(direction{k});
        keep = keep & ((d == 1 & ub(r) > 0) | (d == -1 & lb(r) < 0));
    end
    up = up(keep,:); original = original(keep);
    if isempty(up), error('MinimalDiff:EmptyEligibleBlock1','No eligible Block 1 motifs remain.'); end
    upValues = table2array(up(:,stage1Names));
    inputs = unique(up.m_Y1,'stable');
    names = [{'block1_idx','input_idx'},stage1Names, ...
        {'m_d','m_d1','m_X2','m_M2','m_Y2','m_d2', ...
         'r_S1_2','r_S2_2','r_F2','r_I2','r_R2', ...
         'dir_S1_2','dir_S2_2','dir_F2','dir_I2','dir_R2'}];
else
    error('MinimalDiff:Stage','stage must be 1 or 2.');
end
if isfinite(o.maxInputs), inputs = inputs(1:min(numel(inputs),o.maxInputs)); end
nInputs = numel(inputs);
inputGroups = cell(nInputs,1);
for a = 1:nInputs
    if o.stage == 1
        inputGroups{a} = find(double([paths.anchorMetIdx]).' == inputs(a));
    else
        inputGroups{a} = find(up.m_Y1 == inputs(a));
    end
end

% Joined hits, not raw chains, consume the global budget. Grow in chunks.
capacity = min(1000,o.maxHits);
rows = zeros(capacity,numel(names));
nHits = 0; nDone = 0; rejected = 0; joinChecks = 0; stats = zeros(nInputs,5);
completed = false(nInputs,1);
pool = [];
if o.useParallel && exist('gcp','file') == 2
    pool = gcp('nocreate');
    if isempty(pool) && o.autoPool
        if isempty(o.workers), pool = parpool; else, pool = parpool(o.workers); end
    end
end
parallelMode = ~isempty(pool);
if o.useParallel && ~parallelMode
    warning('MinimalDiff:SerialFallback','No allocated MATLAB pool; using serial search.');
end
clock = tic;
if parallelMode
    % Keep at most NumWorkers inputs in flight. Cancel outstanding work when
    % MAX_HITS is reached; do not submit every input before checking the cap.
    graphConstant = parallel.pool.Constant(graph);
    queue = [];
    if o.showProgress
        queue = parallel.pool.DataQueue;
        afterEach(queue,@(message) fprintf('%s\n',message));
    end
    futures = parallel.FevalFuture.empty(0,1);
    futureInputs = zeros(nInputs,1);
    workCleanup = onCleanup(@cleanWorkers);
    next = 1;
    for k = 1:min(nInputs,pool.NumWorkers)
        submit(next); next = next + 1;
    end
    while nDone < nInputs && nHits < o.maxHits && joinChecks < o.maxJoinChecks
        [slot,candidates,stat] = fetchNext(futures);
        a = futureInputs(slot);
        accept(a,candidates,stat);
        if nHits >= o.maxHits || joinChecks >= o.maxJoinChecks, break; end
        if next <= nInputs, submit(next); next = next + 1; end
    end
    clear workCleanup; % cancel all unused tasks and release shared graph
else
    for a = 1:nInputs
        [candidates,stat] = localSearch(graph,inputs(a), ...
            o.perInputCap,o.maxChecks,o.showProgress,o.stage,[]);
        accept(a,candidates,stat);
        if nHits >= o.maxHits || joinChecks >= o.maxJoinChecks, break; end
    end
end
motifs = array2table(rows(1:nHits,:),'VariableNames',names);
diagnostics = table(inputs,string(model.mets(inputs)),completed, ...
    stats(:,1),stats(:,2),stats(:,3),logical(stats(:,4)),logical(stats(:,5)), ...
    'VariableNames',{'InputIndex','InputID','Completed','CandidateChecks', ...
    'BranchOrChainPairs','CoreCandidates','PerInputHitLimitReached','CheckLimitReached'});
info = struct('topology',"five_reaction_differentiators_v2",'nInputs',nInputs, ...
    'nInputsCompleted',nDone,'nJoinedHits',nHits,'nRejectedJoins',rejected, ...
    'globalHitLimitReached',nHits >= o.maxHits,'parallelUsed',parallelMode, ...
    'joinChecks',joinChecks,'joinCheckLimitReached',joinChecks >= o.maxJoinChecks, ...
    'perInputSearchTruncated',any(stats(:,4:5) ~= 0,'all'), ...
    'inputLimitReached',o.stage == 1 && nInputs < numel(unique([paths.anchorMetIdx])), ...
    'elapsedSeconds',toc(clock));
fprintf('Block %d finished: %d/%d inputs; %d joined hits; %d rejected joins.\n', ...
    o.stage,nDone,nInputs,nHits,rejected);
if info.globalHitLimitReached
    fprintf('Stopped at MAX_HITS.\n');
    if parallelMode, fprintf('Remaining parallel searches were cancelled.\n'); end
end
if info.perInputSearchTruncated
    fprintf('Per-input search was truncated; inspect diagnostics before interpreting zero hits.\n');
end
if info.joinCheckLimitReached
    fprintf('Stopped at MAX_DIFF_JOIN_CHECKS. Rejected joins also count toward this budget.\n');
end

    function submit(a)
        futures(end+1,1) = parfeval(pool,@localSearchConstant,2, ...
            graphConstant,inputs(a),o.perInputCap,o.maxChecks,o.stage,queue);
        futureInputs(numel(futures)) = a;
    end

    function cleanWorkers()
        cancel(futures);
        delete(graphConstant);
    end

    function accept(a,candidates,stat)
        stats(a,:) = stat; completed(a) = true; nDone = nDone + 1;
        members = inputGroups{a};
        for c = 1:size(candidates,1)
            candidate = candidates(c,:);
            for member = members(:).'
                if nHits >= o.maxHits || joinChecks >= o.maxJoinChecks, return; end
                joinChecks = joinChecks+1;
                if o.stage == 1
                    if o.distinctAdapter && any(ismember(candidate(5:9),paths(member).reactionIndices))
                        rejected = rejected + 1; continue;
                    end
                    row = [paths(member).pathIndex,candidate];
                else
                    prefix = upValues(member,:);
                    if o.distinctStageRxns && any(ismember(candidate(4:8),prefix(6:10)))
                        rejected = rejected + 1; continue;
                    end
                    if o.distinctStageMets && any(ismember(candidate(2:3),prefix(2:5)))
                        rejected = rejected + 1; continue;
                    end
                    p = find(pathIds == prefix(1),1);
                    if o.distinctAdapter && any(ismember(candidate(4:8),paths(p).reactionIndices))
                        rejected = rejected + 1; continue;
                    end
                    if o.crossTalk && ~localCascadeConsumption(model.S,prefix,candidate)
                        rejected = rejected + 1; continue;
                    end
                    row = [original(member),a,prefix,prefix(4),prefix(4),candidate(1:3),candidate(3),candidate(4:13)];
                end
                nHits = nHits + 1;
                if nHits > size(rows,1)
                    rows(min(o.maxHits,max(nHits,2*size(rows,1))),numel(names)) = 0;
                end
                rows(nHits,:) = row;
                if nHits >= o.maxHits, break; end
            end
            if nHits >= o.maxHits, break; end
        end
        if o.showProgress
            fprintf('Block %d: %5.1f%% inputs | %d/%d | joined hits %d | %.1f min\n', ...
                o.stage,100*nDone/nInputs,nDone,nInputs,nHits,toc(clock)/60);
        end
    end
end

function [rows,stat] = localSearchConstant(graphConstant,input,cap,checks,stage,queue)
[rows,stat] = localSearch(graphConstant.Value,input,cap,checks,~isempty(queue),stage,queue);
end

function [rows,stat] = localSearch(g,x,cap,maxChecks,show,stage,queue)
if stage == 2
    [rows,stat] = localSearchPositive(g,x,cap,maxChecks,show,queue);
    return;
end
% Numeric columns: X M Y F [S1 S2 F I R] [five directions].
rows = zeros(min(256,cap),14); hits = 0;
stat = zeros(1,5); pulse = tic;
inputConsumers = g.consumers{x};
for slow = inputConsumers
    if ~tick(), finish(); return; end
    for m = g.products{slow}
        if m == x || g.bad(m), continue; end
        for turnover = g.consumers{m}
            if ~tick(), finish(); return; end
            if g.rxn(slow) == g.rxn(turnover), continue; end
            for y = g.products{turnover}
                if y == x || y == m || g.bad(y), continue; end
                if ~tick(), finish(); return; end
                for fast = inputConsumers
                    if g.rxn(fast) == g.rxn(slow) || g.rxn(fast) == g.rxn(turnover), continue; end
                    for f = g.products{fast}
                        if g.bad(f) || ismember(f,[x m y]), continue; end
                        if ~tick(), finish(); return; end
                        stat(2) = stat(2)+1;
                        mets = [x m y f];
                        if ~localColumn(g.S,mets,slow,[1;0;0;0],[0;1;0;0]) || ...
                                ~localColumn(g.S,mets,turnover,[0;1;0;0],[0;0;1;0]) || ...
                                ~localColumn(g.S,mets,fast,[1;0;0;0],[0;0;0;1])
                            continue;
                        end
                        jointOptions = intersect(g.consumers{f},g.consumers{y});
                        jointOptions = localRoleEdges(g,mets,jointOptions,[0;0;1;1],zeros(4,1));
                        removalOptions = localRoleEdges(g,mets,g.consumers{f},[0;0;0;1],zeros(4,1));
                        for joint = jointOptions
                            for removal = removalOptions
                                if ~tick(), finish(); return; end
                                edges = [slow turnover fast joint removal];
                                if numel(unique(g.rxn(edges))) ~= 5, continue; end
                                if ~localCoreConsumption(g.S,mets,edges), continue; end
                                hits = hits+1; stat(3) = hits;
                                if hits > size(rows,1), rows(min(cap,2*size(rows,1)),14) = 0; end
                                rows(hits,:) = [x m y f,g.rxn(edges).',g.dir(edges).'];
                                if hits >= cap, stat(4) = 1; finish(); return; end
                            end
                        end
                    end
                end
            end
        end
    end
end
finish();

    function allowed = tick()
        allowed = stat(1) < maxChecks;
        if ~allowed, stat(5) = 1; return; end
        stat(1) = stat(1)+1;
        if show && toc(pulse) >= 5
            message = sprintf('Block 1 input %d: %d checks, %d five-reaction candidates...',x,stat(1),hits);
            if isempty(queue), fprintf('%s\n',message); else, send(queue,message); end
            pulse = tic;
        end
    end

    function finish()
        rows = rows(1:hits,:);
    end
end

function passes = localCoreConsumption(S,m,r)
consumes = full(S(m,r)) < 0;
intended = logical([1 0 1 0 0;0 1 0 0 0;0 0 0 1 0;0 0 0 1 1]);
passes = isequal(consumes,intended);
% Signal species have only their intended producers inside the algorithm.
% Other metabolites/reactions in the network remain unrestricted.
produces = full(S(m,r)) > 0;
passes = passes && isequal(produces,logical([0 0 0 0 0;1 0 0 0 0;0 1 0 0 0;0 0 1 0 0]));
end

function yes = localColumn(S,mets,edge,consumption,production)
signed = full(S(mets,edge));
yes = isequal(signed < 0,logical(consumption)) && isequal(signed > 0,logical(production));
end

function edges = localRoleEdges(g,mets,edges,consumption,production)
if isempty(edges), return; end
signed = full(g.S(mets,edges));
valid = all((signed < 0) == logical(consumption),1) & all((signed > 0) == logical(production),1);
hasProducts = reshape(cellfun(@(p) ~isempty(p),g.products(edges)),1,[]);
edges = reshape(edges(valid & hasProducts),1,[]);
end

function [rows,stat] = localSearchPositive(g,x,cap,maxChecks,show,queue)
% Numeric columns: X2 M2 Y2 [S1_2 S2_2 F2 I2 R2] [five directions].
% S2_2 converts M2 to Y2. F2 also directly produces Y2 from X2.
rows = zeros(min(256,cap),13); hits = 0;
stat = zeros(1,5); pulse = tic;
slowProducers = g.consumers{x}; fastProducers = slowProducers;
for slow = slowProducers
    if ~tick(), finish(); return; end
    for m = g.products{slow}
        if m == x || g.bad(m), continue; end
        for fast = fastProducers
            if g.rxn(slow) == g.rxn(fast), continue; end
            for y = g.products{fast}
                if y == x || y == m || g.bad(y), continue; end
                if ~tick(), finish(); return; end
                stat(2) = stat(2)+1;
                mets = [x m y];
                if ~isequal(full(g.S(mets,slow)) < 0,logical([1;0;0])) || ...
                        ~isequal(full(g.S(mets,slow)) > 0,logical([0;1;0])) || ...
                        ~isequal(full(g.S(mets,fast)) < 0,logical([1;0;0])) || ...
                        ~isequal(full(g.S(mets,fast)) > 0,logical([0;0;1]))
                    continue;
                end
                slowRemoval = localRoleEdges(g,mets,g.consumers{m},[0;1;0],[0;0;1]);
                sharedRemoval = localRoleEdges(g,mets,intersect(g.consumers{m},g.consumers{y}),[0;1;1],[0;0;0]);
                outputRemoval = localRoleEdges(g,mets,g.consumers{y},[0;0;1],[0;0;0]);
                if isempty(slowRemoval) || isempty(sharedRemoval) || isempty(outputRemoval), continue; end
                for turnover = slowRemoval
                    for joint = sharedRemoval
                        for leak = outputRemoval
                            if ~tick(), finish(); return; end
                            edges = [slow turnover fast joint leak];
                            if numel(unique(g.rxn(edges))) ~= 5, continue; end
                            signed = full(g.S(mets,edges));
                            consumes = logical([1 0 1 0 0;0 1 0 1 0;0 0 0 1 1]);
                            produces = logical([0 0 0 0 0;1 0 0 0 0;0 1 1 0 0]);
                            if ~isequal(signed < 0,consumes) || ~isequal(signed > 0,produces), continue; end
                            hits = hits+1; stat(3) = hits;
                            if hits > size(rows,1), rows(min(cap,2*size(rows,1)),13) = 0; end
                            rows(hits,:) = [x m y,g.rxn(edges).',g.dir(edges).'];
                            if hits >= cap, stat(4) = 1; finish(); return; end
                        end
                    end
                end
            end
        end
    end
end
finish();

    function allowed = tick()
        allowed = stat(1) < maxChecks;
        if ~allowed, stat(5) = 1; return; end
        stat(1) = stat(1)+1;
        if show && toc(pulse) >= 5
            message = sprintf('Block 2 input %d: %d checks, %d five-reaction candidates...',x,stat(1),hits);
            if isempty(queue), fprintf('%s\n',message); else, send(queue,message); end
            pulse = tic;
        end
    end

    function finish()
        rows = rows(1:hits,:);
    end
end

function passes = localCascadeConsumption(S,p,c)
m = [p(2:5),c(2:3)];
r = [p(6:10),c(4:8)]; dirs = [p(11:15),c(9:13)];
signed = full(S(m,r)).*dirs;
consumes = signed < 0;
produces = signed > 0;
% Six physical species and ten reaction roles; only Y1 joins the stages.
schema = minimalCircuitSchema(2);
allowed = schema.consumes;
allowedProduction = schema.produces;
% If roles deliberately reuse a metabolite/reaction, combine the allowed
% consumers for that physical metabolite and directed reaction.
for i = 1:numel(m)
    merged = any(allowed(m == m(i),:),1);
    allowedEdges = [r(merged).',dirs(merged).'];
    if any(consumes(i,:) & ~ismember([r.',dirs.'],allowedEdges,'rows').')
        passes = false; return;
    end
    mergedProducers = any(allowedProduction(m == m(i),:),1);
    producerEdges = [r(mergedProducers).',dirs(mergedProducers).'];
    if any(produces(i,:) & ~ismember([r.',dirs.'],producerEdges,'rows').')
        passes = false; return;
    end
end
passes = true;
end

function localValidateUpstream(up,model,names)
values = table2array(up(:,names));
if ~isnumeric(values) || ~isreal(values) || any(~isfinite(values),'all')
    error('MinimalDiff:InvalidBlock1','Block 1 indices must be finite numbers.');
end
nm = size(model.S,1); nr = size(model.S,2);
m = values(:,2:5); r = values(:,6:10); dirs = values(:,11:15);
if any(m < 1 | m > nm | m ~= floor(m),'all') || ...
        any(r < 1 | r > nr | r ~= floor(r),'all') || any(~ismember(dirs,[-1 1]),'all')
    error('MinimalDiff:ModelMismatch','Block 1 indices/directions do not match the loaded model.');
end
uniqueRows = unique(values(:,2:15),'rows');
for k = 1:size(uniqueRows,1)
    u = uniqueRows(k,:);
    coreS = model.S(u(1:4),u(5:9))*diag(u(10:14));
    jointProducts = model.S(:,u(8))*u(13) > 0;
    removalProducts = model.S(:,u(9))*u(14) > 0;
    if numel(unique(u(1:4))) ~= 4 || numel(unique(u(5:9))) ~= 5 || ...
            ~localCoreConsumption(coreS,1:4,1:5) || ~any(jointProducts) || ~any(removalProducts)
        error('MinimalDiff:ModelMismatch','An upstream motif does not realize the five-reaction Block 1 topology.');
    end
end
end

function text = localGPRText(required)
if required, text = '+GPR'; else, text = ''; end
end
