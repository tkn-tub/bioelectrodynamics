function [motifs,diagnostics,paths,info] = mapMinimalDifferentiator(model,bad,paths,up,o)
%MAPMINIMALDIFFERENTIATOR Three distinct reactions per differentiator.
% Only selected reactions are tested for unintended internal consumption.
% External producers/consumers are not counted or prohibited.
% +1/-1 directions refer to model.S, subject to model.ub/model.lb.
validateattributes(model.S,{'numeric'},{'2d','real','finite'},mfilename,'model.S');
[nm,nr] = size(model.S);
if numel(model.mets) ~= nm || numel(model.rxns) ~= nr
    error('MinimalDiff:ModelSize','S, mets and rxns have inconsistent sizes.');
end
if numel(bad) ~= nm, error('MinimalDiff:CurrencyMask','Incorrect isBadMet length.'); end
bad = logical(bad(:));
for name = {'maxHits','maxInputs','perInputCap','maxChecks'}
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

stage1Names = {'adapterPathIndex','m_in','m_M1','m_Y1', ...
    'r_11','r_12','r_13','dir_11','dir_12','dir_13'};
original = [];
upValues = [];
if o.stage == 1
    inputs = unique(double([paths.anchorMetIdx]).','stable');
    names = stage1Names;
elseif o.stage == 2
    if ~istable(up) || ~all(ismember(stage1Names,up.Properties.VariableNames))
        error('MinimalDiff:OldBlock1', ...
            'Block 1 uses the old topology. Rerun the new Block 1 before Block 2.');
    end
    if isempty(up), error('MinimalDiff:EmptyBlock1','Block 1 has no motifs to extend.'); end
    localValidateUpstream(up,model,stage1Names);
    [up,original] = requireECAnnotatedRows(up,eligible,{'r_11','r_12','r_13'},paths);
    % Recheck direction feasibility, including on a fresh checkpoint load.
    keep = true(height(up),1);
    for k = 1:3
        role = {'r_11','r_12','r_13'}; direction = {'dir_11','dir_12','dir_13'};
        r = up.(role{k}); d = up.(direction{k});
        keep = keep & ((d == 1 & ub(r) > 0) | (d == -1 & lb(r) < 0));
    end
    up = up(keep,:); original = original(keep);
    if isempty(up), error('MinimalDiff:EmptyEligibleBlock1','No eligible Block 1 motifs remain.'); end
    upValues = table2array(up(:,stage1Names));
    inputs = unique(up.m_Y1,'stable');
    names = [{'block1_idx','input_idx'},stage1Names, ...
        {'m_X2','m_M2','m_Y2','m_d','r_21','r_22','r_23','dir_21','dir_22','dir_23'}];
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
nHits = 0; nDone = 0; rejected = 0; stats = zeros(nInputs,5);
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
    while nDone < nInputs && nHits < o.maxHits
        [slot,candidates,stat] = fetchNext(futures);
        a = futureInputs(slot);
        accept(a,candidates,stat);
        if nHits >= o.maxHits, break; end
        if next <= nInputs, submit(next); next = next + 1; end
    end
    clear workCleanup; % cancel all unused tasks and release shared graph
else
    for a = 1:nInputs
        [candidates,stat] = localSearch(graph,inputs(a), ...
            o.perInputCap,o.maxChecks,o.showProgress,o.stage,[]);
        accept(a,candidates,stat);
        if nHits >= o.maxHits, break; end
    end
end
motifs = array2table(rows(1:nHits,:),'VariableNames',names);
diagnostics = table(inputs,string(model.mets(inputs)),completed, ...
    stats(:,1),stats(:,2),stats(:,3),logical(stats(:,4)),logical(stats(:,5)), ...
    'VariableNames',{'InputIndex','InputID','Completed','CandidateChecks', ...
    'ChainPairs','CoreCandidates','PerInputHitLimitReached','CheckLimitReached'});
info = struct('topology',"minimal_three_reaction_v1",'nInputs',nInputs, ...
    'nInputsCompleted',nDone,'nJoinedHits',nHits,'nRejectedJoins',rejected, ...
    'globalHitLimitReached',nHits >= o.maxHits,'parallelUsed',parallelMode, ...
    'perInputSearchTruncated',any(stats(:,4:5) ~= 0,'all'), ...
    'inputLimitReached',o.stage == 1 && nInputs < numel(unique([paths.anchorMetIdx])), ...
    'elapsedSeconds',toc(clock));
fprintf('Block %d finished: %d/%d inputs; %d joined hits; %d rejected joins.\n', ...
    o.stage,nDone,nInputs,nHits,rejected);
if info.globalHitLimitReached
    fprintf('Stopped at MAX_HITS. Remaining parallel searches were cancelled.\n');
end
if info.perInputSearchTruncated
    fprintf('Per-input search was truncated; inspect diagnostics before interpreting zero hits.\n');
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
                if o.stage == 1
                    if o.distinctAdapter && any(ismember(candidate(4:6),paths(member).reactionIndices))
                        rejected = rejected + 1; continue;
                    end
                    row = [paths(member).pathIndex,candidate];
                else
                    prefix = upValues(member,:);
                    if o.distinctStageRxns && any(ismember(candidate(4:6),prefix(5:7)))
                        rejected = rejected + 1; continue;
                    end
                    if o.distinctStageMets && any(ismember(candidate(2:3),prefix(2:4)))
                        rejected = rejected + 1; continue;
                    end
                    p = find(pathIds == prefix(1),1);
                    if o.distinctAdapter && any(ismember(candidate(4:6),paths(p).reactionIndices))
                        rejected = rejected + 1; continue;
                    end
                    if o.crossTalk && ~localCascadeConsumption(model.S,prefix,candidate)
                        rejected = rejected + 1; continue;
                    end
                    row = [original(member),a,prefix,candidate(1:3),candidate(3),candidate(4:9)];
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
% Numeric columns: X M Y rA rB rSink dirA dirB dirSink.
rows = zeros(min(256,cap),9); hits = 0;
stat = zeros(1,5); pulse = tic;
inputConsumers = g.consumers{x};
for eA = inputConsumers
    for m = g.products{eA}
        if m == x || g.bad(m), continue; end
        for eB = g.consumers{m}
            if g.rxn(eA) == g.rxn(eB), continue; end
            for y = g.products{eB}
                if y == x || y == m || g.bad(y), continue; end
                % Count every potential (rA,M,rB,Y), including dead ends.
                if stat(1) >= maxChecks, stat(5) = 1; rows = rows(1:hits,:); return; end
                stat(1) = stat(1)+1; stat(2) = stat(2)+1;
                if show && toc(pulse) >= 5
                    message = sprintf('Block %d input %d: %d checks, %d core candidates...', ...
                        stage,x,stat(1),hits);
                    if isempty(queue), fprintf('%s\n',message); else, send(queue,message); end
                    pulse = tic;
                end
                sinks = intersect(inputConsumers,g.consumers{y});
                for eSink = sinks(:).'
                    if stat(1) >= maxChecks, stat(5) = 1; rows = rows(1:hits,:); return; end
                    stat(1) = stat(1)+1;
                    e = [eA eB eSink];
                    if numel(unique(g.rxn(e))) ~= 3 || isempty(g.products{eSink}), continue; end
                    if ~localCoreConsumption(g.S,[x m y],e), continue; end
                    hits = hits+1;
                    stat(3) = hits;
                    if hits > size(rows,1), rows(min(cap,2*size(rows,1)),9) = 0; end
                    rows(hits,:) = [x m y,g.rxn(e).',g.dir(e).'];
                    if hits >= cap
                        stat(3) = hits; stat(4) = 1; rows = rows(1:hits,:); return;
                    end
                end
            end
        end
    end
end
stat(3) = hits;
rows = rows(1:hits,:);
end

function passes = localCoreConsumption(S,m,r)
consumes = full(S(m,r)) < 0;
intended = logical([1 0 1;0 1 0;0 0 1]);
passes = isequal(consumes,intended);
% Signal species have only their intended producers inside the algorithm.
% Other metabolites/reactions in the network remain unrestricted.
produces = full(S(m,r)) > 0;
passes = passes && isequal(produces,logical([0 0 0;1 0 0;0 1 0]));
end

function passes = localCascadeConsumption(S,p,c)
m = [p(2:4),c(2:3)];
r = [p(5:7),c(4:6)]; dirs = [p(8:10),c(7:9)];
signed = full(S(m,r)).*dirs;
consumes = signed < 0;
produces = signed > 0;
allowed = logical([1 0 1 0 0 0;0 1 0 0 0 0;0 0 1 1 0 1; ...
    0 0 0 0 1 0;0 0 0 0 0 1]);
allowedProduction = logical([0 0 0 0 0 0;1 0 0 0 0 0;0 1 0 0 0 0; ...
    0 0 0 1 0 0;0 0 0 0 1 0]);
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
m = values(:,2:4); r = values(:,5:7); dirs = values(:,8:10);
if any(m < 1 | m > nm | m ~= floor(m),'all') || ...
        any(r < 1 | r > nr | r ~= floor(r),'all') || any(~ismember(dirs,[-1 1]),'all')
    error('MinimalDiff:ModelMismatch','Block 1 indices/directions do not match the loaded model.');
end
uniqueRows = unique(values(:,2:10),'rows');
for k = 1:size(uniqueRows,1)
    u = uniqueRows(k,:);
    coreS = model.S(u(1:3),u(4:6))*diag(u(7:9));
    sinkProducts = model.S(:,u(6))*u(9) > 0;
    if numel(unique(u(1:3))) ~= 3 || numel(unique(u(4:6))) ~= 3 || ...
            ~localCoreConsumption(coreS,1:3,1:3) || ~any(sinkProducts)
        error('MinimalDiff:ModelMismatch','An upstream motif does not realize the three-reaction topology.');
    end
end
end

function text = localGPRText(required)
if required, text = '+GPR'; else, text = ''; end
end
