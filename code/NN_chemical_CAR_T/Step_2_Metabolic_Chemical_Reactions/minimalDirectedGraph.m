function g = minimalDirectedGraph(model,bad,requireGPR)
[nm,nr] = size(model.S);
if numel(bad) ~= nm || numel(model.mets) ~= nm || numel(model.rxns) ~= nr
    error('MinimalGraph:Size','Model and currency mask dimensions do not match.');
end
g.eligible = getReactionECAnnotations(model);
if requireGPR
    [~,hasGPR] = prepareModelGPRCompatibility(model);
    g.eligible = g.eligible & hasGPR;
end
lb = zeros(nr,1); ub = inf(nr,1);
if isfield(model,'lb'), lb = model.lb(:); end
if isfield(model,'ub'), ub = model.ub(:); end
if numel(lb) ~= nr || numel(ub) ~= nr || any(isnan(lb) | isnan(ub) | lb > ub)
    error('MinimalGraph:Bounds','Invalid model bounds.');
end
forward = find(g.eligible & ub > 0); reverse = find(g.eligible & lb < 0);
g.rxn = [forward;reverse]; g.dir = [ones(numel(forward),1);-ones(numel(reverse),1)];
g.S = model.S(:,g.rxn)*spdiags(g.dir,0,numel(g.dir),numel(g.dir));
g.bad = logical(bad(:));
g.producers = cell(nm,1); g.consumers = cell(nm,1);
for m = 1:nm
    p = find(g.S(m,:) > 0); c = find(g.S(m,:) < 0);
    g.producers{m} = p(:).'; g.consumers{m} = c(:).';
end
g.substrates = cell(numel(g.rxn),1); g.products = g.substrates;
for e = 1:numel(g.rxn)
    g.substrates{e} = find(g.S(:,e) < 0).';
    g.products{e} = find(g.S(:,e) > 0).';
end
end
