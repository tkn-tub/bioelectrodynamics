function [ok,reason] = checkMinimalCircuit(model,row,stage)
% Only selected reactions are inspected. External network edges are allowed.
s = minimalCircuitSchema(stage,isfield(row,'r_CsLeak'));
m = zeros(1,numel(s.mets)); r = zeros(1,numel(s.rxns)); d = r;
for k = 1:numel(m), m(k) = row.(s.mets{k}); end
for k = 1:numel(r), r(k) = row.(s.rxns{k}); d(k) = row.(s.dirs{k}); end
if any(~isfinite([m r d])) || any(m < 1 | m > size(model.S,1) | m ~= floor(m)) || ...
        any(r < 1 | r > size(model.S,2) | r ~= floor(r)) || any(~ismember(d,[-1 1]))
    error('MinimalCircuit:ModelIndex','Invalid metabolite/reaction index or direction.');
end
ok = false; reason = "";
if numel(unique(m)) ~= numel(m)
    reason = "Protected signal metabolites overlap"; return;
end
if numel(unique(r)) ~= numel(r)
    reason = "Selected reaction IDs overlap"; return;
end
net = full(model.S(m,r)).*d;
badConsumption = (net < 0) ~= s.consumes;
badProduction = (net > 0) ~= s.produces;
[i,j] = find(badConsumption | badProduction,1);
if ~isempty(i)
    reason = sprintf('%s has an unintended or missing edge in %s',s.mets{i},s.rxns{j});
    return;
end
ok = true;
end
