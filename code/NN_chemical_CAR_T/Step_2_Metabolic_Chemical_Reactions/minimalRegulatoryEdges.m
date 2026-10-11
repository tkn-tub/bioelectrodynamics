function [allSupported,edges] = minimalRegulatoryEdges(row,stage,evidence)
s = minimalCircuitSchema(stage);
edges = zeros(numel(s.targets),5);
keys = evidence{:,{'MetIndex','RxnIndex','Direction','Effect'}};
for k = 1:numel(s.targets)
    key = [row.(s.regulators{k}),row.(s.targets{k}), ...
        row.(strrep(s.targets{k},'r_','dir_')),s.effects(k)];
    supported = ismember(key,keys,'rows');
    edges(k,:) = [key,supported];
end
allSupported = all(edges(:,5) ~= 0);
end
