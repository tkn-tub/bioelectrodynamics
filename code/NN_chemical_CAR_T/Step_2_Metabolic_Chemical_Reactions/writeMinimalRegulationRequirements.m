function requirements = writeMinimalRegulationRequirements(motifs,model,evidence)
% Unique requested interactions, not repeated text for every motif row.
s = minimalCircuitSchema(4); keys = zeros(0,4);
for k = 1:numel(s.targets)
    n = height(motifs);
    keys = [keys;motifs.(s.regulators{k}),motifs.(s.targets{k}), ...
        motifs.(strrep(s.targets{k},'r_','dir_')),repmat(s.effects(k),n,1)]; %#ok<AGROW>
end
keys = unique(keys,'rows');
supported = ismember(keys,evidence{:,{'MetIndex','RxnIndex','Direction','Effect'}},'rows');
effect = repmat("activation",size(keys,1),1); effect(keys(:,4) == -1) = "inhibition";
sources = strings(size(keys,1),1);
for k = 1:size(keys,1)
    match = ismember(evidence{:,{'MetIndex','RxnIndex','Direction','Effect'}},keys(k,:),'rows');
    if any(match), sources(k) = strjoin(evidence.Evidence(match),'; '); end
end
requirements = table(string(model.mets(keys(:,1))),string(model.rxns(keys(:,2))), ...
    keys(:,3),effect,sources,supported,'VariableNames', ...
    {'RegulatorMetID','TargetRxnID','TargetDirection','Effect','Evidence','Supported'});
writetable(requirements,'Regulation_To_Validate.csv');
end
