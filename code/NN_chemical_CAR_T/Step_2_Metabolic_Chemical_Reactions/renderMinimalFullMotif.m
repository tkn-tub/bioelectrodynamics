function text = renderMinimalFullMotif(data,selection)
% Print one full circuit using checkpoint dictionaries, without model/COBRA.
tableNames = {'Motifs_NN_Block_6','Motifs_NN_Block_5','Motifs_NN_Block_4'};
motifs = [];
for k = 1:numel(tableNames)
    if isfield(data,tableNames{k}), motifs = data.(tableNames{k}); break; end
end
if ~istable(motifs), error('MinimalPrint:Table','No full motif table in this checkpoint.'); end
validateattributes(selection,{'numeric'},{'scalar','integer','positive','<=',height(motifs)});
row = table2struct(motifs(selection,:)); schema = minimalCircuitSchema(4);
lines = {sprintf('MINIMAL DIFFERENTIATORS + REGULATED COUNTER + HILL NN\nSelected row: %d\n',selection)};
for field = {'block1_idx','block2_idx','block3_idx','counter_idx','nn_idx', ...
        'motifIndexBeforeMembraneFilter','motifIndexBeforeCrossTalkFilter'}
    if isfield(row,field{1}), lines{end+1} = sprintf('%s: %d\n',field{1},row.(field{1})); end
end
lines{end+1} = sprintf('Individual regulatory links supported by supplied evidence: %d\n',row.RegulatoryEvidenceComplete);
lines{end+1} = sprintf('Kinetics validated: %d\n',row.KineticsValidated);
lines{end+1} = sprintf('\nBLOCK 0: cAMP TO THE FIXED INPUT\n');
p = find([data.InputAdapterPaths.pathIndex] == row.adapterPathIndex,1);
if isempty(p), error('MinimalPrint:Adapter','Missing adapter path %d.',row.adapterPathIndex); end
path = data.InputAdapterPaths(p);
lines{end+1} = sprintf('Adapter path: %d\n',path.pathIndex);
for field = {'metabolitePath','reactionPath','reactionNamePath'}
    if isfield(path,field{1}), lines{end+1} = sprintf('%s: %s\n',field{1},string(path.(field{1}))); end
end
lines{end+1} = sprintf('\nPROTECTED SIGNAL METABOLITES\n');
for k = 1:numel(schema.mets)
    role = schema.mets{k}; idx = row.(role);
    entry = data.Differentiator_Metabolites(data.Differentiator_Metabolites.Index == idx,:);
    if height(entry) ~= 1, error('MinimalPrint:MetIndex','Missing/ambiguous saved %s index %d.',role,idx); end
    lines{end+1} = sprintf('%s | %d | %s | %s | %s\n',role,idx,entry.ID,entry.Name,entry.Formula);
end
lines{end+1} = sprintf('Physical cascade: X2=Y1; Y2=m_d2=m_d. Sole counter signal: m_d2. Counter output: m_c. NN output: m_out.\n');
for k = 1:numel(schema.rxns)
    if k == 1, lines{end+1} = sprintf('\nBLOCK 1: THREE-REACTION DIFFERENTIATOR\n'); end
    if k == 4, lines{end+1} = sprintf('\nBLOCK 2: THREE-REACTION DIFFERENTIATOR\n'); end
    if k == 7, lines{end+1} = sprintf('\nBLOCK 3: COUNTER AND Cs BUFFERING\n'); end
    if k == 12, lines{end+1} = sprintf('\nBLOCK 4: HILL NN PRODUCTION AND REMOVAL\n'); end
    role = schema.rxns{k}; idx = row.(role); direction = row.(schema.dirs{k});
    entry = data.Differentiator_Reactions(data.Differentiator_Reactions.Index == idx,:);
    if height(entry) ~= 1 || ~ismember(direction,[-1 1])
        error('MinimalPrint:RxnIndex','Invalid saved %s index/direction.',role);
    end
    equation = entry.ForwardEquation;
    if direction == -1
        sides = split(equation,' -> '); equation = sides(2)+" -> "+sides(1);
    end
    lines{end+1} = sprintf('\n%s | %d | %s | %s | direction %+d\nEC: %s\nGPR: %s\n%s\n', ...
        role,idx,entry.ID,entry.Name,direction,entry.EC,entry.GPR,equation);
end
lines{end+1} = sprintf('\nNONCONSUMING REGULATORY LINKS\n');
evidence = data.RegulatoryEvidenceInput;
for k = 1:numel(schema.targets)
    target = schema.targets{k}; met = row.(schema.regulators{k});
    rxn = row.(target); direction = row.(strrep(target,'r_','dir_'));
    if schema.effects(k) == 1, effect = 'activation'; else, effect = 'inhibition'; end
    supported = evidence.MetIndex == met & evidence.RxnIndex == rxn & ...
        evidence.Direction == direction & evidence.Effect == schema.effects(k);
    source = 'PENDING: verify this natural interaction or engineer it';
    if any(supported), source = char(strjoin(evidence.Evidence(supported),'; ')); end
    lines{end+1} = sprintf('%s (%s) -> %s (%s): %s | %s\n', ...
        schema.regulators{k},data.EC_ModelMetIDs(met),target,data.EC_ModelRxnIDs(rxn),effect,source);
end
lines{end+1} = sprintf('\nKinetic requirements still to verify/tune:\n');
lines{end+1} = sprintf('Counter gate: 1-h2(m_d2) on Cp and Cr; no Y1-to-counter regulation.\n');
lines{end+1} = sprintf('Cs supply and basal removal approximate k_s*(Cs_baseline-Cs).\n');
lines{end+1} = sprintf('NN: Hill(C), Hill(Ni), Hill(Ni), Hill(weight*(Na+Nb)); independent removal.\n');
lines{end+1} = sprintf('Precursor concentrations, baseline balance, rate parameters and leakage must be checked.\n');
text = strjoin(lines,'');
end
