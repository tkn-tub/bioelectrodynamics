function text = renderFiveReactionDifferentiators(data,selection)
% Print the current Block 2 motif using only its compact checkpoint.
if ~isfield(data,'TopologyVersion') || ...
        ~isequal(string(data.TopologyVersion),"five_reaction_differentiators_v2")
    error('PolarityPrint:Topology','Load the current five-plus-five-reaction Block 2 checkpoint.');
end
motifs = data.Motifs_Diff_Block_2;
validateattributes(selection,{'numeric'},{'scalar','integer','positive','<=',height(motifs)});
row = table2struct(motifs(selection,:));
parts = {sprintf('TEN-REACTION DIFFERENTIATOR CASCADE\nTopology: %s\nMotif %d | Block 1 row %d\n', ...
    data.TopologyVersion,selection,row.block1_idx)};
p = find([data.InputAdapterPaths.pathIndex] == row.adapterPathIndex,1);
if isempty(p), error('PolarityPrint:Adapter','The selected adapter path is unavailable.'); end
path = data.InputAdapterPaths(p);
parts{end+1} = sprintf('\nBLOCK 0: cAMP TO THE FIXED INPUT\nAdapter path index: %d\n',path.pathIndex);
for field = {'metabolitePath','reactionPath','reactionNamePath'}
    if isfield(path,field{1}), parts{end+1} = sprintf('%s: %s\n',field{1},string(path.(field{1}))); end
end
parts{end+1} = sprintf('\nSELECTED METABOLITES\nm_d1 = m_d = Y1; X2 = Y1; m_d2 = Y2 is the counter signal.\n');
for field = {'m_in','m_M1','m_Y1','m_F1_delay','m_M2','m_Y2'}
    index = row.(field{1});
    entry = data.Differentiator_Metabolites(data.Differentiator_Metabolites.Index == index,:);
    if height(entry) ~= 1, error('PolarityPrint:Dictionary','Missing/ambiguous metabolite annotation.'); end
    parts{end+1} = sprintf('%s | %d | %s | %s | %s\n',field{1},index,entry.ID,entry.Name,entry.Formula);
end
roles = {'r_S1','r_S2','r_F','r_I','r_R','r_S1_2','r_S2_2','r_F2','r_I2','r_R2'};
labels = {'X to M1','M1 to Y1','X to F1','joint F1 and Y1 removal','independent F1 removal', ...
    'Y1 to M2','M2 to Y2','Y1 directly to Y2','joint M2 and Y2 removal','independent Y2 removal'};
for k = 1:numel(roles)
    if k == 1, parts{end+1} = sprintf('\nBLOCK 1: FIVE-REACTION SLOW/FAST TOPOLOGY\n'); end
    if k == 6, parts{end+1} = sprintf('\nBLOCK 2: FIVE-REACTION POSITIVE-POLARITY TOPOLOGY\n'); end
    role = roles{k}; index = row.(role); direction = row.(strrep(role,'r_','dir_'));
    entry = data.Differentiator_Reactions(data.Differentiator_Reactions.Index == index,:);
    if height(entry) ~= 1 || ~ismember(direction,[-1 1]), error('PolarityPrint:Dictionary','Invalid reaction annotation/direction.'); end
    equation = entry.ForwardEquation;
    if direction == -1
        sides = split(equation,' -> '); equation = sides(2)+" -> "+sides(1);
    end
    parts{end+1} = sprintf('\n%s (%s) | %d | %s | %s | direction %+d\nEC: %s\nGPR: %s\n%s\n', ...
        role,labels{k},index,entry.ID,entry.Name,direction,entry.EC,entry.GPR,equation);
end
parts{end+1} = sprintf('\nKinetics validated: 0. Reaction rates and saturation parameters still need tuning.\n');
text = strjoin(parts,'');
end
