function text = renderMinimalDifferentiator(data,selection)
% Render one new Block 2 checkpoint row without a model or toolbox.
motifs = data.Motifs_Diff_Block_2;
validateattributes(selection,{'numeric'},{'scalar','integer','positive','<=',height(motifs)});
row = motifs(selection,:);
parts = {sprintf('SIX-REACTION DIFFERENTIATOR CASCADE\nMotif %d | Block 1 row %d\n',selection,row.block1_idx)};
p = find([data.InputAdapterPaths.pathIndex] == row.adapterPathIndex,1);
if isempty(p), error('MinimalPrint:Adapter','Saved adapter path is missing.'); end
path = data.InputAdapterPaths(p);
parts{end+1} = sprintf('\nBLOCK 0: cAMP TO THE FIXED INPUT\nAdapter path index: %d\n',path.pathIndex);
if isfield(path,'metabolitePath'), parts{end+1} = sprintf('%s\n',path.metabolitePath); end
if isfield(path,'reactionPath'), parts{end+1} = sprintf('Reactions: %s\n',path.reactionPath); end
if isfield(path,'reactionNamePath'), parts{end+1} = sprintf('Names: %s\n',path.reactionNamePath); end
parts{end+1} = sprintf('\nSELECTED DIFFERENTIATOR MOTIF\nX2=Y1; final physical output Y2=m_d\n');
for name = {'m_in','m_M1','m_Y1','m_M2','m_Y2'}
    idx = row.(name{1}); entry = data.Differentiator_Metabolites(data.Differentiator_Metabolites.Index == idx,:);
    parts{end+1} = sprintf('%s | %d | %s | %s | %s\n',name{1},idx,entry.ID,entry.Name,entry.Formula);
end
for name = {'r_11','r_12','r_13','r_21','r_22','r_23'}
    idx = row.(name{1}); entry = data.Differentiator_Reactions(data.Differentiator_Reactions.Index == idx,:);
    direction = row.(strrep(name{1},'r_','dir_'));
    equation = entry.ForwardEquation;
    if direction == -1
        sides = split(equation,' -> '); equation = sides(2)+" -> "+sides(1);
    end
    parts{end+1} = sprintf('\n%s | %d | %s | %s | direction %+d\nEC: %s\nGPR: %s\n%s\n', ...
        name{1},idx,entry.ID,entry.Name,direction,entry.EC,entry.GPR,equation);
end
text = strjoin(parts,'');
end
