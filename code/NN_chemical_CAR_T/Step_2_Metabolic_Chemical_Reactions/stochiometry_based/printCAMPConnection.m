function printCAMPConnection(fid,model,motif,path,rankInfo,p,cfg)
%PRINTCAMPCONNECTION Full connector reactions in the circuit report format.
% fid=1 prints to the MATLAB command window; an open file ID prints to a report.
if isfield(rankInfo,'metadata')
    meta = rankInfo.metadata;
else
    [~,meta] = prepareCircuitNetwork(model,cfg,p);
end
fprintf(fid,'\n------------------------------------------------------------\n');
fprintf(fid,'cAMP INPUT CONNECTION (OUTPUT m_in)\n');
fprintf(fid,'------------------------------------------------------------\n\n');
fprintf(fid,'  Status: %s | conversion steps: %g\n',path.status,path.distance);
if ~isfinite(path.distance)
    fprintf(fid,'  No connecting chain was found under the configured filters/depth.\n');
    return;
end
if isempty(path.reactionIndex)
    m = motif.metIndex(strcmp(p.metNames,'m_in'));
    fprintf(fid,'  m_in = cAMP: %s | %s\n',meta.metIDs{m},meta.metNames{m});
    fprintf(fid,'  No connecting reaction is required.\n');
    return;
end
labels = repmat({''},numel(model.mets),1);
for i = 1:numel(p.metNames), labels{motif.metIndex(i)} = p.metNames{i}; end
labels{path.metIndex(1)} = 'm_cAMP';
for step = 1:numel(path.reactionIndex)-1
    labels{path.metIndex(step+1)} = sprintf('m_CAMP_%d',step);
end
labels{path.metIndex(end)} = 'm_in';
fprintf(fid,'  CONNECTOR SIGNAL METABOLITES\n');
for i = 1:numel(path.metIndex)
    m = path.metIndex(i);
    fprintf(fid,'    %-10s | %-22s | %-42s | %s | %s\n', ...
        labels{m},meta.metIDs{m},meta.metNames{m},meta.formulas{m},meta.compartments{m});
end
fprintf(fid,'\n  Path: %s\n',strjoin(meta.metIDs(path.metIndex),' -> '));
fprintf(fid,'  Shared enzyme among connector reactions: %d\n',path.enzymeConflict);
if isfield(path,'unknownCompartments')
    fprintf(fid,'  Unknown compartment in connector: %d\n',path.unknownCompartments);
end
fprintf(fid,'  These connector candidates are reported separately from the mapped circuit.\n\n');
for step = 1:numel(path.reactionIndex)
    r = path.reactionIndex(step); direction = path.direction(step);
    col = direction*model.S(:,r);
    species = find(abs(col) > cfg.TOLERANCE);
    comps = unique(meta.compartments(species)); comps = comps(~cellfun(@isempty,comps));
    fprintf(fid,'  r_CAMP_%d (cAMP-to-input connector step %d)\n',step,step);
    fprintf(fid,'      ID  : %s\n',meta.rxnIDs{r});
    fprintf(fid,'      Name: %s\n',meta.rxnNames{r});
    fprintf(fid,'      GPR : %s\n',meta.gpr{r});
    fprintf(fid,'      EC  : %s\n',meta.ec{r});
    fprintf(fid,'      Enzyme conflict keys: %s\n',strjoin(meta.enzymeTokens{r},'; '));
    fprintf(fid,'      Compartment(s): %s\n',strjoin(comps,', '));
    fprintf(fid,'      Direction: %+d relative to model.S\n',direction);
    fprintf(fid,'      Original bounds: [%g, %g]\n',model.lb(r),model.ub(r));
    fprintf(fid,'      Eqn   : %s\n',equation(col,meta.formulas,meta.metIDs,labels,cfg.TOLERANCE));
    fprintf(fid,'      Names : %s\n',equation(col,meta.metNames,meta.metIDs,labels,cfg.TOLERANCE));
    fprintf(fid,'      IDs   : %s\n\n',equation(col,meta.metIDs,meta.metIDs,labels,cfg.TOLERANCE));
end
end

function text = equation(col,names,fallback,labels,tol)
text = [side(find(col < -tol),col,names,fallback,labels),'  ->  ', ...
    side(find(col > tol),col,names,fallback,labels)];
end
function text = side(indices,col,names,fallback,labels)
parts = cell(numel(indices),1);
for k = 1:numel(indices)
    i = indices(k); name = names{i};
    if isempty(name), name = fallback{i}; end
    if ~isempty(labels{i}), name = [name,'(',labels{i},')']; end
    coefficient = abs(full(col(i)));
    if abs(coefficient-1) > 1e-12, name = [sprintf('%.12g ',coefficient),name]; end
    parts{k} = name;
end
if isempty(parts), text = '0'; else, text = strjoin(parts,' + '); end
end
