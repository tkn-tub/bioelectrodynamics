function writeMetabolicCircuitReport(fileName,model,motifs,info,p,cfg,ranking,paths,rankInfo)
%WRITEMETABOLICCIRCUITREPORT Use the labels and report layout in the attachment.
fid = fopen(fileName,'w','n','UTF-8');
if fid < 0, error('Circuit:Report','Cannot open report: %s',fileName); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
meta = info.metadata;
fprintf(fid,'############################################################\n');
fprintf(fid,'SELECTED CHEMICAL CIRCUIT MOTIFS\n');
fprintf(fid,'Topology: %s\n',p.name);
fprintf(fid,'Model: %s\n',cfg.MODEL_FILE);
fprintf(fid,'Available motifs: %d\n',numel(motifs));
if isfield(info,'oneMotifPerInput')
    fprintf(fid,'One motif per distinct input: %d\n',info.oneMotifPerInput);
    fprintf(fid,'Distinct inputs returned: %d\n',info.distinctInputsReturned);
    fprintf(fid,'Search coverage: %s\n',info.exhaustiveScope);
end
fprintf(fid,'Search status: %s\n',info.status);
fprintf(fid,'Search exhaustive: %d\n',info.exhaustive);
fprintf(fid,'Elapsed search seconds: %.3f\n',info.elapsedSeconds);
fprintf(fid,'Visited search states: %d\n',info.nodes);
fprintf(fid,'Parallel used: %d\n',info.parallelUsed);
fprintf(fid,'Forward only: %d\n',cfg.FORWARD_ONLY);
fprintf(fid,'Enzyme evidence: %s\n',cfg.ENZYME_EVIDENCE);
if isempty(cfg.ENZYME_ID_FIELD)
    fprintf(fid,'Enzyme uniqueness: no shared complete EC number OR GPR gene.\n');
    fprintf(fid,'EC grouping is a conservative activity-class proxy for enzyme identity.\n');
else
    fprintf(fid,'Enzyme uniqueness: no shared ID from %s OR GPR gene.\n',cfg.ENZYME_ID_FIELD);
end
fprintf(fid,'Metabolite identity: %s\n',cfg.METABOLITE_IDENTITY);
fprintf(fid,'Metabolite uniqueness scope: %s\n',cfg.METABOLITE_SCOPE);
fprintf(fid,'Auxiliary reuse exceptions: %s\n',joinSetting(cfg.AUXILIARY_REUSE_ALLOWED));
fprintf(fid,'Excluded role metabolites: %s\n',joinSetting(cfg.EXCLUDED_ROLE_METABOLITES));
fprintf(fid,'Compartment filter: %s\n',cfg.COMPARTMENT);
fprintf(fid,'Eligible enzymatic conversions: %d / %d\n', ...
    meta.eligibleReactions,meta.originalReactions);
fprintf(fid,'Allowed oriented reactions: %d\n',meta.orientedReactions);
fprintf(fid,'Complete EC annotations: %d\n',meta.annotatedECReactions);
fprintf(fid,'GPR annotations: %d\n',meta.annotatedGPRReactions);
fprintf(fid,'############################################################\n\n');
if isempty(motifs)
    fprintf(fid,'NO COMPLETE MOTIF RETURNED\n');
    if info.exhaustive
        fprintf(fid,'No assignment satisfies this exact topology and these annotation-based restrictions.\n');
    else
        fprintf(fid,'The search stopped at a limit. This does not establish that no match exists.\n');
    end
    fprintf(fid,'Initial candidate counts (before global propagation):\n');
    for j = 1:numel(p.rxnNames)
        fprintf(fid,'  %-9s : %d\n',p.rxnNames{j},info.initialReactionCandidates(j));
    end
    fprintf(fid,'\nAn enzymatic removal must convert its rail into other species.\n');
    fprintf(fid,'Unannotated reactions and one-sided exchange/demand reactions are excluded.\n');
end
for k = 1:numel(motifs)
    motif = motifs(k);
    fprintf(fid,'============================================================\nMOTIF %d\n',k);
    fprintf(fid,'============================================================\n\n');
    fprintf(fid,'SELECTED MOTIF: METABOLITES\n');
    for i = 1:numel(p.metNames)
        m = motif.metIndex(i);
        fprintf(fid,'  %-9s | %-22s | %-42s | %s | %s\n', ...
            p.metNames{i},meta.metIDs{m},meta.metNames{m},meta.formulas{m},meta.compartments{m});
    end
    labels = repmat({''},numel(model.mets),1);
    for i = 1:numel(p.metNames), labels{motif.metIndex(i)} = p.metNames{i}; end
    if nargin >= 9
        rankRow = find(ranking.MotifIndex == k,1);
        if ~isempty(rankRow)
            fprintf(fid,'\nInput connection rank: %d\n',ranking.Rank(rankRow));
            printCAMPConnection(fid,model,motif,paths(rankRow),rankInfo,p,cfg);
        end
    end
    for block = 1:numel(p.blocks)
        fprintf(fid,'\n------------------------------------------------------------\n');
        fprintf(fid,'BLOCK %d: %s\n',block,p.blockNames{block});
        fprintf(fid,'------------------------------------------------------------\n\n');
        for j = p.blocks{block}
            r = motif.reactionIndex(j);
            col = motif.direction(j)*model.S(:,r);
            species = find(abs(col) > cfg.TOLERANCE);
            comps = unique(meta.compartments(species)); comps = comps(~cellfun(@isempty,comps));
            fprintf(fid,'  %s (%s)\n',p.rxnNames{j},p.descriptions{j});
            fprintf(fid,'      ID  : %s\n',meta.rxnIDs{r});
            fprintf(fid,'      Name: %s\n',meta.rxnNames{r});
            fprintf(fid,'      GPR : %s\n',meta.gpr{r});
            fprintf(fid,'      EC  : %s\n',meta.ec{r});
            fprintf(fid,'      Enzyme conflict keys: %s\n',strjoin(meta.enzymeTokens{r},'; '));
            fprintf(fid,'      Compartment(s): %s\n',strjoin(comps,', '));
            fprintf(fid,'      Direction: %+d relative to model.S\n',motif.direction(j));
            fprintf(fid,'      Original bounds: [%g, %g]\n',model.lb(r),model.ub(r));
            fprintf(fid,'      Eqn   : %s\n',equation(col,meta.formulas,meta.metIDs,labels,cfg.TOLERANCE));
            fprintf(fid,'      Names : %s\n',equation(col,meta.metNames,meta.metIDs,labels,cfg.TOLERANCE));
            fprintf(fid,'      IDs   : %s\n\n',equation(col,meta.metIDs,meta.metIDs,labels,cfg.TOLERANCE));
        end
    end
    fprintf(fid,'AUXILIARY METABOLITE AUDIT\n');
    shared = false;
    keys = unique(meta.metKeys);
    roleKeys = meta.metKeys(motif.metIndex);
    for a = 1:numel(keys)
        if any(strcmp(roleKeys,keys{a})), continue; end
        rows = find(strcmp(meta.metKeys,keys{a}));
        occurrences = find(full(any(abs(model.S(rows,motif.reactionIndex)) > cfg.TOLERANCE,1)));
        if numel(occurrences) > 1
            shared = true;
            fprintf(fid,'  %s | shared by: %s\n',keys{a},strjoin(p.rxnNames(occurrences),', '));
        end
    end
    if ~shared, fprintf(fid,'  No shared auxiliary compound.\n'); end
    fprintf(fid,'\nASSUMPTIONS AND INTERFACES\n');
    fprintf(fid,'Every diagram connection is direct; no added bridge reactions or input adapter.\n');
    fprintf(fid,'m_in, m_D1, m_D2 and m_c are shared only at the explicitly drawn interfaces.\n');
    fprintf(fid,'m_out is the neuron signal output, produced by r_Nw and r_Nb.\n');
    fprintf(fid,'No selected role participates in an unrequested selected reaction.\n');
    fprintf(fid,'Removal arrows to empty set mean conversion to unlabelled products, not destruction of matter.\n');
    if strcmpi(cfg.NEURON_INPUT_MODE,'substrate')
        fprintf(fid,'consuming_N = true: r_Nw consumes m_c and produces m_out.\n');
    else
        fprintf(fid,'consuming_N = false: m_c activates r_Nw using supplied evidence; not inferred from S.\n');
        fprintf(fid,'r_Nw consumes other substrates, produces m_out, and has zero net stoichiometry for m_c.\n');
        pairs = cfg.ACTIVATION_PAIRS;
        for a = 1:size(pairs,1)
            if strcmp(char(pairs{a,1}),meta.rxnIDs{motif.reactionIndex(16)}) && ...
                    strcmp(char(pairs{a,2}),meta.metIDs{motif.metIndex(10)})
                fprintf(fid,'Supplied activation pair: %s <- %s',char(pairs{a,1}),char(pairs{a,2}));
                if size(pairs,2) >= 3, fprintf(fid,' | %s',char(pairs{a,3})); end
                fprintf(fid,'\n');
            end
        end
    end
    fprintf(fid,'m_I and m_b are externally supplied/buffered assumptions; replenishment is not mapped.\n');
    fprintf(fid,'Fast/slow names identify topology roles; kinetic separation has not been established.\n');
    fprintf(fid,'Connectivity validated: 1. Kinetics validated: 0. Flux feasibility validated: 0.\n');
    fprintf(fid,'The net stoichiometric model does not establish enzyme regulation or the required circuit dynamics.\n\n');
end
end

function text = equation(col,names,fallback,labels,tol)
left = find(col < -tol); right = find(col > tol);
text = [side(left,col,names,fallback,labels),'  ->  ',side(right,col,names,fallback,labels)];
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
function text = joinSetting(values)
if isempty(values), text = '(none)';
else
    if ischar(values) || isstring(values), values = cellstr(values); end
    text = strjoin(values,', ');
end
end
