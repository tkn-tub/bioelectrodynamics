function writeCAMPInputRankingReport(fileName,model,motifs,ranking,paths,rankInfo,p,cfg)
%WRITECAMPINPUTRANKINGREPORT Explain distances and list the connecting reactions.
fid = fopen(fileName,'w','n','UTF-8');
if fid < 0, error('Circuit:CAMPReport','Cannot open report: %s',fileName); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
mets = cellstr(string(model.mets(:)));
fprintf(fid,'cAMP TO CIRCUIT INPUT RANKING\n');
fprintf(fid,'Model: %s\n',cfg.MODEL_FILE);
fprintf(fid,'Returned circuit candidates: %d\n',numel(motifs));
fprintf(fid,'Ranking status: %s\n',rankInfo.status);
fprintf(fid,'Maximum conversion steps: %g\n',cfg.CAMP_MAX_DEPTH);
fprintf(fid,'cAMP source compartment: %s\n',cfg.CAMP_SOURCE_COMPARTMENT);
fprintf(fid,'Transport permitted: %d\n',cfg.CAMP_ALLOW_TRANSPORT);
fprintf(fid,'\nDistance 0: m_in is an identified cAMP model species.\n');
fprintf(fid,'Distance N: a directed substrate-to-product path with N enzymatic reactions.\n');
fprintf(fid,'Inf: not reached under these filters/depth; this does not prove biological disconnection.\n');
fprintf(fid,'Ranks are among supplied motifs, not every possible circuit in the model.\n');
fprintf(fid,'Paths use bound-permitted orientations and the configured enzyme-evidence gate.\n');
fprintf(fid,'Common cofactors/hubs are excluded as path nodes; AMP and adenosine remain available.\n');
fprintf(fid,'Path reactions cannot reuse circuit reactions/enzymes or touch another circuit signal.\n');
fprintf(fid,'Compartments retain separate model species; no implicit compartment jumps.\n');
fprintf(fid,'A flagged shared enzyme within a path still needs review before adding that path to a circuit.\n');
fprintf(fid,'Paths are structural candidates: kinetics, steady-state flux and atom transfer are not established.\n');
fprintf(fid,'Receptor regulation is not inferred from stoichiometry. No adapter is added to the circuit.\n');
fprintf(fid,'\nBIOLOGICAL INPUT RATIONALE\n');
fprintf(fid,'Clayton et al. (2011), DOI 10.4049/jimmunol.1003884:\n');
fprintf(fid,'Cancer-exosome CD39/CD73 generate extracellular adenosine from ATP/AMP.\n');
fprintf(fid,'Adenosine acts through the A2A receptor; adenylyl cyclase generates intracellular cAMP.\n');
fprintf(fid,'This supports a cAMP anchor, but does not establish vesicle-count specificity or preservation of peak timing.\n');
fprintf(fid,'\nIDENTIFIED cAMP SOURCES\n');
if isfield(rankInfo,'sourceMetaboliteIDs') && ~isempty(rankInfo.sourceMetaboliteIDs)
    fprintf(fid,'  %s\n',strjoin(rankInfo.sourceMetaboliteIDs,', '));
else
    fprintf(fid,'  None found in the requested compartment. Set CAMP_SOURCE_METABOLITES to exact model IDs.\n');
    if isfield(rankInfo,'sourceCandidateMetaboliteIDs') && ~isempty(rankInfo.sourceCandidateMetaboliteIDs)
        fprintf(fid,'  Candidates before compartment filtering: %s\n', ...
            strjoin(rankInfo.sourceCandidateMetaboliteIDs,', '));
    end
end
if isfield(rankInfo,'excludedPathMetaboliteIDs')
    fprintf(fid,'\nEXCLUDED PATH NODES\n  %s\n',strjoin(rankInfo.excludedPathMetaboliteIDs,', '));
end
fprintf(fid,'\nRANKED CANDIDATES\n');
for k = 1:height(ranking)
    row = ranking(k,:); motif = motifs(row.MotifIndex); path = paths(k);
    fprintf(fid,'\n------------------------------------------------------------\n');
    fprintf(fid,'Rank %d | original motif %d | distance %g | %s\n', ...
        row.Rank,row.MotifIndex,row.Distance,row.Status{1});
    fprintf(fid,'Input: %s (%s)\n',row.InputMetabolite{1},row.InputName{1});
    fprintf(fid,'Signal roles (current labels):\n');
    for i = 1:numel(p.metNames)
        fprintf(fid,'  %-10s = %s\n',p.metNames{i},mets{motif.metIndex(i)});
    end
    printCAMPConnection(fid,model,motif,path,rankInfo,p,cfg);
end
end
