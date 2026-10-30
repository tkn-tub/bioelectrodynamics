%% A_Master.m - Strict +1 transport, differentiator, counter and two neurons
% ONE FILE: discover exchange transport -> transporter GPR -> D1 candidates.
% Complete motifs use stoichiometric exchange only; protein-path regulation is excluded.
% Requires COBRA Toolbox readCbModel, MATLAB with webread (online enrichment optional).
% Run: matlab -batch "A_Master"
% No checkpoints or network-cache directories are generated.
clearvars; clc;
MODEL_NAME = 'Puniya2021_TNM1055.xml';
OUT_NAME = 'Transport_Stages1to5_BRENDA_Results_RankedTop10';

%% Ranking: distance means the number of cAMP -> m_in reactions.
PRINT_TOP_MOTIFS = 10;
% Exclude complete motifs containing reactions annotated with these EC numbers.
% EC annotations are read from the model reaction fields (rxnECNumbers/rxnEC).
EXCLUDED_MOTIF_ECS = ["2.4.1.23"];
BLOCK0_SOURCE_MET_ID = 'camp[c]';
BLOCK0_MAX_DEPTH = Inf;
BLOCK0_MAX_EXPANSIONS = Inf;
BLOCK0_ALLOW_ZERO_STEP = false;
BLOCK0_DIAGNOSTICS_ONLY = false;
RUN_RANKING_SELF_TESTS = true;
% All regulator/exchange targets are considered. Reaction branching and trial
% caps still bound the combinatorial search; use Inf to remove an individual cap.
D1_MAX_TARGETS = Inf;
D1_MAX_MOTIFS = Inf;
D1_MAX_PER_TARGET = 100;       % Closest retained D1 motifs per target.
D1_MAX_BRANCH = 14;
D1_MAX_RXN_CAND = 60;
D1_MAX_TESTS_PER_TARGET = 40000;
EXCHANGE_D1_MAX_TARGETS = Inf;       % Exchange reaction/counter-substrate pairs to inspect.
EXCHANGE_D1_MAX_CONNECTIONS = Inf;
EXCHANGE_METABOLIC_MAX_DEPTH = 1;          % Prefer m_Sd = counter-substrate (0) or one forward step away (1).
EXCHANGE_MAX_SIGNAL_PRECURSORS = Inf;      % Per exchange pair; set finite only to trade coverage for speed.
EXCHANGE_REQUIRE_GPR = true;               % Keep only transport reactions with a model gene rule.
% Complete-circuit search admits only stoichiometric exchange transport.
% Regulatory protein-path candidates are not passed into the final motif search.
EXCHANGE_ONLY = true;
N8_MAX_BRANCH = 25;
N8_MAX_TESTS = 250000;
% Discovery continues after ten matches to count distinct solutions within
% these limits. Only PRINT_TOP_MOTIFS records are retained for detailed printing.
% Equal-distance motifs may share m_in. The TXT reports found/printed separately.

STAGE5_MAX_HOPS = 2; % reaction steps backward; structural only
STAGE5_MAX_MATCH_ROWS = 250000; % output safety limit
STAGE5_MAX_PATH_ROWS = 100000;
MAX_HOPS = 3;
ONLINE = true;             % false for offline run (local TSV and EC overrides only)
% Strict policy: all searches use +1 relative to the original model.S.
% Backward traversal in Stage 5 traces precursors; it does not reverse chemistry.
TIMEOUT = 30;
MAX_RETRIES = 4;
RETRY_PAUSE = 2; % seconds, exponential backoff
EC_BATCH_SIZE = 15;
% No EC protein cap: every identified regulatory protein is considered.


assert(isscalar(PRINT_TOP_MOTIFS) && isfinite(PRINT_TOP_MOTIFS) && ...
    PRINT_TOP_MOTIFS>=1 && PRINT_TOP_MOTIFS==floor(PRINT_TOP_MOTIFS), ...
    'Ranking:OutputCount','PRINT_TOP_MOTIFS must be a positive integer.');
for capName={'D1_MAX_TARGETS','D1_MAX_MOTIFS','D1_MAX_PER_TARGET', ...
        'D1_MAX_BRANCH','D1_MAX_RXN_CAND','D1_MAX_TESTS_PER_TARGET', ...
        'EXCHANGE_D1_MAX_TARGETS','EXCHANGE_D1_MAX_CONNECTIONS','N8_MAX_BRANCH','N8_MAX_TESTS'}
    capValue=eval(capName{1});
    assert(isnumeric(capValue) && isscalar(capValue) && ~isnan(capValue) && ...
        capValue>=1 && (isinf(capValue) || capValue==floor(capValue)), ...
        'Ranking:SearchCap','%s must be a positive integer or Inf.',capName{1});
end
assert(isscalar(EXCHANGE_METABOLIC_MAX_DEPTH) && isfinite(EXCHANGE_METABOLIC_MAX_DEPTH) && ...
    EXCHANGE_METABOLIC_MAX_DEPTH>=0 && EXCHANGE_METABOLIC_MAX_DEPTH==floor(EXCHANGE_METABOLIC_MAX_DEPTH), ...
    'ExchangeSearch:DepthCap','EXCHANGE_METABOLIC_MAX_DEPTH must be a nonnegative integer.');
assert(isscalar(EXCHANGE_MAX_SIGNAL_PRECURSORS) && ~isnan(EXCHANGE_MAX_SIGNAL_PRECURSORS) && ...
    EXCHANGE_MAX_SIGNAL_PRECURSORS>=1 && ...
    (isinf(EXCHANGE_MAX_SIGNAL_PRECURSORS) || EXCHANGE_MAX_SIGNAL_PRECURSORS==floor(EXCHANGE_MAX_SIGNAL_PRECURSORS)), ...
    'ExchangeSearch:PrecursorCap','EXCHANGE_MAX_SIGNAL_PRECURSORS must be a positive integer or Inf.');
if RUN_RANKING_SELF_TESTS
    selfTestBlock0Connectivity();
    selfTestDistanceRanking();
    selfTestMotifCounts();
end

% Inline equivalent of addingPathParentFolderByName('code').
base=pwd; codeRoot='';
while true
    [up,nm]=fileparts(base);
    if strcmpi(nm,'code'), codeRoot=base; break; end
    if isempty(up) || strcmp(up,base), break; end
    base=up;
end
if ~isempty(codeRoot), addpath(genpath(codeRoot)); end
if isempty(codeRoot), codeRoot=pwd; end
outDir=fullfile(pwd,OUT_NAME);
if ~exist(outDir,'dir'), mkdir(outDir); end
logPath=fullfile(outDir,'stage1to4_report.txt');
logID=fopen(logPath,'w'); assert(logID>0,'Cannot write report');
logGuard=onCleanup(@() fclose(logID));
report(logID,'Cell-network search: import transport, stoichiometric exchange, and regulatory candidates.');
report(logID,['Started: ' datestr(now,31)]);
report(logID,['MATLAB version: ' version]);
report(logID,'Direction policy: STRICT +1 relative to original model.S in every stage.');

% Locate model, including the known Puniya folder structure.
modelPath=which(MODEL_NAME);
if isempty(modelPath)
    p=fullfile(codeRoot,'database','Car_T_Cell_models',MODEL_NAME);
    if isfile(p),modelPath=p;end
end
if isempty(modelPath)
    modelPath=findNamed(codeRoot,MODEL_NAME);
end
assert(~isempty(modelPath),'Model not found: %s. Put it under code/database/Car_T_Cell_models.',MODEL_NAME);
assert(exist('readCbModel','file')~=0,'COBRA Toolbox readCbModel is not on MATLAB path');
model=readCbModel(modelPath);
report(logID,['Model: ' modelPath]);
report(logID,sprintf('Model reactions: %d; metabolites: %d',numel(model.rxns),numel(model.mets)));

%% STAGE 1: find extracellular -> cytosolic transport in stoichiometry.
metIDs=string(model.mets(:)); rxnIDs=string(model.rxns(:));
if isfield(model,'rxnNames'),rxnNames=string(model.rxnNames(:));else,rxnNames=rxnIDs;end
if isfield(model,'lb'), lb=model.lb(:);else,lb=zeros(numel(rxnIDs),1);end
if isfield(model,'ub'), ub=model.ub(:);else,ub=ones(numel(rxnIDs),1)*1000;end
% Only standard [e] and [c] or _e and _c compartment naming is used.
[baseMet,comp]=splitCompartments(metIDs);
% A single forward-only breadth-first search supplies every later motif.
block0Cache=buildBlock0Cache(model,metIDs,comp,ub, ...
    BLOCK0_SOURCE_MET_ID,BLOCK0_MAX_DEPTH,BLOCK0_MAX_EXPANSIONS);
block0InputAllowed=block0Cache.reachable;
if ~BLOCK0_ALLOW_ZERO_STEP
    block0InputAllowed(block0Cache.start)=false;
end
report(logID,sprintf('Block 0: %d eligible m_in species reachable from %s; depth <= %g; expansions = %d; expansion cap reached = %d.', ...
    nnz(block0InputAllowed),BLOCK0_SOURCE_MET_ID,BLOCK0_MAX_DEPTH, ...
    block0Cache.expansions,block0Cache.expansionLimitReached));
if ~any(block0InputAllowed)
    report(logID,'Block 0: no eligible connected inputs; complete-circuit search will return no motifs.');
end
if isfield(model,'metNames') && numel(model.metNames)==numel(model.mets)
    metNames=string(model.metNames(:));
else
    metNames=string(model.mets(:));
end
block0Rows=find(comp=="c");
block0Status=repmat("no_path_within_depth",numel(block0Rows),1);
if isinf(BLOCK0_MAX_DEPTH)
    block0Status(:)="no_path_in_forward_cytosolic_graph";
end
if block0Cache.expansionLimitReached
    block0Status(:)="not_reached_before_expansion_cap";
end
block0Status(block0Cache.reachable(block0Rows))="connected";
if ~BLOCK0_ALLOW_ZERO_STEP
    block0Status(block0Rows==block0Cache.start)="zero_step_not_allowed";
end
block0MetPaths=strings(numel(block0Rows),1);
block0ReactionPaths=strings(numel(block0Rows),1);
for b0row=1:numel(block0Rows)
    b0target=block0Rows(b0row);
    if block0Cache.reachable(b0target)
        [b0mets,b0rx]=block0CachedPath(block0Cache,b0target);
        block0MetPaths(b0row)=join(metIDs(b0mets),' -> ');
        if ~isempty(b0rx)
            block0ReactionPaths(b0row)=join(rxnIDs(b0rx),' -> ');
        end
    end
end
Block0Reachability=table(metIDs(block0Rows),metNames(block0Rows), ...
    block0Cache.reachable(block0Rows),block0InputAllowed(block0Rows), ...
    block0Cache.depth(block0Rows),block0Status,block0MetPaths,block0ReactionPaths, ...
    'VariableNames',{'MetaboliteID','MetaboliteName','Reachable', ...
    'AllowedAsInput','ReactionDepth','Status','PathMetaboliteIDs','PathReactionIDs'});
writetable(Block0Reachability,fullfile(outDir,'00_Block0_Reachable_Inputs.csv'));
report(logID,['Block 0 reachability report: ' fullfile(outDir,'00_Block0_Reachable_Inputs.csv')]);
if BLOCK0_DIAGNOSTICS_ONLY
    report(logID,'Block 0 diagnostics complete; no online enrichment or circuit search was run.');
    clear logGuard;
    return;
end

external=find(comp=="e"); cytosol=find(comp=="c");
rows=cell(0,9);
for k=1:numel(external)
    e=external(k);
    cs=cytosol(baseMet(cytosol)==baseMet(e));
    if isempty(cs),continue;end
    for cc=reshape(cs,1,[])
        % Forward: extracellular consumed and cytosolic produced.
        se=full(model.S(e,:));
        sc=full(model.S(cc,:));
        jj=find(se<0 & sc>0 & ub'>0); %#ok<FNDSB>
        for j=reshape(jj,1,[])
            rows(end+1,:)={char(rxnIDs(j)),char(rxnNames(j)),char(baseMet(e)), ...
                char(metIDs(e)),char(metIDs(cc)),'forward',lb(j),ub(j),j}; %#ok<SAGROW>
        end
    end
end
if isempty(rows)
    T=table('Size',[0 9],'VariableTypes', ...
        {'string','string','string','string','string','string','double','double','double'}, ...
        'VariableNames',{'reaction','reactionName','cargo','extracellular', ...
        'cytosolic','direction','lowerBound','upperBound','modelIndex'});
else
    T=cell2table(rows,'VariableNames',{'reaction','reactionName','cargo','extracellular', ...
        'cytosolic','direction','lowerBound','upperBound','modelIndex'});
end
% Remove repeated rows, preserve one per cargo/direction/reaction.
if ~isempty(T), T=unique(T,'rows','stable');end
assert(all(string(T.direction)=="forward"), ...
    'ForwardOnly:TransportPool','Transporter pool contains a non-forward direction.');
writetable(T,fullfile(outDir,'01_transport_reactions.csv'));
report(logID,sprintf('Stage 1: %d transport rows; %d unique reactions',height(T),numel(unique(string(T.reaction)))));

%% STAGE 1B: identify stoichiometric extracellular/cytosolic exchange candidates.
% An exchange candidate imports one extracellular metabolite while exporting
% a different cytosolic metabolite in the SAME forward reaction. The latter
% can serve as m_Sd_1 if a D1 motif produces it. This is a stoichiometric
% counter-substrate mechanism; it does not add a transcriptional regulation step.
exchangeRows=cell(0,12);
for ti=1:height(T)
    j=double(T.modelIndex(ti));
    col=full(model.S(:,j));
    eCargo=find(metIDs==string(T.extracellular(ti)),1,'first');
    cCargo=find(metIDs==string(T.cytosolic(ti)),1,'first');
    if isempty(eCargo) || isempty(cCargo),continue;end
    importE=-col(eCargo);importC=col(cCargo);
    if importE<=0 || importC<=0 || ...
            abs(importE-importC)>1e-9*max([1,importE,importC])
        continue;
    end
    cargoBase=baseMet(eCargo);
    counterC=find(comp=="c" & col<0);
    counterE=find(comp=="e" & col>0);
    for ci=reshape(counterC,1,[])
        if baseMet(ci)==cargoBase,continue;end
        matchingE=counterE(baseMet(counterE)==baseMet(ci));
        for ei=reshape(matchingE,1,[])
            counterCStoich=-col(ci);counterEStoich=col(ei);
            if counterCStoich<=0 || counterEStoich<=0 || ...
                    abs(counterCStoich-counterEStoich)> ...
                    1e-9*max([1,counterCStoich,counterEStoich])
                continue;
            end
            exchangeRows(end+1,:)={char(rxnIDs(j)),char(rxnNames(j)), ...
                char(T.cargo(ti)),char(T.extracellular(ti)),char(T.cytosolic(ti)), ...
                char(metIDs(ci)),char(metIDs(ei)),'forward',importC, ...
                counterCStoich,ti,j}; %#ok<SAGROW>
        end
    end
end
if isempty(exchangeRows)
    X=table('Size',[0 12],'VariableTypes', ...
        {'string','string','string','string','string','string','string','string', ...
         'double','double','double','double'}, ...
        'VariableNames',{'reaction','reactionName','importedCargo','extracellularCargo', ...
        'cytosolicCargo','counterSubstrate','counterProduct','direction','importStoich', ...
        'counterStoich','transportRow','modelIndex'});
else
    X=cell2table(exchangeRows,'VariableNames',{'reaction','reactionName','importedCargo', ...
        'extracellularCargo','cytosolicCargo','counterSubstrate','counterProduct', ...
        'direction','importStoich','counterStoich','transportRow','modelIndex'});
    X=unique(X,'rows','stable');
end
assert(all(string(X.direction)=="forward"), ...
    'ForwardOnly:ExchangePool','Exchange pool contains a non-forward direction.');
writetable(X,fullfile(outDir,'01b_exchange_transport_candidates.csv'));
report(logID,sprintf('Stage 1B: %d stoichiometric exchange rows; %d distinct reactions', ...
    height(X),numel(unique(string(X.reaction)))));

%% STAGE 2: recover GPR from grRules or rules -> model.genes.
if isfield(model,'genes'),genes=string(model.genes(:));else,genes=strings(0,1);end
if isfield(model,'grRules'),gr=string(model.grRules(:));else,gr=strings(numel(rxnIDs),1);end
if isfield(model,'rules'),rr=string(model.rules(:));else,rr=strings(numel(rxnIDs),1);end
fullGPR=strings(height(T),1); geneList=strings(height(T),1); gprStatus=strings(height(T),1);
for i=1:height(T)
    j=T.modelIndex(i); [gpr,glist,status]=decodeGPR(gr(j),rr(j),genes);
    fullGPR(i)=gpr; geneList(i)=glist; gprStatus(i)=status;
end
T.full_GPR=fullGPR;T.genes=geneList;T.gpr_status=gprStatus;
writetable(T,fullfile(outDir,'02_transport_reactions_GPR.csv'));
% One row per gene and transporter reaction (AND/OR expression kept).
C=table('Size',[0 6],'VariableTypes',repmat({'string'},1,6), ...
    'VariableNames',{'reaction','transporter_gene','full_GPR','cargo','direction','gpr_status'});
for i=1:height(T)
    ids=split(geneList(i),';'); ids=strtrim(ids);
    ids=ids(ids~="" & ~ismissing(ids));
    for k=1:numel(ids)
        C(end+1,:)={string(T.reaction(i)),ids(k),fullGPR(i), ...
            string(T.cargo(i)),string(T.direction(i)),gprStatus(i)}; %#ok<SAGROW>
    end
end
if height(C)>0,C=unique(C,'rows','stable');end
writetable(C,fullfile(outDir,'03_transporter_gene_candidates.csv'));
ids=unique(C.transporter_gene,'stable');
report(logID,sprintf('Stage 2: %d transport rows with GPR; %d candidate gene rows; %d unique gene IDs', ...
    sum(gprStatus=="recovered_rules" | gprStatus=="grRules"),height(C),numel(ids)));
X.full_GPR=strings(height(X),1);X.transporter_genes=strings(height(X),1);
X.gpr_status=strings(height(X),1);
for xi=1:height(X)
    j=double(X.modelIndex(xi));
    [gpr,glist,status]=decodeGPR(gr(j),rr(j),genes);
    X.full_GPR(xi)=gpr;X.transporter_genes(xi)=glist;X.gpr_status(xi)=status;
end
writetable(X,fullfile(outDir,'02_exchange_transport_candidates_GPR.csv'));
report(logID,sprintf('Stage 2: %d exchange rows have a recovered transporter gene rule', ...
    nnz(X.gpr_status=="recovered_rules" | X.gpr_status=="grRules")));

%% STAGE 3A: reviewed human UniProt mapping (scalar-safe).
% Query MyGene for human gene symbols, then UniProt REST for reviewed proteins.
mapFile=fullfile(outDir,'04_gene_to_protein.csv');
gene_id=ids(:); symbol=strings(numel(ids),1);
protein_name=strings(numel(ids),1);uniprot_accession=strings(numel(ids),1);
mapping_status=strings(numel(ids),1);
for i=1:numel(ids)
    mapping_status(i)="not_mapped";
    if ~ONLINE, mapping_status(i)="offline";continue;end
    try
        u=['https://mygene.info/v3/gene/' char(ids(i)) '?fields=symbol,name,taxid'];
        x=webreadRetry(u,weboptions('Timeout',TIMEOUT,'ContentType','json'),MAX_RETRIES,RETRY_PAUSE);
        if ~isfield(x,'taxid') || double(x.taxid)~=9606
            mapping_status(i)="not_verified_human";continue;
        end
        if ~isfield(x,'symbol'),mapping_status(i)="no_symbol";continue;end
        symbol(i)=string(x.symbol);
        if isfield(x,'name'),protein_name(i)=string(x.name);end
        % Build URL using character vectors to avoid MATLAB string-array
        % implicit expansion and incompatible-size comparisons.
        q=['(gene_exact:' char(symbol(i)) ') AND (organism_id:9606) AND (reviewed:true)'];
        u=['https://rest.uniprot.org/uniprotkb/search?query=' urlencodeLocal(q) ...
            '&format=json&size=20'];
        z=webreadRetry(u,weboptions('Timeout',TIMEOUT,'ContentType','json'),MAX_RETRIES,RETRY_PAUSE);
        acc=strings(0,1);
        if isfield(z,'results')
            hits=z.results;
            for h=1:numel(hits)
                hit=hits(h);
                if ~isfield(hit,'primaryAccession'),continue;end
                exact=false;
                if isfield(hit,'genes')
                    for g=1:numel(hit.genes)
                        gn=hit.genes(g);
                        if isfield(gn,'geneName') && isfield(gn.geneName,'value')
                            if strcmpi(char(string(gn.geneName.value)),char(symbol(i)))
                                exact=true;break;
                            end
                        end
                    end
                end
                if exact,acc(end+1,1)=string(hit.primaryAccession);end %#ok<SAGROW>
            end
        end
        if isempty(acc)
            mapping_status(i)="no_reviewed_exact_gene";
        else
            uniprot_accession(i)=join(unique(acc,'stable'),';');
            mapping_status(i)="reviewed_exact_gene";
        end
    catch ME
        mapping_status(i)="request_failed";
        report(logID,['Mapping ' char(ids(i)) ': ' ME.message]);
    end
    if mod(i,10)==0,report(logID,sprintf('Stage 3A progress: %d/%d genes',i,numel(ids)));end
end
MAP=table(gene_id,symbol,protein_name,uniprot_accession,mapping_status);
writetable(MAP,mapFile);
report(logID,sprintf('Stage 3A: %d/%d genes mapped to reviewed UniProt accessions', ...
    nnz(mapping_status=="reviewed_exact_gene"),numel(ids)));

%% STAGE 3B: obtain valid OmniPath human causal interaction TSV.
% SIGNOR is a source within OmniPath, NOT a valid datasets=signor category.
% Prefer an existing local TSV; otherwise download to a temporary file.
% No network-cache directory or persistent download is created.
rawPath=findNamed(codeRoot,'omnipath_human_raw.tsv');
if isempty(rawPath),rawPath=findNamed(codeRoot,'omnipath_signor_raw.tsv');end
if isempty(rawPath),rawPath=findNamed(codeRoot,'download_attempt_3_response.txt');end
if ~isempty(rawPath)
    report(logID,['Using local regulatory TSV: ' rawPath]);
end
rawTempGuard=[];
if isempty(rawPath) && ONLINE
    try
        url=['https://omnipathdb.org/interactions?datasets=omnipath' ...
             '&organisms=9606&fields=references,sources&format=tsv'];
        txt=webreadRetry(url,weboptions('Timeout',120,'ContentType','text'),MAX_RETRIES,RETRY_PAUSE);
        if isstring(txt),txt=char(txt);end
        firstLine=regexp(txt,'[^\r\n]+','match','once');
        assert(~isempty(firstLine) && contains(firstLine,sprintf('\t')) && ...
            contains(firstLine,'source') && contains(firstLine,'target'), ...
            'Response not valid OmniPath TSV');
        rawPath=[tempname '.tsv'];
        fid=fopen(rawPath,'w');assert(fid>0,'Cannot write temporary TSV');
        fwrite(fid,txt,'char');fclose(fid);
        rawTempGuard=onCleanup(@()delete(rawPath));
        report(logID,'Downloaded OmniPath human interactions (temporary file).');
    catch ME
        report(logID,['OmniPath download failed: ' ME.message]);
        rawPath='';
    end
end
E=table('Size',[0 5],'VariableTypes',repmat({'string'},1,5), ...
    'VariableNames',{'source','target','effect','database','reference'});
if isfile(rawPath)
    try
        R=readtable(rawPath,'FileType','text','Delimiter','\t','TextType','string', ...
            'VariableNamingRule','preserve');
        assert(all(ismember({'source','target','is_stimulation','is_inhibition'}, ...
            R.Properties.VariableNames)),'Unexpected SIGNOR columns');
        for i=1:height(R)
            srcCol=R.('source');dstCol=R.('target');
            stimCol=R.('is_stimulation');inhCol=R.('is_inhibition');
            a=string(srcCol(i)); b=string(dstCol(i));
            st=asBool(stimCol(i));inh=asBool(inhCol(i));
            if ~st && ~inh,continue;end
            if st && inh,eff="ambiguous";elseif st,eff="activation";else,eff="inhibition";end
            ref="";if ismember('references',R.Properties.VariableNames),refCol=R.('references');ref=string(refCol(i));end
            E(end+1,:)={a,b,eff,"OmniPath_SIGNOR",ref}; %#ok<SAGROW>
        end
    catch ME,report(logID,['SIGNOR parse failed: ' ME.message]);end
end
if ~isempty(rawTempGuard),clear rawTempGuard;end
% Optional offline curated edge file: source,target,effect,reference
localEdges=findNamed(codeRoot,'regulatory_edges.csv');
if ~isempty(localEdges)
    try
        L=readtable(localEdges,'TextType','string','VariableNamingRule','preserve');
        if all(ismember({'source','target','effect'},L.Properties.VariableNames))
            for i=1:height(L)
                eff=lower(strtrim(string(L.effect(i))));
                if ~ismember(eff,["activation","inhibition"]),continue;end
                ref="";if ismember('reference',L.Properties.VariableNames),ref=string(L.reference(i));end
                E(end+1,:)={string(L.source(i)),string(L.target(i)),eff,"local",ref}; %#ok<SAGROW>
            end
        end
    catch ME,report(logID,['Local regulatory edges: ' ME.message]);end
end
writetable(E,fullfile(outDir,'05_causal_protein_edges.csv'));

% Diagnostic: coverage of mapped reviewed UniProt accessions in regulatory targets.
allAcc=strings(0,1);
for ii=1:height(MAP)
    aa=split(MAP.uniprot_accession(ii),';');aa=strtrim(aa);
    allAcc=[allAcc;aa(aa~="" & ~ismissing(aa))]; %#ok<AGROW>
end
allAcc=unique(allAcc);
report(logID,sprintf('Stage 3B: %d reviewed accessions; %d occur as regulatory targets', ...
    numel(allAcc),nnz(ismember(allAcc,E.target))));

% Diagnostic: overlap between reviewed protein IDs and regulatory targets.
allAcc=strings(0,1);
for mi=1:height(MAP)
    aa=split(string(MAP.uniprot_accession(mi)),';');
    aa=strtrim(aa);aa=aa(aa~="" & ~ismissing(aa));
    allAcc=[allAcc;aa(:)]; %#ok<AGROW>
end
allAcc=unique(allAcc);
matchedTargets=intersect(allAcc,unique(E.target));
report(logID,sprintf('Stage 3B: %d reviewed accessions; %d occur as regulatory targets', ...
    numel(allAcc),numel(matchedTargets)));

%% STAGE 3C: search backward to MAX_HOPS; sign parity tracks net effect.
P=table('Size',[0 8],'VariableTypes',{'string','string','string','double', ...
    'string','string','string','string'},'VariableNames', ...
    {'transporter_gene','transporter_accession','upstream_protein','hops', ...
     'net_effect','protein_path','edge_effects','references'});
for i=1:height(MAP)
    acc=split(string(MAP.uniprot_accession(i)),';');acc=strtrim(acc);
    for k=1:numel(acc)
        target=acc(k);if target=="" || ismissing(target),continue;end
        qNode=target;qPath=target;qSign=1;qDepth=0;qEff="";qRef="";head=1;
        while head<=numel(qNode)
            node=qNode(head);pth=qPath(head);sg=qSign(head);dep=qDepth(head);
            ef=qEff(head);rf=qRef(head);head=head+1;
            if dep>=MAX_HOPS,continue;end
            incoming=find(E.target==node & ismember(E.effect,["activation","inhibition"]));
            for j=reshape(incoming,1,[])
                up=E.source(j);if any(split(pth,'>')==up),continue;end
                sgn=sg;if E.effect(j)=="inhibition",sgn=-sgn;end
                route=up+">"+pth;newEf=E.effect(j)+">"+ef;
                newRf=E.reference(j)+";"+rf;
                qNode(end+1,1)=up;qPath(end+1,1)=route;qSign(end+1,1)=sgn;
                qDepth(end+1,1)=dep+1;qEff(end+1,1)=newEf;qRef(end+1,1)=newRf;
                if sgn>0,net="net_activation";else,net="net_inhibition";end
                P(end+1,:)={string(MAP.gene_id(i)),target,up,dep+1,net,route,newEf,newRf}; %#ok<SAGROW>
            end
        end
    end
end
writetable(P,fullfile(outDir,'06_backward_protein_paths.csv'));
if height(P)>0 && height(C)>0
    J=innerjoin(C,P,'Keys','transporter_gene');
else
    J=outerjoin(C([],:),P([],:),'Keys','transporter_gene','MergeKeys',true,'Type','left');
end
writetable(J,fullfile(outDir,'07_reaction_protein_paths.csv'));
report(logID,sprintf('Stage 3: %d causal edges; %d backward protein paths; %d net activating', ...
    height(E),height(P),sum(P.net_effect=="net_activation")));

%% STAGE 4: map regulatory proteins to EC numbers and BRENDA activators.
% Protein -> EC via reviewed UniProtKB; EC -> activator/inhibitor via local BRENDA.
% EC annotations are enzyme-activity level evidence, NOT protein-specific proof.
% All unique proteins are included; no arbitrary first-800 restriction.
proteinIDs=unique([allAcc(:); P.upstream_protein(:)],'stable');
proteinIDs=proteinIDs(~ismissing(proteinIDs) & proteinIDs~="" & ...
    ~startsWith(proteinIDs,"COMPLEX:") & ~contains(proteinIDs,"_"));
near=unique([allAcc(:); P.upstream_protein(P.hops<=2)],'stable');
near=near(ismember(near,proteinIDs));
proteinIDs=[near; proteinIDs(~ismember(proteinIDs,near))];
report(logID,sprintf('Stage 4: full EC lookup for %d proteins (no cap)',numel(proteinIDs)));
ECMAP=table('Size',[0 3],'VariableTypes',{'string','string','string'}, ...
    'VariableNames',{'protein_accession','ec_number','status'});
% Optional local override: protein_accession,ec_number (one per row).
ecOverride=findNamed(codeRoot,'protein_ec_mapping.csv');
if ~isempty(ecOverride)
    try
        L=readtable(ecOverride,'TextType','string','VariableNamingRule','preserve');
        if all(ismember({'protein_accession','ec_number'},L.Properties.VariableNames))
            for ii=1:height(L)
                ECMAP(end+1,:)={strtrim(string(L.protein_accession(ii))), ...
                    strtrim(string(L.ec_number(ii))),"local_mapping"};
            end
        end
    catch ME,report(logID,['Protein EC override failed: ' ME.message]);end
end
% Query every UniProt EC batch; no persistent network cache.
missingIDs=proteinIDs(~ismember(proteinIDs,ECMAP.protein_accession));
for first=1:EC_BATCH_SIZE:numel(missingIDs)
    block=missingIDs(first:min(first+EC_BATCH_SIZE-1,numel(missingIDs)));
    raw="";
    if ONLINE
        query=strjoin(cellstr("accession:"+block),' OR ');
        url=['https://rest.uniprot.org/uniprotkb/search?query=' ...
            urlencodeLocal(['(' query ') AND (reviewed:true) AND (organism_id:9606)']) ...
            '&format=tsv&fields=accession,ec&size=100'];
        try
            raw=string(webreadRetry(url,weboptions('Timeout',TIMEOUT,'ContentType','text'), ...
                MAX_RETRIES,RETRY_PAUSE));
            if ~contains(raw,sprintf('\t'))
                error('Invalid UniProt TSV response');
            end
        catch ME
            report(logID,sprintf('Stage 4 EC batch %d failed: %s',first,ME.message));
        end
    end
    if raw~=""
        lines=splitlines(raw);
        for li=2:numel(lines)
            fields=split(lines(li),sprintf('\t'));
            if numel(fields)<2,continue;end
            acc=strtrim(fields(1));
            if ~ismember(acc,block),continue;end
            ecs=regexp(char(fields(2)),'\d+\.\d+\.\d+\.[\d-]+','match');
            for ei=1:numel(ecs)
                ECMAP(end+1,:)={acc,string(ecs{ei}),"UniProt_reviewed"}; %#ok<SAGROW>
            end
        end
    end
    if mod(ceil(first/EC_BATCH_SIZE),10)==0 || first+EC_BATCH_SIZE>numel(missingIDs)
        report(logID,sprintf('Stage 4 EC progress: %d/%d proteins', ...
            min(first+EC_BATCH_SIZE-1,numel(missingIDs)),numel(missingIDs)));
    end
end
if height(ECMAP)>0,ECMAP=unique(ECMAP,'rows','stable');end
writetable(ECMAP,fullfile(outDir,'08_protein_EC_mapping.csv'));
report(logID,sprintf('Stage 4: %d proteins with EC mapping; %d distinct EC numbers', ...
    numel(unique(ECMAP.protein_accession)),numel(unique(ECMAP.ec_number))));

brendaPath=findNamed(codeRoot,'brenda_2026_1.txt');
B=table('Size',[0 5],'VariableTypes',repmat({'string'},1,5), ...
    'VariableNames',{'ec_number','metabolite','effect','organism_tag','brenda_raw'});
if isempty(brendaPath)
    report(logID,'Stage 4: brenda_2026_1.txt NOT FOUND (check code path).');
else
    report(logID,['Stage 4: reading local BRENDA: ' brendaPath]);
    wanted=unique(ECMAP.ec_number);
    fid=fopen(brendaPath,'r');
    if fid<0
        report(logID,'Stage 4: cannot open BRENDA file.');
    else
        cleanupBrenda=onCleanup(@()fclose(fid));
        currentEC=""; nlines=0; nact=0;
        while true
            line=fgetl(fid);if ~ischar(line),break;end
            nlines=nlines+1;
            % BRENDA flat file sections start with ID\tEC_NUMBER.
            id=regexp(line,'^ID\s+([0-9]+\.[0-9]+\.[0-9]+\.[0-9-]+)','tokens','once');
            if ~isempty(id),currentEC=string(id{1});continue;end
            if currentEC=="" || ~ismember(currentEC,wanted),continue;end
            % AC = activating compound, IN = inhibitor in BRENDA flat-file notation.
            hit=regexp(line,'^(AC|IN)\s+(.+)$','tokens','once');
            if isempty(hit),continue;end
            typ=string(hit{1}); raw=string(strtrim(hit{2}));
            % Compound is usually before first reference <...> or comment {...}.
            % Organism number (#n#) is retained for organism-specific review.
            org=regexp(char(raw),'#([^#]+)#','tokens','once');
            if isempty(org),orgTag="";else,orgTag=string(org{1});end
            compound=regexprep(char(raw),'^\s*#[^#]+#\s*','');
            compound=regexprep(compound,'\s*[<{].*$','');
            compound=string(strtrim(compound));
            if compound=="",continue;end
            if typ=="AC",effect="activation";else,effect="inhibition";end
            B(end+1,:)={currentEC,compound,effect,orgTag,raw}; %#ok<SAGROW>
            nact=nact+1;
        end
        clear cleanupBrenda;
        report(logID,sprintf('Stage 4: scanned %d BRENDA lines, found %d AC/IN records for mapped ECs',nlines,nact));
    end
end
writetable(B,fullfile(outDir,'09_BRENDA_enzyme_modulators.csv'));
report(logID,sprintf('Stage 4 BRENDA: %d modulator records written',height(B)));

% Relational join replaces expensive nested table-row append loops.
A=table('Size',[0 7],'VariableTypes',repmat({'string'},1,7), ...
    'VariableNames',{'protein_accession','ec_number','metabolite','effect', ...
    'organism_tag','brenda_raw','evidence_class'});
if ~isempty(ECMAP) && ~isempty(B)
    A=innerjoin(ECMAP(:,{'protein_accession','ec_number'}),B,'Keys','ec_number');
    A.evidence_class=repmat("EC_level_hypothesis_not_protein_specific",height(A),1);
    A=movevars(A,'evidence_class','After','brenda_raw');
end
writetable(A,fullfile(outDir,'10_protein_BRENDA_modulator_candidates.csv'));
report(logID,sprintf('Stage 4: wrote %d protein-EC-modulator candidates',height(A)));

% Join modulator evidence to direct and indirect protein paths.
M=table('Size',[0 9],'VariableTypes',repmat({'string'},1,9), ...
    'VariableNames',{'transporter_gene','transporter_accession','upstream_protein', ...
    'protein_path','metabolite','EC','metabolite_effect','overall_effect','evidence_class'});
% Direct mapping: each UniProt accession may correspond to several gene IDs.
D=table('Size',[0 2],'VariableTypes',{'string','string'}, ...
    'VariableNames',{'transporter_gene','protein_accession'});
for ii=1:height(MAP)
    aa=split(string(MAP.uniprot_accession(ii)),';');
    aa=strtrim(aa);aa=aa(aa~="" & ~ismissing(aa));
    if ~isempty(aa)
        D=[D;table(repmat(string(MAP.gene_id(ii)),numel(aa),1),aa(:), ...
            'VariableNames',D.Properties.VariableNames)]; %#ok<AGROW>
    end
end
if ~isempty(D) && ~isempty(A)
    Q=innerjoin(D,A,'Keys','protein_accession');
    MD=table(Q.transporter_gene,Q.protein_accession,Q.protein_accession, ...
        Q.protein_accession,Q.metabolite,Q.ec_number,Q.effect,Q.effect,Q.evidence_class, ...
        'VariableNames',M.Properties.VariableNames);
    M=[M;MD];
end
report(logID,sprintf('Stage 4 direct candidates: %d',height(M)));
if ~isempty(P) && ~isempty(A)
    % Join by the actual regulatory protein, not by transporter accession.
    Q=innerjoin(P,A,'LeftKeys','upstream_protein','RightKeys','protein_accession');
    positive=(Q.net_effect=="net_activation");
    positive(Q.effect=="inhibition")=~positive(Q.effect=="inhibition");
    overall=repmat("net_inhibition",height(Q),1);overall(positive)="net_activation";
    MI=table(Q.transporter_gene,Q.transporter_accession,Q.upstream_protein, ...
        Q.protein_path,Q.metabolite,Q.ec_number,Q.effect,overall,Q.evidence_class, ...
        'VariableNames',M.Properties.VariableNames);
    M=[M;MI];
end
report(logID,sprintf('Stage 4 joined candidates: %d',height(M)));
writetable(M,fullfile(outDir,'11_BRENDA_to_transporter_hypotheses.csv'));
report(logID,sprintf('Stage 4: %d protein-EC-modulator candidates; %d putative regulatory path matches',height(A),height(M)));
report(logID,'WARNING: BRENDA EC-level modulation is NOT verified protein-specific, human/T-cell-specific, or direct transporter regulation.');
report(logID,'WARNING: compound names are not yet matched to Puniya metabolite IDs; no D1 pathway is established.');
report(logID,'Do not interpret missing paths as biological absence.');
report(logID,'Transport direction uses bounds, but physiological import still requires verification.');
report(logID,'GPR membership does not prove a gene product is a direct membrane transporter.');
%% STAGE 5: clean BRENDA names, match Puniya metabolites, rank and trace.
% The results are HYPOTHESES. Enzyme modulation does not imply transport
% stimulation, and reaction connectivity does not imply flux feasibility.
report(logID,'Stage 5: cleaning BRENDA compounds and matching Puniya metabolites.');
if isfield(model,'metNames') && numel(model.metNames)==numel(model.mets)
    metNames=string(model.metNames(:));
else
    metNames=string(model.mets(:));
end
metIDs=string(model.mets(:));
[metBase,metComp]=splitCompartments(metIDs);
metKey=strings(numel(metIDs),1);
for ii=1:numel(metIDs)
    metKey(ii)=normalizeCompound(metNames(ii));
end
% Also accept exact base IDs; never perform fuzzy matches automatically.
baseKey=strings(numel(metIDs),1);
for ii=1:numel(metIDs)
    baseKey(ii)=normalizeCompound(metBase(ii));
end
cleanName=strings(height(B),1); cleanKey=strings(height(B),1);
cleanStatus=strings(height(B),1);
for ii=1:height(B)
    [cleanName(ii),cleanStatus(ii)]=cleanBrendaName(B.metabolite(ii));
    cleanKey(ii)=normalizeCompound(cleanName(ii));
end
BC=B;
BC.clean_compound=cleanName;
BC.normalized_compound=cleanKey;
BC.clean_status=cleanStatus;
writetable(BC,fullfile(outDir,'12_clean_BRENDA_modulators.csv'));
report(logID,sprintf('Stage 5: %d/%d BRENDA records have usable cleaned names', ...
    nnz(cleanStatus=="usable"),height(B)));

% Build a compact one-row-per-unique-name matching table, then join by
% normalized name. Matching the raw 1.2M hypotheses first is wasteful.
keys=unique(cleanKey(cleanStatus=="usable" & cleanKey~=""),'stable');
matchKey=strings(0,1); matchName=strings(0,1);
matchID=strings(0,1); matchComp=strings(0,1); matchIndex=zeros(0,1);
matchMethod=strings(0,1);
for ii=1:numel(keys)
    k=keys(ii);
    hit=find(metKey==k & metKey~="");
    method="exact_metabolite_name";
    if isempty(hit)
        hit=find(baseKey==k & baseKey~="");
        method="exact_metabolite_id";
    end
    for jj=reshape(hit,1,[])
        matchKey(end+1,1)=k; %#ok<SAGROW>
        matchName(end+1,1)=metNames(jj); %#ok<SAGROW>
        matchID(end+1,1)=metIDs(jj); %#ok<SAGROW>
        matchComp(end+1,1)=metComp(jj); %#ok<SAGROW>
        matchIndex(end+1,1)=jj; %#ok<SAGROW>
        matchMethod(end+1,1)=method; %#ok<SAGROW>
    end
end
MATCH=table(matchKey,matchName,matchID,matchComp,matchIndex,matchMethod, ...
    'VariableNames',{'normalized_compound','model_metabolite_name', ...
    'model_metabolite_id','compartment','met_index','match_method'});
writetable(MATCH,fullfile(outDir,'13_Puniya_metabolite_matches.csv'));
report(logID,sprintf('Stage 5: %d compound keys matched to %d model metabolite rows', ...
    numel(unique(MATCH.normalized_compound)),height(MATCH)));

% Keep only BRENDA records with an exact Puniya name/ID match. The
% organism tag is NOT a species label; no human-specific assertion is made.
RANK=table();
if ~isempty(M) && ~isempty(MATCH)
    hKey=strings(height(M),1);
    for ii=1:height(M)
        [nm,st]=cleanBrendaName(M.metabolite(ii));
        if st=="usable",hKey(ii)=normalizeCompound(nm);end
    end
    H=M(:,{'transporter_gene','transporter_accession','upstream_protein', ...
        'protein_path','metabolite','EC','metabolite_effect', ...
        'overall_effect','evidence_class'});
    H.normalized_compound=hKey;
    H=H(H.normalized_compound~="" & H.overall_effect=="net_activation",:);
    if height(H)>STAGE5_MAX_MATCH_ROWS
        % Stable first rows are only a computational subset, not a ranked
        % selection. State the truncation explicitly in the report.
        report(logID,sprintf('Stage 5 WARNING: limiting %d eligible rows to first %d before join', ...
            height(H),STAGE5_MAX_MATCH_ROWS));
        H=H(1:STAGE5_MAX_MATCH_ROWS,:);
    end
    RANK=innerjoin(H,MATCH,'Keys','normalized_compound');
    if ~isempty(RANK)
        RANK.hops=zeros(height(RANK),1);
        for ii=1:height(RANK)
            RANK.hops(ii)=numel(split(RANK.protein_path(ii),'>'))-1;
        end
        RANK.rank_score=100-15*RANK.hops;
        RANK.rank_score(RANK.compartment=="c")=RANK.rank_score(RANK.compartment=="c")+10;
        RANK.rank_score(RANK.match_method=="exact_metabolite_name")= ...
            RANK.rank_score(RANK.match_method=="exact_metabolite_name")+5;
        RANK=sortrows(RANK,{'rank_score','hops'},{'descend','ascend'});
    end
end
if isempty(RANK) && width(RANK)==0
    RANK=table('Size',[0 2],'VariableTypes',{'string','string'}, ...
        'VariableNames',{'model_metabolite_id','status'});
end
writetable(RANK,fullfile(outDir,'14_regulatory_candidates_ranked.csv'));
report(logID,sprintf('Stage 5: %d net-activating EC-level hypotheses matched to model metabolites',height(RANK)));

% Search reaction orientations backwards from each unique matched metabolite.
% A step is admitted only when the model bounds allow that orientation.
% Cofactors and additional reactants are NOT guaranteed available.
PATH=table('Size',[0 7],'VariableTypes', ...
    {'string','string','string','double','string','string','string'}, ...
    'VariableNames',{'target_metabolite','precursor_metabolite', ...
    'reaction_path','steps','reaction_directions','feasibility','target_compartment'});
if ~isempty(RANK)
    targetIdx=unique(RANK.met_index,'stable');
    for tt=1:numel(targetIdx)
        target=targetIdx(tt);
        queueNode=target;queuePath="";queueDirs="";queueDepth=0;
        head=1;seen=string(target)+"|0";
        while head<=numel(queueNode) && height(PATH)<STAGE5_MAX_PATH_ROWS
            curr=queueNode(head); pth=queuePath(head);
            dirs=queueDirs(head); dep=queueDepth(head);head=head+1;
            if dep>=STAGE5_MAX_HOPS,continue;end
            rrIdx=find(model.S(curr,:)~=0);
            for r=reshape(rrIdx,1,[])
                coeff=full(model.S(curr,r));
                if coeff>0 && ub(r)>0
                    orient=1; label="forward";
                else
                    continue;
                end
                subs=find(full(model.S(:,r))*orient < 0);
                for pre=reshape(subs,1,[])
                    if pre==curr,continue;end
                    newPath=string(rxnIDs(r));newDirs=label;
                    if pth~="",newPath=newPath+">"+pth;newDirs=label+">"+dirs;end
                    PATH(end+1,:)={metIDs(target),metIDs(pre),newPath,dep+1, ...
                        newDirs,"structural_only_not_flux_verified",metComp(target)}; %#ok<SAGROW>
                    if dep+1<STAGE5_MAX_HOPS
                        key=string(pre)+"|"+string(dep+1);
                        if ~ismember(key,seen)
                            seen(end+1,1)=key; %#ok<SAGROW>
                            queueNode(end+1,1)=pre;queuePath(end+1,1)=newPath; %#ok<SAGROW>
                            queueDirs(end+1,1)=newDirs;queueDepth(end+1,1)=dep+1; %#ok<SAGROW>
                        end
                    end
                    if height(PATH)>=STAGE5_MAX_PATH_ROWS,break;end
                end
                if height(PATH)>=STAGE5_MAX_PATH_ROWS,break;end
            end
        end
        if mod(tt,25)==0
            report(logID,sprintf('Stage 5 backward progress: %d/%d target metabolites',tt,numel(targetIdx)));
        end
        if height(PATH)>=STAGE5_MAX_PATH_ROWS
            report(logID,'Stage 5 WARNING: structural path output reached safety cap.');break;
        end
    end
end
writetable(PATH,fullfile(outDir,'15_backward_metabolic_paths.csv'));

% Endpoint candidates are metabolites with an incoming allowable reaction.
% They are NOT validated D1 motifs, and no transient kinetics are claimed.
ENDPOINT=table();
if ~isempty(RANK) && ~isempty(PATH)
    validTargets=unique(PATH.target_metabolite);
    ENDPOINT=RANK(ismember(RANK.model_metabolite_id,validTargets),:);
    if ~isempty(ENDPOINT)
        ENDPOINT=unique(ENDPOINT,'rows','stable');
        ENDPOINT.d1_status=repmat("candidate_endpoint_only_D1_not_tested",height(ENDPOINT),1);
    end
end
if isempty(ENDPOINT) && width(ENDPOINT)==0
    ENDPOINT=table('Size',[0 2],'VariableTypes',{'string','string'}, ...
        'VariableNames',{'model_metabolite_id','d1_status'});
end
writetable(ENDPOINT,fullfile(outDir,'16_D1_endpoint_candidates.csv'));
report(logID,sprintf('Stage 5: %d structural backward paths; %d endpoint candidate rows', ...
    height(PATH),height(ENDPOINT)));
report(logID,'Stage 5 WARNING: exact names are candidate mappings, not verified chemical identities.');
report(logID,'Stage 5 WARNING: no flux-feasibility or D1 kinetic differentiation has been established.');

%% STAGE 6: print readable candidate connections (NOT complete D1/counter circuits).
% Uses the same section-oriented format as First_Selected_Full_Connection.txt.
% Selects up to ten distinct gene / regulatory metabolite / cargo combinations.
% The actual D1 motif and downstream counter have NOT been searched here.
report(logID,'Stage 6: printing selected regulatory-to-transport candidate connections.');
N_PRINT_CONNECTIONS=10;
reportPath=fullfile(outDir,'17_Selected_Transport_Connections.txt');
fid6=fopen(reportPath,'w');
assert(fid6>0,'Could not create Stage 6 text report');
close6=onCleanup(@()fclose(fid6));
if EXCHANGE_ONLY
    fprintf(fid6,'REGULATORY PROTEIN-PATH CONNECTIONS EXCLUDED\n');
else
    fprintf(fid6,'SELECTED CANDIDATE METABOLIC CONNECTIONS\n');
    fprintf(fid6,'D1 ENDPOINT -> REGULATORY PROTEIN(S) -> MEMBRANE TRANSPORT\n');
end
fprintf(fid6,'IMPORTANT: D1 motif and downstream counter are NOT yet reconstructed.\n');
fprintf(fid6,'All chemical reaction directions are +1 relative to original model.S.\n');
fprintf(fid6,'BRENDA enzyme activity evidence is EC-level and NOT transporter-specific proof.\n\n');
nPrinted=0;
seenConnections=strings(0,1);
if ~EXCHANGE_ONLY && ~isempty(RANK) && ~isempty(C) && ~isempty(T)
    for si=1:height(RANK)
        if nPrinted>=N_PRINT_CONNECTIONS,break;end
        rec=RANK(si,:);
        if ~ismember('transporter_gene',RANK.Properties.VariableNames),break;end
        candidateC=find(C.transporter_gene==rec.transporter_gene,1,'first');
        if isempty(candidateC),continue;end
        % Prefer a transporter carrying cargo distinct from the regulator.
        choices=find(C.transporter_gene==rec.transporter_gene);
        chosen=0;chosenT=0;
        for ci=reshape(choices,1,[])
            if C.cargo(ci)==rec.normalized_compound,continue;end
            ti=find(string(T.reaction)==C.reaction(ci) & string(T.direction)==C.direction(ci),1,'first');
            if isempty(ti),continue;end
            key=rec.transporter_gene+"|"+rec.model_metabolite_id+"|"+C.reaction(ci)+"|"+C.direction(ci);
            if ismember(key,seenConnections),continue;end
            chosen=ci;chosenT=ti;break;
        end
        if chosen==0,continue;end
        seenConnections(end+1,1)=key; %#ok<SAGROW>
        nPrinted=nPrinted+1;
        fprintf(fid6,'============================================================\n');
        fprintf(fid6,'SELECTED CANDIDATE CONNECTION %d\n',nPrinted);
        fprintf(fid6,'============================================================\n');
        fprintf(fid6,'SELECTION AND EVIDENCE\n');
        fprintf(fid6,'  Rank score: %.1f (heuristic, NOT evidence strength)\n',rec.rank_score);
        fprintf(fid6,'  Status: EC-level hypothesis; biological validation required\n');
        fprintf(fid6,'  Transporter-associated gene: %s\n',char(rec.transporter_gene));
        fprintf(fid6,'  UniProt accession: %s\n',char(rec.transporter_accession));
        fprintf(fid6,'  BRENDA EC: %s\n',char(rec.EC));
        fprintf(fid6,'  BRENDA compound: %s\n',char(rec.metabolite));
        fprintf(fid6,'  Regulatory effect on transport (sign composition): %s\n',char(rec.overall_effect));
        fprintf(fid6,'  Protein pathway: %s\n\n',char(rec.protein_path));
        fprintf(fid6,'COMPLETE SIGNAL-METABOLITE MAP (IDENTIFIED PART)\n');
        fprintf(fid6,'  m_Sd_1 candidate: %s | %s | compartment=%s\n',char(rec.model_metabolite_id),char(rec.model_metabolite_name),char(rec.compartment));
        fprintf(fid6,'  m_c candidate: %s | %s\n',char(T.cytosolic(chosenT)),char(T.cargo(chosenT)));
        fprintf(fid6,'  m_in: NOT IDENTIFIED; D1 differentiator NOT VERIFIED\n');
        fprintf(fid6,'  m_out: NOT IDENTIFIED; counter NOT SEARCHED\n\n');
        fprintf(fid6,'FULL FUNCTIONAL CONNECTION (PARTIAL)\n');
        fprintf(fid6,'  Candidate D1 endpoint -> regulatory pathway: %s -> %s\n',char(rec.model_metabolite_id),char(rec.protein_path));
        fprintf(fid6,'  Transport: %s -> %s via %s\n',char(T.extracellular(chosenT)),char(T.cytosolic(chosenT)),char(T.reaction(chosenT)));
        fprintf(fid6,'  WARNING: the metabolite-to-protein activation and resulting transport flux are NOT verified.\n\n');
        fprintf(fid6,'BACKWARD METABOLIC REACTIONS TO CANDIDATE m_Sd_1\n');
        pidx=find(PATH.target_metabolite==rec.model_metabolite_id,1,'first');
        if isempty(pidx)
            fprintf(fid6,'  No structural backward path retained under the Stage 5 cap.\n');
        else
            fprintf(fid6,'  Precursor: %s\n',char(PATH.precursor_metabolite(pidx)));
            rp=split(PATH.reaction_path(pidx),'>');
            rd=split(PATH.reaction_directions(pidx),'>');
            for kk=1:numel(rp)
                ri=find(rxnIDs==rp(kk),1,'first');
                if isempty(ri),continue;end
                assert(kk<=numel(rd) && rd(kk)=="forward", ...
                    'ForwardOnly:PathDirection','A metabolic path contains a non-forward or missing direction.');
                orient=1;
                printModelReaction(fid6,model,ri,orient,metIDs,metNames,rec.model_metabolite_id,'m_Sd_1');
            end
            fprintf(fid6,'  Feasibility: structural only; cofactors and flux not verified.\n');
        end
        fprintf(fid6,'\nMEMBRANE TRANSPORT REACTION\n');
        ri=double(T.modelIndex(chosenT));
        assert(string(T.direction(chosenT))=="forward", ...
            'ForwardOnly:TransportDirection','A transporter has non-forward direction.');
        orient=1;
        printModelReaction(fid6,model,ri,orient,metIDs,metNames,T.cytosolic(chosenT),'m_c');
        fprintf(fid6,'  GPR: %s\n',char(C.full_GPR(chosen)));
        fprintf(fid6,'  Gene association is NOT proof of direct transport protein identity.\n\n');
        fprintf(fid6,'COUNTER REACTIONS\n');
        fprintf(fid6,'  r_ac: NOT IDENTIFIED (no counter-production reaction selected)\n');
        fprintf(fid6,'  r_cr: NOT IDENTIFIED (no counter-consumption reaction selected)\n');
        fprintf(fid6,'  No complete D1 -> counter chemical circuit is claimed.\n\n');
    end
end
if nPrinted==0
    if EXCHANGE_ONLY
        fprintf(fid6,'Regulatory protein-path connections were excluded by EXCHANGE_ONLY=true.\n');
        fprintf(fid6,'See 18b_Exchange_D1_Connections.txt for the exchange-mediated search.\n');
    else
        fprintf(fid6,'No candidate could be linked to a transport reaction.\n');
    end
end
clear close6;
report(logID,sprintf('Stage 6: wrote %d candidate connections to %s',nPrinted,reportPath));

%% STAGE 7: constrained SEVEN-reaction D1 motif reconstruction.
% Backward from a regulatory metabolite m_Sd_1, then seek a parallel
% fast rail from the SAME m_in. A structural motif is not a kinetic D1.
report(logID,'Stage 7: searching seven-reaction D1 topology with equal-stoichiometry cancellation.');
D1_REQUIRE_CANCELLER_PRODUCTS=true; % true: r_I_1 must form at least one product
D1_SKIP_TRANSPORT=true;       % use intracellular conversions for D1
outD1=fullfile(outDir,'18_Selected_Complete_D1_Transport_Connections.txt');
fid7=fopen(outD1,'w');assert(fid7>0,'Cannot create D1 text report');
close7=onCleanup(@()fclose(fid7));
fprintf(fid7,'REGULATORY-PATH D1 DIAGNOSTIC (DISABLED IN EXCHANGE-ONLY MODE)\n');
fprintf(fid7,'SEVEN D1 REACTIONS: r_Ff_1, r_Fd_1, r_Fr_1, r_Sf_1, r_Sd_1, r_Sr_1, r_I_1\n');
fprintf(fid7,'Direction policy: STRICT +1 relative to original model.S.\n');
fprintf(fid7,'REGULATORY SIGNAL: m_Sd_1 is produced by r_Sd_1, NOT by r_I_1.\n');
fprintf(fid7,'WARNING: topology only; no kinetic differentiation or flux feasibility verified.\n');
fprintf(fid7,'WARNING: BRENDA EC-level effects do not prove transport regulation.\n\n');
% Precompute allowable oriented intracellular conversions once.
[di,dm,do,ds]=d1BuildOriented(model,comp,lb,ub,D1_SKIP_TRANSPORT);
assert(all(dm==1),'ForwardOnly:OrientedPool', ...
    'Intracellular reaction pool contains a non-forward orientation.');
% di = reaction indices; dm = orientations; do = product lists; ds = substrate lists.
% Exclude ubiquitous currency metabolites from the SIGNAL intermediates.
currency=lower(["h2o","h","atp","adp","pi","ppi","nad","nadh", ...
    "nadp","nadph","coa","co2","o2","na1","k","cl"]);
% AMP and NH4 remain eligible as possible signals. Edit currency exclusions as needed.
% m_Sd candidates must be producible in the D1 reaction pool. Block 0
% reachability is applied later to m_in inside d1FindMotifs, not to m_Sd.
d1TargetProduced=false(numel(metIDs),1);
for d1q=1:numel(do)
    d1TargetProduced(do{d1q})=true;
end
% WARNING: no candidate is accepted without a real stoichiometric canceler.
foundD1=0; inspected=0; testedTotal=0; seenD1=strings(0,1);
d1TestLimitedCalls=0;d1BranchLimitedCalls=0;d1RetentionLimitedCalls=0;
discoveredD1Keys=containers.Map('KeyType','char','ValueType','logical');
selectedD1=struct('d',{},'target',{},'transportRow',{},'regulatorRow',{},'exchangeRow',{}, ...
    'exchangePathMets',{},'exchangePathRxns',{});
if ~EXCHANGE_ONLY && ~isempty(RANK) && ~isempty(C) && ~isempty(T)
    for si=1:height(RANK)
        if foundD1>=D1_MAX_MOTIFS || inspected>=D1_MAX_TARGETS,break;end
        rec=RANK(si,:);
        if ~ismember('transporter_gene',RANK.Properties.VariableNames),break;end
        if rec.compartment~="c",continue;end
        target=find(metIDs==rec.model_metabolite_id,1,'first');
        if isempty(target),continue;end
        crows=find(C.transporter_gene==rec.transporter_gene);
        chosenT=0; chosenC=0;
        for cr=reshape(crows,1,[])
            ti=find(string(T.reaction)==C.reaction(cr) & string(T.direction)==C.direction(cr),1);
            if ~isempty(ti) && T.cytosolic(ti)~=rec.model_metabolite_id
                chosenT=ti;chosenC=cr;break;
            end
        end
        if chosenT==0,continue;end
        key=rec.model_metabolite_id+"|"+rec.transporter_gene+"|"+T.reaction(chosenT);
        if ismember(key,seenD1),continue;end
        seenD1(end+1,1)=key;inspected=inspected+1; %#ok<SAGROW>
        [motifs,nTests,d1Stats]=d1FindMotifs(target,di,dm,do,ds,metIDs,comp,currency, ...
            D1_MAX_BRANCH,D1_MAX_RXN_CAND,D1_MAX_TESTS_PER_TARGET,D1_MAX_PER_TARGET, ...
            model.S,D1_REQUIRE_CANCELLER_PRODUCTS,block0InputAllowed,block0Cache.depth);
        testedTotal=testedTotal+nTests;
        d1TestLimitedCalls=d1TestLimitedCalls+d1Stats.testLimitReached;
        d1BranchLimitedCalls=d1BranchLimitedCalls+d1Stats.branchLimitReached;
        d1RetentionLimitedCalls=d1RetentionLimitedCalls+d1Stats.retentionLimitReached;
        mergeChemicalKeys(discoveredD1Keys,d1Stats.chemicalKeys);
        if isempty(motifs),continue;end
        for mi=1:numel(motifs)
            if foundD1>=D1_MAX_MOTIFS,break;end
            d=motifs(mi);
            d1ids=[d.rff,d.rfd,d.rfr,d.rsf,d.rsd,d.rsr,d.ri];
            [~,block0Reactions]=block0CachedPath(block0Cache,d.min);
            motifPrefix=[block0Reactions,d1ids,double(T.modelIndex(chosenT))];
            if reactionHasExcludedEC(model,motifPrefix,EXCLUDED_MOTIF_ECS),continue;end
            assert(block0InputAllowed(d.min),'Ranking:UnconnectedD1','Unconnected input cannot be selected.');
            foundD1=foundD1+1;
            selectedD1(end+1)=struct('d',d,'target',target,'transportRow',chosenT, ...
                'regulatorRow',si,'exchangeRow',0, ...
                'exchangePathMets',zeros(1,0),'exchangePathRxns',zeros(1,0)); %#ok<SAGROW>
            fprintf(fid7,'============================================================\n');
            fprintf(fid7,'STRUCTURAL D1 CANDIDATE %d\n',foundD1);
            fprintf(fid7,'============================================================\n');
            fprintf(fid7,'SELECTED SIGNAL-METABOLITE MAP\n');
            dlabels=["m_in","m_Ff_1","m_Fd_1","m_Sf_1","m_Sd_1"];
            mlist=[d.min,d.mff,d.mfd,d.msf,target];
            for kk=1:5
                fprintf(fid7,'  %s: %s | %s\n',char(dlabels(kk)),char(metIDs(mlist(kk))),char(metNames(mlist(kk))));
            end
            fprintf(fid7,'  m_c: %s\n',char(T.cytosolic(chosenT)));
            fprintf(fid7,'  m_out: NOT IDENTIFIED; downstream counter not searched\n\n');
            fprintf(fid7,'FULL FUNCTIONAL CONNECTION (STRUCTURAL ONLY)\n');
            fprintf(fid7,'  Fast rail: %s -> %s -> %s -> removal\n', ...
                char(metIDs(d.min)),char(metIDs(d.mff)),char(metIDs(d.mfd)));
            fprintf(fid7,'  Slow rail: %s -> %s -> %s -> removal\n', ...
                char(metIDs(d.min)),char(metIDs(d.msf)),char(metIDs(target)));
            fprintf(fid7,'  Cancellation r_I_1 consumes BOTH %s and %s in equal stoichiometry.\n', ...
                char(metIDs(d.mfd)),char(metIDs(target)));
            fprintf(fid7,'  m_Sd_1 is produced by r_Sd_1 and branches to regulation BEFORE cancellation/removal.\n');
            fprintf(fid7,'  Hypothetical regulation (non-stoichiometric): %s -> %s\n',char(metIDs(target)),char(rec.protein_path));
            fprintf(fid7,'  Transport: %s -> %s via %s\n\n', ...
                char(T.extracellular(chosenT)),char(T.cytosolic(chosenT)),char(T.reaction(chosenT)));
            fprintf(fid7,'DIFFERENTIATOR 1 REACTIONS\n');
            labels=["r_Ff_1","r_Fd_1","r_Fr_1","r_Sf_1","r_Sd_1","r_Sr_1","r_I_1"];
            rids=[d.rff,d.rfd,d.rfr,d.rsf,d.rsd,d.rsr,d.ri];
            ror=[d.off,d.ofd,d.ofr,d.osf,d.osd,d.osr,d.oi];
            for kk=1:7
                fprintf(fid7,'%s\n',char(labels(kk)));
                printModelReaction(fid7,model,rids(kk),ror(kk),metIDs,metNames,metIDs(target),'m_Sd_1');
                fprintf(fid7,'  EC: %s\n',char(d1EC(model,rids(kk))));
                [rule,gl,~]=decodeGPR(gr(rids(kk)),rr(rids(kk)),genes);
                fprintf(fid7,'  GPR: %s\n  Resolved genes: %s\n\n',char(rule),char(gl));
            end
            fprintf(fid7,'D1 SIGNAL AND CANCELLER CHECK\n');
            fprintf(fid7,'  r_Sd_1 produces m_Sd_1: %s\n',char(metIDs(target)));
            fprintf(fid7,'  r_I_1 consumes m_Fd_1 and m_Sd_1 equally: %.8g each\n',d.cancelCoeff);
            fprintf(fid7,'  r_Fr_1 consumes m_Fd_1; r_Sr_1 consumes m_Sd_1.\n');
            fprintf(fid7,'  Regulatory signal is m_Sd_1, NOT a canceler product.\n\n');
            fprintf(fid7,'REGULATORY INTERFACE (NOT CHEMICAL STOICHIOMETRY)\n');
            fprintf(fid7,'  BRENDA compound: %s\n  EC: %s\n',char(rec.metabolite),char(rec.EC));
            fprintf(fid7,'  Protein pathway: %s\n  Sign: %s\n',char(rec.protein_path),char(rec.overall_effect));
            fprintf(fid7,'  Status: EC-level hypothesis, not verified for this human protein/cargo.\n\n');
            fprintf(fid7,'MEMBRANE TRANSPORT REACTION\n');
            assert(string(T.direction(chosenT))=="forward", ...
            'ForwardOnly:TransportDirection','A transporter has non-forward direction.');
        orient=1;
            printModelReaction(fid7,model,double(T.modelIndex(chosenT)),orient,metIDs,metNames,T.cytosolic(chosenT),'m_c');
            fprintf(fid7,'  GPR: %s\n\n',char(C.full_GPR(chosenC)));
            fprintf(fid7,'COUNTER REACTIONS\n');
            fprintf(fid7,'  r_ac: NOT IDENTIFIED\n  r_cr: NOT IDENTIFIED\n');
            fprintf(fid7,'  No complete D1-to-counter circuit or transient response is claimed.\n\n');
        end
        if mod(inspected,20)==0
            report(logID,sprintf('Stage 7 progress: %d targets; %d structural motifs; %d trials', ...
                inspected,foundD1,testedTotal));
        end
    end
end
if EXCHANGE_ONLY
    fprintf(fid7,'Regulatory-protein-path D1 search was disabled (EXCHANGE_ONLY=true).\n');
    fprintf(fid7,'Only exchange-mediated transport candidates are considered for complete motifs.\n');
elseif foundD1==0
    fprintf(fid7,'No complete seven-reaction D1 topology WITH a canceler found within configured bounds.\n');
    fprintf(fid7,'This does not prove no D1 exists; increase Stage 7 limits or revise constraints.\n');
end
clear close7;
report(logID,sprintf('Stage 7: inspected %d endpoints; tested %d partial motifs; printed %d seven-reaction motifs', ...
    inspected,testedTotal,foundD1));
report(logID,['Stage 7 file: ' outD1]);

%% STAGE 7B: anchor at exchange transport, then search backward for m_Sd_1.
% Enumerate each forward stoichiometric exchange pair first. Starting from
% its intracellular exported counter-substrate, traverse the forward-reaction
% graph backward to candidate m_Sd_1 species. Search each candidate as the
% D1 endpoint; Block 0 reachability constrains m_in, not m_Sd_1.
% Only direct (0-step) and one-reaction (1-step) m_Sd-to-counter paths are
% searched by default. Transporter regulation/expression is not inferred.
EXCHANGE_D1_MAX_TESTS_PER_TARGET=D1_MAX_TESTS_PER_TARGET;
EXCHANGE_D1_MAX_BRANCH=D1_MAX_BRANCH;
EXCHANGE_D1_MAX_RXN_CAND=D1_MAX_RXN_CAND;
outExchangeD1=fullfile(outDir,'18b_Exchange_D1_Connections.txt');
fid7x=fopen(outExchangeD1,'w');assert(fid7x>0,'Cannot create exchange D1 report');
close7x=onCleanup(@()fclose(fid7x));
fprintf(fid7x,'STOICHIOMETRIC EXCHANGE-TRANSPORT D1 CANDIDATES\n');
fprintf(fid7x,'Search order: exchange pair -> backward precursor search for m_Sd_1 -> D1 -> Block 0 and neurons.\n');
fprintf(fid7x,'No protein-transcription regulation step is inserted in this route.\n');
fprintf(fid7x,'m_Sd-to-counter reaction-distance cap=%d; precursor candidates per exchange pair=%g.\n', ...
    EXCHANGE_METABOLIC_MAX_DEPTH,EXCHANGE_MAX_SIGNAL_PRECURSORS);
fprintf(fid7x,'Require transporter GPR: %d; expression data are not included in the input.\n',EXCHANGE_REQUIRE_GPR);
fprintf(fid7x,'Direction policy: STRICT +1 relative to original model.S.\n');
fprintf(fid7x,'Transporter presence, exchange kinetics, and flux feasibility require validation.\n\n');
exchangeMatchRows=cell(0,15);
exchangeFound=0;exchangeInspected=0;exchangeTests=0;
exchangeMotifCache=containers.Map('KeyType','char','ValueType','any');
if height(X)>0
    % Search each transport/counter pair independently so every transporter
    % receives its own closest successful m_Sd path.
    exchangeTargetRows=(1:height(X))';
    for xt=1:min(numel(exchangeTargetRows),EXCHANGE_D1_MAX_TARGETS)
        xiTarget=exchangeTargetRows(xt);
        counterID=string(X.counterSubstrate(xiTarget));
        counterMet=find(metIDs==counterID,1,'first');
        if isempty(counterMet) || comp(counterMet)~="c",continue;end
        counterBase=regexprep(char(counterID),'\[c\]$','');
        if ismember(lower(string(counterBase)),currency),continue;end
        pairGPR=strtrim(string(X.full_GPR(xiTarget)));
        if EXCHANGE_REQUIRE_GPR && (ismissing(pairGPR) || pairGPR=="" || ...
                string(X.gpr_status(xiTarget))=="no_rule")
            continue;
        end
        exchangeInspected=exchangeInspected+1;
        [signalCandidates,pathMetsBySource,pathRxnsBySource]=findForwardPrecursorPaths( ...
            counterMet,di,do,ds,metIDs,comp,currency,EXCHANGE_METABOLIC_MAX_DEPTH);
        % Keep only metabolites that can actually be products of a D1-pool
        % reaction. Do NOT require m_Sd itself to be cAMP-reachable as m_in.
        keepSignal=d1TargetProduced(signalCandidates);
        signalCandidates=signalCandidates(keepSignal);
        pathMetsBySource=pathMetsBySource(keepSignal);
        pathRxnsBySource=pathRxnsBySource(keepSignal);
        if isempty(signalCandidates),continue;end
        candidateDepth=zeros(numel(signalCandidates),1);
        for sc=1:numel(signalCandidates)
            candidateDepth(sc)=numel(pathRxnsBySource{sc});
        end
        % Nearest chemistry first; metabolite index is only a stable tie-break.
        [~,candidateOrder]=sortrows([candidateDepth,double(signalCandidates(:))],[1 2]);
        signalCandidates=signalCandidates(candidateOrder);
        pathMetsBySource=pathMetsBySource(candidateOrder);
        pathRxnsBySource=pathRxnsBySource(candidateOrder);
        if isfinite(EXCHANGE_MAX_SIGNAL_PRECURSORS)
            keepN=min(numel(signalCandidates),EXCHANGE_MAX_SIGNAL_PRECURSORS);
            signalCandidates=signalCandidates(1:keepN);
            pathMetsBySource=pathMetsBySource(1:keepN);
            pathRxnsBySource=pathRxnsBySource(1:keepN);
        end
        matchedX=xiTarget;
        bestCounterPathDepth=inf;
        for si=1:numel(signalCandidates)
            if exchangeFound>=EXCHANGE_D1_MAX_CONNECTIONS,break;end
            currentPathDepth=numel(pathRxnsBySource{si});
            if currentPathDepth>bestCounterPathDepth,break;end
            signalMet=signalCandidates(si);
            signalKey=char(metIDs(signalMet));
            if ~isKey(exchangeMotifCache,signalKey)
                [signalMotifs,nTests,d1Stats]=d1FindMotifs(signalMet,di,dm,do,ds,metIDs,comp,currency, ...
                    EXCHANGE_D1_MAX_BRANCH,EXCHANGE_D1_MAX_RXN_CAND, ...
                    EXCHANGE_D1_MAX_TESTS_PER_TARGET,D1_MAX_PER_TARGET,model.S, ...
                    D1_REQUIRE_CANCELLER_PRODUCTS,block0InputAllowed,block0Cache.depth);
                exchangeMotifCache(signalKey)=signalMotifs;
                exchangeTests=exchangeTests+nTests;
                d1TestLimitedCalls=d1TestLimitedCalls+d1Stats.testLimitReached;
                d1BranchLimitedCalls=d1BranchLimitedCalls+d1Stats.branchLimitReached;
                d1RetentionLimitedCalls=d1RetentionLimitedCalls+d1Stats.retentionLimitReached;
                mergeChemicalKeys(discoveredD1Keys,d1Stats.chemicalKeys);
            end
            motifs=exchangeMotifCache(signalKey);
            if isempty(motifs),continue;end
            for xi=reshape(matchedX,1,[])
                if exchangeFound>=EXCHANGE_D1_MAX_CONNECTIONS,break;end
                gprText=strtrim(string(X.full_GPR(xi)));
                if EXCHANGE_REQUIRE_GPR && (ismissing(gprText) || gprText=="" || ...
                        string(X.gpr_status(xi))=="no_rule")
                    continue;
                end
                tr=double(X.transportRow(xi));
                for mi=1:numel(motifs)
                    if exchangeFound>=EXCHANGE_D1_MAX_CONNECTIONS,break;end
                    d=motifs(mi);
                    pathMets=pathMetsBySource{si};
                    pathRxns=pathRxnsBySource{si};
                    d1ids=[d.rff,d.rfd,d.rfr,d.rsf,d.rsd,d.rsr,d.ri];
                    [~,block0Reactions]=block0CachedPath(block0Cache,d.min);
                    motifPrefix=[block0Reactions,d1ids,pathRxns,double(T.modelIndex(tr))];
                    if reactionHasExcludedEC(model,motifPrefix,EXCLUDED_MOTIF_ECS),continue;end
                    assert(block0InputAllowed(d.min),'Ranking:UnconnectedExchange','Unconnected input cannot be selected.');
                    exchangeFound=exchangeFound+1;
                    bestCounterPathDepth=currentPathDepth;
                    selectedD1(end+1)=struct('d',d,'target',signalMet,'transportRow',tr, ...
                        'regulatorRow',0,'exchangeRow',xi, ...
                        'exchangePathMets',pathMets,'exchangePathRxns',pathRxns); %#ok<SAGROW>
                    exchangeMatchRows(end+1,:)={exchangeFound,string(metIDs(signalMet)), ...
                        join(metIDs(pathMets),' -> '),join(rxnIDs(pathRxns),'; '), ...
                        numel(pathRxns),counterID, ...
                        string(T.cytosolic(tr)),string(X.counterProduct(xi)),string(X.reaction(xi)), ...
                        double(X.importStoich(xi)),double(X.counterStoich(xi)), ...
                        string(rxnIDs(d.rsd)),string(rxnIDs(d.ri)),string(rxnIDs(d.rsr)), ...
                        string(X.full_GPR(xi))}; %#ok<SAGROW>

            fprintf(fid7x,'============================================================\n');
            fprintf(fid7x,'EXCHANGE D1 CANDIDATE %d\n',exchangeFound);
            fprintf(fid7x,'============================================================\n');
            fprintf(fid7x,'m_Sd_1 pathway source: %s | %s\n', ...
                char(metIDs(signalMet)),char(metNames(signalMet)));
            fprintf(fid7x,'m_Sd-to-counter reaction distance: %d\n',numel(pathRxns));
            fprintf(fid7x,'Exchange counter-substrate: %s | %s\n', ...
                char(counterID),char(metNames(counterMet)));
            fprintf(fid7x,'Forward intracellular path: %s\n',char(join(metIDs(pathMets),' -> ')));
            if isempty(pathRxns)
                fprintf(fid7x,'m_Sd_1 is itself the exported counter-substrate (zero-step path).\n');
            else
                fprintf(fid7x,'Path reactions: %s\n',char(join(rxnIDs(pathRxns),' -> ')));
                for prx=reshape(pathRxns,1,[])
                    printModelReaction(fid7x,model,prx,1,metIDs,metNames, ...
                        metIDs(signalMet),'m_Sd_1');
                    fprintf(fid7x,'  EC: %s\n',char(d1EC(model,prx)));
                end
            end
            fprintf(fid7x,'m_c imported cargo: %s | %s\n', ...
                char(T.cytosolic(tr)),char(T.cargo(tr)));
            fprintf(fid7x,'r_T: %s | %s\n',char(X.reaction(xi)),char(X.reactionName(xi)));
            fprintf(fid7x,'Exchange: %s out; %s in\n', ...
                char(X.counterProduct(xi)),char(T.extracellular(tr)));
            fprintf(fid7x,'Stoichiometry: counter %.8g; imported cargo %.8g\n', ...
                double(X.counterStoich(xi)),double(X.importStoich(xi)));
            fprintf(fid7x,'Transporter GPR: %s\n',char(X.full_GPR(xi)));
            fprintf(fid7x,'Resolved transporter genes: %s\n',char(X.transporter_genes(xi)));
            fprintf(fid7x,'GPR status: %s\n',char(X.gpr_status(xi)));
            fprintf(fid7x,'The counter-substrate is metabolically reachable from m_Sd_1.\n');
            fprintf(fid7x,'The model GPR is reported above. Cell-state expression is not available in this input.\n');
            fprintf(fid7x,'Transporter abundance is held fixed; no transcriptional response is modeled.\n\n');
            fprintf(fid7x,'D1 motif reactions: r_Ff=%s, r_Fd=%s, r_Fr=%s, r_Sf=%s, r_Sd=%s, r_Sr=%s, r_I=%s\n', ...
                char(rxnIDs(d.rff)),char(rxnIDs(d.rfd)),char(rxnIDs(d.rfr)), ...
                char(rxnIDs(d.rsf)),char(rxnIDs(d.rsd)),char(rxnIDs(d.rsr)),char(rxnIDs(d.ri)));
            fprintf(fid7x,'m_in: %s | %s\n',char(metIDs(d.min)),char(metNames(d.min)));
            fprintf(fid7x,'m_Ff: %s | %s\n',char(metIDs(d.mff)),char(metNames(d.mff)));
            fprintf(fid7x,'m_Fd: %s | %s\n',char(metIDs(d.mfd)),char(metNames(d.mfd)));
            fprintf(fid7x,'m_Sf: %s | %s\n',char(metIDs(d.msf)),char(metNames(d.msf)));
            fprintf(fid7x,'m_Sd: %s | %s\n',char(metIDs(signalMet)),char(metNames(signalMet)));
            fprintf(fid7x,'Equal-stoichiometry canceler coefficient: %.8g\n',d.cancelCoeff);
            fprintf(fid7x,'Status: structural match only; no flux or kinetic validation.\n\n');
        end
            end
        end
        if exchangeFound>=EXCHANGE_D1_MAX_CONNECTIONS,break;end
    end
end
if exchangeFound==0
    fprintf(fid7x,'No D1 m_Sd metabolite matched an exchange counter-substrate within depth %g and the configured search caps.\n', ...
        EXCHANGE_METABOLIC_MAX_DEPTH);
    fprintf(fid7x,'This is not proof of biological absence; review the exchange reactions or increase the path-depth/candidate caps.\n');
end
clear close7x;
% EXCHANGE_ONLY keeps regulatory protein-path candidates out of final motifs.
[selectedD1,parentDistances]=rankD1Parents(selectedD1,block0Cache.depth);
assert(all([selectedD1.exchangeRow]>0),'ExchangeOnly:ParentMode', ...
    'Only stoichiometric exchange parents may enter the complete-motif search.');
% Avoid spending the neuron budget again on the same D1/transporter chemistry.
parentChemicalKeys=containers.Map('KeyType','char','ValueType','logical');
keepParent=false(numel(selectedD1),1);
for pi=1:numel(selectedD1)
    pk=[d1ChemicalKey(selectedD1(pi)),'|',sprintf('%d,', ...
        [selectedD1(pi).exchangePathRxns,selectedD1(pi).transportRow])];
    if ~isKey(parentChemicalKeys,pk)
        parentChemicalKeys(pk)=true;keepParent(pi)=true;
    end
end
selectedD1=selectedD1(keepParent);parentDistances=parentDistances(keepParent);
parentRank=(1:numel(selectedD1))';
parentInputIDs=strings(numel(selectedD1),1);
parentTransportIDs=strings(numel(selectedD1),1);
parentCounterDistances=zeros(numel(selectedD1),1);
parentMode=strings(numel(selectedD1),1);
for pr=1:numel(selectedD1)
    parentInputIDs(pr)=metIDs(selectedD1(pr).d.min);
    parentTransportIDs(pr)=string(T.reaction(selectedD1(pr).transportRow));
    parentCounterDistances(pr)=numel(selectedD1(pr).exchangePathRxns);
    parentMode(pr)="stoichiometric_exchange";
end
ParentRanking=table(parentRank,parentCounterDistances,parentDistances,parentInputIDs, ...
    parentTransportIDs,parentMode,'VariableNames', ...
    {'ParentRank','m_SdToCounterReactionDistance','CAMPReactionDistance', ...
    'm_in_ID','TransportReactionID','TransportMode'});
writetable(ParentRanking,fullfile(outDir,'18c_D1_Parents_Ranked_By_Counter_Distance.csv'));
if isempty(exchangeMatchRows)
    exchangeMatches=table('Size',[0 15],'VariableTypes', ...
        {'double','string','string','string','double','string','string','string','string', ...
         'double','double','string','string','string','string'}, ...
        'VariableNames',{'candidate','m_Sd_metabolite','m_Sd_to_counter_metabolites', ...
        'm_Sd_to_counter_reactions','m_Sd_to_counter_reaction_distance', ...
        'exchange_counter_substrate','m_c_imported_cargo', ...
        'counter_exported_metabolite','exchange_reaction','importStoich','counterStoich', ...
        'r_Sd_reaction','r_I_reaction','r_Sr_reaction','transporter_GPR'});
else
    exchangeMatches=cell2table(exchangeMatchRows,'VariableNames', ...
        {'candidate','m_Sd_metabolite','m_Sd_to_counter_metabolites', ...
        'm_Sd_to_counter_reactions','m_Sd_to_counter_reaction_distance', ...
        'exchange_counter_substrate','m_c_imported_cargo', ...
        'counter_exported_metabolite','exchange_reaction','importStoich','counterStoich', ...
        'r_Sd_reaction','r_I_reaction','r_Sr_reaction','transporter_GPR'});
end
if height(exchangeMatches)>0
    exchangeMatches=sortrows(exchangeMatches, ...
        {'m_Sd_to_counter_reaction_distance','candidate'});
end
writetable(exchangeMatches,fullfile(outDir,'18b_Exchange_D1_Connections.csv'));
report(logID,sprintf('Stage 7B: inspected %d exchange pairs; tested %d motif steps; found %d D1/exchange matches', ...
    exchangeInspected,exchangeTests,exchangeFound));
report(logID,['Stage 7B file: ' outExchangeD1]);
%% STAGE 8: STOICHIOMETRIC counter loss and two consuming neurons.
% Neuron 1 consumes m_c; neuron 2 consumes m_Nd. No catalytic input
% is assumed. With EXCHANGE_ONLY=true, Stage 8 receives exchange candidates only.
report(logID,'Stage 8: searching CONSUMING counter -> neuron 1 -> neuron 2 chains.');
outN8=fullfile(outDir,'Step_2_Results.txt');
fid8=fopen(outN8,'w');assert(fid8>0,'Cannot create Stage 8 TXT');
close8=onCleanup(@()fclose(fid8));
fprintf(fid8,'SELECTED BLOCK 0 -> D1 -> TRANSPORT -> COUNTER -> TWO CONSUMING NEURONS\n');
fprintf(fid8,'Direction policy: STRICT +1 relative to original model.S.\n');
fprintf(fid8,'Selection: up to %d complete motifs, ranked by shortest m_Sd-to-counter path, then cAMP-to-m_in path.\n',PRINT_TOP_MOTIFS);
fprintf(fid8,'Transport mechanism filter: stoichiometric exchange only; protein-path regulation is excluded.\n');
fprintf(fid8,'Excluded reaction ECs: %s\n',char(join(EXCLUDED_MOTIF_ECS,', ')));
fprintf(fid8,'Scope: retained D1 parents and configured branch/trial limits; global optimality is not assumed.\n');
fprintf(fid8,'Equal-distance motifs may share m_in. Ties retain discovery order after both path distances.\n');
fprintf(fid8,'Block 0 depth cap=%g; expansion cap=%g; expansion limit reached=%d.\n', ...
    BLOCK0_MAX_DEPTH,BLOCK0_MAX_EXPANSIONS,block0Cache.expansionLimitReached);
fprintf(fid8,'Upstream caps: protein hops=%g; metabolite-match rows=%g; backward metabolic steps=%g; backward path rows=%g.\n', ...
    MAX_HOPS,STAGE5_MAX_MATCH_ROWS,STAGE5_MAX_HOPS,STAGE5_MAX_PATH_ROWS);
fprintf(fid8,'D1 caps: targets=%g; motifs=%g; retained/target=%g; branch=%g; reaction candidates=%g; tests/target=%g.\n', ...
    D1_MAX_TARGETS,D1_MAX_MOTIFS,D1_MAX_PER_TARGET,D1_MAX_BRANCH,D1_MAX_RXN_CAND,D1_MAX_TESTS_PER_TARGET);
fprintf(fid8,'Exchange caps: pairs=%g; connections=%g. Neuron caps: branch=%g; total tests=%g.\n', ...
    EXCHANGE_D1_MAX_TARGETS,EXCHANGE_D1_MAX_CONNECTIONS,N8_MAX_BRANCH,N8_MAX_TESTS);
fprintf(fid8,'m_Sd-to-counter path caps: depth=%g reactions; precursors/exchange pair=%g.\n', ...
    EXCHANGE_METABOLIC_MAX_DEPTH,EXCHANGE_MAX_SIGNAL_PRECURSORS);
fprintf(fid8,'D1 calls reaching trial/branch/retention limits: %d / %d / %d.\n', ...
    d1TestLimitedCalls,d1BranchLimitedCalls,d1RetentionLimitedCalls);
fprintf(fid8,'r_Nf CONSUMES m_c; r_Nf_2 CONSUMES m_Nd.\n');
fprintf(fid8,'m_Sd_1 must be produced by r_Sd_1 and lie on a shortest configured forward path to the exported counter-substrate.\n');
fprintf(fid8,'Stoichiometric exchange topology only; transporter abundance, flux, and kinetics are unverified.\n\n');
% Preindex substrate and product incidence by metabolite to avoid repeated
% full reaction scans and the former arbitrary first-80-reaction selection.
nMet8=numel(metIDs);
cons8=cell(nMet8,1);
for q8=1:numel(di)
    for mm8=reshape(ds{q8},1,[])
        cons8{mm8}(end+1)=q8;
    end
end
% A counter solution is a D1 + transporter r_T + forward r_cr consuming
% imported m_c. Count these BEFORE looking for neurons, so a counter without
% a complete downstream chain is still included in the counter total.
retainedD1Keys=containers.Map('KeyType','char','ValueType','logical');
counterKeys=containers.Map('KeyType','char','ValueType','logical');
counterCoreKeys=containers.Map('KeyType','char','ValueType','logical');
counterRows=cell(0,12);counterRowKeys=cell(0,1);
for cp=1:numel(selectedD1)
    assert(selectedD1(cp).exchangeRow>0,'ExchangeOnly:ParentMode', ...
        'A regulatory protein-path parent reached the complete-motif search.');
    parent=selectedD1(cp);dd=parent.d;tt=parent.transportRow;
    retainedD1Keys(d1ChemicalKey(parent))=true;
    mm=find(metIDs==string(T.cytosolic(tt)),1);
    if isempty(mm) || comp(mm)~="c",continue;end
    candidates=cons8{mm};candidates=candidates(1:min(end,N8_MAX_BRANCH));
    d1reactions=[dd.rff,dd.rfd,dd.rfr,dd.rsf,dd.rsd,dd.rsr,dd.ri];
    for cr=reshape(candidates,1,[])
        if ismember(di(cr),d1reactions),continue;end
        transportReaction=double(T.modelIndex(tt));
        if di(cr)==transportReaction,continue;end
        ck=counterChemicalKey(parent,transportReaction,di(cr),mm);
        if isKey(counterKeys,ck),continue;end
        counterKeys(ck)=true;
        counterCoreKeys(sprintf('%d,',[transportReaction,di(cr),mm]))=true;
        mode="stoichiometric_exchange";
        counterRows(end+1,:)={size(counterRows,1)+1,cp, ...
            numel(parent.exchangePathRxns),block0Cache.depth(dd.min), ...
            metIDs(dd.min),metIDs(parent.target),metIDs(mm), ...
            rxnIDs(transportReaction),rxnIDs(di(cr)),mode,d1ChemicalKey(parent), ...
            join(rxnIDs(parent.exchangePathRxns),'; ')}; %#ok<SAGROW>
        counterRowKeys{end+1,1}=ck; %#ok<SAGROW>
    end
end
report(logID,sprintf('Counters: %d distinct D1-counter solutions; %d distinct transport/loss cores before neuron search.', ...
    counterKeys.Count,counterCoreKeys.Count));

found8=0; tests8=0; examined8=0;
seen8=containers.Map('KeyType','char','ValueType','logical');
countersWithComplete=containers.Map('KeyType','char','ValueType','logical');
completeExchange8=0;completeRegulatory8=0;
best8=struct('parentIndex',{},'mc',{},'nf',{},'nd',{},'nf2',{}, ...
    'mout',{},'chosen',{},'counterDistance',{},'distance',{},'discoveryIndex',{});
rankingRows8=cell(0,12);
for ii=1:numel(selectedD1)
    if tests8>=N8_MAX_TESTS,break;end
    sel=selectedD1(ii);tr=sel.transportRow;
    mc=find(metIDs==string(T.cytosolic(tr)),1);
    if isempty(mc) || comp(mc)~="c",continue;end
    examined8=examined8+1;
    % r_cr and r_Nf are distinct reactions consuming m_c.
    inputRx=cons8{mc};
    if numel(inputRx)<2,continue;end
    inputRx=inputRx(1:min(end,N8_MAX_BRANCH));
    for i1=reshape(inputRx,1,[])
        if tests8>=N8_MAX_TESTS,break;end
        nfCandidates=do{i1};
        nfCandidates=nfCandidates(arrayfun(@(v)d1SignalOK(v,mc,metIDs,comp,currency),nfCandidates));
        nfCandidates=nfCandidates(1:min(end,N8_MAX_BRANCH));
        for nf=reshape(nfCandidates,1,[])
            secondRx=cons8{nf};secondRx=secondRx(1:min(end,N8_MAX_BRANCH));
            for i2=reshape(secondRx,1,[])
                if di(i2)==di(i1),continue;end
                ndCandidates=do{i2};
                ndCandidates=ndCandidates(arrayfun(@(v)d1SignalOK(v,nf,metIDs,comp,currency) && v~=mc,ndCandidates));
                ndCandidates=ndCandidates(1:min(end,N8_MAX_BRANCH));
                for nd=reshape(ndCandidates,1,[])
                    % r_Nr and r_Nf_2 are distinct reactions consuming m_Nd.
                    ndRx=cons8{nd};ndRx=ndRx(1:min(end,N8_MAX_BRANCH));
                    if numel(ndRx)<2,continue;end
                    for i4=reshape(ndRx,1,[])
                        if ismember(di(i4),di([i1,i2])),continue;end
                        nf2Candidates=do{i4};
                        nf2Candidates=nf2Candidates(arrayfun(@(v)d1SignalOK(v,nd,metIDs,comp,currency) && ~ismember(v,[mc,nf]),nf2Candidates));
                        nf2Candidates=nf2Candidates(1:min(end,N8_MAX_BRANCH));
                        for nf2=reshape(nf2Candidates,1,[])
                            fifthRx=cons8{nf2};fifthRx=fifthRx(1:min(end,N8_MAX_BRANCH));
                            for i5=reshape(fifthRx,1,[])
                                if ismember(di(i5),di([i1,i2,i4])),continue;end
                                outCandidates=do{i5};
                                outCandidates=outCandidates(arrayfun(@(v)d1SignalOK(v,nf2,metIDs,comp,currency) && ~ismember(v,[mc,nf,nd]),outCandidates));
                                outCandidates=outCandidates(1:min(end,N8_MAX_BRANCH));
                                for mout=reshape(outCandidates,1,[])
                                    endRx=cons8{mout};endRx=endRx(1:min(end,N8_MAX_BRANCH));
                                    if isempty(endRx),continue;end
                                    % Find independent loss reactions for m_c and m_Nd.
                                    for irc=reshape(inputRx,1,[])
                                        if di(irc)==di(i1),continue;end
                                        for inr=reshape(ndRx,1,[])
                                            if di(inr)==di(i4),continue;end
                                            for inr2=reshape(endRx,1,[])
                                                if tests8>=N8_MAX_TESTS,break;end
                                                tests8=tests8+1;
                                                chosen=[irc,i1,i2,inr,i4,i5,inr2];
                                                if numel(unique(di(chosen)))~=7,continue;end
                                                % Do not use D1 reactions as downstream neuron reactions.
                                                d=sel.d;
                                                d1ids=[d.rff,d.rfd,d.rfr,d.rsf,d.rsd,d.rsr,d.ri];
                                                if any(ismember(di(chosen),d1ids)),continue;end
                                                [~,block0Reactions]=block0CachedPath(block0Cache,d.min);
                                                allMotifReactions=[block0Reactions,d1ids,sel.exchangePathRxns, ...
                                                    double(T.modelIndex(tr)),reshape(di(chosen),1,[])];
                                                if reactionHasExcludedEC(model,allMotifReactions,EXCLUDED_MOTIF_ECS)
                                                    continue;
                                                end
                                                assert(block0InputAllowed(d.min),'Ranking:UnconnectedCircuit', ...
                                                    'Cannot accept a complete motif without a Block 0 path.');
                                                chemicalKey=sprintf('%d,', ...
                                                    [d.min,d.mff,d.mfd,d.msf,sel.target,mc,nf,nd,nf2,mout, ...
                                                     d1ids,sel.exchangePathRxns,double(T.modelIndex(tr)),reshape(di(chosen),1,[])]);
                                                if isKey(seen8,chemicalKey),continue;end
                                                seen8(chemicalKey)=true;
                                                found8=found8+1;
                                                distance=block0Cache.depth(d.min);
                                                counterKey=counterChemicalKey(sel,double(T.modelIndex(tr)),di(irc),mc);
                                                countersWithComplete(counterKey)=true;
                                                if sel.exchangeRow>0
                                                    completeExchange8=completeExchange8+1;
                                                else
                                                    completeRegulatory8=completeRegulatory8+1;
                                                end
                                                item8=struct('parentIndex',ii,'mc',mc,'nf',nf,'nd',nd, ...
                                                    'nf2',nf2,'mout',mout,'chosen',chosen, ...
                                                    'counterDistance',numel(sel.exchangePathRxns), ...
                                                    'distance',distance,'discoveryIndex',found8);
                                                best8(end+1)=item8; %#ok<SAGROW>
                                                best8=keepTopCompleteMotifs(best8,PRINT_TOP_MOTIFS);
                                            end
                                            if tests8>=N8_MAX_TESTS,break;end
                                        end
                                        if tests8>=N8_MAX_TESTS,break;end
                                    end
                                    if tests8>=N8_MAX_TESTS,break;end
                                end
                                if tests8>=N8_MAX_TESTS,break;end
                            end
                            if tests8>=N8_MAX_TESTS,break;end
                        end
                        if tests8>=N8_MAX_TESTS,break;end
                    end
                    if tests8>=N8_MAX_TESTS,break;end
                end
                if tests8>=N8_MAX_TESTS,break;end
            end
            if tests8>=N8_MAX_TESTS,break;end
        end
    end
    if mod(ii,10)==0
        report(logID,sprintf('Stage 8 progress: %d D1 motifs; %d downstream candidates; %d trials',ii,found8,tests8));
    end
end
printed8=numel(best8);
assert(found8==seen8.Count && printed8==min(found8,PRINT_TOP_MOTIFS), ...
    'Counts:Complete','Found/printed complete-motif counts are inconsistent.');
assert(completeExchange8+completeRegulatory8==found8, ...
    'Counts:Mechanism','Transport-mode counts do not sum to the complete count.');
assert(counterKeys.Count==size(counterRows,1), ...
    'Counts:Counter','Counter count and counter table differ.');
completedCounterNames=feval('keys',countersWithComplete); % Stage 5 also has a variable named keys.
for completedIndex=1:numel(completedCounterNames)
    assert(isKey(counterKeys,completedCounterNames{completedIndex}), ...
        'Counts:CounterExtension','A completed counter is missing from the counter catalog.');
end
fprintf(fid8,'\n============================================================\n');
fprintf(fid8,'SEARCH COUNTS: FOUND VERSUS PRINTED\n');
fprintf(fid8,'============================================================\n');
fprintf(fid8,'Distinct D1 motifs discovered in bounded searches: %d\n',discoveredD1Keys.Count);
fprintf(fid8,'Distinct D1 motifs retained for downstream search: %d\n',retainedD1Keys.Count);
fprintf(fid8,'D1/transporter parent connections retained: %d\n',numel(selectedD1));
fprintf(fid8,'Distinct counter cores (r_T, m_c, r_cr): %d\n',counterCoreKeys.Count);
fprintf(fid8,'Distinct D1-counter solutions in retained parent pool: %d\n',counterKeys.Count);
fprintf(fid8,'D1-counter solutions extended to both neurons: %d\n',countersWithComplete.Count);
fprintf(fid8,'Distinct complete motifs found: %d\n',found8);
fprintf(fid8,'  Exchange mode: %d; regulatory-hypothesis mode: %d\n',completeExchange8,completeRegulatory8);
fprintf(fid8,'Complete motifs printed: %d (requested maximum %d)\n',printed8,PRINT_TOP_MOTIFS);
fprintf(fid8,'Count definitions: chemical structures are deduplicated; alternate protein annotations do not add motifs.\n');
fprintf(fid8,'One shortest Block 0 path is used per m_in; alternate connectors do not multiply motif counts.\n');
fprintf(fid8,'Counter total includes counters with no two-neuron solution; complete count requires both neurons.\n');
fprintf(fid8,'D1 parents examined for neurons: %d/%d; trials=%d; neuron trial cap reached=%d.\n', ...
    examined8,numel(selectedD1),tests8,tests8>=N8_MAX_TESTS);
fprintf(fid8,'Counts describe the configured search only; finite branch/trial/retention caps can hide additional solutions.\n');
fprintf(fid8,'Selecting ten for printing does NOT stop discovery/counting at ten.\n\n');

for printedRank8=1:printed8
    item8=best8(printedRank8);
    ii=item8.parentIndex;sel=selectedD1(ii);d=sel.d;tr=sel.transportRow;
    mc=item8.mc;nf=item8.nf;nd=item8.nd;nf2=item8.nf2;mout=item8.mout;
    chosen=item8.chosen;distance=item8.distance;
    [~,block0Reactions]=block0CachedPath(block0Cache,d.min);
    mode="stoichiometric_exchange";
    counterDistance=numel(sel.exchangePathRxns);
    rankingRows8(end+1,:)={printedRank8,counterDistance,distance, ...
        metIDs(d.min),metIDs(sel.target), ...
        metIDs(mc),metIDs(mout),string(T.reaction(tr)),mode,ii, ...
        join(rxnIDs(block0Reactions),' -> '), ...
        join(rxnIDs(di(chosen)),'; ')}; %#ok<SAGROW>
    fprintf(fid8,'============================================================\nCANDIDATE %d\n',printedRank8);
    fprintf(fid8,'m_Sd-to-counter distance: %d reaction(s) | cAMP-to-m_in distance: %d reaction(s) | D1 parent rank: %d\n', ...
        counterDistance,distance,ii);
    fprintf(fid8,'============================================================\n');
    % Named signals: same metabolite ID is labelled consistently in every equation.
    signalIDs=metIDs([d.min,d.mff,d.mfd,d.msf,sel.target,mc,nf,nd,nf2,mout]);
    signalNames=["m_in","m_Ff","m_Fd","m_Sf","m_Sd","m_c", ...
        "m_Nf_1","m_Nd_1","m_Nf_2","m_out"];
    fprintf(fid8,'SIGNAL METABOLITE MAP\n');
    for kk=1:numel(signalIDs)
        ix=find(metIDs==signalIDs(kk),1);
        fprintf(fid8,'  %-9s: %s | %s\n',char(signalNames(kk)), ...
            char(signalIDs(kk)),char(metNames(ix)));
    end
    fprintf(fid8,'  r_Nr_1 and r_Nr_2 are removal REACTIONS; there is no m_Nr state in the supplied dynamics.\n\n');
    fprintf(fid8,'============================================================\n');
    fprintf(fid8,'BLOCK 0: cAMP TO m_in\n');
    fprintf(fid8,'============================================================\n');
    printBlock0(fid8,model,metIDs,metNames,d.min,signalIDs,signalNames,block0Cache);
    fprintf(fid8,'\n============================================================\n');
    fprintf(fid8,'DIFFERENTIATOR (SEVEN REACTIONS)\n');
    fprintf(fid8,'============================================================\n');
    labelsD=["r_Ff","r_Fd","r_Fr","r_Sf","r_Sd","r_Sr","r_I"];
    idsD=[d.rff,d.rfd,d.rfr,d.rsf,d.rsd,d.rsr,d.ri];
    orsD=[d.off,d.ofd,d.ofr,d.osf,d.osd,d.osr,d.oi];
    for kk=1:7
        fprintf(fid8,'%s\n',char(labelsD(kk)));
        printModelReaction(fid8,model,idsD(kk),orsD(kk),metIDs,metNames,"","",signalIDs,signalNames);
        fprintf(fid8,'  EC: %s\n',char(d1EC(model,idsD(kk))));
        [ruleD,genesD,~]=decodeGPR(gr(idsD(kk)),rr(idsD(kk)),genes);
        fprintf(fid8,'  GPR: %s\n  Resolved genes: %s\n\n',char(ruleD),char(genesD));
    end
    fprintf(fid8,'  CHECK: r_Sd produces m_Sd; r_I consumes m_Fd and m_Sd equally (%.8g each).\n',d.cancelCoeff);
    if sel.exchangeRow>0
        ex=sel.exchangeRow;
        fprintf(fid8,'\nTRANSPORT CONTROL: STOICHIOMETRIC EXCHANGE\n');
        fprintf(fid8,'  m_Sd is the source metabolite for the intracellular path to the exported counter-substrate.\n');
        fprintf(fid8,'  m_Sd-to-counter path: %s\n',char(join(metIDs(sel.exchangePathMets),' -> ')));
        if isempty(sel.exchangePathRxns)
            fprintf(fid8,'  Path length: zero; m_Sd is itself the exchanged counter-substrate.\n');
        else
            fprintf(fid8,'  Path reactions: %s\n',char(join(rxnIDs(sel.exchangePathRxns),' -> ')));
            for pathIndex=reshape(sel.exchangePathRxns,1,[])
                printModelReaction(fid8,model,pathIndex,1,metIDs,metNames, ...
                    metIDs(sel.target),'m_Sd');
                fprintf(fid8,'  EC: %s\n',char(d1EC(model,pathIndex)));
            end
        end
        fprintf(fid8,'  Imported cargo: %s -> %s\n', ...
            char(T.extracellular(tr)),char(T.cytosolic(tr)));
        fprintf(fid8,'  Counter-substrate: %s -> %s\n', ...
            char(X.counterSubstrate(ex)),char(X.counterProduct(ex)));
        fprintf(fid8,'  Stoichiometry: import %.8g; counter export %.8g\n', ...
            double(X.importStoich(ex)),double(X.counterStoich(ex)));
        fprintf(fid8,'  Transporter GPR: %s\n',char(X.full_GPR(ex)));
        fprintf(fid8,'  Transporter abundance is treated as pre-existing; no protein-path or transcription step is modeled.\n');
    else
        error('ExchangeOnly:ParentMode','A regulatory protein-path motif reached final output.');
    end
    fprintf(fid8,'\nMEMBRANE TRANSPORT: r_T (production of m_c)\n');
    assert(string(T.direction(tr))=="forward", ...
        'ForwardOnly:TransportDirection','A transporter has non-forward direction.');
    ot=1;
    printModelReaction(fid8,model,double(T.modelIndex(tr)),ot,metIDs,metNames,"","",signalIDs,signalNames);
    fprintf(fid8,'\nCOUNTER LOSS AND TWO CONSUMING NEURONS\n');
    labels8=["r_cr","r_Nf_1","r_Nd_1","r_Nr_1","r_Nf_2","r_Nd_2","r_Nr_2"];
    for kk=1:7
        j=chosen(kk);
        fprintf(fid8,'%s\n',char(labels8(kk)));
        printModelReaction(fid8,model,di(j),dm(j),metIDs,metNames,"","",signalIDs,signalNames);
        fprintf(fid8,'  EC: %s\n\n',char(d1EC(model,di(j))));
    end
    fprintf(fid8,'MODEL REQUIREMENTS: r_Nf consumes m_c; r_Nf_2 consumes m_Nd.\n');
    fprintf(fid8,'WARNING: no flux feasibility, kinetic differentiation, or validated transporter regulation.\n\n');

end
if printed8==0
    fprintf(fid8,'No cAMP-connected complete two-neuron motif passed the exchange-only filter within configured search limits.\n');
    fprintf(fid8,'No regulatory protein-path fallback was used. Check 18b_Exchange_D1_Connections.txt for exchange matches.\n');
    fprintf(fid8,'This is NOT proof of biological absence.\n');
end
rankNames8={'Rank','m_SdToCounterReactionDistance','CAMPReactionDistance', ...
    'm_in_ID','m_Sd_ID','m_c_ID', ...
    'm_out_ID','TransportReactionID','TransportMode','ParentRank', ...
    'Block0ReactionIDs','CounterNeuronReactionIDs'};
if isempty(rankingRows8)
    MotifRanking=table('Size',[0 12],'VariableTypes', ...
        {'double','double','double','string','string','string','string','string','string','double','string','string'}, ...
        'VariableNames',rankNames8);
else
    MotifRanking=cell2table(rankingRows8,'VariableNames',rankNames8);
    if height(MotifRanking)>1
        counterDistanceOrder=MotifRanking.m_SdToCounterReactionDistance;
        campDistanceOrder=MotifRanking.CAMPReactionDistance;
        ordered=counterDistanceOrder(1:end-1)<counterDistanceOrder(2:end) | ...
            (counterDistanceOrder(1:end-1)==counterDistanceOrder(2:end) & ...
            campDistanceOrder(1:end-1)<=campDistanceOrder(2:end));
        assert(all(ordered),'Ranking:OutputOrder', ...
            'Final motif ranking must prioritize m_Sd-to-counter distance, then cAMP distance.');
    end
end
writetable(MotifRanking,fullfile(outDir,'19_Selected_Motifs_Counter_Distance_Ranking.csv'));
counterNames={'CounterIndex','ParentRank','m_SdToCounterReactionDistance', ...
    'CAMPReactionDistance','m_in_ID','m_Sd_ID', ...
    'm_c_ID','TransportReactionID','CounterLossReactionID','TransportMode','D1ChemicalKey', ...
    'ExchangePathReactionIDs'};
if isempty(counterRows)
    CounterSolutions=table('Size',[0 12],'VariableTypes', ...
        {'double','double','double','double','string','string','string','string','string','string','string','string'}, ...
        'VariableNames',counterNames);
else
    CounterSolutions=cell2table(counterRows,'VariableNames',counterNames);
end
hasComplete=false(size(counterRows,1),1);
for ckIndex=1:numel(counterRowKeys)
    hasComplete(ckIndex)=isKey(countersWithComplete,counterRowKeys{ckIndex});
end
CounterSolutions.ExtendedToTwoNeurons=hasComplete;
writetable(CounterSolutions,fullfile(outDir,'18d_All_Counter_Solutions.csv'));
clear close8;
report(logID,sprintf('Stage 8: found %d distinct complete motifs; printed %d; %d distinct D1-counter solutions; %d trials.', ...
    found8,printed8,counterKeys.Count,tests8));
report(logID,['Stage 8 file: ' outN8]);

report(logID,['Finished: ' datestr(now,31)]);
clear logGuard;
fprintf('\nFinished. Outputs: %s\n',outDir);

%% LOCAL HELPERS - no external custom MATLAB files.
function report(fid,msg)
    fprintf('%s\n',char(msg));fprintf(fid,'%s\n',char(msg));
end
function path=findNamed(root,name)
    path=which(name);if ~isempty(path) && isfile(path),return;end
    p=fullfile(pwd,name);if isfile(p),path=p;return;end
    p=fullfile(root,name);if isfile(p),path=p;return;end
    path='';
    % MATLAB recursive dir ** supported on modern versions.
    try
        d=dir(fullfile(root,'**',name));
        if ~isempty(d),path=fullfile(d(1).folder,d(1).name);end
    catch
    end
end
function [base,comp]=splitCompartments(mets)
    base=strings(size(mets));comp=strings(size(mets));
    for ii=1:numel(mets)
        x=char(mets(ii));
        t=regexp(x,'^(.*)\[([a-zA-Z])\]$','tokens','once');
        if isempty(t),t=regexp(x,'^(.*)_([ecm])$','tokens','once');end
        if ~isempty(t),base(ii)=string(t{1});comp(ii)=string(t{2});
        else,base(ii)=string(x);comp(ii)="";end
    end
end
function [gpr,glist,status]=decodeGPR(gr,rule,genes)
    gpr="";glist="";status="no_rule";
    gr=strtrim(string(gr));rule=strtrim(string(rule));
    if ~ismissing(gr) && gr~="" && lower(gr)~="nan"
        gpr=gr;status="grRules";
        % grRules may contain gene IDs rather than x(i). Keep original text.
        % Recover gene membership from exact gene tokens where possible.
        present=strings(0,1);
        for j=1:numel(genes)
            if ~isempty(regexp(char(gr),['(?<![A-Za-z0-9_.-])' regexptranslate('escape',char(genes(j))) '(?![A-Za-z0-9_.-])'],'once'))
                present(end+1,1)=genes(j); %#ok<AGROW>
            end
        end
        glist=join(present,';');
    end
    if (gpr=="" || glist=="") && rule~="" && ~ismissing(rule) && lower(rule)~="nan"
        tok=regexp(char(rule),'x\((\d+)\)','tokens');
        if ~isempty(tok)
            idx=unique(cellfun(@(v)str2double(v{1}),tok),'stable');
            if any(idx<1 | idx>numel(genes) | ~isfinite(idx))
                status="invalid_gene_index";return;
            end
            glist=join(genes(idx),';');
            gpr=rule;
            [st,en,tokens]=regexp(char(rule),'x\((\d+)\)','start','end','tokens');
            chars=char(rule);
            for q=numel(st):-1:1
                ix=str2double(tokens{q}{1});
                chars=[chars(1:st(q)-1) '(' char(genes(ix)) ')' chars(en(q)+1:end)]; %#ok<AGROW>
            end
            gpr=string(chars);status="recovered_rules";
        end
    end
end
function b=asBool(v)
    if isnumeric(v)||islogical(v),b=logical(v);return;end
    b=ismember(lower(strtrim(string(v))),["1","true","t","yes"]);
end

function out=urlencodeLocal(in)
    bytes=unicode2native(char(in),'UTF-8');
    out='';
    for n=1:numel(bytes)
        c=char(bytes(n));
        if isstrprop(c,'alphanum') || any(c=='-_.~')
            out=[out c]; %#ok<AGROW>
        else
            out=[out sprintf('%%%02X',bytes(n))]; %#ok<AGROW>
        end
    end
end

function result=webreadRetry(url,opts,nTries,pauseSeconds)
% Retry transient HTTP/network failures; do not cache failed responses.
    for attempt=1:nTries
        try
            result=webread(url,opts);
            return;
        catch ME
            if attempt==nTries,rethrow(ME);end
            pause(pauseSeconds*2^(attempt-1));
        end
    end
end

function [name,status]=cleanBrendaName(raw)
% Conservative cleaning: remove BRENDA comments and citation fragments.
    x=char(strtrim(string(raw)));
    x=regexprep(x,'^\s*#[^#]*#\s*','');
    x=regexprep(x,'\s*\([^)]*#.*$','');
    x=regexprep(x,'\s*[<{].*$','');
    x=regexprep(x,'\s*\([^)]*(?:activity|inhibit|stimul|protein|enzyme)[^)]*\).*$','', 'ignorecase');
    x=strtrim(x);
    name=string(x);status="usable";
    if strlength(name)<2 || strlength(name)>120 || ...
       contains(name,"#") || contains(name,"(") || contains(name,")") || ...
       any(strcmpi(name,["more","less","unknown","none","protein","n/a","nan"]))
        status="reject_uncertain_or_noncompound";name="";
    end
end

function key=normalizeCompound(name)
    key=lower(strtrim(string(name)));
    key=regexprep(key,'\s+',' ');
    key=regexprep(key,'[‐‑–—]','-');
    if ismissing(key),key="";end
end

function cache=buildBlock0Cache(model,metIDs,comp,ub,sourceID,maxDepth,maxExpansions)
% Explore ALL structurally reachable inputs once, using +1 reactions only.
% A reaction is an edge from each substrate to each product. Other substrates
% are not assumed to be produced by this path; flux feasibility is unverified.
    assert(isnumeric(maxDepth) && isscalar(maxDepth) && ...
        ~isnan(maxDepth) && maxDepth>=0 && ...
        (isinf(maxDepth) || maxDepth==floor(maxDepth)), ...
        'Block0:Depth','BLOCK0_MAX_DEPTH must be a nonnegative integer or Inf.');
    assert(isnumeric(maxExpansions) && isscalar(maxExpansions) && ...
        ~isnan(maxExpansions) && maxExpansions>=0 && ...
        (isinf(maxExpansions) || maxExpansions==floor(maxExpansions)), ...
        'Block0:ExpansionCap','BLOCK0_MAX_EXPANSIONS must be a nonnegative integer or Inf.');
    start=find(metIDs==string(sourceID));
    assert(numel(start)==1,'Block0:Anchor', ...
        'Expected one Block 0 anchor with exact model ID %s.',char(string(sourceID)));
    assert(comp(start)=="c",'Block0:Compartment', ...
        'The Block 0 anchor must be cytosolic.');
    nm=numel(metIDs);
    cache=struct('start',start,'sourceID',string(sourceID), ...
        'maxDepth',maxDepth,'maxExpansions',maxExpansions, ...
        'reachable',false(nm,1),'depth',inf(nm,1), ...
        'parentMet',zeros(nm,1),'parentRx',zeros(nm,1), ...
        'expansions',0,'expansionLimitReached',false);
    cache.reachable(start)=true;cache.depth(start)=0;
    S=model.S;
    tol=1e-12;
    eligible=ub(:)>0;
    eligible=eligible & full(any(S < -tol,1))' & ...
        full(any(S > tol,1))';
    % No extracellular or organelle participant may occur in Block 0.
    outside=comp~="c";
    eligible=eligible & ~full(any(abs(S(outside,:))>tol,1))';
    queue=zeros(nm,1);queue(1)=start;head=1;tail=1;
    while head<=tail
        current=queue(head);head=head+1;
        if cache.depth(current)>=maxDepth,continue;end
        reactions=find(full(S(current,:)) < -tol & eligible');
        for j=reshape(reactions,1,[])
            if cache.expansions>=maxExpansions
                cache.expansionLimitReached=true;
                return;
            end
            cache.expansions=cache.expansions+1;
            products=find(S(:,j)>tol);
            for next=reshape(products,1,[])
                if cache.reachable(next),continue;end
                cache.reachable(next)=true;
                cache.parentMet(next)=current;
                cache.parentRx(next)=j;
                cache.depth(next)=cache.depth(current)+1;
                tail=tail+1;queue(tail)=next;
            end
        end
    end
end

function [mets,rx]=block0CachedPath(cache,target)
% Reconstruct the shortest verified forward path from the BFS tree.
    assert(cache.reachable(target),'Block0:UnconnectedInput', ...
        'Requested m_in is not reachable in the cached Block 0 search.');
    mets=target;rx=zeros(1,0);current=target;
    while current~=cache.start
        r=cache.parentRx(current);
        parent=cache.parentMet(current);
        assert(r>0 && parent>0,'Block0:BrokenPath', ...
            'Cached Block 0 predecessor chain is incomplete.');
        rx=[r,rx];mets=[parent,mets];current=parent; %#ok<AGROW>
    end
    assert(numel(rx)==cache.depth(target),'Block0:PathLength', ...
        'Cached path depth and reaction count differ.');
end

function printBlock0(fid,model,metIDs,metNames,target,signalIDs,signalNames,cache)
% Print precisely the path used to accept this motif; never rerun a search.
    [mets,rx]=block0CachedPath(cache,target);
    fprintf(fid,'  Anchor: %s | %s\n', ...
        char(metIDs(cache.start)),char(metNames(cache.start)));
    fprintf(fid,'  Connection status: CONNECTED (+1 only)\n');
    fprintf(fid,'  Search depth limit: %g; expansion cap reached: %d\n', ...
        cache.maxDepth,cache.expansionLimitReached);
    fprintf(fid,'  Path: %s (%d reaction(s))\n', ...
        char(join(metIDs(mets),' -> ')),numel(rx));
    if isempty(rx)
        fprintf(fid,'  m_in is cAMP; zero reactions explicitly allowed.\n');
        return;
    end
    for k=1:numel(rx)
        assert(model.S(mets(k),rx(k))<0 && ...
            model.S(mets(k+1),rx(k))>0,'Block0:InvalidForwardEdge', ...
            'A cached reaction is inconsistent with +1 model.S direction.');
        fprintf(fid,'r_CAMP_%d\n',k);
        printModelReaction(fid,model,rx(k),1,metIDs,metNames,"","",signalIDs,signalNames);
        fprintf(fid,'  EC: %s\n',char(d1EC(model,rx(k))));
    end
    fprintf(fid,'  WARNING: shortest structural path; additional substrates and flux feasibility unverified.\n');
end

function selfTestBlock0Connectivity()
% Small synthetic regression cases, executed without network or COBRA.
    fixture.S=sparse([-1 -1 0 0 0 0 0; ...
        1 0 -1 0 0 0 0; ...
        0 1 0 0 0 0 0; ...
        0 0 1 1 -1 0 0; ...
        0 0 0 -1 0 0 0; ...
        0 0 0 0 1 -1 0; ...
        0 0 0 0 0 1 -1; ...
        0 0 0 0 0 0 1]);
    ids=["camp[c]";"a[c]";"direct[c]";"b[c]"; ...
        "reverse_only[c]";"outside[e]";"via_outside[c]";"blocked[c]"];
    comps=["c";"c";"c";"c";"c";"e";"c";"c"];
    ub=[1000;1000;1000;1000;1000;1000;0];
    fixture.S(1,8)=-1;fixture.S(8,8)=1;ub(8)=0;
    c=buildBlock0Cache(fixture,ids,comps,ub,'camp[c]',Inf,Inf);
    assert(c.reachable(4) && c.depth(4)==2, ...
        'Block0 self-test: two-step forward path missing.');
    assert(c.reachable(3) && c.depth(3)==1, ...
        'Block0 self-test: direct path missing.');
    assert(~any(c.reachable(5:8)), ...
        'Block0 self-test: reverse, extracellular, or blocked path accepted.');
    [m,r]=block0CachedPath(c,4);
    assert(isequal(m,[1 2 4]) && isequal(r,[1 3]), ...
        'Block0 self-test: path reconstruction failed.');
    c1=buildBlock0Cache(fixture,ids,comps,ub,'camp[c]',1,Inf);
    assert(~c1.reachable(4) && c1.reachable(2) && c1.reachable(3), ...
        'Block0 self-test: depth cap failed.');
    c2=buildBlock0Cache(fixture,ids,comps,ub,'camp[c]',Inf,0);
    assert(c2.expansionLimitReached && nnz(c2.reachable)==1, ...
        'Block0 self-test: expansion truncation not flagged.');
    fixture.S(1,9)=-1;fixture.S(4,9)=1;ub(9)=1000;
    c3=buildBlock0Cache(fixture,ids,comps,ub,'camp[c]',Inf,Inf);
    [m,r]=block0CachedPath(c3,4);
    assert(c3.depth(4)==1 && isequal(m,[1 4]) && isequal(r,9), ...
        'Block0 self-test: shortest path was not preferred.');
    fprintf('PASS: Block 0 connectivity, +1 direction, compartments, limits and shortest paths.\n');
end

function printModelReaction(fid,model,idx,orient,metIDs,metNames,markedID,markedLabel,varargin)
% Print original stoichiometry; never relabel or reverse a stored reaction.
    assert(isscalar(orient) && orient==1, ...
        'ForwardOnly:PrintedDirection', ...
        'Refusing to print reaction %s with direction other than +1.', ...
        char(string(model.rxns(idx))));
    if isfield(model,'ub')
        assert(model.ub(idx)>0,'ForwardOnly:ReactionBound', ...
            'Reaction %s does not permit positive flux.', ...
            char(string(model.rxns(idx))));
    end
    rxn=string(model.rxns(idx));
    rxnName=rxn;
    if isfield(model,'rxnNames') && numel(model.rxnNames)>=idx
        rxnName=string(model.rxnNames(idx));
    end
    signalIDs=strings(0,1);signalNames=strings(0,1);
    if numel(varargin)>=2
        signalIDs=string(varargin{1});signalNames=string(varargin{2});
    end
    coeff=full(model.S(:,idx))*orient;
    left=find(coeff < -1e-12);right=find(coeff > 1e-12);
    lhs=strings(numel(left),1);rhs=strings(numel(right),1);
    for q=1:numel(left)
        lhs(q)=reactionTerm(-coeff(left(q)),left(q),metIDs,metNames,markedID,markedLabel,signalIDs,signalNames);
    end
    for q=1:numel(right)
        rhs(q)=reactionTerm(coeff(right(q)),right(q),metIDs,metNames,markedID,markedLabel,signalIDs,signalNames);
    end
    fprintf(fid,'  ID: %s\n',char(rxn));
    fprintf(fid,'  Name: %s\n',char(rxnName));
    fprintf(fid,'  Direction: %+d relative to model.S\n',orient);
    fprintf(fid,'  Equation: %s  ->  %s\n',char(join(lhs,' + ')),char(join(rhs,' + ')));
    fprintf(fid,'\n');
end
function term=reactionTerm(coef,idx,metIDs,metNames,markedID,markedLabel,signalIDs,signalNames)
    label=metNames(idx);
    if label=="" || ismissing(label),label=metIDs(idx);end
    tags=strings(0,1);
    if nargin>=8 && ~isempty(signalIDs)
        tags=signalNames(signalIDs==metIDs(idx));
    end
    if isempty(tags) && metIDs(idx)==markedID && strlength(string(markedLabel))>0
        tags=string(markedLabel);
    end
    if ~isempty(tags),label=label+" ("+join(unique(tags,'stable'),', ')+")";end
    if abs(coef-1)>1e-10,label=string(sprintf('%.8g ',coef))+label;end
    term=label+" ["+metIDs(idx)+"]";
end

function [rx,sgn,products,substrates]=d1BuildOriented(model,comp,lb,ub,skipTransport)
    rx=zeros(0,1);sgn=zeros(0,1);products=cell(0,1);substrates=cell(0,1);
    for j=1:numel(model.rxns)
        ids=find(model.S(:,j));
        if numel(ids)<2,continue;end
        if skipTransport && (any(comp(ids)~="c")),continue;end
        if ub(j)>0
            coef=full(model.S(ids,j));
            a=ids(coef<0);b=ids(coef>0);
            if ~isempty(a) && ~isempty(b)
                rx(end+1,1)=j;sgn(end+1,1)=1; %#ok<AGROW>
                products{end+1,1}=b;substrates{end+1,1}=a; %#ok<AGROW>
            end
        end
    end
end

function [sources,metPaths,rxnPaths]=findForwardPrecursorPaths(target,rx,products,substrates,metIDs,comp,currency,maxDepth)
% Find cytosolic metabolites that can produce target through forward reactions.
% The returned paths run from each source metabolite toward target.
    nMet=numel(metIDs);
    allowed=comp(:)=="c";
    for mm=1:nMet
        base=regexprep(char(metIDs(mm)),'\[c\]$','');
        if ismember(lower(string(base)),currency),allowed(mm)=false;end
    end
    producers=cell(nMet,1);
    for q=1:numel(rx)
        for product=reshape(products{q},1,[])
            producers{product}(end+1)=q; %#ok<AGROW>
        end
    end
    visited=false(nMet,1);depth=inf(nMet,1);
    nextMet=zeros(nMet,1);nextRx=zeros(nMet,1);
    visited(target)=true;depth(target)=0;
    queue=zeros(nMet,1);queue(1)=target;head=1;tail=1;
    while head<=tail
        downstream=queue(head);head=head+1;
        if depth(downstream)>=maxDepth,continue;end
        for q=reshape(producers{downstream},1,[])
            for upstream=reshape(substrates{q},1,[])
                if ~allowed(upstream) || visited(upstream),continue;end
                visited(upstream)=true;
                depth(upstream)=depth(downstream)+1;
                nextMet(upstream)=downstream;
                nextRx(upstream)=rx(q);
                tail=tail+1;queue(tail)=upstream;
            end
        end
    end
    sources=find(visited & allowed);
    metPaths=cell(numel(sources),1);rxnPaths=cell(numel(sources),1);
    for si=1:numel(sources)
        source=sources(si);current=source;
        pathMets=source;pathRxns=zeros(1,0);
        while current~=target
            r=nextRx(current);downstream=nextMet(current);
            assert(r>0 && downstream>0,'ExchangeSearch:BrokenPath', ...
                'Could not reconstruct an intracellular precursor path.');
            pathRxns(end+1)=r; %#ok<AGROW>
            pathMets(end+1)=downstream; %#ok<AGROW>
            current=downstream;
        end
        metPaths{si}=pathMets;
        rxnPaths{si}=pathRxns;
    end
end

function [out,tests,stats]=d1FindMotifs(target,rx,sgn,products,substrates,metIDs,comp,currency,maxBranch,maxCand,maxTests,maxFound,S,requireProducts,inputAllowed,inputDepth)
% Search exact seven-reaction topology, INCLUDING r_I_1.
% r_Sd_1 must PRODUCE target m_Sd_1. The canceler must CONSUME target
% and m_Fd_1 in equal stoichiometric amounts in the allowed orientation.
% The regulation of transporter is NOT represented as metabolite consumption.
    out=struct([]);tests=0;
    stats=struct('testLimitReached',false,'branchLimitReached',false, ...
        'retentionLimitReached',false);
    stats.chemicalKeys=containers.Map('KeyType','char','ValueType','logical');
    prodTarget=find(cellfun(@(x)any(x==target),products));
    stats.branchLimitReached=stats.branchLimitReached || numel(prodTarget)>maxCand;
    prodTarget=prodTarget(1:min(numel(prodTarget),maxCand));
    srAll=find(cellfun(@(x)any(x==target),substrates));
    stats.branchLimitReached=stats.branchLimitReached || numel(srAll)>maxBranch;
    srAll=srAll(1:min(numel(srAll),maxBranch));
    % Preindex all reactions consuming the target, including possible cancelers.
    cancelAll=find(cellfun(@(x)any(x==target),substrates));
    for a=reshape(prodTarget,1,[])
        for msf=reshape(substrates{a},1,[])
            if ~d1SignalOK(msf,target,metIDs,comp,currency),continue;end
            sfCandidates=find(cellfun(@(x)any(x==msf),products));
            stats.branchLimitReached=stats.branchLimitReached || numel(sfCandidates)>maxBranch;
            sfCandidates=sfCandidates(1:min(numel(sfCandidates),maxBranch));
            for b=reshape(sfCandidates,1,[])
                if rx(b)==rx(a),continue;end
                inputs=substrates{b};inputs=inputs(inputAllowed(inputs));
                [~,inputOrder]=sort(inputDepth(inputs));inputs=inputs(inputOrder);
                for minMet=reshape(inputs,1,[])
                    if ~d1SignalOK(minMet,msf,metIDs,comp,currency) || minMet==target,continue;end
                    ffCandidates=find(cellfun(@(x)any(x==minMet),substrates));
                    stats.branchLimitReached=stats.branchLimitReached || numel(ffCandidates)>maxCand;
                    ffCandidates=ffCandidates(1:min(numel(ffCandidates),maxCand));
                    for c=reshape(ffCandidates,1,[])
                        if ismember(rx(c),[rx(a),rx(b)]),continue;end
                        for mff=reshape(products{c},1,[])
                            if ~d1SignalOK(mff,minMet,metIDs,comp,currency) || ismember(mff,[target,msf]),continue;end
                            fdCandidates=find(cellfun(@(x)any(x==mff),substrates));
                            stats.branchLimitReached=stats.branchLimitReached || numel(fdCandidates)>maxBranch;
                            fdCandidates=fdCandidates(1:min(numel(fdCandidates),maxBranch));
                            for d=reshape(fdCandidates,1,[])
                                if ismember(rx(d),[rx(a),rx(b),rx(c)]),continue;end
                                for mfd=reshape(products{d},1,[])
                                    if ~d1SignalOK(mfd,mff,metIDs,comp,currency) || ismember(mfd,[minMet,msf,target]),continue;end
                                    % Require a real reaction consuming BOTH terminal pools.
                                    cancels=cancelAll(cellfun(@(x)any(x==mfd),substrates(cancelAll)));
                                    goodCancel=zeros(0,1);coeffs=zeros(0,1);
                                    for z=reshape(cancels,1,[])
                                        if ismember(rx(z),[rx(a),rx(b),rx(c),rx(d)]),continue;end
                                        if requireProducts && isempty(products{z}),continue;end
                                        cf=-sgn(z)*full(S(mfd,rx(z)));
                                        cs=-sgn(z)*full(S(target,rx(z)));
                                        if cf<=0 || cs<=0 || abs(cf-cs)>1e-9*max([1,cf,cs]),continue;end
                                        goodCancel(end+1,1)=z;coeffs(end+1,1)=cf; %#ok<AGROW>
                                    end
                                    if isempty(goodCancel),continue;end
                                    frCandidates=find(cellfun(@(x)any(x==mfd),substrates));
                                    stats.branchLimitReached=stats.branchLimitReached || numel(frCandidates)>maxBranch;
                                    frCandidates=frCandidates(1:min(numel(frCandidates),maxBranch));
                                    for e=reshape(frCandidates,1,[])
                                        if ismember(rx(e),[rx(a),rx(b),rx(c),rx(d)]),continue;end
                                        % Removal reaction must NOT itself be the shared canceler.
                                        if any(substrates{e}==target),continue;end
                                        for f=reshape(srAll,1,[])
                                            if ismember(rx(f),[rx(a),rx(b),rx(c),rx(d),rx(e)]),continue;end
                                            if any(substrates{f}==mfd),continue;end
                                            for zi=1:numel(goodCancel)
                                                z=goodCancel(zi);
                                                if tests>=maxTests
                                                    stats.testLimitReached=true;
                                                    out=keepClosestD1Motifs(out,inputDepth,maxFound);
                                                    return;
                                                end
                                                tests=tests+1;
                                                if ismember(rx(z),[rx(a),rx(b),rx(c),rx(d),rx(e),rx(f)]),continue;end
                                                obj=struct('min',minMet,'mff',mff,'mfd',mfd,'msf',msf, ...
                                                    'rff',rx(c),'rfd',rx(d),'rfr',rx(e), ...
                                                    'rsf',rx(b),'rsd',rx(a),'rsr',rx(f),'ri',rx(z), ...
                                                    'off',sgn(c),'ofd',sgn(d),'ofr',sgn(e), ...
                                                    'osf',sgn(b),'osd',sgn(a),'osr',sgn(f), ...
                                                    'oi',sgn(z),'cancelCoeff',coeffs(zi));
                                                localKey=sprintf('%d,',[minMet,mff,mfd,msf,target, ...
                                                    rx(c),rx(d),rx(e),rx(b),rx(a),rx(f),rx(z)]);
                                                if isKey(stats.chemicalKeys,localKey),continue;end
                                                stats.chemicalKeys(localKey)=true;
                                                if isempty(out),out=obj;else,out(end+1)=obj;end %#ok<AGROW>
                                                stats.retentionLimitReached=stats.retentionLimitReached || numel(out)>maxFound;
                                                if isfinite(maxFound) && numel(out)>2*maxFound
                                                    out=keepClosestD1Motifs(out,inputDepth,maxFound);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    out=keepClosestD1Motifs(out,inputDepth,maxFound);
end

function motifs=keepClosestD1Motifs(motifs,depth,limit)
% Keep the closest discovered D1s; do not stop discovery at the first match.
    if isempty(motifs),return;end
    distances=reshape(depth([motifs.min]),[],1);
    serial=(1:numel(motifs))';
    [~,order]=sortrows([distances,serial],[1 2]);
    motifs=motifs(order);
    if isfinite(limit) && numel(motifs)>limit
        motifs=motifs(1:limit);
    end
end

function [parents,distances]=rankD1Parents(parents,depth)
% Prefer the shortest m_Sd-to-counter path, then the shortest cAMP-to-m_in path.
    distances=zeros(numel(parents),1);
    if isempty(parents),return;end
    for k=1:numel(parents)
        distances(k)=depth(parents(k).d.min);
    end
    assert(all(isfinite(distances)),'Ranking:DisconnectedParent', ...
        'Only parents with a cAMP path may enter the final ranking.');
    counterDistances=reshape(arrayfun(@(p)numel(p.exchangePathRxns),parents),[],1);
    nonExchange=reshape([parents.exchangeRow]==0,[],1);
    serial=(1:numel(parents))';
    [~,order]=sortrows([counterDistances,distances,double(nonExchange),serial],[1 2 3 4]);
    parents=parents(order);distances=distances(order);
end

function motifs=keepTopCompleteMotifs(motifs,limit)
    if isempty(motifs),return;end
    counterDistances=reshape([motifs.counterDistance],[],1);
    distances=reshape([motifs.distance],[],1);
    serial=reshape([motifs.discoveryIndex],[],1);
    [~,order]=sortrows([counterDistances,distances,serial],[1 2 3]);
    motifs=motifs(order);
    if numel(motifs)>limit,motifs=motifs(1:limit);end
end


function mergeChemicalKeys(destination,source)
% containers.Map is a handle: merge exact chemical keys across search calls.
    names=keys(source);
    for k=1:numel(names)
        destination(names{k})=true;
    end
end


function key=d1ChemicalKey(parent)
% Different regulatory annotations of the same chemical motif count once.
    d=parent.d;
    key=sprintf('%d,',[d.min,d.mff,d.mfd,d.msf,parent.target, ...
        d.rff,d.rfd,d.rfr,d.rsf,d.rsd,d.rsr,d.ri]);
end

function key=counterChemicalKey(parent,rTransport,rCr,mc)
    pathRxns=zeros(1,0);
    if isfield(parent,'exchangePathRxns'),pathRxns=parent.exchangePathRxns;end
    key=[d1ChemicalKey(parent),'|',sprintf('%d,', ...
        [pathRxns,rTransport,rCr,mc])];
end

function selfTestDistanceRanking()
    depth=[5;1;3;2];
    motifs=struct('min',{1,3,4,2});
    nearest=keepClosestD1Motifs(motifs,depth,2);
    assert(isequal([nearest.min],[2 4]), ...
        'Ranking self-test: a later shorter D1 was discarded.');
    p(1)=struct('d',struct('min',1),'exchangeRow',1,'exchangePathRxns',[1 2]);
    p(2)=struct('d',struct('min',2),'exchangeRow',0,'exchangePathRxns',[]);
    p(3)=struct('d',struct('min',4),'exchangeRow',1,'exchangePathRxns',3);
    [sorted,dist]=rankD1Parents(p,depth);
    assert(isequal(dist,[1;2;5]) && sorted(1).exchangeRow==0 && ...
        numel(sorted(1).exchangePathRxns)==0, ...
        'Ranking self-test: counter-path distance was not prioritized.');
    complete=struct('counterDistance',{1,0,0,1},'distance',{1,5,3,1}, ...
        'discoveryIndex',{1,2,3,4});
    complete=keepTopCompleteMotifs(complete,2);
    assert(isequal([complete.discoveryIndex],[3 2]), ...
        'Ranking self-test: complete motifs did not prioritize counter-path distance.');
    fprintf('PASS: counter-path priority, cAMP-distance ties, and later shorter motifs.\n');
end

function selfTestMotifCounts()
% Different annotations of the same chemistry must not inflate motif counts.
    d=struct('min',1,'mff',2,'mfd',3,'msf',4, ...
        'rff',1,'rfd',2,'rfr',3,'rsf',4,'rsd',5,'rsr',6,'ri',7);
    p=struct('d',d,'target',5,'transportRow',1,'regulatorRow',1,'exchangeRow',0);
    q=p;q.regulatorRow=99;
    assert(strcmp(d1ChemicalKey(p),d1ChemicalKey(q)), ...
        'Count self-test: protein annotations changed chemical identity.');
    c1=counterChemicalKey(p,20,21,6);
    c2=counterChemicalKey(q,20,21,6);
    c3=counterChemicalKey(p,20,22,6);
    counterSet=containers.Map('KeyType','char','ValueType','logical');
    counterSet(c1)=true;counterSet(c2)=true;counterSet(c3)=true;
    assert(counterSet.Count==2, ...
        'Count self-test: counter identity ignored duplicates or loss reaction.');
    A=containers.Map('KeyType','char','ValueType','logical');
    B=containers.Map('KeyType','char','ValueType','logical');
    A(c1)=true;B(c1)=true;B(c3)=true;mergeChemicalKeys(A,B);
    assert(A.Count==2,'Count self-test: merging repeated keys inflated the count.');
    records=struct('counterDistance',num2cell(zeros(1,12)), ...
        'distance',num2cell(12:-1:1),'discoveryIndex',num2cell(1:12));
    printed=keepTopCompleteMotifs(records,10);
    assert(numel(records)==12 && numel(printed)==10 && ...
        isequal([printed.distance],1:10), ...
        'Count self-test: printing cap stopped counting or discarded later shorter motifs.');
    fprintf('PASS: chemical deduplication, counter identity and found-versus-printed counts.\n');
end

function tf=d1SignalOK(m,other,metIDs,comp,currency)
    tf=(m~=other && comp(m)=="c");
    if ~tf,return;end
    base=regexprep(char(metIDs(m)),'\[c\]$','');
    tf=~ismember(lower(string(base)),currency);
end
function ec=d1EC(model,j)
    ec="";
    if isfield(model,'rxnECNumbers') && numel(model.rxnECNumbers)>=j
        ec=string(model.rxnECNumbers(j));
    elseif isfield(model,'rxnEC') && numel(model.rxnEC)>=j
        ec=string(model.rxnEC(j));
    end
    if ismissing(ec) || lower(ec)=="nan",ec="";end
end

function tf=reactionHasExcludedEC(model,reactionIDs,excludedECs)
% True if any reaction annotation matches an excluded EC number.
    tf=false;
    excludedECs=string(excludedECs(:));
    for jj=reshape(unique(reactionIDs),1,[])
        ec=d1EC(model,jj);
        if ec=="",continue;end
        found=string(regexp(char(ec),'\d+\.\d+\.\d+\.\d+','match'));
        if any(ismember(found,excludedECs))
            tf=true;
            return;
        end
    end
end
