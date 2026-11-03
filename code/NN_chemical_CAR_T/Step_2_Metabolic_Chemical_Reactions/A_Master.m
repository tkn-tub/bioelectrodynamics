%% A_Master.m - Strict +1 transport, differentiator, counter and two neurons
% ONE FILE: discover exchange transport -> transporter GPR -> D1 candidates.
% Complete motifs use stoichiometric exchange only; protein-path regulation is excluded.
% Requires COBRA Toolbox readCbModel, MATLAB with webread (online enrichment optional).
% Run: matlab -batch "A_Master"
% No checkpoints or network-cache directories are generated.
clearvars; clc;
MODEL_NAME = 'Puniya2021_TNM1055.xml';
OUT_NAME = 'Transport_Stages1to5_BRENDA_Results_EC_Required';
% Every selected chemical reaction EXCEPT the r_T role requires a complete
% EC annotation. Model fields are merged; missing annotations may be resolved
% by exact reaction-ID links in BiGG or a sourced BRENDA reaction mapping.
EC_TRY_BIGG = true;              % Also requires ONLINE=true.
EC_BIGG_TIMEOUT = 10;             % Seconds per request; failures do not admit reactions.
EC_BIGG_MAX_LOOKUPS = Inf;
EC_BIGG_MAX_CONSECUTIVE_FAILURES = 3;
EC_EXTERNAL_MAPPING_FILE = 'Reaction_EC_Evidence.tsv'; % Optional, found with which().
EC_BRENDA_FILE = 'brenda_2026_1.txt'; % Validates BRENDA-sourced mapping EC entries.
% Optional TSV columns: ReactionID, EC, Source, Reference.
% Source must be BRENDA or BiGG; exact reaction IDs and a source reference
% are required. No mapping is inferred merely from a similar enzyme name.


%% Ranking: distance means the number of cAMP -> m_in reactions.
PRINT_TOP_MOTIFS = 10;
% Exclude complete motifs containing reactions annotated with these EC numbers.
% EC annotations are read from the model reaction fields (rxnECNumbers/rxnEC).
EXCLUDED_MOTIF_ECS = ["2.4.1.23"];
REQUIRE_UNIQUE_COMPLETE_ECS = true; % Each fully specified EC can be assigned to at most one reaction in a complete motif.
BLOCK0_SOURCE_MET_ID = 'camp[c]';
BLOCK0_MAX_DEPTH = Inf;
BLOCK0_MAX_EXPANSIONS = Inf;
BLOCK0_ALLOW_ZERO_STEP = false;
BLOCK0_DIAGNOSTICS_ONLY = false;
% Per constrained path search. Inf removes a budget; capped failures are
% reported as unresolved within the budget, never as biological absence.
BLOCK0_MAX_PATH_STATES = Inf;
% Currency metabolites may occur in equations, but cannot carry a Block 0
% connection or be m_in. Match exact base IDs, without [c] or _c suffixes.
BLOCK0_CURRENCY_BASE_IDS = block0DefaultCurrencyIDs(); % Edit this list/helper if needed.
% AMP is deliberately retained: cAMP -> AMP is a valid possible first step.
RUN_RANKING_SELF_TESTS = true;
% r_I kinetic ranking: structural counter distance ALWAYS remains first.
RI_KCAT_RANKING = true;
RI_KCAT_ORGANISM = 'Homo sapiens';
RI_KCAT_BRENDA_FILE = 'brenda_2026_1.txt'; % located with findNamed(codeRoot,...)
RI_KCAT_CACHE_FILE = ''; % empty -> output directory/rI_brenda_ec_index.mat
RI_KCAT_ALLOW_ANY_ORGANISM = true;
RI_KCAT_ALLOW_ANY_SUBSTRATE = true;
RI_KCAT_ALLOW_SA_MW_ESTIMATE = true;
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
N8_MAX_TESTS = Inf; % Legacy global cap disabled; Stage 8 now uses a fair per-parent budget.
N8_MAX_TESTS_PER_PARENT = 5000; % Deterministic budget for every retained D1 parent.
N8_USE_PARALLEL = true;         % Use parfor when Parallel Computing Toolbox is available.
N8_PARALLEL_WORKERS = 0;        % 0 = scheduler/default pool size; otherwise requested worker count.
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
        'EXCHANGE_D1_MAX_TARGETS','EXCHANGE_D1_MAX_CONNECTIONS','N8_MAX_BRANCH','N8_MAX_TESTS_PER_PARENT'}
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
assert(isnumeric(EC_BIGG_TIMEOUT) && isscalar(EC_BIGG_TIMEOUT) && ...
    isfinite(EC_BIGG_TIMEOUT) && EC_BIGG_TIMEOUT>0, ...
    'ECPolicy:Timeout','EC_BIGG_TIMEOUT must be a positive finite number.');
assert(isnumeric(EC_BIGG_MAX_LOOKUPS) && isscalar(EC_BIGG_MAX_LOOKUPS) && ...
    ~isnan(EC_BIGG_MAX_LOOKUPS) && EC_BIGG_MAX_LOOKUPS>=0 && ...
    (isinf(EC_BIGG_MAX_LOOKUPS) || EC_BIGG_MAX_LOOKUPS==floor(EC_BIGG_MAX_LOOKUPS)), ...
    'ECPolicy:LookupCap','EC_BIGG_MAX_LOOKUPS must be a nonnegative integer or Inf.');
assert(isnumeric(EC_BIGG_MAX_CONSECUTIVE_FAILURES) && ...
    isscalar(EC_BIGG_MAX_CONSECUTIVE_FAILURES) && ...
    isfinite(EC_BIGG_MAX_CONSECUTIVE_FAILURES) && ...
    EC_BIGG_MAX_CONSECUTIVE_FAILURES>=1 && ...
    EC_BIGG_MAX_CONSECUTIVE_FAILURES==floor(EC_BIGG_MAX_CONSECUTIVE_FAILURES), ...
    'ECPolicy:FailureCap','EC_BIGG_MAX_CONSECUTIVE_FAILURES must be a positive integer.');
if RUN_RANKING_SELF_TESTS
    selfTestBlock0Connectivity();
    selfTestBlock0CurrencyPolicy();
    selfTestBlock0Isolation();
    selfTestDistanceRanking();
    selfTestMotifCounts();
    selfTestRequiredECPolicy();
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
% Resolve annotations once, BEFORE Block 0 and all reaction candidate pools.
% Workers receive this immutable annotation table; no parallel network lookups.
ecOptions=struct('online',ONLINE && EC_TRY_BIGG,'timeout',EC_BIGG_TIMEOUT, ...
    'maxLookups',EC_BIGG_MAX_LOOKUPS, ...
    'maxFailures',EC_BIGG_MAX_CONSECUTIVE_FAILURES, ...
    'mappingFile',EC_EXTERNAL_MAPPING_FILE,'brendaFile',EC_BRENDA_FILE);
[model,ECAnnotationAudit,ecResolutionInfo]=resolveRequiredECAnnotations( ...
    model,comp,ub,codeRoot,ecOptions);
writetable(ECAnnotationAudit,fullfile(outDir,'00_Reaction_EC_Annotation_Audit.csv'));
writetable(ECAnnotationAudit(ECAnnotationAudit.CandidateNonTransport & ...
    ~ECAnnotationAudit.HasCompleteEC,:), ...
    fullfile(outDir,'00_Excluded_Missing_EC_Reactions.csv'));
report(logID,sprintf('EC policy: %d/%d forward intracellular reactions have a complete EC; %d excluded as unresolved. r_T alone is exempt.', ...
    nnz(ECAnnotationAudit.CandidateNonTransport & ECAnnotationAudit.HasCompleteEC), ...
    nnz(ECAnnotationAudit.CandidateNonTransport), ...
    nnz(ECAnnotationAudit.CandidateNonTransport & ~ECAnnotationAudit.HasCompleteEC)));
report(logID,ecResolutionInfo);

% A single forward-only breadth-first search supplies every later motif.
block0Cache=buildBlock0Cache(model,metIDs,comp,ub, ...
    BLOCK0_SOURCE_MET_ID,BLOCK0_MAX_DEPTH,BLOCK0_MAX_EXPANSIONS, ...
    BLOCK0_CURRENCY_BASE_IDS);
assert(isscalar(BLOCK0_MAX_PATH_STATES) && isnumeric(BLOCK0_MAX_PATH_STATES) && ...
    ~isnan(BLOCK0_MAX_PATH_STATES) && BLOCK0_MAX_PATH_STATES>=1 && ...
    (isinf(BLOCK0_MAX_PATH_STATES) || BLOCK0_MAX_PATH_STATES==floor(BLOCK0_MAX_PATH_STATES)), ...
    'Block0:StateCap','BLOCK0_MAX_PATH_STATES must be a positive integer or Inf.');
block0Cache.pathOptions=struct('maxStates',BLOCK0_MAX_PATH_STATES, ...
    'excludedECs',EXCLUDED_MOTIF_ECS,'requireUniqueEC',REQUIRE_UNIQUE_COMPLETE_ECS);
block0InputAllowed=block0Cache.reachable;
if ~BLOCK0_ALLOW_ZERO_STEP
    block0InputAllowed(block0Cache.start)=false;
end
report(logID,sprintf('Block 0 structural prefilter: %d candidate m_in species reachable from %s; depth <= %g; expansions = %d; expansion cap reached = %d.', ...
    nnz(block0InputAllowed),BLOCK0_SOURCE_MET_ID,BLOCK0_MAX_DEPTH, ...
    block0Cache.expansions,block0Cache.expansionLimitReached));
report(logID,sprintf('Block 0: %d cytosolic currency metabolites excluded as connectors and m_in; permitted as reaction co-participants.', ...
    nnz(block0Cache.currencyMask & comp=="c")));
report(logID,['Block 0 excluded currency base IDs: ' char(join(block0Cache.currencyBaseIDs,', '))]);
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
    block0Status(:)="no_path_in_forward_cytosolic_noncurrency_graph";
end
if block0Cache.expansionLimitReached
    block0Status(:)="not_reached_before_expansion_cap";
end
block0Status(block0Cache.reachable(block0Rows))="structural_prefilter_only";
block0Status(block0Cache.currencyMask(block0Rows))="excluded_currency_connector";
if ~BLOCK0_ALLOW_ZERO_STEP
    block0Status(block0Rows==block0Cache.start)="zero_step_not_allowed";
end
block0MetPaths=strings(numel(block0Rows),1);
block0ReactionPaths=strings(numel(block0Rows),1);
block0InternalPass=false(numel(block0Rows),1);
block0InternalReason=repmat("not_reached",numel(block0Rows),1);
for b0row=1:numel(block0Rows)
    b0target=block0Rows(b0row);
    if block0Cache.reachable(b0target)
        [b0mets,b0rx]=block0CachedPath(block0Cache,b0target);
        [block0InternalPass(b0row),block0InternalReason(b0row)]= ...
            validateBlock0Isolation(block0Cache,b0mets,b0rx,[],[]);
        block0MetPaths(b0row)=join(metIDs(b0mets),' -> ');
        if ~isempty(b0rx)
            block0ReactionPaths(b0row)=join(rxnIDs(b0rx),' -> ');
        end
    end
end
Block0Reachability=table(metIDs(block0Rows),metNames(block0Rows), ...
    block0Cache.reachable(block0Rows),block0InputAllowed(block0Rows), ...
    block0Cache.depth(block0Rows),block0Status,block0MetPaths,block0ReactionPaths, ...
    block0Cache.currencyMask(block0Rows), ...
    block0InternalPass,block0InternalReason, ...
    'VariableNames',{'MetaboliteID','MetaboliteName','Reachable', ...
    'AllowedAsInput','ReactionDepth','Status','PathMetaboliteIDs','PathReactionIDs', ...
    'IsCurrency','CachedPathInternallyIsolated','CachedPathIsolationStatus'});
writetable(Block0Reachability,fullfile(outDir,'00_Block0_Reachable_Inputs.csv'));
report(logID,['Block 0 reachability report: ' fullfile(outDir,'00_Block0_Reachable_Inputs.csv')]);
report(logID,'Block 0 CSV paths are optimistic prefilters; accepted motifs use separately searched, isolated paths.');
if BLOCK0_DIAGNOSTICS_ONLY
    report(logID,'Block 0 diagnostics complete; EC annotation resolution and Block 0 only. Downstream circuit search was skipped.');
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
                if coeff>0 && ub(r)>0 && ecReactionAllowed(model,r)
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
        printModelReaction(fid6,model,ri,orient,metIDs,metNames,T.cytosolic(chosenT),'m_c',[],[],true);
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
    'exchangePathMets',{},'exchangePathRxns',{},'block0Path',{});
block0Memo=containers.Map('KeyType','char','ValueType','any');
block0ParentRejected=0;block0ParentCapped=0;
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
            parentDownstream=[d1ids,double(T.modelIndex(chosenT))];
            b0=memoizedBlock0Path(model,block0Cache,d.min,parentDownstream, ...
                [d.rff,d.rsf],[],block0Memo);
            if ~b0.found
                block0ParentRejected=block0ParentRejected+1;
                block0ParentCapped=block0ParentCapped+b0.expansionLimitReached+b0.stateLimitReached;
                continue;
            end
            block0Reactions=b0.rxns;
            motifPrefix=[block0Reactions,d1ids,double(T.modelIndex(chosenT))];
            if ~reactionsHaveRequiredEC(model,[block0Reactions,d1ids]),continue;end
            if reactionHasExcludedEC(model,motifPrefix,EXCLUDED_MOTIF_ECS),continue;end
            assert(block0InputAllowed(d.min),'Ranking:UnconnectedD1','Unconnected input cannot be selected.');
            foundD1=foundD1+1;
            selectedD1(end+1)=struct('d',d,'target',target,'transportRow',chosenT, ...
                'regulatorRow',si,'exchangeRow',0, ...
                'exchangePathMets',zeros(1,0),'exchangePathRxns',zeros(1,0), ...
                'block0Path',b0); %#ok<SAGROW>
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
            printModelReaction(fid7,model,double(T.modelIndex(chosenT)),orient,metIDs,metNames,T.cytosolic(chosenT),'m_c',[],[],true);
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
                    parentDownstream=[d1ids,pathRxns,double(T.modelIndex(tr))];
                    b0=memoizedBlock0Path(model,block0Cache,d.min,parentDownstream, ...
                        [d.rff,d.rsf],[],block0Memo);
                    if ~b0.found
                        block0ParentRejected=block0ParentRejected+1;
                        block0ParentCapped=block0ParentCapped+b0.expansionLimitReached+b0.stateLimitReached;
                        continue;
                    end
                    block0Reactions=b0.rxns;
                    motifPrefix=[block0Reactions,d1ids,pathRxns,double(T.modelIndex(tr))];
                    if ~reactionsHaveRequiredEC(model,[block0Reactions,d1ids,pathRxns]),continue;end
                    if reactionHasExcludedEC(model,motifPrefix,EXCLUDED_MOTIF_ECS),continue;end
                    assert(block0InputAllowed(d.min),'Ranking:UnconnectedExchange','Unconnected input cannot be selected.');
                    exchangeFound=exchangeFound+1;
                    bestCounterPathDepth=currentPathDepth;
                    selectedD1(end+1)=struct('d',d,'target',signalMet,'transportRow',tr, ...
                        'regulatorRow',0,'exchangeRow',xi, ...
                        'exchangePathMets',pathMets,'exchangePathRxns',pathRxns, ...
                        'block0Path',b0); %#ok<SAGROW>
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
            fprintf(fid7x,'Isolated Block 0 path (%d reactions): %s\n', ...
                b0.distance,char(join(metIDs(b0.mets),' -> ')));
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
% Enrich each retained D1 parent with canceller (r_I) kcat BEFORE ranking.
% The lookup uses the local BRENDA flat file and caches only ECs actually needed.
% Ranking remains lexicographic: counter distance first, then r_I evidence/kcat,
% then cAMP-to-m_in distance. Missing kinetics never invalidate a motif.
for riq=1:numel(selectedD1)
    selectedD1(riq).rIKcat=NaN;
    selectedD1(riq).rIKcatEvidence=0;
    selectedD1(riq).rIKcatSource="not queried";
end
if RI_KCAT_RANKING && ~isempty(selectedD1)
    riBrenda=findNamed(codeRoot,RI_KCAT_BRENDA_FILE);
    if isempty(riBrenda)
        report(logID,['r_I kinetics: local BRENDA file not found: ' RI_KCAT_BRENDA_FILE '. Ranking falls back to distances.']);
    else
        if isempty(RI_KCAT_CACHE_FILE)
            riCache=fullfile(outDir,'rI_brenda_ec_index.mat');
        else
            riCache=RI_KCAT_CACHE_FILE;
        end
        [selectedD1,RIKinetics]=attachRIKcat(selectedD1,model,metIDs,metNames, ...
            RI_KCAT_ORGANISM,RI_KCAT_ALLOW_ANY_ORGANISM,RI_KCAT_ALLOW_ANY_SUBSTRATE, ...
            RI_KCAT_ALLOW_SA_MW_ESTIMATE,riBrenda,riCache);
        writetable(RIKinetics,fullfile(outDir,'18a_rI_Kinetics.csv'));
        report(logID,sprintf('r_I kinetics: %d/%d parents received a positive kcat; details: 18a_rI_Kinetics.csv', ...
            nnz(isfinite([selectedD1.rIKcat]) & [selectedD1.rIKcat]>0),numel(selectedD1)));
    end
end
[selectedD1,parentDistances]=rankD1Parents(selectedD1,block0Cache.depth);
report(logID,sprintf('Block 0 isolation: %d parent connections rejected; %d budget-limited rejections; %d distinct path queries memoized.', ...
    block0ParentRejected,block0ParentCapped,block0Memo.Count));
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
parentRIReaction=strings(numel(selectedD1),1); parentRIKcat=nan(numel(selectedD1),1);
parentRIEvidence=zeros(numel(selectedD1),1); parentRISource=strings(numel(selectedD1),1);
for pr=1:numel(selectedD1)
    parentInputIDs(pr)=metIDs(selectedD1(pr).d.min);
    parentTransportIDs(pr)=string(T.reaction(selectedD1(pr).transportRow));
    parentCounterDistances(pr)=numel(selectedD1(pr).exchangePathRxns);
    parentMode(pr)="stoichiometric_exchange";
    parentRIReaction(pr)=rxnIDs(selectedD1(pr).d.ri);
    if isfield(selectedD1,'rIKcat'), parentRIKcat(pr)=selectedD1(pr).rIKcat; end
    if isfield(selectedD1,'rIKcatEvidence'), parentRIEvidence(pr)=selectedD1(pr).rIKcatEvidence; end
    if isfield(selectedD1,'rIKcatSource'), parentRISource(pr)=selectedD1(pr).rIKcatSource; end
end
ParentRanking=table(parentRank,parentCounterDistances,parentRIEvidence,parentRIKcat,parentRIReaction, ...
    parentRISource,parentDistances,parentInputIDs,parentTransportIDs,parentMode,'VariableNames', ...
    {'ParentRank','m_SdToCounterReactionDistance','r_I_KcatEvidenceRank','r_I_kcat_s', ...
    'r_I_ReactionID','r_I_KcatSource','CAMPReactionDistance','m_in_ID','TransportReactionID','TransportMode'});
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
fprintf(fid8,'Selection: up to %d complete motifs, ranked by m_Sd-to-counter distance, r_I kinetic evidence/kcat, then isolated cAMP-to-m_in distance.\n',PRINT_TOP_MOTIFS);
fprintf(fid8,'Transport mechanism filter: stoichiometric exchange only; protein-path regulation is excluded.\n');
fprintf(fid8,'EC policy: complete annotation required for Block 0, D1, intracellular connectors, r_cr and every neuron reaction; r_T is exempt.\n');
fprintf(fid8,'Annotation audit: 00_Reaction_EC_Annotation_Audit.csv; unresolved exclusions: 00_Excluded_Missing_EC_Reactions.csv.\n');
fprintf(fid8,'An EC annotation does not guarantee matching KM or kcat measurements.\n');
fprintf(fid8,'Excluded reaction ECs: %s\n',char(join(EXCLUDED_MOTIF_ECS,', ')));
fprintf(fid8,'Scope: retained D1 parents and configured branch/trial limits; global optimality is not assumed.\n');
fprintf(fid8,'Equal-distance motifs may share m_in. Ties retain discovery order after both path distances.\n');
fprintf(fid8,'Block 0 depth cap=%g; expansion cap=%g; expansion limit reached=%d.\n', ...
    BLOCK0_MAX_DEPTH,BLOCK0_MAX_EXPANSIONS,block0Cache.expansionLimitReached);
fprintf(fid8,'Block 0: currency metabolites cannot carry path connections or be m_in; they may remain in reaction equations.\n');
fprintf(fid8,'Block 0 excluded currency base IDs: %s\n',char(join(block0Cache.currencyBaseIDs,', ')));
fprintf(fid8,'Block 0 isolation: all non-currency participants checked; sharing within Block 0 only at adjacent chain steps.\n');
fprintf(fid8,'Block 0 shares only m_in with downstream chemistry, consumed only by r_Ff and r_Sf; reaction reuse is forbidden.\n');
fprintf(fid8,'Alternative Block 0 paths are searched per motif. Per-query expansion cap=%g; partial-path state cap=%g.\n', ...
    BLOCK0_MAX_EXPANSIONS,BLOCK0_MAX_PATH_STATES);
fprintf(fid8,'Upstream caps: protein hops=%g; metabolite-match rows=%g; backward metabolic steps=%g; backward path rows=%g.\n', ...
    MAX_HOPS,STAGE5_MAX_MATCH_ROWS,STAGE5_MAX_HOPS,STAGE5_MAX_PATH_ROWS);
fprintf(fid8,'D1 caps: targets=%g; motifs=%g; retained/target=%g; branch=%g; reaction candidates=%g; tests/target=%g.\n', ...
    D1_MAX_TARGETS,D1_MAX_MOTIFS,D1_MAX_PER_TARGET,D1_MAX_BRANCH,D1_MAX_RXN_CAND,D1_MAX_TESTS_PER_TARGET);
fprintf(fid8,'Exchange caps: pairs=%g; connections=%g. Neuron caps: branch=%g; tests/parent=%g; parallel=%d.\n', ...
    EXCHANGE_D1_MAX_TARGETS,EXCHANGE_D1_MAX_CONNECTIONS,N8_MAX_BRANCH,N8_MAX_TESTS_PER_PARENT,N8_USE_PARALLEL);
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
block0CounterRejected=0;block0CounterCapped=0;
for cp=1:numel(selectedD1)
    assert(selectedD1(cp).exchangeRow>0,'ExchangeOnly:ParentMode', ...
        'A regulatory protein-path parent reached the complete-motif search.');
    parent=selectedD1(cp);dd=parent.d;tt=parent.transportRow;
    selectedD1(cp).counterLinks=struct('reaction',{},'path',{});
    retainedD1Keys(d1ChemicalKey(parent))=true;
    mm=find(metIDs==string(T.cytosolic(tt)),1);
    if isempty(mm) || comp(mm)~="c",continue;end
    candidates=cons8{mm};candidates=candidates(1:min(end,N8_MAX_BRANCH));
    d1reactions=[dd.rff,dd.rfd,dd.rfr,dd.rsf,dd.rsd,dd.rsr,dd.ri];
    for cr=reshape(candidates,1,[])
        if ismember(di(cr),d1reactions),continue;end
        transportReaction=double(T.modelIndex(tt));
        if di(cr)==transportReaction,continue;end
        if ~ecReactionAllowed(model,di(cr)),continue;end
        counterDownstream=[d1reactions,parent.exchangePathRxns,transportReaction,di(cr)];
        counterPath=memoizedBlock0Path(model,block0Cache,dd.min,counterDownstream, ...
            [dd.rff,dd.rsf],parent.block0Path,block0Memo);
        if ~counterPath.found
            block0CounterRejected=block0CounterRejected+1;
            block0CounterCapped=block0CounterCapped+counterPath.expansionLimitReached+counterPath.stateLimitReached;
            continue;
        end
        selectedD1(cp).counterLinks(end+1)=struct('reaction',di(cr),'path',counterPath);
        ck=counterChemicalKey(parent,transportReaction,di(cr),mm);
        if isKey(counterKeys,ck),continue;end
        counterKeys(ck)=true;
        counterCoreKeys(sprintf('%d,',[transportReaction,di(cr),mm]))=true;
        mode="stoichiometric_exchange";
        counterRows(end+1,:)={size(counterRows,1)+1,cp, ...
            numel(parent.exchangePathRxns),counterPath.distance, ...
            metIDs(dd.min),metIDs(parent.target),metIDs(mm), ...
            rxnIDs(transportReaction),rxnIDs(di(cr)),mode,d1ChemicalKey(parent), ...
            join(rxnIDs(parent.exchangePathRxns),'; ')}; %#ok<SAGROW>
        counterRowKeys{end+1,1}=ck; %#ok<SAGROW>
    end
end
report(logID,sprintf('Counters: %d distinct D1-counter solutions; %d distinct transport/loss cores before neuron search.', ...
    counterKeys.Count,counterCoreKeys.Count));
report(logID,sprintf('Block 0 isolation at counter stage: %d rejected; %d budget limited.', ...
    block0CounterRejected,block0CounterCapped));

found8=0; tests8=0; examined8=0;
seen8=containers.Map('KeyType','char','ValueType','logical');
countersWithComplete=containers.Map('KeyType','char','ValueType','logical');
completeExchange8=0;completeRegulatory8=0;
best8=struct('parentIndex',{},'mc',{},'nf',{},'nd',{},'nf2',{}, ...
    'mout',{},'chosen',{},'block0Path',{},'counterDistance',{},'distance',{},'rIKcat',{},'rIKcatEvidence',{},'discoveryIndex',{});
rankingRows8=cell(0,12);

% Stage 8 is embarrassingly parallel by retained D1 parent. Each worker gets
% the SAME deterministic per-parent trial budget. No worker writes files or
% mutates shared maps; all deduplication/counting is done serially below.
nParents8=numel(selectedD1);
parentResults8=cell(nParents8,1);
useParallel8=false;
if N8_USE_PARALLEL && nParents8>1 && license('test','Distrib_Computing_Toolbox')
    try
        pool8=gcp('nocreate');
        if isempty(pool8)
            if N8_PARALLEL_WORKERS>0
                pool8=parpool(N8_PARALLEL_WORKERS);
            else
                pool8=parpool;
            end
        end
        useParallel8=~isempty(pool8);
        report(logID,sprintf('Stage 8 HPC: parallel=%d; workers=%d; per-parent trial cap=%d; parents=%d.', ...
            useParallel8,pool8.NumWorkers,N8_MAX_TESTS_PER_PARENT,nParents8));
    catch ME
        report(logID,['Stage 8 HPC: parallel pool unavailable; serial fallback. ' ME.message]);
        useParallel8=false;
    end
else
    report(logID,sprintf('Stage 8 HPC: parallel=0; per-parent trial cap=%d; parents=%d.', ...
        N8_MAX_TESTS_PER_PARENT,nParents8));
end

if useParallel8
    parfor ii=1:nParents8
        parentResults8{ii}=searchStage8Parent(ii,selectedD1,T,metIDs,comp,currency, ...
            cons8,di,do,model,block0Cache,EXCLUDED_MOTIF_ECS, ...
            REQUIRE_UNIQUE_COMPLETE_ECS,N8_MAX_BRANCH,N8_MAX_TESTS_PER_PARENT);
    end
else
    for ii=1:nParents8
        parentResults8{ii}=searchStage8Parent(ii,selectedD1,T,metIDs,comp,currency, ...
            cons8,di,do,model,block0Cache,EXCLUDED_MOTIF_ECS, ...
            REQUIRE_UNIQUE_COMPLETE_ECS,N8_MAX_BRANCH,N8_MAX_TESTS_PER_PARENT);
        if mod(ii,100)==0
            report(logID,sprintf('Stage 8 serial progress: %d/%d parents.',ii,nParents8));
        end
    end
end

% Serial deterministic merge. Parent order, then each parent's discovery
% order, defines the stable discoveryIndex independently of worker scheduling.
reject8=struct('duplicateReaction',0,'d1Overlap',0,'excludedEC',0, ...
    'repeatedEC',0,'missingEC',0,'block0Isolation',0,'block0PathCapped',0, ...
    'duplicateChemistry',0,'accepted',0,'parentsCapped',0);
for ii=1:nParents8
    R8=parentResults8{ii};
    if isempty(R8),continue;end
    examined8=examined8+R8.examined;
    tests8=tests8+R8.tests;
    reject8.duplicateReaction=reject8.duplicateReaction+R8.rejectDuplicateReaction;
    reject8.d1Overlap=reject8.d1Overlap+R8.rejectD1Overlap;
    reject8.excludedEC=reject8.excludedEC+R8.rejectExcludedEC;
    reject8.repeatedEC=reject8.repeatedEC+R8.rejectRepeatedEC;
    reject8.missingEC=reject8.missingEC+R8.rejectMissingEC;
    reject8.block0Isolation=reject8.block0Isolation+R8.rejectBlock0Isolation;
    reject8.block0PathCapped=reject8.block0PathCapped+R8.block0PathCapped;
    reject8.parentsCapped=reject8.parentsCapped+R8.capReached;
    for jj=1:numel(R8.candidates)
        c8=R8.candidates(jj);
        if isKey(seen8,c8.chemicalKey)
            reject8.duplicateChemistry=reject8.duplicateChemistry+1;
            continue;
        end
        seen8(c8.chemicalKey)=true;
        found8=found8+1;reject8.accepted=reject8.accepted+1;
        countersWithComplete(c8.counterKey)=true;
        sel=selectedD1(ii);
        if sel.exchangeRow>0,completeExchange8=completeExchange8+1;
        else,completeRegulatory8=completeRegulatory8+1;end
        item8=struct('parentIndex',ii,'mc',c8.mc,'nf',c8.nf,'nd',c8.nd, ...
            'nf2',c8.nf2,'mout',c8.mout,'chosen',c8.chosen,'block0Path',c8.block0Path, ...
            'counterDistance',numel(sel.exchangePathRxns), ...
            'distance',c8.block0Path.distance,'rIKcat',sel.rIKcat, ...
            'rIKcatEvidence',sel.rIKcatEvidence,'discoveryIndex',found8);
        best8(end+1)=item8; %#ok<SAGROW>
        best8=keepTopCompleteMotifs(best8,PRINT_TOP_MOTIFS);
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
fprintf(fid8,'Distinct structural D1 motifs discovered before Block 0 isolation: %d\n',discoveredD1Keys.Count);
fprintf(fid8,'Distinct D1 motifs retained for downstream search: %d\n',retainedD1Keys.Count);
fprintf(fid8,'D1/transporter parent connections retained: %d\n',numel(selectedD1));
fprintf(fid8,'Distinct counter cores (r_T, m_c, r_cr): %d\n',counterCoreKeys.Count);
fprintf(fid8,'Distinct D1-counter solutions in retained parent pool: %d\n',counterKeys.Count);
fprintf(fid8,'D1-counter solutions extended to both neurons: %d\n',countersWithComplete.Count);
fprintf(fid8,'Distinct complete motifs found: %d\n',found8);
fprintf(fid8,'  Exchange mode: %d; regulatory-hypothesis mode: %d\n',completeExchange8,completeRegulatory8);
fprintf(fid8,'Complete motifs printed: %d (requested maximum %d)\n',printed8,PRINT_TOP_MOTIFS);
fprintf(fid8,'Count definitions: chemical structures are deduplicated; alternate protein annotations do not add motifs.\n');
fprintf(fid8,'One shortest isolated Block 0 path is used per accepted motif; alternate connectors do not multiply motif counts.\n');
fprintf(fid8,'Parent/counter distances use their prefix constraints; complete-motif distances also include both neurons.\n');
fprintf(fid8,'Block 0 parent rejections=%d; parent path-budget limits=%d; counter rejections=%d; counter path-budget limits=%d.\n', ...
    block0ParentRejected,block0ParentCapped,block0CounterRejected,block0CounterCapped);
fprintf(fid8,'Final-assembly Block 0 isolation rejections=%d; path-budget limits=%d.\n', ...
    reject8.block0Isolation,reject8.block0PathCapped);
fprintf(fid8,'Counter total includes counters with no two-neuron solution; complete count requires both neurons.\n');
fprintf(fid8,'D1 parents examined for neurons: %d/%d; total trials=%d; parents reaching per-parent cap=%d.\n', ...
    examined8,numel(selectedD1),tests8,reject8.parentsCapped);
fprintf(fid8,'Missing EC rejections during final assembly: %d (EC filtering also occurred before candidate enumeration).\n',reject8.missingEC);
fprintf(fid8,'Stage 8 rejection diagnostics: duplicate reaction=%d; D1 overlap=%d; excluded EC=%d; repeated/unassignable EC=%d; duplicate chemistry=%d; accepted=%d.\n', ...
    reject8.duplicateReaction,reject8.d1Overlap,reject8.excludedEC,reject8.repeatedEC, ...
    reject8.duplicateChemistry,reject8.accepted);
fprintf(fid8,'Counts describe the configured search only; finite branch/trial/retention caps can hide additional solutions.\n');
fprintf(fid8,'Budget-limited path failures are unresolved searches, not proof that an isolated connection is absent.\n');
fprintf(fid8,'Selecting ten for printing does NOT stop discovery/counting at ten.\n\n');

block0AuditRows=cell(0,8);
for printedRank8=1:printed8
    item8=best8(printedRank8);
    ii=item8.parentIndex;sel=selectedD1(ii);d=sel.d;tr=sel.transportRow;
    mc=item8.mc;nf=item8.nf;nd=item8.nd;nf2=item8.nf2;mout=item8.mout;
    chosen=item8.chosen;distance=item8.distance;
    block0Reactions=item8.block0Path.rxns;
    finalDownstream=[d.rff,d.rfd,d.rfr,d.rsf,d.rsd,d.rsr,d.ri, ...
        sel.exchangePathRxns,double(T.modelIndex(tr)),reshape(di(chosen),1,[])];
    [isolated,why]=validateBlock0Isolation(block0Cache,item8.block0Path.mets, ...
        block0Reactions,finalDownstream,[d.rff,d.rsf]);
    assert(isolated,'Block0:FinalIsolation','Final Block 0 cross-talk check failed: %s.',char(why));
    participants=find(full(any(abs(model.S(:,block0Reactions))>1e-12,2)) & ~block0Cache.currencyMask);
    for bm=reshape(participants,1,[])
        role="side_participant";
        if bm==d.min,role="m_in";
        elseif bm==block0Cache.start,role="cAMP";
        elseif ismember(bm,item8.block0Path.mets),role="internal_connector";end
        producing=block0Reactions(full(model.S(bm,block0Reactions))>1e-12);
        consuming=block0Reactions(full(model.S(bm,block0Reactions))<-1e-12);
        shared=finalDownstream(abs(full(model.S(bm,finalDownstream)))>1e-12);
        producerText="";consumerText="";sharedText="";
        if ~isempty(producing),producerText=join(rxnIDs(producing),'; ');end
        if ~isempty(consuming),consumerText=join(rxnIDs(consuming),'; ');end
        if ~isempty(shared),sharedText=join(rxnIDs(shared),'; ');end
        block0AuditRows(end+1,:)={printedRank8,metIDs(bm),metNames(bm),role, ...
            producerText,consumerText,sharedText,isempty(shared) || bm==d.min}; %#ok<SAGROW>
    end
    requiredFinalReactions=[block0Reactions,d.rff,d.rfd,d.rfr,d.rsf,d.rsd, ...
        d.rsr,d.ri,sel.exchangePathRxns,reshape(di(chosen),1,[])];
    assert(reactionsHaveRequiredEC(model,requiredFinalReactions), ...
        'ECPolicy:FinalMotifMissingEC', ...
        'A final motif contains a non-r_T reaction without a complete EC annotation.');
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
    printBlock0(fid8,model,metIDs,metNames,d.min,signalIDs,signalNames, ...
        block0Cache,item8.block0Path,finalDownstream,[d.rff,d.rsf]);
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
    printModelReaction(fid8,model,double(T.modelIndex(tr)),ot,metIDs,metNames,"","",signalIDs,signalNames,true);
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
    fprintf(fid8,'Block 0 isolation was required, including all co-substrates/products and all selected downstream reactions.\n');
end
auditNames={'MotifRank','MetaboliteID','MetaboliteName','Role', ...
    'Block0ProducerReactions','Block0ConsumerReactions','DownstreamReactions','IsolationPass'};
if isempty(block0AuditRows)
    Block0IsolationAudit=table('Size',[0 8],'VariableTypes', ...
        {'double','string','string','string','string','string','string','logical'}, ...
        'VariableNames',auditNames);
else
    Block0IsolationAudit=cell2table(block0AuditRows,'VariableNames',auditNames);
end
writetable(Block0IsolationAudit,fullfile(outDir,'19b_Block0_Isolation_Audit.csv'));
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
        expectedOrder=keepTopCompleteMotifs(best8,numel(best8));
        assert(isequal([best8.discoveryIndex],[expectedOrder.discoveryIndex]), ...
            'Ranking:OutputOrder','Final ranking must follow counter distance, r_I evidence/kcat, then isolated cAMP distance.');
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
    % A supplied full path (including an off-path temporary test file) is
    % already a file location. Do not prepend pwd or root to it.
    if isfile(name),path=char(name);return;end
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

function ids=block0DefaultCurrencyIDs()
% Block 0 only: do not change the existing currency policy of other blocks.
% AMP and NH4 intentionally remain available as potential signal metabolites.
    ids=["h","h+","proton","h2o","water","atp","adp","pi","ppi", ...
        "nad","nadh","nadp","nadph","coa","co2","o2", ...
        "na1","na","k","k1","cl","cl1"];
end

function [mask,excluded]=block0CurrencyMask(metIDs,currencyBaseIDs)
% Exact base-ID matching, accepting [c], _c and an optional SBML M_ prefix.
% Never use substring matching: h[c] must not also exclude hco3[c], etc.
    [base,~]=splitCompartments(lower(strtrim(string(metIDs(:)))));
    base=regexprep(base,'^m_','');
    excluded=lower(strtrim(string(currencyBaseIDs(:))));
    assert(all(~ismissing(excluded)) && all(strlength(excluded)>0), ...
        'Block0:CurrencyList','Currency IDs must be nonempty, nonmissing strings.');
    [excluded,~]=splitCompartments(excluded);
    excluded=unique(regexprep(excluded,'^m_',''),'stable');
    mask=ismember(base,excluded);
end

function cache=buildBlock0Cache(model,metIDs,comp,ub,sourceID,maxDepth,maxExpansions,currencyBaseIDs)
% Explore ALL structurally reachable inputs once, using +1 reactions only.
% Only NON-CURRENCY substrates and products can carry path edges. Currency
% co-participants stay in S; reactions using water, protons, ATP, etc. remain
% available when they also link two allowed metabolites. Other substrates
% need not be produced by this path; atom transfer and flux are unverified.
    if nargin<8,currencyBaseIDs=block0DefaultCurrencyIDs();end
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
    [currencyMask,currencyBaseIDs]=block0CurrencyMask(metIDs,currencyBaseIDs);
    assert(~currencyMask(start),'Block0:CurrencyAnchor', ...
        'The Block 0 anchor cannot be an excluded currency metabolite.');
    cache=struct('start',start,'sourceID',string(sourceID), ...
        'maxDepth',maxDepth,'maxExpansions',maxExpansions, ...
        'currencyMask',currencyMask,'currencyBaseIDs',currencyBaseIDs, ...
        'reachable',false(nm,1),'depth',inf(nm,1), ...
        'parentMet',zeros(nm,1),'parentRx',zeros(nm,1), ...
        'expansions',0,'expansionLimitReached',false);
    cache.reachable(start)=true;cache.depth(start)=0;
    S=model.S;
    tol=1e-12;
    eligible=ub(:)>0 & ecSearchMask(model);
    eligible=eligible & full(any(S(~currencyMask,:) < -tol,1))' & ...
        full(any(S(~currencyMask,:) > tol,1))';
    % No extracellular or organelle participant may occur in Block 0.
    outside=comp~="c";
    eligible=eligible & ~full(any(abs(S(outside,:))>tol,1))';
    % Immutable adjacency for motif-specific searches. The original BFS below
    % remains an optimistic distance/reachability prefilter, not an accepted path.
    cache.S=S;cache.eligible=eligible;
    cache.consumers=cell(nm,1);cache.producers=cell(nm,1);
    cache.substrates=cell(size(S,2),1);cache.products=cell(size(S,2),1);
    cache.participants=cell(size(S,2),1);
    for j=reshape(find(eligible),1,[])
        a=find(S(:,j)<-tol & ~currencyMask);
        b=find(S(:,j)>tol & ~currencyMask);
        cache.substrates{j}=a;cache.products{j}=b;
        cache.participants{j}=find(abs(S(:,j))>tol & ~currencyMask);
        for m=reshape(a,1,[]),cache.consumers{m}(end+1)=j;end
        for m=reshape(b,1,[]),cache.producers{m}(end+1)=j;end
    end
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
            products=find(S(:,j)>tol & ~currencyMask);
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

function [ok,reason,blocked]=block0InterfaceMask(cache,target,downstreamRx,inputRx)
% m_in is the only permitted non-currency interface with downstream chemistry,
% and only the two D1 input reactions may consume it.
    downstreamRx=unique(reshape(downstreamRx,1,[]));
    inputRx=unique(reshape(inputRx,1,[]));
    blocked=false(size(cache.S,1),1);ok=false;reason="";
    if ~isempty(downstreamRx)
        if isempty(inputRx) || ~all(ismember(inputRx,downstreamRx)) || ...
                any(full(cache.S(target,inputRx))>=-1e-12)
            reason="invalid_D1_input_interface";return;
        end
        other=setdiff(downstreamRx,inputRx);
        if any(abs(full(cache.S(target,other)))>1e-12)
            reason="m_in_used_outside_D1_input_reactions";return;
        end
        blocked=full(any(abs(cache.S(:,downstreamRx))>1e-12,2)) & ~cache.currencyMask;
        blocked(target)=false;
    end
    ok=true;reason="PASS";
end

function [ok,reason]=validateBlock0Isolation(cache,mets,rx,downstreamRx,inputRx)
% Every repeated non-currency participant must be the intended intermediate
% between two adjacent steps. This checks ALL substrates/products, not only
% the metabolite chosen to draw each graph edge.
    mets=reshape(mets,1,[]);rx=reshape(rx,1,[]);
    ok=false;reason="invalid_path_shape";
    if isempty(mets) || numel(mets)~=numel(rx)+1 || mets(1)~=cache.start || ...
            numel(unique(mets))~=numel(mets) || numel(unique(rx))~=numel(rx),return;end
    if any(cache.currencyMask(mets)),reason="currency_connector";return;end
    if any(~cache.eligible(rx)),reason="ineligible_reaction";return;end
    if any(ismember(rx,downstreamRx)),reason="reaction_reused_between_blocks";return;end
    [interfaceOK,reason,blocked]=block0InterfaceMask(cache,mets(end),downstreamRx,inputRx);
    if ~interfaceOK,return;end
    seen=false(size(cache.S,1),1);
    for k=1:numel(rx)
        if cache.S(mets(k),rx(k))>=-1e-12 || cache.S(mets(k+1),rx(k))<=1e-12
            reason="invalid_forward_edge";return;
        end
        part=cache.participants{rx(k)};
        if any(blocked(part)),reason="noncurrency_overlap_with_downstream";return;end
        shared=part(seen(part));
        if k==1
            if ~isempty(shared),reason="repeated_nonadjacent_participant";return;end
        elseif numel(shared)~=1 || shared~=mets(k)
            reason="repeated_nonadjacent_or_side_participant";return;
        end
        seen(part)=true;
    end
    ok=true;reason="PASS";
end

function result=emptyBlock0Connection(target)
    result=struct('found',false,'mets',zeros(1,0),'rxns',zeros(1,0), ...
        'distance',Inf,'target',target,'downstreamRx',zeros(1,0), ...
        'inputRx',zeros(1,0),'shortestVerified',false,'expansions',0, ...
        'states',0,'expansionLimitReached',false,'stateLimitReached',false, ...
        'status',"no_isolated_path");
end

function result=findIsolatedBlock0Path(model,cache,target,downstreamRx,inputRx,preferred)
% Breadth-first search over PATH STATES. A visited-metabolite-only BFS is
% insufficient: two prefixes at the same metabolite can have different side
% participants and therefore different legal continuations.
    if nargin<6,preferred=[];end
    result=emptyBlock0Connection(target);
    downstreamRx=unique(reshape(downstreamRx,1,[]));
    inputRx=unique(reshape(inputRx,1,[]));
    result.downstreamRx=downstreamRx;result.inputRx=inputRx;
    options=struct('maxStates',Inf,'excludedECs',strings(0,1),'requireUniqueEC',false);
    if isfield(cache,'pathOptions'),options=cache.pathOptions;end
    [interfaceOK,why,blocked]=block0InterfaceMask(cache,target,downstreamRx,inputRx);
    if ~interfaceOK,result.status=why;return;end
    if options.requireUniqueEC && ~reactionsHaveUniqueECAssignment(model,downstreamRx)
        result.status="downstream_EC_assignment_conflict";return;
    end
    if ~cache.reachable(target)
        result.status="not_reached_by_structural_prefilter";return;
    end
    if isempty(preferred)
        [m,r]=block0CachedPath(cache,target);
        preferred=emptyBlock0Connection(target);preferred.mets=m;preferred.rxns=r;
    end
    [valid,~]=validateBlock0Isolation(cache,preferred.mets,preferred.rxns,downstreamRx,inputRx);
    inheritedProof=preferred.shortestVerified && preferred.target==target && ...
        all(ismember(preferred.downstreamRx,downstreamRx)) && ...
        isequal(preferred.inputRx,inputRx);
    if valid && preferred.mets(end)==target && numel(preferred.rxns)<=cache.maxDepth && ...
            (numel(preferred.rxns)==cache.depth(target) || inheritedProof) && ...
            ~reactionHasExcludedEC(model,[downstreamRx,preferred.rxns],options.excludedECs) && ...
            (~options.requireUniqueEC || reactionsHaveUniqueECAssignment(model,[downstreamRx,preferred.rxns]))
        result.found=true;result.mets=preferred.mets;result.rxns=preferred.rxns;
        result.distance=numel(preferred.rxns);result.shortestVerified=true;
        result.status="isolated_shortest_path";return;
    end
    eligible=cache.eligible;
    eligible(downstreamRx)=false;
    eligible=eligible & ~full(any(abs(cache.S(blocked,:))>1e-12,1))';
    % The final output cannot be a co-substrate of an upstream step, and the
    % source cannot be regenerated later in the connector.
    eligible=eligible & full(cache.S(target,:)>=-1e-12)' & ...
        full(cache.S(cache.start,:)<=1e-12)';
    for j=reshape(find(eligible),1,[])
        if reactionHasExcludedEC(model,j,options.excludedECs) || ...
                (options.requireUniqueEC && ~reactionsHaveUniqueECAssignment(model,[downstreamRx,j]))
            eligible(j)=false;
        end
    end
    if reactionHasExcludedEC(model,downstreamRx,options.excludedECs)
        result.status="downstream_excluded_EC";return;
    end
    % Reverse graph distances are only lower bounds; they prune impossible
    % branches without merging distinct path states.
    remaining=inf(size(cache.S,1),1);remaining(target)=0;
    q=target;head=1;
    while head<=numel(q)
        v=q(head);head=head+1;
        for j=reshape(cache.producers{v},1,[])
            if ~eligible(j),continue;end
            for pre=reshape(cache.substrates{j},1,[])
                if isfinite(remaining(pre)),continue;end
                remaining(pre)=remaining(v)+1;q(end+1)=pre; %#ok<AGROW>
            end
        end
    end
    if ~isfinite(remaining(cache.start)) || remaining(cache.start)>cache.maxDepth
        result.status="no_path_in_allowed_graph_within_depth";return;
    end
    stateMets={cache.start};stateRx={zeros(1,0)};stateUsed={zeros(1,0)};
    head=1;tail=1;result.states=1;
    while head<=tail
        mets=stateMets{head};rx=stateRx{head};used=stateUsed{head};
        stateMets{head}=[];stateRx{head}=[];stateUsed{head}=[];head=head+1;
        v=mets(end);depth=numel(rx);
        if depth>=cache.maxDepth,continue;end
        for j=reshape(cache.consumers{v},1,[])
            if ~eligible(j),continue;end
            if result.expansions>=cache.maxExpansions
                result.expansionLimitReached=true;result.status="path_expansion_budget_reached";return;
            end
            result.expansions=result.expansions+1;
            part=reshape(cache.participants{j},1,[]);
            shared=intersect(part,used);
            if depth>0 && (numel(shared)~=1 || shared~=v),continue;end
            newRx=[rx,j];
            if options.requireUniqueEC && ~reactionsHaveUniqueECAssignment(model,[downstreamRx,newRx]),continue;end
            nexts=reshape(cache.products{j},1,[]);
            if any(nexts==target),nexts=target;end
            for next=nexts
                if ismember(next,mets) || ~isfinite(remaining(next)) || ...
                        depth+1+remaining(next)>cache.maxDepth,continue;end
                newMets=[mets,next];
                if next==target
                    [valid,why]=validateBlock0Isolation(cache,newMets,newRx,downstreamRx,inputRx);
                    assert(valid,'Block0:IsolationSearchBug','Invalid constrained path: %s.',char(why));
                    result.found=true;result.mets=newMets;result.rxns=newRx;
                    result.distance=numel(newRx);result.shortestVerified=true;
                    result.status="isolated_shortest_path";return;
                end
                if tail>=options.maxStates
                    result.stateLimitReached=true;result.status="path_state_budget_reached";return;
                end
                tail=tail+1;stateMets{tail}=newMets;stateRx{tail}=newRx; %#ok<AGROW>
                stateUsed{tail}=union(used,part);result.states=tail; %#ok<AGROW>
            end
        end
    end
end

function result=memoizedBlock0Path(model,cache,target,downstreamRx,inputRx,preferred,memo)
% One result per exact constraint set; alternate connectors do not multiply
% motif counts. Workers create their own memo rather than sharing handles.
    key=sprintf('%d,',[target,-1,sort(unique(downstreamRx)),-2,sort(unique(inputRx))]);
    if isKey(memo,key),result=memo(key);return;end
    result=findIsolatedBlock0Path(model,cache,target,downstreamRx,inputRx,preferred);
    memo(key)=result;
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
    assert(~any(cache.currencyMask(mets)),'Block0:CurrencyPath', ...
        'Refusing a Block 0 path containing an excluded currency metabolite.');
end

function printBlock0(fid,model,metIDs,metNames,target,signalIDs,signalNames,cache,connection,downstreamRx,inputRx)
% Print precisely the path used to accept this motif; never rerun a search.
    assert(connection.found && connection.shortestVerified && connection.target==target, ...
        'Block0:UnverifiedPrintedPath','Block 0 output requires an accepted isolated path.');
    mets=connection.mets;rx=connection.rxns;
    [valid,why]=validateBlock0Isolation(cache,mets,rx,downstreamRx,inputRx);
    assert(valid,'Block0:PrintedCrossTalk','Refusing a Block 0 path with cross-talk: %s.',char(why));
    fprintf(fid,'  Anchor: %s | %s\n', ...
        char(metIDs(cache.start)),char(metNames(cache.start)));
    fprintf(fid,'  Connection status: CONNECTED (+1 only)\n');
    fprintf(fid,'  Currency check: PASS (no currency connector or m_in; co-participants retained).\n');
    fprintf(fid,'  Internal isolation: PASS (non-currency sharing only between adjacent path steps).\n');
    fprintf(fid,'  Downstream isolation: PASS (only m_in shared, consumed by r_Ff and r_Sf).\n');
    fprintf(fid,'  No Block 0 reaction is reused in another block.\n');
    if ~isempty(rx),fprintf(fid,'  m_in is produced only by the last Block 0 step and is never consumed within Block 0.\n');end
    fprintf(fid,'  Shared currency co-participants remain permitted; their buffering is not verified.\n');
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
    fprintf(fid,'  WARNING: shortest isolated structural path within the selected motif; atom transfer, additional substrates and flux feasibility unverified.\n');
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

function selfTestBlock0CurrencyPolicy()
% A false two-step cAMP -> proton -> f6p shortcut must lose to the real
% three-step route. Water/cofactor-only routes must remain unreachable.
    ids=["camp[c]";"amp[c]";"h[c]";"h2o[c]";"atp[c]";"adp[c]"; ...
        "nad[c]";"nadh[c]";"f6p[c]";"water_only[c]"; ...
        "atp_only[c]";"nad_only[c]";"real_mid[c]"];
    comps=repmat("c",numel(ids),1);ub=1000*ones(9,1);
    fixture.S=sparse(numel(ids),9);
    fixture.S([1 4 2 3],1)=[-1;-1;1;1]; % camp + water -> amp + proton
    fixture.S([3 9],2)=[-1;1];          % proton -> f6p (false link)
    fixture.S([4 10],3)=[-1;1];         % water-only link
    fixture.S([2 5 7 13 6 8],4)=[-1;-1;-1;1;1;1];
    fixture.S([6 11],5)=[-1;1];         % ADP-only link
    fixture.S([8 12],6)=[-1;1];         % NADH-only link
    fixture.S([13 3 9 4],7)=[-1;-1;1;1]; % real_mid + proton -> f6p + water
    fixture.S([13 5],8)=[-1;1];         % ATP must not become a connector
    fixture.S([13 7],9)=[-1;1];         % NAD must not become a connector
    c=buildBlock0Cache(fixture,ids,comps,ub,'camp[c]',Inf,Inf);
    assert(~any(c.reachable(3:8)) && all(c.currencyMask(3:8)), ...
        'Block0:SelfTestCurrencyNodes','A currency metabolite became reachable.');
    assert(~any(c.reachable(10:12)), ...
        'Block0:SelfTestCurrencyShortcut','A currency-only connection was accepted.');
    assert(c.reachable(2) && c.depth(2)==1, ...
        'Block0:SelfTestAMP','The valid cAMP -> AMP first step was removed.');
    [m,r]=block0CachedPath(c,9);
    assert(c.depth(9)==3 && isequal(m,[1 2 13 9]) && isequal(r,[1 4 7]), ...
        'Block0:SelfTestCurrencyDistance', ...
        'Shortest distance did not use the allowed route with currency co-participants.');
    unfiltered=buildBlock0Cache(fixture,ids,comps,ub,'camp[c]',Inf,Inf,strings(0,1));
    assert(unfiltered.depth(9)==2 && all(unfiltered.reachable(10:12)), ...
        'Block0:SelfTestFixture','The fixture does not reproduce currency shortcuts.');
    custom=buildBlock0Cache(fixture,ids,comps,ub,'camp[c]',Inf,Inf, ...
        [block0DefaultCurrencyIDs(),"amp"]);
    assert(nnz(custom.reachable)==1, ...
        'Block0:SelfTestCurrencyConfig','A configured extra currency was ignored.');
    aliases=block0CurrencyMask(["M_H2O_c";"H[c]";"ATP_c";"amp[c]"; ...
        "hco3[c]";"water[c]";"f6p[c]"],block0DefaultCurrencyIDs());
    assert(isequal(aliases,[true;true;true;false;false;true;false]), ...
        'Block0:SelfTestCurrencyIDs','Currency alias matching is too broad or incomplete.');
    rejected=false;
    try
        buildBlock0Cache(fixture,ids,comps,ub,'h[c]',Inf,Inf);
    catch ME
        rejected=strcmp(ME.identifier,'Block0:CurrencyAnchor');
    end
    assert(rejected,'Block0:SelfTestCurrencyAnchor','A currency anchor was accepted.');
    corrupted=c;corrupted.currencyMask(2)=true;rejected=false;
    try
        block0CachedPath(corrupted,9);
    catch ME
        rejected=strcmp(ME.identifier,'Block0:CurrencyPath');
    end
    assert(rejected,'Block0:SelfTestCurrencyGuard','Cached-path currency guard failed.');
    fprintf('PASS: Block 0 currency exclusion, retained AMP/co-participants, filtered distances and path guard.\n');
end

function selfTestBlock0Isolation()
% Covers the reported f6p co-substrate/recycling problem, hidden side-product
% sharing, downstream leaks, longer alternatives and path-state BFS behavior.
    ids=["camp[c]";"a[c]";"b[c]";"f6p[c]";"u[c]";"v[c]";"w[c]"; ...
        "x[c]";"fast[c]";"slow[c]";"h[c]";"h2o[c]";"side[c]"];
    fixture.S=sparse(numel(ids),13);fixture.rxns=cellstr("test_"+string((1:13)'));
    fixture.S([1 4 2],1)=[-1;-1;1]; % invalid: m_in consumed before it is made
    fixture.S([2 4],2)=[-1;1];
    fixture.S([1 12 3 11],3)=[-1;-1;1;1];
    fixture.S([3 5],4)=[-1;1];fixture.S([5 4],5)=[-1;1];
    fixture.S([4 9],6)=[-1;1];fixture.S([4 10],7)=[-1;1];
    fixture.S([1 6],8)=[-1;1];fixture.S([6 7],9)=[-1;1];
    fixture.S([7 8],10)=[-1;1];fixture.S([8 4],11)=[-1;1];
    fixture.S([3 13],12)=[-1;1];fixture.S([13 9],13)=[-1;1];
    comp=repmat("c",numel(ids),1);ub=1000*ones(13,1);
    cache=buildBlock0Cache(fixture,ids,comp,ub,'camp[c]',Inf,Inf);
    [rawM,rawR]=block0CachedPath(cache,4);
    [ok,~]=validateBlock0Isolation(cache,rawM,rawR,[6 7],[6 7]);
    assert(cache.depth(4)==2 && ~ok,'Block0:SelfTestPrematureInput', ...
        'The f6p recycling shortcut was not reproduced/rejected.');
    first=findIsolatedBlock0Path(fixture,cache,4,[6 7],[6 7],[]);
    assert(first.found && first.distance==3 && isequal(first.rxns,[3 4 5]), ...
        'Block0:SelfTestLongerAlternative','A valid longer isolated route was missed.');
    second=findIsolatedBlock0Path(fixture,cache,4,[6 7 12],[6 7],first);
    assert(second.found && second.distance==4 && isequal(second.rxns,[8 9 10 11]), ...
        'Block0:SelfTestDownstreamReroute','A downstream leak did not trigger an alternative path.');
    [ok,~]=validateBlock0Isolation(cache,first.mets,first.rxns,[3 6 7],[6 7]);
    assert(~ok,'Block0:SelfTestSharedReaction','A reaction reused across blocks was accepted.');
    leaking=cache;leaking.S(4,12)=1;
    rejected=findIsolatedBlock0Path(fixture,leaking,4,[6 7 12],[6 7],first);
    assert(~rejected.found && rejected.status=="m_in_used_outside_D1_input_reactions", ...
        'Block0:SelfTestInputLeak','A downstream producer of m_in was accepted.');
    limited=cache;limited.maxDepth=2;
    rejected=findIsolatedBlock0Path(fixture,limited,4,[6 7],[6 7],[]);
    assert(~rejected.found,'Block0:SelfTestIsolationDepth','The isolated path exceeded the depth cap.');
    limited=cache;limited.maxExpansions=0;
    rejected=findIsolatedBlock0Path(fixture,limited,4,[6 7],[6 7],[]);
    assert(~rejected.found && rejected.expansionLimitReached, ...
        'Block0:SelfTestIsolationBudget','Path expansion-budget failure was not reported.');
    limited=cache;limited.pathOptions=struct('maxStates',1, ...
        'excludedECs',strings(0,1),'requireUniqueEC',false);
    rejected=findIsolatedBlock0Path(fixture,limited,4,[6 7],[6 7],[]);
    assert(~rejected.found && rejected.stateLimitReached, ...
        'Block0:SelfTestIsolationStates','Path state-budget failure was not reported.');
    side=fixture;side.S=sparse(numel(ids),2);side.rxns={'side_1';'side_2'};
    side.S([1 13 2],1)=[-1;-1;1];side.S([2 4 13],2)=[-1;1;1];
    sc=buildBlock0Cache(side,ids,comp,1000*ones(2,1),'camp[c]',Inf,Inf);
    rejected=findIsolatedBlock0Path(side,sc,4,[],[],[]);
    assert(~rejected.found,'Block0:SelfTestSideMetabolite', ...
        'Cross-talk through a non-chain co-substrate/product was accepted.');
    % Two prefixes reach x. Only the second leaves side available for the last
    % reaction. A single global visited(x) flag would discard that solution.
    prefix=fixture;prefix.S=sparse(numel(ids),5);
    prefix.rxns=cellstr("prefix_"+string((1:5)'));
    prefix.S([1 13 2],1)=[-1;-1;1];prefix.S([2 8],2)=[-1;1];
    prefix.S([1 3],3)=[-1;1];prefix.S([3 8],4)=[-1;1];
    prefix.S([8 13 4],5)=[-1;-1;1];
    pc=buildBlock0Cache(prefix,ids,comp,1000*ones(5,1),'camp[c]',Inf,Inf);
    alternative=findIsolatedBlock0Path(prefix,pc,4,[],[],[]);
    assert(alternative.found && isequal(alternative.rxns,[3 4 5]), ...
        'Block0:SelfTestDistinctPrefixes','Distinct path states were incorrectly merged.');
    fprintf('PASS: Block 0 recycling, side-participant and downstream isolation; alternative paths, budgets and distinct prefixes.\n');
end

function printModelReaction(fid,model,idx,orient,metIDs,metNames,markedID,markedLabel,varargin)
% Print original stoichiometry; never relabel or reverse a stored reaction.
    isTransportRole=false;
    if numel(varargin)>=3,isTransportRole=isequal(varargin{3},true);end
    if isfield(model,'circuitECPolicyActive') && model.circuitECPolicyActive && ~isTransportRole
        assert(ecReactionAllowed(model,idx),'ECPolicy:MissingAnnotation', ...
            'Reaction %s has no complete EC annotation; only the r_T role is exempt.', ...
            char(string(model.rxns(idx))));
    end
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
    if isfield(model,'circuitResolvedEC')
        fprintf(fid,'  Resolved EC: %s\n',char(d1EC(model,idx)));
        fprintf(fid,'  EC annotation source: %s\n',char(model.circuitECSource(idx)));
        fprintf(fid,'  EC annotation reference: %s\n',char(model.circuitECReference(idx)));
        if isTransportRole && ~ecReactionAllowed(model,idx)
            fprintf(fid,'  EC requirement: exempt for r_T only; transporter GPR is reported separately.\n');
        end
    end
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
    ecAllowed=ecSearchMask(model);
    for j=1:numel(model.rxns)
        if ~ecAllowed(j),continue;end
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
% Rank: shortest m_Sd-to-counter path; best r_I evidence; largest r_I kcat; then cAMP distance.
    distances=zeros(numel(parents),1);
    if isempty(parents),return;end
    for k=1:numel(parents)
        if isfield(parents,'block0Path'),distances(k)=parents(k).block0Path.distance;
        else,distances(k)=depth(parents(k).d.min);end
    end
    assert(all(isfinite(distances)),'Ranking:DisconnectedParent','Only cAMP-connected parents may be ranked.');
    counterDistances=reshape(arrayfun(@(p)numel(p.exchangePathRxns),parents),[],1);
    evidence=zeros(numel(parents),1); kc=nan(numel(parents),1);
    if isfield(parents,'rIKcatEvidence'), evidence=reshape([parents.rIKcatEvidence],[],1); end
    if isfield(parents,'rIKcat'), kc=reshape([parents.rIKcat],[],1); end
    known=isfinite(kc)&kc>0; ksort=kc; ksort(~known)=-Inf;
    serial=(1:numel(parents))';
    [~,order]=sortrows([counterDistances,-evidence,-double(known),-ksort,distances,serial],[1 2 3 4 5 6]);
    parents=parents(order);distances=distances(order);
end

function motifs=keepTopCompleteMotifs(motifs,limit)
    if isempty(motifs),return;end
    counterDistances=reshape([motifs.counterDistance],[],1); distances=reshape([motifs.distance],[],1);
    evidence=zeros(numel(motifs),1); kc=nan(numel(motifs),1);
    if isfield(motifs,'rIKcatEvidence'), evidence=reshape([motifs.rIKcatEvidence],[],1); end
    if isfield(motifs,'rIKcat'), kc=reshape([motifs.rIKcat],[],1); end
    known=isfinite(kc)&kc>0; ksort=kc; ksort(~known)=-Inf; serial=reshape([motifs.discoveryIndex],[],1);
    [~,order]=sortrows([counterDistances,-evidence,-double(known),-ksort,distances,serial],[1 2 3 4 5 6]);
    motifs=motifs(order); if numel(motifs)>limit,motifs=motifs(1:limit);end
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
function [parents,K]=attachRIKcat(parents,model,metIDs,metNames,organism,allowAnyOrg,allowAnySub,allowSAMW,flatFile,cacheFile)
% Local BRENDA kinetics for r_I. Evidence rank: 3 substrate+organism matched;
% 2 substrate matched with organism fallback / any-substrate same organism;
% 1 EC-level turnover or SA*MW estimate; 0 unresolved.
    n=numel(parents); rxn=strings(n,1); ec=strings(n,1); substrate=strings(n,1);
    for i=1:n
        j=parents(i).d.ri; rxn(i)=string(model.rxns{j}); ec(i)=d1EC(model,j);
        subs=find(full(model.S(:,j))<0);
        % r_I consumes m_Fd and m_Sd; use m_Sd as the primary tagged canceller substrate.
        target=parents(i).target;
        if any(subs==target), substrate(i)=metNames(target); elseif ~isempty(subs), substrate(i)=metNames(subs(1)); end
    end
    fullEC=strings(0,1);
    for i=1:n, fullEC=[fullEC;string(regexp(char(ec(i)),'\d+\.\d+\.\d+\.\d+','match'))']; end %#ok<AGROW>
    fullEC=unique(fullEC,'stable');
    db=loadRIFlat(flatFile,fullEC,cacheFile);
    kval=nan(n,1); evidence=zeros(n,1); source=strings(n,1); usedEC=strings(n,1); usedSub=strings(n,1); usedOrg=strings(n,1);
    for i=1:n
        codes=string(regexp(char(ec(i)),'\d+\.\d+\.\d+\.\d+','match'));
        best=[-Inf -Inf]; bestVal=NaN; bestSource=""; bestEC=""; bestSub=""; bestOrg="";
        for code=codes
            if ~isKey(db,char(code)),continue;end
            d=db(char(code)); tn=d.TN;
            if ~isempty(tn)
                vals=[tn.value]'; ok=isfinite(vals)&vals>0; tn=tn(ok);
                if ~isempty(tn)
                    org=string({tn.organism})'; sub=string({tn.substrate})';
                    exactSub=normRI(sub)==normRI(substrate(i)); exactOrg=normRI(org)==normRI(string(organism));
                    masks={exactSub&exactOrg, exactSub, exactOrg, true(size(exactSub))}; ranks=[3 2 2 1]; labels=["BRENDA turnover: substrate + organism matched","BRENDA turnover: substrate matched; organism fallback","BRENDA turnover: organism matched; substrate fallback","BRENDA turnover: EC-level any substrate/organism"];
                    for q=1:4
                        if q==2 && ~allowAnyOrg,continue;end
                        if q>=3 && ~allowAnySub,continue;end
                        m=masks{q}; if any(m)
                            v=median([tn(m).value]); score=[ranks(q),v];
                            if score(1)>best(1) || (score(1)==best(1)&&score(2)>best(2))
                                best=score;bestVal=v;bestSource=labels(q);bestEC=code;bestSub=strjoin(unique(sub(m),'stable'),'; ');bestOrg=strjoin(unique(org(m),'stable'),'; ');
                            end
                            break
                        end
                    end
                end
            end
            if ~isfinite(bestVal) && allowSAMW
                sa=[d.SA.value]; mw=[d.MW.value]; sa=sa(isfinite(sa)&sa>0); mw=mw(isfinite(mw)&mw>0);
                if ~isempty(sa)&&~isempty(mw)
                    v=median(sa)*median(mw)/60/1000; % SA umol/min/mg * MW Da -> s^-1
                    best=[1 v];bestVal=v;bestSource="ESTIMATED SA*MW/60; EC-level, mass basis requires review";bestEC=code;
                end
            end
        end
        kval(i)=bestVal;evidence(i)=max(0,best(1));source(i)=bestSource;usedEC(i)=bestEC;usedSub(i)=bestSub;usedOrg(i)=bestOrg;
        parents(i).rIKcat=kval(i);parents(i).rIKcatEvidence=evidence(i);parents(i).rIKcatSource=source(i);
    end
    K=table(rxn,ec,substrate,kval,evidence,source,usedEC,usedSub,usedOrg,'VariableNames', ...
        {'r_I_ReactionID','ModelEC','PrimarySubstrate','kcat_s','EvidenceRank','Source','ECused','SubstrateUsed','OrganismUsed'});
end

function db=loadRIFlat(flatFile,wantedECs,cacheFile)
% Small independent BRENDA flat-file index for the r_I ECs only.
    wantedECs=unique(string(wantedECs)); wantedECs=wantedECs(wantedECs~="");
    db=containers.Map('KeyType','char','ValueType','any'); sig=dir(flatFile); signature=[sig.bytes sig.datenum];
    if isfile(cacheFile)
        try, S=load(cacheFile,'dbCache','signatureCache'); if isequal(S.signatureCache,signature),db=S.dbCache;end, catch, end
    end
    missing=wantedECs(~cellfun(@(x)isKey(db,x),cellstr(wantedECs))); if isempty(missing),return;end
    want=containers.Map(cellstr(missing),num2cell(true(size(missing))));
    fid=fopen(flatFile,'r','n','UTF-8'); assert(fid>0,'Cannot open BRENDA flat file: %s',flatFile); c=onCleanup(@()fclose(fid));
    current=""; keep=false; orgMap=containers.Map('KeyType','char','ValueType','char');
    while true
        line=fgetl(fid); if ~ischar(line),break;end
        id=regexp(line,'^ID\s+([0-9]+\.[0-9]+\.[0-9]+\.[0-9-]+)','tokens','once');
        if ~isempty(id), current=string(id{1});keep=isKey(want,char(current));orgMap=containers.Map('KeyType','char','ValueType','char'); if keep&&~isKey(db,char(current)),db(char(current))=struct('TN',emptyRIRec(),'SA',emptyRIRec(),'MW',emptyRIRec());end;continue;end
        if ~keep,continue;end
        pr=regexp(line,'^PR\s+#(\d+)#\s*([^<{#]*)','tokens','once'); if ~isempty(pr),orgMap(strtrim(pr{1}))=strtrim(pr{2});continue;end
        hit=regexp(line,'^(TN|SA|MW)\s+(.+)$','tokens','once'); if isempty(hit),continue;end
        tag=hit{1}; body=hit{2}; pm=regexp(body,'^#([\d,\s]*)#','tokens','once'); protein=''; if ~isempty(pm),pp=strsplit(strtrim(pm{1}),',');protein=strtrim(pp{1});end
        rest=regexprep(body,'^#[^#]*#\s*',''); vm=regexp(rest,'^\s*(-?[\d.]+(?:[eE][+-]?\d+)?)','tokens','once'); if isempty(vm),continue;end
        v=str2double(vm{1}); if ~isfinite(v)||v<=0,continue;end
        sm=regexp(rest,'\{([^}]*)\}','tokens','once'); sub='';if ~isempty(sm),sub=strtrim(sm{1});end
        org='';if ~isempty(protein)&&isKey(orgMap,protein),org=orgMap(protein);end
        d=db(char(current)); rec=struct('value',v,'substrate',sub,'organism',org);d.(tag)(end+1)=rec;db(char(current))=d;
    end
    dbCache=db;signatureCache=signature; %#ok<NASGU>
    try,save(cacheFile,'dbCache','signatureCache','-v7');catch,end
end
function s=emptyRIRec(),s=struct('value',{},'substrate',{},'organism',{});end
function s=normRI(s),s=lower(strtrim(string(s)));s=regexprep(s,'\s+',' ');end


function mappings=readECMappingTable(mappingPath)
% EC codes are identifiers, never numeric values or dates. Force every column
% to string BEFORE reading; TextType='string' alone still permits type inference.
    [~,~,extension]=fileparts(mappingPath);delimiter=sprintf('\t');
    if strcmpi(extension,'.csv'),delimiter=',';end
    opts=detectImportOptions(mappingPath,'FileType','text','Delimiter',delimiter, ...
        'VariableNamingRule','preserve');
    opts.VariableNamesLine=1;
    opts.DataLines=[2 Inf];
    opts=setvartype(opts,opts.VariableNames,'string');
    mappings=readtable(mappingPath,opts);
    required={'ReactionID','EC','Source','Reference'};
    assert(all(ismember(required,mappings.Properties.VariableNames)), ...
        'ECPolicy:MappingColumns', ...
        'EC evidence file needs ReactionID, EC, Source, and Reference columns.');
    for k=1:numel(required)
        mappings.(required{k})=strtrim(mappings.(required{k}));
    end
end

function [model,audit,info]=resolveRequiredECAnnotations(model,comp,ub,root,options)
% Resolve exact reaction annotations once. No enzyme-name similarity matching.
    n=numel(model.rxns);ids=string(model.rxns(:));names=ids;
    if isfield(model,'rxnNames'),names=string(model.rxnNames(:));end
    originalEC=strings(n,1);resolved=strings(n,1);source=strings(n,1);
    reference=strings(n,1);status=repmat("unresolved",n,1);
    for j=1:n
        [tokens,fields]=modelECTokens(model,j);
        if ~isempty(tokens)
            originalEC(j)=join(tokens,'; ');
            resolved(j)=originalEC(j);source(j)=join(fields,'; ');
            reference(j)="model reaction "+ids(j);status(j)="model_annotation";
        end
    end
    S=model.S;tol=1e-12;
    candidate=ub(:)>0 & full(any(S < -tol,1))' & full(any(S > tol,1))';
    candidate=candidate & ~full(any(abs(S(comp~="c",:))>tol,1))';
    mappingPath=which(options.mappingFile);
    if isempty(mappingPath),mappingPath=findNamed(root,options.mappingFile);end
    mappings=table();
    if ~isempty(mappingPath)
        mappings=readECMappingTable(mappingPath);
    end
    brendaPath=which(options.brendaFile);
    if isempty(brendaPath),brendaPath=findNamed(root,options.brendaFile);end
    wanted=strings(0,1);
    if ~isempty(mappings)
        for k=1:height(mappings)
            if strcmpi(strtrim(string(mappings.Source(k))),'BRENDA')
                wanted=[wanted;completeECTokens(mappings.EC(k))]; %#ok<AGROW>
            end
        end
    end
    brendaEntries=brendaECEntries(brendaPath,unique(wanted));
    queries=0;failures=0;stopped=false;
    if ~isempty(mappings)
        for k=1:height(mappings)
            j=find(ids==strtrim(string(mappings.ReactionID(k))),1);
            if isempty(j) || resolved(j)~="",continue;end
            codes=completeECTokens(mappings.EC(k));
            evidenceSource=strtrim(string(mappings.Source(k)));
            ref=strtrim(string(mappings.Reference(k)));
            if isempty(codes) || ismissing(ref) || ref=="",continue;end
            if strcmpi(evidenceSource,'BRENDA')
                valid=false(size(codes));
                for q=1:numel(codes),valid(q)=isKey(brendaEntries,char(codes(q)));end
                codes=codes(valid);
                if isempty(codes),status(j)="BRENDA_mapping_not_confirmed_in_local_file";continue;end
                source(j)="BRENDA curated reaction mapping; local ID entry confirmed";
                reference(j)=ref+" | "+string(brendaPath);
            elseif strcmpi(evidenceSource,'BiGG')
                % A user-supplied exact reaction-ID mapping is evidence, not a
                % name guess. Require a BiGG source URL for auditability.
                if isempty(regexpi(char(ref),'^https?://(www\.)?(bigg\.bio|bigg\.ucsd\.edu|biggdata\.ucsd\.edu)/','once'))
                    status(j)="BiGG_mapping_missing_BiGG_source_URL";continue;
                end
                if ~options.online || stopped || queries>=options.maxLookups
                    status(j)="BiGG_mapping_could_not_be_verified";continue;
                end
                queries=queries+1;
                [confirmed,canonicalRef,responded]=lookupBiGGEC(ids(j),options.timeout);
                if responded,failures=0;else,failures=failures+1;end
                if failures>=options.maxFailures,stopped=true;end
                codes=intersect(codes,confirmed,'stable');
                if isempty(codes)
                    status(j)="BiGG_mapping_not_confirmed_in_database_response";continue;
                end
                source(j)="BiGG exact reaction-ID annotation; supplied mapping confirmed";
                reference(j)=canonicalRef+" | "+ref;
            else
                status(j)="mapping_source_not_BiGG_or_BRENDA";continue;
            end
            resolved(j)=join(codes,'; ');status(j)="external_sourced_mapping";
        end
    end
    for j=reshape(find(candidate & resolved==""),1,[])
        if ~options.online
            status(j)=status(j)+";online_lookup_disabled";continue;
        end
        if queries>=options.maxLookups
            status(j)="unresolved_lookup_budget";continue;
        end
        if stopped
            status(j)="unresolved_BiGG_service_unavailable";continue;
        end
        queries=queries+1;
        [codes,ref,serviceResponded]=lookupBiGGEC(ids(j),options.timeout);
        if serviceResponded,failures=0;else,failures=failures+1;end
        if ~isempty(codes)
            resolved(j)=join(codes,'; ');source(j)="BiGG exact reaction-ID annotation";
            reference(j)=ref;status(j)="BiGG_annotation";
        elseif serviceResponded
            status(j)="unresolved_BiGG_no_complete_EC";
        else
            status(j)="unresolved_BiGG_lookup_failed";
        end
        if failures>=options.maxFailures,stopped=true;end
        if mod(queries,50)==0
            fprintf('EC annotation resolution: %d BiGG lookups; %d annotated intracellular reactions.\n', ...
                queries,nnz(candidate & resolved~=""));
        end
    end
    model.circuitResolvedEC=resolved;
    model.circuitECSource=source;
    model.circuitECReference=reference;
    model.circuitECEligible=resolved~="";
    model.circuitECPolicyActive=true;
    audit=table(ids,names,originalEC,resolved,source,reference,status, ...
        model.circuitECEligible,candidate, ...
        'VariableNames',{'ReactionID','ReactionName','ModelEC','ResolvedEC', ...
        'AnnotationSource','Reference','ResolutionStatus','HasCompleteEC','CandidateNonTransport'});
    info=sprintf('EC resolution: %d model-annotated reactions; %d BiGG requests; circuit breaker=%d. Unresolved means excluded, not proven absent from every database.', ...
        nnz(originalEC~=""),queries,stopped);
end

function [codes,fields]=modelECTokens(model,j)
% Union all known EC fields, including a populated fallback after an empty one.
    codes=strings(0,1);fields=strings(0,1);
    names={'rxnECNumbers','rxnEC','eccodes','ecCodes','ecNumbers','rxnECNumber', ...
        'rxnECcodes','rxnEnzymeClassification','circuitSBMLEC'};
    for q=1:numel(names)
        if ~isfield(model,names{q}),continue;end
        value=ecAnnotationAt(model.(names{q}),j,numel(model.rxns));
        found=completeECTokens(value);
        if ~isempty(found)
            codes=[codes;found];fields(end+1,1)="model."+string(names{q}); %#ok<AGROW>
        end
    end
    if isfield(model,'rxnMiriams')
        node=ecAnnotationAt(model.rxnMiriams,j,numel(model.rxns));
        found=namedECAnnotations(node);
        if ~isempty(found)
            codes=[codes;found];fields(end+1,1)="model.rxnMiriams"; %#ok<AGROW>
        end
    end
    codes=unique(codes,'stable');fields=unique(fields,'stable');
end

function value=ecAnnotationAt(values,j,n)
    value=[];
    if iscell(values)
        if numel(values)==n,value=values{j};end
    elseif isstring(values) || isstruct(values)
        if numel(values)==n,value=values(j);end
    elseif ischar(values)
        if n==1,value=values;
        elseif size(values,1)==n,value=values(j,:);end
    end
end

function codes=completeECTokens(value)
% Only full four-field ECs qualify. Partial, blank and non-EC labels do not.
    codes=strings(0,1);
    if iscell(value)
        for k=1:numel(value),codes=[codes;completeECTokens(value{k})];end %#ok<AGROW>
    elseif isstruct(value)
        fields=fieldnames(value);
        for k=1:numel(value)
            for q=1:numel(fields)
                codes=[codes;completeECTokens(value(k).(fields{q}))]; %#ok<AGROW>
            end
        end
    elseif ischar(value) || isstring(value)
        text=string(value);text=text(~ismissing(text));
        if isempty(text),return;end
        text=char(join(text(:),' '));
        hit=regexp(text,'(?<![0-9.])[1-7]\.[0-9]+\.[0-9]+\.[0-9]+(?![0-9.\-])','match');
        codes=string(hit(:));
    end
    codes=unique(codes,'stable');
end

function codes=namedECAnnotations(node)
% Parse only EC-labelled namespaces, keys, or explicit EC-code URIs.
% A bare number elsewhere in a database payload is NOT an EC annotation.
    codes=strings(0,1);
    if iscell(node)
        for k=1:numel(node),codes=[codes;namedECAnnotations(node{k})];end %#ok<AGROW>
    elseif isstruct(node)
        names=fieldnames(node);
        for k=1:numel(node)
            % COBRA MIRIAM: name='ec-code', ids={'1.2.3.4'}.
            for nsName={'name','namespace','database','data_source'}
                if isfield(node,nsName{1})
                    ns=node(k).(nsName{1});
                    if (ischar(ns) || isstring(ns)) && isECNamespace(ns)
                        for idName={'ids','id','identifier','identifiers','accession','value'}
                            if isfield(node,idName{1})
                                codes=[codes;completeECTokens(node(k).(idName{1}))]; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
            for q=1:numel(names)
                val=node(k).(names{q});
                if isECNamespace(names{q})
                    codes=[codes;completeECTokens(val)]; %#ok<AGROW>
                else
                    codes=[codes;namedECAnnotations(val)]; %#ok<AGROW>
                end
            end
        end
    elseif ischar(node) || isstring(node)
        vals=string(node);vals=vals(~ismissing(vals));
        for k=1:numel(vals)
            text=char(vals(k));
            if ~isempty(regexpi(text,'(ec-code|ec_code|enzyme-commission|(^|\s)EC:)[/: ]','once'))
                codes=[codes;completeECTokens(text)]; %#ok<AGROW>
            end
        end
    end
    codes=unique(codes,'stable');
end

function tf=isECNamespace(value)
    text=lower(string(value));text=regexprep(text,'[^a-z0-9]','');
    tf=isscalar(text) && any(text==["ec","eccode","eccodes","ecnumber","ecnumbers", ...
        "enzymecommission","enzymecommissionnumber","enzymeclassification", ...
        "ecx0x2dcode","ec0x2dcode"]);
end

function entries=brendaECEntries(path,wanted)
% Validate explicitly sourced BRENDA mappings against ID records, not names.
    entries=containers.Map('KeyType','char','ValueType','logical');
    if isempty(path) || isempty(wanted),return;end
    fid=fopen(path,'r','n','UTF-8');
    if fid<0,return;end
    closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
    wanted=unique(wanted);
    while true
        line=fgetl(fid);if ~ischar(line),break;end
        hit=regexp(line,'^ID\s+([1-7]\.[0-9]+\.[0-9]+\.[0-9]+)(?![0-9.\-])','tokens','once');
        if ~isempty(hit) && any(wanted==string(hit{1}))
            entries(hit{1})=true;
            if entries.Count==numel(wanted),break;end
        end
    end
end

function [codes,reference,responded]=lookupBiGGEC(reactionID,timeout)
% Exact ID requests only. Prefer the current v3 object API, then legacy v2.
% Missing EC fields, errors and unknown payload shapes all fail closed.
    codes=strings(0,1);reference="";responded=false;
    id=char(string(reactionID));
    if isempty(regexp(id,'^[A-Za-z0-9_\-]+$','once')),return;end
    try
        endpoint='https://bigg.bio/api/v3/objects';
        record=webwrite(endpoint,struct('type','universalreaction','id',id), ...
            weboptions('Timeout',timeout,'MediaType','application/json','ContentType','json'));
        if isfield(record,'data') && isstruct(record.data),record=record.data;end
        if isstruct(record) && isscalar(record) && isfield(record,'bigg_id') && ...
                strcmp(string(record.bigg_id),string(reactionID))
            responded=true;codes=namedECAnnotations(record);
            if ~isempty(codes)
                reference="https://bigg.bio/universal/reactions/"+string(reactionID);return;
            end
        end
    catch ME
        if contains(ME.identifier,'HTTP404') || contains(ME.message,'404')
            responded=true; % Unknown reaction is not a network outage.
        end
        % Try the legacy exact-ID endpoint; never infer an EC from an error.
    end
    try
        endpoint=['https://bigg.ucsd.edu/api/v2/universal/reactions/',id];
        record=webread(endpoint,weboptions('Timeout',timeout,'ContentType','json'));
        if isstruct(record) && isscalar(record) && isfield(record,'bigg_id') && ...
                strcmp(string(record.bigg_id),string(reactionID))
            responded=true;codes=namedECAnnotations(record);
            if ~isempty(codes),reference=string(endpoint);end
        end
    catch ME
        if contains(ME.identifier,'HTTP404') || contains(ME.message,'404')
            responded=true;
        end
    end
end

function mask=ecSearchMask(model)
% Synthetic legacy fixtures without an active policy retain their old behavior.
    if isfield(model,'S'),n=size(model.S,2);else,n=numel(model.rxns);end
    mask=true(n,1);
    if isfield(model,'circuitECPolicyActive') && model.circuitECPolicyActive
        mask=logical(model.circuitECEligible(:));
    end
end

function tf=ecReactionAllowed(model,j)
    mask=ecSearchMask(model);tf=all(mask(j));
end

function tf=reactionsHaveRequiredEC(model,reactionIDs)
    tf=true;
    for j=reshape(unique(reactionIDs),1,[])
        if isempty(completeECTokens(d1EC(model,j))),tf=false;return;end
    end
end

function selfTestRequiredECPolicy()
% No external calls. In particular, a missing r_cr must never be admitted.
    fixture.rxns={'annotated';'fallback';'partial';'MCLACCYSR';'transporter'};
    fixture.rxnECNumbers={'1.1.1.1';'';'2.7.1.-';'';''};
    fixture.rxnEC={'';'2.7.1.1';'';'';''};
    [codes,~]=modelECTokens(fixture,2);
    assert(isequal(codes,"2.7.1.1"),'EC self-test: empty primary field hid a valid fallback.');
    assert(isempty(completeECTokens('2.7.1.-; NaN')), ...
        'EC self-test: partial EC was admitted.');
    record=struct('database_links',struct('ec_code', ...
        struct('id','1.2.3.4','link','https://identifiers.org/ec-code/1.2.3.4')));
    assert(isequal(namedECAnnotations(record),"1.2.3.4"), ...
        'EC self-test: BiGG EC namespace was not parsed.');
    miriam=struct('name','ec-code','ids',{{'2.3.4.5'}});
    assert(isequal(namedECAnnotations(miriam),"2.3.4.5"), ...
        'EC self-test: model MIRIAM EC namespace was not parsed.');
    assert(isempty(namedECAnnotations(struct('unrelated_id','1.2.3.4'))), ...
        'EC self-test: unrelated numeric identifier was treated as an EC.');
    fixture.circuitResolvedEC=["1.1.1.1";"2.7.1.1";"";"";""];
    fixture.circuitECEligible=fixture.circuitResolvedEC~="";
    fixture.circuitECPolicyActive=true;
    assert(~ecReactionAllowed(fixture,4) && ~reactionsHaveRequiredEC(fixture,[1 4]), ...
        'EC self-test: missing MCLACCYSR EC passed the required-EC filter.');
    assert(reactionsHaveRequiredEC(fixture,[1 2]), ...
        'EC self-test: annotated nontransport reactions were rejected.');
    % The required list excludes the r_T role only; using that unannotated
    % reaction in any other role is still forbidden.
    assert(~reactionsHaveRequiredEC(fixture,[1 2 5]), ...
        'EC self-test: exemption leaked from a role to a reaction ID.');
    folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
    fixture.rxns{6,1}='external_test';
    fixture.rxnECNumbers{6,1}='';fixture.rxnEC{6,1}='';
    fixture.S=sparse([-ones(1,6);ones(1,6)]);
    fixture.rxnNames=fixture.rxns;
    brenda=fullfile(folder,'brenda_test.txt');fid=fopen(brenda,'w');
    assert(fid>0,'EC self-test: cannot create temporary BRENDA fixture.');
    fprintf(fid,'ID\t1.2.3.4\nRN\tSynthetic enzyme\n');fclose(fid);
    mapping=fullfile(folder,'mapping.tsv');fid=fopen(mapping,'w');
    assert(fid>0,'EC self-test: cannot create temporary mapping fixture.');
    fprintf(fid,'ReactionID\tEC\tSource\tReference\n');
    fprintf(fid,'external_test\t1.2.3.4\tBRENDA\tSynthetic test source\n');fclose(fid);
    foundMapping=findNamed(folder,mapping);foundBrenda=findNamed(folder,brenda);
    assert(~isempty(foundMapping) && isfile(foundMapping), ...
        'ECPolicy:SelfTestMappingPath','Explicit temporary mapping path was not resolved.');
    assert(~isempty(foundBrenda) && isfile(foundBrenda), ...
        'ECPolicy:SelfTestBRENDAPath','Explicit temporary BRENDA path was not resolved.');
    imported=readECMappingTable(foundMapping);
    assert(height(imported)==1 && isstring(imported.EC) && ...
        imported.ReactionID(1)=="external_test" && imported.EC(1)=="1.2.3.4" && ...
        imported.Source(1)=="BRENDA" && imported.Reference(1)=="Synthetic test source", ...
        'ECPolicy:SelfTestMappingText','TSV mapping did not preserve the original reaction/EC/source text.');
    csvMapping=fullfile(folder,'mapping with spaces.csv');fid=fopen(csvMapping,'w');
    assert(fid>0,'EC self-test: cannot create CSV mapping fixture.');
    fprintf(fid,'ReactionID,EC,Source,Reference\n');
    fprintf(fid,'external_test,"1.2.3.4",BRENDA,"Synthetic, test source"\n');fclose(fid);
    importedCSV=readECMappingTable(csvMapping);
    assert(height(importedCSV)==1 && importedCSV.EC(1)=="1.2.3.4" && ...
        importedCSV.Reference(1)=="Synthetic, test source", ...
        'ECPolicy:SelfTestCSVMapping','CSV mapping did not preserve quoted EC/reference text.');
    options=struct('online',false,'timeout',1,'maxLookups',0,'maxFailures',1, ...
        'mappingFile',mapping,'brendaFile',brenda);
    [resolved,audit,~]=resolveRequiredECAnnotations(fixture,["c";"c"],ones(6,1),folder,options);
    assert(resolved.circuitECEligible(1),'ECPolicy:SelfTestPrimaryEC', ...
        'Model EC was lost: status=%s; EC=%s.', ...
        char(audit.ResolutionStatus(1)),char(audit.ResolvedEC(1)));
    assert(resolved.circuitECEligible(2),'ECPolicy:SelfTestFallbackEC', ...
        'Fallback model EC was lost: status=%s; EC=%s.', ...
        char(audit.ResolutionStatus(2)),char(audit.ResolvedEC(2)));
    assert(resolved.circuitECEligible(6),'ECPolicy:SelfTestBRENDAImport', ...
        'Synthetic BRENDA mapping was not loaded: status=%s; EC=%s; mapping=%s; BRENDA=%s.', ...
        char(audit.ResolutionStatus(6)),char(audit.ResolvedEC(6)),mapping,brenda);
    assert(~resolved.circuitECEligible(4),'ECPolicy:SelfTestUnknownReaction', ...
        'Unannotated MCLACCYSR was incorrectly accepted: status=%s; EC=%s.', ...
        char(audit.ResolutionStatus(4)),char(audit.ResolvedEC(4)));
    assert(audit.ResolvedEC(6)=="1.2.3.4" && contains(audit.AnnotationSource(6),"BRENDA"), ...
        'EC self-test: BRENDA mapping source was not recorded.');
    [pool,~,~,~]=d1BuildOriented(resolved,["c";"c"],zeros(6,1),ones(6,1),true);
    assert(isequal(pool,[1;2;6]),'EC self-test: unannotated reaction survived pool filtering.');
    % Printing provides an independent guard. The SAME unannotated model
    % reaction is forbidden in r_cr and allowed only in an explicit r_T role.
    probe=fopen(fullfile(folder,'reaction_role_test.txt'),'w');
    assert(probe>0,'ECPolicy:SelfTestReport','Cannot create temporary role report.');
    closeProbe=onCleanup(@()fclose(probe)); %#ok<NASGU>
    rejected=false;
    try
        printModelReaction(probe,resolved,4,1,["camp[c]";"amp[c]"], ...
            ["cAMP";"AMP"],"","");
    catch ME
        rejected=strcmp(ME.identifier,'ECPolicy:MissingAnnotation');
    end
    assert(rejected,'ECPolicy:SelfTestCounterPrint', ...
        'An unannotated r_cr reaction passed the final printing guard.');
    printModelReaction(probe,resolved,4,1,["camp[c]";"amp[c]"], ...
        ["cAMP";"AMP"],"","",[],[],true);
    clear closeProbe;
    % A missing-EC direct connector must not hide an annotated longer path.
    connector.S=sparse([-1 -1 0;1 0 -1;0 1 1]);
    connector.rxns={'first';'unannotated_shortcut';'second'};
    connector.rxnECNumbers={'3.1.4.17';'';'3.5.4.6'};
    connector.circuitResolvedEC=["3.1.4.17";"";"3.5.4.6"];
    connector.circuitECEligible=[true;false;true];
    connector.circuitECPolicyActive=true;
    c=buildBlock0Cache(connector,["camp[c]";"amp[c]";"imp[c]"], ...
        repmat("c",3,1),ones(3,1),'camp[c]',Inf,Inf);
    [~,connectorRx]=block0CachedPath(c,3);
    assert(c.depth(3)==2 && isequal(connectorRx,[1 3]), ...
        'ECPolicy:SelfTestBlock0','An unannotated shortcut was selected in Block 0.');
    fprintf('PASS: complete EC presence, fallback annotations, namespace parsing and role-only transport exemption.\n');
end

function ec=d1EC(model,j)
    if isfield(model,'circuitResolvedEC')
        ec=string(model.circuitResolvedEC(j));
    else
        [ecs,~]=modelECTokens(model,j);
        ec="";
        if ~isempty(ecs),ec=join(ecs,'; ');end
    end
    if ismissing(ec),ec="";end
end

function tf=reactionsHaveUniqueECAssignment(model,reactionIDs)
% True when all reactions carrying complete EC annotations can be assigned
% distinct EC numbers. This is a bipartite matching problem because a model
% reaction may carry several EC annotations. Repeated reaction indices are
% collapsed: one chemical reaction appearing twice in the assembled list is
% still one reaction, not two enzymes.
% EC presence is enforced separately, even if this uniqueness check is off.
% Blank annotations may occur here only for the independently exempt r_T role.
    reactionIDs=unique(reshape(reactionIDs,1,[]),'stable');
    ecSets=cell(numel(reactionIDs),1);
    constrained=false(numel(reactionIDs),1);
    for k=1:numel(reactionIDs)
        raw=d1EC(model,reactionIDs(k));
        if raw=="",continue;end
        % Only complete four-field ECs participate. Partial codes such as
        % 2.6.1.- and blank annotations do not prove enzyme identity.
        ecs=string(regexp(char(raw),'(?<![\d.])\d+\.\d+\.\d+\.\d+(?![\d.])','match'));
        ecs=unique(ecs,'stable');
        if isempty(ecs),continue;end
        ecSets{k}=ecs(:)';constrained(k)=true;
    end
    ecSets=ecSets(constrained);
    if numel(ecSets)<=1,tf=true;return;end
    % Fail fast when there are fewer distinct EC labels than constrained reactions.
    allEC=unique([ecSets{:}]);
    if numel(allEC)<numel(ecSets),tf=false;return;end
    % Hardest reactions first improves the small backtracking search.
    [~,ord]=sort(cellfun(@numel,ecSets));ecSets=ecSets(ord);
    tf=assignDistinctEC(ecSets,1,strings(0,1));
end

function tf=assignDistinctEC(ecSets,k,used)
% Depth-first exact matching; motif reaction counts are small, so this is
% deterministic and cheaper than introducing graph-toolbox dependencies.
    if k>numel(ecSets),tf=true;return;end
    choices=ecSets{k};
    choices=choices(~ismember(choices,used));
    tf=false;
    for q=1:numel(choices)
        if assignDistinctEC(ecSets,k+1,[used;choices(q)])
            tf=true;return;
        end
    end
end

function R=searchStage8Parent(ii,selectedD1,T,metIDs,comp,currency,cons8,di,do, ...
        model,block0Cache,excludedECs,requireUniqueEC,maxBranch,maxTests)
% Pure per-parent Stage 8 search. Safe for parfor: no file I/O, no shared maps.
    R=struct('examined',0,'tests',0,'capReached',0, ...
        'rejectDuplicateReaction',0,'rejectD1Overlap',0,'rejectExcludedEC',0, ...
        'rejectRepeatedEC',0,'rejectMissingEC',0,'rejectBlock0Isolation',0, ...
        'block0PathCapped',0,'candidates',struct('mc',{},'nf',{},'nd',{}, ...
        'nf2',{},'mout',{},'chosen',{},'block0Path',{},'chemicalKey',{},'counterKey',{}));
    sel=selectedD1(ii);tr=sel.transportRow;
    pathMemo=containers.Map('KeyType','char','ValueType','any');
    if isempty(sel.counterLinks),return;end
    mc=find(metIDs==string(T.cytosolic(tr)),1);
    if isempty(mc) || comp(mc)~="c",return;end
    R.examined=1;
    inputRx=cons8{mc};
    if numel(inputRx)<2,return;end
    inputRx=inputRx(1:min(end,maxBranch));
    stop=false;
    for i1=reshape(inputRx,1,[])
        if stop,break;end
        nfCandidates=do{i1};
        nfCandidates=nfCandidates(arrayfun(@(v)d1SignalOK(v,mc,metIDs,comp,currency),nfCandidates));
        nfCandidates=nfCandidates(1:min(end,maxBranch));
        for nf=reshape(nfCandidates,1,[])
            if stop,break;end
            secondRx=cons8{nf};secondRx=secondRx(1:min(end,maxBranch));
            for i2=reshape(secondRx,1,[])
                if stop,break;end
                if di(i2)==di(i1),continue;end
                ndCandidates=do{i2};
                ndCandidates=ndCandidates(arrayfun(@(v)d1SignalOK(v,nf,metIDs,comp,currency) && v~=mc,ndCandidates));
                ndCandidates=ndCandidates(1:min(end,maxBranch));
                for nd=reshape(ndCandidates,1,[])
                    if stop,break;end
                    ndRx=cons8{nd};ndRx=ndRx(1:min(end,maxBranch));
                    if numel(ndRx)<2,continue;end
                    for i4=reshape(ndRx,1,[])
                        if stop,break;end
                        if ismember(di(i4),di([i1,i2])),continue;end
                        nf2Candidates=do{i4};
                        nf2Candidates=nf2Candidates(arrayfun(@(v)d1SignalOK(v,nd,metIDs,comp,currency) && ~ismember(v,[mc,nf]),nf2Candidates));
                        nf2Candidates=nf2Candidates(1:min(end,maxBranch));
                        for nf2=reshape(nf2Candidates,1,[])
                            if stop,break;end
                            fifthRx=cons8{nf2};fifthRx=fifthRx(1:min(end,maxBranch));
                            for i5=reshape(fifthRx,1,[])
                                if stop,break;end
                                if ismember(di(i5),di([i1,i2,i4])),continue;end
                                outCandidates=do{i5};
                                outCandidates=outCandidates(arrayfun(@(v)d1SignalOK(v,nf2,metIDs,comp,currency) && ~ismember(v,[mc,nf,nd]),outCandidates));
                                outCandidates=outCandidates(1:min(end,maxBranch));
                                for mout=reshape(outCandidates,1,[])
                                    if stop,break;end
                                    endRx=cons8{mout};endRx=endRx(1:min(end,maxBranch));
                                    if isempty(endRx),continue;end
                                    for irc=reshape(inputRx,1,[])
                                        if stop,break;end
                                        if di(irc)==di(i1),continue;end
                                        counterLink=find([sel.counterLinks.reaction]==di(irc),1);
                                        if isempty(counterLink),continue;end
                                        for inr=reshape(ndRx,1,[])
                                            if stop,break;end
                                            if di(inr)==di(i4),continue;end
                                            for inr2=reshape(endRx,1,[])
                                                if R.tests>=maxTests,R.capReached=1;stop=true;break;end
                                                R.tests=R.tests+1;
                                                chosen=[irc,i1,i2,inr,i4,i5,inr2];
                                                if numel(unique(di(chosen)))~=7
                                                    R.rejectDuplicateReaction=R.rejectDuplicateReaction+1;continue;
                                                end
                                                d=sel.d;d1ids=[d.rff,d.rfd,d.rfr,d.rsf,d.rsd,d.rsr,d.ri];
                                                if any(ismember(di(chosen),d1ids))
                                                    R.rejectD1Overlap=R.rejectD1Overlap+1;continue;
                                                end
                                                downstreamRx=[d1ids,sel.exchangePathRxns, ...
                                                    double(T.modelIndex(tr)),reshape(di(chosen),1,[])];
                                                if reactionHasExcludedEC(model,downstreamRx,excludedECs)
                                                    R.rejectExcludedEC=R.rejectExcludedEC+1;continue;
                                                end
                                                if requireUniqueEC && ~reactionsHaveUniqueECAssignment(model,downstreamRx)
                                                    R.rejectRepeatedEC=R.rejectRepeatedEC+1;continue;
                                                end
                                                b0=memoizedBlock0Path(model,block0Cache,d.min,downstreamRx, ...
                                                    [d.rff,d.rsf],sel.counterLinks(counterLink).path,pathMemo);
                                                if ~b0.found
                                                    R.rejectBlock0Isolation=R.rejectBlock0Isolation+1;
                                                    R.block0PathCapped=R.block0PathCapped+b0.expansionLimitReached+b0.stateLimitReached;
                                                    continue;
                                                end
                                                block0Reactions=b0.rxns;
                                                allRx=[block0Reactions,d1ids,sel.exchangePathRxns, ...
                                                    double(T.modelIndex(tr)),reshape(di(chosen),1,[])];
                                                requiredECReactions=[block0Reactions,d1ids,sel.exchangePathRxns, ...
                                                    reshape(di(chosen),1,[])];
                                                if ~reactionsHaveRequiredEC(model,requiredECReactions)
                                                    R.rejectMissingEC=R.rejectMissingEC+1;continue;
                                                end
                                                if reactionHasExcludedEC(model,allRx,excludedECs)
                                                    R.rejectExcludedEC=R.rejectExcludedEC+1;continue;
                                                end
                                                if requireUniqueEC && ~reactionsHaveUniqueECAssignment(model,allRx)
                                                    R.rejectRepeatedEC=R.rejectRepeatedEC+1;continue;
                                                end
                                                chemicalKey=sprintf('%d,',[d.min,d.mff,d.mfd,d.msf,sel.target, ...
                                                    mc,nf,nd,nf2,mout,d1ids,sel.exchangePathRxns, ...
                                                    double(T.modelIndex(tr)),reshape(di(chosen),1,[])]);
                                                counterKey=counterChemicalKey(sel,double(T.modelIndex(tr)),di(irc),mc);
                                                c=struct('mc',mc,'nf',nf,'nd',nd,'nf2',nf2,'mout',mout, ...
                                                    'chosen',chosen,'block0Path',b0, ...
                                                    'chemicalKey',chemicalKey,'counterKey',counterKey);
                                                R.candidates(end+1)=c; %#ok<AGROW>
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

function tf=reactionHasExcludedEC(model,reactionIDs,excludedECs)
% True if any reaction annotation matches an excluded EC number.
    tf=false;
    excludedECs=string(excludedECs(:));
    if isempty(excludedECs),return;end
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
