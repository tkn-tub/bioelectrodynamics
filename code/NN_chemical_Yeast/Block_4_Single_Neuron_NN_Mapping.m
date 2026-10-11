%% Block 4: map the single-neuron NN to yeast-GEM
%
% The counter output is the NN input (t_peak in the diagram):
%
%   m_c --r_omega--> m_out (= m_dhat)
%                       ^
%                       |
%             buffered substrates --r_b--
%
% Both reactions produce the SAME output metabolite. Their kinetic rates
% implement the affine neuron
%
%   d_hat = omega*t_peak + b.
%
% omega and b are kinetic parameters, not metabolites. This script finds
% candidate reactions r_omega and r_b; their enzyme concentration, kcat and
% Km values are assigned later to realize the desired numerical omega and b.
%
% This version produces the same type of full connection table as Block 3.
% It saves:
%
%   1) Motifs_NN_Block_4: full Block 1 + Block 2 + counter + NN rows;
%   2) a compact list of NN candidates indexed by the counter output m_c;
%   3) a link table mapping every Block 3 motif to its NN input group.
%
% The full table performs only the connection through the shared m_c. No
% cross-block reaction, metabolite, gene or crosstalk filters are applied here;
% those operations are intentionally deferred to a later filtering script.
%
% No separate m_out-removal reaction is included because the supplied NN
% diagram contains exactly these two reaction roles. Output clearance can be
% added later if a steady-state rather than an accumulated readout is needed.
%
% Expected master-file sequence:
%   load Block_3_counter.mat
%   Block_4_Single_Neuron_NN_Mapping
%
% Required workspace variables:
%   model, model.S (or S), isBadMet
% Block 3 must provide either:
%   results_counter             struct array, or
%   Motifs_Counter_Block_3      table

%% Validate and prepare the workspace
if ~exist('model','var') || ~isfield(model,'S')
    error('Block4:MissingModel', ...
          'The COBRA model must be loaded before running Block 4.');
end

if ~exist('S','var'), S = model.S; end
[nMets,nRxns] = size(S);

if ~exist('isBadMet','var') || numel(isBadMet) ~= nMets
    error('Block4:MissingCurrencyMask', ...
        ['isBadMet is missing or has the wrong length. Run the currency-', ...
         'metabolite setup in the master file before Block 4.']);
end
badMetIdx = find(isBadMet);

if ~exist('hasGPR','var') || numel(hasGPR) ~= nRxns
    hasGPR = false(nRxns,1);
    for r = 1:nRxns
        rule = strtrim(string(model.grRules{r}));
        hasGPR(r) = strlength(rule) > 0 && ~strcmpi(rule,"s0001");
    end
end

if ~exist('substrateSet','var') || numel(substrateSet) ~= nRxns || ...
        ~exist('productSet','var') || numel(productSet) ~= nRxns
    substrateSet = cell(nRxns,1);
    productSet = cell(nRxns,1);
    for r = 1:nRxns
        substrateSet{r} = find(S(:,r) < 0);
        productSet{r} = find(S(:,r) > 0);
    end
end

if ~exist('nProducers','var') || numel(nProducers) ~= nMets
    nProducers = full(sum(S > 0,2));
end

if exist('results_counter','var') && isstruct(results_counter)
    block3Results = results_counter(:);
elseif exist('Motifs_Counter_Block_3','var') && ...
        istable(Motifs_Counter_Block_3)
    block3Results = table2struct(Motifs_Counter_Block_3);
    block3Results = block3Results(:);
else
    error('Block4:MissingBlock3', ...
        ['Load Block_3_counter.mat before running Block 4. Expected ', ...
         'results_counter or Motifs_Counter_Block_3.']);
end

if ~isfield(block3Results,'m_c')
    error('Block4:IncompleteBlock3', ...
          'Block 3 results do not contain the counter output field m_c.');
end

nBlock3 = numel(block3Results);
if nBlock3 == 0
    error('Block4:EmptyBlock3', ...
          'Block 3 contains no counter motifs to connect to the NN.');
end

%% Topology-search controls
if ~exist('REQUIRE_GPR_NN','var')
    if exist('REQUIRE_GPR','var')
        REQUIRE_GPR_NN = REQUIRE_GPR;
    else
        REQUIRE_GPR_NN = true;
    end
end

% The bias must be input-independent: it produces m_out but does not consume
% m_c. Other cross-block constraints are intentionally postponed.
if ~exist('STRICT_BIAS_REACTION','var')
    STRICT_BIAS_REACTION = true;
end

% Optional local filters. They default to broad discovery because the user
% requested that filtering be performed later.
if ~exist('REQUIRE_BIAS_NONCURRENCY_SUBSTRATE','var')
    REQUIRE_BIAS_NONCURRENCY_SUBSTRATE = false;
end
if ~exist('REQUIRE_NN_GENE_DISJOINT','var')
    REQUIRE_NN_GENE_DISJOINT = false;
end

if ~exist('MAX_PRODUCERS_MOUT','var')
    if exist('MAX_PRODUCERS_MNC','var')
        MAX_PRODUCERS_MOUT = MAX_PRODUCERS_MNC;
    else
        MAX_PRODUCERS_MOUT = inf;
    end
end

if ~exist('MAX_HITS_NN','var')
    if exist('MAX_HITS','var')
        MAX_HITS_NN = MAX_HITS;
    else
        MAX_HITS_NN = inf;
    end
end

% Maximum number of rows in the flattened end-to-end table. A finite default
% prevents the Cartesian expansion from exhausting memory. Set a larger value
% when required. An unbounded request is protected by a safety cap unless
% ALLOW_UNBOUNDED_FULL_NN is explicitly enabled.
if ~exist('MAX_FULL_MOTIFS_NN','var')
    MAX_FULL_MOTIFS_NN = 1e5;
end
if ~exist('FULL_TABLE_SAFETY_CAP','var')
    FULL_TABLE_SAFETY_CAP = 1e5;
end
if ~exist('ALLOW_UNBOUNDED_FULL_NN','var')
    ALLOW_UNBOUNDED_FULL_NN = false;
end

if ~exist('SHOW_PROGRESS','var'), SHOW_PROGRESS = true; end
if ~exist('PROGRESS_STEP_PERCENT','var'), PROGRESS_STEP_PERCENT = 1; end
PROGRESS_STEP_PERCENT = max(1,min(100,round(PROGRESS_STEP_PERCENT)));

% With a finite hit cap, serial is the default because it stops exactly at
% MAX_HITS_NN. Set USE_PARALLEL_NN=true on the HPC cluster to scan distinct
% m_c inputs concurrently. The number of distinct inputs is often small, so
% parallel speedup depends on how many different m_c metabolites are present.
if ~exist('USE_PARALLEL_NN','var')
    USE_PARALLEL_NN = ~isfinite(MAX_HITS_NN);
end
if ~exist('AUTO_START_PARPOOL','var'), AUTO_START_PARPOOL = true; end
if ~exist('PARPOOL_WORKERS','var'), PARPOOL_WORKERS = []; end

reactionGenes = cell(nRxns,1);
if REQUIRE_NN_GENE_DISJOINT
    for r = 1:nRxns
        gpr = string(model.grRules{r});
        tokens = string(regexp(char(gpr), ...
            '[A-Za-z0-9_:\-\.]+','match'));
        tokens = tokens(~ismember(lower(tokens),["and","or"]));
        reactionGenes{r} = unique(tokens);
    end
end

%% Build the compact Block 3 -> NN input connection
allMC = [block3Results.m_c].';
[uniqueMC,~,inputOfBlock3] = unique(allMC,'stable');
nInputs = numel(uniqueMC);

fprintf('\nBlock 4 single-neuron input preparation:\n');
fprintf('  Differentiator + counter motifs: %d\n',nBlock3);
fprintf('  Distinct m_c inputs             : %d\n',nInputs);
fprintf('  Cross-block filters             : deferred\n');

% Every Block 3 row points to exactly one NN input group. This compact mapping
% is retained alongside the requested flattened table for traceability and for
% later filtering without needing to materialize additional combinations.
Counter_NN_Link_Block_4 = table( ...
    (1:nBlock3).',inputOfBlock3,allMC, ...
    'VariableNames',{'block3_idx','nn_input_idx','m_c'});

NN_Input_Lookup_Block_4 = table( ...
    (1:nInputs).',uniqueMC, ...
    'VariableNames',{'nn_input_idx','m_c'});

candidateTemplate = struct( ...
    'nn_candidate_idx',[], ...
    'nn_input_idx',[], ...
    'm_c',[], ...
    'm_out',[], ...
    'm_dhat',[], ...
    'r_omega',[], ...
    'r_b',[], ...
    'r_Nomega',[], ...
    'r_Nb',[]);

results_nn_candidates = repmat(candidateTemplate,0,1);
nnResultsByInput = cell(nInputs,1);

%% Select serial or parallel topology search
runInParallel = false;
poolObject = [];

if USE_PARALLEL_NN
    try
        hasParallelToolbox = ~isempty(ver('parallel')) && ...
                             license('test','Distrib_Computing_Toolbox');
        if hasParallelToolbox
            poolObject = gcp('nocreate');
            if isempty(poolObject) && AUTO_START_PARPOOL
                if isempty(PARPOOL_WORKERS)
                    poolObject = parpool;
                else
                    poolObject = parpool(PARPOOL_WORKERS);
                end
            end
            runInParallel = ~isempty(poolObject);
        end
    catch parallelError
        warning('Block4:ParallelFallback', ...
            'Could not start/reuse a parallel pool (%s). Using serial search.', ...
            parallelError.message);
    end
end

if USE_PARALLEL_NN && ~runInParallel
    warning('Block4:NoParallelPool', ...
        ['Parallel NN search was requested, but no pool is available. ', ...
         'The search will run serially.']);
end

searchTimer = tic;

if runInParallel
    fprintf('\nRunning compact NN search on %d MATLAB workers.\n', ...
            poolObject.NumWorkers);

    if isfinite(MAX_HITS_NN)
        warning('Block4:ParallelMaxHits', ...
            ['In parallel mode MAX_HITS_NN is applied after all m_c ', ...
             'workers finish.']);
        maxParallelInputSolutions = MAX_HITS_NN;
    else
        maxParallelInputSolutions = inf;
    end

    if SHOW_PROGRESS
        progressQueue = parallel.pool.DataQueue;
        nnParallelProgress('init',0,nInputs,PROGRESS_STEP_PERCENT);
        progressListener = afterEach(progressQueue, ...
            @(localHits) nnParallelProgress( ...
                'tick',localHits,nInputs,PROGRESS_STEP_PERCENT));
    else
        progressQueue = [];
        progressListener = [];
    end

    parfor inputIdx = 1:nInputs
        localResults = searchSingleNeuronForInput( ...
            uniqueMC(inputIdx),S,substrateSet,productSet,badMetIdx, ...
            hasGPR,nProducers,reactionGenes,REQUIRE_GPR_NN, ...
            STRICT_BIAS_REACTION,REQUIRE_BIAS_NONCURRENCY_SUBSTRATE, ...
            REQUIRE_NN_GENE_DISJOINT,MAX_PRODUCERS_MOUT, ...
            maxParallelInputSolutions);
        nnResultsByInput{inputIdx} = localResults;
        if SHOW_PROGRESS
            send(progressQueue,numel(localResults));
        end
    end

    if SHOW_PROGRESS
        nnParallelProgress('finish',0,nInputs,PROGRESS_STEP_PERCENT);
        delete(progressListener);
    end

else
    lastProgressPercent = 0;
    if SHOW_PROGRESS
        fprintf(['Compact NN search:   0%% | m_c inputs scanned: 0/%d', ...
                 ' | NN pairs: 0 | elapsed: 0.0 min'],nInputs);
    end

    serialRawHits = 0;
    for inputIdx = 1:nInputs
        progressPercent = floor(100*(inputIdx-1)/max(nInputs,1));
        if SHOW_PROGRESS && ...
                progressPercent >= lastProgressPercent + PROGRESS_STEP_PERCENT
            fprintf(['\rCompact NN search: %3d%% | m_c inputs scanned: %d/%d', ...
                     ' | NN pairs: %d | elapsed: %.1f min'], ...
                    progressPercent,inputIdx-1,nInputs,serialRawHits, ...
                    toc(searchTimer)/60);
            lastProgressPercent = progressPercent;
            drawnow;
        end

        if isfinite(MAX_HITS_NN)
            maxInputSolutions = MAX_HITS_NN - serialRawHits;
        else
            maxInputSolutions = inf;
        end
        if maxInputSolutions <= 0, break; end

        localResults = searchSingleNeuronForInput( ...
            uniqueMC(inputIdx),S,substrateSet,productSet,badMetIdx, ...
            hasGPR,nProducers,reactionGenes,REQUIRE_GPR_NN, ...
            STRICT_BIAS_REACTION,REQUIRE_BIAS_NONCURRENCY_SUBSTRATE, ...
            REQUIRE_NN_GENE_DISJOINT,MAX_PRODUCERS_MOUT, ...
            maxInputSolutions);
        nnResultsByInput{inputIdx} = localResults;
        serialRawHits = serialRawHits + numel(localResults);
    end

    if SHOW_PROGRESS
        fprintf(['\rCompact NN search: 100%% | m_c inputs scanned: %d/%d', ...
                 ' | NN pairs: %d | elapsed: %.1f min\n'], ...
                min(inputIdx,nInputs),nInputs,serialRawHits, ...
                toc(searchTimer)/60);
    end
end

%% Combine compact candidates in deterministic m_c order
nHitsNN = 0;
nInputsWithSolutions = 0;
stopCombine = false;

for inputIdx = 1:nInputs
    localResults = nnResultsByInput{inputIdx};
    if isempty(localResults), continue; end
    nInputsWithSolutions = nInputsWithSolutions + 1;

    for qN = 1:numel(localResults)
        nHitsNN = nHitsNN + 1;
        local = localResults(qN);

        rec = candidateTemplate;
        rec.nn_candidate_idx = nHitsNN;
        rec.nn_input_idx = inputIdx;
        rec.m_c = uniqueMC(inputIdx);
        rec.m_out = local.m_out;
        rec.m_dhat = local.m_out;
        rec.r_omega = local.r_omega;
        rec.r_b = local.r_b;
        rec.r_Nomega = local.r_omega;
        rec.r_Nb = local.r_b;
        results_nn_candidates(nHitsNN,1) = rec;

        if nHitsNN >= MAX_HITS_NN
            stopCombine = true;
            break;
        end
    end
    if stopCombine, break; end
end

elapsedMinutes = toc(searchTimer)/60;

%% Materialize the full Block 3 + NN table without applying filters
Motifs_NN_Candidates_Block_4 = struct2table(results_nn_candidates);

if nHitsNN > 0
    candidateInputIdx = [results_nn_candidates.nn_input_idx].';
    candidatesPerInput = accumarray(candidateInputIdx,1,[nInputs,1]);
else
    candidateInputIdx = zeros(0,1);
    candidatesPerInput = zeros(nInputs,1);
end
counterMotifsPerInput = accumarray(inputOfBlock3,1,[nInputs,1]);
totalFullCombinations = sum( ...
    double(counterMotifsPerInput).*double(candidatesPerInput));

effectiveFullLimit = MAX_FULL_MOTIFS_NN;
if ~isfinite(effectiveFullLimit) && ~ALLOW_UNBOUNDED_FULL_NN && ...
        totalFullCombinations > FULL_TABLE_SAFETY_CAP
    warning('Block4:FullTableSafetyCap', ...
        ['The unfiltered full table contains %.0f possible rows. ', ...
         'This run is limited to FULL_TABLE_SAFETY_CAP=%.0f rows to ', ...
         'protect memory. Set a finite MAX_FULL_MOTIFS_NN, or set ', ...
         'ALLOW_UNBOUNDED_FULL_NN=true to materialize every row.'], ...
        totalFullCombinations,FULL_TABLE_SAFETY_CAP);
    effectiveFullLimit = FULL_TABLE_SAFETY_CAP;
end

if isfinite(effectiveFullLimit)
    nFullTarget = min(totalFullCombinations, ...
                      max(0,floor(effectiveFullLimit)));
else
    nFullTarget = totalFullCombinations;
end

jointTemplate = block3Results(1);
jointTemplate.block3_idx = [];
jointTemplate.nn_input_idx = [];
jointTemplate.nn_candidate_idx = [];
jointTemplate.m_out = [];
jointTemplate.m_dhat = [];
jointTemplate.r_omega = [];
jointTemplate.r_b = [];
jointTemplate.r_Nomega = [];
jointTemplate.r_Nb = [];

results_nn = repmat(jointTemplate,nFullTarget,1);
candidateIndicesByInput = cell(nInputs,1);
for inputIdx = 1:nInputs
    candidateIndicesByInput{inputIdx} = find(candidateInputIdx == inputIdx);
end

fprintf('\nFull Block 4 table construction:\n');
fprintf('  Possible unfiltered combinations: %.0f\n',totalFullCombinations);
fprintf('  Rows to materialize             : %.0f\n',nFullTarget);

nFull = 0;
lastFullProgress = 0;
fullTableTimer = tic;
if SHOW_PROGRESS
    fprintf(['Full Block 4 table:   0%% | rows created: 0/%.0f', ...
             ' | elapsed: 0.0 min'],nFullTarget);
end

for block3Idx = 1:nBlock3
    inputIdx = inputOfBlock3(block3Idx);
    candidateIndices = candidateIndicesByInput{inputIdx};
    if isempty(candidateIndices), continue; end

    upstream = block3Results(block3Idx);
    for qN = 1:numel(candidateIndices)
        if nFull >= nFullTarget, break; end

        candidateIdx = candidateIndices(qN);
        candidate = results_nn_candidates(candidateIdx);
        nFull = nFull + 1;

        joint = upstream;
        joint.block3_idx = block3Idx;
        joint.nn_input_idx = inputIdx;
        joint.nn_candidate_idx = candidate.nn_candidate_idx;
        joint.m_out = candidate.m_out;
        joint.m_dhat = candidate.m_dhat;
        joint.r_omega = candidate.r_omega;
        joint.r_b = candidate.r_b;
        joint.r_Nomega = candidate.r_Nomega;
        joint.r_Nb = candidate.r_Nb;
        results_nn(nFull,1) = joint;

        fullProgress = floor(100*nFull/max(nFullTarget,1));
        if SHOW_PROGRESS && ...
                (fullProgress >= lastFullProgress + PROGRESS_STEP_PERCENT || ...
                 nFull == nFullTarget)
            fprintf(['\rFull Block 4 table: %3d%% | rows created: %.0f/%.0f', ...
                     ' | elapsed: %.1f min'], ...
                    fullProgress,nFull,nFullTarget,toc(fullTableTimer)/60);
            lastFullProgress = fullProgress;
            drawnow;
        end
    end
    if nFull >= nFullTarget, break; end
end

if nFull < nFullTarget
    results_nn = results_nn(1:nFull);
end
if SHOW_PROGRESS
    fprintf('\n');
end

Motifs_NN_Block_4 = struct2table(results_nn);

%% Save full table, compact candidates and connection indices
scriptFolder = fileparts(mfilename('fullpath'));
outputFile = fullfile(scriptFolder,'Block_4_single_neuron_NN.mat');
save(outputFile, ...
    'Motifs_NN_Block_4','results_nn', ...
    'Motifs_NN_Candidates_Block_4','results_nn_candidates', ...
    'Counter_NN_Link_Block_4','NN_Input_Lookup_Block_4');

fprintf('\nSingle-neuron search diagnostics:\n');
fprintf('  Block 3 motifs linked       : %d\n',nBlock3);
fprintf('  Distinct m_c input groups   : %d\n',nInputs);
fprintf('  Input groups with solutions : %d\n',nInputsWithSolutions);
fprintf('  Compact NN candidates       : %d\n',nHitsNN);
fprintf('  Search wall time            : %.1f min\n',elapsedMinutes);
fprintf('  Possible full combinations  : %.0f\n',totalFullCombinations);
fprintf('  Full rows materialized      : %d\n',nFull);
fprintf('  Full table columns          : %d\n',width(Motifs_NN_Block_4));
fprintf('Saved Block 4 results to:\n%s\n',outputFile);

fprintf('\nTo retrieve the materialized full rows for Block 3 motif k:\n');
fprintf('  rows = Motifs_NN_Block_4.block3_idx == k;\n');
fprintf('  fullSolutions = Motifs_NN_Block_4(rows,:);\n');

if nFull > 0
    fprintf('\n=== First full differentiator + counter + NN rows ===\n');
    for k = 1:min(5,nFull)
        rec = results_nn(k);
        fprintf('\n--- Full candidate %d (Block 3 motif %d) ---\n', ...
                k,rec.block3_idx);
        fprintf('  Input m_c     = %s | %s\n', ...
                model.mets{rec.m_c},model.metNames{rec.m_c});
        fprintf('  Output m_dhat = %s | %s\n', ...
                model.mets{rec.m_out},model.metNames{rec.m_out});

        labels = {'r_omega','r_b'};
        rxns = [rec.r_omega,rec.r_b];
        for q = 1:numel(rxns)
            fprintf('  %-7s = %s | %s\n',labels{q}, ...
                    model.rxns{rxns(q)},model.rxnNames{rxns(q)});
            fprintf('            GPR: %s\n',model.grRules{rxns(q)});
        end
    end
end

%% Local single-neuron search
function nnResults = searchSingleNeuronForInput( ...
        m_c,S,substrateSet,productSet,badMetIdx,hasGPR,nProducers, ...
        reactionGenes,REQUIRE_GPR_NN,STRICT_BIAS_REACTION, ...
        REQUIRE_BIAS_NONCURRENCY_SUBSTRATE,REQUIRE_NN_GENE_DISJOINT, ...
        MAX_PRODUCERS_MOUT,maxHits)

nnTemplate = struct('m_out',[],'r_omega',[],'r_b',[]);
nnResults = repmat(nnTemplate,0,1);
nLocal = 0;

% Multiplier reaction: m_c is a substrate and m_out is a product.
rOmegaCandidates = find(S(m_c,:) < 0);
if REQUIRE_GPR_NN
    rOmegaCandidates = rOmegaCandidates(hasGPR(rOmegaCandidates));
end

for qW = 1:numel(rOmegaCandidates)
    r_omega = rOmegaCandidates(qW);
    outputCandidates = setdiff(productSet{r_omega},badMetIdx);
    outputCandidates = setdiff(outputCandidates,m_c);

    for qO = 1:numel(outputCandidates)
        m_out = outputCandidates(qO);

        % Two different producers are necessary: r_omega and r_b.
        if nProducers(m_out) < 2 || ...
                nProducers(m_out) > MAX_PRODUCERS_MOUT
            continue;
        end

        rBCandidates = find(S(m_out,:) > 0);
        if REQUIRE_GPR_NN
            rBCandidates = rBCandidates(hasGPR(rBCandidates));
        end
        rBCandidates = setdiff(rBCandidates,r_omega);
        if STRICT_BIAS_REACTION
            rBCandidates = rBCandidates(S(m_c,rBCandidates) >= 0);
        end

        for qB = 1:numel(rBCandidates)
            r_b = rBCandidates(qB);

            if REQUIRE_BIAS_NONCURRENCY_SUBSTRATE && ...
                    isempty(setdiff(substrateSet{r_b},badMetIdx))
                continue;
            end
            if REQUIRE_NN_GENE_DISJOINT && ...
                    reactionsShareGenesNN( ...
                        [r_omega,r_b],reactionGenes)
                continue;
            end

            nLocal = nLocal + 1;
            rec = nnTemplate;
            rec.m_out = m_out;
            rec.r_omega = r_omega;
            rec.r_b = r_b;
            nnResults(nLocal,1) = rec;

            if nLocal >= maxHits, break; end
        end
        if nLocal >= maxHits, break; end
    end
    if nLocal >= maxHits, break; end
end

if ~isempty(nnResults)
    key = [[nnResults.m_out].',[nnResults.r_omega].',[nnResults.r_b].'];
    [~,uniqueIdx] = unique(key,'rows','stable');
    nnResults = nnResults(uniqueIdx);
end
end

function tf = reactionsShareGenesNN(rxns,reactionGenes)
% True when any distinct reaction pair shares at least one GPR token.
rxns = unique(rxns(:).','stable');
tf = false;
for i = 1:numel(rxns)
    for j = i+1:numel(rxns)
        if ~isempty(intersect(reactionGenes{rxns(i)}, ...
                              reactionGenes{rxns(j)}))
            tf = true;
            return;
        end
    end
end
end

function nnParallelProgress(command,localHits,nInputs,stepPercent)
% Client-side progress display for DataQueue messages from parfor workers.
persistent completedInputs totalHits lastPercent progressTimer

switch command
    case 'init'
        completedInputs = 0;
        totalHits = 0;
        lastPercent = 0;
        progressTimer = tic;
        fprintf(['Compact NN parallel search:   0%%', ...
                 ' | m_c inputs finished: 0/%d', ...
                 ' | NN pairs: 0 | elapsed: 0.0 min'],nInputs);

    case 'tick'
        completedInputs = completedInputs + 1;
        totalHits = totalHits + localHits;
        percent = floor(100*completedInputs/max(nInputs,1));
        if percent >= lastPercent + stepPercent || completedInputs == nInputs
            fprintf(['\rCompact NN parallel search: %3d%%', ...
                     ' | m_c inputs finished: %d/%d', ...
                     ' | NN pairs: %d | elapsed: %.1f min'], ...
                    percent,completedInputs,nInputs,totalHits, ...
                    toc(progressTimer)/60);
            lastPercent = percent;
            drawnow;
        end

    case 'finish'
        fprintf('\n');
end
end
