%% Block 3: map the three-reaction chemical counter to yeast-GEM
%
% Counter topology from the block diagram:
%
%   m_d --r_C--> m_1 (= m_c, counter output)
%
%   m_1 + m_s --r_CI1--> ignored products
%   m_d + m_s --r_CI2--> ignored products
%
% m_d is the bottom-branch output of the second differentiator. m_s is the
% shared counter species. The names r_CI1/r_CI2 avoid a collision with r_I2,
% which is already the cancellation reaction of the second differentiator.
%
% Compatibility aliases used by the earlier counter/NN code are exported:
%   m_c  = m_1,       m_Cs = m_s
%   r_Cp = r_C,       r_Cr = r_CI1,       r_Cg = r_CI2
% The diagram labels are also exported without colliding with the two
% differentiators: r_I1_counter = r_CI1 and r_I2_counter = r_CI2.
%
% Expected master-file sequence:
%   load Block_2_second_differentiator.mat
%   Block_3_Counter_Mapping
%
% Required workspace variables:
%   model, model.S (or S), isBadMet
% Block 2 must provide either:
%   results_diff2             struct array, or
%   Motifs_Diff_Block_2       table

%% Validate and prepare the workspace
if ~exist('model','var') || ~isfield(model,'S')
    error('Block3:MissingModel', ...
          'The COBRA model must be loaded before running Block 3.');
end

if ~exist('S','var'), S = model.S; end
[nMets,nRxns] = size(S);

if ~exist('isBadMet','var') || numel(isBadMet) ~= nMets
    error('Block3:MissingCurrencyMask', ...
        ['isBadMet is missing or has the wrong length. Run the currency-', ...
         'metabolite setup in the master file before Block 3.']);
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

if exist('results_diff2','var') && isstruct(results_diff2)
    block2Results = results_diff2(:);
elseif exist('Motifs_Diff_Block_2','var') && istable(Motifs_Diff_Block_2)
    block2Results = table2struct(Motifs_Diff_Block_2);
    block2Results = block2Results(:);
else
    error('Block3:MissingBlock2', ...
        ['Load Block_2_second_differentiator.mat before running Block 3. ', ...
         'Expected results_diff2 or Motifs_Diff_Block_2.']);
end

requiredBlock2Fields = { ...
    'm_in','m_S1','m_F1','m_Sg1','m_Fg1', ...
    'm_S2','m_F2','m_Sg2','m_Fg2','m_d', ...
    'r_S1','r_Sg1','r_Ss1','r_F1','r_Fg1','r_Fs1','r_I1', ...
    'r_S2','r_Sg2','r_Ss2','r_F2','r_Fg2','r_Fs2','r_I2'};
missingFields = setdiff(requiredBlock2Fields,fieldnames(block2Results));
if ~isempty(missingFields)
    error('Block3:IncompleteBlock2', ...
        'Block 2 results are missing fields: %s',strjoin(missingFields,', '));
end

nBlock2 = numel(block2Results);
if nBlock2 == 0
    error('Block3:EmptyBlock2', ...
          'Block 2 contains no motifs to extend with the counter.');
end

%% Search controls
if ~exist('REQUIRE_GPR','var'), REQUIRE_GPR = true; end
if ~exist('STRICT_COUNTER_EDGES','var')
    % r_CI1 may consume m_1+m_s, but not m_d; r_CI2 may consume m_d+m_s,
    % but not m_1. This prevents unintended shortcuts inside the counter.
    STRICT_COUNTER_EDGES = true;
end
if ~exist('ALLOW_MD_AS_MS','var'), ALLOW_MD_AS_MS = false; end
if ~exist('ALLOW_CURRENCY_MS','var')
    % m_s is a buffered/shared counter species.  When enabled, ATP, NAD(P),
    % AMP and other currency metabolites may fill only this role.  The
    % counter input m_d and output m_1 remain non-currency metabolites.
    ALLOW_CURRENCY_MS = false;
end

if ~exist('MAX_PRODUCERS_M1','var')
    if exist('MAX_PRODUCERS_MC','var')
        MAX_PRODUCERS_M1 = MAX_PRODUCERS_MC;
    else
        MAX_PRODUCERS_M1 = inf;
    end
end
if ~exist('MAX_PRODUCERS_MS','var')
    if exist('MAX_PRODUCERS_MCS','var')
        MAX_PRODUCERS_MS = MAX_PRODUCERS_MCS;
    else
        MAX_PRODUCERS_MS = inf;
    end
end

% Broad discovery is the default; these can be enabled for a stricter join.
if ~exist('REQUIRE_COUNTER_DISJOINT_RXNS','var')
    REQUIRE_COUNTER_DISJOINT_RXNS = false;
end
if ~exist('REQUIRE_COUNTER_DISJOINT_METS','var')
    REQUIRE_COUNTER_DISJOINT_METS = false;
end

% The attached earlier counter mapper rejected counter reactions that share
% a gene. Keep that behavior by default. The cross-block form is optional,
% because requiring disjoint genes across all 17 reactions is much stricter
% and can remove every otherwise valid metabolic realization.
if ~exist('REQUIRE_COUNTER_GENE_DISJOINT','var')
    REQUIRE_COUNTER_GENE_DISJOINT = true;
end
if ~exist('REQUIRE_COUNTER_UPSTREAM_GENE_DISJOINT','var')
    REQUIRE_COUNTER_UPSTREAM_GENE_DISJOINT = false;
end

reactionGenes = cell(nRxns,1);
if REQUIRE_COUNTER_GENE_DISJOINT || ...
        REQUIRE_COUNTER_UPSTREAM_GENE_DISJOINT
    for r = 1:nRxns
        gpr = string(model.grRules{r});
        tokens = string(regexp(char(gpr), ...
            '[A-Za-z0-9_:\-\.]+','match'));
        tokens = tokens(~ismember(lower(tokens),["and","or"]));
        reactionGenes{r} = unique(tokens);
    end
end

if ~exist('MAX_HITS_COUNTER','var')
    if exist('MAX_HITS','var')
        MAX_HITS_COUNTER = MAX_HITS;
    else
        MAX_HITS_COUNTER = inf;
    end
end

if ~exist('SHOW_PROGRESS','var'), SHOW_PROGRESS = true; end
if ~exist('PROGRESS_STEP_PERCENT','var'), PROGRESS_STEP_PERCENT = 1; end
PROGRESS_STEP_PERCENT = max(1,min(100,round(PROGRESS_STEP_PERCENT)));

% Parallelize independent m_d searches when requested. With a finite hit
% cap, serial is the default because it can stop the nested search exactly
% when MAX_HITS_COUNTER is reached. Set USE_PARALLEL_COUNTER=true on the HPC
% cluster to scan all distinct m_d inputs concurrently; the final saved
% output is still capped exactly, but workers already in flight cannot be
% stopped by that cap.
if ~exist('USE_PARALLEL_COUNTER','var')
    USE_PARALLEL_COUNTER = ~isfinite(MAX_HITS_COUNTER);
end
if ~exist('AUTO_START_PARPOOL','var'), AUTO_START_PARPOOL = true; end
if ~exist('PARPOOL_WORKERS','var'), PARPOOL_WORKERS = []; end

%% Search each distinct second-differentiator output only once
allMD = [block2Results.m_d].';
[uniqueMD,~,inputOfBlock2] = unique(allMD,'stable');
nInputs = numel(uniqueMD);

fprintf('\nBlock 3 counter input preparation:\n');
fprintf('  Two-differentiator motifs : %d\n',nBlock2);
fprintf('  Distinct m_d inputs       : %d\n',nInputs);

% Extend the complete Block 2 record so every final row remains traceable.
jointTemplate = block2Results(1);
jointTemplate.block2_idx = [];
jointTemplate.counter_input_idx = [];
jointTemplate.m_1 = [];
jointTemplate.m_s = [];
jointTemplate.m_c = [];
jointTemplate.m_Cs = [];
jointTemplate.r_C = [];
jointTemplate.r_CI1 = [];
jointTemplate.r_CI2 = [];
jointTemplate.r_I1_counter = [];
jointTemplate.r_I2_counter = [];
jointTemplate.r_Cp = [];
jointTemplate.r_Cr = [];
jointTemplate.r_Cg = [];

results_counter = repmat(jointTemplate,0,1);
nHitsCounter = 0;
stopSearch = false;
nInputsWithSolutions = 0;
nCounterMotifs = 0;

%% Optionally search the distinct m_d inputs in parallel
runInParallel = false;
poolObject = [];
counterResultsByInput = cell(nInputs,1);
counterDiagnosticsByInput = cell(nInputs,1);

if USE_PARALLEL_COUNTER
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
        warning('Block3:ParallelFallback', ...
            'Could not start/reuse a parallel pool (%s). Using serial search.', ...
            parallelError.message);
    end
end

if USE_PARALLEL_COUNTER && ~runInParallel
    warning('Block3:NoParallelPool', ...
        ['Parallel counter search was requested, but no pool is available. ', ...
         'The search will run serially.']);
end

if runInParallel
    fprintf('\nRunning counter search on %d MATLAB workers.\n', ...
            poolObject.NumWorkers);
    if isfinite(MAX_HITS_COUNTER)
        warning('Block3:ParallelMaxHits', ...
            ['In parallel mode MAX_HITS_COUNTER is an exact output cap, ', ...
             'not a worker early-stop condition.']);
        maxParallelInputSolutions = MAX_HITS_COUNTER;
    else
        maxParallelInputSolutions = inf;
    end

    if SHOW_PROGRESS
        progressQueue = parallel.pool.DataQueue;
        counterParallelProgress('init',0,nInputs,PROGRESS_STEP_PERCENT);
        progressListener = afterEach(progressQueue, ...
            @(localHits) counterParallelProgress( ...
                'tick',localHits,nInputs,PROGRESS_STEP_PERCENT));
    else
        progressQueue = [];
        progressListener = [];
    end

    parfor parallelInputIdx = 1:nInputs
        localMD = uniqueMD(parallelInputIdx);
        [localResults,localDiagnostics] = searchCounterForInput( ...
            localMD,S,substrateSet,productSet,badMetIdx,hasGPR, ...
            nProducers,reactionGenes, ...
            REQUIRE_GPR,STRICT_COUNTER_EDGES,ALLOW_MD_AS_MS, ...
            ALLOW_CURRENCY_MS, ...
            REQUIRE_COUNTER_GENE_DISJOINT, ...
            MAX_PRODUCERS_M1,MAX_PRODUCERS_MS, ...
            maxParallelInputSolutions);
        counterResultsByInput{parallelInputIdx} = localResults;
        counterDiagnosticsByInput{parallelInputIdx} = localDiagnostics;
        if SHOW_PROGRESS
            send(progressQueue,numel(localResults));
        end
    end

    if SHOW_PROGRESS
        counterParallelProgress('finish',0,nInputs,PROGRESS_STEP_PERCENT);
        delete(progressListener);
    end
end

searchTimer = tic;
lastProgressPercent = 0;
if runInParallel
    progressLabel = 'Counter join';
else
    progressLabel = 'Counter search';
end
if SHOW_PROGRESS
    fprintf(['%s:   0%% | m_d inputs scanned: 0/%d', ...
             ' | joint candidates: 0 | elapsed: 0.0 min'], ...
             progressLabel,nInputs);
end

for inputIdx = 1:nInputs
    progressPercent = floor(100*(inputIdx-1)/max(nInputs,1));
    if SHOW_PROGRESS && ...
            progressPercent >= lastProgressPercent + PROGRESS_STEP_PERCENT
        fprintf(['\r%s: %3d%% | m_d inputs scanned: %d/%d', ...
                 ' | joint candidates: %d | elapsed: %.1f min'], ...
                progressLabel,progressPercent,inputIdx-1,nInputs,nHitsCounter, ...
                toc(searchTimer)/60);
        lastProgressPercent = progressPercent;
        drawnow;
    end

    m_d = uniqueMD(inputIdx);
    if runInParallel
        counterResults = counterResultsByInput{inputIdx};
        localDiagnostics = counterDiagnosticsByInput{inputIdx};
    else
        if isfinite(MAX_HITS_COUNTER)
            maxInputSolutions = MAX_HITS_COUNTER - nHitsCounter;
        else
            maxInputSolutions = inf;
        end
        if maxInputSolutions <= 0
            stopSearch = true;
            break;
        end

        [counterResults,localDiagnostics] = searchCounterForInput( ...
            m_d,S,substrateSet,productSet,badMetIdx,hasGPR,nProducers, ...
            reactionGenes, ...
            REQUIRE_GPR,STRICT_COUNTER_EDGES,ALLOW_MD_AS_MS, ...
            ALLOW_CURRENCY_MS, ...
            REQUIRE_COUNTER_GENE_DISJOINT, ...
            MAX_PRODUCERS_M1,MAX_PRODUCERS_MS,maxInputSolutions);
        counterDiagnosticsByInput{inputIdx} = localDiagnostics;
    end

    if isempty(counterResults), continue; end
    nInputsWithSolutions = nInputsWithSolutions + 1;
    nCounterMotifs = nCounterMotifs + numel(counterResults);

    block2Indices = find(inputOfBlock2 == inputIdx);
    for q2 = 1:numel(block2Indices)
        block2Idx = block2Indices(q2);
        upstream = block2Results(block2Idx);

        upstreamRxns = [ ...
            upstream.r_S1,upstream.r_Sg1,upstream.r_Ss1, ...
            upstream.r_F1,upstream.r_Fg1,upstream.r_Fs1,upstream.r_I1, ...
            upstream.r_S2,upstream.r_Sg2,upstream.r_Ss2, ...
            upstream.r_F2,upstream.r_Fg2,upstream.r_Fs2,upstream.r_I2];

        % All fields beginning with m_ are computational metabolites. The
        % explicit list avoids accidentally treating numeric indices as mets.
        upstreamMets = [ ...
            upstream.m_in,upstream.m_S1,upstream.m_F1, ...
            upstream.m_Sg1,upstream.m_Fg1, ...
            upstream.m_S2,upstream.m_F2, ...
            upstream.m_Sg2,upstream.m_Fg2];

        for qC = 1:numel(counterResults)
            counter = counterResults(qC);
            counterRxns = [counter.r_C,counter.r_CI1,counter.r_CI2];
            counterNewMets = [counter.m_1,counter.m_s];

            if REQUIRE_COUNTER_DISJOINT_RXNS && ...
                    ~isempty(intersect(upstreamRxns,counterRxns))
                continue;
            end
            if REQUIRE_COUNTER_DISJOINT_METS && ...
                    ~isempty(intersect(upstreamMets,counterNewMets))
                continue;
            end
            if REQUIRE_COUNTER_UPSTREAM_GENE_DISJOINT && ...
                    reactionsShareGenes(upstreamRxns,counterRxns, ...
                                        reactionGenes)
                continue;
            end

            nHitsCounter = nHitsCounter + 1;
            joint = upstream;
            joint.block2_idx = block2Idx;
            joint.counter_input_idx = inputIdx;

            joint.m_1 = counter.m_1;
            joint.m_s = counter.m_s;
            joint.m_c = counter.m_1;
            joint.m_Cs = counter.m_s;

            joint.r_C = counter.r_C;
            joint.r_CI1 = counter.r_CI1;
            joint.r_CI2 = counter.r_CI2;
            joint.r_I1_counter = counter.r_CI1;
            joint.r_I2_counter = counter.r_CI2;
            joint.r_Cp = counter.r_C;
            joint.r_Cr = counter.r_CI1;
            joint.r_Cg = counter.r_CI2;
            results_counter(nHitsCounter,1) = joint;

            if nHitsCounter >= MAX_HITS_COUNTER
                stopSearch = true;
                break;
            end
        end
        if stopSearch, break; end
    end
    if stopSearch, break; end
end

elapsedMinutes = toc(searchTimer)/60;
if SHOW_PROGRESS
    if stopSearch
        progressPercent = floor(100*(inputIdx-1)/max(nInputs,1));
        fprintf(['\r%s stopped at %3d%% | m_d inputs scanned: %d/%d', ...
                 ' | joint candidates: %d | elapsed: %.1f min', ...
                 ' (MAX_HITS_COUNTER reached)\n'], ...
                progressLabel,progressPercent,inputIdx-1,nInputs, ...
                nHitsCounter,elapsedMinutes);
    else
        fprintf(['\r%s: 100%% | m_d inputs scanned: %d/%d', ...
                 ' | joint candidates: %d | elapsed: %.1f min\n'], ...
                progressLabel,nInputs,nInputs,nHitsCounter,elapsedMinutes);
    end
end

%% Assemble per-input stage diagnostics
validDiagnostics = ~cellfun(@isempty,counterDiagnosticsByInput);
if any(validDiagnostics)
    diagnosticRecords = vertcat( ...
        counterDiagnosticsByInput{validDiagnostics});
    Counter_Search_Diagnostics_Block_3 = struct2table(diagnosticRecords);
    Counter_Search_Diagnostics_Block_3.m_d_ID = ...
        string(model.mets(Counter_Search_Diagnostics_Block_3.m_d));
    Counter_Search_Diagnostics_Block_3.m_d_Name = ...
        string(model.metNames(Counter_Search_Diagnostics_Block_3.m_d));
    Counter_Search_Diagnostics_Block_3 = movevars( ...
        Counter_Search_Diagnostics_Block_3,{'m_d_ID','m_d_Name'}, ...
        'After','m_d');
else
    Counter_Search_Diagnostics_Block_3 = table;
end

fprintf('\nCounter stage diagnostics by m_d:\n');
disp(Counter_Search_Diagnostics_Block_3);

%% Save beside this script and print representative solutions
Motifs_Counter_Block_3 = struct2table(results_counter);
scriptFolder = fileparts(mfilename('fullpath'));
outputFile = fullfile(scriptFolder,'Block_3_counter.mat');
save(outputFile,'Motifs_Counter_Block_3','results_counter', ...
    'Counter_Search_Diagnostics_Block_3');

fprintf('\nCounter search diagnostics:\n');
fprintf('  Distinct m_d inputs examined : %d / %d\n', ...
        min(inputIdx,nInputs),nInputs);
fprintf('  Inputs with counter solutions: %d\n',nInputsWithSolutions);
fprintf('  Counter motifs before join   : %d\n',nCounterMotifs);
fprintf('  Joint motifs                 : %d\n',nHitsCounter);
fprintf('Saved Block 3 results to:\n%s\n',outputFile);

if nHitsCounter > 0
    fprintf('\n=== First joint differentiator + counter candidates ===\n');
    for k = 1:min(5,nHitsCounter)
        rec = results_counter(k);
        fprintf('\n--- Candidate %d (Block 2 motif %d) ---\n', ...
                k,rec.block2_idx);
        fprintf('  m_d = %s | %s\n', ...
                model.mets{rec.m_d},model.metNames{rec.m_d});
        fprintf('  m_1 = m_c = %s | %s\n', ...
                model.mets{rec.m_1},model.metNames{rec.m_1});
        fprintf('  m_s = m_Cs = %s | %s\n', ...
                model.mets{rec.m_s},model.metNames{rec.m_s});

        labels = {'r_C','r_CI1','r_CI2'};
        rxns = [rec.r_C,rec.r_CI1,rec.r_CI2];
        for q = 1:numel(rxns)
            fprintf('  %-6s = %s | %s\n',labels{q}, ...
                    model.rxns{rxns(q)},model.rxnNames{rxns(q)});
            fprintf('           GPR: %s\n',model.grRules{rxns(q)});
        end
    end
end

%% Local counter search
function [counterResults,diagnostics] = searchCounterForInput( ...
        m_d,S,substrateSet,productSet,badMetIdx,hasGPR,nProducers, ...
        reactionGenes, ...
        REQUIRE_GPR,STRICT_COUNTER_EDGES,ALLOW_MD_AS_MS, ...
        ALLOW_CURRENCY_MS, ...
        REQUIRE_COUNTER_GENE_DISJOINT, ...
        MAX_PRODUCERS_M1,MAX_PRODUCERS_MS,maxHits)

counterTemplate = struct( ...
    'm_1',[],'m_s',[],'r_C',[],'r_CI1',[],'r_CI2',[]);
counterResults = repmat(counterTemplate,0,1);
nLocal = 0;

diagnostics = struct( ...
    'm_d',m_d, ...
    'nMdConsumers',0, ...
    'nM1Candidates',0, ...
    'nM1Accepted',0, ...
    'nRCI1Candidates',0, ...
    'nRCI1AfterStrict',0, ...
    'nMsCandidates',0, ...
    'nMsAccepted',0, ...
    'nRCI2BeforeStrict',0, ...
    'nRCI2AfterStrict',0, ...
    'nCompleteBeforeGene',0, ...
    'nCompleteMotifs',0);

% r_C consumes m_d and produces the counter output m_1. Reuse the same
% consumer list later for r_CI2 instead of rescanning the sparse row.
mdConsumers = find(S(m_d,:) < 0);
if REQUIRE_GPR
    mdConsumers = mdConsumers(hasGPR(mdConsumers));
end
rCCandidates = mdConsumers;
diagnostics.nMdConsumers = numel(mdConsumers);

for qC = 1:numel(rCCandidates)
    r_C = rCCandidates(qC);
    m1Candidates = setdiff(productSet{r_C},badMetIdx);
    m1Candidates = setdiff(m1Candidates,m_d);
    diagnostics.nM1Candidates = diagnostics.nM1Candidates + ...
        numel(m1Candidates);

    for q1 = 1:numel(m1Candidates)
        m_1 = m1Candidates(q1);
        if nProducers(m_1) < 1 || nProducers(m_1) > MAX_PRODUCERS_M1
            continue;
        end
        diagnostics.nM1Accepted = diagnostics.nM1Accepted + 1;

        % r_CI1 consumes m_1 together with the shared species m_s.
        rCI1Candidates = find(S(m_1,:) < 0);
        if REQUIRE_GPR
            rCI1Candidates = rCI1Candidates(hasGPR(rCI1Candidates));
        end
        rCI1Candidates = setdiff(rCI1Candidates,r_C);
        diagnostics.nRCI1Candidates = diagnostics.nRCI1Candidates + ...
            numel(rCI1Candidates);

        for qI1 = 1:numel(rCI1Candidates)
            r_CI1 = rCI1Candidates(qI1);
            if STRICT_COUNTER_EDGES && S(m_d,r_CI1) < 0
                continue;
            end
            diagnostics.nRCI1AfterStrict = ...
                diagnostics.nRCI1AfterStrict + 1;

            if ALLOW_CURRENCY_MS
                msCandidates = substrateSet{r_CI1};
            else
                msCandidates = setdiff(substrateSet{r_CI1},badMetIdx);
            end
            msCandidates = setdiff(msCandidates,m_1);
            if ~ALLOW_MD_AS_MS
                msCandidates = setdiff(msCandidates,m_d);
            end
            diagnostics.nMsCandidates = diagnostics.nMsCandidates + ...
                numel(msCandidates);

            for qS = 1:numel(msCandidates)
                m_s = msCandidates(qS);
                if nProducers(m_s) < 1 || ...
                        nProducers(m_s) > MAX_PRODUCERS_MS
                    continue;
                end
                diagnostics.nMsAccepted = diagnostics.nMsAccepted + 1;

                % r_CI2 consumes the same m_s together with m_d.
                rCI2Candidates = intersect( ...
                    mdConsumers,find(S(m_s,:) < 0));
                rCI2Candidates = setdiff(rCI2Candidates,[r_C,r_CI1]);
                diagnostics.nRCI2BeforeStrict = ...
                    diagnostics.nRCI2BeforeStrict + ...
                    numel(rCI2Candidates);
                if STRICT_COUNTER_EDGES
                    rCI2Candidates = rCI2Candidates( ...
                        S(m_1,rCI2Candidates) >= 0);
                end
                diagnostics.nRCI2AfterStrict = ...
                    diagnostics.nRCI2AfterStrict + ...
                    numel(rCI2Candidates);

                for qI2 = 1:numel(rCI2Candidates)
                    r_CI2 = rCI2Candidates(qI2);
                    if numel(unique([r_C,r_CI1,r_CI2])) < 3
                        continue;
                    end
                    diagnostics.nCompleteBeforeGene = ...
                        diagnostics.nCompleteBeforeGene + 1;
                    if REQUIRE_COUNTER_GENE_DISJOINT && ...
                            reactionsShareGenes( ...
                                [r_C,r_CI1,r_CI2],[],reactionGenes)
                        continue;
                    end

                    nLocal = nLocal + 1;
                    diagnostics.nCompleteMotifs = ...
                        diagnostics.nCompleteMotifs + 1;
                    rec = counterTemplate;
                    rec.m_1 = m_1;
                    rec.m_s = m_s;
                    rec.r_C = r_C;
                    rec.r_CI1 = r_CI1;
                    rec.r_CI2 = r_CI2;
                    counterResults(nLocal,1) = rec;

                    if nLocal >= maxHits, break; end
                end
                if nLocal >= maxHits, break; end
            end
            if nLocal >= maxHits, break; end
        end
        if nLocal >= maxHits, break; end
    end
    if nLocal >= maxHits, break; end
end

if ~isempty(counterResults)
    key = [[counterResults.m_1].',[counterResults.m_s].', ...
           [counterResults.r_C].',[counterResults.r_CI1].', ...
           [counterResults.r_CI2].'];
    [~,uniqueIdx] = unique(key,'rows','stable');
    counterResults = counterResults(uniqueIdx);
end
end

function tf = reactionsShareGenes(rxnsA,rxnsB,reactionGenes)
% True when any distinct reaction pair shares at least one GPR token.
if isempty(rxnsB)
    rxns = unique(rxnsA(:).','stable');
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
else
    rxnsA = unique(rxnsA(:).','stable');
    rxnsB = unique(rxnsB(:).','stable');
    tf = false;
    for i = 1:numel(rxnsA)
        for j = 1:numel(rxnsB)
            if ~isempty(intersect(reactionGenes{rxnsA(i)}, ...
                                  reactionGenes{rxnsB(j)}))
                tf = true;
                return;
            end
        end
    end
end
end

function counterParallelProgress(command,localHits,nInputs,stepPercent)
% Client-side progress display for DataQueue messages from parfor workers.
persistent completedInputs totalHits lastPercent progressTimer

switch command
    case 'init'
        completedInputs = 0;
        totalHits = 0;
        lastPercent = 0;
        progressTimer = tic;
        fprintf(['Counter parallel search:   0%% | m_d inputs finished: 0/%d', ...
                 ' | counter motifs: 0 | elapsed: 0.0 min'],nInputs);

    case 'tick'
        completedInputs = completedInputs + 1;
        totalHits = totalHits + localHits;
        percent = floor(100*completedInputs/max(nInputs,1));
        if percent >= lastPercent + stepPercent || completedInputs == nInputs
            fprintf(['\rCounter parallel search: %3d%%', ...
                     ' | m_d inputs finished: %d/%d', ...
                     ' | counter motifs: %d | elapsed: %.1f min'], ...
                    percent,completedInputs,nInputs,totalHits, ...
                    toc(progressTimer)/60);
            lastPercent = percent;
            drawnow;
        end

    case 'finish'
        fprintf('\n');
end
end
