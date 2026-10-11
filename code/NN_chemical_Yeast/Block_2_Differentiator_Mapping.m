%% Block 2: map the second chemical differentiator to a COBRA model
%
% This script extends the solutions saved by Block 1 with the second
% differentiator shown in the cascade diagram:
%
%   m_Sg1 --r_S2--> m_S2 --r_Sg2--> m_Sg2 (= m_I2)
%                              \--r_Ss2--> ignored products
%
%   m_Fg1 --r_F2--> m_F2 --r_Fg2--> m_Fg2 (= m_d, output)
%                              \--r_Fs2--> ignored products
%
%   m_Sg2 + m_Fg2 --r_I2--> ignored products
%
% The two inputs are fixed by each Block 1 solution. The search is performed
% once per unique (m_Sg1,m_Fg1) pair, then joined back to every compatible
% Block 1 motif. This avoids repeating the same metabolic search when many
% Block 1 motifs expose the same pair of reporter metabolites.
%
% Expected master-file sequence:
%   load Block_1_differentiator.mat
%   Block_2_Second_Differentiator_Mapping
%
% Required workspace variables:
%   model, model.S (or S), isBadMet
% Block 1 must provide either:
%   results_diff              struct array, or
%   Motifs_Diff_Block_1       table

%% Validate and prepare the workspace
if ~exist('model','var') || ~isfield(model,'S')
    error('Block2:MissingModel', ...
          'The COBRA model must be loaded before running Block 2.');
end

if ~exist('S','var'), S = model.S; end
[nMets,nRxns] = size(S); %#ok<ASGLU>

if ~exist('isBadMet','var') || numel(isBadMet) ~= nMets
    error('Block2:MissingCurrencyMask', ...
        ['isBadMet is missing or has the wrong length. Run the currency-', ...
         'metabolite setup in the master file before Block 2.']);
end
badMetIdx = find(isBadMet);

if ~exist('hasGPR','var') || numel(hasGPR) ~= nRxns
    if exist('prepareModelGPRCompatibility','file') == 2
        [model,hasGPR] = prepareModelGPRCompatibility(model);
    elseif isfield(model,'grRules') && numel(model.grRules) == nRxns
        hasGPR = false(nRxns,1);
        for r = 1:nRxns
            rule = strtrim(string(model.grRules{r}));
            hasGPR(r) = strlength(rule) > 0 && ...
                ~strcmpi(rule,"s0001") && ...
                ~strcmpi(rule,"spontaneous");
        end
    else
        error('Block2:MissingGPRData', ...
            ['hasGPR and model.grRules are unavailable. Run ', ...
             'prepareModelGPRCompatibility after loading the model.']);
    end
end

block2RuleText = strings(nRxns,1);
if isfield(model,'grRules') && numel(model.grRules) == nRxns
    block2RuleText = string(model.grRules(:));
elseif isfield(model,'rules') && numel(model.rules) == nRxns
    block2RuleText = string(model.rules(:));
elseif isfield(model,'rxnGeneMat') && size(model.rxnGeneMat,1) == nRxns
    block2RuleText(logical(hasGPR(:))) = ...
        "gene association present in rxnGeneMat";
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

if exist('results_diff','var') && isstruct(results_diff)
    block1Results = results_diff(:);
elseif exist('Motifs_Diff_Block_1','var') && istable(Motifs_Diff_Block_1)
    block1Results = table2struct(Motifs_Diff_Block_1);
    block1Results = block1Results(:);
else
    error('Block2:MissingBlock1', ...
        ['Load Block_1_differentiator.mat before running Block 2. ', ...
         'Expected results_diff or Motifs_Diff_Block_1.']);
end

requiredBlock1Fields = { ...
    'm_in','m_S','m_F','m_Sg','m_Fg', ...
    'r_S','r_Sg','r_Ss','r_F','r_Fg','r_Fs','r_I'};
missingFields = setdiff(requiredBlock1Fields, fieldnames(block1Results));
if ~isempty(missingFields)
    error('Block2:IncompleteBlock1', ...
        'Block 1 results are missing fields: %s', strjoin(missingFields,', '));
end

nBlock1 = numel(block1Results);
if nBlock1 == 0
    error('Block2:EmptyBlock1', ...
          'Block 1 contains no motifs to extend.');
end

%% Search controls
if ~exist('REQUIRE_GPR','var'), REQUIRE_GPR = true; end
if ~exist('REQUIRE_DISTINCT_STAGE2_METS','var')
    REQUIRE_DISTINCT_STAGE2_METS = true;
end
if ~exist('STRICT_STAGE2_INPUTS','var')
    % r_S2 may consume m_Sg1, but not m_Fg1; conversely for r_F2.
    STRICT_STAGE2_INPUTS = true;
end

% These optional filters are left false during broad motif discovery. They
% can be enabled later to reject reuse across the two differentiators.
if ~exist('REQUIRE_STAGE_DISJOINT_RXNS','var')
    REQUIRE_STAGE_DISJOINT_RXNS = false;
end
if ~exist('REQUIRE_STAGE_DISJOINT_METS','var')
    REQUIRE_STAGE_DISJOINT_METS = false;
end

if ~exist('MAX_REMOVERS_MS2','var'), MAX_REMOVERS_MS2 = inf; end
if ~exist('MAX_REMOVERS_MF2','var'), MAX_REMOVERS_MF2 = inf; end

if ~exist('MAX_HITS_BLOCK2','var')
    if exist('MAX_HITS','var')
        MAX_HITS_BLOCK2 = MAX_HITS;
    else
        MAX_HITS_BLOCK2 = inf;
    end
end

if ~exist('SHOW_PROGRESS','var'), SHOW_PROGRESS = true; end
if ~exist('PROGRESS_STEP_PERCENT','var'), PROGRESS_STEP_PERCENT = 1; end
PROGRESS_STEP_PERCENT = max(1,min(100,round(PROGRESS_STEP_PERCENT)));

% Parallelize only across the independent unique input pairs. An existing
% HPC pool is reused; no nested parallel loops are created.
if ~exist('USE_PARALLEL_BLOCK2','var')
    USE_PARALLEL_BLOCK2 = true;
end
if ~exist('AUTO_START_PARPOOL','var')
    AUTO_START_PARPOOL = true;
end
if ~exist('PARPOOL_WORKERS','var')
    PARPOOL_WORKERS = [];
end

% In parallel mode a global hit cap cannot interrupt workers already running.
% This separate cap protects each pair search. With a finite global limit, it
% defaults to that same value. For an otherwise unbounded search, a visible
% safety default is used; set it to inf explicitly for exhaustive enumeration.
if ~exist('MAX_STAGE2_PER_PAIR','var')
    if isfinite(MAX_HITS_BLOCK2)
        MAX_STAGE2_PER_PAIR = MAX_HITS_BLOCK2;
    else
        MAX_STAGE2_PER_PAIR = 1e4;
        warning('Block2:DefaultPerPairCap', ...
            ['MAX_HITS_BLOCK2 is infinite. MAX_STAGE2_PER_PAIR therefore ', ...
             'defaults to 10000; set it to inf for an exhaustive search.']);
    end
end

% Building slowOptions/fastOptions happens before stage-2 motifs are counted.
% Without these two budgets, a pathological input pair can remain inside the
% option builder or test millions of option pairs while the progress display
% appears frozen on the last few percent.
if ~exist('MAX_BRANCH_OPTIONS_PER_INPUT','var')
    if isfinite(MAX_STAGE2_PER_PAIR)
        MAX_BRANCH_OPTIONS_PER_INPUT = min(5000,MAX_STAGE2_PER_PAIR);
    else
        MAX_BRANCH_OPTIONS_PER_INPUT = 5000;
    end
end
if ~exist('MAX_STAGE2_OPTION_PAIRS_TESTED','var')
    MAX_STAGE2_OPTION_PAIRS_TESTED = 2e5;
end

if ~isscalar(MAX_HITS_BLOCK2) || isnan(MAX_HITS_BLOCK2) || ...
        MAX_HITS_BLOCK2 < 0
    error('Block2:InvalidMaxHits', ...
        'MAX_HITS_BLOCK2 must be a nonnegative scalar or inf.');
end
if ~isscalar(MAX_STAGE2_PER_PAIR) || isnan(MAX_STAGE2_PER_PAIR) || ...
        MAX_STAGE2_PER_PAIR < 0
    error('Block2:InvalidPerPairCap', ...
        'MAX_STAGE2_PER_PAIR must be a nonnegative scalar or inf.');
end
if ~isscalar(MAX_BRANCH_OPTIONS_PER_INPUT) || ...
        isnan(MAX_BRANCH_OPTIONS_PER_INPUT) || ...
        MAX_BRANCH_OPTIONS_PER_INPUT < 0
    error('Block2:InvalidBranchOptionCap', ...
        'MAX_BRANCH_OPTIONS_PER_INPUT must be nonnegative or inf.');
end
if ~isscalar(MAX_STAGE2_OPTION_PAIRS_TESTED) || ...
        isnan(MAX_STAGE2_OPTION_PAIRS_TESTED) || ...
        MAX_STAGE2_OPTION_PAIRS_TESTED < 0
    error('Block2:InvalidOptionPairBudget', ...
        'MAX_STAGE2_OPTION_PAIRS_TESTED must be nonnegative or inf.');
end

if isfinite(MAX_HITS_BLOCK2)
    MAX_HITS_BLOCK2 = floor(MAX_HITS_BLOCK2);
end
if isfinite(MAX_STAGE2_PER_PAIR)
    MAX_STAGE2_PER_PAIR = floor(MAX_STAGE2_PER_PAIR);
end
if isfinite(MAX_BRANCH_OPTIONS_PER_INPUT)
    MAX_BRANCH_OPTIONS_PER_INPUT = floor(MAX_BRANCH_OPTIONS_PER_INPUT);
end
if isfinite(MAX_STAGE2_OPTION_PAIRS_TESTED)
    MAX_STAGE2_OPTION_PAIRS_TESTED = ...
        floor(MAX_STAGE2_OPTION_PAIRS_TESTED);
end

%% Identify unique Block 1 output pairs
inputPairs = [[block1Results.m_Sg].', [block1Results.m_Fg].'];
[uniqueInputPairs,~,pairOfBlock1] = unique(inputPairs,'rows','stable');
nInputPairs = size(uniqueInputPairs,1);

% A corrected currency mask should have excluded ATP, ADP, protons and other
% high-degree currency metabolites before Block 1 was generated. Detect stale
% Block 1 results here instead of launching a combinatorially explosive search.
stage1SignalMets = [[block1Results.m_S].', [block1Results.m_F].', ...
                    [block1Results.m_Sg].', [block1Results.m_Fg].'];
rowsWithCurrencySignals = any(reshape( ...
    isBadMet(stage1SignalMets(:)),size(stage1SignalMets)),2);
if any(rowsWithCurrencySignals)
    firstBadRow = find(rowsWithCurrencySignals,1,'first');
    badInFirstRow = stage1SignalMets(firstBadRow, ...
        isBadMet(stage1SignalMets(firstBadRow,:)));
    badIds = strjoin(string(model.mets(badInFirstRow)),', ');
    error('Block2:CurrencySignalsInBlock1', ...
        ['%d Block 1 motifs contain currency metabolites in m_S, m_F, ', ...
         'm_Sg or m_Fg. Example row %d: %s. Correct isBadMet and rerun ', ...
         'Block 1 before Block 2.'], ...
        sum(rowsWithCurrencySignals),firstBadRow,badIds);
end

fprintf('\nBlock 2 input preparation:\n');
fprintf('  Block 1 motifs                  : %d\n',nBlock1);
fprintf('  Unique (m_Sg1,m_Fg1) pairs     : %d\n',nInputPairs);
fprintf('  Maximum stage-2 motifs per pair : %.0f\n',MAX_STAGE2_PER_PAIR);
fprintf('  Maximum branch options/input    : %.0f\n', ...
    MAX_BRANCH_OPTIONS_PER_INPUT);
fprintf('  Maximum option pairs tested     : %.0f\n', ...
    MAX_STAGE2_OPTION_PAIRS_TESTED);
fprintf('  Maximum final joined motifs     : %.0f\n',MAX_HITS_BLOCK2);

%% Joint Block 1 + Block 2 result template
jointTemplate = struct( ...
    'block1_idx',[], 'input_pair_idx',[], ...
    'm_in',[], ...
    'm_S1',[], 'm_F1',[], 'm_Sg1',[], 'm_Fg1',[], ...
    'r_S1',[], 'r_Sg1',[], 'r_Ss1',[], ...
    'r_F1',[], 'r_Fg1',[], 'r_Fs1',[], 'r_I1',[], ...
    'm_S2',[], 'm_F2',[], 'm_Sg2',[], 'm_Fg2',[], ...
    'm_I2',[], 'm_d',[], ...
    'r_S2',[], 'r_Sg2',[], 'r_Ss2',[], ...
    'r_F2',[], 'r_Fg2',[], 'r_Fs2',[], 'r_I2',[]);

if isfinite(MAX_HITS_BLOCK2)
    % Preallocate the final join so 100000 motifs do not trigger 100000
    % increasingly expensive struct-array reallocations.
    results_diff2 = repmat(jointTemplate,MAX_HITS_BLOCK2,1);
else
    results_diff2 = repmat(jointTemplate,0,1);
end
nHits2 = 0;
stopSearch = false;

%% Select parallel or serial pair search
runInParallel = false;
poolObject = [];

if USE_PARALLEL_BLOCK2 && nInputPairs > 1 && MAX_STAGE2_PER_PAIR > 0
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
        warning('Block2:ParallelFallback', ...
            'Could not start/reuse a parallel pool (%s). Using serial search.', ...
            parallelError.message);
    end
end

if USE_PARALLEL_BLOCK2 && nInputPairs > 1 && ~runInParallel
    warning('Block2:NoParallelPool', ...
        ['Parallel Block 2 search was requested, but no pool is available. ', ...
         'The input-pair search will run serially.']);
end

%% Search each unique input pair independently
stage2ResultsByPair = cell(nInputPairs,1);
pairSearchTimer = tic;

if MAX_HITS_BLOCK2 == 0 || MAX_STAGE2_PER_PAIR == 0
    fprintf('\nBlock 2 search skipped because a hit limit is zero.\n');

elseif runInParallel
    fprintf('\nRunning Block 2 pair search on %d MATLAB workers.\n', ...
        poolObject.NumWorkers);

    if SHOW_PROGRESS
        progressQueue = parallel.pool.DataQueue;
        block2ParallelProgress('init',[],nInputPairs, ...
            PROGRESS_STEP_PERCENT);
        progressListener = afterEach(progressQueue, ...
            @(pairSummary) block2ParallelProgress( ...
                'tick',pairSummary,nInputPairs,PROGRESS_STEP_PERCENT));
    else
        progressQueue = [];
        progressListener = [];
    end

    parfor pairIdx = 1:nInputPairs
        m_Sg1 = uniqueInputPairs(pairIdx,1);
        m_Fg1 = uniqueInputPairs(pairIdx,2);

        localResults = searchSecondDifferentiatorPair( ...
            m_Sg1,m_Fg1,S,substrateSet,productSet,badMetIdx,hasGPR, ...
            REQUIRE_GPR,REQUIRE_DISTINCT_STAGE2_METS, ...
            STRICT_STAGE2_INPUTS,MAX_REMOVERS_MS2,MAX_REMOVERS_MF2, ...
            MAX_STAGE2_PER_PAIR,MAX_BRANCH_OPTIONS_PER_INPUT, ...
            MAX_STAGE2_OPTION_PAIRS_TESTED);

        stage2ResultsByPair{pairIdx} = localResults;
        if SHOW_PROGRESS
            send(progressQueue,[pairIdx,numel(localResults)]);
        end
    end

    if SHOW_PROGRESS
        block2ParallelProgress('finish',[],nInputPairs, ...
            PROGRESS_STEP_PERCENT);
        delete(progressListener);
    end

else
    fprintf('\nRunning Block 2 input-pair search serially.\n');
    if SHOW_PROGRESS
        fprintf(['Second differentiator pair search:   0%%', ...
                 ' | pairs finished: 0/%d | stage-2 motifs: 0', ...
                 ' | elapsed: 0.0 min'],nInputPairs);
    end

    totalStage2SoFar = 0;
    lastProgressPercent = 0;
    for pairIdx = 1:nInputPairs
        m_Sg1 = uniqueInputPairs(pairIdx,1);
        m_Fg1 = uniqueInputPairs(pairIdx,2);

        stage2ResultsByPair{pairIdx} = searchSecondDifferentiatorPair( ...
            m_Sg1,m_Fg1,S,substrateSet,productSet,badMetIdx,hasGPR, ...
            REQUIRE_GPR,REQUIRE_DISTINCT_STAGE2_METS, ...
            STRICT_STAGE2_INPUTS,MAX_REMOVERS_MS2,MAX_REMOVERS_MF2, ...
            MAX_STAGE2_PER_PAIR,MAX_BRANCH_OPTIONS_PER_INPUT, ...
            MAX_STAGE2_OPTION_PAIRS_TESTED);

        totalStage2SoFar = totalStage2SoFar + ...
            numel(stage2ResultsByPair{pairIdx});
        progressPercent = floor(100*pairIdx/max(nInputPairs,1));
        if SHOW_PROGRESS && ...
                (progressPercent >= lastProgressPercent + ...
                 PROGRESS_STEP_PERCENT || pairIdx == nInputPairs)
            fprintf(['\rSecond differentiator pair search: %3d%%', ...
                     ' | pairs finished: %d/%d | stage-2 motifs: %d', ...
                     ' | elapsed: %.1f min'], ...
                    progressPercent,pairIdx,nInputPairs,totalStage2SoFar, ...
                    toc(pairSearchTimer)/60);
            lastProgressPercent = progressPercent;
            drawnow;
        end
    end
    if SHOW_PROGRESS
        fprintf('\n');
    end
end

pairSearchMinutes = toc(pairSearchTimer)/60;
pairSolutionCounts = cellfun(@numel,stage2ResultsByPair);
nPairsWithSolutions = sum(pairSolutionCounts > 0);
nStage2Solutions = sum(pairSolutionCounts);

if isfinite(MAX_STAGE2_PER_PAIR)
    potentiallyCappedPairs = find( ...
        pairSolutionCounts >= MAX_STAGE2_PER_PAIR & ...
        pairSolutionCounts > 0);
else
    potentiallyCappedPairs = zeros(0,1);
end

if ~isempty(potentiallyCappedPairs)
    warning('Block2:PerPairCapReached', ...
        ['%d input pair(s) returned MAX_STAGE2_PER_PAIR=%d motifs and ', ...
         'may have additional solutions. Increase the cap or set it to ', ...
         'inf for exhaustive enumeration.'], ...
        numel(potentiallyCappedPairs),MAX_STAGE2_PER_PAIR);
end

%% Join the compact pair solutions back to the full Block 1 table
joinTimer = tic;
lastJoinProgress = 0;
if SHOW_PROGRESS
    fprintf(['Full Block 1 + Block 2 join:   0%%', ...
             ' | input groups joined: 0/%d | full motifs: 0', ...
             ' | elapsed: 0.0 min'],nInputPairs);
end

for pairIdx = 1:nInputPairs
    stage2Results = stage2ResultsByPair{pairIdx};

    if ~isempty(stage2Results)
        block1Indices = find(pairOfBlock1 == pairIdx);
        for q1 = 1:numel(block1Indices)
            block1Idx = block1Indices(q1);
            rec1 = block1Results(block1Idx);

            stage1Rxns = [rec1.r_S,rec1.r_Sg,rec1.r_Ss, ...
                          rec1.r_F,rec1.r_Fg,rec1.r_Fs,rec1.r_I];
            stage1Mets = [rec1.m_in,rec1.m_S,rec1.m_F, ...
                          rec1.m_Sg,rec1.m_Fg];

            for q2 = 1:numel(stage2Results)
                rec2 = stage2Results(q2);

                stage2Rxns = [rec2.r_S2,rec2.r_Sg2,rec2.r_Ss2, ...
                              rec2.r_F2,rec2.r_Fg2,rec2.r_Fs2, ...
                              rec2.r_I2];
                stage2NewMets = [rec2.m_S2,rec2.m_F2, ...
                                 rec2.m_Sg2,rec2.m_Fg2];

                if REQUIRE_STAGE_DISJOINT_RXNS && ...
                        ~isempty(intersect(stage1Rxns,stage2Rxns))
                    continue;
                end
                if REQUIRE_STAGE_DISJOINT_METS && ...
                        ~isempty(intersect(stage1Mets,stage2NewMets))
                    continue;
                end

                nHits2 = nHits2 + 1;
                joint = jointTemplate;
                joint.block1_idx = block1Idx;
                joint.input_pair_idx = pairIdx;

                joint.m_in = rec1.m_in;
                joint.m_S1 = rec1.m_S;
                joint.m_F1 = rec1.m_F;
                joint.m_Sg1 = rec1.m_Sg;
                joint.m_Fg1 = rec1.m_Fg;
                joint.r_S1 = rec1.r_S;
                joint.r_Sg1 = rec1.r_Sg;
                joint.r_Ss1 = rec1.r_Ss;
                joint.r_F1 = rec1.r_F;
                joint.r_Fg1 = rec1.r_Fg;
                joint.r_Fs1 = rec1.r_Fs;
                joint.r_I1 = rec1.r_I;

                joint.m_S2 = rec2.m_S2;
                joint.m_F2 = rec2.m_F2;
                joint.m_Sg2 = rec2.m_Sg2;
                joint.m_Fg2 = rec2.m_Fg2;
                joint.m_I2 = rec2.m_Sg2;
                joint.m_d = rec2.m_Fg2;
                joint.r_S2 = rec2.r_S2;
                joint.r_Sg2 = rec2.r_Sg2;
                joint.r_Ss2 = rec2.r_Ss2;
                joint.r_F2 = rec2.r_F2;
                joint.r_Fg2 = rec2.r_Fg2;
                joint.r_Fs2 = rec2.r_Fs2;
                joint.r_I2 = rec2.r_I2;
                results_diff2(nHits2,1) = joint;

                if nHits2 >= MAX_HITS_BLOCK2
                    stopSearch = true;
                    break;
                end
            end
            if stopSearch, break; end
        end
    end

    joinProgress = floor(100*pairIdx/max(nInputPairs,1));
    if SHOW_PROGRESS && ...
            (joinProgress >= lastJoinProgress + PROGRESS_STEP_PERCENT || ...
             pairIdx == nInputPairs || stopSearch)
        fprintf(['\rFull Block 1 + Block 2 join: %3d%%', ...
                 ' | input groups joined: %d/%d | full motifs: %d', ...
                 ' | elapsed: %.1f min'], ...
                joinProgress,pairIdx,nInputPairs,nHits2, ...
                toc(joinTimer)/60);
        lastJoinProgress = joinProgress;
        drawnow;
    end

    if stopSearch, break; end
end

joinMinutes = toc(joinTimer)/60;
if SHOW_PROGRESS
    if stopSearch
        fprintf(' (MAX_HITS_BLOCK2 reached)\n');
    else
        fprintf('\n');
    end
end

%% Save and print
results_diff2 = results_diff2(1:nHits2);
Motifs_Diff_Block_2 = struct2table(results_diff2);

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end
outputFile = fullfile(scriptFolder,'Block_2_differentiator.mat');

save(outputFile, ...
     'Motifs_Diff_Block_2','results_diff2', ...
     'stage2ResultsByPair','uniqueInputPairs','pairOfBlock1', ...
     'pairSolutionCounts','potentiallyCappedPairs','runInParallel', ...
     'MAX_HITS_BLOCK2','MAX_STAGE2_PER_PAIR', ...
     'MAX_BRANCH_OPTIONS_PER_INPUT','MAX_STAGE2_OPTION_PAIRS_TESTED', ...
     '-v7.3');

fprintf('\nSecond differentiator diagnostics:\n');
if runInParallel
    executionMode = 'parallel';
else
    executionMode = 'serial';
end
fprintf('  Execution mode                    : %s\n',executionMode);
fprintf('  Unique input pairs examined       : %d / %d\n', ...
        nInputPairs,nInputPairs);
fprintf('  Input pairs with solutions        : %d\n',nPairsWithSolutions);
fprintf('  Stage 2 motifs before Block 1 join: %d\n',nStage2Solutions);
fprintf('  Joint two-differentiator motifs   : %d\n',nHits2);
fprintf('  Pair-search wall time             : %.1f min\n',pairSearchMinutes);
fprintf('  Full-join wall time               : %.1f min\n',joinMinutes);
fprintf('Saved Block 2 results to:\n%s\n',outputFile);

if nHits2 > 0
    fprintf('\n=== First joint two-differentiator candidates ===\n');
    for k = 1:min(5,nHits2)
        rec = results_diff2(k);
        fprintf('\n--- Candidate %d (Block 1 motif %d) ---\n', ...
                k,rec.block1_idx);
        fprintf('  Input pair: m_Sg1=%s | m_Fg1=%s\n', ...
                model.mets{rec.m_Sg1},model.mets{rec.m_Fg1});
        fprintf('  Output m_d = %s | %s\n', ...
                model.mets{rec.m_d},model.metNames{rec.m_d});

        labels = {'r_S2','r_Sg2','r_Ss2','r_F2','r_Fg2','r_Fs2','r_I2'};
        rxns = [rec.r_S2,rec.r_Sg2,rec.r_Ss2, ...
                rec.r_F2,rec.r_Fg2,rec.r_Fs2,rec.r_I2];
        for q = 1:numel(rxns)
            fprintf('  %-5s = %s | %s\n',labels{q}, ...
                    model.rxns{rxns(q)},model.rxnNames{rxns(q)});
            fprintf('          GPR: %s\n',block2RuleText(rxns(q)));
        end
    end
end

%% Local search functions
function stage2Results = searchSecondDifferentiatorPair( ...
        m_Sg1,m_Fg1,S,substrateSet,productSet,badMetIdx,hasGPR, ...
        REQUIRE_GPR,REQUIRE_DISTINCT_STAGE2_METS,STRICT_STAGE2_INPUTS, ...
        MAX_REMOVERS_MS2,MAX_REMOVERS_MF2,maxHits,maxBranchOptions, ...
        maxOptionPairsTested)

optionTemplate = struct('m_tracker',[],'m_reporter',[], ...
                        'r_produce',[],'r_report',[],'r_remove',[]);

slowOptions = buildSecondStageBranch( ...
    m_Sg1,m_Fg1,S,substrateSet,productSet,badMetIdx,hasGPR, ...
    REQUIRE_GPR,STRICT_STAGE2_INPUTS,MAX_REMOVERS_MS2, ...
    maxBranchOptions,optionTemplate);

fastOptions = buildSecondStageBranch( ...
    m_Fg1,m_Sg1,S,substrateSet,productSet,badMetIdx,hasGPR, ...
    REQUIRE_GPR,STRICT_STAGE2_INPUTS,MAX_REMOVERS_MF2, ...
    maxBranchOptions,optionTemplate);

stageTemplate = struct( ...
    'm_S2',[],'m_F2',[],'m_Sg2',[],'m_Fg2',[],'m_d',[], ...
    'r_S2',[],'r_Sg2',[],'r_Ss2',[], ...
    'r_F2',[],'r_Fg2',[],'r_Fs2',[],'r_I2',[]);
nLocal = 0;
nOptionPairsTested = 0;

if maxHits <= 0 || maxOptionPairsTested <= 0
    stage2Results = repmat(stageTemplate,0,1);
    return;
end

% Preallocation avoids repeated struct-array growth for every motif.
if isfinite(maxHits)
    stage2Results = repmat(stageTemplate,maxHits,1);
else
    stage2Results = repmat(stageTemplate,0,1);
end

for iS = 1:numel(slowOptions)
    slow = slowOptions(iS);
    for iF = 1:numel(fastOptions)
        fast = fastOptions(iF);
        nOptionPairsTested = nOptionPairsTested + 1;
        if nOptionPairsTested > maxOptionPairsTested
            stage2Results = stage2Results(1:nLocal);
            stage2Results = deduplicateStage2Results(stage2Results);
            return;
        end

        if REQUIRE_DISTINCT_STAGE2_METS && ...
                numel(unique([m_Sg1,m_Fg1,slow.m_tracker,fast.m_tracker, ...
                              slow.m_reporter,fast.m_reporter])) < 6
            continue;
        end

        selectedBranchRxns = [slow.r_produce,slow.r_report,slow.r_remove, ...
                              fast.r_produce,fast.r_report,fast.r_remove];
        if numel(unique(selectedBranchRxns)) < 6, continue; end

        rI2Candidates = intersect( ...
            find(S(slow.m_reporter,:) < 0), ...
            find(S(fast.m_reporter,:) < 0));
        if REQUIRE_GPR
            rI2Candidates = rI2Candidates(hasGPR(rI2Candidates));
        end
        rI2Candidates = setdiff(rI2Candidates,selectedBranchRxns);

        for qI = 1:numel(rI2Candidates)
            nLocal = nLocal + 1;
            rec = stageTemplate;
            rec.m_S2 = slow.m_tracker;
            rec.m_F2 = fast.m_tracker;
            rec.m_Sg2 = slow.m_reporter;
            rec.m_Fg2 = fast.m_reporter;
            rec.m_d = fast.m_reporter;
            rec.r_S2 = slow.r_produce;
            rec.r_Sg2 = slow.r_report;
            rec.r_Ss2 = slow.r_remove;
            rec.r_F2 = fast.r_produce;
            rec.r_Fg2 = fast.r_report;
            rec.r_Fs2 = fast.r_remove;
            rec.r_I2 = rI2Candidates(qI);
            stage2Results(nLocal,1) = rec;

            % Immediate per-pair stop at the actual candidate-creation point.
            if nLocal >= maxHits
                stage2Results = stage2Results(1:nLocal);
                stage2Results = deduplicateStage2Results(stage2Results);
                return;
            end
        end
    end
end

stage2Results = stage2Results(1:nLocal);
stage2Results = deduplicateStage2Results(stage2Results);
end

function options = buildSecondStageBranch( ...
        inputMet,otherInputMet,S,substrateSet,productSet,badMetIdx,hasGPR, ...
        REQUIRE_GPR,STRICT_STAGE2_INPUTS,MAX_REMOVERS,maxOptions, ...
        optionTemplate)

nOptions = 0;

if maxOptions <= 0
    options = repmat(optionTemplate,0,1);
    return;
elseif isfinite(maxOptions)
    options = repmat(optionTemplate,maxOptions,1);
else
    options = repmat(optionTemplate,0,1);
end

rProduceCandidates = find(S(inputMet,:) < 0);
if REQUIRE_GPR
    rProduceCandidates = rProduceCandidates(hasGPR(rProduceCandidates));
end
if STRICT_STAGE2_INPUTS
    rProduceCandidates = rProduceCandidates( ...
        S(otherInputMet,rProduceCandidates) >= 0);
end

for qP = 1:numel(rProduceCandidates)
    rProduce = rProduceCandidates(qP);
    trackerCandidates = setdiff(productSet{rProduce},badMetIdx);
    trackerCandidates = setdiff(trackerCandidates,[inputMet;otherInputMet]);

    for qT = 1:numel(trackerCandidates)
        tracker = trackerCandidates(qT);
        trackerConsumers = find(S(tracker,:) < 0);
        if REQUIRE_GPR
            trackerConsumers = trackerConsumers(hasGPR(trackerConsumers));
        end
        trackerConsumers = setdiff(trackerConsumers,rProduce);
        if numel(trackerConsumers) < 2, continue; end

        for qR = 1:numel(trackerConsumers)
            rReport = trackerConsumers(qR);
            reporterCandidates = setdiff(productSet{rReport},badMetIdx);
            reporterCandidates = setdiff(reporterCandidates, ...
                                         [inputMet;otherInputMet;tracker]);
            if isempty(reporterCandidates), continue; end

            rRemoveCandidates = setdiff(trackerConsumers,rReport);
            if isfinite(MAX_REMOVERS) && ...
                    numel(rRemoveCandidates) > MAX_REMOVERS
                continue;
            end

            for qD = 1:numel(rRemoveCandidates)
                rRemove = rRemoveCandidates(qD);
                if numel(unique([rProduce,rReport,rRemove])) < 3
                    continue;
                end

                for qM = 1:numel(reporterCandidates)
                    nOptions = nOptions + 1;
                    option = optionTemplate;
                    option.m_tracker = tracker;
                    option.m_reporter = reporterCandidates(qM);
                    option.r_produce = rProduce;
                    option.r_report = rReport;
                    option.r_remove = rRemove;
                    options(nOptions,1) = option;

                    % Stop during option construction, not after it finishes.
                    if nOptions >= maxOptions
                        options = options(1:nOptions);
                        options = deduplicateBranchOptions(options);
                        return;
                    end
                end
            end
        end
    end
end

options = options(1:nOptions);
options = deduplicateBranchOptions(options);
end

function stage2Results = deduplicateStage2Results(stage2Results)
if isempty(stage2Results), return; end
key = [[stage2Results.m_S2].',[stage2Results.m_F2].', ...
       [stage2Results.m_Sg2].',[stage2Results.m_Fg2].', ...
       [stage2Results.r_S2].',[stage2Results.r_Sg2].', ...
       [stage2Results.r_Ss2].',[stage2Results.r_F2].', ...
       [stage2Results.r_Fg2].',[stage2Results.r_Fs2].', ...
       [stage2Results.r_I2].'];
[~,uniqueIdx] = unique(key,'rows','stable');
stage2Results = stage2Results(uniqueIdx);
end

function options = deduplicateBranchOptions(options)
if isempty(options), return; end
key = [[options.m_tracker].',[options.m_reporter].', ...
       [options.r_produce].',[options.r_report].', ...
       [options.r_remove].'];
[~,uniqueIdx] = unique(key,'rows','stable');
options = options(uniqueIdx);
end

function block2ParallelProgress(action,pairSummary,nPairs,stepPercent)
% Client-side progress display for DataQueue messages from parfor workers.
persistent completedPairs totalMotifs lastPercent progressTimer

switch action
    case 'init'
        completedPairs = 0;
        totalMotifs = 0;
        lastPercent = 0;
        progressTimer = tic;
        fprintf(['Block 2 parallel pair search:   0%%', ...
                 ' | pairs finished: 0/%d | stage-2 motifs: 0', ...
                 ' | elapsed: 0.0 min'],nPairs);

    case 'tick'
        completedPairs = completedPairs + 1;
        totalMotifs = totalMotifs + pairSummary(2);
        percent = floor(100*completedPairs/max(nPairs,1));
        if percent >= lastPercent + stepPercent || completedPairs == nPairs
            fprintf(['\rBlock 2 parallel pair search: %3d%%', ...
                     ' | pairs finished: %d/%d | stage-2 motifs: %d', ...
                     ' | latest pair: %d | elapsed: %.1f min'], ...
                    percent,completedPairs,nPairs,totalMotifs, ...
                    pairSummary(1),toc(progressTimer)/60);
            lastPercent = percent;
            drawnow;
        end

    case 'finish'
        fprintf('\n');
end
end
