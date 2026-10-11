%% Block 1: map the first differentiator from Block 0 input adapters
%
% Complete connection:
%
%   camp[c] --Block 0 path--> m_A
%
%   m_A --r_S--> m_S --r_Sg--> m_Sg (= m_d, positive output)
%                       \--r_Ss--> ignored products
%
%   m_A --r_F--> m_F --r_Fg--> m_Fg
%                       \--r_Fs--> ignored products
%
%   m_Sg + m_Fg --r_I--> ignored products
%
% The seven-reaction search is performed once per unique m_A from Block 0.
% Each result row retains its complete Block 0 metabolite/reaction path.

%% Search controls
if ~exist('MAX_PRODUCERS_MDS','var'), MAX_PRODUCERS_MDS = inf; end
if ~exist('MAX_PRODUCERS_MDF','var'), MAX_PRODUCERS_MDF = inf; end
if ~exist('MAX_PRODUCERS_MD','var'),  MAX_PRODUCERS_MD  = inf; end
if ~exist('MAX_REMOVERS_MS','var'),   MAX_REMOVERS_MS   = inf; end
if ~exist('MAX_REMOVERS_MF','var'),   MAX_REMOVERS_MF   = inf; end
if ~exist('MAX_HITS','var'),          MAX_HITS          = inf; end
if ~exist('REQUIRE_GPR','var'),       REQUIRE_GPR       = true; end
if ~exist('REQUIRE_DISTINCT_DIFF_METS','var')
    REQUIRE_DISTINCT_DIFF_METS = true;
end
if ~exist('REQUIRE_DISTINCT_ADAPTER_RXNS','var')
    REQUIRE_DISTINCT_ADAPTER_RXNS = true;
end
if ~exist('MAX_INPUT_ANCHORS','var'), MAX_INPUT_ANCHORS = inf; end

if ~(isscalar(MAX_HITS) && MAX_HITS > 0 && ...
        (isinf(MAX_HITS) || MAX_HITS == floor(MAX_HITS)))
    error('Block1:InvalidMaxHits', ...
        'MAX_HITS must be a positive integer or inf.');
end
if ~(isscalar(MAX_INPUT_ANCHORS) && MAX_INPUT_ANCHORS > 0 && ...
        (isinf(MAX_INPUT_ANCHORS) || ...
         MAX_INPUT_ANCHORS == floor(MAX_INPUT_ANCHORS)))
    error('Block1:InvalidAnchorLimit', ...
        'MAX_INPUT_ANCHORS must be a positive integer or inf.');
end

if ~exist('SHOW_PROGRESS','var'), SHOW_PROGRESS = true; end
if ~exist('PROGRESS_STEP_PERCENT','var'), PROGRESS_STEP_PERCENT = 1; end
PROGRESS_STEP_PERCENT = max(1,min(100,round(PROGRESS_STEP_PERCENT)));

if ~exist('USE_PARALLEL','var'),       USE_PARALLEL = true; end
if ~exist('AUTO_START_PARPOOL','var'), AUTO_START_PARPOOL = true; end
if ~exist('PARPOOL_WORKERS','var'),    PARPOOL_WORKERS = []; end

%% Validate/reconstruct common model variables
if ~exist('model','var') || ~isfield(model,'S') || ...
        ~isfield(model,'mets') || ~isfield(model,'rxns')
    error('Block1:MissingModel', ...
        'model.S, model.mets and model.rxns are required.');
end

if ~exist('S','var') || ~isequal(size(S),size(model.S))
    S = model.S;
end
[nMets,nRxns] = size(S);
modelMetIds = string(model.mets(:));
modelRxnIds = string(model.rxns(:));

if numel(modelMetIds) ~= nMets || numel(modelRxnIds) ~= nRxns
    error('Block1:ModelSizeMismatch', ...
        'model.mets/model.rxns do not agree with model.S.');
end

if isfield(model,'metNames') && numel(model.metNames) == nMets
    modelMetNames = string(model.metNames(:));
else
    modelMetNames = modelMetIds;
end
if isfield(model,'rxnNames') && numel(model.rxnNames) == nRxns
    modelRxnNames = string(model.rxnNames(:));
else
    modelRxnNames = modelRxnIds;
end

% Robust GPR handling for the Puniya model and models without grRules.
if ~exist('hasGPR','var') || numel(hasGPR) ~= nRxns
    if isfield(model,'grRules') && numel(model.grRules) == nRxns
        hasGPR = strlength(strtrim(string(model.grRules(:)))) > 0;
    elseif isfield(model,'rules') && numel(model.rules) == nRxns
        hasGPR = strlength(strtrim(string(model.rules(:)))) > 0;
    elseif isfield(model,'rxnGeneMat') && size(model.rxnGeneMat,1) == nRxns
        hasGPR = full(any(model.rxnGeneMat ~= 0,2));
    else
        warning('Block1:NoGPRField', ...
            ['No usable grRules, rules or rxnGeneMat was found. ', ...
             'All reactions will pass the GPR test.']);
        hasGPR = true(nRxns,1);
    end
else
    hasGPR = logical(hasGPR(:));
end

rxnRuleText = strings(nRxns,1);
if isfield(model,'grRules') && numel(model.grRules) == nRxns
    rxnRuleText = string(model.grRules(:));
elseif isfield(model,'rules') && numel(model.rules) == nRxns
    rxnRuleText = string(model.rules(:));
elseif isfield(model,'rxnGeneMat') && size(model.rxnGeneMat,1) == nRxns
    rxnRuleText(hasGPR) = "gene association present in rxnGeneMat";
end

if ~exist('isBadMet','var') || numel(isBadMet) ~= nMets
    warning('Block1:BuildingCurrencyMask', ...
        'isBadMet is unavailable; constructing the standard mask.');
    currencyBaseIds = ["h2o","atp","adp","amp","pi","ppi","h", ...
        "nad","nadh","nadp","nadph","co2","o2","nh3","nh4", ...
        "coa","fad","fadh","fadh2"];
    metBaseIds = lower(regexprep(modelMetIds,'\[[^\]]+\]$',''));
    metBaseIds = regexprep(metBaseIds,'_([a-z0-9]+)$','');
    isBadMet = ismember(metBaseIds,currencyBaseIds);
else
    isBadMet = logical(isBadMet(:));
end

if ~exist('substrateSet','var') || numel(substrateSet) ~= nRxns
    substrateSet = cell(nRxns,1);
    for r = 1:nRxns
        substrateSet{r} = find(S(:,r) < 0);
    end
end
if ~exist('nProducers','var') || numel(nProducers) ~= nMets
    nProducers = full(sum(S > 0,2));
else
    nProducers = nProducers(:);
end

%% Load Block 0 results
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder), scriptFolder = pwd; end

haveBlock0 = exist('InputAdapterPaths','var') && ...
             exist('Block0AnchorMetIdx','var');
BLOCK0_RESULTS_FILE_USED = "workspace variables";

if ~haveBlock0
    if exist('BLOCK0_RESULTS_FILE','var') && ...
            strlength(string(BLOCK0_RESULTS_FILE)) > 0
        block0File = char(string(BLOCK0_RESULTS_FILE));
    else
        block0File = '';
        candidates = {fullfile(scriptFolder,'Block_0_Input_Adapter.mat'), ...
                      fullfile(pwd,'Block_0_Input_Adapter.mat')};
        for k = 1:numel(candidates)
            if isfile(candidates{k})
                block0File = candidates{k};
                break;
            end
        end
    end

    if isempty(block0File) || ~isfile(block0File)
        error('Block1:MissingBlock0Results', ...
            ['Block_0_Input_Adapter.mat was not found. Run Block 0 first ', ...
             'or set BLOCK0_RESULTS_FILE to its full path.']);
    end

    B0 = load(block0File);
    if ~isfield(B0,'InputAdapterPaths') || ...
            ~isfield(B0,'Block0AnchorMetIdx')
        error('Block1:InvalidBlock0Results', ...
            'The Block 0 MAT file lacks paths or anchor indices.');
    end

    InputAdapterPaths = B0.InputAdapterPaths;
    Block0AnchorMetIdx = B0.Block0AnchorMetIdx;
    if isfield(B0,'InputAdapterTable')
        InputAdapterTable = B0.InputAdapterTable;
    else
        InputAdapterTable = struct2table(InputAdapterPaths);
    end
    if isfield(B0,'Block0AnchorMetIDs')
        Block0AnchorMetIDs = B0.Block0AnchorMetIDs;
    else
        Block0AnchorMetIDs = modelMetIds(Block0AnchorMetIdx);
    end
    if isfield(B0,'Block0Parameters')
        Block0Parameters = B0.Block0Parameters;
    else
        Block0Parameters = struct;
    end
    BLOCK0_RESULTS_FILE_USED = string(block0File);
end

if isempty(InputAdapterPaths) || isempty(Block0AnchorMetIdx)
    error('Block1:NoBlock0Anchors', ...
        'Block 0 contains no paths or no anchors to extend.');
end

% Also support a master file that loaded only the two required Block 0
% variables instead of every variable saved in the MAT file.
if ~exist('InputAdapterTable','var')
    InputAdapterTable = struct2table(InputAdapterPaths);
end
if ~exist('Block0AnchorMetIDs','var')
    Block0AnchorMetIDs = modelMetIds(Block0AnchorMetIdx);
end
if ~exist('Block0Parameters','var')
    Block0Parameters = struct;
end

neededPathFields = {'anchorMetIdx','reactionIndices','inputMetIdx', ...
    'inputMetID','anchorMetID','nReactions','metabolitePath', ...
    'reactionPath','reactionNamePath'};
for k = 1:numel(neededPathFields)
    if ~isfield(InputAdapterPaths,neededPathFields{k})
        error('Block1:InvalidAdapterPaths', ...
            'InputAdapterPaths.%s is missing.',neededPathFields{k});
    end
end

inputAnchorIdx = unique(double(Block0AnchorMetIdx(:)),'stable');
if any(~isfinite(inputAnchorIdx) | inputAnchorIdx < 1 | ...
        inputAnchorIdx > nMets | inputAnchorIdx ~= floor(inputAnchorIdx))
    error('Block1:InvalidAnchorIndex', ...
        'Block0AnchorMetIdx contains an invalid model index.');
end
if isfinite(MAX_INPUT_ANCHORS) && numel(inputAnchorIdx) > MAX_INPUT_ANCHORS
    inputAnchorIdx = inputAnchorIdx(1:MAX_INPUT_ANCHORS);
    warning('Block1:AnchorLimit', ...
        'Searching only the first MAX_INPUT_ANCHORS=%d anchors.', ...
        MAX_INPUT_ANCHORS);
end

nInputAnchors = numel(inputAnchorIdx);
inputAnchorIDs = modelMetIds(inputAnchorIdx);
inputAnchorNames = modelMetNames(inputAnchorIdx);

%% Precompute branch masks for every m_A
baseBadMetIdx = find(isBadMet);
anchorBadMetIdx = cell(nInputAnchors,1);
inputConsumerRxnsByAnchor = cell(nInputAnchors,1);
inputConsumerMaskByAnchor = cell(nInputAnchors,1);
inputProductMaskByAnchor = cell(nInputAnchors,1);
nEligibleConsumers = zeros(nInputAnchors,1);
nPossibleBranchProducts = zeros(nInputAnchors,1);

fprintf('\n============================================================\n');
fprintf('BLOCK 1 - DIFFERENTIATOR FROM BLOCK 0 ANCHORS\n');
fprintf('============================================================\n');
fprintf('Block 0 source: %s\n',BLOCK0_RESULTS_FILE_USED);
fprintf('Unique m_A anchors: %d\n\n',nInputAnchors);

for a = 1:nInputAnchors
    m_A = inputAnchorIdx(a);
    localBadMetIdx = baseBadMetIdx;
    if isBadMet(m_A)
        warning('Block1:AnchorMarkedAsCurrency', ...
            ['Anchor %s is marked as currency. It is allowed as m_in ', ...
             'only, and remains excluded from every other role.'], ...
            inputAnchorIDs(a));
        localBadMetIdx = setdiff(localBadMetIdx,m_A);
    end
    anchorBadMetIdx{a} = localBadMetIdx;

    consumers = find(S(m_A,:) < 0).';
    if REQUIRE_GPR, consumers = consumers(hasGPR(consumers)); end

    consumerMask = false(nRxns,1);
    consumerMask(consumers) = true;
    productMask = false(nMets,1);
    if ~isempty(consumers)
        productMask = full(any(S(:,consumers) > 0,2));
    end
    productMask(localBadMetIdx) = false;
    productMask(m_A) = false;

    inputConsumerRxnsByAnchor{a} = consumers;
    inputConsumerMaskByAnchor{a} = consumerMask;
    inputProductMaskByAnchor{a} = productMask;
    nEligibleConsumers(a) = numel(consumers);
    nPossibleBranchProducts(a) = sum(productMask);

    fprintf('Anchor %4d/%4d: %6d  %-24s | consumers %4d | products %4d\n', ...
        a,nInputAnchors,m_A,inputAnchorIDs(a), ...
        nEligibleConsumers(a),nPossibleBranchProducts(a));
end

%% Seven-reaction search: one task = one (m_A,r_I) pair
resultTemplate = struct( ...
    'm_in',[], 'm_S',[], 'm_F',[], 'm_Sg',[], 'm_Fg',[], ...
    'r_S',[], 'r_Sg',[], 'r_Ss',[], ...
    'r_F',[], 'r_Fg',[], 'r_Fs',[], 'r_I',[], ...
    'm_Ds',[], 'm_Df',[], 'm_d',[], ...
    'r_DS',[], 'r_DF',[], 'r_Dp',[], 'r_Di',[]);

nTasks = nInputAnchors*nRxns;
resultsByTask = cell(nTasks,1);
diagnosticsByTask = zeros(nTasks,4);

runInParallel = false;
poolObject = [];
if USE_PARALLEL
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
        warning('Block1:ParallelFallback', ...
            'Parallel pool unavailable (%s); using serial search.', ...
            parallelError.message);
    end
end
if USE_PARALLEL && ~runInParallel
    warning('Block1:NoParallelPool', ...
        'Parallel execution was requested, but the search is serial.');
end

searchTimer = tic;

if runInParallel
    fprintf('\nRunning %d tasks on %d MATLAB workers.\n', ...
        nTasks,poolObject.NumWorkers);
    if isfinite(MAX_HITS)
        fprintf(['Important: parallel MAX_HITS is applied after workers ', ...
                 'finish. Use USE_PARALLEL=false for an immediate stop.\n']);
    end

    if SHOW_PROGRESS
        progressQueue = parallel.pool.DataQueue;
        block1ParallelProgress('init',0,nTasks,PROGRESS_STEP_PERCENT);
        progressListener = afterEach(progressQueue, ...
            @(localHits) block1ParallelProgress( ...
                'tick',localHits,nTasks,PROGRESS_STEP_PERCENT));
    else
        progressQueue = [];
        progressListener = [];
    end

    parfor taskIdx = 1:nTasks
        anchorOrdinal = floor((taskIdx-1)/nRxns) + 1;
        r_I = mod(taskIdx-1,nRxns) + 1;
        [localResults,localDiagnostics] = searchOneDifferentiatorAnchor( ...
            r_I,S,substrateSet,anchorBadMetIdx{anchorOrdinal}, ...
            hasGPR,nProducers,inputAnchorIdx(anchorOrdinal), ...
            inputConsumerMaskByAnchor{anchorOrdinal}, ...
            inputProductMaskByAnchor{anchorOrdinal}, ...
            REQUIRE_GPR,MAX_PRODUCERS_MD,MAX_PRODUCERS_MDS, ...
            MAX_PRODUCERS_MDF,MAX_REMOVERS_MS,MAX_REMOVERS_MF, ...
            REQUIRE_DISTINCT_DIFF_METS,MAX_HITS,resultTemplate);
        resultsByTask{taskIdx} = localResults;
        diagnosticsByTask(taskIdx,:) = localDiagnostics;
        if SHOW_PROGRESS, send(progressQueue,numel(localResults)); end
    end

    rawHitCount = sum(cellfun(@numel,resultsByTask));
    if SHOW_PROGRESS
        block1ParallelProgress('finish',rawHitCount,nTasks, ...
            PROGRESS_STEP_PERCENT);
        delete(progressListener);
    end
else
    fprintf('\nRunning %d tasks serially.\n',nTasks);
    accumulatedHits = 0;
    completedTasks = 0;
    lastPercent = -PROGRESS_STEP_PERCENT;

    for taskIdx = 1:nTasks
        anchorOrdinal = floor((taskIdx-1)/nRxns) + 1;
        r_I = mod(taskIdx-1,nRxns) + 1;
        if isfinite(MAX_HITS)
            remainingHits = max(0,MAX_HITS-accumulatedHits);
        else
            remainingHits = inf;
        end
        if remainingHits == 0, break; end

        [localResults,localDiagnostics] = searchOneDifferentiatorAnchor( ...
            r_I,S,substrateSet,anchorBadMetIdx{anchorOrdinal}, ...
            hasGPR,nProducers,inputAnchorIdx(anchorOrdinal), ...
            inputConsumerMaskByAnchor{anchorOrdinal}, ...
            inputProductMaskByAnchor{anchorOrdinal}, ...
            REQUIRE_GPR,MAX_PRODUCERS_MD,MAX_PRODUCERS_MDS, ...
            MAX_PRODUCERS_MDF,MAX_REMOVERS_MS,MAX_REMOVERS_MF, ...
            REQUIRE_DISTINCT_DIFF_METS,remainingHits,resultTemplate);
        resultsByTask{taskIdx} = localResults;
        diagnosticsByTask(taskIdx,:) = localDiagnostics;
        accumulatedHits = accumulatedHits + numel(localResults);
        completedTasks = taskIdx;

        percent = floor(100*taskIdx/max(nTasks,1));
        if SHOW_PROGRESS && (percent >= lastPercent+PROGRESS_STEP_PERCENT || ...
                taskIdx == nTasks)
            fprintf(['Block 1 progress: %3d%% | tasks %d/%d | ', ...
                     'raw candidates %d | elapsed %.1f min\n'], ...
                percent,taskIdx,nTasks,accumulatedHits,toc(searchTimer)/60);
            lastPercent = percent;
            drawnow;
        end
    end

    if isfinite(MAX_HITS) && accumulatedHits >= MAX_HITS
        fprintf(['Serial search stopped immediately at MAX_HITS=%d ', ...
                 '(task %d/%d).\n'],MAX_HITS,completedTasks,nTasks);
    end
end

nonEmptyTasks = ~cellfun(@isempty,resultsByTask);
if any(nonEmptyTasks)
    rawResults = vertcat(resultsByTask{nonEmptyTasks});
else
    rawResults = repmat(resultTemplate,0,1);
end

if isfinite(MAX_HITS) && numel(rawResults) > MAX_HITS
    rawResults = rawResults(1:MAX_HITS);
    warning('Block1:ParallelMaxHits', ...
        'Combined output was truncated to MAX_HITS=%d.',MAX_HITS);
end

% Exact differentiator deduplication before attaching adapter paths.
if ~isempty(rawResults)
    rawKey = [[rawResults.m_in].',[rawResults.m_S].', ...
        [rawResults.m_F].',[rawResults.m_Sg].',[rawResults.m_Fg].', ...
        [rawResults.r_S].',[rawResults.r_Sg].',[rawResults.r_Ss].', ...
        [rawResults.r_F].',[rawResults.r_Fg].',[rawResults.r_Fs].', ...
        [rawResults.r_I].'];
    [~,uniqueRawIdx] = unique(rawKey,'rows','stable');
    rawResults = rawResults(uniqueRawIdx);
end

%% Attach compatible Block 0 paths to Block 1 assignments
fullTemplate = resultTemplate;
fullTemplate.adapterPathIndex = [];
fullTemplate.adapterInputMetIdx = [];
fullTemplate.adapterInputMetID = "";
fullTemplate.adapterAnchorMetIdx = [];
fullTemplate.adapterAnchorMetID = "";
fullTemplate.adapterNReactions = [];
fullTemplate.adapterMetabolitePath = "";
fullTemplate.adapterReactionPath = "";
fullTemplate.adapterReactionNamePath = "";

results_diff = repmat(fullTemplate,0,1);
nRejectedReactionOverlap = 0;
adapterAnchorIdx = double([InputAdapterPaths.anchorMetIdx]).';
diffFields = fieldnames(resultTemplate);

for i = 1:numel(rawResults)
    rec = rawResults(i);
    compatiblePaths = find(adapterAnchorIdx == rec.m_in);
    diffRxns = [rec.r_S,rec.r_Sg,rec.r_Ss, ...
                rec.r_F,rec.r_Fg,rec.r_Fs,rec.r_I];

    for p = compatiblePaths(:).'
        adapterRxns = double(InputAdapterPaths(p).reactionIndices(:)).';
        if REQUIRE_DISTINCT_ADAPTER_RXNS && ...
                any(ismember(adapterRxns,diffRxns))
            nRejectedReactionOverlap = nRejectedReactionOverlap + 1;
            continue;
        end

        fullRec = fullTemplate;
        for f = 1:numel(diffFields)
            fullRec.(diffFields{f}) = rec.(diffFields{f});
        end
        if isfield(InputAdapterPaths,'pathIndex') && ...
                ~isempty(InputAdapterPaths(p).pathIndex)
            fullRec.adapterPathIndex = InputAdapterPaths(p).pathIndex;
        else
            fullRec.adapterPathIndex = p;
        end
        fullRec.adapterInputMetIdx = InputAdapterPaths(p).inputMetIdx;
        fullRec.adapterInputMetID = string(InputAdapterPaths(p).inputMetID);
        fullRec.adapterAnchorMetIdx = InputAdapterPaths(p).anchorMetIdx;
        fullRec.adapterAnchorMetID = string(InputAdapterPaths(p).anchorMetID);
        fullRec.adapterNReactions = InputAdapterPaths(p).nReactions;
        fullRec.adapterMetabolitePath = ...
            string(InputAdapterPaths(p).metabolitePath);
        fullRec.adapterReactionPath = ...
            string(InputAdapterPaths(p).reactionPath);
        fullRec.adapterReactionNamePath = ...
            string(InputAdapterPaths(p).reactionNamePath);
        results_diff(end+1,1) = fullRec; %#ok<SAGROW>

        if isfinite(MAX_HITS) && numel(results_diff) >= MAX_HITS, break; end
    end
    if isfinite(MAX_HITS) && numel(results_diff) >= MAX_HITS, break; end
end

nHits = numel(results_diff);

%% Diagnostics and full connection table
diagnosticsByAnchor = zeros(nInputAnchors,4);
for a = 1:nInputAnchors
    taskRange = (a-1)*nRxns + (1:nRxns);
    diagnosticsByAnchor(a,:) = sum(diagnosticsByTask(taskRange,:),1);
end

Differentiator_Input_Anchor_Block_1 = table( ...
    inputAnchorIdx,inputAnchorIDs,inputAnchorNames, ...
    nEligibleConsumers,nPossibleBranchProducts, ...
    diagnosticsByAnchor(:,1),diagnosticsByAnchor(:,2), ...
    diagnosticsByAnchor(:,3),diagnosticsByAnchor(:,4), ...
    'VariableNames',{'metaboliteIndex','metaboliteId','metaboliteName', ...
    'nEligibleConsumers','nPossibleBranchProducts','nRIExamined', ...
    'nReporterPairs','nProductionPairs','nBothRemovalAssignments'});

Motifs_Diff_Block_1 = struct2table(results_diff);

fprintf('\n############################################################\n');
fprintf('BLOCK 1 RESULTS\n');
fprintf('############################################################\n');
fprintf('Block 0 anchors searched            : %d\n',nInputAnchors);
fprintf('Raw differentiator assignments      : %d\n',numel(rawResults));
fprintf('Adapter/reaction overlaps rejected  : %d\n', ...
    nRejectedReactionOverlap);
fprintf('Complete connected motifs retained  : %d\n',nHits);
fprintf('Wall time                            : %.1f min\n', ...
    toc(searchTimer)/60);

if nHits > 0
    fprintf('\n=== First complete connected candidates ===\n');
    for k = 1:min(5,nHits)
        rec = results_diff(k);
        fprintf('\n--- Candidate %d ---\n',k);
        fprintf('  Adapter metabolites: %s\n',rec.adapterMetabolitePath);
        fprintf('  Adapter reactions  : %s\n',rec.adapterReactionPath);

        metLabels = {'m_in (= m_A)','m_S','m_F','m_Sg (= m_d)','m_Fg'};
        metIdx = [rec.m_in,rec.m_S,rec.m_F,rec.m_Sg,rec.m_Fg];
        for q = 1:numel(metIdx)
            fprintf('  %-14s = %s | %s\n',metLabels{q}, ...
                modelMetIds(metIdx(q)),modelMetNames(metIdx(q)));
        end

        rxnLabels = {'r_S','r_Sg','r_Ss','r_F','r_Fg','r_Fs','r_I'};
        rxnIdx = [rec.r_S,rec.r_Sg,rec.r_Ss, ...
                  rec.r_F,rec.r_Fg,rec.r_Fs,rec.r_I];
        for q = 1:numel(rxnIdx)
            fprintf('  %-4s = %s | %s\n',rxnLabels{q}, ...
                modelRxnIds(rxnIdx(q)),modelRxnNames(rxnIdx(q)));
            fprintf('       GPR: %s\n',rxnRuleText(rxnIdx(q)));
        end
    end
else
    warning('Block1:NoConnectedMotifs', ...
        ['No complete connected adapter+differentiator motif was found. ', ...
         'Inspect Differentiator_Input_Anchor_Block_1 by anchor.']);
end

%% Save
if ~exist('BLOCK1_OUTPUT_DIR','var') || isempty(BLOCK1_OUTPUT_DIR)
    BLOCK1_OUTPUT_DIR = scriptFolder;
end
BLOCK1_OUTPUT_DIR = char(string(BLOCK1_OUTPUT_DIR));
if ~isfolder(BLOCK1_OUTPUT_DIR), mkdir(BLOCK1_OUTPUT_DIR); end
outputFile = fullfile(BLOCK1_OUTPUT_DIR,'Block_1_differentiator.mat');

save(outputFile,'Motifs_Diff_Block_1','results_diff','rawResults', ...
    'Differentiator_Input_Anchor_Block_1','InputAdapterPaths', ...
    'InputAdapterTable','Block0AnchorMetIdx','Block0AnchorMetIDs', ...
    'Block0Parameters','BLOCK0_RESULTS_FILE_USED','inputAnchorIdx', ...
    'inputAnchorIDs','inputConsumerRxnsByAnchor', ...
    'REQUIRE_DISTINCT_ADAPTER_RXNS','-v7.3');

fprintf('\nSaved connected Block 1 results to:\n%s\n',outputFile);
fprintf('############################################################\n');

%% Local search function
function [localResults,diagnostics] = searchOneDifferentiatorAnchor( ...
        r_I,S,substrateSet,badMetIdx,hasGPR,nProducers, ...
        inputMetIdx,inputConsumerMask,inputProductMask, ...
        REQUIRE_GPR,MAX_PRODUCERS_MD,MAX_PRODUCERS_MDS, ...
        MAX_PRODUCERS_MDF,MAX_REMOVERS_MS,MAX_REMOVERS_MF, ...
        REQUIRE_DISTINCT_DIFF_METS,localHitLimit,resultTemplate)

localResults = repmat(resultTemplate,0,1);
localHits = 0;
% [r_I examined, reporter pairs, production pairs, remover assignments]
diagnostics = zeros(1,4);

if REQUIRE_GPR && ~hasGPR(r_I), return; end
reporterSubs = setdiff(substrateSet{r_I},badMetIdx);
if numel(reporterSubs) < 2, return; end
diagnostics(1) = 1;

for iSg = 1:numel(reporterSubs)
    m_Sg = reporterSubs(iSg);
    if nProducers(m_Sg) < 1 || nProducers(m_Sg) > MAX_PRODUCERS_MD
        continue;
    end
    rSgCandidates = find(S(m_Sg,:) > 0);
    if REQUIRE_GPR, rSgCandidates = rSgCandidates(hasGPR(rSgCandidates)); end
    rSgCandidates(rSgCandidates == r_I) = [];

    for iFg = 1:numel(reporterSubs)
        if iFg == iSg, continue; end
        m_Fg = reporterSubs(iFg);
        rFgCandidates = find(S(m_Fg,:) > 0);
        if REQUIRE_GPR
            rFgCandidates = rFgCandidates(hasGPR(rFgCandidates));
        end
        rFgCandidates(rFgCandidates == r_I) = [];

        for a = 1:numel(rSgCandidates)
            r_Sg = rSgCandidates(a);
            mSCandidates = setdiff(substrateSet{r_Sg},badMetIdx);
            mSCandidates = setdiff(mSCandidates,[m_Sg;m_Fg]);
            mSCandidates = mSCandidates(inputProductMask(mSCandidates));

            for b = 1:numel(rFgCandidates)
                r_Fg = rFgCandidates(b);
                if r_Fg == r_Sg, continue; end
                mFCandidates = setdiff(substrateSet{r_Fg},badMetIdx);
                mFCandidates = setdiff(mFCandidates,[m_Sg;m_Fg]);
                mFCandidates = mFCandidates(inputProductMask(mFCandidates));

                for iS = 1:numel(mSCandidates)
                    m_S = mSCandidates(iS);
                    if nProducers(m_S) < 1 || ...
                            nProducers(m_S) > MAX_PRODUCERS_MDS
                        continue;
                    end
                    rSCandidates = find(S(m_S,:) > 0);
                    if REQUIRE_GPR
                        rSCandidates = rSCandidates(hasGPR(rSCandidates));
                    end
                    rSCandidates = setdiff(rSCandidates,[r_I,r_Sg,r_Fg]);
                    rSCandidates = rSCandidates(inputConsumerMask(rSCandidates));
                    if isempty(rSCandidates), continue; end

                    for iF = 1:numel(mFCandidates)
                        m_F = mFCandidates(iF);
                        if m_F == m_S, continue; end
                        if nProducers(m_F) < 1 || ...
                                nProducers(m_F) > MAX_PRODUCERS_MDF
                            continue;
                        end
                        if REQUIRE_DISTINCT_DIFF_METS && ...
                                numel(unique([m_S,m_F,m_Sg,m_Fg])) < 4
                            continue;
                        end
                        diagnostics(2) = diagnostics(2) + 1;

                        rFCandidates = find(S(m_F,:) > 0);
                        if REQUIRE_GPR
                            rFCandidates = rFCandidates(hasGPR(rFCandidates));
                        end
                        rFCandidates = setdiff(rFCandidates,[r_I,r_Sg,r_Fg]);
                        rFCandidates = rFCandidates(inputConsumerMask(rFCandidates));
                        if isempty(rFCandidates), continue; end

                        rSsCandidates = find(S(m_S,:) < 0);
                        if REQUIRE_GPR
                            rSsCandidates = rSsCandidates(hasGPR(rSsCandidates));
                        end
                        rSsCandidates = setdiff(rSsCandidates,[r_I,r_Sg,r_Fg]);
                        if isfinite(MAX_REMOVERS_MS) && ...
                                numel(rSsCandidates) > MAX_REMOVERS_MS
                            continue;
                        end

                        rFsCandidates = find(S(m_F,:) < 0);
                        if REQUIRE_GPR
                            rFsCandidates = rFsCandidates(hasGPR(rFsCandidates));
                        end
                        rFsCandidates = setdiff(rFsCandidates,[r_I,r_Sg,r_Fg]);
                        if isfinite(MAX_REMOVERS_MF) && ...
                                numel(rFsCandidates) > MAX_REMOVERS_MF
                            continue;
                        end
                        if isempty(rSsCandidates) || isempty(rFsCandidates)
                            continue;
                        end

                        for c = 1:numel(rSCandidates)
                            r_S = rSCandidates(c);
                            for d = 1:numel(rFCandidates)
                                r_F = rFCandidates(d);
                                if r_F == r_S, continue; end
                                diagnostics(3) = diagnostics(3) + 1;

                                for e = 1:numel(rSsCandidates)
                                    r_Ss = rSsCandidates(e);
                                    if ismember(r_Ss,[r_S,r_F,r_Sg,r_Fg,r_I])
                                        continue;
                                    end
                                    for f = 1:numel(rFsCandidates)
                                        r_Fs = rFsCandidates(f);
                                        selectedRxns = [r_S,r_Sg,r_Ss, ...
                                                        r_F,r_Fg,r_Fs,r_I];
                                        if numel(unique(selectedRxns)) < 7
                                            continue;
                                        end
                                        diagnostics(4) = diagnostics(4) + 1;

                                        m_in = inputMetIdx;
                                        if REQUIRE_DISTINCT_DIFF_METS && ...
                                                numel(unique([m_in,m_S,m_F, ...
                                                            m_Sg,m_Fg])) < 5
                                            continue;
                                        end

                                        localHits = localHits + 1;
                                        rec = resultTemplate;
                                        rec.m_in = m_in;
                                        rec.m_S = m_S;
                                        rec.m_F = m_F;
                                        rec.m_Sg = m_Sg;
                                        rec.m_Fg = m_Fg;
                                        rec.r_S = r_S;
                                        rec.r_Sg = r_Sg;
                                        rec.r_Ss = r_Ss;
                                        rec.r_F = r_F;
                                        rec.r_Fg = r_Fg;
                                        rec.r_Fs = r_Fs;
                                        rec.r_I = r_I;
                                        rec.m_Ds = m_S;
                                        rec.m_Df = m_F;
                                        rec.m_d = m_Sg;
                                        rec.r_DS = r_S;
                                        rec.r_DF = r_F;
                                        rec.r_Dp = r_Sg;
                                        rec.r_Di = r_I;
                                        localResults(localHits,1) = rec;
                                        if localHits >= localHitLimit, return; end
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

function block1ParallelProgress(action,localHits,nTasks,stepPercent)
persistent completedTasks accumulatedHits progressTimer lastPercent
switch action
    case 'init'
        completedTasks = 0;
        accumulatedHits = 0;
        progressTimer = tic;
        lastPercent = 0;
        fprintf(['Block 1 progress:   0%% | tasks 0/%d | ', ...
                 'raw candidates 0 | elapsed 0.0 min\n'],nTasks);
    case 'tick'
        completedTasks = completedTasks + 1;
        accumulatedHits = accumulatedHits + localHits;
        percent = floor(100*completedTasks/max(nTasks,1));
        if percent >= lastPercent+stepPercent || completedTasks == nTasks
            fprintf(['Block 1 progress: %3d%% | tasks %d/%d | ', ...
                     'raw candidates %d | elapsed %.1f min\n'], ...
                percent,completedTasks,nTasks,accumulatedHits, ...
                toc(progressTimer)/60);
            lastPercent = percent;
            drawnow;
        end
    case 'finish'
        fprintf(['Block 1 progress: 100%% | tasks %d/%d | ', ...
                 'raw candidates %d | elapsed %.1f min\n'], ...
            nTasks,nTasks,localHits,toc(progressTimer)/60);
end
end
