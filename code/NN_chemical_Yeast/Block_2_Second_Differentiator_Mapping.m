%% Block 2: map the second chemical differentiator to yeast-GEM
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

%% Identify unique Block 1 output pairs
inputPairs = [[block1Results.m_Sg].', [block1Results.m_Fg].'];
[uniqueInputPairs,~,pairOfBlock1] = unique(inputPairs,'rows','stable');
nInputPairs = size(uniqueInputPairs,1);

fprintf('\nBlock 2 input preparation:\n');
fprintf('  Block 1 motifs                  : %d\n',nBlock1);
fprintf('  Unique (m_Sg1,m_Fg1) pairs     : %d\n',nInputPairs);

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

results_diff2 = repmat(jointTemplate,0,1);
nHits2 = 0;
stopSearch = false;
nPairsWithSolutions = 0;
nStage2Solutions = 0;

%% Serial search with exact MAX_HITS_BLOCK2 early stopping
searchTimer = tic;
lastProgressPercent = 0;
if SHOW_PROGRESS
    fprintf(['Second differentiator:   0%% | input pairs scanned: 0/%d', ...
             ' | joint candidates: 0 | elapsed: 0.0 min'],nInputPairs);
end

for pairIdx = 1:nInputPairs
    progressPercent = floor(100*(pairIdx-1)/max(nInputPairs,1));
    if SHOW_PROGRESS && ...
            progressPercent >= lastProgressPercent + PROGRESS_STEP_PERCENT
        fprintf(['\rSecond differentiator: %3d%% | input pairs scanned: %d/%d', ...
                 ' | joint candidates: %d | elapsed: %.1f min'], ...
                progressPercent,pairIdx-1,nInputPairs,nHits2, ...
                toc(searchTimer)/60);
        lastProgressPercent = progressPercent;
        drawnow;
    end

    m_Sg1 = uniqueInputPairs(pairIdx,1);
    m_Fg1 = uniqueInputPairs(pairIdx,2);

    if isfinite(MAX_HITS_BLOCK2)
        maxPairSolutions = MAX_HITS_BLOCK2 - nHits2;
    else
        maxPairSolutions = inf;
    end
    if maxPairSolutions <= 0
        stopSearch = true;
        break;
    end

    stage2Results = searchSecondDifferentiatorPair( ...
        m_Sg1,m_Fg1,S,substrateSet,productSet,badMetIdx,hasGPR, ...
        REQUIRE_GPR,REQUIRE_DISTINCT_STAGE2_METS,STRICT_STAGE2_INPUTS, ...
        MAX_REMOVERS_MS2,MAX_REMOVERS_MF2,maxPairSolutions);

    if isempty(stage2Results), continue; end
    nPairsWithSolutions = nPairsWithSolutions + 1;
    nStage2Solutions = nStage2Solutions + numel(stage2Results);

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
                          rec2.r_F2,rec2.r_Fg2,rec2.r_Fs2,rec2.r_I2];
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
    if stopSearch, break; end
end

elapsedMinutes = toc(searchTimer)/60;
if SHOW_PROGRESS
    if stopSearch
        progressPercent = floor(100*(pairIdx-1)/max(nInputPairs,1));
        fprintf(['\rSecond differentiator stopped at %3d%% | input pairs scanned: %d/%d', ...
                 ' | joint candidates: %d | elapsed: %.1f min', ...
                 ' (MAX_HITS_BLOCK2 reached)\n'], ...
                progressPercent,pairIdx-1,nInputPairs,nHits2,elapsedMinutes);
    else
        fprintf(['\rSecond differentiator: 100%% | input pairs scanned: %d/%d', ...
                 ' | joint candidates: %d | elapsed: %.1f min\n'], ...
                nInputPairs,nInputPairs,nHits2,elapsedMinutes);
    end
end

%% Save and print
Motifs_Diff_Block_2 = struct2table(results_diff2);
save('Block_2_second_differentiator.mat', ...
     'Motifs_Diff_Block_2','results_diff2');

fprintf('\nSecond differentiator diagnostics:\n');
fprintf('  Unique input pairs examined       : %d / %d\n', ...
        min(pairIdx,nInputPairs),nInputPairs);
fprintf('  Input pairs with solutions        : %d\n',nPairsWithSolutions);
fprintf('  Stage 2 motifs before Block 1 join: %d\n',nStage2Solutions);
fprintf('  Joint two-differentiator motifs   : %d\n',nHits2);
fprintf('Saved Block_2_second_differentiator.mat\n');

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
            fprintf('          GPR: %s\n',model.grRules{rxns(q)});
        end
    end
end

%% Local search functions
function stage2Results = searchSecondDifferentiatorPair( ...
        m_Sg1,m_Fg1,S,substrateSet,productSet,badMetIdx,hasGPR, ...
        REQUIRE_GPR,REQUIRE_DISTINCT_STAGE2_METS,STRICT_STAGE2_INPUTS, ...
        MAX_REMOVERS_MS2,MAX_REMOVERS_MF2,maxHits)

optionTemplate = struct('m_tracker',[],'m_reporter',[], ...
                        'r_produce',[],'r_report',[],'r_remove',[]);

slowOptions = buildSecondStageBranch( ...
    m_Sg1,m_Fg1,S,substrateSet,productSet,badMetIdx,hasGPR, ...
    REQUIRE_GPR,STRICT_STAGE2_INPUTS,MAX_REMOVERS_MS2,optionTemplate);

fastOptions = buildSecondStageBranch( ...
    m_Fg1,m_Sg1,S,substrateSet,productSet,badMetIdx,hasGPR, ...
    REQUIRE_GPR,STRICT_STAGE2_INPUTS,MAX_REMOVERS_MF2,optionTemplate);

stageTemplate = struct( ...
    'm_S2',[],'m_F2',[],'m_Sg2',[],'m_Fg2',[],'m_d',[], ...
    'r_S2',[],'r_Sg2',[],'r_Ss2',[], ...
    'r_F2',[],'r_Fg2',[],'r_Fs2',[],'r_I2',[]);
stage2Results = repmat(stageTemplate,0,1);
nLocal = 0;

for iS = 1:numel(slowOptions)
    slow = slowOptions(iS);
    for iF = 1:numel(fastOptions)
        fast = fastOptions(iF);

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

            if nLocal >= maxHits, break; end
        end
        if nLocal >= maxHits, break; end
    end
    if nLocal >= maxHits, break; end
end

% Remove exact duplicates caused by reactions with multiple valid products.
if ~isempty(stage2Results)
    key = [[stage2Results.m_S2].',[stage2Results.m_F2].', ...
           [stage2Results.m_Sg2].',[stage2Results.m_Fg2].', ...
           [stage2Results.r_S2].',[stage2Results.r_Sg2].', ...
           [stage2Results.r_Ss2].',[stage2Results.r_F2].', ...
           [stage2Results.r_Fg2].',[stage2Results.r_Fs2].', ...
           [stage2Results.r_I2].'];
    [~,uniqueIdx] = unique(key,'rows','stable');
    stage2Results = stage2Results(uniqueIdx);
end
end

function options = buildSecondStageBranch( ...
        inputMet,otherInputMet,S,substrateSet,productSet,badMetIdx,hasGPR, ...
        REQUIRE_GPR,STRICT_STAGE2_INPUTS,MAX_REMOVERS,optionTemplate)

options = repmat(optionTemplate,0,1);
nOptions = 0;

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
                end
            end
        end
    end
end

if ~isempty(options)
    key = [[options.m_tracker].',[options.m_reporter].', ...
           [options.r_produce].',[options.r_report].', ...
           [options.r_remove].'];
    [~,uniqueIdx] = unique(key,'rows','stable');
    options = options(uniqueIdx);
end
end
