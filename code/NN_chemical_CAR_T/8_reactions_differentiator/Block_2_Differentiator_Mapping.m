%% Second three-reaction differentiator: X2=Y1 -> M2 -> Y2; Y1+Y2 -> waste.
% Searches once per distinct Y1, then connects to its original Block 1 rows.
% m_d is an alias of the PHYSICAL concentration Y2, not a branch difference.
if ~exist('model','var'), error('Block2:MissingModel','Load model first.'); end
if ~exist('isBadMet','var') || numel(isBadMet) ~= size(model.S,1)
    error('Block2:CurrencyMask','Prepare isBadMet for this model first.');
end
% Prefer a table over a potentially stale compatibility struct.
if exist('Motifs_Diff_Block_1','var') && istable(Motifs_Diff_Block_1)
    b2Upstream = Motifs_Diff_Block_1;
elseif exist('results_diff','var') && isstruct(results_diff)
    b2Upstream = struct2table(results_diff);
else
    b2Loaded = load('Block_1_differentiator.mat');
    if ~isfield(b2Loaded,'Motifs_Diff_Block_1')
        error('Block2:OldCheckpoint','Rerun the new three-reaction Block 1.');
    end
    if ~isfield(b2Loaded,'EC_ModelMetIDs') || ~isfield(b2Loaded,'EC_ModelRxnIDs') || ...
            ~isequal(b2Loaded.EC_ModelMetIDs,string(model.mets(:))) || ...
            ~isequal(b2Loaded.EC_ModelRxnIDs,string(model.rxns(:)))
        error('Block2:ModelMismatch','Block 1 was saved with a different model or index ordering.');
    end
    b2Upstream = b2Loaded.Motifs_Diff_Block_1;
    if isfield(b2Loaded,'InputAdapterPaths'), InputAdapterPaths = b2Loaded.InputAdapterPaths; end
    clear b2Loaded;
end
if ~exist('InputAdapterPaths','var')
    b2Adapter = load('Block_0_Input_Adapter.mat','InputAdapterPaths');
    InputAdapterPaths = b2Adapter.InputAdapterPaths;
    clear b2Adapter;
end
if ~exist('BLOCK0_REQUIRE_EC','var'), BLOCK0_REQUIRE_EC = false; end
if ~exist('REQUIRE_GPR','var'), REQUIRE_GPR = true; end
if ~exist('MAX_HITS_BLOCK2','var'), MAX_HITS_BLOCK2 = 1e5; end
if ~exist('MAX_STAGE2_PER_INPUT','var')
    if exist('MAX_STAGE2_PER_PAIR','var')
        MAX_STAGE2_PER_INPUT = MAX_STAGE2_PER_PAIR; % old budget remains usable
    else
        MAX_STAGE2_PER_INPUT = 10000;
    end
end
if ~exist('MAX_DIFF_CANDIDATE_CHECKS','var'), MAX_DIFF_CANDIDATE_CHECKS = 2e5; end
if ~exist('USE_PARALLEL_BLOCK2','var'), USE_PARALLEL_BLOCK2 = false; end
if ~exist('AUTO_START_PARPOOL','var'), AUTO_START_PARPOOL = false; end
if ~exist('PARPOOL_WORKERS','var'), PARPOOL_WORKERS = []; end
if ~exist('SHOW_PROGRESS','var'), SHOW_PROGRESS = true; end
if ~exist('FILTER_INTERNAL_CROSSTALK_DURING_SEARCH','var'), FILTER_INTERNAL_CROSSTALK_DURING_SEARCH = true; end
if ~exist('REQUIRE_STAGE_DISJOINT_RXNS','var'), REQUIRE_STAGE_DISJOINT_RXNS = false; end
if ~exist('REQUIRE_STAGE_DISJOINT_METS','var'), REQUIRE_STAGE_DISJOINT_METS = false; end
if ~exist('REQUIRE_DISTINCT_ADAPTER_RXNS','var'), REQUIRE_DISTINCT_ADAPTER_RXNS = true; end
if ~exist('BLOCK2_OUTPUT_DIR','var') || isempty(BLOCK2_OUTPUT_DIR), BLOCK2_OUTPUT_DIR = pwd; end
b2Settings = struct('stage',2,'maxHits',MAX_HITS_BLOCK2, ...
    'maxInputs',inf,'perInputCap',MAX_STAGE2_PER_INPUT, ...
    'maxChecks',MAX_DIFF_CANDIDATE_CHECKS,'requireGPR',REQUIRE_GPR, ...
    'requireAdapterEC',BLOCK0_REQUIRE_EC,'showProgress',SHOW_PROGRESS, ...
    'useParallel',USE_PARALLEL_BLOCK2,'autoPool',AUTO_START_PARPOOL, ...
    'workers',PARPOOL_WORKERS,'distinctAdapter',REQUIRE_DISTINCT_ADAPTER_RXNS, ...
    'distinctStageRxns',REQUIRE_STAGE_DISJOINT_RXNS, ...
    'distinctStageMets',REQUIRE_STAGE_DISJOINT_METS, ...
    'crossTalk',FILTER_INTERNAL_CROSSTALK_DURING_SEARCH);
[Motifs_Diff_Block_2,Differentiator_Input_Anchor_Block_2, ...
    InputAdapterPaths,b2Info] = mapMinimalDifferentiator( ...
    model,isBadMet,InputAdapterPaths,b2Upstream,b2Settings);
results_diff2 = table2struct(Motifs_Diff_Block_2);
b2Checkpoint = minimalDifferentiatorCheckpoint(model,Motifs_Diff_Block_2, ...
    Differentiator_Input_Anchor_Block_2,InputAdapterPaths,b2Settings,b2Info);
b2Checkpoint.Motifs_Diff_Block_2 = Motifs_Diff_Block_2;
b2Checkpoint.Differentiator_Input_Anchor_Block_2 = Differentiator_Input_Anchor_Block_2;
saveCompactCheckpoint(fullfile(BLOCK2_OUTPUT_DIR,'Block_2_differentiator.mat'),b2Checkpoint);
fprintf('Block 2: %d connected six-reaction cascades. Output: m_Y2 (=m_d).\n',height(Motifs_Diff_Block_2));
clear b2Checkpoint b2Settings b2Info b2Upstream;
