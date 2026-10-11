%% Three-reaction differentiator: X1 -> M1 -> Y1; X1 + Y1 -> waste.
% Run Minimal_Differentiator_Settings BEFORE Block 0, then this script.
% A topology match is a candidate, not a validation of kinetic parameters.
if ~exist('model','var'), error('Block1:MissingModel','Load model first.'); end
if ~exist('isBadMet','var') || numel(isBadMet) ~= size(model.S,1)
    error('Block1:CurrencyMask','Prepare isBadMet for this model first.');
end
if ~exist('InputAdapterPaths','var')
    b1Adapter = load('Block_0_Input_Adapter.mat','InputAdapterPaths');
    InputAdapterPaths = b1Adapter.InputAdapterPaths;
    clear b1Adapter;
end
if ~exist('BLOCK0_REQUIRE_EC','var'), BLOCK0_REQUIRE_EC = false; end
if ~exist('REQUIRE_GPR','var'), REQUIRE_GPR = true; end
if ~exist('MAX_HITS','var'), MAX_HITS = 1e5; end
if ~exist('MAX_INPUT_ANCHORS','var'), MAX_INPUT_ANCHORS = inf; end
if ~exist('MAX_STAGE1_PER_INPUT','var'), MAX_STAGE1_PER_INPUT = 10000; end
if ~exist('MAX_DIFF_CANDIDATE_CHECKS','var'), MAX_DIFF_CANDIDATE_CHECKS = 2e5; end
if ~exist('SHOW_PROGRESS','var'), SHOW_PROGRESS = true; end
if ~exist('USE_PARALLEL','var'), USE_PARALLEL = false; end
if ~exist('AUTO_START_PARPOOL','var'), AUTO_START_PARPOOL = false; end
if ~exist('PARPOOL_WORKERS','var'), PARPOOL_WORKERS = []; end
if ~exist('REQUIRE_DISTINCT_ADAPTER_RXNS','var'), REQUIRE_DISTINCT_ADAPTER_RXNS = true; end
if ~exist('BLOCK1_OUTPUT_DIR','var') || isempty(BLOCK1_OUTPUT_DIR), BLOCK1_OUTPUT_DIR = pwd; end

b1Settings = struct('stage',1,'maxHits',MAX_HITS, ...
    'maxInputs',MAX_INPUT_ANCHORS,'perInputCap',MAX_STAGE1_PER_INPUT, ...
    'maxChecks',MAX_DIFF_CANDIDATE_CHECKS,'requireGPR',REQUIRE_GPR, ...
    'requireAdapterEC',BLOCK0_REQUIRE_EC,'showProgress',SHOW_PROGRESS, ...
    'useParallel',USE_PARALLEL,'autoPool',AUTO_START_PARPOOL, ...
    'workers',PARPOOL_WORKERS,'distinctAdapter',REQUIRE_DISTINCT_ADAPTER_RXNS, ...
    'distinctStageRxns',false,'distinctStageMets',false,'crossTalk',true);
[Motifs_Diff_Block_1,Differentiator_Input_Anchor_Block_1, ...
    InputAdapterPaths,b1Info] = mapMinimalDifferentiator( ...
    model,isBadMet,InputAdapterPaths,table(),b1Settings);
results_diff = table2struct(Motifs_Diff_Block_1);
Block0AnchorMetIdx = unique([InputAdapterPaths.anchorMetIdx].','stable');
Block0AnchorMetIDs = string(model.mets(Block0AnchorMetIdx));
b1Checkpoint = minimalDifferentiatorCheckpoint(model,Motifs_Diff_Block_1, ...
    Differentiator_Input_Anchor_Block_1,InputAdapterPaths,b1Settings,b1Info);
b1Checkpoint.Motifs_Diff_Block_1 = Motifs_Diff_Block_1;
b1Checkpoint.Differentiator_Input_Anchor_Block_1 = Differentiator_Input_Anchor_Block_1;
saveCompactCheckpoint(fullfile(BLOCK1_OUTPUT_DIR,'Block_1_differentiator.mat'),b1Checkpoint);
fprintf('Block 1: %d connected three-reaction motifs.\n',height(Motifs_Diff_Block_1));
clear b1Checkpoint b1Settings b1Info;
