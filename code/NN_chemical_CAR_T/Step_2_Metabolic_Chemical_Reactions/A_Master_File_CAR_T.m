%% CAR-T master: full minimal differentiator/counter/Hill-NN structural search.
% This is a separate entry point; the original full master is preserved.
clearvars;
clc;
addingPathParentFolderByName('code');
masterFolder = fileparts(mfilename('fullpath'));
addpath(masterFolder);

%% Configuration (all result files go to the current folder)
MODEL_FILE = 'Puniya2021_TNM1055.xml';
INITIALIZE_COBRA_TOOLBOX = true; % false if your HPC job initialized it
RUN_BLOCK0 = true;             % false loads the local adapter checkpoint
RUN_BLOCK1 = true;             % false loads NEW three-reaction results
RUN_BLOCK2 = true;             % false loads NEW six-reaction results
EXCLUDE_CURRENCY = true;
REQUIRE_GPR = true;
BLOCK0_REQUIRE_EC = false;
BLOCK0_MAX_DEPTH = 10;
BLOCK0_MAX_PATHS = 1e5;
BLOCK0_MAX_SHORTEST_PATHS_PER_MET = 3;
BLOCK0_CONTINUE_NONCURRENCY = true;
BLOCK0_REQUIRE_BRANCHABLE_ANCHOR = false; % no old two-branch test
BLOCK0_OUTPUT_DIR = pwd;
BLOCK0_SAVE_BESIDE_SCRIPT = false;
MAX_HITS = 1e5;
MAX_HITS_BLOCK2 = 1e5;
MAX_INPUT_ANCHORS = inf;
MAX_STAGE1_PER_INPUT = 10000;
MAX_STAGE2_PER_INPUT = 10000;
MAX_DIFF_CANDIDATE_CHECKS = 2e5;
USE_PARALLEL = false;
USE_PARALLEL_BLOCK2 = false;   % true reuses your allocated HPC pool
AUTO_START_PARPOOL = false;
PARPOOL_WORKERS = [];
REQUIRE_DISTINCT_ADAPTER_RXNS = true;
REQUIRE_STAGE_DISJOINT_RXNS = true;
REQUIRE_STAGE_DISJOINT_METS = true; % Y1=X2 is still permitted
FILTER_INTERNAL_CROSSTALK_DURING_SEARCH = true;
SHOW_PROGRESS = true;
PROGRESS_STEP_PERCENT = 1;
BLOCK1_OUTPUT_DIR = pwd;
BLOCK2_OUTPUT_DIR = pwd;
PRINT_SELECTED_DIFFERENTIATOR = false; % printer defaults to motif 1

%% Downstream policy and budgets (see README for regulatory evidence mode)
REGULATION_MODE = "design";
REGULATION_FILE = '';
RUN_BLOCK3 = true;
RUN_BLOCK4 = true;
RUN_BLOCK5 = true;
RUN_BLOCK6 = true;
MAX_HITS_COUNTER = 1e5;
MAX_HITS_NN = 1e5;
MAX_COUNTER_POOL = 2000;
MAX_NN_POOL = 2000;
MAX_DOWNSTREAM_CANDIDATE_CHECKS = 2e6;
MAX_DOWNSTREAM_JOIN_CHECKS = 5e6;
MAX_OPTIONS_PER_SPECIES = 50;
MAX_NN_NODE_OPTIONS = 2000;
NN_NODE_OPTIONS_PER_ROLE = 40;
ALLOW_CURRENCY_COUNTER_CS = false;
FILTER_ADAPTER_MEMBRANE = true;
MOTIF_TO_PRINT = 1;
PRINT_SELECTED_MOTIF = true;

if INITIALIZE_COBRA_TOOLBOX
    if exist('initCobraToolbox','file') ~= 2
        error('Master:MissingCOBRA', ...
            ['initCobraToolbox was not found. Add the COBRA Toolbox to ', ...
             'the MATLAB path or set INITIALIZE_COBRA_TOOLBOX=false ', ...
             'after initializing it in the HPC job.']);
    end
    initCobraToolbox(false);
end

%% ================================================================
%% Load the Puniya T-cell metabolic model
%% ================================================================

modelFilePath = '';
modelCandidates = {MODEL_FILE,fullfile(masterFolder,MODEL_FILE), ...
                   which(MODEL_FILE)};
for k = 1:numel(modelCandidates)
    candidate = modelCandidates{k};
    if ~isempty(candidate) && isfile(candidate)
        modelFilePath = candidate;
        break;
    end
end

if isempty(modelFilePath)
    error('Master:MissingModelFile', ...
        ['Could not find %s. Put it beside this master file, use an ', ...
         'absolute MODEL_FILE path, or add its directory to MATLAB path.'], ...
        MODEL_FILE);
end

fprintf('\nLoading T-cell model:\n  %s\n',modelFilePath);
model = readCbModel(modelFilePath);

if ~isfield(model,'S') || ~isfield(model,'mets') || ...
        ~isfield(model,'rxns')
    error('Master:InvalidModel', ...
        'The loaded COBRA model must contain S, mets and rxns.');
end

S = model.S;
[nMets,nRxns] = size(S);
fprintf('T-cell model loaded: %d metabolites, %d reactions.\n', ...
    nMets,nRxns);

%% ================================================================
%% Normalize GPR information
%% ================================================================

if exist('prepareModelGPRCompatibility','file') == 2
    [model,hasGPR,gprSource] = prepareModelGPRCompatibility(model);
else
    warning('Master:MissingGPRHelper', ...
        ['prepareModelGPRCompatibility.m was not found. Using available ', ...
         'model GPR fields directly.']);

    if isfield(model,'grRules') && numel(model.grRules) == nRxns
        model.grRules = cellstr(string(model.grRules(:)));
    elseif isfield(model,'rules') && numel(model.rules) == nRxns
        model.grRules = cellstr(string(model.rules(:)));
    else
        model.grRules = repmat({''},nRxns,1);
    end

    hasGPR = strlength(strtrim(string(model.grRules(:)))) > 0;
    gprSource = "master-file fallback";
end

hasGPR = logical(hasGPR(:));
fprintf('GPR source retained by master: %s\n',gprSource);
[hasEC,EC_Codes,EC_AnnotationFields] = getReactionECAnnotations(model);
fprintf('Complete EC-annotated reactions: %d / %d\n',sum(hasEC),nRxns);
if ~any(hasEC)
    error('Master:NoECAnnotations', ...
        ['No complete EC codes were found in eccodes, rxnECNumbers or ecNumbers. ', ...
         'Check the loaded model annotations before searching.']);
end

%% ================================================================
%% Currency-metabolite mask
%% ================================================================

currencyBaseIDs = [ ...
    "h2o","atp","adp","amp","pi","ppi","h", ...
    "nad","nadh","nadp","nadph","co2","o2", ...
    "nh3","nh4","coa","fad","fadh","fadh2", ...
    "gtp","gdp","dpcoa"];

modelMetIDs = lower(string(model.mets(:)));
metBaseIDs = regexprep(modelMetIDs,'\[[^\]]+\]$','');
metBaseIDs = regexprep(metBaseIDs, ...
    '_(c|e|m|n|r|g|x|l|v|p|u|ce|pm)$','');

isBadMet = false(nMets,1);
if EXCLUDE_CURRENCY
    isBadMet = ismember(metBaseIDs,currencyBaseIDs);
end

fprintf('Currency metabolites excluded from signal roles: %d / %d\n', ...
    sum(isBadMet),nMets);

%% ================================================================
%% Resolve the cAMP input and AMP adapter bridge
%% ================================================================

inputCandidates = ["camp_c","camp[c]"];
inputMatch = find(ismember(modelMetIDs,lower(inputCandidates)));
if numel(inputMatch) ~= 1
    error('Master:CAMPResolution', ...
        ['Expected exactly one cytoplasmic cAMP metabolite matching ', ...
         'camp_c or camp[c], but found %d.'],numel(inputMatch));
end
BLOCK0_INPUT_MET_ID = string(model.mets(inputMatch));

bridgeCandidates = ["amp_c","amp[c]"];
bridgeMatch = find(ismember(modelMetIDs,lower(bridgeCandidates)));
if numel(bridgeMatch) ~= 1
    error('Master:AMPResolution', ...
        ['Expected exactly one cytoplasmic AMP metabolite matching ', ...
         'amp_c or amp[c], but found %d.'],numel(bridgeMatch));
end
BLOCK0_BRIDGE_MET_IDS = string(model.mets(bridgeMatch));

fprintf('Block 0 input anchor : %s\n',BLOCK0_INPUT_MET_ID);
fprintf('Block 0 AMP bridge   : %s\n',BLOCK0_BRIDGE_MET_IDS);



%% Run or resume the three blocks
if RUN_BLOCK0
    Block_0_Input_Adapter;
else
    adapterSaved = load('Block_0_Input_Adapter.mat','InputAdapterPaths');
    InputAdapterPaths = adapterSaved.InputAdapterPaths;
    clear adapterSaved;
end

if isempty(InputAdapterPaths)
    fprintf('Block 0 has no adapter paths; no differentiator search started.\n');
    return;
end

if RUN_BLOCK1
    Block_1_Differentiator_Mapping;
else
    b1Saved = localLoadMinimal('Block_1_differentiator.mat',model);
    Motifs_Diff_Block_1 = b1Saved.Motifs_Diff_Block_1;
    InputAdapterPaths = b1Saved.InputAdapterPaths;
    results_diff = table2struct(Motifs_Diff_Block_1);
    clear b1Saved;
end
if isempty(Motifs_Diff_Block_1)
    fprintf('Block 1 found no three-reaction motif; inspect its input diagnostics.\n');
    return;
end

if RUN_BLOCK2
    Block_2_Differentiator_Mapping;
else
    b2Saved = localLoadMinimal('Block_2_differentiator.mat',model);
    Motifs_Diff_Block_2 = b2Saved.Motifs_Diff_Block_2;
    results_diff2 = table2struct(Motifs_Diff_Block_2);
    clear b2Saved;
end

Continue_Minimal_Pipeline;

function saved = localLoadMinimal(filename,model)
saved = load(filename);
if ~isfield(saved,'TopologyVersion') || saved.TopologyVersion ~= "minimal_three_reaction_v1"
    error('MinimalMaster:OldFile','%s uses the old topology; regenerate it.',filename);
end
if ~isfield(saved,'EC_ModelMetIDs') || ~isfield(saved,'EC_ModelRxnIDs') || ...
        ~isequal(saved.EC_ModelMetIDs,string(model.mets(:))) || ...
        ~isequal(saved.EC_ModelRxnIDs,string(model.rxns(:)))
    error('MinimalMaster:ModelMismatch','%s belongs to a different model/index order.',filename);
end
end



