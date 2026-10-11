%% MASTER FILE - cAMP-ANCHORED T-CELL METABOLIC COMPUTATION
%
% Pipeline:
%   Block 0 : cAMP input adapter
%   Block 1 : first chemical differentiator
%   Block 2 : second chemical differentiator
%   Block 3 : counter
%   Block 4 : single-neuron NN
%   Block 5 : remove complete motifs containing outer-cell-membrane reactions
%   Block 6 : remove motifs containing internal cross-talk
%   Block 6B: verify complete EC annotations, including the input adapter
%   Block 7 : print one selected cAMP-to-NN motif
%
% Put this master file and the eight Block_*.m files in the same directory.
% The Puniya2021_TNM1055.xml model may be in that directory or anywhere on
% the MATLAB path. Each block saves its own MAT file beside its script.

clearvars;
clc;
addingPathParentFolderByName('code');

%% ================================================================
%% User configuration
%% ================================================================

MODEL_FILE = 'Puniya2021_TNM1055.xml';

% Set false if the COBRA Toolbox has already been initialized by the HPC job.
INITIALIZE_COBRA_TOOLBOX = true;

% Execution controls. If a block is false, its previously saved MAT file is
% loaded so that the following block receives the same workspace variables.
RUN_BLOCK0 = true; % Regenerate adapters under the configured passive-input policy.
RUN_BLOCK1 = true;
RUN_BLOCK2 = true;
RUN_BLOCK3 = true;
RUN_BLOCK4 = true;
RUN_BLOCK5 = true;
RUN_BLOCK6 = true;
RUN_BLOCK7 = true;

% General search controls
EXCLUDE_CURRENCY = true;
REQUIRE_GPR = true;
% Complete EC annotations are required in Blocks 1-4, independently
% of REQUIRE_GPR. Rebuild Blocks 0-4 once to obtain the new search coverage.
SHOW_PROGRESS = true;
PROGRESS_STEP_PERCENT = 1;

% Reject internal cross-talk while Blocks 2-4 assemble candidates. The
% MAX_HITS limits then count accepted motifs instead of raw motifs that would
% all be removed later by Block 6. Reactions outside the selected algorithm
% are not tested, and extra network producers remain allowed.
FILTER_INTERNAL_CROSSTALK_DURING_SEARCH = true;
MAX_CROSSTALK_TESTS_BLOCK3 = 1e7;
MAX_CROSSTALK_TESTS_BLOCK4 = 1e7;

% Parallel pool controls. On the HPC system, configure the MATLAB cluster
% profile before executing this file. Leave PARPOOL_WORKERS empty to use the
% worker count assigned by the active profile/job.
AUTO_START_PARPOOL = true;
PARPOOL_WORKERS = [];
USE_PARALLEL = false;         % Block 1
USE_PARALLEL_BLOCK2 = true;   % Block 2
USE_PARALLEL_COUNTER = true;  % Block 3
USE_PARALLEL_NN = true;       % Block 4

% Block 0: cAMP adapter
% Passive adapter reactions need no EC annotation for enzyme control.
BLOCK0_REQUIRE_EC = false;
BLOCK0_MAX_DEPTH = 10;
BLOCK0_MAX_PATHS = 1e5;
BLOCK0_MAX_SHORTEST_PATHS_PER_MET = 3;
BLOCK0_CONTINUE_NONCURRENCY = true;
BLOCK0_REQUIRE_BRANCHABLE_ANCHOR = true;
BLOCK0_SAVE_BESIDE_SCRIPT = true;

% Block 1: first differentiator
MAX_HITS = 1e5;
MAX_INPUT_ANCHORS = inf;
MAX_PRODUCERS_MD = inf;
MAX_PRODUCERS_MDF = inf;
MAX_PRODUCERS_MDS = inf;
MAX_REMOVERS_MS = inf;
MAX_REMOVERS_MF = inf;
REQUIRE_DISTINCT_DIFF_METS = true;
REQUIRE_DISTINCT_ADAPTER_RXNS = true;

% Block 2: second differentiator. The three finite inner limits prevent a
% small number of pathological metabolite pairs from holding the search at
% 93-99 percent indefinitely.
MAX_HITS_BLOCK2 = 1e5;
MAX_STAGE2_PER_PAIR = 2000;
MAX_BRANCH_OPTIONS_PER_INPUT = 2000;
MAX_STAGE2_OPTION_PAIRS_TESTED = 2e5;
MAX_REMOVERS_MS2 = inf;
MAX_REMOVERS_MF2 = inf;
STRICT_STAGE2_INPUTS = true;
REQUIRE_DISTINCT_STAGE2_METS = true;
REQUIRE_STAGE_DISJOINT_RXNS = false;
REQUIRE_STAGE_DISJOINT_METS = false;

% Block 3: counter
MAX_HITS_COUNTER = 1e5;
MAX_COUNTER_CANDIDATES_PER_INPUT = 1e5;
MAX_PRODUCERS_M1 = inf;
MAX_PRODUCERS_MS = inf;
STRICT_COUNTER_EDGES = true;
ALLOW_MD_AS_MS = false;
ALLOW_CURRENCY_MS = false;
REQUIRE_COUNTER_DISJOINT_RXNS = false;
REQUIRE_COUNTER_DISJOINT_METS = false;
REQUIRE_COUNTER_GENE_DISJOINT = false;
REQUIRE_COUNTER_UPSTREAM_GENE_DISJOINT = false;

% Block 4: single-neuron NN
MAX_HITS_NN = 1e5;
MAX_NN_CANDIDATES_PER_INPUT = 1e5;
MAX_FULL_MOTIFS_NN = 1e5;
FULL_TABLE_SAFETY_CAP = 1e5;
ALLOW_UNBOUNDED_FULL_NN = false;
MAX_PRODUCERS_MOUT = inf;
STRICT_BIAS_REACTION = true;
REQUIRE_BIAS_NONCURRENCY_SUBSTRATE = false;
REQUIRE_NN_GENE_DISJOINT = false;

% Block 5: outer-cell-membrane filter
SHOW_MEMBRANE_FILTER_PROGRESS = true;
MEMBRANE_PROGRESS_STEP_PERCENT = 10;
PRINT_REJECTED_MOTIFS = false;
MAX_REJECTED_TO_PRINT = 20;
PRINT_RETAINED_INDICES = false;

% Block 6: internal cross-talk filter
SHOW_CROSSTALK_FILTER_PROGRESS = true;
CROSSTALK_PROGRESS_STEP_PERCENT = 10;
PRINT_REJECTED_CROSSTALK_MOTIFS = false;
MAX_CROSSTALK_REJECTED_TO_PRINT = 20;

% Block 7: selected-motif printer
MOTIF_TO_PRINT = 1;
ADAPTER_MATCH_TO_PRINT = 1;
% "auto" prints from Block 6 when strict survivors exist and otherwise
% falls back to Block 5. Use "block6" to forbid that fallback.
MOTIF_SOURCE = "block6";

%% ================================================================
%% Locate scripts and initialize COBRA
%% ================================================================

masterFolder = fileparts(mfilename('fullpath'));
if isempty(masterFolder)
    masterFolder = pwd;
end

BLOCK0_SCRIPT = fullfile(masterFolder,'Block_0_Input_Adapter.m');
BLOCK1_SCRIPT = fullfile(masterFolder,'Block_1_Differentiator_Mapping.m');
BLOCK2_SCRIPT = fullfile(masterFolder,'Block_2_Differentiator_Mapping.m');
BLOCK3_SCRIPT = fullfile(masterFolder,'Block_3_Counter_Mapping.m');
BLOCK4_SCRIPT = fullfile(masterFolder,'Block_4_Single_Neuron_NN_Mapping.m');
BLOCK5_SCRIPT = fullfile(masterFolder, ...
    'Block_5_Remove_Cell_Membrane_Reactions.m');
BLOCK6_SCRIPT = fullfile(masterFolder, ...
    'Block_6_Remove_Internal_CrossTalk.m');
BLOCK6B_SCRIPT = fullfile(masterFolder,'Block_6B_Require_EC_Annotations.m');
BLOCK7_SCRIPT = fullfile(masterFolder, ...
    'Block_7_Print_Selected_Motif.m');
CROSSTALK_HELPER = fullfile(masterFolder,'passesInternalCrossTalk.m');
CHECKPOINT_HELPER = fullfile(masterFolder,'saveCompactCheckpoint.m');

blockScripts = {BLOCK0_SCRIPT,BLOCK1_SCRIPT,BLOCK2_SCRIPT, ...
                BLOCK3_SCRIPT,BLOCK4_SCRIPT,BLOCK5_SCRIPT, ...
                BLOCK6_SCRIPT,BLOCK6B_SCRIPT,BLOCK7_SCRIPT};
for k = 1:numel(blockScripts)
    if ~isfile(blockScripts{k})
        error('Master:MissingBlockScript', ...
            'Required block script was not found: %s',blockScripts{k});
    end
end

if FILTER_INTERNAL_CROSSTALK_DURING_SEARCH && ~isfile(CROSSTALK_HELPER)
    error('Master:MissingCrossTalkHelper', ...
        'Required incremental cross-talk helper was not found: %s', ...
        CROSSTALK_HELPER);
end

addpath(masterFolder);
ecHelpers = {'getReactionECAnnotations.m','requireECAnnotatedRows.m', ...
    'requireECAnnotatedAdapters.m'};
for k = 1:numel(ecHelpers)
    if ~isfile(fullfile(masterFolder,ecHelpers{k}))
        error('Master:MissingECHelper','Required EC helper is missing: %s.',ecHelpers{k});
    end
end
if exist('saveCompactCheckpoint','file') ~= 2
    error('Master:MissingCheckpointHelper', ...
        'Required compact checkpoint helper was not found: %s', ...
        CHECKPOINT_HELPER);
end

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

%% ================================================================
%% Common sparse-model caches used by Blocks 1-4
%% ================================================================

nProducers = full(sum(S > 0,2));
substrateSet = cell(nRxns,1);
productSet = cell(nRxns,1);
for r = 1:nRxns
    substrateSet{r} = find(S(:,r) < 0);
    productSet{r} = find(S(:,r) > 0);
end

%% ================================================================
%% BLOCK 0 - cAMP input adapter
%% ================================================================

if RUN_BLOCK0
    run(BLOCK0_SCRIPT);
else
    block0File = fullfile(fileparts(BLOCK0_SCRIPT), ...
        'Block_0_Input_Adapter.mat');
    assert(isfile(block0File), ...
        'Master:MissingBlock0Results', ...
        'Block 0 result file not found: %s',block0File);
    load(block0File);
end

% Apply the same gate to resumed adapter checkpoints without renumbering paths.
InputAdapterPaths = requireECAnnotatedAdapters(model,InputAdapterPaths,hasEC,BLOCK0_REQUIRE_EC);
InputAdapterTable = struct2table(InputAdapterPaths);
Block0AnchorMetIdx = unique(double([InputAdapterPaths.anchorMetIdx]).','stable');
Block0AnchorMetIDs = string(model.mets(Block0AnchorMetIdx));

%% ================================================================
%% BLOCK 1 - first differentiator
%% ================================================================

if RUN_BLOCK1
    run(BLOCK1_SCRIPT);
else
    block1File = fullfile(fileparts(BLOCK1_SCRIPT), ...
        'Block_1_differentiator.mat');
    assert(isfile(block1File), ...
        'Master:MissingBlock1Results', ...
        'Block 1 result file not found: %s',block1File);
    [Motifs_Diff_Block_1,results_diff] = localLoadMotifCheckpoint( ...
        block1File,'Motifs_Diff_Block_1','results_diff');
end

%% ================================================================
%% BLOCK 2 - second differentiator
%% ================================================================

if RUN_BLOCK2
    run(BLOCK2_SCRIPT);
else
    block2File = fullfile(fileparts(BLOCK2_SCRIPT), ...
        'Block_2_differentiator.mat');
    assert(isfile(block2File), ...
        'Master:MissingBlock2Results', ...
        'Block 2 result file not found: %s',block2File);
    [Motifs_Diff_Block_2,results_diff2] = localLoadMotifCheckpoint( ...
        block2File,'Motifs_Diff_Block_2','results_diff2');
end

%% ================================================================
%% BLOCK 3 - counter
%% ================================================================

if RUN_BLOCK3
    run(BLOCK3_SCRIPT);
else
    block3File = fullfile(fileparts(BLOCK3_SCRIPT), ...
        'Block_3_counter.mat');
    assert(isfile(block3File), ...
        'Master:MissingBlock3Results', ...
        'Block 3 result file not found: %s',block3File);
    [Motifs_Counter_Block_3,results_counter] = localLoadMotifCheckpoint( ...
        block3File,'Motifs_Counter_Block_3','results_counter');
end

if ~exist('results_counter','var') || isempty(results_counter)
    error('Master:NoCounterMotifs', ...
        ['Block 3 found no counter motif. Inspect ', ...
         'Counter_Search_Diagnostics_Block_3 before running Block 4.']);
end

%% ================================================================
%% BLOCK 4 - single-neuron NN
%% ================================================================

if RUN_BLOCK4
    run(BLOCK4_SCRIPT);
else
    block4File = fullfile(fileparts(BLOCK4_SCRIPT), ...
        'Block_4_single_neuron_NN.mat');
    assert(isfile(block4File), ...
        'Master:MissingBlock4Results', ...
        'Block 4 result file not found: %s',block4File);
    [Motifs_NN_Block_4,results_nn] = localLoadMotifCheckpoint( ...
        block4File,'Motifs_NN_Block_4','results_nn');
end

%% ================================================================
%% BLOCK 5 - remove motifs containing cell-membrane reactions
%% ================================================================
%
% A complete row is rejected when any canonical reaction in either
% differentiator, the counter, or the single-neuron NN:
%   1) explicitly involves cell envelope/plasma membrane, or
%   2) directly connects cytoplasm/cytosol and extracellular space.
%
% Internal mitochondrial, ER, Golgi and vacuolar membranes are retained.

if RUN_BLOCK5
    model = localEnsureModelCompartments(model);
    run(BLOCK5_SCRIPT);
else
    block5File = fullfile(fileparts(BLOCK5_SCRIPT), ...
        'Block_5_no_cell_membrane.mat');
    assert(isfile(block5File), ...
        'Master:MissingBlock5Results', ...
        'Block 5 result file not found: %s',block5File);
    [Motifs_NN_Block_5,results_nn_no_membrane] = localLoadMotifCheckpoint( ...
        block5File,'Motifs_NN_Block_5','results_nn_no_membrane');
    fullJoinRows = results_nn_no_membrane;
    nF = numel(fullJoinRows);
end

%% ================================================================
%% BLOCK 6 - remove motifs containing internal cross-talk
%% ================================================================

if RUN_BLOCK6
    run(BLOCK6_SCRIPT);
else
    block6File = fullfile(fileparts(BLOCK6_SCRIPT), ...
        'Block_6_no_internal_crosstalk.mat');
    assert(isfile(block6File), ...
        'Master:MissingBlock6Results', ...
        'Block 6 result file not found: %s',block6File);
    [Motifs_NN_Block_6,results_nn_no_crosstalk] = localLoadMotifCheckpoint( ...
        block6File,'Motifs_NN_Block_6','results_nn_no_crosstalk');
    fullJoinRows = results_nn_no_crosstalk;
    nF = numel(fullJoinRows);
end

%% ================================================================
%% BLOCK 6B - final complete-EC check for current or resumed checkpoints
%% ================================================================

nMotifsAfterCrossTalk = numel(results_nn_no_crosstalk);
run(BLOCK6B_SCRIPT);

%% ================================================================
%% BLOCK 7 - print one selected motif
%% ================================================================

if RUN_BLOCK7 && nF > 0
    run(BLOCK7_SCRIPT);
elseif RUN_BLOCK7
    fprintf('\nNo EC-annotated motif survived; Block 7 was skipped.\n');
end

%% ================================================================
%% Final master summary
%% ================================================================

fprintf('\n============================================================\n');
fprintf('MASTER PIPELINE COMPLETE\n');
fprintf('============================================================\n');
fprintf('Model                         : %s\n',modelFilePath);
fprintf('Block 0 adapter paths          : %d\n', ...
    localCountRows(InputAdapterPaths));
fprintf('Block 1 complete motifs        : %d\n', ...
    localCountRows(results_diff));
fprintf('Block 2 complete motifs        : %d\n', ...
    localCountRows(results_diff2));
fprintf('Block 3 complete motifs        : %d\n', ...
    localCountRows(results_counter));
fprintf('Block 4 complete motifs        : %d\n', ...
    localCountRows(results_nn));
fprintf('After cell-membrane filter     : %d\n', ...
    localCountRows(results_nn_no_membrane));
fprintf('After internal cross-talk filter: %d\n', ...
    nMotifsAfterCrossTalk);
fprintf('After complete EC check        : %d\n',nF);
if RUN_BLOCK7 && nF > 0
    fprintf('Printed motif                  : %d\n',MOTIF_TO_PRINT);
    fprintf('Requested print source         : %s\n',string(MOTIF_SOURCE));
end
fprintf('Final table variable           : Motifs_NN_Block_6\n');
fprintf('Final struct-array variable    : results_nn_no_crosstalk\n');
fprintf('Compatibility struct variable : fullJoinRows\n');
fprintf('============================================================\n');

%% Local helper
function [motifTable,rows] = localLoadMotifCheckpoint( ...
        filePath,tableName,legacyStructName)
% Accept compact and old checkpoints, preferring the one table. Selective
% loading avoids bringing duplicate structs/search caches into memory or
% accidentally overriding current search settings with saved settings.
variables = whos('-file',filePath);
names = {variables.name};
if ismember(tableName,names)
    saved = load(filePath,tableName);
    motifTable = saved.(tableName);
    if ~istable(motifTable)
        error('Master:InvalidCheckpointTable', ...
            '%s in %s must be a table.',tableName,filePath);
    end
    rows = table2struct(motifTable);
elseif ismember(legacyStructName,names)
    saved = load(filePath,legacyStructName);
    rows = saved.(legacyStructName);
    if ~isstruct(rows)
        error('Master:InvalidCheckpointStruct', ...
            '%s in %s must be a struct array.',legacyStructName,filePath);
    end
    motifTable = struct2table(rows(:));
else
    error('Master:InvalidCheckpoint', ...
        '%s contains neither %s nor %s.', ...
        filePath,tableName,legacyStructName);
end
rows = rows(:);
fprintf('Loaded %d motif rows from:\n%s\n',height(motifTable),filePath);
end

function n = localCountRows(value)
if istable(value)
    n = height(value);
else
    n = numel(value);
end
end

function model = localEnsureModelCompartments(model)
% Reconstruct model.metComps from metabolite ID suffixes when an SBML
% importer provides model.comps/compNames but omits the per-metabolite map.

nMetsLocal = size(model.S,1);

if isfield(model,'metComps') && numel(model.metComps) == nMetsLocal
    model.metComps = double(model.metComps(:));
else
    metCompSymbols = strings(nMetsLocal,1);

    for m = 1:nMetsLocal
        metID = char(string(model.mets(m)));
        token = regexp(metID,'\[([^\]]+)\]$','tokens','once');

        if isempty(token)
            token = regexp(metID,'_([^_]+)$','tokens','once');
        end

        if isempty(token)
            error('Master:UnknownMetCompartment', ...
                'Could not determine the compartment of metabolite %s.', ...
                metID);
        end

        metCompSymbols(m) = lower(strtrim(string(token{1})));
    end

    if isfield(model,'comps') && ~isempty(model.comps)
        rawComps = model.comps;
        if ischar(rawComps)
            rawComps = cellstr(rawComps);
        end
        modelCompSymbols = lower(strtrim(string(rawComps(:))));
    else
        modelCompSymbols = unique(metCompSymbols,'stable');
        model.comps = cellstr(modelCompSymbols);
    end

    [compartmentFound,metCompIdx] = ismember( ...
        metCompSymbols,modelCompSymbols);

    if any(~compartmentFound)
        unknownSymbols = unique(metCompSymbols(~compartmentFound));
        error('Master:UnknownCompartmentSymbols', ...
            'Unknown metabolite compartment symbols: %s', ...
            strjoin(unknownSymbols,', '));
    end

    model.metComps = double(metCompIdx(:));
    fprintf('Reconstructed model.metComps for %d metabolites.\n', ...
        nMetsLocal);
end

% Block 5 indexes compNames using model.metComps, so compNames must follow
% the same ordering as model.comps.
if ~isfield(model,'comps') || isempty(model.comps)
    error('Master:MissingModelCompartments', ...
        ['model.comps is missing. The master cannot align the reconstructed ', ...
         'model.metComps indices with compartment names.']);
end

rawComps = model.comps;
if ischar(rawComps)
    rawComps = cellstr(rawComps);
end
modelCompSymbols = lower(strtrim(string(rawComps(:))));

validCompNames = false;
if isfield(model,'compNames') && ~isempty(model.compNames)
    rawCompNames = model.compNames;
    if ischar(rawCompNames)
        rawCompNames = cellstr(rawCompNames);
    end
    if numel(rawCompNames) == numel(modelCompSymbols)
        model.compNames = cellstr(string(rawCompNames(:)));
        validCompNames = true;
    end
end

if ~validCompNames
    compNames = strings(numel(modelCompSymbols),1);
    for c = 1:numel(modelCompSymbols)
        switch modelCompSymbols(c)
            case "c"
                compNames(c) = "cytoplasm";
            case "e"
                compNames(c) = "extracellular";
            case "m"
                compNames(c) = "mitochondrion";
            case "n"
                compNames(c) = "nucleus";
            case "r"
                compNames(c) = "endoplasmic reticulum";
            case "g"
                compNames(c) = "golgi";
            case "x"
                compNames(c) = "peroxisome";
            case "p"
                compNames(c) = "peroxisome";
            case "l"
                compNames(c) = "lysosome";
            case "v"
                compNames(c) = "vacuole";
            case "ce"
                compNames(c) = "cell envelope";
            case "pm"
                compNames(c) = "plasma membrane";
            otherwise
                compNames(c) = modelCompSymbols(c);
        end
    end
    model.compNames = cellstr(compNames);
end
end
