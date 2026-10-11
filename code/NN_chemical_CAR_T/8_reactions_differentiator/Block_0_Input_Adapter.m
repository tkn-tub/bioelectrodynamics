%% ================================================================
%% BLOCK 0 - cAMP INPUT ADAPTER
%% ================================================================
%
% Purpose
% -------
% Find directed intracellular reaction paths that start at camp[c], cross
% only explicitly permitted currencies, and continue through non-currency
% metabolites to a bounded depth. Every branch-capable non-currency
% metabolite is exported as a possible m_A for Block 1.
%
% Default topology:
%
%   camp[c] --PDE1--> amp[c] -- ... --> m_A
%
% amp[c] is allowed only as an adapter bridge.  It remains excluded as a
% computational signal in Blocks 1-4.
%
% Required workspace variable
% ---------------------------
%   model       COBRA model containing S, mets and rxns
%
% Optional workspace variables
% ----------------------------
%   isBadMet               logical vector made by the master file
%   REQUIRE_GPR            default: true
%   BLOCK0_REQUIRE_EC      default: false (passive adapter)
%   BLOCK0_INPUT_MET_ID    default: "camp[c]"
%   BLOCK0_BRIDGE_MET_IDS  default: "amp[c]"
%   BLOCK0_MAX_DEPTH       default: 6 reactions, including PDE1
%   BLOCK0_MAX_PATHS       default: 100000
%   BLOCK0_MAX_SHORTEST_PATHS_PER_MET default: 3
%   BLOCK0_CONTINUE_NONCURRENCY      default: true
%   BLOCK0_REQUIRE_BRANCHABLE_ANCHOR default: true
%   BLOCK0_EXTRA_CURRENCY_BASE_IDS   extra T-cell currency IDs
%   BLOCK0_SAVE_BESIDE_SCRIPT         default: true
%   BLOCK0_OUTPUT_DIR      used only when SAVE_BESIDE_SCRIPT=false
%
% Saved output
% ------------
%   Block_0_Input_Adapter.mat
%       InputAdapterPaths       complete paths, one structure per path
%       InputAdapterTable       printable full-path table
%       Block0AnchorMetIdx      unique candidate m_A model indices
%       Block0AnchorMetIDs      unique candidate m_A identifiers
%       Block0ReachableTable    every reached metabolite and its depth
%       Block0Parameters        settings used for the search

%% ================================================================
%% Configuration
%% ================================================================

if ~exist('BLOCK0_REQUIRE_EC','var'), BLOCK0_REQUIRE_EC = false; end

if ~exist('BLOCK0_INPUT_MET_ID','var') || isempty(BLOCK0_INPUT_MET_ID)
    BLOCK0_INPUT_MET_ID = "camp[c]";
end

if ~exist('BLOCK0_BRIDGE_MET_IDS','var') || isempty(BLOCK0_BRIDGE_MET_IDS)
    % Keep this list deliberately short.  Do not add h[c], h2o[c] or pi[c]:
    % expanding through those ubiquitous metabolites causes a graph explosion.
    BLOCK0_BRIDGE_MET_IDS = "amp[c]";
end

if ~exist('BLOCK0_MAX_DEPTH','var') || isempty(BLOCK0_MAX_DEPTH)
    BLOCK0_MAX_DEPTH = 6;
end

if ~exist('BLOCK0_MAX_PATHS','var') || isempty(BLOCK0_MAX_PATHS)
    BLOCK0_MAX_PATHS = 100000;
end

if ~exist('BLOCK0_MAX_SHORTEST_PATHS_PER_MET','var') || ...
        isempty(BLOCK0_MAX_SHORTEST_PATHS_PER_MET)
    BLOCK0_MAX_SHORTEST_PATHS_PER_MET = 3;
end

if ~exist('BLOCK0_MAX_PRINT_PATHS','var') || isempty(BLOCK0_MAX_PRINT_PATHS)
    BLOCK0_MAX_PRINT_PATHS = 50;
end

if ~exist('BLOCK0_CONTINUE_NONCURRENCY','var') || ...
        isempty(BLOCK0_CONTINUE_NONCURRENCY)
    BLOCK0_CONTINUE_NONCURRENCY = true;
end

% A Block 1 anchor must at least have two eligible consuming reactions and
% two possible non-currency branch products. This is necessary, although not
% sufficient, for the complete differentiator topology.
if ~exist('BLOCK0_REQUIRE_BRANCHABLE_ANCHOR','var') || ...
        isempty(BLOCK0_REQUIRE_BRANCHABLE_ANCHOR)
    BLOCK0_REQUIRE_BRANCHABLE_ANCHOR = true;
end

if ~exist('BLOCK0_EXTRA_CURRENCY_BASE_IDS','var') || ...
        isempty(BLOCK0_EXTRA_CURRENCY_BASE_IDS)
    % GDP/GTP are energy-coupling currencies. dpCoA is a cofactor precursor.
    BLOCK0_EXTRA_CURRENCY_BASE_IDS = ["gtp","gdp","dpcoa"];
end

if ~exist('REQUIRE_GPR','var') || isempty(REQUIRE_GPR)
    REQUIRE_GPR = true;
end

if ~exist('BLOCK0_SAVE_BESIDE_SCRIPT','var') || ...
        isempty(BLOCK0_SAVE_BESIDE_SCRIPT)
    BLOCK0_SAVE_BESIDE_SCRIPT = true;
end

block0_scriptFolder = fileparts(mfilename('fullpath'));
if isempty(block0_scriptFolder)
    block0_scriptFolder = pwd;
end

if BLOCK0_SAVE_BESIDE_SCRIPT
    % Deliberately overwrite any stale output-directory variable left by a
    % previous run. This is the safest behavior for master files/HPC jobs.
    BLOCK0_OUTPUT_DIR = block0_scriptFolder;
elseif ~exist('BLOCK0_OUTPUT_DIR','var') || isempty(BLOCK0_OUTPUT_DIR)
    BLOCK0_OUTPUT_DIR = pwd;
end

BLOCK0_INPUT_MET_ID   = string(BLOCK0_INPUT_MET_ID);
BLOCK0_BRIDGE_MET_IDS = string(BLOCK0_BRIDGE_MET_IDS(:));
BLOCK0_EXTRA_CURRENCY_BASE_IDS = ...
    lower(string(BLOCK0_EXTRA_CURRENCY_BASE_IDS(:)));
BLOCK0_OUTPUT_DIR     = char(string(BLOCK0_OUTPUT_DIR));

if ~isscalar(BLOCK0_MAX_DEPTH) || BLOCK0_MAX_DEPTH < 1 || ...
        BLOCK0_MAX_DEPTH ~= floor(BLOCK0_MAX_DEPTH)
    error('BLOCK0_MAX_DEPTH must be a positive integer.');
end

if ~isscalar(BLOCK0_MAX_PATHS) || BLOCK0_MAX_PATHS < 1 || ...
        BLOCK0_MAX_PATHS ~= floor(BLOCK0_MAX_PATHS)
    error('BLOCK0_MAX_PATHS must be a positive integer.');
end

if ~isscalar(BLOCK0_MAX_SHORTEST_PATHS_PER_MET) || ...
        BLOCK0_MAX_SHORTEST_PATHS_PER_MET < 1 || ...
        BLOCK0_MAX_SHORTEST_PATHS_PER_MET ~= ...
        floor(BLOCK0_MAX_SHORTEST_PATHS_PER_MET)
    error('BLOCK0_MAX_SHORTEST_PATHS_PER_MET must be a positive integer.');
end

%% ================================================================
%% Validate the model
%% ================================================================

requiredFields = {'S','mets','rxns'};
for block0_f = 1:numel(requiredFields)
    if ~isfield(model,requiredFields{block0_f})
        error('model.%s is missing.',requiredFields{block0_f});
    end
end

block0_S     = model.S;
block0_nMets = size(block0_S,1);
block0_nRxns = size(block0_S,2);

if numel(model.mets) ~= block0_nMets || numel(model.rxns) ~= block0_nRxns
    error('The dimensions of model.S do not agree with model.mets/model.rxns.');
end

block0_metIDs = string(model.mets(:));
block0_rxnIDs = string(model.rxns(:));

block0_inputIdx = find(strcmpi(block0_metIDs,BLOCK0_INPUT_MET_ID));
if isempty(block0_inputIdx)
    error('Input metabolite %s was not found in model.mets.',BLOCK0_INPUT_MET_ID);
elseif numel(block0_inputIdx) > 1
    error('Input metabolite %s occurs more than once in model.mets.',BLOCK0_INPUT_MET_ID);
end

%% ================================================================
%% Metabolite names and reaction names
%% ================================================================

if isfield(model,'metNames') && numel(model.metNames) == block0_nMets
    block0_metNames = string(model.metNames(:));
else
    block0_metNames = block0_metIDs;
end

if isfield(model,'rxnNames') && numel(model.rxnNames) == block0_nRxns
    block0_rxnNames = string(model.rxnNames(:));
else
    block0_rxnNames = block0_rxnIDs;
end

%% ================================================================
%% Currency mask
%% ================================================================

block0_currencyBaseIDs = [ ...
    "h2o","atp","adp","amp","pi","ppi","h", ...
    "nad","nadh","nadp","nadph","co2","o2", ...
    "nh3","nh4","coa","fad","fadh","fadh2"];

if exist('isBadMet','var') && numel(isBadMet) == block0_nMets
    block0_isCurrency = logical(isBadMet(:));
else
    warning(['isBadMet was not found (or has the wrong length). ', ...
             'Block 0 is constructing its own currency mask.']);

    block0_isCurrency = false(block0_nMets,1);
end

% Convert amp[c] or amp_c to the compartment-free identifier amp. This is
% needed even when the master file supplied isBadMet, because Block 0 adds a
% small T-cell-specific extension without modifying the master's mask.
block0_metBaseIDs = lower(regexprep(block0_metIDs,'\[[^\]]+\]$',''));
block0_metBaseIDs = regexprep(block0_metBaseIDs,'_([a-z0-9]+)$','');

block0_isCurrency = block0_isCurrency | ...
    ismember(block0_metBaseIDs,block0_currencyBaseIDs);
block0_isCurrency = block0_isCurrency | ...
    ismember(block0_metBaseIDs,BLOCK0_EXTRA_CURRENCY_BASE_IDS);

% Currency metabolites that may be crossed inside the adapter.  This does
% not change isBadMet and therefore does not relax later blocks.
block0_isPermittedBridge = ...
    ismember(lower(block0_metIDs),lower(BLOCK0_BRIDGE_MET_IDS));

block0_missingBridgeIDs = BLOCK0_BRIDGE_MET_IDS(~ismember( ...
    lower(BLOCK0_BRIDGE_MET_IDS),lower(block0_metIDs)));

if ~isempty(block0_missingBridgeIDs)
    warning('These requested bridge metabolites were not found: %s', ...
        strjoin(block0_missingBridgeIDs,', '));
end

%% ================================================================
%% Robust GPR mask
%% ================================================================

block0_hasGPR = true(block0_nRxns,1);

if REQUIRE_GPR
    if isfield(model,'grRules') && numel(model.grRules) == block0_nRxns
        block0_hasGPR = strlength(strtrim(string(model.grRules(:)))) > 0;

    elseif isfield(model,'rules') && numel(model.rules) == block0_nRxns
        block0_hasGPR = strlength(strtrim(string(model.rules(:)))) > 0;

    elseif isfield(model,'rxnGeneMat') && ...
            size(model.rxnGeneMat,1) == block0_nRxns
        block0_hasGPR = full(any(model.rxnGeneMat ~= 0,2));

    else
        warning(['REQUIRE_GPR=true, but the model has no usable grRules, ', ...
                 'rules or rxnGeneMat field. All reactions will be allowed.']);
        block0_hasGPR = true(block0_nRxns,1);
    end
end

%% ================================================================
%% Direction bounds
%% ================================================================

[block0_hasEC,block0_ECCodes] = getReactionECAnnotations(model);
fprintf('Block 0 EC-annotated reactions: %d / %d\n', ...
    sum(block0_hasEC),block0_nRxns);
fprintf('Block 0 requires EC annotations: %d\n',BLOCK0_REQUIRE_EC);

if isfield(model,'lb') && numel(model.lb) == block0_nRxns
    block0_lb = full(model.lb(:));
else
    block0_lb = zeros(block0_nRxns,1);
end

if isfield(model,'ub') && numel(model.ub) == block0_nRxns
    block0_ub = full(model.ub(:));
else
    block0_ub = ones(block0_nRxns,1);
end

%% ================================================================
%% Detect extracellular/cell-envelope metabolites and reactions
%% ================================================================

block0_isCytMet      = false(block0_nMets,1);
block0_isExtraMet    = false(block0_nMets,1);
block0_isEnvelopeMet = false(block0_nMets,1);

if isfield(model,'compNames') && isfield(model,'metComps') && ...
        numel(model.metComps) == block0_nMets

    block0_compNames = lower(strtrim(string(model.compNames(:))));
    block0_metComps  = double(model.metComps(:));

    block0_idxCyt = find(block0_compNames == "cytoplasm" | ...
                         block0_compNames == "cytosol");
    block0_idxExtra = find(block0_compNames == "extracellular" | ...
                           block0_compNames == "extracellular space");
    block0_idxEnvelope = find(block0_compNames == "cell envelope" | ...
                              block0_compNames == "plasma membrane");

    block0_isCytMet      = ismember(block0_metComps,block0_idxCyt);
    block0_isExtraMet    = ismember(block0_metComps,block0_idxExtra);
    block0_isEnvelopeMet = ismember(block0_metComps,block0_idxEnvelope);
else
    % Fallback for models whose metabolite IDs use compartment suffixes.
    block0_lowerMetIDs   = lower(block0_metIDs);
    block0_isCytMet      = endsWith(block0_lowerMetIDs,"[c]");
    block0_isExtraMet    = endsWith(block0_lowerMetIDs,"[e]");
    block0_isEnvelopeMet = endsWith(block0_lowerMetIDs,"[ce]") | ...
                           endsWith(block0_lowerMetIDs,"[pm]");
end

block0_forbiddenRxn = false(block0_nRxns,1);

for block0_r = 1:block0_nRxns
    if BLOCK0_REQUIRE_EC && ~block0_hasEC(block0_r)
        block0_forbiddenRxn(block0_r) = true;
        continue;
    end
    block0_participants = find(block0_S(:,block0_r) ~= 0);

    if isempty(block0_participants)
        block0_forbiddenRxn(block0_r) = true;
        continue;
    end

    block0_usesEnvelope = any(block0_isEnvelopeMet(block0_participants));
    block0_crossesCellBoundary = ...
        any(block0_isCytMet(block0_participants)) && ...
        any(block0_isExtraMet(block0_participants));

    % Also remove conventional exchange, demand and sink reactions.
    block0_isBoundaryID = ~isempty(regexpi(char(block0_rxnIDs(block0_r)), ...
        '^(EX_|DM_|SK_|SINK_)','once'));

    block0_forbiddenRxn(block0_r) = ...
        block0_usesEnvelope || block0_crossesCellBoundary || block0_isBoundaryID;
end

%% ================================================================
%% Cheap Block 1 anchor prefilter
%% ================================================================
%
% This does not attempt the seven-reaction search. It only prevents Block 1
% from receiving metabolites that cannot possibly form two initial branches.

block0_anchorEligible = ~block0_isCurrency & ...
    ~block0_isExtraMet & ~block0_isEnvelopeMet;
block0_anchorConsumerCount = zeros(block0_nMets,1);
block0_anchorProductCount = zeros(block0_nMets,1);

if BLOCK0_REQUIRE_BRANCHABLE_ANCHOR
    for block0_m = 1:block0_nMets
        if ~block0_anchorEligible(block0_m), continue; end

        block0_consumers = find(block0_S(block0_m,:) < 0 & ...
            block0_ub.' > 0).';
        block0_consumers = block0_consumers( ...
            ~block0_forbiddenRxn(block0_consumers));
        % The computational branches of Block 1 still require EC codes.
        block0_consumers = block0_consumers(block0_hasEC(block0_consumers));
        if REQUIRE_GPR
            block0_consumers = block0_consumers( ...
                block0_hasGPR(block0_consumers));
        end

        block0_anchorConsumerCount(block0_m) = numel(block0_consumers);
        if numel(block0_consumers) < 2
            block0_anchorEligible(block0_m) = false;
            continue;
        end

        block0_possibleProducts = find(full(any( ...
            block0_S(:,block0_consumers) > 0,2)));
        block0_possibleProducts = block0_possibleProducts( ...
            ~block0_isCurrency(block0_possibleProducts) & ...
            ~block0_isExtraMet(block0_possibleProducts) & ...
            ~block0_isEnvelopeMet(block0_possibleProducts) & ...
            block0_possibleProducts ~= block0_m);

        block0_anchorProductCount(block0_m) = ...
            numel(unique(block0_possibleProducts));
        if block0_anchorProductCount(block0_m) < 2
            block0_anchorEligible(block0_m) = false;
        end
    end
end

%% ================================================================
%% Breadth-first search
%% ================================================================
%
% Permitted currency bridges and non-currency metabolites can be expanded.
% A metabolite is expanded only on its first (therefore shortest) BFS path.
% Eligible non-currency products are also recorded as candidate m_A anchors.

block0_stateTemplate = struct( ...
    'currentMetIdx',[], ...
    'metPath',[], ...
    'rxnPath',[], ...
    'directionPath',[]);

block0_states = block0_stateTemplate;
block0_states.currentMetIdx = block0_inputIdx;
block0_states.metPath = block0_inputIdx;
block0_states.rxnPath = zeros(1,0);
block0_states.directionPath = zeros(1,0);

% Retain only a bounded number of shortest expansion paths per metabolite.
block0_bestDepth = inf(block0_nMets,1);
block0_bestDepth(block0_inputIdx) = 0;
block0_pathsAtBestDepth = zeros(block0_nMets,1);
block0_pathsAtBestDepth(block0_inputIdx) = 1;

block0_pathTemplate = struct( ...
    'pathIndex',[], ...
    'inputMetIdx',[], ...
    'inputMetID',"", ...
    'anchorMetIdx',[], ...
    'anchorMetID',"", ...
    'anchorMetName',"", ...
    'anchorConsumerCount',[], ...
    'anchorProductCount',[], ...
    'nReactions',[], ...
    'metaboliteIndices',[], ...
    'reactionIndices',[], ...
    'reactionDirections',[], ...
    'metabolitePath',"", ...
    'reactionPath',"", ...
    'reactionNamePath',"");

InputAdapterPaths = repmat(block0_pathTemplate,0,1);
block0_stopBecauseMaxPaths = false;
block0_examinedDirectedRxns = 0;

fprintf('\n============================================================\n');
fprintf('BLOCK 0 - cAMP INPUT ADAPTER\n');
fprintf('============================================================\n');
fprintf('Input metabolite       : %s (index %d)\n', ...
    block0_metIDs(block0_inputIdx),block0_inputIdx);
fprintf('Permitted bridge mets  : %s\n', ...
    strjoin(BLOCK0_BRIDGE_MET_IDS,', '));
fprintf('Maximum reaction depth : %d\n',BLOCK0_MAX_DEPTH);
fprintf('Shortest paths/met     : %d\n', ...
    BLOCK0_MAX_SHORTEST_PATHS_PER_MET);
fprintf('Continue non-currency  : %d\n',BLOCK0_CONTINUE_NONCURRENCY);
fprintf('Branchable anchors only: %d\n',BLOCK0_REQUIRE_BRANCHABLE_ANCHOR);
fprintf('GPR required           : %d\n',REQUIRE_GPR);
fprintf('Cell-boundary reactions: excluded\n\n');

for block0_depth = 1:BLOCK0_MAX_DEPTH

    if isempty(block0_states)
        fprintf('No traversable states remain at depth %d; search is complete.\n', ...
            block0_depth);
        break;
    end

    block0_nextStates = repmat(block0_stateTemplate,0,1);
    block0_nStatesThisDepth = numel(block0_states);
    block0_lastPercent = -1;

    for block0_s = 1:block0_nStatesThisDepth

        block0_levelPercent = floor(100*block0_s/block0_nStatesThisDepth);
        if block0_levelPercent >= block0_lastPercent + 10 || ...
                block0_s == block0_nStatesThisDepth
            block0_overallPercent = 100*((block0_depth-1) + ...
                block0_s/block0_nStatesThisDepth)/BLOCK0_MAX_DEPTH;
            fprintf(['Block 0 progress: %5.1f%% | depth %d/%d | ', ...
                     'states %d/%d | paths %d\n'], ...
                block0_overallPercent,block0_depth,BLOCK0_MAX_DEPTH, ...
                block0_s,block0_nStatesThisDepth,numel(InputAdapterPaths));
            block0_lastPercent = block0_levelPercent;
        end

        block0_currentMet = block0_states(block0_s).currentMetIdx;

        % Forward use: the current metabolite is on the substrate side.
        block0_forwardRxns = find( ...
            block0_S(block0_currentMet,:) < 0 & block0_ub.' > 0);

        % Reverse use: the current metabolite is on the product side and
        % the reaction has a negative lower bound.
        block0_reverseRxns = find( ...
            block0_S(block0_currentMet,:) > 0 & block0_lb.' < 0);

        block0_candidateRxns = [block0_forwardRxns(:);block0_reverseRxns(:)];
        block0_candidateDirs = [ones(numel(block0_forwardRxns),1); ...
                               -ones(numel(block0_reverseRxns),1)];

        for block0_q = 1:numel(block0_candidateRxns)

            block0_r   = block0_candidateRxns(block0_q);
            block0_dir = block0_candidateDirs(block0_q);

            if block0_forbiddenRxn(block0_r) || ~block0_hasGPR(block0_r) || ...
                    (BLOCK0_REQUIRE_EC && ~block0_hasEC(block0_r))
                continue;
            end

            % Do not reuse a reaction inside one adapter path.
            if any(block0_states(block0_s).rxnPath == block0_r)
                continue;
            end

            block0_examinedDirectedRxns = block0_examinedDirectedRxns + 1;

            if block0_dir > 0
                block0_products = find(block0_S(:,block0_r) > 0);
            else
                block0_products = find(block0_S(:,block0_r) < 0);
            end

            for block0_p = 1:numel(block0_products)

                block0_nextMet = block0_products(block0_p);

                % Prevent cycles and exclude extracellular/envelope products.
                if any(block0_states(block0_s).metPath == block0_nextMet) || ...
                        block0_isExtraMet(block0_nextMet) || ...
                        block0_isEnvelopeMet(block0_nextMet)
                    continue;
                end

                block0_newMetPath = ...
                    [block0_states(block0_s).metPath,block0_nextMet];
                block0_newRxnPath = ...
                    [block0_states(block0_s).rxnPath,block0_r];
                block0_newDirPath = ...
                    [block0_states(block0_s).directionPath,block0_dir];

                % Unlisted currencies are neither anchors nor traversal
                % nodes. amp[c] remains traversable through the explicit
                % permitted-bridge list.
                block0_isNonCurrency = ...
                    ~block0_isCurrency(block0_nextMet);
                block0_canTraverse = block0_isNonCurrency || ...
                    block0_isPermittedBridge(block0_nextMet);

                if ~block0_canTraverse
                    continue;
                end

                % Breadth-first order guarantees that the first visit is a
                % shortest path. Ignore equal/longer rediscoveries.
                if block0_depth > block0_bestDepth(block0_nextMet)
                    continue;
                end
                if block0_depth < block0_bestDepth(block0_nextMet)
                    block0_bestDepth(block0_nextMet) = block0_depth;
                    block0_pathsAtBestDepth(block0_nextMet) = 0;
                end
                if block0_pathsAtBestDepth(block0_nextMet) >= ...
                        BLOCK0_MAX_SHORTEST_PATHS_PER_MET
                    continue;
                end
                block0_pathsAtBestDepth(block0_nextMet) = ...
                    block0_pathsAtBestDepth(block0_nextMet) + 1;

                if block0_isNonCurrency && ...
                        block0_anchorEligible(block0_nextMet)
                    % Record this branch-capable non-currency endpoint as m_A.
                    block0_newPath = block0_pathTemplate;
                    block0_newPath.pathIndex = numel(InputAdapterPaths) + 1;
                    block0_newPath.inputMetIdx = block0_inputIdx;
                    block0_newPath.inputMetID = block0_metIDs(block0_inputIdx);
                    block0_newPath.anchorMetIdx = block0_nextMet;
                    block0_newPath.anchorMetID = block0_metIDs(block0_nextMet);
                    block0_newPath.anchorMetName = block0_metNames(block0_nextMet);
                    block0_newPath.anchorConsumerCount = ...
                        block0_anchorConsumerCount(block0_nextMet);
                    block0_newPath.anchorProductCount = ...
                        block0_anchorProductCount(block0_nextMet);
                    block0_newPath.nReactions = numel(block0_newRxnPath);
                    block0_newPath.metaboliteIndices = block0_newMetPath;
                    block0_newPath.reactionIndices = block0_newRxnPath;
                    block0_newPath.reactionDirections = block0_newDirPath;
                    block0_newPath.metabolitePath = strjoin( ...
                        block0_metIDs(block0_newMetPath)," -> ");

                    block0_signedRxnIDs = block0_rxnIDs(block0_newRxnPath);
                    block0_signedRxnNames = block0_rxnNames(block0_newRxnPath);
                    block0_reverseMask = block0_newDirPath(:) < 0;
                    block0_signedRxnIDs(block0_reverseMask) = ...
                        block0_signedRxnIDs(block0_reverseMask) + " (reverse)";
                    block0_signedRxnNames(block0_reverseMask) = ...
                        block0_signedRxnNames(block0_reverseMask) + " (reverse)";

                    block0_newPath.reactionPath = ...
                        strjoin(block0_signedRxnIDs," -> ");
                    block0_newPath.reactionNamePath = ...
                        strjoin(block0_signedRxnNames," -> ");

                    InputAdapterPaths(end+1,1) = block0_newPath; %#ok<SAGROW>

                    if numel(InputAdapterPaths) >= BLOCK0_MAX_PATHS
                        block0_stopBecauseMaxPaths = true;
                        break;
                    end
                end

                % Unlike the original version, a non-currency endpoint can
                % also become an internal adapter node at the next depth.
                block0_expandThisMet = ...
                    block0_isPermittedBridge(block0_nextMet) || ...
                    (BLOCK0_CONTINUE_NONCURRENCY && block0_isNonCurrency);

                if block0_expandThisMet && ...
                        block0_depth < BLOCK0_MAX_DEPTH
                    block0_newState = block0_stateTemplate;
                    block0_newState.currentMetIdx = block0_nextMet;
                    block0_newState.metPath = block0_newMetPath;
                    block0_newState.rxnPath = block0_newRxnPath;
                    block0_newState.directionPath = block0_newDirPath;
                    block0_nextStates(end+1,1) = block0_newState; %#ok<SAGROW>
                end
            end

            if block0_stopBecauseMaxPaths
                break;
            end
        end

        if block0_stopBecauseMaxPaths
            break;
        end
    end

    if block0_stopBecauseMaxPaths
        break;
    end

    block0_states = block0_nextStates;
end

%% ================================================================
%% Remove exact duplicate paths and construct outputs
%% ================================================================

fprintf('Block 0 progress: 100.0%% | search complete | raw paths %d\n', ...
    numel(InputAdapterPaths));

Block0ReachableMetIdx = find(isfinite(block0_bestDepth) & ...
    (1:block0_nMets).' ~= block0_inputIdx & ...
    ~block0_isExtraMet & ~block0_isEnvelopeMet);
Block0ReachableTable = table( ...
    Block0ReachableMetIdx, ...
    block0_metIDs(Block0ReachableMetIdx), ...
    block0_metNames(Block0ReachableMetIdx), ...
    block0_bestDepth(Block0ReachableMetIdx), ...
    block0_isCurrency(Block0ReachableMetIdx), ...
    block0_anchorEligible(Block0ReachableMetIdx), ...
    block0_anchorConsumerCount(Block0ReachableMetIdx), ...
    block0_anchorProductCount(Block0ReachableMetIdx), ...
    'VariableNames',{ ...
        'metaboliteIndex','metaboliteId','metaboliteName','shortestDepth', ...
        'isCurrency','isBranchableAnchor','nEligibleConsumers', ...
        'nPossibleBranchProducts'});

if isempty(InputAdapterPaths)
    InputAdapterTable = struct2table(InputAdapterPaths);
    Block0AnchorMetIdx = zeros(0,1);
    Block0AnchorMetIDs = strings(0,1);
else
    block0_pathKeys = strings(numel(InputAdapterPaths),1);
    for block0_k = 1:numel(InputAdapterPaths)
        block0_pathKeys(block0_k) = sprintf('%d,', ...
            [InputAdapterPaths(block0_k).reactionIndices, ...
             InputAdapterPaths(block0_k).reactionDirections, ...
             InputAdapterPaths(block0_k).anchorMetIdx]);
    end

    [~,block0_uniquePathIdx] = unique(block0_pathKeys,'stable');
    InputAdapterPaths = InputAdapterPaths(block0_uniquePathIdx);

    for block0_k = 1:numel(InputAdapterPaths)
        InputAdapterPaths(block0_k).pathIndex = block0_k;
    end

    InputAdapterTable = struct2table(InputAdapterPaths);

    block0_allAnchorIdx = [InputAdapterPaths.anchorMetIdx].';
    Block0AnchorMetIdx = unique(block0_allAnchorIdx,'stable');
    Block0AnchorMetIDs = block0_metIDs(Block0AnchorMetIdx);
end

Block0Parameters = struct;
Block0Parameters.inputMetID = BLOCK0_INPUT_MET_ID;
Block0Parameters.inputMetIdx = block0_inputIdx;
Block0Parameters.bridgeMetIDs = BLOCK0_BRIDGE_MET_IDS;
Block0Parameters.extraCurrencyBaseIDs = BLOCK0_EXTRA_CURRENCY_BASE_IDS;
Block0Parameters.maxDepth = BLOCK0_MAX_DEPTH;
Block0Parameters.maxPaths = BLOCK0_MAX_PATHS;
Block0Parameters.maxShortestPathsPerMetabolite = ...
    BLOCK0_MAX_SHORTEST_PATHS_PER_MET;
Block0Parameters.continueNonCurrency = BLOCK0_CONTINUE_NONCURRENCY;
Block0Parameters.requireBranchableAnchor = ...
    BLOCK0_REQUIRE_BRANCHABLE_ANCHOR;
Block0Parameters.requireGPR = REQUIRE_GPR;
Block0Parameters.saveBesideScript = BLOCK0_SAVE_BESIDE_SCRIPT;
Block0Parameters.requireCompleteEC = logical(BLOCK0_REQUIRE_EC);
Block0Parameters.outputDirectory = string(BLOCK0_OUTPUT_DIR);
Block0Parameters.currencyBridgeDoesNotModifyIsBadMet = true;
Block0Parameters.cellBoundaryReactionsExcluded = true;
Block0Parameters.examinedDirectedReactions = block0_examinedDirectedRxns;
Block0Parameters.reachedMetabolites = sum(isfinite(block0_bestDepth));
Block0Parameters.stoppedAtMaxPaths = block0_stopBecauseMaxPaths;

%% ================================================================
%% Summary and save
%% ================================================================

fprintf('\n############################################################\n');
fprintf('BLOCK 0 INPUT-ADAPTER RESULTS\n');
fprintf('############################################################\n');
fprintf('Complete adapter paths : %d\n',numel(InputAdapterPaths));
fprintf('Unique anchor mets m_A : %d\n',numel(Block0AnchorMetIdx));
fprintf('Metabolites reached     : %d\n',sum(isfinite(block0_bestDepth)));
fprintf('Directed reactions read: %d\n',block0_examinedDirectedRxns);

if block0_stopBecauseMaxPaths
    warning(['BLOCK0_MAX_PATHS=%d was reached. The saved results are ', ...
             'partial; increase the limit only after inspecting them.'], ...
             BLOCK0_MAX_PATHS);
end

if isempty(InputAdapterPaths)
    warning(['No branch-capable intracellular non-currency anchor was ', ...
             'reached from %s ', ...
             'within %d reactions using bridge set {%s}.'], ...
             BLOCK0_INPUT_MET_ID,BLOCK0_MAX_DEPTH, ...
             strjoin(BLOCK0_BRIDGE_MET_IDS,', '));
else
    fprintf('\nCandidate anchors:\n');
    for block0_k = 1:numel(Block0AnchorMetIdx)
        fprintf('  %6d  %-24s  %s\n',Block0AnchorMetIdx(block0_k), ...
            Block0AnchorMetIDs(block0_k), ...
            block0_metNames(Block0AnchorMetIdx(block0_k)));
    end

    block0_nPrintPaths = min(numel(InputAdapterPaths),BLOCK0_MAX_PRINT_PATHS);
    fprintf('\nFirst %d complete paths:\n',block0_nPrintPaths);
    for block0_k = 1:block0_nPrintPaths
        fprintf('  Path %d: %s\n',block0_k, ...
            InputAdapterPaths(block0_k).metabolitePath);
        fprintf('          %s\n',InputAdapterPaths(block0_k).reactionPath);
    end
    if numel(InputAdapterPaths) > block0_nPrintPaths
        fprintf('  ... %d additional paths are stored in InputAdapterTable.\n', ...
            numel(InputAdapterPaths)-block0_nPrintPaths);
    end
end

if ~isfolder(BLOCK0_OUTPUT_DIR)
    [block0_mkdirOK,block0_mkdirMessage] = mkdir(BLOCK0_OUTPUT_DIR);
    if ~block0_mkdirOK
        error('Block0:CannotCreateOutputDirectory', ...
            'Cannot create output directory %s: %s', ...
            BLOCK0_OUTPUT_DIR,block0_mkdirMessage);
    end
end

block0_saveFile = fullfile(BLOCK0_OUTPUT_DIR,'Block_0_Input_Adapter.mat');
Block0OutputFile = string(block0_saveFile);

fprintf('\nSaving Block 0 results to:\n  %s\n',block0_saveFile);

try
    save(block0_saveFile, ...
        'InputAdapterPaths', ...
        'InputAdapterTable', ...
        'Block0AnchorMetIdx', ...
        'Block0AnchorMetIDs', ...
        'Block0ReachableMetIdx', ...
        'Block0ReachableTable', ...
        'Block0Parameters', ...
        'Block0OutputFile', ...
        '-v7.3');
catch block0_saveError
    error('Block0:SaveFailed', ...
        'Could not save %s: %s',block0_saveFile,block0_saveError.message);
end

if ~isfile(block0_saveFile)
    error('Block0:SaveVerificationFailed', ...
        'MATLAB returned from save, but the MAT file is not present: %s', ...
        block0_saveFile);
end

block0_savedFileInfo = dir(block0_saveFile);

fprintf('Verified MAT file: %.2f MB\n', ...
    block0_savedFileInfo.bytes/(1024^2));
fprintf('You can load it with:\n  load(''%s'')\n',block0_saveFile);
fprintf('############################################################\n');
