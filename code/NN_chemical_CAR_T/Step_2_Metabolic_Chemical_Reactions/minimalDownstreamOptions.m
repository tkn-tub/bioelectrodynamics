%% Shared defaults. Existing explicit master values take priority.
if ~exist('REGULATION_MODE','var'), REGULATION_MODE = "design"; end
if ~exist('REGULATION_FILE','var'), REGULATION_FILE = ''; end
if ~exist('REQUIRE_GPR','var'), REQUIRE_GPR = true; end
if ~exist('BLOCK0_REQUIRE_EC','var'), BLOCK0_REQUIRE_EC = false; end
if ~exist('MAX_HITS_COUNTER','var'), MAX_HITS_COUNTER = 1e5; end
if ~exist('MAX_HITS_NN','var'), MAX_HITS_NN = 1e5; end
if ~exist('MAX_COUNTER_POOL','var'), MAX_COUNTER_POOL = 2000; end
if ~exist('MAX_NN_POOL','var'), MAX_NN_POOL = 2000; end
if ~exist('MAX_DOWNSTREAM_CANDIDATE_CHECKS','var'), MAX_DOWNSTREAM_CANDIDATE_CHECKS = 2e6; end
if ~exist('MAX_DOWNSTREAM_JOIN_CHECKS','var'), MAX_DOWNSTREAM_JOIN_CHECKS = 5e6; end
if ~exist('MAX_OPTIONS_PER_SPECIES','var'), MAX_OPTIONS_PER_SPECIES = 50; end
if ~exist('MAX_NN_NODE_OPTIONS','var'), MAX_NN_NODE_OPTIONS = 2000; end
if ~exist('NN_NODE_OPTIONS_PER_ROLE','var'), NN_NODE_OPTIONS_PER_ROLE = 40; end
if ~exist('ALLOW_CURRENCY_COUNTER_CS','var'), ALLOW_CURRENCY_COUNTER_CS = false; end
if ~exist('FILTER_INTERNAL_CROSSTALK_DURING_SEARCH','var'), FILTER_INTERNAL_CROSSTALK_DURING_SEARCH = true; end
if ~exist('REQUIRE_DISTINCT_ADAPTER_RXNS','var'), REQUIRE_DISTINCT_ADAPTER_RXNS = true; end
if ~exist('SHOW_PROGRESS','var'), SHOW_PROGRESS = true; end
if ~exist('MAP_COUNTER_CS_LEAK','var'), MAP_COUNTER_CS_LEAK = false; end
if ~isscalar(MAP_COUNTER_CS_LEAK) || ~ismember(MAP_COUNTER_CS_LEAK,[0 1])
    error('MinimalOptions:BufferMode','MAP_COUNTER_CS_LEAK must be true or false.');
end
downstreamOptions = struct('requireGPR',REQUIRE_GPR,'requireAdapterEC',BLOCK0_REQUIRE_EC, ...
    'regulationMode',string(REGULATION_MODE),'regulationFile',REGULATION_FILE, ...
    'maxHits',MAX_HITS_COUNTER,'poolCap',MAX_COUNTER_POOL, ...
    'maxChecks',MAX_DOWNSTREAM_CANDIDATE_CHECKS,'maxJoinChecks',MAX_DOWNSTREAM_JOIN_CHECKS, ...
    'optionsCap',MAX_OPTIONS_PER_SPECIES,'nodeCap',MAX_NN_NODE_OPTIONS, ...
    'nodesPerRole',NN_NODE_OPTIONS_PER_ROLE,'allowCurrencyCs',ALLOW_CURRENCY_COUNTER_CS, ...
    'distinctAdapter',REQUIRE_DISTINCT_ADAPTER_RXNS,'mapCsLeak',logical(MAP_COUNTER_CS_LEAK), ...
    'crossTalk',FILTER_INTERNAL_CROSSTALK_DURING_SEARCH,'showProgress',SHOW_PROGRESS);
