%% Default candidate search: explicit regulation to validate or engineer later.
REGULATION_MODE = "design"; % "evidence" requires all supplied individual regulatory links
REGULATION_FILE = '';       % e.g. 'Regulatory_Interactions.csv'; empty = no evidence
MAX_HITS_COUNTER = 1e5;
MAX_HITS_NN = 1e5;
MAX_COUNTER_POOL = 2000;
MAX_NN_POOL = 2000;
MAX_DOWNSTREAM_CANDIDATE_CHECKS = 2e6;
MAX_OPTIONS_PER_SPECIES = 50;
MAX_NN_NODE_OPTIONS = 2000;
NN_NODE_OPTIONS_PER_ROLE = 40;
ALLOW_CURRENCY_COUNTER_CS = false;
MAP_COUNTER_CS_LEAK = false; % complete Cs buffer term is provided externally
% true: positive external supply, with an additional mapped basal Cs removal.
FILTER_INTERNAL_CROSSTALK_DURING_SEARCH = true;
REQUIRE_DISTINCT_ADAPTER_RXNS = true;
SHOW_PROGRESS = true;
MOTIF_TO_PRINT = 1;
PRINT_SELECTED_MOTIF = true;
FILTER_ADAPTER_MEMBRANE = true;
