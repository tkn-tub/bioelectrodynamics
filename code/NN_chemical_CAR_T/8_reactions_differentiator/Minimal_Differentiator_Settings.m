%% Put this call AFTER the master search settings and BEFORE Block 0.
% Minimal_Differentiator_Settings
BLOCK0_REQUIRE_BRANCHABLE_ANCHOR = false; % discard the obsolete two-branch prefilter
BLOCK0_REQUIRE_EC = false;               % passive adapter exception
MAX_HITS = 1e5;
MAX_HITS_BLOCK2 = 1e5;
MAX_STAGE1_PER_INPUT = 10000;
MAX_STAGE2_PER_INPUT = 10000;
MAX_DIFF_CANDIDATE_CHECKS = 2e5;          % inf for exhaustive per-input enumeration
USE_PARALLEL = false;
USE_PARALLEL_BLOCK2 = false;             % true to use your allocated HPC pool
AUTO_START_PARPOOL = false;
FILTER_INTERNAL_CROSSTALK_DURING_SEARCH = true;
REQUIRE_DISTINCT_ADAPTER_RXNS = true;
REQUIRE_STAGE_DISJOINT_RXNS = true;
REQUIRE_STAGE_DISJOINT_METS = true;       % shared Y1=X2 is always permitted
SHOW_PROGRESS = true;
BLOCK1_OUTPUT_DIR = pwd;
BLOCK2_OUTPUT_DIR = pwd;
% In the old master, insert return DIRECTLY after its Block 2 section.
% RUN_BLOCK3=false loads cached results in that master; it does not skip it.
% Alternatively use A_Master_File_CAR_T_Minimal_Differentiators from this bundle.
