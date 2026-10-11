%% Optional settings: call after your master settings and BEFORE Block 0.
% Keep cAMP and your existing Block 0 adapter search unchanged.
% Its old two-branch endpoint prefilter does not apply to this topology.
BLOCK0_REQUIRE_BRANCHABLE_ANCHOR = false;
BLOCK0_REQUIRE_EC = false; % passive adapter; controlled differentiator reactions require EC
MAX_HITS = 1e5;
MAX_HITS_BLOCK2 = 1e5;
MAX_STAGE1_PER_INPUT = 10000;
MAX_STAGE2_PER_INPUT = 10000;
MAX_DIFF_CANDIDATE_CHECKS = 2e5; % per distinct input; inf for exhaustive search
MAX_DIFF_JOIN_CHECKS = 2e6;      % accepted AND rejected connection attempts
USE_PARALLEL = false;
USE_PARALLEL_BLOCK2 = false;    % true to use your allocated MATLAB/HPC pool
AUTO_START_PARPOOL = false;
REQUIRE_GPR = true;
FILTER_INTERNAL_CROSSTALK_DURING_SEARCH = true;
REQUIRE_DISTINCT_ADAPTER_RXNS = true;
REQUIRE_STAGE_DISJOINT_RXNS = true;
REQUIRE_STAGE_DISJOINT_METS = true;
SHOW_PROGRESS = true;
BLOCK1_OUTPUT_DIR = pwd;
BLOCK2_OUTPUT_DIR = pwd;
