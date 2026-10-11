FIVE-REACTION DIFFERENTIATORS: OCTOBER 6 KINETICS DESIGN
=====================================================

Install and rerun
-----------------

Extract ALL files into your active MATLAB code folder, replacing files with
the same names. Helpers are required: copying only the two Block scripts is
not sufficient. Your existing model, isBadMet and cAMP/Block 0 setup are used.
No absolute user-machine paths or COBRA initialization are added.

With model and isBadMet loaded and InputAdapterPaths available, run:

  Block_1_Differentiator_Mapping
  Block_2_Differentiator_Mapping

If InputAdapterPaths is absent, Block 1 loads it from the local
Block_0_Input_Adapter.mat. Existing valid Block 0 results may be reused.
Old Block 1 and Block 2 results have different chemistry and must be rerun.
The scripts overwrite Block_1_differentiator.mat and Block_2_differentiator.mat
in the current folder, unless BLOCK1_OUTPUT_DIR/BLOCK2_OUTPUT_DIR are set.

Alternatively run:
  Run_Five_Reaction_Differentiators
This invokes both scripts, skips Block 2 if Block 1 has zero matches, and
prints the first connected differentiator motif locally. It preserves your
explicit master settings. No full network search for Blocks 3-6 is launched.

Exact chemical roles
--------------------
Block 1 has the separate fast pool from the latest kinetics discussion:

  r_S1    m_in -> m_M1
  r_S2    m_M1 -> m_d1
  r_F     m_in -> m_F1_delay
  r_I     m_F1_delay + m_d1 -> products
  r_R     m_F1_delay -> products

Output aliases: m_d1=m_Y1=m_d. The independent r_R removes the FAST pool,
not m_d1. The coupled r_I consumes both the fast pool and the output.

Block 2 starts from that m_d1 output:

  r_S1_2  m_d1 -> m_M2
  r_S2_2  m_M2 -> m_d2
  r_F2    m_d1 -> m_d2
  r_I2    m_M2 + m_d2 -> products
  r_R2    m_d2 -> products

Here the slow conversion r_S2_2 DOES produce m_d2. The fast reaction r_F2
also produces m_d2: both are intended producers inside the algorithm.
The earlier October 5 mapping, where r_S2_2 removed M2 without producing
Y2, is explicitly rejected. m_d2=m_Y2, and m_X2=m_d1.

The first block stores a 17-column table (five directed reaction roles, four
physical species, adapter provenance and output aliases). The second stores
a 33-column table containing the full ten-reaction association. Shared IDs,
names, equations, EC and GPR descriptions are stored once as dictionaries,
without model.S or duplicated result struct arrays in the MAT checkpoint.
Compatibility struct arrays results_diff/results_diff2 are still created in
the workspace for inspection.

Search rules
------------
Each block requires five distinct reaction IDs and distinct species within
that block. Additional co-substrates/products are allowed if they do not add
edges involving the protected signal species. Matching uses signed net
stoichiometry in model.S and feasible directions from model.lb/model.ub.
Stoichiometric magnitudes and kinetic parameters are not fitted by this
search. Pure exchange/demand reactions with no reaction products do not
implement the mapped removal roles.

All differentiator roles require annotated EC codes BEFORE hit counting.
REQUIRE_GPR defaults to true. BLOCK0_REQUIRE_EC defaults to false: the
input adapter need not be enzyme-controlled, as previously requested.
Currency metabolites may participate as co-substrates/co-products, but
new signal species are chosen from the noncurrency pool. cAMP anchoring
remains in Block 0, independent of identifier conventions in yeast/CAR-T.

Block 2 searches once for each DISTINCT first output, then joins candidates
to the corresponding original Block 1 rows. External network reactions may
produce or consume those metabolites. With
  FILTER_INTERNAL_CROSSTALK_DURING_SEARCH=true
only the selected algorithm reactions must have the intended signal edges.
Their intended producers include BOTH r_S2_2 and r_F2 for m_d2.

The recommended distinct-stage settings, defaulted by Block 2 if absent, are:
  REQUIRE_STAGE_DISJOINT_RXNS=true
  REQUIRE_STAGE_DISJOINT_METS=true
Only the shared first output/second input is exempt. Explicit values already
in your master take precedence. The downstream full-motif filter requires
distinct physical signal species and reaction IDs.

Budgets, progress and HPC
------------------------
Existing settings are preserved. Default joined-hit limits are:
  MAX_HITS=1e5              (Block 1)
  MAX_HITS_BLOCK2=1e5       (Block 2)
Default per-input motif limits are MAX_STAGE1_PER_INPUT=10000 and
MAX_STAGE2_PER_INPUT=10000. A legacy MAX_STAGE2_PER_PAIR is accepted as
the latter value when MAX_STAGE2_PER_INPUT is absent.

MAX_DIFF_CANDIDATE_CHECKS (default 2e5) limits each distinct-input search.
MAX_DIFF_JOIN_CHECKS (default 2e6) includes accepted AND rejected joins.
SHOW_PROGRESS defaults to true.
Percentage is the fraction of distinct inputs completed; within a long
input search the code prints checks and candidate counts. Saving follows
the search, so search completion is not a serialization progress estimate.

For an already allocated MATLAB/HPC pool set:
  USE_PARALLEL=true
  USE_PARALLEL_BLOCK2=true
The code keeps at most NumWorkers inputs in flight, shares the graph via
parallel.pool.Constant and cancels pending tasks after the joined-hit cap.
AUTO_START_PARPOOL defaults to false. Set it true only if pool creation is
part of your MATLAB/HPC setup; PARPOOL_WORKERS controls its requested size.
No pool/toolbox available -> serial fallback. Parallel completion order may
change which candidates are retained at a finite global cap. A single
distinct input is handled by one worker and gains little from this scheme.

Inspect saved SearchInfo and diagnostics before interpreting zero results:
a capped search with zero matches is not proof of absence. Raising budgets
is different from changing Block 0's adapter depth.

Optional settings script
------------------------
Five_Reaction_Differentiator_Settings supplies example settings and disables
the older Block 0 branchability prefilter. Run it BEFORE Block 0 only if you
want to reset your settings to that example. It need not be called to rerun
Blocks 1-2 from an existing adapter checkpoint.

Printing and continuation
-------------------------
  Print_Five_Reaction_Differentiators
loads only the compact Block 2 MAT file and writes Differentiator_Motif_1.txt
locally, including the adapter path, species and all ten reactions. Set
POLARITY_MOTIF_TO_PRINT before invoking it to choose another available row.
No model load or COBRA is needed for this focused report.

Downstream compatibility files are included so all ten reaction roles are
checked, saved and printed consistently. After finding new Block 2 motifs:
  Resume_Five_Reaction_Pipeline
recomputes Blocks 3-6 and prints Block 7. In a master that still calls
Resume_Polarity_Pipeline, replace that call with Resume_Five_Reaction_Pipeline.
Continue_Minimal_Pipeline now also accepts the current five-reaction format.

Current differentiator TopologyVersion: five_reaction_differentiators_v2
Current downstream TopologyVersion: five_reaction_regulated_full_v2
Old checkpoints are rejected, never silently relabelled. An old print-model
cache is usable only if its model ID ordering still agrees with the results.

Downstream counter roles remain r_Cp,r_Cg,r_Cr with only m_d2 as its regulator.
Cs supply/removal support is external by default (the COMPLETE net buffer
term). MAP_COUNTER_CS_LEAK=true optionally maps an additional basal removal
when only positive Cs supply is provided externally. The full mapped motif
has 21 reactions by default (10 differentiator + 3 counter + 8 NN), or 22
with CsLeak. No CsSupply is mapped.

REGULATION_MODE="design" produces structural candidates with regulatory
requirements pending validation. "evidence" requires the supplied individual
links. Structural matching alone cannot establish a slow/fast timescale,
derivative polarity, adaptation, stability, or biological regulation.
KineticsValidated remains false in every checkpoint.

Validation
----------
Independent numerical checks passed for both five-reaction templates,
old-topology rejection, each required role's eligibility, reversed reaction
directions, the complete ten-reaction cascade, all column/alias positions,
the 21/22-reaction downstream schemas and report label alignment. Basic
MATLAB block/bracket balance was also checked.

MATLAB/Octave are not installed in the development environment. The included
native tests were not executed here. In MATLAB run:
  test_five_reaction_differentiators
  test_five_reaction_pipeline
Optional allocated-pool hit-limit test:
  test_five_reaction_differentiators(true)
These use synthetic models and temporary folders, require no COBRA, and
cover the actual search, MAT roundtrips, continuation, filters and printing.
Your real model and its candidate availability were not tested here.
