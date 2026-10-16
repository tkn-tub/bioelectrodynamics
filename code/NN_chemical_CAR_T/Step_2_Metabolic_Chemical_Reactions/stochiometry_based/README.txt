METABOLIC CIRCUIT: SIGNAL LABELS AND cAMP INPUT RANKING

DISTINCT INPUT SEARCH AND SMALLER RESULT FILES

ONE_MOTIF_PER_INPUT = true returns one circuit for each distinct m_in.
The search stops exploring an input once it finds that input's first valid
circuit. Input groups, including compartment aliases, are assigned to
disjoint parallel workers. MAX_MOTIFS = 20 therefore means up to 20 different
inputs, subject to MAX_SECONDS and model feasibility. Input identity follows
METABOLITE_IDENTITY: the default 'chemical' groups compartment copies together;
'model_id' treats individual species IDs as distinct inputs.

Set ONE_MOTIF_PER_INPUT = false to enumerate multiple circuits per input.
With a fixed m_in, distinct-input mode returns at most one circuit.
The report's exhaustive status refers to coverage of input representatives in
this mode, not enumeration of every circuit for each input. The representative
is the first circuit found for that input; it is not optimized for its cAMP
connector. The cAMP ranking evaluates these returned representatives afterward.

New search and ranking MAT files use compressed '-v7', retaining the complete
model and all result variables. This avoids the HDF5 container overhead of
'-v7.3' for cell-heavy model data. Actual file-size reduction depends on the
dataset. MAT v7 has a 2 GiB per-variable limit.

To compact an existing result file without repeating the search:
  Compact_Circuit_Search_Results('/full/path/Circuit_Search_Results.mat')
This writes Circuit_Search_Results_compact.mat and keeps the original file.
The function reloads the compact file and verifies every variable's value
before finalizing it. It prints original size, compact size and percentage
saved. The compact file works with Rank_Existing_Circuit_Results.

Put these .m files together in your existing code tree. The TNM1055 model and
COBRA Toolbox are external dependencies and are not included in this archive.

RANK YOUR EXISTING SOLUTIONS WITHOUT REPEATING THE SEARCH

Call the provided MATLAB function with your saved results file:
  Rank_Existing_Circuit_Results('/full/path/Circuit_Search_Results.mat')

With no arguments it selects the latest results under Circuit_Search_Results
beside the function, or a results file already on the MATLAB path:
  Rank_Existing_Circuit_Results

Optional settings can be supplied as a second argument:
  opts = struct('CAMP_MAX_DEPTH',8);
  [ranking,paths,rankInfo,outputFolder] = Rank_Existing_Circuit_Results('/full/path/Circuit_Search_Results.mat',opts);

The script uses the model and motifs saved with your original search, keeps
their indices, applies the new labels, and writes a CSV table, a text path
report, a relabelled circuit report, and Circuit_CAMP_Ranking.mat.

NEW SEARCHES

Run A_Master_File. It now ranks the returned solutions automatically and saves
campRanking, campPaths and campRankInfo beside the existing results variables.
Both reports print complete connector reactions using the circuit format:
ID, name, GPR, EC, enzyme keys, compartments, direction, bounds, and formulas,
chemical names and metabolite IDs with full stoichiometric coefficients.
Every motif's connection is embedded in Selected_Circuit_Motifs.txt before
its ordinary circuit blocks. The closest finite connection is also printed
in the MATLAB command window (PRINT_CAMP_CHAIN_TO_CONSOLE controls this in
the master). Intermediate connector labels are m_CAMP_1, m_CAMP_2, etc.,
matching r_CAMP_1, r_CAMP_2, etc.; the source is m_cAMP and the target is m_in.
The search still allows all shared non-signal metabolites (scope 'roles').
Existing enzyme exclusivity, exact selected-signal connectivity and direction
checks remain active. Metabolite isolation applies to the selected reactions,
not to every unrelated reaction elsewhere in the metabolic model.

SIGNAL LABELS (stored index order is preserved)

Index   Previous label   Current label   Producing role
1       m_in             m_in            external circuit input
2       m_F              m_Ff            r_Ff
3       m_S              m_Sf            r_Sf
4       m_P1             m_Fd            r_Fd
5       m_d1             m_D1            r_Sd (differentiator 1 output)
6       m_F2             m_Ff_2          r_Ff_2
7       m_S2             m_Sf_2          r_Sf_2
8       m_d2             m_D2            r_Fd_2 (differentiator 2 output)
9       m_N2             m_Sd_2          r_Sd_2
10      m_c              m_c             r_c (counter output)
11      m_I              m_I             buffered counter input
12      m_b              m_b             buffered bias input
13      m_out            m_out           r_Nw and r_Nb (neuron output)

Labels do not change the chemical model's biological metabolite names or IDs.
Legacy FIXED_METABOLITES role labels remain accepted by the circuit search.
Older result files can be relabelled when their reaction roles and topology
match; the saved-results runner checks this before using their indices.

WHAT THE RANKING MEANS

The target is m_in, not an internal signal or m_out. Distance counts the
smallest number of allowed directed substrate-to-product reaction steps
from an identified intracellular cAMP species to that input under the
configured filters. It ranks only the supplied/discovered motifs.

Distance 0: m_in is cAMP itself (the same exact model species).
Distance 1: one connecting enzymatic conversion.
Distance N: N connecting conversions, each listed with its direction.
Distance Inf: no path found with these settings. A finite depth limit does
not establish that longer paths or unmodelled regulatory links do not exist.

Rows with equal distances share a rank and retain their original motif order.
MotifIndex refers to the original motifs array, which is never reordered.
The table and campPaths have the same sorted row order.

CAMP SETTINGS

CAMP_SOURCE_METABOLITES = {} automatically recognizes exact conventional
3',5'-cAMP IDs/names. If automatic detection misses your model annotation, set
exact model IDs, e.g. {'camp[c]'}. The report records the sources actually used.
Explicit source IDs are your designation of cAMP; opaque IDs are not guessed.
Species explicitly annotated as 2',3'-cAMP are not accepted as 3',5'-cAMP.

CAMP_SOURCE_COMPARTMENT = 'c' selects cytosolic sources. Change this to your
model's cytosol code if necessary. '' accepts sources from any compartment.
CAMP_ALLOW_TRANSPORT = false excludes reactions spanning known compartments.
true permits explicitly modelled transports; it never adds implicit jumps.
Paths containing species with unknown compartments are flagged in the report.
CAMP_MAX_DEPTH = 6 limits reaction steps. 0 tests exact cAMP inputs only;
Inf explores the full filtered graph.
CAMP_EXCLUDED_PATH_METABOLITES adds exclusions by exact ID or base ID.

Common cofactors/hubs (ATP, ADP, water, protons, phosphate, etc.) cannot serve
as path intermediates or input targets for this ranking. This exclusion is
only for the cAMP connector ranking and does not reject a discovered circuit.
AMP and adenosine are retained for possible cAMP breakdown pathways.

Connecting reactions cannot reuse selected circuit reactions or share an
EC/GPR resource with the circuit. They cannot touch another selected signal,
including a compartment alias. They may arrive at m_in, but cannot use it as
an intermediate. Other co-substrates are assumed available, as in your model.

BridgeEnzymeConflict flags enzyme reuse between connecting steps. Such a row
remains a shortest structural connection; review it before installing that
connector as additional independently tunable reactions. A longer compatible
connector may exist. The ranking does not search for that alternative.

KINETIC INTERPRETATION

This is a stoichiometric graph projection. It does not establish steady-state
flux feasibility, atom transfer, reaction rate separation, or transmission of
the input waveform. It does not infer regulation from substrate/product edges
and does not add the connecting path to the mapped circuit.

Clayton et al. (2011), DOI 10.4049/jimmunol.1003884, support this anchor:
exosome-associated CD39/CD73 generate extracellular adenosine from ATP/AMP;
adenosine receptor A2A signaling leads to intracellular cAMP production through
adenylyl cyclase. Adenosine is a regulator of this pathway, not a stoichiometric
substrate converted directly into cAMP. The paper measured a cAMP response;
it did not establish preservation of vesicle-count peak timing or tumor-vesicle
specificity. That transfer function remains to be evaluated for the application.

VALIDATION

On your MATLAB machine, run:
  selfTestCircuitSearch
  selfTestCAMPRanking
  selfTestDistinctInputs

Optional parallel uniqueness checks:
  selfTestDistinctInputs(true)

The synthetic tests cover circuit constraints and directed cAMP ranking,
including bounds, hub filtering, signal/enzymatic conflicts, compartment
transport, path reconstruction and stable motif indices. The package received
static checks and code review here; MATLAB/Octave was unavailable, so these
runtime tests have not been executed in this workspace.
