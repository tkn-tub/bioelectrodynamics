function results = A_Master(modelInput, userOptions)
% A_MASTER Connect branches, Anderson canceller and counter.
% Usage:
%   A_Master
%   R = A_Master('model_aux.xml');
%   R = A_Master(model); % import + a seed-dependent autocatalytic module
%   R = A_Master(model, struct('productionMode',"single_reaction"));
%   R = A_Master(model, struct('productionMinimumOnly',true));
%   A_Master('selftest') % no COBRA/model/network needed
%   R = A_Master(model, struct('block0MaxAlternatives',10,'block0MaxDepth',20));
% Block 0 enumerates alternative balanced cAMP-to-m_in modules. It permits
% nonunit stoichiometry and shared INTERNAL intermediates, while excluding
% every downstream noncurrency participant except the intended m_in output.
% Alternative enumeration is bounded; a zero count is not an impossibility proof.
%
% Target architecture:
%   S controls a membrane protein (evidence separate from stoichiometry)
%   that imports a CO-SUBSTRATE X of the production module:
%   X[outside] -(S)-> X[inside]; effective production module:
%   I[inside] + X -> 2 I[inside] + W.
%   S MUST gate the X importer, never the I importer. A gated I import
%   contributes k*S*I[outside], i.e. a SOURCE term in dI/dt; only a
%   reaction in which S-controlled X meets I can realise k*S*I[inside].
%   requireIImport=false (DEFAULT): no I import reaction is searched at
%   all. Single-shot operation assumes I(0) > 0 is already present in the
%   vessel, so eps = 0 and there is no Block A. I pools then come from
%   cancellerI instead of from the transporter list. Set it true to
%   restore the earlier import-paired search.
%   With eps = 0 the two arms must be matched at rest or I decays
%   exponentially while it waits; that is a calibration requirement the
%   structural search does NOT check.
%   Stages 06a/06b search and report the S-gated X importers; stages
%   03/05 keep the earlier I-import pairing for reference.
%   One OR several existing reactions may implement the production module.
%   Actual I consumption and production must sum to 1 I -> 2 I after
%   normalization; non-I internal intermediates must cancel. X -> I alone
%   is not an autocatalyst. Single-reaction matching remains optional.
%   F may gate import of a co-substrate L used in I + L -> products.
%   Exact gross 2 I -> I + W and alternative 2 I -> W are catalogued separately.
%   An alternative pair-loss event removes TWO I; its ODE contribution is -2v.
%   Stoichiometric coefficients do NOT establish reaction orders or k*I^2.
% X and W are vectors of feed substrates and waste products, not one species.
% All gross/core coefficients and reaction flux multipliers are retained.
%
% IMPORTANT: this is a structural candidate search, NOT a proof of dI=k*S*I.
% Import feeds an intracellular autocatalytic core; I[e] and I[c] are distinct
% state variables. No export/recycling loop, closure, or kinetic fit is added.
% GPR means gene association, not S-dependent gating. Positive protein edges
% are reported as hypotheses, NEVER reclassified as opening evidence.
%
% Adapted from the supplied A_Master(20261113-121519).m: model discovery,
% direction checks, GPR decoding, currency IDs and EC evidence tools.
% No delay-line or enzyme-activator search.
% All canceller roles use the same precise intracellular I species index.
% Branch extension: common m_in -> m_Ff -> m_Fd plus r_Fr removal, and
% common m_in -> m_Sf -> m_Sd plus r_Sr removal. Every branch reaction needs
% a complete EC annotation. Fast/slow are DESIGN ROLES, not measured speeds.
% m_Sd is the proposed S gate signal; m_Fd is the proposed F gate signal.
% Unidentified gate interfaces remain provisional, never biochemical proof.
% The attached simulation instead currently uses an I/Q annihilation ODE;
% this search preserves the requested THREE-REACTION Anderson architecture.
% Counter default: a primitive single-copy conversion I + X -> C + W,
% with a distinct noncurrency product C. This ADDS a loss term to dI/dt.
% Formation-only accumulators and counters with a separate C-removal step
% are reported separately. C-removal must not regenerate I.
% Optional counterMode="membrane" retains the earlier I-gated cargo search.
% Step 2 default: the neuron consumes m_c directly in its inside compartment.
% Search 2..neuronChainMaxDepth conversion steps plus an independent output
% removal. All selected reactions and intermediate metabolites are reported.
% Rejected consuming chains are reported as DIAGNOSTICS ONLY, separately
% from isolated neuron candidates; I-feedback checks are never bypassed.
% neuronPlacement="outside" retains the optional adapter/export search.
% Up to neuronAdapterMaxDepth intracellular conversions may first turn m_c
% into an exportable m_exp. Adapter paths retain every co-participant and
% need EC evidence, distinct controlled-role ECs and full isolation from I.
% A zero-step adapter retains the original direct-export possibility.
% Extracellular I is allowed only when it has no selected upstream pool
% overlap; external GSH feeding the existing I importer is rejected.
% Only model-recorded export and extracellular reactions are searched;
% no extracellular species, transporter or enzyme relocation is invented.
% One consuming neuron (m_Nin_1 -> Nf1 -> Nd1, plus Nd1 removal),
% with output m_out = m_Nd_1. A constrained Block 0 starts at
% exact intracellular cAMP and ends at branch input m_in. Every added
% metabolic step requires an EC. Reports distinguish partial from complete.
% Step_2_Results.txt follows the supplied signal-map/reaction-block format.
%
% QUICK START (COBRA Toolbox already initialized):
%   A_Master('selftest');
%   R = A_Master;
% Linux/HPC, with no display or figures:
%   matlab -batch "A_Master('selftest')"
%   matlab -batch "A_Master"
% The no-argument run loads model_aux.xml. Keep it in the current folder or
% under the code directory. The search settings are included in this file.
% Only console, TXT, CSV and MAT outputs are produced; no graphics are used.
% SBML metadata extraction is native MATLAB and also works without a JVM.
% COBRA model loading remains unchanged.
% Protein discovery additionally needs organismTaxID and gene annotations:
%   opt = struct('onlineProteinDiscovery',true);
%   R = A_Master('model_aux.xml',opt);
% Positive protein edges remain hypotheses, not membrane-opening evidence.
%
% DEFAULTS:
%   Forward and reverse orientations where the model bounds permit them.
%   STRUCTURAL SCREEN: EC and transporter GPR requirements are OFF.
%   Missing EC/GPR annotations are reported, never inferred from stoichiometry.
%   Compartments come from SBML/model metadata. model_aux searches u and x
%   separately; the two blocks are never merged or treated as proteins.
%   cellFreeMerge=true (default): ONE reaction vessel, no membranes. Inside
%   each block, organelle copies of a species (cytosol, mitochondrion,
%   peroxisome, ...) become one species; transports with no net chemistry
%   are removed and identical isoenzyme copies become one reaction. The
%   outside compartment and the S-gated I import are kept. Set
%   cellFreeMerge=false to recover the compartmentalised cell model.
%   cancellerI="pyr" (default): the canceller is searched ONLY for this I
%   (BaseID such as "pyr", or an exact ID such as "AUXOSCpyr[u]"; a string
%   array allows several). Set cancellerI="" to search every imported I.
%   counterUseEmbeddediMM904EC=true (default): the counter, neuron and
%   Block 0 stages may use the iMM904 authors' EC annotations (760
%   reactions, exact reaction ID plus matching chemistry). The fast/slow
%   branch stage keeps its own EC rules.
%   counterRequireFormationEC=false, counterRequireLeakEC=false (default):
%   the counter only has to accumulate C linearly with the I it consumes,
%   so its reactions need no EC. A counter reaction WITH an EC still joins
%   the distinct-EC checks; one without an EC is simply left out of them.
%   Step 2 speed: step2ShareBlock0Search=true runs the cAMP -> m_in search
%   once per (branch module, canceller) and only re-checks each counter
%   variant against it (the counter can only forbid more, never allow more).
%   step2MaxSeconds=21600 (6 h) stops Step 2 cleanly and still writes all
%   outputs. Progress is written to Step2_Progress.txt in outputDir.
%   step2SearchBlock0=false (default): Step 2 searches the neuron only; the
%   cAMP -> m_in connector (Block 0) is skipped and parents with a passing
%   neuron are reported as "neuron_found_Block0_not_searched". Set it true
%   to add Block 0. block0RequireEC=false (default): Block 0 reactions need
%   no EC and stay out of the distinct-EC assignment.
%   excludedReactionIDs: lumped pseudo-reactions that no single enzyme
%   catalyses (default: thiazole-phosphate synthesis THZPSN1_SC/2_SC).
%   Original currency exclusions are retained. No EC numbers are excluded.
%   productionMode="multi_reaction" is the default, with up to 10 existing
%   steps for an I-dependent module implementing effective I+X->2I+W.
%   Modules are searched in layers of 1, 2, 3, ... distinct model reactions.
%   Shortest modules for each imported I are written separately. Longer
%   modules remain available for downstream wiring unless productionMinimumOnly
%   is true. Minimum certificates apply only to the simple-cycle search class
%   and require every smaller layer to finish without solver-unresolved cases.
%   productionMode="single_reaction" optionally requires one explicit
%   primitive I+X->2I+W model reaction and applies the direct atom filter.
%   Module atom balance is reported separately. Missing/unsupported formulas
%   remain explicitly unresolved. No reaction or hidden I catalyst is invented.
%   Each I has a 20,000-state budget; reached limits appear in coverage CSV.
%   Protein signals S only when validating supplied gating evidence.
%   Network retrieval is OFF; only explicitly named protein CSVs are reused.
%   Core/removal screening imposes no unique-EC constraint.
%   Removal connection is ON. Unknown F controls remain explicit placeholders.
%   By default, connected roles cannot reuse reactions or share noncurrency
%   species except I. Partial loss components are still written to CSV.
%   Biomass-sized loss reactions (>16 participants) are audited, not used.
%   Branches are ON; all six branch reactions require EC annotations.
%   An embedded, source-linked six-reaction EC seed map is included. It is
%   a partial annotation set, not full yeast-model EC coverage. The map is
%   applied only to explicit SBML BiGG references with matching chemistry.
%   Six distinct EC assignments and no shared known branch genes are required.
%   Counter search is ON in consume_I mode. Formation and any selected
%   removal need ECs distinct from the assigned branch ECs and each other.
%   Formation-only accumulators are retained with removal explicitly absent.
%   All stoichiometric co-substrates, products and cross-talk are retained.
%
% SELECTED OPTIONS (pass a struct; see defaultOptions for the full list):
%   maxCoreReactions=12                    % permit longer production modules
%   productionMode="single_reaction"      % optional direct-reaction filter
%   allowReverse=false                   % restrict to forward orientations
%   requireCoreEC=true, requireTransportGPR=true % annotated candidates only
%   onlineEC=true                        % exact-ID annotation enrichment
%   requireGateEvidence=true             % retain only evidence-linked pairs
%   requireProteinS=false                % also admit nonprotein gating ligands
%   gateEvidenceFile='/absolute/path/Gate_Evidence.tsv'
%   outputDir='/absolute/path/results'
%   connectRemovalTerms=false             % production-only comparison
%   allowPairLossAlternative=false        % demand exact gross 2I -> I
%
% MAIN OUTPUTS: 03_Production_Candidates.csv; Anderson_Production_Report.txt;
% Anderson_Production_Results.mat. Additional CSVs audit ECs, gross vs net
% stoichiometry, importers, cores, search coverage and protein hypotheses.
% DIRECT PRODUCTION: 02_Direct_Production_Filter_Audit.csv;
% 03_Import_Autocatalytic_Pairs.csv; Eq_Block_1_Import_Production_Report.txt.
% MINIMUM PRODUCTION: 02_Production_Depth_Coverage.csv;
% 02_Minimum_Production_Modules.csv; 03_Minimum_Import_Production_Pairs.csv;
% Minimum_Production_Report.txt. The importer is not counted as a core step.
% EXTERNAL NEURON: 14_Neuron_Export_Audit.csv;
% 15_Transported_Neuron_Candidates.csv; 16_Neuron_Compartment_Audit.csv.
% Export and neuron-only matches are separate from complete Block0 circuits.
% DIRECT CONSUMPTION: 15_Consuming_Neuron_Chains.csv records actual model
% paths, including clearly flagged feedback diagnostics. Only
% 15_Isolated_Neuron_Candidates.csv contributes isolated neuron candidates.
% 15_Neuron_Chain_Coverage.csv records separate clean/diagnostic budgets.
% ADAPTER: 14_Counter_Adapter_Audit.csv; 14_Counter_Adapter_Coverage.csv;
% 14_Counter_Adapter_Export_Links.csv. m_A1, m_A2, ... label adapter
% intermediates, m_exp is the final intracellular export cargo, and r_A_N1_XX
% identify the connecting reactions. No I-return relaxation is applied.
% Each pair reports transporter ID, I_out/I_in, R1 ID/direction, gross
% equation, module equation, X substrates, W products, atom/gate status.
% Step_2_Results.txt includes upstream counts and the first empty stage.
% NEW: 05_Linear_Removal_Candidates.csv; 05_Quadratic_Removal_Candidates.csv;
% 05_Removal_Coverage.csv; 05_Removal_Filter_Audit.csv;
% 06_Connected_Anderson_Candidates.csv; Anderson_Three_Reaction_Report.txt.
% BRANCH OUTPUTS: 07_Branch_EC_Audit.csv; 07_Branch_Modules.csv;
% 07_Branch_Reaction_Assignments.csv; 08_Branches_Canceller_Wiring.csv;
% 08_Matched_Control_Interfaces.csv; Anderson_Fast_Slow_Branch_Report.txt.
% Branch_EC_Evidence.tsv can extend annotations using columns
% SourceReactionID, EC, Reference. Unknown annotations stay excluded.
% COUNTER OUTPUTS: 09_Counter_Components.csv; 09_Counter_EC_Audit.csv;
% 10_Canceller_Counter_Connections.csv; 10_Counter_Kinetics_Template.csv;
% 11_Branches_Canceller_Counter_Wiring.csv; Anderson_Counter_Report.txt.
% Counter_EC_Evidence.tsv extends removal annotations; Counter_Gate_Evidence.tsv
% supplies exact I/transporter/direction/cargo/side mappings and references.
% In membrane mode, rows without gate evidence are STRUCTURAL HYPOTHESES.
% In default consume_I mode, the gate is not part of the counter. Read
% 10_Consuming_Counter_Connection_Audit.csv for missing EC/reuse/overlap
% rejections. counterRequireRemoval=true demands a separate C-removal step.
% Gate_Evidence_TEMPLATE.tsv is written if there is no supplied evidence file.
% Fill it only from mechanism-specific evidence. It uses exact reaction ID,
% direction (+1/-1), intracellular cargo ID and model gene ID. Effect must be
% 'opening'; Mechanism: ligand_gating/protein_gating/pore_assembly;
% S_NotConsumed: true; Organism must match organismName. A supplied citation
% is recorded, not automatically verified or treated as a measured rate law.
%
% LIMITS: bounded simple-cycle search; no general branched-network search.
% Currency co-substrates are assumed supplied. Local seed dependence is a
% qualitative AND-substrate reachability check, not a quantitative transient
% or a demonstration that the intact cell has no alternative source of I.
% One imported copy plus one normalized R1 turnover gives net TWO internal
% I molecules from I_out plus X. The intracellular seed I is still required
% for R1 to fire. The paired equation does NOT assert matched fluxes.
% S-gated I import controls entry/seeding, not R1 itself. S-dependent stopping
% requires additional production/X-supply control and kinetic evidence.
% linprog or a configured COBRA LP solver is needed only if a candidate's
% flux weights are underdetermined. Such cases are reported if unsolved.

% EDIT HERE to choose the model for a no-argument Linux/HPC batch run.
if nargin < 1, modelInput = 'model_aux.xml'; end
if nargin < 2, userOptions = struct(); end
if (ischar(modelInput) || isstring(modelInput)) && strcmpi(modelInput,'selftest')
    selfTestAnderson(); results = struct('selfTestsPassed',true); return;
end
cfg = defaultOptions();
keys = fieldnames(userOptions);
for k=1:numel(keys)
    assert(isfield(cfg,keys{k}),'Anderson:Option','Unknown option: %s',keys{k});
    cfg.(keys{k}) = userOptions.(keys{k});
end
validateOptions(cfg);
codeRoot = addCodeParent();
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
fid = fopen(fullfile(cfg.outputDir,'Anderson_Production_Report.txt'),'w');
assert(fid>0,'Anderson:Report','Cannot write output report.');
guard = onCleanup(@() fclose(fid)); %#ok<NASGU>
report(fid,'ANDERSON: S-GATED IMPORT + INTRACELLULAR AUTOCATALYTIC PRODUCTION');
report(fid,'Reader revision: native MATLAB SBML metadata; no Java required, 2026-11-13.');
report(fid,['Search file: ' mfilename('fullpath')]);
report(fid,['Started: ' datestr(now,31)]);
report(fid,'Scope: production plus optional same-I linear/quadratic removal connections.');
report(fid,'Candidates are structural; S gating, kinetics and spatial closure are separate checks.');
report(fid,['Production mode: ' char(cfg.productionMode)]);
report(fid,sprintf('Depth <= %d; path-state budget per imported I = %g.', ...
    cfg.maxCoreReactions,cfg.maxStatesPerI));

modelPath = '';
if isstruct(modelInput)
    model = modelInput;
else
    modelPath = findNamed(codeRoot,char(modelInput));
    assert(~isempty(modelPath),'Anderson:Model', ...
        'Model not found: %s. Put it in the current folder or under the code directory, or pass its full path.',modelInput);
    assert(exist('readCbModel','file')~=0,'Anderson:COBRA', ...
        'COBRA Toolbox readCbModel must be on the MATLAB path.');
    model = readCbModel(modelPath);
end
model = normalizeModel(model);
report(fid,sprintf('Model: %s; %d reactions; %d metabolites.', ...
    modelPath,numel(model.rxns),numel(model.mets)));
[baseMet,comp,layout,nativeMetadata] = readModelLayout(model,modelPath);
model.circuitCompartment=comp;
model.circuitBaseIDs=baseMet;
cfg=resolveModelProfile(cfg,layout,comp);
report(fid,['Model identity: ' char(layout.modelName)]);
report(fid,['Inside compartments (searched separately): ' char(joinText(cfg.insideCompartment,'; '))]);
report(fid,['Outside compartment: ' char(cfg.outsideCompartment)]);
for label=reshape(unique(comp),1,[])
    report(fid,sprintf('  Compartment %s: %d metabolites.',label,nnz(comp==label)));
end
if layout.communityIDs
    report(fid,'AUXOS/PROTS prefixes identify separate model blocks, not protein synthesis.');
    if cfg.cellFreeMerge
        report(fid,'cellFreeMerge=1: organelle copies of each species are merged inside each block (details below).');
    else
        report(fid,'Original organelle labels remain in species IDs; no metabolite states are merged.');
    end
end
report(fid,sprintf('Annotation filters: requireCoreEC=%d; requireTransportGPR=%d.', ...
    cfg.requireCoreEC,cfg.requireTransportGPR));
if ~cfg.requireCoreEC || ~cfg.requireTransportGPR
    report(fid,'STRUCTURAL SCREEN: missing enzyme/protein evidence does not exclude topology candidates.');
    if isempty(cfg.excludedECs)
        report(fid,'No EC exclusions are applied.');
    else
        report(fid,'Configured EC exclusions apply only to resolved EC annotations.');
    end
end
report(fid,sprintf('Reactions with a nonempty gene rule: %d.', ...
    nnz(strlength(strtrim(string(model.grRules)))>0 | strlength(strtrim(string(model.rules)))>0)));
currency = ismember(lower(regexprep(baseMet,'^M_','')),lower(cfg.currencyBaseIDs));
metNames=string(model.mets);
if isfield(model,'metNames') && numel(model.metNames)==numel(model.mets)
    metNames=string(model.metNames(:));
end
metCatalog=table((1:numel(model.mets))',string(model.mets),metNames,baseMet,comp,currency,layout.source, ...
    'VariableNames',{'ModelIndex','MetaboliteID','MetaboliteName','BaseID','Compartment','ExcludedAsSignal','CompartmentSource'});
writetable(metCatalog,fullfile(cfg.outputDir,'00_Metabolite_Catalog.csv'));
% Reuse EC policy from the original code; map the chosen inside compartment
% to its original c-only policy without changing species IDs or stoichiometry.
ecComp = repmat("other",size(comp)); ecComp(ismember(comp,cfg.insideCompartment))="c";
[Rgross,Pgross,grossKnown,grossAudit,sbmlEC,sbmlBiGG,sbmlInfo] = readGrossStoichiometry(model,modelPath,nativeMetadata);
model.circuitSBMLEC = sbmlEC;
model.circuitSourceBiGG = sbmlBiGG;
report(fid,sprintf('SBML metadata: %d reaction nodes; %d matched model reactions; %d explicit BiGG IDs; %d gross equations recovered.', ...
    sbmlInfo.parsedReactionCount,sbmlInfo.matchedReactionCount,sbmlInfo.sourceBiGGCount,nnz(grossKnown)));
if sbmlInfo.error~="",report(fid,['SBML metadata error: ' char(sbmlInfo.error)]);end
ecOpt = struct('online',cfg.onlineEC,'timeout',10,'maxLookups',cfg.maxECLookups, ...
    'maxFailures',3,'mappingFile',cfg.ecMappingFile,'brendaFile',cfg.brendaFile, ...
    'modelOnly',~cfg.requireCoreEC && isempty(cfg.excludedECs) && ~cfg.onlineEC);
[model,ecAudit,ecInfo] = resolveRequiredECAnnotations(model,ecComp,model.ub,codeRoot,ecOpt);
writetable(ecAudit,fullfile(cfg.outputDir,'00_EC_Annotation_Audit.csv'));
writetable(grossAudit,fullfile(cfg.outputDir,'00_Gross_Stoichiometry_Audit.csv'));
report(fid,ecInfo);

% CELL-FREE MODE: merge organelle copies BEFORE every search stage, so that
% all downstream blocks see one pot with one pool per chemical species.
cellFreeAudit=table();cellFreeMetAudit=table();cellFreeInfo=struct();
if cfg.cellFreeMerge
    [model,baseMet,comp,currency,Rgross,Pgross,grossKnown,ecComp,metSource,cellFreeAudit, ...
        cellFreeMetAudit,cellFreeInfo]=mergeCellFreeCompartments(model,baseMet,comp,currency, ...
        Rgross,Pgross,grossKnown,ecComp,layout.source,layout,cfg); %#ok<ASGLU>
    layout.source=metSource;
    writetable(cellFreeAudit,fullfile(cfg.outputDir,'00_CellFree_Reaction_Merge.csv'));
    writetable(cellFreeMetAudit,fullfile(cfg.outputDir,'00_CellFree_Metabolite_Merge.csv'));
    report(fid,'CELL-FREE MODE: one vessel, no membranes; organelle copies merged inside each block.');
    report(fid,sprintf('  Metabolites %d -> %d (%d species had several organelle copies).', ...
        cellFreeInfo.metabolitesBefore,cellFreeInfo.metabolitesAfter,cellFreeInfo.mergedSpecies));
    report(fid,sprintf('  Reactions %d -> %d: removed %d intracellular transports with no net chemistry and %d with a currency-only net change.', ...
        cellFreeInfo.reactionsBefore,cellFreeInfo.reactionsAfter, ...
        cellFreeInfo.removedNoNetChemistry,cellFreeInfo.removedCurrencyOnly));
    report(fid,sprintf('  Organelle-spanning reactions kept with shared species cancelled: %d; identical copies merged: %d; representatives with widened bounds: %d.', ...
        cellFreeInfo.spanningKept,cellFreeInfo.duplicatesMerged,cellFreeInfo.boundsWidened));
    report(fid,'  The outside compartment and the I import transporters are NOT merged: the S-gated import stays in the topology.');
    report(fid,'  Inspect 00_CellFree_Reaction_Merge.csv and 00_CellFree_Metabolite_Merge.csv.');
    metNames=string(model.mets);
    if isfield(model,'metNames') && numel(model.metNames)==numel(model.mets)
        metNames=string(model.metNames(:));
    end
    metCatalog=table((1:numel(model.mets))',string(model.mets),metNames,baseMet,comp,currency,layout.source, ...
        model.circuitMergedMetIDs, ...
        'VariableNames',{'ModelIndex','MetaboliteID','MetaboliteName','BaseID','Compartment', ...
        'ExcludedAsSignal','CompartmentSource','CellFreeMergedIDs'});
    writetable(metCatalog,fullfile(cfg.outputDir,'00_Metabolite_Catalog.csv'));
end

[oriented,transport,poolAudit] = buildPools(model,baseMet,comp,currency,Rgross,Pgross,grossKnown,cfg);
writetable(poolAudit,fullfile(cfg.outputDir,'00_Reaction_Filter_Audit.csv'));
writetable(transport,fullfile(cfg.outputDir,'01_Import_Transporters.csv'));
report(fid,sprintf('Eligible intracellular orientations: %d; import rows: %d.', ...
    numel(oriented.rxn),height(transport)));
for label=reshape(cfg.insideCompartment,1,[])
    report(fid,sprintf('  Import rows into %s: %d.',label, ...
        nnz(comp(transport.InsideIndex)==label)));
end
if isempty(oriented.rxn) || isempty(transport)
    report(fid,'SEARCH INPUT EMPTY: a zero result here cannot establish absence of autocatalysis.');
    report(fid,'Inspect 00_Metabolite_Catalog.csv and 00_Reaction_Filter_Audit.csv before increasing search limits.');
end
cores = emptyCores();statsRows=cell(0,12);directProductionAudit=directProductionAuditTable(cell(0,15));
productionDepthCoverage=productionDepthAuditTable(cell(0,8));
wantedI=lower(strtrim(string(cfg.cancellerI(:))));wantedI=wantedI(wantedI~="");
if cfg.requireIImport || isempty(wantedI)
    if ~cfg.requireIImport
        report(fid,'requireIImport=false but cancellerI is empty: falling back to the imported-I pool list.');
    end
    signals = unique(transport.InsideIndex,'stable');
else
    % No Block A: I(0) is already in the vessel, so an importer is not
    % required and the I pools come from cancellerI directly.
    signals = find(ismember(comp,cfg.insideCompartment) & ~currency);
    report(fid,sprintf(['requireIImport=false: I pools taken from cancellerI, not from the ' ...
        'transporter list (%d intracellular noncurrency pools before filtering).'],numel(signals)));
end
signals=reshape(signals,[],1);
if ~isempty(wantedI)
    keepI=ismember(lower(baseMet(signals)),wantedI)|ismember(lower(string(model.mets(signals))),wantedI);
    report(fid,sprintf('cancellerI=%s: %d of %d imported I pools kept (%s).',joinText(wantedI,','), ...
        nnz(keepI),numel(signals),joinText(string(model.mets(signals(keepI))),'; ')));
    if ~any(keepI)
        report(fid,'cancellerI matched no imported I pool: check the BaseID in 00_Metabolite_Catalog.csv and 01_Import_Transporters.csv.');
    end
    signals=signals(keepI);
end
for q=1:numel(signals)
    seed = signals(q);
    [found,stats,oneAudit,oneDepth] = findCores(seed,oriented,currency,cfg,model);
    directProductionAudit=[directProductionAudit;oneAudit]; %#ok<AGROW>
    productionDepthCoverage=[productionDepthCoverage;oneDepth]; %#ok<AGROW>
    cores = [cores found]; %#ok<AGROW>
    statsRows(end+1,:) = {string(model.mets{seed}),stats.states,stats.closures, ...
        numel(found),stats.capped,stats.unsolved,stats.smallestFound,stats.minimumProven, ...
        stats.minimumLayerComplete,stats.completedDepth,stats.minimumStatus,stats.minimumScope}; %#ok<AGROW>
    report(fid,sprintf('I = %s: %d cores; %d path states; capped=%d; solver-unresolved=%d.', ...
        model.mets{seed},numel(found),stats.states,stats.capped,stats.unsolved));
    report(fid,sprintf('  Smallest found=%g; minimum proven within search class=%d; completed depth=%d; %s.', ...
        stats.smallestFound,stats.minimumProven,stats.completedDepth,stats.minimumStatus));
end
statsTable = cellTable(statsRows, ...
    {'I_ID','PathStates','ClosedPathsTested','CoreCount','StateBudgetReached','SolverUnresolved', ...
    'SmallestFoundReactionCount','MinimumProvenWithinSearchClass','MinimumLayerComplete', ...
    'CompletedReactionDepth','MinimumStatus','MinimumScope'}, ...
    {'string','double','double','double','logical','double','double','logical','logical','double','string','string'});
writetable(statsTable,fullfile(cfg.outputDir,'02_Search_Coverage.csv'));
coreTable = describeCores(cores,oriented,model);
coreTable=tagProductionMinima(coreTable,statsTable,'I_ID','CoreReactionCount');
minimumCoreTable=coreTable(coreTable.IsShortestFoundForI,:);
writetable(coreTable,fullfile(cfg.outputDir,'02_Autocatalytic_Cores.csv'));
writetable(minimumCoreTable,fullfile(cfg.outputDir,'02_Minimum_Production_Modules.csv'));
writetable(productionDepthCoverage,fullfile(cfg.outputDir,'02_Production_Depth_Coverage.csv'));
writetable(directProductionAudit,fullfile(cfg.outputDir,'02_Direct_Production_Filter_Audit.csv'));

% No artificial S assignment. Optional exact reaction/cargo/gene evidence.
gatePath = findNamed(codeRoot,cfg.gateEvidenceFile);
[gates,gateAudit] = readGateEvidence(gatePath,transport,cfg);
writetable(gateAudit,fullfile(cfg.outputDir,'03_Gate_Evidence_Audit.csv'));
if isempty(gatePath)
    writeGateTemplate(fullfile(cfg.outputDir,'Gate_Evidence_TEMPLATE.tsv'));
end
[candidates,details] = assembleCandidates(cores,oriented,transport,gates,model,cfg);
writetable(candidates,fullfile(cfg.outputDir,'03_Production_Candidates.csv'));
productionPairs=describeImportProductionPairs(candidates,details,cores,oriented,transport,model);
productionPairs=tagProductionMinima(productionPairs,statsTable,'I_Inside','ProductionReactionCount');
minimumProductionPairs=productionPairs(productionPairs.IsShortestFoundForI,:);
writetable(productionPairs,fullfile(cfg.outputDir,'03_Import_Autocatalytic_Pairs.csv'));
writetable(minimumProductionPairs,fullfile(cfg.outputDir,'03_Minimum_Import_Production_Pairs.csv'));
writeMinimumProductionReport(fullfile(cfg.outputDir,'Minimum_Production_Report.txt'), ...
    statsTable,productionDepthCoverage,minimumProductionPairs,cfg);
writeImportProductionReport(fullfile(cfg.outputDir,'Eq_Block_1_Import_Production_Report.txt'), ...
    productionPairs,directProductionAudit,cfg);

% Reuse protein-level discovery from the supplied code, only for transporters
% that actually connect to a core. Never query BRENDA enzyme activators here.
usedTransport = unique(details.TransportRow);usedTransport(usedTransport<1)=[];
[proteinMap,proteinEdges,regulators] = discoverRegulators( ...
    transport(usedTransport,:),codeRoot,cfg,fid);
writetable(proteinMap,fullfile(cfg.outputDir,'04_Transporter_Protein_Map.csv'));
writetable(proteinEdges,fullfile(cfg.outputDir,'04_Protein_Edges.csv'));
writetable(regulators,fullfile(cfg.outputDir,'04_Protein_Regulator_Hypotheses.csv'));

% Connect the other two formal terms without inventing kinetic orders or F.
removalSignals=signals;if ~cfg.connectRemovalTerms,removalSignals=[];end
[linearRemoval,quadraticRemoval,removalCoverage,removalAudit]=findRemovalTerms( ...
    removalSignals,oriented,transport,model,comp,currency,cfg);
[connected,connectionStats]=connectAndersonTerms(candidates,details,cores,oriented, ...
    transport,linearRemoval,quadraticRemoval,model,comp,currency,Rgross,Pgross,cfg);
writetable(linearRemoval,fullfile(cfg.outputDir,'05_Linear_Removal_Candidates.csv'));
writetable(quadraticRemoval,fullfile(cfg.outputDir,'05_Quadratic_Removal_Candidates.csv'));
writetable(removalCoverage,fullfile(cfg.outputDir,'05_Removal_Coverage.csv'));
writetable(removalAudit,fullfile(cfg.outputDir,'05_Removal_Filter_Audit.csv'));
writetable(connected,fullfile(cfg.outputDir,'06_Connected_Anderson_Candidates.csv'));
writeThreeReactionReport(fullfile(cfg.outputDir,'Anderson_Three_Reaction_Report.txt'), ...
    candidates,details,cores,oriented,transport,linearRemoval,quadraticRemoval, ...
    connected,connectionStats,model,cfg);

% S gate on the CO-SUBSTRATE importer, not on the I importer. The I import
% found above stays in the architecture as the basal seed term only.
[sGateImports,sGateStats]=findSGatedCoreImports(cores,oriented,transport,model,currency,cfg);
[gatedAnderson,gatedStats]=attachSGates(connected,details,sGateImports,cores,oriented, ...
    transport,linearRemoval,quadraticRemoval,model,currency,Rgross,Pgross,cfg);
writetable(sGateImports,fullfile(cfg.outputDir,'06a_S_Gated_Cosubstrate_Importers.csv'));
writetable(gatedAnderson,fullfile(cfg.outputDir,'06b_Gated_Anderson_Candidates.csv'));
writeSGateReport(fullfile(cfg.outputDir,'Gated_Anderson_Report.txt'), ...
    sGateImports,gatedAnderson,sGateStats,gatedStats,cfg);
report(fid,sprintf(['S gate: %d cores, %d feeds tested, %d importers, %d cores gated; ' ...
    'crossed %d, reuse-rejected %d, crosstalk-rejected %d, %d connected rows had no ' ...
    'S-gated feed -> %d gated candidates.'], ...
    sGateStats.cores,sGateStats.feedsTested,sGateStats.importersFound, ...
    sGateStats.coresWithSGate,gatedStats.tested,gatedStats.reuseRejected, ...
    gatedStats.crossTalkRejected,gatedStats.candidatesWithoutSGate,height(gatedAnderson)));

% Require ECs for the new rails without discarding structural cancellers.
[branchEC,branchECAudit]=resolveBranchEC(model,baseMet,codeRoot,cfg);
writetable(branchECAudit,fullfile(cfg.outputDir,'07_Branch_EC_Audit.csv'));
[branchRails,branchModules,branchAssignments,branchStats]=findBranchModules( ...
    oriented,branchEC,model,comp,currency,cfg);
[branchWiring,branchWiringStats]=wireBranchesToCancellers(branchModules,branchRails, ...
    connected,candidates,details,cores,oriented,transport,linearRemoval,quadraticRemoval, ...
    model,currency,Rgross,Pgross,cfg);
writetable(branchModules,fullfile(cfg.outputDir,'07_Branch_Modules.csv'));
writetable(branchAssignments,fullfile(cfg.outputDir,'07_Branch_Reaction_Assignments.csv'));
writetable(branchWiring,fullfile(cfg.outputDir,'08_Branches_Canceller_Wiring.csv'));
writetable(branchWiring(branchWiring.BothControlInterfacesMatched,:), ...
    fullfile(cfg.outputDir,'08_Matched_Control_Interfaces.csv'));
writeBranchReport(fullfile(cfg.outputDir,'Anderson_Fast_Slow_Branch_Report.txt'), ...
    branchECAudit,branchModules,branchAssignments,branchWiring,branchStats,branchWiringStats,sbmlInfo,cfg);

% Default counter consumes I; optional membrane mode retains the old search.
[counterEC,counterECAudit]=resolveCounterEC(model,baseMet,branchEC,codeRoot,cfg);
if cfg.counterMode=="consume_I"
    [counterComponents,counterParts,counterAudit,counterStats]=findConsumingCounterFormations( ...
        connected,oriented,counterEC,model,baseMet,comp,currency,cfg);
    [counterConnections,counterConnectionParts,counterConnectionAudit,counterConnectionStats]= ...
        connectConsumingCounters(counterComponents,counterParts,connected,details,cores, ...
        oriented,transport,linearRemoval,quadraticRemoval,counterEC,model,baseMet,currency,Rgross,Pgross,cfg);
    [counterWiring,counterWiringStats]=wireConsumingCounters(counterConnections,counterComponents, ...
        counterConnectionParts,branchWiring,branchModules,branchRails,oriented,model,currency,cfg);
    counterGates=table();counterGateAudit=table();
    writetable(counterConnectionAudit,fullfile(cfg.outputDir,'10_Consuming_Counter_Connection_Audit.csv'));
    writeConsumingCounterOutputs(counterComponents,counterConnections,counterWiring, ...
        counterAudit,counterConnectionAudit,counterStats,counterConnectionStats,counterWiringStats, ...
        counterECAudit,model,cfg);
else
[counterComponents,counterParts,counterAudit,counterStats]=findCounterComponents( ...
    oriented,transport,counterEC,model,baseMet,comp,currency,Rgross,Pgross,cfg);
[counterGates,counterGateAudit]=readCounterGateEvidence( ...
    findNamed(codeRoot,cfg.counterGateEvidenceFile),counterComponents,connected,cfg);
[counterConnections,counterConnectionStats]=connectCounters(counterComponents,counterParts, ...
    counterGates,connected,details,cores,oriented,transport,linearRemoval,quadraticRemoval, ...
    model,baseMet,currency,Rgross,Pgross,cfg);
[counterWiring,counterWiringStats]=wireCountersToBranches(counterConnections,counterComponents, ...
    counterParts,branchWiring,branchModules,branchRails,model,currency,cfg);
writetable(counterECAudit,fullfile(cfg.outputDir,'09_Counter_EC_Audit.csv'));
writetable(counterComponents,fullfile(cfg.outputDir,'09_Counter_Components.csv'));
writetable(counterAudit,fullfile(cfg.outputDir,'09_Counter_Filter_Audit.csv'));
writetable(counterConnections,fullfile(cfg.outputDir,'10_Canceller_Counter_Connections.csv'));
writetable(counterConnections(counterConnections.GateEvidenceRows>0,:), ...
    fullfile(cfg.outputDir,'10_Counter_Evidence_Linked_Connections.csv'));
writetable(counterGateAudit,fullfile(cfg.outputDir,'10_Counter_Gate_Evidence_Audit.csv'));
writetable(counterWiring,fullfile(cfg.outputDir,'11_Branches_Canceller_Counter_Wiring.csv'));
writeCounterTemplates(counterComponents,counterConnections,model,cfg);
writeCounterReport(fullfile(cfg.outputDir,'Anderson_Counter_Report.txt'),counterComponents, ...
    counterConnections,counterWiring,counterStats,counterConnectionStats,counterWiringStats,cfg);
counterConnectionParts=cell(0,1);counterConnectionAudit=table();
end

% Extend current consuming-counter parents; preserve existing search stages.
[step2EC,step2ECAudit]=resolveStep2EC(model,baseMet,counterEC,codeRoot,cfg);
circuitContext=struct('model',model,'base',baseMet,'comp',comp,'currency',currency, ...
    'R',Rgross,'P',Pgross,'grossKnown',grossKnown,'O',oriented,'T',transport,'C',connected,'D',details, ...
    'cores',cores,'L',linearRemoval,'Q',quadraticRemoval,'M',branchModules, ...
    'rails',branchRails,'X',counterConnections,'F',counterComponents,'Y',counterWiring, ...
    'XP',{counterConnectionParts},'E',step2EC);
upstreamCoverage=cellTable({ ...
    'Intracellular reaction orientations',numel(oriented.rxn),true,'00_Reaction_Filter_Audit.csv'; ...
    'I import transporters',height(transport),true,'01_Import_Transporters.csv'; ...
    'Imported I pools',numel(signals),true,'02_Search_Coverage.csv'; ...
    'Seed-dependent production modules',numel(cores),true,'02_Autocatalytic_Cores.csv'; ...
    'Import plus production pairs',height(candidates),true,'03_Production_Candidates.csv'; ...
    'Linear I-removal candidates',height(linearRemoval),cfg.connectRemovalTerms,'05_Linear_Removal_Candidates.csv'; ...
    'Pair-loss candidates',height(quadraticRemoval),cfg.connectRemovalTerms,'05_Quadratic_Removal_Candidates.csv'; ...
    'Connected Anderson cancellers',height(connected),cfg.connectRemovalTerms,'06_Connected_Anderson_Candidates.csv'; ...
    'S-gated co-substrate importers',height(sGateImports),true,'06a_S_Gated_Cosubstrate_Importers.csv'; ...
    'Gated Anderson candidates',height(gatedAnderson),cfg.connectRemovalTerms,'06b_Gated_Anderson_Candidates.csv'; ...
    'EC-qualified fast/slow modules',height(branchModules),cfg.connectBranches,'07_Branch_Modules.csv'; ...
    'Branch plus canceller connections',height(branchWiring),cfg.connectBranches,'08_Branches_Canceller_Wiring.csv'; ...
    'Counter components',height(counterComponents),cfg.connectCounter,'09_Counter_Components.csv'; ...
    'Canceller plus counter connections',height(counterConnections),cfg.connectCounter,'10_Canceller_Counter_Connections.csv'; ...
    'Branch-canceller-counter parents',height(counterWiring),cfg.connectCounter,'11_Branches_Canceller_Counter_Wiring.csv'}, ...
    {'Stage','Count','Enabled','DiagnosticFile'},{'string','double','logical','string'});
circuitContext.upstreamCoverage=upstreamCoverage;
circuitContext.productionCoverage=statsTable;
circuitContext.directProductionAudit=directProductionAudit;
[step2,circuitContext]=searchStep2Circuit(circuitContext,cfg);
writetable(upstreamCoverage,fullfile(cfg.outputDir,'12_Upstream_Search_Coverage.csv'));
writetable(step2ECAudit,fullfile(cfg.outputDir,'12_Step2_EC_Annotation_Audit.csv'));
writetable(step2.parents,fullfile(cfg.outputDir,'13_Step2_Parent_Coverage.csv'));
writetable(step2.entryAudit,fullfile(cfg.outputDir,'14_Neuron_Entry_Audit.csv'));
writetable(step2.transportAudit,fullfile(cfg.outputDir,'14_Neuron_Export_Audit.csv'));
writetable(step2.adapterAudit,fullfile(cfg.outputDir,'14_Counter_Adapter_Audit.csv'));
writetable(step2.adapterCoverage,fullfile(cfg.outputDir,'14_Counter_Adapter_Coverage.csv'));
writetable(step2.adapterLinks,fullfile(cfg.outputDir,'14_Counter_Adapter_Export_Links.csv'));
writetable(step2.neuronCandidates,fullfile(cfg.outputDir,'15_Transported_Neuron_Candidates.csv'));
writetable(step2.neuronCandidates,fullfile(cfg.outputDir,'15_Isolated_Neuron_Candidates.csv'));
writetable(step2.consumingChains,fullfile(cfg.outputDir,'15_Consuming_Neuron_Chains.csv'));
writetable(step2.chainCoverage,fullfile(cfg.outputDir,'15_Neuron_Chain_Coverage.csv'));
writetable(step2.compartmentAudit,fullfile(cfg.outputDir,'16_Neuron_Compartment_Audit.csv'));
writetable(step2.block0Alternatives,fullfile(cfg.outputDir,'16_Block0_Alternatives.csv'));
writetable(step2.block0Coverage,fullfile(cfg.outputDir,'16_Block0_Coverage.csv'));
writetable(step2.block0FilterAudit,fullfile(cfg.outputDir,'16_Block0_Filter_Audit.csv'));
writetable(step2.complete,fullfile(cfg.outputDir,'15_Complete_Step2_Circuits.csv'));
writeStep2Report(fullfile(cfg.outputDir,'Step_2_Results.txt'),step2,circuitContext,cfg);

report(fid,sprintf('\nDistinct intracellular cores: %d; transport/core combinations: %d.', ...
    numel(cores),height(candidates)));
report(fid,sprintf('Combinations with supplied direct-gating evidence: %d.', ...
    nnz(candidates.GateEvidenceRows>0)));
report(fid,sprintf('Combinations with all core ECs and a transporter gene rule present: %d.', ...
    nnz(candidates.AllCoreECResolved & candidates.TransportHasGPR)));
report(fid,sprintf('Connected structural triples: %d; exact formal quadratic stoichiometry: %d.', ...
    height(connected),nnz(connected.ExactFormalQuadratic)));
report(fid,'Connected triples do not establish S/F gating, quadratic kinetics or equal rate constants.');
report(fid,sprintf('EC-qualified fast/slow branch modules: %d; provisional wiring rows: %d; matched control interfaces: %d.', ...
    height(branchModules),height(branchWiring),nnz(branchWiring.BothControlInterfacesMatched)));
report(fid,sprintf('Counter mode=%s; components=%d; canceller/counter rows=%d; branch/canceller/counter rows=%d.', ...
    cfg.counterMode,height(counterComponents),height(counterConnections),height(counterWiring)));
report(fid,sprintf('Step 2: parents checked=%d; single-neuron rails=%d; complete cAMP-to-single-neuron circuits=%d. Report: Step_2_Results.txt.', ...
    step2.stats.parentsTested,step2.stats.neuron1Rails,height(step2.complete)));
if isempty(signals)
    report(fid,'RESULT STATUS: no imported I could be tested; the pathway loop was not entered.');
elseif isempty(oriented.rxn)
    report(fid,'RESULT STATUS: no eligible intracellular steps; no pathway could be tested.');
else
    if cfg.requireIImport
        report(fid,sprintf('RESULT STATUS: %d distinct imported I pools searched.',numel(signals)));
    else
        report(fid,sprintf('RESULT STATUS: %d distinct I pools searched (no import required).',numel(signals)));
    end
end
report(fid,'No S is inferred from stoichiometry, GPR membership or a positive regulatory edge.');
report(fid,'Seed checks use ALL substrates in the selected core, not graph connectivity alone.');
report(fid,'Currency co-participants are assumed available; their buffering is not demonstrated.');
report(fid,'Only selected-core seed dependence is checked; native-cell bypasses remain possible.');
report(fid,'No genome-scale FBA feasibility, thermodynamic or kinetic validation is claimed.');
report(fid,'A finite depth or reached state budget can hide other solutions.');
report(fid,'A solver-unresolved closure is untested, not a failed biological candidate.');
if cfg.requireIImport
    report(fid,'Spatial warning: external and internal I are distinct pools; import adds its own flux.');
    report(fid,'For each block: dI_inside/dt = J_T(S,I_outside,I_inside) + v_core(I_inside,X) + other terms.');
else
    report(fid,'requireIImport=false: no import term. dI/dt = v_core(I,X(S)) + other terms, with I(0) > 0 given.');
    report(fid,'Single-shot only: with eps = 0 the node fires once, and a resting arm mismatch kills it while it waits.');
end
report(fid,'This is not yet the single-pool Anderson equation +k*S*I.');
for n=1:min(height(candidates),cfg.printTop)
    d=details(n,:); c=cores(d.CoreIndex);
    if d.TransportRow>0
        t=transport(d.TransportRow,:);
        report(fid,sprintf('\nCANDIDATE %d | %s -> %s | core size %d', ...
            n,t.OutsideMetaboliteID,t.InsideMetaboliteID,numel(c.path)));
    else
        report(fid,sprintf('\nCANDIDATE %d | no I import searched | I = %s | core size %d', ...
            n,string(model.mets{c.seed}),numel(c.path)));
    end
    report(fid,['  Gate status: ' char(candidates.GateStatus(n))]);
    report(fid,['  Annotation status: ' char(candidates.AnnotationStatus(n))]);
    report(fid,['  S from supplied evidence: ' char(candidates.S_IDs(n))]);
    if d.TransportRow>0
        report(fid,['  r_T: ' char(t.ReactionID) ' | ' char(t.Equation)]);
        report(fid,['  Transporter GPR: ' char(t.GPR)]);
    else
        report(fid,'  r_T: none searched; I(0) assumed already present in the vessel (eps = 0).');
    end
    for k=1:numel(c.path)
        o=c.path(k); j=oriented.rxn(o);
        report(fid,sprintf('  core[%g] %s (%+d): %s | EC %s',c.flux(k), ...
            model.rxns{j},oriented.direction(o), ...
            equation(Rgross(:,j)*(oriented.direction(o)==1)+Pgross(:,j)*(oriented.direction(o)==-1), ...
            Pgross(:,j)*(oriented.direction(o)==1)+Rgross(:,j)*(oriented.direction(o)==-1),model), ...
            model.circuitResolvedEC(j)));
    end
    report(fid,['  Gross I accounting: 1 consumed -> 2 produced; net +1.']);
    report(fid,['  Core net equation: ' char(candidates.CoreNetEquation(n))]);
    if d.TransportRow>0
        report(fid,['  Full transport + core net: ' char(candidates.ModuleNetEquation(n))]);
    else
        report(fid,['  Module net (no transport): ' char(candidates.ModuleNetEquation(n))]);
    end
end
results = struct('options',cfg,'modelPath',modelPath,'modelLayout',layout,'poolAudit',poolAudit,'transport',transport, ...
    'cores',cores,'coreTable',coreTable,'orientedCorePool',oriented, ...
    'candidates',candidates,'details',details, ...
    'coverage',statsTable,'gateEvidence',gates,'gateAudit',gateAudit, ...
    'proteinMap',proteinMap,'proteinEdges',proteinEdges,'regulatorHypotheses',regulators, ...
    'linearRemoval',linearRemoval,'quadraticRemoval',quadraticRemoval, ...
    'removalCoverage',removalCoverage,'removalAudit',removalAudit, ...
    'connectedCandidates',connected,'connectionStats',connectionStats);
results.sGateImports=sGateImports;results.sGateStats=sGateStats;
results.gatedAndersonCandidates=gatedAnderson;results.gatedAndersonStats=gatedStats;
results.branchEC=branchEC;results.branchECAudit=branchECAudit;results.branchRails=branchRails;
results.branchModules=branchModules;results.branchReactionAssignments=branchAssignments;
results.branchWiring=branchWiring;results.branchSearchStats=branchStats;
results.branchWiringStats=branchWiringStats;
results.counterEC=counterEC;results.counterECAudit=counterECAudit;
results.counterComponents=counterComponents;results.counterComponentParts=counterParts;
results.counterFilterAudit=counterAudit;results.counterGates=counterGates;
results.counterGateAudit=counterGateAudit;results.counterConnections=counterConnections;
results.counterWiring=counterWiring;results.counterSearchStats=counterStats;
results.counterConnectionStats=counterConnectionStats;results.counterWiringStats=counterWiringStats;
results.counterConnectionParts=counterConnectionParts;results.counterConnectionAudit=counterConnectionAudit;
results.step2EC=step2EC;results.step2ECAudit=step2ECAudit;results.step2=step2;
results.upstreamCoverage=upstreamCoverage;
results.productionPairs=productionPairs;results.directProductionAudit=directProductionAudit;
results.minimumCoreTable=minimumCoreTable;results.minimumProductionPairs=minimumProductionPairs;
results.productionDepthCoverage=productionDepthCoverage;
results.sbmlMetadata=sbmlInfo;
results.cellFreeReactionMerge=cellFreeAudit;results.cellFreeMetaboliteMerge=cellFreeMetAudit;
results.cellFreeSummary=cellFreeInfo;
save(fullfile(cfg.outputDir,'Anderson_Production_Results.mat'),'results');
report(fid,['Finished. Outputs: ' cfg.outputDir]);
end

function cfg = defaultOptions()
cfg = struct('outputDir',fullfile(pwd,'Anderson_Transport_Production_Results'), ...
    'insideCompartment',"auto",'outsideCompartment',"e", ...
    'allowReverse',true,'requireTransportGPR',false,'requireCoreEC',false, ...
    'excludedECs',strings(0,1),'productionMode',"multi_reaction", ...
    'maxCoreReactions',10,'maxStatesPerI',20000,'maxProductionParticipants',16, ...
    'productionMinimumOnly',false, ...
    'rejectAtomImbalancedProduction',true,'requireResolvedProductionAtomBalance',false, ...
    'printTop',20,'onlineEC',false,'maxECLookups',100, ...
    'ecMappingFile','Reaction_EC_Evidence.tsv','brendaFile','brenda_2026_1.txt', ...
    'gateEvidenceFile','Gate_Evidence.tsv','requireGateEvidence',false, ...
    'requireProteinS',true,'onlineProteinDiscovery',false, ...
    'proteinMapFile','Transporter_Gene_To_Protein.tsv','proteinEdgesFile','Transporter_Protein_Edges.tsv', ...
    'organismTaxID',0,'organismName',"auto",'networkTimeout',15, ...
    'connectRemovalTerms',true,'requireFControlTopology',true,'requireRemovalIsolation',true, ...
    'allowPairLossAlternative',true,'maxRemovalParticipants',16, ...
    'maxTripletTests',200000,'maxConnectedCandidates',10000, ...
    'requireSGateIrreversible',false,'requireSGateIsolation',true,'maxSGateTests',200000, ...
    'requireIImport',false, ...
    'connectBranches',true,'branchUseVerifiedSeedMap',true, ...
    'branchECFile','Branch_EC_Evidence.tsv','branchInputID',"", ...
    'branchRequireDistinctECs',true,'branchRejectSharedKnownGenes',true, ...
    'branchRequireMatchedInterfaces',false,'branchMaxParticipants',16, ...
    'branchMaxRails',5000,'branchMaxPairTests',200000,'branchMaxModules',500, ...
    'branchMaxWiringTests',100000,'branchMaxWiring',10000, ...
    'connectCounter',true,'counterMode',"consume_I",'counterRequireFormationEC',false, ...
    'counterRequireRemoval',false,'counterRequireLeakEC',false,'counterUseVerifiedSeedMap',true, ...
    'counterECFile','Counter_EC_Evidence.tsv','counterGateEvidenceFile','Counter_Gate_Evidence.tsv', ...
    'counterCargoID',"",'counterRequireGateEvidence',false,'counterRequireIsolation',true, ...
    'counterDistinctECFromBranches',true,'counterRequireBoundaryFeed',false, ...
    'counterMaxParticipants',16,'counterMaxComponentTests',100000,'counterMaxComponents',5000, ...
    'counterMaxConnectionTests',200000,'counterMaxConnections',10000, ...
    'counterMaxWiringTests',200000,'counterMaxWiring',20000, ...
    'connectStep2',true,'step2UseVerifiedSeedMap',true,'step2ECFile','Step2_EC_Evidence.tsv', ...
    'neuronPlacement',"inside",'neuronCompartment',"auto", ...
    'neuronChainMinDepth',2,'neuronChainMaxDepth',6,'neuronChainMaxStates',2000, ...
    'neuronDiagnoseRejectedChains',true,'neuronDiagnosticMaxStates',500, ...
    'neuronDiagnosticMaxTests',5000,'neuronChainPrintTop',3, ...
    'neuronTransportRequireGPR',false,'neuronTransportMaxTests',100000, ...
    'neuronTransportMaxLinksPerParent',1000, ...
    'neuronAdapterMaxDepth',3,'neuronAdapterMaxStates',2000,'neuronAdapterMaxExpansions',20000, ...
    'step2DistinctECs',true,'step2PrintTop',10,'step2PrintPartialWhenEmpty',true, ...
    'step2MaxParents',1000,'step2MaxSeconds',21600,'step2ShareBlock0Search',true, ...
    'step2SearchBlock0',false,'block0RequireEC',false, ...
    'step2MaxParticipants',16,'neuronMaxTestsPerRail',50000, ...
    'neuronMaxRails',2000,'step2MaxNeuronTests',100000,'step2MaxCompleteMotifs',1000, ...
    'block0MaxDepth',20,'block0MaxStates',50000,'block0MaxExpansions',100000, ...
    'block0MaxAlternatives',10,'block0MaxParticipants',16, ...
    'cellFreeMerge',true,'cancellerI',"pyr",'counterUseEmbeddediMM904EC',true, ...
    'excludedReactionIDs',["THZPSN1_SC";"THZPSN2_SC"], ...
    'currencyBaseIDs',block0DefaultCurrencyIDs());
end

function validateOptions(c)
for key={'maxCoreReactions','maxStatesPerI','printTop','maxECLookups','maxProductionParticipants', ...
        'maxRemovalParticipants','maxTripletTests','maxConnectedCandidates','maxSGateTests', ...
        'branchMaxParticipants','branchMaxRails','branchMaxPairTests','branchMaxModules', ...
        'branchMaxWiringTests','branchMaxWiring','counterMaxParticipants','counterMaxComponentTests', ...
        'counterMaxComponents','counterMaxConnectionTests','counterMaxConnections', ...
        'counterMaxWiringTests','counterMaxWiring','step2PrintTop','step2MaxParents', ...
        'step2MaxParticipants','neuronMaxTestsPerRail','neuronMaxRails','step2MaxNeuronTests', ...
        'step2MaxCompleteMotifs','block0MaxDepth','block0MaxStates','block0MaxExpansions','step2MaxSeconds', ...
        'block0MaxAlternatives','block0MaxParticipants', ...
        'neuronTransportMaxTests','neuronTransportMaxLinksPerParent', ...
        'neuronAdapterMaxStates','neuronAdapterMaxExpansions','neuronChainMinDepth', ...
        'neuronChainMaxDepth','neuronChainMaxStates','neuronDiagnosticMaxStates', ...
        'neuronDiagnosticMaxTests','neuronChainPrintTop'}
    x=c.(key{1});
    assert(isnumeric(x) && isscalar(x) && x>=1 && (isinf(x)||x==floor(x)), ...
        'Anderson:Option','%s must be a positive integer or Inf.',key{1});
end
assert(isfinite(c.maxCoreReactions),'Anderson:Depth','maxCoreReactions must be finite.');
assert(isfinite(c.block0MaxDepth),'Anderson:Block0Depth','block0MaxDepth must be finite.');
assert(isfinite(c.neuronChainMaxDepth)&&c.neuronChainMinDepth>=2 ...
    &&c.neuronChainMinDepth<=c.neuronChainMaxDepth,'Anderson:NeuronDepth', ...
    'Neuron conversion depths must satisfy 2 <= neuronChainMinDepth <= finite neuronChainMaxDepth.');
assert(isnumeric(c.neuronAdapterMaxDepth)&&isscalar(c.neuronAdapterMaxDepth) ...
    &&isfinite(c.neuronAdapterMaxDepth)&&c.neuronAdapterMaxDepth>=0 ...
    &&c.neuronAdapterMaxDepth==floor(c.neuronAdapterMaxDepth), ...
    'Anderson:AdapterDepth','neuronAdapterMaxDepth must be a nonnegative finite integer.');
assert(~isempty(string(c.insideCompartment)) && isscalar(string(c.outsideCompartment)) ...
    && all(string(c.insideCompartment(:))~=string(c.outsideCompartment)),'Anderson:Compartments', ...
    'Inside labels must be nonempty and distinct from the outside compartment.');
for key={'allowReverse','requireTransportGPR','requireCoreEC','onlineEC', ...
        'requireGateEvidence','requireProteinS','onlineProteinDiscovery', ...
        'connectRemovalTerms','requireFControlTopology','requireRemovalIsolation','allowPairLossAlternative', ...
        'requireSGateIrreversible','requireSGateIsolation','requireIImport', ...
        'connectBranches','branchUseVerifiedSeedMap','branchRequireDistinctECs', ...
        'branchRejectSharedKnownGenes','branchRequireMatchedInterfaces','connectCounter', ...
        'counterRequireLeakEC','counterUseVerifiedSeedMap','counterRequireGateEvidence', ...
        'counterRequireIsolation','counterDistinctECFromBranches','counterRequireBoundaryFeed', ...
        'counterRequireFormationEC','counterRequireRemoval','connectStep2','step2UseVerifiedSeedMap', ...
        'step2DistinctECs','step2PrintPartialWhenEmpty','rejectAtomImbalancedProduction', ...
        'requireResolvedProductionAtomBalance','productionMinimumOnly','neuronTransportRequireGPR', ...
        'neuronDiagnoseRejectedChains','cellFreeMerge','counterUseEmbeddediMM904EC','step2ShareBlock0Search', ...
        'step2SearchBlock0','block0RequireEC'}
    assert(islogical(c.(key{1})) && isscalar(c.(key{1})), ...
        'Anderson:Flag','%s must be a logical scalar.',key{1});
end
assert(isstring(string(c.cancellerI))&&~any(ismissing(string(c.cancellerI))), ...
    'Anderson:CancellerI','cancellerI must be "" or a string (array) of BaseIDs or metabolite IDs.');
assert(isscalar(string(c.counterMode))&&ismember(string(c.counterMode),["consume_I","membrane"]), ...
    'Anderson:CounterMode','counterMode must be consume_I or membrane.');
assert(isscalar(string(c.productionMode))&&ismember(string(c.productionMode),["single_reaction","multi_reaction"]), ...
    'Anderson:ProductionMode','productionMode must be single_reaction or multi_reaction.');
assert(isscalar(string(c.neuronPlacement))&&ismember(string(c.neuronPlacement),["outside","inside"]), ...
    'Anderson:NeuronPlacement','neuronPlacement must be outside or inside.');
assert(isscalar(string(c.neuronCompartment))&&strlength(string(c.neuronCompartment))>0, ...
    'Anderson:NeuronCompartment','neuronCompartment must be one compartment label or auto.');
if c.neuronPlacement=="outside"
    target=step2NeuronCompartment(c);
    assert(~ismember(target,string(c.insideCompartment)), ...
        'Anderson:NeuronCompartment','The outside neuron compartment must differ from every inside compartment.');
end
end

function root = addCodeParent()
% Inline version of the attached addingPathParentFolderByName helper.
base=pwd; root=pwd;
while true
    [up,nm]=fileparts(base);
    if strcmpi(nm,'code'), root=base; addpath(genpath(root)); break; end
    if isempty(up)||strcmp(up,base), break; end
    base=up;
end
end

function model = normalizeModel(model)
assert(isfield(model,'S') && isfield(model,'mets') && isfield(model,'rxns'), ...
    'Anderson:Model','Model must contain S, mets and rxns.');
model.mets=cellstr(string(model.mets(:))); model.rxns=cellstr(string(model.rxns(:)));
[m,n]=size(model.S);
assert(m==numel(model.mets) && n==numel(model.rxns),'Anderson:Dimensions','S dimensions mismatch.');
assert(numel(unique(model.mets))==m && numel(unique(model.rxns))==n, ...
    'Anderson:IDs','Model metabolite and reaction IDs must be unique.');
assert(all(isfinite(nonzeros(model.S))),'Anderson:Stoichiometry','Nonfinite stoichiometry.');
model.S=sparse(model.S);
assert(isfield(model,'lb') && isfield(model,'ub'),'Anderson:Bounds','Model needs explicit lb and ub.');
model.lb=model.lb(:); model.ub=model.ub(:);
assert(numel(model.lb)==n && numel(model.ub)==n && all(model.lb<=model.ub) ...
    && ~any(isnan(model.lb)|isnan(model.ub)),'Anderson:Bounds','Invalid reaction bounds.');
for key={'rxnNames','grRules','rules'}
    if ~isfield(model,key{1}), model.(key{1})=repmat({''},n,1); end
    assert(numel(model.(key{1}))==n,'Anderson:Annotation','%s size mismatch.',key{1});
    model.(key{1})=cellstr(string(model.(key{1})(:)));
end
if ~isfield(model,'genes'), model.genes=cell(0,1); end
end

function [base,comp,layout,metadata] = readModelLayout(model,path)
% Preserve model species; canonical BaseID is used ONLY for cargo matching
% and currency filtering. Different compartments never share a state index.
ids=string(model.mets(:)); [~,comp]=splitCompartments(ids);metadata=[];
source=repmat("identifier_suffix",numel(ids),1);source(comp=="")="unresolved";
modelName="unspecified";
for key={'modelName','description'}
    if isfield(model,key{1})
        value=string(model.(key{1}));
        if isscalar(value) && ~ismissing(value) && value~="",modelName=value;break;end
    end
end
if isfield(model,'metComps') && isfield(model,'comps') && numel(model.metComps)==numel(ids)
    labels=string(model.comps(:));ix=full(model.metComps(:));
    valid=isfinite(ix)&ix==floor(ix)&ix>=1&ix<=numel(labels);
    comp(valid)=labels(ix(valid));source(valid)="model.metComps/model.comps";
end
if ~isempty(path) && endsWith(lower(string(path)),[".xml",".sbml"])
    try
        metadata=readSBMLNative(path);
        if metadata.modelName~="",modelName=metadata.modelName;end
        for k=1:numel(metadata.speciesIDs)
            i=uniqueIDMatch(metadata.speciesIDs(k),ids,'M_');
            if isempty(i),continue;end
            value=metadata.speciesCompartments(k);
            if value~="",comp(i)=value;source(i)="SBML_species_compartment";end
        end
    catch ME
        warning('Anderson:Layout','SBML compartment reading unavailable: %s.',ME.message);
    end
end
raw=regexprep(ids,'^M_','');base=raw;
for i=1:numel(ids)
    if comp(i)=="",continue;end
    label=regexptranslate('escape',char(comp(i)));
    base(i)=regexprep(base(i),['(?:\[' label '\]|_' label ')$'],'');
end
% Exact wrapper convention found in the attached model_aux.xml: AUXOS and
% PROTS blocks plus an original organelle code C/M/X/N/R/V/G. Strip only
% those wrappers for chemical identity matching, not arbitrary ID prefixes.
isU=comp=="u";isX=comp=="x";
communityIDs=any(isU)&&any(isX)&&all(startsWith(base(isU),"AUXOS")) ...
    && all(startsWith(base(isX),"PROTS"));
if communityIDs
    base(isU)=regexprep(base(isU),'^AUXOS[CMXNRVG]','');
    base(isX)=regexprep(base(isX),'^PROTS[CMXNRVG]','');
end
layout=struct('modelName',modelName,'communityIDs',communityIDs,'source',source);
end

function cfg = resolveModelProfile(cfg,layout,comp)
available=unique(comp(comp~=""),'stable');
inside=string(cfg.insideCompartment);inside=inside(:);outside=string(cfg.outsideCompartment);
if isscalar(inside) && inside=="auto"
    if any(available=="c")
        inside="c";
    elseif layout.communityIDs
        inside=["u";"x"];
    else
        inside=setdiff(available,outside,'stable');
        assert(numel(inside)==1,'Anderson:ChooseCompartments', ...
            'Set insideCompartment in defaultOptions. Available compartments: %s.',joinText(available,'; '));
    end
end
assert(all(ismember(inside,available)) && any(available==outside), ...
    'Anderson:CompartmentMismatch','Configured inside/outside labels are absent. Available: %s.', ...
    joinText(available,'; '));
cfg.insideCompartment=inside;cfg.outsideCompartment=outside;
if string(cfg.organismName)=="auto",cfg.organismName=layout.modelName;end
assert(~cfg.onlineProteinDiscovery || cfg.organismTaxID>0,'Anderson:Taxon', ...
    'Online protein discovery requires an explicitly configured organismTaxID.');
end

function [R,P,known,audit,sbmlEC,sbmlBiGG,info] = readGrossStoichiometry(model,path,metadata)
% model.S alone cancels I in I+X -> 2I+W. NEVER invent that hidden catalyst.
% Exact R/P can alternatively be supplied in a model struct under these fields.
R=max(-model.S,0); P=max(model.S,0); n=numel(model.rxns);
if nargin<3,metadata=[];end % reuse the native parse already performed for layout
known=false(n,1); status=repmat("net_only_no_hidden_catalysts_inferred",n,1);
sbmlEC=strings(n,1);sbmlBiGG=strings(n,1);rawIDs=strings(n,1);
info=struct('parsedReactionCount',0,'matchedReactionCount',0,'sourceBiGGCount',0, ...
    'unmatchedReactionExamples',strings(0,1),'error',"",'reader',"native_MATLAB_no_Java");
if isfield(model,'circuitSourceBiGG') && numel(model.circuitSourceBiGG)==n
    sbmlBiGG=string(model.circuitSourceBiGG(:));
end
if isfield(model,'andersonGrossReactants') && isfield(model,'andersonGrossProducts')
    a=sparse(model.andersonGrossReactants); b=sparse(model.andersonGrossProducts);
    assert(isequal(size(a),size(R)) && isequal(size(b),size(P)) && ...
        all(nonzeros(a)>=0) && all(nonzeros(b)>=0) && ...
        all(isfinite(nonzeros(a))) && all(isfinite(nonzeros(b))) && ...
        norm((b-a)-model.S,Inf)<1e-9,'Anderson:Gross','Invalid supplied gross stoichiometry.');
    R=a; P=b; known(:)=true; status(:)="supplied_gross_matrices_net_consistent";
elseif ~isempty(path) && endsWith(lower(string(path)),[".xml",".sbml"])
    try
        if isempty(metadata),metadata=readSBMLNative(path);end
        info.parsedReactionCount=numel(metadata.reactions);
        assert(~isempty(metadata.reactions),'Anderson:XMLReactionNodes', ...
            'The XML reader found no reaction elements; no EC/source-ID search can run.');
        for k=1:numel(metadata.reactions)
            reaction=metadata.reactions(k);id=reaction.id;
            j=uniqueIDMatch(id,string(model.rxns),'R_');
            if isempty(j)
                if numel(info.unmatchedReactionExamples)<5
                    info.unmatchedReactionExamples(end+1,1)=id;
                end
                continue;
            end
            info.matchedReactionCount=info.matchedReactionCount+1;rawIDs(j)=id;
            sbmlEC(j)=joinText(unique(reaction.ec,'stable'),'; ');
            sourceIDs=unique(reaction.bigg,'stable');
            if numel(sourceIDs)==1,sbmlBiGG(j)=sourceIDs(1);end
            [a,ok1]=nativeSBMLSide(reaction.reactantIDs,reaction.reactantCoefficients,model.mets);
            [b,ok2]=nativeSBMLSide(reaction.productIDs,reaction.productCoefficients,model.mets);
            if reaction.grossSupported && ok1 && ok2 && norm((b-a)-model.S(:,j),Inf)<1e-9
                R(:,j)=a; P(:,j)=b; known(j)=true; status(j)="SBML_gross_net_consistent";
            else
                status(j)="SBML_"+reaction.grossStatus+"_unmapped_or_net_mismatch_using_net_only";
            end
        end
    catch ME
        info.error=string(ME.message);
        warning('Anderson:GrossRead','Gross SBML parsing unavailable: %s. Using net topology only.',ME.message);
    end
end
info.sourceBiGGCount=nnz(sbmlBiGG~="");
audit=table(string(model.rxns),known,status,rawIDs,sbmlBiGG, ...
    'VariableNames',{'ReactionID','GrossKnown','Status','RawSBMLReactionID','SourceReactionID'});
end

function [v,ok] = nativeSBMLSide(ids,coefficients,mets)
v=sparse(numel(mets),1);ok=true;
for k=1:numel(ids)
    i=uniqueIDMatch(ids(k),string(mets),'M_');coef=coefficients(k);
    if isempty(i)||~isfinite(coef)||coef<0,ok=false;return;end
    v(i)=v(i)+coef;
end
end

function metadata = readSBMLNative(path)
% Focused SBML metadata reader using MATLAB text operations only. No JVM,
% xmlread, Python process, extra toolbox or network request is required.
% A quote-aware XML tag stream and an element stack keep RDF, sides and
% reaction scopes separate. The COBRA loader remains responsible for model
% loading. Every recovered gross equation must also agree with model.S.
text=fileread(path);
text=regexprep(text,'(?s)<!--.*?-->','');
text=regexprep(text,'(?s)<!\[CDATA\[.*?\]\]>','');
text=regexprep(text,'(?s)<\?.*?\?>','');
assert(isempty(regexpi(text,'<!DOCTYPE\b','once')),'Anderson:NativeXML', ...
    'DOCTYPE/entity declarations are not supported by the native SBML metadata reader.');
tagPattern='<(?:"[^"]*"|''[^'']*''|[^''">])*>';
tags=regexp(text,tagPattern,'match');
template=struct('id',"",'reactantIDs',strings(0,1),'reactantCoefficients',zeros(0,1), ...
    'productIDs',strings(0,1),'productCoefficients',zeros(0,1), ...
    'ec',strings(0,1),'bigg',strings(0,1),'referenceIDs',strings(0,1), ...
    'grossSupported',true,'grossStatus',"constant_explicit");
metadata=struct('modelName',"",'speciesIDs',strings(0,1), ...
    'speciesCompartments',strings(0,1),'reactions',template([]), ...
    'reader',"native_MATLAB_no_Java");
stack=strings(0,1);qualifiedStack=strings(0,1);
current=0;reactionDepth=0;side=0;sideDepth=0;annotationDepth=0;
assignedSymbols=strings(0,1);
for t=1:numel(tags)
    tag=tags{t};head=regexp(tag,'^<\s*/?\s*([A-Za-z_][A-Za-z0-9_.:-]*)','tokens','once');
    if isempty(head),error('Anderson:NativeXML','Unsupported XML markup encountered.');end
    qualified=string(head{1});name=regexprep(qualified,'^.*:','');
    closing=~isempty(regexp(tag,'^<\s*/','once'));
    if closing
        assert(~isempty(stack)&&qualifiedStack(end)==qualified,'Anderson:NativeXML', ...
            'Unbalanced XML element: %s.',qualified);
        depth=numel(stack);
        if current>0&&depth==reactionDepth
            current=0;reactionDepth=0;side=0;sideDepth=0;annotationDepth=0;
        elseif depth==sideDepth
            side=0;sideDepth=0;
        elseif depth==annotationDepth
            annotationDepth=0;
        end
        stack(end)=[];qualifiedStack(end)=[];continue;
    end
    depth=numel(stack)+1;parent="";if ~isempty(stack),parent=stack(end);end
    emptyTag=~isempty(regexp(tag,'/\s*>$','once'));
    if name=="model"&&parent=="sbml"&&numel(stack)==1
        attrs=nativeSBMLAttributes(tag);value=nativeSBMLAttribute(attrs,'name');
        if value=="",value=nativeSBMLAttribute(attrs,'id');end
        metadata.modelName=value;
    elseif name=="species"&&parent=="listOfSpecies"&&numel(stack)==3&&stack(2)=="model"
        attrs=nativeSBMLAttributes(tag);
        metadata.speciesIDs(end+1,1)=nativeSBMLAttribute(attrs,'id');
        metadata.speciesCompartments(end+1,1)=nativeSBMLAttribute(attrs,'compartment');
    elseif name=="reaction"&&parent=="listOfReactions"&&numel(stack)==3&&stack(2)=="model"
        attrs=nativeSBMLAttributes(tag);
        metadata.reactions(end+1,1)=template;current=numel(metadata.reactions);
        metadata.reactions(current).id=nativeSBMLAttribute(attrs,'id');
        reactionDepth=depth;side=0;sideDepth=0;annotationDepth=0;
    elseif current>0&&numel(stack)==reactionDepth ...
            &&(name=="listOfReactants"||name=="listOfProducts")
        side=1;if name=="listOfProducts",side=2;end
        sideDepth=depth;
    elseif current>0&&name=="annotation"&&numel(stack)==reactionDepth
        annotationDepth=depth;
    elseif current>0&&side>0&&name=="speciesReference"&&numel(stack)==sideDepth
        attrs=nativeSBMLAttributes(tag);id=nativeSBMLAttribute(attrs,'species');
        raw=nativeSBMLAttribute(attrs,'stoichiometry');coef=1;
        if raw~="",coef=str2double(raw);end
        constant=lower(nativeSBMLAttribute(attrs,'constant'));
        denominator=nativeSBMLAttribute(attrs,'denominator');
        if id==""||~isfinite(coef)||coef<0||ismember(constant,["false","0"]) ...
                ||(denominator~=""&&str2double(denominator)~=1)
            metadata.reactions(current).grossSupported=false;
            metadata.reactions(current).grossStatus="unsupported_or_nonconstant_speciesReference";
        end
        referenceID=nativeSBMLAttribute(attrs,'id');
        if referenceID~=""
            metadata.reactions(current).referenceIDs(end+1,1)=referenceID;
        end
        if side==1
            metadata.reactions(current).reactantIDs(end+1,1)=id;
            metadata.reactions(current).reactantCoefficients(end+1,1)=coef;
        else
            metadata.reactions(current).productIDs(end+1,1)=id;
            metadata.reactions(current).productCoefficients(end+1,1)=coef;
        end
    elseif current>0&&side>0&&name=="stoichiometryMath"
        metadata.reactions(current).grossSupported=false;
        metadata.reactions(current).grossStatus="stoichiometryMath_not_evaluated";
    elseif current>0&&annotationDepth>0&&name=="li"
        attrs=nativeSBMLAttributes(tag);uri=char(nativeSBMLAttribute(attrs,'resource'));
        hit=regexp(uri,'(?:ec-code|ec_code|enzyme-commission)[/:]([1-7]\.[0-9]+\.[0-9]+\.[0-9]+)(?![0-9.\-])','tokens','once');
        if ~isempty(hit),metadata.reactions(current).ec(end+1,1)=string(hit{1});end
        hit=regexp(uri,'bigg\.reaction[/:]([^/?#\s]+)','tokens','once');
        if ~isempty(hit),metadata.reactions(current).bigg(end+1,1)=string(hit{1});end
    elseif name=="initialAssignment"||name=="assignmentRule"||name=="rateRule"
        attrs=nativeSBMLAttributes(tag);symbol=nativeSBMLAttribute(attrs,'symbol');
        if symbol=="",symbol=nativeSBMLAttribute(attrs,'variable');end
        if symbol~="",assignedSymbols(end+1,1)=symbol;end
    end
    if emptyTag
        if current>0&&depth==reactionDepth
            current=0;reactionDepth=0;side=0;sideDepth=0;annotationDepth=0;
        elseif depth==sideDepth
            side=0;sideDepth=0;
        elseif depth==annotationDepth
            annotationDepth=0;
        end
    else
        stack(end+1,1)=name;qualifiedStack(end+1,1)=qualified;
    end
end
assert(isempty(stack),'Anderson:NativeXML','An XML element was not closed.');
for k=1:numel(metadata.reactions)
    metadata.reactions(k).ec=unique(metadata.reactions(k).ec,'stable');
    metadata.reactions(k).bigg=unique(metadata.reactions(k).bigg,'stable');
    if any(ismember(metadata.reactions(k).referenceIDs,assignedSymbols))
        metadata.reactions(k).grossSupported=false;
        metadata.reactions(k).grossStatus="assigned_stoichiometry_not_evaluated";
    end
end
end

function attrs = nativeSBMLAttributes(tag)
% Only quoted XML attribute values are accepted; > inside a quoted value
% does not terminate the tag. Decode entities exactly once.
pattern='(?s)([A-Za-z_][A-Za-z0-9_.:-]*)\s*=\s*(["''])(.*?)\2';
parts=regexp(tag,pattern,'tokens');names=strings(numel(parts),1);values=names;
for k=1:numel(parts)
    p=parts{k};names(k)=string(p{1});raw='';
    if numel(p)>=3,raw=p{3};end
    values(k)=string(nativeXMLDecode(raw));
end
attrs=struct('names',names,'localNames',regexprep(names,'^.*:',''),'values',values);
end

function value = nativeSBMLAttribute(attrs,name)
index=find(attrs.names==string(name),1);
if isempty(index),index=find(attrs.localNames==string(name),1);end
value="";if ~isempty(index),value=attrs.values(index);end
end

function text = nativeXMLDecode(text)
[start,stop,tokens]=regexp(text,'&([^;\s]+);','start','end','tokens');
for k=numel(start):-1:1
    entity=tokens{k}{1};
    switch entity
        case 'amp',value='&';
        case 'lt',value='<';
        case 'gt',value='>';
        case 'quot',value='"';
        case 'apos',value='''';
        otherwise
            code=NaN;
            if startsWith(entity,'#x'),code=hex2dec(entity(3:end));
            elseif startsWith(entity,'#'),code=str2double(entity(2:end));end
            assert(isfinite(code)&&code>=1&&code<=1114111&&code==floor(code) ...
                &&~(code>=55296&&code<=57343),'Anderson:NativeXML', ...
                'Unsupported XML entity: &%s;',entity);
            if code<=65535
                value=char(code);
            else
                code=code-65536;value=char([55296+floor(code/1024),56320+mod(code,1024)]);
            end
    end
    text=[text(1:start(k)-1),value,text(stop(k)+1:end)]; %#ok<AGROW>
end
end

function ix = uniqueIDMatch(id,ids,prefix)
ix=find(ids==id); if numel(ix)==1,return;end
id=regexprep(id,['^' prefix],'');
id=replace(id,["_LSQBKT_","_RSQBKT_","_LPAREN_","_RPAREN_","_DASH_"], ...
    ["[","]","(",")","-"]);
clean=regexprep(ids,['^' prefix],''); ix=find(clean==id);
if numel(ix)==1,return;end
% COBRA converts raw SBML metabolite suffixes _u/_x/_e to [u]/[x]/[e].
% Match that spelling only for metabolites; never merge compartment pools.
if strcmp(prefix,'M_')
    clean=replace(clean,["_LSQBKT_","_RSQBKT_"],["[","]"]);
    clean=regexprep(clean,'\[([A-Za-z_][A-Za-z0-9_]*)\]$','_$1');
    id=regexprep(id,'\[([A-Za-z_][A-Za-z0-9_]*)\]$','_$1');
    ix=find(clean==id);
end
if numel(ix)~=1,ix=[];end
end

function [O,T,A] = buildPools(model,base,comp,currency,R,P,known,cfg)
n=numel(model.rxns); m=numel(model.mets);
O=struct('rxn',zeros(1,0),'direction',zeros(1,0),'R',sparse(m,0), ...
    'P',sparse(m,0),'N',sparse(m,0),'consumers',{cell(m,1)},'grossKnown',false(1,0));
rows=cell(0,17);auditRows=cell(0,9);
% Lumped multi-enzyme pseudo-reactions listed in cfg.excludedReactionIDs
% (matched after removing the AUXOS/PROTS wrapper) never enter any pool.
excludedRx=false(n,1);
if isfield(cfg,'excludedReactionIDs')&&~isempty(cfg.excludedReactionIDs)
    unwrapped=regexprep(regexprep(string(model.rxns(:)),'^R_',''),'^(AUXOS|PROTS)','');
    excludedRx=ismember(lower(unwrapped),lower(string(cfg.excludedReactionIDs(:))));
end
for j=1:n
    if excludedRx(j)
        auditRows(end+1,:)={string(model.rxns{j}),0,"excluded_by_excludedReactionIDs",false, ...
            false,false,false,0,0}; %#ok<AGROW>
        continue;
    end
    dirs=[]; if model.ub(j)>0,dirs=1;end
    if cfg.allowReverse && model.lb(j)<0,dirs=[dirs -1];end %#ok<AGROW>
    [gpr,glist,~]=decodeGPR(string(model.grRules{j}),string(model.rules{j}),string(model.genes));
    ec=completeECTokens(model.circuitResolvedEC(j));
    for direction=dirs
        a=R(:,j); b=P(:,j); if direction<0,tmp=a;a=b;b=tmp;end
        net=b-a; participants=find(a+b>1e-12);
        % Each reaction must belong to ONE selected inside block. A reaction
        % spanning u and x is not admitted as an intracellular core step.
        insideOnly=~isempty(participants)&&ismember(comp(participants(1)),cfg.insideCompartment) ...
            &&all(comp(participants)==comp(participants(1)));
        coreAdmitted=false;importMatches=0;importAdmitted=0;
        if insideOnly && any(a>1e-12) && any(b>1e-12)
            ecOK=~cfg.requireCoreEC||~isempty(ec);
            if ecOK && ~any(ismember(ec,cfg.excludedECs))
                O.rxn(end+1)=j; O.direction(end+1)=direction;
                O.R(:,end+1)=a; O.P(:,end+1)=b; O.N(:,end+1)=net;
                O.grossKnown(end+1)=known(j);
                coreAdmitted=true;
            end
        end
        outside=find(comp==cfg.outsideCompartment & net < -1e-12 & ~currency);
        for e=reshape(outside,1,[])
            ins=find(ismember(comp,cfg.insideCompartment) & base==base(e) & net>1e-12 & ~currency);
            for c=reshape(ins,1,[])
                if abs(net(c)+net(e))>1e-9,continue;end % equal cargo copies
                importMatches=importMatches+1;
                if cfg.requireTransportGPR && strlength(gpr)==0,continue;end
                importAdmitted=importAdmitted+1;
                noncargo=setdiff(participants,[e c]);
                coupled=noncargo(~currency(noncargo));
                mechanism="same_species_import";
                if ~isempty(coupled),mechanism="coupled_import_not_gating_proof";end
                rows(end+1,:)={string(model.rxns{j}),string(model.rxnNames{j}),j,direction, ...
                    string(model.mets{e}),string(model.mets{c}),e,c,net(c), ...
                    gpr,glist,mechanism,equation(a,b,model), ...
                    joinText(string(model.mets(coupled)),';'),model.lb(j),model.ub(j),known(j)}; %#ok<AGROW>
            end
        end
        auditRows(end+1,:)={string(model.rxns{j}),direction, ...
            joinText(unique(comp(participants)),';'),insideOnly,~isempty(ec), ...
            strlength(gpr)>0,coreAdmitted,importMatches,importAdmitted}; %#ok<AGROW>
    end
end
for i=1:m,O.consumers{i}=find(O.R(i,:)>1e-12);end
T=cellTable(rows,{'ReactionID','ReactionName','ModelIndex','Direction', ...
    'OutsideMetaboliteID','InsideMetaboliteID','OutsideIndex','InsideIndex','ImportCoefficient', ...
    'GPR','Genes','TransportTopology','Equation','OtherNoncurrencyParticipants', ...
    'LowerBound','UpperBound','GrossKnown'}, ...
    {'string','string','double','double','string','string','double','double','double', ...
    'string','string','string','string','string','double','double','logical'});
A=cellTable(auditRows,{'ReactionID','Direction','Compartments','SingleInsideBlock', ...
    'HasCompleteEC','HasGPR','CoreAdmitted','ImportCargoMatches','ImportRowsAdmitted'}, ...
    {'string','double','string','logical','logical','logical','logical','double','double'});
end

function out = emptyCores()
out=struct('seed',{},'path',{},'flux',{},'net',{},'feed',{},'waste',{},'buffered',{},'grossExact',{});
end

function A = directProductionAuditTable(rows)
A=cellTable(rows,{'I_ID','ReactionID','Direction','GrossKnown','PrimitiveScale','I_GrossConsumed', ...
    'I_GrossProduced','X_SubstrateIDs','X_Coefficients','W_ProductIDs','W_Coefficients', ...
    'AtomBalanceStatus','AtomImbalance','PrimitiveEquation','Status'}, ...
    {'string','string','double','logical','double','double','double','string','string','string', ...
    'string','string','string','string','string'});
end

function [out,stats,A] = findSingleReactionCores(seed,O,currency,cfg,model)
% Do not turn net X->I into I+X->2I by adding an invented catalyst.
% Primitive normalization removes arbitrary column scaling but does not
% accept a primitive 2I+X->4I event as a one-copy reaction.
out=emptyCores();rows=cell(0,15);tol=1e-8;
stats=struct('states',0,'closures',0,'capped',false,'unsolved',0);
if isempty(model)
    model=struct('mets',{cellstr("met_"+string((1:size(O.R,1))'))}, ...
        'rxns',{cellstr("reaction_"+string((1:max([0 O.rxn]))'))});
end
if isfield(model,'circuitBaseIDs')&&numel(model.circuitBaseIDs)==numel(model.mets)
    base=string(model.circuitBaseIDs(:));
else
    [base,~]=splitCompartments(string(model.mets));
end
for o=reshape(O.consumers{seed},1,[])
    if stats.states>=cfg.maxStatesPerI,stats.capped=true;break;end
    stats.states=stats.states+1;j=O.rxn(o);a=O.R(:,o);b=O.P(:,o);parts=find(a+b>tol);
    scale=NaN;ci=NaN;pi=NaN;feed=[];waste=[];net=sparse(size(a,1),1);
    atomStatus="not_checked";imbalance="";eq="";status="accepted_single_reaction_I_plus_X_to_2I_plus_W";
    if ~O.grossKnown(o),status="gross_I_participation_not_known";
    elseif numel(parts)>cfg.maxProductionParticipants,status="production_participant_limit";
    else
        [scale,ok]=primitiveReactionScale(a,b);
        if ~ok,status="primitive_normalization_unresolved";
        else
            ci=full(a(seed))*scale;pi=full(b(seed))*scale;net=(b-a)*scale;
            feed=find(net < -tol);waste=find(net>tol);waste(waste==seed)=[];
            eq=equation(a*scale,b*scale,model);
            if b(seed)>tol,stats.closures=stats.closures+1;end
            if abs(ci-1)>tol||abs(pi-2)>tol,status="not_primitive_1I_consumed_2I_produced";
            elseif isempty(feed),status="no_additional_consumed_X_substrate";
            elseif any(base(parts)==base(seed)&parts~=seed)
                status="another_I_pool_cannot_be_counted_as_X";
            else
                [atomStatus,imbalance]=productionAtomBalance(net,model);
                if cfg.rejectAtomImbalancedProduction&&atomStatus=="unbalanced_atoms"
                    status="production_atom_imbalance";
                elseif cfg.requireResolvedProductionAtomBalance&&atomStatus~="balanced_atoms"
                    status="production_atom_balance_unresolved";
                else
                    [core,valid,unresolved]=validateCore(seed,o,O,currency);
                    stats.unsolved=stats.unsolved+unresolved;
                    if valid,out(end+1)=core; %#ok<AGROW>
                    else,status="additional_substrate_or_seed_dependence_check_failed";end
                end
            end
        end
    end
    rows(end+1,:)={string(model.mets{seed}),string(model.rxns{j}),O.direction(o),O.grossKnown(o), ...
        scale,ci,pi,joinText(string(model.mets(feed)),';'),joinText(string(full(-net(feed))),';'), ...
        joinText(string(model.mets(waste)),';'),joinText(string(full(net(waste))),';'), ...
        atomStatus,imbalance,eq,status}; %#ok<AGROW>
end
A=directProductionAuditTable(rows);
end

function [status,detail] = productionAtomBalance(net,model)
% Conservative atom ledger. Missing, polymeric, grouped or unsupported
% formulas are UNRESOLVED, never silently treated as balanced. Charge,
% energy, thermodynamic and kinetic feasibility are not established here.
status="unresolved_no_metabolite_formulas";detail="";
if ~isfield(model,'metFormulas')||numel(model.metFormulas)~=numel(model.mets),return;end
if ~iscell(model.metFormulas)&&~isstring(model.metFormulas)
    status="unresolved_unsupported_formula_field";return;
end
elements=split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og");
balance=zeros(numel(elements),1);parts=find(abs(net)>1e-10);
for i=reshape(parts,1,[])
    if iscell(model.metFormulas)
        value=model.metFormulas{i};
        if ~(ischar(value)||isstring(value))
            status="unresolved_unsupported_formula_entry";return;
        end
        raw=string(value);
    else
        raw=model.metFormulas(i);
    end
    if isempty(raw),status="unresolved_missing_formula";return;end
    if ~isscalar(raw),status="unresolved_unsupported_formula_entry";return;end
    if ismissing(raw)||raw=="",status="unresolved_missing_formula";return;end
    formula=char(strtrim(raw));
    if isempty(regexp(formula,'^(?:[A-Z][a-z]?[0-9]*(?:\.[0-9]+)?)+$','once'))
        status="unresolved_unsupported_formula";return;
    end
    tokens=regexp(formula,'([A-Z][a-z]?)([0-9]*(?:\.[0-9]+)?)','tokens');
    for k=1:numel(tokens)
        e=find(elements==string(tokens{k}{1}),1);
        if isempty(e),status="unresolved_non_element_or_generic_group";return;end
        amount=1;if numel(tokens{k})>=2&&~isempty(tokens{k}{2}),amount=str2double(tokens{k}{2});end
        if ~isfinite(amount)||amount<=0,status="unresolved_atom_count";return;end
        balance(e)=balance(e)+full(net(i))*amount;
    end
end
bad=find(abs(balance)>1e-7);
status="balanced_atoms";
if ~isempty(bad)
    status="unbalanced_atoms";pieces=strings(numel(bad),1);
    for k=1:numel(bad),pieces(k)=elements(bad(k))+":"+string(sprintf('%+.12g',balance(bad(k))));end
    detail=joinText(pieces,';');
end
end

function P = describeImportProductionPairs(C,D,cores,O,T,model)
rows=cell(0,23);
for p=1:height(C)
    core=cores(D.CoreIndex(p));single=numel(core.path)==1;
    if D.TransportRow(p)>0
        t=T(D.TransportRow(p),:);
        tID=t.ReactionID;tDir=t.Direction;tOut=t.OutsideMetaboliteID;
        tIn=t.InsideMetaboliteID;tEq=t.Equation;
    else
        tID="none_no_I_import";tDir=0;tOut="";
        tIn=string(model.mets{core.seed});
        tEq="no import searched; I(0) assumed present in the vessel";
    end
    reactionID=joinText(string(model.rxns(O.rxn(core.path))),';');
    direction=joinText(string(O.direction(core.path)),';');gross="not_applicable_multi_reaction_see_R1_EffectiveModuleEquation";
    if single,gross=equation(O.R(:,core.path)*core.flux(:),O.P(:,core.path)*core.flux(:),model);end
    [atomStatus,imbalance]=productionAtomBalance(core.net,model);
    rows(end+1,:)={C.CandidateID(p),tID,tDir,tOut,tIn, ...
        tEq,reactionID,direction,single,gross, ...
        joinText(string(model.mets(core.feed)),';'),joinText(string(full(-core.net(core.feed))),';'), ...
        joinText(string(model.mets(core.waste)),';'),joinText(string(full(core.net(core.waste))),';'), ...
        atomStatus,imbalance,C.GateStatus(p),C.S_IDs(p), ...
        "I0_in_vessel_no_import_S_gate_on_X_importer_kSI_not_established", ...
        productionModuleEquation(core,O,model),numel(core.path), ...
        full(O.R(core.seed,core.path)*core.flux(:)), ...
        full(O.P(core.seed,core.path)*core.flux(:))}; %#ok<AGROW>
end
P=cellTable(rows,{'CandidateID','TransporterID','TransportDirection','I_Outside','I_Inside','ImportEquation', ...
    'R1_ReactionID','R1_Direction','SingleProductionReaction','R1_PrimitiveGrossEquation','X_SubstrateIDs', ...
    'X_Coefficients','W_ProductIDs','W_Coefficients','AtomBalanceStatus','AtomImbalance','ImportGateStatus', ...
    'Supplied_S_IDs','ProductionControlStatus','R1_EffectiveModuleEquation', ...
    'ProductionReactionCount','I_ConsumedPerModule','I_ProducedPerModule'}, ...
    {'double','string','double','string','string','string','string','string','logical','string', ...
    'string','string','string','string','string','string','string','string','string', ...
    'string','double','double','double'});
end

function value = productionModuleEquation(core,O,model)
% Cancel non-I intermediates but retain the I consumption/production
% actually present in the selected reactions; never invent I on both sides.
a=max(-core.net,0);b=max(core.net,0);
a(core.seed)=full(O.R(core.seed,core.path)*core.flux(:));
b(core.seed)=full(O.P(core.seed,core.path)*core.flux(:));
value=equation(a,b,model);
end

function writeImportProductionReport(path,P,A,cfg)
fid=fopen(path,'w');assert(fid>0,'Cannot write Eq. Block 1 report.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
report(fid,'EQ. BLOCK 1: I_OUT IMPORT + I_IN + X -> 2 I_IN + W');
report(fid,['Production search mode: ' char(cfg.productionMode)]);
report(fid,sprintf('Transport/production pairs=%d; production depth bound=%d.',height(P),cfg.maxCoreReactions));
if cfg.productionMode=="single_reaction"
    report(fid,sprintf('Direct-only mode: %d production orientations inspected; one explicit primitive 1I->2I reaction required.',height(A)));
else
    report(fid,'Module mode: one or several recorded reactions may implement effective I+X->2I+W.');
    report(fid,'The direct-only audit is not applicable in module mode; see 02_Search_Coverage.csv and 02_Autocatalytic_Cores.csv.');
end
report(fid,'Actual selected-reaction I consumption=1 and production=2 are required. Net X->I is never relabelled as autocatalysis.');
report(fid,'All noncurrency internal intermediates must balance; all selected steps must be reachable with I and supplied feeds, but cannot produce I without its seed.');
report(fid,'X denotes all additional net-consumed substrates; W denotes net-produced co-products. Full gross equations retain catalysts and co-participants.');
report(fid,'The exact imported I_inside state must be the I consumed and produced by R1. Other I pools cannot serve as an extra X.');
report(fid,'The optional direct-only mode rejects known atom imbalance by default. Module mode retains the original structural checks and reports its net atom ledger.');
report(fid,'Unsupported/missing formulas are unresolved. The direct-only strict formula option applies only in single_reaction mode.');
report(fid,'Atom balance is separate from charge/energy balance, flux feasibility, thermodynamic direction and kinetics.');
report(fid,'S denotes the proposed controller of I import; supplied gate evidence is reported, not invented from GPR or stoichiometry.');
report(fid,'Direct I import seeds/adds I. It does NOT by itself turn off R1 when S disappears. Production/X-supply control is still required for that behavior.');
report(fid,'The paired normalization imports one I and executes one normalized R1 turnover; it does not assert matched physiological fluxes.');
report(fid,'Mass balances: dI_inside/dt = J_import + v_R1_net - other losses; dX/dt includes consumption by R1.');
report(fid,'No k*S*I production law or measured reaction order is claimed by this structural search.');
for status=reshape(unique(A.Status,'stable'),1,[])
    report(fid,sprintf('  Production filter %s: %d',status,nnz(A.Status==status)));
end
for p=1:min(height(P),cfg.printTop)
    report(fid,sprintf('\nPAIR %d | imported seed I=%s',P.CandidateID(p),P.I_Inside(p)));
    report(fid,['  Transporter: ' char(P.TransporterID(p)) ' | ' char(P.ImportEquation(p))]);
    report(fid,sprintf('  R1 module: %d existing reaction(s); actual I accounting %.12g consumed -> %.12g produced.', ...
        P.ProductionReactionCount(p),P.I_ConsumedPerModule(p),P.I_ProducedPerModule(p)));
    report(fid,['  R1 reaction IDs: ' char(P.R1_ReactionID(p))]);
    report(fid,['  R1 effective equation: ' char(P.R1_EffectiveModuleEquation(p))]);
    if P.SingleProductionReaction(p)
        report(fid,['  R1 primitive gross equation: ' char(P.R1_PrimitiveGrossEquation(p))]);
    end
    report(fid,['  X substrates: ' char(P.X_SubstrateIDs(p)) ' | coefficients: ' char(P.X_Coefficients(p))]);
    report(fid,['  W products: ' char(P.W_ProductIDs(p)) ' | coefficients: ' char(P.W_Coefficients(p))]);
    report(fid,['  Atom balance: ' char(P.AtomBalanceStatus(p)) ' | ' char(P.AtomImbalance(p))]);
    report(fid,['  S-import gate status: ' char(P.ImportGateStatus(p))]);
end
end

function [out,stats,audit,depthAudit] = findCores(seed,O,currency,cfg,model)
if nargin<5,model=[];end
audit=directProductionAuditTable(cell(0,15));
seedID="met_"+string(seed);
if ~isempty(model),seedID=string(model.mets{seed});end
if cfg.productionMode=="single_reaction"
    [out,stats,audit]=findSingleReactionCores(seed,O,currency,cfg,model);
    complete=~stats.capped;
    depthAudit=productionDepthAuditTable({seedID,1,stats.states,stats.closures, ...
        numel(out),stats.unsolved,complete,stats.capped});
    stats=finishProductionMinimum(stats,out,complete,stats.unsolved, ...
        "direct_reactions_for_this_I_under_current_filters");return;
end
% Iterative deepening by number of DISTINCT model reactions. A whole layer
% is searched before any larger module is tested. Only simple metabolite
% cycles are enumerated: this is not a search of arbitrary branched modules.
out=emptyCores();stats=struct('states',0,'closures',0,'capped',false,'unsolved',0);
rows=cell(0,8);layerComplete=false(0,1);layerUnresolved=zeros(0,1);
seen=containers.Map('KeyType','char','ValueType','logical');
for depth=1:cfg.maxCoreReactions
    beforeStates=stats.states;beforeClosures=stats.closures;
    beforeFound=numel(out);beforeUnsolved=stats.unsolved;
    stackNode={seed};stackPath={[]};stackMets={seed};
    while ~isempty(stackNode)
        if stats.states>=cfg.maxStatesPerI,stats.capped=true;break;end
        node=stackNode{end};path=stackPath{end};chain=stackMets{end};
        stackNode(end)=[];stackPath(end)=[];stackMets(end)=[];
        stats.states=stats.states+1;
        for o=reshape(O.consumers{node},1,[])
            if any(O.rxn(path)==O.rxn(o)),continue;end
            chosen=[path o];
            if numel(chosen)==depth && O.P(seed,o)>1e-12
                key=sprintf('%d,',sort(chosen));
                if ~isKey(seen,key)
                    seen(key)=true;stats.closures=stats.closures+1;
                    [c,valid,unresolved]=validateCore(seed,chosen,O,currency);
                    stats.unsolved=stats.unsolved+unresolved;
                    if valid,out(end+1)=c;end %#ok<AGROW>
                end
            end
            if numel(chosen)>=depth,continue;end
            products=find(O.P(:,o)>1e-12 & ~currency);
            products=products(~ismember(products,chain));
            for p=reshape(products,1,[])
                stackNode{end+1}=p;stackPath{end+1}=chosen;stackMets{end+1}=[chain p]; %#ok<AGROW>
            end
        end
    end
    complete=isempty(stackNode);
    layerComplete(depth,1)=complete;
    layerUnresolved(depth,1)=stats.unsolved-beforeUnsolved;
    rows(end+1,:)={seedID,depth,stats.states-beforeStates, ...
        stats.closures-beforeClosures,numel(out)-beforeFound, ...
        layerUnresolved(depth),complete,stats.capped}; %#ok<AGROW>
    if stats.capped || (cfg.productionMinimumOnly && numel(out)>beforeFound),break;end
end
depthAudit=productionDepthAuditTable(rows);
stats=finishProductionMinimum(stats,out,layerComplete,layerUnresolved, ...
    "simple_cycle_modules_for_this_I_under_current_filters");
end

function T = productionDepthAuditTable(rows)
T=cellTable(rows,{'I_ID','ReactionCount','PathStates','ClosedModulesTested', ...
    'ValidModulesFound','SolverUnresolved','LayerComplete','StateBudgetReached'}, ...
    {'string','double','double','double','double','double','logical','logical'});
end

function stats = finishProductionMinimum(stats,cores,complete,unresolved,scope)
% A valid size-k module is minimum in this search class once ALL smaller
% layers are exhausted and resolved. Unfinished size-k ties do not invalidate
% that size bound, but are reported separately as MinimumLayerComplete=false.
complete=logical(complete(:));unresolved=unresolved(:);
stats.completedDepth=0;
for k=1:numel(complete)
    if ~complete(k),break;end
    stats.completedDepth=k;
end
stats.smallestFound=NaN;stats.minimumProven=false;stats.minimumLayerComplete=false;
stats.minimumScope=scope;
if isempty(cores)
    stats.minimumStatus="no_match_search_incomplete";
    if ~isempty(complete)&&all(complete)&&all(unresolved==0),stats.minimumStatus="no_match_through_tested_depth";end
    return;
end
sizes=arrayfun(@(c)numel(c.path),cores);k=min(sizes);stats.smallestFound=k;
stats.minimumProven=numel(complete)>=k-1 && all(complete(1:k-1)) && all(unresolved(1:k-1)==0);
stats.minimumLayerComplete=numel(complete)>=k && complete(k) && unresolved(k)==0;
stats.minimumStatus="smallest_found_shorter_search_incomplete";
if stats.minimumProven,stats.minimumStatus="minimum_proven_within_search_class";end
end

function T = tagProductionMinima(T,coverage,iField,countField)
T.SmallestFoundReactionCount=nan(height(T),1);
T.IsShortestFoundForI=false(height(T),1);T.MinimumProvenWithinSearchClass=false(height(T),1);
T.MinimumLayerComplete=false(height(T),1);T.MinimumStatus=strings(height(T),1);
T.MinimumScope=strings(height(T),1);ids=string(T.(iField));counts=T.(countField);
for k=1:height(T)
    r=find(coverage.I_ID==ids(k),1);if isempty(r),continue;end
    T.SmallestFoundReactionCount(k)=coverage.SmallestFoundReactionCount(r);
    T.IsShortestFoundForI(k)=counts(k)==coverage.SmallestFoundReactionCount(r);
    T.MinimumScope(k)=coverage.MinimumScope(r);
    T.MinimumStatus(k)="longer_module_retained_for_downstream_wiring";
    if T.IsShortestFoundForI(k)
        T.MinimumProvenWithinSearchClass(k)=coverage.MinimumProvenWithinSearchClass(r);
        T.MinimumLayerComplete(k)=coverage.MinimumLayerComplete(r);
        T.MinimumStatus(k)=coverage.MinimumStatus(r);
    end
end
end

function writeMinimumProductionReport(path,C,D,P,cfg)
fid=fopen(path,'w');assert(fid>0,'Anderson:MinimumReport','Cannot write minimum-production report.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
report(fid,'MINIMUM NUMBER OF PRODUCTION REACTIONS FOR EACH IMPORTED I');
report(fid,'Objective: minimize the number of distinct model reactions in the intracellular I+X->2I+W module.');
report(fid,'The I importer is separate and adds one reaction to an import/production pair. Reaction flux multipliers do not count as extra reactions.');
report(fid,'Search order: complete reaction-count layers 1, 2, 3, ... with the existing direction, compartment, currency and seed-dependence checks.');
report(fid,'Scope: simple metabolite cycles for each I under the current filters and numerical validation, not arbitrary branched reaction sets or the entire circuit.');
report(fid,'MinimumProvenWithinSearchClass=true requires a valid size-k module and exhaustive, solver-resolved searches of every smaller layer.');
report(fid,'MinimumLayerComplete separately records whether all size-k alternatives were examined and resolved.');
report(fid,'A budget hit or a smaller solver-unresolved case gives only a smallest-found result. No-match results are limited to the tested layers and search class.');
report(fid,sprintf('Mode=%s; maximum production steps=%d; shared per-I state budget=%g; productionMinimumOnly=%d.', ...
    cfg.productionMode,cfg.maxCoreReactions,cfg.maxStatesPerI,cfg.productionMinimumOnly));
report(fid,'Repeated prefix expansions at successive depths count toward the same per-I state budget. See 02_Production_Depth_Coverage.csv.');
if ~cfg.productionMinimumOnly
    report(fid,'Longer modules are retained for downstream wiring; the shortest isolated production module may not connect to the removal/counter blocks.');
else
    report(fid,'Only the first layer containing valid modules is retained per I; longer downstream-compatible alternatives are not searched.');
end
report(fid,sprintf('I pools searched=%d; I pools with candidates=%d; minimum-size certificates=%d; depth rows=%d.', ...
    height(C),nnz(~isnan(C.SmallestFoundReactionCount)),nnz(C.MinimumProvenWithinSearchClass),height(D)));
for k=1:height(C)
    report(fid,sprintf('I=%s | smallest=%g | minimum proven=%d | ties complete=%d | completed depth=%d | %s', ...
        C.I_ID(k),C.SmallestFoundReactionCount(k),C.MinimumProvenWithinSearchClass(k), ...
        C.MinimumLayerComplete(k),C.CompletedReactionDepth(k),C.MinimumStatus(k)));
end
if ~isempty(P)
    [~,order]=sortrows([P.ProductionReactionCount P.CandidateID],[1 2]);P=P(order,:);
end
for k=1:min(height(P),cfg.printTop)
    report(fid,sprintf('\nPAIR %d | I=%s | production steps=%d | minimum proven=%d', ...
        P.CandidateID(k),P.I_Inside(k),P.ProductionReactionCount(k),P.MinimumProvenWithinSearchClass(k)));
    report(fid,['  Importer: ' char(P.TransporterID(k))]);
    report(fid,['  Production reactions: ' char(P.R1_ReactionID(k))]);
    report(fid,['  Effective module: ' char(P.R1_EffectiveModuleEquation(k))]);
    report(fid,['  Scope/status: ' char(P.MinimumScope(k)) ' | ' char(P.MinimumStatus(k))]);
end
report(fid,'All shortest-found modules are in 02_Minimum_Production_Modules.csv; paired importers are in 03_Minimum_Import_Production_Pairs.csv.');
report(fid,'S-dependent gating, reaction orders and whole-cell feasibility are not established by a reaction-count minimum.');
end

function [c,valid,unresolved] = validateCore(seed,path,O,currency)
c=[];valid=false;unresolved=false;tol=1e-8;
a=O.R(:,path);b=O.P(:,path);N=b-a;
% Every noncurrency species both consumed and produced is an internal
% intermediate, except I. It must cancel in the overall core equation.
internal=find(any(a>tol,2)&any(b>tol,2)&~currency);
internal(internal==seed)=[];
A=[full(N(internal,:));full(a(seed,:));full(b(seed,:))];
target=[zeros(numel(internal),1);1;2];
[v,solved]=positiveFlux(A,target);
if ~solved,unresolved=true;return;end
if isempty(v),return;end
net=N*v;net(abs(net)<tol)=0;
if abs(net(seed)-1)>tol,return;end
feed=find(net < -tol); waste=find(net > tol); waste(waste==seed)=[];
if isempty(feed),return;end
% Currency metabolites participating catalytically can be externally buffered.
buffered=find(currency & any(a+b>tol,2));
available=false(size(N,1),1); available([feed;buffered])=true;
noSeed=substrateClosure(available,a,b);
if noSeed(seed),return;end % selectable chemistry can self-start from feeds
available(seed)=true;
[withSeed,fired]=substrateClosure(available,a,b);
if ~withSeed(seed)||~all(fired),return;end % unsupported co-substrate/catalyst
% A one-reaction autocatalyst requires uncancelled evidence, not a net guess.
if numel(path)==1 && ~O.grossKnown(path),return;end
c=struct('seed',seed,'path',path,'flux',v(:)', ...
    'net',net,'feed',feed,'waste',waste,'buffered',buffered, ...
    'grossExact',all(O.grossKnown(path)));
valid=true;
end

function [v,solved] = positiveFlux(A,b)
% Solve the caller's exact positive-flux accounting constraints A*v=b.
% Production uses gross I consumption=1/production=2; Block 0 instead
% fixes cAMP consumption=1 and balances its shared internal intermediates.
% Use the unique linear solution where possible, else LP. Flux multipliers
% are stoichiometric weights, NOT physiological rates or fitted constants.
tol=1e-8; n=size(A,2); v=[];solved=true;
if rank(A)==n
    candidate=A\b;
    if all(candidate>1e-9) && norm(A*candidate-b,Inf)<tol,v=candidate;end
    return;
end
if rank([A b])>rank(A),return;end
try
    if exist('linprog','file')~=0
        options=optimoptions('linprog','Display','none');
        [candidate,~,flag]=linprog(ones(n,1),[],[],A,b,ones(n,1)*1e-7,[],options);
        if flag==-2,return;end
        if flag<=0,solved=false;return;end
    elseif exist('solveCobraLP','file')~=0
        problem=struct('A',sparse(A),'b',b,'c',ones(n,1),'lb',ones(n,1)*1e-7, ...
            'ub',inf(n,1),'csense',repmat('E',size(A,1),1),'osense',1);
        answer=solveCobraLP(problem);
        if answer.stat==0,return;end
        if answer.stat~=1,solved=false;return;end
        candidate=answer.full;
    else
        solved=false;return;
    end
    if numel(candidate)==n && all(candidate>1e-9) && norm(A*candidate-b,Inf)<tol
        v=candidate;
    else
        solved=false;
    end
catch
    solved=false;
end
end

function [available,fired] = substrateClosure(available,a,b)
fired=false(size(a,2),1);changed=true;
while changed
    changed=false;
    for j=1:size(a,2)
        if ~fired(j) && all(available(a(:,j)>1e-12))
            fired(j)=true;available(b(:,j)>1e-12)=true;changed=true;
        end
    end
end
end

function T = describeCores(cores,O,model)
rows=cell(0,11);
for k=1:numel(cores)
    c=cores(k);
    rows(end+1,:)={k,string(model.mets{c.seed}),numel(c.path), ...
        joinText(string(model.rxns(O.rxn(c.path))),';'),joinText(string(O.direction(c.path)),';'), ...
        joinText(string(c.flux),';'),equation(max(-c.net,0),max(c.net,0),model), ...
        joinText(string(model.mets(c.feed)),';'),joinText(string(model.mets(c.waste)),';'), ...
        joinText(string(model.mets(c.buffered)),';'),c.grossExact}; %#ok<AGROW>
end
T=cellTable(rows,{'CoreIndex','I_ID','CoreReactionCount','ReactionIDs','Directions','FluxMultipliers', ...
    'NetEquation','FeedMetabolites','WasteMetabolites','AssumedBufferedCurrency','AllGrossKnown'}, ...
    {'double','string','double','string','string','string','string','string','string','string','logical'});
T.AllCoreECResolved=false(height(T),1);
for k=1:height(T)
    T.AllCoreECResolved(k)=all(strlength(strtrim(model.circuitResolvedEC(O.rxn(cores(k).path))))>0);
end
end

function names = gateColumns()
names={'TransportReactionID','Direction','ImportedMetaboliteID','TargetGene', ...
    'S_ID','S_Type','Effect','Mechanism','S_NotConsumed','Organism','Reference'};
end

function writeGateTemplate(path)
fid=fopen(path,'w');assert(fid>0,'Cannot write gate template');g=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',strjoin(gateColumns(),sprintf('\t')));
end

function [G,audit] = readGateEvidence(path,T,cfg)
names=gateColumns();G=cellTable(cell(0,numel(names)),names,repmat({'string'},1,numel(names)));
audit=G;audit.Accepted=false(0,1);audit.Status=strings(0,1);
if isempty(path),return;end
G=readTextTable(path);
assert(all(ismember(names,G.Properties.VariableNames)),'Anderson:GateColumns', ...
    'Gate evidence must have columns: %s',strjoin(names,', '));
G=G(:,names);audit=G;audit.Accepted=false(height(G),1);audit.Status=strings(height(G),1);
for i=1:height(G)
    row=G(i,:);status="accepted_supplied_direct_gating_evidence";
    tr=find(T.ReactionID==row.TransportReactionID & T.Direction==str2double(row.Direction) ...
        & T.InsideMetaboliteID==row.ImportedMetaboliteID);
    if isempty(tr),status="no_exact_transport_direction_cargo_match";
    elseif row.S_ID=="" || row.TargetGene=="" || row.Reference==""
        status="missing_signal_target_or_reference";
    elseif cfg.requireProteinS && lower(row.S_Type)~="protein"
        status="S_not_protein_under_current_policy";
    elseif ~strcmpi(row.Organism,cfg.organismName)
        status="organism_mismatch";
    elseif lower(row.Effect)~="opening" || ~ismember(lower(row.Mechanism), ...
            ["ligand_gating","protein_gating","pore_assembly"])
        status="not_direct_opening_mechanism";
    elseif ~asBool(row.S_NotConsumed)
        status="S_not_declared_unconsumed";
    else
        member=false;
        for t=reshape(tr,1,[])
            member=member||any(split(T.Genes(t),';')==row.TargetGene);
        end
        if ~member,status="target_gene_not_in_exact_transport_GPR";end
    end
    audit.Status(i)=status;
    audit.Accepted(i)=status=="accepted_supplied_direct_gating_evidence";
end
G=G(audit.Accepted,:);
end

function [T,D] = assembleCandidates(cores,O,transport,gates,model,cfg)
rows=cell(0,12);drows=cell(0,2);
for k=1:numel(cores)
    c=cores(k);
    if ~cfg.requireIImport
        % No Block A. One candidate per core, with no import reaction.
        coreEq=equation(max(-c.net,0),max(c.net,0),model);
        rows(end+1,:)={heightFromRows(rows)+1,"none_no_I_import",0, ...
            "",string(model.mets{c.seed}),k,numel(c.path), ...
            0,"","no_I_import_searched_S_gate_is_on_the_X_importer", ...
            coreEq,coreEq}; %#ok<AGROW>
        drows(end+1,:)={0,k}; %#ok<AGROW>
        continue;
    end
    for t=reshape(find(transport.InsideIndex==c.seed),1,[])
        tr=transport(t,:);j=tr.ModelIndex;
        if any(O.rxn(c.path)==j),continue;end
        ge=find(gates.TransportReactionID==tr.ReactionID & ...
            str2double(gates.Direction)==tr.Direction & gates.ImportedMetaboliteID==tr.InsideMetaboliteID);
        if cfg.requireGateEvidence && isempty(ge),continue;end
        status="S_unassigned_gating_unverified";
        if ~isempty(ge),status="supplied_direct_gating_evidence_kinetics_unverified";end
        % One imported molecule is illustrative normalization only. It does
        % not assert equality of transport and core fluxes in the ODE.
        moduleNet=c.net+tr.Direction*model.S(:,j)/tr.ImportCoefficient;
        rows(end+1,:)={heightFromRows(rows)+1,tr.ReactionID,tr.Direction, ...
            tr.OutsideMetaboliteID,tr.InsideMetaboliteID,k,numel(c.path), ...
            numel(ge),joinText(unique(gates.S_ID(ge)),';'),status, ...
            equation(max(-c.net,0),max(c.net,0),model), ...
            equation(max(-moduleNet,0),max(moduleNet,0),model)}; %#ok<AGROW>
        drows(end+1,:)={t,k}; %#ok<AGROW>
    end
end
T=cellTable(rows,{'CandidateID','TransportReactionID','Direction','I_Outside','I_Inside', ...
    'CoreIndex','CoreReactionCount','GateEvidenceRows','S_IDs','GateStatus', ...
    'CoreNetEquation','ModuleNetEquation'}, ...
    {'double','string','double','string','string','double','double','double','string','string','string','string'});
D=cellTable(drows,{'TransportRow','CoreIndex'},{'double','double'});
if ~isempty(T)
    [~,order]=sortrows([T.CoreReactionCount,-double(T.GateEvidenceRows>0),T.CandidateID],[1 2 3]);
    T=T(order,:);D=D(order,:);T.CandidateID=(1:height(T))';
end
T.TransportHasGPR=false(height(T),1);T.AllCoreECResolved=false(height(T),1);
T.AnnotationStatus=strings(height(T),1);T.InsideCompartment=strings(height(T),1);
if isfield(model,'circuitCompartment')
    modelComp=model.circuitCompartment;
else
    [~,modelComp]=splitCompartments(string(model.mets));
end
for k=1:height(T)
    c=cores(D.CoreIndex(k));
    T.AllCoreECResolved(k)=all(strlength(strtrim(model.circuitResolvedEC(O.rxn(c.path))))>0);
    T.InsideCompartment(k)=modelComp(c.seed);
    if D.TransportRow(k)>0
        tr=transport(D.TransportRow(k),:);
        T.TransportHasGPR(k)=strlength(strtrim(tr.GPR))>0;
        if T.TransportHasGPR(k) && T.AllCoreECResolved(k)
            T.AnnotationStatus(k)="EC_and_GPR_present_gating_and_kinetics_unverified";
        else
            T.AnnotationStatus(k)="structural_only_missing_EC_or_GPR";
        end
    else
        T.TransportHasGPR(k)=false;
        if T.AllCoreECResolved(k)
            T.AnnotationStatus(k)="core_EC_present_no_I_import_searched";
        else
            T.AnnotationStatus(k)="structural_only_missing_core_EC";
        end
    end
end
end

function [SG,stats] = findSGatedCoreImports(cores,O,T,model,currency,cfg)
% S-GATED CO-SUBSTRATE IMPORT for the autocatalytic production core.
% Mirrors the F_gated_cosubstrate_import branch of findRemovalTerms.
% S gates a transporter importing a FEED metabolite X of the core. X then
% meets I inside that core, so the realised production flux carries the
% product S*I. Gating the I importer itself CANNOT do this: that flux is
% proportional to I[outside] and enters dI/dt as a source term, never as
% k*S*I[inside]. Transporter identity, gating and rate laws stay unverified.
rows=cell(0,16);
stats=struct('cores',numel(cores),'feedsTested',0,'importersFound',0, ...
    'coresWithSGate',0,'rejectedCoreReuse',0,'rejectedReversible',0);
for k=1:numel(cores)
    c=cores(k);found=false;
    for x=reshape(c.feed,1,[])
        if currency(x) || x==c.seed,continue;end
        stats.feedsTested=stats.feedsTested+1;
        for t=reshape(find(T.InsideIndex==x),1,[])
            if any(O.rxn(c.path)==T.ModelIndex(t))
                stats.rejectedCoreReuse=stats.rejectedCoreReuse+1;continue;
            end
            % A freely reversible carrier equilibrates X across the membrane
            % in BOTH directions, so the gate carries no sign. Reported
            % always; rejected only when requireSGateIrreversible is set.
            reversible=T.LowerBound(t)<0 && T.UpperBound(t)>0;
            if cfg.requireSGateIrreversible && reversible
                stats.rejectedReversible=stats.rejectedReversible+1;continue;
            end
            found=true;
            rows(end+1,:)={size(rows,1)+1,k,numel(c.path),c.seed, ...
                string(model.mets{c.seed}),x,string(model.mets{x}), ...
                full(-c.net(x)),t,T.ReactionID(t),T.Direction(t), ...
                T.OutsideMetaboliteID(t),T.Equation(t),T.GPR(t),reversible, ...
                "S_identity_and_gating_unverified"}; %#ok<AGROW>
        end
    end
    if found,stats.coresWithSGate=stats.coresWithSGate+1;end
end
SG=cellTable(rows,{'SGateID','CoreIndex','CoreReactionCount','I_Index','I_ID', ...
    'X_Index','X_ID','X_PerModule','STransportRow','STransportReactionID', ...
    'STransportDirection','X_Outside','STransportEquation','STransportGPR', ...
    'STransportReversible','SStatus'}, ...
    {'double','double','double','double','string','double','string','double','double', ...
    'string','double','string','string','string','logical','string'});
end

function [G,stats] = attachSGates(C,D,SG,cores,O,T,L,Q,model,currency,Rgross,Pgross,cfg)
% Cross every connected Anderson candidate with the S-gated co-substrate
% importers available for ITS core, reusing the reaction-reuse and
% isolation checks already applied to the F transporter. A connected row
% with no S-gated feed is counted, never silently promoted.
rows=cell(0,17);
stats=struct('tested',0,'reuseRejected',0,'crossTalkRejected',0, ...
    'candidatesWithoutSGate',0,'testLimitReached',false);
tol=1e-8;stop=false;
for i=1:height(C)
    sgs=find(SG.CoreIndex==C.CoreIndex(i));
    if isempty(sgs)
        stats.candidatesWithoutSGate=stats.candidatesWithoutSGate+1;continue;
    end
    core=cores(C.CoreIndex(i));seed=core.seed;
    iRow=D.TransportRow(C.ProductionCandidateID(i));
    prodRx=O.rxn(core.path);
    prodParts=find(any(O.R(:,core.path)+O.P(:,core.path)>tol,2));
    if iRow>0
        iImport=T.ModelIndex(iRow);prodRx=unique([prodRx,iImport]);
        prodParts=union(prodParts,find(Rgross(:,iImport)+Pgross(:,iImport)>tol));
    end
    l=find(L.LinearID==C.LinearID(i),1);q=find(Q.QuadraticID==C.QuadraticID(i),1);
    if isempty(l)||isempty(q),continue;end
    linRx=L.ReactionIndex(l);
    linParts=find(O.R(:,L.OrientationIndex(l))+O.P(:,L.OrientationIndex(l))>tol);
    if L.FTransportRow(l)>0
        ft=T.ModelIndex(L.FTransportRow(l));linRx=[linRx ft]; %#ok<AGROW>
        linParts=union(linParts,find(Rgross(:,ft)+Pgross(:,ft)>tol));
    end
    quadRx=Q.ReactionIndex(q);
    quadParts=find(O.R(:,Q.OrientationIndex(q))+O.P(:,Q.OrientationIndex(q))>tol);
    for sg=reshape(sgs,1,[])
        if stats.tested>=cfg.maxSGateTests
            stats.testLimitReached=true;stop=true;break;
        end
        stats.tested=stats.tested+1;
        st=T.ModelIndex(SG.STransportRow(sg));
        if ismember(st,[prodRx linRx quadRx])
            stats.reuseRejected=stats.reuseRejected+1;continue;
        end
        sParts=find(Rgross(:,st)+Pgross(:,st)>tol);
        % I and X are shared with the core BY DESIGN: that shared X is
        % precisely what makes the realised rate carry S*I. Everything
        % else shared with the removal terms is cross-talk.
        exempt=[seed;SG.X_Index(sg)];
        shared=unique([setdiff(intersect(sParts,prodParts),exempt); ...
                       setdiff(intersect(sParts,union(linParts,quadParts)),exempt)]);
        shared(currency(shared))=[];
        isolated=isempty(shared);
        if cfg.requireSGateIsolation && ~isolated
            stats.crossTalkRejected=stats.crossTalkRejected+1;continue;
        end
        rows(end+1,:)={size(rows,1)+1,C.ConnectedID(i),SG.SGateID(sg),C.CoreIndex(i), ...
            C.CoreReactionCount(i),C.I_ID(i),SG.X_ID(sg),SG.X_PerModule(sg), ...
            SG.STransportReactionID(sg),SG.STransportReversible(sg), ...
            C.FMode(i),C.FTransportReactionID(i),C.FControlCargo(i), ...
            C.PairClass(i),isolated,joinText(string(model.mets(shared)),';'), ...
            "S_gates_X_import_F_gates_L_import_rate_laws_unverified"}; %#ok<AGROW>
    end
    if stop,break;end
end
G=cellTable(rows,{'GatedID','ConnectedID','SGateID','CoreIndex','CoreReactionCount', ...
    'I_ID','S_Gated_X_ID','X_PerModule','S_TransportReactionID','S_TransportReversible', ...
    'FMode','F_TransportReactionID','F_ControlCargo','PairClass', ...
    'NoSharedNoncurrencyOutsideIandX','SharedNoncurrency','ArchitectureStatus'}, ...
    {'double','double','double','double','double','string','string','double','string', ...
    'logical','string','string','string','string','logical','string','string'});
if ~isempty(G)
    [~,order]=sortrows([G.CoreReactionCount, ...
        double(G.S_TransportReversible),G.GatedID],[1 2 3]);
    G=G(order,:);G.GatedID=(1:height(G))';
end
end

function writeSGateReport(path,SG,G,sstats,gstats,cfg)
fid=fopen(path,'w');assert(fid>0,'Anderson:Report','Cannot write S-gate report.');
guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
report(fid,'ANDERSON: S-GATED CO-SUBSTRATE IMPORT ON THE AUTOCATALYTIC CORE');
report(fid,sprintf('Started: %s.',datestr(now,31)));
report(fid,sprintf('Cores examined: %d; noncurrency feeds tested: %d; importers found: %d.', ...
    sstats.cores,sstats.feedsTested,sstats.importersFound));
report(fid,sprintf('Cores with at least one S-gated feed: %d; core-reuse rejections: %d.', ...
    sstats.coresWithSGate,sstats.rejectedCoreReuse));
report(fid,sprintf('Reversible-carrier rejections: %d (requireSGateIrreversible=%d).', ...
    sstats.rejectedReversible,cfg.requireSGateIrreversible));
report(fid,sprintf('Gated candidates: %d; crossings tested: %d.',height(G),gstats.tested));
report(fid,sprintf('Reaction-reuse rejections: %d; cross-talk rejections: %d (isolation required=%d).', ...
    gstats.reuseRejected,gstats.crossTalkRejected,cfg.requireSGateIsolation));
report(fid,sprintf('Connected rows with no S-gated feed: %d; test limit reached=%d.', ...
    gstats.candidatesWithoutSGate,gstats.testLimitReached));
if height(SG)>0
    report(fid,sprintf('Reversible S-carrier rows retained: %d of %d.', ...
        nnz(SG.STransportReversible),height(SG)));
end
report(fid,'S gates the X importer. A gated I importer yields a SOURCE term, not k*S*I.');
report(fid,'The I importer found in stage 03 is retained ONLY as the basal seed eps.');
report(fid,'A reversible S carrier equilibrates X both ways, so its gate carries no sign.');
report(fid,'Shared X between importer and core is intended; it is exempt from isolation.');
report(fid,'Gate identity, membrane control by S and all rate laws remain UNVERIFIED here.');
report(fid,'Quasi-steady-state collapse of the core to a single I state is NOT checked here.');
end

function [L,Q,C,A] = findRemovalTerms(signals,O,T,model,comp,currency,cfg)
% These are stoichiometric component candidates, NOT fitted loss-rate laws.
linearRows=cell(0,21);quadRows=cell(0,17);coverageRows=cell(0,8);auditRows=cell(0,9);
tol=1e-8;
for seed=reshape(signals,1,[])
    lstart=size(linearRows,1);qstart=size(quadRows,1);lossCount=0;oversized=0;unresolved=0;
    for o=reshape(O.consumers{seed},1,[])
        if O.N(seed,o)>=-tol,continue;end
        lossCount=lossCount+1;j=O.rxn(o);a=O.R(:,o);b=O.P(:,o);
        participants=find(a+b>tol);status="not_single_copy_or_pair_loss";
        if numel(participants)>cfg.maxRemovalParticipants
            oversized=oversized+1;
            auditRows(end+1,:)={string(model.mets{seed}),string(model.rxns{j}), ...
                O.direction(o),numel(participants),"skipped_participant_limit",NaN,NaN,NaN,"not_searched"}; %#ok<AGROW>
            continue;
        end
        [scale,ok]=primitiveReactionScale(a,b);
        if ~ok
            unresolved=unresolved+1;
            auditRows(end+1,:)={string(model.mets{seed}),string(model.rxns{j}), ...
                O.direction(o),numel(participants),"normalization_unresolved",NaN,NaN,NaN,"not_searched"}; %#ok<AGROW>
            continue;
        end
        ar=a*scale;br=b*scale;ai=full(ar(seed));bi=full(br(seed));loss=ai-bi;
        extra=find((br-ar)<-tol);extra(extra==seed)=[];
        products=find(br>tol);
        common={seed,string(model.mets{seed}),comp(seed),j,string(model.rxns{j}), ...
            O.direction(o),o,scale,loss,equation(a,b,model),equation(ar,br,model), ...
            joinText(string(model.mets(extra)),';'),joinText(string(model.mets(products)),';'),O.grossKnown(o)};
        if abs(ai-1)<tol && abs(bi)<tol && any(br>tol)
            status="single_copy_sink_kinetics_unverified";
            % Always retain a partial sink even when F cannot yet be linked.
            linearRows(end+1,:)=[{size(linearRows,1)+1},common, ...
                {"F_unassigned","",0,"","","F_control_missing"}]; %#ok<AGROW>
            % F can be an explicit regenerated noncurrency participant, but
            % only gross reactant/product evidence can reveal that role.
            if O.grossKnown(o)
                catalysts=find(ar>tol & br>tol & abs(br-ar)<tol & ~currency);
                catalysts(catalysts==seed)=[];
                for f=reshape(catalysts,1,[])
                    linearRows(end+1,:)=[{size(linearRows,1)+1},common, ...
                        {"explicit_regenerated_F",string(model.mets{f}),0,"","", ...
                        "gross_F_regeneration_present_kinetics_unverified"}]; %#ok<AGROW>
                end
            end
            % Preferred membrane architecture: F gates a transporter importing
            % co-substrate L; L is consumed, while the controller F is NOT L.
            cargos=extra(~currency(extra));
            for cargo=reshape(cargos,1,[])
                for t=reshape(find(T.InsideIndex==cargo),1,[])
                    if T.ModelIndex(t)==j,continue;end
                    linearRows(end+1,:)=[{size(linearRows,1)+1},common, ...
                        {"F_gated_cosubstrate_import","",t,T.ReactionID(t), ...
                        string(model.mets{cargo}),"F_identity_and_gating_unverified"}]; %#ok<AGROW>
                end
            end
        end
        if abs(ai-2)<tol && (abs(bi)<tol || abs(bi-1)<tol) && any(br>tol)
            exact=abs(bi-1)<tol && O.grossKnown(o);
            if exact
                pairClass="exact_gross_2I_to_I";
            elseif abs(bi)<tol
                pairClass="alternative_2I_to_products_two_copy_loss";
            else
                pairClass="gross_evidence_missing";
            end
            if pairClass~="gross_evidence_missing"
                quadRows(end+1,:)=[{size(quadRows,1)+1},common,{pairClass,exact}]; %#ok<AGROW>
                status=pairClass+"_kinetic_order_unverified";
            end
        end
        auditRows(end+1,:)={string(model.mets{seed}),string(model.rxns{j}), ...
            O.direction(o),numel(participants),"primitive_normalized",scale,ai,bi,status}; %#ok<AGROW>
    end
    theseL=linearRows(lstart+1:end,:);theseQ=quadRows(qstart+1:end,:);
    linkedF=0;exactQ=0;
    if ~isempty(theseL),modes=vertcat(theseL{:,16});linkedF=nnz(modes~="F_unassigned");end
    if ~isempty(theseQ),exactQ=nnz(cell2mat(theseQ(:,17)));end
    coverageRows(end+1,:)={string(model.mets{seed}),lossCount,size(theseL,1),linkedF, ...
        size(theseQ,1),exactQ,oversized,unresolved}; %#ok<AGROW>
end
commonNames={'I_Index','I_ID','InsideCompartment','ReactionIndex','ReactionID','Direction', ...
    'OrientationIndex','PrimitiveScale','LossPerEvent','OriginalEquation','PrimitiveEquation', ...
    'OtherNetConsumed','Products','GrossKnown'};
commonTypes={'double','string','string','double','string','double','double','double', ...
    'double','string','string','string','string','logical'};
L=cellTable(linearRows,[{'LinearID'},commonNames,{'FMode','FSpeciesID','FTransportRow', ...
    'FTransportReactionID','FControlCargo','FStatus'}], ...
    [{'double'},commonTypes,{'string','string','double','string','string','string'}]);
Q=cellTable(quadRows,[{'QuadraticID'},commonNames,{'PairClass','ExactFormalStoichiometry'}], ...
    [{'double'},commonTypes,{'string','logical'}]);
C=cellTable(coverageRows,{'I_ID','ConsumingOrientations','LinearRows','FLinkedRows', ...
    'PairLossRows','ExactGross2IToI','SkippedParticipantLimit','NormalizationUnresolved'}, ...
    {'string','double','double','double','double','double','double','double'});
A=cellTable(auditRows,{'I_ID','ReactionID','Direction','ParticipantCount', ...
    'NormalizationStatus','PrimitiveScale','PrimitiveIConsumed','PrimitiveIProduced','Status'}, ...
    {'string','string','double','double','string','double','double','double','string'});
end

function [scale,ok] = primitiveReactionScale(a,b)
% Remove arbitrary whole-column scaling WITHOUT cancelling gross catalysts.
% E.g. 2I+2X->2P is primitive I+X->P, not evidence of an I-pair event.
values=full(nonzeros([a;b]));scale=NaN;ok=false;
if isempty(values)||any(~isfinite(values)|values<=0),return;end
[num,den]=rat(values,1e-9);common=1;
for k=1:numel(den)
    common=lcm(common,den(k));
    if ~isfinite(common)||common>1e6,return;end
end
integers=round(num.*(common./den));divisor=integers(1);
for k=2:numel(integers),divisor=gcd(divisor,integers(k));end
if divisor<=0 || any(abs(integers)>1e12),return;end
scale=common/divisor;
ok=all(abs(values*scale-round(values*scale))<1e-7);
end

function [C,stats] = connectAndersonTerms(P,D,cores,O,T,L,Q,model,comp,currency,Rgross,Pgross,cfg)
rows=cell(0,18);stats=struct('tested',0,'reactionReuseRejected',0, ...
    'crossTalkRejected',0,'testLimitReached',false,'candidateLimitReached',false);
stop=false;tol=1e-8;
for p=1:height(P)
    core=cores(D.CoreIndex(p));seed=core.seed;iRow=D.TransportRow(p);
    prodRx=O.rxn(core.path);
    prodParts=find(any(O.R(:,core.path)+O.P(:,core.path)>tol,2));
    if iRow>0
        s=T.ModelIndex(iRow);prodRx=unique([prodRx,s]);
        prodParts=union(prodParts,find(Rgross(:,s)+Pgross(:,s)>tol));
    end
    ls=find(L.I_Index==seed);
    if cfg.requireFControlTopology,ls=ls(L.FMode(ls)~="F_unassigned");end
    qs=find(Q.I_Index==seed);
    if ~cfg.allowPairLossAlternative,qs=qs(Q.ExactFormalStoichiometry(qs));end
    for l=reshape(ls,1,[])
        linRx=L.ReactionIndex(l);linParts=find(O.R(:,L.OrientationIndex(l))+O.P(:,L.OrientationIndex(l))>tol);
        if L.FTransportRow(l)>0
            ft=T.ModelIndex(L.FTransportRow(l));linRx=[linRx ft]; %#ok<AGROW>
            linParts=union(linParts,find(Rgross(:,ft)+Pgross(:,ft)>tol));
        end
        for q=reshape(qs,1,[])
            if stats.tested>=cfg.maxTripletTests
                stats.testLimitReached=true;stop=true;break;
            end
            stats.tested=stats.tested+1;quadRx=Q.ReactionIndex(q);
            if any(ismember(linRx,prodRx)) || ismember(quadRx,[prodRx linRx]) ...
                    || numel(unique(linRx))<numel(linRx)
                stats.reactionReuseRejected=stats.reactionReuseRejected+1;continue;
            end
            quadParts=find(O.R(:,Q.OrientationIndex(q))+O.P(:,Q.OrientationIndex(q))>tol);
            shared=unique([intersect(prodParts,linParts);intersect(prodParts,quadParts); ...
                intersect(linParts,quadParts)]);
            shared(shared==seed | currency(shared))=[];
            isolated=isempty(shared);
            if cfg.requireRemovalIsolation && ~isolated
                stats.crossTalkRejected=stats.crossTalkRejected+1;continue;
            end
            rows(end+1,:)={size(rows,1)+1,p,D.CoreIndex(p),numel(core.path),seed, ...
                string(model.mets{seed}),comp(seed),L.LinearID(l),Q.QuadraticID(q), ...
                L.FMode(l),L.FTransportReactionID(l),L.FControlCargo(l),Q.PairClass(q), ...
                Q.LossPerEvent(q),Q.ExactFormalStoichiometry(q),isolated, ...
                joinText(string(model.mets(shared)),';'), ...
                "S_F_control_and_rate_laws_unverified"}; %#ok<AGROW>
            if size(rows,1)>=cfg.maxConnectedCandidates
                stats.candidateLimitReached=true;stop=true;break;
            end
        end
        if stop,break;end
    end
    if stop,break;end
end
C=cellTable(rows,{'ConnectedID','ProductionCandidateID','CoreIndex','CoreReactionCount', ...
    'I_Index','I_ID','InsideCompartment','LinearID','QuadraticID','FMode', ...
    'FTransportReactionID','FControlCargo','PairClass','QuadraticLossPerEvent', ...
    'ExactFormalQuadratic','NoSharedNoncurrencyOutsideI','SharedNoncurrencyOutsideI','RateStatus'}, ...
    {'double','double','double','double','double','string','string','double','double', ...
    'string','string','string','string','double','logical','logical','string','string'});
if ~isempty(C)
    [~,order]=sortrows([C.CoreReactionCount,-double(C.ExactFormalQuadratic),C.ConnectedID],[1 2 3]);
    C=C(order,:);C.ConnectedID=(1:height(C))';
end
end

function writeThreeReactionReport(path,P,D,cores,O,T,L,Q,C,stats,model,cfg)
fid=fopen(path,'w');assert(fid>0,'Anderson:Report','Cannot write three-reaction report.');
guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
report(fid,'ANDERSON: CONNECTED PRODUCTION, F-CONTROLLED LOSS, AND PAIR-LOSS CANDIDATES');
report(fid,sprintf('Started: %s; removal stage enabled=%d.',datestr(now,31),cfg.connectRemovalTerms));
report(fid,sprintf('Production combinations: %d; linear component rows: %d; F-linked rows: %d.', ...
    height(P),height(L),nnz(L.FMode~="F_unassigned")));
report(fid,sprintf('Pair-loss component rows: %d; exact gross 2I->I rows: %d.', ...
    height(Q),nnz(Q.ExactFormalStoichiometry)));
report(fid,sprintf('Connected structural triples: %d; exact formal quadratic matches: %d.', ...
    height(C),nnz(C.ExactFormalQuadratic)));
report(fid,sprintf('Triples tested: %d; reaction-reuse rejections: %d; cross-talk rejections: %d.', ...
    stats.tested,stats.reactionReuseRejected,stats.crossTalkRejected));
report(fid,sprintf('Test limit reached=%d; candidate limit reached=%d; isolation required=%d.', ...
    stats.testLimitReached,stats.candidateLimitReached,cfg.requireRemovalIsolation));
report(fid,'All three roles use the SAME model species index; e/u/x pools are never identified.');
report(fid,'F-gated co-substrate import is a topology hypothesis: F itself is NOT that consumed cargo.');
report(fid,'A linear sink lacking a linked F mechanism remains a partial component.');
report(fid,'No kinetic orders are inferred from stoichiometric coefficients or normalized column scales.');
report(fid,'Exact 2I->I events lose one copy; alternative 2I->products events lose TWO copies.');
report(fid,'For loss v_q=a*I^2, dI/dt includes -nu*a*I^2. Matching -k*I^2 needs a=k/nu.');
report(fid,'Enzymatic pair consumption alone does NOT demonstrate v_q proportional to I^2.');
report(fid,'No example is certified as a complete Anderson kinetic implementation.');
report(fid,'Reaction directions follow model bounds; thermodynamic and whole-cell flux feasibility are unverified.');
report(fid,'The current transport-fed production design still includes an explicit import flux.');
report(fid,'A flux-ratio, membrane-opening law, feed buffering and equal effective coefficients need validation.');
report(fid,'When F gates co-substrate L, a rate proportional to F*I requires an appropriate transport/kinetic reduction.');
report(fid,'Components above the participant limit and unresolved rational normalizations are in the filter audit.');
report(fid,'Production depth/state limits and pair/single-reaction removal scope can hide alternatives.');
if isempty(C)
    report(fid,'No connected structural triple passed this configuration; inspect partial component CSVs.');
end
for n=1:min(height(C),cfg.printTop)
    row=C(n,:);p=row.ProductionCandidateID;core=cores(D.CoreIndex(p));
    l=L(row.LinearID,:);q=Q(row.QuadraticID,:);
    report(fid,sprintf('\nCONNECTED CANDIDATE %d | I=%s | block=%s | core size=%d', ...
        n,row.I_ID,row.InsideCompartment,numel(core.path)));
    if D.TransportRow(p)>0
        tr=T(D.TransportRow(p),:);
        report(fid,['  I import (not the S gate): ' char(tr.ReactionID) ' | ' char(tr.Equation)]);
    else
        report(fid,'  I import: none searched; I(0) assumed present in the vessel (eps = 0).');
    end
    report(fid,['  Production core net: ' char(P.CoreNetEquation(p))]);
    for k=1:numel(core.path)
        o=core.path(k);
        report(fid,sprintf('    core[%g] %s (%+d): %s',core.flux(k), ...
            model.rxns{O.rxn(o)},O.direction(o),equation(O.R(:,o),O.P(:,o),model)));
    end
    report(fid,['  F topology: ' char(l.FMode) ' | F species: ' char(l.FSpeciesID)]);
    if l.FTransportRow>0
        ft=T(l.FTransportRow,:);
        report(fid,['  F-controlled co-substrate import: ' char(ft.ReactionID) ' | ' char(ft.Equation)]);
        report(fid,['  Consumed co-substrate (not F): ' char(l.FControlCargo)]);
    end
    report(fid,['  Linear sink: ' char(l.ReactionID) ' | ' char(l.OriginalEquation)]);
    report(fid,['  F status: ' char(l.FStatus)]);
    report(fid,['  Pair-loss candidate: ' char(q.ReactionID) ' | ' char(q.OriginalEquation)]);
    report(fid,sprintf('  Pair class: %s | primitive scale=%g | loss/event=%g | exact 2I->I=%d', ...
        q.PairClass,q.PrimitiveScale,q.LossPerEvent,q.ExactFormalStoichiometry));
    report(fid,['  Pair primitive equation: ' char(q.PrimitiveEquation)]);
    report(fid,['  Additional pair-reaction substrates (buffering unverified): ' char(q.OtherNetConsumed)]);
    report(fid,['  Rates: ' char(row.RateStatus)]);
end
end

function [E,A] = resolveBranchEC(model,base,root,cfg)
% Annotation scope is ONLY the fast/slow rails. Never guess from a name.
n=numel(model.rxns);sourceID=strings(n,1);
if isfield(model,'circuitSourceBiGG'),sourceID=string(model.circuitSourceBiGG(:));end
E=struct('codes',strings(n,1),'source',strings(n,1),'reference',strings(n,1));
status=repmat("unresolved_excluded_from_branches",n,1);
for j=1:n
    codes=modelECTokens(model,j);resolvedFallback=false;
    if isempty(codes)&&isfield(model,'circuitResolvedEC')
        codes=completeECTokens(model.circuitResolvedEC(j));resolvedFallback=~isempty(codes);
    end
    if ~isempty(codes)
        E.codes(j)=joinText(codes,';');E.source(j)="model_EC_annotation";
        E.reference(j)="model reaction "+string(model.rxns{j});status(j)="model_annotated";
        if resolvedFallback
            if isfield(model,'circuitECSource'),E.source(j)=model.circuitECSource(j);end
            if isfield(model,'circuitECReference'),E.reference(j)=model.circuitECReference(j);end
            status(j)="existing_resolved_EC_evidence";
        end
    end
end
file=findNamed(root,cfg.branchECFile);
if ~isempty(file)
    supplied=readTextTable(file);
    assert(all(ismember({'SourceReactionID','EC','Reference'},supplied.Properties.VariableNames)), ...
        'Anderson:BranchECFile','Branch EC evidence needs SourceReactionID, EC, Reference columns.');
    for k=1:height(supplied)
        code=completeECTokens(supplied.EC(k));
        if isempty(code)||supplied.Reference(k)=="",continue;end
        matches=find(sourceID==supplied.SourceReactionID(k) & sourceID~="");
        if ismember('ReactionID',supplied.Properties.VariableNames) && supplied.ReactionID(k)~=""
            matches=find(string(model.rxns)==supplied.ReactionID(k));
        end
        for j=reshape(matches,1,[])
            old=completeECTokens(E.codes(j));E.codes(j)=joinText(unique([old;code],'stable'),';');
            E.source(j)="supplied_branch_EC_evidence";
            refs=[E.reference(j);supplied.Reference(k)];refs=refs(refs~="");
            E.reference(j)=joinText(unique(refs,'stable'),'; ');status(j)="supplied_reference_not_independently_checked";
        end
    end
end
if cfg.branchUseVerifiedSeedMap
    seed=verifiedBranchECSeeds();
    for j=reshape(find(E.codes=="" & sourceID~=""),1,[])
        k=cellFreeSeedMatch(seed.SourceReactionID,sourceID(j),model,j);
        if isempty(k),status(j)="source_ID_outside_partial_EC_seed_map";continue;end
        if ~branchSignatureMatches(model.S(:,j),base,seed.ChemicalSignature(k))
            status(j)="explicit_BiGG_link_but_seed_chemistry_mismatch";continue;
        end
        E.codes(j)=seed.EC(k);E.source(j)="verified_partial_BiGG_seed_map_2026_11_13";
        E.reference(j)=seed.Reference(k);status(j)="source_ID_and_chemical_signature_matched";
    end
    missing=E.codes=="" & sourceID=="";
    status(missing)="missing_explicit_source_ID_for_EC_seed_lookup";
end
A=table(string(model.rxns),sourceID,E.codes,E.source,E.reference,status,E.codes~="", ...
    'VariableNames',{'ReactionID','SourceReactionID','EC','AnnotationSource','Reference','Status','BranchECEligible'});
template=fullfile(cfg.outputDir,'Branch_EC_Evidence_TEMPLATE.tsv');
if ~exist(template,'file')
    blank=table(strings(0,1),strings(0,1),strings(0,1), ...
        'VariableNames',{'SourceReactionID','EC','Reference'});
    writetable(blank,template,'FileType','text','Delimiter','\t');
end
end

function T = verifiedBranchECSeeds()
% Primary-source pages checked 2026-11-13. This is deliberately a PARTIAL
% starting set. G6PDH2 is explicitly listed among G6PDH2r's old identifiers.
% The full reaction signature is checked before transferring any annotation.
rows={ ...
    "PGI","5.3.1.9","g6p:-1;f6p:1","https://bigg.ucsd.edu/universal/reactions/PGI"; ...
    "G6PDH2","1.1.1.49;1.1.1.363","g6p:-1;nadp:-1;6pgl:1;h:1;nadph:1", ...
        "https://bigg.ucsd.edu/universal/reactions/G6PDH2r"; ...
    "G6PDH2r","1.1.1.49;1.1.1.363","g6p:-1;nadp:-1;6pgl:1;h:1;nadph:1", ...
        "https://bigg.ucsd.edu/universal/reactions/G6PDH2r"; ...
    "PGL","3.1.1.31","6pgl:-1;h2o:-1;6pgc:1;h:1","https://bigg.ucsd.edu/universal/reactions/PGL"; ...
    "PFK","2.7.1.11","atp:-1;f6p:-1;adp:1;fdp:1;h:1","https://bigg.ucsd.edu/universal/reactions/PFK"; ...
    "FBA","4.1.2.13","fdp:-1;dhap:1;g3p:1","https://bigg.ucsd.edu/universal/reactions/FBA"; ...
    "GAPD","1.2.1.12;1.2.1.59","g3p:-1;nad:-1;pi:-1;13dpg:1;h:1;nadh:1", ...
        "https://bigg.ucsd.edu/universal/reactions/GAPD"};
T=cellTable(rows,{'SourceReactionID','EC','ChemicalSignature','Reference'},repmat({'string'},1,4));
end

function ok = branchSignatureMatches(net,base,signature)
% Compare exact chemistry, permitting only a uniform orientation/scale change.
ok=false;parts=find(abs(net)>1e-10);names=base(parts);
if isempty(parts)||numel(unique(names))~=numel(names),return;end
terms=split(signature,';');expected=strings(numel(terms),1);coeff=zeros(numel(terms),1);
for k=1:numel(terms)
    pair=split(terms(k),':');if numel(pair)~=2,return;end
    expected(k)=pair(1);coeff(k)=str2double(pair(2));
end
if numel(parts)~=numel(expected)||any(~isfinite(coeff)|coeff==0),return;end
[found,index]=ismember(names,expected);if ~all(found),return;end
actual=full(net(parts));ratios=actual./coeff(index);
ok=abs(ratios(1))>1e-12&&max(abs(ratios-ratios(1)))<1e-8;
end

function [rails,M,A,stats] = findBranchModules(O,E,model,comp,currency,cfg)
rails=struct('input',{},'mid',{},'output',{},'o',{},'rx',{},'parts',{});
moduleRows=cell(0,20);assignmentRows=cell(0,20);tol=1e-8;
stats=struct('eligibleOrientations',0,'rails',0,'railLimitReached',false, ...
    'pairsTested',0,'pairLimitReached',false,'moduleLimitReached',false, ...
    'reactionReuseRejected',0,'crossTalkRejected',0,'ECRejected',0,'knownGeneOverlapRejected',0);
eligible=false(1,numel(O.rxn));scale=nan(1,numel(O.rxn));
edgeO=[];edgeFrom=[];edgeTo=[];outEdges=cell(numel(model.mets),1);
if cfg.connectBranches
    for o=1:numel(O.rxn)
        if E.codes(O.rxn(o))==""||nnz(O.R(:,o)+O.P(:,o))>cfg.branchMaxParticipants,continue;end
        [scale(o),ok]=primitiveReactionScale(O.R(:,o),O.P(:,o));if ~ok,continue;end
        eligible(o)=true;a=O.R(:,o)*scale(o);b=O.P(:,o)*scale(o);
        from=find(abs(a-1)<tol & b<tol & ~currency);
        to=find(abs(b-1)<tol & a<tol & ~currency);
        for u=reshape(from,1,[])
            for v=reshape(to,1,[])
                if u==v,continue;end
                edgeO(end+1)=o;edgeFrom(end+1)=u;edgeTo(end+1)=v; %#ok<AGROW>
                outEdges{u}(end+1)=numel(edgeO);
            end
        end
    end
end
stats.eligibleOrientations=nnz(eligible);stop=false;
for e1=1:numel(edgeO)
    u=edgeFrom(e1);mid=edgeTo(e1);o1=edgeO(e1);
    if string(cfg.branchInputID)~=""&&string(model.mets{u})~=string(cfg.branchInputID),continue;end
    for e2=reshape(outEdges{mid},1,[])
        out=edgeTo(e2);o2=edgeO(e2);
        if out==u||O.rxn(o1)==O.rxn(o2),continue;end
        p1=find(O.R(:,o1)+O.P(:,o1)>tol);p2=find(O.R(:,o2)+O.P(:,o2)>tol);
        shared=intersect(p1,p2);shared(shared==mid|currency(shared))=[];
        if ~isempty(shared)||ismember(out,p1)||ismember(u,p2),continue;end
        for o3=reshape(O.consumers{out},1,[])
            if ~eligible(o3)||ismember(O.rxn(o3),O.rxn([o1 o2])),continue;end
            ar=O.R(:,o3)*scale(o3);br=O.P(:,o3)*scale(o3);
            if abs(full(ar(out))-1)>tol||br(out)>tol||~any(br>tol),continue;end
            p3=find(ar+br>tol);
            s13=intersect(p1,p3);s13(currency(s13))=[];
            s23=intersect(p2,p3);s23(s23==out|currency(s23))=[];
            if ~isempty(s13)||~isempty(s23),continue;end
            rails(end+1)=struct('input',u,'mid',mid,'output',out,'o',[o1 o2 o3], ...
                'rx',O.rxn([o1 o2 o3]),'parts',unique([p1;p2;p3])); %#ok<AGROW>
            if numel(rails)>=cfg.branchMaxRails,stats.railLimitReached=true;stop=true;break;end
        end
        if stop,break;end
    end
    if stop,break;end
end
stats.rails=numel(rails);stop=false;
for f=1:numel(rails)
    for s=f+1:numel(rails)
        if rails(f).input~=rails(s).input,continue;end
        if stats.pairsTested>=cfg.branchMaxPairTests,stats.pairLimitReached=true;stop=true;break;end
        stats.pairsTested=stats.pairsTested+1;u=rails(f).input;
        rx=[rails(f).rx,rails(s).rx];
        if numel(unique(rx))~=6,stats.reactionReuseRejected=stats.reactionReuseRejected+1;continue;end
        shared=intersect(rails(f).parts,rails(s).parts);shared(shared==u|currency(shared))=[];
        if ~isempty(shared),stats.crossTalkRejected=stats.crossTalkRejected+1;continue;end
        sets=cell(1,6);geneSets=cell(1,6);allGPR=true;
        for k=1:6
            sets{k}=completeECTokens(E.codes(rx(k)));
            [gpr,genes,~]=decodeGPR(string(model.grRules{rx(k)}),string(model.rules{rx(k)}),string(model.genes));
            geneSets{k}=split(genes,';');geneSets{k}=geneSets{k}(geneSets{k}~="");
            allGPR=allGPR&&strlength(gpr)>0&&~isempty(geneSets{k});
        end
        [assigned,ecOK]=distinctBranchECs(sets,cfg.branchRequireDistinctECs);
        if ~ecOK,stats.ECRejected=stats.ECRejected+1;continue;end
        overlap=strings(0,1);
        for i=1:6,for j=i+1:6,overlap=union(overlap,intersect(geneSets{i},geneSets{j}));end;end
        if cfg.branchRejectSharedKnownGenes&&~isempty(overlap)
            stats.knownGeneOverlapRejected=stats.knownGeneOverlapRejected+1;continue;
        end
        % Both assignments are designs: EC labels do not measure which rail
        % is faster. The user will choose enzyme levels after kinetic lookup.
        for swap=0:1
            fi=f;si=s;ec=assigned;
            if swap,fi=s;si=f;ec=assigned([4 5 6 1 2 3]);end
            id=size(moduleRows,1)+1;fr=rails(fi);sr=rails(si);
            moduleRows(end+1,:)={id,u,string(model.mets{u}),fr.mid,string(model.mets{fr.mid}), ...
                fr.output,string(model.mets{fr.output}),sr.mid,string(model.mets{sr.mid}), ...
                sr.output,string(model.mets{sr.output}),fi,si,comp(u),joinText(ec,';'), ...
                true,numel(unique(ec))==6,allGPR,joinText(overlap,';'), ...
                "EC_qualified_enzyme_independence_and_rates_unverified"}; %#ok<AGROW>
            ro=[fr.o sr.o];input=[u fr.mid fr.output u sr.mid sr.output];
            output=[fr.mid fr.output 0 sr.mid sr.output 0];
            roles=["r_Ff","r_Fd","r_Fr","r_Sf","r_Sd","r_Sr"];
            symbols=["k_Ff","k_Fd","2*k_Fr","k_Sf","k_Sd","k_Sr"];
            for k=1:6
                o=ro(k);j=O.rxn(o);outID="";if output(k)>0,outID=string(model.mets{output(k)});end
                assignmentRows(end+1,:)={id,roles(k),j,string(model.rxns{j}),O.direction(o), ...
                    branchSourceID(model,j),ec(k),E.codes(j),E.source(j),E.reference(j), ...
                    string(model.mets{input(k)}),outID,equation(O.R(:,o),O.P(:,o),model), ...
                    full(-O.N(input(k),o)),symbols(k),NaN,NaN,NaN,NaN, ...
                    "needs_substrate_specific_KM_kcat_and_linear_regime"}; %#ok<AGROW>
            end
            if size(moduleRows,1)>=cfg.branchMaxModules,stats.moduleLimitReached=true;stop=true;break;end
        end
        if stop,break;end
    end
    if stop,break;end
end
M=cellTable(moduleRows,{'ModuleID','InputIndex','InputID','FastMidIndex','FastMidID', ...
    'FastSignalIndex','FastSignalID','SlowMidIndex','SlowMidID','SlowSignalIndex','SlowSignalID', ...
    'FastRailIndex','SlowRailIndex','InsideCompartment','AssignedECs','AllBranchECsPresent', ...
    'DistinctECAssignment','AllGPRPresent','KnownGeneOverlap','ControlStatus'}, ...
    {'double','double','string','double','string','double','string','double','string', ...
    'double','string','double','double','string','string','logical','logical','logical','string','string'});
A=cellTable(assignmentRows,{'ModuleID','Role','ReactionIndex','ReactionID','Direction','SourceReactionID', ...
    'AssignedEC','AllECs','AnnotationSource','Reference','InputSignalID','OutputSignalID','Equation', ...
    'SignalConsumedPerRawEvent','TargetRateSymbol','TargetRate_per_s','KM_mM','Kcat_per_s', ...
    'EnzymeConcentration_mM','RateStatus'}, ...
    {'double','string','double','string','double','string','string','string','string','string', ...
    'string','string','string','double','string','double','double','double','double','string'});
end

function id = branchSourceID(model,j)
id="";if isfield(model,'circuitSourceBiGG'),id=model.circuitSourceBiGG(j);end
end

function [assignment,ok] = distinctBranchECs(sets,requireDistinct)
assignment=strings(1,numel(sets));ok=all(~cellfun(@isempty,sets));if ~ok,return;end
if ~requireDistinct
    for k=1:numel(sets),assignment(k)=sets{k}(1);end
    return;
end
[~,order]=sort(cellfun(@numel,sets));
[assignment,ok]=assignBranchECRecursive(sets,order,1,assignment,strings(0,1));
end

function [assignment,ok] = assignBranchECRecursive(sets,order,k,assignment,used)
if k>numel(order),ok=true;return;end
index=order(k);ok=false;
for code=reshape(sets{index},1,[])
    if any(used==code),continue;end
    trial=assignment;trial(index)=code;
    [trial,valid]=assignBranchECRecursive(sets,order,k+1,trial,[used;code]);
    if valid,assignment=trial;ok=true;return;end
end
end

function [W,stats] = wireBranchesToCancellers(M,rails,C,P,D,cores,O,T,L,Q,model,currency,Rgross,Pgross,cfg)
rows=cell(0,14);tol=1e-8;stop=false;
stats=struct('tested',0,'reactionReuseRejected',0,'crossTalkRejected',0, ...
    'unmatchedInterfaceRejected',0,'testLimitReached',false,'wiringLimitReached',false);
for b=1:height(M)
    fast=rails(M.FastRailIndex(b));slow=rails(M.SlowRailIndex(b));
    branchRx=[fast.rx slow.rx];branchParts=union(fast.parts,slow.parts);
    for c=reshape(find(C.InsideCompartment==M.InsideCompartment(b)),1,[])
        if stats.tested>=cfg.branchMaxWiringTests,stats.testLimitReached=true;stop=true;break;end
        stats.tested=stats.tested+1;prod=C.ProductionCandidateID(c);core=cores(D.CoreIndex(prod));
        iRow=D.TransportRow(prod);l=L(C.LinearID(c),:);q=Q(C.QuadraticID(c),:);
        cancelRx=[O.rxn(core.path),l.ReactionIndex,q.ReactionIndex];
        cancelParts=find(any(O.R(:,core.path)+O.P(:,core.path)>tol,2));
        if iRow>0
            st=T.ModelIndex(iRow);cancelRx(end+1)=st; %#ok<AGROW>
            cancelParts=union(cancelParts,find(Rgross(:,st)+Pgross(:,st)>tol));
        end
        cancelParts=union(cancelParts,find(O.R(:,l.OrientationIndex)+O.P(:,l.OrientationIndex)>tol));
        cancelParts=union(cancelParts,find(O.R(:,q.OrientationIndex)+O.P(:,q.OrientationIndex)>tol));
        if l.FTransportRow>0
            ft=T.ModelIndex(l.FTransportRow);cancelRx(end+1)=ft; %#ok<AGROW>
            cancelParts=union(cancelParts,find(Rgross(:,ft)+Pgross(:,ft)>tol));
        end
        if any(ismember(branchRx,cancelRx))
            stats.reactionReuseRejected=stats.reactionReuseRejected+1;continue;
        end
        % Cargo I/formate is NEVER silently relabelled as gate controller S/F.
        sMatched=any(split(P.S_IDs(prod),';')==M.SlowSignalID(b));
        fMatched=l.FMode=="explicit_regenerated_F" && l.FSpeciesID==M.FastSignalID(b);
        both=sMatched&&fMatched;allowed=[];
        if sMatched&&slow.output~=C.I_Index(c),allowed(end+1)=slow.output;end %#ok<AGROW>
        if fMatched&&fast.output~=C.I_Index(c),allowed(end+1)=fast.output;end %#ok<AGROW>
        shared=intersect(branchParts,cancelParts);shared(currency(shared)|ismember(shared,allowed))=[];
        if ~isempty(shared),stats.crossTalkRejected=stats.crossTalkRejected+1;continue;end
        if cfg.branchRequireMatchedInterfaces&&~both
            stats.unmatchedInterfaceRejected=stats.unmatchedInterfaceRejected+1;continue;
        end
        status="provisional_SF_gate_interfaces_and_kinetics_unverified";
        if both,status="control_species_matched_rates_and_enzyme_independence_unverified";end
        rows(end+1,:)={size(rows,1)+1,M.ModuleID(b),C.ConnectedID(c),M.InputID(b), ...
            M.FastSignalID(b),M.SlowSignalID(b),M.InsideCompartment(b),sMatched,fMatched,both, ...
            joinText(string(model.mets(shared)),';'),true,true,status}; %#ok<AGROW>
        if size(rows,1)>=cfg.branchMaxWiring,stats.wiringLimitReached=true;stop=true;break;end
    end
    if stop,break;end
end
W=cellTable(rows,{'WiringID','ModuleID','ConnectedID','InputID','FastSignalID','SlowSignalID', ...
    'InsideCompartment','S_InterfaceMatched','F_InterfaceMatched','BothControlInterfacesMatched', ...
    'SharedNoncurrencyOutsideMatchedInterfaces','NoReactionReuse','NoncurrencyIsolationPassed','InterfaceStatus'}, ...
    {'double','double','double','string','string','string','string','logical','logical','logical', ...
    'string','logical','logical','string'});
end

function writeBranchReport(path,E,M,A,W,stats,wstats,sbmlInfo,cfg)
fid=fopen(path,'w');assert(fid>0,'Anderson:Report','Cannot write branch report.');
guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
report(fid,'FAST/SLOW ENZYMATIC RAILS AND ANDERSON CANCELLER INTERFACES');
report(fid,'Reader revision: native MATLAB SBML metadata; no Java required, 2026-11-13.');
report(fid,['Search file: ' mfilename('fullpath')]);
report(fid,sprintf('Started: %s; branch stage enabled=%d.',datestr(now,31),cfg.connectBranches));
report(fid,'Topology: m_in --r_Ff--> m_Ff --r_Fd--> m_Fd --r_Fr--> waste.');
report(fid,'Topology: m_in --r_Sf--> m_Sf --r_Sd--> m_Sd --r_Sr--> waste.');
report(fid,'Connection hypothesis: m_Sd controls the S gate; m_Fd controls the F gate.');
report(fid,'All six rail reactions require complete EC annotations. Cancellers retain structural EC-free screening.');
report(fid,sprintf('SBML metadata: %d reaction nodes; %d matched model reactions; %d explicit BiGG IDs.', ...
    sbmlInfo.parsedReactionCount,sbmlInfo.matchedReactionCount,sbmlInfo.sourceBiGGCount));
if sbmlInfo.error~="",report(fid,['SBML metadata error: ' char(sbmlInfo.error)]);end
if ~isempty(sbmlInfo.unmatchedReactionExamples)
    report(fid,['Unmatched raw reaction examples: ' char(joinText(sbmlInfo.unmatchedReactionExamples,'; '))]);
end
seed=verifiedBranchECSeeds();
report(fid,sprintf('Embedded EC map: %d source-ID matches; %d chemistry matches; %d chemistry mismatches.', ...
    nnz(ismember(E.SourceReactionID,seed.SourceReactionID)), ...
    nnz(E.Status=="source_ID_and_chemical_signature_matched"), ...
    nnz(E.Status=="explicit_BiGG_link_but_seed_chemistry_mismatch")));
if nnz(E.BranchECEligible)==0
    report(fid,'BRANCH SEARCH NOT STARTED: zero EC-qualified reactions. This is an annotation/metadata diagnostic, not absence of branches.');
    report(fid,'Inspect RawSBMLReactionID/SourceReactionID in 00_Gross_Stoichiometry_Audit.csv and Status in 07_Branch_EC_Audit.csv.');
end
report(fid,sprintf('EC-annotated model reactions: %d/%d; eligible rail orientations: %d.', ...
    nnz(E.BranchECEligible),height(E),stats.eligibleOrientations));
report(fid,'The embedded BiGG starting map covers SIX reaction classes, not the complete model.');
report(fid,'Only exact annotated source IDs with matching chemical signatures receive embedded EC annotations.');
report(fid,'Add sourced SourceReactionID/EC/Reference rows to Branch_EC_Evidence.tsv for broader coverage.');
report(fid,sprintf('Single rails: %d; pairs tested: %d; branch modules: %d; wiring rows: %d; matched interfaces: %d.', ...
    stats.rails,stats.pairsTested,height(M),height(W),nnz(W.BothControlInterfacesMatched)));
report(fid,sprintf('Pair rejections: reaction reuse=%d; shared noncurrency=%d; EC assignment=%d; shared known genes=%d.', ...
    stats.reactionReuseRejected,stats.crossTalkRejected,stats.ECRejected,stats.knownGeneOverlapRejected));
report(fid,sprintf('Limits reached: rails=%d; pairs=%d; modules=%d; wiring tests=%d; wiring rows=%d.', ...
    stats.railLimitReached,stats.pairLimitReached,stats.moduleLimitReached, ...
    wstats.testLimitReached,wstats.wiringLimitReached));
report(fid,sprintf('Wiring rejections: reused reactions=%d; shared noncurrency=%d; required interface mismatch=%d.', ...
    wstats.reactionReuseRejected,wstats.crossTalkRejected,wstats.unmatchedInterfaceRejected));
report(fid,'Fast/slow are assigned design roles. EC identifiers do not measure rates or prove independent proteins.');
report(fid,'Distinct EC assignments are a screening constraint; gene/protein identity and specificity still need validation.');
report(fid,'Missing GPRs remain unknown. Shared known branch genes are rejected under the default policy.');
report(fid,'Under an unsaturated single-signal approximation, k_eff=nu*kcat*[E]/KM.');
report(fid,'Thus [E]=k_eff*KM/(nu*kcat), with direction/substrate-specific kinetics and other substrates buffered.');
report(fid,'Assignment CSV: nu=SignalConsumedPerRawEvent. TargetRate, KM, kcat and [E] are deliberately NaN.');
report(fid,'The attached simulation uses 2*k_Fr for the fast endpoint loss; the target-rate symbols retain that factor.');
report(fid,'Its active lines J_I=alpha*m_Sd and J_Q=alpha*m_Fd feed an I+Q annihilator, not Anderson SI-FI-I^2.');
report(fid,'This network search preserves Anderson three-reaction candidates; no simulation ODE is silently replaced.');
report(fid,'No gate ligand is inferred from a transported or consumed cargo. Unmatched S/F interfaces are provisional.');
report(fid,'08_Matched_Control_Interfaces.csv contains only rows matching BOTH known control species.');
report(fid,'Even matched interfaces do not establish membrane opening, quadratic kinetics, speed ordering or thermal feasibility.');
for n=1:min(height(M),cfg.printTop)
    row=M(n,:);report(fid,sprintf('\nBRANCH MODULE %d | block=%s | input=%s', ...
        row.ModuleID,row.InsideCompartment,row.InputID));
    report(fid,['  Fast: ' char(row.FastMidID) ' -> ' char(row.FastSignalID)]);
    report(fid,['  Slow: ' char(row.SlowMidID) ' -> ' char(row.SlowSignalID)]);
    report(fid,sprintf('  Complete ECs=%d; distinct assignment=%d; all GPRs present=%d.', ...
        row.AllBranchECsPresent,row.DistinctECAssignment,row.AllGPRPresent));
    for k=reshape(find(A.ModuleID==row.ModuleID),1,[])
        report(fid,sprintf('    %s | %s (%+d) | EC %s | %s', ...
            A.Role(k),A.ReactionID(k),A.Direction(k),A.AssignedEC(k),A.Equation(k)));
        report(fid,['      EC source: ' char(A.Reference(k))]);
    end
end
end

function [E,A] = resolveCounterEC(model,base,branchEC,root,cfg)
% Separate annotation scope: the working six-class branch map is unchanged.
E=branchEC;n=numel(model.rxns);sourceID=strings(n,1);
if isfield(model,'circuitSourceBiGG'),sourceID=string(model.circuitSourceBiGG(:));end
status=repmat("unresolved_excluded_from_counter_removal",n,1);
status(E.codes~="")="inherited_existing_EC_evidence";
file=findNamed(root,cfg.counterECFile);
if ~isempty(file)
    supplied=readTextTable(file);
    assert(all(ismember({'SourceReactionID','EC','Reference'},supplied.Properties.VariableNames)), ...
        'Anderson:CounterECFile','Counter EC evidence needs SourceReactionID, EC, Reference.');
    for k=1:height(supplied)
        codes=completeECTokens(supplied.EC(k));
        if isempty(codes)||supplied.Reference(k)=="",continue;end
        matches=find(sourceID==supplied.SourceReactionID(k)&sourceID~="");
        if ismember('ReactionID',supplied.Properties.VariableNames)&&supplied.ReactionID(k)~=""
            matches=find(string(model.rxns)==supplied.ReactionID(k));
        end
        for j=reshape(matches,1,[])
            E.codes(j)=joinText(unique([completeECTokens(E.codes(j));codes],'stable'),';');
            E.source(j)="supplied_counter_EC_evidence";
            refs=[E.reference(j);supplied.Reference(k)];refs=refs(refs~="");
            E.reference(j)=joinText(unique(refs,'stable'),'; ');
            status(j)="supplied_reference_not_independently_checked";
        end
    end
end
if cfg.counterUseVerifiedSeedMap
    seeds=verifiedCounterECSeeds();
    for j=reshape(find(E.codes==""&sourceID~=""),1,[])
        k=cellFreeSeedMatch(seeds.SourceReactionID,sourceID(j),model,j);if isempty(k),continue;end
        if ~branchSignatureMatches(model.S(:,j),base,seeds.ChemicalSignature(k))
            status(j)="explicit_BiGG_link_but_counter_seed_chemistry_mismatch";continue;
        end
        E.codes(j)=seeds.EC(k);E.source(j)="verified_partial_counter_BiGG_seed_2026_11_13";
        E.reference(j)=seeds.Reference(k);status(j)="source_ID_and_chemical_signature_matched";
    end
end
if cfg.counterUseEmbeddediMM904EC
    [E,status]=applyEmbeddediMM904EC(E,status,model,base);
end
A=table(string(model.rxns),sourceID,E.codes,E.source,E.reference,status,E.codes~="", ...
    'VariableNames',{'ReactionID','SourceReactionID','EC','AnnotationSource','Reference','Status','CounterECEligible'});
end

function T = verifiedCounterECSeeds()
% Primary BiGG pages checked 2026-11-13. Complete signatures are required.
% SFGTH is listed as an old identifier on the linked SFGTHi model page.
% These annotations establish enzyme classes, not rates or thermodynamics.
rows={ ...
    "GALK","2.7.1.6","atp:-1;gal:-1;adp:1;gal1p:1;h:1","https://bigg.ucsd.edu/universal/reactions/GALKr"; ...
    "GALKr","2.7.1.6","atp:-1;gal:-1;adp:1;gal1p:1;h:1","https://bigg.ucsd.edu/universal/reactions/GALKr"; ...
    "GTMLT","2.3.2.2","ala_L:-1;gthrd:-1;cgly:1;gluala:1","https://bigg.ucsd.edu/universal/reactions/GTMLT"; ...
    "LGTHL","4.4.1.5","gthrd:-1;mthgxl:-1;lgt_S:1","https://bigg.ucsd.edu/universal/reactions/LGTHL"; ...
    "FALDH","1.1.1.284;1.2.1.1","fald:-1;gthrd:-1;nad:-1;Sfglutth:1;h:1;nadh:1", ...
        "https://bigg.ucsd.edu/universal/reactions/FALDH"; ...
    "SFGTH","3.1.2.12","Sfglutth:-1;h2o:-1;for:1;gthrd:1;h:1", ...
        "https://bigg.ucsd.edu/models/iCHOv1/reactions/SFGTHi"; ...
    "SFGTHi","3.1.2.12","Sfglutth:-1;h2o:-1;for:1;gthrd:1;h:1", ...
        "https://bigg.ucsd.edu/universal/reactions/SFGTHi"};
T=cellTable(rows,{'SourceReactionID','EC','ChemicalSignature','Reference'},repmat({'string'},1,4));
end

function [F,support,A,stats] = findConsumingCounterFormations(C,O,E,model,base,comp,currency,cfg)
% Catalog chemistry before EC/isolation filtering, including unresolved ECs.
rows=cell(0,17);support=cell(0,1);audit=cell(0,7);tol=1e-8;stop=false;
stats=struct('tested',0,'testLimitReached',false,'componentLimitReached',false);
signals=unique(C.I_Index,'stable');if ~cfg.connectCounter,signals=[];end
for seed=reshape(signals,1,[])
    for o=reshape(O.consumers{seed},1,[])
        if stats.tested>=cfg.counterMaxComponentTests,stats.testLimitReached=true;stop=true;break;end
        stats.tested=stats.tested+1;j=O.rxn(o);a=O.R(:,o);b=O.P(:,o);parts=find(a+b>tol);
        [scale,ok]=primitiveReactionScale(a,b);products=[];status="formation_chemistry_catalogued";
        if ~O.grossKnown(o),status="gross_stoichiometry_unknown_check_metabolite_ID_mapping";
        elseif numel(parts)>cfg.counterMaxParticipants,status="participant_limit";
        elseif ~ok||abs(full(a(seed))*scale-1)>tol||b(seed)>tol,status="not_primitive_single_I_consumption";
        elseif any(base(parts)==base(seed)&parts~=seed),status="another_pool_of_I_in_formation_chemistry";
        else
            products=find(b>tol&a<tol&~currency&comp==comp(seed)&base~=base(seed));
            if cfg.counterCargoID~=""
                products=products(string(model.mets(products))==cfg.counterCargoID|base(products)==cfg.counterCargoID);
            end
            if isempty(products),status="no_distinct_noncurrency_C_product";end
        end
        audit(end+1,:)={string(model.mets{seed}),string(model.rxns{j}),O.direction(o), ...
            O.grossKnown(o),E.codes(j)~="",numel(products),status}; %#ok<AGROW>
        for c=reshape(products,1,[])
            rows(end+1,:)={size(rows,1)+1,seed,string(model.mets{seed}),c,string(model.mets{c}), ...
                comp(seed),o,j,string(model.rxns{j}),O.direction(o),scale, ...
                full(a(seed)-b(seed))*scale,full(b(c)-a(c))*scale,E.codes(j),E.reference(j), ...
                equation(a,b,model),E.codes(j)~=""}; %#ok<AGROW>
            support{end+1,1}=parts; %#ok<AGROW>
            if size(rows,1)>=cfg.counterMaxComponents,stats.componentLimitReached=true;stop=true;break;end
        end
        if stop,break;end
    end
    if stop,break;end
end
F=cellTable(rows,{'FormationID','I_Index','I_ID','C_Index','C_ID','InsideCompartment', ...
    'OrientationIndex','ReactionIndex','ReactionID','Direction','PrimitiveScale', ...
    'I_LossPerEvent','C_GainPerEvent','ECs','ECReference','Equation','ECResolved'}, ...
    {'double','double','string','double','string','string','double','double','string','double', ...
    'double','double','double','string','string','string','logical'});
A=cellTable(audit,{'I_ID','ReactionID','Direction','GrossKnown','ECResolved','EligibleProductCount','Status'}, ...
    {'string','string','double','logical','logical','double','string'});
end

function [X,support,A,stats] = connectConsumingCounters(F,fp,C,D,cores,O,T,L,Q,E,model,base,currency,R,P,cfg)
rows=cell(0,16);support=cell(0,1);audit=cell(0,10);tol=1e-8;stop=false;
stats=struct('tested',0,'leakTests',0,'testLimitReached',false,'leakTestLimitReached',false, ...
    'connectionLimitReached',false);
for c=1:height(C)
    [cr,cp]=counterCancellerSupport(c,C,D,cores,O,T,L,Q,R,P);
    for f=reshape(find(F.I_Index==C.I_Index(c)),1,[])
        if stats.tested>=cfg.counterMaxConnectionTests,stats.testLimitReached=true;stop=true;break;end
        stats.tested=stats.tested+1;seed=F.I_Index(f);cargo=F.C_Index(f);fr=F.ReactionIndex(f);
        shared=intersect(fp{f},cp);shared(shared==seed|currency(shared))=[];
        isolated=isempty(shared);notReused=~ismember(fr,cr);
        ecOK=F.ECResolved(f)||~cfg.counterRequireFormationEC;status="formation_passed";
        if ~ecOK,status="formation_EC_unresolved";
        elseif ~notReused,status="formation_reuses_canceller_reaction";
        elseif cfg.counterRequireIsolation&&~isolated,status="formation_shares_noncurrency_outside_I";
        end
        audit(end+1,:)={C.ConnectedID(c),F.FormationID(f),"formation",F.ReactionID(f),F.C_ID(f), ...
            ecOK,notReused,joinText(string(model.mets(shared)),';'),isolated,status}; %#ok<AGROW>
        if status~="formation_passed",continue;end
        variants=reshape(O.consumers{cargo},1,[]);
        if ~cfg.counterRequireRemoval,variants=[0 variants];end
        for o=variants
            leakID="";leakEC="";direction=0;scale=NaN;parts=fp{f};allShared=shared;
            type="accumulator_without_selected_removal";
            if o>0
                if stats.leakTests>=cfg.counterMaxComponentTests
                    stats.leakTestLimitReached=true;break;
                end
                stats.leakTests=stats.leakTests+1;lr=O.rxn(o);a=O.R(:,o);b=O.P(:,o);
                lp=find(a+b>tol);[scale,ok]=primitiveReactionScale(a,b);
                directShared=intersect(fp{f},lp);directShared(directShared==cargo|currency(directShared))=[];
                cancelShared=intersect(lp,cp);cancelShared(currency(cancelShared))=[];
                allShared=union(shared,union(directShared,cancelShared));
                leakEC=E.codes(lr);ecOK=leakEC~=""||~cfg.counterRequireLeakEC;
                notReused=~ismember(lr,[cr fr]);status="removal_passed";
                if ~O.grossKnown(o),status="removal_gross_stoichiometry_unknown";
                elseif numel(lp)>cfg.counterMaxParticipants,status="removal_participant_limit";
                elseif ~ok||abs(full(a(cargo))*scale-1)>tol||b(cargo)>tol,status="not_primitive_single_C_removal";
                elseif any(base(lp)==base(seed)),status="removal_contains_or_regenerates_I";
                elseif ~ecOK,status="removal_EC_unresolved";
                elseif ~notReused,status="removal_reuses_formation_or_canceller_reaction";
                elseif cfg.counterRequireIsolation&&~isempty(allShared),status="removal_noncurrency_overlap";
                end
                audit(end+1,:)={C.ConnectedID(c),F.FormationID(f),"removal",string(model.rxns{lr}), ...
                    F.C_ID(f),ecOK,notReused,joinText(string(model.mets(allShared)),';'),isempty(allShared),status}; %#ok<AGROW>
                if status~="removal_passed",continue;end
                leakID=string(model.rxns{lr});direction=O.direction(o);parts=union(parts,lp);
                type="leaky_counter_with_selected_removal";
            end
            rows(end+1,:)={size(rows,1)+1,C.ConnectedID(c),F.FormationID(f),F.I_ID(f),F.C_ID(f), ...
                F.InsideCompartment(f),F.ECs(f),o,leakID,direction,leakEC,scale,type, ...
                isempty(allShared),joinText(string(model.mets(allShared)),';'), ...
                "counter_consumes_I_adds_canceller_loss_kinetics_unverified"}; %#ok<AGROW>
            support{end+1,1}=parts; %#ok<AGROW>
            if size(rows,1)>=cfg.counterMaxConnections,stats.connectionLimitReached=true;stop=true;break;end
        end
        if stop,break;end
    end
    if stop,break;end
end
X=cellTable(rows,{'CounterConnectionID','ConnectedID','FormationID','I_ID','C_ID','InsideCompartment', ...
    'FormationECs','LeakOrientationIndex','LeakReactionID','LeakDirection','LeakECs','LeakPrimitiveScale', ...
    'CounterType','NoncurrencyIsolationPassed','SharedNoncurrencyOutsideI','Status'}, ...
    {'double','double','double','string','string','string','string','double','string','double', ...
    'string','double','string','logical','string','string'});
A=cellTable(audit,{'ConnectedID','FormationID','Stage','ReactionID','C_ID','ECEligible', ...
    'NoReactionReuse','SharedNoncurrencyOutsideI','IsolationPassed','Status'}, ...
    {'double','double','string','string','string','logical','logical','string','logical','string'});
end

function [Y,stats] = wireConsumingCounters(X,F,support,W,M,rails,O,model,currency,cfg)
rows=cell(0,14);stop=false;
stats=struct('tested',0,'reactionReuseRejected',0,'crossTalkRejected',0,'ECRejected',0, ...
    'knownGeneOverlapRejected',0,'testLimitReached',false,'wiringLimitReached',false);
for x=1:height(X)
    f=X.FormationID(x);rx=F.ReactionIndex(f);
    if X.LeakOrientationIndex(x)>0,rx(end+1)=O.rxn(X.LeakOrientationIndex(x));end %#ok<AGROW>
    for w=reshape(find(W.ConnectedID==X.ConnectedID(x)),1,[])
        if stats.tested>=cfg.counterMaxWiringTests,stats.testLimitReached=true;stop=true;break;end
        stats.tested=stats.tested+1;b=find(M.ModuleID==W.ModuleID(w),1);
        fast=rails(M.FastRailIndex(b));slow=rails(M.SlowRailIndex(b));brx=[fast.rx slow.rx];
        if any(ismember(rx,brx)),stats.reactionReuseRejected=stats.reactionReuseRejected+1;continue;end
        shared=intersect(support{x},union(fast.parts,slow.parts));shared(currency(shared))=[];
        if cfg.counterRequireIsolation&&~isempty(shared),stats.crossTalkRejected=stats.crossTalkRejected+1;continue;end
        assigned=completeECTokens(M.AssignedECs(b));sets=cell(1,numel(assigned));
        for k=1:numel(assigned),sets{k}=assigned(k);end
        % A counter reaction without an EC joins the distinct-EC check only
        % when an EC is required for it (counterRequireFormationEC/LeakEC).
        formationSet=completeECTokens(F.ECs(f));
        hasFormationEC=~isempty(formationSet)||cfg.counterRequireFormationEC;
        if hasFormationEC,sets{end+1}=formationSet;end
        hasLeakEC=false;
        if X.LeakOrientationIndex(x)>0
            leakSet=completeECTokens(X.LeakECs(x));
            hasLeakEC=~isempty(leakSet)||cfg.counterRequireLeakEC;
            if hasLeakEC,sets{end+1}=leakSet;end
        end
        [codes,ecOK]=distinctBranchECs(sets,cfg.counterDistinctECFromBranches);
        if ~ecOK,stats.ECRejected=stats.ECRejected+1;continue;end
        genes=strings(0,1);branchGenes=strings(0,1);
        for j=reshape(rx,1,[])
            [~,g,~]=decodeGPR(string(model.grRules{j}),string(model.rules{j}),string(model.genes));genes=[genes;split(g,';')]; %#ok<AGROW>
        end
        for j=reshape(brx,1,[])
            [~,g,~]=decodeGPR(string(model.grRules{j}),string(model.rules{j}),string(model.genes));branchGenes=[branchGenes;split(g,';')]; %#ok<AGROW>
        end
        knownOverlap=intersect(genes(genes~=""),branchGenes(branchGenes~=""));
        if cfg.branchRejectSharedKnownGenes&&~isempty(knownOverlap)
            stats.knownGeneOverlapRejected=stats.knownGeneOverlapRejected+1;continue;
        end
        formationEC="";leakEC="";
        if hasFormationEC,formationEC=codes(numel(assigned)+1);end
        if hasLeakEC,leakEC=codes(end);end
        rows(end+1,:)={size(rows,1)+1,W.WiringID(w),W.ModuleID(w),X.ConnectedID(x), ...
            X.CounterConnectionID(x),F.FormationID(f),X.I_ID(x),X.C_ID(x),X.CounterType(x), ...
            formationEC,leakEC,numel(unique(codes))==numel(codes),W.BothControlInterfacesMatched(w), ...
            "I_is_consumed_upstream_SF_interfaces_and_all_rates_remain_unverified"}; %#ok<AGROW>
        if size(rows,1)>=cfg.counterMaxWiring,stats.wiringLimitReached=true;stop=true;break;end
    end
    if stop,break;end
end
Y=cellTable(rows,{'CounterWiringID','BranchWiringID','ModuleID','ConnectedID','CounterConnectionID', ...
    'FormationID','I_ID','C_ID','CounterType','AssignedFormationEC','AssignedLeakEC', ...
    'DistinctECs','UpstreamSFInterfacesMatched','Status'}, ...
    {'double','double','double','double','double','double','string','string','string','string','string', ...
    'logical','logical','string'});
end

function writeConsumingCounterOutputs(F,X,Y,FA,A,stats,cstats,wstats,EA,model,cfg)
writetable(EA,fullfile(cfg.outputDir,'09_Counter_EC_Audit.csv'));
writetable(F,fullfile(cfg.outputDir,'09_Counter_Components.csv'));
writetable(FA,fullfile(cfg.outputDir,'09_Counter_Filter_Audit.csv'));
writetable(X,fullfile(cfg.outputDir,'10_Canceller_Counter_Connections.csv'));
writetable(Y,fullfile(cfg.outputDir,'11_Branches_Canceller_Counter_Wiring.csv'));
rows=cell(0,13);
for x=1:height(X)
    f=X.FormationID(x);
    rows(end+1,:)={X.CounterConnectionID(x),X.I_ID(x),X.C_ID(x),F.ReactionID(f),F.ECs(f), ...
        X.LeakReactionID(x),X.LeakECs(x),F.I_LossPerEvent(f),F.C_GainPerEvent(f), ...
        NaN,NaN,X.CounterType(x),"unmeasured_include_counter_loss_in_I_dynamics"}; %#ok<AGROW>
end
KT=cellTable(rows,{'CounterConnectionID','I_ID','C_ID','FormationReactionID','FormationECs', ...
    'LeakReactionID','LeakECs','I_LossPerPrimitiveEvent','C_GainPerPrimitiveEvent', ...
    'EffectiveFormationK_per_s','EffectiveLeakK_per_s','CounterType','KineticsStatus'}, ...
    {'double','string','string','string','string','string','string','double','double','double','double','string','string'});
writetable(KT,fullfile(cfg.outputDir,'10_Counter_Kinetics_Template.csv'));
fid=fopen(fullfile(cfg.outputDir,'Anderson_Counter_Report.txt'),'w');assert(fid>0,'Cannot write consuming counter report.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
report(fid,'ANDERSON CANCELLER -> COUNTER BY DIRECT CONSUMPTION OF I');
report(fid,sprintf('Mode=%s; direct I-to-C formations=%d; EC-resolved=%d; canceller/counter=%d; with branches=%d.', ...
    cfg.counterMode,height(F),nnz(F.ECResolved),height(X),height(Y)));
report(fid,sprintf('Formation-only accumulators=%d; counters with a selected separate removal=%d; require removal=%d.', ...
    nnz(X.LeakOrientationIndex==0),nnz(X.LeakOrientationIndex>0),cfg.counterRequireRemoval));
report(fid,sprintf('Tested formation orientations with recovered gross stoichiometry=%d of %d.',nnz(FA.GrossKnown),height(FA)));
report(fid,'No counter membrane or I-opening evidence is required in this mode. Upstream S/F controls remain separate and unverified.');
report(fid,'Target: I + X -> nu_C C + W. X and W denote full co-substrate/product vectors, retained in the equation CSV.');
report(fid,'Counter formation consumes one I per primitive event; C is a distinct chemical state, not I in another compartment.');
report(fid,'Exact mass balances: dI/dt = f_canceller - (I_LossPerEvent/PrimitiveScale)*v_formation;');
report(fid,'                     dC/dt = (C_GainPerEvent/PrimitiveScale)*v_formation - v_removal/LeakPrimitiveScale.');
report(fid,'For rows without selected removal the selected-module removal flux is zero; this does not establish absence of native-cell C loss.');
report(fid,'Only in a justified linear regime with buffered co-substrates: dI/dt = f_canceller - k_counter*I; dC/dt = nu_C*k_counter*I - kcr*C.');
report(fid,'This changes the canceller dynamics and generally its steady-state response. Weak loading may reduce the disturbance but also slows accumulation.');
report(fid,'For the ideal single-pool Anderson form, the linear I coefficient becomes k*(S-F)-k_counter; a positive root at zero background requires S-F > k_counter/k.');
report(fid,'All selected C-removal reactions are screened against containing or regenerating I; feedback paths are excluded and audited.');
report(fid,'No reaction reuse or noncurrency overlap outside the intended I/C interfaces is admitted under default isolation.');
report(fid,'Formation chemistry is catalogued before EC filtering. Separate EC, gross-stoichiometry, reuse and overlap rejections are written to CSV.');
report(fid,'EC annotations are from model/supplied evidence or exact source-ID AND complete-signature matches in a partial primary-source seed map.');
report(fid,'EC separation does not establish independent proteins, kinetic constants, flux feasibility or absence of endogenous bypasses.');
report(fid,'Counter values integrate the consumed I flux; proportional pulse counts require standardized pulse areas and a suitable retention window.');
report(fid,'A formation-only accumulator is not a complete leaky counter. No leakage kinetics or steady-state capacity has been inferred.');
report(fid,'Read 09_Counter_Components.csv, 10_Canceller_Counter_Connections.csv, 10_Consuming_Counter_Connection_Audit.csv and 11_Branches_Canceller_Counter_Wiring.csv for this run.');
report(fid,'Any membrane-evidence CSVs left from an older run are not outputs of consume_I mode.');
report(fid,sprintf('Formation tests=%d; test cap=%d; component cap=%d.',stats.tested,stats.testLimitReached,stats.componentLimitReached));
report(fid,sprintf('Connection tests=%d; removal tests=%d; connection test cap=%d; removal test cap=%d; connection cap=%d.', ...
    cstats.tested,cstats.leakTests,cstats.testLimitReached,cstats.leakTestLimitReached,cstats.connectionLimitReached));
report(fid,sprintf('Branch joins=%d; reuse=%d; overlap=%d; EC assignment rejected=%d; known gene overlap=%d; test cap=%d; wiring cap=%d.', ...
    wstats.tested,wstats.reactionReuseRejected,wstats.crossTalkRejected,wstats.ECRejected, ...
    wstats.knownGeneOverlapRejected,wstats.testLimitReached,wstats.wiringLimitReached));
report(fid,'\nFORMATION FILTER COUNTS');
for s=reshape(unique(FA.Status,'stable'),1,[]),report(fid,sprintf('  %s: %d',s,nnz(FA.Status==s)));end
report(fid,'CONNECTION/REMOVAL FILTER COUNTS');
for s=reshape(unique(A.Status,'stable'),1,[]),report(fid,sprintf('  %s: %d',s,nnz(A.Status==s)));end
for x=1:min(height(X),cfg.printTop)
    f=X.FormationID(x);
    report(fid,sprintf('\nCOUNTER %d | canceller=%d | I=%s -> C=%s | %s', ...
        X.CounterConnectionID(x),X.ConnectedID(x),X.I_ID(x),X.C_ID(x),X.CounterType(x)));
    report(fid,['  Formation: ' char(F.ReactionID(f)) ' | EC ' char(F.ECs(f)) ' | ' char(F.Equation(f))]);
    if X.LeakOrientationIndex(x)>0
        report(fid,['  Separate removal: ' char(X.LeakReactionID(x)) ' | EC ' char(X.LeakECs(x))]);
    else
        report(fid,'  Separate C-removal: not selected; accumulation only.');
    end
end
end

function [E,A] = resolveStep2EC(model,base,previous,root,cfg)
% New annotations are confined to Step 2; previous successful stages persist.
E=previous;n=numel(model.rxns);source=strings(n,1);
if isfield(model,'circuitSourceBiGG'),source=string(model.circuitSourceBiGG(:));end
status=repmat("unresolved",n,1);status(E.codes~="")="inherited_EC_evidence";
file=findNamed(root,cfg.step2ECFile);
if ~isempty(file)
    supplied=readTextTable(file);
    assert(all(ismember({'SourceReactionID','EC','Reference'},supplied.Properties.VariableNames)), ...
        'Anderson:Step2EC','Step2 EC evidence needs SourceReactionID, EC, Reference.');
    for k=1:height(supplied)
        codes=completeECTokens(supplied.EC(k));if isempty(codes)||supplied.Reference(k)=="",continue;end
        matches=find(source==supplied.SourceReactionID(k)&source~="");
        if ismember('ReactionID',supplied.Properties.VariableNames)&&supplied.ReactionID(k)~=""
            matches=find(string(model.rxns)==supplied.ReactionID(k));
        end
        for j=reshape(matches,1,[])
            E.codes(j)=joinText(unique([completeECTokens(E.codes(j));codes],'stable'),';');
            E.source(j)="supplied_Step2_EC_reference";E.reference(j)=E.reference(j)+"; "+supplied.Reference(k);
            status(j)="supplied_reference_not_independently_checked";
        end
    end
end
if cfg.step2UseVerifiedSeedMap
    seeds=verifiedStep2ECSeeds();
    for j=reshape(find(E.codes==""&source~=""),1,[])
        k=cellFreeSeedMatch(seeds.SourceReactionID,source(j),model,j);if isempty(k),continue;end
        if ~branchSignatureMatches(model.S(:,j),base,seeds.ChemicalSignature(k))
            status(j)="Step2_seed_chemical_signature_mismatch";continue;
        end
        E.codes(j)=seeds.EC(k);E.source(j)="partial_primary_source_Step2_seed_2026_11_13";
        E.reference(j)=seeds.Reference(k);status(j)="explicit_source_ID_and_full_chemistry_matched";
    end
end
A=table(string(model.rxns),source,E.codes,E.source,E.reference,status,E.codes~="", ...
    'VariableNames',{'ReactionID','SourceReactionID','EC','AnnotationSource','Reference','Status','ECResolved'});
blank=table(strings(0,1),strings(0,1),strings(0,1),'VariableNames',{'SourceReactionID','EC','Reference'});
writetable(blank,fullfile(cfg.outputDir,'Step2_EC_Evidence_TEMPLATE.tsv'),'FileType','text','Delimiter','\t');
end

function T = verifiedStep2ECSeeds()
% Primary BiGG/IUBMB pages checked 2026-11-13. This is a PARTIAL map.
% IUBMB entries support chemistry-based EC classes for the exact SBML source
% IDs below; no gene, protein identity, rate or preferred direction is inferred.
rows={ ...
    "PDE1","3.1.4.17;3.1.4.53","camp:-1;h2o:-1;amp:1;h:1","https://bigg.ucsd.edu/models/ic_1306/reactions/PDE1"; ...
    "AMPN","3.2.2.4","amp:-1;h2o:-1;ade:1;r5p:1","https://iubmb.qmul.ac.uk/enzyme/EC3/2/2/4.html"; ...
    "TKT1","2.2.1.1","r5p:-1;xu5p_D:-1;g3p:1;s7p:1","https://bigg.ucsd.edu/models/e_coli_core/reactions/TKT1"; ...
    "TKT2","2.2.1.1","e4p:-1;xu5p_D:-1;f6p:1;g3p:1","https://bigg.ucsd.edu/universal/reactions/TKT2"; ...
    "TALA","2.2.1.2","g3p:-1;s7p:-1;e4p:1;f6p:1","https://bigg.ucsd.edu/universal/reactions/TALA"; ...
    "RPI","5.3.1.6","r5p:-1;ru5p_D:1","https://bigg.ucsd.edu/universal/reactions/RPI"; ...
    "RPE","5.1.3.1","ru5p_D:-1;xu5p_D:1","https://bigg.ucsd.edu/universal/reactions/RPE"; ...
    "GLYOX","3.1.2.6","h2o:-1;lgt_S:-1;gthrd:1;h:1;lac_D:1","https://bigg.ucsd.edu/universal/reactions/GLYOX"; ...
    "D_LACDm","1.1.2.4","ficytc:-2;lac_D:-1;focytc:2;pyr:1","https://iubmb.qmul.ac.uk/enzyme/EC1/1/2/4.html"; ...
    "ALATA_L","2.6.1.2","akg:-1;ala_L:-1;glu_L:1;pyr:1","https://iubmb.qmul.ac.uk/enzyme/EC2/6/1/2.html"; ...
    "ALATA_Lm","2.6.1.2","akg:-1;ala_L:-1;glu_L:1;pyr:1","https://iubmb.qmul.ac.uk/enzyme/EC2/6/1/2.html"; ...
    "GLUCYS","6.3.2.2","atp:-1;cys_L:-1;glu_L:-1;adp:1;glucys:1;h:1;pi:1","https://bigg.ucsd.edu/universal/reactions/GLUCYS"; ...
    "GTHS","6.3.2.3","atp:-1;glucys:-1;gly:-1;adp:1;gthrd:1;h:1;pi:1","https://bigg.ucsd.edu/universal/reactions/GTHS"; ...
    "VALTA","2.6.1.42;2.6.1.6","akg:-1;val_L:-1;3mob:1;glu_L:1","https://bigg.ucsd.edu/universal/reactions/VALTA"; ...
    "GTHP","1.11.1.9","gthrd:-2;h2o2:-1;gthox:1;h2o:2","https://iubmb.qmul.ac.uk/enzyme/EC1/11/1/9.html"};
T=cellTable(rows,{'SourceReactionID','EC','ChemicalSignature','Reference'},repmat({'string'},1,4));
end

function target = step2NeuronCompartment(cfg)
target=string(cfg.neuronCompartment);
if target=="auto",target=string(cfg.outsideCompartment);end
end

function [ctx,A] = prepareStep2NeuronCompartment(ctx,cfg)
% Append actual outside-only orientations. Preserve every upstream O index.
rows=cell(0,5);
if cfg.neuronPlacement=="outside"
    target=step2NeuronCompartment(cfg);tol=1e-8;
    assert(~ismember(target,string(cfg.insideCompartment)), ...
        'Anderson:NeuronCompartment','The neuron must be outside the upstream inside compartments.');
    for j=1:numel(ctx.model.rxns)
        parts=find(ctx.R(:,j)+ctx.P(:,j)>tol);
        if ~any(ctx.comp(parts)==target),continue;end
        dirs=[];if ctx.model.ub(j)>0,dirs=1;end
        if cfg.allowReverse&&ctx.model.lb(j)<0,dirs=[dirs -1];end %#ok<AGROW>
        for direction=dirs
            a=ctx.R(:,j);b=ctx.P(:,j);if direction<0,tmp=a;a=b;b=tmp;end
            status="outside_reaction_crosses_compartments";
            if all(ctx.comp(parts)==target)
                status="outside_reaction_has_no_conversion_products";
                if any(a>tol)&&any(b>tol)
                    status="outside_orientation_admitted";
                    if ~any(ctx.O.rxn==j&ctx.O.direction==direction)
                        o=numel(ctx.O.rxn)+1;
                        ctx.O.rxn(o)=j;ctx.O.direction(o)=direction;
                        ctx.O.R(:,o)=a;ctx.O.P(:,o)=b;ctx.O.N(:,o)=b-a;
                        known=false;
                        if isfield(ctx,'grossKnown'),known=ctx.grossKnown(j);end
                        ctx.O.grossKnown(o)=known;
                        for i=reshape(find(a>tol),1,[]),ctx.O.consumers{i}(end+1)=o;end
                    end
                end
            end
            rows(end+1,:)={string(ctx.model.rxns{j}),direction,target, ...
                joinText(unique(ctx.comp(parts)),';'),status}; %#ok<AGROW>
        end
    end
end
A=cellTable(rows,{'ReactionID','Direction','NeuronCompartment','ParticipantCompartments','Status'}, ...
    {'string','double','string','string','string'});
end

function link = makeStep2ExportLink(rx,direction,counter,input,parts,coefficient,reverseAllowed)
link=struct('rx',rx,'direction',direction,'counter',counter,'input',input,'parts',parts, ...
    'cargoCoefficient',coefficient,'reverseAllowed',reverseAllowed,'cargo',counter, ...
    'adapterO',[],'adapterMets',counter,'adapterRx',[],'adapterFlux',[],'adapterParts',[]);
end

function links = emptyStep2Exports()
link=makeStep2ExportLink(0,0,0,0,[],0,false);links=link([]);
end

function A = step2ExportAuditTable(rows)
A=cellTable(rows,{'ParentID','ExportCargoID','TargetCompartment','ReactionID','Direction', ...
    'NeuronInputID','HasGPR','ReverseAllowed','ExportCoefficient','ImportCoefficient', ...
    'SharedNoncurrency','Status','Equation'}, ...
    {'double','string','string','string','double','string','logical','logical', ...
    'double','double','string','string','string'});
end

function [links,A,stats] = findStep2Exports(parent,counter,seed,prefixRx,prefixParts,ctx,cfg)
% True export: same canonical chemical, different existing compartment state.
% This function never relabels a conversion or an exchange sink as transport.
links=emptyStep2Exports();rows=cell(0,13);tol=1e-8;
stats=struct('tests',0,'testLimitReached',false,'linkLimitReached',false);
target=step2NeuronCompartment(cfg);outside=find(ctx.base==ctx.base(counter)&ctx.comp==target);
if isempty(outside)
    rows={parent,string(ctx.model.mets{counter}),target,"",0,"",false,false,NaN,NaN,"", ...
        "no_matching_counter_state_in_target_compartment",""};
    A=step2ExportAuditTable(rows);return;
end
stop=false;
for j=reshape(find(ctx.R(counter,:)+ctx.P(counter,:)>tol),1,[])
    dirs=[];if ctx.model.ub(j)>0,dirs=1;end
    if cfg.allowReverse&&ctx.model.lb(j)<0,dirs=[dirs -1];end %#ok<AGROW>
    for direction=dirs
        if stats.tests>=cfg.neuronTransportMaxTests,stats.testLimitReached=true;stop=true;break;end
        stats.tests=stats.tests+1;a=ctx.R(:,j);b=ctx.P(:,j);
        if direction<0,tmp=a;a=b;b=tmp;end
        net=b-a;parts=find(a+b>tol);cargo=outside(net(outside)>tol);
        input=0;ci=NaN;if numel(cargo)==1,input=cargo;ci=full(net(input));end
        ce=-full(net(counter));shared=step2Overlap(parts,prefixParts,counter,ctx.currency);
        [gpr,~,~]=decodeGPR(string(ctx.model.grRules{j}),string(ctx.model.rules{j}),string(ctx.model.genes));
        reverseAllowed=(direction==1&&cfg.allowReverse&&ctx.model.lb(j)<0) ...
            ||(direction==-1&&ctx.model.ub(j)>0);
        status="accepted_counter_export";
        if ismember(j,prefixRx),status="export_reuses_upstream_reaction";
        elseif ce<=tol||b(counter)>tol,status="reaction_does_not_export_counter";
        elseif numel(cargo)~=1,status="no_unique_same_chemical_external_product";
        elseif abs(ce-ci)>tol*max([1 ce ci])||a(input)>tol,status="cargo_copies_not_conserved";
        elseif any(~ismember(ctx.comp(parts),[ctx.comp(counter);target]))
            status="export_crosses_an_additional_compartment";
        elseif any(ctx.base(parts)==ctx.base(counter)&~ismember(parts,[counter;input]))
            status="export_has_additional_counter_pools";
        elseif any(ctx.base(parts)==ctx.base(seed)),status="export_itself_couples_to_I";
        elseif ~isempty(shared),status="export_shares_other_upstream_noncurrency_states";
        elseif cfg.neuronTransportRequireGPR&&strlength(gpr)==0,status="export_GPR_required_but_missing";
        elseif ~step2TransportGenesOK(j,prefixRx,ctx,cfg),status="export_shares_known_upstream_genes";
        end
        inputID="";if input>0,inputID=string(ctx.model.mets{input});end
        rows(end+1,:)={parent,string(ctx.model.mets{counter}),target,string(ctx.model.rxns{j}), ...
            direction,inputID,strlength(gpr)>0,reverseAllowed,ce,ci, ...
            joinText(string(ctx.model.mets(shared)),';'),status,equation(a,b,ctx.model)}; %#ok<AGROW>
        if status~="accepted_counter_export",continue;end
        links(end+1)=makeStep2ExportLink(j,direction,counter,input,parts,ce,reverseAllowed); %#ok<AGROW>
        if numel(links)>=cfg.neuronTransportMaxLinksPerParent,stats.linkLimitReached=true;stop=true;break;end
    end
    if stop,break;end
end
if isempty(rows)
    rows={parent,string(ctx.model.mets{counter}),target,"",0,"",false,false,NaN,NaN,"", ...
        "no_bound_allowed_counter_reaction",""};
end
A=step2ExportAuditTable(rows);
end

function A = step2AdapterAuditTable(rows)
A=cellTable(rows,{'ParentID','AdapterDepth','FromID','ToID','ReactionID','Direction', ...
    'ECResolved','SharedUpstream','SharedAdapter','IStateIDs','Status'}, ...
    {'double','double','string','string','string','double','logical','string','string','string','string'});
end

function [links,TA,AA,stats] = findStep2AdapterExports(parent,counter,seed,prefixRx,prefixParts,control,ctx,G,cfg)
% Breadth-first path STATES preserve full side-participant history. Each
% adapter stays in the counter compartment; only the existing exporter crosses
% the membrane. The exported species need not be the original m_c chemical.
links=emptyStep2Exports();TA=step2ExportAuditTable(cell(0,13));rows=cell(0,11);
stats=struct('tests',0,'testLimitReached',false,'linkLimitReached',false,'states',0, ...
    'expansions',0,'stateLimitReached',false,'expansionLimitReached',false, ...
    'depthLimitReached',false,'maxDepthReached',0);
queue=struct('node',counter,'path',[],'mets',counter,'parts',[]);head=1;tol=1e-8;
seen=containers.Map('KeyType','char','ValueType','logical');stop=false;O=ctx.O;
while head<=numel(queue)
    if stats.states>=cfg.neuronAdapterMaxStates,stats.stateLimitReached=true;break;end
    state=queue(head);head=head+1;stats.states=stats.states+1;depth=numel(state.path);
    stats.maxDepthReached=max(stats.maxDepthReached,depth);
    adapterRx=O.rxn(state.path);fullRx=unique([prefixRx adapterRx],'stable');
    fullParts=union(prefixParts,state.parts);
    remaining=cfg.neuronTransportMaxTests-stats.tests;
    if remaining<=0,stats.testLimitReached=true;break;end
    localCfg=cfg;localCfg.neuronTransportMaxTests=remaining;
    localCfg.neuronTransportMaxLinksPerParent=cfg.neuronTransportMaxLinksPerParent-numel(links);
    [found,a,ts]=findStep2Exports(parent,state.node,seed,fullRx,fullParts,ctx,localCfg);
    TA=[TA;a];stats.tests=stats.tests+ts.tests; %#ok<AGROW>
    stats.testLimitReached=stats.testLimitReached||ts.testLimitReached;
    for k=1:numel(found)
        found(k).counter=counter;found(k).cargo=state.node;
        found(k).adapterO=state.path;found(k).adapterMets=state.mets;
        found(k).adapterRx=adapterRx;found(k).adapterFlux=G.scale(state.path);
        found(k).adapterParts=state.parts;found(k).parts=union(found(k).parts,state.parts);
        links(end+1)=found(k); %#ok<AGROW>
    end
    if numel(links)>=cfg.neuronTransportMaxLinksPerParent,stats.linkLimitReached=true;break;end
    if stats.testLimitReached,break;end
    if depth>=cfg.neuronAdapterMaxDepth
        stats.depthLimitReached=stats.depthLimitReached||~isempty(O.consumers{state.node});continue;
    end
    allowed=[];if depth==0,allowed=counter;end
    for o=reshape(O.consumers{state.node},1,[])
        if stats.expansions>=cfg.neuronAdapterMaxExpansions
            stats.expansionLimitReached=true;stop=true;break;
        end
        stats.expansions=stats.expansions+1;j=O.rxn(o);parts=G.parts{o};
        sharedUp=step2Overlap(parts,prefixParts,allowed,ctx.currency);
        sharedPath=step2Overlap(parts,state.parts,state.node,ctx.currency);
        iParts=parts(ctx.base(parts)==ctx.base(seed));
        edges=G.byFrom{state.node};edges=edges(G.o(edges)==o);
        status="accepted_adapter_step";
        if ismember(j,fullRx),status="adapter_reuses_upstream_or_adapter_reaction";
        elseif any(ctx.comp(parts)~=ctx.comp(counter)),status="adapter_leaves_counter_compartment";
        elseif step2NeuronIForbidden(parts,seed,prefixParts,ctx,cfg),status="adapter_couples_to_I_pool";
        elseif ~isempty(sharedUp),status="adapter_shares_other_upstream_noncurrency_states";
        elseif ~isempty(sharedPath),status="adapter_reuses_nonadjacent_noncurrency_state";
        elseif ~G.eligible(o),status="adapter_gross_primitive_or_size_policy_failed";
        elseif ~G.ecKnown(o),status="adapter_EC_unresolved";
        elseif isempty(edges),status="adapter_not_a_unit_consuming_conversion";
        else
            [~,ok]=step2AssignEC([control adapterRx j],ctx,cfg);
            if ~ok,status="adapter_EC_or_known_gene_assignment_conflict";end
        end
        if status~="accepted_adapter_step"
            rows(end+1,:)={parent,depth+1,string(ctx.model.mets{state.node}),"", ...
                string(ctx.model.rxns{j}),O.direction(o),G.ecKnown(o), ...
                joinText(string(ctx.model.mets(sharedUp)),';'),joinText(string(ctx.model.mets(sharedPath)),';'), ...
                joinText(string(ctx.model.mets(iParts)),';'),status}; %#ok<AGROW>
            continue;
        end
        for e=reshape(edges,1,[])
            node=G.to(e);chosen=[state.path o];key=sprintf('%d,',[node sort(chosen)]);
            reason="accepted_adapter_step";
            if isKey(seen,key),reason="duplicate_adapter_state";
            elseif numel(queue)>=cfg.neuronAdapterMaxStates
                reason="adapter_state_retention_budget_reached";stats.stateLimitReached=true;
            else
                seen(key)=true;queue(end+1)=struct('node',node,'path',chosen, ...
                    'mets',[state.mets node],'parts',union(state.parts,parts)); %#ok<AGROW>
            end
            rows(end+1,:)={parent,depth+1,string(ctx.model.mets{state.node}),string(ctx.model.mets{node}), ...
                string(ctx.model.rxns{j}),O.direction(o),G.ecKnown(o),"","", ...
                joinText(string(ctx.model.mets(iParts)),';'),reason}; %#ok<AGROW>
        end
    end
    if stop,break;end
end
AA=step2AdapterAuditTable(rows);
end

function T = describeStep2AdapterLinks(links,parents,ctx)
rows=cell(0,12);
for k=1:numel(links)
    v=links(k);rows(end+1,:)={k,parents(k),string(ctx.model.mets{v.counter}), ...
        string(ctx.model.mets{v.cargo}),string(ctx.model.mets{v.input}),numel(v.adapterO), ...
        joinText(string(ctx.model.rxns(v.adapterRx)),';'),joinText(string(ctx.O.direction(v.adapterO)),';'), ...
        joinText(string(ctx.model.mets(v.adapterMets)),' -> '),string(ctx.model.rxns{v.rx}), ...
        v.direction,v.reverseAllowed}; %#ok<AGROW>
end
T=cellTable(rows,{'InterfaceID','ParentID','CounterID','ExportCargoID','NeuronInputID', ...
    'AdapterReactionCount','AdapterReactionIDs','AdapterDirections','AdapterMetabolitePath', ...
    'ExportReactionID','ExportDirection','ExportReverseAllowed'}, ...
    {'double','double','string','string','string','double','string','string','string','string','double','logical'});
end

function ok = step2TransportGenesOK(j,otherRx,ctx,cfg)
ok=true;if j==0||~cfg.branchRejectSharedKnownGenes,return;end
[~,value,~]=decodeGPR(string(ctx.model.grRules{j}),string(ctx.model.rules{j}),string(ctx.model.genes));
genes=split(value,';');genes=genes(genes~="");if isempty(genes),return;end
for k=reshape(unique(otherRx),1,[])
    [~,value,~]=decodeGPR(string(ctx.model.grRules{k}),string(ctx.model.rules{k}),string(ctx.model.genes));
    others=split(value,';');if ~isempty(intersect(genes,others(others~=""))),ok=false;return;end
end
end

function blocked = step2NeuronIForbidden(parts,seed,prefixParts,ctx,cfg)
same=ctx.base(parts)==ctx.base(seed);
if cfg.neuronPlacement=="outside"
    % Distinct compartments can isolate I. Another spelling of the same
    % chemical in the same outside compartment cannot bypass an I importer.
    upstreamI=prefixParts(ctx.base(prefixParts)==ctx.base(seed));
    same=same&ismember(ctx.comp(parts),unique(ctx.comp(upstreamI)));
end
blocked=any(same);
end

function G = buildStep2Graph(ctx,cfg)
O=ctx.O;m=numel(ctx.model.mets);n=numel(O.rxn);tol=1e-8;
G=struct('o',[],'from',[],'to',[],'byFrom',{cell(m,1)},'byTo',{cell(m,1)}, ...
    'parts',{cell(n,1)},'scale',nan(1,n),'eligible',false(1,n),'ecKnown',false(1,n));
for o=1:n
    parts=find(O.R(:,o)+O.P(:,o)>tol);G.parts{o}=parts;
    codes=completeECTokens(ctx.E.codes(O.rxn(o)));G.ecKnown(o)=~isempty(codes);
    if ~O.grossKnown(o)||numel(parts)>cfg.step2MaxParticipants||any(ismember(codes,cfg.excludedECs)),continue;end
    [scale,ok]=primitiveReactionScale(O.R(:,o),O.P(:,o));if ~ok,continue;end
    G.scale(o)=scale;G.eligible(o)=true;a=O.R(:,o)*scale;b=O.P(:,o)*scale;
    from=find(abs(a-1)<tol&b<tol&~ctx.currency);to=find(abs(b-1)<tol&a<tol&~ctx.currency);
    for u=reshape(from,1,[])
        for v=reshape(to,1,[])
            e=numel(G.o)+1;G.o(e)=o;G.from(e)=u;G.to(e)=v;
            G.byFrom{u}(end+1)=e;G.byTo{v}(end+1)=e;
        end
    end
end
end

function shared = step2Overlap(a,b,allowed,currency)
shared=intersect(a,b);shared(currency(shared)|ismember(shared,allowed))=[];
end

function [rx,parts,control,b,x,c] = step2Prefix(ctx,w)
b=find(ctx.M.ModuleID==ctx.Y.ModuleID(w),1);x=find(ctx.X.CounterConnectionID==ctx.Y.CounterConnectionID(w),1);
c=find(ctx.C.ConnectedID==ctx.Y.ConnectedID(w),1);
[rx,parts]=counterCancellerSupport(c,ctx.C,ctx.D,ctx.cores,ctx.O,ctx.T,ctx.L,ctx.Q,ctx.R,ctx.P);
fast=ctx.rails(ctx.M.FastRailIndex(b));slow=ctx.rails(ctx.M.SlowRailIndex(b));
control=[fast.rx slow.rx];rx=[rx control];parts=union(parts,union(fast.parts,slow.parts));
% Counter reactions are always part of the circuit (rx, parts), but take
% part in the distinct-EC assignment only when they carry an EC: the counter
% may use EC-less chemistry (counterRequireFormationEC/LeakEC=false).
f=ctx.X.FormationID(x);fr=ctx.F.ReactionIndex(f);rx(end+1)=fr;
if ctx.E.codes(fr)~="",control(end+1)=fr;end
parts=union(parts,ctx.XP{x});
if ctx.X.LeakOrientationIndex(x)>0
    lr=ctx.O.rxn(ctx.X.LeakOrientationIndex(x));rx(end+1)=lr;
    if ctx.E.codes(lr)~="",control(end+1)=lr;end
end
rx=unique(rx,'stable');
end

function [codes,ok] = step2AssignEC(rx,ctx,cfg)
% ctx.ecTokenSets/geneTokenSets (built once in searchStep2Circuit) avoid
% re-parsing EC text and gene rules inside the Block 0/neuron inner loops.
cached=isfield(ctx,'ecTokenSets')&&numel(ctx.ecTokenSets)==numel(ctx.model.rxns);
sets=cell(1,numel(rx));
for k=1:numel(rx)
    if cached,sets{k}=ctx.ecTokenSets{rx(k)};else,sets{k}=completeECTokens(ctx.E.codes(rx(k)));end
end
[codes,ok]=distinctBranchECs(sets,cfg.step2DistinctECs);if ~ok,return;end
if cfg.branchRejectSharedKnownGenes
    seen=strings(0,1);
    for j=reshape(rx,1,[])
        if cached
            genes=ctx.geneTokenSets{j};
        else
            [~,g,~]=decodeGPR(string(ctx.model.grRules{j}),string(ctx.model.rules{j}),string(ctx.model.genes));
            genes=split(g,';');genes=genes(genes~="");
        end
        if ~isempty(intersect(seen,genes)),ok=false;return;end
        seen=union(seen,genes);
    end
end
end

function A = step2EntryRows(parent,input,seed,prefixRx,prefixParts,ctx,G,cfg)
rows=cell(0,10);os=ctx.O.consumers{input};
if isempty(os)
    rows={parent,string(ctx.model.mets{input}),"",0,false,false,"", ...
        "no_consuming_reaction_in_neuron_compartment",ctx.comp(input),""};
end
for o=reshape(os,1,[])
    j=ctx.O.rxn(o);parts=G.parts{o};shared=step2Overlap(parts,prefixParts,input,ctx.currency);
    madeI=find(ctx.base==ctx.base(seed)&ctx.O.P(:,o)>1e-8);
    regenerated=~isempty(madeI);status="eligible_entry_for_neuron_search";
    if cfg.neuronPlacement=="outside"&&step2NeuronIForbidden(madeI,seed,prefixParts,ctx,cfg)
        status="entry_external_I_feeds_existing_upstream_pool";
    elseif step2NeuronIForbidden(parts,seed,prefixParts,ctx,cfg)
        if regenerated,status="entry_would_regenerate_I";
        else,status="I_participates_in_neuron_chemistry";end
    elseif ismember(j,prefixRx),status="entry_reuses_upstream_reaction";
    elseif ~isempty(shared),status="entry_shares_noncurrency_outside_neuron_input";
    elseif ~G.eligible(o),status="entry_gross_primitive_or_size_policy_failed";
    elseif ~G.ecKnown(o),status="entry_EC_unresolved";
    end
    rows(end+1,:)={parent,string(ctx.model.mets{input}),string(ctx.model.rxns{j}),ctx.O.direction(o), ...
        G.ecKnown(o),regenerated,joinText(string(ctx.model.mets(shared)),';'),status, ...
        ctx.comp(input),joinText(string(ctx.model.mets(madeI)),';')}; %#ok<AGROW>
end
A=rows;
end

function rails = emptyStep2NeuronChains()
rails=struct('input',{},'mid',{},'output',{},'o',{},'rx',{},'parts',{}, ...
    'mets',{},'conversionO',{},'removalO',{});
end

function T = step2ConsumingChainTable(rows)
T=cellTable(rows,{'InputID','OutputID','NeuronCompartment','ConversionCount','ReactionIDs', ...
    'Directions','OrientationIndices','MetabolitePath','NetInputChange','NetIChange', ...
    'IInteraction','ReusesUpstreamReaction','SharedUpstreamIDs','NonadjacentSharedIDs', ...
    'ECGeneAssignmentPassed','IndependentRemovalCount','Classification','NetEquation','SearchPass'}, ...
    {'string','string','string','double','string','string','string','string','double','double', ...
    'logical','logical','string','string','logical','double','string','string','string'});
end

function [rails,stats,T] = findStep2NeuronRails(input,seed,prefixRx,prefixParts,ctx,G,cfg,control)
if nargin<8,control=[];end
[rails,stats,T]=scanStep2NeuronChains(input,seed,prefixRx,prefixParts,control,ctx,G,cfg,true);
stats.diagnosticStates=0;stats.diagnosticTests=0;stats.diagnosticLimitReached=false;
if isempty(rails)&&cfg.neuronDiagnoseRejectedChains
    % A separate budget prevents feedback routes from displacing a valid
    % isolated-neuron search. Diagnostics never enter rails or full motifs.
    diagnosticCfg=cfg;diagnosticCfg.neuronChainMaxStates=cfg.neuronDiagnosticMaxStates;
    diagnosticCfg.neuronMaxTestsPerRail=cfg.neuronDiagnosticMaxTests;
    [~,d,DT]=scanStep2NeuronChains(input,seed,prefixRx,prefixParts,control,ctx,G,diagnosticCfg,false);
    T=[T;DT];stats.diagnosticStates=d.states;stats.diagnosticTests=d.tests;
    stats.diagnosticLimitReached=d.stateLimitReached||d.testLimitReached;
end
end

function [rails,stats,T] = scanStep2NeuronChains(input,seed,prefixRx,prefixParts,control,ctx,G,cfg,strict)
% Paths stay in the input compartment and consume the exact input state.
% Every substrate/product contributes to overlap checks. I can be a side
% product in a DIAGNOSTIC, but is never used as a neuron path connector.
rails=emptyStep2NeuronChains();rows=cell(0,19);O=ctx.O;tol=1e-8;
stats=struct('tests',0,'states',0,'testLimitReached',false,'railLimitReached',false, ...
    'stateLimitReached',false,'depthLimitReached',false);
queue=struct('node',input,'path',[],'mets',input,'parts',[], ...
    'sharedUp',[],'sharedHistory',[]);head=1;stop=false;
seen=containers.Map('KeyType','char','ValueType','logical');
while head<=numel(queue)
    if stats.states>=cfg.neuronChainMaxStates,stats.stateLimitReached=true;break;end
    state=queue(head);head=head+1;stats.states=stats.states+1;depth=numel(state.path);
    rx=O.rxn(state.path);iConflict=step2NeuronIForbidden(state.parts,seed,prefixParts,ctx,cfg);
    reuse=any(ismember(rx,prefixRx));[~,ecOK]=step2AssignEC([control rx],ctx,cfg);
    clean=~iConflict&&~reuse&&isempty(state.sharedUp)&&isempty(state.sharedHistory)&&ecOK;
    removals=0;
    if depth>0
        weights=reshape(G.scale(state.path),[],1);net=O.N(:,state.path)*weights;
        if strict&&clean&&depth>=cfg.neuronChainMinDepth&&net(input)<-tol
            for o=reshape(O.consumers{state.node},1,[])
                if stats.tests>=cfg.neuronMaxTestsPerRail,stats.testLimitReached=true;stop=true;break;end
                stats.tests=stats.tests+1;parts=G.parts{o};j=O.rxn(o);
                if ~G.eligible(o)||~G.ecKnown(o)||ismember(j,[prefixRx rx]) ...
                    ||any(ctx.comp(parts)~=ctx.comp(input)) ...
                    ||abs(full(O.R(state.node,o))*G.scale(o)-1)>tol||O.P(state.node,o)>tol ...
                    ||step2NeuronIForbidden(parts,seed,prefixParts,ctx,cfg) ...
                    ||~isempty(step2Overlap(parts,prefixParts,[],ctx.currency)) ...
                    ||~isempty(step2Overlap(parts,state.parts,state.node,ctx.currency)),continue;end
                [~,ok]=step2AssignEC([control rx j],ctx,cfg);if ~ok,continue;end
                removals=removals+1;
                rails(end+1)=struct('input',input,'mid',state.mets(2),'output',state.node, ...
                    'o',[state.path o],'rx',[rx j],'parts',union(state.parts,parts), ...
                    'mets',state.mets,'conversionO',state.path,'removalO',o); %#ok<AGROW>
                if numel(rails)>=cfg.neuronMaxRails,stats.railLimitReached=true;stop=true;break;end
            end
        end
        classification="diagnostic_feedback_overlap_or_control_conflict";
        if clean
            classification="isolated_chain_without_independent_output_removal";
            if depth<cfg.neuronChainMinDepth,classification="consuming_prefix_below_minimum_neuron_depth";end
            if removals>0,classification="isolated_neuron_with_independent_output_removal";end
        end
        if clean&&~strict,classification="diagnostic_clean_chain_removal_not_evaluated";end
        if strict&&clean&&removals==0&&stats.testLimitReached
            classification="isolated_chain_removal_search_budget_reached";
        end
        if net(input)>=-tol,classification="diagnostic_input_recycled_or_not_net_consumed";end
        pass="diagnostic_only";if strict,pass="isolated_search";end
        rows(end+1,:)={string(ctx.model.mets{input}),string(ctx.model.mets{state.node}),ctx.comp(input), ...
            depth,joinText(string(ctx.model.rxns(rx)),';'),joinText(string(O.direction(state.path)),';'), ...
            joinText(string(state.path),';'),joinText(string(ctx.model.mets(state.mets)),' -> '), ...
            full(net(input)),full(sum(net(ctx.base==ctx.base(seed)))),iConflict,reuse, ...
            joinText(string(ctx.model.mets(state.sharedUp)),';'), ...
            joinText(string(ctx.model.mets(state.sharedHistory)),';'),ecOK,removals,classification, ...
            equation(max(-net,0),max(net,0),ctx.model),pass}; %#ok<AGROW>
    end
    if stop,break;end
    if depth>=cfg.neuronChainMaxDepth
        stats.depthLimitReached=stats.depthLimitReached||~isempty(O.consumers{state.node});continue;
    end
    allowed=[];if depth==0,allowed=input;end
    for o=reshape(O.consumers{state.node},1,[])
        if stats.tests>=cfg.neuronMaxTestsPerRail,stats.testLimitReached=true;stop=true;break;end
        stats.tests=stats.tests+1;parts=G.parts{o};j=O.rxn(o);
        if ~G.eligible(o)||~G.ecKnown(o)||ismember(j,rx)||any(ctx.comp(parts)~=ctx.comp(input)),continue;end
        sharedUp=union(state.sharedUp,step2Overlap(parts,prefixParts,allowed,ctx.currency));
        sharedHistory=union(state.sharedHistory,step2Overlap(parts,state.parts,state.node,ctx.currency));
        chosen=[state.path o];
        if strict
            if ismember(j,prefixRx)||step2NeuronIForbidden(parts,seed,prefixParts,ctx,cfg) ...
                ||~isempty(sharedUp)||~isempty(sharedHistory),continue;end
            [~,ok]=step2AssignEC([control O.rxn(chosen)],ctx,cfg);if ~ok,continue;end
        end
        edges=G.byFrom{state.node};edges=edges(G.o(edges)==o);
        for e=reshape(edges,1,[])
            node=G.to(e);
            if ismember(node,state.mets)||ctx.base(node)==ctx.base(seed),continue;end
            key=sprintf('%d,',[node chosen]);if isKey(seen,key),continue;end
            if numel(queue)>=cfg.neuronChainMaxStates,stats.stateLimitReached=true;continue;end
            seen(key)=true;queue(end+1)=struct('node',node,'path',chosen, ...
                'mets',[state.mets node],'parts',union(state.parts,parts), ...
                'sharedUp',sharedUp,'sharedHistory',sharedHistory); %#ok<AGROW>
        end
    end
    if stop,break;end
end
T=step2ConsumingChainTable(rows);
end

function G = buildStep2Block0Graph(ctx,cfg)
% A separate graph leaves the neuron's primitive one-copy policy unchanged.
% Edges identify net-consumed/net-produced connectors; full gross sides are
% retained for isolation, positive-flux balance and substrate-startup checks.
O=ctx.O;m=numel(ctx.model.mets);n=numel(O.rxn);tol=1e-8;
G=struct('o',[],'from',[],'to',[],'byFrom',{cell(m,1)},'byTo',{cell(m,1)}, ...
    'parts',{cell(n,1)},'eligible',false(1,n),'ecKnown',false(1,n));
for o=1:n
    parts=find(O.R(:,o)+O.P(:,o)>tol);G.parts{o}=parts;
    codes=completeECTokens(ctx.E.codes(O.rxn(o)));G.ecKnown(o)=~isempty(codes);
    if ~O.grossKnown(o)||numel(parts)>cfg.block0MaxParticipants ...
            ||any(ismember(codes,cfg.excludedECs)),continue;end
    net=O.P(:,o)-O.R(:,o);G.eligible(o)=true;
    from=find(net < -tol&~ctx.currency);to=find(net > tol&~ctx.currency);
    for u=reshape(from,1,[])
        for v=reshape(to,1,[])
            if ctx.comp(u)~=ctx.comp(v),continue;end
            e=numel(G.o)+1;G.o(e)=o;G.from(e)=u;G.to(e)=v;
            G.byFrom{u}(end+1)=e;G.byTo{v}(end+1)=e;
        end
    end
end
end

function result = emptyStep2Path(target)
result=struct('found',false,'mets',[],'o',[],'target',target,'distance',Inf, ...
    'status',"no_isolated_path_within_configured_search",'states',0,'expansions',0, ...
    'stateLimitReached',false,'expansionLimitReached',false,'depthLimitReached',false,'ECComplete',false, ...
    'alternativeLimitReached',false,'topologyReachable',false,'searchStatus',"not_started", ...
    'moduleTests',0,'balanceRejected',0,'startupRejected',0,'seedIndependentRejected',0,'solverUnresolved',0, ...
    'flux',[],'net',[],'feed',[],'waste',[],'internal',[],'buffered',[]);
end

function [result,alternatives,audit] = findStep2Block0(target,prefixRx,prefixParts,controlRx,ctx,G,cfg,requireEC)
% Enumerate bounded simple connector traces, then balance their reaction
% supports. Internal side intermediates may cancel. This is not exhaustive
% enumeration of arbitrary branched hypernetworks or a whole-cell flux test.
result=emptyStep2Path(target);alternatives=repmat(result,1,0);O=ctx.O;
audit=cellTable(cell(0,4),{'Reason','OrientationCount','SharedStateIDs','ExampleReactionIDs'}, ...
    {'string','double','string','string'});
starts=find(ctx.base=="camp"&ctx.comp==ctx.comp(target)&~ctx.currency);
if isempty(starts),result.status="intracellular_cAMP_anchor_missing";return;end
if ctx.currency(target),result.status="m_in_is_currency";return;end
reasons=strings(1,numel(O.rxn));sharedIDs=cell(1,numel(O.rxn));
for o=1:numel(O.rxn)
    shared=step2Overlap(G.parts{o},prefixParts,target,ctx.currency);sharedIDs{o}=shared;
    if ~G.eligible(o),reasons(o)="gross_unknown_participant_cap_or_excluded_EC";
    elseif ismember(O.rxn(o),prefixRx),reasons(o)="reuses_downstream_reaction";
    elseif ~isempty(shared),reasons(o)="shares_downstream_noncurrency";
    elseif any(ctx.comp(G.parts{o}(~ctx.currency(G.parts{o})))~=ctx.comp(target))
        reasons(o)="noncurrency_participant_in_other_compartment";
    elseif O.R(target,o)>1e-8,reasons(o)="consumes_m_in";
    elseif ~G.ecKnown(o)&&requireEC,reasons(o)="missing_complete_EC";
    end
end
for reason=reshape(unique(reasons(reasons~=""),'stable'),1,[])
    os=find(reasons==reason);shared=[];
    for o=reshape(os,1,[]),shared=union(shared,sharedIDs{o});end
    examples=unique(O.rxn(os),'stable');examples=examples(1:min(5,numel(examples)));
    audit(end+1,:)={reason,numel(os),joinText(string(ctx.model.mets(shared)),';'), ...
        joinText(string(ctx.model.rxns(examples)),';')}; %#ok<AGROW>
end
allowed=reasons=="";
if requireEC
    [~,ok]=step2AssignEC(controlRx,ctx,cfg);
    if ~ok,result.status="downstream_control_EC_or_known_gene_conflict";return;end
end
seenSupports=containers.Map('KeyType','char','ValueType','logical');stop=false;
for start=reshape(starts,1,[])
    if start==target||ismember(start,prefixParts),continue;end
    permitted=allowed&full(O.P(start,:)<1e-8);
    remaining=inf(numel(ctx.model.mets),1);remaining(target)=0;q=target;head=1;
    while head<=numel(q)
        v=q(head);head=head+1;
        for e=reshape(G.byTo{v},1,[])
            if ~permitted(G.o(e)),continue;end
            u=G.from(e);if isfinite(remaining(u)),continue;end
            remaining(u)=remaining(v)+1;q(end+1)=u; %#ok<AGROW>
        end
    end
    if ~isfinite(remaining(start)),continue;end
    result.topologyReachable=true;
    if remaining(start)>cfg.block0MaxDepth,result.depthLimitReached=true;continue;end
    if result.states>=cfg.block0MaxStates,result.stateLimitReached=true;break;end
    stateM={start};stateO={[]};head=1;tail=1;result.states=result.states+1;
    while head<=tail
        mets=stateM{head};path=stateO{head};
        stateM{head}=[];stateO{head}=[];head=head+1;
        u=mets(end);depth=numel(path);
        if depth>=cfg.block0MaxDepth,result.depthLimitReached=true;continue;end
        for e=reshape(G.byFrom{u},1,[])
            o=G.o(e);v=G.to(e);
            if ~permitted(o)||ismember(v,mets)||~isfinite(remaining(v)) ...
                    ||any(O.rxn(path)==O.rxn(o)),continue;end
            if depth+1+remaining(v)>cfg.block0MaxDepth,result.depthLimitReached=true;continue;end
            if result.expansions>=cfg.block0MaxExpansions
                result.expansionLimitReached=true;stop=true;break;
            end
            result.expansions=result.expansions+1;next=[path o];
            if requireEC
                [~,ok]=step2AssignEC([controlRx O.rxn(next)],ctx,cfg);if ~ok,continue;end
            end
            if v==target
                key=sprintf('%d,',[start sort(next)]);if isKey(seenSupports,key),continue;end
                seenSupports(key)=true;result.moduleTests=result.moduleTests+1;
                [candidate,reason]=validateStep2Block0Module(start,target,next,ctx);
                if reason~="accepted"
                    result.(char(reason))=result.(char(reason))+1;continue;
                end
                candidate.found=true;candidate.mets=[mets v];candidate.o=next;candidate.distance=numel(next);
                candidate.ECComplete=all(G.ecKnown(next));
                candidate.status="isolated_balanced_Block0_EC_qualified";
                if ~requireEC,candidate.status="isolated_structural_path_not_EC_qualified";end
                alternatives(end+1)=candidate; %#ok<AGROW>
                if numel(alternatives)>=cfg.block0MaxAlternatives
                    result.alternativeLimitReached=true;stop=true;break;
                end
                continue;
            end
            if result.states>=cfg.block0MaxStates
                result.stateLimitReached=true;stop=true;break;
            end
            tail=tail+1;stateM{tail}=[mets v];stateO{tail}=next; %#ok<AGROW>
            result.states=result.states+1;
        end
        if stop,break;end
    end
    if stop,break;end
end
result.searchStatus="bounded_connector_trace_enumeration_completed";
if ~result.topologyReachable,result.searchStatus="disconnected_after_downstream_and_annotation_filters";end
if result.depthLimitReached,result.searchStatus="depth_bound_reached";end
if result.solverUnresolved>0,result.searchStatus="positive_flux_solver_unresolved";end
if result.stateLimitReached,result.searchStatus="path_state_budget_reached";end
if result.expansionLimitReached,result.searchStatus="expansion_budget_reached";end
if result.alternativeLimitReached,result.searchStatus="alternative_retention_limit_reached";end
if isempty(alternatives),result.status="no_accepted_Block0_"+result.searchStatus;return;end
[~,order]=sort([alternatives.distance]);alternatives=alternatives(order);
fields={'states','expansions','stateLimitReached','expansionLimitReached','depthLimitReached', ...
    'alternativeLimitReached','topologyReachable','searchStatus','moduleTests','balanceRejected', ...
    'startupRejected','seedIndependentRejected','solverUnresolved'};
for k=1:numel(alternatives)
    for j=1:numel(fields),alternatives(k).(fields{j})=result.(fields{j});end
end
result=alternatives(1);
end

function [candidate,reason] = validateStep2Block0Module(start,target,path,ctx)
candidate=emptyStep2Path(target);reason="balanceRejected";tol=1e-8;
a=ctx.O.R(:,path);b=ctx.O.P(:,path);N=b-a;
internal=find(~ctx.currency&any(a>tol,2)&any(b>tol,2));
internal(ismember(internal,[start target]))=[];
% Normalize to one cAMP consumed. All shared internal states must cancel.
A=full([N(internal,:);N(start,:)]);rhs=[zeros(numel(internal),1);-1];
[v,solved]=positiveFlux(A,rhs);
if ~solved,reason="solverUnresolved";return;end
if isempty(v),return;end
net=N*v;if net(target)<=tol||abs(net(start)+1)>tol,return;end
feed=find(net < -tol&~ctx.currency);feed(feed==start)=[];
waste=find(net > tol&~ctx.currency);waste(waste==target)=[];
buffered=find(ctx.currency&any(a+b>tol,2));available=false(size(N,1),1);
available([feed;buffered])=true;without=substrateClosure(available,a,b);
if without(target),reason="seedIndependentRejected";return;end
available(start)=true;[with,fired]=substrateClosure(available,a,b);
if ~with(target)||~all(fired),reason="startupRejected";return;end
candidate.flux=v(:)';candidate.net=net;candidate.feed=feed;candidate.waste=waste;
candidate.internal=internal;candidate.buffered=buffered;reason="accepted";
end

function T = step2Block0AlternativeTable(rows)
T=cellTable(rows,{'ParentID','AlternativeID','ReactionCount','MetaboliteTrace','ReactionIDs', ...
    'Directions','FluxMultipliers','NetEquation','AdditionalFeedIDs','WasteIDs','InternalIDs', ...
    'ECQualified','Status'}, ...
    {'double','double','double','string','string','string','string','string','string','string','string','logical','string'});
end

function rows = step2Block0AlternativeRows(parent,alternatives,ctx)
rows=cell(0,13);
for k=1:numel(alternatives)
    p=alternatives(k);rows(end+1,:)={parent,k,p.distance, ...
        joinText(string(ctx.model.mets(p.mets)),' -> '),joinText(string(ctx.model.rxns(ctx.O.rxn(p.o))),';'), ...
        joinText(string(ctx.O.direction(p.o)),';'),joinText(compose('%.9g',p.flux),';'), ...
        equation(max(-p.net,0),max(p.net,0),ctx.model),joinText(string(ctx.model.mets(p.feed)),';'), ...
        joinText(string(ctx.model.mets(p.waste)),';'),joinText(string(ctx.model.mets(p.internal)),';'), ...
        p.ECComplete&&~contains(p.status,"not_EC"),p.status}; %#ok<AGROW>
end
end

function path = step2DiagnosticGraphPath(target,ctx,G)
% A diagnostic only: this ignores side-participant isolation and EC uniqueness.
path=emptyStep2Path(target);starts=find(ctx.base=="camp"&ctx.comp==ctx.comp(target));
for start=reshape(starts,1,[])
    seen=false(numel(ctx.model.mets),1);seen(start)=true;q=start;head=1;prev=zeros(size(seen));edge=prev;
    while head<=numel(q)
        u=q(head);head=head+1;
        for e=reshape(G.byFrom{u},1,[])
            v=G.to(e);if seen(v),continue;end
            seen(v)=true;prev(v)=u;edge(v)=G.o(e);q(end+1)=v; %#ok<AGROW>
        end
        if seen(target),break;end
    end
    if ~seen(target),continue;end
    mets=target;os=[];v=target;
    while v~=start,os=[edge(v) os];v=prev(v);mets=[v mets];end %#ok<AGROW>
    path.found=true;path.mets=mets;path.o=os;path.distance=numel(os);
    path.ECComplete=all(G.ecKnown(os));path.status="diagnostic_graph_path_NOT_an_isolated_solution";return;
end
end

function [step2,ctx] = searchStep2Circuit(ctx,cfg)
[ctx,compartmentAudit]=prepareStep2NeuronCompartment(ctx,cfg);
G=buildStep2Graph(ctx,cfg);G0=buildStep2Block0Graph(ctx,cfg);
parentRows=cell(0,12);entryRows=cell(0,10);
block0Rows=cell(0,13);block0CoverageRows=cell(0,16);block0AuditRows=cell(0,6);
completeRows=cell(0,9);neuronRows=cell(0,17);transportAudit=step2ExportAuditTable(cell(0,13));
chainInventory=step2ConsumingChainTable(cell(0,19));chainInventory.ParentID=zeros(0,1);chainCoverageRows=cell(0,10);
adapterAudit=step2AdapterAuditTable(cell(0,11));adapterCoverageRows=cell(0,10);
allInterfaces=emptyStep2Exports();interfaceParents=zeros(0,1);
details=struct('id',{},'w',{},'x',{},'b',{},'c',{},'rx',{},'parts',{},'controlRx',{}, ...
    'block0',{},'block0Alternatives',{},'diagnostic',{},'neuronStatus',{},'previewNeuron',{},'previewTransport',{});
motifs=struct('parentID',{},'n1',{},'transport',{},'block0',{},'controlRx',{},'assignedECs',{});
stats=struct('parentsTested',0,'parentLimitReached',false,'neuron1Rails',0, ...
    'neuronCandidatesTested',0,'neuronBudgetLimitedCalls',0,'neuronTestLimitReached',false,'completeLimitReached',false, ...
    'parentBlock0Found',0,'parentBlock0StructuralOnly',0,'block0BudgetLimitedCalls',0,'block0Queries',0, ...
    'transportTests',0,'exportLinks',0,'exportBudgetLimitedCalls',0,'qualifiedNeuronCandidates',0, ...
    'adapterStates',0,'adapterExpansions',0,'adapterBudgetLimitedCalls',0, ...
    'neuronChainStates',0,'neuronDiagnosticStates',0,'neuronDiagnosticLimitCalls',0, ...
    'timeLimitReached',false,'sharedBlock0Searches',0,'sharedBlock0Reused',0,'block0FullReruns',0,'elapsedSeconds',0);
cache=containers.Map('KeyType','char','ValueType','any');stop=false;
shared=containers.Map('KeyType','char','ValueType','any');
% Parse every reaction's EC and gene sets once (used in all inner loops).
nR=numel(ctx.model.rxns);ctx.ecTokenSets=cell(nR,1);ctx.geneTokenSets=cell(nR,1);
for j=1:nR
    ctx.ecTokenSets{j}=completeECTokens(ctx.E.codes(j));
    [~,g,~]=decodeGPR(string(ctx.model.grRules{j}),string(ctx.model.rules{j}),string(ctx.model.genes));
    g=split(string(g),';');ctx.geneTokenSets{j}=g(g~="");
end
progressPath=fullfile(cfg.outputDir,'Step2_Progress.txt');clock0=tic;
step2Progress(progressPath,sprintf('Step 2 started %s: %d parent circuits available, cap %g, time limit %g s.', ...
    datestr(now,31),height(ctx.Y),cfg.step2MaxParents,cfg.step2MaxSeconds),true);
if cfg.connectStep2&&cfg.counterMode=="consume_I"
% Visit the least expensive parents first (fewest canceller reactions), so
% that parent, neuron-test and time caps cut the expensive ones.
parentCost=zeros(height(ctx.Y),1);
for w=1:height(ctx.Y)
    parentCost(w)=step2CancellerReactionCount(ctx,find(ctx.C.ConnectedID==ctx.Y.ConnectedID(w),1));
end
[~,parentOrder]=sortrows([parentCost (1:height(ctx.Y))']);
for w=reshape(parentOrder,1,[])
    if stats.parentsTested>=cfg.step2MaxParents,stats.parentLimitReached=true;break;end
    if toc(clock0)>cfg.step2MaxSeconds
        stats.timeLimitReached=true;
        step2Progress(progressPath,sprintf('Time limit reached after %d parents; writing outputs.',stats.parentsTested),false);
        break;
    end
    stats.parentsTested=stats.parentsTested+1;
    [rx,parts,control,b,x,c]=step2Prefix(ctx,w);target=ctx.M.InputIndex(b);
    key=sprintf('%d,',[target sort(rx) 0 sort(control)]);
    if ~cfg.step2SearchBlock0
        b0=emptyStep2Path(target);b0.status="Block0_not_searched_step2SearchBlock0_false";
        b0alts=repmat(b0,1,0);passes=struct('mode',{},'result',{},'retained',{},'audit',{});
    elseif isKey(cache,key),saved=cache(key);b0=saved.best;b0alts=saved.alternatives;passes=saved.passes;
    elseif ~cfg.block0RequireEC
        % Block 0 without EC requirement: one structural pass only.
        stats.block0Queries=stats.block0Queries+1;
        [baseRx,baseParts]=step2PrefixWithoutCounter(ctx,w);
        [b0,b0alts,b0audit,how]=step2SharedBlock0(target,rx,parts,[],baseRx,baseParts,[], ...
            ctx,G0,cfg,false,shared);
        stats=step2CountShared(stats,how);
        passes=struct('mode',"structural_only",'result',b0,'retained',numel(b0alts),'audit',b0audit);
        saved=struct('best',b0,'alternatives',b0alts,'passes',passes);cache(key)=saved;
    else
        stats.block0Queries=stats.block0Queries+1;
        [baseRx,baseParts,baseControl]=step2PrefixWithoutCounter(ctx,w);
        [b0,b0alts,b0audit,how]=step2SharedBlock0(target,rx,parts,control,baseRx,baseParts,baseControl, ...
            ctx,G0,cfg,true,shared);
        stats=step2CountShared(stats,how);
        passes=struct('mode',"EC_qualified",'result',b0,'retained',numel(b0alts),'audit',b0audit);
        if ~b0.found
            [raw,rawalts,rawaudit,how]=step2SharedBlock0(target,rx,parts,[],baseRx,baseParts,[], ...
                ctx,G0,cfg,false,shared);
            stats=step2CountShared(stats,how);
            passes(end+1)=struct('mode',"structural_only",'result',raw,'retained',numel(rawalts),'audit',rawaudit);
            if raw.found,b0=raw;b0alts=rawalts;end
        end
        saved=struct('best',b0,'alternatives',b0alts,'passes',passes);cache(key)=saved;
    end
    parentID=ctx.Y.CounterWiringID(w);
    block0Rows=[block0Rows;step2Block0AlternativeRows(parentID,b0alts,ctx)]; %#ok<AGROW>
    budgetLimited=false;
    for pass=reshape(passes,1,[])
        p=pass.result;status=p.searchStatus;if status=="not_started",status=p.status;end
        block0CoverageRows(end+1,:)={parentID,pass.mode,status,p.states,p.expansions,p.topologyReachable, ...
            pass.retained,p.moduleTests,p.balanceRejected,p.startupRejected,p.seedIndependentRejected,p.solverUnresolved, ...
            p.stateLimitReached,p.expansionLimitReached,p.depthLimitReached,p.alternativeLimitReached}; %#ok<AGROW>
        a=table2cell(pass.audit);
        block0AuditRows=[block0AuditRows;[repmat({parentID,pass.mode},size(a,1),1) a]]; %#ok<AGROW>
        budgetLimited=budgetLimited||p.stateLimitReached||p.expansionLimitReached;
    end
    stats.parentBlock0Found=stats.parentBlock0Found+double(b0.found&&b0.ECComplete&&~contains(b0.status,"not_EC"));
    stats.parentBlock0StructuralOnly=stats.parentBlock0StructuralOnly+double(b0.found&&contains(b0.status,"not_EC"));
    stats.block0BudgetLimitedCalls=stats.block0BudgetLimitedCalls+double(budgetLimited);
    f=ctx.X.FormationID(x);counter=ctx.F.C_Index(f);seed=ctx.F.I_Index(f);
    if cfg.neuronPlacement=="outside"
        [links,ta,aa,ts]=findStep2AdapterExports(ctx.Y.CounterWiringID(w),counter,seed,rx,parts,control,ctx,G,cfg);
        transportAudit=[transportAudit;ta];adapterAudit=[adapterAudit;aa]; %#ok<AGROW>
        allInterfaces=[allInterfaces links]; %#ok<AGROW>
        interfaceParents=[interfaceParents;repmat(ctx.Y.CounterWiringID(w),numel(links),1)]; %#ok<AGROW>
        stats.adapterStates=stats.adapterStates+ts.states;stats.adapterExpansions=stats.adapterExpansions+ts.expansions;
        stats.adapterBudgetLimitedCalls=stats.adapterBudgetLimitedCalls+double(ts.stateLimitReached||ts.expansionLimitReached);
        adapterCoverageRows(end+1,:)={ctx.Y.CounterWiringID(w),ts.states,ts.expansions,numel(links), ...
            ts.maxDepthReached,ts.stateLimitReached,ts.expansionLimitReached,ts.depthLimitReached, ...
            ts.testLimitReached,ts.linkLimitReached}; %#ok<AGROW>
        stats.transportTests=stats.transportTests+ts.tests;stats.exportLinks=stats.exportLinks+numel(links);
        stats.exportBudgetLimitedCalls=stats.exportBudgetLimitedCalls+double(ts.testLimitReached||ts.linkLimitReached);
    else
        links=makeStep2ExportLink(0,0,counter,counter,[],0,false);
    end
    neuronTests=0;parentRails=0;preview=[];previewTransport=[];
    if ~isempty(links),previewTransport=links(1);end
    neuronStatus="no_exported_input_passed_current_policies";
    if ~isempty(links),neuronStatus="no_single_neuron_rail_passed_current_policies";end
    for t=1:numel(links)
        link=links(t);transportRx=link.rx(link.rx>0);prefix=unique([rx link.adapterRx transportRx],'stable');
        prefixParts=union(parts,link.parts);
        entryRows=[entryRows;step2EntryRows(ctx.Y.CounterWiringID(w),link.input,seed,prefix,prefixParts,ctx,G,cfg)]; %#ok<AGROW>
        [n1,s1,chains]=findStep2NeuronRails(link.input,seed,prefix,prefixParts,ctx,G,cfg,[control link.adapterRx]);
        chains.ParentID=repmat(ctx.Y.CounterWiringID(w),height(chains),1);chainInventory=[chainInventory;chains]; %#ok<AGROW>
        stats.neuronChainStates=stats.neuronChainStates+s1.states;
        stats.neuronDiagnosticStates=stats.neuronDiagnosticStates+s1.diagnosticStates;
        stats.neuronDiagnosticLimitCalls=stats.neuronDiagnosticLimitCalls+double(s1.diagnosticLimitReached);
        chainCoverageRows(end+1,:)={ctx.Y.CounterWiringID(w),string(ctx.model.mets{link.input}),s1.states, ...
            s1.tests,s1.stateLimitReached,s1.testLimitReached,s1.depthLimitReached, ...
            s1.diagnosticStates,s1.diagnosticTests,s1.diagnosticLimitReached}; %#ok<AGROW>
        stats.neuron1Rails=stats.neuron1Rails+numel(n1);parentRails=parentRails+numel(n1);
        stats.neuronBudgetLimitedCalls=stats.neuronBudgetLimitedCalls+double(s1.testLimitReached||s1.railLimitReached||s1.stateLimitReached);
        for a=1:numel(n1)
            if stats.neuronCandidatesTested>=cfg.step2MaxNeuronTests
                stats.neuronTestLimitReached=true;stop=true;break;
            end
            stats.neuronCandidatesTested=stats.neuronCandidatesTested+1;neuronTests=neuronTests+1;
            ctrl=[control link.adapterRx n1(a).rx];[~,ok]=step2AssignEC(ctrl,ctx,cfg);
            if ~ok||~step2TransportGenesOK(link.rx,ctrl,ctx,cfg),continue;end
            stats.qualifiedNeuronCandidates=stats.qualifiedNeuronCandidates+1;
            if isempty(preview),preview=n1(a);previewTransport=link;end
            if ~startsWith(neuronStatus,"complete_")
                neuronStatus="neuron_found_Block0_or_full_circuit_pending";
                if ~cfg.step2SearchBlock0,neuronStatus="neuron_found_Block0_not_searched";end
            end
            exportID="";if link.rx>0,exportID=string(ctx.model.rxns{link.rx});end
            neuronRows(end+1,:)={size(neuronRows,1)+1,ctx.Y.CounterWiringID(w), ...
                string(ctx.model.mets{counter}),string(ctx.model.mets{link.input}),ctx.comp(link.input), ...
                exportID,link.direction,link.reverseAllowed, ...
                joinText(string(ctx.model.rxns(n1(a).rx)),';'),string(ctx.model.mets{n1(a).mid}), ...
                string(ctx.model.mets{n1(a).output}),"neuron_only_candidate_not_a_complete_circuit", ...
                string(ctx.model.mets{link.cargo}),numel(link.adapterO),joinText(string(ctx.model.rxns(link.adapterRx)),';'), ...
                numel(n1(a).conversionO),joinText(string(ctx.model.mets(n1(a).mets)),' -> ')}; %#ok<AGROW>
            if ~cfg.step2SearchBlock0,continue;end % neuron-only mode: no per-neuron Block 0
            fullRx=unique([prefix n1(a).rx],'stable');fullParts=union(prefixParts,n1(a).parts);
            blockControl=[];if cfg.block0RequireEC,blockControl=ctrl;end
            [full0,fullAlternatives]=findStep2Block0(target,fullRx,fullParts,blockControl,ctx,G0,cfg,cfg.block0RequireEC);
            if ~full0.found,continue;end
            % A transport-gene conflict in the first retained alternative
            % must not hide a later compatible Block 0 support. Without
            % block0RequireEC, Block 0 reactions stay out of the EC assignment.
            ok=false;
            for z=1:numel(fullAlternatives)
                chosen=fullAlternatives(z);controlled=[ctrl ctx.O.rxn(chosen.o)];
                ecSet=ctrl;if cfg.block0RequireEC,ecSet=controlled;end
                [assigned,compatible]=step2AssignEC(ecSet,ctx,cfg);
                if compatible&&step2TransportGenesOK(link.rx,controlled,ctx,cfg)
                    full0=chosen;ok=true;break;
                end
            end
            if ~ok,continue;end
            neuronStatus="complete_single_neuron_circuit_found_controls_and_kinetics_unverified";
            id=numel(motifs)+1;motifs(id)=struct('parentID',numel(details)+1,'n1',n1(a),'transport',link, ...
                'block0',full0,'controlRx',ctrl,'assignedECs',assigned); %#ok<AGROW>
            completeRows(end+1,:)={id,ctx.Y.CounterWiringID(w),ctx.Y.ModuleID(w),ctx.Y.ConnectedID(w), ...
                ctx.X.CounterConnectionID(x),full0.distance,string(ctx.model.mets{n1(a).output}), ...
                ctx.Y.UpstreamSFInterfacesMatched(w),neuronStatus}; %#ok<AGROW>
            if numel(motifs)>=cfg.step2MaxCompleteMotifs,stats.completeLimitReached=true;stop=true;break;end
        end
        if stop,break;end
    end
    step2Progress(progressPath,sprintf('parent %d/%d (wiring %d): neuron rails %d, tested %d -> %s | Block 0 %s | %.0f s', ...
        stats.parentsTested,min(height(ctx.Y),cfg.step2MaxParents),ctx.Y.CounterWiringID(w), ...
        parentRails,neuronTests,neuronStatus,b0.status,toc(clock0)),false);
    id=numel(details)+1;details(id)=struct('id',ctx.Y.CounterWiringID(w),'w',w,'x',x,'b',b,'c',c, ...
        'rx',rx,'parts',parts,'controlRx',control,'block0',b0,'block0Alternatives',b0alts, ...
        'diagnostic',step2DiagnosticGraphPath(target,ctx,G0),'neuronStatus',neuronStatus, ...
        'previewNeuron',preview,'previewTransport',previewTransport); %#ok<AGROW>
    parentRows(end+1,:)={ctx.Y.CounterWiringID(w),ctx.X.CounterConnectionID(x),ctx.Y.ModuleID(w), ...
        ctx.Y.ConnectedID(w),ctx.M.InsideCompartment(b),ctx.X.C_ID(x),ctx.M.InputID(b), ...
        b0.status,b0.distance,parentRails,neuronTests,neuronStatus}; %#ok<AGROW>
    if stop,break;end
end
end
stats.elapsedSeconds=toc(clock0);
step2Progress(progressPath,sprintf('Step 2 finished: %d parents in %.0f s (time limit reached=%d).', ...
    stats.parentsTested,stats.elapsedSeconds,stats.timeLimitReached),false);
step2=struct('stats',stats,'details',details,'motifs',motifs,'reactionPool',ctx.O);
step2.block0Alternatives=step2Block0AlternativeTable(block0Rows);
step2.block0Coverage=cellTable(block0CoverageRows,{'ParentID','Pass','Status','PathStates','Expansions', ...
    'TopologyReachable','RetainedAlternatives','ModuleTests','BalanceRejected','StartupRejected', ...
    'SeedIndependentRejected','SolverUnresolved','StateLimitReached','ExpansionLimitReached', ...
    'DepthBoundReached','AlternativeLimitReached'}, ...
    {'double','string','string','double','double','logical','double','double','double','double','double','double', ...
    'logical','logical','logical','logical'});
step2.block0FilterAudit=cellTable(block0AuditRows,{'ParentID','Pass','Reason','OrientationCount', ...
    'SharedStateIDs','ExampleReactionIDs'},{'double','string','string','double','string','string'});
step2.parents=cellTable(parentRows,{'ParentID','CounterConnectionID','ModuleID','ConnectedID','InsideCompartment', ...
    'CounterMetaboliteID','InputID','Block0Status','Block0Distance','Neuron1RailCount','NeuronCandidateCount','Status'}, ...
    {'double','double','double','double','string','string','string','string','double','double','double','string'});
step2.entryAudit=cellTable(entryRows,{'ParentID','InputID','ReactionID','Direction','ECResolved','WouldRegenerateI', ...
    'SharedNoncurrencyOutsideInput','Status','NeuronCompartment','IProductIDs'}, ...
    {'double','string','string','double','logical','logical','string','string','string','string'});
step2.transportAudit=transportAudit;step2.compartmentAudit=compartmentAudit;
ctx.currentExportAudit=transportAudit;ctx.currentAdapterAudit=adapterAudit;
ctx.currentChainInventory=chainInventory;step2.consumingChains=chainInventory;
step2.chainCoverage=cellTable(chainCoverageRows,{'ParentID','InputID','PathStates','Tests', ...
    'StateLimitReached','TestLimitReached','DepthBoundReached','DiagnosticStates','DiagnosticTests','DiagnosticLimitReached'}, ...
    {'double','string','double','double','logical','logical','logical','double','double','logical'});
step2.adapterAudit=adapterAudit;step2.adapterLinks=describeStep2AdapterLinks(allInterfaces,interfaceParents,ctx);
step2.adapterCoverage=cellTable(adapterCoverageRows,{'ParentID','PathStates','ReactionExpansions','ExportLinks', ...
    'MaxAdapterDepthReached','StateLimitReached','ExpansionLimitReached','DepthBoundReached', ...
    'TransportBudgetReached','LinkLimitReached'}, ...
    {'double','double','double','double','double','logical','logical','logical','logical','logical'});
step2.neuronCandidates=cellTable(neuronRows,{'NeuronCandidateID','ParentID','CounterID','NeuronInputID', ...
    'NeuronCompartment','ExportReactionID','ExportDirection','ExportReverseAllowed','NeuronReactionIDs', ...
    'NeuronMidID','NeuronOutputID','Status','ExportCargoID','AdapterReactionCount','AdapterReactionIDs', ...
    'NeuronConversionCount','NeuronMetabolitePath'}, ...
    {'double','double','string','string','string','string','double','logical','string','string','string','string', ...
    'string','double','string','double','string'});
step2.complete=cellTable(completeRows,{'CandidateID','ParentID','ModuleID','ConnectedID','CounterConnectionID', ...
    'Block0Distance','OutputID','UpstreamSFInterfacesMatched','Status'}, ...
    {'double','double','double','double','double','double','string','logical','string'});
end

function writeUpstreamCoverage(fid,ctx,cfg)
step2Section(fid,'UPSTREAM SEARCH DIAGNOSTICS');
report(fid,sprintf('Production mode=%s; maximum module reactions=%d; I-state budget=%g.', ...
    cfg.productionMode,cfg.maxCoreReactions,cfg.maxStatesPerI));
report(fid,sprintf('Stages enabled: removal=%d; branches=%d; counter=%d; Step2=%d. Required supplied import-gate evidence=%d.', ...
    cfg.connectRemovalTerms,cfg.connectBranches,cfg.connectCounter,cfg.connectStep2,cfg.requireGateEvidence));
if ~isfield(ctx,'upstreamCoverage')
    report(fid,'Upstream stage counts are unavailable in this context.');return;
end
U=ctx.upstreamCoverage;
for k=1:height(U)
    report(fid,sprintf('  %s: %d; enabled=%d; inspect %s.',U.Stage(k),U.Count(k),U.Enabled(k),U.DiagnosticFile(k)));
end
first=find(U.Enabled & U.Count==0,1);
if ~isempty(first)
    report(fid,['First empty required upstream stage: ' char(U.Stage(first)) '.']);
    report(fid,['Inspect ' char(U.DiagnosticFile(first)) ' and its preceding stage before diagnosing neuron connectivity.']);
end
if isfield(ctx,'productionCoverage')
    P=ctx.productionCoverage;
    report(fid,sprintf('Production coverage: budget reached for %d I pools; solver-unresolved closures=%d. See 02_Search_Coverage.csv.', ...
        nnz(P.StateBudgetReached),sum(P.SolverUnresolved)));
    if ismember('MinimumProvenWithinSearchClass',P.Properties.VariableNames)
        report(fid,sprintf('Production minimum sizes proven within the enumerated search class for %d I pools. See Minimum_Production_Report.txt.', ...
            nnz(P.MinimumProvenWithinSearchClass)));
        report(fid,'A minimum applies to production alone. Longer modules may be needed to satisfy downstream isolation and wiring.');
    end
end
if cfg.productionMode=="single_reaction" && isfield(ctx,'directProductionAudit')
    A=ctx.directProductionAudit;
    report(fid,'Direct-only production filter: 02_Direct_Production_Filter_Audit.csv.');
    for status=reshape(unique(A.Status,'stable'),1,[])
        report(fid,sprintf('  %s: %d',status,nnz(A.Status==status)));
    end
    report(fid,'Direct-only mode excludes multi-step realizations. Set productionMode="multi_reaction" to search reaction modules.');
end
if ~cfg.connectStep2
    report(fid,'Step 2 was disabled by connectStep2=false.');
elseif cfg.counterMode~="consume_I"
    report(fid,'Step 2 requires counterMode="consume_I".');
elseif isempty(ctx.Y)
    report(fid,'Neuron and Block 0 candidate searches were not entered: there are no upstream branch-canceller-counter parents.');
    report(fid,'This report therefore does not establish a neuron-connectivity failure.');
end
end

function writeStep2Report(path,S,ctx,cfg)
fid=fopen(path,'w');assert(fid>0,'Anderson:Step2Report','Cannot write Step_2_Results.txt.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
report(fid,'BLOCK 0 -> FAST/SLOW BRANCHES -> ANDERSON CANCELLER -> CONSUMING COUNTER -> ONE CONSUMING NEURON');
report(fid,'Code revision: direct_intracellular_consuming_neuron_chains_2026_11_16; cell-free merge and cancellerI 2026-11-19.');
report(fid,sprintf('cancellerI=%s; counter/neuron/Block 0 may use embedded iMM904 ECs=%d.', ...
    joinText(string(cfg.cancellerI),','),cfg.counterUseEmbeddediMM904EC));
if cfg.cellFreeMerge
    report(fid,'CELL-FREE MODE: one vessel; organelle copies merged per block; intracellular transports removed; identical isoenzymes merged.');
end
if ~isfield(ctx,'currentChainInventory'),ctx.currentChainInventory=S.consumingChains;end
report(fid,'Direction policy: model-bound orientations; reverse steps included only when allowReverse and the lower bound permit them.');
report(fid,sprintf('allowReverse=%d; counterMode=%s; EC map is partial and source/signature checked.',cfg.allowReverse,cfg.counterMode));
if cfg.neuronPlacement=="outside"
    report(fid,['Neuron compartment: ' char(step2NeuronCompartment(cfg)) '. Existing extracellular reactions only; no reaction relocation is assumed.']);
    report(fid,'Neuron layout: m_c[inside] -> zero or more adapter steps -> m_exp[inside] --r_T_N1 export--> m_Nin_1[outside] -> m_Nf_1[outside] -> m_Nd_1[outside], plus r_Nr_1 removal.');
    report(fid,'m_exp and m_Nin_1 are the same chemical in different states; m_c may be a different chemical. The first adapter consumes m_c; a zero-step adapter exports m_c directly.');
    report(fid,sprintf('Adapter limits per parent: conversions=%g; path states=%g; reaction expansions=%g. Shorter paths are examined first; all co-participants are checked.', ...
        cfg.neuronAdapterMaxDepth,cfg.neuronAdapterMaxStates,cfg.neuronAdapterMaxExpansions));
    report(fid,'Every adapter needs EC evidence and distinct controlled-role EC/gene assignments. Any adapter participation in an upstream I pool is rejected, including GSH regeneration in its first reaction.');
    report(fid,'Extracellular GSH is not intracellular I, but GSH shared with the selected upstream importer is rejected as feedback.');
    report(fid,'Transport does not require an EC; optional GPR requirements and all allowed directions are audited. A reversible exporter is explicitly flagged.');
    ids=unique(S.compartmentAudit.ReactionID(S.compartmentAudit.Status=="outside_orientation_admitted"));
    [~,j]=ismember(ids,string(ctx.model.rxns));j=j(j>0);
    report(fid,sprintf('Recorded extracellular conversion reaction IDs=%d; with complete EC evidence=%d; target-compartment species=%d.', ...
        numel(j),nnz(ctx.E.codes(j)~=""),nnz(ctx.comp==step2NeuronCompartment(cfg))));
    if ~any(ctx.E.codes(j)~="")
        report(fid,'EXTERNAL NEURON BOTTLENECK: no recorded extracellular conversion has complete EC evidence. Adapter/export interfaces can still be found, but an EC-qualified external neuron cannot pass the present annotation policy.');
    end
else
    report(fid,'Neuron layout: m_c = m_Nin_1 is consumed directly inside the counter compartment; no exporter or external input pool is required.');
    report(fid,'The first product is m_Nf_1; optional intermediate products are m_Nmid_1_XX; the last product is m_Nd_1 = m_out. An independent r_Nr_1 removes that output.');
end
report(fid,sprintf('Neuron conversion-chain depth=%g..%g; state cap=%g; test cap=%g; independent output removal is required.', ...
    cfg.neuronChainMinDepth,cfg.neuronChainMaxDepth,cfg.neuronChainMaxStates,cfg.neuronMaxTestsPerRail));
report(fid,'Search scope: one compartment, existing bound-permitted gross reactions, EC-annotated steps, noncurrency connectors with one primitive consumed/produced copy; all side participants are retained.');
report(fid,'Chemical consuming paths can exist but fail isolation. I interaction, upstream reaction reuse, side-participant overlap and EC/gene conflicts remain excluded from valid neuron candidates.');
report(fid,sprintf('Rejected-chain diagnostics enabled=%d; separate diagnostic state/test caps=%g/%g. Diagnostic chains NEVER count as isolated neurons or complete circuits.', ...
    cfg.neuronDiagnoseRejectedChains,cfg.neuronDiagnosticMaxStates,cfg.neuronDiagnosticMaxTests));
report(fid,'m_out = m_Nd_1. Adapter/export interfaces and neuron-only candidates are retained independently of Block 0 completion.');
report(fid,'The neuron input is consumed as a substrate; transporter or enzyme activation is not inferred.');
report(fid,'Block 0 always starts at exact intracellular cAMP and stops at the exact m_in used by both branches.');
report(fid,'All substrates and products are checked for cross-talk; currency may be buffered but cannot carry a path connection.');
report(fid,'Block 0 shares only m_in with downstream noncurrency chemistry and reuses no downstream reaction. G3P and every other selected downstream participant remain forbidden, even if net-cancelled.');
report(fid,'Block 0 uses a separate net-consumption/production graph that permits nonunit coefficients and balanced sharing of its own internal side intermediates. Neuron matching is unchanged.');
report(fid,'Each Block 0 module has positive stoichiometric weights normalized to one cAMP consumed; shared internal noncurrency intermediates cancel. All selected reactions must be substrate-reachable from cAMP, additional feeds and buffered currencies.');
report(fid,'A Block 0 output reachable without cAMP in the selected module is rejected. These are local stoichiometry/reachability checks, not kinetic or atom-tracing validation.');
report(fid,'Scope: bounded simple connector traces with balanced supports; arbitrary branched reaction sets are not exhaustively enumerated.');
report(fid,'16_Block0_Alternatives.csv retains parent-stage alternatives even without a neuron; full circuits recheck Block 0 against the selected neuron and transport.');
report(fid,'16_Block0_Coverage.csv separates qualified/structural passes, module rejections, unresolved flux solves and every cap. 16_Block0_Filter_Audit.csv summarizes first-match reaction exclusions with example IDs.');
report(fid,'Accepted neurons share only their intended input with preceding selected chemistry; return to selected upstream I pools is excluded. Rejected diagnostic paths are labelled separately.');
report(fid,'New Block 0, adapter and neuron roles require complete ECs. Existing canceller roles retain any unresolved annotations explicitly shown.');
report(fid,sprintf('Distinct EC assignment for controlled branch/counter/neuron/Block0 roles=%d; known shared-gene rejection=%d.',cfg.step2DistinctECs,cfg.branchRejectSharedKnownGenes));
report(fid,'An EC does not establish a protein identity, KM, kcat, first-order kinetics or independent biochemical control.');
report(fid,'Upstream S/F gate interfaces remain hypotheses unless their supplied evidence matches; a structural extension does not validate those controls.');
report(fid,sprintf('Block 0 caps: depth=%g; path states=%g; expansions=%g. Neuron caps: tests/call=%g; rails/call=%g; candidate tests=%g.', ...
    cfg.block0MaxDepth,cfg.block0MaxStates,cfg.block0MaxExpansions,cfg.neuronMaxTestsPerRail,cfg.neuronMaxRails,cfg.step2MaxNeuronTests));
report(fid,sprintf('Block 0 alternatives retained per search=%g; maximum gross participants/reaction=%g.', ...
    cfg.block0MaxAlternatives,cfg.block0MaxParticipants));
report(fid,'A reached budget is unresolved search coverage, not evidence of biochemical impossibility.');
writeUpstreamCoverage(fid,ctx,cfg);
step2Section(fid,'SEARCH COUNTS: FOUND VERSUS PRINTED');
report(fid,sprintf('Branch-canceller-counter parents available: %d; examined: %d.',height(ctx.Y),S.stats.parentsTested));
if cfg.neuronPlacement=="outside"
report(fid,sprintf('Counter-export orientations tested: %d; accepted export links: %d; budget-limited calls: %d.', ...
    S.stats.transportTests,S.stats.exportLinks,S.stats.exportBudgetLimitedCalls));
report(fid,sprintf('Intracellular adapter path states=%d; reaction expansions=%d; budget-limited calls=%d; conversion/export interfaces retained=%d.', ...
    S.stats.adapterStates,S.stats.adapterExpansions,S.stats.adapterBudgetLimitedCalls,height(S.adapterLinks)));
report(fid,'See 14_Counter_Adapter_Export_Links.csv for interfaces even when no external neuron or Block 0 is available.');
end
report(fid,sprintf('EC-qualified isolated neuron candidates, without requiring Block 0: %d. See 15_Transported_Neuron_Candidates.csv.', ...
    S.stats.qualifiedNeuronCandidates));
report(fid,sprintf('Isolated-search path states=%d; diagnostic path states=%d; diagnostic budget-limited calls=%d; consuming-path inventory rows=%d.', ...
    S.stats.neuronChainStates,S.stats.neuronDiagnosticStates,S.stats.neuronDiagnosticLimitCalls,height(S.consumingChains)));
report(fid,'15_Consuming_Neuron_Chains.csv separates clean prefixes from flagged feedback/overlap diagnostics; 15_Isolated_Neuron_Candidates.csv contains only qualified neurons.');
report(fid,'15_Neuron_Chain_Coverage.csv records both search passes. A reached state, depth or test limit is incomplete coverage.');
report(fid,sprintf('Parents with an EC-qualified isolated cAMP-to-m_in connector: %d.',S.stats.parentBlock0Found));
report(fid,sprintf('Parents with an isolated structural connector awaiting EC qualification: %d.',S.stats.parentBlock0StructuralOnly));
report(fid,sprintf('Single-neuron rails passing current policies: %d; neuron candidates tested: %d.', ...
    S.stats.neuron1Rails,S.stats.neuronCandidatesTested));
report(fid,sprintf('Distinct complete motifs found: %d; complete motifs printed: %d.',numel(S.motifs),min(numel(S.motifs),cfg.step2PrintTop)));
report(fid,sprintf('Limits reached: parent=%d; neuron calls=%d; neuron candidate tests=%d; complete retention=%d; parent Block0 budget-limited calls=%d.', ...
    S.stats.parentLimitReached,S.stats.neuronBudgetLimitedCalls,S.stats.neuronTestLimitReached, ...
    S.stats.completeLimitReached,S.stats.block0BudgetLimitedCalls));
if isfield(S.stats,'timeLimitReached')
    report(fid,sprintf('Step 2 time limit reached=%d (step2MaxSeconds=%g); elapsed %.0f s.', ...
        S.stats.timeLimitReached,cfg.step2MaxSeconds,S.stats.elapsedSeconds));
    report(fid,sprintf('Block 0: shared upstream searches=%d; counter variants answered from a shared search=%d; full per-parent re-searches=%d.', ...
        S.stats.sharedBlock0Searches,S.stats.sharedBlock0Reused,S.stats.block0FullReruns));
end
report(fid,'Complete counts require one neuron and a Block 0 connector rechecked against all downstream reactions and EC assignments.');
report(fid,sprintf('Block 0 searched=%d; Block 0 requires EC=%d. With Block 0 off, look for neuronStatus neuron_found_Block0_not_searched.', ...
    cfg.step2SearchBlock0,cfg.block0RequireEC));
if isfield(S,'neuronCandidates')
    report(fid,sprintf('Parents with at least one passing neuron: %d; passing neuron candidates: %d.', ...
        numel(unique(S.neuronCandidates.ParentID)),height(S.neuronCandidates)));
end
report(fid,'Partial parents and rejected diagnostic paths never contribute to the complete count. Printing does not determine search limits.');
report(fid,'No whole-cell flux, atom-transfer, thermodynamic or kinetic validation is claimed. Additional co-substrates remain required.');
if cfg.neuronPlacement=="outside"
step2Section(fid,'COUNTER EXPORT REJECTIONS / ACCEPTED LINKS');
for status=reshape(unique(S.transportAudit.Status,'stable'),1,[])
    report(fid,sprintf('  %s: %d',status,nnz(S.transportAudit.Status==status)));
end
step2Section(fid,'INTRACELLULAR ADAPTER REJECTIONS / ACCEPTED STEPS');
for status=reshape(unique(S.adapterAudit.Status,'stable'),1,[])
    report(fid,sprintf('  %s: %d',status,nnz(S.adapterAudit.Status==status)));
end
report(fid,'See 14_Counter_Adapter_Audit.csv and 14_Counter_Adapter_Coverage.csv for path decisions and reached limits.');
report(fid,'See 14_Neuron_Export_Audit.csv for full transporter equations and 16_Neuron_Compartment_Audit.csv for available extracellular chemistry.');
end
step2Section(fid,'NEURON ENTRY REJECTIONS');
for status=reshape(unique(S.entryAudit.Status,'stable'),1,[])
    report(fid,sprintf('  %s: %d',status,nnz(S.entryAudit.Status==status)));
end
if ~isempty(S.motifs)
    % Least expensive first: fewest canceller reactions, then fewest Block 0
    % reactions, then fewest neuron reactions; ties keep the search order.
    nMotif=numel(S.motifs);cost=zeros(nMotif,3);
    for k=1:nMotif
        motif=S.motifs(k);detail=S.details(motif.parentID);
        cost(k,:)=[step2CancellerReactionCount(ctx,detail.c) motif.block0.distance numel(motif.n1.rx)];
    end
    [~,printOrder]=sortrows([cost (1:nMotif)']);
    report(fid,sprintf(['\nComplete motifs are printed least expensive first: canceller reactions, then Block 0 ' ...
        'reactions, then neuron reactions. Cheapest: %d canceller, %g Block 0, %d neuron reactions.'], ...
        cost(printOrder(1),1),cost(printOrder(1),2),cost(printOrder(1),3)));
    for k=1:min(nMotif,cfg.step2PrintTop)
        motif=S.motifs(printOrder(k));detail=S.details(motif.parentID);
        report(fid,sprintf('\nRANK %d | canceller reactions=%d | Block 0 reactions=%g | neuron reactions=%d', ...
            k,cost(printOrder(k),1),cost(printOrder(k),2),cost(printOrder(k),3)));
        writeStep2Candidate(fid,k,true,detail,motif,ctx,cfg);
    end
elseif cfg.step2PrintPartialWhenEmpty
    % Most advanced parents first (a passing neuron), then the least
    % expensive: fewest canceller reactions, then fewest Block 0 reactions
    % (Inf when none was found or searched); ties keep the search order.
    nDetail=numel(S.details);neuronStatus=strings(nDetail,1);cost=zeros(nDetail,2);
    for k=1:nDetail
        neuronStatus(k)=string(S.details(k).neuronStatus);
        cost(k,:)=[step2CancellerReactionCount(ctx,S.details(k).c) S.details(k).block0.distance];
    end
    printRank=2*ones(nDetail,1);printRank(startsWith(neuronStatus,"neuron_found"))=1;
    [~,printOrder]=sortrows([printRank cost (1:nDetail)']);
    report(fid,sprintf('\nNo complete motif passed. Printing %d of %d PARTIAL parents with all known reactions and explicit missing roles.', ...
        min(numel(S.details),cfg.step2PrintTop),numel(S.details)));
    report(fid,sprintf('Parents with a passing neuron are printed first: %d of %d parents have one.', ...
        nnz(printRank==1),numel(S.details)));
    report(fid,'Within each group: fewest canceller reactions first, then fewest Block 0 reactions.');
    for k=1:min(numel(S.details),cfg.step2PrintTop)
        report(fid,sprintf('\nRANK %d | neuron found=%d | canceller reactions=%d | Block 0 reactions=%g', ...
            k,printRank(printOrder(k))==1,cost(printOrder(k),1),cost(printOrder(k),2)));
        writeStep2Candidate(fid,k,false,S.details(printOrder(k)),[],ctx,cfg);
    end
end
end

function n = step2CancellerReactionCount(ctx,c)
% Distinct model reactions of one Anderson canceller: production core, I
% import, linear sink (plus its F-gated co-substrate import) and pair loss.
rx=counterCancellerSupport(c,ctx.C,ctx.D,ctx.cores,ctx.O,ctx.T,ctx.L,ctx.Q,ctx.R,ctx.P);
n=numel(unique(rx));
end

function writeStep2Candidate(fid,k,complete,d,motif,ctx,cfg)
b=d.b;x=d.x;c=d.c;f=ctx.X.FormationID(x);O=ctx.O;model=ctx.model;
signalNames={'m_in','m_Ff','m_Fd','m_Sf','m_Sd','I','m_c'};
signalIdx=[ctx.M.InputIndex(b),ctx.M.FastMidIndex(b),ctx.M.FastSignalIndex(b), ...
    ctx.M.SlowMidIndex(b),ctx.M.SlowSignalIndex(b),ctx.C.I_Index(c),ctx.F.C_Index(f)];
ctx.signalIndices=signalIdx;ctx.signalLabels=string(signalNames);
selectedNeuron=[];selectedTransport=[];
if complete
    selectedNeuron=motif.n1;selectedTransport=motif.transport;
else
    if isfield(d,'previewNeuron'),selectedNeuron=d.previewNeuron;end
    if isfield(d,'previewTransport'),selectedTransport=d.previewTransport;end
end
if ~isempty(selectedTransport)
    ctx.signalIndices(end+1)=selectedTransport.input;ctx.signalLabels(end+1)="m_Nin_1";
    if cfg.neuronPlacement=="outside"
        ctx.signalIndices(end+1)=selectedTransport.cargo;ctx.signalLabels(end+1)="m_exp";
        for n=2:numel(selectedTransport.adapterMets)-1
            ctx.signalIndices(end+1)=selectedTransport.adapterMets(n);
            ctx.signalLabels(end+1)="m_A"+string(n-1);
        end
    end
elseif cfg.neuronPlacement=="inside"
    ctx.signalIndices(end+1)=ctx.F.C_Index(f);ctx.signalLabels(end+1)="m_Nin_1";
end
if ~isempty(selectedNeuron)
    ctx.signalIndices=[ctx.signalIndices selectedNeuron.mid selectedNeuron.output selectedNeuron.output];
    ctx.signalLabels=[ctx.signalLabels "m_Nf_1" "m_Nd_1" "m_out"];
    for n=3:numel(selectedNeuron.mets)-1
        ctx.signalIndices(end+1)=selectedNeuron.mets(n);
        ctx.signalLabels(end+1)="m_Nmid_1_"+compose('%02d',n-2);
    end
end
title=sprintf('CANDIDATE %d',k);if ~complete,title=[title ' | PARTIAL PARENT, NOT A COMPLETE CIRCUIT'];end
step2Section(fid,title);
report(fid,sprintf('Parent=%d; module=%d; canceller=%d; counter=%d; block=%s.', ...
    d.id,ctx.M.ModuleID(b),ctx.C.ConnectedID(c),ctx.X.CounterConnectionID(x),ctx.M.InsideCompartment(b)));
report(fid,['Neuron status: ' char(d.neuronStatus)]);
report(fid,sprintf('S/F control interfaces matched: %d.',ctx.Y.UpstreamSFInterfacesMatched(d.w)));
report(fid,'SIGNAL METABOLITE MAP');
for n=1:numel(signalIdx),step2Signal(fid,signalNames{n},signalIdx(n),model);end
if ~isempty(selectedTransport)
    step2Signal(fid,'m_Nin_1',selectedTransport.input,model);
    if cfg.neuronPlacement=="outside"
        step2Signal(fid,'m_exp',selectedTransport.cargo,model);
        for n=2:numel(selectedTransport.adapterMets)-1
            step2Signal(fid,char("m_A"+string(n-1)),selectedTransport.adapterMets(n),model);
        end
    end
elseif cfg.neuronPlacement=="inside"
    step2Signal(fid,'m_Nin_1',ctx.F.C_Index(f),model);
else
    report(fid,'  m_exp, m_Nin_1: UNASSIGNED (no accepted conversion/export interface).');
end
if ~isempty(selectedNeuron)
    step2Signal(fid,'m_Nf_1',selectedNeuron.mid,model);step2Signal(fid,'m_Nd_1',selectedNeuron.output,model);
    step2Signal(fid,'m_out',selectedNeuron.output,model);
    for n=3:numel(selectedNeuron.mets)-1
        step2Signal(fid,sprintf('m_Nmid_1_%02d',n-2),selectedNeuron.mets(n),model);
    end
    if ~complete,report(fid,'  Neuron roles below are a NEURON-ONLY candidate; Block 0/full-circuit completion is pending.');end
else
    report(fid,'  m_Nf_1, m_Nd_1, m_out: UNASSIGNED (no accepted single-neuron extension).');
end
report(fid,'  r_Nr_1 is the neuron output-removal reaction; it does not introduce an m_Nr state.');
step2Section(fid,'BLOCK 0: cAMP TO m_in');
path=d.block0;if complete,path=motif.block0;end
report(fid,['  Connection status: ' char(path.status)]);
report(fid,['  Search termination: ' char(path.searchStatus)]);
report(fid,sprintf('  Limits reached: state=%d; expansion=%d; depth=%d; alternative retention=%d; unresolved flux solves=%d.', ...
    path.stateLimitReached,path.expansionLimitReached,path.depthLimitReached, ...
    path.alternativeLimitReached,path.solverUnresolved));
if path.found
    step2PrintBlock0Path(fid,path,ctx,'r_CAMP');
else
    report(fid,'  No isolated cAMP-to-m_in connector accepted within the configured search.');
    if d.diagnostic.found
        diagnostic=d.diagnostic;parts=[];
        for o=reshape(diagnostic.o,1,[]),parts=union(parts,find(O.R(:,o)+O.P(:,o)>1e-8));end
        shared=step2Overlap(parts,d.parts,ctx.M.InputIndex(b),ctx.currency);
        report(fid,'  REJECTED SHORTEST GRAPH ROUTE (diagnostic only; not a Block 0 solution):');
        report(fid,['  Path: ' char(joinText(string(model.mets(diagnostic.mets)),' -> '))]);
        report(fid,['  Shared downstream noncurrency species outside m_in: ' char(joinText(string(model.mets(shared)),'; '))]);
        report(fid,'  This route is diagnostic only. Cross-block overlap, failed internal balance/startup or EC/gene conflicts prevent acceptance.');
        for n=1:numel(diagnostic.o),step2PrintOrientation(fid,sprintf('r_CAMP_DIAGNOSTIC_%d',n),diagnostic.o(n),ctx);end
    end
end
report(fid,sprintf('  Parent-stage Block 0 alternatives retained: %d (before neuron rechecking).',numel(d.block0Alternatives)));
for z=1:numel(d.block0Alternatives)
    alternative=d.block0Alternatives(z);
    report(fid,sprintf('  Alternative %d | %s | %d reactions | %s',z,alternative.status,alternative.distance, ...
        joinText(string(model.rxns(O.rxn(alternative.o))),'; ')));
    report(fid,['    Connector trace: ' char(joinText(string(model.mets(alternative.mets)),' -> '))]);
end
report(fid,'  Full alternative reaction/direction/flux/feed/internal-state lists: 16_Block0_Alternatives.csv.');
step2Section(fid,'FAST AND SLOW BRANCHES (SIX REACTIONS)');
fast=ctx.rails(ctx.M.FastRailIndex(b));slow=ctx.rails(ctx.M.SlowRailIndex(b));
roles={'r_Ff','r_Fd','r_Fr','r_Sf','r_Sd','r_Sr'};os=[fast.o slow.o];
for n=1:6,step2PrintOrientation(fid,roles{n},os(n),ctx);end
step2Section(fid,'ANDERSON CANCELLER: PRODUCTION CORE AND TWO REMOVAL TERMS');
prod=ctx.C.ProductionCandidateID(c);core=ctx.cores(ctx.D.CoreIndex(prod));iRow=ctx.D.TransportRow(prod);
if iRow>0
    step2Section(fid,'EQ. BLOCK 1: I_OUT IMPORT + I + X -> 2 I + W');
else
    step2Section(fid,'EQ. BLOCK 1: I + X -> 2 I + W   (no I import; I(0) present in the vessel)');
end
report(fid,['  Proposed S signal: m_Sd = ' char(ctx.M.SlowSignalID(b)) '. The transporter-control link is not inferred from stoichiometry.']);
if iRow>0
    st=ctx.T(iRow,:);
    step2PrintReaction(fid,'r_T_S',st.ModelIndex,st.Direction,ctx);
    report(fid,'  S-dependent opening: supplied-interface status remains separate from cargo transport.');
else
    report(fid,'  No I import searched (requireIImport=false): eps = 0 and I(0) is assumed already present.');
    report(fid,'  S gates the co-substrate importer instead; see 06a/06b and Gated_Anderson_Report.txt.');
end
productionCtx=ctx;
for n=1:numel(core.feed)
    label="X_"+string(n);productionCtx.signalIndices(end+1)=core.feed(n);productionCtx.signalLabels(end+1)=label;
    step2Signal(fid,char(label),core.feed(n),model);
end
for n=1:numel(core.waste)
    label="W_"+string(n);productionCtx.signalIndices(end+1)=core.waste(n);productionCtx.signalLabels(end+1)=label;
    step2Signal(fid,char(label),core.waste(n),model);
end
effectiveR=max(-core.net,0);effectiveP=max(core.net,0);
effectiveR(core.seed)=full(O.R(core.seed,core.path)*core.flux(:));
effectiveP(core.seed)=full(O.P(core.seed,core.path)*core.flux(:));
report(fid,['  R1 effective module equation (non-I intermediates cancelled): ' ...
    char(step2NamedSide(effectiveR,productionCtx)+" -> "+step2NamedSide(effectiveP,productionCtx))]);
report(fid,sprintf('  Actual I accounting: %.12g consumed -> %.12g produced across %d reaction(s).', ...
    effectiveR(core.seed),effectiveP(core.seed),numel(core.path)));
if isfield(ctx,'productionCoverage') && ismember('SmallestFoundReactionCount',ctx.productionCoverage.Properties.VariableNames)
    coverage=ctx.productionCoverage;r=find(coverage.I_ID==string(model.mets{core.seed}),1);
    if ~isempty(r)
        report(fid,sprintf('  Smallest production module found for this I: %g steps; minimum proven within search class=%d; minimum-layer alternatives complete=%d.', ...
            coverage.SmallestFoundReactionCount(r),coverage.MinimumProvenWithinSearchClass(r),coverage.MinimumLayerComplete(r)));
        report(fid,['  Production-only minimum status: ' char(coverage.MinimumStatus(r)) '. Import adds one separate reaction.']);
    end
end
if numel(core.path)==1
    o=core.path;step2PrintOrientation(fid,'R1',o,productionCtx);
    report(fid,['  Primitive normalized R1: ' char(step2NamedSide(O.R(:,o)*core.flux,productionCtx)+" -> "+step2NamedSide(O.P(:,o)*core.flux,productionCtx))]);
else
    report(fid,'  R1 is the effective production role, implemented by the following recorded reactions.');
    for n=1:numel(core.path)
        report(fid,sprintf('  Eq. Block 1 / R1 module step %d of %d:',n,numel(core.path)));
        step2PrintOrientation(fid,sprintf('r_I_production_%02d',n),core.path(n),productionCtx);
        report(fid,sprintf('  Core flux multiplier per normalized turnover: %.12g',core.flux(n)));
    end
end
report(fid,'  X comprises all net-consumed precursors/energy substrates; W comprises other net-produced products.');
[atomStatus,imbalance]=productionAtomBalance(core.net,model);
report(fid,['  R1 module atom balance: ' char(atomStatus) ' | ' char(imbalance)]);
report(fid,'  S-controlled I import supplies seed/additional I. S-dependent stopping of R1 or a k*S*I law has not been established.');
lin=ctx.L(ctx.C.LinearID(c),:);quad=ctx.Q(ctx.C.QuadraticID(c),:);
step2Section(fid,'EQ. BLOCK 2: F-LINKED I REMOVAL');
report(fid,['  Proposed F signal: m_Fd = ' char(ctx.M.FastSignalID(b)) '. A consumed co-substrate is not F.']);
step2PrintOrientation(fid,'r_I_linear_removal',lin.OrientationIndex,ctx);
if lin.FTransportRow>0
    ft=ctx.T(lin.FTransportRow,:);step2PrintReaction(fid,'r_T_F',ft.ModelIndex,ft.Direction,ctx);
    report(fid,'  Imported co-substrate is consumed; it is not automatically the F controller.');
end
step2Section(fid,'EQ. BLOCK 3: I PAIR LOSS');
step2PrintOrientation(fid,'r_I_quadratic_removal',quad.OrientationIndex,ctx);
report(fid,sprintf('  I loss per primitive quadratic event: %g; exact formal 2I->I=%d; quadratic kinetics unverified.', ...
    quad.LossPerEvent,quad.ExactFormalStoichiometry));
step2Section(fid,'COUNTER AND ONE CONSUMING NEURON');
step2PrintOrientation(fid,'r_c',ctx.F.OrientationIndex(f),ctx);
report(fid,sprintf('  Counter consumes I; primitive I loss=%g; C gain=%g.',ctx.F.I_LossPerEvent(f),ctx.F.C_GainPerEvent(f)));
if ctx.X.LeakOrientationIndex(x)>0,step2PrintOrientation(fid,'r_cr',ctx.X.LeakOrientationIndex(x),ctx);
else,report(fid,'r_cr');report(fid,'  UNASSIGNED: formation-only counter has no selected independent C-removal reaction.');end
if cfg.neuronPlacement=="outside"
    if ~isempty(selectedTransport)
        report(fid,'  Counter-to-neuron interface:');
        report(fid,['  Intracellular adapter path: ' char(joinText(string(model.mets(selectedTransport.adapterMets)),' -> '))]);
        if isempty(selectedTransport.adapterO)
            report(fid,'  Zero-step adapter: m_exp = m_c.');
        else
            for n=1:numel(selectedTransport.adapterO)
                step2PrintOrientation(fid,sprintf('r_A_N1_%02d',n),selectedTransport.adapterO(n),ctx);
                report(fid,sprintf('  Adapter primitive flux multiplier: %.12g',selectedTransport.adapterFlux(n)));
            end
        end
        step2PrintReaction(fid,'r_T_N1',selectedTransport.rx,selectedTransport.direction,ctx);
        report(fid,sprintf('  Exported copies per model event=%g; reverse orientation allowed=%d.', ...
            selectedTransport.cargoCoefficient,selectedTransport.reverseAllowed));
    else
        report(fid,'r_T_N1');report(fid,'  UNASSIGNED: no isolated adapter plus bound-permitted cargo exporter was found within the configured limits.');
    end
    report(fid,'  The adapter consumes m_c and makes m_exp. Export transfers m_exp; the external neuron consumes m_Nin_1. These roles must not be merged.');
else
    report(fid,'  r_Nf_1 consumes m_c = m_Nin_1 and supplies a counter loss flux.');
end
if ~isempty(selectedNeuron)
    report(fid,['  Neuron consuming path: ' char(joinText(string(model.mets(selectedNeuron.mets)),' -> '))]);
    for n=1:numel(selectedNeuron.conversionO)
        role='r_Nf_1';
        if n==numel(selectedNeuron.conversionO),role='r_Nd_1';
        elseif n>1,role=sprintf('r_Nmid_1_%02d',n-1);end
        step2PrintOrientation(fid,role,selectedNeuron.conversionO(n),ctx);
    end
    step2PrintOrientation(fid,'r_Nr_1',selectedNeuron.removalO,ctx);
    if complete
        report(fid,['Assigned distinct controlled-role ECs (branches, counter, neuron, then Block0): ' char(joinText(motif.assignedECs,'; '))]);
    else
        report(fid,'  The neuron conversion chain, independent removal and prefix EC assignments passed, but this parent is not a complete circuit.');
    end
else
    for role=["r_Nf_1","r_Nd_1","r_Nr_1"]
        report(fid,char(role));report(fid,'  UNASSIGNED: no accepted single-neuron extension.');
    end
    if cfg.neuronPlacement=="outside"
        audit=ctx.currentExportAudit;
        audit=audit(audit.ParentID==d.id,:);
        for n=1:min(height(audit),cfg.step2PrintTop)
            report(fid,sprintf('  Export audit: %s (%+d) | input=%s | %s', ...
                audit.ReactionID(n),audit.Direction(n),audit.NeuronInputID(n),audit.Status(n)));
        end
        aa=ctx.currentAdapterAudit;aa=aa(aa.ParentID==d.id,:);
        for n=1:min(height(aa),cfg.step2PrintTop)
            report(fid,sprintf('  Adapter audit: depth=%d | %s (%+d) | %s -> %s | %s | I states=%s', ...
                aa.AdapterDepth(n),aa.ReactionID(n),aa.Direction(n),aa.FromID(n),aa.ToID(n),aa.Status(n),aa.IStateIDs(n)));
        end
        report(fid,'  See 14_Neuron_Export_Audit.csv and 14_Neuron_Entry_Audit.csv; external GSH sharing the selected upstream import pool remains disallowed.');
    else
        rows=ctx.O.consumers{ctx.F.C_Index(f)};
        for o=reshape(rows,1,[]),step2PrintOrientation(fid,'REJECTED_r_Nf_1_ENTRY',o,ctx);end
    end
end
if isempty(selectedNeuron),writeStep2ConsumingDiagnostics(fid,d.id,ctx,cfg);end
report(fid,'MODEL REQUIREMENTS: r_Nf_1 consumes m_Nin_1; m_out is m_Nd_1; r_Nr_1 removes m_out; Block0 ends at m_in and shares no other noncurrency state.');
report(fid,'STATUS: structural topology only; S/F control, buffered feeds, driving forces, fluxes and rate constants require validation.');
end

function writeStep2ConsumingDiagnostics(fid,parent,ctx,cfg)
if ~isfield(ctx,'currentChainInventory'),return;end
T=ctx.currentChainInventory;T=T(T.ParentID==parent,:);if isempty(T),return;end
step2Section(fid,'CONSUMING PATHS: DIAGNOSTICS, NOT ACCEPTED NEURONS');
report(fid,'These paths show model-recorded consumption. Feedback/overlap/control conflicts or missing output removal prevent their use as isolated neurons.');
for k=1:min(height(T),cfg.neuronChainPrintTop)
    report(fid,sprintf('  Path %d: %s | %s',k,T.MetabolitePath(k),T.Classification(k)));
    report(fid,sprintf('  Conversions=%d; net m_c/input change=%g; net I change=%g; I interaction=%d; upstream reaction reuse=%d; EC/gene assignment passed=%d.', ...
        T.ConversionCount(k),T.NetInputChange(k),T.NetIChange(k),T.IInteraction(k), ...
        T.ReusesUpstreamReaction(k),T.ECGeneAssignmentPassed(k)));
    report(fid,['  Shared upstream states: ' char(T.SharedUpstreamIDs(k))]);
    report(fid,['  Repeated nonadjacent states: ' char(T.NonadjacentSharedIDs(k))]);
    report(fid,['  Net path equation: ' char(T.NetEquation(k))]);
    diagnosticCtx=ctx;ids=split(T.MetabolitePath(k),' -> ');
    for n=2:numel(ids)
        i=find(string(ctx.model.mets)==ids(n),1);if isempty(i),continue;end
        diagnosticCtx.signalIndices(end+1)=i;
        diagnosticCtx.signalLabels(end+1)="m_diag_"+compose('%02d',n-1);
    end
    os=str2double(split(T.OrientationIndices(k),';'));
    for n=1:numel(os)
        step2PrintOrientation(fid,sprintf('DIAGNOSTIC_r_Npath_%02d',n),os(n),diagnosticCtx);
    end
end
report(fid,'Full path inventory and flags: 15_Consuming_Neuron_Chains.csv. These rows never contribute to the complete-circuit count.');
end

function step2PrintBlock0Path(fid,path,ctx,prefix)
report(fid,['  Connector trace: ' char(joinText(string(ctx.model.mets(path.mets)),' -> '))]);
report(fid,sprintf('  Reactions=%d; complete EC evidence=%d; states=%d; expansions=%d.', ...
    path.distance,path.ECComplete,path.states,path.expansions));
idx=path.mets(2:end-1);labels="m_B0_"+compose('%02d',1:numel(idx));
side=setdiff(path.internal,path.mets,'stable');idx=[reshape(idx,1,[]) reshape(side,1,[])];
labels=[reshape(labels,1,[]) "m_B0_side_"+compose('%02d',1:numel(side))];
ctx.signalIndices=[ctx.signalIndices idx];ctx.signalLabels=[ctx.signalLabels labels];
step2Signal(fid,'cAMP',path.mets(1),ctx.model);
for j=1:numel(idx),step2Signal(fid,char(labels(j)),idx(j),ctx.model);end
report(fid,['  Additional supplied feeds: ' char(joinText(string(ctx.model.mets(path.feed)),'; '))]);
report(fid,['  Waste products: ' char(joinText(string(ctx.model.mets(path.waste)),'; '))]);
report(fid,['  Buffered currency participants: ' char(joinText(string(ctx.model.mets(path.buffered)),'; '))]);
report(fid,['  Module net equation: ' char(step2NamedSide(max(-path.net,0),ctx)+" -> "+step2NamedSide(max(path.net,0),ctx))]);
report(fid,'  Stoichiometric weights below normalize one cAMP consumed; they are not kinetic rate constants.');
for n=1:numel(path.o)
    report(fid,sprintf('  %s_%d module flux multiplier: %.9g',prefix,n,path.flux(n)));
    step2PrintOrientation(fid,sprintf('%s_%d',prefix,n),path.o(n),ctx);
end
end

function step2PrintOrientation(fid,role,o,ctx)
step2PrintReaction(fid,role,ctx.O.rxn(o),ctx.O.direction(o),ctx);
end

function step2PrintReaction(fid,role,j,direction,ctx)
model=ctx.model;a=ctx.R(:,j);b=ctx.P(:,j);if direction<0,tmp=a;a=b;b=tmp;end
report(fid,char(role));report(fid,['  ID: ' model.rxns{j}]);report(fid,['  Name: ' model.rxnNames{j}]);
if isfield(model,'circuitMergedRxnIDs')&&numel(model.circuitMergedRxnIDs)==numel(model.rxns) ...
        &&contains(model.circuitMergedRxnIDs(j),";")
    report(fid,['  Cell-free merged copies: ' char(model.circuitMergedRxnIDs(j))]);
end
report(fid,sprintf('  Direction: %+d relative to model.S',direction));
ec=ctx.E.codes(j);if ec=="",ec="UNRESOLVED (no annotation inferred)";end
report(fid,['  Resolved EC: ' char(ec)]);
report(fid,['  EC annotation source: ' char(ctx.E.source(j))]);report(fid,['  EC annotation reference: ' char(ctx.E.reference(j))]);
report(fid,['  Equation: ' char(step2NamedSide(a,ctx)+" -> "+step2NamedSide(b,ctx))]);
[gpr,genes,~]=decodeGPR(string(model.grRules{j}),string(model.rules{j}),string(model.genes));
report(fid,['  GPR: ' char(gpr)]);report(fid,['  Resolved genes: ' char(genes)]);report(fid,'');
end

function s = step2NamedSide(v,ctx)
model=ctx.model;
idx=find(v>1e-8);pieces=strings(numel(idx),1);
for k=1:numel(idx)
    i=idx(k);name=step2MetName(model,i);coefficient=full(v(i));prefix="";
    if abs(coefficient-1)>1e-8,prefix=string(sprintf('%.12g ',coefficient));end
    if isfield(ctx,'signalIndices')&&isfield(ctx,'signalLabels')
        labels=string(ctx.signalLabels(ctx.signalIndices==i));
        labels=unique(labels(labels~=""),'stable');
        if ~isempty(labels),name=name+" ("+joinText(labels,' / ')+")";end
    end
    pieces(k)=prefix+name+" ["+string(model.mets{i})+"]";
end
s=joinText(pieces,' + ');if isempty(idx),s="0";end
end

function name = step2MetName(model,i)
name=string(model.mets{i});
if isfield(model,'metNames')&&numel(model.metNames)==numel(model.mets)
    candidate=string(model.metNames(i));if ~ismissing(candidate)&&candidate~="",name=candidate;end
end
end

function step2Signal(fid,role,i,model)
report(fid,sprintf('  %-9s: %s | %s',role,string(model.mets{i}),step2MetName(model,i)));
end

function step2Section(fid,title)
report(fid,sprintf('\n============================================================\n%s\n============================================================',title));
end

function [K,support,A,stats] = findCounterComponents(O,T,E,model,base,comp,currency,R,P,cfg)
% Pair import of exact C pool with a primitive single-C enzymatic loss.
rows=cell(0,22);support=cell(0,1);audit=cell(0,4);tol=1e-8;stop=false;
stats=struct('tested',0,'testLimitReached',false,'componentLimitReached',false);
if cfg.connectCounter
    for t=1:height(T)
        c=T.InsideIndex(t);tr=T.ModelIndex(t);
        if cfg.counterCargoID~=""&&string(model.mets{c})~=cfg.counterCargoID&&base(c)~=cfg.counterCargoID,continue;end
        tp=find(R(:,tr)+P(:,tr)>tol);
        [feed,feedStatus]=counterBoundaryFeed(T.OutsideIndex(t),model,R,P);
        for o=reshape(O.consumers{c},1,[])
            if stats.tested>=cfg.counterMaxComponentTests,stats.testLimitReached=true;stop=true;break;end
            stats.tested=stats.tested+1;lr=O.rxn(o);lp=find(O.R(:,o)+O.P(:,o)>tol);
            status="accepted_structural_counter_component";
            [scale,primitiveOK]=primitiveReactionScale(O.R(:,o),O.P(:,o));
            shared=intersect(tp,lp);shared(shared==c|currency(shared))=[];
            if ~T.GrossKnown(t)||~O.grossKnown(o),status="gross_stoichiometry_unknown";
            elseif lr==tr,status="transport_and_removal_reuse_reaction";
            elseif numel(tp)>cfg.counterMaxParticipants||numel(lp)>cfg.counterMaxParticipants
                status="participant_limit";
            elseif ~primitiveOK||abs(full(O.R(c,o))*scale-1)>tol||abs(full(O.P(c,o)))>tol
                status="not_primitive_single_C_loss";
            elseif E.codes(lr)==""&&cfg.counterRequireLeakEC,status="removal_EC_unresolved";
            elseif cfg.counterRequireIsolation&&~isempty(shared),status="transport_removal_noncurrency_overlap";
            elseif cfg.counterRequireBoundaryFeed&&~feed,status="outside_feed_not_allowed_by_boundary_bounds";
            end
            audit(end+1,:)={T.ReactionID(t),string(model.mets{c}),string(model.rxns{lr}),status}; %#ok<AGROW>
            if status~="accepted_structural_counter_component",continue;end
            rows(end+1,:)={size(rows,1)+1,comp(c),c,string(model.mets{c}),T.OutsideIndex(t), ...
                T.OutsideMetaboliteID(t),t,T.ReactionID(t),T.Direction(t),T.ImportCoefficient(t), ...
                o,lr,string(model.rxns{lr}),O.direction(o),scale, ...
                full(O.R(c,o)-O.P(c,o))*scale,E.codes(lr),E.reference(lr),T.Equation(t), ...
                equation(O.R(:,o),O.P(:,o),model),feed,feedStatus}; %#ok<AGROW>
            support{end+1,1}=union(tp,lp); %#ok<AGROW>
            if size(rows,1)>=cfg.counterMaxComponents,stats.componentLimitReached=true;stop=true;break;end
        end
        if stop,break;end
    end
end
K=cellTable(rows,{'ComponentID','InsideCompartment','C_Index','C_ID','C_OutsideIndex','C_OutsideID', ...
    'TransportRow','TransportReactionID','TransportDirection','ImportCoefficient', ...
    'LeakOrientationIndex','LeakReactionIndex','LeakReactionID','LeakDirection','LeakPrimitiveScale', ...
    'LeakLossPerEvent','LeakECs','LeakECReference','TransportEquation','LeakEquation', ...
    'OutsideFeedAllowed','BoundaryFeedStatus'}, ...
    {'double','string','double','string','double','string','double','string','double','double', ...
    'double','double','string','double','double','double','string','string','string','string','logical','string'});
A=cellTable(audit,{'TransportReactionID','C_ID','LeakReactionID','Status'},repmat({'string'},1,4));
end

function [allowed,status] = counterBoundaryFeed(e,model,R,P)
% Local boundary check only. No FBA, no changes to exchanges or the medium.
allowed=false;found=false;
for j=reshape(find(R(e,:)+P(e,:)>1e-10),1,[])
    if nnz(R(:,j)+P(:,j)>1e-10)~=1,continue;end
    found=true;net=full(P(e,j)-R(e,j));
    if (net>1e-10&&model.ub(j)>0)||(net < -1e-10&&model.lb(j)<0),allowed=true;end
end
status="no_single_species_boundary_supply_found";
if found,status="single_species_boundary_supply_closed";end
if allowed,status="single_species_boundary_supply_allowed_not_FBA_validated";end
end

function names = counterGateColumns()
names={'I_ID','TransportReactionID','Direction','CounterMetaboliteID','SignalSide', ...
    'Effect','Mechanism','I_NotConsumed','TransporterProteinID','TargetMappingReference','Organism','Reference'};
end

function [G,A] = readCounterGateEvidence(path,K,C,cfg)
names=counterGateColumns();G=cellTable(cell(0,numel(names)),names,repmat({'string'},1,numel(names)));
A=G;A.Accepted=false(0,1);A.Status=strings(0,1);if isempty(path),return;end
G=readTextTable(path);
assert(all(ismember(names,G.Properties.VariableNames)),'Anderson:CounterGateColumns', ...
    'Counter gate evidence needs: %s',strjoin(names,', '));
G=G(:,names);A=G;A.Accepted=false(height(G),1);A.Status=strings(height(G),1);
for i=1:height(G)
    row=G(i,:);status="supplied_gate_reference_linked_not_independently_verified";
    k=find(K.TransportReactionID==row.TransportReactionID & K.TransportDirection==str2double(row.Direction) ...
        & K.C_ID==row.CounterMetaboliteID);
    if isempty(k),status="no_exact_transport_direction_C_match";
    elseif ~any(C.I_ID==row.I_ID & ismember(C.InsideCompartment,K.InsideCompartment(k)))
        status="I_pool_or_signal_side_compartment_mismatch";
    elseif lower(row.SignalSide)~="inside",status="intracellular_I_requires_inside_access";
    elseif lower(row.Effect)~="opening"||~ismember(lower(row.Mechanism), ...
            ["reversible_binding","nonconsuming_protein_interaction"])
        status="nonconsuming_opening_mechanism_not_declared";
    elseif ~asBool(row.I_NotConsumed),status="I_not_declared_unconsumed";
    elseif row.TransporterProteinID==""||row.TargetMappingReference==""||row.Reference==""
        status="missing_protein_mapping_or_mechanism_reference";
    elseif ~strcmpi(row.Organism,cfg.organismName),status="organism_mismatch";
    end
    A.Status(i)=status;A.Accepted(i)=status=="supplied_gate_reference_linked_not_independently_verified";
end
G=G(A.Accepted,:);
end

function [rx,parts] = counterCancellerSupport(c,C,D,cores,O,T,L,Q,R,P)
prod=C.ProductionCandidateID(c);core=cores(D.CoreIndex(prod));
iRow=D.TransportRow(prod);l=L(C.LinearID(c),:);q=Q(C.QuadraticID(c),:);
rx=[O.rxn(core.path),l.ReactionIndex,q.ReactionIndex];
parts=find(any(O.R(:,core.path)+O.P(:,core.path)>1e-8,2));
if iRow>0
    st=T.ModelIndex(iRow);rx(end+1)=st;
    parts=union(parts,find(R(:,st)+P(:,st)>1e-8));
end
parts=union(parts,find(O.R(:,l.OrientationIndex)+O.P(:,l.OrientationIndex)>1e-8));
parts=union(parts,find(O.R(:,q.OrientationIndex)+O.P(:,q.OrientationIndex)>1e-8));
if l.FTransportRow>0
    ft=T.ModelIndex(l.FTransportRow);rx(end+1)=ft;
    parts=union(parts,find(R(:,ft)+P(:,ft)>1e-8));
end
end

function [X,stats] = connectCounters(K,support,G,C,D,cores,O,T,L,Q,model,base,currency,R,P,cfg)
rows=cell(0,15);stop=false;
stats=struct('tested',0,'signalInChemistryRejected',0,'reactionReuseRejected',0, ...
    'crossTalkRejected',0,'missingGateEvidenceRejected',0,'testLimitReached',false,'connectionLimitReached',false);
for c=1:height(C)
    [cr,cp]=counterCancellerSupport(c,C,D,cores,O,T,L,Q,R,P);
    for k=reshape(find(K.InsideCompartment==C.InsideCompartment(c)),1,[])
        if stats.tested>=cfg.counterMaxConnectionTests,stats.testLimitReached=true;stop=true;break;end
        stats.tested=stats.tested+1;kp=support{k};
        % Exclude I from ALL gross chemistry, including outside I or a
        % regenerated I catalyst. I acts only in the proposed gate mechanism.
        if any(base(kp)==base(C.I_Index(c)))
            stats.signalInChemistryRejected=stats.signalInChemistryRejected+1;continue;
        end
        kr=[T.ModelIndex(K.TransportRow(k)),K.LeakReactionIndex(k)];
        if any(ismember(kr,cr)),stats.reactionReuseRejected=stats.reactionReuseRejected+1;continue;end
        shared=intersect(kp,cp);shared(currency(shared))=[];isolated=isempty(shared);
        if cfg.counterRequireIsolation&&~isolated,stats.crossTalkRejected=stats.crossTalkRejected+1;continue;end
        matches=find(G.I_ID==C.I_ID(c)&G.TransportReactionID==K.TransportReactionID(k) ...
            &str2double(G.Direction)==K.TransportDirection(k)&G.CounterMetaboliteID==K.C_ID(k));
        if cfg.counterRequireGateEvidence&&isempty(matches)
            stats.missingGateEvidenceRejected=stats.missingGateEvidenceRejected+1;continue;
        end
        gateStatus="I_dependent_opening_unverified";
        if ~isempty(matches),gateStatus="supplied_reference_linked_not_independently_verified";end
        rows(end+1,:)={size(rows,1)+1,C.ConnectedID(c),K.ComponentID(k),C.I_Index(c),C.I_ID(c), ...
            K.C_ID(k),K.InsideCompartment(k),true,isolated,joinText(string(model.mets(shared)),';'), ...
            numel(matches),joinText(G.Reference(matches),'; '),gateStatus,K.OutsideFeedAllowed(k), ...
            "provisional_gating_transport_removal_kinetics_and_upstream_interfaces"}; %#ok<AGROW>
        if size(rows,1)>=cfg.counterMaxConnections,stats.connectionLimitReached=true;stop=true;break;end
    end
    if stop,break;end
end
X=cellTable(rows,{'CounterConnectionID','ConnectedID','ComponentID','I_Index','I_ID','C_ID', ...
    'InsideCompartment','I_AbsentFromCounterChemistry','NoncurrencyIsolationPassed', ...
    'SharedNoncurrency','GateEvidenceRows','GateReferences','GateStatus','OutsideFeedAllowed','Status'}, ...
    {'double','double','double','double','string','string','string','logical','logical', ...
    'string','double','string','string','logical','string'});
end

function [Y,stats] = wireCountersToBranches(X,K,support,W,M,rails,model,currency,cfg)
rows=cell(0,15);stop=false;
stats=struct('tested',0,'reactionReuseRejected',0,'crossTalkRejected',0,'ECOverlapRejected',0, ...
    'testLimitReached',false,'wiringLimitReached',false);
for x=1:height(X)
    k=X.ComponentID(x);kr=[modelIndexForCounter(K,k,model),K.LeakReactionIndex(k)];kp=support{k};
    for w=reshape(find(W.ConnectedID==X.ConnectedID(x)),1,[])
        if stats.tested>=cfg.counterMaxWiringTests,stats.testLimitReached=true;stop=true;break;end
        stats.tested=stats.tested+1;b=find(M.ModuleID==W.ModuleID(w),1);
        f=rails(M.FastRailIndex(b));s=rails(M.SlowRailIndex(b));
        if any(ismember(kr,[f.rx s.rx])),stats.reactionReuseRejected=stats.reactionReuseRejected+1;continue;end
        shared=intersect(kp,union(f.parts,s.parts));shared(currency(shared))=[];
        isolated=isempty(shared);
        if cfg.counterRequireIsolation&&~isolated,stats.crossTalkRejected=stats.crossTalkRejected+1;continue;end
        codes=completeECTokens(K.LeakECs(k));assigned=completeECTokens(M.AssignedECs(b));
        available=codes(~ismember(codes,assigned));distinct=~isempty(available);
        if cfg.counterDistinctECFromBranches&&~distinct,stats.ECOverlapRejected=stats.ECOverlapRejected+1;continue;end
        selected="";if distinct,selected=available(1);elseif ~isempty(codes),selected=codes(1);end
        rows(end+1,:)={size(rows,1)+1,W.WiringID(w),W.ModuleID(w),X.ConnectedID(x), ...
            X.CounterConnectionID(x),K.ComponentID(k),X.I_ID(x),X.C_ID(x), ...
            selected,distinct,isolated,joinText(string(model.mets(shared)),';'), ...
            W.BothControlInterfacesMatched(w),X.GateEvidenceRows(x)>0, ...
            "provisional_all_gate_mechanisms_rates_and_enzyme_independence_require_validation"}; %#ok<AGROW>
        if size(rows,1)>=cfg.counterMaxWiring,stats.wiringLimitReached=true;stop=true;break;end
    end
    if stop,break;end
end
Y=cellTable(rows,{'CounterWiringID','BranchWiringID','ModuleID','ConnectedID','CounterConnectionID', ...
    'ComponentID','I_ID','C_ID','AssignedCounterLeakEC','DistinctFromAssignedBranchECs', ...
    'NoncurrencyIsolationPassed','SharedNoncurrency','UpstreamSFInterfacesMatched', ...
    'CounterGateReferenceLinked','Status'}, ...
    {'double','double','double','double','double','double','string','string','string','logical', ...
    'logical','string','logical','logical','string'});
end

function j = modelIndexForCounter(K,k,model)
j=find(string(model.rxns)==K.TransportReactionID(k),1);
assert(~isempty(j),'Anderson:CounterTransport','Counter transporter is absent from model.');
end

function writeCounterTemplates(K,X,model,cfg)
% Leave all kinetic constants unmeasured. ECs are not kinetic parameters.
path=fullfile(cfg.outputDir,'Counter_EC_Evidence_TEMPLATE.tsv');
if ~exist(path,'file')
    blank=table(strings(0,1),strings(0,1),strings(0,1), ...
        'VariableNames',{'SourceReactionID','EC','Reference'});
    writetable(blank,path,'FileType','text','Delimiter','\t');
end
path=fullfile(cfg.outputDir,'Counter_Gate_Evidence_TEMPLATE.tsv');
if ~exist(path,'file')
    names=counterGateColumns();blank=cellTable(cell(0,numel(names)),names,repmat({'string'},1,numel(names)));
    writetable(blank,path,'FileType','text','Delimiter','\t');
end
rows=cell(0,16);
for x=1:height(X)
    k=X.ComponentID(x);j=K.LeakReactionIndex(k);
    rows(end+1,:)={X.CounterConnectionID(x),X.I_ID(x),X.C_ID(x),K.TransportReactionID(k), ...
        K.LeakReactionID(k),K.LeakECs(k),string(model.grRules{j}), ...
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,"unmeasured_effective_leak_requires_low_C_regime"}; %#ok<AGROW>
end
KT=cellTable(rows,{'CounterConnectionID','I_ID','C_ID','TransportReactionID','LeakReactionID','LeakECs', ...
    'LeakGPR','GateKd_concentration','GateCapacity_concentration','TransportVmax_concentration_per_s', ...
    'EffectiveAlpha_per_s','EffectiveKcr_per_s','LeakKm_concentration','LeakKcat_per_s', ...
    'LeakEnzyme_concentration','KineticsStatus'}, ...
    {'double','string','string','string','string','string','string','double','double','double', ...
    'double','double','double','double','double','string'});
writetable(KT,fullfile(cfg.outputDir,'10_Counter_Kinetics_Template.csv'));
end

function writeCounterReport(path,K,X,Y,stats,cstats,wstats,cfg)
fid=fopen(path,'w');assert(fid>0,'Anderson:Report','Cannot write counter report.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
report(fid,'ANDERSON CANCELLER -> I-CONTROLLED MEMBRANE -> COUNTER');
report(fid,sprintf('Counter enabled=%d; components=%d; canceller/counter=%d; branch/canceller/counter=%d.', ...
    cfg.connectCounter,height(K),height(X),height(Y)));
report(fid,sprintf('Connections with supplied I-opening references=%d (references not independently verified).',nnz(X.GateEvidenceRows>0)));
report(fid,sprintf('Components with outside supply allowed by single-species boundary bounds=%d of %d.',nnz(K.OutsideFeedAllowed),height(K)));
report(fid,'Closed supply is retained as an explicit external-reservoir design assumption by default. No bounds were changed.');
report(fid,'This is a partial EC-based structural screen, not a complete yeast transporter/enzyme survey.');
report(fid,'The separate GALK seed uses an explicit BiGG source ID AND the complete chemical signature; it does not change branch annotations.');
report(fid,'Counter chemistry neither consumes, regenerates nor transports I: all gross participants with the same canonical chemical ID are excluded.');
report(fid,'Proposed regulatory coupling: G + I <-> GI; GI + C_out -> GI + C_in (other transport co-substrates retained in CSV).');
report(fid,'G is a gate state, not an inferred metabolic enzyme. Binding is hypothetical until mechanism-specific evidence is supplied.');
report(fid,'Binding alone conserves I_total = I_free + GI, but transiently sequesters free I. Model both pools if gate occupancy is appreciable.');
report(fid,'Do not add I as an extra substrate/product in an ordinary COBRA S column to encode activation; that cancels to zero and supplies no kinetic control.');
report(fid,'Counter mass balance: dC/dt = ImportCoefficient*v_transport - v_leak/LeakPrimitiveScale.');
report(fid,'Proposed effective gate occupancy: GI = G_total*I_free/(Kd+I_free), under rapid reversible binding.');
report(fid,'Transport must also depend on cargo availability, the membrane driving force and any coupled species.');
report(fid,'Only with I_free << Kd, stable driving/feed conditions and unsaturated removal: dC/dt approximately alpha*I_free - kcr*C.');
report(fid,'A metabolic removal coefficient of one does not establish first-order removal. For Michaelis-Menten loss, kcr approximately kcat*E/Km only when C << Km.');
report(fid,'This is a leaky integral of I pulse area. It approximates pulse count only for standardized pulse areas and slow leakage over the counting window.');
report(fid,'EC separation enables separate parameter roles in the design; it does not prove distinct proteins, independent control or absent cellular bypasses.');
report(fid,'No I-gated transporter, whole-cell feasibility, membrane implementation or kinetic simulation has been established by these rows.');
report(fid,'Upstream S/F interface status is retained in the final wiring CSV; adding a counter does not resolve it.');
report(fid,sprintf('Component tests=%d; test cap=%d; component cap=%d.',stats.tested,stats.testLimitReached,stats.componentLimitReached));
report(fid,sprintf('Connection tests=%d; I-in-chemistry rejected=%d; reaction reuse=%d; overlap=%d; missing evidence=%d; test cap=%d; connection cap=%d.', ...
    cstats.tested,cstats.signalInChemistryRejected,cstats.reactionReuseRejected,cstats.crossTalkRejected, ...
    cstats.missingGateEvidenceRejected,cstats.testLimitReached,cstats.connectionLimitReached));
report(fid,sprintf('Branch joins tested=%d; reaction reuse=%d; overlap=%d; EC overlap=%d; test cap=%d; wiring cap=%d.', ...
    wstats.tested,wstats.reactionReuseRejected,wstats.crossTalkRejected,wstats.ECOverlapRejected,wstats.testLimitReached,wstats.wiringLimitReached));
for k=1:min(height(K),cfg.printTop)
    report(fid,sprintf('\nCOMPONENT %d | C=%s | %s',K.ComponentID(k),K.C_ID(k),K.InsideCompartment(k)));
    report(fid,['  Import: ' char(K.TransportReactionID(k)) ' | ' char(K.TransportEquation(k))]);
    report(fid,['  Removal: ' char(K.LeakReactionID(k)) ' | EC ' char(K.LeakECs(k)) ' | ' char(K.LeakEquation(k))]);
    report(fid,['  Boundary: ' char(K.BoundaryFeedStatus(k))]);
    report(fid,'  I-dependent opening and all kinetic constants: unverified/unmeasured.');
end
end

function n = heightFromRows(rows)
n=size(rows,1);
end

function [MAP,E,H] = discoverRegulators(T,root,cfg,fid)
% Direct positive causal edges only, retaining the full GPR elsewhere.
MAP=cellTable(cell(0,3),{'gene_id','uniprot_accession','mapping_status'}, ...
    {'string','string','string'});
E=cellTable(cell(0,5),{'source','target','effect','database','reference'},repmat({'string'},1,5));
H=cellTable(cell(0,8),{'TransportReactionID','ImportedMetaboliteID','TransporterGene', ...
    'TransporterAccession','S_ProteinHypothesis','EvidenceClass','Reference','Source'},repmat({'string'},1,8));
if isempty(T),return;end
ids=unique(split(joinText(T.Genes,';'),';')); ids=ids(ids~=""&~ismissing(ids));
if isempty(ids),return;end % no transporter genes to map in an unannotated model
mapPath=findNamed(root,cfg.proteinMapFile);
if ~isempty(mapPath)
    L=readTextTable(mapPath);
    assert(all(ismember({'gene_id','uniprot_accession'},L.Properties.VariableNames)), ...
        'Anderson:ProteinMap','Protein mapping needs gene_id and uniprot_accession.');
    for k=1:height(L)
        if any(ids==L.gene_id(k)) && L.uniprot_accession(k)~=""
            MAP(end+1,:)={L.gene_id(k),L.uniprot_accession(k),"supplied_exact_gene_mapping"}; %#ok<AGROW>
        end
    end
end
if cfg.onlineProteinDiscovery
    for k=1:numel(ids)
        if any(MAP.gene_id==ids(k)),continue;end
        acc="";status="mapping_failed";
        try
            gene=webread(['https://mygene.info/v3/gene/' urlencodeLocal(ids(k)) ...
                '?fields=symbol,taxid'],weboptions('Timeout',cfg.networkTimeout));
            assert(isfield(gene,'taxid') && double(gene.taxid)==cfg.organismTaxID,'Taxon mismatch');
            query=['(gene_exact:' char(gene.symbol) ') AND (organism_id:' ...
                num2str(cfg.organismTaxID) ') AND (reviewed:true)'];
            data=webread(['https://rest.uniprot.org/uniprotkb/search?query=' urlencodeLocal(query) ...
                '&format=json&size=25'],weboptions('Timeout',cfg.networkTimeout,'ContentType','json'));
            accessions=strings(0,1);
            if isfield(data,'results')
                for h=1:numel(data.results)
                    hit=data.results(h);exact=false;
                    if isfield(hit,'genes')
                        for g=1:numel(hit.genes)
                            gn=hit.genes(g);
                            if isfield(gn,'geneName') && strcmpi(gn.geneName.value,gene.symbol),exact=true;end
                        end
                    end
                    if exact,accessions(end+1,1)=string(hit.primaryAccession);end %#ok<AGROW>
                end
            end
            acc=joinText(unique(accessions),';');status="reviewed_exact_gene";
            if acc=="",status="no_reviewed_exact_gene";end
        catch ME
            report(fid,['Protein mapping ' char(ids(k)) ': ' ME.message]);
        end
        MAP(end+1,:)={ids(k),acc,status}; %#ok<AGROW>
    end
end
edgePath=findNamed(root,cfg.proteinEdgesFile);
if ~isempty(edgePath)
    L=readTextTable(edgePath);
    if all(ismember({'source','target','effect'},L.Properties.VariableNames))
        for k=1:height(L)
            ref="";db="supplied";
            if ismember('reference',L.Properties.VariableNames),ref=L.reference(k);end
            if ismember('database',L.Properties.VariableNames),db=L.database(k);end
            E(end+1,:)={L.source(k),L.target(k),lower(L.effect(k)),db,ref}; %#ok<AGROW>
        end
    else
        E=parseOmniPath(L);
    end
elseif cfg.onlineProteinDiscovery
    try
        url=['https://omnipathdb.org/interactions?datasets=omnipath&organisms=' ...
            num2str(cfg.organismTaxID) '&fields=references,sources&format=tsv'];
        txt=webread(url,weboptions('Timeout',cfg.networkTimeout,'ContentType','text'));
        temp=[tempname '.tsv'];f=fopen(temp,'w');assert(f>0);fwrite(f,txt,'char');fclose(f);
        cleanup=onCleanup(@()delete(temp)); %#ok<NASGU>
        E=parseOmniPath(readTextTable(temp));
    catch ME
        report(fid,['Protein edge retrieval unavailable: ' ME.message]);
    end
end
E=unique(E,'rows','stable');
for t=1:height(T)
    genes=split(T.Genes(t),';');
    for g=reshape(find(ismember(MAP.gene_id,genes)),1,[])
        acc=split(MAP.uniprot_accession(g),';');acc=acc(acc~="");
        for e=reshape(find(ismember(E.target,acc)&E.effect=="activation"&E.reference~=""),1,[])
            H(end+1,:)={T.ReactionID(t),T.InsideMetaboliteID(t),MAP.gene_id(g), ...
                E.target(e),E.source(e),"positive_protein_edge_NOT_opening_proof", ...
                E.reference(e),E.database(e)}; %#ok<AGROW>
        end
    end
end
H=unique(H,'rows','stable');
report(fid,sprintf('Protein discovery: %d mappings; %d edges; %d hypotheses (not gating evidence).', ...
    height(MAP),height(E),height(H)));
end

function E = parseOmniPath(L)
E=cellTable(cell(0,5),{'source','target','effect','database','reference'},repmat({'string'},1,5));
assert(all(ismember({'source','target','is_stimulation','is_inhibition'},L.Properties.VariableNames)), ...
    'Anderson:Edges','Unrecognized causal-edge table.');
for k=1:height(L)
    st=asBool(L.is_stimulation(k));inh=asBool(L.is_inhibition(k));
    if ~st&&~inh,continue;end
    eff="ambiguous";if st&&~inh,eff="activation";elseif inh&&~st,eff="inhibition";end
    ref="";if ismember('references',L.Properties.VariableNames),ref=L.references(k);end
    E(end+1,:)={L.source(k),L.target(k),eff,"OmniPath",ref}; %#ok<AGROW>
end
end

function T = readTextTable(path)
[~,~,ext]=fileparts(path); delim=sprintf('\t');if strcmpi(ext,'.csv'),delim=',';end
opts=detectImportOptions(path,'FileType','text','Delimiter',delim,'VariableNamingRule','preserve');
opts=setvartype(opts,opts.VariableNames,'string');T=readtable(path,opts);
for k=1:width(T)
    nm=T.Properties.VariableNames{k};v=strtrim(string(T.(nm)));v(ismissing(v))="";T.(nm)=v;
end
end

function T = cellTable(rows,names,types)
T=table('Size',[0 numel(names)],'VariableTypes',types,'VariableNames',names);
for k=1:size(rows,1),T(end+1,:)=rows(k,:);end %#ok<AGROW>
end

function s = equation(a,b,model)
s=side(a,model)+" -> "+side(b,model);
end

function s = side(v,model)
ids=find(v>1e-8);pieces=strings(numel(ids),1);
for k=1:numel(ids)
    i=ids(k);
    % Indexing a sparse matrix returns a sparse scalar. sprintf requires a
    % full numeric input. Convert only this scalar, retaining sparse matrices.
    coefficient=full(v(i));
    prefix="";
    if abs(coefficient-1)>1e-8
        prefix=string(sprintf('%.9g ',coefficient));
    end
    pieces(k)=prefix+string(model.mets{i});
end
s=joinText(pieces,' + ');if isempty(ids),s="0";end
end

function selfTestAnderson()
% Meaningful chemistry regressions; no real model/network/toolbox required.
cfg=defaultOptions();cfg.requireIImport=true;
selfTestBlock0Alternatives();
assert(cfg.productionMode=="multi_reaction"&&cfg.maxCoreReactions==10, ...
    'Default search must permit the original multi-step production modules.');
cfg.maxCoreReactions=3;cfg.maxStatesPerI=Inf;
cfg.insideCompartment="c";cfg.requireCoreEC=true;cfg.requireTransportGPR=true;
m=struct('mets',{{'i[e]';'i[c]';'x[c]';'y[c]';'w[c]';'q[c]'}}, ...
    'rxns',{{'T';'A';'B'}},'S',sparse([-1 0 0;1 -1 2;0 -1 0;0 1 -1;0 0 1;0 0 0]), ...
    'lb',[0;0;0],'ub',[100;100;100],'genes',{{'gT';'gA';'gB'}}, ...
    'grRules',{{'gT';'gA';'gB'}},'rxnECNumbers',{{'';'1.1.1.1';'2.2.2.2'}});
m=normalizeModel(m);m.circuitResolvedEC=string(m.rxnECNumbers);
% Sparse reporting regression: non-unit coefficients reach sprintf as full
% scalars. Check integer, fractional, unit and empty equation sides.
assert(side(sparse([0;2;0;0;0;0]),m)=="2 i[c]", ...
    'Sparse integer coefficient formatting failed.');
assert(side(sparse([0;0.5;0;0;0;0]),m)=="0.5 i[c]", ...
    'Sparse fractional coefficient formatting failed.');
assert(side(sparse([0;1;0;0;0;0]),m)=="i[c]", ...
    'Unit coefficient formatting failed.');
assert(side(sparse(6,1),m)=="0",'Empty equation-side formatting failed.');
[base,comp]=splitCompartments(string(m.mets));cur=false(6,1);
R=max(-m.S,0);P=max(m.S,0);
[O,T]=buildPools(m,base,comp,cur,R,P,true(3,1),cfg);
assert(height(T)==1 && T.InsideIndex==2,'Forward importer missing.');
[found,stats]=findCores(2,O,cur,cfg);
assert(numel(found)==1 && ~stats.capped && stats.unsolved==0,'Two-step core missing.');
assert(norm(found.net-[0;1;-1;0;1;0],Inf)<1e-9,'Wrong core net equation.');
moduleGates=cellTable(cell(0,11),gateColumns(),repmat({'string'},1,11));
[moduleCandidates,moduleDetails]=assembleCandidates(found,O,T,moduleGates,m,cfg);
modulePairs=describeImportProductionPairs(moduleCandidates,moduleDetails,found,O,T,m);
assert(height(modulePairs)==1&&~modulePairs.SingleProductionReaction&&modulePairs.ProductionReactionCount==2, ...
    'Default import/production matching discarded a multi-step autocatalytic realization.');
assert(modulePairs.R1_EffectiveModuleEquation=="i[c] + x[c] -> 2 i[c] + w[c]" ...
    &&modulePairs.I_ConsumedPerModule==1&&modulePairs.I_ProducedPerModule==2, ...
    'Effective R1 did not cancel its internal intermediate or retain actual I consumption/production.');
% Net-zero recycling must not be mistaken for duplication.
bad=m;bad.S(2,3)=1;Rb=max(-bad.S,0);Pb=max(bad.S,0);
[O2,~]=buildPools(bad,base,comp,cur,Rb,Pb,true(3,1),cfg);
f=findCores(2,O2,cur,cfg);assert(isempty(f),'Net-zero loop accepted.');
% Direct hidden catalyst requires gross evidence, not an invented I term.
direct=m;direct.S(:,2)=[0;1;-1;0;1;0];direct.ub(3)=0;
Rd=max(-direct.S,0);Pd=max(direct.S,0);Rd(2,2)=1;Pd(2,2)=2;
[Od,~]=buildPools(direct,base,comp,cur,Rd,Pd,true(3,1),cfg);
f=findCores(2,Od,cur,cfg);assert(numel(f)==1,'Gross one-step core missing.');
[On,~]=buildPools(direct,base,comp,cur,max(-direct.S,0),max(direct.S,0),false(3,1),cfg);
f=findCores(2,On,cur,cfg);assert(isempty(f),'Net producer falsely made autocatalytic.');
% An unprovided noncurrency catalyst prevents a core from operating.
Rq=R;Pq=P;Rq(6,3)=1;Pq(6,3)=1;
[Oq,~]=buildPools(m,base,comp,cur,Rq,Pq,true(3,1),cfg);
f=findCores(2,Oq,cur,cfg);assert(isempty(f),'Unavailable co-substrate ignored.');
% Conversely a declared buffered currency cofactor may be supplied.
curQ=cur;curQ(6)=true;
[Oq,~]=buildPools(m,base,comp,curQ,Rq,Pq,true(3,1),cfg);
f=findCores(2,Oq,curQ,cfg);assert(numel(f)==1,'Buffered cofactor policy broken.');
% A feed-driven reaction makes the selected set seed-independent.
a=sparse([0 0;1 0;1 1;0 1;0 0;0 0]);b=sparse([0 0;0 3;0 0;1 0;0 0;0 0]);
available=false(6,1);available(3)=true;
reach=substrateClosure(available,a,b);assert(~reach(2),'AND-substrate check broken.');
% Bound and reverse-policy regression.
rev=m;rev.S(:,1)=-rev.S(:,1);rev.lb(1)=-100;rev.ub(1)=0;
cfg.allowReverse=false;
[~,Tr]=buildPools(rev,base,comp,cur,max(-rev.S,0),max(rev.S,0),true(3,1),cfg);
assert(isempty(Tr),'Reverse importer bypassed strict +1 policy.');
cfg.allowReverse=true;
[~,Tr]=buildPools(rev,base,comp,cur,max(-rev.S,0),max(rev.S,0),true(3,1),cfg);
assert(height(Tr)==1 && Tr.Direction==-1,'Explicit reverse option failed.');
cfg.allowReverse=false;
blocked=m;blocked.ub(1)=0;
[~,Tb]=buildPools(blocked,base,comp,cur,R,P,true(3,1),cfg);
assert(isempty(Tb),'Blocked transporter accepted.');
% Explicit cycle-depth and state caps must be visible.
cfg.maxCoreReactions=1;f=findCores(2,O,cur,cfg);assert(isempty(f),'Depth cap ignored.');
cfg.maxCoreReactions=3;cfg.maxStatesPerI=1;
[~,st]=findCores(2,O,cur,cfg);assert(st.capped,'State cap not reported.');
% Gene AND rule and EC fallback are preserved from original helpers.
[g,genes,~]=decodeGPR("","x(1) & x(2)",["ga";"gb"]);
assert(contains(g,'&') && genes=="ga;gb",'GPR AND rule lost.');
assert(isempty(completeECTokens('1.2.3.-')),'Partial EC admitted.');
% Candidate module accounting includes imported I, so total inside gain=2.
emptyGate=cellTable(cell(0,11),gateColumns(),repmat({'string'},1,11));
[CT,~]=assembleCandidates(found,O,T,emptyGate,m,setfield(defaultOptions(),'requireIImport',true));
assert(height(CT)==1 && CT.GateEvidenceRows==0 && CT.S_IDs=="",'S was invented.');
assert(contains(CT.ModuleNetEquation,'2 i[c]'),'Import/core module accounting wrong.');
% Gate evidence is matched to reaction, direction, cargo AND GPR gene.
gateFile=[tempname '.tsv'];cleanup=onCleanup(@()delete(gateFile)); %#ok<NASGU>
gate=cellTable(cell(0,11),gateColumns(),repmat({'string'},1,11));
gate(end+1,:)={"T","1","i[c]","gT","S_fixture","protein","opening", ...
    "protein_gating","true","Homo sapiens","SYNTHETIC_TEST_ONLY"};
gate(end+1,:)=gate(1,:);gate.Effect(2)="activation";
gate(end+1,:)=gate(1,:);gate.S_NotConsumed(3)="false";
gate(end+1,:)=gate(1,:);gate.TargetGene(4)="wrong_gene";
writetable(gate,gateFile,'FileType','text','Delimiter','\t');
gateCfg=defaultOptions();gateCfg.requireIImport=true;gateCfg.organismName="Homo sapiens";
[accepted,audit]=readGateEvidence(gateFile,T,gateCfg);
assert(height(accepted)==1 && nnz(audit.Accepted)==1, ...
    'Generic activation, S consumption or wrong-gene evidence admitted as gating.');
% Missing annotations affect evidence status, not a relaxed structural screen.
unannotated=m;unannotated.grRules=repmat({''},3,1);unannotated.rules=repmat({''},3,1);
unannotated.circuitResolvedEC=strings(3,1);
strictCfg=defaultOptions();strictCfg.requireIImport=true;strictCfg.insideCompartment="c";strictCfg.productionMode="multi_reaction";strictCfg.maxCoreReactions=3;
strictCfg.requireCoreEC=true;strictCfg.requireTransportGPR=true;
[Os,Ts]=buildPools(unannotated,base,comp,cur,R,P,true(3,1),strictCfg);
assert(isempty(Os.rxn)&&isempty(Ts),'Strict annotation filters were bypassed.');
relaxedCfg=strictCfg;relaxedCfg.requireCoreEC=false;relaxedCfg.requireTransportGPR=false;
[Ou,Tu]=buildPools(unannotated,base,comp,cur,R,P,true(3,1),relaxedCfg);
assert(numel(Ou.rxn)==2&&height(Tu)==1,'Unannotated structural topology was excluded.');
unannotatedCores=findCores(2,Ou,cur,relaxedCfg);
[unannotatedCandidates,~]=assembleCandidates(unannotatedCores,Ou,Tu,emptyGate,unannotated,relaxedCfg);
assert(height(unannotatedCandidates)==1 && ~unannotatedCandidates.AllCoreECResolved ...
    && ~unannotatedCandidates.TransportHasGPR ...
    && unannotatedCandidates.AnnotationStatus=="structural_only_missing_EC_or_GPR", ...
    'Missing enzyme/protein evidence was not exposed on the structural candidate.');
% Community IDs and metComps must recover both imports without merging pools.
community=struct('mets',{{'M_2hb_e';'M_AUXOSC2hb_u';'M_h_e';'M_AUXOSCh_u'; ...
    'M_PROTSC2hb_x';'M_PROTSCh_x'}},'rxns',{{'T_u';'T_x';'cross_block'}}, ...
    'S',sparse([-1 -1 0;1 0 -1;-1 -1 0;1 0 0;0 1 1;0 1 0]), ...
    'lb',[0;0;0],'ub',[100;100;100],'comps',{{'e';'u';'x'}},'metComps',[1;2;1;2;3;3]);
community=normalizeModel(community);community.circuitResolvedEC=strings(3,1);
[cb,cc,cl]=readModelLayout(community,'');communityCfg=resolveModelProfile(setfield(defaultOptions(),'requireIImport',true),cl,cc);
assert(cl.communityIDs&&cb(1)==cb(2)&&cb(1)==cb(5),'Community cargo identities were not recovered.');
assert(cc(2)=="u"&&cc(5)=="x",'Independent community pools were merged.');
currencyCommunity=ismember(cb,communityCfg.currencyBaseIDs);
[Oc,Tc]=buildPools(community,cb,cc,currencyCommunity,max(-community.S,0), ...
    max(community.S,0),true(3,1),communityCfg);
assert(isempty(Oc.rxn)&&height(Tc)==2,'Community import/currency handling failed.');
assert(all(sort(Tc.InsideIndex)==[2;5]),'Import matched the wrong block or a currency species.');
selfTestNativeSBML(direct);
selfTestThreeReactionExtension();
selfTestSingleReactionProduction();
selfTestMinimumProduction();
selfTestCellFreeMerge();
selfTestEmbeddediMM904EC();
fprintf('PASS: embedded iMM904 EC transfer with chemistry guard.\n');
fprintf('PASS: cell-free merge of organelle copies, transport removal and isoenzyme collapse.\n');
fprintf('PASS: import, gross/net distinction, two-step amplification, net-zero rejection,\n');
fprintf('      co-substrate availability, buffering, bounds, reverse policy, depth/budgets,\n');
fprintf('      GPR/EC parsing, gating-evidence filtering and transport/core accounting.\n');
fprintf('      Missing-annotation screening, metadata compartments and separate community pools.\n');
fprintf('      Same-I removal joins, primitive scaling, pair-loss factors, reaction reuse and cross-talk.\n');
fprintf('      EC-qualified rails, distinct EC matching, shared-gene rejection and provisional gate interfaces.\n');
fprintf('      Native SBML without Java, namespace prefixes, entities and nonconstant-stoichiometry rejection.\n');
end

function selfTestSingleReactionProduction()
% Optional synthetic regressions only; no real model or network is queried.
cfg=defaultOptions();cfg.requireIImport=true;cfg.insideCompartment="c";cfg.maxStatesPerI=Inf;
cfg.productionMode="single_reaction"; % optional strict mode is tested explicitly
R=sparse(4,2);P=sparse(4,2);R(1,1)=1;P(2,1)=1;
R([2 3],2)=1;P(2,2)=2;P(4,2)=1;
m=struct('mets',{{'i[e]';'i[c]';'x[c]';'w[c]'}},'rxns',{{'T_seed';'R1'}}, ...
    'S',P-R,'lb',[0;0],'ub',[100;100],'metFormulas',{{'C';'C';'CH2O';'H2O'}});
m=normalizeModel(m);m.circuitResolvedEC=["";"1.1.1.1"];
[base,comp]=splitCompartments(string(m.mets));m.circuitBaseIDs=base;m.circuitCompartment=comp;
currency=false(4,1);[O,T]=buildPools(m,base,comp,currency,R,P,true(2,1),cfg);
[cores,~,audit]=findCores(2,O,currency,cfg,m);
assert(numel(cores)==1&&numel(cores.path)==1&&audit.I_GrossConsumed==1&&audit.I_GrossProduced==2, ...
    'The direct I+X->2I+W reaction was not recovered as one production reaction.');
assert(audit.AtomBalanceStatus=="balanced_atoms"&&cores.net(2)==1);
gates=cellTable(cell(0,11),gateColumns(),repmat({'string'},1,11));
[C,D]=assembleCandidates(cores,O,T,gates,m,cfg);
pairs=describeImportProductionPairs(C,D,cores,O,T,m);
assert(height(pairs)==1&&pairs.TransporterID=="T_seed"&&pairs.R1_ReactionID=="R1" ...
    &&pairs.X_SubstrateIDs=="x[c]"&&pairs.W_ProductIDs=="w[c]"&&pairs.SingleProductionReaction, ...
    'The importer, production reaction, X and W were not reported as separate roles.');
assert(contains(pairs.ImportGateStatus,"unverified")&&contains(pairs.ProductionControlStatus,"not_established"), ...
    'Missing S-import control or S-dependent R1 stopping was promoted to evidence.');
% requireIImport=false: no Block A. The core must still yield one candidate,
% carrying TransportRow 0, and every downstream join must tolerate it.
noImp=cfg;noImp.requireIImport=false;
[Cn,Dn]=assembleCandidates(cores,O,T,gates,m,noImp);
assert(height(Cn)==1&&Dn.TransportRow==0&&Cn.TransportReactionID=="none_no_I_import", ...
    'With requireIImport=false a core must give exactly one import-free candidate.');
assert(Cn.CoreNetEquation==Cn.ModuleNetEquation, ...
    'Without an import the module net equation must equal the core net equation.');
pairsN=describeImportProductionPairs(Cn,Dn,cores,O,T,m);
assert(height(pairsN)==1&&pairsN.TransporterID=="none_no_I_import"&&pairsN.TransportDirection==0 ...
    &&pairsN.R1_ReactionID=="R1"&&pairsN.X_SubstrateIDs=="x[c]", ...
    'The import-free pair lost its production roles.');
assert(contains(pairsN.ProductionControlStatus,"not_established"), ...
    'The import-free candidate must still disclaim the k*S*I law.');
[SGn,sstN]=findSGatedCoreImports(cores,O,T,m,currency,noImp);
assert(isstruct(sstN)&&sstN.cores==numel(cores), ...
    'The S-gate stage must still run when no I import is searched.');
assert(all(SGn.CoreIndex>=1), 'S-gate rows must reference a real core.');
% Same net equation, different gross chemistry: X->I cannot be relabelled.
netOnly=O;netOnly.R(2,1)=0;netOnly.P(2,1)=1;netOnly.N=netOnly.P-netOnly.R;
[none,~,a]=findCores(2,netOnly,currency,cfg,m);
assert(isempty(none)&&all(a.Status=="not_primitive_1I_consumed_2I_produced"), ...
    'A net X->I reaction was accepted by inventing an I catalyst.');
unknown=O;unknown.grossKnown(:)=false;none=findCores(2,unknown,currency,cfg,m);
assert(isempty(none),'Unknown gross participation was accepted as explicit autocatalysis.');
pair=O;pair.R(2,1)=2;pair.P(2,1)=4;pair.N=pair.P-pair.R;
none=findCores(2,pair,currency,cfg,m);
assert(isempty(none),'A primitive two-copy event was rescaled into one-copy autocatalysis.');
scaled=O;scaled.R=2*O.R;scaled.P=2*O.P;scaled.N=scaled.P-scaled.R;
[scaledCore,~,a]=findCores(2,scaled,currency,cfg,m);
assert(numel(scaledCore)==1&&abs(scaledCore.flux-0.5)<1e-10&&a.PrimitiveScale==0.5, ...
    'Uniform model-column scaling changed the primitive reaction identity.');
noX=O;noX.R(3,1)=0;noX.N=noX.P-noX.R;
[none,~,a]=findCores(2,noX,currency,cfg,m);
assert(isempty(none)&&a.Status=="no_additional_consumed_X_substrate",'I duplication without a consumed precursor was accepted.');
otherI=m;otherI.circuitBaseIDs(3)="i";
[none,~,a]=findCores(2,O,currency,cfg,otherI);
assert(isempty(none)&&a.Status=="another_I_pool_cannot_be_counted_as_X", ...
    'Another I pool was treated as a new precursor chemical.');
unbalanced=m;unbalanced.metFormulas{3}='C2H2O';
[none,~,a]=findCores(2,O,currency,cfg,unbalanced);
assert(isempty(none)&&a.AtomBalanceStatus=="unbalanced_atoms",'Known atom imbalance was accepted.');
missing=m;missing.metFormulas{3}='';
[partial,~,a]=findCores(2,O,currency,cfg,missing);
assert(numel(partial)==1&&contains(a.AtomBalanceStatus,"unresolved"), ...
    'Missing formula evidence was silently called balanced or lost from the audit.');
strict=cfg;strict.requireResolvedProductionAtomBalance=true;
none=findCores(2,O,currency,strict,missing);assert(isempty(none));
wrong=T;wrong.InsideIndex(:)=3;
[none,~]=assembleCandidates(cores,O,wrong,gates,m,cfg);
assert(isempty(none),'An importer feeding a different intracellular state was joined.');
rev=m;rev.lb(2)=-100;rev.ub(2)=0;rr=R;pp=P;rr(:,2)=P(:,2);pp(:,2)=R(:,2);rev.S=pp-rr;
[Or,~]=buildPools(rev,base,comp,currency,rr,pp,true(2,1),cfg);
reversed=findCores(2,Or,currency,cfg,rev);
assert(numel(reversed)==1&&Or.direction(reversed.path)==-1,'An allowed reverse production orientation was missed.');
strict=cfg;strict.allowReverse=false;
[Of,~]=buildPools(rev,base,comp,currency,rr,pp,true(2,1),strict);
none=findCores(2,Of,currency,strict,rev);assert(isempty(none),'A reverse reaction bypassed the forward-only policy.');
end

function selfTestMinimumProduction()
% Two known realizations: I+X->Y->2I+W (2 reactions), and
% I+X->Z->T->2I+W (3 reactions). Tests use only this synthetic pool.
cfg=defaultOptions();cfg.requireIImport=true;cfg.productionMode="multi_reaction";cfg.maxCoreReactions=3;
cfg.maxStatesPerI=Inf;cfg.productionMinimumOnly=false;
R=sparse(7,5);P=sparse(7,5);
R([1 2],1)=1;P(3,1)=1;R(3,2)=1;P(1,2)=2;P(7,2)=1;
R([1 2],3)=1;P(4,3)=1;R(4,4)=1;P(5,4)=1;
R(5,5)=1;P(1,5)=2;P(7,5)=1;
O=struct('rxn',1:5,'direction',ones(1,5),'R',R,'P',P,'N',P-R, ...
    'consumers',{cell(7,1)},'grossKnown',true(1,5));
for k=1:7,O.consumers{k}=find(R(k,:)>0);end
currency=false(7,1);
[cores,stats,~,depths]=findCores(1,O,currency,cfg);
assert(isequal(arrayfun(@(c)numel(c.path),cores),[2 3]), ...
    'The search did not enumerate the shorter module before its longer alternative.');
assert(stats.smallestFound==2&&stats.minimumProven&&stats.minimumLayerComplete ...
    &&stats.completedDepth==3&&isequal(depths.ReactionCount,[1;2;3]), ...
    'A completed smaller-layer search failed to establish the known two-reaction minimum.');
cfg.productionMinimumOnly=true;
[shortest,stats,~,depths]=findCores(1,O,currency,cfg);
assert(numel(shortest)==1&&numel(shortest.path)==2&&stats.minimumProven ...
    &&height(depths)==2,'Minimum-only mode searched or retained a larger layer.');
cfg.productionMinimumOnly=false;cfg.maxStatesPerI=4;
[capped,stats]=findCores(1,O,currency,cfg);
assert(numel(capped)==1&&numel(capped.path)==2&&stats.capped&&stats.minimumProven ...
    &&stats.minimumLayerComplete&&stats.completedDepth==2, ...
    'A later-layer budget hit incorrectly invalidated the already established smaller minimum.');
cfg.maxStatesPerI=2;
[none,stats]=findCores(1,O,currency,cfg);
assert(isempty(none)&&stats.capped&&~stats.minimumProven, ...
    'An unfinished small-layer search was labelled a proven minimum.');
baseStats=struct('states',0,'closures',0,'capped',false,'unsolved',1);
uncertain=finishProductionMinimum(baseStats,cores(1),[true;true],[1;0],"synthetic_scope");
assert(~uncertain.minimumProven&&uncertain.minimumStatus=="smallest_found_shorter_search_incomplete", ...
    'A solver-unresolved smaller module was treated as ruled out.');
ties=finishProductionMinimum(baseStats,cores(1),[true;false],[0;0],"synthetic_scope");
assert(ties.minimumProven&&~ties.minimumLayerComplete, ...
    'Minimum-size proof was confused with enumeration of all equal-size alternatives.');
end

function selfTestNativeSBML(direct)
assert(uniqueIDMatch("M_x_u",["i[e]";"x[u]";"x[x]"],'M_')==2, ...
    'SBML underscore to COBRA bracket compartment matching failed.');
assert(isempty(uniqueIDMatch("M_x_u",["x[c]";"x[x]"],'M_')), ...
    'Compartment matching merged distinct metabolite pools.');
% These tests run unconditionally, including a MATLAB -nojvm session.
rawFile=[tempname '.xml'];rawCleanup=onCleanup(@()delete(rawFile)); %#ok<NASGU>
xml=['<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core">' ...
    '<model name="Synthetic &amp; control > fixture">' ...
    '<!-- <reaction id="COMMENT_ONLY"/> -->' ...
    '<notes><![CDATA[<reaction id="CDATA_ONLY"/>]]></notes>' ...
    '<listOfSpecies><species id="M_i_LSQBKT_c_RSQBKT_" compartment="c"/>' ...
    '<species id="M_x_LSQBKT_c_RSQBKT_" compartment="c"/>' ...
    '<species id="M_w_LSQBKT_c_RSQBKT_" compartment="c"/></listOfSpecies>' ...
    '<listOfReactions><reaction id="R_A"><annotation>' ...
    '<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" ' ...
    'xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">' ...
    '<rdf:Description><bqbiol:is><rdf:Bag>' ...
    '<rdf:li rdf:resource="http://identifiers.org/bigg.reaction/A"/>' ...
    '<rdf:li rdf:resource="http://identifiers.org/ec-code/1.1.1.1"/>' ...
    '</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation><listOfReactants>' ...
    '<speciesReference species="M_i_LSQBKT_c_RSQBKT_" stoichiometry="1"/>' ...
    '<speciesReference species="M_x_LSQBKT_c_RSQBKT_" stoichiometry="1"/>' ...
    '</listOfReactants><listOfProducts>' ...
    '<speciesReference species="M_i_LSQBKT_c_RSQBKT_" stoichiometry="2"/>' ...
    '<speciesReference species="M_w_LSQBKT_c_RSQBKT_" stoichiometry="1"/>' ...
    '</listOfProducts></reaction></listOfReactions></model></sbml>'];
prefixed=regexprep(xml,'<(/?)(sbml|model|notes|listOfSpecies|species|listOfReactions|reaction|annotation|listOfReactants|listOfProducts|speciesReference)(?=[\s>/])','<$1s:$2');
prefixed=strrep(prefixed,'<s:sbml xmlns="','<s:sbml xmlns:s="');
for variant=1:3
    payload=xml;if variant==2,payload=prefixed;end
    if variant==3,payload=strrep(xml,'"',char(39));end
    rf=fopen(rawFile,'w');assert(rf>0);fwrite(rf,payload,'char');fclose(rf);
    meta=readSBMLNative(rawFile);
    assert(meta.modelName=="Synthetic & control > fixture"&&numel(meta.speciesIDs)==3 ...
        &&numel(meta.reactions)==1,'Native XML scope/quote/entity parsing regression.');
    [ar,ap,known,~,testEC,testBiGG,info]=readGrossStoichiometry(direct,rawFile,meta);
    assert(known(2)&&ar(2,2)==1&&ap(2,2)==2&&testBiGG(2)=="A" ...
        &&testEC(2)=="1.1.1.1"&&info.parsedReactionCount==1 ...
        &&info.matchedReactionCount==1&&info.error=="", ...
        'Native SBML gross stoichiometry or RDF identifier recovery failed.');
end
assert(strcmp(nativeXMLDecode('&amp;lt;'),'&lt;'),'XML entities were decoded more than once.');
assert(strcmp(nativeXMLDecode('A&#32;&#x42;'),'A B'),'Numeric XML entity decoding failed.');
payload=strrep(xml,'stoichiometry="1"','stoichiometry="1" constant="false"');
rf=fopen(rawFile,'w');assert(rf>0);fwrite(rf,payload,'char');fclose(rf);
[~,~,known,~,~,ids,info]=readGrossStoichiometry(direct,rawFile);
assert(~known(2)&&ids(2)=="A"&&info.error=="", ...
    'Nonconstant stoichiometry was accepted as exact, or its source ID was lost.');
payload=strrep(xml,'<speciesReference species="M_x_LSQBKT_c_RSQBKT_" stoichiometry="1"/>', ...
    ['<speciesReference species="M_x_LSQBKT_c_RSQBKT_" stoichiometry="1">' ...
    '<stoichiometryMath><math><cn>1</cn></math></stoichiometryMath></speciesReference>']);
rf=fopen(rawFile,'w');assert(rf>0);fwrite(rf,payload,'char');fclose(rf);
[~,~,known,~,~,~,info]=readGrossStoichiometry(direct,rawFile);
assert(~known(2)&&info.error=="",'Unevaluated stoichiometryMath was promoted to exact gross stoichiometry.');
payload=strrep(xml,'species="M_x_LSQBKT_c_RSQBKT_"','id="sr_test" species="M_x_LSQBKT_c_RSQBKT_"');
payload=strrep(payload,'<listOfReactions>', ...
    ['<listOfInitialAssignments><initialAssignment symbol="sr_test">' ...
    '<math><cn>1</cn></math></initialAssignment></listOfInitialAssignments><listOfReactions>']);
rf=fopen(rawFile,'w');assert(rf>0);fwrite(rf,payload,'char');fclose(rf);
[~,~,known,~,~,~,info]=readGrossStoichiometry(direct,rawFile);
assert(~known(2)&&info.error=="",'Assigned stoichiometry was promoted to exact gross stoichiometry.');
payload=strrep(xml,'<listOfReactions>','<listOfReactions><reaction id="UNMATCHED_EMPTY"/>');
rf=fopen(rawFile,'w');assert(rf>0);fwrite(rf,payload,'char');fclose(rf);
[~,~,known,~,~,~,info]=readGrossStoichiometry(direct,rawFile);
assert(known(2)&&info.parsedReactionCount==2&&info.matchedReactionCount==1, ...
    'A self-closing reaction corrupted the following reaction scope.');
end

function selfTestThreeReactionExtension()
cfg=defaultOptions();cfg.requireIImport=true;cfg.insideCompartment="c";cfg.maxCoreReactions=3;cfg.maxStatesPerI=Inf;
cfg.productionMode="multi_reaction";
m=struct('mets',{{'i[e]';'i[c]';'x[c]';'y[c]';'w[c]';'l[e]';'l[c]'; ...
    'p[c]';'z[c]';'q[c]';'r[c]';'i[m]'}}, ...
    'rxns',{{'Ts';'A';'B';'Tf';'L';'Q2';'scaled';'Q1'}}, ...
    'lb',zeros(8,1),'ub',ones(8,1)*100);
R=sparse(12,8);P=sparse(12,8);
R(1,1)=1;P(2,1)=1; % S-controlled import topology
R([2 3],2)=1;P(4,2)=1; % I + X -> Y
R(4,3)=1;P(2,3)=2;P(5,3)=1; % Y -> 2I + W
R(6,4)=1;P(7,4)=1; % F-controlled co-substrate import topology
R([2 7],5)=1;P(8,5)=1; % I + L -> P
R(2,6)=2;R(10,6)=1;P(9,6)=1; % alternative 2I + Q -> Z
R([2 3],7)=2;P(11,7)=2; % scaled I+X->R, not a pair-loss event
R(2,8)=2;R(10,8)=1;P([2 9],8)=1; % exact gross 2I+Q -> I+Z
m.S=P-R;m=normalizeModel(m);m.circuitResolvedEC=strings(8,1);
[base,comp]=splitCompartments(string(m.mets));m.circuitCompartment=comp;currency=false(12,1);
[O,T]=buildPools(m,base,comp,currency,R,P,true(8,1),cfg);
[cores,~]=findCores(2,O,currency,cfg);
assert(numel(cores)==1,'Synthetic production core was lost.');
g=cellTable(cell(0,11),gateColumns(),repmat({'string'},1,11));
[prod,details]=assembleCandidates(cores,O,T,g,m,cfg);
[L,Q,~,~]=findRemovalTerms(2,O,T,m,comp,currency,cfg);
assert(height(Q)==2 && nnz(Q.ExactFormalStoichiometry)==1,'Exact/alternative pair losses were conflated.');
assert(~any(Q.ReactionID=="scaled"),'An arbitrary doubled equation was accepted as a pair event.');
[scale,ok]=primitiveReactionScale(R(:,6)/2,P(:,6)/2);
assert(ok&&abs(scale-2)<1e-9,'Fractionally scaled pair equation was not recovered.');
[C,~]=connectAndersonTerms(prod,details,cores,O,T,L,Q,m,comp,currency,R,P,cfg);
assert(height(C)==2&&nnz(C.ExactFormalQuadratic)==1,'The three same-I components were not joined.');
assert(all(sort(C.QuadraticLossPerEvent)==[1;2]),'The factor of two for alternative pair loss was lost.');
assert(all(contains(C.RateStatus,"unverified")),'Unverified kinetics were promoted to a solution.');
assert(all(L.FSpeciesID(L.FMode=="F_gated_cosubstrate_import")==""), ...
    'The consumed co-substrate was incorrectly assigned as regenerated F.');
exactCfg=cfg;exactCfg.allowPairLossAlternative=false;
[Ce,~]=connectAndersonTerms(prod,details,cores,O,T,L,Q,m,comp,currency,R,P,exactCfg);
assert(height(Ce)==1&&Ce.ExactFormalQuadratic,'Exact-only mode admitted the alternative pair loss.');
wrongPool=L;wrongPool.I_Index(:)=12;
[Cw,~]=connectAndersonTerms(prod,details,cores,O,T,wrongPool,Q,m,comp,currency,R,P,cfg);
assert(isempty(Cw),'A different intracellular I pool was joined.');
reused=L;linked=find(reused.FMode=="F_gated_cosubstrate_import");
reused.FTransportRow(linked)=find(T.ModelIndex==1,1);
[Cr,sr]=connectAndersonTerms(prod,details,cores,O,T,reused,Q,m,comp,currency,R,P,cfg);
assert(isempty(Cr)&&sr.reactionReuseRejected==2,'Reused production/removal transport was admitted.');
crossO=O;linO=L.OrientationIndex(linked(1));crossO.P(3,linO)=1;crossO.N=crossO.P-crossO.R;
[Cx,sx]=connectAndersonTerms(prod,details,cores,crossO,T,L,Q,m,comp,currency,R,P,cfg);
assert(isempty(Cx)&&sx.crossTalkRejected==2,'A sink feeding the production core bypassed isolation.');
capCfg=cfg;capCfg.maxTripletTests=1;
[Cc,sc]=connectAndersonTerms(prod,details,cores,O,T,L,Q,m,comp,currency,R,P,capCfg);
assert(height(Cc)==1&&sc.testLimitReached,'The connection budget was not reported.');
selfTestBranchExtension(m,R,P);
end

function selfTestBranchExtension(cancelModel,cancelR,cancelP)
% Two independent three-reaction rails, plus negative EC/GPR/interface cases.
cfg=defaultOptions();cfg.requireIImport=true;cfg.insideCompartment="c";cfg.maxCoreReactions=3;cfg.maxStatesPerI=Inf;
cfg.productionMode="multi_reaction";
b=struct('mets',{{'u[c]';'fm[c]';'fs[c]';'fw[c]';'sm[c]';'ss[c]';'sw[c]'}}, ...
    'rxns',{{'Ff';'Fd';'Fr';'Sf';'Sd';'Sr'}},'lb',zeros(6,1),'ub',ones(6,1)*100);
R=sparse(7,6);P=sparse(7,6);
R(1,1)=1;P(2,1)=1;R(2,2)=1;P(3,2)=1;R(3,3)=1;P(4,3)=1;
R(1,4)=1;P(5,4)=1;R(5,5)=1;P(6,5)=1;R(6,6)=1;P(7,6)=1;
b.S=P-R;b=normalizeModel(b);b.circuitResolvedEC=strings(6,1);
[base,comp]=splitCompartments(string(b.mets));currency=false(7,1);
[O,~]=buildPools(b,base,comp,currency,R,P,true(6,1),cfg);
E=struct('codes',["1.1.1.1";"1.1.1.2";"1.1.1.3";"1.1.1.4";"1.1.1.5";"1.1.1.6"], ...
    'source',repmat("synthetic_fixture",6,1),'reference',repmat("TEST_ONLY",6,1));
[~,M,A,~]=findBranchModules(O,E,b,comp,currency,cfg);
assert(height(M)==2&&height(A)==12&&all(M.DistinctECAssignment), ...
    'The two EC-qualified fast/slow assignments were not recovered.');
assert(~any(M.AllGPRPresent)&&all(isnan(A.Kcat_per_s)), ...
    'EC annotations invented enzyme identities or kinetic measurements.');
missing=E;missing.codes(6)="";
[~,Mm,~,~]=findBranchModules(O,missing,b,comp,currency,cfg);
assert(isempty(Mm),'A branch reaction without an EC was accepted.');
sameEC=E;sameEC.codes(:)="1.1.1.1";
[~,Ms,~,~]=findBranchModules(O,sameEC,b,comp,currency,cfg);
assert(isempty(Ms),'The six-role distinct EC requirement was bypassed.');
sets={ ["1.1.1.1";"1.1.1.2"],"1.1.1.1","1.1.1.3","1.1.1.4","1.1.1.5","1.1.1.6" };
[assignment,ok]=distinctBranchECs(sets,true);
assert(ok&&assignment(1)=="1.1.1.2"&&numel(unique(assignment))==6, ...
    'EC matching failed when a nongreedy assignment was required.');
sharedGene=b;sharedGene.genes={'ga';'gb';'gc';'gd';'ge';'gf'};
sharedGene.grRules={'ga';'gb';'gc';'ga';'ge';'gf'};
[~,Mg,~,sg]=findBranchModules(O,E,sharedGene,comp,currency,cfg);
assert(isempty(Mg)&&sg.knownGeneOverlapRejected>0,'Shared known branch genes were not rejected.');
assert(branchSignatureMatches(sparse([-1;1;0]),["g6p";"f6p";"extra"],"g6p:-1;f6p:1"), ...
    'An exact source reaction signature failed.');
assert(~branchSignatureMatches(sparse([-1;1;1]),["g6p";"f6p";"extra"],"g6p:-1;f6p:1"), ...
    'A modified reaction incorrectly inherited a source EC annotation.');
% Join this module with the independent canceller fixture. Unknown S/F
% controllers must remain provisional despite complete rail EC annotations.
allR=blkdiag(cancelR,R);allP=blkdiag(cancelP,P);
fullModel=struct('mets',{[cancelModel.mets;b.mets]},'rxns',{[cancelModel.rxns;b.rxns]}, ...
    'S',allP-allR,'lb',[cancelModel.lb;b.lb],'ub',[cancelModel.ub;b.ub]);
fullModel=normalizeModel(fullModel);fullModel.circuitResolvedEC=strings(numel(fullModel.rxns),1);
[fullBase,fullComp]=splitCompartments(string(fullModel.mets));fullModel.circuitCompartment=fullComp;
fullCurrency=false(numel(fullModel.mets),1);
[fullO,T]=buildPools(fullModel,fullBase,fullComp,fullCurrency,allR,allP,true(size(allR,2),1),cfg);
[cores,~]=findCores(2,fullO,fullCurrency,cfg);
g=cellTable(cell(0,11),gateColumns(),repmat({'string'},1,11));
[prod,D]=assembleCandidates(cores,fullO,T,g,fullModel,cfg);
[L,Q,~,~]=findRemovalTerms(2,fullO,T,fullModel,fullComp,fullCurrency,cfg);
[C,~]=connectAndersonTerms(prod,D,cores,fullO,T,L,Q,fullModel,fullComp,fullCurrency,allR,allP,cfg);
fullE=struct('codes',[strings(numel(cancelModel.rxns),1);E.codes], ...
    'source',[strings(numel(cancelModel.rxns),1);E.source], ...
    'reference',[strings(numel(cancelModel.rxns),1);E.reference]);
[rails,modules,~,~]=findBranchModules(fullO,fullE,fullModel,fullComp,fullCurrency,cfg);
[W,~]=wireBranchesToCancellers(modules,rails,C,prod,D,cores,fullO,T,L,Q, ...
    fullModel,fullCurrency,allR,allP,cfg);
assert(height(W)==4&&~any(W.BothControlInterfacesMatched), ...
    'Provisional rail wiring was lost or unknown control interfaces were promoted.');
strict=cfg;strict.branchRequireMatchedInterfaces=true;
[Ws,~]=wireBranchesToCancellers(modules,rails,C,prod,D,cores,fullO,T,L,Q, ...
    fullModel,fullCurrency,allR,allP,strict);
assert(isempty(Ws),'Matched-interface-only mode admitted unverified S/F connections.');
selfTestCounterExtension(fullModel,allR,allP,fullE);
end

function selfTestCounterExtension(previous,oldR,oldP,oldE)
% Add independent counter cargo; test loading, overlap and evidence policy.
cfg=defaultOptions();cfg.requireIImport=true;cfg.insideCompartment="c";cfg.maxCoreReactions=3;cfg.maxStatesPerI=Inf;
cfg.productionMode="multi_reaction";
cfg.counterMode="membrane";
cfg.organismName="synthetic_test_organism";
nm=numel(previous.mets);nr=numel(previous.rxns);e=nm+1;c=nm+2;w=nm+3;
R=blkdiag(oldR,sparse([1 0;0 1;0 0]));
P=blkdiag(oldP,sparse([0 0;1 0;0 1]));
m=struct('mets',{[previous.mets;{'counter[e]';'counter[c]';'counter_product[c]'}]}, ...
    'rxns',{[previous.rxns;{'Tcounter';'Lcounter'}]},'S',P-R, ...
    'lb',[previous.lb;0;0],'ub',[previous.ub;100;100]);
m=normalizeModel(m);m.circuitResolvedEC=strings(nr+2,1);
[base,comp]=splitCompartments(string(m.mets));m.circuitCompartment=comp;currency=false(nm+3,1);
[O,T]=buildPools(m,base,comp,currency,R,P,true(nr+2,1),cfg);
E=struct('codes',[oldE.codes;"";"2.7.1.6"], ...
    'source',[oldE.source;"";"synthetic_fixture"],'reference',[oldE.reference;"";"TEST_ONLY"]);
[K,support,~,~]=findCounterComponents(O,T,E,m,base,comp,currency,R,P,cfg);
assert(height(K)==1&&K.C_Index==c&&K.LeakLossPerEvent==1, ...
    'Independent C import/removal was not recovered with the exact pool.');
assert(~K.OutsideFeedAllowed,'Unspecified extracellular counter feed was invented.');
noEC=E;noEC.codes(end)="";
[missing,~,~,~]=findCounterComponents(O,T,noEC,m,base,comp,currency,R,P,cfg);
assert(isempty(missing),'Counter removal without a complete EC was accepted.');
strictFeed=cfg;strictFeed.counterRequireBoundaryFeed=true;
[closed,~,~,~]=findCounterComponents(O,T,E,m,base,comp,currency,R,P,strictFeed);
assert(isempty(closed),'A closed or absent counter feed passed strict boundary mode.');
fed=m;fed.rxns{end+1}='EX_counter';fed.lb(end+1)=-100;fed.ub(end+1)=100;
feedR=R;feedP=P;feedR(e,end+1)=1;feedP(:,end+1)=0;
[allowed,~]=counterBoundaryFeed(e,fed,feedR,feedP);
assert(allowed,'Allowed uptake through a negative exchange flux was lost.');
[cores,~]=findCores(2,O,currency,cfg);
g=cellTable(cell(0,11),gateColumns(),repmat({'string'},1,11));
[prod,D]=assembleCandidates(cores,O,T,g,m,cfg);
[L,Q,~,~]=findRemovalTerms(2,O,T,m,comp,currency,cfg);
[C,~]=connectAndersonTerms(prod,D,cores,O,T,L,Q,m,comp,currency,R,P,cfg);
[G,~]=readCounterGateEvidence('',K,C,cfg);
[X,~]=connectCounters(K,support,G,C,D,cores,O,T,L,Q,m,base,currency,R,P,cfg);
assert(height(X)==height(C)&&all(X.I_AbsentFromCounterChemistry)&&~any(X.GateEvidenceRows), ...
    'Counter joining lost same-I cancellers or promoted an unsupported gate.');
strict=cfg;strict.counterRequireGateEvidence=true;
[none,~]=connectCounters(K,support,G,C,D,cores,O,T,L,Q,m,base,currency,R,P,strict);
assert(isempty(none),'Gate-evidence-only mode admitted an unsupported I gate.');
wrongPool=K;wrongPool.InsideCompartment(:)="m";
[none,~]=connectCounters(wrongPool,support,G,C,D,cores,O,T,L,Q,m,base,currency,R,P,cfg);
assert(isempty(none),'Counter and I from different inside blocks were joined.');
badParts=support;badParts{1}=union(badParts{1},2);
[none,rejected]=connectCounters(K,badParts,G,C,D,cores,O,T,L,Q,m,base,currency,R,P,cfg);
assert(isempty(none)&&rejected.signalInChemistryRejected>0,'Counter chemistry was allowed to contain I.');
badBase=base;badBase(e)=base(2);
[none,rejected]=connectCounters(K,support,G,C,D,cores,O,T,L,Q,m,badBase,currency,R,P,cfg);
assert(isempty(none)&&rejected.signalInChemistryRejected>0,'Another compartment of I was accepted as counter cargo.');
overlap=support;overlap{1}=union(overlap{1},3);
[none,rejected]=connectCounters(K,overlap,G,C,D,cores,O,T,L,Q,m,base,currency,R,P,cfg);
assert(isempty(none)&&rejected.crossTalkRejected>0,'Counter feeding a core intermediate passed isolation.');
reuse=K;reuse.LeakReactionIndex(:)=1;
[none,rejected]=connectCounters(reuse,support,G,C,D,cores,O,T,L,Q,m,base,currency,R,P,cfg);
assert(isempty(none)&&rejected.reactionReuseRejected>0,'Counter/canceller reaction reuse was accepted.');
[rails,M,~,~]=findBranchModules(O,E,m,comp,currency,cfg);
[W,~]=wireBranchesToCancellers(M,rails,C,prod,D,cores,O,T,L,Q,m,currency,R,P,cfg);
[Y,~]=wireCountersToBranches(X,K,support,W,M,rails,m,currency,cfg);
assert(height(Y)==height(W)&&~any(Y.UpstreamSFInterfacesMatched)&&~any(Y.CounterGateReferenceLinked), ...
    'Full counter wiring was lost or unsupported interfaces were promoted.');
sameEC=K;sameEC.LeakECs(:)="1.1.1.1";
[none,rejected]=wireCountersToBranches(X,sameEC,support,W,M,rails,m,currency,cfg);
assert(isempty(none)&&rejected.ECOverlapRejected>0,'Counter EC overlap with an assigned branch EC was accepted.');
% Evidence rows can be linked, but must not become independently verified.
path=[tempname '.tsv'];cleanup=onCleanup(@()delete(path)); %#ok<NASGU>
names=counterGateColumns();row={C.I_ID(1),K.TransportReactionID(1),string(K.TransportDirection(1)), ...
    K.C_ID(1),"inside","opening","reversible_binding","true","TEST_PROTEIN", ...
    "TEST_MAPPING_ONLY",cfg.organismName,"TEST_MECHANISM_ONLY"};
evidence=cellTable(row,names,repmat({'string'},1,numel(names)));
writetable(evidence,path,'FileType','text','Delimiter','\t');
[linked,audit]=readCounterGateEvidence(path,K,C,cfg);
assert(height(linked)==1&&contains(audit.Status,"not_independently_verified"), ...
    'Supplied reference linkage was lost or upgraded to independently verified.');
evidence.I_NotConsumed(:)="false";writetable(evidence,path,'FileType','text','Delimiter','\t');
[linked,~]=readCounterGateEvidence(path,K,C,cfg);
assert(isempty(linked),'An I-consuming gate passed the evidence filter.');
evidence.I_NotConsumed(:)="true";evidence.SignalSide(:)="outside";
writetable(evidence,path,'FileType','text','Delimiter','\t');
[linked,~]=readCounterGateEvidence(path,K,C,cfg);
assert(isempty(linked),'An outside-only gate was linked to intracellular I.');
assert(w==numel(m.mets),'Counter fixture size mismatch.');
selfTestConsumingCounterExtension(previous,oldR,oldP,oldE);
end

function selfTestConsumingCounterExtension(previous,oldR,oldP,oldE)
cfg=defaultOptions();cfg.requireIImport=true;cfg.insideCompartment="c";cfg.maxCoreReactions=3;cfg.maxStatesPerI=Inf;
cfg.productionMode="multi_reaction";
nm=numel(previous.mets);nr=numel(previous.rxns);cargo=nm+1;
R=blkdiag(oldR,sparse(2,2));P=blkdiag(oldP,sparse(2,2));
R(2,nr+1)=1;P(cargo,nr+1)=1;R(cargo,nr+2)=1;P(nm+2,nr+2)=1;
m=struct('mets',{[previous.mets;{'counter_total[c]';'counter_waste[c]'}]}, ...
    'rxns',{[previous.rxns;{'I_to_counter';'counter_removal'}]},'S',P-R, ...
    'lb',[previous.lb;0;0],'ub',[previous.ub;100;100]);
m=normalizeModel(m);m.circuitResolvedEC=strings(nr+2,1);
[base,comp]=splitCompartments(string(m.mets));m.circuitCompartment=comp;currency=false(nm+2,1);
[O,T]=buildPools(m,base,comp,currency,R,P,true(nr+2,1),cfg);
E=struct('codes',[oldE.codes;"2.7.1.6";"3.1.3.1"], ...
    'source',[oldE.source;"synthetic_fixture";"synthetic_fixture"], ...
    'reference',[oldE.reference;"TEST_ONLY";"TEST_ONLY"]);
[cores,~]=findCores(2,O,currency,cfg);
g=cellTable(cell(0,11),gateColumns(),repmat({'string'},1,11));
[prod,D]=assembleCandidates(cores,O,T,g,m,cfg);
[L,Q,~,~]=findRemovalTerms(2,O,T,m,comp,currency,cfg);
[C,~]=connectAndersonTerms(prod,D,cores,O,T,L,Q,m,comp,currency,R,P,cfg);
[F,fp,~,~]=findConsumingCounterFormations(C,O,E,m,base,comp,currency,cfg);
f=find(F.ReactionID=="I_to_counter");assert(numel(f)==1&&F.I_LossPerEvent(f)==1&&F.C_GainPerEvent(f)==1);
[X,xp,~,~]=connectConsumingCounters(F,fp,C,D,cores,O,T,L,Q,E,m,base,currency,R,P,cfg);
assert(height(X)==2*height(C)&&nnz(X.LeakOrientationIndex==0)==height(C), ...
    'Accumulation-only and leaky consuming-counter variants were not distinguished.');
[rails,M,~,~]=findBranchModules(O,E,m,comp,currency,cfg);
[W,~]=wireBranchesToCancellers(M,rails,C,prod,D,cores,O,T,L,Q,m,currency,R,P,cfg);
[Y,~]=wireConsumingCounters(X,F,xp,W,M,rails,O,m,currency,cfg);
assert(height(Y)==2*height(W)&&~any(Y.UpstreamSFInterfacesMatched), ...
    'Consuming-counter branch joining lost variants or promoted unknown interfaces.');
strict=cfg;strict.counterRequireRemoval=true;
[onlyLeaky,~,~,~]=connectConsumingCounters(F,fp,C,D,cores,O,T,L,Q,E,m,base,currency,R,P,strict);
assert(height(onlyLeaky)==height(C)&&all(onlyLeaky.LeakOrientationIndex>0));
noEC=F;noEC.ECResolved(f)=false;noEC.ECs(f)="";
needEC=cfg;needEC.counterRequireFormationEC=true;
[none,~,audit,~]=connectConsumingCounters(noEC,fp,C,D,cores,O,T,L,Q,E,m,base,currency,R,P,needEC);
assert(isempty(none)&&any(audit.Status=="formation_EC_unresolved"),'Missing formation EC was silently accepted when required.');
[relaxed,relaxedParts,~,~]=connectConsumingCounters(noEC,fp,C,D,cores,O,T,L,Q,E,m,base,currency,R,P,cfg);
assert(height(relaxed)==height(X),'A counter without EC was rejected although no counter EC is required.');
[relaxedY,~]=wireConsumingCounters(relaxed,noEC,relaxedParts,W,M,rails,O,m,currency,cfg);
assert(height(relaxedY)==height(Y)&&all(relaxedY.AssignedFormationEC==""), ...
    'Branch wiring rejected or mislabelled an EC-less counter formation.');
bad=fp;bad{f}=union(bad{f},3);
[none,~,audit,~]=connectConsumingCounters(F,bad,C,D,cores,O,T,L,Q,E,m,base,currency,R,P,cfg);
assert(isempty(none)&&any(audit.Status=="formation_shares_noncurrency_outside_I"),'Formation cross-talk was accepted.');
reused=F;reused.ReactionIndex(f)=1;
[none,~,audit,~]=connectConsumingCounters(reused,fp,C,D,cores,O,T,L,Q,E,m,base,currency,R,P,cfg);
assert(isempty(none)&&any(audit.Status=="formation_reuses_canceller_reaction"),'Formation reaction reuse was accepted.');
leakO=find(O.rxn==nr+2&O.direction==1,1);feedback=O;feedback.P(2,leakO)=1;feedback.N=feedback.P-feedback.R;
[none,~,audit,~]=connectConsumingCounters(F,fp,C,D,cores,feedback,T,L,Q,E,m,base,currency,R,P,strict);
assert(isempty(none)&&any(audit.Status=="removal_contains_or_regenerates_I"),'A C sink regenerating I passed isolation.');
overlapEC=X;overlapEC.LeakECs(overlapEC.LeakOrientationIndex>0)="2.7.1.6";
[partial,ws]=wireConsumingCounters(overlapEC,F,xp,W,M,rails,O,m,currency,cfg);
assert(height(partial)==height(W)&&ws.ECRejected>0,'Formation/removal EC independence was bypassed.');
scaledO=O;fo=F.OrientationIndex(f);scaledO.R(:,fo)=2*scaledO.R(:,fo);scaledO.P(:,fo)=2*scaledO.P(:,fo);
[scaled,~,~,~]=findConsumingCounterFormations(C,scaledO,E,m,base,comp,currency,cfg);
sf=find(scaled.ReactionID=="I_to_counter");
assert(scaled.PrimitiveScale(sf)==0.5&&scaled.I_LossPerEvent(sf)/scaled.PrimitiveScale(sf)==2, ...
    'Formation stoichiometry lost the per-model-flux loss after primitive normalization.');
selfTestStep2Circuit(m,R,P,E,F,X,xp,Y,cfg);
end

function selfTestBlock0Alternatives()
% Toy reaction supports only: these are not proposed biochemical reactions.
cfg=defaultOptions();cfg.requireIImport=true;cfg.insideCompartment="c";cfg.block0MaxDepth=4;
cfg.block0MaxAlternatives=Inf;cfg.block0MaxStates=1000;cfg.block0MaxExpansions=10000;
m=struct('mets',{{'camp[c]';'blocked_F[c]';'goal[c]';'a[c]';'b[c]';'x[c]';'w[c]';'q[c]'}}, ...
    'rxns',{{'bad_first';'bad_last';'nonunit_first';'nonunit_last'; ...
    'alt_first';'alt_middle';'alt_last';'parallel_first';'parallel_last'}}, ...
    'S',sparse(8,9),'lb',zeros(9,1),'ub',ones(9,1)*100);
R=sparse(8,9);P=R;
R(1,1)=1;P(2,1)=1;R(2,2)=1;P(3,2)=1;
R([1 6],3)=1;P(4,3)=2;R(4,4)=1;P([3 7],4)=1;
R(1,5)=1;P(5,5)=1;R(5,6)=1;P(8,6)=1;R(8,7)=1;P(3,7)=1;
R([1 6],8)=1;P([4 5],8)=1;R([4 5],9)=1;P([3 7],9)=1;
m.S=P-R;m=normalizeModel(m);m.circuitResolvedEC="1.1.1."+string((1:9)');
[base,comp]=splitCompartments(string(m.mets));m.circuitCompartment=comp;
[O,T]=buildPools(m,base,comp,false(8,1),R,P,true(9,1),cfg);
E=struct('codes',m.circuitResolvedEC,'source',repmat("SYNTHETIC_ONLY",9,1), ...
    'reference',repmat("NOT_BIOLOGICAL_EVIDENCE",9,1));
ctx=struct('model',m,'base',base,'comp',comp,'currency',false(8,1),'O',O,'T',T,'E',E);
G=buildStep2Block0Graph(ctx,cfg);
[best,alternatives,audit]=findStep2Block0(3,[],[2 3],[],ctx,G,cfg,true);
assert(best.found&&numel(alternatives)>=3,'Alternative Block 0 enumeration failed.');
nonunit=false;parallel=false;threeStep=false;
for k=1:numel(alternatives)
    p=alternatives(k);parts=find(any(O.R(:,p.o)+O.P(:,p.o)>1e-8,2));
    assert(~ismember(2,parts),'A forbidden downstream signal survived as a side participant.');
    assert(all(p.flux>0)&&abs(p.net(1)+1)<1e-8&&all(abs(p.net(p.internal))<1e-8));
    if isequal(sort(O.rxn(p.o)),[3 4])
        nonunit=true;assert(norm(p.flux-[1 2],Inf)<1e-8&&abs(p.net(3)-2)<1e-8);
    end
    if isequal(sort(O.rxn(p.o)),[8 9])
        parallel=true;assert(all(ismember([4 5],p.internal)));
    end
    if isequal(sort(O.rxn(p.o)),[5 6 7]),threeStep=true;end
end
assert(nonunit&&parallel&&threeStep,'Nonunit, shared-internal or longer alternatives were missed.');
assert(any(audit.Reason=="shares_downstream_noncurrency"),'Cross-talk exclusions were not audited.');
blocked=findStep2Block0(3,[],[2 3 4 5 8],[],ctx,G,cfg,true);
assert(~blocked.found&&~blocked.topologyReachable,'Disconnected-after-exclusions status is incorrect.');
small=cfg;small.block0MaxAlternatives=1;
limited=findStep2Block0(3,[],[2 3],[],ctx,G,small,true);
assert(limited.found&&limited.alternativeLimitReached,'Alternative retention coverage was lost.');
small=cfg;small.block0MaxDepth=1;
limited=findStep2Block0(3,[],[2 3],[],ctx,G,small,true);
assert(~limited.found&&limited.depthLimitReached,'A depth-limited search was reported as exhaustive.');
conflict=ctx;conflict.E.codes(4)=conflict.E.codes(3);
[~,qualified]=findStep2Block0(3,[],[2 3],[],conflict,G,cfg,true);
for k=1:numel(qualified),assert(~all(ismember([3 4],O.rxn(qualified(k).o))));end
[raw,structural]=findStep2Block0(3,[],[2 3],[],conflict,G,cfg,false);
assert(raw.found&&all(contains([structural.status],"not_EC")), ...
    'Structural-only alternatives were mislabeled as independently EC-qualified.');
unsupported=ctx;first=find(O.rxn==3,1);last=find(O.rxn==4,1);
unsupported.O.R(8,first)=1;unsupported.O.P(8,first)=1;
[~,reason]=validateStep2Block0Module(1,3,[first last],unsupported);
assert(reason=="startupRejected",'An unseeded internal catalyst was accepted.');
end

function selfTestStep2Circuit(previous,oldR,oldP,oldE,F,X,xp,Y,cfg)
% Synthetic single-neuron/report regressions, run only by explicit selftest.
cfg.neuronPlacement="inside"; % preserve the earlier direct-neuron regressions
cfg.step2SearchBlock0=true;cfg.block0RequireEC=true; % original complete-circuit regressions
nm=numel(previous.mets);nr=numel(previous.rxns);a=nm+(1:5);minIndex=find(string(previous.mets)=="u[c]",1);
assert(~isempty(minIndex));input=find(string(previous.mets)=="counter_total[c]",1);
R=blkdiag(oldR,sparse(5,5));P=blkdiag(oldP,sparse(5,5));
R(input,nr+1)=1;P(a(1),nr+1)=1;R(a(1),nr+2)=1;P(a(2),nr+2)=1;
R(a(2),nr+3)=1;P(a(3),nr+3)=1;
R(a(4),nr+4)=1;P(a(5),nr+4)=1;R(a(5),nr+5)=1;P(minIndex,nr+5)=1;
m=struct('mets',{[previous.mets;{'nf1[c]';'nd1[c]';'nw1[c]';'camp[c]';'camp_bridge[c]'}]}, ...
    'rxns',{[previous.rxns;{'Nf1';'Nd1';'Nr1';'CAMP_first';'CAMP_last'}]}, ...
    'S',P-R,'lb',[previous.lb;zeros(5,1)],'ub',[previous.ub;100*ones(5,1)]);
m=normalizeModel(m);m.circuitResolvedEC=strings(nr+5,1);[base,comp]=splitCompartments(string(m.mets));
m.circuitCompartment=comp;currency=false(nm+5,1);[O,T]=buildPools(m,base,comp,currency,R,P,true(nr+5,1),cfg);
E=struct('codes',[oldE.codes;"2.1.1.1";"2.1.1.2";"2.1.1.3";"3.1.1.1";"5.3.1.1"], ...
    'source',[oldE.source;repmat("synthetic_fixture",5,1)],'reference',[oldE.reference;repmat("TEST_ONLY",5,1)]);
[cores,~]=findCores(2,O,currency,cfg);g=cellTable(cell(0,11),gateColumns(),repmat({'string'},1,11));
[prod,D]=assembleCandidates(cores,O,T,g,m,cfg);[L,Q,~,~]=findRemovalTerms(2,O,T,m,comp,currency,cfg);
[C,~]=connectAndersonTerms(prod,D,cores,O,T,L,Q,m,comp,currency,R,P,cfg);
[rails,M,~,~]=findBranchModules(O,E,m,comp,currency,cfg);
% Earlier parent indices remain valid: append-only construction preserves O.
ctx=struct('model',m,'base',base,'comp',comp,'currency',currency,'R',R,'P',P,'O',O,'T',T, ...
    'C',C,'D',D,'cores',cores,'L',L,'Q',Q,'M',M,'rails',rails,'X',X,'F',F,'Y',Y,'XP',{xp},'E',E);
result=searchStep2Circuit(ctx,cfg);
assert(~isempty(result.complete)&&all(result.complete.Block0Distance==2), ...
    'The complete single-neuron circuit with exact cAMP-to-input Block 0 was not found.');
neuronOnly=cfg;neuronOnly.step2SearchBlock0=false;only=searchStep2Circuit(ctx,neuronOnly);
assert(isempty(only.complete)&&~isempty(only.neuronCandidates)&&only.stats.block0Queries==0 ...
    &&any([only.details.neuronStatus]=="neuron_found_Block0_not_searched"), ...
    'Neuron-only mode searched Block 0 or lost the neuron.');
noEC=cfg;noEC.block0RequireEC=false;relaxed=searchStep2Circuit(ctx,noEC);
assert(~isempty(relaxed.complete)&&all(relaxed.complete.Block0Distance==2), ...
    'Block 0 without EC requirement lost the complete circuit.');
for k=1:numel(result.motifs)
    motif=result.motifs(k);assert(motif.n1.input==input&&motif.n1.output==a(2));
    assert(motif.block0.mets(1)==a(4)&&motif.block0.mets(end)==minIndex);
    assert(numel(unique(motif.assignedECs))==numel(motif.assignedECs));
end
assert(all(result.complete.OutputID==string(m.mets{a(2)})),'The single-neuron output must be Nd1.');
missing=ctx;missing.base(a(4))="not_cAMP";bad=searchStep2Circuit(missing,cfg);
assert(isempty(bad.complete),'A full circuit was accepted without the exact cAMP source.');
feedback=ctx;fo=find(O.rxn==nr+1,1);feedback.O.P(2,fo)=1;feedback.O.N=feedback.O.P-feedback.O.R;
bad=searchStep2Circuit(feedback,cfg);
assert(isempty(bad.complete)&&any(bad.entryAudit.WouldRegenerateI),'A neuron regenerating canceller I was accepted.');
cross=ctx;last=find(O.rxn==nr+5,1);cross.O.P(M.FastMidIndex(1),last)=1;cross.O.N=cross.O.P-cross.O.R;
bad=searchStep2Circuit(cross,cfg);
assert(isempty(bad.complete),'A Block 0 byproduct shared with a branch passed isolation.');
path=[tempname '.txt'];cleanup=onCleanup(@()delete(path)); %#ok<NASGU>
writeStep2Report(path,result,ctx,cfg);text=fileread(path);
assert(contains(text,'SIGNAL METABOLITE MAP')&&contains(text,'r_CAMP_1')&&contains(text,'r_Nf_1') ...
    &&contains(text,'Resolved EC:')&&contains(text,'Equation:'),'The requested Step_2 report reaction sections were lost.');
assert(contains(text,' (m_in) [u[c]]')&&contains(text,' (m_Fd) [') ...
    &&contains(text,' (m_c / m_Nin_1) [counter_total[c]]') ...
    &&contains(text,' (m_Nd_1 / m_out) [nd1[c]]'), ...
    'Main signal names or same-metabolite aliases are missing from reaction equations.');
assert(~contains(text,'r_Nf_2')&&~contains(text,'r_Nd_2')&&~contains(text,'r_Nr_2'), ...
    'The single-neuron report still requires a second-neuron role.');
selfTestDirectConsumingChains(ctx,cfg,input,nr+2,a(2));
selfTestExternalNeuron(ctx,cfg,input);
end

function selfTestDirectConsumingChains(ctx,cfg,input,secondStep,oldOutput)
% Split the two-step toy neuron into THREE conversions, then removal.
% This distinguishes the new variable-length search from the old fixed motif.
defaults=defaultOptions();defaults.requireIImport=true;assert(defaults.neuronPlacement=="inside",'The default neuron must consume m_c inside.');
cfg.neuronPlacement="inside";cfg.neuronChainMinDepth=3;cfg.neuronChainMaxDepth=3;
nm=numel(ctx.model.mets);nr=numel(ctx.model.rxns);mid=nm+1;j=nr+1;
R=blkdiag(ctx.R,sparse(1,1));P=blkdiag(ctx.P,sparse(1,1));
P(oldOutput,secondStep)=0;P(mid,secondStep)=1;R(mid,j)=1;P(oldOutput,j)=1;
m=ctx.model;m.mets=[m.mets;{'neuron_bridge[c]'}];m.rxns=[m.rxns;{'Neuron_bridge_step'}];
for key={'rxnNames','grRules','rules'},m.(key{1})=[m.(key{1});{''}];end
m.S=P-R;m.lb=[m.lb;0];m.ub=[m.ub;100];[base,comp]=splitCompartments(string(m.mets));
m.circuitCompartment=comp;t=ctx;t.model=m;t.base=base;t.comp=comp;t.R=R;t.P=P;
t.currency=[ctx.currency;false];t.grossKnown=true(nr+1,1);
t.E.codes=[ctx.E.codes;"5.4.99.101"];t.E.source=[ctx.E.source;"SYNTHETIC_CHAIN_TEST_ONLY"];
t.E.reference=[ctx.E.reference;"TEST_ONLY_NOT_BIOLOGICAL_EVIDENCE"];
t.model.circuitResolvedEC=t.E.codes;
[t.O,~]=buildPools(t.model,base,comp,t.currency,R,P,t.grossKnown,cfg);
[result,outctx]=searchStep2Circuit(t,cfg);
assert(~isempty(result.complete)&&all(result.neuronCandidates.NeuronConversionCount==3), ...
    'The direct consuming search failed to find a three-conversion neuron.');
for k=1:numel(result.motifs)
    v=result.motifs(k);assert(v.n1.input==input&&numel(v.n1.conversionO)==3 ...
        &&numel(v.n1.o)==4&&v.transport.rx==0);
end
short=cfg;short.neuronChainMinDepth=2;short.neuronChainMaxDepth=2;
bounded=searchStep2Circuit(t,short);
assert(all(bounded.neuronCandidates.NeuronConversionCount<=2),'The neuron conversion-depth bound was bypassed.');
feedback=t;first=find(t.O.rxn==secondStep-1,1);assert(~isempty(first));seed=ctx.F.I_Index(1);
feedback.O.P(seed,first)=1;feedback.O.N=feedback.O.P-feedback.O.R;
bad=searchStep2Circuit(feedback,cfg);
assert(isempty(bad.complete)&&isempty(bad.neuronCandidates) ...
    &&any(bad.consumingChains.IInteraction&bad.consumingChains.SearchPass=="diagnostic_only"), ...
    'Feedback consuming chains were either hidden or counted as accepted neurons.');
quiet=cfg;quiet.neuronDiagnoseRejectedChains=false;noDiagnostics=searchStep2Circuit(feedback,quiet);
assert(~any(noDiagnostics.consumingChains.SearchPass=="diagnostic_only"), ...
    'Rejected-chain diagnostics ignored their off switch.');
path=[tempname '.txt'];cleanup=onCleanup(@()delete(path)); %#ok<NASGU>
writeStep2Report(path,result,outctx,cfg);value=fileread(path);
assert(contains(value,'r_Nmid_1_01')&&contains(value,' (m_Nmid_1_01) [neuron_bridge[c]]'), ...
    'Longer-neuron reactions or intermediate signals were left unannotated.');
end

function selfTestExternalNeuron(ctx,cfg,input)
% Synthetic outside exporter and three-neuron-reaction fixture only.
% Never executed during a normal search; no real-model search is run here.
cfg.neuronPlacement="outside";cfg.neuronCompartment="e";
nm=numel(ctx.model.mets);nr=numel(ctx.model.rxns);a=nm+(1:4);
R=blkdiag(ctx.R,sparse(4,4));P=blkdiag(ctx.P,sparse(4,4));
R(input,nr+1)=1;P(a(1),nr+1)=1;
R(a(1),nr+2)=1;P(a(2),nr+2)=1;
R(a(2),nr+3)=1;P(a(3),nr+3)=1;
R(a(3),nr+4)=1;P(a(4),nr+4)=1;
m=ctx.model;m.mets=[m.mets;cellstr([ctx.base(input)+"[e]";"external_nf[e]";"external_nd[e]";"external_w[e]"])];
m.rxns=[m.rxns;{'T_counter_export';'External_Nf';'External_Nd';'External_Nr'}];
for key={'rxnNames','grRules','rules'},m.(key{1})=[m.(key{1});repmat({''},4,1)];end
m.S=P-R;m.lb=[m.lb;zeros(4,1)];m.ub=[m.ub;100*ones(4,1)];m=normalizeModel(m);
[base,comp]=splitCompartments(string(m.mets));m.circuitCompartment=comp;
t=ctx;t.model=m;t.base=base;t.comp=comp;t.currency=[ctx.currency;false(4,1)];
t.R=R;t.P=P;t.grossKnown=true(nr+4,1);
t.O.R=[ctx.O.R;sparse(4,numel(ctx.O.rxn))];
t.O.P=[ctx.O.P;sparse(4,numel(ctx.O.rxn))];t.O.N=t.O.P-t.O.R;
t.O.consumers=[ctx.O.consumers;cell(4,1)];
t.E.codes=[ctx.E.codes;"";"7.1.1.101";"7.1.1.102";"7.1.1.103"];
t.E.source=[ctx.E.source;repmat("SYNTHETIC_TEST_ONLY",4,1)];
t.E.reference=[ctx.E.reference;repmat("TEST_ONLY_NOT_BIOLOGICAL_EVIDENCE",4,1)];
% The old pool remains append-only; preparation adds only outside reactions.
[result,outctx]=searchStep2Circuit(t,cfg);
assert(~isempty(result.complete)&&~isempty(result.neuronCandidates)&&result.stats.exportLinks>0, ...
    'An existing exporter and isolated external neuron were not connected.');
for k=1:numel(result.motifs)
    v=result.motifs(k);assert(v.transport.counter==input&&v.transport.input==a(1)&&v.n1.input==a(1));
    assert(all(comp(v.n1.parts)=="e")&&v.transport.rx==nr+1);
end
assert(isequal(outctx.O.rxn(1:numel(ctx.O.rxn)),ctx.O.rxn),'Outside preparation changed an upstream orientation index.');
noExport=t;noExport.model.ub(nr+1)=0;bad=searchStep2Circuit(noExport,cfg);
assert(isempty(bad.complete)&&isempty(bad.neuronCandidates),'An external neuron bypassed its mandatory exporter.');
wrong=t;wrong.P(a(1),nr+1)=2;wrong.model.S=wrong.P-wrong.R;
bad=searchStep2Circuit(wrong,cfg);
assert(isempty(bad.complete)&&any(bad.transportAudit.Status=="cargo_copies_not_conserved"), ...
    'Unequal counter-copy export was accepted as transport.');
feedback=t;outsideI=find(base==base(ctx.F.I_Index(1))&comp=="e",1);assert(~isempty(outsideI));
feedback.P(outsideI,nr+2)=1;feedback.model.S=feedback.P-feedback.R;
bad=searchStep2Circuit(feedback,cfg);
assert(isempty(bad.complete)&&any(bad.entryAudit.Status=="entry_external_I_feeds_existing_upstream_pool"), ...
    'Outside I that can feed the existing upstream import pool escaped the feedback check.');
noBlock=t;noBlock.base(noBlock.base=="camp")="not_cAMP";
bad=searchStep2Circuit(noBlock,cfg);
assert(isempty(bad.complete)&&~isempty(bad.neuronCandidates), ...
    'Missing Block 0 concealed otherwise valid external-neuron candidates.');
path=[tempname '.txt'];cleanup=onCleanup(@()delete(path)); %#ok<NASGU>
writeStep2Report(path,result,outctx,cfg);value=fileread(path);
assert(contains(value,'r_T_N1')&&contains(value,' (m_c / m_exp) [counter_total[c]]') ...
    &&contains(value,' (m_Nin_1) [counter_total[e]]'), ...
    'The exported counter and external neuron input were merged or left unannotated.');
selfTestNeuronAdapter(t,cfg,input,nr+1,nr+2);
end

function selfTestNeuronAdapter(ctx,cfg,counter,oldExport,externalEntry)
% A two-reaction adapter is REQUIRED to reach the only enabled exporter.
% All annotations here are synthetic fixture data, never biological evidence.
nm=numel(ctx.model.mets);nr=numel(ctx.model.rxns);v=nm+(1:3);j=nr+(1:3);
R=blkdiag(ctx.R,sparse(3,3));P=blkdiag(ctx.P,sparse(3,3));
R(counter,j(1))=1;P(v(1),j(1))=1;
R(v(1),j(2))=1;P(v(2),j(2))=1;
R(v(2),j(3))=1;P(v(3),j(3))=1;
R(:,externalEntry)=0;R(v(3),externalEntry)=1;
m=ctx.model;m.mets=[m.mets;{'adapter_mid[c]';'export_cargo[c]';'export_cargo[e]'}];
m.rxns=[m.rxns;{'Adapter_first';'Adapter_last';'T_adapter_cargo'}];
for key={'rxnNames','grRules','rules'},m.(key{1})=[m.(key{1});repmat({''},3,1)];end
m.S=P-R;m.lb=[m.lb;zeros(3,1)];m.ub=[m.ub;100*ones(3,1)];m.ub(oldExport)=0;
[base,comp]=splitCompartments(string(m.mets));m.circuitCompartment=comp;
t=ctx;t.model=m;t.base=base;t.comp=comp;t.R=R;t.P=P;
t.currency=[ctx.currency;false(3,1)];t.grossKnown=true(nr+3,1);
t.E.codes=[ctx.E.codes;"6.3.1.201";"6.3.1.202";""];
t.E.source=[ctx.E.source;repmat("SYNTHETIC_ADAPTER_TEST_ONLY",3,1)];
t.E.reference=[ctx.E.reference;repmat("TEST_ONLY_NOT_BIOLOGICAL_EVIDENCE",3,1)];
t.model.circuitResolvedEC=t.E.codes;
[t.O,~]=buildPools(t.model,base,comp,t.currency,R,P,t.grossKnown,cfg);
cfg.neuronAdapterMaxDepth=2;
[result,outctx]=searchStep2Circuit(t,cfg);
assert(~isempty(result.complete)&&~isempty(result.adapterLinks) ...
    &&all(result.adapterLinks.AdapterReactionCount==2), ...
    'The two-reaction adapter failed to reach its existing exporter.');
for k=1:numel(result.motifs)
    link=result.motifs(k).transport;
    assert(link.counter==counter&&link.cargo==v(2)&&link.input==v(3) ...
        &&isequal(link.adapterRx,j(1:2))&&result.motifs(k).n1.input==v(3));
end
short=cfg;short.neuronAdapterMaxDepth=1;bad=searchStep2Circuit(t,short);
assert(isempty(bad.complete)&&isempty(bad.adapterLinks),'The adapter depth bound was bypassed.');
feedback=t;seed=ctx.F.I_Index(1);feedback.P(seed,j(1))=1;feedback.model.S=feedback.P-feedback.R;
[feedback.O,~]=buildPools(feedback.model,base,comp,feedback.currency,feedback.R,feedback.P,feedback.grossKnown,cfg);
bad=searchStep2Circuit(feedback,cfg);
assert(isempty(bad.adapterLinks)&&any(bad.adapterAudit.Status=="adapter_couples_to_I_pool"), ...
    'The adapter was allowed to regenerate the upstream canceller I.');
missing=t;missing.E.codes(externalEntry)="";bad=searchStep2Circuit(missing,cfg);
assert(isempty(bad.neuronCandidates)&&~isempty(bad.adapterLinks), ...
    'An unannotated external neuron concealed valid adapter/export interfaces.');
path=[tempname '.txt'];cleanup=onCleanup(@()delete(path)); %#ok<NASGU>
writeStep2Report(path,result,outctx,cfg);value=fileread(path);
assert(contains(value,'r_A_N1_01')&&contains(value,'r_A_N1_02') ...
    &&contains(value,' (m_A1) [adapter_mid[c]]') ...
    &&contains(value,' (m_exp) [export_cargo[c]]') ...
    &&contains(value,' (m_Nin_1) [export_cargo[e]]'), ...
    'The adapter intermediates, export cargo and external input were not annotated separately.');
end

function [model,base,comp,currency,R,P,known,ecComp,metSource,A,MA,info] = mergeCellFreeCompartments( ...
    model,base,comp,currency,R,P,known,ecComp,metSource,layout,cfg)
% CELL-FREE MODE (cfg.cellFreeMerge): one reaction vessel, every component
% supplied, no membranes. Inside each block, all organelle copies of one
% chemical species (same BaseID in cytosol, mitochondrion, peroxisome, ...)
% become ONE species. The outside compartment is never merged, so the
% S-gated I import and the rest of the Anderson topology are unchanged.
% Gross and net stoichiometry are transformed together:
%  * a reaction that only moved species between organelles has no net
%    chemistry left and is removed (carriers, antiporters, ATP-driven
%    pumps, proton carriers, ATP synthase);
%  * an organelle-spanning reaction with real chemistry is kept, with any
%    species that appeared on both sides cancelled;
%  * isoenzyme copies with identical merged chemistry become ONE reaction
%    (the same chemistry in one pot is one enzyme activity). Their bounds
%    are united and their ECs, gene rules and source IDs are pooled.
tol=1e-9;[m,n]=size(model.S);ids=string(model.mets(:));rxIDs=string(model.rxns(:));
base=string(base(:));comp=string(comp(:));
inside=string(cfg.insideCompartment(:));outside=string(cfg.outsideCompartment);
organelle=comp;block=comp;
if layout.communityIDs
    for i=1:m
        t=regexp(char(ids(i)),'^(?:M_)?(?:AUXOS|PROTS)([CMXNRVG])','tokens','once');
        if ~isempty(t),organelle(i)=string(t{1});end
    end
    cytosol="C";
else
    assert(isscalar(inside),'Anderson:CellFree', ...
        'cellFreeMerge on a model without AUXOS/PROTS blocks needs exactly one insideCompartment.');
    move=comp~=""&comp~=outside;block(move)=inside;cytosol=inside;
end
mergeable=ismember(block,inside)&comp~="";

% ---- metabolites: one state per (block, BaseID); cytosolic ID kept -------
key="#"+string((1:m)');key(mergeable)=block(mergeable)+"|"+base(mergeable);
[~,~,group]=unique(key,'stable');group=group(:);nNew=max(group);
repIdx=zeros(nNew,1);copies=zeros(nNew,1);mergedMetIDs=strings(nNew,1);
for g=1:nNew
    members=find(group==g);pick=members(organelle(members)==cytosol);
    if isempty(pick),pick=members;end
    repIdx(g)=pick(1);copies(g)=numel(members);mergedMetIDs(g)=joinText(ids(members),';');
end
G=sparse(group,(1:m)',1,nNew,m);
Rm=G*R;Pm=G*P;newCurrency=currency(repIdx);newCurrency=newCurrency(:);

% ---- reactions that span organelles of one block -------------------------
spanning=false(n,1);cytosolic=false(n,1);
for j=1:n
    parts=find(((R(:,j)+P(:,j))>tol)&mergeable);
    cytosolic(j)=~isempty(parts)&&all(organelle(parts)==cytosol);
    if numel(parts)<2,continue;end
    blocks=block(parts);
    for b=reshape(unique(blocks),1,[])
        if numel(unique(organelle(parts(blocks==b))))>1,spanning(j)=true;break;end
    end
end
action=repmat("kept_unchanged",n,1);keepRx=true(n,1);
for j=reshape(find(spanning),1,[])
    net=full(Pm(:,j)-Rm(:,j));net(abs(net)<tol)=0;
    Rm(:,j)=sparse(max(-net,0));Pm(:,j)=sparse(max(net,0));known(j)=false;
    changed=find(net);
    if isempty(changed)
        keepRx(j)=false;action(j)="removed_intra_cell_transport_no_net_chemistry";
    elseif all(newCurrency(changed))
        keepRx(j)=false;action(j)="removed_intra_cell_transport_currency_only";
    else
        action(j)="kept_organelle_spanning_chemistry_shared_species_cancelled";
    end
end

% ---- identical merged chemistry: one reaction ----------------------------
ec=strings(n,1);
if isfield(model,'circuitResolvedEC'),ec=string(model.circuitResolvedEC(:));end
canon=strings(n,1);sense=ones(n,1);
for j=reshape(find(keepRx),1,[])
    f=cellFreeReactionKey(Rm(:,j),Pm(:,j));r=cellFreeReactionKey(Pm(:,j),Rm(:,j));
    ordered=sort([f;r]);
    if ordered(1)~=f,sense(j)=-1;end
    canon(j)=ordered(1);
end
memberOf=(1:n)';orientation=ones(n,1);duplicates=0;
kept=find(keepRx);
if ~isempty(kept)
    [~,~,dg]=unique(canon(kept),'stable');dg=dg(:);
    for g=1:max(dg)
        members=kept(dg==g);
        if numel(members)<2,continue;end
        % Representative: annotated EC first, then all-cytosolic, then order.
        [~,order]=sortrows([-double(ec(members)~="") -double(cytosolic(members)) members]);
        rep=members(order(1));
        for j=reshape(members,1,[])
            memberOf(j)=rep;orientation(j)=sense(j)*sense(rep);
            if j~=rep
                keepRx(j)=false;action(j)="merged_into_identical_reaction";duplicates=duplicates+1;
            end
        end
        if action(rep)=="kept_unchanged",action(rep)="kept_representative_of_identical_reactions";end
    end
end
lb=model.lb(:);ub=model.ub(:);
for j=reshape(find(memberOf~=(1:n)'),1,[])
    rep=memberOf(j);
    if orientation(j)>0,l=model.lb(j);u=model.ub(j);else,l=-model.ub(j);u=-model.lb(j);end
    lb(rep)=min(lb(rep),l);ub(rep)=max(ub(rep),u);
end

% ---- audits (equations use the ORIGINAL species IDs) ----------------------
oldEq=strings(n,1);
for j=1:n,oldEq(j)=equation(R(:,j),P(:,j),model);end
repID=rxIDs(memberOf);repID(~keepRx & memberOf==(1:n)')="";
A=table(rxIDs,action,repID,orientation,model.lb(:),model.ub(:),spanning,oldEq, ...
    'VariableNames',{'OriginalReactionID','CellFreeAction','RepresentativeID', ...
    'OrientationVsRepresentative','OriginalLB','OriginalUB','SpannedOrganelles','OriginalEquation'});
MA=table(ids(repIdx),base(repIdx),block(repIdx),copies,mergedMetIDs, ...
    'VariableNames',{'CellFreeMetaboliteID','BaseID','Block','OrganelleCopies','MergedOriginalIDs'});

% ---- build the merged model ----------------------------------------------
old=model;keep=find(keepRx);nKeep=numel(keep);
model=cellFreeSubsetRows(model,repIdx,keep,m,n);
model.mets=cellstr(ids(repIdx));
R=Rm(:,keep);P=Pm(:,keep);model.S=P-R;
model.lb=lb(keep);model.ub=ub(keep);known=known(:);known=known(keep);
if isfield(model,'rev')&&numel(model.rev)==nKeep,model.rev=double(model.lb<0);end
mergedRxnIDs=strings(nKeep,1);mergedSource=strings(nKeep,1);
hasSource=isfield(old,'circuitSourceBiGG')&&numel(old.circuitSourceBiGG)==n;
for k=1:nKeep
    rep=keep(k);members=find(memberOf==rep);
    mergedRxnIDs(k)=joinText(rxIDs(members),';');
    if hasSource
        values=string(old.circuitSourceBiGG(members));
        mergedSource(k)=joinText(unique(values(values~=""),'stable'),';');
    end
    if numel(members)<2,continue;end
    model.grRules{k}=char(cellFreeOrJoin(string(old.grRules(members)),' or '));
    model.rules{k}=char(cellFreeOrJoin(string(old.rules(members)),' | '));
    if isfield(old,'rxnGeneMat')&&size(old.rxnGeneMat,1)==n&&size(model.rxnGeneMat,1)==nKeep
        model.rxnGeneMat(k,:)=any(old.rxnGeneMat(members,:)~=0,1);
    end
    for f={'circuitResolvedEC','circuitSBMLEC'}
        if isfield(old,f{1})&&numel(old.(f{1}))==n
            codes=completeECTokens(string(old.(f{1})(members)));
            if ~isempty(codes),model.(f{1})(k)=join(codes,'; ');end
        end
    end
    for f={'circuitECSource','circuitECReference'}
        if isfield(old,f{1})&&numel(old.(f{1}))==n
            values=string(old.(f{1})(members));values=unique(values(values~=""),'stable');
            model.(f{1})(k)=joinText(values,'; ');
        end
    end
    if isfield(old,'circuitSourceBiGG')&&numel(old.circuitSourceBiGG)==n&&model.circuitSourceBiGG(k)==""
        values=string(old.circuitSourceBiGG(members));values=values(values~="");
        if ~isempty(values),model.circuitSourceBiGG(k)=values(1);end
    end
end
if isfield(model,'circuitResolvedEC')&&isfield(model,'circuitECEligible')
    model.circuitECEligible=string(model.circuitResolvedEC(:))~="";
end
base=base(repIdx);comp=comp(repIdx);
if ~layout.communityIDs,comp(mergeable(repIdx))=inside;end
currency=newCurrency;ecComp=ecComp(repIdx);metSource=metSource(repIdx);
model.circuitCompartment=comp;model.circuitBaseIDs=base;
model.circuitMergedMetIDs=mergedMetIDs;model.circuitMergedRxnIDs=mergedRxnIDs;
model.circuitMergedSourceBiGG=mergedSource;
info=struct('metabolitesBefore',m,'metabolitesAfter',nNew,'mergedSpecies',nnz(copies>1), ...
    'reactionsBefore',n,'reactionsAfter',nKeep, ...
    'removedNoNetChemistry',nnz(action=="removed_intra_cell_transport_no_net_chemistry"), ...
    'removedCurrencyOnly',nnz(action=="removed_intra_cell_transport_currency_only"), ...
    'spanningKept',nnz(action=="kept_organelle_spanning_chemistry_shared_species_cancelled"), ...
    'duplicatesMerged',duplicates, ...
    'boundsWidened',nnz(model.lb<old.lb(keep)-tol|model.ub>old.ub(keep)+tol));
end

function k = cellFreeSeedMatch(seedIDs,sourceID,model,j)
% EC seed row for reaction j: its own source ID first, then the source IDs
% of the isoenzyme copies merged into it (cell-free mode only).
k=find(seedIDs==sourceID,1);
if ~isempty(k)||~isfield(model,'circuitMergedSourceBiGG') ...
        ||numel(model.circuitMergedSourceBiGG)~=numel(model.rxns),return;end
alternatives=split(string(model.circuitMergedSourceBiGG(j)),";");
alternatives=alternatives(alternatives~="");
for q=1:numel(alternatives)
    k=find(seedIDs==alternatives(q),1);
    if ~isempty(k),return;end
end
end

function k = cellFreeReactionKey(a,b)
% Exact text key of one gross reaction: substrate side | product side.
ia=find(a);ib=find(b);
k=string(sprintf('%d:%.9g,',[ia(:)';full(a(ia))']))+"|"+ ...
    string(sprintf('%d:%.9g,',[ib(:)';full(b(ib))']));
end

function s = cellFreeOrJoin(values,op)
% Isoenzymes are alternatives: OR their gene rules, dropping empty ones.
values=strtrim(values(:));values=values(~ismissing(values));
values=unique(values(values~=""),'stable');
if isempty(values)
    s="";
elseif isscalar(values)
    s=values;
else
    s="("+join(values,")"+string(op)+"(")+")";
end
end

function model = cellFreeSubsetRows(model,metRows,rxnRows,m,n)
% Keep every per-metabolite and per-reaction model field consistent.
skip={'S','genes','geneNames','comps','compNames','modelName','modelID', ...
    'description','proteins','modelVersion','osenseStr'};
metFields={'mets','b','csense','circuitCompartment','circuitBaseIDs'};
rxnFields={'rxns','lb','ub','c','rev','grRules','rules','subSystems','eccodes','ecCodes', ...
    'ecNumbers','circuitSBMLEC','circuitSourceBiGG','circuitResolvedEC','circuitECSource', ...
    'circuitECReference','circuitECEligible'};
names=fieldnames(model);
for q=1:numel(names)
    f=names{q};
    if any(strcmp(f,skip)),continue;end
    x=model.(f);
    isMet=startsWith(f,'met')||any(strcmp(f,metFields));
    isRxn=startsWith(f,'rxn')||any(strcmp(f,rxnFields));
    if m==n&&~isMet&&~isRxn,continue;end % ambiguous size: leave untouched
    if size(x,1)==m&&(isMet||(~isRxn&&m~=n))
        model.(f)=x(metRows,:);
    elseif size(x,1)==n&&(isRxn||(~isMet&&m~=n))
        model.(f)=x(rxnRows,:);
    elseif isvector(x)&&numel(x)==m&&isMet
        model.(f)=x(metRows);
    elseif isvector(x)&&numel(x)==n&&isRxn
        model.(f)=x(rxnRows);
    end
end
end

function selfTestCellFreeMerge()
% Organelle copies merge per block; transport without chemistry disappears;
% identical isoenzymes collapse with united bounds, ECs and gene rules.
ids={'M_a_e';'M_AUXOSCa_u';'M_AUXOSMa_u';'M_AUXOSCb_u';'M_AUXOSMb_u';'M_AUXOSCc_u'; ...
    'M_AUXOSMd_u';'M_AUXOSCatp_u';'M_AUXOSCadp_u';'M_AUXOSCh_u';'M_AUXOSMh_u'; ...
    'M_PROTSCa_x';'M_PROTSMa_x'};
rx={'E1c';'E1m';'E1rev';'Tcm';'Pump';'Hleak';'Span';'SpanH';'Imp1';'Imp2';'Xblock'};
R=zeros(13,11);P=zeros(13,11);
R(2,1)=1;P(4,1)=1;            % E1c   a[C] -> b[C]
R(3,2)=1;P(5,2)=1;            % E1m   a[M] -> b[M]
R(5,3)=1;P(3,3)=1;            % E1rev b[M] -> a[M]  (written backwards)
R(2,4)=1;P(3,4)=1;            % Tcm   a[C] -> a[M]
R([2 8],5)=1;P([3 9],5)=1;    % Pump  a[C] + atp -> a[M] + adp
R(10,6)=1;P(11,6)=1;          % Hleak h[C] -> h[M]
R([2 6],7)=1;P(7,7)=1;        % Span  a[C] + c[C] -> d[M]
R([2 6 11],8)=1;P([7 10],8)=1;% SpanH a[C] + c[C] + h[M] -> d[M] + h[C]
R(1,9)=1;P(2,9)=1;            % Imp1  a[e] -> a[C]
R(1,10)=1;P(3,10)=1;          % Imp2  a[e] -> a[M]
R(12,11)=1;P(13,11)=1;        % Xblock a[C,x] -> a[M,x]
m=struct('mets',{ids},'rxns',{rx},'S',sparse(P-R),'lb',[0;0;0;0;0;0;0;0;0;0;0], ...
    'ub',[100;100;20;100;100;100;100;100;100;100;100],'comps',{{'e';'u';'x'}}, ...
    'metComps',[1;2;2;2;2;2;2;2;2;2;2;3;3],'genes',{{'g1';'g2';'g3'}}, ...
    'grRules',{{'g1';'g2';'g3';'';'';'';'';'';'';'';''}});
m=normalizeModel(m);
m.circuitResolvedEC=["1.1.1.1";"1.1.1.2";"";"";"";"";"";"";"";"";""];
m.circuitSBMLEC=strings(11,1);m.circuitSourceBiGG=["";"E1m";"";"";"";"";"";"";"";"";""];
[b,c,layout]=readModelLayout(m,'');cfg=resolveModelProfile(setfield(defaultOptions(),'requireIImport',true),layout,c);
cur=ismember(lower(b),lower(cfg.currencyBaseIDs));
[mm,bb,cc,cu,RR,PP,kk,~,~,A,MA,info]=mergeCellFreeCompartments(m,b,c,cur,sparse(R),sparse(P), ...
    true(11,1),repmat("other",13,1),repmat("test",13,1),layout,cfg);
assert(numel(mm.mets)==9&&numel(bb)==9&&numel(cc)==9&&numel(cu)==9&&height(MA)==9, ...
    'Cell-free merge: wrong number of merged species.');
assert(any(string(mm.mets)=="M_AUXOSCa_u")&&~any(string(mm.mets)=="M_AUXOSMa_u") ...
    &&any(string(mm.mets)=="M_AUXOSMd_u")&&any(string(mm.mets)=="M_PROTSCa_x") ...
    &&any(string(mm.mets)=="M_a_e"),'Cell-free merge: wrong representative species or block mixing.');
assert(isequal(string(mm.rxns),["E1c";"Span";"Imp1"])&&size(mm.S,2)==3&&numel(kk)==3, ...
    'Cell-free merge: transports were not removed or identical reactions were not collapsed.');
assert(norm(full(mm.S-(PP-RR)),Inf)<1e-12,'Cell-free merge: gross and net stoichiometry disagree.');
assert(mm.lb(1)==-20&&mm.ub(1)==100,'Cell-free merge: bounds of a reversed isoenzyme were not united.');
codes=completeECTokens(mm.circuitResolvedEC(1));
assert(all(ismember(["1.1.1.1";"1.1.1.2"],codes))&&mm.circuitSourceBiGG(1)=="E1m" ...
    &&mm.circuitMergedSourceBiGG(1)=="E1m"&&mm.circuitMergedRxnIDs(1)=="E1c;E1m;E1rev", ...
    'Cell-free merge: EC numbers or source IDs of isoenzymes were not pooled.');
assert(cellFreeSeedMatch(["X";"E1m"],"",mm,1)==2&&isempty(cellFreeSeedMatch("X","",mm,1)), ...
    'Cell-free merge: EC seed lookup does not try merged source IDs.');
assert(contains(string(mm.grRules{1}),"g1")&&contains(string(mm.grRules{1}),"g3") ...
    &&contains(string(mm.grRules{1})," or "),'Cell-free merge: gene rules were not joined by OR.');
assert(A.CellFreeAction(4)=="removed_intra_cell_transport_no_net_chemistry" ...
    &&A.CellFreeAction(5)=="removed_intra_cell_transport_currency_only" ...
    &&A.CellFreeAction(11)=="removed_intra_cell_transport_no_net_chemistry" ...
    &&A.RepresentativeID(3)=="E1c"&&A.OrientationVsRepresentative(3)==-1, ...
    'Cell-free merge: audit actions are wrong.');
assert(info.removedNoNetChemistry==3&&info.removedCurrencyOnly==1&&info.duplicatesMerged==4, ...
    'Cell-free merge: summary counts are wrong.');
% A merged pool must still expose the import and the merged core step.
[O,T]=buildPools(mm,bb,cc,cu,RR,PP,kk,cfg);
assert(height(T)==1&&T.InsideMetaboliteID=="M_AUXOSCa_u"&&any(O.rxn==1), ...
    'Cell-free merge: import or core admission failed on the merged model.');
end

function [E,status,count] = applyEmbeddediMM904EC(E,status,model,base)
% Counter/neuron/Block 0 EC evidence: the iMM904 authors' own annotation
% (PROTEIN_CLASS), transferred only on an exact iMM904 reaction ID (after
% removing the AUXOS/PROTS wrapper, including the IDs of isoenzyme copies
% merged in cell-free mode) AND a matching net chemistry. Fills only empty
% rows, so model fields, supplied files and verified seeds keep priority.
T=embeddediMM904ECTable();count=0;n=numel(model.rxns);
keys=iMM904ReactionKey(T.SourceReactionID);
map=containers.Map(cellstr(keys),num2cell((1:height(T))'));
ids=string(model.rxns(:));alt=strings(n,1);
if isfield(model,'circuitMergedRxnIDs')&&numel(model.circuitMergedRxnIDs)==n
    alt=string(model.circuitMergedRxnIDs(:));
end
names=lower(regexprep(string(base(:)),'[^A-Za-z0-9]+','_'));
names=regexprep(names,'^_+|_+$','');
for j=reshape(find(E.codes==""),1,[])
    candidates=[ids(j);split(alt(j),";")];candidates=candidates(candidates~="");
    candidates=iMM904ReactionKey(regexprep(regexprep(candidates,'^R_',''),'^(AUXOS|PROTS)',''));
    row=0;mismatch=false;
    for q=1:numel(candidates)
        if ~isKey(map,char(candidates(q))),continue;end
        k=map(char(candidates(q)));
        if iMM904SignatureMatches(model.S(:,j),names,T.ChemicalSignature(k)),row=k;break;end
        mismatch=true;
    end
    if row==0
        if mismatch,status(j)="embedded_iMM904_ID_matched_chemistry_mismatch_EC_withheld";end
        continue;
    end
    E.codes(j)=T.EC(row);E.source(j)="embedded_iMM904_PROTEIN_CLASS_2009";
    E.reference(j)="iMM904 SBML reaction notes PROTEIN_CLASS for "+T.SourceReactionID(row)+ ...
        " (Mo, Palsson & Herrgard 2009, BMC Syst Biol 3:37)";
    status(j)="embedded_iMM904_ID_and_chemistry_matched";count=count+1;
end
end

function key = iMM904ReactionKey(x)
key=regexprep(string(x),'[^A-Za-z0-9]+','_');
key=regexprep(key,'^_+|_+$','');
end

function ok = iMM904SignatureMatches(net,names,signature)
% Same chemistry up to one uniform scale/orientation. Organelle tags in the
% table are ignored, so merged (cell-free) species still match.
ok=false;parts=find(abs(net)>1e-10);
if isempty(parts),return;end
[actualNames,~,g]=unique(names(parts));actual=accumarray(g(:),full(net(parts)));
terms=split(string(signature),";");expectedNames=strings(numel(terms),1);coefficient=zeros(numel(terms),1);
for k=1:numel(terms)
    pair=split(terms(k),":");if numel(pair)~=2,return;end
    expectedNames(k)=regexprep(pair(1),'@[a-z]$','');coefficient(k)=str2double(pair(2));
end
[uniqueExpected,~,h]=unique(expectedNames);expected=accumarray(h(:),coefficient);
keepA=abs(actual)>1e-10;keepE=abs(expected)>1e-10;
actualNames=actualNames(keepA);actual=actual(keepA);
uniqueExpected=uniqueExpected(keepE);expected=expected(keepE);
if isempty(actual)||numel(actualNames)~=numel(uniqueExpected)||~all(ismember(actualNames,uniqueExpected)),return;end
[~,loc]=ismember(actualNames,uniqueExpected);ratio=actual./expected(loc);
ok=all(isfinite(ratio))&&abs(ratio(1))>1e-12&&max(abs(ratio-ratio(1)))<1e-8;
end

function selfTestEmbeddediMM904EC()
% A wrapped and a merged iMM904 ID resolve; wrong chemistry is withheld.
m=struct('mets',{{'AUXOSCg6p[u]';'AUXOSCf6p[u]';'AUXOSCfdp[u]';'AUXOSCpyr[u]';'AUXOSCala_L[u]'; ...
    'AUXOSCakg[u]';'AUXOSCglu_L[u]'}},'rxns',{{'AUXOSPGI';'AUXOSPFK';'AUXOSALT_merged'}}, ...
    'S',sparse([-1 0 0;1 -1 0;0 1 0;0 0 1;0 0 -1;0 0 -1;0 0 1]),'lb',[-10;0;-10],'ub',[10;10;10]);
m=normalizeModel(m);m.circuitMergedRxnIDs=["AUXOSPGI";"AUXOSPFK";"AUXOSALT_merged;AUXOSALATA_Lm;AUXOSALATA_L"];
E=struct('codes',strings(3,1),'source',strings(3,1),'reference',strings(3,1));
[E,status,count]=applyEmbeddediMM904EC(E,repmat("unresolved",3,1),m, ...
    ["g6p";"f6p";"fdp";"pyr";"ala_L";"akg";"glu_L"]);
assert(count==2&&E.codes(1)=="5.3.1.9"&&E.codes(3)=="2.6.1.2"&&E.codes(2)=="" ...
    &&status(2)=="embedded_iMM904_ID_matched_chemistry_mismatch_EC_withheld", ...
    'Embedded iMM904 EC transfer or its chemistry guard failed.');
end

function T = embeddediMM904ECTable()
% EC numbers exactly as annotated by the iMM904 authors (reaction notes field
% PROTEIN_CLASS; Mo, Palsson & Herrgard 2009, BMC Syst Biol 3:37), read from
% sbml/iMM904.xml in github.com/opencobra/m_model_collection on 2026-11-13.
% Row format: iMM904 reaction ID | EC | net stoichiometry as
% metabolite@organelle:coefficient. These are 2009 annotations: some EC
% entries have since been transferred, so check them before a BRENDA query.
raw=[
"13BGH|3.2.1.58|13bdglcn@c:-1;h2o@c:-1;glc_d@c:1"
"13BGHe|3.2.1.58|13bdglcn@e:-1;h2o@e:-1;glc_d@e:1"
"13GS|2.4.1.34|udpg@c:-1;13bdglcn@c:1;h@c:1;udp@c:1"
"2DOXG6PP|3.1.3.68|2doxg6p@c:-1;h2o@c:-1;2dglc@c:1;pi@c:1"
"2HBO|1.1.1.27|2hb@c:-1;nad@c:-1;2obut@c:1;h@c:1;nadh@c:1"
"34HPPOR|1.13.11.27|34hpp@c:-1;o2@c:-1;co2@c:1;hgentis@c:1"
"3DSPHR|1.1.1.102|3dsphgn@c:-1;h@c:-1;nadph@c:-1;nadp@c:1;sphgn@c:1"
"3HAO|1.13.11.6|3hanthrn@c:-1;o2@c:-1;cmusa@c:1;h@c:1"
"3IPM3MT|2.1.1.145|3c2hmp@c:-1;amet@c:-1;3ipmmest@c:1;ahcys@c:1"
"4HGLSDm|1.5.1.12|4hglusa@m:-1;h2o@m:-1;nad@m:-1;e4hglu@m:1;h@m:2;nadh@m:1"
"AABTN|3.5.1.63|4aabutn@c:-1;h2o@c:-1;4abut@c:1;ac@c:1"
"AASAD1|1.2.1.31|l2aadp@c:-1;atp@c:-1;h@c:-1;nadph@c:-1;l2aadp6sa@c:1;amp@c:1;nadp@c:1;ppi@c:1"
"AASAD2|1.2.1.31|l2aadp@c:-1;atp@c:-1;h@c:-1;nadh@c:-1;l2aadp6sa@c:1;amp@c:1;nad@c:1;ppi@c:1"
"AATA|2.6.1.39|2oxoadp@c:-1;glu_l@c:-1;l2aadp@c:1;akg@c:1"
"ABTA|2.6.1.19|4abut@c:-1;akg@c:-1;glu_l@c:1;sucsal@c:1"
"ABUTDm|1.2.1.19|4abutn@m:-1;h2o@m:-1;nad@m:-1;4abut@m:1;h@m:2;nadh@m:1"
"ACACT1|2.3.1.9|accoa@c:-2;aacoa@c:1;coa@c:1"
"ACACT1m|2.3.1.9|accoa@m:-2;aacoa@m:1;coa@m:1"
"ACACT4p|2.3.1.16|3odcoa@x:-1;coa@x:-1;accoa@x:1;occoa@x:1"
"ACACT5p|2.3.1.16|3oddcoa@x:-1;coa@x:-1;accoa@x:1;dcacoa@x:1"
"ACACT6p|2.3.1.16|3otdcoa@x:-1;coa@x:-1;accoa@x:1;ddcacoa@x:1"
"ACACT7p|2.3.1.16|3ohdcoa@x:-1;coa@x:-1;accoa@x:1;tdcoa@x:1"
"ACACT8p|2.3.1.16|3ohodcoa@x:-1;coa@x:-1;accoa@x:1;pmtcoa@x:1"
"ACACT9p|2.3.1.16|3ohxccoa@x:-1;coa@x:-1;accoa@x:1;ttccoa@x:1"
"ACCOACr|6.4.1.2|accoa@c:-1;atp@c:-1;hco3@c:-1;adp@c:1;h@c:1;malcoa@c:1;pi@c:1"
"ACCOACrm|6.4.1.2|accoa@m:-1;atp@m:-1;hco3@m:-1;adp@m:1;h@m:1;malcoa@m:1;pi@m:1"
"ACGAM6PS|2.3.1.4|accoa@c:-1;gam6p@c:-1;acgam6p@c:1;coa@c:1;h@c:1"
"ACGAMPM|5.4.2.3|acgam6p@c:-1;acgam1p@c:1"
"ACGKm|2.7.2.8|acglu@m:-1;atp@m:-1;acg5p@m:1;adp@m:1"
"ACGSm|2.3.1.1|accoa@m:-1;glu_l@m:-1;acglu@m:1;coa@m:1;h@m:1"
"ACHBSm|4.1.3.18|2obut@m:-1;h@m:-1;pyr@m:-1;2ahbut@m:1;co2@m:1"
"ACLSm|4.1.3.18|h@m:-1;pyr@m:-2;alac_s@m:1;co2@m:1"
"ACOAHim|3.1.2.1|accoa@m:-1;h2o@m:-1;ac@m:1;coa@m:1;h@m:1"
"ACOAO4p|1.3.3.6|dcacoa@x:-1;o2@x:-1;dc2coa@x:1;h2o2@x:1"
"ACOAO5p|1.3.3.6|ddcacoa@x:-1;o2@x:-1;dd2coa@x:1;h2o2@x:1"
"ACOAO6p|1.3.3.6|o2@x:-1;tdcoa@x:-1;h2o2@x:1;td2coa@x:1"
"ACOAO7p|1.3.3.6|o2@x:-1;pmtcoa@x:-1;h2o2@x:1;hdd2coa@x:1"
"ACOAO8p|1.3.3.6|o2@x:-1;stcoa@x:-1;h2o2@x:1;od2coa@x:1"
"ACOAO9p|1.3.3.6|hexccoa@x:-1;o2@x:-1;h2o2@x:1;hxc2coa@x:1"
"ACOATA|2.3.1.38|acp@c:-1;accoa@c:-1;acacp@c:1;coa@c:1"
"ACOATAm|2.3.1.41|acp@m:-1;accoa@m:-1;acacp@m:1;coa@m:1"
"ACON3MT|2.1.1.145|acon_t@c:-1;amet@c:-1;acon5m@c:1;ahcys@c:1"
"ACONT|4.2.1.3|cit@c:-1;icit@c:1"
"ACONT3m|4.2.1.3|2mcit@m:-1;micit@m:1"
"ACONTm|4.2.1.3|cit@m:-1;icit@m:1"
"ACOTAim|2.6.1.11|acg5sa@m:-1;glu_l@m:-1;acorn@m:1;akg@m:1"
"ACP1e|3.1.3.2|fmn@e:-1;h2o@e:-1;pi@e:1;ribflv@e:1"
"ACS|6.2.1.1|ac@c:-1;atp@c:-1;coa@c:-1;accoa@c:1;amp@c:1;ppi@c:1"
"ACSm|6.2.1.1|ac@m:-1;atp@m:-1;coa@m:-1;accoa@m:1;amp@m:1;ppi@m:1"
"ACSp|6.2.1.1|ac@x:-1;atp@x:-1;coa@x:-1;accoa@x:1;amp@x:1;ppi@x:1"
"ADA|3.5.4.4|adn@c:-1;h@c:-1;h2o@c:-1;ins@c:1;nh4@c:1"
"ADD|3.5.4.2|ade@c:-1;h@c:-1;h2o@c:-1;hxan@c:1;nh4@c:1"
"ADK1|2.7.4.3|amp@c:-1;atp@c:-1;adp@c:2"
"ADK1m|2.7.4.3|amp@m:-1;atp@m:-1;adp@m:2"
"ADK4m|2.7.4.3|amp@m:-1;itp@m:-1;adp@m:1;idp@m:1"
"ADMDC|4.1.1.50|amet@c:-1;h@c:-1;ametam@c:1;co2@c:1"
"ADNCYC|4.6.1.1|atp@c:-1;camp@c:1;ppi@c:1"
"ADNK1|2.7.1.20|adn@c:-1;atp@c:-1;adp@c:1;amp@c:1;h@c:1"
"ADNUC|3.2.2.8|adn@c:-1;h2o@c:-1;ade@c:1;rib_d@c:1"
"ADPT|2.4.2.7|ade@c:-1;prpp@c:-1;amp@c:1;ppi@c:1"
"ADSK|2.7.1.25|aps@c:-1;atp@c:-1;adp@c:1;h@c:1;paps@c:1"
"ADSL1r|4.3.2.2|dcamp@c:-1;amp@c:1;fum@c:1"
"ADSL2r|4.3.2.2|25aics@c:-1;aicar@c:1;fum@c:1"
"ADSS|6.3.4.4|asp_l@c:-1;gtp@c:-1;imp@c:-1;dcamp@c:1;gdp@c:1;h@c:2;pi@c:1"
"AGAT_SC|2.3.1.51|1ag3p_sc@c:-0.01;dcacoa@c:-0.02;ddcacoa@c:-0.06;hdcoa@c:-0.17;ocdycacoa@c:-0.09;odecoa@c:-0.24;pmtcoa@c:-0.27;stcoa@c:-0.05;tdcoa@c:-0.1;coa@c:1;pa_sc@c:0.01"
"AGPRim|1.2.1.38|acg5p@m:-1;h@m:-1;nadph@m:-1;acg5sa@m:1;nadp@m:1;pi@m:1"
"AGTi|2.6.1.44|ala_l@c:-1;glx@c:-1;gly@c:1;pyr@c:1"
"AHCi|3.3.1.1|ahcys@c:-1;h2o@c:-1;adn@c:1;hcys_l@c:1"
"AHSERL|4.2.99.10|achms@c:-1;ch4s@c:-1;ac@c:1;h@c:1;met_l@c:1"
"AHSERL2|4.2.99.8|achms@c:-1;h2s@c:-1;ac@c:1;h@c:1;hcys_l@c:1"
"AICART|2.1.2.3|10fthf@c:-1;aicar@c:-1;fprica@c:1;thf@c:1"
"AIRCr|4.1.1.21|air@c:-1;co2@c:-1;5aizc@c:1;h@c:1"
"AKGDam|1.2.4.2|akg@m:-1;h@m:-1;lpam@m:-1;co2@m:1;sdhlam@m:1"
"AKGDbm|2.3.1.61|coa@m:-1;sdhlam@m:-1;dhlam@m:1;succoa@m:1"
"AKP1|3.1.3.1|ahdt@c:-1;h2o@c:-3;dhnpt@c:1;h@c:2;pi@c:3"
"ALASm|2.3.1.37|gly@m:-1;h@m:-1;succoa@m:-1;5aop@m:1;co2@m:1;coa@m:1"
"ALATA_L|2.6.1.2|akg@c:-1;ala_l@c:-1;glu_l@c:1;pyr@c:1"
"ALATA_Lm|2.6.1.2|akg@m:-1;ala_l@m:-1;glu_l@m:1;pyr@m:1"
"ALATRS|6.1.1.7|ala_l@c:-1;atp@c:-1;trnaala@c:-1;alatrna@c:1;amp@c:1;ppi@c:1"
"ALCD2if|1.1.1.1|etoh@c:-1;nad@c:-1;acald@c:1;h@c:1;nadh@c:1"
"ALCD2ir|1.1.1.1|acald@c:-1;h@c:-1;nadh@c:-1;etoh@c:1;nad@c:1"
"ALCD2irm|1.1.1.1|acald@m:-1;h@m:-1;nadh@m:-1;etoh@m:1;nad@m:1"
"ALCD2x|1.1.1.1|etoh@c:-1;nad@c:-1;acald@c:1;h@c:1;nadh@c:1"
"ALDD19x|1.2.1.39|h2o@c:-1;nad@c:-1;pacald@c:-1;h@c:2;nadh@c:1;pac@c:1"
"ALDD22x|1.2.1.39|aproa@c:-1;h2o@c:-1;nad@c:-1;ala_b@c:1;h@c:2;nadh@c:1"
"ALDD2xm|1.2.1.3|acald@m:-1;h2o@m:-1;nad@m:-1;ac@m:1;h@m:2;nadh@m:1"
"ALDD2y|1.2.1.4|acald@c:-1;h2o@c:-1;nadp@c:-1;ac@c:1;h@c:2;nadph@c:1"
"ALDD2ym|1.2.1.4|acald@m:-1;h2o@m:-1;nadp@m:-1;ac@m:1;h@m:2;nadph@m:1"
"ALLTAHr|3.5.3.4|alltt@c:-1;h2o@c:-1;urdglyc@c:1;urea@c:1"
"ALLTNr|3.5.2.5|alltn@c:-1;h2o@c:-1;alltt@c:1;h@c:1"
"ALPHNH|3.5.1.54|allphn@c:-1;h@c:-3;h2o@c:-1;co2@c:2;nh4@c:2"
"AMAOTr|2.6.1.62|8aonn@c:-1;amet@c:-1;amob@c:1;dann@c:1"
"AMID|3.5.1.4|4gudbd@c:-1;h2o@c:-1;4gudbutn@c:1;nh4@c:1"
"AMID2|3.5.1.4|h2o@c:-1;pad@c:-1;nh4@c:1;pac@c:1"
"AMID3|3.5.1.4|h2o@c:-1;iad@c:-1;ind3ac@c:1;nh4@c:1"
"AMPN|3.2.2.4|amp@c:-1;h2o@c:-1;ade@c:1;r5p@c:1"
"ANPRT|2.4.2.18|anth@c:-1;prpp@c:-1;ppi@c:1;pran@c:1"
"ANS|4.1.3.27|chor@c:-1;gln_l@c:-1;anth@c:1;glu_l@c:1;h@c:1;pyr@c:1"
"AP4AHr|3.6.1.41|ap4a@c:-1;h2o@c:-1;adp@c:2;h@c:2"
"APRTO2|1.4.3.4|aprut@c:-1;h2o@c:-1;o2@c:-1;h2o2@c:1;n4abutn@c:1;nh4@c:1"
"ARAB1D2|1.1.1.117|arab_d@c:-1;nadp@c:-1;dara14lac@c:1;h@c:1;nadph@c:1"
"ARABR|1.1.1.21|arab_l@c:-1;h@c:-1;nadph@c:-1;abt@c:1;nadp@c:1"
"ARGN|3.5.3.1|arg_l@c:-1;h2o@c:-1;orn@c:1;urea@c:1"
"ARGSL|4.3.2.1|argsuc@c:-1;arg_l@c:1;fum@c:1"
"ARGSSr|6.3.4.5|asp_l@c:-1;atp@c:-1;citr_l@c:-1;amp@c:1;argsuc@c:1;h@c:1;ppi@c:1"
"ARGTRS|6.1.1.19|arg_l@c:-1;atp@c:-1;trnaarg@c:-1;amp@c:1;argtrna@c:1;ppi@c:1"
"ARGTRSm|6.1.1.19|arg_l@m:-1;atp@m:-1;trnaarg@m:-1;amp@m:1;argtrna@m:1;ppi@m:1"
"ASADi|1.2.1.11|4pasp@c:-1;h@c:-1;nadph@c:-1;aspsa@c:1;nadp@c:1;pi@c:1"
"ASNN|3.5.1.1|asn_l@c:-1;h2o@c:-1;asp_l@c:1;nh4@c:1"
"ASNNe|3.5.1.1|asn_l@e:-1;h2o@e:-1;asp_l@e:1;nh4@e:1"
"ASNS1|6.3.5.4|asp_l@c:-1;atp@c:-1;gln_l@c:-1;h2o@c:-1;amp@c:1;asn_l@c:1;glu_l@c:1;h@c:1;ppi@c:1"
"ASNTRS|6.1.1.22|asn_l@c:-1;atp@c:-1;trnaasn@c:-1;amp@c:1;asntrna@c:1;ppi@c:1"
"ASNTRSm|6.1.1.22|asn_l@m:-1;atp@m:-1;trnaasn@m:-1;amp@m:1;asntrna@m:1;ppi@m:1"
"ASPCT|2.1.3.2|asp_l@c:-1;cbp@c:-1;cbasp@c:1;h@c:1;pi@c:1"
"ASPKi|2.7.2.4|asp_l@c:-1;atp@c:-1;4pasp@c:1;adp@c:1"
"ASPTA|2.6.1.1|akg@c:-1;asp_l@c:-1;glu_l@c:1;oaa@c:1"
"ASPTAm|2.6.1.1|akg@m:-1;asp_l@m:-1;glu_l@m:1;oaa@m:1"
"ASPTAp|2.6.1.1|akg@x:-1;asp_l@x:-1;glu_l@x:1;oaa@x:1"
"ASPTRS|6.1.1.12|asp_l@c:-1;atp@c:-1;trnaasp@c:-1;amp@c:1;asptrna@c:1;ppi@c:1"
"ASPTRSm|6.1.1.12|asp_l@m:-1;atp@m:-1;trnaasp@m:-1;amp@m:1;asptrna@m:1;ppi@m:1"
"ATPATF1|2.7.7.53|adp@c:-1;atp@c:-1;h@c:-1;ap4a@c:1;pi@c:1"
"ATPATF2|2.7.7.53|adp@c:-1;gtp@c:-1;h@c:-1;ap4g@c:1;pi@c:1"
"ATPATF3|2.7.7.53|gdp@c:-1;gtp@c:-1;h@c:-1;gp4g@c:1;pi@c:1"
"ATPH1|3.6.1.5|atp@c:-1;h2o@c:-2;amp@c:1;h@c:2;pi@c:2"
"ATPPRT|2.4.2.17|atp@c:-1;prpp@c:-1;ppi@c:1;prbatp@c:1"
"ATPS|3.6.3.6|atp@c:-1;h2o@c:-1;adp@c:1;pi@c:1;h@e:1"
"ATPS3m|3.6.3.14|h@c:-3;adp@m:-1;pi@m:-1;atp@m:1;h@m:2;h2o@m:1"
"ATPS3v|3.6.3.14|h@c:-3;adp@v:-1;pi@v:-1;atp@v:1;h@v:2;h2o@v:1"
"BACCL|6.3.4.15|atp@c:-1;btn@c:-1;h@c:-1;btamp@c:1;ppi@c:1"
"BPNT|3.1.3.7|h2o@c:-1;pap@c:-1;amp@c:1;pi@c:1"
"BTDD_RR|1.1.1.4|btd_rr@c:-1;nad@c:-1;actn_r@c:1;h@c:1;nadh@c:1"
"BTSr|2.8.1.6|dtbt@c:-1;s@c:-1;btn@c:1;h@c:2"
"C3STDH1|5.3.3.1|4mzym_int1@c:-1;nad@c:-1;4mzym_int2@c:1;co2@c:1;h@c:1;nadh@c:1"
"C3STDH2|5.3.3.1|nad@c:-1;zym_int1@c:-1;co2@c:1;h@c:1;nadh@c:1;zym_int2@c:1"
"CAT|1.11.1.6|h2o2@c:-2;h2o@c:2;o2@c:1"
"CATp|1.11.1.6|h2o2@x:-2;h2o@x:2;o2@x:1"
"CBPS|6.3.5.5|atp@c:-2;gln_l@c:-1;h2o@c:-1;hco3@c:-1;adp@c:2;cbp@c:1;glu_l@c:1;h@c:2;pi@c:1"
"CDPDGPm_SC|2.7.8.8|cdpdag_sc@m:-0.01;glyc3p@m:-1;cmp@m:1;h@m:1;pgp_sc@m:0.01"
"CHLPCTD|2.7.7.15|cholp@c:-1;ctp@c:-1;h@c:-1;cdpchol@c:1;ppi@c:1"
"CHLSTI|5.3.3.5|amet@c:-1;o2@c:-1;zymst@c:-1;ahcys@c:1;ergtetrol@c:1;h@c:1;h2o@c:2"
"CHOLK|2.7.1.32|atp@c:-1;chol@c:-1;adp@c:1;cholp@c:1;h@c:1"
"CHORM|5.4.99.5|chor@c:-1;pphn@c:1"
"CHORS|4.2.3.5|3psme@c:-1;chor@c:1;pi@c:1"
"CHTNDA|3.5.1.41|chitin@c:-1;h2o@c:-1;ac@c:1;chitos@c:1;h@c:1"
"CHTNS|2.4.1.16|udpacgal@c:-1;chitin@c:1;h@c:1;udp@c:1"
"CLPNSm_SC|2.7.8.5|cdpdag_sc@m:-0.01;pg_sc@m:-0.01;clpn_sc@m:0.01;cmp@m:1;h@m:1"
"CMPN|3.2.2.10|cmp@c:-1;h2o@c:-1;csn@c:1;r5p@c:1"
"CPPPGO|1.3.3.3|cpppg3@c:-1;h@c:-2;o2@c:-1;co2@c:2;h2o@c:2;pppg9@c:1"
"CSNATm|2.3.1.7|acrn@m:-1;coa@m:-1;accoa@m:1;crn@m:1"
"CSNATp|2.3.1.7|accoa@x:-1;crn@x:-1;acrn@x:1;coa@x:1"
"CSNATr|2.3.1.7|accoa@c:-1;crn@c:-1;acrn@c:1;coa@c:1"
"CSND|3.5.4.1|csn@c:-1;h@c:-1;h2o@c:-1;nh4@c:1;ura@c:1"
"CTPS1|6.3.4.2|atp@c:-1;nh4@c:-1;utp@c:-1;adp@c:1;ctp@c:1;h@c:2;pi@c:1"
"CTPS2|6.3.4.2|atp@c:-1;gln_l@c:-1;h2o@c:-1;utp@c:-1;adp@c:1;ctp@c:1;glu_l@c:1;h@c:2;pi@c:1"
"CYOOm|1.9.3.1|focytc@m:-4;h@m:-6;o2@m:-1;h@c:6;ficytc@m:4;h2o@m:2"
"CYOR_u6m|1.10.2.2|ficytc@m:-2;h@m:-1.5;q6h2@m:-1;h@c:1.5;focytc@m:2;q6@m:1"
"CYSTGL|4.4.1.1|cyst_l@c:-1;h2o@c:-1;2obut@c:1;cys_l@c:1;nh4@c:1"
"CYSTL|4.4.1.8|cyst_l@c:-1;h2o@c:-1;hcys_l@c:1;nh4@c:1;pyr@c:1"
"CYSTLp|4.4.1.8|cyst_l@x:-1;h2o@x:-1;hcys_l@x:1;nh4@x:1;pyr@x:1"
"CYSTRS|6.1.1.16|atp@c:-1;cys_l@c:-1;trnacys@c:-1;amp@c:1;cystrna@c:1;ppi@c:1"
"CYSTS|4.2.1.22|hcys_l@c:-1;ser_l@c:-1;cyst_l@c:1;h2o@c:1"
"CYTD|3.5.4.5|cytd@c:-1;h@c:-1;h2o@c:-1;nh4@c:1;uri@c:1"
"CYTK1|2.7.4.14|atp@c:-1;cmp@c:-1;adp@c:1;cdp@c:1"
"CYTK2|2.7.4.14|atp@c:-1;dcmp@c:-1;adp@c:1;dcdp@c:1"
"DADK|2.7.4.11|atp@c:-1;damp@c:-1;adp@c:1;dadp@c:1"
"DAGCPT_SC|2.7.8.2|12dgr_sc@c:-0.01;cdpchol@c:-1;cmp@c:1;h@c:1;pc_sc@c:0.01"
"DAGPYP_SC|3.1.3.4|h2o@c:-1;pa_sc@c:-0.01;12dgr_sc@c:0.01;pi@c:1"
"DASYN_SC|2.7.7.41|ctp@c:-1;h@c:-1;pa_sc@c:-0.01;cdpdag_sc@c:0.01;ppi@c:1"
"DASYNm_SC|2.7.7.41|ctp@m:-1;h@m:-1;pa_sc@m:-0.01;cdpdag_sc@m:0.01;ppi@m:1"
"DBTSr|6.3.3.3|atp@c:-1;co2@c:-1;dann@c:-1;adp@c:1;dtbt@c:1;h@c:3;pi@c:1"
"DCMPDA|3.5.4.12|dcmp@c:-1;h@c:-1;h2o@c:-1;dump@c:1;nh4@c:1"
"DCTPD|3.5.4.13|dctp@c:-1;h@c:-1;h2o@c:-1;dutp@c:1;nh4@c:1"
"DCYTD|3.5.4.14|dcyt@c:-1;h@c:-1;h2o@c:-1;duri@c:1;nh4@c:1"
"DDPA|4.1.2.15|e4p@c:-1;h2o@c:-1;pep@c:-1;2dda7p@c:1;pi@c:1"
"DDPAm|4.1.2.15|e4p@m:-1;h2o@m:-1;pep@m:-1;2dda7p@m:1;pi@m:1"
"DESAT14|1.14.19.1|h@c:-1;nadph@c:-1;o2@c:-1;tdcoa@c:-1;h2o@c:2;nadp@c:1;tdecoa@c:1"
"DESAT16|1.14.19.1|h@c:-1;nadph@c:-1;o2@c:-1;pmtcoa@c:-1;h2o@c:2;hdcoa@c:1;nadp@c:1"
"DESAT18|1.14.19.1|h@c:-1;nadph@c:-1;o2@c:-1;stcoa@c:-1;h2o@c:2;nadp@c:1;odecoa@c:1"
"DESAT18_2|1.14.19.1|h@c:-1;nadph@c:-1;o2@c:-1;odecoa@c:-1;h2o@c:2;nadp@c:1;ocdycacoa@c:1"
"DGK1|2.7.4.8|atp@c:-1;dgmp@c:-1;adp@c:1;dgdp@c:1"
"DHAD1m|4.2.1.9|23dhmb@m:-1;3mob@m:1;h2o@m:1"
"DHAK|2.7.1.29|atp@c:-1;dha@c:-1;adp@c:1;dhap@c:1;h@c:1"
"DHFRi|1.5.1.3|dhf@c:-1;h@c:-1;nadph@c:-1;nadp@c:1;thf@c:1"
"DHFRim|1.5.1.3|dhf@m:-1;h@m:-1;nadph@m:-1;nadp@m:1;thf@m:1"
"DHFS|6.3.2.12|atp@c:-1;dhpt@c:-1;glu_l@c:-1;adp@c:1;dhf@c:1;h@c:1;pi@c:1"
"DHNPAm|4.1.2.25|dhnpt@m:-1;2ahhmp@m:1;gcald@m:1;h@m:1"
"DHORD4i|1.3.3.1|dhor_s@c:-1;q6@m:-1;orot@c:1;q6h2@m:1"
"DHORDi|1.3.3.1|dhor_s@c:-1;o2@c:-1;h2o2@c:1;orot@c:1"
"DHORTS|3.5.2.3|dhor_s@c:-1;h2o@c:-1;cbasp@c:1;h@c:1"
"DHPSm|2.5.1.15|2ahhmp@m:-1;4abz@m:-1;dhpt@m:1;h2o@m:1"
"DHQS|4.2.3.4|2dda7p@c:-1;3dhq@c:1;pi@c:1"
"DHQTi|4.2.1.10|3dhq@c:-1;3dhsk@c:1;h2o@c:1"
"DIAT|2.6.1.29|accoa@c:-1;sprm@c:-1;n1sprm@c:1;coa@c:1;h@c:1"
"DIPS|2.1.1.98|amet@c:-1;caphis@c:-1;ahcys@c:1;cmaphis@c:1;h@c:1"
"DMATT|2.5.1.1|dmpp@c:-1;ipdp@c:-1;grdp@c:1;ppi@c:1"
"DOLK|2.7.1.108|ctp@c:-1;dolichol@c:-1;cdp@c:1;dolp@c:1;h@c:1"
"DOLPMMer|2.4.1.109|dolmanp@r:-1;dolp@r:1;h@r:1;mannan@r:1"
"DOLPMTcer|2.4.1.83|dolp@c:-1;gdpmann@c:-1;gdp@c:1;dolmanp@r:1"
"DPCOAK|2.7.1.24|atp@c:-1;dpcoa@c:-1;adp@c:1;coa@c:1;h@c:1"
"DPCOAKm|2.7.1.24|atp@m:-1;dpcoa@m:-1;adp@m:1;coa@m:1;h@m:1"
"DPMVD|4.1.1.33|5dpmev@c:-1;atp@c:-1;adp@c:1;co2@c:1;ipdp@c:1;pi@c:1"
"DPR|1.1.1.169|2dhp@c:-1;h@c:-1;nadph@c:-1;nadp@c:1;pant_r@c:1"
"DPRm|1.1.1.169|2dhp@m:-1;h@m:-1;nadph@m:-1;nadp@m:1;pant_r@m:1"
"DRBK|2.7.1.15|atp@c:-1;drib@c:-1;2dr5p@c:1;adp@c:1;h@c:1"
"DTMPK|2.7.4.9|atp@c:-1;dtmp@c:-1;adp@c:1;dtdp@c:1"
"DUTPDP|3.6.1.23|dutp@c:-1;h2o@c:-1;dump@c:1;h@c:1;ppi@c:1"
"DXHPScm|1.5.99.6|h2o@c:-1;spmd@c:-1;q6@m:-1;13dampp@c:1;4abutn@c:1;q6h2@m:1"
"D_LACDcm|1.1.2.4|lac_d@c:-1;ficytc@m:-2;pyr@c:1;focytc@m:2"
"D_LACDm|1.1.2.4|ficytc@m:-2;lac_d@m:-1;focytc@m:2;pyr@m:1"
"ECOAH11p|4.2.1.17|h2o@x:-1;hxc2coa@x:-1;3hxccoa@x:1"
"ECOAH4p|4.2.1.17|3hdcoa@x:-1;dc2coa@x:1;h2o@x:1"
"ECOAH5p|4.2.1.17|3hddcoa@x:-1;dd2coa@x:1;h2o@x:1"
"ECOAH6p|4.2.1.17|3htdcoa@x:-1;h2o@x:1;td2coa@x:1"
"ECOAH7p|4.2.1.17|3hhdcoa@x:-1;h2o@x:1;hdd2coa@x:1"
"ECOAH8p|4.2.1.17|3hodcoa@x:-1;h2o@x:1;od2coa@x:1"
"EHGLAT|2.6.1.1|akg@c:-1;e4hglu@c:-1;4h2oglt@c:1;glu_l@c:1"
"EHGLATm|2.6.1.1|akg@m:-1;e4hglu@m:-1;4h2oglt@m:1;glu_l@m:1"
"EHGLATp|2.6.1.1|akg@x:-1;e4hglu@x:-1;4h2oglt@x:1;glu_l@x:1"
"ENO|4.2.1.11|2pg@c:-1;h2o@c:1;pep@c:1"
"EPMGH|3.2.1.22|epm@c:-1;h2o@c:-1;gal@c:1;man@c:1"
"ETHAK|2.7.1.82|atp@c:-1;etha@c:-1;adp@c:1;ethamp@c:1;h@c:1"
"ETHAPT_SC|2.7.8.1|12dgr_sc@c:-0.01;cdpea@c:-1;cmp@c:1;h@c:1;pe_sc@c:0.01"
"FA120ACPH|3.1.2.14|ddcaacp@c:-1;h2o@c:-1;acp@c:1;ddca@c:1;h@c:1"
"FA140ACPH|3.1.2.14|h2o@c:-1;myrsacp@c:-1;acp@c:1;h@c:1;ttdca@c:1"
"FA141ACPH|3.1.2.14|h2o@c:-1;tdeacp@c:-1;acp@c:1;h@c:1;ttdcea@c:1"
"FA160ACPH|3.1.2.14|h2o@c:-1;palmacp@c:-1;acp@c:1;h@c:1;hdca@c:1"
"FA161ACPH|3.1.2.14|h2o@c:-1;hdeacp@c:-1;acp@c:1;h@c:1;hdcea@c:1"
"FA180ACPH|3.1.2.14|h2o@c:-1;ocdcaacp@c:-1;acp@c:1;h@c:1;ocdca@c:1"
"FA181ACPH|3.1.2.14|h2o@c:-1;octeacp@c:-1;acp@c:1;h@c:1;ocdcea@c:1"
"FA182ACPH|3.1.2.14|h2o@c:-1;ocdcyaacp@c:-1;acp@c:1;h@c:1;ocdcya@c:1"
"FACOAL100p|6.2.1.3|atp@x:-1;coa@x:-1;dca@x:-1;amp@x:1;dcacoa@x:1;ppi@x:1"
"FACOAL120p|6.2.1.3|atp@x:-1;coa@x:-1;ddca@x:-1;amp@x:1;ddcacoa@x:1;ppi@x:1"
"FACOAL140|6.2.1.3|atp@c:-1;coa@c:-1;ttdca@c:-1;amp@c:1;ppi@c:1;tdcoa@c:1"
"FACOAL140p|6.2.1.3|atp@x:-1;coa@x:-1;ttdca@x:-1;amp@x:1;ppi@x:1;tdcoa@x:1"
"FACOAL141|6.2.1.3|atp@c:-1;coa@c:-1;ttdcea@c:-1;amp@c:1;ppi@c:1;tdecoa@c:1"
"FACOAL141p|6.2.1.3|atp@x:-1;coa@x:-1;ttdcea@x:-1;amp@x:1;ppi@x:1;tdecoa@x:1"
"FACOAL160|6.2.1.3|atp@c:-1;coa@c:-1;hdca@c:-1;amp@c:1;pmtcoa@c:1;ppi@c:1"
"FACOAL160p|6.2.1.3|atp@x:-1;coa@x:-1;hdca@x:-1;amp@x:1;pmtcoa@x:1;ppi@x:1"
"FACOAL161|6.2.1.3|atp@c:-1;coa@c:-1;hdcea@c:-1;amp@c:1;hdcoa@c:1;ppi@c:1"
"FACOAL161p|6.2.1.3|atp@x:-1;coa@x:-1;hdcea@x:-1;amp@x:1;hdcoa@x:1;ppi@x:1"
"FACOAL180|6.2.1.3|atp@c:-1;coa@c:-1;ocdca@c:-1;amp@c:1;ppi@c:1;stcoa@c:1"
"FACOAL181|6.2.1.3|atp@c:-1;coa@c:-1;ocdcea@c:-1;amp@c:1;odecoa@c:1;ppi@c:1"
"FACOAL182|6.2.1.3|atp@c:-1;coa@c:-1;ocdcya@c:-1;amp@c:1;ocdycacoa@c:1;ppi@c:1"
"FACOAL240p|6.2.1.3|atp@x:-1;coa@x:-1;ttc@x:-1;amp@x:1;ppi@x:1;ttccoa@x:1"
"FACOAL260p|6.2.1.3|atp@x:-1;coa@x:-1;hexc@x:-1;amp@x:1;hexccoa@x:1;ppi@x:1"
"FACOAL80p|6.2.1.2|atp@x:-1;coa@x:-1;octa@x:-1;amp@x:1;occoa@x:1;ppi@x:1"
"FALDH|1.2.1.1|fald@c:-1;gthrd@c:-1;nad@c:-1;sfglutth@c:1;h@c:1;nadh@c:1"
"FAS100|2.3.1.85|h@c:-3;malcoa@c:-1;nadph@c:-2;octa@c:-1;co2@c:1;coa@c:1;dca@c:1;h2o@c:1;nadp@c:2"
"FAS100COA|2.3.1.86|h@c:-3;malcoa@c:-1;nadph@c:-2;occoa@c:-1;co2@c:1;coa@c:1;dcacoa@c:1;h2o@c:1;nadp@c:2"
"FAS120|2.3.1.85|dca@c:-1;h@c:-3;malcoa@c:-1;nadph@c:-2;co2@c:1;coa@c:1;ddca@c:1;h2o@c:1;nadp@c:2"
"FAS120COA|2.3.1.86|dcacoa@c:-1;h@c:-3;malcoa@c:-1;nadph@c:-2;co2@c:1;coa@c:1;ddcacoa@c:1;h2o@c:1;nadp@c:2"
"FAS140|2.3.1.85|ddca@c:-1;h@c:-3;malcoa@c:-1;nadph@c:-2;co2@c:1;coa@c:1;h2o@c:1;nadp@c:2;ttdca@c:1"
"FAS140COA|2.3.1.86|ddcacoa@c:-1;h@c:-3;malcoa@c:-1;nadph@c:-2;co2@c:1;coa@c:1;h2o@c:1;nadp@c:2;tdcoa@c:1"
"FAS141|2.3.1.85|ddca@c:-1;h@c:-4;malcoa@c:-1;nadph@c:-3;o2@c:-1;co2@c:1;coa@c:1;h2o@c:3;nadp@c:3;ttdcea@c:1"
"FAS160|2.3.1.85|h@c:-3;malcoa@c:-1;nadph@c:-2;ttdca@c:-1;co2@c:1;coa@c:1;h2o@c:1;hdca@c:1;nadp@c:2"
"FAS160COA|2.3.1.86|h@c:-3;malcoa@c:-1;nadph@c:-2;tdcoa@c:-1;co2@c:1;coa@c:1;h2o@c:1;nadp@c:2;pmtcoa@c:1"
"FAS161|2.3.1.85|h@c:-4;malcoa@c:-1;nadph@c:-3;o2@c:-1;ttdca@c:-1;co2@c:1;coa@c:1;h2o@c:3;hdcea@c:1;nadp@c:3"
"FAS180|2.3.1.85|h@c:-3;hdca@c:-1;malcoa@c:-1;nadph@c:-2;co2@c:1;coa@c:1;h2o@c:1;nadp@c:2;ocdca@c:1"
"FAS180COA|2.3.1.86|h@c:-3;malcoa@c:-1;nadph@c:-2;pmtcoa@c:-1;co2@c:1;coa@c:1;h2o@c:1;nadp@c:2;stcoa@c:1"
"FAS181|2.3.1.85|h@c:-4;hdca@c:-1;malcoa@c:-1;nadph@c:-3;o2@c:-1;co2@c:1;coa@c:1;h2o@c:3;nadp@c:3;ocdcea@c:1"
"FAS240_L|2.3.1.85|h@c:-9;malcoa@c:-3;nadph@c:-6;ocdca@c:-1;co2@c:3;coa@c:3;h2o@c:3;nadp@c:6;ttc@c:1"
"FAS260|2.3.1.85|h@c:-3;malcoa@c:-1;nadph@c:-2;ttc@c:-1;co2@c:1;coa@c:1;h2o@c:1;hexc@c:1;nadp@c:2"
"FAS80COA_L|2.3.1.86|accoa@c:-1;h@c:-9;malcoa@c:-3;nadph@c:-6;co2@c:3;coa@c:3;h2o@c:3;nadp@c:6;occoa@c:1"
"FAS80_L|2.3.1.85|accoa@c:-1;h@c:-8;malcoa@c:-3;nadph@c:-6;co2@c:3;coa@c:4;h2o@c:2;nadp@c:6;octa@c:1"
"FBA|4.1.2.13|fdp@c:-1;dhap@c:1;g3p@c:1"
"FBA2|4.1.2.13|f1p@c:-1;dhap@c:1;glyald@c:1"
"FBA3|4.1.2.13|s17bp@c:-1;dhap@c:1;e4p@c:1"
"FBP|3.1.3.11|fdp@c:-1;h2o@c:-1;f6p@c:1;pi@c:1"
"FBP26|3.1.3.46|f26bp@c:-1;h2o@c:-1;f6p@c:1;pi@c:1"
"FCLTm|4.99.1.1|fe2@m:-1;ppp9@m:-1;h@m:2;pheme@m:1"
"FDH|1.2.1.2|for@c:-1;nad@c:-1;co2@c:1;nadh@c:1"
"FDNG|1.2.2.1|for@c:-1;h@c:-1;q6@m:-1;co2@c:1;q6h2@m:1"
"FKYNH|3.5.1.9|lfmkynr@c:-1;h2o@c:-1;lkynr@c:1;for@c:1;h@c:1"
"FMETTRSm|2.1.2.9|10fthf@m:-1;mettrna@m:-1;fmettrna@m:1;h@m:1;thf@m:1"
"FMNAT|2.7.7.2|atp@c:-1;fmn@c:-1;h@c:-1;fad@c:1;ppi@c:1"
"FMNATm|2.7.7.2|atp@m:-1;fmn@m:-1;h@m:-1;fad@m:1;ppi@m:1"
"FOLD3m|2.5.1.15|2ahhmd@m:-1;4abz@m:-1;dhpt@m:1;ppi@m:1"
"FRDcm|1.3.99.1|fum@c:-1;fadh2@m:-1;succ@c:1;fad@m:1"
"FRDm|1.3.99.1|fadh2@m:-1;fum@m:-1;fad@m:1;succ@m:1"
"FRTT|2.5.1.29|frdp@c:-1;ipdp@c:-1;ggdp@c:1;ppi@c:1"
"FRUK|2.7.1.56|atp@c:-1;f1p@c:-1;adp@c:1;fdp@c:1;h@c:1"
"FTHFCL|6.3.3.2|5fthf@c:-1;atp@c:-1;adp@c:1;methf@c:1;pi@c:1"
"FTHFCLm|6.3.3.2|5fthf@m:-1;atp@m:-1;adp@m:1;methf@m:1;pi@m:1"
"FTHFLi|6.3.4.3|atp@c:-1;for@c:-1;thf@c:-1;10fthf@c:1;adp@c:1;pi@c:1"
"FTHFLmi|6.3.4.3|atp@m:-1;for@m:-1;thf@m:-1;10fthf@m:1;adp@m:1;pi@m:1"
"FUM|4.2.1.2|fum@c:-1;h2o@c:-1;mal_l@c:1"
"FUMAC|3.7.1.2|4fumacac@c:-1;h2o@c:-1;acac@c:1;fum@c:1;h@c:1"
"FUMm|4.2.1.2|fum@m:-1;h2o@m:-1;mal_l@m:1"
"G12MT1g|2.4.1.131|gdpmann@g:-1;m2macchitppdol@g:-1;gdp@g:1;h@g:1;m3macchitppdol@g:1"
"G12MT2g|2.4.1.131|gdpmann@g:-1;m3macchitppdol@g:-1;gdp@g:1;h@g:1;m4macchitppdol@g:1"
"G13MTg|2.4.1.132|gdpmann@g:-1;macchitppdol@g:-1;gdp@g:1;h@g:1;m1macchitppdol@g:1"
"G3PD1ir|1.1.1.8|dhap@c:-1;h@c:-1;nadh@c:-1;glyc3p@c:1;nad@c:1"
"G3PD1irm|1.1.1.8|dhap@m:-1;h@m:-1;nadh@m:-1;glyc3p@m:1;nad@m:1"
"G3PDm|1.1.99.5|fad@m:-1;glyc3p@m:-1;dhap@m:1;fadh2@m:1"
"G5SD|1.2.1.41|glu5p@c:-1;h@c:-1;nadph@c:-1;glu5sa@c:1;nadp@c:1;pi@c:1"
"G5SD2|1.2.1.41|glu5p@c:-1;h@c:-1;nadh@c:-1;glu5sa@c:1;nad@c:1;pi@c:1"
"G6PDA|3.5.99.6|gam6p@c:-1;h2o@c:-1;f6p@c:1;nh4@c:1"
"G6PDH2|1.1.1.49|g6p@c:-1;nadp@c:-1;6pgl@c:1;h@c:1;nadph@c:1"
"G6PI|5.3.1.9|g6p@c:-1;g6p_b@c:1"
"G6PI3|5.3.1.9|g6p_b@c:-1;f6p@c:1"
"GALIGH|3.2.1.22|1dgali@c:-1;h2o@c:-1;gal@c:1;inost@c:1"
"GALK|2.7.1.6|atp@c:-1;gal@c:-1;adp@c:1;gal1p@c:1;h@c:1"
"GALS3|3.2.1.22|h2o@c:-1;melib@c:-1;gal@c:1;glc_d@c:1"
"GALT|2.7.7.10|gal1p@c:-1;h@c:-1;utp@c:-1;ppi@c:1;udpgal@c:1"
"GALU|2.7.7.9|g1p@c:-1;h@c:-1;utp@c:-1;ppi@c:1;udpg@c:1"
"GAPD|1.2.1.12|g3p@c:-1;nad@c:-1;pi@c:-1;13dpg@c:1;h@c:1;nadh@c:1"
"GARFTi|2.1.2.2|10fthf@c:-1;gar@c:-1;fgam@c:1;h@c:1;thf@c:1"
"GBEZ|2.4.1.18|14glun@c:-1;glycogen@c:1;h2o@c:1"
"GCC2am|1.4.4.2|gly@m:-1;h@m:-1;lpam@m:-1;alpam@m:1;co2@m:1"
"GCC2bim|2.1.2.10|alpam@m:-1;thf@m:-1;dhlam@m:1;mlthf@m:1;nh4@m:1"
"GCC2cm|1.8.1.4|dhlam@m:-1;nad@m:-1;h@m:1;lpam@m:1;nadh@m:1"
"GCCam|1.4.4.2|gly@m:-1;h@m:-1;lpro@m:-1;alpro@m:1;co2@m:1"
"GCCbim|2.1.2.10|alpro@m:-1;thf@m:-1;dhlpro@m:1;mlthf@m:1;nh4@m:1"
"GCCcm|1.8.1.4|dhlpro@m:-1;nad@m:-1;h@m:1;lpro@m:1;nadh@m:1"
"GDPMANNtg|2.4.1.83|gdpmann@c:-1;gmp@g:-1;gmp@c:1;gdpmann@g:1"
"GF6PTA|2.6.1.16|f6p@c:-1;gln_l@c:-1;gam6p@c:1;glu_l@c:1"
"GGLGH|3.2.1.22|ggl@c:-1;h2o@c:-1;gal@c:1;glyc@c:1"
"GHMT2r|2.1.2.1|ser_l@c:-1;thf@c:-1;gly@c:1;h2o@c:1;mlthf@c:1"
"GHMT2rm|2.1.2.1|ser_l@m:-1;thf@m:-1;gly@m:1;h2o@m:1;mlthf@m:1"
"GK1|2.7.4.8|atp@c:-1;gmp@c:-1;adp@c:1;gdp@c:1"
"GK2|2.7.4.8|datp@c:-1;gmp@c:-1;dadp@c:1;gdp@c:1"
"GLCGSD|3.2.1.3|glycogen@c:-1;h2o@c:-1;glc_d@c:1"
"GLCGSDv|3.2.1.3|glycogen@v:-1;h2o@v:-1;glc_d@v:1"
"GLCP|2.4.1.1|glycogen@c:-1;pi@c:-1;g1p@c:1"
"GLCS2|2.4.1.11|udpg@c:-1;glycogen@c:1;h@c:1;udp@c:1"
"GLNS|6.3.1.2|atp@c:-1;glu_l@c:-1;nh4@c:-1;adp@c:1;gln_l@c:1;h@c:1;pi@c:1"
"GLNTRS|6.1.1.18|atp@c:-1;gln_l@c:-1;trnagln@c:-1;amp@c:1;glntrna@c:1;ppi@c:1"
"GLPT|2.3.1.97|glp@c:-1;tdcoa@c:-1;coa@c:1;tglp@c:1"
"GLU5K|2.7.2.11|atp@c:-1;glu_l@c:-1;adp@c:1;glu5p@c:1"
"GLUCYS|6.3.2.2|atp@c:-1;cys_l@c:-1;glu_l@c:-1;adp@c:1;glucys@c:1;h@c:1;pi@c:1"
"GLUDC|4.1.1.15|glu_l@c:-1;h@c:-1;4abut@c:1;co2@c:1"
"GLUDxi|1.4.1.2|glu_l@c:-1;h2o@c:-1;nad@c:-1;akg@c:1;h@c:1;nadh@c:1;nh4@c:1"
"GLUDyi|1.4.1.4|akg@c:-1;h@c:-1;nadph@c:-1;nh4@c:-1;glu_l@c:1;h2o@c:1;nadp@c:1"
"GLUK|2.7.1.2|atp@c:-1;glc_d@c:-1;adp@c:1;g6p_b@c:1;h@c:1"
"GLUN|3.5.1.2|gln_l@c:-1;h2o@c:-1;glu_l@c:1;nh4@c:1"
"GLUPRT|2.4.2.14|gln_l@c:-1;h2o@c:-1;prpp@c:-1;glu_l@c:1;ppi@c:1;pram@c:1"
"GLUSx|1.4.1.14|akg@c:-1;gln_l@c:-1;h@c:-1;nadh@c:-1;glu_l@c:2;nad@c:1"
"GLUTRS|6.1.1.17|atp@c:-1;glu_l@c:-1;trnaglu@c:-1;amp@c:1;glutrna@c:1;ppi@c:1"
"GLUTRSm|6.1.1.17|atp@m:-1;glu_l@m:-1;trnaglu@m:-1;amp@m:1;glutrna@m:1;ppi@m:1"
"GLYATi|2.3.1.29|2aobut@c:-1;coa@c:-1;accoa@c:1;gly@c:1"
"GLYCDy|1.1.1.6|glyc@c:-1;nadp@c:-1;dha@c:1;h@c:1;nadph@c:1"
"GLYCLm|2.1.2.10|gly@m:-1;nad@m:-1;thf@m:-1;co2@m:1;mlthf@m:1;nadh@m:1;nh4@m:1"
"GLYGS|2.4.1.11|h2o@c:-1;udpg@c:-1;14glun@c:1;h@c:1;udp@c:1"
"GLYK|2.7.1.30|atp@c:-1;glyc@c:-1;adp@c:1;glyc3p@c:1;h@c:1"
"GLYOX|3.1.2.6|h2o@c:-1;lgt_s@c:-1;gthrd@c:1;h@c:1;lac_d@c:1"
"GLYOXm|3.1.2.6|h2o@m:-1;lgt_s@m:-1;gthrd@m:1;h@m:1;lac_d@m:1"
"GLYTRS|6.1.1.14|atp@c:-1;gly@c:-1;trnagly@c:-1;amp@c:1;glytrna@c:1;ppi@c:1"
"GMPS2|6.3.5.2|atp@c:-1;gln_l@c:-1;h2o@c:-1;xmp@c:-1;amp@c:1;glu_l@c:1;gmp@c:1;h@c:2;ppi@c:1"
"GND|1.1.1.44|6pgc@c:-1;nadp@c:-1;co2@c:1;nadph@c:1;ru5p_d@c:1"
"GPDDA1|3.1.4.46|g3pc@c:-1;h2o@c:-1;chol@c:1;glyc3p@c:1;h@c:1"
"GRTT|2.5.1.10|grdp@c:-1;ipdp@c:-1;frdp@c:1;ppi@c:1"
"GTHO|1.8.1.7|gthox@c:-1;h@c:-1;nadph@c:-1;gthrd@c:2;nadp@c:1"
"GTHOm|1.8.1.7|gthox@m:-1;h@m:-1;nadph@m:-1;gthrd@m:2;nadp@m:1"
"GTHP|1.11.1.9|gthrd@c:-2;h2o2@c:-1;gthox@c:1;h2o@c:2"
"GTHPm|1.11.1.9|gthrd@m:-2;h2o2@m:-1;gthox@m:1;h2o@m:2"
"GTHS|6.3.2.3|atp@c:-1;glucys@c:-1;gly@c:-1;adp@c:1;gthrd@c:1;h@c:1;pi@c:1"
"GTMLT|2.3.2.2|ala_l@c:-1;gthrd@c:-1;cgly@c:1;gluala@c:1"
"GTPCI|3.5.4.16|gtp@c:-1;h2o@c:-1;ahdt@c:1;for@c:1;h@c:1"
"GTPCII|3.5.4.25|gtp@c:-1;h2o@c:-3;25dhpp@c:1;for@c:1;h@c:2;ppi@c:1"
"GTPH1|3.6.1.5|gtp@c:-1;h2o@c:-2;gmp@c:1;h@c:2;pi@c:2"
"GUAD|3.5.4.3|gua@c:-1;h@c:-1;h2o@c:-1;nh4@c:1;xan@c:1"
"GUAPRT|2.4.2.8|gua@c:-1;prpp@c:-1;gmp@c:1;ppi@c:1"
"HACD10p|1.1.1.35|3hxccoa@x:-1;nad@x:-1;3ohxccoa@x:1;h@x:1;nadh@x:1"
"HACD4p|1.1.1.35|3odcoa@x:-1;h@x:-1;nadh@x:-1;3hdcoa@x:1;nad@x:1"
"HACD5p|1.1.1.35|3oddcoa@x:-1;h@x:-1;nadh@x:-1;3hddcoa@x:1;nad@x:1"
"HACD6p|1.1.1.35|3otdcoa@x:-1;h@x:-1;nadh@x:-1;3htdcoa@x:1;nad@x:1"
"HACD7p|1.1.1.35|3ohdcoa@x:-1;h@x:-1;nadh@x:-1;3hhdcoa@x:1;nad@x:1"
"HACD8p|1.1.1.35|3ohodcoa@x:-1;h@x:-1;nadh@x:-1;3hodcoa@x:1;nad@x:1"
"HACNHm|4.2.1.36|b124tc@m:-1;h2o@m:-1;hicit@m:1"
"HCITSm|4.1.3.21|accoa@m:-1;akg@m:-1;h2o@m:-1;coa@m:1;h@m:1;hcit@m:1"
"HCITSn|4.1.3.21|accoa@n:-1;akg@n:-1;h2o@n:-1;coa@n:1;h@n:1;hcit@n:1"
"HCO3E|4.2.1.1|co2@c:-1;h2o@c:-1;h@c:1;hco3@c:1"
"HCYSMT|2.1.1.10|amet@c:-1;hcys_l@c:-1;ahcys@c:1;h@c:1;met_l@c:1"
"HETZK|2.7.1.50|4mhetz@c:-1;atp@c:-1;4mpetz@c:1;adp@c:1;h@c:1"
"HEX1|2.7.1.1|atp@c:-1;glc_d@c:-1;adp@c:1;g6p@c:1;h@c:1"
"HEX4|2.7.1.7|atp@c:-1;man@c:-1;adp@c:1;h@c:1;man6p@c:1"
"HEX7|2.7.1.1|atp@c:-1;fru@c:-1;adp@c:1;f6p@c:1;h@c:1"
"HGNTOR|1.13.11.5|hgentis@c:-1;o2@c:-1;4mlacac@c:1;h@c:1"
"HICITDm|1.1.1.155|hicit@m:-1;nad@m:-1;h@m:1;nadh@m:1;oxag@m:1"
"HISTD|1.1.1.23|h2o@c:-1;histd@c:-1;nad@c:-2;h@c:3;his_l@c:1;nadh@c:2"
"HISTP|3.1.3.15|h2o@c:-1;hisp@c:-1;histd@c:1;pi@c:1"
"HISTRS|6.1.1.21|atp@c:-1;his_l@c:-1;trnahis@c:-1;amp@c:1;histrna@c:1;ppi@c:1"
"HISTRSm|6.1.1.21|atp@m:-1;his_l@m:-1;trnahis@m:-1;amp@m:1;histrna@m:1;ppi@m:1"
"HKYNH|3.7.1.3|h2o@c:-1;hlkynr@c:-1;3hanthrn@c:1;ala_l@c:1"
"HMBS|4.3.1.8|h2o@c:-1;ppbng@c:-4;hmbil@c:1;nh4@c:4"
"HMGCOAR|1.1.1.34|coa@c:-1;mev_r@c:-1;nadp@c:-2;h@c:2;hmgcoa@c:1;nadph@c:2"
"HMGCOAS|4.1.3.5|coa@c:-1;h@c:-1;hmgcoa@c:-1;aacoa@c:1;accoa@c:1;h2o@c:1"
"HMGCOASm|4.1.3.5|coa@m:-1;h@m:-1;hmgcoa@m:-1;aacoa@m:1;accoa@m:1;h2o@m:1"
"HMPK1|2.7.1.49|4ahmmp@c:-1;atp@c:-1;4ampm@c:1;adp@c:1;h@c:1"
"HPPKm|2.7.6.3|2ahhmp@m:-1;atp@m:-1;2ahhmd@m:1;amp@m:1;h@m:1"
"HPROa|1.5.1.2|1p3h5c@c:-1;h@c:-2;nadh@c:-1;4hpro_lt@c:1;nad@c:1"
"HPROb|1.5.1.2|1p3h5c@c:-1;h@c:-2;nadph@c:-1;4hpro_lt@c:1;nadp@c:1"
"HPROxm|1.5.1.12|4hpro_lt@m:-1;nad@m:-1;1p3h5c@m:1;h@m:2;nadh@m:1"
"HPROym|1.5.1.12|4hpro_lt@m:-1;nadp@m:-1;1p3h5c@m:1;h@m:2;nadph@m:1"
"HSDxi|1.1.1.3|aspsa@c:-1;h@c:-1;nadh@c:-1;hom_l@c:1;nad@c:1"
"HSDyi|1.1.1.3|aspsa@c:-1;h@c:-1;nadph@c:-1;hom_l@c:1;nadp@c:1"
"HSERTA|2.3.1.31|accoa@c:-1;hom_l@c:-1;achms@c:1;coa@c:1"
"HSK|2.7.1.39|atp@c:-1;hom_l@c:-1;adp@c:1;h@c:1;phom@c:1"
"HSTPT|2.6.1.9|glu_l@c:-1;imacp@c:-1;akg@c:1;hisp@c:1"
"HXPRT|2.4.2.8|hxan@c:-1;prpp@c:-1;imp@c:1;ppi@c:1"
"ICDHxm|1.1.1.41|icit@m:-1;nad@m:-1;akg@m:1;co2@m:1;nadh@m:1"
"ICDHy|1.1.1.42|icit@c:-1;nadp@c:-1;akg@c:1;co2@c:1;nadph@c:1"
"ICDHym|1.1.1.42|icit@m:-1;nadp@m:-1;akg@m:1;co2@m:1;nadph@m:1"
"ICDHyp|1.1.1.42|icit@x:-1;nadp@x:-1;akg@x:1;co2@x:1;nadph@x:1"
"ICL|4.1.3.1|icit@c:-1;glx@c:1;succ@c:1"
"IGPDH|4.2.1.19|eig3p@c:-1;h2o@c:1;imacp@c:1"
"IGPS|4.1.1.48|2cpr5p@c:-1;h@c:-1;3ig3p@c:1;co2@c:1;h2o@c:1"
"ILETA|2.6.1.42|akg@c:-1;ile_l@c:-1;3mop@c:1;glu_l@c:1"
"ILETAm|2.6.1.42|akg@m:-1;ile_l@m:-1;3mop@m:1;glu_l@m:1"
"ILETRS|6.1.1.5|atp@c:-1;ile_l@c:-1;trnaile@c:-1;amp@c:1;iletrna@c:1;ppi@c:1"
"ILETRSm|6.1.1.5|atp@m:-1;ile_l@m:-1;trnaile@m:-1;amp@m:1;iletrna@m:1;ppi@m:1"
"IMPC|3.5.4.10|h2o@c:-1;imp@c:-1;fprica@c:1"
"IMPD|1.1.1.205|h2o@c:-1;imp@c:-1;nad@c:-1;h@c:1;nadh@c:1;xmp@c:1"
"INDPYRD|4.1.1.74|h@c:-1;indpyr@c:-1;co2@c:1;id3acald@c:1"
"INSK|2.7.1.73|atp@c:-1;ins@c:-1;adp@c:1;h@c:1;imp@c:1"
"IPDDI|5.3.3.2|ipdp@c:-1;dmpp@c:1"
"IPMD|1.1.1.85|3c2hmp@c:-1;nad@c:-1;3c4mop@c:1;h@c:1;nadh@c:1"
"IPPMIa|4.2.1.33|3c2hmp@c:-1;2ippm@c:1;h2o@c:1"
"IPPMIb|4.2.1.33|2ippm@c:-1;h2o@c:-1;3c3hmp@c:1"
"IPPS|4.1.3.12|3mob@c:-1;accoa@c:-1;h2o@c:-1;3c3hmp@c:1;coa@c:1;h@c:1"
"IPPSm|4.1.3.12|3mob@m:-1;accoa@m:-1;h2o@m:-1;3c3hmp@m:1;coa@m:1;h@m:1"
"ITCOALm|6.2.1.5|atp@m:-1;coa@m:-1;itacon@m:-1;adp@m:1;itaccoa@m:1;pi@m:1"
"KARA1im|1.1.1.86|alac_s@m:-1;h@m:-1;nadph@m:-1;23dhmb@m:1;nadp@m:1"
"KARA2im|1.1.1.86|2ahbut@m:-1;h@m:-1;nadph@m:-1;23dhmp@m:1;nadp@m:1"
"KYN|3.7.1.3|lkynr@c:-1;h2o@c:-1;ala_l@c:1;anth@c:1;h@c:1"
"KYN3OX|1.14.13.9|lkynr@c:-1;h@c:-1;nadph@c:-1;o2@c:-1;h2o@c:1;hlkynr@c:1;nadp@c:1"
"LALDO3|1.1.1.21|h@c:-1;mthgxl@c:-1;nadph@c:-1;lald_l@c:1;nadp@c:1"
"LCADi|1.2.1.21|h2o@c:-1;lald_l@c:-1;nad@c:-1;h@c:2;lac_l@c:1;nadh@c:1"
"LEUTA|2.6.1.6|akg@c:-1;leu_l@c:-1;4mop@c:1;glu_l@c:1"
"LEUTAm|2.6.1.6|akg@m:-1;leu_l@m:-1;4mop@m:1;glu_l@m:1"
"LEUTRS|6.1.1.4|atp@c:-1;leu_l@c:-1;trnaleu@c:-1;amp@c:1;leutrna@c:1;ppi@c:1"
"LEUTRSm|6.1.1.4|atp@m:-1;leu_l@m:-1;trnaleu@m:-1;amp@m:1;leutrna@m:1;ppi@m:1"
"LGTHL|4.4.1.5|gthrd@c:-1;mthgxl@c:-1;lgt_s@c:1"
"LNS14DM|1.14.14.1|h@c:-2;lanost@c:-1;nadph@c:-3;o2@c:-3;44mctr@c:1;for@c:1;h2o@c:4;nadp@c:3"
"LNS14DMx|1.14.14.1|h@c:-2;lanost@c:-1;nadh@c:-3;o2@c:-3;44mctr@c:1;for@c:1;h2o@c:4;nad@c:3"
"LNSTLS|5.4.99.7|ssq23epx@c:-1;lanost@c:1"
"LPCAT_SC|2.3.1.23|1agpc_sc@c:-0.01;dcacoa@c:-0.02;ddcacoa@c:-0.06;hdcoa@c:-0.17;ocdycacoa@c:-0.09;odecoa@c:-0.24;pmtcoa@c:-0.27;stcoa@c:-0.05;tdcoa@c:-0.1;coa@c:1;pc_sc@c:0.01"
"LYSTRS|6.1.1.6|atp@c:-1;lys_l@c:-1;trnalys@c:-1;amp@c:1;lystrna@c:1;ppi@c:1"
"LYSTRSm|6.1.1.6|atp@m:-1;lys_l@m:-1;trnalys@m:-1;amp@m:1;lystrna@m:1;ppi@m:1"
"L_LACD2cm|1.1.2.3|lac_l@c:-1;ficytc@m:-2;pyr@c:1;focytc@m:2"
"MACACI|5.2.1.2|4mlacac@c:-1;4fumacac@c:1"
"MALS|4.1.3.2|accoa@c:-1;glx@c:-1;h2o@c:-1;coa@c:1;h@c:1;mal_l@c:1"
"MALSp|4.1.3.2|accoa@x:-1;glx@x:-1;h2o@x:-1;coa@x:1;h@x:1;mal_l@x:1"
"MALT|3.2.1.20|h2o@c:-1;malt@c:-1;glc_d@c:2"
"MAN1PT|2.7.7.13|gtp@c:-1;h@c:-1;man1p@c:-1;gdpmann@c:1;ppi@c:1"
"MAN6PI|5.3.1.8|man6p@c:-1;f6p@c:1"
"MCITDm|4.2.1.79|hcit@m:-1;b124tc@m:1;h2o@m:1"
"MCITL2m|4.1.3.30|micit@m:-1;pyr@m:1;succ@m:1"
"MCITSm|4.1.3.31|h2o@m:-1;oaa@m:-1;ppcoa@m:-1;2mcit@m:1;coa@m:1;h@m:1"
"MCOATA|2.3.1.39|acp@c:-1;malcoa@c:-1;coa@c:1;malacp@c:1"
"MCOATAm|2.3.1.39|acp@m:-1;malcoa@m:-1;coa@m:1;malacp@m:1"
"MDH|1.1.1.37|mal_l@c:-1;nad@c:-1;h@c:1;nadh@c:1;oaa@c:1"
"MDHm|1.1.1.37|mal_l@m:-1;nad@m:-1;h@m:1;nadh@m:1;oaa@m:1"
"MDHp|1.1.1.37|mal_l@x:-1;nad@x:-1;h@x:1;nadh@x:1;oaa@x:1"
"ME1m|1.1.1.38|mal_l@m:-1;nad@m:-1;co2@m:1;nadh@m:1;pyr@m:1"
"ME2m|1.1.1.40|mal_l@m:-1;nadp@m:-1;co2@m:1;nadph@m:1;pyr@m:1"
"MELTGH|3.2.1.22|h2o@c:-1;melt@c:-1;gal@c:1;sbt_d@c:1"
"METAT|2.5.1.6|atp@c:-1;h2o@c:-1;met_l@c:-1;amet@c:1;pi@c:1;ppi@c:1"
"METS|2.1.1.13|5mthf@c:-1;hcys_l@c:-1;h@c:1;met_l@c:1;thf@c:1"
"METTRS|6.1.1.10|atp@c:-1;met_l@c:-1;trnamet@c:-1;amp@c:1;mettrna@c:1;ppi@c:1"
"METTRSm|6.1.1.10|atp@m:-1;met_l@m:-1;trnamet@m:-1;amp@m:1;mettrna@m:1;ppi@m:1"
"MEVK1|2.7.1.36|atp@c:-1;mev_r@c:-1;5pmev@c:1;adp@c:1;h@c:1"
"MEVK2|2.7.1.36|ctp@c:-1;mev_r@c:-1;5pmev@c:1;cdp@c:1;h@c:1"
"MEVK3|2.7.1.36|gtp@c:-1;mev_r@c:-1;5pmev@c:1;gdp@c:1;h@c:1"
"MEVK4|2.7.1.36|mev_r@c:-1;utp@c:-1;5pmev@c:1;h@c:1;udp@c:1"
"MFAPS_SC|2.1.1.16|amet@c:-1;ptdmeeta_sc@c:-0.01;ahcys@c:1;h@c:1;ptd2meeta_sc@c:0.01"
"MGSA|4.2.3.3|dhap@c:-1;mthgxl@c:1;pi@c:1"
"MHPGLUT|2.1.1.14|hcys_l@c:-1;mhpglu@c:-1;hpglu@c:1;met_l@c:1"
"MI1456PKn|2.7.1.151|atp@n:-1;mi1456p@n:-1;adp@n:1;h@n:1;mi13456p@n:1"
"MI145P6Kn|2.7.1.151|atp@n:-1;mi145p@n:-1;adp@n:1;h@n:1;mi1456p@n:1"
"MI145PKn|2.7.1.127|atp@n:-1;mi145p@n:-1;adp@n:1;h@n:1;mi1345p@n:1"
"MI1PP|3.1.3.25|h2o@c:-1;mi1p_d@c:-1;inost@c:1;pi@c:1"
"MI1PS|5.5.1.4|g6p@c:-1;mi1p_d@c:1"
"MOHMT|2.1.2.11|3mob@c:-1;h2o@c:-1;mlthf@c:-1;2dhp@c:1;thf@c:1"
"MTAP|2.4.2.28|5mta@c:-1;pi@c:-1;5mdr1p@c:1;ade@c:1"
"MTHFC|3.5.4.9|h2o@c:-1;methf@c:-1;10fthf@c:1;h@c:1"
"MTHFCm|3.5.4.9|h2o@m:-1;methf@m:-1;10fthf@m:1;h@m:1"
"MTHFD|1.5.1.5|mlthf@c:-1;nadp@c:-1;methf@c:1;nadph@c:1"
"MTHFD2i|1.5.1.5|mlthf@c:-1;nad@c:-1;methf@c:1;nadh@c:1"
"MTHFDm|1.5.1.5|mlthf@m:-1;nadp@m:-1;methf@m:1;nadph@m:1"
"MTHFR3|1.5.1.20|h@c:-2;mlthf@c:-1;nadph@c:-1;5mthf@c:1;nadp@c:1"
"MTRI|5.3.1.23|5mdr1p@c:-1;5mdru1p@c:1"
"NABTNO|1.2.1.3|h2o@c:-1;n4abutn@c:-1;nad@c:-1;4aabutn@c:1;h@c:2;nadh@c:1"
"NADDPp|3.6.1.22|h2o@x:-1;nad@x:-1;amp@x:1;h@x:2;nmn@x:1"
"NADH2_u6cm|1.6.99.3|h@c:-1;nadh@c:-1;q6@m:-1;nad@c:1;q6h2@m:1"
"NADH2_u6m|1.6.99.3|h@m:-1;nadh@m:-1;q6@m:-1;nad@m:1;q6h2@m:1"
"NADK|2.7.1.23|atp@c:-1;nad@c:-1;adp@c:1;h@c:1;nadp@c:1"
"NADKm|2.7.1.23|atp@m:-1;nad@m:-1;adp@m:1;h@m:1;nadp@m:1"
"NADN|3.2.2.5|h2o@c:-1;nad@c:-1;adprib@c:1;h@c:1;ncam@c:1"
"NADNm|3.2.2.5|h2o@m:-1;nad@m:-1;adprib@m:1;h@m:1;ncam@m:1"
"NADS1|6.3.1.5|atp@c:-1;dnad@c:-1;nh4@c:-1;amp@c:1;h@c:1;nad@c:1;ppi@c:1"
"NADS1n|6.3.1.5|atp@n:-1;dnad@n:-1;nh4@n:-1;amp@n:1;h@n:1;nad@n:1;ppi@n:1"
"NAPRT|2.4.2.11|h@c:-1;nac@c:-1;prpp@c:-1;nicrnt@c:1;ppi@c:1"
"NAPRTm|2.4.2.11|h@m:-1;nac@m:-1;prpp@m:-1;nicrnt@m:1;ppi@m:1"
"NDP1|3.6.1.6|adp@c:-1;h2o@c:-1;amp@c:1;h@c:1;pi@c:1"
"NDP3|3.6.1.6|gdp@c:-1;h2o@c:-1;gmp@c:1;h@c:1;pi@c:1"
"NDP3g|3.6.1.6|gdp@g:-1;h2o@g:-1;gmp@g:1;h@g:1;pi@g:1"
"NDP4|3.6.1.6|dgdp@c:-1;h2o@c:-1;dgmp@c:1;h@c:1;pi@c:1"
"NDP7|3.6.1.6|h2o@c:-1;udp@c:-1;h@c:1;pi@c:1;ump@c:1"
"NDPK1|2.7.4.6|atp@c:-1;gdp@c:-1;adp@c:1;gtp@c:1"
"NDPK2|2.7.4.6|atp@c:-1;udp@c:-1;adp@c:1;utp@c:1"
"NDPK3|2.7.4.6|atp@c:-1;cdp@c:-1;adp@c:1;ctp@c:1"
"NDPK4|2.7.4.6|atp@c:-1;dtdp@c:-1;adp@c:1;dttp@c:1"
"NDPK5|2.7.4.6|atp@c:-1;dgdp@c:-1;adp@c:1;dgtp@c:1"
"NDPK6|2.7.4.6|atp@c:-1;dudp@c:-1;adp@c:1;dutp@c:1"
"NDPK7|2.7.4.6|atp@c:-1;dcdp@c:-1;adp@c:1;dctp@c:1"
"NDPK8|2.7.4.6|atp@c:-1;dadp@c:-1;adp@c:1;datp@c:1"
"NDPK9|2.7.4.6|atp@c:-1;idp@c:-1;adp@c:1;itp@c:1"
"NMNAT|2.7.7.1|atp@c:-1;h@c:-1;nmn@c:-1;nad@c:1;ppi@c:1"
"NMNATm|2.7.7.1|atp@m:-1;h@m:-1;nmn@m:-1;nad@m:1;ppi@m:1"
"NNAMr|3.5.1.19|h2o@c:-1;ncam@c:-1;nac@c:1;nh4@c:1"
"NNAMrm|3.5.1.19|h2o@m:-1;ncam@m:-1;nac@m:1;nh4@m:1"
"NNAT|2.7.7.18|atp@c:-1;h@c:-1;nicrnt@c:-1;dnad@c:1;ppi@c:1"
"NNATm|2.7.7.18|atp@m:-1;h@m:-1;nicrnt@m:-1;dnad@m:1;ppi@m:1"
"NNDPR|2.4.2.19|h@c:-2;prpp@c:-1;quln@c:-1;co2@c:1;nicrnt@c:1;ppi@c:1"
"NNDPRm|2.4.2.19|h@m:-2;prpp@m:-1;quln@m:-1;co2@m:1;nicrnt@m:1;ppi@m:1"
"NNMT|2.1.1.1|amet@c:-1;ncam@c:-1;1mncam@c:1;ahcys@c:1"
"NTD1|3.1.3.5|dump@c:-1;h2o@c:-1;duri@c:1;pi@c:1"
"NTD10|3.1.3.5|h2o@c:-1;xmp@c:-1;pi@c:1;xtsn@c:1"
"NTD11|3.1.3.5|h2o@c:-1;imp@c:-1;ins@c:1;pi@c:1"
"NTD2|3.1.3.5|h2o@c:-1;ump@c:-1;pi@c:1;uri@c:1"
"NTD3|3.1.3.5|dcmp@c:-1;h2o@c:-1;dcyt@c:1;pi@c:1"
"NTD4|3.1.3.5|cmp@c:-1;h2o@c:-1;cytd@c:1;pi@c:1"
"NTD5|3.1.3.5|dtmp@c:-1;h2o@c:-1;pi@c:1;thymd@c:1"
"NTD6|3.1.3.5|damp@c:-1;h2o@c:-1;dad_2@c:1;pi@c:1"
"NTD7|3.1.3.5|amp@c:-1;h2o@c:-1;adn@c:1;pi@c:1"
"NTD8|3.1.3.5|dgmp@c:-1;h2o@c:-1;dgsn@c:1;pi@c:1"
"NTD9|3.1.3.5|gmp@c:-1;h2o@c:-1;gsn@c:1;pi@c:1"
"NTP3|3.6.1.15|gtp@c:-1;h2o@c:-1;gdp@c:1;h@c:1;pi@c:1"
"NTP4|3.6.1.15|dgtp@c:-1;h2o@c:-1;dgdp@c:1;h@c:1;pi@c:1"
"NTRLASE|3.5.5.1|h2o@c:-2;ind3acnl@c:-1;ind3ac@c:1;nh4@c:1"
"NTRLASE2|3.5.5.1|aprop@c:-1;h2o@c:-2;ala_l@c:1;nh4@c:1"
"NTRLASE3|3.5.5.1|acybut@c:-1;h2o@c:-2;glu_l@c:1;nh4@c:1"
"OBDHm|1.2.7.2|2obut@m:-1;coa@m:-1;nad@m:-1;co2@m:1;nadh@m:1;ppcoa@m:1"
"OCBTi|2.1.3.3|cbp@c:-1;orn@c:-1;citr_l@c:1;h@c:1;pi@c:1"
"OHACT1|2.3.1.84|accoa@c:-1;etoh@c:-1;aces@c:1;coa@c:1"
"OHACT2|2.3.1.84|accoa@c:-1;iamoh@c:-1;coa@c:1;iamac@c:1"
"OHACT3|2.3.1.84|2mbtoh@c:-1;accoa@c:-1;2mbac@c:1;coa@c:1"
"OHACT4|2.3.1.84|accoa@c:-1;ibutoh@c:-1;coa@c:1;ibutac@c:1"
"OHACT5|2.3.1.84|2phetoh@c:-1;accoa@c:-1;coa@c:1;pheac@c:1"
"OHPBAT|2.6.1.52|glu_l@c:-1;ohpb@c:-1;akg@c:1;phthr@c:1"
"OMPDC|4.1.1.23|h@c:-1;orot5p@c:-1;co2@c:1;ump@c:1"
"ORNDC|4.1.1.17|h@c:-1;orn@c:-1;co2@c:1;ptrc@c:1"
"ORNTA|2.6.1.13|akg@c:-1;orn@c:-1;glu_l@c:1;glu5sa@c:1"
"ORNTACim|2.3.1.35|acorn@m:-1;glu_l@m:-1;acglu@m:1;orn@m:1"
"ORPT|2.4.2.10|orot5p@c:-1;ppi@c:-1;orot@c:1;prpp@c:1"
"P5CDm|1.5.1.12|1pyr5c@m:-1;h2o@m:-2;nad@m:-1;glu_l@m:1;h@m:1;nadh@m:1"
"P5CR|1.5.1.2|1pyr5c@c:-1;h@c:-2;nadph@c:-1;nadp@c:1;pro_l@c:1"
"PANTS|6.3.2.1|ala_b@c:-1;atp@c:-1;pant_r@c:-1;amp@c:1;h@c:1;pnto_r@c:1;ppi@c:1"
"PAPSR|1.8.4.8|paps@c:-1;trdrd@c:-1;h@c:2;pap@c:1;so3@c:1;trdox@c:1"
"PC|6.4.1.1|atp@c:-1;hco3@c:-1;pyr@c:-1;adp@c:1;h@c:1;oaa@c:1;pi@c:1"
"PCDAGAT|2.3.1.43|12dgr_sc@c:-1;pc_sc@c:-1;1agpc_sc@c:1;triglyc_sc@c:1"
"PDE1|3.1.4.17|camp@c:-1;h2o@c:-1;amp@c:1;h@c:1"
"PDE2|3.1.4.17|35cdamp@c:-1;h2o@c:-1;damp@c:1;h@c:1"
"PDE3|3.1.4.17|35cimp@c:-1;h2o@c:-1;h@c:1;imp@c:1"
"PDE4|3.1.4.17|35cgmp@c:-1;h2o@c:-1;gmp@c:1;h@c:1"
"PDE5|3.1.4.17|35ccmp@c:-1;h2o@c:-1;cmp@c:1;h@c:1"
"PDHcm|1.8.1.4|dhlam@m:-1;nad@m:-1;h@m:1;lpam@m:1;nadh@m:1"
"PDHm|1.2.4.1|coa@m:-1;nad@m:-1;pyr@m:-1;accoa@m:1;co2@m:1;nadh@m:1"
"PDX5PO|1.4.3.5|o2@c:-1;pdx5p@c:-1;h2o2@c:1;pydx5p@c:1"
"PEPAT|2.3.1.88|accoa@c:-1;pepd@c:-1;apep@c:1;coa@c:1;h@c:1"
"PETHCT|2.7.7.14|ctp@c:-1;ethamp@c:-1;h@c:-1;cdpea@c:1;ppi@c:1"
"PETOHM_SC|2.1.1.17|amet@c:-1;pe_sc@c:-0.01;ahcys@c:1;h@c:1;ptdmeeta_sc@c:0.01"
"PFK|2.7.1.11|atp@c:-1;f6p@c:-1;adp@c:1;fdp@c:1;h@c:1"
"PFK26|2.7.1.105|atp@c:-1;f6p@c:-1;adp@c:1;f26bp@c:1;h@c:1"
"PFK_2|2.7.1.11|atp@c:-1;tag6p_d@c:-1;adp@c:1;h@c:1;tagdp_d@c:1"
"PFK_3|2.7.1.11|atp@c:-1;s7p@c:-1;adp@c:1;h@c:1;s17bp@c:1"
"PGAMT|5.4.2.10|gam1p@c:-1;gam6p@c:1"
"PGCD|1.1.1.95|3pg@c:-1;nad@c:-1;3php@c:1;h@c:1;nadh@c:1"
"PGI|5.3.1.9|g6p@c:-1;f6p@c:1"
"PGK|2.7.2.3|3pg@c:-1;atp@c:-1;13dpg@c:1;adp@c:1"
"PGL|3.1.1.31|6pgl@c:-1;h2o@c:-1;6pgc@c:1;h@c:1"
"PGM|5.4.2.1|2pg@c:-1;3pg@c:1"
"PGMT|5.4.2.2|g1p@c:-1;g6p@c:1"
"PGPPAm_SC|3.1.3.27|h2o@m:-1;pgp_sc@m:-0.01;pg_sc@m:0.01;pi@m:1"
"PHCDm|1.5.1.12|1p3h5c@m:-1;h2o@m:-2;nad@m:-1;e4hglu@m:1;h@m:1;nadh@m:1"
"PHETA1|2.6.1.58|akg@c:-1;phe_l@c:-1;glu_l@c:1;phpyr@c:1"
"PHETRS|6.1.1.20|atp@c:-1;phe_l@c:-1;trnaphe@c:-1;amp@c:1;phetrna@c:1;ppi@c:1"
"PHETRSm|6.1.1.20|atp@m:-1;phe_l@m:-1;trnaphe@m:-1;amp@m:1;phetrna@m:1;ppi@m:1"
"PI45BPP_SC|3.1.4.11|h2o@c:-1;ptd145bp_sc@c:-0.01;12dgr_sc@c:0.01;h@c:1;mi145p@c:1"
"PI4P5K_SC|2.7.1.68|atp@c:-1;ptd4ino_sc@c:-0.01;adp@c:1;h@c:1;ptd145bp_sc@c:0.01"
"PIN3K_SC|2.7.1.137|atp@c:-1;ptd1ino_sc@c:-0.01;adp@c:1;h@c:1;ptd3ino_sc@c:0.01"
"PIN4K_SC|2.7.1.67|atp@c:-1;ptd1ino_sc@c:-0.01;adp@c:1;h@c:1;ptd4ino_sc@c:0.01"
"PIN4Kn_SC|2.7.1.67|atp@n:-1;ptd1ino_sc@n:-0.01;adp@n:1;h@n:1;ptd4ino_sc@n:0.01"
"PINOS_SC|2.7.8.11|cdpdag_sc@c:-0.01;inost@c:-1;cmp@c:1;h@c:1;ptd1ino_sc@c:0.01"
"PMANM|5.4.2.8|man1p@c:-1;man6p@c:1"
"PMETM_SC|2.1.1.16|amet@c:-1;ptd2meeta_sc@c:-0.01;ahcys@c:1;h@c:1;pc_sc@c:0.01"
"PMEVK|2.7.4.2|5pmev@c:-1;atp@c:-1;5dpmev@c:1;adp@c:1"
"PMI12346PH|3.6.1.52|h2o@c:-1;ppmi12346p@c:-1;h@c:1;minohp@c:1;pi@c:1"
"PMPK|2.7.4.7|4ampm@c:-1;atp@c:-1;2mahmp@c:1;adp@c:1"
"PNP|2.4.2.1|pi@c:-1;rnam@c:-1;h@c:1;ncam@c:1;r1p@c:1"
"PNTK|2.7.1.33|atp@c:-1;pnto_r@c:-1;4ppan@c:1;adp@c:1;h@c:1"
"POLYAO|1.5.3.11|n1aspmd@c:-1;h2o@c:-1;o2@c:-1;aproa@c:1;aprut@c:1;h2o2@c:1"
"POLYAO2|1.5.3.11|n1sprm@c:-1;h2o@c:-1;o2@c:-1;n1aspmd@c:1;aproa@c:1;h2o2@c:1"
"POLYAO3|1.5.3.11|h2o@c:-1;o2@c:-1;sprm@c:-1;aproa@c:1;h2o2@c:1;spmd@c:1"
"PPA|3.6.1.1|h2o@c:-1;ppi@c:-1;h@c:1;pi@c:2"
"PPAm|3.6.1.1|h2o@m:-1;ppi@m:-1;h@m:1;pi@m:2"
"PPBNGS|4.2.1.24|5aop@c:-2;h@c:1;h2o@c:2;ppbng@c:1"
"PPCDC|4.1.1.36|4ppcys@c:-1;h@c:-1;co2@c:1;pan4p@c:1"
"PPCK|4.1.1.49|atp@c:-1;oaa@c:-1;adp@c:1;co2@c:1;pep@c:1"
"PPM|5.4.2.7|r1p@c:-1;r5p@c:1"
"PPNCL2|6.3.2.5|4ppan@c:-1;ctp@c:-1;cys_l@c:-1;4ppcys@c:1;cmp@c:1;h@c:1;ppi@c:1"
"PPND|1.3.1.12|nad@c:-1;pphn@c:-1;34hpp@c:1;co2@c:1;nadh@c:1"
"PPND2|1.3.1.13|nadp@c:-1;pphn@c:-1;34hpp@c:1;co2@c:1;nadph@c:1"
"PPNDH|4.2.1.51|h@c:-1;pphn@c:-1;co2@c:1;h2o@c:1;phpyr@c:1"
"PPPGOm|1.3.3.4|o2@m:-3;pppg9@m:-2;h2o@m:6;ppp9@m:2"
"PPTTm|2.5.1.33|ipdp@m:-1;pendp@m:-1;hexdp@m:1;ppi@m:1"
"PRAGSr|6.3.4.13|atp@c:-1;gly@c:-1;pram@c:-1;adp@c:1;gar@c:1;h@c:1;pi@c:1"
"PRAIS|6.3.3.1|atp@c:-1;fpram@c:-1;adp@c:1;air@c:1;h@c:2;pi@c:1"
"PRAIi|5.3.1.24|pran@c:-1;2cpr5p@c:1"
"PRAMPC|3.5.4.19|h2o@c:-1;prbamp@c:-1;prfp@c:1"
"PRASCS|6.3.2.6|5aizc@c:-1;asp_l@c:-1;atp@c:-1;25aics@c:1;adp@c:1;h@c:1;pi@c:1"
"PRATPP|3.6.1.31|h2o@c:-1;prbatp@c:-1;h@c:1;ppi@c:1;prbamp@c:1"
"PRFGS|6.3.5.3|atp@c:-1;fgam@c:-1;gln_l@c:-1;h2o@c:-1;adp@c:1;fpram@c:1;glu_l@c:1;h@c:1;pi@c:1"
"PRMICIi|5.3.1.16|prfp@c:-1;prlp@c:1"
"PROTRS|6.1.1.15|atp@c:-1;pro_l@c:-1;trnapro@c:-1;amp@c:1;ppi@c:1;protrna@c:1"
"PRPPS|2.7.6.1|atp@c:-1;r5p@c:-1;amp@c:1;h@c:1;prpp@c:1"
"PSCVTi|2.5.1.19|pep@c:-1;skm5p@c:-1;3psme@c:1;pi@c:1"
"PSERDg_SC|4.1.1.65|h@g:-1;ps_sc@g:-0.01;co2@g:1;pe_sc@g:0.01"
"PSERDm_SC|4.1.1.65|h@m:-1;ps_sc@m:-0.01;co2@m:1;pe_sc@m:0.01"
"PSERDv_SC|4.1.1.65|h@v:-1;ps_sc@v:-0.01;co2@v:1;pe_sc@v:0.01"
"PSERS_SC|2.7.8.8|cdpdag_sc@c:-0.01;ser_l@c:-1;cmp@c:1;h@c:1;ps_sc@c:0.01"
"PSERSm_SC|2.7.8.8|cdpdag_sc@m:-0.01;ser_l@m:-1;cmp@m:1;h@m:1;ps_sc@m:0.01"
"PSERT|2.6.1.52|3php@c:-1;glu_l@c:-1;akg@c:1;pser_l@c:1"
"PSP_L|3.1.3.3|h2o@c:-1;pser_l@c:-1;pi@c:1;ser_l@c:1"
"PTPATi|2.7.7.3|atp@c:-1;h@c:-1;pan4p@c:-1;dpcoa@c:1;ppi@c:1"
"PTPATim|2.7.7.3|atp@m:-1;h@m:-1;pan4p@m:-1;dpcoa@m:1;ppi@m:1"
"PUNP1|2.4.2.1|adn@c:-1;pi@c:-1;ade@c:1;r1p@c:1"
"PUNP1m|2.4.2.1|adn@m:-1;pi@m:-1;ade@m:1;r1p@m:1"
"PUNP2|2.4.2.1|dad_2@c:-1;pi@c:-1;2dr1p@c:1;ade@c:1"
"PUNP3|2.4.2.1|gsn@c:-1;pi@c:-1;gua@c:1;r1p@c:1"
"PUNP3m|2.4.2.1|gsn@m:-1;pi@m:-1;gua@m:1;r1p@m:1"
"PUNP4|2.4.2.1|dgsn@c:-1;pi@c:-1;2dr1p@c:1;gua@c:1"
"PUNP5|2.4.2.1|ins@c:-1;pi@c:-1;hxan@c:1;r1p@c:1"
"PUNP6|2.4.2.1|din@c:-1;pi@c:-1;2dr1p@c:1;hxan@c:1"
"PUNP7|2.4.2.1|pi@c:-1;xtsn@c:-1;r1p@c:1;xan@c:1"
"PYAM5PO|1.4.3.5|h2o@c:-1;o2@c:-1;pyam5p@c:-1;h2o2@c:1;nh4@c:1;pydx5p@c:1"
"PYDAMK|2.7.1.35|atp@c:-1;pydam@c:-1;adp@c:1;h@c:1;pyam5p@c:1"
"PYDXK|2.7.1.35|atp@c:-1;pydx@c:-1;adp@c:1;h@c:1;pydx5p@c:1"
"PYDXNK|2.7.1.35|atp@c:-1;pydxn@c:-1;adp@c:1;h@c:1;pdx5p@c:1"
"PYDXNO|1.4.3.5|o2@c:-1;pydxn@c:-1;h2o2@c:1;pydx@c:1"
"PYDXO|1.4.3.5|h2o@c:-2;nh4@c:-1;o2@c:-0.5;pydx@c:-1;h2o2@c:2;pydam@c:1"
"PYK|2.7.1.40|adp@c:-1;h@c:-1;pep@c:-1;atp@c:1;pyr@c:1"
"PYNP2r|2.4.2.2|pi@c:-1;uri@c:-1;r1p@c:1;ura@c:1"
"PYR5CDm|1.5.1.12|glu5sa@m:-1;h2o@m:-1;nadp@m:-1;glu_l@m:1;h@m:2;nadph@m:1"
"PYRDC|4.1.1.1|h@c:-1;pyr@c:-1;acald@c:1;co2@c:1"
"PYRDC2|4.1.1.1|acald@c:-1;h@c:-1;pyr@c:-1;actn_r@c:1;co2@c:1"
"RAFGH|3.2.1.22|h2o@c:-1;raffin@c:-1;gal@c:1;sucr@c:1"
"RBFK|2.7.1.26|atp@c:-1;ribflv@c:-1;adp@c:1;fmn@c:1;h@c:1"
"RBFKm|2.7.1.26|atp@m:-1;ribflv@m:-1;adp@m:1;fmn@m:1;h@m:1"
"RBFSa|2.5.1.9|4r5au@c:-1;db4p@c:-1;dmlz@c:1;h2o@c:2;pi@c:1"
"RBFSb|2.5.1.9|dmlz@c:-2;4r5au@c:1;ribflv@c:1"
"RBK|2.7.1.15|atp@c:-1;rib_d@c:-1;adp@c:1;h@c:1;r5p@c:1"
"RNDR1|1.17.4.1|adp@c:-1;trdrd@c:-1;dadp@c:1;h2o@c:1;trdox@c:1"
"RNDR1n|1.17.4.1|adp@n:-1;trdrd@n:-1;dadp@n:1;h2o@n:1;trdox@n:1"
"RNDR2|1.17.4.1|gdp@c:-1;trdrd@c:-1;dgdp@c:1;h2o@c:1;trdox@c:1"
"RNDR2n|1.17.4.1|gdp@n:-1;trdrd@n:-1;dgdp@n:1;h2o@n:1;trdox@n:1"
"RNDR3|1.17.4.1|cdp@c:-1;trdrd@c:-1;dcdp@c:1;h2o@c:1;trdox@c:1"
"RNDR3n|1.17.4.1|cdp@n:-1;trdrd@n:-1;dcdp@n:1;h2o@n:1;trdox@n:1"
"RNDR4|1.17.4.1|trdrd@c:-1;udp@c:-1;dudp@c:1;h2o@c:1;trdox@c:1"
"RNDR4n|1.17.4.1|trdrd@n:-1;udp@n:-1;dudp@n:1;h2o@n:1;trdox@n:1"
"RNMK|2.7.1.22|atp@c:-1;rnam@c:-1;adp@c:1;h@c:1;nmn@c:1"
"RNTR1|1.17.4.2|atp@c:-1;trdrd@c:-1;datp@c:1;h2o@c:1;trdox@c:1"
"RNTR2|1.17.4.2|gtp@c:-1;trdrd@c:-1;dgtp@c:1;h2o@c:1;trdox@c:1"
"RNTR3|1.17.4.2|ctp@c:-1;trdrd@c:-1;dctp@c:1;h2o@c:1;trdox@c:1"
"RNTR4|1.17.4.2|trdrd@c:-1;utp@c:-1;dutp@c:1;h2o@c:1;trdox@c:1"
"RPE|5.1.3.1|ru5p_d@c:-1;xu5p_d@c:1"
"RPI|5.3.1.6|r5p@c:-1;ru5p_d@c:1"
"SACCD1|1.5.1.10|l2aadp6sa@c:-1;glu_l@c:-1;h@c:-1;nadph@c:-1;h2o@c:1;nadp@c:1;saccrp_l@c:1"
"SACCD2|1.5.1.7|h2o@c:-1;nad@c:-1;saccrp_l@c:-1;akg@c:1;h@c:1;lys_l@c:1;nadh@c:1"
"SADT|2.7.7.4|atp@c:-1;h@c:-1;so4@c:-1;aps@c:1;ppi@c:1"
"SAM24MT|2.1.1.41|amet@c:-1;zymst@c:-1;ahcys@c:1;fecost@c:1;h@c:1"
"SBTD_D2|1.1.1.14|nad@c:-1;sbt_d@c:-1;fru@c:1;h@c:1;nadh@c:1"
"SBTD_L|1.1.1.14|nad@c:-1;sbt_l@c:-1;h@c:1;nadh@c:1;srb_l@c:1"
"SBTR|1.1.1.21|glc_d@c:-1;h@c:-1;nadph@c:-1;nadp@c:1;sbt_d@c:1"
"SERATi|2.3.1.30|accoa@c:-1;ser_l@c:-1;acser@c:1;coa@c:1"
"SERD_L|4.3.1.17|ser_l@c:-1;nh4@c:1;pyr@c:1"
"SERPT|2.3.1.50|h@c:-1;pmtcoa@c:-1;ser_l@c:-1;3dsphgn@c:1;co2@c:1;coa@c:1"
"SERTRS|6.1.1.11|atp@c:-1;ser_l@c:-1;trnaser@c:-1;amp@c:1;ppi@c:1;sertrna@c:1"
"SFGTH|3.1.2.12|sfglutth@c:-1;h2o@c:-1;for@c:1;gthrd@c:1;h@c:1"
"SHK3D|1.1.1.25|3dhsk@c:-1;h@c:-1;nadph@c:-1;nadp@c:1;skm@c:1"
"SHKK|2.7.1.71|atp@c:-1;skm@c:-1;adp@c:1;h@c:1;skm5p@c:1"
"SHSL1|4.2.99.9|cys_l@c:-1;suchms@c:-1;cyst_l@c:1;h@c:1;succ@c:1"
"SHSL4r|4.2.99.9|h2o@c:-1;suchms@c:-1;2obut@c:1;h@c:1;nh4@c:1;succ@c:1"
"SLFAT|2.7.7.5|adp@c:-1;h@c:-1;so4@c:-1;aps@c:1;pi@c:1"
"SPMDAT1|2.3.1.57|accoa@c:-1;spmd@c:-1;n1aspmd@c:1;coa@c:1;h@c:1"
"SPMS|2.5.1.16|ametam@c:-1;ptrc@c:-1;5mta@c:1;h@c:1;spmd@c:1"
"SPRMS|2.5.1.22|ametam@c:-1;spmd@c:-1;5mta@c:1;h@c:1;sprm@c:1"
"SQLEr|1.14.99.7|h@r:-1;nadph@r:-1;o2@r:-1;sql@r:-1;ssq23epx@r:1;h2o@r:1;nadp@r:1"
"SQLErx|1.14.99.7|h@r:-1;nadh@r:-1;o2@r:-1;sql@r:-1;ssq23epx@r:1;h2o@r:1;nad@r:1"
"SQLS|2.5.1.21|frdp@c:-2;h@c:-1;nadph@c:-1;nadp@c:1;ppi@c:2;sql@c:1"
"SSALy|1.2.1.16|h2o@c:-1;nadp@c:-1;sucsal@c:-1;h@c:2;nadph@c:1;succ@c:1"
"SUCD1m|1.3.99.1|fad@m:-1;succ@m:-1;fadh2@m:1;fum@m:1"
"SUCD2_u6m|1.3.5.1|q6@m:-1;succ@m:-1;fum@m:1;q6h2@m:1"
"SUCD3_u6m|1.3.5.1|fadh2@m:-1;q6@m:-1;fad@m:1;q6h2@m:1"
"SUCOASm|6.2.1.5|atp@m:-1;coa@m:-1;succ@m:-1;adp@m:1;pi@m:1;succoa@m:1"
"SUCRe|3.2.1.26|h2o@e:-1;sucr@e:-1;fru@e:1;glc_d@e:1"
"SULR|1.8.2.2|h2o@c:-3;h2s@c:-1;nadp@c:-3;h@c:5;nadph@c:3;so3@c:1"
"TALA|2.2.1.2|g3p@c:-1;s7p@c:-1;e4p@c:1;f6p@c:1"
"TAUDO|1.14.11.17|akg@c:-1;o2@c:-1;taur@c:-1;aacald@c:1;co2@c:1;h@c:1;so3@c:1;succ@c:1"
"THFATm|2.1.2.10|h2o@m:-1;methf@m:-1;5fthf@m:1;h@m:1"
"THFGLUS|6.3.2.17|atp@c:-1;glu_l@c:-1;thf@c:-1;adp@c:1;h@c:1;pi@c:1;thfglu@c:1"
"THMDPe|3.1.3.2|h2o@e:-2;thmpp@e:-1;h@e:1;pi@e:2;thm@e:1"
"THMPe|3.1.3.2|h2o@e:-1;thmmp@e:-1;pi@e:1;thm@e:1"
"THRA2i|4.1.2.5|athr_l@c:-1;acald@c:1;gly@c:1"
"THRAi|4.1.2.5|thr_l@c:-1;acald@c:1;gly@c:1"
"THRD_L|4.3.1.19|thr_l@c:-1;2obut@c:1;nh4@c:1"
"THRD_Lm|4.3.1.19|thr_l@m:-1;2obut@m:1;nh4@m:1"
"THRS|4.2.3.1|h2o@c:-1;phom@c:-1;pi@c:1;thr_l@c:1"
"THRTRS|6.1.1.3|atp@c:-1;thr_l@c:-1;trnathr@c:-1;amp@c:1;ppi@c:1;thrtrna@c:1"
"THRTRSm|6.1.1.3|atp@m:-1;thr_l@m:-1;trnathr@m:-1;amp@m:1;ppi@m:1;thrtrna@m:1"
"TKT1|2.2.1.1|r5p@c:-1;xu5p_d@c:-1;g3p@c:1;s7p@c:1"
"TKT2|2.2.1.1|e4p@c:-1;xu5p_d@c:-1;f6p@c:1;g3p@c:1"
"TMDK1|2.7.1.21|atp@c:-1;thymd@c:-1;adp@c:1;dtmp@c:1;h@c:1"
"TMDPK|2.7.6.2|atp@c:-1;thm@c:-1;amp@c:1;h@c:1;thmpp@c:1"
"TMDPP|2.4.2.4|pi@c:-1;thymd@c:-1;2dr1p@c:1;thym@c:1"
"TMDPPK|2.7.4.15|atp@c:-1;thmpp@c:-1;adp@c:1;thmtp@c:1"
"TMDS|2.1.1.45|dump@c:-1;mlthf@c:-1;dhf@c:1;dtmp@c:1"
"TMN|3.5.99.2|h2o@c:-1;thm@c:-1;4ahmmp@c:1;4mhetz@c:1;h@c:1"
"TMPKr|2.7.4.16|atp@c:-1;thmmp@c:-1;adp@c:1;thmpp@c:1"
"TMPPP|2.5.1.3|2mahmp@c:-1;4mpetz@c:-1;h@c:-1;ppi@c:1;thmmp@c:1"
"TPI|5.3.1.1|dhap@c:-1;g3p@c:1"
"TRDR|1.8.1.9|h@c:-1;nadph@c:-1;trdox@c:-1;nadp@c:1;trdrd@c:1"
"TRDRm|1.8.1.9|h@m:-1;nadph@m:-1;trdox@m:-1;nadp@m:1;trdrd@m:1"
"TRE6PP|3.1.3.12|h2o@c:-1;tre6p@c:-1;pi@c:1;tre@c:1"
"TRE6PS|2.4.1.15|g6p@c:-1;udpg@c:-1;h@c:1;tre6p@c:1;udp@c:1"
"TREH|3.2.1.28|h2o@c:-1;tre@c:-1;glc_d@c:2"
"TREHv|3.2.1.28|h2o@v:-1;tre@v:-1;glc_d@v:2"
"TRPO2|1.13.11.11|o2@c:-1;trp_l@c:-1;lfmkynr@c:1"
"TRPS1|4.2.1.20|3ig3p@c:-1;ser_l@c:-1;g3p@c:1;h2o@c:1;trp_l@c:1"
"TRPTRS|6.1.1.2|atp@c:-1;trnatrp@c:-1;trp_l@c:-1;amp@c:1;ppi@c:1;trptrna@c:1"
"TRPTRSm|6.1.1.2|atp@m:-1;trnatrp@m:-1;trp_l@m:-1;amp@m:1;ppi@m:1;trptrna@m:1"
"TYRNFT|2.1.2.3|10fthf@c:-1;tyr_l@c:-1;nfortyr@c:1;h@c:1;thf@c:1"
"TYRTAi|2.6.1.5|34hpp@c:-1;glu_l@c:-1;akg@c:1;tyr_l@c:1"
"TYRTAim|2.6.1.5|34hpp@m:-1;glu_l@m:-1;akg@m:1;tyr_l@m:1"
"TYRTAip|2.6.1.5|34hpp@x:-1;glu_l@x:-1;akg@x:1;tyr_l@x:1"
"TYRTRS|6.1.1.1|atp@c:-1;trnatyr@c:-1;tyr_l@c:-1;amp@c:1;ppi@c:1;tyrtrna@c:1"
"TYRTRSm|6.1.1.1|atp@m:-1;trnatyr@m:-1;tyr_l@m:-1;amp@m:1;ppi@m:1;tyrtrna@m:1"
"UDPACGLP|2.7.7.23|acgam1p@c:-1;h@c:-1;utp@c:-1;ppi@c:1;udpacgal@c:1"
"UDPG4E|5.1.3.2|udpg@c:-1;udpgal@c:1"
"UGLT|2.7.7.12|gal1p@c:-1;udpg@c:-1;g1p@c:1;udpgal@c:1"
"UGLYCHr|3.5.3.19|h@c:-2;h2o@c:-1;urdglyc@c:-1;co2@c:1;glx@c:1;nh4@c:2"
"UMPK|2.7.4.14|atp@c:-1;ump@c:-1;adp@c:1;udp@c:1"
"UPP3MT|2.1.1.107|amet@c:-2;uppg3@c:-1;ahcys@c:2;dscl@c:1;h@c:1"
"UPP3S|4.2.1.75|hmbil@c:-1;h2o@c:1;uppg3@c:1"
"UPPDC1|4.1.1.37|h@c:-4;uppg3@c:-1;co2@c:4;cpppg3@c:1"
"UPPRT|2.4.2.9|prpp@c:-1;ura@c:-1;ppi@c:1;ump@c:1"
"UREASE|6.3.4.6|atp@c:-1;hco3@c:-1;urea@c:-1;adp@c:1;allphn@c:1;h@c:1;pi@c:1"
"URIK1|2.7.1.48|atp@c:-1;uri@c:-1;adp@c:1;h@c:1;ump@c:1"
"URIK2|2.7.1.48|gtp@c:-1;uri@c:-1;gdp@c:1;h@c:1;ump@c:1"
"UTPH1|3.6.1.5|h2o@c:-2;utp@c:-1;h@c:2;pi@c:2;ump@c:1"
"VALTA|2.6.1.42|akg@c:-1;val_l@c:-1;3mob@c:1;glu_l@c:1"
"VALTAim|2.6.1.42|akg@m:-1;val_l@m:-1;3mob@m:1;glu_l@m:1"
"VALTRS|6.1.1.9|atp@c:-1;trnaval@c:-1;val_l@c:-1;amp@c:1;ppi@c:1;valtrna@c:1"
"VALTRSm|6.1.1.9|atp@m:-1;trnaval@m:-1;val_l@m:-1;amp@m:1;ppi@m:1;valtrna@m:1"
"XPPT|2.4.2.22|prpp@c:-1;xan@c:-1;ppi@c:1;xmp@c:1"
"XYLK|2.7.1.17|atp@c:-1;xylu_d@c:-1;adp@c:1;h@c:1;xu5p_d@c:1"
"XYLR|1.1.1.21|h@c:-1;nadph@c:-1;xyl_d@c:-1;nadp@c:1;xylt@c:1"
"XYLTD_D|1.1.1.9|nad@c:-1;xylt@c:-1;h@c:1;nadh@c:1;xylu_d@c:1"
"YUMPS|4.2.1.70|r5p@c:-1;ura@c:-1;h2o@c:1;psd5p@c:1"
];
cells=split(raw,"|");
T=table(cells(:,1),cells(:,2),cells(:,3), ...
    'VariableNames',{'SourceReactionID','EC','ChemicalSignature'});
end

function step2Progress(path,msg,restart)
% Console plus a small text file that survives a closed MATLAB session.
fprintf('%s\n',msg);
mode='a';if restart,mode='w';end
fid=fopen(path,mode);if fid<0,return;end
fprintf(fid,'%s\n',msg);fclose(fid);
end

function [rx,parts,control] = step2PrefixWithoutCounter(ctx,w)
% Branch module + canceller support of parent w, without its counter.
b=find(ctx.M.ModuleID==ctx.Y.ModuleID(w),1);c=find(ctx.C.ConnectedID==ctx.Y.ConnectedID(w),1);
[rx,parts]=counterCancellerSupport(c,ctx.C,ctx.D,ctx.cores,ctx.O,ctx.T,ctx.L,ctx.Q,ctx.R,ctx.P);
fast=ctx.rails(ctx.M.FastRailIndex(b));slow=ctx.rails(ctx.M.SlowRailIndex(b));
control=[fast.rx slow.rx];rx=[rx control];parts=union(parts,union(fast.parts,slow.parts));
end

function [result,alternatives,audit,how] = step2SharedBlock0(target,rx,parts,control, ...
    baseRx,baseParts,baseControl,ctx,G,cfg,requireEC,shared)
% Adding a counter only forbids more reactions/species and adds controlled
% ECs, so every Block 0 valid for the full parent is also valid without
% the counter. Search once without the counter, then keep the alternatives
% that also pass with it. Re-search the full parent only when the shared
% search found alternatives but the counter rules out all of them.
if ~cfg.step2ShareBlock0Search
    [result,alternatives,audit]=findStep2Block0(target,rx,parts,control,ctx,G,cfg,requireEC);
    how="full";return;
end
key=sprintf('%d,',[double(requireEC) target sort(baseRx) 0 sort(baseControl)]);
how="reused";
if isKey(shared,key)
    base=shared(key);
else
    [r,a,au]=findStep2Block0(target,baseRx,baseParts,baseControl,ctx,G,cfg,requireEC);
    base=struct('result',r,'alternatives',a,'audit',au);shared(key)=base;how="shared_new";
end
if ~base.result.found
    % Nothing without the counter means nothing with it (within the same
    % budgets; a budget-limited shared search is reported as such).
    result=base.result;alternatives=base.alternatives;audit=base.audit;return;
end
keep=false(1,numel(base.alternatives));
for k=1:numel(base.alternatives)
    ors=base.alternatives(k).o;ok=~any(ismember(ctx.O.rxn(ors),rx));
    for o=reshape(ors,1,[])
        if ~ok,break;end
        ok=isempty(step2Overlap(G.parts{o},parts,target,ctx.currency));
    end
    if ok&&requireEC,[~,ok]=step2AssignEC([control ctx.O.rxn(ors)],ctx,cfg);end
    keep(k)=ok;
end
if any(keep)
    alternatives=base.alternatives(keep);result=alternatives(1);audit=base.audit;return;
end
[result,alternatives,audit]=findStep2Block0(target,rx,parts,control,ctx,G,cfg,requireEC);
how="full_rerun";
end

function stats = step2CountShared(stats,how)
switch how
    case "shared_new",stats.sharedBlock0Searches=stats.sharedBlock0Searches+1;
    case "reused",stats.sharedBlock0Reused=stats.sharedBlock0Reused+1;
    case "full_rerun",stats.block0FullReruns=stats.block0FullReruns+1;
end
end

function s = joinText(values,delimiter)
% Guarantee a scalar empty string for zero-row joins (including empty tables).
values=string(values);values=values(~ismissing(values));
if isempty(values),s="";else,s=join(values(:),delimiter);end
end

% Helpers retained from the supplied A_Master (EC provenance and GPR policy).

function report(fid,msg)
    fprintf('%s\n',char(msg));fprintf(fid,'%s\n',char(msg));
end

function path=findNamed(root,name)
    path=which(name);if ~isempty(path) && isfile(path),return;end
    % A supplied full path (including an off-path temporary test file) is
    % already a file location. Do not prepend pwd or root to it.
    if isfile(name),path=char(name);return;end
    p=fullfile(pwd,name);if isfile(p),path=p;return;end
    p=fullfile(root,name);if isfile(p),path=p;return;end
    path='';
    % MATLAB recursive dir ** supported on modern versions.
    try
        d=dir(fullfile(root,'**',name));
        if ~isempty(d),path=fullfile(d(1).folder,d(1).name);end
    catch
    end
end

function [base,comp]=splitCompartments(mets)
    base=strings(size(mets));comp=strings(size(mets));
    for ii=1:numel(mets)
        x=char(mets(ii));
        t=regexp(x,'^(.*)\[([a-zA-Z])\]$','tokens','once');
        if isempty(t),t=regexp(x,'^(.*)_([ecm])$','tokens','once');end
        if ~isempty(t),base(ii)=string(t{1});comp(ii)=string(t{2});
        else,base(ii)=string(x);comp(ii)="";end
    end
end

function [gpr,glist,status]=decodeGPR(gr,rule,genes)
    gpr="";glist="";status="no_rule";
    gr=strtrim(string(gr));rule=strtrim(string(rule));
    if ~ismissing(gr) && gr~="" && lower(gr)~="nan"
        gpr=gr;status="grRules";
        % grRules may contain gene IDs rather than x(i). Keep original text.
        % Recover gene membership from exact gene tokens where possible.
        present=strings(0,1);
        for j=1:numel(genes)
            if ~isempty(regexp(char(gr),['(?<![A-Za-z0-9_.-])' regexptranslate('escape',char(genes(j))) '(?![A-Za-z0-9_.-])'],'once'))
                present(end+1,1)=genes(j); %#ok<AGROW>
            end
        end
        glist=join(present,';');
    end
    if (gpr=="" || glist=="") && rule~="" && ~ismissing(rule) && lower(rule)~="nan"
        tok=regexp(char(rule),'x\((\d+)\)','tokens');
        if ~isempty(tok)
            idx=unique(cellfun(@(v)str2double(v{1}),tok),'stable');
            if any(idx<1 | idx>numel(genes) | ~isfinite(idx))
                status="invalid_gene_index";return;
            end
            glist=join(genes(idx),';');
            gpr=rule;
            [st,en,tokens]=regexp(char(rule),'x\((\d+)\)','start','end','tokens');
            chars=char(rule);
            for q=numel(st):-1:1
                ix=str2double(tokens{q}{1});
                chars=[chars(1:st(q)-1) '(' char(genes(ix)) ')' chars(en(q)+1:end)]; %#ok<AGROW>
            end
            gpr=string(chars);status="recovered_rules";
        end
    end
end

function b=asBool(v)
    if isnumeric(v)||islogical(v),b=logical(v);return;end
    b=ismember(lower(strtrim(string(v))),["1","true","t","yes"]);
end

function out=urlencodeLocal(in)
    bytes=unicode2native(char(in),'UTF-8');
    out='';
    for n=1:numel(bytes)
        c=char(bytes(n));
        if isstrprop(c,'alphanum') || any(c=='-_.~')
            out=[out c]; %#ok<AGROW>
        else
            out=[out sprintf('%%%02X',bytes(n))]; %#ok<AGROW>
        end
    end
end

function ids=block0DefaultCurrencyIDs()
% Block 0 only: do not change the existing currency policy of other blocks.
% AMP and NH4 intentionally remain available as potential signal metabolites.
    ids=["h","h+","proton","h2o","water","atp","adp","pi","ppi", ...
        "nad","nadh","nadp","nadph","coa","co2","o2", ...
        "na1","na","k","k1","cl","cl1"];
end

function mappings=readECMappingTable(mappingPath)
% EC codes are identifiers, never numeric values or dates. Force every column
% to string BEFORE reading; TextType='string' alone still permits type inference.
    [~,~,extension]=fileparts(mappingPath);delimiter=sprintf('\t');
    if strcmpi(extension,'.csv'),delimiter=',';end
    opts=detectImportOptions(mappingPath,'FileType','text','Delimiter',delimiter, ...
        'VariableNamingRule','preserve');
    opts.VariableNamesLine=1;
    opts.DataLines=[2 Inf];
    opts=setvartype(opts,opts.VariableNames,'string');
    mappings=readtable(mappingPath,opts);
    required={'ReactionID','EC','Source','Reference'};
    assert(all(ismember(required,mappings.Properties.VariableNames)), ...
        'ECPolicy:MappingColumns', ...
        'EC evidence file needs ReactionID, EC, Source, and Reference columns.');
    for k=1:numel(required)
        mappings.(required{k})=strtrim(mappings.(required{k}));
    end
end

function [model,audit,info]=resolveRequiredECAnnotations(model,comp,ub,root,options)
% Resolve exact reaction annotations once. No enzyme-name similarity matching.
    n=numel(model.rxns);ids=string(model.rxns(:));names=ids;
    if isfield(model,'rxnNames'),names=string(model.rxnNames(:));end
    originalEC=strings(n,1);resolved=strings(n,1);source=strings(n,1);
    reference=strings(n,1);status=repmat("unresolved",n,1);
    for j=1:n
        [tokens,fields]=modelECTokens(model,j);
        if ~isempty(tokens)
            originalEC(j)=join(tokens,'; ');
            resolved(j)=originalEC(j);source(j)=join(fields,'; ');
            reference(j)="model reaction "+ids(j);status(j)="model_annotation";
        end
    end
    S=model.S;tol=1e-12;
    candidate=ub(:)>0 & full(any(S < -tol,1))' & full(any(S > tol,1))';
    candidate=candidate & ~full(any(abs(S(comp~="c",:))>tol,1))';
    modelOnly=isfield(options,'modelOnly')&&options.modelOnly;
    mappingPath='';
    if ~modelOnly
        mappingPath=which(options.mappingFile);
        if isempty(mappingPath),mappingPath=findNamed(root,options.mappingFile);end
    end
    mappings=table();
    if ~isempty(mappingPath)
        mappings=readECMappingTable(mappingPath);
    end
    brendaPath='';
    if ~modelOnly
        brendaPath=which(options.brendaFile);
        if isempty(brendaPath),brendaPath=findNamed(root,options.brendaFile);end
    end
    wanted=strings(0,1);
    if ~isempty(mappings)
        for k=1:height(mappings)
            if strcmpi(strtrim(string(mappings.Source(k))),'BRENDA')
                wanted=[wanted;completeECTokens(mappings.EC(k))]; %#ok<AGROW>
            end
        end
    end
    brendaEntries=brendaECEntries(brendaPath,unique(wanted));
    queries=0;failures=0;stopped=false;
    if ~isempty(mappings)
        for k=1:height(mappings)
            j=find(ids==strtrim(string(mappings.ReactionID(k))),1);
            if isempty(j) || resolved(j)~="",continue;end
            codes=completeECTokens(mappings.EC(k));
            evidenceSource=strtrim(string(mappings.Source(k)));
            ref=strtrim(string(mappings.Reference(k)));
            if isempty(codes) || ismissing(ref) || ref=="",continue;end
            if strcmpi(evidenceSource,'BRENDA')
                valid=false(size(codes));
                for q=1:numel(codes),valid(q)=isKey(brendaEntries,char(codes(q)));end
                codes=codes(valid);
                if isempty(codes),status(j)="BRENDA_mapping_not_confirmed_in_local_file";continue;end
                source(j)="BRENDA curated reaction mapping; local ID entry confirmed";
                reference(j)=ref+" | "+string(brendaPath);
            elseif strcmpi(evidenceSource,'BiGG')
                % A user-supplied exact reaction-ID mapping is evidence, not a
                % name guess. Require a BiGG source URL for auditability.
                if isempty(regexpi(char(ref),'^https?://(www\.)?(bigg\.bio|bigg\.ucsd\.edu|biggdata\.ucsd\.edu)/','once'))
                    status(j)="BiGG_mapping_missing_BiGG_source_URL";continue;
                end
                if ~options.online || stopped || queries>=options.maxLookups
                    status(j)="BiGG_mapping_could_not_be_verified";continue;
                end
                queries=queries+1;
                [confirmed,canonicalRef,responded]=lookupBiGGEC(ids(j),options.timeout);
                if responded,failures=0;else,failures=failures+1;end
                if failures>=options.maxFailures,stopped=true;end
                codes=intersect(codes,confirmed,'stable');
                if isempty(codes)
                    status(j)="BiGG_mapping_not_confirmed_in_database_response";continue;
                end
                source(j)="BiGG exact reaction-ID annotation; supplied mapping confirmed";
                reference(j)=canonicalRef+" | "+ref;
            else
                status(j)="mapping_source_not_BiGG_or_BRENDA";continue;
            end
            resolved(j)=join(codes,'; ');status(j)="external_sourced_mapping";
        end
    end
    for j=reshape(find(candidate & resolved==""),1,[])
        if ~options.online
            status(j)=status(j)+";online_lookup_disabled";continue;
        end
        if queries>=options.maxLookups
            status(j)="unresolved_lookup_budget";continue;
        end
        if stopped
            status(j)="unresolved_BiGG_service_unavailable";continue;
        end
        queries=queries+1;
        [codes,ref,serviceResponded]=lookupBiGGEC(ids(j),options.timeout);
        if serviceResponded,failures=0;else,failures=failures+1;end
        if ~isempty(codes)
            resolved(j)=join(codes,'; ');source(j)="BiGG exact reaction-ID annotation";
            reference(j)=ref;status(j)="BiGG_annotation";
        elseif serviceResponded
            status(j)="unresolved_BiGG_no_complete_EC";
        else
            status(j)="unresolved_BiGG_lookup_failed";
        end
        if failures>=options.maxFailures,stopped=true;end
        if mod(queries,50)==0
            fprintf('EC annotation resolution: %d BiGG lookups; %d annotated intracellular reactions.\n', ...
                queries,nnz(candidate & resolved~=""));
        end
    end
    model.circuitResolvedEC=resolved;
    model.circuitECSource=source;
    model.circuitECReference=reference;
    model.circuitECEligible=resolved~="";
    model.circuitECPolicyActive=true;
    audit=table(ids,names,originalEC,resolved,source,reference,status, ...
        model.circuitECEligible,candidate, ...
        'VariableNames',{'ReactionID','ReactionName','ModelEC','ResolvedEC', ...
        'AnnotationSource','Reference','ResolutionStatus','HasCompleteEC','CandidateNonTransport'});
    info=sprintf('EC resolution: %d model-annotated reactions; %d BiGG requests; circuit breaker=%d. Unresolved ECs exclude steps only when requireCoreEC is enabled.', ...
        nnz(originalEC~=""),queries,stopped);
end

function [codes,fields]=modelECTokens(model,j)
% Union all known EC fields, including a populated fallback after an empty one.
    codes=strings(0,1);fields=strings(0,1);
    names={'rxnECNumbers','rxnEC','eccodes','ecCodes','ecNumbers','rxnECNumber', ...
        'rxnECcodes','rxnEnzymeClassification','circuitSBMLEC'};
    for q=1:numel(names)
        if ~isfield(model,names{q}),continue;end
        value=ecAnnotationAt(model.(names{q}),j,numel(model.rxns));
        found=completeECTokens(value);
        if ~isempty(found)
            codes=[codes;found];fields(end+1,1)="model."+string(names{q}); %#ok<AGROW>
        end
    end
    if isfield(model,'rxnMiriams')
        node=ecAnnotationAt(model.rxnMiriams,j,numel(model.rxns));
        found=namedECAnnotations(node);
        if ~isempty(found)
            codes=[codes;found];fields(end+1,1)="model.rxnMiriams"; %#ok<AGROW>
        end
    end
    codes=unique(codes,'stable');fields=unique(fields,'stable');
end

function value=ecAnnotationAt(values,j,n)
    value=[];
    if iscell(values)
        if numel(values)==n,value=values{j};end
    elseif isstring(values) || isstruct(values)
        if numel(values)==n,value=values(j);end
    elseif ischar(values)
        if n==1,value=values;
        elseif size(values,1)==n,value=values(j,:);end
    end
end

function codes=completeECTokens(value)
% Only full four-field ECs qualify. Partial, blank and non-EC labels do not.
    codes=strings(0,1);
    if iscell(value)
        for k=1:numel(value),codes=[codes;completeECTokens(value{k})];end %#ok<AGROW>
    elseif isstruct(value)
        fields=fieldnames(value);
        for k=1:numel(value)
            for q=1:numel(fields)
                codes=[codes;completeECTokens(value(k).(fields{q}))]; %#ok<AGROW>
            end
        end
    elseif ischar(value) || isstring(value)
        text=string(value);text=text(~ismissing(text));
        if isempty(text),return;end
        text=char(join(text(:),' '));
        hit=regexp(text,'(?<![0-9.])[1-7]\.[0-9]+\.[0-9]+\.[0-9]+(?![0-9.\-])','match');
        codes=string(hit(:));
    end
    codes=unique(codes,'stable');
end

function codes=namedECAnnotations(node)
% Parse only EC-labelled namespaces, keys, or explicit EC-code URIs.
% A bare number elsewhere in a database payload is NOT an EC annotation.
    codes=strings(0,1);
    if iscell(node)
        for k=1:numel(node),codes=[codes;namedECAnnotations(node{k})];end %#ok<AGROW>
    elseif isstruct(node)
        names=fieldnames(node);
        for k=1:numel(node)
            % COBRA MIRIAM: name='ec-code', ids={'1.2.3.4'}.
            for nsName={'name','namespace','database','data_source'}
                if isfield(node,nsName{1})
                    ns=node(k).(nsName{1});
                    if (ischar(ns) || isstring(ns)) && isECNamespace(ns)
                        for idName={'ids','id','identifier','identifiers','accession','value'}
                            if isfield(node,idName{1})
                                codes=[codes;completeECTokens(node(k).(idName{1}))]; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
            for q=1:numel(names)
                val=node(k).(names{q});
                if isECNamespace(names{q})
                    codes=[codes;completeECTokens(val)]; %#ok<AGROW>
                else
                    codes=[codes;namedECAnnotations(val)]; %#ok<AGROW>
                end
            end
        end
    elseif ischar(node) || isstring(node)
        vals=string(node);vals=vals(~ismissing(vals));
        for k=1:numel(vals)
            text=char(vals(k));
            if ~isempty(regexpi(text,'(ec-code|ec_code|enzyme-commission|(^|\s)EC:)[/: ]','once'))
                codes=[codes;completeECTokens(text)]; %#ok<AGROW>
            end
        end
    end
    codes=unique(codes,'stable');
end

function tf=isECNamespace(value)
    text=lower(string(value));text=regexprep(text,'[^a-z0-9]','');
    tf=isscalar(text) && any(text==["ec","eccode","eccodes","ecnumber","ecnumbers", ...
        "enzymecommission","enzymecommissionnumber","enzymeclassification", ...
        "ecx0x2dcode","ec0x2dcode"]);
end

function entries=brendaECEntries(path,wanted)
% Validate explicitly sourced BRENDA mappings against ID records, not names.
    entries=containers.Map('KeyType','char','ValueType','logical');
    if isempty(path) || isempty(wanted),return;end
    fid=fopen(path,'r','n','UTF-8');
    if fid<0,return;end
    closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
    wanted=unique(wanted);
    while true
        line=fgetl(fid);if ~ischar(line),break;end
        hit=regexp(line,'^ID\s+([1-7]\.[0-9]+\.[0-9]+\.[0-9]+)(?![0-9.\-])','tokens','once');
        if ~isempty(hit) && any(wanted==string(hit{1}))
            entries(hit{1})=true;
            if entries.Count==numel(wanted),break;end
        end
    end
end

function [codes,reference,responded]=lookupBiGGEC(reactionID,timeout)
% Exact ID requests only. Prefer the current v3 object API, then legacy v2.
% Missing EC fields, errors and unknown payload shapes all fail closed.
    codes=strings(0,1);reference="";responded=false;
    id=char(string(reactionID));
    if isempty(regexp(id,'^[A-Za-z0-9_\-]+$','once')),return;end
    try
        endpoint='https://bigg.bio/api/v3/objects';
        record=webwrite(endpoint,struct('type','universalreaction','id',id), ...
            weboptions('Timeout',timeout,'MediaType','application/json','ContentType','json'));
        if isfield(record,'data') && isstruct(record.data),record=record.data;end
        if isstruct(record) && isscalar(record) && isfield(record,'bigg_id') && ...
                strcmp(string(record.bigg_id),string(reactionID))
            responded=true;codes=namedECAnnotations(record);
            if ~isempty(codes)
                reference="https://bigg.bio/universal/reactions/"+string(reactionID);return;
            end
        end
    catch ME
        if contains(ME.identifier,'HTTP404') || contains(ME.message,'404')
            responded=true; % Unknown reaction is not a network outage.
        end
        % Try the legacy exact-ID endpoint; never infer an EC from an error.
    end
    try
        endpoint=['https://bigg.ucsd.edu/api/v2/universal/reactions/',id];
        record=webread(endpoint,weboptions('Timeout',timeout,'ContentType','json'));
        if isstruct(record) && isscalar(record) && isfield(record,'bigg_id') && ...
                strcmp(string(record.bigg_id),string(reactionID))
            responded=true;codes=namedECAnnotations(record);
            if ~isempty(codes),reference=string(endpoint);end
        end
    catch ME
        if contains(ME.identifier,'HTTP404') || contains(ME.message,'404')
            responded=true;
        end
    end
end
