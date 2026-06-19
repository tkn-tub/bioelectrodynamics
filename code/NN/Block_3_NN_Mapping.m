%\[text] # Step_3_NN_Mapping.m
%\[text] Step 3 of the metabolic mapping: find 2-neuron NN implementations in iMM904
%\[text] that can be wired directly to the chemical counter's output, i.e. two
%\[text] reactions r\_Na and r\_Nb that both consume the SAME counter-output
%\[text] metabolite m\_c identified in Step\_2\_Counter\_Mapping.m's validated
%\[text] joint solutions (`joint_passed` from its Step 5).
%\[text] NN CRN topology (Fig. 1, right-hand box): r\_Ni distributes m\_c -> m\_Na, m\_Nb;
%\[text] r\_Na consumes m\_Na producing m\_Nc; r\_Nb consumes m\_Nb producing m\_Nc;
%\[text] r\_No consumes m\_Nc (from both branches) -> m\_out.
%\[text] RESTRICTED VARIANT used here: r\_Ni is taken to be the identity/no-op
%\[text] distributor, i.e. the two neuron reactions r\_Na and r\_Nb consume the
%\[text] counter output m\_c DIRECTLY (m\_Na = m\_Nb = m\_c). This is the
%\[text] "shared counter-output input" wiring you asked for.
%\[text] ALGORITHM (anchored on m\_c FIRST, as requested)
%\[text] - Step 3.0: build, once, a cache of gene lists per reaction and a
%\[text]   per-m\_c lookup of the upstream (differentiator+counter) reactions
%\[text]   and genes that must stay disjoint. This used to be recomputed
%\[text]   inside every inner loop iteration, which is what made the script
%\[text]   slow; now it is computed exactly once per distinct m\_c.
%\[text] - Step 3.1: for each distinct m\_c, find every gene-disjoint pair of
%\[text]   reactions (r\_Na, r\_Nb) that consume m\_c (excluding the upstream
%\[text]   chain reactions, which already consume m\_c via r\_Cr).
%\[text] - Step 3.2: for each pair, look for a metabolite m\_Nc that BOTH r\_Na
%\[text]   and r\_Nb produce (the convergence node, mirroring Step 1 of
%\[text]   2\_NN\_Restricted.m but starting from the shared input side).
%\[text] - Step 3.3: for each (r\_Na, r\_Nb, m\_Nc), find a gene-disjoint merge
%\[text]   reaction r\_No that consumes m\_Nc, completing the 2-neuron NN.
%\[text] ASSUMES: model, S, nRxns, hasGPR, isBadMet, nProducers in workspace
%\[text] (set by A\_Master\_File.m), and the validated joint
%\[text] differentiator+counter solutions from Step\_2\_Counter\_Mapping.m's
%\[text] STEP 5, i.e. `joint_passed` (gene-disjoint across all 7 reactions:
%\[text] r\_DF, r\_DS, r\_Dp, r\_Di, r\_Cp, r\_Cg, r\_Cr). We anchor on this set
%\[text] rather than the raw counter-only `step4`, since `joint_passed` is the
%\[text] set that has actually been confirmed compatible with a differentiator.

if ~exist('joint_passed','var') || isempty(joint_passed) %[output:group:36a7eac9]
    error(['No validated joint differentiator+counter solutions in workspace. ', ... %[output:81fe3904]
           'Run Step_2_Counter_Mapping.m first (need variable joint_passed ', ... %[output:81fe3904]
           'from its Step 5).']); %[output:81fe3904]
end %[output:group:36a7eac9]

% Local alias so the rest of this script can stay close to its previous form
step5 = joint_passed;

%\[text] #### Restrictions
MAX_PRODUCERS_MNC     = inf;  % cap on producers of the convergence metabolite m_Nc (inf = no cap)
REQUIRE_GPR_NN         = true; % enzyme-mediated reactions only
EXCLUDE_CURRENCY_NN     = true; % reuse isBadMet computed in the master file
MAX_HITS_NN             = inf; % safety cap on number of motifs collected per stage
MAX_CONSUMERS_PER_MC    = 60;  % cap on candidate r_Na/r_Nb pool size per m_c, to bound nchoosek blow-up
                                 % (a metabolite with more consumers than this is almost certainly a
                                 % currency-like hub and is reported, not silently skipped)

tStart = tic;

%\[text] ## =========================================================
%\[text] STEP 3.0: pre-compute once: (a) a gene-list cache for every reaction,
%\[text] (b) per-distinct-m_c upstream reactions/genes from joint_passed.
%\[text] This replaces the previous design, which recomputed both of these
%\[text] from scratch inside the innermost loops -- the actual cause of the
%\[text] slowdown.
%\[text] =========================================================
geneCache = cell(nRxns,1);
for r = 1:nRxns
    geneCache{r} = localGeneList(model, r);
end
geneListOf = @(rxnIdx) geneCache{rxnIdx};

m_c_list = unique([step5.m_c]);
nMc = numel(m_c_list);
fprintf('Step 3: anchoring on %d distinct counter-output metabolite(s) m_c ', nMc);
fprintf('(from %d validated joint differentiator+counter solutions)\n', numel(step5));

% Map: m_c (double) -> indices into step5 sharing that m_c
mc_to_step5idx = containers.Map('KeyType','double','ValueType','any');
for kk = 1:numel(step5)
    key = double(step5(kk).m_c);
    if isKey(mc_to_step5idx, key)
        mc_to_step5idx(key) = [mc_to_step5idx(key), kk];
    else
        mc_to_step5idx(key) = kk;
    end
end

% Map: m_c (double) -> struct('rxns', upstream reaction list, 'genes', upstream gene list)
mc_upstream = containers.Map('KeyType','double','ValueType','any');
for ii = 1:nMc
    s_mc = double(m_c_list(ii));
    idxs = mc_to_step5idx(s_mc);
    rxnsHere = [];
    for kk = idxs
        sol = step5(kk);
        rxnsHere = [rxnsHere, sol.r_DF, sol.r_DS, sol.r_Dp, sol.r_Di, sol.r_Cp, sol.r_Cg, sol.r_Cr]; %#ok<AGROW>
    end
    rxnsHere = unique(rxnsHere);
    genesHere = strings(0,1);
    for rr = rxnsHere
        genesHere = [genesHere; geneListOf(rr)]; %#ok<AGROW>
    end
    genesHere = unique(genesHere);
    mc_upstream(s_mc) = struct('rxns', rxnsHere, 'genes', genesHere);
end
fprintf('Step 3.0: gene cache (%d reactions) and per-m_c upstream lookup built in %.2f s\n', ...
        nRxns, toc(tStart));

%\[text] ## =========================================================
%\[text] STEP 3.1: for each m_c, gene-disjoint reaction pairs (r_Na, r_Nb)
%\[text] that both consume that same m_c
%\[text] =========================================================
nnStep1 = struct('m_c',{},'r_Na',{},'r_Nb',{});
n1 = 0;
mcReport = strings(0,1);

for ii = 1:nMc
    s_mc = double(m_c_list(ii));

    consumers = find(S(s_mc,:) < 0);
    if REQUIRE_GPR_NN
        consumers = consumers(hasGPR(consumers));
    end

    up = mc_upstream(s_mc);
    consumers = setdiff(consumers, up.rxns);

    if numel(consumers) > MAX_CONSUMERS_PER_MC
        mcReport(end+1) = sprintf('  m_c = %s (%s): %d candidate consumers > cap (%d) -- truncated', ...
            model.mets{s_mc}, model.metNames{s_mc}, numel(consumers), MAX_CONSUMERS_PER_MC); %#ok<AGROW>
        consumers = consumers(1:MAX_CONSUMERS_PER_MC);
    end

    if numel(consumers) < 2, continue; end

    pairs = nchoosek(consumers, 2);
    for p = 1:size(pairs,1)
        r_Na = pairs(p,1);
        r_Nb = pairs(p,2);

        genesA = geneListOf(r_Na);
        genesB = geneListOf(r_Nb);
        if isempty(genesA) || isempty(genesB), continue; end
        if isequal(genesA, genesB), continue; end
        if ~isempty(intersect(genesA, genesB)), continue; end
        % Keep the NN reactions gene-disjoint from the upstream
        % differentiator+counter chain that produced this m_c
        if ~isempty(intersect(genesA, up.genes)), continue; end
        if ~isempty(intersect(genesB, up.genes)), continue; end

        n1 = n1 + 1;
        nnStep1(n1).m_c  = s_mc;
        nnStep1(n1).r_Na = r_Na;
        nnStep1(n1).r_Nb = r_Nb;
        if n1 >= MAX_HITS_NN, break; end
    end
    if n1 >= MAX_HITS_NN, break; end
end
fprintf('Step 3.1: %d (m_c, r_Na, r_Nb) candidate pairs sharing the same counter-output input (%.2f s elapsed)\n', ...
        n1, toc(tStart));
if ~isempty(mcReport)
    fprintf('Note: some m_c had more raw consumers than MAX_CONSUMERS_PER_MC and were truncated:\n');
    fprintf('%s\n', strjoin(mcReport, '\n'));
end

if n1 == 0
    fprintf('\nNo NN candidates found at Step 3.1. Diagnostic:\n');
    for ii = 1:nMc
        s_mc = double(m_c_list(ii));
        consumers = find(S(s_mc,:) < 0);
        fprintf('  m_c = %s (%s): %d raw consumers (before GPR/exclusion filters)\n', ...
                model.mets{s_mc}, model.metNames{s_mc}, numel(consumers));
    end
    return;
end

%\[text] ## =========================================================
%\[text] STEP 3.2: for each pair, find a shared output metabolite m_Nc
%\[text] (i.e. a metabolite produced by BOTH r_Na and r_Nb)
%\[text] =========================================================
nnStep2 = struct('m_c',{},'r_Na',{},'r_Nb',{},'m_Nc',{});
n2 = 0;
badMetIdx = find(isBadMet);

for k = 1:numel(nnStep1)
    rec = nnStep1(k);

    prodA = find(S(:,rec.r_Na) > 0);
    prodB = find(S(:,rec.r_Nb) > 0);
    if EXCLUDE_CURRENCY_NN
        prodA = setdiff(prodA, badMetIdx);
        prodB = setdiff(prodB, badMetIdx);
    end

    common = intersect(prodA, prodB);
    common = setdiff(common, rec.m_c);   % avoid trivially reusing the input as output

    for c = 1:numel(common)
        s_mNc = common(c);
        if nProducers(s_mNc) > MAX_PRODUCERS_MNC, continue; end

        n2 = n2 + 1;
        nnStep2(n2).m_c  = rec.m_c;
        nnStep2(n2).r_Na = rec.r_Na;
        nnStep2(n2).r_Nb = rec.r_Nb;
        nnStep2(n2).m_Nc = s_mNc;
        if n2 >= MAX_HITS_NN, break; end
    end
    if n2 >= MAX_HITS_NN, break; end
end
fprintf('Step 3.2: %d candidates with a shared convergence metabolite m_Nc (%.2f s elapsed)\n', ...
        n2, toc(tStart));

%\[text] ## =========================================================
%\[text] STEP 3.3: find a gene-disjoint merge reaction r_No consuming m_Nc
%\[text] (completes the 2-neuron NN: r_Na, r_Nb, r_No)
%\[text] =========================================================
nnStep3 = struct('m_c',{},'r_Na',{},'r_Nb',{},'m_Nc',{},'r_No',{});
n3 = 0;

if n2 > 0
    for k = 1:numel(nnStep2)
        rec = nnStep2(k);
        up = mc_upstream(double(rec.m_c));   % cached lookup, not recomputed

        consumersNc = find(S(rec.m_Nc,:) < 0);
        if REQUIRE_GPR_NN
            consumersNc = consumersNc(hasGPR(consumersNc));
        end
        consumersNc = setdiff(consumersNc, [rec.r_Na, rec.r_Nb, up.rxns]);

        genesA = geneListOf(rec.r_Na);
        genesB = geneListOf(rec.r_Nb);
        genesAB = union(genesA, genesB);

        for u = 1:numel(consumersNc)
            r_No = consumersNc(u);
            genesO = geneListOf(r_No);
            if isempty(genesO), continue; end
            if ~isempty(intersect(genesO, genesAB)), continue; end
            if ~isempty(intersect(genesO, up.genes)), continue; end

            n3 = n3 + 1;
            nnStep3(n3).m_c  = rec.m_c;
            nnStep3(n3).r_Na = rec.r_Na;
            nnStep3(n3).r_Nb = rec.r_Nb;
            nnStep3(n3).m_Nc = rec.m_Nc;
            nnStep3(n3).r_No = r_No;
            if n3 >= MAX_HITS_NN, break; end
        end
        if n3 >= MAX_HITS_NN, break; end
    end
end
fprintf('Step 3.3: %d full NN candidates (r_Na, r_Nb, r_No) with a merge reaction found (%.2f s elapsed)\n', ...
        n3, toc(tStart));

%\[text] ## =========================================================
%\[text] Build tables for inspection
%\[text] =========================================================
if n1 > 0
    NN_Pairs = struct2table(nnStep1);
else
    NN_Pairs = table();
end

if n2 > 0
    NN_Convergent = struct2table(nnStep2);
else
    NN_Convergent = table();
end

if n3 > 0
    NN_Full = struct2table(nnStep3);
else
    NN_Full = table();
end

%\[text] #### Join NN pairs with their parent joint differentiator+counter
%\[text] solution(s) on m_c (full chain: differentiator -> counter -> NN)
joinRows = struct('diff_idx',{},'cnt_idx',{}, ...
                   'm_in',{},'m_Df',{},'m_d',{},'m_Ds',{},'m_c',{},'m_Cs',{}, ...
                   'r_DF',{},'r_DS',{},'r_Dp',{},'r_Di',{}, ...
                   'r_Cp',{},'r_Cg',{},'r_Cr',{}, ...
                   'r_Na',{},'r_Nb',{});
nJ = 0;

% Group nnStep1 hits by m_c once, instead of scanning nnStep1 fully for every step5 row
nn1ByMc = containers.Map('KeyType','double','ValueType','any');
for h = 1:numel(nnStep1)
    key = double(nnStep1(h).m_c);
    if isKey(nn1ByMc, key)
        nn1ByMc(key) = [nn1ByMc(key), h];
    else
        nn1ByMc(key) = h;
    end
end

for k = 1:numel(step5)
    sol = step5(k);
    key = double(sol.m_c);
    if ~isKey(nn1ByMc, key), continue; end
    hitIdx = nn1ByMc(key);
    for h = hitIdx
        nJ = nJ + 1;
        joinRows(nJ).diff_idx = sol.diff_idx;
        joinRows(nJ).cnt_idx  = sol.cnt_idx;
        joinRows(nJ).m_in  = sol.m_in;
        joinRows(nJ).m_Df  = sol.m_Df;
        joinRows(nJ).m_d   = sol.m_d;
        joinRows(nJ).m_Ds  = sol.m_Ds;
        joinRows(nJ).m_c   = sol.m_c;
        joinRows(nJ).m_Cs  = sol.m_Cs;
        joinRows(nJ).r_DF  = sol.r_DF;
        joinRows(nJ).r_DS  = sol.r_DS;
        joinRows(nJ).r_Dp  = sol.r_Dp;
        joinRows(nJ).r_Di  = sol.r_Di;
        joinRows(nJ).r_Cp  = sol.r_Cp;
        joinRows(nJ).r_Cg  = sol.r_Cg;
        joinRows(nJ).r_Cr  = sol.r_Cr;
        joinRows(nJ).r_Na  = nnStep1(h).r_Na;
        joinRows(nJ).r_Nb  = nnStep1(h).r_Nb;
    end
end

% if nJ > 0
%     Counter_NN_Joint = struct2table(joinRows);
%     fprintf('\n=== Differentiator + counter + NN joint table (first %d rows) ===\n', min(nJ,nJ));
%     disp(Counter_NN_Joint(1:min(nJ,nJ),:));
% else
%     Counter_NN_Joint = table();
% end

%\[text] #### Printing a few full candidates with detail (if any were found)
if n3 > 0
    fprintf('\n=== Full 2-neuron NN candidates anchored on counter output m_c ===\n');
    for i = 1:min(5, n3)
        rec = nnStep3(i);
        fprintf('\n--- NN candidate %d ---\n', i);
        fprintf('  m_c  = %s | %s\n', model.mets{rec.m_c},  model.metNames{rec.m_c});
        fprintf('  r_Na = %s | %s\n', model.rxns{rec.r_Na}, model.rxnNames{rec.r_Na});
        fprintf('  r_Nb = %s | %s\n', model.rxns{rec.r_Nb}, model.rxnNames{rec.r_Nb});
        fprintf('  m_Nc = %s | %s\n', model.mets{rec.m_Nc}, model.metNames{rec.m_Nc});
        fprintf('  r_No = %s | %s\n', model.rxns{rec.r_No}, model.rxnNames{rec.r_No});
        fprintf('  GPR(r_Na) = %s\n', model.grRules{rec.r_Na});
        fprintf('  GPR(r_Nb) = %s\n', model.grRules{rec.r_Nb});
        fprintf('  GPR(r_No) = %s\n', model.grRules{rec.r_No});
    end
elseif n1 > 0
    fprintf('\n=== NN input-sharing candidates (no convergence/merge reaction found yet) ===\n');
    for i = 1:min(5, n1)
        rec = nnStep1(i);
        fprintf('\n--- NN candidate %d ---\n', i);
        fprintf('  m_c  = %s | %s\n', model.mets{rec.m_c},  model.metNames{rec.m_c});
        fprintf('  r_Na = %s | %s\n', model.rxns{rec.r_Na}, model.rxnNames{rec.r_Na});
        fprintf('  r_Nb = %s | %s\n', model.rxns{rec.r_Nb}, model.rxnNames{rec.r_Nb});
        fprintf('  GPR(r_Na) = %s\n', model.grRules{rec.r_Na});
        fprintf('  GPR(r_Nb) = %s\n', model.grRules{rec.r_Nb});
    end
end



%\[text] #### Build and save the full diff+counter+NN join, ranked by gene complexity
%\[text] This was previously only ever pasted ad hoc into a separate print
%\[text] script and never saved here -- that's why Block_3_nn_full_motifs.mat
%\[text] kept ending up incomplete. It now lives in Block_3 itself.

if n3 > 0
    sourceNN  = nnStep3;
    hasMerge  = true;
else
    sourceNN  = nnStep1;
    hasMerge  = false;
end

mc_to_step5 = containers.Map('KeyType','double','ValueType','any');
for k = 1:numel(step5)
    key = double(step5(k).m_c);
    if isKey(mc_to_step5, key)
        mc_to_step5(key) = [mc_to_step5(key), k];
    else
        mc_to_step5(key) = k;
    end
end

fullJoinRows = struct('diff_idx',{},'cnt_idx',{}, ...
                       'm_in',{},'m_Df',{},'m_d',{},'m_Ds',{},'m_c',{},'m_Cs',{}, ...
                       'r_DF',{},'r_DS',{},'r_Dp',{},'r_Di',{}, ...
                       'r_Cp',{},'r_Cg',{},'r_Cr',{}, ...
                       'r_Na',{},'r_Nb',{},'m_Nc',{},'r_No',{},'nGenes',{});
nF = 0;
MAX_CANDIDATES_FOR_RANKING = 2000;
nBuild = min(numel(sourceNN), MAX_CANDIDATES_FOR_RANKING);

for h = 1:nBuild
    key = double(sourceNN(h).m_c);
    if ~isKey(mc_to_step5, key), continue; end
    step5idxs = mc_to_step5(key);
    sol = step5(step5idxs(1));

    nF = nF + 1;
    fullJoinRows(nF).diff_idx = sol.diff_idx;
    fullJoinRows(nF).cnt_idx  = sol.cnt_idx;
    fullJoinRows(nF).m_in  = sol.m_in;
    fullJoinRows(nF).m_Df  = sol.m_Df;
    fullJoinRows(nF).m_d   = sol.m_d;
    fullJoinRows(nF).m_Ds  = sol.m_Ds;
    fullJoinRows(nF).m_c   = sol.m_c;
    fullJoinRows(nF).m_Cs  = sol.m_Cs;
    fullJoinRows(nF).r_DF  = sol.r_DF;
    fullJoinRows(nF).r_DS  = sol.r_DS;
    fullJoinRows(nF).r_Dp  = sol.r_Dp;
    fullJoinRows(nF).r_Di  = sol.r_Di;
    fullJoinRows(nF).r_Cp  = sol.r_Cp;
    fullJoinRows(nF).r_Cg  = sol.r_Cg;
    fullJoinRows(nF).r_Cr  = sol.r_Cr;
    fullJoinRows(nF).r_Na  = sourceNN(h).r_Na;
    fullJoinRows(nF).r_Nb  = sourceNN(h).r_Nb;
    if hasMerge
        fullJoinRows(nF).m_Nc = sourceNN(h).m_Nc;
        fullJoinRows(nF).r_No = sourceNN(h).r_No;
    else
        fullJoinRows(nF).m_Nc = NaN;
        fullJoinRows(nF).r_No = NaN;
    end

    rxnList = [fullJoinRows(nF).r_DF, fullJoinRows(nF).r_DS, fullJoinRows(nF).r_Dp, fullJoinRows(nF).r_Di, ...
               fullJoinRows(nF).r_Cp, fullJoinRows(nF).r_Cg, fullJoinRows(nF).r_Cr, ...
               fullJoinRows(nF).r_Na, fullJoinRows(nF).r_Nb];
    if hasMerge
        rxnList = [rxnList, fullJoinRows(nF).r_No];
    end
    [fullJoinRows(nF).nGenes, hadOverlap] = motifGeneCount(model, rxnList);
    if hadOverlap
        fprintf('WARNING: motif %d (diff_idx=%d, cnt_idx=%d) has reactions sharing a gene.\n', ...
                nF, fullJoinRows(nF).diff_idx, fullJoinRows(nF).cnt_idx);
    end
end

if nF > 0
    [~, sortIdx] = sort([fullJoinRows.nGenes], 'ascend');
    fullJoinRows = fullJoinRows(sortIdx);
    FullMotifTable = struct2table(fullJoinRows);
    save('Block_3_nn_full_motifs.mat', 'FullMotifTable', 'fullJoinRows', 'hasMerge', 'NN_Pairs', 'NN_Convergent', 'NN_Full', ...
     'nnStep1', 'nnStep2', 'nnStep3','model');
    fprintf('Saved Block_3_nn_full_motifs.mat (%d full motifs, ranked by gene complexity)\n', nF);
else
    fprintf('No full motifs to save (nF = 0).\n');
end

fprintf('\nSaved nn_candidates.mat. Total Step 3 runtime: %.2f s\n', toc(tStart));

%\[text] ## ---------------- Helper functions ----------------
function genes = localGeneList(model, rxnIdx)
% Returns the unique set of gene tokens in a reaction's GPR rule,
% same convention used in Step_1 and Step_2 (AND/OR tokens stripped).
gpr = string(model.grRules{rxnIdx});
if strlength(strtrim(gpr)) == 0
    genes = strings(0,1);
    return;
end
toks = regexp(gpr, '[A-Za-z0-9_:\-\.]+', 'match');
toks = toks(~ismember(lower(toks), {'and','or'}));
genes = unique(string(toks));
genes = genes(:);   % force column vector so vertcat with strings(0,1) always works
end

function [n, hasOverlap] = motifGeneCount(model, rxnList)
% Returns the total gene count across a list of reactions. Since every
% motif here has already passed gene-disjointness checks (Step 1, Step 2's
% keepJ, Step 3's pairwise/upstream checks), summing per-reaction gene
% counts should equal the count of the union -- no two reactions should
% share a gene. hasOverlap flags it (without erroring) if that assumption
% is ever violated, which would indicate a bug upstream in the pipeline.
sumCount = 0;
allGenes = strings(0,1);
for i = 1:numel(rxnList)
    gpr = string(model.grRules{rxnList(i)});
    if strlength(strtrim(gpr)) == 0, continue; end
    toks = regexp(gpr, '[A-Za-z0-9_:\-\.]+', 'match');
    toks = toks(~ismember(lower(toks), {'and','or'}));
    genesHere = unique(string(toks));
    sumCount = sumCount + numel(genesHere);
    allGenes = [allGenes; genesHere(:)]; %#ok<AGROW>
end
n = sumCount;
hasOverlap = (numel(unique(allGenes)) ~= sumCount);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:81fe3904]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"No validated joint differentiator+counter solutions in workspace. Run Step_2_Counter_Mapping.m first (need variable joint_passed from its Step 5)."}}
%---
