%% Counter_Mapping.m
%
% Step 2 of the metabolic mapping: find chemical-counter implementations in
% iMM904, then match each result to the differentiator candidates from
% Differentiator_Mapping_anchorRDi_4rxn.m via the shared metabolite m_d.
%
% Counter CRN topology (3 reactions):
%   r_Cp:  m_d   -> m_c                  (m_d converted to counter output)
%   r_Cg:  m_d   + m_Cs -> sink          (m_d sequesters m_Cs to gate reset off)
%   r_Cr:  m_c   + m_Cs -> sink          (m_c reset by m_Cs)
%
% The S-pattern, in the same convention as the differentiator search:
%   row m_d  : negative entries in r_Cp and r_Cg
%   row m_c  : positive entry in r_Cp; negative entry in r_Cr
%   row m_Cs : negative entries in r_Cr and r_Cg (cofactor co-consumed in both)
%
% STRATEGY (Variant 2): global counter search, then match to differentiator.
%   STEP 1. Enumerate every yeast reaction r_Cp that has at least one substrate
%           (= m_d candidate) and produces at least one M_internal-style m_c.
%   STEP 2. For each (r_Cp, m_d, m_c) found, look for r_Cg consuming both m_d
%           and some m_Cs (any cofactor with relaxed-producer constraint).
%   STEP 3. For the same m_Cs, look for r_Cr consuming both m_c and m_Cs.
%   STEP 4. Apply enzyme-mediation, gene-disjointness, currency filters.
%   STEP 5. Match each counter (by its m_d) to differentiator candidates
%           whose m_d is the same yeast metabolite. Produce joint shortlist.
%
% ASSUMES: model loaded; results from differentiator search in `results_diff`
% (or pass an empty array to skip matching and just report all counters).

clc; clearvars -except model results_diff;

% Map differentiator results into a known variable name for this script
if exist('results_diff','var') && ~isempty(results)
    results_diff = results;
else
    error(['No differentiator results in workspace. Run ', ...
           'Differentiator_Mapping_anchorRDi_4rxn.m first.']);
end

%% Configuration
MAX_PRODUCERS_MC   = 1;     % m_c controllability: single-producer in iMM904
MAX_PRODUCERS_MCS  = inf;   % m_Cs is a buffered cofactor; many producers OK
EXCLUDE_CURRENCY   = true;  % drop water/protons/etc from substrate pairings
REQUIRE_GPR        = true;  % enzyme-mediated reactions only
ALLOW_MD_AS_MCS    = false; % m_d ≠ m_Cs (different roles)
MAX_HITS           = 1000;

%% Setup -- assumes model already loaded
S = model.S;
[nMets, nRxns] = size(S);
fprintf('Counter search on iMM904: %d mets, %d rxns.\n', nMets, nRxns);

% Currency-met filter
badMetPatterns = ["h[", "h2o[", "atp[", "adp[", "amp[", "pi[", "ppi[", ...
                  "nad[", "nadh[", "nadp[", "nadph[", ...
                  "co2[", "o2[", "nh3[", "nh4[", "coa[", ...
                  "fad[", "fadh[", "fadh2["];
isBadMet = false(nMets,1);
if EXCLUDE_CURRENCY
    for k = 1:numel(badMetPatterns)
        isBadMet = isBadMet | contains(lower(string(model.mets)), badMetPatterns(k));
    end
end

% GPR flag
hasGPR = false(nRxns,1);
for r = 1:nRxns
    rule = strtrim(string(model.grRules{r}));
    if strlength(rule) > 0 && ~strcmpi(rule, "s0001")
        hasGPR(r) = true;
    end
end

% Producer counts
nProducers = full(sum(S > 0, 2));

% Substrate sets
substrateSet = cell(nRxns,1);
productSet   = cell(nRxns,1);
for r = 1:nRxns
    substrateSet{r} = find(S(:,r) < 0);
    productSet{r}   = find(S(:,r) > 0);
end

%% =========================================================
% STEP 1: enumerate (r_Cp, m_d, m_c) triples
% =========================================================
% r_Cp is any enzyme-mediated reaction that has:
%   - at least one non-currency substrate s_md  (m_d candidate)
%   - at least one non-currency product s_mc with single-producer in iMM904
% Crucially we do NOT require r_Cp to be bi-substrate; it can be a simple
% conversion. m_d does not have to be the only substrate.

step1 = struct('m_d',{},'m_c',{},'r_Cp',{});
n1 = 0;

for r_Cp = 1:nRxns
    if REQUIRE_GPR && ~hasGPR(r_Cp), continue; end
    subs = setdiff(substrateSet{r_Cp}, find(isBadMet));
    if isempty(subs), continue; end
    prods = setdiff(productSet{r_Cp}, find(isBadMet));
    if isempty(prods), continue; end

    for i = 1:numel(subs)
        s_md = subs(i);
        for j = 1:numel(prods)
            s_mc = prods(j);
            if s_mc == s_md, continue; end
            if nProducers(s_mc) < 1 || nProducers(s_mc) > MAX_PRODUCERS_MC, continue; end
            n1 = n1 + 1;
            step1(n1).m_d  = s_md;
            step1(n1).m_c  = s_mc;
            step1(n1).r_Cp = r_Cp;
            if n1 >= MAX_HITS, break; end
        end
        if n1 >= MAX_HITS, break; end
    end
    if n1 >= MAX_HITS, break; end
end
fprintf('Step 1: %d (r_Cp, m_d, m_c) triples\n', n1);

%% =========================================================
% STEP 2: for each (r_Cp, m_d, m_c), find r_Cg co-consuming m_d + m_Cs
% =========================================================
step2 = struct('m_d',{},'m_c',{},'m_Cs',{},'r_Cp',{},'r_Cg',{});
n2 = 0;

for k = 1:numel(step1)
    rec = step1(k);
    s_md = rec.m_d;

    % Reactions that consume m_d
    md_consumers = find(S(s_md,:) < 0);
    if REQUIRE_GPR
        md_consumers = md_consumers(hasGPR(md_consumers));
    end
    md_consumers = setdiff(md_consumers, rec.r_Cp);   % r_Cg ≠ r_Cp

    for u = 1:numel(md_consumers)
        r_Cg = md_consumers(u);
        % Other substrates of r_Cg = candidates for m_Cs
        cg_subs = setdiff(substrateSet{r_Cg}, [find(isBadMet); s_md]);
        for v = 1:numel(cg_subs)
            s_mCs = cg_subs(v);
            if ~ALLOW_MD_AS_MCS && s_mCs == s_md, continue; end
            if s_mCs == rec.m_c, continue; end
            if nProducers(s_mCs) < 1 || nProducers(s_mCs) > MAX_PRODUCERS_MCS, continue; end
            n2 = n2 + 1;
            step2(n2).m_d  = s_md;
            step2(n2).m_c  = rec.m_c;
            step2(n2).m_Cs = s_mCs;
            step2(n2).r_Cp = rec.r_Cp;
            step2(n2).r_Cg = r_Cg;
            if n2 >= MAX_HITS, break; end
        end
        if n2 >= MAX_HITS, break; end
    end
    if n2 >= MAX_HITS, break; end
end
fprintf('Step 2: %d (r_Cp, r_Cg, m_d, m_c, m_Cs) tuples\n', n2);

%% =========================================================
% STEP 3: find r_Cr co-consuming m_c + m_Cs
% =========================================================
step3 = struct('m_d',{},'m_c',{},'m_Cs',{}, ...
               'r_Cp',{},'r_Cg',{},'r_Cr',{});
n3 = 0;

for k = 1:numel(step2)
    rec = step2(k);
    % consumers of m_c
    mc_consumers = find(S(rec.m_c,:) < 0);
    if REQUIRE_GPR
        mc_consumers = mc_consumers(hasGPR(mc_consumers));
    end
    mc_consumers = setdiff(mc_consumers, [rec.r_Cp, rec.r_Cg]);

    for u = 1:numel(mc_consumers)
        r_Cr = mc_consumers(u);
        % Require r_Cr to also consume m_Cs
        if S(rec.m_Cs, r_Cr) >= 0, continue; end
        n3 = n3 + 1;
        step3(n3).m_d  = rec.m_d;
        step3(n3).m_c  = rec.m_c;
        step3(n3).m_Cs = rec.m_Cs;
        step3(n3).r_Cp = rec.r_Cp;
        step3(n3).r_Cg = rec.r_Cg;
        step3(n3).r_Cr = r_Cr;
        if n3 >= MAX_HITS, break; end
    end
    if n3 >= MAX_HITS, break; end
end
fprintf('Step 3: %d full counter triples (r_Cp, r_Cg, r_Cr)\n', n3);

if n3 == 0
    fprintf('\nNo counters found. Diagnostic:\n');
    % Where did Step 2 die?
    nWith_mDs_consumer = 0;
    for k = 1:numel(step1)
        if any(S(step1(k).m_d,:) < 0)
            nWith_mDs_consumer = nWith_mDs_consumer + 1;
        end
    end
    fprintf('  Of %d step-1 triples, %d had at least one m_d consumer (other than r_Cp).\n', ...
            numel(step1), nWith_mDs_consumer);
    return;
end

%% =========================================================
% STEP 4: gene-disjointness filter on (r_Cp, r_Cg, r_Cr)
% =========================================================
keep4 = false(numel(step3),1);
overlapInfo = strings(numel(step3),1);

for k = 1:numel(step3)
    rec = step3(k);
    rxns = [rec.r_Cp, rec.r_Cg, rec.r_Cr];
    geneSets = cell(3,1);
    for i = 1:3
        gpr = string(model.grRules{rxns(i)});
        if strlength(gpr) == 0
            geneSets{i} = strings(0,1);
        else
            toks = regexp(gpr, '[A-Za-z0-9_:\-\.]+', 'match');
            toks = toks(~ismember(lower(toks), {'and','or'}));
            geneSets{i} = unique(string(toks));
        end
    end
    isDisj = true;
    msg = "";
    labels = {'r_Cp','r_Cg','r_Cr'};
    for i = 1:3
        for j = i+1:3
            common = intersect(geneSets{i}, geneSets{j});
            if ~isempty(common)
                isDisj = false;
                msg = msg + sprintf("%s<>%s share %s; ", labels{i}, labels{j}, ...
                                    strjoin(common,','));
            end
        end
    end
    keep4(k) = isDisj;
    overlapInfo(k) = msg;
end
step4 = step3(keep4);
fprintf('Step 4: %d counters after gene-disjointness filter\n', numel(step4));

%% =========================================================
% STEP 5: match counters to differentiator candidates by shared m_d
% =========================================================
joint = struct('diff_idx',{},'cnt_idx',{}, ...
               'm_in',{},'m_Df',{},'m_d',{},'m_Ds',{},'m_c',{},'m_Cs',{}, ...
               'r_DF',{},'r_DS',{},'r_Dp',{},'r_Di',{}, ...
               'r_Cp',{},'r_Cg',{},'r_Cr',{});
nJoint = 0;

% Collect differentiator m_d values
diff_md = arrayfun(@(s) s.m_d, results_diff);

for c = 1:numel(step4)
    cnt = step4(c);
    matchingDiffs = find(diff_md == cnt.m_d);
    for u = 1:numel(matchingDiffs)
        d_idx = matchingDiffs(u);
        d = results_diff(d_idx);
        nJoint = nJoint + 1;
        joint(nJoint).diff_idx = d_idx;
        joint(nJoint).cnt_idx  = c;
        joint(nJoint).m_in     = d.m_in;
        joint(nJoint).m_Df     = d.m_Df;
        joint(nJoint).m_d      = d.m_d;
        joint(nJoint).m_Ds     = d.m_Ds;
        joint(nJoint).m_c      = cnt.m_c;
        joint(nJoint).m_Cs     = cnt.m_Cs;
        joint(nJoint).r_DF     = d.r_DF;
        joint(nJoint).r_DS     = d.r_DS;
        joint(nJoint).r_Dp     = d.r_Dp;
        joint(nJoint).r_Di     = d.r_Di;
        joint(nJoint).r_Cp     = cnt.r_Cp;
        joint(nJoint).r_Cg     = cnt.r_Cg;
        joint(nJoint).r_Cr     = cnt.r_Cr;
    end
end
fprintf('Step 5: %d joint (differentiator + counter) candidates\n', nJoint);

%% =========================================================
% Joint-level filters: ensure differentiator and counter reactions
% don't overlap or share genes.
% =========================================================
keepJ = false(nJoint,1);
jointOverlapInfo = strings(nJoint,1);

for k = 1:nJoint
    j = joint(k);
    allRxns = [j.r_DF, j.r_DS, j.r_Dp, j.r_Di, j.r_Cp, j.r_Cg, j.r_Cr];

    % No reaction reused across blocks
    if numel(unique(allRxns)) ~= 7
        jointOverlapInfo(k) = "reaction reused across blocks";
        continue;
    end

    % Gene-disjointness across all 7
    geneSets = cell(7,1);
    for i = 1:7
        gpr = string(model.grRules{allRxns(i)});
        if strlength(gpr) == 0
            geneSets{i} = strings(0,1);
        else
            toks = regexp(gpr, '[A-Za-z0-9_:\-\.]+', 'match');
            toks = toks(~ismember(lower(toks), {'and','or'}));
            geneSets{i} = unique(string(toks));
        end
    end
    isDisj = true;
    msg = "";
    labels = {'r_DF','r_DS','r_Dp','r_Di','r_Cp','r_Cg','r_Cr'};
    for i = 1:7
        for jj = i+1:7
            common = intersect(geneSets{i}, geneSets{jj});
            if ~isempty(common)
                isDisj = false;
                msg = msg + sprintf("%s<>%s share %s; ", labels{i}, labels{jj}, ...
                                    strjoin(common,','));
            end
        end
    end
    keepJ(k) = isDisj;
    if ~isDisj
        jointOverlapInfo(k) = msg;
    end
end

joint_passed = joint(keepJ);
nJpass = numel(joint_passed);
fprintf('Joint candidates passing 7-reaction gene-disjointness: %d / %d\n', nJpass, nJoint);

%% =========================================================
% Build a JointTable for inspection
% =========================================================
if nJpass > 0
    rows = repmat(struct(),nJpass,1);
    for k = 1:nJpass
        j = joint_passed(k);
        rows(k).diff_idx = j.diff_idx;
        rows(k).cnt_idx  = j.cnt_idx;
        rows(k).m_in   = string(model.mets{j.m_in});
        rows(k).m_Df   = string(model.mets{j.m_Df});
        rows(k).m_d    = string(model.mets{j.m_d});
        rows(k).m_Ds   = string(model.mets{j.m_Ds});
        rows(k).m_c    = string(model.mets{j.m_c});
        rows(k).m_Cs   = string(model.mets{j.m_Cs});
        rows(k).r_DF   = string(model.rxns{j.r_DF});
        rows(k).r_DS   = string(model.rxns{j.r_DS});
        rows(k).r_Dp   = string(model.rxns{j.r_Dp});
        rows(k).r_Di   = string(model.rxns{j.r_Di});
        rows(k).r_Cp   = string(model.rxns{j.r_Cp});
        rows(k).r_Cg   = string(model.rxns{j.r_Cg});
        rows(k).r_Cr   = string(model.rxns{j.r_Cr});
    end
    JointTable = struct2table(rows);
    fprintf('\n=== Joint differentiator + counter table (first %d rows) ===\n', ...
            min(10,nJpass));
    disp(JointTable(1:min(10,nJpass),:));
    save('joint_diff_counter_candidates.mat', 'JointTable', 'joint_passed', ...
         'step3', 'step4', 'joint');
    fprintf('\nSaved joint_diff_counter_candidates.mat\n');
else
    fprintf(['\nNo joint candidate passed gene-disjointness across all 7 reactions.\n', ...
             'Try inspecting `joint` directly to see what overlaps exist.\n']);
end
