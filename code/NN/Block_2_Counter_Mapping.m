%[text] # Counter\_Mapping.m
%[text] Step 2 of the metabolic mapping: find chemical-counter implementations in iMM904, then match each result to the differentiator candidates from Differentiator\_Mapping\_anchorRDi\_4rxn.m via the shared metabolite m\_d.
%[text] Counter CRN topology (3 reactions): r\_Cp: m\_d -\> m\_c (m\_d converted to counter output) r\_Cg: m\_d + m\_Cs -\> sink (m\_d sequesters m\_Cs to gate reset off) r\_Cr: m\_c + m\_Cs -\> sink (m\_c reset by m\_Cs)
%[text] The S-pattern, in the same convention as the differentiator search: row m\_d : negative entries in r\_Cp and r\_Cg row m\_c : positive entry in r\_Cp; negative entry in r\_Cr row m\_Cs : negative entries in r\_Cr and r\_Cg (cofactor co-consumed in both)
%[text] STRATEGY (Variant 2): global counter search, then match to differentiator. STEP 1. Enumerate every yeast reaction r\_Cp that has at least one substrate (= m\_d candidate) and produces at least one M\_internal-style m\_c. STEP 2. For each (r\_Cp, m\_d, m\_c) found, look for r\_Cg consuming both m\_d and some m\_Cs (any cofactor with relaxed-producer constraint). STEP 3. For the same m\_Cs, look for r\_Cr consuming both m\_c and m\_Cs. STEP 4. Apply enzyme-mediation, gene-disjointness, currency filters. STEP 5. Match each counter (by its m\_d) to differentiator candidates whose m\_d is the same yeast metabolite. Produce joint shortlist.
%[text] ASSUMES: model loaded; results from differentiator search in `results\_diff` (or pass an empty array to skip matching and just report all counters).
% clc; clearvars -except model results_diff;

%load differentiator_candidates_anchorRDi_4rxn.mat

% Map differentiator results into a known variable name for this script
if exist('results_diff','var') && ~isempty(results_diff) %[output:group:428e0322]
    fprintf('Loaded the differentiator') %[output:762ffa34]

else
    error(['No differentiator results in workspace. Run ', ...
           'Differentiator_Mapping_anchorRDi_4rxn.m first.']);
end %[output:group:428e0322]
%[text] ## =========================================================
%[text] STEP 1: enumerate (r\_Cp, m\_d, m\_c) triples ========================================================= r\_Cp is any enzyme-mediated reaction that has: - at least one non-currency substrate s\_md (m\_d candidate) - at least one non-currency product s\_mc with single-producer in iMM904 Crucially we do NOT require r\_Cp to be bi-substrate; it can be a simple conversion. m\_d does not have to be the only substrate.
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
fprintf('Step 1: %d (r_Cp, m_d, m_c) triples\n', n1); %[output:10078da4]
%[text] ## =========================================================
%[text] STEP 2: for each (r\_Cp, m\_d, m\_c), find r\_Cg co-consuming m\_d + m\_Cs =========================================================
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
fprintf('Step 2: %d (r_Cp, r_Cg, m_d, m_c, m_Cs) tuples\n', n2); %[output:2268009a]
%[text] ## =========================================================
%[text] STEP 3: find r\_Cr co-consuming m\_c + m\_Cs =========================================================
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
fprintf('Step 3: %d full counter triples (r_Cp, r_Cg, r_Cr)\n', n3); %[output:92f15b5e]

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
%[text] ## =========================================================
%[text] STEP 4: gene-disjointness filter on (r\_Cp, r\_Cg, r\_Cr) =========================================================
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
fprintf('Step 4: %d counters after gene-disjointness filter\n', numel(step4)); %[output:3564bbb3]

%% Printing solutions at Step 4
% if numel(step4) > 0
% 
%     fprintf('\n=== First candidates (4-reaction differentiator) ===\n');
%     for i = 1:numel(step4) %min(5, nHits)
%         rec = step4(i);
%         fprintf('\n--- Candidate %d ---\n', i);
%         fprintf('  m_d   = %s | %s\n', model.mets{rec.m_d},   model.metNames{rec.m_d});
%         fprintf('  m_c  = %s | %s\n', model.mets{rec.m_c},  model.metNames{rec.m_c});
%         fprintf('  m_Cs  = %s | %s\n', model.rxns{rec.m_Cs},  model.rxnNames{rec.m_Cs});
%         fprintf('  GPR(r_DF) = %s\n', model.grRules{rec.m_d});
%         fprintf('  GPR(r_DS) = %s\n', model.grRules{rec.m_c});
%         fprintf('  GPR(r_Dp) = %s\n', model.grRules{rec.m_Cs});
%     end
%     save('differentiator_candidates_anchorRDi_4rxn.mat', 'Motifs_Diff_4rxn', 'results_diff');
% end
%[text] ## =========================================================
%[text] STEP 5: match counters to differentiator candidates by shared m\_d =========================================================
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
fprintf('Step 5: %d joint (differentiator + counter) candidates\n', nJoint); %[output:00c2e1ba]
%[text] ## =========================================================
%[text] Joint-level filters: ensure differentiator and counter reactions don't overlap or share genes. =========================================================
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
fprintf('Joint candidates passing 7-reaction gene-disjointness: %d / %d\n', nJpass, nJoint); %[output:890139b4]
%[text] ## =========================================================
%[text] Build a JointTable for inspection =========================================================
if nJpass > 0 %[output:group:88f84177]
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
    fprintf('\n=== Joint differentiator + counter table (first %d rows) ===\n', ... %[output:08093c42]
            min(nJpass,nJpass)); %[output:08093c42]
    disp(JointTable(1:min(nJpass,nJpass),:)); %[output:948bb2ed]
    save('Block_2_joint_diff_counter_candidates.mat', 'JointTable', 'joint_passed', ...
         'step3', 'step4', 'joint');
    fprintf('\nSaved joint_diff_counter_candidates.mat\n'); %[output:606d1ac5]
else
    fprintf(['\nNo joint candidate passed gene-disjointness across all 7 reactions.\n', ...
             'Try inspecting `joint` directly to see what overlaps exist.\n']);
end %[output:group:88f84177]

if numel(step4) > 0 %[output:group:8d59eaae]
    Motifs_Cnt = struct2table(step4);
    fprintf('\n=== Counter candidates (3-reaction CRN) ===\n'); %[output:7654d8dc]
    for i = 1:numel(step4)    % use 1:min(5, numel(step4)) if you only want first few
        rec = step4(i);
        fprintf('\n--- Counter candidate %d ---\n', i); %[output:596067ef] %[output:12e70884] %[output:727af5ac] %[output:243092d5] %[output:2cd283dc] %[output:963f3445] %[output:7f450d27] %[output:8d70a876] %[output:1074bcf1] %[output:72225947] %[output:7cae1de4] %[output:91792c1a] %[output:99229002] %[output:5febab0d] %[output:73fc6839] %[output:92cee311] %[output:32b3c01d] %[output:7127c606] %[output:5053a5b1] %[output:7ddb315d] %[output:5f932f1d] %[output:95fbf5f6] %[output:8369850b] %[output:65c9fedc] %[output:7d1d72fe] %[output:304b9fa1] %[output:09b18094] %[output:7e562fbe] %[output:21c461f3] %[output:6f9fd200] %[output:5bd01679] %[output:50f104cb] %[output:88c79437] %[output:8327b81c] %[output:12cdf022] %[output:354a9704] %[output:1e6f07e9] %[output:9a811f1d] %[output:8628159f] %[output:72953939] %[output:06f1ec41] %[output:38b6a466] %[output:88d71ed5] %[output:579c7420] %[output:348c6870] %[output:5fe43a3a] %[output:6f46adf2] %[output:8d737617] %[output:1cbc65dd] %[output:8f894226] %[output:21a82635] %[output:30b82c99] %[output:6057c254] %[output:688d4779] %[output:74718a55] %[output:870f125e] %[output:81b22fde] %[output:4323dbc6] %[output:72b63095] %[output:3f78c1d8] %[output:1e9aff1d] %[output:0daea8c6] %[output:512d87b4] %[output:105479e6] %[output:11088e9f] %[output:9882674e] %[output:4c43160c] %[output:3a1713f4] %[output:386b408f] %[output:86ab39cd] %[output:73838060] %[output:33bae055] %[output:7dee8915] %[output:1b418a9b] %[output:8f0eab80] %[output:6ef810fe] %[output:19021326] %[output:8ded4164] %[output:257e13db] %[output:37cd056e] %[output:6a1050ee] %[output:264bc238] %[output:8914c546] %[output:8695b7e1] %[output:6526e319] %[output:342073c7] %[output:98865c68] %[output:0c5c47fc] %[output:2e2bca33] %[output:1a9d7d9a] %[output:70367429] %[output:596e1c5d] %[output:455bc3c0] %[output:84a2d9ca] %[output:220a7d77] %[output:988c3d90] %[output:05792a32] %[output:720c2e61] %[output:785e4e5f] %[output:11b95d13] %[output:883569b1] %[output:886c7f11] %[output:95a18892] %[output:94daa8cf] %[output:9916dc31] %[output:246819d5] %[output:603a87d4] %[output:4fdba30c] %[output:1c16a5e6] %[output:7af38729] %[output:36ff2443] %[output:10d76a8c] %[output:183ac662] %[output:6dc480a7] %[output:0c152fd4] %[output:56da8b1b] %[output:778f141c] %[output:0016392f] %[output:8cf2de9b] %[output:21aa9f25] %[output:8bc75d47] %[output:56573bd7] %[output:5e632d0d] %[output:72fcbc73] %[output:4a5fe891] %[output:5bce843c] %[output:8d649804] %[output:0bcf98cf] %[output:2a1f0d84] %[output:471da9e7] %[output:378d728b] %[output:0d3c9f4a] %[output:1af9645d] %[output:3d49a259] %[output:23a95609] %[output:27b732e5] %[output:7a61a555] %[output:28d5d04e] %[output:8e03de47] %[output:0135cfa5] %[output:11ae98d3] %[output:8dc715b0] %[output:0f531385] %[output:04ffd46a] %[output:5c567304] %[output:9fa04f30] %[output:2b6c86ab] %[output:180cbeee] %[output:7e21387e] %[output:48834740] %[output:7005970f] %[output:37300a7f] %[output:4512eaaa] %[output:04164d86] %[output:7fc3f50e] %[output:2d97170b] %[output:78d38b2e] %[output:9f854ea8] %[output:5ab8b04c] %[output:478b990f] %[output:5ae71839] %[output:9189419d] %[output:5df7cb6e] %[output:63976968] %[output:5a7d35a3] %[output:054d7b2b] %[output:845da39a] %[output:3984af47] %[output:6e3d93e2] %[output:5ff23c50] %[output:3c08ef07] %[output:79a32d79] %[output:08a743d1] %[output:73156168] %[output:729fa55b] %[output:138cb531] %[output:1141e437] %[output:05b99f0e] %[output:5bf75a91] %[output:5a946bb4] %[output:81b0fada] %[output:09a71e6a] %[output:69b9ce04] %[output:88052660] %[output:2b51d49a] %[output:6b48f107] %[output:84dfadc1] %[output:87aca4ae] %[output:45e4caa8] %[output:915c3a8e] %[output:36e116f0] %[output:8f761848] %[output:595a949c] %[output:1e2eece8] %[output:7b3feb59] %[output:3cf2d74e] %[output:80514bbf] %[output:209484ba] %[output:3b7ea8b3] %[output:4dcdb699] %[output:64b03aab] %[output:0ea97bcf] %[output:4899ffa2] %[output:0c35caaf] %[output:54cae95f] %[output:260655bb] %[output:0b4e31bb] %[output:7a8ae9cb] %[output:007d59da] %[output:08815d7c] %[output:6fb16ea8] %[output:71feb355] %[output:02c4464b] %[output:948cbbf3] %[output:0a02b2be] %[output:2f0a73a7] %[output:027b3cad] %[output:64a3fcd6] %[output:39ede1e1] %[output:070cfd19] %[output:34d647ca] %[output:9f040b25] %[output:3a933269] %[output:350d437b] %[output:66e382ce] %[output:96de8d38] %[output:12bbc1d8] %[output:732e7a90] %[output:67f54263] %[output:11b49691] %[output:451d6d8c] %[output:3e290cbb] %[output:8d47466e] %[output:1bd6de1e] %[output:722a5704] %[output:044dbf8f] %[output:2429dbb9] %[output:24c45d42] %[output:229c0d4a] %[output:862f941b] %[output:027c93da] %[output:52b0dcb4] %[output:65d6f6fd] %[output:2b321227] %[output:4818b790] %[output:8a399894] %[output:21039420] %[output:8deab526] %[output:1c7e9725] %[output:447c6b55] %[output:2bf1d9ef] %[output:2098807d] %[output:71afa7d4] %[output:52961f30] %[output:044cc5b4] %[output:75bd85b0] %[output:3a5566a9] %[output:19ec6b1f] %[output:1f8ec22b] %[output:4618b100] %[output:34c49aa0] %[output:85038a64] %[output:36306094] %[output:2ea76327] %[output:3f2d921d] %[output:4695f474] %[output:5adde840] %[output:0506b13c] %[output:65069091] %[output:4f824830] %[output:4965bf67] %[output:56e7ed0b] %[output:19141129] %[output:55e383fb] %[output:00378585] %[output:316c594d] %[output:9f99e26a] %[output:36e44eca] %[output:8291cf97] %[output:4c967884] %[output:85e4b623] %[output:54901dc5] %[output:72596dc1] %[output:8b1a0cb4] %[output:31eb9918] %[output:747867af] %[output:6c8ba680] %[output:32bdb8a7] %[output:3ad8bbd5] %[output:76d013a6] %[output:9959c62d] %[output:7d03650c] %[output:78e3325b] %[output:9bcf1135] %[output:8fe64110] %[output:71654441] %[output:571c50e7] %[output:76acd579] %[output:1b898cb1] %[output:18bd9a5e] %[output:828a9c40] %[output:033087f2] %[output:1c059e6b] %[output:237c5433] %[output:5ea6ac57] %[output:0208c275] %[output:00155077] %[output:0b8bba99] %[output:2693bec6] %[output:1f10c3f7] %[output:7e9f5eb5] %[output:81e95f02] %[output:5fa65f8d] %[output:4b2686c3] %[output:79107700] %[output:50aecc58] %[output:2795b153] %[output:8ce780c4] %[output:9e31f3aa] %[output:1160e836] %[output:78a8f627] %[output:226fafcc] %[output:6eef967b] %[output:2fd9cf95] %[output:27404310] %[output:56dafdaf] %[output:7901937a] %[output:135ccc49] %[output:11e7d7ac] %[output:7173e99e] %[output:7a0322b0] %[output:23b229ce] %[output:23dd2c31] %[output:02db2a31] %[output:20f7be31] %[output:481897eb] %[output:8e8dfbae] %[output:8a6cf90a] %[output:506a8dff] %[output:17e01612] %[output:1bad0bca] %[output:9ed31a80] %[output:6a77370a] %[output:1dfc9a15] %[output:09d166fe] %[output:98e412e3] %[output:2ce5a054] %[output:207adbd3] %[output:16bf7568] %[output:5cc41b5d] %[output:44710edc] %[output:017bd402] %[output:70cc05a1] %[output:01bd9336] %[output:93b76927] %[output:9276cf8a] %[output:866fc1ab] %[output:16e04d58] %[output:6ebc0b3b] %[output:6f944047] %[output:9936dae0] %[output:118a6e79] %[output:1652a8b1] %[output:7c5776c6] %[output:89baa316] %[output:2a4ffe6e] %[output:319447d7] %[output:780c4a95] %[output:6b7c2229] %[output:1e5190ab] %[output:86b18adf] %[output:092dbd8c] %[output:98fab286] %[output:9599c3d3] %[output:578ad184] %[output:2fe5e793] %[output:22d22eba] %[output:2e5bcdf3] %[output:032282c0] %[output:42f622d5] %[output:4f304123] %[output:5e40e727] %[output:12faa176] %[output:9b2bbe5d] %[output:25e5e4d6] %[output:99418a8d] %[output:0919b219] %[output:2eb01979] %[output:133a26ce] %[output:30a5a747] %[output:4c7a7e8d] %[output:95579a84] %[output:3c59d843] %[output:776cce4e] %[output:2f9150b2] %[output:051092cc] %[output:5770aa16] %[output:7bcfb166] %[output:5f3d7ae6] %[output:0f1311d0] %[output:760f65b1] %[output:7a82b0af] %[output:2c59e2fb] %[output:2278d136] %[output:5109bf7a] %[output:842d6512] %[output:5fd3b80a] %[output:64ae93eb] %[output:5d3be5ea] %[output:40d89321] %[output:2b1b62dd] %[output:781cf433] %[output:3fbe1002] %[output:23e2b252] %[output:12adb112] %[output:5c31c3fb] %[output:5c7ebe62] %[output:6531f418] %[output:133458a8] %[output:6c299f58] %[output:8d152646] %[output:0b107991] %[output:08905c89] %[output:4dc81ee0] %[output:1c2f23cf] %[output:174b7e3a] %[output:21b07c49] %[output:42a67e1d] %[output:6289aaa3] %[output:8e98dea8] %[output:6e089475] %[output:0fdcc375] %[output:01754cc2] %[output:9f52a311] %[output:7d1b6684] %[output:8b5b8f49] %[output:5ef6b4d8] %[output:3f537d17] %[output:2a579d89] %[output:943f8d29] %[output:63575f60] %[output:0b32938f] %[output:05f5471d] %[output:64d96122] %[output:1fa3274c] %[output:1f15f735] %[output:92d6b4ca] %[output:3689ad23] %[output:54818bec] %[output:220ba4ad] %[output:7495be47] %[output:884bcca6] %[output:6da8979f] %[output:19a880b3] %[output:98852bd5] %[output:5e921ffb] %[output:01b32f91] %[output:4e4622f1] %[output:953901c1] %[output:7b9c3205] %[output:57cb16b1] %[output:03be723b] %[output:3ddfe554] %[output:332ff50c] %[output:36283bc6] %[output:68e1387f] %[output:0c619111] %[output:780f2186] %[output:08a0e8bd] %[output:3ea3fa0a] %[output:0d3f4c24] %[output:0b24a315] %[output:4b16be22] %[output:62bd98ec] %[output:243378da] %[output:583a9124] %[output:87e2f1ba] %[output:6c095e46] %[output:6968500c] %[output:51988360] %[output:330eb5d3] %[output:6255bf45] %[output:472f9b40] %[output:810040b0] %[output:9b7472e0] %[output:1b3e1aa2] %[output:43dda632] %[output:69ebd86a] %[output:73a80d0f] %[output:54445450] %[output:08a8ddda] %[output:425b7e9d] %[output:7ddbbc7c] %[output:174914bc] %[output:8cd72412] %[output:84b0e411] %[output:65155b44] %[output:15efc94f] %[output:4d5e3bda]
        fprintf('  m_d  = %s | %s\n', model.mets{rec.m_d},  model.metNames{rec.m_d}); %[output:3d67aa74] %[output:2be6c316] %[output:6d05f825] %[output:938b7b8f] %[output:1cdd234c] %[output:77c477b4] %[output:779e6a55] %[output:4b597c27] %[output:206d0b47] %[output:83728a93] %[output:065b66b7] %[output:8010da10] %[output:2d363e7e] %[output:79fd2239] %[output:541bf84c] %[output:2a5a9fc3] %[output:72d000c2] %[output:3ef9b704] %[output:748f34b6] %[output:273374e9] %[output:8302329b] %[output:5ec285c4] %[output:452c666a] %[output:33dbaf23] %[output:9f1cd97b] %[output:978d4116] %[output:5ba43c97] %[output:02a0d147] %[output:378973be] %[output:3970db3f] %[output:67b272b9] %[output:15c3e391] %[output:1133a651] %[output:233411b7] %[output:74429697] %[output:211fdbc9] %[output:2a5c20c4] %[output:8dd3ee05] %[output:1a66b34f] %[output:549ad6c7] %[output:24778d34] %[output:75a7a70e] %[output:02a8ecda] %[output:8206a718] %[output:1e27c755] %[output:218eea8a] %[output:5db47e85] %[output:838866da] %[output:0f09f978] %[output:6871322b] %[output:21e39684] %[output:4732a548] %[output:1dcb4f15] %[output:668be0eb] %[output:7f0ac505] %[output:23e2da42] %[output:6495fb90] %[output:2600ac11] %[output:542a99b2] %[output:38e7ad6d] %[output:5676361d] %[output:7fbeecac] %[output:523b4e35] %[output:08acdf06] %[output:07afa78c] %[output:43d4a80b] %[output:9dd3adbc] %[output:3e1d8caa] %[output:73ac57eb] %[output:076e2ccc] %[output:4a8eb149] %[output:8b81c326] %[output:69efc685] %[output:2083147d] %[output:190f5cd8] %[output:3c80f1a9] %[output:860ea298] %[output:9ba4b911] %[output:444bccdd] %[output:315a7e33] %[output:34b9ea6b] %[output:596e285d] %[output:29da7539] %[output:9537147e] %[output:23a66deb] %[output:256d806b] %[output:286838a7] %[output:77cfd510] %[output:27a4dbd6] %[output:219edde0] %[output:27ca9601] %[output:976d7355] %[output:461eb6c7] %[output:00d6f576] %[output:7ce1d45b] %[output:28f22c9d] %[output:0c008d96] %[output:655f29d4] %[output:7c402e6e] %[output:7d973e83] %[output:46eaeaee] %[output:33fc7e28] %[output:8fc7c3c8] %[output:2faf9595] %[output:908f7b41] %[output:03b8749f] %[output:392c7b01] %[output:7bc2e00a] %[output:8bcf63c6] %[output:3feaa8fe] %[output:73e88950] %[output:43f6bc88] %[output:5e904ec9] %[output:3d5c32e6] %[output:94c89ef9] %[output:239407ba] %[output:7ae98b04] %[output:4d78d608] %[output:28c9364a] %[output:06f2243e] %[output:0a71171b] %[output:0231fcff] %[output:0f50bba1] %[output:5aa43ac1] %[output:7764ee68] %[output:6a31871b] %[output:405baf28] %[output:6c396e26] %[output:55a23363] %[output:600e874b] %[output:0121f727] %[output:751a463a] %[output:1ce17b05] %[output:20487279] %[output:3fa929de] %[output:28f7dd98] %[output:8f2d1a4d] %[output:9cdb9b35] %[output:7222ddf7] %[output:60c1be90] %[output:686eb1e9] %[output:6c957c80] %[output:342fd214] %[output:88158109] %[output:03c29dd4] %[output:0bf3ea6f] %[output:2805c352] %[output:9c124ea5] %[output:1f74a554] %[output:21f73b03] %[output:9ab76b2f] %[output:5234572c] %[output:265ff359] %[output:5f8ebce4] %[output:624dd1ae] %[output:91083374] %[output:1332a1fe] %[output:89511d57] %[output:30e0b0e7] %[output:5359be80] %[output:18c7b257] %[output:46dc9760] %[output:0b0d7033] %[output:471dad89] %[output:79057ff9] %[output:01e898d0] %[output:9090cee2] %[output:160539c4] %[output:96b732a2] %[output:4a408d6a] %[output:0eef1021] %[output:75149875] %[output:3298229e] %[output:2423cf4d] %[output:0c3fc659] %[output:1c880e2f] %[output:38d9c5fe] %[output:9ce498ff] %[output:69ed582e] %[output:902d4510] %[output:18977cf3] %[output:9ac35571] %[output:6bfcc808] %[output:83419f72] %[output:13afe2b8] %[output:2f122382] %[output:0b37a037] %[output:68429bd9] %[output:3dfda2bb] %[output:7d2a73b7] %[output:279524c9] %[output:26f6a218] %[output:9acbd008] %[output:1d1bafda] %[output:3622a6d0] %[output:9e7afa63] %[output:3c09b8cc] %[output:9e242f9a] %[output:3dcf9971] %[output:5d641a51] %[output:41eaf15c] %[output:8ce975b9] %[output:06250398] %[output:4c352c78] %[output:7acfccb3] %[output:367ed87d] %[output:8423a821] %[output:8183cda4] %[output:0395e44e] %[output:8168cb53] %[output:73f31171] %[output:743df2e0] %[output:60f78edd] %[output:054f3d5f] %[output:51ed4cff] %[output:007769bf] %[output:3acca036] %[output:63c22942] %[output:506b98d4] %[output:5618dd7f] %[output:41658d3d] %[output:2ac971f8] %[output:9184e244] %[output:94ba7200] %[output:22158923] %[output:33b402e4] %[output:83db84a0] %[output:9eb7508c] %[output:1fbfa608] %[output:1b385bd9] %[output:49796456] %[output:27a5fe7e] %[output:334d9f73] %[output:919949b5] %[output:32dcb09a] %[output:3b235e53] %[output:4d96c85b] %[output:8dd05024] %[output:2f09f4cb] %[output:2e945266] %[output:4d7c6bef] %[output:820b6611] %[output:7bc4937a] %[output:83344567] %[output:09ce14ab] %[output:75036819] %[output:0bcdbd07] %[output:2c085798] %[output:99644a19] %[output:11496f55] %[output:066721ff] %[output:01bd6126] %[output:3526ea44] %[output:42e79b6e] %[output:2edb9eaa] %[output:7466ce34] %[output:465eebd5] %[output:7716e15d] %[output:9aa0d725] %[output:277cc47a] %[output:5c609d61] %[output:23048b03] %[output:2c88a0c7] %[output:37e8cbf0] %[output:855a252b] %[output:2c2524e5] %[output:57ef575b] %[output:4b82a77b] %[output:09ee34a5] %[output:7cde0059] %[output:219a5da9] %[output:3fadaa40] %[output:441e34ba] %[output:10ea8f2c] %[output:18efc3e7] %[output:664d170c] %[output:7373ca0c] %[output:05c15ecc] %[output:37a897a4] %[output:87cc2884] %[output:52ab0fa5] %[output:6dbe1eb1] %[output:9a56a51c] %[output:51f001f1] %[output:9a303278] %[output:2a6c8654] %[output:6b4ea7f0] %[output:57406199] %[output:73e509e8] %[output:7fe90879] %[output:1ce62b41] %[output:6507f527] %[output:0671b825] %[output:0d59532b] %[output:6101ff12] %[output:751d3b39] %[output:051fec21] %[output:86ab6c9d] %[output:4ab0bb6e] %[output:3ea3026a] %[output:25039871] %[output:0329d734] %[output:09297d9a] %[output:1b35c56c] %[output:1782f8b6] %[output:91fe8fc0] %[output:5f2a4aca] %[output:7fbede02] %[output:89afb399] %[output:602b20f5] %[output:2c312842] %[output:4fad01e4] %[output:6c4fe67b] %[output:26686d6d] %[output:829fe5e0] %[output:58f9baf8] %[output:28fa7788] %[output:79d2eff3] %[output:257ecfcd] %[output:8e25644d] %[output:28a657c4] %[output:3eca4361] %[output:487974f7] %[output:5fc375b6] %[output:9401a413] %[output:3562f533] %[output:1abd9d73] %[output:2b5b22dd] %[output:75ab5cd0] %[output:2dcc2f06] %[output:8c2dac52] %[output:875511ae] %[output:9fc2812b] %[output:26decc7c] %[output:42d7c83a] %[output:89421343] %[output:48bb7b10] %[output:127b2b69] %[output:9128712b] %[output:8eedb8dd] %[output:5117c54a] %[output:58284af8] %[output:991721a1] %[output:67be94a0] %[output:15621b79] %[output:4e8dca6b] %[output:5493ff6c] %[output:82682ee8] %[output:9893cf01] %[output:345872d1] %[output:82e0a147] %[output:7724996e] %[output:482387fc] %[output:1b85ac05] %[output:55128a5d] %[output:095d27ff] %[output:2f527a98] %[output:584b9bdb] %[output:47e3099e] %[output:3c046561] %[output:8d8b385c] %[output:697a15f5] %[output:67aaa049] %[output:06d400aa] %[output:18e0fa50] %[output:4051cc8f] %[output:09d9fc48] %[output:10f83ef6] %[output:10858acd] %[output:7c722c12] %[output:0826466c] %[output:7e76e449] %[output:81bf6c2b] %[output:810b4fba] %[output:81cab6c4] %[output:713d68a9] %[output:1ea73661] %[output:17f22f49] %[output:175d72ea] %[output:0a9899d7] %[output:4eca035b] %[output:43238c70] %[output:2a48ab8d] %[output:0b98c157] %[output:22b6ab77] %[output:4f5bc31c] %[output:3a00fee6] %[output:0a4ef3b8] %[output:092762a1] %[output:6e184154] %[output:9e68d627] %[output:08f51f72] %[output:59e9a6af] %[output:26253693] %[output:7ecd3b36] %[output:0fbba460] %[output:81d5e56f] %[output:69ff599c] %[output:3f38afc3] %[output:5c1b939f] %[output:613597ba] %[output:111a9ea6] %[output:85ad1006] %[output:3534668e] %[output:3dca75c5] %[output:3c5a2db0] %[output:7262880d] %[output:58fcc7c1] %[output:323112ca] %[output:2bba4fa7] %[output:92eb782f] %[output:58fb12ed] %[output:59173ff9] %[output:49bccf6e] %[output:35cd31db] %[output:6385a270] %[output:57d863e4] %[output:858a4990] %[output:5fc9d260] %[output:34f94b64] %[output:49ffb4d8] %[output:8d7c480d] %[output:012da4b0] %[output:6dd12d38] %[output:9738338c] %[output:643166df] %[output:3e18db6f] %[output:4dc9bf56] %[output:0d6f8b35] %[output:63f835ca] %[output:496f46fd] %[output:99af7427] %[output:2dfbb5c1] %[output:72d43335] %[output:6b36550b] %[output:9d1f8744] %[output:68ad3163] %[output:92b99444] %[output:63fa5cb7] %[output:726cc0a7] %[output:2985f875] %[output:16932fa5] %[output:4b8c1648] %[output:67f4ba6c] %[output:2d43b53a] %[output:02a73106] %[output:97f4d66c] %[output:73643e93] %[output:7a38413b] %[output:44459907] %[output:316afea3] %[output:717da4ae] %[output:71ea91f8] %[output:292bfb06] %[output:0cc8463a] %[output:51bfc35f] %[output:1383aee1] %[output:658d793f] %[output:562cef7e] %[output:8e4cc1ed] %[output:88bd5d89] %[output:15371d92] %[output:3758b38f] %[output:84cf394f] %[output:3f12bcb3] %[output:003b351b] %[output:60fec791] %[output:1fe0333b] %[output:7b3f5128] %[output:07d121f2] %[output:4a7e00bb] %[output:2417f3ea] %[output:8afbecf8] %[output:72e449aa] %[output:0e1bb0f9] %[output:7b4f3964] %[output:0b524a00] %[output:0621487a] %[output:9f37109b] %[output:4e770d23] %[output:8c69ac28] %[output:1e4b9ec3] %[output:12e01a0c] %[output:8eda4cc1] %[output:5f3ba98c] %[output:9123b4ec] %[output:9706c1c1] %[output:7c5cf69b] %[output:1baea7ba] %[output:4c0d1604] %[output:2f2298ee] %[output:88375cbf] %[output:02275730] %[output:00489b40] %[output:648ef5fd] %[output:143ae5fb] %[output:7bf730ad] %[output:58461595] %[output:2255eded]
        fprintf('  m_c  = %s | %s\n', model.mets{rec.m_c},  model.metNames{rec.m_c}); %[output:5deabfec] %[output:30e55f88] %[output:810aa9e7] %[output:7246e0e8] %[output:273f522a] %[output:33d215d8] %[output:327c6217] %[output:98141f95] %[output:8587b5dc] %[output:8d30f0c4] %[output:9a4373b9] %[output:69666f37] %[output:30ba1ed4] %[output:0b706432] %[output:8500c489] %[output:4195e092] %[output:6d10fd32] %[output:99b3537d] %[output:740695ea] %[output:85352222] %[output:495bcb17] %[output:202025a4] %[output:682a1171] %[output:3a24b888] %[output:47ce84b0] %[output:1b93d727] %[output:335e57fa] %[output:2452a7fb] %[output:764ef4ba] %[output:59beae70] %[output:804b0e82] %[output:5daf201a] %[output:8d76d34b] %[output:9b2ec8f7] %[output:3d79d93d] %[output:2fcb32d7] %[output:20fb71ce] %[output:7dd4258d] %[output:5af6b794] %[output:1b0765c8] %[output:3d10d268] %[output:8dc15161] %[output:3e4d24a5] %[output:949f7f4c] %[output:5ab15f68] %[output:21fd2028] %[output:8c2552b9] %[output:7866ac24] %[output:7348f770] %[output:4a88527e] %[output:31127ad6] %[output:8f5f6cd7] %[output:04f9fbb6] %[output:11fd5494] %[output:1f7d5b3c] %[output:515b01b6] %[output:0e5ed95d] %[output:8593dec4] %[output:044a0c3f] %[output:404cb5f4] %[output:335af21f] %[output:707f3b4c] %[output:22d06d92] %[output:436cc9ca] %[output:43fdb9ed] %[output:8ff8636f] %[output:81084beb] %[output:24176cf9] %[output:125bdce3] %[output:611cc880] %[output:7a7c2bb9] %[output:5ef4139d] %[output:44107a08] %[output:85922ebf] %[output:20c8a213] %[output:28030707] %[output:2259fed6] %[output:185722eb] %[output:96ab9118] %[output:92591168] %[output:77934d24] %[output:2f20d0c0] %[output:227653a6] %[output:451d46dd] %[output:8dea331a] %[output:4df0008d] %[output:97a1b478] %[output:69542f99] %[output:8f49e488] %[output:0cadc856] %[output:8f5016b8] %[output:9ebe2ae2] %[output:6bd20915] %[output:9a601da9] %[output:55eb7a0e] %[output:01045405] %[output:37f41a9c] %[output:4338fa29] %[output:84e9c6aa] %[output:4552d605] %[output:7ed69ab3] %[output:3a400476] %[output:0cb36107] %[output:7d9d983b] %[output:3cd231a8] %[output:11a09786] %[output:404f7317] %[output:8ada2ad1] %[output:63ed2dfa] %[output:5c1d11e1] %[output:20e037a8] %[output:44860991] %[output:6b8b6016] %[output:20a41a27] %[output:1cd31c42] %[output:07e4ed86] %[output:13921e26] %[output:6cbe00d2] %[output:40fe822f] %[output:144d024c] %[output:70effe05] %[output:7d58aaf3] %[output:11df0a6c] %[output:2d83e7ff] %[output:58ac33f1] %[output:0762c4da] %[output:26cb82e0] %[output:13cdfb40] %[output:390950ef] %[output:2912060a] %[output:3d14f584] %[output:405c604f] %[output:4032b353] %[output:4fa7a25b] %[output:2b855b80] %[output:53e7d15f] %[output:5ec2e67e] %[output:003aa0b8] %[output:13b9e6aa] %[output:283c6725] %[output:7e048f63] %[output:18e62070] %[output:1af896cb] %[output:834516a0] %[output:51ef09ce] %[output:029f99bb] %[output:0e0bd2e7] %[output:8f7b2aff] %[output:916fc741] %[output:9a216205] %[output:4ffdbbca] %[output:54c29f01] %[output:293c0fe3] %[output:5d67851c] %[output:6a4ae9ce] %[output:3d8d5523] %[output:47351e70] %[output:8da8479d] %[output:3a3e83fc] %[output:500cd33a] %[output:8a2a0dbd] %[output:6b901616] %[output:426f067b] %[output:123a3273] %[output:7bde413c] %[output:00aa4729] %[output:64ea656a] %[output:170bbd7d] %[output:13f6f81e] %[output:0dd60b59] %[output:06d71b32] %[output:6b93a191] %[output:16c03c5e] %[output:0939d0d1] %[output:7fa9d06e] %[output:57c407d5] %[output:22096909] %[output:832ad6a2] %[output:6e7066dc] %[output:9f1e73fd] %[output:5e154e2e] %[output:9dc5e0c8] %[output:5b4dddb9] %[output:7e9bde79] %[output:22f0e9e2] %[output:89fcb61f] %[output:78264476] %[output:4e475476] %[output:3426c9e3] %[output:73878c87] %[output:1302d4ec] %[output:6ee7f019] %[output:10a91027] %[output:1131e20c] %[output:36cc0546] %[output:3b83bcf6] %[output:8bf0fa99] %[output:565102f9] %[output:0231ebcc] %[output:1cb0c2ca] %[output:6dff560a] %[output:3753a8c6] %[output:54335c98] %[output:7f499ecb] %[output:193b09b0] %[output:86281533] %[output:0145ed57] %[output:74794c1f] %[output:526cc4a4] %[output:197d80cf] %[output:3384609c] %[output:8af491aa] %[output:5db9c53b] %[output:49aa337c] %[output:871e95fd] %[output:21a5997a] %[output:8504ebb1] %[output:5f0c043b] %[output:2c46c402] %[output:4f3a9654] %[output:32875c2b] %[output:8a9f24f7] %[output:398e665a] %[output:7f0cf9a9] %[output:8a92ad6c] %[output:7fba6f7f] %[output:2bbd7f2b] %[output:0ab61f85] %[output:105cdd2f] %[output:7e03e276] %[output:82c86b80] %[output:5908a2de] %[output:6451264a] %[output:8da1554e] %[output:942be068] %[output:3750789b] %[output:0226794e] %[output:40cbc811] %[output:55145175] %[output:7114c91c] %[output:5ad03dc2] %[output:7a04deb5] %[output:7561bb21] %[output:38df06b3] %[output:04c6cc1e] %[output:59cdde8c] %[output:8db969f7] %[output:03aa8b1d] %[output:604fa63f] %[output:44d9363c] %[output:17cf9d66] %[output:2197c516] %[output:1f6cddd2] %[output:43677fa4] %[output:47360468] %[output:8ab440b3] %[output:0e032535] %[output:360f628f] %[output:3196e42f] %[output:5a409d35] %[output:766bdfa0] %[output:32ad3be7] %[output:67aaec13] %[output:21043f07] %[output:269dd0ad] %[output:0a755934] %[output:885f4191] %[output:01516819] %[output:66424d7a] %[output:9a5651e0] %[output:2c16f4d2] %[output:43574a1e] %[output:2d170537] %[output:269b68ce] %[output:67e9a9bf] %[output:58427691] %[output:96e763d8] %[output:364eac83] %[output:12c4fd54] %[output:3c706868] %[output:3f219f2b] %[output:74a9f862] %[output:3ad4bc92] %[output:268a96d5] %[output:7f2b5e1d] %[output:2b361cb5] %[output:6fd19446] %[output:65e997b4] %[output:631614df] %[output:2c1bb521] %[output:85402f60] %[output:484152f1] %[output:9393079a] %[output:10be65f5] %[output:0ae5ddda] %[output:9902be2d] %[output:4caf4487] %[output:2f890ad2] %[output:8d531fa2] %[output:2c35f629] %[output:35c567ed] %[output:431e29e7] %[output:8d94946e] %[output:48c4f104] %[output:5c365f69] %[output:71ea00ee] %[output:0fb5e2b8] %[output:5f4cdce3] %[output:2d3ad30f] %[output:1c71a002] %[output:9aa73e54] %[output:9b388da7] %[output:490e9777] %[output:35819416] %[output:298bd542] %[output:1b7e458d] %[output:2e1e36a5] %[output:1b2f2692] %[output:7adccd10] %[output:3ac1a829] %[output:657485a2] %[output:5c24c4f2] %[output:40994209] %[output:958078bb] %[output:2c8a80f3] %[output:0d76f144] %[output:662e877f] %[output:75536170] %[output:2245e65d] %[output:71136f15] %[output:66015d65] %[output:88a72a8d] %[output:263bf04a] %[output:14578637] %[output:7b3e3450] %[output:0315d6bd] %[output:68dd3624] %[output:4644af50] %[output:868713a8] %[output:49f9b334] %[output:5e5e9f5d] %[output:1211c64f] %[output:3abbbd5a] %[output:705f0ba7] %[output:37f283e6] %[output:841b7545] %[output:5b94d4fb] %[output:3a2f7533] %[output:986e44d4] %[output:3dd90dc4] %[output:5773a26c] %[output:384b0936] %[output:3aef273f] %[output:57d1a236] %[output:9ffe48b9] %[output:7d9dcebf] %[output:1ff24d29] %[output:39104e4e] %[output:39f37174] %[output:19df4136] %[output:7b4acf07] %[output:5d63cb7f] %[output:66641eda] %[output:4f62a34e] %[output:62f9b160] %[output:7e399add] %[output:9cd5c4c7] %[output:09c34be3] %[output:590f0e08] %[output:06a0dbc7] %[output:868a019c] %[output:288dc84d] %[output:651f8d84] %[output:9f188243] %[output:2dda2c32] %[output:7c827be7] %[output:84c6948e] %[output:6ce6199b] %[output:9ab496f0] %[output:106c1849] %[output:0c0b9532] %[output:8eb6c0c3] %[output:978fbe13] %[output:44e63ea3] %[output:47329d18] %[output:20983280] %[output:598a18f8] %[output:5e198369] %[output:096e4178] %[output:761859b0] %[output:5879d413] %[output:0684b5f3] %[output:3f396eb9] %[output:0411fba4] %[output:6284c3e0] %[output:2d6f8393] %[output:043f699a] %[output:04f70217] %[output:6abb8865] %[output:0614f41d] %[output:248f905b] %[output:5a16052a] %[output:8458e8b1] %[output:266764be] %[output:84d95b2d] %[output:15e3319a] %[output:0181af2c] %[output:64955de5] %[output:468c3142] %[output:8280a161] %[output:5e5d6bff] %[output:5b8110b5] %[output:0b65bad8] %[output:5a26d20b] %[output:49cccf24] %[output:9724f269] %[output:2c762a4c] %[output:72f8a101] %[output:9041b75e] %[output:7698677a] %[output:3e6aa175] %[output:5956a609] %[output:2f1ad408] %[output:904174e3] %[output:61801f2b] %[output:5f22683c] %[output:77a75075] %[output:4f6a3ed0] %[output:5cc5b89b] %[output:36f6ca17] %[output:7e5c19f2] %[output:0058c71e] %[output:48c1b2e3] %[output:8347a9e3] %[output:45668607] %[output:2e2a8db2] %[output:43f27157] %[output:849f6db8] %[output:199b24db] %[output:95e51498] %[output:530a4e15] %[output:81b4a461] %[output:83263e76] %[output:7d4f373a] %[output:9c8bbe9d] %[output:8d411675] %[output:9673bd5d] %[output:295c53f3] %[output:1955e1c5] %[output:86c8352d] %[output:3cb8c61e] %[output:590598d1] %[output:79c1c6b0] %[output:43c867b3] %[output:363c7ea9] %[output:950136a6] %[output:07dacc93] %[output:9143408d] %[output:9fae1e7e] %[output:73078948] %[output:422914af] %[output:8b1d1a12] %[output:98823ca3] %[output:3ab5bc3e] %[output:7a92033a] %[output:0bc506f1] %[output:174cc008] %[output:36a8e301] %[output:5edfa16c] %[output:2d7ec5e4] %[output:712fbd5e] %[output:20266a4c] %[output:86ac97e4] %[output:2f5dd8df] %[output:65279907] %[output:2fdfcd67] %[output:3f2a9e4f] %[output:28d971ce] %[output:1f16d540] %[output:491b7153] %[output:9b01c973] %[output:6bf45461] %[output:4f2d2b9a] %[output:38f4b3c0] %[output:84e1832c] %[output:35107c89] %[output:9fe0cdc8] %[output:86f6ad8f] %[output:4be950d3] %[output:8baded26] %[output:908b7497] %[output:53442daf] %[output:90118232] %[output:310f1237] %[output:60c9e291] %[output:27387850] %[output:68dcc750] %[output:4b20dd77] %[output:9bf85d19]
        fprintf('  m_Cs = %s | %s\n', model.mets{rec.m_Cs}, model.metNames{rec.m_Cs}); %[output:3e7be64d] %[output:5a089de4] %[output:4a6eeeeb] %[output:6e7ea9f2] %[output:98145ac0] %[output:6f3eefd6] %[output:32a6c7bb] %[output:4bb0aeff] %[output:226c96fb] %[output:41a4e91c] %[output:93d77437] %[output:7af96c1a] %[output:267c8192] %[output:7fb6df0b] %[output:6a865d87] %[output:282f98c0] %[output:48b1f511] %[output:2db086c5] %[output:2245b1fa] %[output:1ff5b3ea] %[output:7bd12cae] %[output:1535371d] %[output:3a957dab] %[output:4f79a254] %[output:27f87fd4] %[output:9601c806] %[output:2751c774] %[output:83d40db2] %[output:3d95d754] %[output:11d0bc76] %[output:90f5d5c9] %[output:0d5a2c42] %[output:87fea622] %[output:838af062] %[output:35fc1789] %[output:76f526b1] %[output:05198aaf] %[output:0c626fe5] %[output:148f3827] %[output:1ce853d4] %[output:11725d22] %[output:69e06d17] %[output:8c685e90] %[output:02c2eb4e] %[output:9c96292e] %[output:4d1f71db] %[output:89edef7c] %[output:1d3f5f17] %[output:7d72e72d] %[output:4ba8096f] %[output:2bfede77] %[output:00dc73b1] %[output:49b554f5] %[output:17305af5] %[output:8ea16393] %[output:6a86667d] %[output:459fffaa] %[output:481bc0a9] %[output:9561625e] %[output:4c4e4284] %[output:6ef21499] %[output:075f2968] %[output:76f04d96] %[output:4ec5a83c] %[output:7d5308ad] %[output:2e247695] %[output:1c90138a] %[output:8d8575a0] %[output:5267237b] %[output:4d9057d8] %[output:21d8045d] %[output:6deb9326] %[output:45466aa0] %[output:37e403bb] %[output:2a15e637] %[output:87f0b6a7] %[output:3511014f] %[output:261f5f29] %[output:45918ad4] %[output:4a09d5d4] %[output:13731fda] %[output:7a29d57b] %[output:1367f115] %[output:3aeb453e] %[output:6d020ae1] %[output:56fc8e54] %[output:948be976] %[output:5aef3b5b] %[output:1d882eb7] %[output:0bdd6800] %[output:71013814] %[output:779b7fcc] %[output:9c75d8dc] %[output:7d475c72] %[output:506f5e20] %[output:111f1f92] %[output:57452e86] %[output:8e83a2d9] %[output:11da2b1f] %[output:5614f277] %[output:1277d0d9] %[output:9f1bdd92] %[output:78f6b265] %[output:97773ed6] %[output:8cc07c19] %[output:1ddf1508] %[output:8f7cf411] %[output:36715477] %[output:69193989] %[output:3bacc4ed] %[output:711560f6] %[output:75ae5515] %[output:6aa2f724] %[output:7b48b4da] %[output:627d0339] %[output:824d7f67] %[output:6961588c] %[output:3f689b34] %[output:85fe7a5f] %[output:4899a452] %[output:2b16f06a] %[output:5a65f208] %[output:4ec43ee0] %[output:40a9b9ef] %[output:525ce569] %[output:92404a45] %[output:7f159dab] %[output:98bac089] %[output:660005c7] %[output:81dc7e06] %[output:22a9111d] %[output:618f745a] %[output:5d04807d] %[output:4286e7b8] %[output:50af68d5] %[output:4878316f] %[output:06fd1d37] %[output:279b2da5] %[output:05e37167] %[output:24a7360c] %[output:6ee28fe1] %[output:213570d2] %[output:4f5898b2] %[output:2241d781] %[output:5e49815a] %[output:1082281b] %[output:2cae2caf] %[output:468c7c8f] %[output:7fba1d13] %[output:8313ad35] %[output:34c6b01d] %[output:535f0f06] %[output:84607e40] %[output:90f07bef] %[output:5c4db51a] %[output:472a0eb1] %[output:6a92e7ae] %[output:54b2fe40] %[output:67202c7c] %[output:48dda778] %[output:927372c7] %[output:3af9c347] %[output:1a47bf67] %[output:00db8793] %[output:0526ce0f] %[output:59b40a9f] %[output:2123174c] %[output:7156d3ed] %[output:016d8d1e] %[output:0f8cf826] %[output:2c05d9b9] %[output:14e5502d] %[output:23092ccd] %[output:8fb3215e] %[output:867d98d9] %[output:0b772f08] %[output:562276ad] %[output:95c17e42] %[output:32e078a3] %[output:63d0dc05] %[output:2a0869b7] %[output:2f6edbdc] %[output:57a78e81] %[output:42edc181] %[output:57609bcc] %[output:1c0af54d] %[output:109cb368] %[output:10a8c6a6] %[output:84b24152] %[output:2c17a8f7] %[output:650935b3] %[output:4270bbc0] %[output:6918d3f2] %[output:40242774] %[output:25298172] %[output:85dd35f6] %[output:8f27a628] %[output:0a1fa2b7] %[output:41970077] %[output:751ff263] %[output:5aa84a6b] %[output:50395dbd] %[output:2fcd94ef] %[output:2631f8f5] %[output:767c3c77] %[output:8f7270de] %[output:20d4fb33] %[output:338c38bc] %[output:8e42e9de] %[output:157011a3] %[output:90e7d3b4] %[output:84a9d554] %[output:84848b80] %[output:449f86d6] %[output:19f705cd] %[output:77c20905] %[output:9b4e4bbe] %[output:297d6391] %[output:9347105f] %[output:4a08f445] %[output:241d70f7] %[output:561d00cf] %[output:0261b8b2] %[output:32fe4291] %[output:27467a5e] %[output:9cd34213] %[output:043644fd] %[output:175105f3] %[output:5fded4d4] %[output:47d4ec2b] %[output:4c92c148] %[output:60659c6d] %[output:0e14c5c2] %[output:1f7c0c9f] %[output:361e2ea9] %[output:8dbf0ebc] %[output:4a725329] %[output:537a6777] %[output:736dca4d] %[output:953373d3] %[output:65570dfd] %[output:049a1148] %[output:58b9a769] %[output:935d4911] %[output:4b313b59] %[output:4601cd9e] %[output:8a62b961] %[output:46c12be6] %[output:429b50b5] %[output:4b498d28] %[output:7c3dd49b] %[output:73e00b0c] %[output:20b4b6d1] %[output:69b111c6] %[output:297c7d35] %[output:4f4888a9] %[output:6c7e01c7] %[output:9cbf7ccd] %[output:0e9ec1ed] %[output:137345a4] %[output:95430e00] %[output:68502e0a] %[output:9900d055] %[output:04633914] %[output:5a6a5b27] %[output:74b09fa8] %[output:6da9730e] %[output:53220182] %[output:4b7e3288] %[output:1f760055] %[output:81b3d865] %[output:76147e6c] %[output:5cd59be2] %[output:830ba580] %[output:746c2cc8] %[output:7932bab8] %[output:43cdcec2] %[output:5de07ba7] %[output:79693faf] %[output:177b0974] %[output:14f956e7] %[output:7b923018] %[output:3de323c6] %[output:6b2fa59c] %[output:7201ee03] %[output:669efd3d] %[output:98aeed7b] %[output:888743e0] %[output:7d57a355] %[output:993c6931] %[output:991be7dc] %[output:66e6e3b2] %[output:0af6b221] %[output:4449cf59] %[output:5d6a8c74] %[output:3f793f4c] %[output:6341270d] %[output:0bfc3e5e] %[output:8839ac44] %[output:222bf5e9] %[output:9400462c] %[output:6915dd6e] %[output:1043f5ad] %[output:96039493] %[output:3d1b0428] %[output:21f89dae] %[output:1b9e9d14] %[output:541863d0] %[output:318355ee] %[output:7f66ef14] %[output:9b569f9c] %[output:4ddbe073] %[output:2e5b9d19] %[output:47519b27] %[output:1b4fc589] %[output:82caaf66] %[output:8cab2090] %[output:7ab46c37] %[output:4717d24d] %[output:0b182de3] %[output:401e9827] %[output:6549e7e4] %[output:34b651d7] %[output:2e634c56] %[output:2cd8cc14] %[output:3d338645] %[output:5fc8009d] %[output:424ca84e] %[output:3ff0463d] %[output:47c2bd0d] %[output:3000cd9b] %[output:5bcf23b0] %[output:09552365] %[output:757db702] %[output:0fd5ff6a] %[output:0d94dcba] %[output:2af1c3bf] %[output:98a342c8] %[output:033128ec] %[output:38134e62] %[output:6b2a982b] %[output:7357db6b] %[output:35391af0] %[output:7fe09a6f] %[output:5d1821c1] %[output:63f01a3c] %[output:7a69d12b] %[output:177c0085] %[output:7c6961f3] %[output:8150ca60] %[output:4362622f] %[output:1eab916f] %[output:59608aa1] %[output:2268fa9f] %[output:66839862] %[output:256173a4] %[output:8c9ad7a0] %[output:2878aa7b] %[output:06c312fd] %[output:11510ea1] %[output:7778fa4b] %[output:7cb31358] %[output:6a9a85f0] %[output:392f5e13] %[output:0931bdad] %[output:8099d5d5] %[output:6e43ae40] %[output:8eabf12f] %[output:052f0b5a] %[output:4bba9f25] %[output:2949901d] %[output:00deabcb] %[output:18bdac07] %[output:469037c9] %[output:60fc8806] %[output:6b08695d] %[output:58606d62] %[output:68f68bf5] %[output:4561ce6b] %[output:622408b6] %[output:9f09b0ab] %[output:68d29d23] %[output:59deb958] %[output:1be086de] %[output:6048f989] %[output:6f42dac9] %[output:8ff2fb1a] %[output:416a46e4] %[output:2eca8274] %[output:932c13d3] %[output:940de893] %[output:6b622c14] %[output:3e534b3c] %[output:4207ffc4] %[output:228a9b3d] %[output:92d3c7d6] %[output:4e7d4204] %[output:791d3519] %[output:0e3974be] %[output:6eec7c0b] %[output:4012bb70] %[output:5d79d1f0] %[output:0a66dd83] %[output:2ab46e6d] %[output:2b5a186a] %[output:76a991f2] %[output:3e14b8bb] %[output:100e3cf9] %[output:9a25b726] %[output:35e55fa0] %[output:3968f335] %[output:8a499b09] %[output:09fd7090] %[output:64bd7e32] %[output:15a04c13] %[output:17240c2c] %[output:14400900] %[output:6a19ddac] %[output:6daf0749] %[output:30153d8b] %[output:56822d5a] %[output:2ae93401] %[output:8e070853] %[output:696b0094] %[output:0aaf942e] %[output:798d14ec] %[output:524928c7] %[output:3457c9b7] %[output:360c42be] %[output:48a09627] %[output:9695a6a1] %[output:7a4f7c5e] %[output:81f608f7] %[output:3584ae6b] %[output:2433509d] %[output:566d1878] %[output:9c3012c0] %[output:134c7722] %[output:97855979] %[output:9ba08091] %[output:968dacab] %[output:7af00410] %[output:3b24bcec] %[output:37d24675] %[output:2be9bda1] %[output:34c1bf65] %[output:4b14734a] %[output:3e91ab8f] %[output:392aace5] %[output:2d771a8f] %[output:1852dc51] %[output:035436b8] %[output:50edf1c4] %[output:17e50f4e] %[output:10f0dedd] %[output:7f298104] %[output:5127a442] %[output:83430433] %[output:2e884947] %[output:4aa657e1] %[output:819230de] %[output:5656626e] %[output:44005adb] %[output:330726ad] %[output:0be398f7] %[output:69652af1] %[output:1670a534] %[output:50e2ae76] %[output:9f1d1085] %[output:6de593d6] %[output:64835c8f] %[output:70fb325d] %[output:2517942f] %[output:643f4110] %[output:880a7d57] %[output:80707306] %[output:502c8bcb] %[output:2155a040] %[output:3e469d47] %[output:5d2717a6] %[output:2dd4e662] %[output:527afa14] %[output:3028218a] %[output:1df2df2c] %[output:62ddff27] %[output:573843d4] %[output:7bea5766] %[output:567e9ce5] %[output:2e61e727] %[output:02f53d13] %[output:95bf4f63] %[output:9e99586e] %[output:2af296de] %[output:3aebab70] %[output:8011f8e6] %[output:0ef6f547] %[output:2faaeec1] %[output:2edaf821] %[output:62e6790b]
        fprintf('  r_Cp = %s | %s\n', model.rxns{rec.r_Cp}, model.rxnNames{rec.r_Cp}); %[output:5d868d56] %[output:2c37fb9b] %[output:10c86049] %[output:83d93a25] %[output:74c23a4c] %[output:6bb59d60] %[output:805e182c] %[output:1c70dfc0] %[output:46681df3] %[output:953cad0b] %[output:1c9255a0] %[output:82799960] %[output:23016969] %[output:0d8aa916] %[output:729c198d] %[output:623f96d6] %[output:0a31d30b] %[output:14914d87] %[output:58c0bc82] %[output:9c66502f] %[output:6faafaf8] %[output:2bd8b4ba] %[output:3da6859f] %[output:98877294] %[output:52808868] %[output:8139b4a3] %[output:9044c80c] %[output:6021a4ac] %[output:24302ac7] %[output:9fc0fb5b] %[output:9a098457] %[output:1a38401f] %[output:36d533f1] %[output:93fdc256] %[output:9c0b327a] %[output:7fc28e82] %[output:8f9eae16] %[output:0642c059] %[output:6fa907a1] %[output:3b64aad0] %[output:01800888] %[output:8ed571c6] %[output:323bf1a7] %[output:885f09e0] %[output:5e473573] %[output:78dce281] %[output:2057fefb] %[output:63eafe62] %[output:7bfa8b5f] %[output:836a99ee] %[output:142c1056] %[output:9609eb5d] %[output:53818287] %[output:8eff27c0] %[output:9d5523b8] %[output:02915393] %[output:1f2b4757] %[output:3c6a4f64] %[output:6c58a366] %[output:77e9283f] %[output:08069b8d] %[output:4c0fd898] %[output:568d3ce9] %[output:2d009b09] %[output:171832e5] %[output:497ff4da] %[output:8b589b0b] %[output:41e5354d] %[output:05384ee2] %[output:8f1d328d] %[output:805128fa] %[output:86297da8] %[output:5b872a36] %[output:4a6ef407] %[output:1d02f969] %[output:10b608d8] %[output:40837d18] %[output:22975d00] %[output:46beff76] %[output:3deeec57] %[output:16160485] %[output:0753b3cb] %[output:69d08b6c] %[output:3ba1e1c0] %[output:48f95032] %[output:70bf81f7] %[output:70f55f0b] %[output:620c2b5d] %[output:96e33254] %[output:8f4c371c] %[output:6426d08d] %[output:67fabede] %[output:3a10161f] %[output:4045b817] %[output:46f0905b] %[output:3a62318f] %[output:68f34b90] %[output:3762adac] %[output:5b4e8c53] %[output:39d02512] %[output:2ddcca1d] %[output:776d4ce8] %[output:640a330d] %[output:08f7fcb9] %[output:4f9ea2dd] %[output:9483c541] %[output:65e450f1] %[output:9d8ae395] %[output:5adade27] %[output:83098b0b] %[output:96b2c83f] %[output:967f5b47] %[output:0ef05ada] %[output:84e85b97] %[output:0c1705f1] %[output:1ababc29] %[output:07233b92] %[output:6260723e] %[output:48b8fb45] %[output:338108f2] %[output:759b9a39] %[output:64f7a24b] %[output:051fbd17] %[output:588fabf7] %[output:0edb288f] %[output:4fa7efa9] %[output:47104afc] %[output:94b6bc00] %[output:2d60fbaa] %[output:4a396962] %[output:826a8bd0] %[output:1f300e9e] %[output:406573cf] %[output:51d443e6] %[output:20bb0d78] %[output:3632c1f9] %[output:71d9f7c6] %[output:503796a3] %[output:503fd89f] %[output:9a4e29ef] %[output:7869d09c] %[output:91e4b6c4] %[output:534dfec2] %[output:0a0ba531] %[output:85713784] %[output:88f60adf] %[output:6d9f9735] %[output:84a0fa97] %[output:79ad2897] %[output:16276320] %[output:7768d038] %[output:21a9c03d] %[output:322344aa] %[output:77943dd0] %[output:281739ac] %[output:917d2c2b] %[output:565ef731] %[output:4168022c] %[output:1396148c] %[output:12d05fe8] %[output:7c96d9b1] %[output:481fa8bb] %[output:8d92b5fd] %[output:04fe9ba5] %[output:7069e94d] %[output:249e8863] %[output:009a76c1] %[output:1c49c597] %[output:88cae9d8] %[output:0aa3c219] %[output:53b312f9] %[output:7358f93f] %[output:26b621c6] %[output:7d72769f] %[output:47877371] %[output:3ebbbf2f] %[output:37abd0f4] %[output:63982503] %[output:5515e694] %[output:23a5c94d] %[output:73085a5a] %[output:9b0ac543] %[output:49731d86] %[output:9e2fa71e] %[output:554102e7] %[output:966ba4f0] %[output:1924e140] %[output:7cefcbd3] %[output:059c16d7] %[output:94de5bb2] %[output:9734b4b6] %[output:7bfe67c5] %[output:990d89cf] %[output:2e871347] %[output:5afc5320] %[output:47f50b11] %[output:86bdbe83] %[output:35420804] %[output:699e62c0] %[output:30f6ed94] %[output:1bfff1ee] %[output:4658e8b2] %[output:8e05d32b] %[output:81c24b18] %[output:08d42a6c] %[output:006b7d49] %[output:47311a3f] %[output:4b7c3e6f] %[output:90e47644] %[output:073e934b] %[output:9fcc6df9] %[output:29392050] %[output:782d8fe7] %[output:27371426] %[output:0347b85c] %[output:3d7fb2f7] %[output:25ad883f] %[output:871c31af] %[output:6c4b1e3c] %[output:53e94600] %[output:5447d571] %[output:9a97b80f] %[output:93a03a5a] %[output:952247ae] %[output:1477760e] %[output:7df9862d] %[output:1bb33e93] %[output:320a0002] %[output:1e84c977] %[output:0d64717c] %[output:0641c1b4] %[output:5101d99b] %[output:3b7f03f2] %[output:5fd98075] %[output:5447bf87] %[output:9397cf4d] %[output:3815c0ac] %[output:2e6f265f] %[output:389fa0e7] %[output:2328d1b4] %[output:2548f1b6] %[output:79ca8bff] %[output:107433aa] %[output:154c1ba8] %[output:8b548875] %[output:75c0eeb8] %[output:3597a13b] %[output:7559cdbe] %[output:95d70951] %[output:75d9a855] %[output:753f5a2e] %[output:6f12363c] %[output:4357f316] %[output:727e7828] %[output:31eff13c] %[output:596b991a] %[output:6636682e] %[output:6488ad7a] %[output:4daa28e6] %[output:441bdd92] %[output:581b3617] %[output:3bae1260] %[output:18a16659] %[output:70e3cfb1] %[output:229c86d4] %[output:677bcfc6] %[output:57b6daf3] %[output:4dbf1b8f] %[output:1ba823a0] %[output:7133743f] %[output:77a7a494] %[output:7bb69dab] %[output:0da47a7e] %[output:508e9a4c] %[output:8ce17f1d] %[output:9a42a3e1] %[output:9d50b394] %[output:7d9e7f10] %[output:28ff98c7] %[output:8beb7463] %[output:14ca1706] %[output:1e45f293] %[output:77c3a373] %[output:37f6d7bc] %[output:728f58fe] %[output:3ed5f22b] %[output:7e5698df] %[output:9853daa4] %[output:259103e0] %[output:9b05435f] %[output:1a62bdd3] %[output:239de070] %[output:84737b11] %[output:25df58dd] %[output:7994c4da] %[output:37e93bf2] %[output:08bcf4c9] %[output:770b3df2] %[output:48f55f49] %[output:00c2472e] %[output:0fd8ded1] %[output:2262280a] %[output:3cf17ebb] %[output:1d8af207] %[output:4329d3bc] %[output:24039353] %[output:50701ce3] %[output:5f3227a5] %[output:1a9cd918] %[output:84e180d3] %[output:60bf9001] %[output:585bacb8] %[output:79366e7b] %[output:53ef2afa] %[output:8e60c6d6] %[output:1a1439b5] %[output:77447882] %[output:972acfd4] %[output:27e6e061] %[output:2632669c] %[output:900138e6] %[output:8c1ffc83] %[output:9a2cbfed] %[output:63e2288f] %[output:27540b19] %[output:1e9b57f6] %[output:2462d50d] %[output:9b99c78e] %[output:06d84711] %[output:4458e01a] %[output:2e3d8f95] %[output:0ef0eed4] %[output:290b8489] %[output:622c9bcd] %[output:2c7c8354] %[output:1ec704f7] %[output:212b8189] %[output:99285de2] %[output:05e37099] %[output:02d1aea6] %[output:42733788] %[output:2ad2d9d7] %[output:2acc8372] %[output:7c87238f] %[output:673b8f00] %[output:3abb4eb2] %[output:112f9527] %[output:1126bdb0] %[output:3216c2af] %[output:2f053480] %[output:5502b76e] %[output:968c1a89] %[output:41dbb4df] %[output:82847765] %[output:3d29f4c3] %[output:7629875f] %[output:3b3381fd] %[output:4a6bb5eb] %[output:4a28c77e] %[output:36379202] %[output:08857e07] %[output:52f9e7dd] %[output:2a8cb659] %[output:6a5c5964] %[output:6e74992f] %[output:60ecbd04] %[output:2258bc7a] %[output:935cdf45] %[output:4698e948] %[output:74d66d85] %[output:120fa9c6] %[output:004c6b37] %[output:7161fea0] %[output:94e56575] %[output:2a2ea5f7] %[output:408783f0] %[output:7548aefb] %[output:02e5d532] %[output:8507bd02] %[output:43f16d75] %[output:13e84767] %[output:88f5a937] %[output:2a31f5dc] %[output:29fdabc3] %[output:68a7db27] %[output:09bf4a41] %[output:204f654b] %[output:2618260c] %[output:164db8b8] %[output:669fda31] %[output:3ef2e1eb] %[output:581cb38d] %[output:28afc39a] %[output:78a227fc] %[output:11807403] %[output:90c72559] %[output:0f54cee8] %[output:919a356d] %[output:8a7f913b] %[output:6b1f8df8] %[output:8297abdc] %[output:4a395970] %[output:6383fb17] %[output:828b775a] %[output:082e92cf] %[output:17b9b2f5] %[output:30eedcb2] %[output:586740a8] %[output:273c3548] %[output:56ebf098] %[output:0227cfa0] %[output:18367dbd] %[output:4b2d4ddb] %[output:19458ab2] %[output:15a09977] %[output:6db91702] %[output:4f1ac81a] %[output:465c56b0] %[output:855bfac5] %[output:1c449c5f] %[output:430b36f2] %[output:9ddd4e05] %[output:3795ec72] %[output:3b110b09] %[output:7e8da264] %[output:4d9235a3] %[output:3746709a] %[output:10c3b6fe] %[output:781f3705] %[output:325fd73f] %[output:2fcea368] %[output:59ce63d2] %[output:97514f4c] %[output:6a2a3885] %[output:9a0ba7e5] %[output:6ec1226a] %[output:7b34277a] %[output:1f02be7c] %[output:4970f1de] %[output:587a456e] %[output:8775cf0d] %[output:0af60b05] %[output:9f988f96] %[output:5aa64547] %[output:1684026f] %[output:19ebba3d] %[output:23934dd4] %[output:7f835af4] %[output:8069be55] %[output:19759201] %[output:2780284d] %[output:09a52fec] %[output:2d2c0617] %[output:52911ca6] %[output:8fc9992f] %[output:3ac70cf6] %[output:960be5cc] %[output:17334f2a] %[output:7e6f88bd] %[output:7aec911d] %[output:15a401e0] %[output:85a69172] %[output:41868ece] %[output:8c79d1bb] %[output:88cb7b6b] %[output:9728b0d8] %[output:7f9ea80b] %[output:378e85db] %[output:5cab25ad] %[output:153f9d57] %[output:01f8cfe2] %[output:376a4837] %[output:0d245b0e] %[output:258cbbf0] %[output:441f962a] %[output:1c1766d5] %[output:5c92e797] %[output:028ecf05] %[output:64b08b37] %[output:4abf8608] %[output:6f834bae] %[output:65beeb3a] %[output:3f847acd] %[output:418d10ef] %[output:33ef2632] %[output:2964fc92] %[output:032865f6] %[output:4bf1d927] %[output:4de659ee] %[output:757be0bb] %[output:2c3a67c0] %[output:4a020190] %[output:94a5a67a] %[output:4903feb4] %[output:798e14fe] %[output:43275654] %[output:00536189] %[output:8a525ca7] %[output:2cfe2231]
        fprintf('  r_Cg = %s | %s\n', model.rxns{rec.r_Cg}, model.rxnNames{rec.r_Cg}); %[output:45ef7353] %[output:80c13a38] %[output:54db9230] %[output:50d800b8] %[output:1af78b63] %[output:50137529] %[output:119a2049] %[output:3c95432a] %[output:4144e535] %[output:047cd987] %[output:18b69a13] %[output:169bf3de] %[output:64cd8aac] %[output:4dc7bab7] %[output:036ca5cd] %[output:6f172576] %[output:015fe111] %[output:8b8a800f] %[output:1e6b1a11] %[output:3857d23c] %[output:283577aa] %[output:8fcda8c3] %[output:0d7ecf9f] %[output:6e9cd27f] %[output:3e9fbb43] %[output:0be9b620] %[output:92bfc140] %[output:78c7940c] %[output:6f0fb172] %[output:9d721458] %[output:9af451b7] %[output:000386ba] %[output:98977b75] %[output:5adc6733] %[output:24daad86] %[output:2fbac2eb] %[output:43cab3d2] %[output:2ad06965] %[output:781b5302] %[output:3a0743e8] %[output:2c2000f8] %[output:090c880e] %[output:0589bfbf] %[output:6020ebfc] %[output:0e70b173] %[output:33b9cdf0] %[output:22027e7d] %[output:2ca0cc5b] %[output:30a3a05b] %[output:300d4473] %[output:0ae79aa6] %[output:48f48955] %[output:79d2818d] %[output:2db01ad6] %[output:040ae057] %[output:9ab3128e] %[output:30802ac1] %[output:05949f02] %[output:989e660d] %[output:9e9bbdaa] %[output:17d08faa] %[output:09ac476b] %[output:38ed6bec] %[output:2fe6f433] %[output:1a6362fe] %[output:833b07a3] %[output:937a73df] %[output:7e6c0b84] %[output:0c53e0cf] %[output:41da4fe7] %[output:2c5d8580] %[output:0d969fcd] %[output:1e6609c6] %[output:12f698ad] %[output:5c0fd54b] %[output:3237ff38] %[output:95d957d1] %[output:7183ac6e] %[output:3cf5fc30] %[output:0e10988b] %[output:7c472d58] %[output:87826620] %[output:8f837eea] %[output:877d76c4] %[output:534174dd] %[output:1df0e12c] %[output:964cf3e3] %[output:983fd48c] %[output:258a5fe9] %[output:1b123b4c] %[output:692f98e2] %[output:17954f18] %[output:50bdabc9] %[output:5411a8ed] %[output:037dd31b] %[output:068b4afe] %[output:652ea116] %[output:30d9e91c] %[output:3ce5e57e] %[output:42df7012] %[output:2ffd719a] %[output:85213e0f] %[output:54a2a669] %[output:26afedb1] %[output:7e116d47] %[output:0f696751] %[output:3db83209] %[output:318d2f2b] %[output:9fa7a0e1] %[output:95e4090a] %[output:093e61b1] %[output:388d1c7a] %[output:35b312f7] %[output:80181e6c] %[output:9556506a] %[output:181b4f92] %[output:3acfa02f] %[output:485abbf1] %[output:80ef7498] %[output:7e3b8a9a] %[output:50de5680] %[output:7b71028e] %[output:2466f9be] %[output:3c88cfbd] %[output:167f55af] %[output:2ccf9583] %[output:4c05dd78] %[output:6ec14811] %[output:6351792b] %[output:07ed936c] %[output:8948bd1b] %[output:51d3bd4d] %[output:08d3b8c0] %[output:8242ac53] %[output:4d4bcf50] %[output:6d15cff7] %[output:9f560ab1] %[output:027fc281] %[output:34c50c3d] %[output:96642862] %[output:63e78da1] %[output:42ec25fa] %[output:177d0406] %[output:0708eca4] %[output:81699372] %[output:607b4d8e] %[output:0b93f9a6] %[output:95b85f43] %[output:43c22e8a] %[output:6af7a312] %[output:8afdbdd9] %[output:2c33d748] %[output:259d3a7e] %[output:10e1cd2a] %[output:995cd326] %[output:4bd64dd4] %[output:3fefacc5] %[output:83f16f97] %[output:279f892a] %[output:7ee089fe] %[output:23c33b76] %[output:57bf3c32] %[output:830b1451] %[output:2f6bb058] %[output:8fa81526] %[output:5f8ac1e9] %[output:8de37d26] %[output:0c8e7bc7] %[output:32ca5770] %[output:2099f60f] %[output:2758241f] %[output:1fe6cd0e] %[output:83c33ece] %[output:20500f90] %[output:13d193c0] %[output:036d74a3] %[output:468ed43f] %[output:99ccb90f] %[output:44b70560] %[output:28df34b2] %[output:1a981335] %[output:7a4d1176] %[output:0b630542] %[output:84107c08] %[output:3428de36] %[output:286cbbb7] %[output:53f84a0f] %[output:5de10db1] %[output:136c1571] %[output:65f83dfc] %[output:87ffb88d] %[output:33d08b1e] %[output:971e0ee6] %[output:6904afd8] %[output:4ea6d06a] %[output:9215f574] %[output:70b8bb3a] %[output:102c9e2c] %[output:166aab94] %[output:1041882b] %[output:290c9410] %[output:3b9bef58] %[output:16fab257] %[output:7845a723] %[output:020c01c1] %[output:3b0bedae] %[output:10735288] %[output:7ed63650] %[output:391d768b] %[output:5fd4b3ef] %[output:17f44f21] %[output:57cb4f7d] %[output:0d10e97c] %[output:6b430ea7] %[output:8f773133] %[output:4317d4ff] %[output:5f7623ae] %[output:8bbaeaa7] %[output:5f8ad696] %[output:2653a610] %[output:47427404] %[output:08f4db37] %[output:2c595fd3] %[output:579ee851] %[output:64d222b7] %[output:66ea54a2] %[output:8eabf50f] %[output:83b542b6] %[output:3f51d8ef] %[output:1ba6f504] %[output:32f0b3d5] %[output:00588c55] %[output:9dfa21b8] %[output:168b08a2] %[output:0007af8a] %[output:7fe28aa9] %[output:63d63ac8] %[output:12c51791] %[output:6768c6d0] %[output:1fab0bd8] %[output:1a8aa862] %[output:0277c903] %[output:0ab1a0ed] %[output:0e8ad17a] %[output:75c6e2e3] %[output:519fda88] %[output:4581891b] %[output:45dc3a0e] %[output:79cbb026] %[output:600ae140] %[output:97569e86] %[output:486fa5f8] %[output:6a7fb119] %[output:4c054e60] %[output:8be7ff00] %[output:55a75d8e] %[output:6ad7645b] %[output:74cb3ea5] %[output:3adea439] %[output:30bf1acc] %[output:55ce52e9] %[output:63e060bd] %[output:292860a7] %[output:829fa097] %[output:14dcaa36] %[output:22149dcb] %[output:888d01d8] %[output:2439bf58] %[output:82773c5e] %[output:57ef70b9] %[output:22b35873] %[output:310faf1d] %[output:5740fd4a] %[output:7bfceca4] %[output:41ed669b] %[output:24871ea0] %[output:61b1cd35] %[output:212d5369] %[output:9dbd0074] %[output:485e16f8] %[output:053bf4a2] %[output:25036e7c] %[output:8cec23a0] %[output:8fb98f52] %[output:84b91028] %[output:718de745] %[output:6d34af46] %[output:6687740c] %[output:4a51aef5] %[output:695eb2da] %[output:65cf6bcb] %[output:40c96c1c] %[output:2eb3b57f] %[output:4e5104f2] %[output:0a989154] %[output:0f324cf3] %[output:41d3b162] %[output:4a498d17] %[output:57fe7f47] %[output:23f500d9] %[output:0c8785ae] %[output:6d57129a] %[output:893431a7] %[output:25925afd] %[output:1f569faa] %[output:35c4906c] %[output:0eda4363] %[output:1906a715] %[output:23b1e420] %[output:4748dcab] %[output:6c0690ec] %[output:59c84627] %[output:0aaf5ac3] %[output:7b6be62f] %[output:4066bd48] %[output:3e9a85b1] %[output:00e4abee] %[output:0c8514f4] %[output:14d597a3] %[output:65282246] %[output:8cab9302] %[output:56969891] %[output:7a5a7375] %[output:6030bbda] %[output:298084ae] %[output:820ba7c9] %[output:3666d0d9] %[output:4f549e01] %[output:5fefd0d6] %[output:9e1a534b] %[output:47c98a4f] %[output:09ff383e] %[output:89fa86fb] %[output:140975fc] %[output:8543b208] %[output:4b1ad3e2] %[output:235cc1f2] %[output:920156cd] %[output:8f4cba31] %[output:57164c77] %[output:8ce864be] %[output:4d6ab2bb] %[output:5051f439] %[output:3db48a68] %[output:8fdbd353] %[output:4440b5c3] %[output:6b200643] %[output:7debd99d] %[output:7f59b140] %[output:5143d211] %[output:763ee520] %[output:25452f2a] %[output:8895254e] %[output:57bdc2e7] %[output:9c0d1e53] %[output:3f591147] %[output:86bb41dd] %[output:21e32ea1] %[output:08ee3074] %[output:56fa33f9] %[output:4c2ecfc1] %[output:184ca560] %[output:1c34d055] %[output:1e43415a] %[output:9537ee5b] %[output:84ca9418] %[output:5fedf2da] %[output:9e92e0ac] %[output:95ea053f] %[output:8115d309] %[output:8f18c823] %[output:839d1460] %[output:8b08cd40] %[output:832fd064] %[output:051f844b] %[output:2b2c052b] %[output:70e562c9] %[output:4301c189] %[output:1fb1506a] %[output:146cdf71] %[output:8f14604d] %[output:6bcd72cf] %[output:80e5921b] %[output:069a9115] %[output:711d701c] %[output:8980b18f] %[output:3b79f521] %[output:6a418232] %[output:3151f38a] %[output:2d3209cb] %[output:47182a1c] %[output:2dc4940b] %[output:44f25d52] %[output:8450ecbb] %[output:3257f413] %[output:03f7ea12] %[output:796baf51] %[output:8cd78cde] %[output:24d2f831] %[output:87b769b9] %[output:78984899] %[output:9d19d009] %[output:9086769d] %[output:849c23e0] %[output:1eae5ec6] %[output:9d7701f4] %[output:1170fc24] %[output:74214652] %[output:565b37e7] %[output:5ca92c1b] %[output:37346636] %[output:161deacc] %[output:98657666] %[output:26a36c10] %[output:0aace77d] %[output:81f95f63] %[output:1045a877] %[output:173d3658] %[output:923ba5a6] %[output:0fd47cf7] %[output:56f30c1a] %[output:37956aad] %[output:099ccb1e] %[output:439c47fc] %[output:9435a3e3] %[output:979e40b5] %[output:76b08a33] %[output:1adecb58] %[output:539c39d0] %[output:7676bec5] %[output:3baeede7] %[output:88ac9df4] %[output:6885054a] %[output:5ca6ce01] %[output:7e0714c4] %[output:405004bc] %[output:64e087e6] %[output:36834873] %[output:68a31698] %[output:91f2fa98] %[output:36e48e48] %[output:85d3febc] %[output:21a84636] %[output:52b35b6f] %[output:56ad8462] %[output:269dbc6e] %[output:76ad496e] %[output:6a59fba6] %[output:8d70a260] %[output:67fb38ad] %[output:537a3514] %[output:397cea1d] %[output:1b08affb] %[output:2fbc05dc] %[output:4f62d1b0] %[output:3398b722] %[output:4f3458c3] %[output:67fb9f16] %[output:3f94d695] %[output:0e28c1f4] %[output:0d35b411] %[output:5d516970] %[output:2d1f4aba] %[output:11ea2aa5] %[output:59aaf96b] %[output:25725416] %[output:56f62df7] %[output:73c80221] %[output:126b9bc7] %[output:206fc9ea] %[output:293723b6] %[output:0a1f48dd] %[output:519a2448] %[output:1a2e6670] %[output:5bbb3fd7] %[output:523f105e] %[output:25114548] %[output:62270d09] %[output:37372cda] %[output:0797dd81] %[output:81282cb9] %[output:07aebe69] %[output:1dde7715] %[output:9f7fd7f5] %[output:53e809aa] %[output:4cf0a0a6] %[output:39acd1bc] %[output:86bce8bb] %[output:236d02fe] %[output:70ac9825] %[output:83aba52a] %[output:4440235a] %[output:93164d40] %[output:2212e609] %[output:46b65a9f] %[output:1a9ad448] %[output:56bc6827] %[output:27e13ecd] %[output:51172f35]
        fprintf('  r_Cr = %s | %s\n', model.rxns{rec.r_Cr}, model.rxnNames{rec.r_Cr}); %[output:4a5015c6] %[output:50e2e668] %[output:8e06b1ca] %[output:3877461e] %[output:4dca5dcc] %[output:46c98303] %[output:21632340] %[output:86b0163c] %[output:8643b86a] %[output:9f715819] %[output:99353125] %[output:260cea50] %[output:8ede528c] %[output:1f7650b8] %[output:9b13ae61] %[output:08bdd325] %[output:133e3197] %[output:77abd6ee] %[output:87c7b0c6] %[output:9dfce923] %[output:5fc3bc7f] %[output:8dfc82dc] %[output:2587ca2d] %[output:8182b134] %[output:46fc9dbd] %[output:40155f3a] %[output:78dc74ef] %[output:311198f1] %[output:2b33e2bf] %[output:14207070] %[output:3c3c921a] %[output:7f39ac57] %[output:630c5390] %[output:2e814bf3] %[output:24b1f45a] %[output:2caadc2c] %[output:945966f7] %[output:01131da1] %[output:858dce96] %[output:8176f9c1] %[output:4b05a013] %[output:213de271] %[output:02fa9b6f] %[output:2fbd38dc] %[output:54967435] %[output:844ecb80] %[output:9b78823b] %[output:0ccee55c] %[output:345a17a5] %[output:7d1e9ed0] %[output:03709de8] %[output:348e8ecf] %[output:9e21e67b] %[output:05554e00] %[output:8da79ab8] %[output:49b9bc67] %[output:29995571] %[output:11f991ee] %[output:16547a72] %[output:1f028dfa] %[output:26c26823] %[output:90b48607] %[output:282f6159] %[output:7a9d3043] %[output:6c5322e0] %[output:83c3cb38] %[output:12a0af94] %[output:2ea0d6cf] %[output:0ac8c280] %[output:8096f4df] %[output:24962699] %[output:6fc579c3] %[output:8fc21a4a] %[output:4c7366fd] %[output:43681c19] %[output:6c68bd5f] %[output:3941f27f] %[output:0d22f074] %[output:26f25dea] %[output:9f1b43f0] %[output:8a6a0f01] %[output:4ca6a836] %[output:7267daf0] %[output:1444d51a] %[output:393e16c7] %[output:9aa967a7] %[output:9b26eae3] %[output:5c22cc6e] %[output:2814787b] %[output:9076c365] %[output:28a33a38] %[output:997d31cf] %[output:2effec92] %[output:107b7cac] %[output:37b96f83] %[output:8d642964] %[output:5725ac1f] %[output:516f6eb8] %[output:372ef0bf] %[output:139be4e6] %[output:5ba9bb8e] %[output:58da7ff4] %[output:012e3820] %[output:58478739] %[output:6ae9e99f] %[output:1c987f05] %[output:4b8fffd9] %[output:8cf1ae62] %[output:3a9f448a] %[output:18be25e1] %[output:47a660b4] %[output:158a4976] %[output:0d70ed5c] %[output:529a38b8] %[output:233acdd1] %[output:69a3341c] %[output:372fac1a] %[output:6fdb5c4f] %[output:5fd70fe7] %[output:72629a22] %[output:13db9f98] %[output:945f99fc] %[output:8bbdcdd6] %[output:913c8fcf] %[output:447f3f3c] %[output:495b52df] %[output:81100de6] %[output:56813306] %[output:2d4d5962] %[output:3daebb5f] %[output:8c074799] %[output:74f1cf1c] %[output:46823e7a] %[output:3867f29e] %[output:38e5acae] %[output:5164cc0f] %[output:957a767b] %[output:0956a49f] %[output:6c344d21] %[output:3027bfca] %[output:29a4bdea] %[output:70ffb2a1] %[output:95a548fb] %[output:780ab89e] %[output:50a9187a] %[output:92bc6272] %[output:7dda08ff] %[output:12267585] %[output:6c88e182] %[output:725ad433] %[output:4748c767] %[output:479717c8] %[output:8adc94b5] %[output:0b7bb37c] %[output:3fa0a74d] %[output:7fff1f9b] %[output:0c0d5a08] %[output:5b14822d] %[output:3639bc89] %[output:35fd4a75] %[output:936d0b6e] %[output:01c3a757] %[output:3d21a863] %[output:97ad8b79] %[output:87c886d6] %[output:25cc4849] %[output:9c9e1e0a] %[output:2bc64d27] %[output:92f9db14] %[output:65041980] %[output:2590c6ac] %[output:9a207f34] %[output:3287329f] %[output:58bfeba3] %[output:52828764] %[output:4ce433a7] %[output:52914fb6] %[output:67bc4cb1] %[output:594ab597] %[output:3bf4a424] %[output:94f5b290] %[output:0761ec6e] %[output:50fc8e4d] %[output:6dcdc2da] %[output:4ada5c0c] %[output:63b93866] %[output:24219de3] %[output:5924b182] %[output:8a44b282] %[output:7870893e] %[output:4439fee4] %[output:4625e3ec] %[output:9d3043a5] %[output:3971af63] %[output:25f8c1e3] %[output:6d760c0c] %[output:8ec6fddc] %[output:67860d09] %[output:52b50af5] %[output:191331fc] %[output:97596f5f] %[output:384bba9e] %[output:1a78529d] %[output:5df81b06] %[output:258cbcc2] %[output:3f5185c7] %[output:8703710c] %[output:606c70d8] %[output:88cf8abd] %[output:40e6ee65] %[output:4bdc7d8c] %[output:09c579c4] %[output:0e68209e] %[output:1454c3df] %[output:1c0eede7] %[output:11b95e0e] %[output:6a7c2897] %[output:9a353137] %[output:50e64157] %[output:333766ca] %[output:07ef2ee7] %[output:5dde3e93] %[output:7543527c] %[output:52c2d37b] %[output:9fead3ec] %[output:10e6482f] %[output:2b6c5bbd] %[output:24556e5b] %[output:306cf3d8] %[output:1099755a] %[output:05fa5d50] %[output:841b5f74] %[output:02dfd280] %[output:3209a1e2] %[output:05d1cfe4] %[output:86448e52] %[output:8a09c5e5] %[output:9b4246af] %[output:5c7b6a5b] %[output:0a223522] %[output:85a4671a] %[output:5132f348] %[output:53f78135] %[output:7054bc5f] %[output:6d5b8a68] %[output:86f92199] %[output:89dedec3] %[output:887b999a] %[output:24b0ee7d] %[output:4fe9691f] %[output:7871e7a8] %[output:313352d1] %[output:6dadfe07] %[output:0cb924db] %[output:0404a6c6] %[output:9e44e5a2] %[output:8f6947db] %[output:799dfef4] %[output:8ec14757] %[output:732589f3] %[output:9244e18e] %[output:87f8bbaf] %[output:1e37b630] %[output:626ac72d] %[output:6f02a1aa] %[output:5740c392] %[output:6da99e38] %[output:7558289f] %[output:76352dde] %[output:1029dbfc] %[output:53bbbb36] %[output:50c3e67c] %[output:8493a6f2] %[output:7ed870c2] %[output:1b85f343] %[output:31642efb] %[output:879ae7b6] %[output:89b0a21b] %[output:039e3cdc] %[output:59395256] %[output:65032e78] %[output:4e40afdd] %[output:65e4fb30] %[output:20a94f4b] %[output:155bdb14] %[output:863e6489] %[output:49574648] %[output:5a23797a] %[output:03fa98a1] %[output:09f7b1d6] %[output:7f0afe56] %[output:60f7fb02] %[output:0442126c] %[output:57298c03] %[output:5bdea67e] %[output:10be0c01] %[output:1c1fa716] %[output:98bcffe6] %[output:28e879be] %[output:8b23d1b7] %[output:4584bb74] %[output:84fc702d] %[output:61433de5] %[output:8458d661] %[output:7d6cbc36] %[output:33254083] %[output:015e541f] %[output:8b55b02c] %[output:53a4f48e] %[output:3ea7d2f7] %[output:8303e6e6] %[output:70e25223] %[output:85761337] %[output:5c2ab548] %[output:46b5436d] %[output:00796948] %[output:3d7c3567] %[output:79158a91] %[output:3f55bbb3] %[output:0868bef2] %[output:74326eae] %[output:582efaa7] %[output:77bcf499] %[output:4b8ba446] %[output:2866f2cd] %[output:7b97c7eb] %[output:11b3f321] %[output:6760a959] %[output:31901581] %[output:7030d79f] %[output:13e203b0] %[output:85492a15] %[output:903da677] %[output:11de238b] %[output:4248617a] %[output:15a6aabb] %[output:922b7bb1] %[output:3c3a274d] %[output:5fb155a7] %[output:7a06bd5e] %[output:3f557202] %[output:18e5e103] %[output:9ff957b8] %[output:69fac5ad] %[output:2f92820f] %[output:62384f37] %[output:9a71e4b1] %[output:56134f30] %[output:38c0a4ab] %[output:12e63ea8] %[output:49e85832] %[output:06567046] %[output:68f0050c] %[output:708f4385] %[output:3788cecf] %[output:4f079df8] %[output:911c203e] %[output:588a2950] %[output:1e0435d9] %[output:2b06f57a] %[output:3ed76987] %[output:08278dd3] %[output:7447593e] %[output:7ab3bdd6] %[output:42569ef8] %[output:85307bdd] %[output:30d0fe70] %[output:5cf30e1a] %[output:900fb9ea] %[output:79e43206] %[output:3d6625eb] %[output:7052205d] %[output:921fadfa] %[output:4d398825] %[output:9787e9dc] %[output:63f60293] %[output:1892b727] %[output:50c41805] %[output:9027eeaf] %[output:684f19e7] %[output:0827b0c0] %[output:85ce146a] %[output:242e70e0] %[output:7007eef9] %[output:5a9fb54d] %[output:26a401c5] %[output:67352dc9] %[output:4098276c] %[output:6316c82d] %[output:0db25b35] %[output:8d02eb64] %[output:419fb1df] %[output:473a8cb6] %[output:94582730] %[output:0680536f] %[output:09a7ebfc] %[output:0fb20694] %[output:9544a24c] %[output:71436b10] %[output:88b7751a] %[output:8b769c98] %[output:351f8354] %[output:18cb1518] %[output:31a38038] %[output:53a6697f] %[output:6f9eddaa] %[output:8bd430e5] %[output:183deec7] %[output:22e83f2d] %[output:225d55c1] %[output:10d9813f] %[output:1c0eda00] %[output:8e76e35c] %[output:0bb19a05] %[output:323066dc] %[output:89e8d927] %[output:5abb16f8] %[output:165fd7a8] %[output:27cf1add] %[output:758d7805] %[output:6ba9b090] %[output:214a97e9] %[output:62cbc513] %[output:976ca751] %[output:115fb835] %[output:78530366] %[output:6bb7d67f] %[output:2ba914d8] %[output:41917a46] %[output:42be268f] %[output:1b864846] %[output:58d325f3] %[output:67d226f5] %[output:8f05db7c] %[output:83c1071b] %[output:414617aa] %[output:93a7c866] %[output:87734737] %[output:47a366aa] %[output:937360cd] %[output:85d98ea6] %[output:59849214] %[output:2a346737] %[output:41279d4f] %[output:2bf051b1] %[output:836ae16e] %[output:2ae17def] %[output:9bcedccb] %[output:83d6b4aa] %[output:7e622bad] %[output:1805106c] %[output:200860a3] %[output:50c5abc0] %[output:1b25e15a] %[output:62d63bee] %[output:8d713a11] %[output:315c5757] %[output:79644e85] %[output:6f2a8f53] %[output:0356a9cb] %[output:92abf472] %[output:4224c0c4] %[output:227d68d9] %[output:62bfdabe] %[output:115ae1dd] %[output:652d079f] %[output:43ba1976] %[output:936fd1e6] %[output:1761c807] %[output:8013802b] %[output:3483b6eb] %[output:65011042] %[output:3f93902f] %[output:5a54fc63] %[output:46ffef91] %[output:586c7b34] %[output:95b5635e] %[output:098ee27a] %[output:23e1e94b] %[output:88baeed7] %[output:8af7ec3d] %[output:0c89f0a3] %[output:2188a004] %[output:4a235452] %[output:42e876d9] %[output:02c59696] %[output:3cdd45bb] %[output:60a1cf34] %[output:3401325a] %[output:965aff62] %[output:0a1b6ec6] %[output:0263f737] %[output:56be215f] %[output:2069d09e] %[output:0cc76f3f] %[output:3e57e121] %[output:5769b848] %[output:6287271d] %[output:821ae56b]
        fprintf('  GPR(r_Cp) = %s\n', model.grRules{rec.r_Cp}); %[output:726e5166] %[output:3a1565df] %[output:233fd097] %[output:82e35756] %[output:403c7dbc] %[output:9943508a] %[output:228cf17d] %[output:8ebea8d3] %[output:39c0c7b7] %[output:4e0c83f1] %[output:420a2d92] %[output:830f323f] %[output:5c87ac3d] %[output:837d9784] %[output:9cc1ee9e] %[output:92e352f5] %[output:7a58bd80] %[output:61e7f1f5] %[output:7d0dd091] %[output:8d85614f] %[output:200a1415] %[output:9f88d6a3] %[output:02962b4e] %[output:1be49d2c] %[output:84e6b666] %[output:2367b67e] %[output:0fcd6821] %[output:9a06b37b] %[output:5672ddee] %[output:251bef3c] %[output:8e0e36de] %[output:76f7e102] %[output:54166a80] %[output:7666cbbe] %[output:273afeb0] %[output:711e02d0] %[output:730e22aa] %[output:705f987d] %[output:3683cd23] %[output:77860761] %[output:59e35237] %[output:2536ed0b] %[output:0b026df8] %[output:23d0594d] %[output:1d2ed9da] %[output:85aabcb4] %[output:9d4f69c5] %[output:371af19b] %[output:42931f50] %[output:4c38ae8c] %[output:1c2f9965] %[output:624bad03] %[output:0120bb78] %[output:3dcdcdc3] %[output:19ae58a6] %[output:31a9c62d] %[output:14520272] %[output:15d64e0d] %[output:5843df3c] %[output:210ecfa3] %[output:7bfd58eb] %[output:6ac09de6] %[output:0ca6383f] %[output:6677f733] %[output:05be9d26] %[output:44496883] %[output:1bae90a4] %[output:5954ed84] %[output:2e84a2c8] %[output:97c5066c] %[output:5f19e928] %[output:795481a1] %[output:34cd3712] %[output:0722df1f] %[output:8f0aa1b2] %[output:1616ddc3] %[output:77efd6bc] %[output:5ea1b95d] %[output:2adb8aba] %[output:120fd3b2] %[output:5266a2a0] %[output:771a6871] %[output:4df45d68] %[output:72800c69] %[output:9a67d2ba] %[output:5885d9d5] %[output:0f5a2cd0] %[output:5f826a31] %[output:916f6269] %[output:67751617] %[output:8cb79b28] %[output:90d5f09c] %[output:52b53f42] %[output:0a6c7bbb] %[output:76d56cb0] %[output:0ae95044] %[output:4ceb9254] %[output:95162f29] %[output:2d9181cb] %[output:2e4232c9] %[output:413c03a5] %[output:90eb24e8] %[output:6cdd017a] %[output:98665350] %[output:0a298316] %[output:08551844] %[output:40547158] %[output:88f2347c] %[output:8f32ead8] %[output:77c3d769] %[output:29a7f084] %[output:90b943c1] %[output:38274cba] %[output:04d17245] %[output:29786336] %[output:7d934f0c] %[output:100080d7] %[output:514dcba3] %[output:3e0b1d5f] %[output:0525186b] %[output:614967e8] %[output:92582d81] %[output:4a07aa46] %[output:786206c1] %[output:0c03371b] %[output:64079c6a] %[output:5c8e7f44] %[output:999ce18e] %[output:3a0b637b] %[output:2b5eff96] %[output:13877bb4] %[output:66104d9f] %[output:01e8a824] %[output:1d2df148] %[output:6371e5ac] %[output:39c87fc0] %[output:06d7751a] %[output:4c1343da] %[output:6839a340] %[output:6aa8bbda] %[output:2a0a6a8c] %[output:8fb61c6b] %[output:7b8587dd] %[output:04b937c5] %[output:93b94068] %[output:07629fcb] %[output:6212879c] %[output:2959d481] %[output:664ee7e4] %[output:86a75a4f] %[output:13c2fd8c] %[output:70a38319] %[output:239a18e6] %[output:410a0ef6] %[output:96770e16] %[output:16727cde] %[output:2f33f6ba] %[output:322e2456] %[output:80d0aef9] %[output:64b05e85] %[output:678c3783] %[output:427be006] %[output:667d2f2a] %[output:4aea5c9a] %[output:3affd9c8] %[output:36d0b64e] %[output:650b397a] %[output:5ec245e4] %[output:6f94daff] %[output:2b280fe6] %[output:219a1f1b] %[output:2e393d7d] %[output:1675ca50] %[output:7cd80102] %[output:9fcb5e5e] %[output:3f1144a6] %[output:5043eae5] %[output:3ad182c5] %[output:3898a663] %[output:006e6e12] %[output:1b385ee4] %[output:068dae8b] %[output:2917d814] %[output:06f338f3] %[output:83fffd2f] %[output:93f87385] %[output:35d444c4] %[output:68b8dd0f] %[output:25543ae0] %[output:74bf4e98] %[output:14aac98e] %[output:82467f82] %[output:6ce2876e] %[output:10cee372] %[output:3244234c] %[output:87732624] %[output:1d7d4954] %[output:2983952e] %[output:7d03a921] %[output:945bee65] %[output:2d7caaf3] %[output:732be933] %[output:314995a7] %[output:67d0bb4d] %[output:899b1066] %[output:5439e6df] %[output:348be378] %[output:1a3a84fc] %[output:9ec97ca4] %[output:081686c0] %[output:6326e146] %[output:02ce0112] %[output:5ac05008] %[output:2103ba0e] %[output:61f8061f] %[output:978e9da5] %[output:48a21072] %[output:091c9f97] %[output:9427462b] %[output:26d939c8] %[output:75a1ace0] %[output:93734db7] %[output:495719ec] %[output:077eb4b8] %[output:74a2f02e] %[output:24b8465d] %[output:461fde69] %[output:69710649] %[output:7569c32e] %[output:60907120] %[output:0c394c2e] %[output:95883223] %[output:3acb8156] %[output:66693ded] %[output:91f557e0] %[output:3838aee6] %[output:2997f0a2] %[output:40c2bc08] %[output:436f916e] %[output:3e938fbd] %[output:00e545f8] %[output:84c54b88] %[output:6a72c8c8] %[output:788d2e16] %[output:269a64ca] %[output:510c1363] %[output:3db89dd9] %[output:5339764a] %[output:276e9073] %[output:59a72174] %[output:595eb7c9] %[output:2d75c292] %[output:1881e162] %[output:00591a2b] %[output:761ff593] %[output:46a4765b] %[output:76570dc2] %[output:2ac85343] %[output:852af08b] %[output:31880522] %[output:74924819] %[output:13f2c075] %[output:32e550c1] %[output:2842cedf] %[output:85c286a4] %[output:5518494c] %[output:9a3006da] %[output:66efc5a8] %[output:3532d8f4] %[output:990570f2] %[output:486e50df] %[output:45cda27d] %[output:1ff1fb64] %[output:682c5025] %[output:32239d46] %[output:6a7ea151] %[output:61b1f72f] %[output:2dba0354] %[output:79460cae] %[output:89879ff9] %[output:31f8eec0] %[output:8d6b83de] %[output:0f292aa1] %[output:6085c791] %[output:94750650] %[output:91ed34c2] %[output:6021bd15] %[output:528e8955] %[output:1eb0b2a5] %[output:41a7d503] %[output:26cba808] %[output:2a12cfcf] %[output:00d630db] %[output:8db21f6e] %[output:905d46b6] %[output:868c9757] %[output:00c791a6] %[output:8ebfb962] %[output:1cd47f78] %[output:959ab479] %[output:9179d770] %[output:02b45064] %[output:9b7399e2] %[output:6b485795] %[output:01543104] %[output:50d2133c] %[output:381c8fbb] %[output:118e5bad] %[output:7e49ac2f] %[output:218fc89e] %[output:05b4ae21] %[output:4d42c96e] %[output:3561fc24] %[output:75c55551] %[output:33a4d49e] %[output:29d4f4ee] %[output:1f4f4185] %[output:5dcc9b1e] %[output:403e15f1] %[output:8470db03] %[output:4c6a361c] %[output:852c370f] %[output:77461bde] %[output:5861f8d1] %[output:52f1b7e0] %[output:13a4e321] %[output:21853e51] %[output:1b656929] %[output:6c9bc5e9] %[output:09c74c0f] %[output:74d5cca5] %[output:7524fa96] %[output:65ebd627] %[output:3dee9eb7] %[output:0642cb04] %[output:10691113] %[output:23dcab47] %[output:8bb6d2df] %[output:02ec6fa3] %[output:95515d2b] %[output:54700101] %[output:7eff6d06] %[output:93be624e] %[output:2f9890c3] %[output:0779382b] %[output:81506b21] %[output:25a1fece] %[output:03176f71] %[output:4a5d026a] %[output:325d40d1] %[output:72295523] %[output:34904b6c] %[output:481657da] %[output:0a2dd6d6] %[output:4d9e75d7] %[output:86c35573] %[output:52161fee] %[output:3668fd60] %[output:47a1bacb] %[output:3b41965c] %[output:9b12a899] %[output:4e347fec] %[output:62f90b2a] %[output:1a27260d] %[output:4d20d4d2] %[output:6745ff59] %[output:7afcbba0] %[output:99d8b6ba] %[output:8ed88977] %[output:9cf8c192] %[output:5b43d3ae] %[output:669248ef] %[output:0211e6ee] %[output:7824e02a] %[output:41f9ecc5] %[output:571dfb46] %[output:39d9e969] %[output:63246bef] %[output:392ecaa3] %[output:25713139] %[output:1e01c3d7] %[output:9b5dab2c] %[output:34f3359e] %[output:71b32469] %[output:02e04187] %[output:3284150b] %[output:9c2f76f5] %[output:828717eb] %[output:8c2ff8ac] %[output:70f8373b] %[output:80af16c7] %[output:2700a183] %[output:1fd5513c] %[output:7c5feb85] %[output:8af1d860] %[output:60d4e98d] %[output:8385661f] %[output:38d9b9e4] %[output:26506441] %[output:7cb4f358] %[output:40eae86a] %[output:6904c29f] %[output:597c8b53] %[output:00e19fb8] %[output:21bcf00a] %[output:87148d16] %[output:62a101f5] %[output:2568e299] %[output:38390145] %[output:19c0028f] %[output:244575a2] %[output:2d8dd0ab] %[output:0b18302a] %[output:532fd258] %[output:9193f3ab] %[output:1a572571] %[output:095f138a] %[output:3427801a] %[output:4d958338] %[output:315435cc] %[output:0273481a] %[output:74149295] %[output:4087f66c] %[output:7c0eeb09] %[output:9c3c6b35] %[output:9bfa801e] %[output:099b943f] %[output:8c6512cd] %[output:3c9b2e6e] %[output:02eb43e1] %[output:11a18b72] %[output:92b7cb2d] %[output:6298fdc0] %[output:20c9297c] %[output:44a2937f] %[output:4b2e22b1] %[output:4554cfde] %[output:369dbd47] %[output:7d67123b] %[output:565ab4b3] %[output:952a3249] %[output:8e8aff58] %[output:4732bcc2] %[output:67735e6c] %[output:67786db2] %[output:48d7db13] %[output:41324801] %[output:519172cb] %[output:0cae8718] %[output:8c84c6db] %[output:237afc80] %[output:920b3d10] %[output:57608b69] %[output:9f4d1525] %[output:8c7a09b4] %[output:61c0290c] %[output:50e02fcc] %[output:5629abae] %[output:8a5af7da] %[output:1853aefb] %[output:13c89cfd] %[output:410b3828] %[output:52c83729] %[output:9b1ccd13] %[output:0f10e896] %[output:1970c947] %[output:03d5607f] %[output:7d35eee1] %[output:1e050781] %[output:1d78db85] %[output:2c0e22f2] %[output:80329be6] %[output:0ec1d2e7] %[output:84ddc6e4] %[output:38c0b7bc] %[output:184af97e] %[output:8f276e1e] %[output:711a9150] %[output:3004c46f] %[output:27f1071c] %[output:13d8528b] %[output:581d631f] %[output:860c0295] %[output:3ae3d33e] %[output:80f809e8] %[output:972ee56c] %[output:7806015a] %[output:0376d8b4] %[output:834f0ab9] %[output:28592058] %[output:0d18da98] %[output:09d3bc8d] %[output:46eb5f70] %[output:424345d9] %[output:169885ee] %[output:4882c043] %[output:75f08f3f] %[output:422fc5c0] %[output:85b96464]
        fprintf('  GPR(r_Cg) = %s\n', model.grRules{rec.r_Cg}); %[output:1345a5f2] %[output:9fd0ddf9] %[output:0c06388c] %[output:78bdea82] %[output:19567586] %[output:183feaf0] %[output:8b239da4] %[output:75a89082] %[output:97d0d789] %[output:9c5b5c6e] %[output:99ff968b] %[output:07ea3d3e] %[output:7765bf22] %[output:02b630d0] %[output:3ba8a70d] %[output:608cd035] %[output:11169a6d] %[output:8f6db872] %[output:0e48ae17] %[output:5fceadd8] %[output:55f6f2f7] %[output:8bf83285] %[output:52c2885b] %[output:4e7928f8] %[output:43294562] %[output:8777d901] %[output:264add6d] %[output:74d76316] %[output:478ee891] %[output:2044631e] %[output:348a4b7a] %[output:733ea912] %[output:42d805b0] %[output:27af2963] %[output:4beedca0] %[output:8e36c44b] %[output:3dcd5b93] %[output:3b25fff0] %[output:47644492] %[output:3c22f4dd] %[output:58444083] %[output:673c322a] %[output:83199545] %[output:7cd5425a] %[output:310205b9] %[output:6fcfc9ea] %[output:2e3669e7] %[output:5e3fca52] %[output:9be95330] %[output:9f05f547] %[output:3b0fff0c] %[output:96f70bb9] %[output:389e05b5] %[output:9440e4d1] %[output:030a5af8] %[output:79dfa8fe] %[output:8f0cd67a] %[output:8847162a] %[output:7345f361] %[output:82ad03b9] %[output:98493506] %[output:00b3b220] %[output:7ff76ffa] %[output:859a29cb] %[output:1d6e85b3] %[output:54b34184] %[output:27ea1218] %[output:099e6741] %[output:9e2d3ae8] %[output:3f7b9b5e] %[output:4ab27dca] %[output:6b3e8dff] %[output:350709b5] %[output:961c8e17] %[output:16b91550] %[output:96cfae11] %[output:2b62cb6e] %[output:0c1a7d85] %[output:2eec927d] %[output:9ee127f6] %[output:104d6523] %[output:53a2560e] %[output:4ec682ed] %[output:5bce85d6] %[output:8e497efc] %[output:2105e2a6] %[output:1d921896] %[output:3ee3b600] %[output:3603d841] %[output:40bf1ba0] %[output:775ffa98] %[output:096f70cf] %[output:2332b963] %[output:69ddfa1e] %[output:9ac2d027] %[output:8f197131] %[output:89b193de] %[output:1433932e] %[output:81f1b08c] %[output:6fc3155d] %[output:4aadeea7] %[output:4dd1b8a6] %[output:00fcb936] %[output:21cfa3ca] %[output:79930d8c] %[output:0b8651eb] %[output:30c85285] %[output:8d925d84] %[output:60da3345] %[output:47147243] %[output:9c30cace] %[output:708d27fd] %[output:80f6ba13] %[output:12e17974] %[output:08a5430f] %[output:73a4a6c7] %[output:826dc833] %[output:6796654e] %[output:25b62f4f] %[output:2de287d1] %[output:758fa92d] %[output:4fcebf17] %[output:2ff81612] %[output:54b75497] %[output:2967e275] %[output:4030a1e1] %[output:5cc0fa92] %[output:0452e2ce] %[output:66bc6551] %[output:4d23e3b3] %[output:0f13f662] %[output:8e85394b] %[output:843d1706] %[output:06ca38aa] %[output:87387e17] %[output:2f6cd200] %[output:63ef1eb4] %[output:666e8ed9] %[output:02367a86] %[output:1fde282c] %[output:5307b551] %[output:7c2b3b8b] %[output:19c9b091] %[output:9bd110e0] %[output:0750c667] %[output:3a791cb1] %[output:893f2cc1] %[output:543c078f] %[output:6e19a930] %[output:929b206d] %[output:13286603] %[output:0b6e135a] %[output:52713ad8] %[output:7c67dd3b] %[output:4ffcba06] %[output:590c26a3] %[output:9f08c9b9] %[output:3be0b2ea] %[output:3b7de89a] %[output:3763ac8a] %[output:3fa4e212] %[output:62c62de3] %[output:8aa3582a] %[output:8a232de9] %[output:0a551857] %[output:0efc91d6] %[output:1a33c583] %[output:69bd9140] %[output:5c41e5c2] %[output:87133656] %[output:75da3f01] %[output:1ceb472d] %[output:8ed36f58] %[output:342a7dd9] %[output:396c0629] %[output:92000377] %[output:4690ea17] %[output:20d8e011] %[output:45528a47] %[output:92c9808f] %[output:0b210460] %[output:5325157d] %[output:783c98b5] %[output:7cde2a09] %[output:238cceca] %[output:947abd2f] %[output:5d347ab8] %[output:87712d69] %[output:2c7ade8b] %[output:7150b0c6] %[output:4e880cf7] %[output:90fe4a6f] %[output:7413ba04] %[output:5a684d5c] %[output:652a09c8] %[output:9a4f5e14] %[output:24be4082] %[output:0e5dcdd4] %[output:5932ddd0] %[output:35477083] %[output:2f12dbe1] %[output:57fd5763] %[output:909cd4d6] %[output:0fa7473a] %[output:70ae1f0f] %[output:975909f3] %[output:7ebef416] %[output:67a86635] %[output:13debe81] %[output:46800862] %[output:9cbdea23] %[output:92b66632] %[output:4875149e] %[output:3cb43c10] %[output:4ea3bdc6] %[output:1a3d1abc] %[output:0229da21] %[output:41c835f7] %[output:6b29a911] %[output:0e2b5f3d] %[output:6ab0424c] %[output:7afded06] %[output:7a13b6bd] %[output:1fa7b38f] %[output:1c2a2fd6] %[output:7357c516] %[output:53295cd7] %[output:106ef868] %[output:3273910f] %[output:5686a756] %[output:484d12ef] %[output:4240a343] %[output:35f42ce5] %[output:89a441ed] %[output:4d0f0609] %[output:2a00d7e5] %[output:02954e8d] %[output:61cf1ff0] %[output:6d70141e] %[output:2e6b5dfe] %[output:5fb6aec2] %[output:0075a7cb] %[output:3028b632] %[output:8f11cef0] %[output:99653d7d] %[output:8b5434d1] %[output:0b1f33fe] %[output:36a475f8] %[output:61f1fe27] %[output:0ac1f6e3] %[output:9b82c576] %[output:433e3e64] %[output:20d96a9f] %[output:9c77f2de] %[output:78d1080b] %[output:6c309ab1] %[output:0fc143ab] %[output:50ec109e] %[output:376917ff] %[output:8766bdce] %[output:191b67dc] %[output:93027442] %[output:06b3bc71] %[output:5d8afce2] %[output:490fc131] %[output:827a4710] %[output:896c3c03] %[output:08df83bd] %[output:7a96fd08] %[output:3c914305] %[output:86f42587] %[output:2a33ffee] %[output:717685a7] %[output:4b1072b2] %[output:43c439b4] %[output:40768988] %[output:9bcf3b0f] %[output:9368a586] %[output:7c0c3b1a] %[output:033bfef1] %[output:65509623] %[output:9015c6a0] %[output:5fcbfc03] %[output:13b8818c] %[output:4e116197] %[output:3e04e210] %[output:82f2ca75] %[output:36293012] %[output:80a25b87] %[output:65499dda] %[output:8c71e29a] %[output:40afbd65] %[output:7d4e00c0] %[output:1515de37] %[output:65922c78] %[output:6711d9a0] %[output:13888a8b] %[output:67dfd722] %[output:296a6525] %[output:8c3ade37] %[output:6d7a8c28] %[output:0ac696e3] %[output:106ed596] %[output:969711a9] %[output:42943ed3] %[output:55b2380c] %[output:41042c5a] %[output:9ae1a1c5] %[output:8c5cc371] %[output:3282d665] %[output:33678a0f] %[output:73701fd9] %[output:61b7255e] %[output:36683f27] %[output:25ad980d] %[output:851ebd44] %[output:3ad80065] %[output:2b433293] %[output:5c3f3d2e] %[output:4f8961bf] %[output:6d2328f1] %[output:4cbf7584] %[output:3a6427c1] %[output:22a2e1c9] %[output:072eedfa] %[output:7d97bf50] %[output:6cb17c15] %[output:88836c64] %[output:9b8aa9c0] %[output:5cf7cdfd] %[output:98583d82] %[output:9bd2fcca] %[output:14af12df] %[output:32684f25] %[output:526df449] %[output:7e6e3814] %[output:277ff934] %[output:7aea7290] %[output:6721bb78] %[output:45c664b3] %[output:17b0fc0c] %[output:890b0ded] %[output:61ba9f4a] %[output:194a240d] %[output:86d8fca5] %[output:1c180562] %[output:23bd3a78] %[output:1e119f98] %[output:3e5e8f25] %[output:5094a13c] %[output:3fd46c54] %[output:3dd188b0] %[output:3e4bc72b] %[output:7d6b94c6] %[output:56671865] %[output:79165323] %[output:70f71b19] %[output:32059e71] %[output:3fd28e93] %[output:157119a4] %[output:0e58077c] %[output:3c44993e] %[output:6ed4ea29] %[output:36267989] %[output:223c18ac] %[output:98ead2ff] %[output:951ec690] %[output:917c863d] %[output:945efed4] %[output:3250e04f] %[output:2ff34beb] %[output:2410ef07] %[output:6ea84fe2] %[output:49a8812e] %[output:6e4700b6] %[output:9a6daeed] %[output:7f2824c9] %[output:2c4b54c6] %[output:3b52104f] %[output:07444ff8] %[output:0385df11] %[output:41988b51] %[output:8c0b7f1c] %[output:319b9f10] %[output:1002af2a] %[output:165fa37b] %[output:85f424c9] %[output:36a3f934] %[output:0bfc58fc] %[output:04cba7db] %[output:0f74106f] %[output:87235cb0] %[output:07bf16fc] %[output:3e819d54] %[output:7a7853bd] %[output:0aadbef7] %[output:5e706948] %[output:5ebe84f4] %[output:75b8978c] %[output:9cc3db8d] %[output:9180cb67] %[output:86bf3394] %[output:47f4a7b7] %[output:60d8a7c5] %[output:8a12f862] %[output:34b551c9] %[output:84e9f59e] %[output:27427522] %[output:92dbd16f] %[output:95207f77] %[output:5d00d2b8] %[output:9b798597] %[output:85204701] %[output:2346fb6d] %[output:2a661c52] %[output:471fda00] %[output:3d5b47b9] %[output:0d0ed8ba] %[output:687293f5] %[output:92bf71cd] %[output:79a40a34] %[output:4059712f] %[output:69191571] %[output:97142c2e] %[output:28e55fdc] %[output:55cc439b] %[output:9c0f3a65] %[output:4d6395c4] %[output:0b38dcf3] %[output:5ce40190] %[output:28db91bd] %[output:044435f6] %[output:2ba32851] %[output:6f9f0bef] %[output:9c40ac4d] %[output:710ad041] %[output:4ce6d6a3] %[output:28dd2a23] %[output:14c9fee3] %[output:3cc911ba] %[output:0e0e5c54] %[output:7661a4a2] %[output:0f3d8361] %[output:25dc5c5c] %[output:8047e4fd] %[output:991d30e2] %[output:826688b5] %[output:7ad6af3d] %[output:93f83329] %[output:3de47713] %[output:8964840f] %[output:755bc37c] %[output:01c09aa1] %[output:6a4ef797] %[output:03d278bd] %[output:26c1a1c3] %[output:7e559e7a] %[output:81c21872] %[output:5bfdc40e] %[output:4c2a9173] %[output:363a9df5] %[output:9f1c778f] %[output:943b0097] %[output:1e74e0a1] %[output:6b3f8dbc] %[output:058e3ec4] %[output:1a6e9eb6] %[output:0f7c2bbf] %[output:79763e59] %[output:2c51715a] %[output:9c2549da] %[output:4c603382] %[output:7531af15] %[output:91228fd5] %[output:8d55f5b5] %[output:4f4450c5] %[output:4c1a6b53] %[output:3bf895e4] %[output:66ec09dc] %[output:460961bb] %[output:447ceef3] %[output:782071b1] %[output:134be620] %[output:4ea9f130] %[output:854c1077] %[output:130a281f] %[output:3fe53062] %[output:44ea66f7] %[output:8aa17401] %[output:033e397f] %[output:9467d873] %[output:6b70f22b] %[output:71b64e80] %[output:02b2cb6a] %[output:22516c38] %[output:42ed2005] %[output:4efd6008] %[output:002fd19b] %[output:9fd412ad]
        fprintf('  GPR(r_Cr) = %s\n', model.grRules{rec.r_Cr}); %[output:001dc101] %[output:9e6bda0c] %[output:2ba43c37] %[output:36ce42bf] %[output:13b09989] %[output:1af3f671] %[output:5f961668] %[output:3825a8e7] %[output:179a9e38] %[output:4009c4a6] %[output:0433f374] %[output:937c5e74] %[output:0d23fd0a] %[output:151e1a5e] %[output:8196dc19] %[output:3db5c5dc] %[output:0bc155db] %[output:9251808a] %[output:79fa7eb7] %[output:69add51a] %[output:753403b7] %[output:015f17fb] %[output:350f5e67] %[output:0077baf8] %[output:154b9fe0] %[output:1a589bc6] %[output:6e333b20] %[output:0c9c65f7] %[output:42b040e9] %[output:4da84241] %[output:0db67d5c] %[output:571d92d2] %[output:70c7e91d] %[output:48528592] %[output:5a978e21] %[output:55650436] %[output:234536c8] %[output:99096ea8] %[output:87e37c46] %[output:37e2412e] %[output:7bde32f2] %[output:62ce6f62] %[output:70e8dded] %[output:82bb5015] %[output:547cfab5] %[output:33ba8656] %[output:3be48109] %[output:139f79c3] %[output:7d95ae1e] %[output:09b3103b] %[output:09a8d452] %[output:47f23de8] %[output:4d7771bc] %[output:26693c13] %[output:4fbd17d5] %[output:946a362f] %[output:076a57d1] %[output:313fb5df] %[output:3fe74138] %[output:11a13fe5] %[output:6c8cfe1e] %[output:6240da1f] %[output:8efe4676] %[output:89e018ae] %[output:00d88bff] %[output:54c653e1] %[output:5ba5b171] %[output:4b9795ac] %[output:3f5b3c58] %[output:8e84d23d] %[output:5599eb24] %[output:3b5b02a3] %[output:858ae2bf] %[output:1522517b] %[output:1f7a8948] %[output:4e9fd9d2] %[output:2ce6cca2] %[output:8b2d5afd] %[output:5811dc69] %[output:261d2f55] %[output:4422d326] %[output:62a73a98] %[output:6a80b85c] %[output:755ca243] %[output:0688f676] %[output:44fa83f7] %[output:57819ba6] %[output:94afb47f] %[output:850003d5] %[output:65915348] %[output:4909b378] %[output:44f7cd07] %[output:77555e26] %[output:46e1e9ac] %[output:0b563091] %[output:957d22b4] %[output:2cee2f7a] %[output:3f5cef0d] %[output:336dfc54] %[output:357ab227] %[output:28c97839] %[output:7b52ec80] %[output:1f524b9a] %[output:40f00c1d] %[output:882a4c5c] %[output:9819a2c6] %[output:455806f9] %[output:85991f6d] %[output:9dbbc68b] %[output:221c5f3c] %[output:9be71a1b] %[output:0892d630] %[output:5d4d5748] %[output:76aa1fb1] %[output:642471da] %[output:9ddefcba] %[output:686712dc] %[output:364f0131] %[output:553cc7b6] %[output:60dd97e7] %[output:0bb498e2] %[output:20198b64] %[output:5995a79e] %[output:81f76006] %[output:3e0538dc] %[output:35d713b3] %[output:1b997a2e] %[output:6be53022] %[output:33dbcd67] %[output:319e763e] %[output:41c95c35] %[output:497b6fc9] %[output:4abd6434] %[output:919cf212] %[output:1b7becac] %[output:3c8a3553] %[output:74f65834] %[output:806476d2] %[output:8b4ef88f] %[output:1beeff45] %[output:943a5b40] %[output:5d97ceb1] %[output:483825c1] %[output:6133c2d6] %[output:9dda78e5] %[output:0d0f5177] %[output:6f6c8f47] %[output:402929b1] %[output:72a4d798] %[output:30df20bd] %[output:487fa04f] %[output:13281d57] %[output:72c16e0e] %[output:59b71001] %[output:187539b8] %[output:5d0e4e1c] %[output:5f821b0c] %[output:03a2008f] %[output:1ddd6416] %[output:0ae5e06a] %[output:6e5730e1] %[output:9c9046d9] %[output:57e1a99e] %[output:128a63a2] %[output:6eeff0cb] %[output:2248a662] %[output:4bb2b288] %[output:8c972825] %[output:6d08639e] %[output:7b382c9f] %[output:063a5f74] %[output:12ad0643] %[output:5d5d2678] %[output:919cfe23] %[output:7bc58ff7] %[output:4dfdb819] %[output:634b3a07] %[output:202bfea7] %[output:90295226] %[output:03b73189] %[output:86123863] %[output:21deaf5b] %[output:1d3a1af5] %[output:411558f4] %[output:3582f492] %[output:40b58781] %[output:092e6233] %[output:91566213] %[output:104c45a0] %[output:6e29a630] %[output:0a51965a] %[output:48596f90] %[output:51b26e0b] %[output:2f660c06] %[output:100fc03b] %[output:48fa1c77] %[output:1ba0293e] %[output:3d34e4f4] %[output:43280c3c] %[output:77c05889] %[output:6421a41b] %[output:38efcd4f] %[output:7983b2e5] %[output:4e368eb4] %[output:76097afc] %[output:7036f211] %[output:7cbd5c20] %[output:64a72d47] %[output:6206e013] %[output:0b643257] %[output:4bb51e16] %[output:6879a94a] %[output:493fc05f] %[output:8f1f4db0] %[output:2cdc4950] %[output:9fe93c72] %[output:9557dfb9] %[output:193a01cf] %[output:4406585f] %[output:202723f8] %[output:7f67e915] %[output:76da6b24] %[output:84c17723] %[output:7f6a4a09] %[output:35e6239f] %[output:86b6aab9] %[output:9daff6a9] %[output:96cc7f0a] %[output:82589d21] %[output:9d910b19] %[output:82db50ad] %[output:31e0bd12] %[output:1641d38f] %[output:6ecd2396] %[output:1c566f66] %[output:935e243b] %[output:66862dc7] %[output:2e5cdb7c] %[output:7e42ee3f] %[output:0e10553a] %[output:3ff4df54] %[output:22dc32aa] %[output:77a4aecc] %[output:4d7dab8b] %[output:504a78ca] %[output:3bcc5608] %[output:073fd259] %[output:9d012e51] %[output:6d1b0638] %[output:4f564a5f] %[output:8e81b7df] %[output:71dbf121] %[output:37221105] %[output:66b56ad1] %[output:47a0700d] %[output:5192ed95] %[output:87f83e45] %[output:3b2cd744] %[output:02e4273f] %[output:0de23352] %[output:689cd198] %[output:8c1b4e08] %[output:46679fd4] %[output:5aa22d1d] %[output:0829f021] %[output:15382885] %[output:56c0d2e1] %[output:03bbf093] %[output:378da015] %[output:9958a29e] %[output:0ceb4ea8] %[output:86ac9c32] %[output:9a0bd714] %[output:25df37af] %[output:9390cda0] %[output:8728f819] %[output:1c043c9c] %[output:04067779] %[output:1f3c623e] %[output:1fb871db] %[output:44046a68] %[output:1ac8d568] %[output:4bb8f540] %[output:0917f6fe] %[output:4d4fd206] %[output:973a7e36] %[output:27bf3395] %[output:83762fb2] %[output:38a1e089] %[output:91888963] %[output:4d899d69] %[output:6006a7d1] %[output:6221f475] %[output:22ca3c6c] %[output:4d7f7cdc] %[output:321ffbdf] %[output:2a4c6ddd] %[output:4f2d21d7] %[output:51eb82cf] %[output:5fd11695] %[output:4b6cb5ce] %[output:8cce911f] %[output:39356964] %[output:04cfde83] %[output:55e20174] %[output:4caf7f09] %[output:6f7a6273] %[output:03494e22] %[output:608cdfe2] %[output:0da20a89] %[output:56576c3f] %[output:85d17995] %[output:00fb54de] %[output:797e2cc4] %[output:07504aad] %[output:8086d95d] %[output:5864f34b] %[output:2ea870b9] %[output:165fb11f] %[output:99826d4e] %[output:7c8cad34] %[output:6febcc7e] %[output:7db09fdc] %[output:413a725d] %[output:305f8159] %[output:93d7b3da] %[output:345b4929] %[output:1ae2c2c0] %[output:8f82df58] %[output:228cc991] %[output:6101aaf0] %[output:9f2a5703] %[output:86dec0b2] %[output:2788ac7c] %[output:5d155a74] %[output:5db1a9c4] %[output:08679e98] %[output:69fcdf80] %[output:24d71fb6] %[output:3ea9feab] %[output:0186f782] %[output:46243af1] %[output:1fe391f9] %[output:98a81625] %[output:95bb4093] %[output:0cdce167] %[output:7d663a68] %[output:78a7d36f] %[output:41ef011f] %[output:6be39167] %[output:6707f046] %[output:0467a03a] %[output:99d5bd06] %[output:0cad9cf5] %[output:00443f90] %[output:784b19b1] %[output:6c8a069d] %[output:822f5d23] %[output:1234bcb4] %[output:0d452ecd] %[output:44734c37] %[output:189f1913] %[output:0a4397a1] %[output:4ed80847] %[output:20cc46a2] %[output:7cf5ce99] %[output:6b51a974] %[output:8c7454e4] %[output:82c127b4] %[output:54480bb3] %[output:48584ea8] %[output:4f858b59] %[output:1bfe81f5] %[output:1c5c27a4] %[output:24730319] %[output:0b0e1090] %[output:52e989ee] %[output:8432d15c] %[output:52c9d934] %[output:4e7e3558] %[output:189e472f] %[output:55253c78] %[output:6e033684] %[output:3c292ff1] %[output:74ed9d9b] %[output:2a0e1bfe] %[output:73c742ef] %[output:18145e65] %[output:98e1a147] %[output:0aa9dfc0] %[output:1b8f2037] %[output:2518fe01] %[output:433db78e] %[output:4c4f7fcc] %[output:5a66f0a1] %[output:715d7728] %[output:32635763] %[output:997e70a0] %[output:3ae9a81b] %[output:7809fec5] %[output:78f00328] %[output:4a2ced72] %[output:16949b7a] %[output:55e314ae] %[output:85cd2d8e] %[output:8f109f19] %[output:4fb290d9] %[output:6a22151a] %[output:14116a56] %[output:19626f5a] %[output:57c0c463] %[output:562fe160] %[output:1844dffe] %[output:9a76a14b] %[output:9449e096] %[output:6b67bfda] %[output:5f0b414b] %[output:2d001e90] %[output:0577262a] %[output:507fe250] %[output:72aa28bb] %[output:55d8284d] %[output:8d1e4816] %[output:8be616a4] %[output:2931141d] %[output:5b00e852] %[output:28db04a6] %[output:215f5b78] %[output:9f85ce8e] %[output:0a01730e] %[output:814ab8e3] %[output:777f4e51] %[output:2de45dc7] %[output:5b353517] %[output:5ae8aba0] %[output:14e0a57e] %[output:7de4f265] %[output:5021fc8e] %[output:37b1271b] %[output:8c671179] %[output:18c5326a] %[output:1319d506] %[output:9790f03d] %[output:18aca330] %[output:02ec9989] %[output:4eebb475] %[output:296d1501] %[output:96b1571c] %[output:5e75bd20] %[output:542c1d14] %[output:7755283a] %[output:9f9d03a6] %[output:14a4a490] %[output:38d3fd67] %[output:9aaa57b7] %[output:0ba3b505] %[output:59e81f23] %[output:0eb7d931] %[output:29e839f9] %[output:40ba01b0] %[output:9bea26c5] %[output:37db6da6] %[output:9815095a] %[output:302d1ff7] %[output:6727688d] %[output:208912a7] %[output:813b4909] %[output:4269db5d] %[output:44484d79] %[output:4f53d0ca] %[output:7c9fbf3d] %[output:00391f9b] %[output:0efe9db2] %[output:9e1c051a] %[output:331747a5] %[output:6c2f0fda] %[output:49d92f23] %[output:40313976] %[output:1f313628] %[output:758d4ff5] %[output:2d710021] %[output:4d706f09] %[output:21843349] %[output:8b3f9c1b] %[output:58a76f72] %[output:3ca11807] %[output:91ba1cc3] %[output:85466d34] %[output:96019eb4] %[output:067abfa9] %[output:79f68bb1] %[output:3c15bbbe] %[output:3ef7aae4] %[output:4a02d1a9] %[output:233869b8] %[output:39f6d939] %[output:5e78be45] %[output:4c27308d] %[output:8f0ab58b]
    end
    save('Block_2_joint_diff_counter_candidates.mat', 'Motifs_Cnt', 'step4');
end %[output:group:8d59eaae]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":45}
%---
%[output:762ffa34]
%   data: {"dataType":"text","outputData":{"text":"Loaded the differentiator","truncated":false}}
%---
%[output:10078da4]
%   data: {"dataType":"text","outputData":{"text":"Step 1: 9301 (r_Cp, m_d, m_c) triples\n","truncated":false}}
%---
%[output:2268009a]
%   data: {"dataType":"text","outputData":{"text":"Step 2: 127373 (r_Cp, r_Cg, m_d, m_c, m_Cs) tuples\n","truncated":false}}
%---
%[output:92f15b5e]
%   data: {"dataType":"text","outputData":{"text":"Step 3: 5048 full counter triples (r_Cp, r_Cg, r_Cr)\n","truncated":false}}
%---
%[output:3564bbb3]
%   data: {"dataType":"text","outputData":{"text":"Step 4: 3548 counters after gene-disjointness filter\n","truncated":false}}
%---
%[output:00c2e1ba]
%   data: {"dataType":"text","outputData":{"text":"Step 5: 143687 joint (differentiator + counter) candidates\n","truncated":false}}
%---
%[output:890139b4]
%   data: {"dataType":"text","outputData":{"text":"Joint candidates passing 7-reaction gene-disjointness: 72568 \/ 143687\n","truncated":false}}
%---
%[output:08093c42]
%   data: {"dataType":"text","outputData":{"text":"\n=== Joint differentiator + counter table (first 72568 rows) ===\n","truncated":false}}
%---
%[output:948bb2ed]
%   data: {"dataType":"text","outputData":{"text":"    <strong>diff_idx<\/strong>    <strong>cnt_idx<\/strong>      <strong>m_in<\/strong>        <strong>m_Df<\/strong>        <strong>m_d<\/strong>         <strong>m_Ds<\/strong>        <strong>m_c<\/strong>         <strong>m_Cs<\/strong>        <strong>r_DF<\/strong>        <strong>r_DS<\/strong>        <strong>r_Dp<\/strong>        <strong>r_Di<\/strong>        <strong>r_Cp<\/strong>        <strong>r_Cg<\/strong>        <strong>r_Cr<\/strong>  \n    <strong>________<\/strong>    <strong>_______<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>\n\n       699         11      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       700         11      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       701         11      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       702         11      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       703         11      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       704         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       705         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       706         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       707         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       708         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       709         11      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       710         11      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       711         11      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       712         11      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       713         11      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       714         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       715         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       716         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       717         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       718         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       719         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       720         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       721         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       722         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       723         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       724         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       725         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       726         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       727         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       728         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       729         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       730         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       731         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       732         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       733         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       734         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       735         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       736         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       737         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       738         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       739         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       740         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       741         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       749         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       750         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       751         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       752         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       753         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       754         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       755         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       756         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       757         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       758         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       759         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       760         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       761         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       762         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       763         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       764         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       765         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       766         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       767         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       768         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       769         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       770         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       771         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       772         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       773         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       774         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       775         11      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0714\"    \"r_2132\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       776         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       777         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       778         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       779         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       780         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       781         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       782         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       783         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       784         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       785         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       786         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       787         11      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       788         11      \"s_0066\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_2132\"    \"r_0714\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       789         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       790         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       791         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       792         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       793         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       794         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       795         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       796         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       797         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       798         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       799         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       800         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       801         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       802         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       803         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       804         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       805         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       806         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       807         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       808         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       809         11      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       810         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       811         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       812         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       813         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       814         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       815         11      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0714\"    \"r_2132\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       816         11      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0459\"    \"r_0307\"\n       694         12      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0018\"    \"r_0307\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       695         12      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0026\"    \"r_0307\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       696         12      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0538\"    \"r_0307\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       697         12      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0918\"    \"r_0307\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       698         12      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_1063\"    \"r_0307\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       704         12      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_0991\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       705         12      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_0991\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       706         12      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_0991\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       707         12      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_0991\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       708         12      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_0991\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       709         12      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0018\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       710         12      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0026\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       711         12      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0538\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       712         12      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0918\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       713         12      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_1063\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       770         12      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0018\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       771         12      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0026\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       772         12      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0538\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       773         12      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0918\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       774         12      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_1063\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       775         12      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0714\"    \"r_2132\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       783         12      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       784         12      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       785         12      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       786         12      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       787         12      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       788         12      \"s_0066\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_2132\"    \"r_0714\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       810         12      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0018\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       811         12      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0026\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       812         12      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0538\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       813         12      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0918\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       814         12      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_1063\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       815         12      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0714\"    \"r_2132\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0469\"    \"r_0476\"\n       699         14      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       700         14      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       701         14      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       702         14      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       703         14      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       704         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       705         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       706         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       707         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       708         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       709         14      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       710         14      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       711         14      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       712         14      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       713         14      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       714         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       715         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       716         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       717         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       718         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       719         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       720         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       721         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       722         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       723         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       724         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       725         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       726         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       727         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       728         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       729         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       730         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       731         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       732         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       733         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       734         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       735         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       736         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       737         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       738         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       739         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       740         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       741         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       749         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       750         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       751         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       752         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       753         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       754         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       755         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       756         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       757         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       758         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       759         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       760         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       761         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       762         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       763         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       764         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       765         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       766         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       767         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       768         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       769         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       770         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       771         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       772         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       773         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       774         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       775         14      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0714\"    \"r_2132\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       776         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       777         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       778         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       779         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       780         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       781         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       782         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       783         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       784         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       785         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       786         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       787         14      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       788         14      \"s_0066\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_2132\"    \"r_0714\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       789         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       790         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       791         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       792         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       793         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       794         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       795         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       796         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       797         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       798         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       799         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       800         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       801         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       802         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       803         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       804         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       805         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       806         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       807         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       808         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       809         14      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       810         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       811         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       812         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       813         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       814         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       815         14      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0714\"    \"r_2132\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       816         14      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0973\"    \"r_0307\"\n       694         15      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0018\"    \"r_0307\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       695         15      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0026\"    \"r_0307\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       696         15      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0538\"    \"r_0307\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       697         15      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0918\"    \"r_0307\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       698         15      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_1063\"    \"r_0307\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       704         15      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_0991\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       705         15      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_0991\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       706         15      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_0991\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       707         15      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_0991\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       708         15      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_0991\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       709         15      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0018\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       710         15      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0026\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       711         15      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0538\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       712         15      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0918\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       713         15      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_1063\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       770         15      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0018\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       771         15      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0026\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       772         15      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0538\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       773         15      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0918\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       774         15      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_1063\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       775         15      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0714\"    \"r_2132\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       783         15      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       784         15      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       785         15      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       786         15      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       787         15      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       788         15      \"s_0066\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_0991\"    \"r_2132\"    \"r_0714\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       810         15      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0018\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       811         15      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0026\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       812         15      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0538\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       813         15      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_0918\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       814         15      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0470\"    \"r_1063\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       815         15      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_0991\"    \"r_0714\"    \"r_2132\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_0989\"    \"r_0476\"\n       699         16      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       700         16      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       701         16      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       702         16      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       703         16      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       704         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       705         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       706         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       707         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       708         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       709         16      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       710         16      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       711         16      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       712         16      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       713         16      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       714         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       715         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       716         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       717         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       718         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       719         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       720         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       721         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       722         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       723         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       724         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       725         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       726         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       727         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       728         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       729         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       730         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       731         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       732         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       733         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       734         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       735         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       736         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       737         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       738         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       739         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       740         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       741         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       749         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       750         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       751         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       752         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       753         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       754         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       755         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       756         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       757         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       758         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       759         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       760         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       761         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       762         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       763         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       764         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       765         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       766         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       767         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       768         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       769         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       770         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       771         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       772         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       773         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       774         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       775         16      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0714\"    \"r_2132\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       776         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       777         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       778         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       779         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       780         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       781         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       782         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       783         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       784         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       785         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       786         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       787         16      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       788         16      \"s_0066\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_2132\"    \"r_0714\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       789         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       790         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       791         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       792         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       793         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       794         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       795         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       796         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       797         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       798         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       799         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       800         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       801         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       802         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       803         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       804         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       805         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       806         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       807         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       808         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       809         16      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       810         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       811         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       812         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       813         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       814         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       815         16      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0714\"    \"r_2132\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       816         16      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1069\"    \"r_0307\"\n       699         17      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       700         17      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       701         17      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       702         17      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       703         17      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       704         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       705         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       706         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       707         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       708         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       709         17      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       710         17      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       711         17      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       712         17      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       713         17      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       714         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       715         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       716         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       717         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       718         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       719         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       720         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       721         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       722         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       723         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       724         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       725         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       726         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       727         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       728         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       729         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       730         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       731         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       732         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       733         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       734         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       735         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       736         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       737         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       738         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       739         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       740         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       741         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       749         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       750         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       751         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       752         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       753         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       754         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       755         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       756         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       757         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       758         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       759         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       760         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       761         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       762         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       763         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       764         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       765         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       766         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       767         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       768         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       769         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       770         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       771         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       772         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       773         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       774         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       775         17      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0714\"    \"r_2132\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       776         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       777         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       778         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       779         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       780         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       781         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       782         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       783         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       784         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       785         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       786         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       787         17      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       788         17      \"s_0066\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_2132\"    \"r_0714\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       789         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       790         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       791         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       792         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       793         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       794         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       795         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       796         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       797         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       798         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       799         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       800         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       801         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       802         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       803         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       804         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       805         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       806         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       807         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       808         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       809         17      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       810         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       811         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       812         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       813         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       814         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       815         17      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0714\"    \"r_2132\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       816         17      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1084\"    \"r_0307\"\n       699         18      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       700         18      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       701         18      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       702         18      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       703         18      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       704         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       705         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       706         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       707         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       708         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       709         18      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       710         18      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       711         18      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       712         18      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       713         18      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       714         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       715         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       716         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       717         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       718         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       719         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       720         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       721         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       722         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       723         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       724         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       725         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       726         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       727         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       728         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       729         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       730         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       731         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       732         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       733         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       734         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       735         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       736         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       737         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       738         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       739         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       740         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       741         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       749         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       750         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       751         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       752         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       753         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       754         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       755         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0514\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       756         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       757         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       758         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       759         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       760         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       761         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       762         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0563\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       763         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       764         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       765         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       766         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       767         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       768         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       769         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0768\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       770         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       771         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       772         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       773         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       774         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       775         18      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0714\"    \"r_2132\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       776         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_0771\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       777         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       778         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       779         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       780         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       781         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       782         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0476\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       783         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       784         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       785         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       786         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       787         18      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       788         18      \"s_0066\"    \"s_0180\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_2132\"    \"r_0714\"    \"r_1029\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       789         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       790         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       791         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       792         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       793         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       794         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       795         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_1197\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       796         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       797         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       798         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       799         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       800         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       801         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       802         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4211\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       803         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       804         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       805         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       806         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       807         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       808         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       809         18      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_4212\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       810         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       811         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       812         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       813         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       814         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       815         18      \"s_0066\"    \"s_1203\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0714\"    \"r_2132\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       816         18      \"s_0991\"    \"s_1203\"    \"s_0794\"    \"s_0999\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0476\"    \"r_4254\"    \"r_0472\"    \"r_0014\"    \"r_1275\"    \"r_0307\"\n       699         19      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       700         19      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       701         19      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       702         19      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       703         19      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_0476\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       704         19      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0018\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       705         19      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0026\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       706         19      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0538\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       707         19      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_0918\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       708         19      \"s_0991\"    \"s_0180\"    \"s_0794\"    \"s_0419\"    \"s_0419\"    \"s_1559\"    \"r_1063\"    \"r_0470\"    \"r_1029\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       709         19      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0018\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       710         19      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0026\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       711         19      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0538\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       712         19      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_0918\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       713         19      \"s_0991\"    \"s_0419\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0470\"    \"r_1063\"    \"r_4278\"    \"r_0471\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       714         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       715         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       716         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       717         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       718         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       719         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       720         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0079\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       721         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       722         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       723         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       724         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       725         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       726         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       727         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0203\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       728         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       729         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       730         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       731         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       732         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       733         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_1063\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       734         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_1203\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0211\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       735         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0018\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       736         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0026\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       737         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0470\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       738         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0538\"    \"r_0250\"    \"r_0472\"    \"r_0014\"    \"r_1276\"    \"r_0307\"\n       739         19      \"s_0991\"    \"s_0999\"    \"s_0794\"    \"s_0180\"    \"s_0419\"    \"s_1559\"    \"r_0476\"    \"r_0918\"    \"r_0250\"    \"r_0472\"    \"r","truncated":true}}
%---
%[output:606d1ac5]
%   data: {"dataType":"text","outputData":{"text":"\nSaved joint_diff_counter_candidates.mat\n","truncated":false}}
%---
%[output:7654d8dc]
%   data: {"dataType":"text","outputData":{"text":"\n=== Counter candidates (3-reaction CRN) ===\n","truncated":false}}
%---
%[output:596067ef]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 1 ---\n","truncated":false}}
%---
%[output:3d67aa74]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0025 | (R)-lactate\n","truncated":false}}
%---
%[output:5deabfec]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:3e7be64d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5d868d56]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0001 | (R)-lactate:ferricytochrome-c 2-oxidoreductase\n","truncated":false}}
%---
%[output:45ef7353]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4236 | (R)-lactate hydro-lyase\n","truncated":false}}
%---
%[output:4a5015c6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0959 | pyruvate decarboxylase\n","truncated":false}}
%---
%[output:726e5166]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YDL174C and YEL039C) or (YDL174C and YJR048W) or (YEL039C and YEL071W) or (YEL071W and YJR048W)\n","truncated":false}}
%---
%[output:1345a5f2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDR533C or YMR322C or YOR391C or YPL280W\n","truncated":false}}
%---
%[output:001dc101]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:12e70884]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 2 ---\n","truncated":false}}
%---
%[output:727af5ac]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 3 ---\n","truncated":false}}
%---
%[output:243092d5]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 4 ---\n","truncated":false}}
%---
%[output:2cd283dc]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 5 ---\n","truncated":false}}
%---
%[output:963f3445]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 6 ---\n","truncated":false}}
%---
%[output:7f450d27]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 7 ---\n","truncated":false}}
%---
%[output:8d70a876]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 8 ---\n","truncated":false}}
%---
%[output:1074bcf1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 9 ---\n","truncated":false}}
%---
%[output:72225947]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 10 ---\n","truncated":false}}
%---
%[output:9e6bda0c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:9fd0ddf9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDR533C or YMR322C or YOR391C or YPL280W\n","truncated":false}}
%---
%[output:3a1565df]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YDL174C and YEL039C) or (YDL174C and YJR048W) or (YEL039C and YEL071W) or (YEL071W and YJR048W)\n","truncated":false}}
%---
%[output:50e2e668]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:80c13a38]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4236 | (R)-lactate hydro-lyase\n","truncated":false}}
%---
%[output:2c37fb9b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0001 | (R)-lactate:ferricytochrome-c 2-oxidoreductase\n","truncated":false}}
%---
%[output:5a089de4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:30e55f88]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:2be6c316]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0025 | (R)-lactate\n","truncated":false}}
%---
%[output:7cae1de4]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 11 ---\n","truncated":false}}
%---
%[output:2ba43c37]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YGL080W and YGR243W) or (YGL080W and YHR162W)\n","truncated":false}}
%---
%[output:91792c1a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 12 ---\n","truncated":false}}
%---
%[output:0c06388c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDR533C or YMR322C or YOR391C or YPL280W\n","truncated":false}}
%---
%[output:99229002]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 13 ---\n","truncated":false}}
%---
%[output:233fd097]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YDL174C and YEL039C) or (YDL174C and YJR048W) or (YEL039C and YEL071W) or (YEL071W and YJR048W)\n","truncated":false}}
%---
%[output:5febab0d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 14 ---\n","truncated":false}}
%---
%[output:8e06b1ca]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2034 | pyruvate transport\n","truncated":false}}
%---
%[output:73fc6839]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 15 ---\n","truncated":false}}
%---
%[output:54db9230]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4236 | (R)-lactate hydro-lyase\n","truncated":false}}
%---
%[output:92cee311]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 16 ---\n","truncated":false}}
%---
%[output:36ce42bf]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:10c86049]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0001 | (R)-lactate:ferricytochrome-c 2-oxidoreductase\n","truncated":false}}
%---
%[output:32b3c01d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 17 ---\n","truncated":false}}
%---
%[output:78bdea82]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDR533C or YMR322C or YOR391C or YPL280W\n","truncated":false}}
%---
%[output:4a6eeeeb]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:7127c606]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 18 ---\n","truncated":false}}
%---
%[output:82e35756]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YDL174C and YEL039C) or (YDL174C and YJR048W) or (YEL039C and YEL071W) or (YEL071W and YJR048W)\n","truncated":false}}
%---
%[output:810aa9e7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:5053a5b1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 19 ---\n","truncated":false}}
%---
%[output:3877461e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:13b09989]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:6d05f825]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0025 | (R)-lactate\n","truncated":false}}
%---
%[output:7ddb315d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 20 ---\n","truncated":false}}
%---
%[output:50d800b8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4236 | (R)-lactate hydro-lyase\n","truncated":false}}
%---
%[output:19567586]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:5f932f1d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 21 ---\n","truncated":false}}
%---
%[output:95fbf5f6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 22 ---\n","truncated":false}}
%---
%[output:83d93a25]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0001 | (R)-lactate:ferricytochrome-c 2-oxidoreductase\n","truncated":false}}
%---
%[output:403c7dbc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YDL174C and YEL039C) or (YDL174C and YJR048W) or (YEL039C and YEL071W) or (YEL071W and YJR048W)\n","truncated":false}}
%---
%[output:8369850b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 23 ---\n","truncated":false}}
%---
%[output:1af3f671]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:65c9fedc]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 24 ---\n","truncated":false}}
%---
%[output:6e7ea9f2]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4dca5dcc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0959 | pyruvate decarboxylase\n","truncated":false}}
%---
%[output:7d1d72fe]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 25 ---\n","truncated":false}}
%---
%[output:183feaf0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:304b9fa1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 26 ---\n","truncated":false}}
%---
%[output:7246e0e8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:1af78b63]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:09b18094]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 27 ---\n","truncated":false}}
%---
%[output:5f961668]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YGL080W and YGR243W) or (YGL080W and YHR162W)\n","truncated":false}}
%---
%[output:9943508a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YDL174C and YEL039C) or (YDL174C and YJR048W) or (YEL039C and YEL071W) or (YEL071W and YJR048W)\n","truncated":false}}
%---
%[output:7e562fbe]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 28 ---\n","truncated":false}}
%---
%[output:938b7b8f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0025 | (R)-lactate\n","truncated":false}}
%---
%[output:74c23a4c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0001 | (R)-lactate:ferricytochrome-c 2-oxidoreductase\n","truncated":false}}
%---
%[output:21c461f3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 29 ---\n","truncated":false}}
%---
%[output:8b239da4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:46c98303]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:6f9fd200]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 30 ---\n","truncated":false}}
%---
%[output:5bd01679]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 31 ---\n","truncated":false}}
%---
%[output:3825a8e7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:98145ac0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:50f104cb]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 32 ---\n","truncated":false}}
%---
%[output:228cf17d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YDL174C and YEL039C) or (YDL174C and YJR048W) or (YEL039C and YEL071W) or (YEL071W and YJR048W)\n","truncated":false}}
%---
%[output:50137529]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:88c79437]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 33 ---\n","truncated":false}}
%---
%[output:8327b81c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 34 ---\n","truncated":false}}
%---
%[output:75a89082]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:273f522a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:12cdf022]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 35 ---\n","truncated":false}}
%---
%[output:179a9e38]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:21632340]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2034 | pyruvate transport\n","truncated":false}}
%---
%[output:6bb59d60]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0001 | (R)-lactate:ferricytochrome-c 2-oxidoreductase\n","truncated":false}}
%---
%[output:354a9704]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 36 ---\n","truncated":false}}
%---
%[output:1e6f07e9]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 37 ---\n","truncated":false}}
%---
%[output:8ebea8d3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YEL039C and YML054C) or (YJR048W and YML054C)\n","truncated":false}}
%---
%[output:1cdd234c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0025 | (R)-lactate\n","truncated":false}}
%---
%[output:9a811f1d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 38 ---\n","truncated":false}}
%---
%[output:97d0d789]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:119a2049]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:4009c4a6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YGL080W and YGR243W) or (YGL080W and YHR162W)\n","truncated":false}}
%---
%[output:6f3eefd6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8628159f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 39 ---\n","truncated":false}}
%---
%[output:72953939]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 40 ---\n","truncated":false}}
%---
%[output:86b0163c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0959 | pyruvate decarboxylase\n","truncated":false}}
%---
%[output:06f1ec41]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 41 ---\n","truncated":false}}
%---
%[output:38b6a466]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 42 ---\n","truncated":false}}
%---
%[output:39c0c7b7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YEL039C and YML054C) or (YJR048W and YML054C)\n","truncated":false}}
%---
%[output:805e182c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0001 | (R)-lactate:ferricytochrome-c 2-oxidoreductase\n","truncated":false}}
%---
%[output:9c5b5c6e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:0433f374]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:33d215d8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:88d71ed5]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 43 ---\n","truncated":false}}
%---
%[output:579c7420]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 44 ---\n","truncated":false}}
%---
%[output:3c95432a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:348c6870]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 45 ---\n","truncated":false}}
%---
%[output:5fe43a3a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 46 ---\n","truncated":false}}
%---
%[output:8643b86a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:32a6c7bb]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4e0c83f1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YEL039C and YML054C) or (YJR048W and YML054C)\n","truncated":false}}
%---
%[output:937c5e74]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPR035W\n","truncated":false}}
%---
%[output:99ff968b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR018C\n","truncated":false}}
%---
%[output:77c477b4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0025 | (R)-lactate\n","truncated":false}}
%---
%[output:6f46adf2]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 47 ---\n","truncated":false}}
%---
%[output:8d737617]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 48 ---\n","truncated":false}}
%---
%[output:1c70dfc0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0004 | (S)-lactate:ferricytochrome-c 2-oxidoreductase\n","truncated":false}}
%---
%[output:1cbc65dd]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 49 ---\n","truncated":false}}
%---
%[output:8f894226]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 50 ---\n","truncated":false}}
%---
%[output:4144e535]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:327c6217]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:0d23fd0a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:9f715819]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2034 | pyruvate transport\n","truncated":false}}
%---
%[output:07ea3d3e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YMR250W\n","truncated":false}}
%---
%[output:420a2d92]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOL066C\n","truncated":false}}
%---
%[output:21a82635]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 51 ---\n","truncated":false}}
%---
%[output:30b82c99]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 52 ---\n","truncated":false}}
%---
%[output:6057c254]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 53 ---\n","truncated":false}}
%---
%[output:4bb0aeff]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:688d4779]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 54 ---\n","truncated":false}}
%---
%[output:74718a55]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 55 ---\n","truncated":false}}
%---
%[output:151e1a5e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:46681df3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0004 | (S)-lactate:ferricytochrome-c 2-oxidoreductase\n","truncated":false}}
%---
%[output:779e6a55]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0025 | (R)-lactate\n","truncated":false}}
%---
%[output:7765bf22]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:047cd987]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:830f323f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOL066C\n","truncated":false}}
%---
%[output:99353125]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:870f125e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 56 ---\n","truncated":false}}
%---
%[output:81b22fde]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 57 ---\n","truncated":false}}
%---
%[output:4323dbc6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 58 ---\n","truncated":false}}
%---
%[output:8196dc19]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPR035W\n","truncated":false}}
%---
%[output:98141f95]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:72b63095]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 59 ---\n","truncated":false}}
%---
%[output:3f78c1d8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 60 ---\n","truncated":false}}
%---
%[output:02b630d0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:226c96fb]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1e9aff1d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 61 ---\n","truncated":false}}
%---
%[output:5c87ac3d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOL066C\n","truncated":false}}
%---
%[output:953cad0b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0004 | (S)-lactate:ferricytochrome-c 2-oxidoreductase\n","truncated":false}}
%---
%[output:260cea50]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0476 | glutamine synthetase\n","truncated":false}}
%---
%[output:3db5c5dc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:18b69a13]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0459 | galactose-1-phosphate uridylyltransferase\n","truncated":false}}
%---
%[output:0daea8c6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 62 ---\n","truncated":false}}
%---
%[output:512d87b4]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 63 ---\n","truncated":false}}
%---
%[output:105479e6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 64 ---\n","truncated":false}}
%---
%[output:3ba8a70d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNR050C\n","truncated":false}}
%---
%[output:4b597c27]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0063 | (S)-lactate\n","truncated":false}}
%---
%[output:11088e9f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 65 ---\n","truncated":false}}
%---
%[output:9882674e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 66 ---\n","truncated":false}}
%---
%[output:837d9784]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOL066C\n","truncated":false}}
%---
%[output:0bc155db]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:8587b5dc]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:4c43160c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 67 ---\n","truncated":false}}
%---
%[output:8ede528c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:41a4e91c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:169bf3de]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:608cd035]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL103C\n","truncated":false}}
%---
%[output:1c9255a0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0014 | 2,5-diamino-6-ribitylamino-4(3H)-pyrimidinone 5'-phosphate deaminase\n","truncated":false}}
%---
%[output:3a1713f4]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 68 ---\n","truncated":false}}
%---
%[output:386b408f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 69 ---\n","truncated":false}}
%---
%[output:9251808a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:86ab39cd]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 70 ---\n","truncated":false}}
%---
%[output:9cc1ee9e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOL066C\n","truncated":false}}
%---
%[output:73838060]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 71 ---\n","truncated":false}}
%---
%[output:33bae055]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 72 ---\n","truncated":false}}
%---
%[output:7dee8915]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 73 ---\n","truncated":false}}
%---
%[output:1f7650b8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:11169a6d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YHL012W or YKL035W\n","truncated":false}}
%---
%[output:206d0b47]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0063 | (S)-lactate\n","truncated":false}}
%---
%[output:1b418a9b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 74 ---\n","truncated":false}}
%---
%[output:79fa7eb7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:64cd8aac]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:8d30f0c4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:82799960]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0014 | 2,5-diamino-6-ribitylamino-4(3H)-pyrimidinone 5'-phosphate deaminase\n","truncated":false}}
%---
%[output:92e352f5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOL066C\n","truncated":false}}
%---
%[output:93d77437]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:8f0eab80]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 75 ---\n","truncated":false}}
%---
%[output:6ef810fe]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 76 ---\n","truncated":false}}
%---
%[output:8f6db872]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR192W\n","truncated":false}}
%---
%[output:19021326]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 77 ---\n","truncated":false}}
%---
%[output:69add51a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR204W\n","truncated":false}}
%---
%[output:9b13ae61]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0476 | glutamine synthetase\n","truncated":false}}
%---
%[output:8ded4164]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 78 ---\n","truncated":false}}
%---
%[output:257e13db]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 79 ---\n","truncated":false}}
%---
%[output:37cd056e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 80 ---\n","truncated":false}}
%---
%[output:4dc7bab7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0973 | ribonucleoside-triphosphate reductase (UTP)\n","truncated":false}}
%---
%[output:7a58bd80]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOL066C\n","truncated":false}}
%---
%[output:6a1050ee]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 81 ---\n","truncated":false}}
%---
%[output:264bc238]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 82 ---\n","truncated":false}}
%---
%[output:0e48ae17]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR192W\n","truncated":false}}
%---
%[output:753403b7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:23016969]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0014 | 2,5-diamino-6-ribitylamino-4(3H)-pyrimidinone 5'-phosphate deaminase\n","truncated":false}}
%---
%[output:83728a93]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0063 | (S)-lactate\n","truncated":false}}
%---
%[output:7af96c1a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:08bdd325]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:9a4373b9]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:8914c546]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 83 ---\n","truncated":false}}
%---
%[output:8695b7e1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 84 ---\n","truncated":false}}
%---
%[output:61e7f1f5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOL066C\n","truncated":false}}
%---
%[output:6526e319]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 85 ---\n","truncated":false}}
%---
%[output:015f17fb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR166C\n","truncated":false}}
%---
%[output:5fceadd8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:036ca5cd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:342073c7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 86 ---\n","truncated":false}}
%---
%[output:98865c68]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 87 ---\n","truncated":false}}
%---
%[output:0c5c47fc]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 88 ---\n","truncated":false}}
%---
%[output:0d8aa916]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0014 | 2,5-diamino-6-ribitylamino-4(3H)-pyrimidinone 5'-phosphate deaminase\n","truncated":false}}
%---
%[output:133e3197]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:2e2bca33]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 89 ---\n","truncated":false}}
%---
%[output:1a9d7d9a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 90 ---\n","truncated":false}}
%---
%[output:350f5e67]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:7d0dd091]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOL066C\n","truncated":false}}
%---
%[output:55f6f2f7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:267c8192]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:70367429]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 91 ---\n","truncated":false}}
%---
%[output:69666f37]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:6f172576]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1069 | UDP-N-acetylglucosamine diphosphorylase\n","truncated":false}}
%---
%[output:065b66b7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:596e1c5d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 92 ---\n","truncated":false}}
%---
%[output:455bc3c0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 93 ---\n","truncated":false}}
%---
%[output:0077baf8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:77abd6ee]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:84a2d9ca]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 94 ---\n","truncated":false}}
%---
%[output:8bf83285]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNL316C\n","truncated":false}}
%---
%[output:8d85614f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR153W\n","truncated":false}}
%---
%[output:729c198d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0014 | 2,5-diamino-6-ribitylamino-4(3H)-pyrimidinone 5'-phosphate deaminase\n","truncated":false}}
%---
%[output:220a7d77]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 95 ---\n","truncated":false}}
%---
%[output:988c3d90]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 96 ---\n","truncated":false}}
%---
%[output:05792a32]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 97 ---\n","truncated":false}}
%---
%[output:7fb6df0b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:154b9fe0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:015fe111]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1084 | UTP-glucose-1-phosphate uridylyltransferase\n","truncated":false}}
%---
%[output:720c2e61]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 98 ---\n","truncated":false}}
%---
%[output:785e4e5f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 99 ---\n","truncated":false}}
%---
%[output:52c2885b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:87c7b0c6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:200a1415]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR153W\n","truncated":false}}
%---
%[output:30ba1ed4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:11b95d13]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 100 ---\n","truncated":false}}
%---
%[output:8010da10]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:1a589bc6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:623f96d6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0014 | 2,5-diamino-6-ribitylamino-4(3H)-pyrimidinone 5'-phosphate deaminase\n","truncated":false}}
%---
%[output:883569b1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 101 ---\n","truncated":false}}
%---
%[output:886c7f11]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 102 ---\n","truncated":false}}
%---
%[output:95a18892]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 103 ---\n","truncated":false}}
%---
%[output:4e7928f8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:8b8a800f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1275 | UTP transport\n","truncated":false}}
%---
%[output:94daa8cf]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 104 ---\n","truncated":false}}
%---
%[output:9f88d6a3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR153W\n","truncated":false}}
%---
%[output:9dfce923]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0732 | methylenetetrahydrofolate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6e333b20]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:6a865d87]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:9916dc31]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 105 ---\n","truncated":false}}
%---
%[output:246819d5]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 106 ---\n","truncated":false}}
%---
%[output:603a87d4]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 107 ---\n","truncated":false}}
%---
%[output:0b706432]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:43294562]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YMR250W\n","truncated":false}}
%---
%[output:0a31d30b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0014 | 2,5-diamino-6-ribitylamino-4(3H)-pyrimidinone 5'-phosphate deaminase\n","truncated":false}}
%---
%[output:4fdba30c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 108 ---\n","truncated":false}}
%---
%[output:1c16a5e6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 109 ---\n","truncated":false}}
%---
%[output:0c9c65f7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:02962b4e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR153W\n","truncated":false}}
%---
%[output:1e6b1a11]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1276 | UTP\/UMP antiport\n","truncated":false}}
%---
%[output:5fc3bc7f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:2d363e7e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:7af38729]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 110 ---\n","truncated":false}}
%---
%[output:36ff2443]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 111 ---\n","truncated":false}}
%---
%[output:8777d901]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YMR250W\n","truncated":false}}
%---
%[output:282f98c0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:10d76a8c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 112 ---\n","truncated":false}}
%---
%[output:42b040e9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:183ac662]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 113 ---\n","truncated":false}}
%---
%[output:6dc480a7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 114 ---\n","truncated":false}}
%---
%[output:1be49d2c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR153W\n","truncated":false}}
%---
%[output:14914d87]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0014 | 2,5-diamino-6-ribitylamino-4(3H)-pyrimidinone 5'-phosphate deaminase\n","truncated":false}}
%---
%[output:0c152fd4]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 115 ---\n","truncated":false}}
%---
%[output:8dfc82dc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0939 | prephenate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:3857d23c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:264add6d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YPR035W\n","truncated":false}}
%---
%[output:8500c489]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:4da84241]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR115C and YGL154C\n","truncated":false}}
%---
%[output:56da8b1b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 116 ---\n","truncated":false}}
%---
%[output:778f141c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 117 ---\n","truncated":false}}
%---
%[output:0016392f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 118 ---\n","truncated":false}}
%---
%[output:79fd2239]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:84e6b666]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER152C or YGL202W or YJL060W\n","truncated":false}}
%---
%[output:48b1f511]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:8cf2de9b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 119 ---\n","truncated":false}}
%---
%[output:21aa9f25]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 120 ---\n","truncated":false}}
%---
%[output:74d76316]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNR050C\n","truncated":false}}
%---
%[output:0db67d5c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR115C and YGL154C\n","truncated":false}}
%---
%[output:2587ca2d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:58c0bc82]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0014 | 2,5-diamino-6-ribitylamino-4(3H)-pyrimidinone 5'-phosphate deaminase\n","truncated":false}}
%---
%[output:283577aa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:8bc75d47]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 121 ---\n","truncated":false}}
%---
%[output:56573bd7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 122 ---\n","truncated":false}}
%---
%[output:5e632d0d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 123 ---\n","truncated":false}}
%---
%[output:2367b67e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER152C or YGL202W or YJL060W\n","truncated":false}}
%---
%[output:4195e092]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:72fcbc73]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 124 ---\n","truncated":false}}
%---
%[output:571d92d2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR204W\n","truncated":false}}
%---
%[output:478ee891]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNR050C\n","truncated":false}}
%---
%[output:4a5fe891]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 125 ---\n","truncated":false}}
%---
%[output:5bce843c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 126 ---\n","truncated":false}}
%---
%[output:8182b134]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:2db086c5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:8d649804]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 127 ---\n","truncated":false}}
%---
%[output:8fcda8c3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0938 | prephenate dehydratase\n","truncated":false}}
%---
%[output:9c66502f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0015 | 2,5-diamino-6-ribosylamino-4(3H)-pyrimidinone 5'-phosphate reductase (NADPH)\n","truncated":false}}
%---
%[output:0fcd6821]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER152C or YGL202W or YJL060W\n","truncated":false}}
%---
%[output:70c7e91d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:541bf84c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:2044631e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YMR250W\n","truncated":false}}
%---
%[output:0bcf98cf]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 128 ---\n","truncated":false}}
%---
%[output:2a1f0d84]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 129 ---\n","truncated":false}}
%---
%[output:471da9e7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 130 ---\n","truncated":false}}
%---
%[output:378d728b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 131 ---\n","truncated":false}}
%---
%[output:46fc9dbd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6d10fd32]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:0d3c9f4a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 132 ---\n","truncated":false}}
%---
%[output:48528592]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR166C\n","truncated":false}}
%---
%[output:1af9645d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 133 ---\n","truncated":false}}
%---
%[output:9a06b37b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER152C or YGL202W or YJL060W\n","truncated":false}}
%---
%[output:348a4b7a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNR050C\n","truncated":false}}
%---
%[output:0d7ecf9f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:2245b1fa]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:6faafaf8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0015 | 2,5-diamino-6-ribosylamino-4(3H)-pyrimidinone 5'-phosphate reductase (NADPH)\n","truncated":false}}
%---
%[output:3d49a259]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 134 ---\n","truncated":false}}
%---
%[output:23a95609]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 135 ---\n","truncated":false}}
%---
%[output:27b732e5]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 136 ---\n","truncated":false}}
%---
%[output:5a978e21]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:40155f3a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:2a5a9fc3]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:7a61a555]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 137 ---\n","truncated":false}}
%---
%[output:733ea912]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:5672ddee]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER152C or YGL202W or YJL060W\n","truncated":false}}
%---
%[output:28d5d04e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 138 ---\n","truncated":false}}
%---
%[output:8e03de47]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 139 ---\n","truncated":false}}
%---
%[output:6e9cd27f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:99b3537d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:55650436]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:0135cfa5]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 140 ---\n","truncated":false}}
%---
%[output:2bd8b4ba]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0015 | 2,5-diamino-6-ribosylamino-4(3H)-pyrimidinone 5'-phosphate reductase (NADPH)\n","truncated":false}}
%---
%[output:1ff5b3ea]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0306 | 5,10-methylenetetrahydrofolate\n","truncated":false}}
%---
%[output:78dc74ef]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:42d805b0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:11ae98d3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 141 ---\n","truncated":false}}
%---
%[output:251bef3c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER152C or YGL202W or YJL060W\n","truncated":false}}
%---
%[output:8dc715b0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 142 ---\n","truncated":false}}
%---
%[output:0f531385]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 143 ---\n","truncated":false}}
%---
%[output:234536c8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:04ffd46a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 144 ---\n","truncated":false}}
%---
%[output:5c567304]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 145 ---\n","truncated":false}}
%---
%[output:3e9fbb43]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:72d000c2]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:9fa04f30]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 146 ---\n","truncated":false}}
%---
%[output:27af2963]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNL316C\n","truncated":false}}
%---
%[output:2b6c86ab]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 147 ---\n","truncated":false}}
%---
%[output:311198f1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:8e0e36de]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER152C or YGL202W or YJL060W\n","truncated":false}}
%---
%[output:99096ea8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:3da6859f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0015 | 2,5-diamino-6-ribosylamino-4(3H)-pyrimidinone 5'-phosphate reductase (NADPH)\n","truncated":false}}
%---
%[output:740695ea]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:7bd12cae]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:180cbeee]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 148 ---\n","truncated":false}}
%---
%[output:7e21387e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 149 ---\n","truncated":false}}
%---
%[output:48834740]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 150 ---\n","truncated":false}}
%---
%[output:4beedca0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:0be9b620]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:7005970f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 151 ---\n","truncated":false}}
%---
%[output:87e37c46]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:37300a7f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 152 ---\n","truncated":false}}
%---
%[output:76f7e102]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YHR063C\n","truncated":false}}
%---
%[output:2b33e2bf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:4512eaaa]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 153 ---\n","truncated":false}}
%---
%[output:04164d86]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 154 ---\n","truncated":false}}
%---
%[output:98877294]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0015 | 2,5-diamino-6-ribosylamino-4(3H)-pyrimidinone 5'-phosphate reductase (NADPH)\n","truncated":false}}
%---
%[output:3ef9b704]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:8e36c44b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:7fc3f50e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 155 ---\n","truncated":false}}
%---
%[output:37e2412e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:1535371d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1377 | prephenate\n","truncated":false}}
%---
%[output:85352222]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:92bfc140]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0476 | glutamine synthetase\n","truncated":false}}
%---
%[output:54166a80]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YHR063C\n","truncated":false}}
%---
%[output:2d97170b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 156 ---\n","truncated":false}}
%---
%[output:14207070]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0678 | L-aminoadipate-semialdehyde dehydrogenase (NADPH)\n","truncated":false}}
%---
%[output:78d38b2e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 157 ---\n","truncated":false}}
%---
%[output:9f854ea8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 158 ---\n","truncated":false}}
%---
%[output:3dcd5b93]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:7bde32f2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:5ab8b04c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 159 ---\n","truncated":false}}
%---
%[output:478b990f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 160 ---\n","truncated":false}}
%---
%[output:52808868]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0018 | 2-aminoadipate transaminase\n","truncated":false}}
%---
%[output:5ae71839]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 161 ---\n","truncated":false}}
%---
%[output:9189419d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 162 ---\n","truncated":false}}
%---
%[output:7666cbbe]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YHR063C\n","truncated":false}}
%---
%[output:5df7cb6e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 163 ---\n","truncated":false}}
%---
%[output:78c7940c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:3c3c921a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0678 | L-aminoadipate-semialdehyde dehydrogenase (NADPH)\n","truncated":false}}
%---
%[output:62ce6f62]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:3b25fff0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:3a957dab]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:748f34b6]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:495bcb17]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:63976968]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 164 ---\n","truncated":false}}
%---
%[output:5a7d35a3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 165 ---\n","truncated":false}}
%---
%[output:054d7b2b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 166 ---\n","truncated":false}}
%---
%[output:273afeb0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YHR063C\n","truncated":false}}
%---
%[output:8139b4a3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0018 | 2-aminoadipate transaminase\n","truncated":false}}
%---
%[output:70e8dded]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:845da39a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 167 ---\n","truncated":false}}
%---
%[output:47644492]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:3984af47]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 168 ---\n","truncated":false}}
%---
%[output:7f39ac57]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0732 | methylenetetrahydrofolate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:6f0fb172]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:6e3d93e2]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 169 ---\n","truncated":false}}
%---
%[output:5ff23c50]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 170 ---\n","truncated":false}}
%---
%[output:4f79a254]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:3c08ef07]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 171 ---\n","truncated":false}}
%---
%[output:82bb5015]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:711e02d0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YHR063C\n","truncated":false}}
%---
%[output:79a32d79]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 172 ---\n","truncated":false}}
%---
%[output:3c22f4dd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:202025a4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:273374e9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:9044c80c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0018 | 2-aminoadipate transaminase\n","truncated":false}}
%---
%[output:630c5390]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:08a743d1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 173 ---\n","truncated":false}}
%---
%[output:9d721458]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:547cfab5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YDR148C and YFL018C and YIL125W and YFR049W) or (YDR148C and YFL018C and YIL125W)\n","truncated":false}}
%---
%[output:73156168]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 174 ---\n","truncated":false}}
%---
%[output:729fa55b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 175 ---\n","truncated":false}}
%---
%[output:730e22aa]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL104C or YOR108W\n","truncated":false}}
%---
%[output:58444083]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:138cb531]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 176 ---\n","truncated":false}}
%---
%[output:1141e437]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 177 ---\n","truncated":false}}
%---
%[output:27f87fd4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:05b99f0e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 178 ---\n","truncated":false}}
%---
%[output:5bf75a91]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 179 ---\n","truncated":false}}
%---
%[output:33ba8656]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:2e814bf3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0939 | prephenate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:5a946bb4]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 180 ---\n","truncated":false}}
%---
%[output:6021a4ac]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0018 | 2-aminoadipate transaminase\n","truncated":false}}
%---
%[output:9af451b7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:673c322a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:705f987d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL104C or YOR108W\n","truncated":false}}
%---
%[output:682a1171]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:81b0fada]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 181 ---\n","truncated":false}}
%---
%[output:8302329b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:3be48109]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:09a71e6a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 182 ---\n","truncated":false}}
%---
%[output:69b9ce04]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 183 ---\n","truncated":false}}
%---
%[output:88052660]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 184 ---\n","truncated":false}}
%---
%[output:24b1f45a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:9601c806]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:83199545]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:2b51d49a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 185 ---\n","truncated":false}}
%---
%[output:3683cd23]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL104C or YOR108W\n","truncated":false}}
%---
%[output:6b48f107]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 186 ---\n","truncated":false}}
%---
%[output:139f79c3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:000386ba]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:24302ac7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0018 | 2-aminoadipate transaminase\n","truncated":false}}
%---
%[output:84dfadc1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 187 ---\n","truncated":false}}
%---
%[output:87aca4ae]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 188 ---\n","truncated":false}}
%---
%[output:3a24b888]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:45e4caa8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 189 ---\n","truncated":false}}
%---
%[output:7cd5425a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:2caadc2c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:915c3a8e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 190 ---\n","truncated":false}}
%---
%[output:7d95ae1e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:77860761]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL104C or YOR108W\n","truncated":false}}
%---
%[output:5ec285c4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:36e116f0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 191 ---\n","truncated":false}}
%---
%[output:2751c774]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0419 | ammonium\n","truncated":false}}
%---
%[output:98977b75]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:8f761848]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 192 ---\n","truncated":false}}
%---
%[output:9fc0fb5b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0018 | 2-aminoadipate transaminase\n","truncated":false}}
%---
%[output:310205b9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL131W or YDL182W\n","truncated":false}}
%---
%[output:595a949c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 193 ---\n","truncated":false}}
%---
%[output:09b3103b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:1e2eece8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 194 ---\n","truncated":false}}
%---
%[output:945966f7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:59e35237]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL104C or YOR108W\n","truncated":false}}
%---
%[output:7b3feb59]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 195 ---\n","truncated":false}}
%---
%[output:3cf2d74e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 196 ---\n","truncated":false}}
%---
%[output:47ce84b0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:80514bbf]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 197 ---\n","truncated":false}}
%---
%[output:209484ba]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 198 ---\n","truncated":false}}
%---
%[output:6fcfc9ea]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YMR250W\n","truncated":false}}
%---
%[output:09a8d452]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:5adc6733]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0938 | prephenate dehydratase\n","truncated":false}}
%---
%[output:3b7ea8b3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 199 ---\n","truncated":false}}
%---
%[output:83d40db2]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:9a098457]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0018 | 2-aminoadipate transaminase\n","truncated":false}}
%---
%[output:2536ed0b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL104C or YOR108W\n","truncated":false}}
%---
%[output:01131da1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:452c666a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:4dcdb699]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 200 ---\n","truncated":false}}
%---
%[output:64b03aab]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 201 ---\n","truncated":false}}
%---
%[output:47f23de8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:2e3669e7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YMR250W\n","truncated":false}}
%---
%[output:0ea97bcf]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 202 ---\n","truncated":false}}
%---
%[output:4899ffa2]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 203 ---\n","truncated":false}}
%---
%[output:0c35caaf]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 204 ---\n","truncated":false}}
%---
%[output:24daad86]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:1b93d727]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:0b026df8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL104C or YOR108W\n","truncated":false}}
%---
%[output:54cae95f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 205 ---\n","truncated":false}}
%---
%[output:858dce96]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:4d7771bc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:260655bb]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 206 ---\n","truncated":false}}
%---
%[output:5e3fca52]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YPR035W\n","truncated":false}}
%---
%[output:1a38401f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0019 | 2-dehydropantoate 2-reductase\n","truncated":false}}
%---
%[output:3d95d754]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0b4e31bb]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 207 ---\n","truncated":false}}
%---
%[output:7a8ae9cb]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 208 ---\n","truncated":false}}
%---
%[output:33dbaf23]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:007d59da]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 209 ---\n","truncated":false}}
%---
%[output:23d0594d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL104C or YOR108W\n","truncated":false}}
%---
%[output:26693c13]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:2fbac2eb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:08815d7c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 210 ---\n","truncated":false}}
%---
%[output:9be95330]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNR050C\n","truncated":false}}
%---
%[output:8176f9c1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:6fb16ea8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 211 ---\n","truncated":false}}
%---
%[output:71feb355]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 212 ---\n","truncated":false}}
%---
%[output:335e57fa]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:36d533f1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0019 | 2-dehydropantoate 2-reductase\n","truncated":false}}
%---
%[output:02c4464b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 213 ---\n","truncated":false}}
%---
%[output:4fbd17d5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:11d0bc76]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1d2ed9da]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL104C\n","truncated":false}}
%---
%[output:948cbbf3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 214 ---\n","truncated":false}}
%---
%[output:9f05f547]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNR050C\n","truncated":false}}
%---
%[output:0a02b2be]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 215 ---\n","truncated":false}}
%---
%[output:43cab3d2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:4b05a013]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:2f0a73a7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 216 ---\n","truncated":false}}
%---
%[output:027b3cad]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 217 ---\n","truncated":false}}
%---
%[output:946a362f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:9f1cd97b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:64a3fcd6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 218 ---\n","truncated":false}}
%---
%[output:39ede1e1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 219 ---\n","truncated":false}}
%---
%[output:85aabcb4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YHR208W or YJR148W or YGL202W or YHR137W\n","truncated":false}}
%---
%[output:3b0fff0c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YFL059W or YNL333W\n","truncated":false}}
%---
%[output:93fdc256]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0019 | 2-dehydropantoate 2-reductase\n","truncated":false}}
%---
%[output:070cfd19]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 220 ---\n","truncated":false}}
%---
%[output:2452a7fb]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:90f5d5c9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:076a57d1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:213de271]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:2ad06965]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:34d647ca]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 221 ---\n","truncated":false}}
%---
%[output:9f040b25]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 222 ---\n","truncated":false}}
%---
%[output:3a933269]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 223 ---\n","truncated":false}}
%---
%[output:96f70bb9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL038W or YOR347C\n","truncated":false}}
%---
%[output:9d4f69c5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YHR208W or YJR148W or YGL202W or YHR137W\n","truncated":false}}
%---
%[output:350d437b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 224 ---\n","truncated":false}}
%---
%[output:66e382ce]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 225 ---\n","truncated":false}}
%---
%[output:313fb5df]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR111C\n","truncated":false}}
%---
%[output:96de8d38]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 226 ---\n","truncated":false}}
%---
%[output:9c0b327a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0019 | 2-dehydropantoate 2-reductase\n","truncated":false}}
%---
%[output:978d4116]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:02fa9b6f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:12bbc1d8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 227 ---\n","truncated":false}}
%---
%[output:781b5302]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:389e05b5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL038W or YOR347C\n","truncated":false}}
%---
%[output:732e7a90]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 228 ---\n","truncated":false}}
%---
%[output:371af19b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YHR208W or YJR148W or YGL202W or YHR137W\n","truncated":false}}
%---
%[output:3fe74138]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOL052C\n","truncated":false}}
%---
%[output:0d5a2c42]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0306 | 5,10-methylenetetrahydrofolate\n","truncated":false}}
%---
%[output:764ef4ba]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:67f54263]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 229 ---\n","truncated":false}}
%---
%[output:11b49691]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 230 ---\n","truncated":false}}
%---
%[output:451d6d8c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 231 ---\n","truncated":false}}
%---
%[output:3e290cbb]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 232 ---\n","truncated":false}}
%---
%[output:2fbd38dc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:9440e4d1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL038W or YOR347C\n","truncated":false}}
%---
%[output:7fc28e82]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0019 | 2-dehydropantoate 2-reductase\n","truncated":false}}
%---
%[output:11a13fe5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:8d47466e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 233 ---\n","truncated":false}}
%---
%[output:42931f50]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YHR208W or YJR148W or YGL202W or YHR137W\n","truncated":false}}
%---
%[output:3a0743e8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:1bd6de1e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 234 ---\n","truncated":false}}
%---
%[output:722a5704]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 235 ---\n","truncated":false}}
%---
%[output:5ba43c97]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:87fea622]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:044dbf8f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 236 ---\n","truncated":false}}
%---
%[output:030a5af8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL038W or YOR347C\n","truncated":false}}
%---
%[output:6c8cfe1e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:59beae70]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0953 | L-2-aminoadipate\n","truncated":false}}
%---
%[output:54967435]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0832 | oxoglutarate dehydrogenase (lipoamide)\n","truncated":false}}
%---
%[output:2429dbb9]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 237 ---\n","truncated":false}}
%---
%[output:4c38ae8c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YHR208W or YJR148W or YGL202W or YHR137W\n","truncated":false}}
%---
%[output:24c45d42]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 238 ---\n","truncated":false}}
%---
%[output:8f9eae16]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0024 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:2c2000f8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:229c0d4a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 239 ---\n","truncated":false}}
%---
%[output:862f941b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 240 ---\n","truncated":false}}
%---
%[output:6240da1f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:79dfa8fe]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YJL060W\n","truncated":false}}
%---
%[output:027c93da]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 241 ---\n","truncated":false}}
%---
%[output:52b0dcb4]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 242 ---\n","truncated":false}}
%---
%[output:65d6f6fd]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 243 ---\n","truncated":false}}
%---
%[output:844ecb80]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:1c2f9965]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR487C\n","truncated":false}}
%---
%[output:838af062]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1377 | prephenate\n","truncated":false}}
%---
%[output:2b321227]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 244 ---\n","truncated":false}}
%---
%[output:02a0d147]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:8efe4676]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:804b0e82]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0953 | L-2-aminoadipate\n","truncated":false}}
%---
%[output:8f0cd67a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YJL060W\n","truncated":false}}
%---
%[output:090c880e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:0642c059]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0024 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:4818b790]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 245 ---\n","truncated":false}}
%---
%[output:8a399894]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 246 ---\n","truncated":false}}
%---
%[output:21039420]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 247 ---\n","truncated":false}}
%---
%[output:624bad03]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR249C or YDR035W\n","truncated":false}}
%---
%[output:9b78823b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:89e018ae]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:8deab526]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 248 ---\n","truncated":false}}
%---
%[output:1c7e9725]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 249 ---\n","truncated":false}}
%---
%[output:8847162a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YJL060W\n","truncated":false}}
%---
%[output:447c6b55]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 250 ---\n","truncated":false}}
%---
%[output:35fc1789]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:2bf1d9ef]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 251 ---\n","truncated":false}}
%---
%[output:0589bfbf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:2098807d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 252 ---\n","truncated":false}}
%---
%[output:6fa907a1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0024 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:00d88bff]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:0120bb78]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR249C or YDR035W\n","truncated":false}}
%---
%[output:71afa7d4]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 253 ---\n","truncated":false}}
%---
%[output:0ccee55c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:7345f361]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YER175C\n","truncated":false}}
%---
%[output:5daf201a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:378973be]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:52961f30]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 254 ---\n","truncated":false}}
%---
%[output:044cc5b4]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 255 ---\n","truncated":false}}
%---
%[output:75bd85b0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 256 ---\n","truncated":false}}
%---
%[output:54c653e1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:3a5566a9]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 257 ---\n","truncated":false}}
%---
%[output:6020ebfc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:3dcdcdc3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR249C or YDR035W\n","truncated":false}}
%---
%[output:76f526b1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:82ad03b9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR011W\n","truncated":false}}
%---
%[output:19ec6b1f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 258 ---\n","truncated":false}}
%---
%[output:345a17a5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:3b64aad0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0024 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:1f8ec22b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 259 ---\n","truncated":false}}
%---
%[output:5ba5b171]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:4618b100]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 260 ---\n","truncated":false}}
%---
%[output:34c49aa0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 261 ---\n","truncated":false}}
%---
%[output:8d76d34b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:85038a64]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 262 ---\n","truncated":false}}
%---
%[output:19ae58a6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR249C or YDR035W\n","truncated":false}}
%---
%[output:98493506]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:3970db3f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:0e70b173]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1838 | homocitrate synthase\n","truncated":false}}
%---
%[output:36306094]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 263 ---\n","truncated":false}}
%---
%[output:4b9795ac]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7d1e9ed0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:2ea76327]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 264 ---\n","truncated":false}}
%---
%[output:05198aaf]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:01800888]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0024 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:3f2d921d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 265 ---\n","truncated":false}}
%---
%[output:4695f474]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 266 ---\n","truncated":false}}
%---
%[output:00b3b220]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:31a9c62d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YLR231C\n","truncated":false}}
%---
%[output:5adde840]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 267 ---\n","truncated":false}}
%---
%[output:3f5b3c58]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:0506b13c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 268 ---\n","truncated":false}}
%---
%[output:65069091]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 269 ---\n","truncated":false}}
%---
%[output:33b9cdf0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:03709de8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:9b2ec8f7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:4f824830]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 270 ---\n","truncated":false}}
%---
%[output:4965bf67]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 271 ---\n","truncated":false}}
%---
%[output:7ff76ffa]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YJR139C\n","truncated":false}}
%---
%[output:67b272b9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:8e84d23d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:14520272]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YLR231C\n","truncated":false}}
%---
%[output:8ed571c6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0024 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:0c626fe5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:56e7ed0b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 272 ---\n","truncated":false}}
%---
%[output:19141129]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 273 ---\n","truncated":false}}
%---
%[output:55e383fb]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 274 ---\n","truncated":false}}
%---
%[output:348e8ecf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:22027e7d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:859a29cb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:5599eb24]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:00378585]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 275 ---\n","truncated":false}}
%---
%[output:316c594d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 276 ---\n","truncated":false}}
%---
%[output:15d64e0d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YLR231C\n","truncated":false}}
%---
%[output:9f99e26a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 277 ---\n","truncated":false}}
%---
%[output:3d79d93d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:36e44eca]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 278 ---\n","truncated":false}}
%---
%[output:323bf1a7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0024 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:8291cf97]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 279 ---\n","truncated":false}}
%---
%[output:148f3827]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3b5b02a3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:1d6e85b3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL061W\n","truncated":false}}
%---
%[output:9e21e67b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:4c967884]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 280 ---\n","truncated":false}}
%---
%[output:2ca0cc5b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0476 | glutamine synthetase\n","truncated":false}}
%---
%[output:5843df3c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCL018W\n","truncated":false}}
%---
%[output:15c3e391]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:85e4b623]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 281 ---\n","truncated":false}}
%---
%[output:54901dc5]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 282 ---\n","truncated":false}}
%---
%[output:72596dc1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 283 ---\n","truncated":false}}
%---
%[output:858ae2bf]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:8b1a0cb4]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 284 ---\n","truncated":false}}
%---
%[output:54b34184]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR011W\n","truncated":false}}
%---
%[output:31eb9918]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 285 ---\n","truncated":false}}
%---
%[output:885f09e0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0024 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:05554e00]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:2fcb32d7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:210ecfa3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:747867af]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 286 ---\n","truncated":false}}
%---
%[output:30a3a05b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:1522517b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:1ce853d4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:6c8ba680]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 287 ---\n","truncated":false}}
%---
%[output:27ea1218]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:32bdb8a7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 288 ---\n","truncated":false}}
%---
%[output:3ad8bbd5]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 289 ---\n","truncated":false}}
%---
%[output:1133a651]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:76d013a6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 290 ---\n","truncated":false}}
%---
%[output:8da79ab8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:7bfd58eb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:1f7a8948]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:9959c62d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 291 ---\n","truncated":false}}
%---
%[output:5e473573]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0025 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:7d03650c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 292 ---\n","truncated":false}}
%---
%[output:099e6741]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:300d4473]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:78e3325b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 293 ---\n","truncated":false}}
%---
%[output:20fb71ce]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:11725d22]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:9bcf1135]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 294 ---\n","truncated":false}}
%---
%[output:4e9fd9d2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YMR250W\n","truncated":false}}
%---
%[output:8fe64110]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 295 ---\n","truncated":false}}
%---
%[output:6ac09de6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:49b9bc67]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:71654441]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 296 ---\n","truncated":false}}
%---
%[output:9e2d3ae8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YJR139C\n","truncated":false}}
%---
%[output:571c50e7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 297 ---\n","truncated":false}}
%---
%[output:76acd579]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 298 ---\n","truncated":false}}
%---
%[output:78dce281]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0026 | 2-keto-4-methylthiobutyrate transamination\n","truncated":false}}
%---
%[output:0ae79aa6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4212 | D-ribulose 5-phosphate,D-glyceraldehyde 3-phosphate pyridoxal 5-phosphate-lyase\n","truncated":false}}
%---
%[output:2ce6cca2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNR050C\n","truncated":false}}
%---
%[output:233411b7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:1b898cb1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 299 ---\n","truncated":false}}
%---
%[output:18bd9a5e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 300 ---\n","truncated":false}}
%---
%[output:0ca6383f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:828a9c40]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 301 ---\n","truncated":false}}
%---
%[output:3f7b9b5e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:29995571]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:69e06d17]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:7dd4258d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:8b2d5afd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPR035W\n","truncated":false}}
%---
%[output:033087f2]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 302 ---\n","truncated":false}}
%---
%[output:1c059e6b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 303 ---\n","truncated":false}}
%---
%[output:237c5433]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 304 ---\n","truncated":false}}
%---
%[output:48f48955]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0962 | pyruvate kinase\n","truncated":false}}
%---
%[output:2057fefb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0026 | 2-keto-4-methylthiobutyrate transamination\n","truncated":false}}
%---
%[output:6677f733]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:4ab27dca]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL061W\n","truncated":false}}
%---
%[output:5ea6ac57]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 305 ---\n","truncated":false}}
%---
%[output:0208c275]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 306 ---\n","truncated":false}}
%---
%[output:5811dc69]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YMR250W\n","truncated":false}}
%---
%[output:11f991ee]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4226 | L-Alanine:2-oxoglutarate aminotransferase\n","truncated":false}}
%---
%[output:00155077]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 307 ---\n","truncated":false}}
%---
%[output:74429697]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:0b8bba99]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 308 ---\n","truncated":false}}
%---
%[output:8c685e90]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2693bec6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 309 ---\n","truncated":false}}
%---
%[output:5af6b794]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:6b3e8dff]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL038W or YOR347C\n","truncated":false}}
%---
%[output:05be9d26]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:261d2f55]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNR050C\n","truncated":false}}
%---
%[output:79d2818d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0962 | pyruvate kinase\n","truncated":false}}
%---
%[output:1f10c3f7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 310 ---\n","truncated":false}}
%---
%[output:63eafe62]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0026 | 2-keto-4-methylthiobutyrate transamination\n","truncated":false}}
%---
%[output:16547a72]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0145 | adenosylmethionine decarboxylase\n","truncated":false}}
%---
%[output:7e9f5eb5]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 311 ---\n","truncated":false}}
%---
%[output:81e95f02]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 312 ---\n","truncated":false}}
%---
%[output:5fa65f8d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 313 ---\n","truncated":false}}
%---
%[output:4b2686c3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 314 ---\n","truncated":false}}
%---
%[output:350709b5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL038W or YOR347C\n","truncated":false}}
%---
%[output:4422d326]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YMR250W\n","truncated":false}}
%---
%[output:79107700]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 315 ---\n","truncated":false}}
%---
%[output:44496883]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C or YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:50aecc58]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 316 ---\n","truncated":false}}
%---
%[output:02c2eb4e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2db01ad6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0962 | pyruvate kinase\n","truncated":false}}
%---
%[output:211fdbc9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:1f028dfa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:2795b153]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 317 ---\n","truncated":false}}
%---
%[output:7bfa8b5f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0026 | 2-keto-4-methylthiobutyrate transamination\n","truncated":false}}
%---
%[output:62a73a98]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNR050C\n","truncated":false}}
%---
%[output:961c8e17]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL038W or YOR347C\n","truncated":false}}
%---
%[output:1b0765c8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:8ce780c4]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 318 ---\n","truncated":false}}
%---
%[output:1bae90a4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C or YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:9e31f3aa]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 319 ---\n","truncated":false}}
%---
%[output:1160e836]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 320 ---\n","truncated":false}}
%---
%[output:78a8f627]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 321 ---\n","truncated":false}}
%---
%[output:226fafcc]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 322 ---\n","truncated":false}}
%---
%[output:040ae057]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0962 | pyruvate kinase\n","truncated":false}}
%---
%[output:6a80b85c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:26c26823]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:16b91550]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL038W or YOR347C\n","truncated":false}}
%---
%[output:6eef967b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 323 ---\n","truncated":false}}
%---
%[output:9c96292e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0182 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:2fd9cf95]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 324 ---\n","truncated":false}}
%---
%[output:5954ed84]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C or YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:836a99ee]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0026 | 2-keto-4-methylthiobutyrate transamination\n","truncated":false}}
%---
%[output:27404310]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 325 ---\n","truncated":false}}
%---
%[output:2a5c20c4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:755ca243]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:3d10d268]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:56dafdaf]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 326 ---\n","truncated":false}}
%---
%[output:96cfae11]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:7901937a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 327 ---\n","truncated":false}}
%---
%[output:90b48607]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:9ab3128e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4761 | 3-hydroxykynurenine aminotransferase (2-oxoglutarate)\n","truncated":false}}
%---
%[output:135ccc49]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 328 ---\n","truncated":false}}
%---
%[output:2e84a2c8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C or YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:11e7d7ac]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 329 ---\n","truncated":false}}
%---
%[output:0688f676]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7173e99e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 330 ---\n","truncated":false}}
%---
%[output:4d1f71db]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:142c1056]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0038 | 3,4-dihydroxy-2-butanone-4-phosphate synthase\n","truncated":false}}
%---
%[output:2b62cb6e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:7a0322b0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 331 ---\n","truncated":false}}
%---
%[output:23b229ce]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 332 ---\n","truncated":false}}
%---
%[output:23dd2c31]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 333 ---\n","truncated":false}}
%---
%[output:282f6159]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:02db2a31]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 334 ---\n","truncated":false}}
%---
%[output:44fa83f7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:97c5066c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C or YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:30802ac1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4761 | 3-hydroxykynurenine aminotransferase (2-oxoglutarate)\n","truncated":false}}
%---
%[output:8dc15161]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:8dd3ee05]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:0c1a7d85]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:20f7be31]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 335 ---\n","truncated":false}}
%---
%[output:481897eb]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 336 ---\n","truncated":false}}
%---
%[output:8e8dfbae]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 337 ---\n","truncated":false}}
%---
%[output:9609eb5d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0042 | 3-deoxy-D-arabino-heptulosonate 7-phosphate synthetase\n","truncated":false}}
%---
%[output:57819ba6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:89edef7c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:7a9d3043]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:5f19e928]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C or YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:8a6cf90a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 338 ---\n","truncated":false}}
%---
%[output:506a8dff]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 339 ---\n","truncated":false}}
%---
%[output:2eec927d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:05949f02]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4761 | 3-hydroxykynurenine aminotransferase (2-oxoglutarate)\n","truncated":false}}
%---
%[output:17e01612]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 340 ---\n","truncated":false}}
%---
%[output:1bad0bca]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 341 ---\n","truncated":false}}
%---
%[output:94afb47f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:9ed31a80]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 342 ---\n","truncated":false}}
%---
%[output:3e4d24a5]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:6a77370a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 343 ---\n","truncated":false}}
%---
%[output:1a66b34f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:795481a1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR127W\n","truncated":false}}
%---
%[output:6c5322e0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:9ee127f6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:53818287]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0042 | 3-deoxy-D-arabino-heptulosonate 7-phosphate synthetase\n","truncated":false}}
%---
%[output:1dfc9a15]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 344 ---\n","truncated":false}}
%---
%[output:850003d5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:1d3f5f17]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0419 | ammonium\n","truncated":false}}
%---
%[output:989e660d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0059 | 3-isopropylmalate 3-methyltransferase\n","truncated":false}}
%---
%[output:09d166fe]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 345 ---\n","truncated":false}}
%---
%[output:98e412e3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 346 ---\n","truncated":false}}
%---
%[output:2ce5a054]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 347 ---\n","truncated":false}}
%---
%[output:207adbd3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 348 ---\n","truncated":false}}
%---
%[output:34cd3712]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR127W\n","truncated":false}}
%---
%[output:104d6523]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:16bf7568]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 349 ---\n","truncated":false}}
%---
%[output:65915348]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNL141W\n","truncated":false}}
%---
%[output:83c3cb38]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:5cc41b5d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 350 ---\n","truncated":false}}
%---
%[output:949f7f4c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:8eff27c0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0042 | 3-deoxy-D-arabino-heptulosonate 7-phosphate synthetase\n","truncated":false}}
%---
%[output:44710edc]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 351 ---\n","truncated":false}}
%---
%[output:9e9bbdaa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:017bd402]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 352 ---\n","truncated":false}}
%---
%[output:7d72e72d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:53a2560e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:4909b378]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNL141W\n","truncated":false}}
%---
%[output:0722df1f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR127W\n","truncated":false}}
%---
%[output:549ad6c7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:70cc05a1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 353 ---\n","truncated":false}}
%---
%[output:12a0af94]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:01bd9336]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 354 ---\n","truncated":false}}
%---
%[output:93b76927]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 355 ---\n","truncated":false}}
%---
%[output:9276cf8a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 356 ---\n","truncated":false}}
%---
%[output:866fc1ab]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 357 ---\n","truncated":false}}
%---
%[output:9d5523b8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0042 | 3-deoxy-D-arabino-heptulosonate 7-phosphate synthetase\n","truncated":false}}
%---
%[output:44f7cd07]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNL141W\n","truncated":false}}
%---
%[output:4ec682ed]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR011W\n","truncated":false}}
%---
%[output:17d08faa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:8f0aa1b2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR127W\n","truncated":false}}
%---
%[output:16e04d58]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 358 ---\n","truncated":false}}
%---
%[output:5ab15f68]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:6ebc0b3b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 359 ---\n","truncated":false}}
%---
%[output:2ea0d6cf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:4ba8096f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6f944047]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 360 ---\n","truncated":false}}
%---
%[output:77555e26]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNL141W\n","truncated":false}}
%---
%[output:9936dae0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 361 ---\n","truncated":false}}
%---
%[output:5bce85d6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:24778d34]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:118a6e79]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 362 ---\n","truncated":false}}
%---
%[output:1616ddc3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNR033W\n","truncated":false}}
%---
%[output:1652a8b1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 363 ---\n","truncated":false}}
%---
%[output:09ac476b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:02915393]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0045 | 3-hydroxy-L-kynurenine hydrolase\n","truncated":false}}
%---
%[output:7c5776c6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 364 ---\n","truncated":false}}
%---
%[output:46e1e9ac]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR018C\n","truncated":false}}
%---
%[output:0ac8c280]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:89baa316]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 365 ---\n","truncated":false}}
%---
%[output:8e497efc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:2a4ffe6e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 366 ---\n","truncated":false}}
%---
%[output:21fd2028]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:2bfede77]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:77efd6bc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNR033W\n","truncated":false}}
%---
%[output:319447d7]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 367 ---\n","truncated":false}}
%---
%[output:780c4a95]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 368 ---\n","truncated":false}}
%---
%[output:0b563091]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR209C\n","truncated":false}}
%---
%[output:6b7c2229]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 369 ---\n","truncated":false}}
%---
%[output:38ed6bec]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:1e5190ab]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 370 ---\n","truncated":false}}
%---
%[output:2105e2a6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YJR139C\n","truncated":false}}
%---
%[output:8096f4df]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:1f2b4757]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0045 | 3-hydroxy-L-kynurenine hydrolase\n","truncated":false}}
%---
%[output:75a7a70e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:86b18adf]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 371 ---\n","truncated":false}}
%---
%[output:5ea1b95d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR019W\n","truncated":false}}
%---
%[output:957d22b4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL103C\n","truncated":false}}
%---
%[output:092dbd8c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 372 ---\n","truncated":false}}
%---
%[output:98fab286]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 373 ---\n","truncated":false}}
%---
%[output:9599c3d3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 374 ---\n","truncated":false}}
%---
%[output:00dc73b1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1d921896]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:8c2552b9]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:2fe6f433]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:24962699]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:578ad184]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 375 ---\n","truncated":false}}
%---
%[output:2cee2f7a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YHL012W or YKL035W\n","truncated":false}}
%---
%[output:2fe5e793]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 376 ---\n","truncated":false}}
%---
%[output:2adb8aba]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR019W\n","truncated":false}}
%---
%[output:3c6a4f64]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0045 | 3-hydroxy-L-kynurenine hydrolase\n","truncated":false}}
%---
%[output:22d22eba]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 377 ---\n","truncated":false}}
%---
%[output:2e5bcdf3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 378 ---\n","truncated":false}}
%---
%[output:3ee3b600]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YAL061W\n","truncated":false}}
%---
%[output:032282c0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 379 ---\n","truncated":false}}
%---
%[output:02a8ecda]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:42f622d5]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 380 ---\n","truncated":false}}
%---
%[output:3f5cef0d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR192W\n","truncated":false}}
%---
%[output:4f304123]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 381 ---\n","truncated":false}}
%---
%[output:6fc579c3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:1a6362fe]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:120fd3b2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR019W\n","truncated":false}}
%---
%[output:49b554f5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5e40e727]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 382 ---\n","truncated":false}}
%---
%[output:3603d841]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDR017C\n","truncated":false}}
%---
%[output:7866ac24]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:6c58a366]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0061 | 3-isopropylmalate dehydrogenase\n","truncated":false}}
%---
%[output:336dfc54]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR192W\n","truncated":false}}
%---
%[output:12faa176]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 383 ---\n","truncated":false}}
%---
%[output:9b2bbe5d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 384 ---\n","truncated":false}}
%---
%[output:25e5e4d6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 385 ---\n","truncated":false}}
%---
%[output:99418a8d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 386 ---\n","truncated":false}}
%---
%[output:8fc21a4a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:5266a2a0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR019W\n","truncated":false}}
%---
%[output:0919b219]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 387 ---\n","truncated":false}}
%---
%[output:40bf1ba0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:833b07a3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:357ab227]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:2eb01979]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 388 ---\n","truncated":false}}
%---
%[output:8206a718]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:17305af5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:133a26ce]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 389 ---\n","truncated":false}}
%---
%[output:77e9283f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0062 | 3-methyl-2-oxobutanoate decarboxylase\n","truncated":false}}
%---
%[output:30a5a747]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 390 ---\n","truncated":false}}
%---
%[output:7348f770]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:771a6871]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR019W\n","truncated":false}}
%---
%[output:775ffa98]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:28c97839]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:4c7366fd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:4c7a7e8d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 391 ---\n","truncated":false}}
%---
%[output:95579a84]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 392 ---\n","truncated":false}}
%---
%[output:937a73df]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:3c59d843]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 393 ---\n","truncated":false}}
%---
%[output:776cce4e]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 394 ---\n","truncated":false}}
%---
%[output:2f9150b2]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 395 ---\n","truncated":false}}
%---
%[output:051092cc]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 396 ---\n","truncated":false}}
%---
%[output:8ea16393]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:7b52ec80]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:096f70cf]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDR017C\n","truncated":false}}
%---
%[output:4df45d68]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C\n","truncated":false}}
%---
%[output:08069b8d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0062 | 3-methyl-2-oxobutanoate decarboxylase\n","truncated":false}}
%---
%[output:43681c19]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:5770aa16]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 397 ---\n","truncated":false}}
%---
%[output:1e27c755]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0376 | acetyl-CoA\n","truncated":false}}
%---
%[output:7bcfb166]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 398 ---\n","truncated":false}}
%---
%[output:7e6c0b84]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:4a88527e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:1f524b9a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:5f3d7ae6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 399 ---\n","truncated":false}}
%---
%[output:2332b963]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:0f1311d0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 400 ---\n","truncated":false}}
%---
%[output:72800c69]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C\n","truncated":false}}
%---
%[output:760f65b1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 401 ---\n","truncated":false}}
%---
%[output:7a82b0af]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 402 ---\n","truncated":false}}
%---
%[output:6c68bd5f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:2c59e2fb]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 403 ---\n","truncated":false}}
%---
%[output:4c0fd898]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0062 | 3-methyl-2-oxobutanoate decarboxylase\n","truncated":false}}
%---
%[output:40f00c1d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:6a86667d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:2278d136]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 404 ---\n","truncated":false}}
%---
%[output:69ddfa1e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:0c53e0cf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:5109bf7a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 405 ---\n","truncated":false}}
%---
%[output:9a67d2ba]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C\n","truncated":false}}
%---
%[output:842d6512]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 406 ---\n","truncated":false}}
%---
%[output:218eea8a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:31127ad6]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:882a4c5c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:3941f27f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:5fd3b80a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 407 ---\n","truncated":false}}
%---
%[output:64ae93eb]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 408 ---\n","truncated":false}}
%---
%[output:9ac2d027]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:5d3be5ea]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 409 ---\n","truncated":false}}
%---
%[output:568d3ce9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0062 | 3-methyl-2-oxobutanoate decarboxylase\n","truncated":false}}
%---
%[output:40d89321]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 410 ---\n","truncated":false}}
%---
%[output:5885d9d5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C\n","truncated":false}}
%---
%[output:41da4fe7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:9819a2c6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:459fffaa]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:2b1b62dd]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 411 ---\n","truncated":false}}
%---
%[output:781cf433]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 412 ---\n","truncated":false}}
%---
%[output:0d22f074]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0476 | glutamine synthetase\n","truncated":false}}
%---
%[output:8f197131]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:3fbe1002]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 413 ---\n","truncated":false}}
%---
%[output:23e2b252]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 414 ---\n","truncated":false}}
%---
%[output:12adb112]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 415 ---\n","truncated":false}}
%---
%[output:8f5f6cd7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:455806f9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:0f5a2cd0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C\n","truncated":false}}
%---
%[output:5db47e85]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:2d009b09]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0062 | 3-methyl-2-oxobutanoate decarboxylase\n","truncated":false}}
%---
%[output:2c5d8580]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:5c31c3fb]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 416 ---\n","truncated":false}}
%---
%[output:89b193de]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:5c7ebe62]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 417 ---\n","truncated":false}}
%---
%[output:26f25dea]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:481bc0a9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:85991f6d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:6531f418]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 418 ---\n","truncated":false}}
%---
%[output:133458a8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 419 ---\n","truncated":false}}
%---
%[output:5f826a31]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL080C\n","truncated":false}}
%---
%[output:6c299f58]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 420 ---\n","truncated":false}}
%---
%[output:8d152646]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 421 ---\n","truncated":false}}
%---
%[output:0b107991]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 422 ---\n","truncated":false}}
%---
%[output:1433932e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:08905c89]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 423 ---\n","truncated":false}}
%---
%[output:0d969fcd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0962 | pyruvate kinase\n","truncated":false}}
%---
%[output:9dbbc68b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:171832e5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0062 | 3-methyl-2-oxobutanoate decarboxylase\n","truncated":false}}
%---
%[output:9f1b43f0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:04f9fbb6]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:4dc81ee0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 424 ---\n","truncated":false}}
%---
%[output:916f6269]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YLR410W or YOR163W\n","truncated":false}}
%---
%[output:838866da]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:9561625e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1416 | S-adenosyl-L-methionine\n","truncated":false}}
%---
%[output:81f1b08c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:1c2f23cf]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 425 ---\n","truncated":false}}
%---
%[output:221c5f3c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:174b7e3a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 426 ---\n","truncated":false}}
%---
%[output:21b07c49]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 427 ---\n","truncated":false}}
%---
%[output:42a67e1d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 428 ---\n","truncated":false}}
%---
%[output:1e6609c6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0962 | pyruvate kinase\n","truncated":false}}
%---
%[output:8a6a0f01]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:6289aaa3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 429 ---\n","truncated":false}}
%---
%[output:67751617]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YLR017W\n","truncated":false}}
%---
%[output:497ff4da]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0064 | 3-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:6fc3155d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:9be71a1b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:8e98dea8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 430 ---\n","truncated":false}}
%---
%[output:6e089475]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 431 ---\n","truncated":false}}
%---
%[output:11fd5494]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:0fdcc375]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 432 ---\n","truncated":false}}
%---
%[output:4c4e4284]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:01754cc2]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 433 ---\n","truncated":false}}
%---
%[output:0f09f978]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:4ca6a836]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:8cb79b28]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YLR017W\n","truncated":false}}
%---
%[output:0892d630]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:4aadeea7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:12f698ad]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0962 | pyruvate kinase\n","truncated":false}}
%---
%[output:9f52a311]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 434 ---\n","truncated":false}}
%---
%[output:7d1b6684]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 435 ---\n","truncated":false}}
%---
%[output:8b589b0b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0064 | 3-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:8b5b8f49]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 436 ---\n","truncated":false}}
%---
%[output:5ef6b4d8]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 437 ---\n","truncated":false}}
%---
%[output:3f537d17]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 438 ---\n","truncated":false}}
%---
%[output:2a579d89]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 439 ---\n","truncated":false}}
%---
%[output:5d4d5748]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:1f7d5b3c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:4dd1b8a6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:90d5f09c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YLR017W\n","truncated":false}}
%---
%[output:7267daf0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:6ef21499]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:5c0fd54b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0962 | pyruvate kinase\n","truncated":false}}
%---
%[output:943f8d29]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 440 ---\n","truncated":false}}
%---
%[output:63575f60]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 441 ---\n","truncated":false}}
%---
%[output:0b32938f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 442 ---\n","truncated":false}}
%---
%[output:76aa1fb1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:41e5354d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0064 | 3-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:6871322b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:00fcb936]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:05f5471d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 443 ---\n","truncated":false}}
%---
%[output:52b53f42]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YLR017W\n","truncated":false}}
%---
%[output:64d96122]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 444 ---\n","truncated":false}}
%---
%[output:1444d51a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:1fa3274c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 445 ---\n","truncated":false}}
%---
%[output:1f15f735]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 446 ---\n","truncated":false}}
%---
%[output:642471da]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR209C\n","truncated":false}}
%---
%[output:3237ff38]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:92d6b4ca]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 447 ---\n","truncated":false}}
%---
%[output:075f2968]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:21cfa3ca]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:515b01b6]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:3689ad23]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 448 ---\n","truncated":false}}
%---
%[output:0a6c7bbb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:05384ee2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0064 | 3-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:54818bec]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 449 ---\n","truncated":false}}
%---
%[output:9ddefcba]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YFR047C\n","truncated":false}}
%---
%[output:393e16c7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:220ba4ad]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 450 ---\n","truncated":false}}
%---
%[output:7495be47]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 451 ---\n","truncated":false}}
%---
%[output:21e39684]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0577 | D-ribulose 5-phosphate\n","truncated":false}}
%---
%[output:79930d8c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:95d957d1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:884bcca6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 452 ---\n","truncated":false}}
%---
%[output:6da8979f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 453 ---\n","truncated":false}}
%---
%[output:76d56cb0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:686712dc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YMR250W\n","truncated":false}}
%---
%[output:19a880b3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 454 ---\n","truncated":false}}
%---
%[output:76f04d96]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:98852bd5]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 455 ---\n","truncated":false}}
%---
%[output:9aa967a7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:8f1d328d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0064 | 3-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:0b8651eb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:0e5ed95d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:5e921ffb]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 456 ---\n","truncated":false}}
%---
%[output:01b32f91]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 457 ---\n","truncated":false}}
%---
%[output:364f0131]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNR050C\n","truncated":false}}
%---
%[output:7183ac6e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:0ae95044]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:4e4622f1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 458 ---\n","truncated":false}}
%---
%[output:953901c1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 459 ---\n","truncated":false}}
%---
%[output:7b9c3205]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 460 ---\n","truncated":false}}
%---
%[output:4732a548]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1360 | phosphoenolpyruvate\n","truncated":false}}
%---
%[output:30c85285]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:9b26eae3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:57cb16b1]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 461 ---\n","truncated":false}}
%---
%[output:553cc7b6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:4ec5a83c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:805128fa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0064 | 3-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:03be723b]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 462 ---\n","truncated":false}}
%---
%[output:4ceb9254]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:3ddfe554]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 463 ---\n","truncated":false}}
%---
%[output:3cf5fc30]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:8593dec4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0955 | L-alanine\n","truncated":false}}
%---
%[output:8d925d84]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:332ff50c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 464 ---\n","truncated":false}}
%---
%[output:60dd97e7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:36283bc6]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 465 ---\n","truncated":false}}
%---
%[output:5c22cc6e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:68e1387f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 466 ---\n","truncated":false}}
%---
%[output:0c619111]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 467 ---\n","truncated":false}}
%---
%[output:780f2186]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 468 ---\n","truncated":false}}
%---
%[output:95162f29]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:08a0e8bd]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 469 ---\n","truncated":false}}
%---
%[output:86297da8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0065 | 3-phosphoshikimate 1-carboxyvinyltransferase\n","truncated":false}}
%---
%[output:60da3345]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:0bb498e2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:7d5308ad]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:0e10988b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:1dcb4f15]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1360 | phosphoenolpyruvate\n","truncated":false}}
%---
%[output:3ea3fa0a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 470 ---\n","truncated":false}}
%---
%[output:2814787b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:0d3f4c24]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 471 ---\n","truncated":false}}
%---
%[output:044a0c3f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:2d9181cb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:0b24a315]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 472 ---\n","truncated":false}}
%---
%[output:20198b64]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:47147243]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:4b16be22]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 473 ---\n","truncated":false}}
%---
%[output:62bd98ec]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 474 ---\n","truncated":false}}
%---
%[output:243378da]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 475 ---\n","truncated":false}}
%---
%[output:5b872a36]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0065 | 3-phosphoshikimate 1-carboxyvinyltransferase\n","truncated":false}}
%---
%[output:7c472d58]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:583a9124]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 476 ---\n","truncated":false}}
%---
%[output:9076c365]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:2e247695]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:5995a79e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR192C or YJL052W or YJR009C\n","truncated":false}}
%---
%[output:2e4232c9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:9c30cace]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:87e2f1ba]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 477 ---\n","truncated":false}}
%---
%[output:6c095e46]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 478 ---\n","truncated":false}}
%---
%[output:668be0eb]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1360 | phosphoenolpyruvate\n","truncated":false}}
%---
%[output:6968500c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 479 ---\n","truncated":false}}
%---
%[output:404cb5f4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0937 | isobutyraldehyde\n","truncated":false}}
%---
%[output:51988360]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 480 ---\n","truncated":false}}
%---
%[output:330eb5d3]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 481 ---\n","truncated":false}}
%---
%[output:81f76006]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR192C or YJL052W or YJR009C\n","truncated":false}}
%---
%[output:87826620]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:28a33a38]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:708d27fd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:413c03a5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:4a6ef407]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0065 | 3-phosphoshikimate 1-carboxyvinyltransferase\n","truncated":false}}
%---
%[output:6255bf45]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 482 ---\n","truncated":false}}
%---
%[output:472f9b40]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 483 ---\n","truncated":false}}
%---
%[output:1c90138a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:810040b0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 484 ---\n","truncated":false}}
%---
%[output:3e0538dc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:9b7472e0]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 485 ---\n","truncated":false}}
%---
%[output:1b3e1aa2]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 486 ---\n","truncated":false}}
%---
%[output:43dda632]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 487 ---\n","truncated":false}}
%---
%[output:80f6ba13]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:7f0ac505]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1360 | phosphoenolpyruvate\n","truncated":false}}
%---
%[output:90eb24e8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:997d31cf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:8f837eea]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:335af21f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0937 | isobutyraldehyde\n","truncated":false}}
%---
%[output:35d713b3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR166C\n","truncated":false}}
%---
%[output:1d02f969]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0065 | 3-phosphoshikimate 1-carboxyvinyltransferase\n","truncated":false}}
%---
%[output:69ebd86a]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 488 ---\n","truncated":false}}
%---
%[output:73a80d0f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 489 ---\n","truncated":false}}
%---
%[output:54445450]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 490 ---\n","truncated":false}}
%---
%[output:12e17974]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:8d8575a0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:08a8ddda]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 491 ---\n","truncated":false}}
%---
%[output:6cdd017a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:425b7e9d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 492 ---\n","truncated":false}}
%---
%[output:1b997a2e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:2effec92]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:7ddbbc7c]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 493 ---\n","truncated":false}}
%---
%[output:877d76c4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:174914bc]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 494 ---\n","truncated":false}}
%---
%[output:8cd72412]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 495 ---\n","truncated":false}}
%---
%[output:08a5430f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YMR300C\n","truncated":false}}
%---
%[output:10b608d8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0066 | 4-amino-4-deoxychorismate synthase\n","truncated":false}}
%---
%[output:84b0e411]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 496 ---\n","truncated":false}}
%---
%[output:707f3b4c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0937 | isobutyraldehyde\n","truncated":false}}
%---
%[output:6be53022]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:98665350]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:23e2da42]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0222 | 3-hydroxy-L-kynurenine\n","truncated":false}}
%---
%[output:65155b44]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 497 ---\n","truncated":false}}
%---
%[output:107b7cac]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0459 | galactose-1-phosphate uridylyltransferase\n","truncated":false}}
%---
%[output:5267237b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:15efc94f]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 498 ---\n","truncated":false}}
%---
%[output:73a4a6c7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YMR300C\n","truncated":false}}
%---
%[output:534174dd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:4d5e3bda]
%   data: {"dataType":"text","outputData":{"text":"\n--- Counter candidate 499 ---\n","truncated":false}}
%---
%[output:33dbcd67]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR221C and YER178W and YFL018C and YGR193C and YNL071W\n","truncated":false}}
%---
%[output:6495fb90]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0222 | 3-hydroxy-L-kynurenine\n","truncated":false}}
%---
%[output:2600ac11]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0222 | 3-hydroxy-L-kynurenine\n","truncated":false}}
%---
%[output:0a298316]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:40837d18]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0066 | 4-amino-4-deoxychorismate synthase\n","truncated":false}}
%---
%[output:542a99b2]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0009 | (2R,3S)-3-isopropylmalate\n","truncated":false}}
%---
%[output:38e7ad6d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:37b96f83]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0973 | ribonucleoside-triphosphate reductase (UTP)\n","truncated":false}}
%---
%[output:826dc833]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:5676361d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:319e763e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR221C and YER178W and YFL018C and YGR193C and YNL071W\n","truncated":false}}
%---
%[output:22d06d92]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0937 | isobutyraldehyde\n","truncated":false}}
%---
%[output:7fbeecac]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:1df0e12c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:4d9057d8]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:08551844]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:523b4e35]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:08acdf06]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:07afa78c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:6796654e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:41c95c35]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL054C\n","truncated":false}}
%---
%[output:22975d00]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0068 | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:8d642964]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1069 | UDP-N-acetylglucosamine diphosphorylase\n","truncated":false}}
%---
%[output:43d4a80b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:9dd3adbc]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:3e1d8caa]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:73ac57eb]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:40547158]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:964cf3e3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:076e2ccc]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:497b6fc9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:25b62f4f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:436cc9ca]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0937 | isobutyraldehyde\n","truncated":false}}
%---
%[output:21d8045d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:4a8eb149]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:5725ac1f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1084 | UTP-glucose-1-phosphate uridylyltransferase\n","truncated":false}}
%---
%[output:8b81c326]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1360 | phosphoenolpyruvate\n","truncated":false}}
%---
%[output:46beff76]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0068 | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:69efc685]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1360 | phosphoenolpyruvate\n","truncated":false}}
%---
%[output:88f2347c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:4abd6434]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:2083147d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1360 | phosphoenolpyruvate\n","truncated":false}}
%---
%[output:2de287d1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:190f5cd8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1360 | phosphoenolpyruvate\n","truncated":false}}
%---
%[output:983fd48c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:3c80f1a9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:860ea298]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:9ba4b911]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:516f6eb8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1275 | UTP transport\n","truncated":false}}
%---
%[output:444bccdd]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:919cf212]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:6deb9326]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8f32ead8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:758fa92d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:43fdb9ed]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0937 | isobutyraldehyde\n","truncated":false}}
%---
%[output:3deeec57]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0068 | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:315a7e33]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:34b9ea6b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:258a5fe9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0074 | 4PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:596e285d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:1b7becac]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:29da7539]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:372ef0bf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1276 | UTP\/UMP antiport\n","truncated":false}}
%---
%[output:9537147e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:4fcebf17]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:77c3d769]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:23a66deb]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:256d806b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:286838a7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:45466aa0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:3c8a3553]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:16160485]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0068 | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:77cfd510]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:1b123b4c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:8ff8636f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0166 | 2-methylbutanal\n","truncated":false}}
%---
%[output:2ff81612]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YFL060C or YNL334C\n","truncated":false}}
%---
%[output:139be4e6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:29a7f084]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:27a4dbd6]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0279 | 4-diphospho-1D-myo-inositol pentakisphosphate\n","truncated":false}}
%---
%[output:219edde0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:74f65834]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:27ca9601]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:976d7355]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:461eb6c7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:00d6f576]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:7ce1d45b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:54b75497]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YFL059W or YNL333W\n","truncated":false}}
%---
%[output:0753b3cb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0068 | 4-aminobutyrate transaminase\n","truncated":false}}
%---
%[output:692f98e2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:90b943c1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:806476d2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:5ba9bb8e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:37e403bb]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:28f22c9d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:0c008d96]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:81084beb]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0166 | 2-methylbutanal\n","truncated":false}}
%---
%[output:655f29d4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:2967e275]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7c402e6e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:7d973e83]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:8b4ef88f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:46eaeaee]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:38274cba]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:33fc7e28]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:58da7ff4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:17954f18]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:69d08b6c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0072 | 4-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:8fc7c3c8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:4030a1e1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNL316C\n","truncated":false}}
%---
%[output:2a15e637]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1beeff45]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:2faf9595]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:908f7b41]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:03b8749f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:04d17245]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:24176cf9]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0166 | 2-methylbutanal\n","truncated":false}}
%---
%[output:392c7b01]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:012e3820]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:7bc2e00a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:5cc0fa92]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:943a5b40]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:50bdabc9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:8bcf63c6]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:3ba1e1c0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0072 | 4-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:3feaa8fe]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:73e88950]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:29786336]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:87f0b6a7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:43f6bc88]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:5e904ec9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:5d97ceb1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER073W or YOR374W\n","truncated":false}}
%---
%[output:0452e2ce]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:58478739]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:3d5c32e6]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:94c89ef9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:5411a8ed]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0306 | CTP synthase (glutamine)\n","truncated":false}}
%---
%[output:125bdce3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0166 | 2-methylbutanal\n","truncated":false}}
%---
%[output:239407ba]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:7d934f0c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:48f95032]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0072 | 4-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:483825c1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YDR148C and YFL018C and YIL125W and YFR049W) or (YDR148C and YFL018C and YIL125W)\n","truncated":false}}
%---
%[output:7ae98b04]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:66bc6551]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = (YCL009C and YMR108W) or YMR108W\n","truncated":false}}
%---
%[output:4d78d608]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:28c9364a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:6ae9e99f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:3511014f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:06f2243e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:0a71171b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:037dd31b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0306 | CTP synthase (glutamine)\n","truncated":false}}
%---
%[output:6133c2d6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:100080d7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:0231fcff]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:4d23e3b3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = (YCL009C and YMR108W) or YMR108W\n","truncated":false}}
%---
%[output:0f50bba1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:5aa43ac1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:70bf81f7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0072 | 4-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:611cc880]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0166 | 2-methylbutanal\n","truncated":false}}
%---
%[output:1c987f05]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:7764ee68]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:9dda78e5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:6a31871b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:405baf28]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:514dcba3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:0f13f662]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBL015W\n","truncated":false}}
%---
%[output:261f5f29]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0419 | ammonium\n","truncated":false}}
%---
%[output:068b4afe]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0306 | CTP synthase (glutamine)\n","truncated":false}}
%---
%[output:6c396e26]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:55a23363]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0799 | H+\n","truncated":false}}
%---
%[output:600e874b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0799 | H+\n","truncated":false}}
%---
%[output:0d0f5177]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:0121f727]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1464 | succinyl-CoA\n","truncated":false}}
%---
%[output:4b8fffd9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:70f55f0b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0072 | 4-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:751a463a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0318 | 5-diphospho-1D-myo-inositol pentakisphosphate\n","truncated":false}}
%---
%[output:8e85394b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:3e0b1d5f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:1ce17b05]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0318 | 5-diphospho-1D-myo-inositol pentakisphosphate\n","truncated":false}}
%---
%[output:7a7c2bb9]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0166 | 2-methylbutanal\n","truncated":false}}
%---
%[output:20487279]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0318 | 5-diphospho-1D-myo-inositol pentakisphosphate\n","truncated":false}}
%---
%[output:6f6c8f47]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:652ea116]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0306 | CTP synthase (glutamine)\n","truncated":false}}
%---
%[output:3fa929de]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0318 | 5-diphospho-1D-myo-inositol pentakisphosphate\n","truncated":false}}
%---
%[output:45918ad4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:28f7dd98]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0318 | 5-diphospho-1D-myo-inositol pentakisphosphate\n","truncated":false}}
%---
%[output:8cf1ae62]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:843d1706]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:8f2d1a4d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0318 | 5-diphospho-1D-myo-inositol pentakisphosphate\n","truncated":false}}
%---
%[output:0525186b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:9cdb9b35]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0318 | 5-diphospho-1D-myo-inositol pentakisphosphate\n","truncated":false}}
%---
%[output:402929b1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:620c2b5d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0072 | 4-methyl-2-oxopentanoate decarboxylase\n","truncated":false}}
%---
%[output:7222ddf7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0318 | 5-diphospho-1D-myo-inositol pentakisphosphate\n","truncated":false}}
%---
%[output:60c1be90]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0557 | D-fructose 6-phosphate\n","truncated":false}}
%---
%[output:686eb1e9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0333 | 6-diphospho-1D-myo-inositol pentakisphosphate\n","truncated":false}}
%---
%[output:30d9e91c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0306 | CTP synthase (glutamine)\n","truncated":false}}
%---
%[output:6c957c80]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0799 | H+\n","truncated":false}}
%---
%[output:06ca38aa]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:5ef4139d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:3a9f448a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:72a4d798]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:614967e8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:342fd214]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0376 | acetyl-CoA\n","truncated":false}}
%---
%[output:4a09d5d4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:88158109]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:03c29dd4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:96e33254]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0073 | 4PP-IP5 depyrophosphorylation to IP6\n","truncated":false}}
%---
%[output:0bf3ea6f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:87387e17]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:2805c352]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:30df20bd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:3ce5e57e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0306 | CTP synthase (glutamine)\n","truncated":false}}
%---
%[output:9c124ea5]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:92582d81]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:18be25e1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:1f74a554]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:21f73b03]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:9ab76b2f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:44107a08]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:2f6cd200]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:487fa04f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:13731fda]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5234572c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:8f4c371c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0075 | 5'-methylthioadenosine phosphorylase\n","truncated":false}}
%---
%[output:265ff359]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:4a07aa46]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:42df7012]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:47a660b4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:5f8ebce4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:624dd1ae]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:13281d57]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:63ef1eb4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:91083374]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:1332a1fe]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:89511d57]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:30e0b0e7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:5359be80]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:786206c1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YGR061C\n","truncated":false}}
%---
%[output:7a29d57b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6426d08d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0075 | 5'-methylthioadenosine phosphorylase\n","truncated":false}}
%---
%[output:72c16e0e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:158a4976]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4652 | aldehyde dehydrogenase (1-propanol, NAD)\n","truncated":false}}
%---
%[output:666e8ed9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:2ffd719a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:85922ebf]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:18c7b257]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0445 | bicarbonate\n","truncated":false}}
%---
%[output:46dc9760]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0445 | bicarbonate\n","truncated":false}}
%---
%[output:0b0d7033]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0445 | bicarbonate\n","truncated":false}}
%---
%[output:471dad89]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0445 | bicarbonate\n","truncated":false}}
%---
%[output:0c03371b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:59b71001]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:79057ff9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0445 | bicarbonate\n","truncated":false}}
%---
%[output:01e898d0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0376 | acetyl-CoA\n","truncated":false}}
%---
%[output:02367a86]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:9090cee2]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0362 | acetate\n","truncated":false}}
%---
%[output:0d70ed5c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4661 | aldehyde dehydrogenase (methionol, NAD)\n","truncated":false}}
%---
%[output:160539c4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0362 | acetate\n","truncated":false}}
%---
%[output:85213e0f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:67fabede]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0075 | 5'-methylthioadenosine phosphorylase\n","truncated":false}}
%---
%[output:1367f115]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:187539b8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:96b732a2]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0993 | L-glutamate\n","truncated":false}}
%---
%[output:64079c6a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:20c8a213]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:1fde282c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL104C\n","truncated":false}}
%---
%[output:4a408d6a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0532 | coenzyme A\n","truncated":false}}
%---
%[output:0eef1021]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0532 | coenzyme A\n","truncated":false}}
%---
%[output:75149875]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0383 | adenine\n","truncated":false}}
%---
%[output:529a38b8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4670 | aldehyde dehydrogenase (tyrosol, NAD)\n","truncated":false}}
%---
%[output:3298229e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:5d0e4e1c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:2423cf4d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:54a2a669]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:0c3fc659]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:5c8e7f44]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:5307b551]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDR017C\n","truncated":false}}
%---
%[output:3a10161f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0075 | 5'-methylthioadenosine phosphorylase\n","truncated":false}}
%---
%[output:1c880e2f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:3aeb453e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:38d9c5fe]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:5f821b0c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:9ce498ff]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:233acdd1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0774 | NAPRtase\n","truncated":false}}
%---
%[output:28030707]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:69ed582e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:902d4510]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:7c2b3b8b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:999ce18e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:26afedb1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:18977cf3]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:03a2008f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:9ac35571]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:4045b817]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:6bfcc808]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:83419f72]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:69a3341c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0786 | nicotinate-nucleotide diphosphorylase (carboxylating)\n","truncated":false}}
%---
%[output:6d020ae1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:19c9b091]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL131W or YDL182W\n","truncated":false}}
%---
%[output:13afe2b8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:3a0b637b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR232W\n","truncated":false}}
%---
%[output:1ddd6416]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:2f122382]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:0b37a037]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:7e116d47]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:2259fed6]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:68429bd9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:3dfda2bb]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:46f0905b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:9bd110e0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNL104C or YOR108W\n","truncated":false}}
%---
%[output:372fac1a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:0ae5e06a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:7d2a73b7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:2b5eff96]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR232W\n","truncated":false}}
%---
%[output:279524c9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:26f6a218]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:56fc8e54]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:9acbd008]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:0f696751]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:1d1bafda]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:0750c667]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:6e5730e1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:3622a6d0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:9e7afa63]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:6fdb5c4f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:13877bb4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR232W\n","truncated":false}}
%---
%[output:185722eb]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:3a62318f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:3c09b8cc]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0386 | adenosine\n","truncated":false}}
%---
%[output:9e242f9a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0386 | adenosine\n","truncated":false}}
%---
%[output:3dcf9971]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0386 | adenosine\n","truncated":false}}
%---
%[output:9c9046d9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:3a791cb1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:5d641a51]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0386 | adenosine\n","truncated":false}}
%---
%[output:3db83209]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:948be976]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:41eaf15c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:66104d9f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR163W\n","truncated":false}}
%---
%[output:5fd70fe7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:8ce975b9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:06250398]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:57e1a99e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:4c352c78]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:893f2cc1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:68f34b90]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:7acfccb3]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:96ab9118]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:367ed87d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:318d2f2b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:01e8a824]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR163W\n","truncated":false}}
%---
%[output:8423a821]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:128a63a2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YGL080W and YGR243W) or (YGL080W and YHR162W)\n","truncated":false}}
%---
%[output:72629a22]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:8183cda4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:543c078f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:5aef3b5b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:0395e44e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:8168cb53]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:73f31171]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:3762adac]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:743df2e0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:6eeff0cb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:1d2df148]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR163W\n","truncated":false}}
%---
%[output:60f78edd]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:9fa7a0e1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:6e19a930]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:13db9f98]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:054f3d5f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:92591168]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:51ed4cff]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:007769bf]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:2248a662]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YDR148C and YFL018C and YIL125W and YFR049W) or (YDR148C and YFL018C and YIL125W)\n","truncated":false}}
%---
%[output:1d882eb7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:3acca036]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:6371e5ac]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR163W\n","truncated":false}}
%---
%[output:63c22942]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:929b206d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:5b4e8c53]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:506b98d4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:945f99fc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:95e4090a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:4bb2b288]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:5618dd7f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:41658d3d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:2ac971f8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:9184e244]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:39c87fc0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR163W\n","truncated":false}}
%---
%[output:13286603]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:77934d24]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:94ba7200]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:0bdd6800]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8c972825]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:22158923]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:8bbdcdd6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0486 | glyceraldehyde-3-phosphate dehydrogenase\n","truncated":false}}
%---
%[output:39d02512]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:093e61b1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:33b402e4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:83db84a0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:0b6e135a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:06d7751a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR163W\n","truncated":false}}
%---
%[output:9eb7508c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:6d08639e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL131W or YDL182W\n","truncated":false}}
%---
%[output:1fbfa608]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:1b385bd9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:49796456]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:27a5fe7e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:913c8fcf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0486 | glyceraldehyde-3-phosphate dehydrogenase\n","truncated":false}}
%---
%[output:2f20d0c0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:71013814]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:52713ad8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:388d1c7a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:7b382c9f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YCL009C and YMR108W) or YMR108W\n","truncated":false}}
%---
%[output:4c1343da]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR163W\n","truncated":false}}
%---
%[output:2ddcca1d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:334d9f73]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:919949b5]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:32dcb09a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:3b235e53]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:4d96c85b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:447f3f3c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:7c67dd3b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:063a5f74]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YCL009C and YMR108W) or YMR108W\n","truncated":false}}
%---
%[output:8dd05024]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:2f09f4cb]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:6839a340]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR163W\n","truncated":false}}
%---
%[output:2e945266]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:35b312f7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:4d7c6bef]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:776d4ce8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:779b7fcc]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:227653a6]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0234 | 3-methylbutanal\n","truncated":false}}
%---
%[output:12ad0643]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR399W\n","truncated":false}}
%---
%[output:4ffcba06]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:820b6611]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:495b52df]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0939 | prephenate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:7bc4937a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:6aa8bbda]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YIL107C or YOL136C\n","truncated":false}}
%---
%[output:83344567]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:09ce14ab]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:75036819]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:80181e6c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:5d5d2678]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:0bcdbd07]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:590c26a3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:2c085798]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0785 | GTP\n","truncated":false}}
%---
%[output:640a330d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:99644a19]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0973 | L-aspartate\n","truncated":false}}
%---
%[output:81100de6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:2a0a6a8c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YLR410W or YOR163W\n","truncated":false}}
%---
%[output:9c75d8dc]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:11496f55]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0973 | L-aspartate\n","truncated":false}}
%---
%[output:919cfe23]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPR035W\n","truncated":false}}
%---
%[output:451d46dd]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0234 | 3-methylbutanal\n","truncated":false}}
%---
%[output:066721ff]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0973 | L-aspartate\n","truncated":false}}
%---
%[output:9f08c9b9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:01bd6126]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0298 | 5'-adenylyl sulfate\n","truncated":false}}
%---
%[output:9556506a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0915 | phosphoribosylpyrophosphate amidotransferase\n","truncated":false}}
%---
%[output:3526ea44]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0298 | 5'-adenylyl sulfate\n","truncated":false}}
%---
%[output:42e79b6e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0298 | 5'-adenylyl sulfate\n","truncated":false}}
%---
%[output:2edb9eaa]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0298 | 5'-adenylyl sulfate\n","truncated":false}}
%---
%[output:8fb61c6b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YLR355C\n","truncated":false}}
%---
%[output:7bc58ff7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:56813306]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:08f7fcb9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:7466ce34]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3be0b2ea]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:465eebd5]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:7d475c72]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:7716e15d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:9aa0d725]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:181b4f92]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0915 | phosphoribosylpyrophosphate amidotransferase\n","truncated":false}}
%---
%[output:4dfdb819]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:8dea331a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0234 | 3-methylbutanal\n","truncated":false}}
%---
%[output:7b8587dd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR207C\n","truncated":false}}
%---
%[output:277cc47a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2d4d5962]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0961 | pyruvate dehydrogenase\n","truncated":false}}
%---
%[output:3b7de89a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:5c609d61]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:23048b03]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:4f9ea2dd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:2c88a0c7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:634b3a07]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPR035W\n","truncated":false}}
%---
%[output:37e8cbf0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:855a252b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:506f5e20]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:04b937c5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:3acfa02f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:3763ac8a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:2c2524e5]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3daebb5f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0961 | pyruvate dehydrogenase\n","truncated":false}}
%---
%[output:57ef575b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:202bfea7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:4b82a77b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:4df0008d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0234 | 3-methylbutanal\n","truncated":false}}
%---
%[output:09ee34a5]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:9483c541]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:7cde0059]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:93b94068]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:3fa4e212]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YJL130C or (YJR109C and YOR303W)\n","truncated":false}}
%---
%[output:219a5da9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3fadaa40]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:90295226]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:485abbf1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:8c074799]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0113 | acetyl-CoA synthetase\n","truncated":false}}
%---
%[output:441e34ba]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:111f1f92]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:10ea8f2c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:18efc3e7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:664d170c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:62c62de3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR218C or YGL062W\n","truncated":false}}
%---
%[output:07629fcb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:03b73189]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:7373ca0c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:65e450f1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:97a1b478]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0234 | 3-methylbutanal\n","truncated":false}}
%---
%[output:05c15ecc]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:74f1cf1c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:80ef7498]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:37a897a4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:87cc2884]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:8aa3582a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR218C or YGL062W\n","truncated":false}}
%---
%[output:86123863]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:52ab0fa5]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:6212879c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:57452e86]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:6dbe1eb1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:9a56a51c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:51f001f1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:9d8ae395]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:46823e7a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:9a303278]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:21deaf5b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR399W\n","truncated":false}}
%---
%[output:8a232de9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR218C or YGL062W\n","truncated":false}}
%---
%[output:7e3b8a9a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:2a6c8654]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2959d481]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:69542f99]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0234 | 3-methylbutanal\n","truncated":false}}
%---
%[output:6b4ea7f0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:57406199]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:73e509e8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:8e83a2d9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:1d3a1af5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR399W\n","truncated":false}}
%---
%[output:7fe90879]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:0a551857]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR218C or YGL062W\n","truncated":false}}
%---
%[output:3867f29e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:1ce62b41]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:5adade27]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:664ee7e4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:50de5680]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:6507f527]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:0671b825]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:411558f4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:0d59532b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:6101ff12]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:0efc91d6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL131W or YDL182W\n","truncated":false}}
%---
%[output:8f49e488]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:751d3b39]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:38e5acae]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:051fec21]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:86a75a4f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:11da2b1f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:3582f492]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPR035W\n","truncated":false}}
%---
%[output:86ab6c9d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:7b71028e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:83098b0b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:1a33c583]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:4ab0bb6e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3ea3026a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:25039871]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:0329d734]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:5164cc0f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:40b58781]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:13c2fd8c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:09297d9a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:1b35c56c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:0cadc856]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0383 | adenine\n","truncated":false}}
%---
%[output:69bd9140]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:1782f8b6]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2466f9be]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4211 | D-ribose 5-phosphate,D-glyceraldehyde 3-phosphate pyridoxal 5-phosphate-lyase (glutamine-hydrolyzing)\n","truncated":false}}
%---
%[output:5614f277]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:96b2c83f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:092e6233]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:91fe8fc0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:5f2a4aca]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:70a38319]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:957a767b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:7fbede02]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:5c41e5c2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YJL071W or YMR062C\n","truncated":false}}
%---
%[output:89afb399]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:602b20f5]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2c312842]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:91566213]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPR035W\n","truncated":false}}
%---
%[output:4fad01e4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3c88cfbd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4212 | D-ribulose 5-phosphate,D-glyceraldehyde 3-phosphate pyridoxal 5-phosphate-lyase\n","truncated":false}}
%---
%[output:6c4fe67b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:8f5016b8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0383 | adenine\n","truncated":false}}
%---
%[output:239a18e6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:967f5b47]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:87133656]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR221C and YER178W and YFL018C and YGR193C and YNL071W\n","truncated":false}}
%---
%[output:0956a49f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:1277d0d9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:104c45a0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:26686d6d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:829fe5e0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:58f9baf8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:28fa7788]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:79d2eff3]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:167f55af]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:410a0ef6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:75da3f01]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR221C and YER178W and YFL018C and YGR193C and YNL071W\n","truncated":false}}
%---
%[output:257ecfcd]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:6e29a630]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:8e25644d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:6c344d21]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:28a657c4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:0ef05ada]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:3eca4361]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:9f1bdd92]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:9ebe2ae2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0383 | adenine\n","truncated":false}}
%---
%[output:487974f7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:1ceb472d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDR441C or YML022W\n","truncated":false}}
%---
%[output:0a51965a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:96770e16]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:5fc375b6]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2ccf9583]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0938 | prephenate dehydratase\n","truncated":false}}
%---
%[output:9401a413]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3027bfca]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:3562f533]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:1abd9d73]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2b5b22dd]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:84e85b97]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:48596f90]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:8ed36f58]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR018C\n","truncated":false}}
%---
%[output:75ab5cd0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:16727cde]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:2dcc2f06]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:78f6b265]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:8c2dac52]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:4c05dd78]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:29a4bdea]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:6bd20915]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0383 | adenine\n","truncated":false}}
%---
%[output:51b26e0b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR209C\n","truncated":false}}
%---
%[output:875511ae]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:342a7dd9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YMR250W\n","truncated":false}}
%---
%[output:9fc2812b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:26decc7c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2f33f6ba]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:42d7c83a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:0c1705f1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:89421343]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:48bb7b10]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2f660c06]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR209C\n","truncated":false}}
%---
%[output:127b2b69]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:70ffb2a1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0178 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:396c0629]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:6ec14811]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:97773ed6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:9128712b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:322e2456]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:8eedb8dd]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:9a601da9]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:100fc03b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR209C\n","truncated":false}}
%---
%[output:5117c54a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:58284af8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:1ababc29]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:92000377]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:991721a1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:95a548fb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0832 | oxoglutarate dehydrogenase (lipoamide)\n","truncated":false}}
%---
%[output:67be94a0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:6351792b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0016 | 2-aceto-2-hydroxybutanoate synthase\n","truncated":false}}
%---
%[output:80d0aef9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:48fa1c77]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR209C\n","truncated":false}}
%---
%[output:15621b79]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:4e8dca6b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:8cc07c19]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:5493ff6c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:4690ea17]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNR050C\n","truncated":false}}
%---
%[output:82682ee8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:9893cf01]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:55eb7a0e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:780ab89e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0062 | 3-methyl-2-oxobutanoate decarboxylase\n","truncated":false}}
%---
%[output:1ba0293e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:07233b92]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:64b05e85]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:345872d1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:07ed936c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0097 | acetolactate synthase\n","truncated":false}}
%---
%[output:82e0a147]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:20d8e011]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL103C\n","truncated":false}}
%---
%[output:7724996e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:482387fc]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:1b85ac05]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3d34e4f4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:1ddf1508]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:55128a5d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:50a9187a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:678c3783]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:095d27ff]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2f527a98]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:45528a47]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YHL012W or YKL035W\n","truncated":false}}
%---
%[output:6260723e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:8948bd1b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4039 | succinyl-CoA:acetate CoA transferase\n","truncated":false}}
%---
%[output:43280c3c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:584b9bdb]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:01045405]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:47e3099e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3c046561]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:8d8b385c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:427be006]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:92bc6272]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:92c9808f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR192W\n","truncated":false}}
%---
%[output:697a15f5]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:77c05889]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:8f7cf411]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:67aaa049]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:06d400aa]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:51d3bd4d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:48b8fb45]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:18e0fa50]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:4051cc8f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:667d2f2a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:0b210460]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR192W\n","truncated":false}}
%---
%[output:6421a41b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:09d9fc48]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:7dda08ff]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:37f41a9c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:10f83ef6]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:10858acd]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:7c722c12]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:36715477]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:08d3b8c0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:0826466c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:38efcd4f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPR035W\n","truncated":false}}
%---
%[output:5325157d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR209C\n","truncated":false}}
%---
%[output:4aea5c9a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:338108f2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:7e76e449]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:12267585]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:81bf6c2b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:810b4fba]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:81cab6c4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:713d68a9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:7983b2e5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:4338fa29]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:783c98b5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YFR047C\n","truncated":false}}
%---
%[output:1ea73661]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3affd9c8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDL141W and YNR016C\n","truncated":false}}
%---
%[output:8242ac53]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:17f22f49]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:69193989]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:6c88e182]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:759b9a39]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:4e368eb4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:175d72ea]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:0a9899d7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:7cde2a09]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR018C\n","truncated":false}}
%---
%[output:4eca035b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:43238c70]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:36d0b64e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBL015W\n","truncated":false}}
%---
%[output:2a48ab8d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:0b98c157]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:4d4bcf50]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:76097afc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPR035W\n","truncated":false}}
%---
%[output:22b6ab77]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:725ad433]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:84e9c6aa]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:238cceca]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YMR250W\n","truncated":false}}
%---
%[output:4f5bc31c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:64f7a24b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:3bacc4ed]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:650b397a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YAL054C or YLR153C\n","truncated":false}}
%---
%[output:3a00fee6]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:7036f211]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:0a4ef3b8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:092762a1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:6e184154]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:6d15cff7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:947abd2f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:4748c767]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:9e68d627]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:08f51f72]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:59e9a6af]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:7cbd5c20]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:5ec245e4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YAL054C or YLR153C\n","truncated":false}}
%---
%[output:26253693]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:051fbd17]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:4552d605]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:711560f6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:5d347ab8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:7ecd3b36]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:0fbba460]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:479717c8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:64a72d47]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:9f560ab1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:81d5e56f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:6f94daff]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOL140W\n","truncated":false}}
%---
%[output:69ff599c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3f38afc3]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:5c1b939f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:87712d69]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNR050C\n","truncated":false}}
%---
%[output:613597ba]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:588fabf7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0079 | 5'-phosphoribosylformyl glycinamidine synthetase\n","truncated":false}}
%---
%[output:6206e013]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBL039C or YJR103W\n","truncated":false}}
%---
%[output:111a9ea6]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:85ad1006]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:8adc94b5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:75ae5515]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:2b280fe6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL148C\n","truncated":false}}
%---
%[output:027fc281]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:7ed69ab3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:2c7ade8b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL103C\n","truncated":false}}
%---
%[output:3534668e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:0b643257]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNL220W\n","truncated":false}}
%---
%[output:3dca75c5]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:3c5a2db0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:7262880d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:58fcc7c1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:0edb288f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:0b7bb37c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:219a1f1b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL148C\n","truncated":false}}
%---
%[output:323112ca]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:7150b0c6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YHL012W or YKL035W\n","truncated":false}}
%---
%[output:4bb51e16]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNL220W\n","truncated":false}}
%---
%[output:2bba4fa7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:34c50c3d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:92eb782f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:6aa2f724]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:58fb12ed]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:3a400476]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:59173ff9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:49bccf6e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:2e393d7d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:6879a94a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNL220W\n","truncated":false}}
%---
%[output:4e880cf7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR192W\n","truncated":false}}
%---
%[output:3fa0a74d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:35cd31db]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:4fa7efa9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:6385a270]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:96642862]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0477 | glutamine-fructose-6-phosphate transaminase\n","truncated":false}}
%---
%[output:57d863e4]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:858a4990]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:5fc9d260]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:493fc05f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNL220W\n","truncated":false}}
%---
%[output:7b48b4da]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:90fe4a6f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR192W\n","truncated":false}}
%---
%[output:1675ca50]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:34f94b64]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:7fff1f9b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:49ffb4d8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:0cb36107]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:8d7c480d]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:47104afc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:8f1f4db0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNL220W\n","truncated":false}}
%---
%[output:63e78da1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0093 | 6PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:012da4b0]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:7413ba04]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:6dd12d38]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:7cd80102]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:9738338c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:643166df]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:0c0d5a08]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:3e18db6f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:2cdc4950]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER070W or YGR180C or YIL066C or YJL026W\n","truncated":false}}
%---
%[output:627d0339]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1386 | PRPP\n","truncated":false}}
%---
%[output:4dc9bf56]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:0d6f8b35]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:5a684d5c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR410W\n","truncated":false}}
%---
%[output:63f835ca]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:42ec25fa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:9fcb5e5e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:94b6bc00]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:7d9d983b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:9fe93c72]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER070W or YGR180C or YIL066C or YJL026W\n","truncated":false}}
%---
%[output:496f46fd]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:5b14822d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:99af7427]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:2dfbb5c1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:652a09c8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDR017C\n","truncated":false}}
%---
%[output:72d43335]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:6b36550b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:824d7f67]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1386 | PRPP\n","truncated":false}}
%---
%[output:3f1144a6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:9557dfb9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YCL050C\n","truncated":false}}
%---
%[output:9d1f8744]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:177d0406]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1838 | homocitrate synthase\n","truncated":false}}
%---
%[output:68ad3163]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:2d60fbaa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0081 | 5-aminolevulinate synthase\n","truncated":false}}
%---
%[output:3639bc89]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:9a4f5e14]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:92b99444]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:63fa5cb7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:3cd231a8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:193a01cf]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER070W or YGR180C or YIL066C or YJL026W\n","truncated":false}}
%---
%[output:726cc0a7]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:5043eae5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:2985f875]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:16932fa5]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:4b8c1648]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:0708eca4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0024 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:24be4082]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR209C\n","truncated":false}}
%---
%[output:6961588c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:35fd4a75]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:4406585f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YCL050C\n","truncated":false}}
%---
%[output:67f4ba6c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:4a396962]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0081 | 5-aminolevulinate synthase\n","truncated":false}}
%---
%[output:2d43b53a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0361 | acetaldehyde\n","truncated":false}}
%---
%[output:3ad182c5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:02a73106]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:97f4d66c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:73643e93]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:0e5dcdd4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR209C\n","truncated":false}}
%---
%[output:11a09786]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:202723f8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER070W or YGR180C or YIL066C or YJL026W\n","truncated":false}}
%---
%[output:7a38413b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:81699372]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:936d0b6e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:44459907]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:316afea3]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:3898a663]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:3f689b34]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:826a8bd0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0081 | 5-aminolevulinate synthase\n","truncated":false}}
%---
%[output:5932ddd0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR209C\n","truncated":false}}
%---
%[output:7f67e915]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER070W or YGR180C or YIL066C or YJL026W\n","truncated":false}}
%---
%[output:717da4ae]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:71ea91f8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:292bfb06]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:0cc8463a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:51bfc35f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:01c3a757]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0959 | pyruvate decarboxylase\n","truncated":false}}
%---
%[output:607b4d8e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:006e6e12]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:1383aee1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:76da6b24]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YCL050C\n","truncated":false}}
%---
%[output:35477083]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR209C\n","truncated":false}}
%---
%[output:404f7317]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:658d793f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:562cef7e]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:1f300e9e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0082 | 5-diphosphoinositol-1,2,3,4,6-pentakisphosphate diphosphohydrolase\n","truncated":false}}
%---
%[output:85fe7a5f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8e4cc1ed]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:88bd5d89]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:3d21a863]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:84c17723]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER070W or YGR180C or YIL066C or YJL026W\n","truncated":false}}
%---
%[output:1b385ee4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:2f12dbe1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR018C\n","truncated":false}}
%---
%[output:15371d92]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:0b93f9a6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:3758b38f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:84cf394f]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:3f12bcb3]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:003b351b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:8ada2ad1]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:7f6a4a09]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER070W or YGR180C or YIL066C or YJL026W\n","truncated":false}}
%---
%[output:406573cf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0082 | 5-diphosphoinositol-1,2,3,4,6-pentakisphosphate diphosphohydrolase\n","truncated":false}}
%---
%[output:60fec791]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:57fd5763]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YMR250W\n","truncated":false}}
%---
%[output:068dae8b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:97ad8b79]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2034 | pyruvate transport\n","truncated":false}}
%---
%[output:4899a452]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1fe0333b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:95b85f43]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:7b3f5128]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:35e6239f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR245C\n","truncated":false}}
%---
%[output:07d121f2]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:4a7e00bb]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:2417f3ea]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:909cd4d6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL171C\n","truncated":false}}
%---
%[output:8afbecf8]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:2917d814]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:72e449aa]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:87c886d6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:51d443e6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0082 | 5-diphosphoinositol-1,2,3,4,6-pentakisphosphate diphosphohydrolase\n","truncated":false}}
%---
%[output:86b6aab9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER070W or YGR180C or YIL066C or YJL026W\n","truncated":false}}
%---
%[output:0e1bb0f9]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:63ed2dfa]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:43c22e8a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:2b16f06a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0fa7473a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:7b4f3964]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:0b524a00]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:06f338f3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:0621487a]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:9daff6a9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER070W or YGR180C or YIL066C or YJL026W\n","truncated":false}}
%---
%[output:9f37109b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:25cc4849]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0832 | oxoglutarate dehydrogenase (lipoamide)\n","truncated":false}}
%---
%[output:4e770d23]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:8c69ac28]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:20bb0d78]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0082 | 5-diphosphoinositol-1,2,3,4,6-pentakisphosphate diphosphohydrolase\n","truncated":false}}
%---
%[output:70ae1f0f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNR050C\n","truncated":false}}
%---
%[output:1e4b9ec3]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:6af7a312]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:12e01a0c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_0794 | H+\n","truncated":false}}
%---
%[output:96cc7f0a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER070W or YGR180C or YIL066C or YJL026W\n","truncated":false}}
%---
%[output:83fffd2f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:8eda4cc1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:5a65f208]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5c1d11e1]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:9c9e1e0a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:5f3ba98c]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:975909f3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL103C\n","truncated":false}}
%---
%[output:9123b4ec]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:9706c1c1]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:82589d21]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:7c5cf69b]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:3632c1f9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0082 | 5-diphosphoinositol-1,2,3,4,6-pentakisphosphate diphosphohydrolase\n","truncated":false}}
%---
%[output:93f87385]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:8afdbdd9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:1baea7ba]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:4c0d1604]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:2f2298ee]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:7ebef416]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YHL012W or YKL035W\n","truncated":false}}
%---
%[output:2bc64d27]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:9d910b19]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:88375cbf]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:4ec43ee0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0764 | glyceraldehyde 3-phosphate\n","truncated":false}}
%---
%[output:02275730]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:20e037a8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:35d444c4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:00489b40]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:648ef5fd]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:2c33d748]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:67a86635]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR192W\n","truncated":false}}
%---
%[output:82db50ad]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:71d9f7c6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0082 | 5-diphosphoinositol-1,2,3,4,6-pentakisphosphate diphosphohydrolase\n","truncated":false}}
%---
%[output:143ae5fb]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:92f9db14]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1838 | homocitrate synthase\n","truncated":false}}
%---
%[output:7bf730ad]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:58461595]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:2255eded]
%   data: {"dataType":"text","outputData":{"text":"  m_d  = s_1318 | phenylacetaldehyde\n","truncated":false}}
%---
%[output:68b8dd0f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:44860991]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:40a9b9ef]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0764 | glyceraldehyde 3-phosphate\n","truncated":false}}
%---
%[output:31e0bd12]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:13debe81]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR192W\n","truncated":false}}
%---
%[output:6b8b6016]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:20a41a27]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:259d3a7e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:1cd31c42]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:65041980]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0016 | 2-aceto-2-hydroxybutanoate synthase\n","truncated":false}}
%---
%[output:503796a3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0082 | 5-diphosphoinositol-1,2,3,4,6-pentakisphosphate diphosphohydrolase\n","truncated":false}}
%---
%[output:07e4ed86]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:25543ae0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:1641d38f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:13921e26]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:46800862]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:6cbe00d2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:40fe822f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:144d024c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:70effe05]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:525ce569]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:10e1cd2a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:2590c6ac]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0097 | acetolactate synthase\n","truncated":false}}
%---
%[output:6ecd2396]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:7d58aaf3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:74bf4e98]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:9cbdea23]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:11df0a6c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:503fd89f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0082 | 5-diphosphoinositol-1,2,3,4,6-pentakisphosphate diphosphohydrolase\n","truncated":false}}
%---
%[output:2d83e7ff]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:58ac33f1]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:0762c4da]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:26cb82e0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:1c566f66]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:13cdfb40]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:390950ef]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0532 | coenzyme A\n","truncated":false}}
%---
%[output:9a207f34]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0562 | hypoxanthine phosphoribosyltransferase (Hypoxanthine)\n","truncated":false}}
%---
%[output:92b66632]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL055C\n","truncated":false}}
%---
%[output:14aac98e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:995cd326]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:2912060a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0532 | coenzyme A\n","truncated":false}}
%---
%[output:92404a45]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1377 | prephenate\n","truncated":false}}
%---
%[output:3d14f584]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0532 | coenzyme A\n","truncated":false}}
%---
%[output:935e243b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:9a4e29ef]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0090 | 6-phosphofructo-2-kinase\n","truncated":false}}
%---
%[output:405c604f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:4032b353]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:4fa7a25b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:4875149e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:2b855b80]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:82467f82]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:3287329f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:53e7d15f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:66862dc7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:4bd64dd4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:5ec2e67e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:003aa0b8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:13b9e6aa]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:7f159dab]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:3cb43c10]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL198C\n","truncated":false}}
%---
%[output:7869d09c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0092 | 6PP-IP5 depyrophosphorylation to IP6\n","truncated":false}}
%---
%[output:283c6725]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:6ce2876e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:2e5cdb7c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:7e048f63]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:58bfeba3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0476 | glutamine synthetase\n","truncated":false}}
%---
%[output:18e62070]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1210 | NADP(+)\n","truncated":false}}
%---
%[output:1af896cb]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:3fefacc5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:834516a0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:4ea3bdc6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:51ef09ce]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:029f99bb]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:7e42ee3f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:0e0bd2e7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:10cee372]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:8f7b2aff]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:91e4b6c4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0096 | acetohydroxy acid isomeroreductase\n","truncated":false}}
%---
%[output:52828764]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:98bac089]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:916fc741]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:1a3d1abc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:9a216205]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:0e10553a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:83f16f97]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:4ffdbbca]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:54c29f01]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:3244234c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:293c0fe3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:5d67851c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:6a4ae9ce]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:4ce433a7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:0229da21]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL055C\n","truncated":false}}
%---
%[output:3ff4df54]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:3d8d5523]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:534dfec2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0108 | acetyl-Coa carboxylase\n","truncated":false}}
%---
%[output:47351e70]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:660005c7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1401 | pyruvate\n","truncated":false}}
%---
%[output:279f892a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:87732624]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL141W\n","truncated":false}}
%---
%[output:8da8479d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:3a3e83fc]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:500cd33a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:22dc32aa]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:41c835f7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL055C\n","truncated":false}}
%---
%[output:8a2a0dbd]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:52914fb6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0476 | glutamine synthetase\n","truncated":false}}
%---
%[output:6b901616]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:426f067b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:123a3273]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:0a0ba531]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:1d7d4954]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YJR105W\n","truncated":false}}
%---
%[output:7bde413c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:77a4aecc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:7ee089fe]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:6b29a911]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:00aa4729]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:81dc7e06]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1401 | pyruvate\n","truncated":false}}
%---
%[output:64ea656a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:67bc4cb1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:170bbd7d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:13f6f81e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0182 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:0dd60b59]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:4d7dab8b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:2983952e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YJR105W\n","truncated":false}}
%---
%[output:06d71b32]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:0e2b5f3d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:6b93a191]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0843 | hypoxanthine\n","truncated":false}}
%---
%[output:85713784]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:23c33b76]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0250 | carbamoyl-phosphate synthase (glutamine-hydrolysing)\n","truncated":false}}
%---
%[output:16c03c5e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:0939d0d1]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:594ab597]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:504a78ca]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:7fa9d06e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:22a9111d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0365 | acetate\n","truncated":false}}
%---
%[output:7d03a921]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YJR105W\n","truncated":false}}
%---
%[output:6ab0424c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:57c407d5]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:22096909]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:832ad6a2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:6e7066dc]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:9f1e73fd]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:3bcc5608]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:57bf3c32]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0958 | pyruvate carboxylase\n","truncated":false}}
%---
%[output:88f60adf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:3bf4a424]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:5e154e2e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:7afded06]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL198C\n","truncated":false}}
%---
%[output:945bee65]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YJR105W\n","truncated":false}}
%---
%[output:9dc5e0c8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0843 | hypoxanthine\n","truncated":false}}
%---
%[output:5b4dddb9]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0843 | hypoxanthine\n","truncated":false}}
%---
%[output:7e9bde79]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:073fd259]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:618f745a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:22f0e9e2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:89fcb61f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:78264476]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:830b1451]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0958 | pyruvate carboxylase\n","truncated":false}}
%---
%[output:7a13b6bd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL198C\n","truncated":false}}
%---
%[output:94f5b290]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:2d7caaf3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:4e475476]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:9d012e51]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:6d9f9735]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:3426c9e3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:73878c87]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:1302d4ec]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:6ee7f019]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:10a91027]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0856 | inosine\n","truncated":false}}
%---
%[output:1fa7b38f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:5d04807d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:1131e20c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0856 | inosine\n","truncated":false}}
%---
%[output:6d1b0638]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:732be933]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:0761ec6e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0562 | hypoxanthine phosphoribosyltransferase (Hypoxanthine)\n","truncated":false}}
%---
%[output:2f6bb058]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0958 | pyruvate carboxylase\n","truncated":false}}
%---
%[output:36cc0546]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0856 | inosine\n","truncated":false}}
%---
%[output:3b83bcf6]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0856 | inosine\n","truncated":false}}
%---
%[output:84a0fa97]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:8bf0fa99]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:1c2a2fd6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNR012W\n","truncated":false}}
%---
%[output:565102f9]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:4f564a5f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL062W or YOR375C\n","truncated":false}}
%---
%[output:0231ebcc]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:1cb0c2ca]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:314995a7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:6dff560a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:50fc8e4d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0562 | hypoxanthine phosphoribosyltransferase (Hypoxanthine)\n","truncated":false}}
%---
%[output:3753a8c6]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:8fa81526]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0958 | pyruvate carboxylase\n","truncated":false}}
%---
%[output:4286e7b8]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:7357c516]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:8e81b7df]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:54335c98]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:7f499ecb]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:79ad2897]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:193b09b0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:67d0bb4d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:86281533]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:0145ed57]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:6dcdc2da]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:74794c1f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:71dbf121]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:53295cd7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:526cc4a4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0419 | ammonium\n","truncated":false}}
%---
%[output:5f8ac1e9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1838 | homocitrate synthase\n","truncated":false}}
%---
%[output:197d80cf]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0849 | IMP\n","truncated":false}}
%---
%[output:3384609c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0849 | IMP\n","truncated":false}}
%---
%[output:50af68d5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:899b1066]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:8af491aa]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0849 | IMP\n","truncated":false}}
%---
%[output:16276320]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:37221105]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR410W\n","truncated":false}}
%---
%[output:5db9c53b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0849 | IMP\n","truncated":false}}
%---
%[output:106ef868]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:4ada5c0c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0476 | glutamine synthetase\n","truncated":false}}
%---
%[output:49aa337c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0849 | IMP\n","truncated":false}}
%---
%[output:871e95fd]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0739 | GDP\n","truncated":false}}
%---
%[output:21a5997a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0739 | GDP\n","truncated":false}}
%---
%[output:8de37d26]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4717 | benzyl-acetate esterase, c\n","truncated":false}}
%---
%[output:8504ebb1]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0739 | GDP\n","truncated":false}}
%---
%[output:5439e6df]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:66b56ad1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR017C\n","truncated":false}}
%---
%[output:5f0c043b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0739 | GDP\n","truncated":false}}
%---
%[output:2c46c402]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0739 | GDP\n","truncated":false}}
%---
%[output:3273910f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:4f3a9654]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0739 | GDP\n","truncated":false}}
%---
%[output:4878316f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:63b93866]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:7768d038]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:32875c2b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0739 | GDP\n","truncated":false}}
%---
%[output:8a9f24f7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0739 | GDP\n","truncated":false}}
%---
%[output:47a0700d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YER053C or YJR077C\n","truncated":false}}
%---
%[output:398e665a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0739 | GDP\n","truncated":false}}
%---
%[output:348be378]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:0c8e7bc7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4717 | benzyl-acetate esterase, c\n","truncated":false}}
%---
%[output:5686a756]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:7f0cf9a9]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:8a92ad6c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:7fba6f7f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:2bbd7f2b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:24219de3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:5192ed95]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:0ab61f85]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:105cdd2f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:21a9c03d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:1a3a84fc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:484d12ef]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:06fd1d37]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:7e03e276]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:32ca5770]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0761 | N-acteylglutamate synthase\n","truncated":false}}
%---
%[output:82c86b80]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:87f83e45]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:5908a2de]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:6451264a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:5924b182]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0476 | glutamine synthetase\n","truncated":false}}
%---
%[output:8da1554e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:942be068]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:4240a343]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:9ec97ca4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:3750789b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:0226794e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:3b2cd744]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:322344aa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:40cbc811]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:2099f60f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0961 | pyruvate dehydrogenase\n","truncated":false}}
%---
%[output:279b2da5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:55145175]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:8a44b282]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:35f42ce5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:7114c91c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:081686c0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:02e4273f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:5ad03dc2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:7a04deb5]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:7561bb21]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:38df06b3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:04c6cc1e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:77943dd0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:2758241f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0961 | pyruvate dehydrogenase\n","truncated":false}}
%---
%[output:89a441ed]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:59cdde8c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:0de23352]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:7870893e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:6326e146]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:8db969f7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:05e37167]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:03aa8b1d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1322 | phosphate\n","truncated":false}}
%---
%[output:604fa63f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:44d9363c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:17cf9d66]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:4d0f0609]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:689cd198]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:2197c516]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:1f6cddd2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:1fe6cd0e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0139 | adenine phosphoribosyltransferase\n","truncated":false}}
%---
%[output:02ce0112]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:4439fee4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:281739ac]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:43677fa4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:47360468]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:8ab440b3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:8c1b4e08]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:2a00d7e5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:24a7360c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:0e032535]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:360f628f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:3196e42f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:5ac05008]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:5a409d35]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:4625e3ec]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:83c33ece]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0459 | galactose-1-phosphate uridylyltransferase\n","truncated":false}}
%---
%[output:46679fd4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:766bdfa0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:02954e8d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL055C\n","truncated":false}}
%---
%[output:917d2c2b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:32ad3be7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:67aaec13]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:21043f07]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:269dd0ad]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:2103ba0e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YML035C\n","truncated":false}}
%---
%[output:6ee28fe1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5aa22d1d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:0a755934]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:9d3043a5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0951 | inosine phosphorylase\n","truncated":false}}
%---
%[output:61cf1ff0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL055C\n","truncated":false}}
%---
%[output:885f4191]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:20500f90]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:01516819]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:66424d7a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:565ef731]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:9a5651e0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:0829f021]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:61f8061f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:2c16f4d2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:43574a1e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:6d70141e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL055C\n","truncated":false}}
%---
%[output:2d170537]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:3971af63]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0951 | inosine phosphorylase\n","truncated":false}}
%---
%[output:269b68ce]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:213570d2]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0852 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:13d193c0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:15382885]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:67e9a9bf]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:58427691]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:978e9da5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:96e763d8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:2e6b5dfe]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL055C\n","truncated":false}}
%---
%[output:4168022c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:364eac83]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:12c4fd54]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:25f8c1e3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0951 | inosine phosphorylase\n","truncated":false}}
%---
%[output:56c0d2e1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:3c706868]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:3f219f2b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:74a9f862]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:036d74a3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0973 | ribonucleoside-triphosphate reductase (UTP)\n","truncated":false}}
%---
%[output:48a21072]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:5fb6aec2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:3ad4bc92]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:4f5898b2]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0182 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:268a96d5]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:03bbf093]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:7f2b5e1d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:1396148c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:6d760c0c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0951 | inosine phosphorylase\n","truncated":false}}
%---
%[output:2b361cb5]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:6fd19446]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:65e997b4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0075a7cb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:091c9f97]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:631614df]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:378da015]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:468ed43f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:2c1bb521]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:85402f60]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:484152f1]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:2241d781]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0232 | 3-methyl-2-oxobutanoate\n","truncated":false}}
%---
%[output:8ec6fddc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:9393079a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:3028b632]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:12d05fe8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:9958a29e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:9427462b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:10be65f5]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0ae5ddda]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:9902be2d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:99ccb90f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1069 | UDP-N-acetylglucosamine diphosphorylase\n","truncated":false}}
%---
%[output:4caf4487]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:2f890ad2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:8d531fa2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:8f11cef0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR209C\n","truncated":false}}
%---
%[output:0ceb4ea8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:67860d09]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:2c35f629]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:26d939c8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:35c567ed]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:5e49815a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:7c96d9b1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:431e29e7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:8d94946e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:44b70560]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1084 | UTP-glucose-1-phosphate uridylyltransferase\n","truncated":false}}
%---
%[output:86ac9c32]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:99653d7d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL198C\n","truncated":false}}
%---
%[output:48c4f104]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:5c365f69]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:52b50af5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:75a1ace0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:71ea00ee]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0fb5e2b8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:5f4cdce3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:2d3ad30f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:9a0bd714]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:1c71a002]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:8b5434d1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL198C\n","truncated":false}}
%---
%[output:481fa8bb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:1082281b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:28df34b2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1275 | UTP transport\n","truncated":false}}
%---
%[output:9aa73e54]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:93734db7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:191331fc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:9b388da7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:25df37af]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:490e9777]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:35819416]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0b1f33fe]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL198C\n","truncated":false}}
%---
%[output:298bd542]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:1b7e458d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:2e1e36a5]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:1b2f2692]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:8d92b5fd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:495719ec]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:9390cda0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR533C or YMR322C or YOR391C or YPL280W\n","truncated":false}}
%---
%[output:1a981335]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1276 | UTP\/UMP antiport\n","truncated":false}}
%---
%[output:97596f5f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:7adccd10]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:36a475f8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YDL198C\n","truncated":false}}
%---
%[output:2cae2caf]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:3ac1a829]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:657485a2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:5c24c4f2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:40994209]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:8728f819]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:958078bb]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:077eb4b8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:2c8a80f3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0d76f144]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:61f1fe27]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR124W or YPR145W\n","truncated":false}}
%---
%[output:384bba9e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0476 | glutamine synthetase\n","truncated":false}}
%---
%[output:7a4d1176]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0774 | NAPRtase\n","truncated":false}}
%---
%[output:04fe9ba5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:662e877f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:1c043c9c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:75536170]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:468c7c8f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:2245e65d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:74a2f02e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:71136f15]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0ac1f6e3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR027C\n","truncated":false}}
%---
%[output:66015d65]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:88a72a8d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:1a78529d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:04067779]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YGL080W and YGR243W) or (YGL080W and YHR162W)\n","truncated":false}}
%---
%[output:263bf04a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0b630542]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0786 | nicotinate-nucleotide diphosphorylase (carboxylating)\n","truncated":false}}
%---
%[output:14578637]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:7069e94d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0109 | acetyl-CoA carboxylase, reaction\n","truncated":false}}
%---
%[output:7b3e3450]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:24b8465d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:9b82c576]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YLR027C\n","truncated":false}}
%---
%[output:0315d6bd]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:68dd3624]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:1f3c623e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7fba1d13]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:4644af50]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:5df81b06]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:868713a8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:49f9b334]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:84107c08]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0459 | galactose-1-phosphate uridylyltransferase\n","truncated":false}}
%---
%[output:5e5e9f5d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:433e3e64]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:461fde69]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:1fb871db]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:1211c64f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:249e8863]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0111 | acetyl-CoA hydrolase\n","truncated":false}}
%---
%[output:3abbbd5a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:705f0ba7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:37f283e6]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:258cbcc2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0476 | glutamine synthetase\n","truncated":false}}
%---
%[output:841b7545]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:8313ad35]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:20d96a9f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:44046a68]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:5b94d4fb]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:69710649]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:3428de36]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:3a2f7533]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:986e44d4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:3dd90dc4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:009a76c1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0112 | acetyl-CoA synthetase\n","truncated":false}}
%---
%[output:5773a26c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:3f5185c7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:1ac8d568]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:9c77f2de]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:384b0936]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:3aef273f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:7569c32e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:57d1a236]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:9ffe48b9]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:286cbbb7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:34c6b01d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:7d9dcebf]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:4bb8f540]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:1ff24d29]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:78d1080b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCL050C\n","truncated":false}}
%---
%[output:39104e4e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:8703710c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:1c49c597]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0112 | acetyl-CoA synthetase\n","truncated":false}}
%---
%[output:60907120]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:39f37174]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:19df4136]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:7b4acf07]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0917f6fe]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:5d63cb7f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:53f84a0f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0973 | ribonucleoside-triphosphate reductase (UTP)\n","truncated":false}}
%---
%[output:6c309ab1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:66641eda]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:4f62a34e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:535f0f06]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0680 | ethanol\n","truncated":false}}
%---
%[output:606c70d8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:0c394c2e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:62f9b160]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:4d4fd206]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:7e399add]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:88cae9d8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0118 | acteylornithine transaminase\n","truncated":false}}
%---
%[output:9cd5c4c7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0fc143ab]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:09c34be3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:590f0e08]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:5de10db1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:06a0dbc7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:868a019c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:973a7e36]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:95883223]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:88cf8abd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0307 | CTP synthase (NH3)\n","truncated":false}}
%---
%[output:288dc84d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:651f8d84]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:50ec109e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:84607e40]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:9f188243]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0aa3c219]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0119 | acyl carrier protein synthase\n","truncated":false}}
%---
%[output:2dda2c32]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:27bf3395]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:7c827be7]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:136c1571]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1069 | UDP-N-acetylglucosamine diphosphorylase\n","truncated":false}}
%---
%[output:3acb8156]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:84c6948e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:40e6ee65]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:376917ff]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:6ce6199b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:9ab496f0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:106c1849]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:83762fb2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:0c0b9532]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:8eb6c0c3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:90f07bef]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:53b312f9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0119 | acyl carrier protein synthase\n","truncated":false}}
%---
%[output:66693ded]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:978fbe13]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:8766bdce]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:65f83dfc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1084 | UTP-glucose-1-phosphate uridylyltransferase\n","truncated":false}}
%---
%[output:4bdc7d8c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:38a1e089]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:44e63ea3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:47329d18]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:20983280]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:598a18f8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:5e198369]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:096e4178]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:91f557e0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:191b67dc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:761859b0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:91888963]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:5879d413]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:7358f93f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:09c579c4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:87ffb88d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1275 | UTP transport\n","truncated":false}}
%---
%[output:5c4db51a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0684b5f3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:3f396eb9]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0411fba4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:93027442]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:4d899d69]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:3838aee6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:6284c3e0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:2d6f8393]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:043f699a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:04f70217]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:0e68209e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:6abb8865]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:33d08b1e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1276 | UTP\/UMP antiport\n","truncated":false}}
%---
%[output:26b621c6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:6006a7d1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:06b3bc71]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:0614f41d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:2997f0a2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:472a0eb1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:248f905b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:5a16052a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:8458e8b1]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:266764be]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0169 | 2-methylbutanol\n","truncated":false}}
%---
%[output:1454c3df]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:6221f475]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:84d95b2d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0169 | 2-methylbutanol\n","truncated":false}}
%---
%[output:5d8afce2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:15e3319a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0169 | 2-methylbutanol\n","truncated":false}}
%---
%[output:971e0ee6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:40c2bc08]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:0181af2c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0169 | 2-methylbutanol\n","truncated":false}}
%---
%[output:7d72769f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:64955de5]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:468c3142]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:22ca3c6c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:6a92e7ae]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8280a161]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:490fc131]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:1c0eede7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0978 | ribonucleotide reductase\n","truncated":false}}
%---
%[output:5e5d6bff]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:5b8110b5]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:436f916e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:0b65bad8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0185 | 2-phenylethanol\n","truncated":false}}
%---
%[output:6904afd8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:4d7f7cdc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:5a26d20b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0185 | 2-phenylethanol\n","truncated":false}}
%---
%[output:49cccf24]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0185 | 2-phenylethanol\n","truncated":false}}
%---
%[output:47877371]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:827a4710]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:9724f269]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0185 | 2-phenylethanol\n","truncated":false}}
%---
%[output:2c762a4c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:11b95e0e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0978 | ribonucleotide reductase\n","truncated":false}}
%---
%[output:72f8a101]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:3e938fbd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:321ffbdf]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:54b2fe40]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:9041b75e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:7698677a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:4ea6d06a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:896c3c03]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:3e6aa175]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:5956a609]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0362 | acetate\n","truncated":false}}
%---
%[output:2f1ad408]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0362 | acetate\n","truncated":false}}
%---
%[output:3ebbbf2f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:2a4c6ddd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:6a7c2897]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0224 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:00e545f8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:904174e3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:61801f2b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:5f22683c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:08df83bd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:77a75075]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:67202c7c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:9215f574]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:4f2d21d7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:4f6a3ed0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:5cc5b89b]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:36f6ca17]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:84c54b88]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:9a353137]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0978 | ribonucleotide reductase\n","truncated":false}}
%---
%[output:7e5c19f2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:7a96fd08]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:37abd0f4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:0058c71e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:51eb82cf]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:48c1b2e3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:8347a9e3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:45668607]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:70b8bb3a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4749 | adenosine phosphorylase\n","truncated":false}}
%---
%[output:2e2a8db2]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:6a72c8c8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:48dda778]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:3c914305]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:50e64157]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0224 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:5fd11695]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:43f27157]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:849f6db8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1205 | NADH\n","truncated":false}}
%---
%[output:199b24db]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1205 | NADH\n","truncated":false}}
%---
%[output:63982503]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:95e51498]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1205 | NADH\n","truncated":false}}
%---
%[output:530a4e15]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1205 | NADH\n","truncated":false}}
%---
%[output:81b4a461]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:788d2e16]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:86f42587]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:4b6cb5ce]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:102c9e2c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4749 | adenosine phosphorylase\n","truncated":false}}
%---
%[output:83263e76]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:333766ca]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0978 | ribonucleotide reductase\n","truncated":false}}
%---
%[output:7d4f373a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:9c8bbe9d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:927372c7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:8d411675]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:9673bd5d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:5515e694]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:8cce911f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:2a33ffee]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:269a64ca]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:295c53f3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:1955e1c5]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:166aab94]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4749 | adenosine phosphorylase\n","truncated":false}}
%---
%[output:07ef2ee7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0978 | ribonucleotide reductase\n","truncated":false}}
%---
%[output:86c8352d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1205 | NADH\n","truncated":false}}
%---
%[output:3cb8c61e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1205 | NADH\n","truncated":false}}
%---
%[output:590598d1]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1205 | NADH\n","truncated":false}}
%---
%[output:39356964]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:79c1c6b0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1205 | NADH\n","truncated":false}}
%---
%[output:717685a7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:43c867b3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:510c1363]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:3af9c347]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:363c7ea9]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:23a5c94d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:950136a6]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:5dde3e93]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0224 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:04cfde83]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:1041882b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4749 | adenosine phosphorylase\n","truncated":false}}
%---
%[output:07dacc93]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:4b1072b2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:9143408d]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:9fae1e7e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:3db89dd9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:73078948]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:422914af]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:8b1d1a12]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:55e20174]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:98823ca3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:1a47bf67]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:7543527c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0978 | ribonucleotide reductase\n","truncated":false}}
%---
%[output:43c439b4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:73085a5a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:290c9410]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0459 | galactose-1-phosphate uridylyltransferase\n","truncated":false}}
%---
%[output:3ab5bc3e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:5339764a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:7a92033a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0799 | H+\n","truncated":false}}
%---
%[output:4caf7f09]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:0bc506f1]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0929 | isoamylol\n","truncated":false}}
%---
%[output:174cc008]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0929 | isoamylol\n","truncated":false}}
%---
%[output:36a8e301]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0929 | isoamylol\n","truncated":false}}
%---
%[output:5edfa16c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0929 | isoamylol\n","truncated":false}}
%---
%[output:40768988]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:2d7ec5e4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:52c2d37b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0978 | ribonucleotide reductase\n","truncated":false}}
%---
%[output:712fbd5e]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:20266a4c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:6f7a6273]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:276e9073]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:3b9bef58]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0469 | glutamate decarboxylase\n","truncated":false}}
%---
%[output:9b0ac543]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:00db8793]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:86ac97e4]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:9bcf3b0f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:2f5dd8df]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:65279907]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0932 | isobutanol\n","truncated":false}}
%---
%[output:2fdfcd67]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0932 | isobutanol\n","truncated":false}}
%---
%[output:03494e22]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:9fead3ec]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0314 | cytidine deaminase\n","truncated":false}}
%---
%[output:3f2a9e4f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0932 | isobutanol\n","truncated":false}}
%---
%[output:59a72174]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:28d971ce]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0932 | isobutanol\n","truncated":false}}
%---
%[output:1f16d540]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:16fab257]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:9368a586]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:491b7153]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:49731d86]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0138 | adenine deaminase\n","truncated":false}}
%---
%[output:608cdfe2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:9b01c973]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:0526ce0f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:6bf45461]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:10e6482f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0974 | ribonucleotide reductase\n","truncated":false}}
%---
%[output:595eb7c9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YNL220W\n","truncated":false}}
%---
%[output:4f2d2b9a]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1207 | NADP(+)\n","truncated":false}}
%---
%[output:38f4b3c0]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:7c0c3b1a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:84e1832c]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:0da20a89]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:35107c89]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:7845a723]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0973 | ribonucleoside-triphosphate reductase (UTP)\n","truncated":false}}
%---
%[output:9fe0cdc8]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:86f6ad8f]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:9e2fa71e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:4be950d3]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_0794 | H+\n","truncated":false}}
%---
%[output:2d75c292]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YKL001C\n","truncated":false}}
%---
%[output:2b6c5bbd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0976 | ribonucleotide reductase\n","truncated":false}}
%---
%[output:033bfef1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:56576c3f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:8baded26]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:59b40a9f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0182 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:908b7497]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:53442daf]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:90118232]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:020c01c1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0989 | saccharopine dehydrogenase (NADP, L-glutamate forming)\n","truncated":false}}
%---
%[output:310f1237]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:60c9e291]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:1881e162]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YKL001C\n","truncated":false}}
%---
%[output:85d17995]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:65509623]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:27387850]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:24556e5b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0978 | ribonucleotide reductase\n","truncated":false}}
%---
%[output:554102e7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:68dcc750]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:4b20dd77]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:9bf85d19]
%   data: {"dataType":"text","outputData":{"text":"  m_c  = s_1203 | NADH\n","truncated":false}}
%---
%[output:2123174c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:7156d3ed]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:00fb54de]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:3b0bedae]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1069 | UDP-N-acetylglucosamine diphosphorylase\n","truncated":false}}
%---
%[output:9015c6a0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:00591a2b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YKL001C\n","truncated":false}}
%---
%[output:016d8d1e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0376 | acetyl-CoA\n","truncated":false}}
%---
%[output:0f8cf826]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1401 | pyruvate\n","truncated":false}}
%---
%[output:306cf3d8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:2c05d9b9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1401 | pyruvate\n","truncated":false}}
%---
%[output:14e5502d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1386 | PRPP\n","truncated":false}}
%---
%[output:966ba4f0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:797e2cc4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:23092ccd]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:8fb3215e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:5fcbfc03]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:867d98d9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:761ff593]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YKL001C\n","truncated":false}}
%---
%[output:0b772f08]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:10735288]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1084 | UTP-glucose-1-phosphate uridylyltransferase\n","truncated":false}}
%---
%[output:562276ad]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:1099755a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:07504aad]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:95c17e42]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:32e078a3]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:63d0dc05]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:13b8818c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:2a0869b7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:1924e140]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:46a4765b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:2f6edbdc]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1386 | PRPP\n","truncated":false}}
%---
%[output:57a78e81]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1386 | PRPP\n","truncated":false}}
%---
%[output:8086d95d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:42edc181]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:7ed63650]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1275 | UTP transport\n","truncated":false}}
%---
%[output:05fa5d50]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:57609bcc]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:4e116197]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:1c0af54d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:109cb368]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:10a8c6a6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:76570dc2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:5864f34b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:84b24152]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:2c17a8f7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:7cefcbd3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:650935b3]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:4270bbc0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:3e04e210]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:841b5f74]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:391d768b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1276 | UTP\/UMP antiport\n","truncated":false}}
%---
%[output:6918d3f2]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:2ea870b9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:40242774]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:2ac85343]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:25298172]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:85dd35f6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:8f27a628]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:0a1fa2b7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:82f2ca75]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:41970077]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:059c16d7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:165fb11f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:02dfd280]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:751ff263]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:5fd4b3ef]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0223 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:852af08b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:5aa84a6b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:50395dbd]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:2fcd94ef]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:36293012]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:2631f8f5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:99826d4e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:767c3c77]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0991 | L-glutamate\n","truncated":false}}
%---
%[output:8f7270de]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:20d4fb33]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:3209a1e2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:338c38bc]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:31880522]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:94de5bb2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:17f44f21]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0224 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:80a25b87]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:7c8cad34]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:8e42e9de]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1559 | UTP\n","truncated":false}}
%---
%[output:157011a3]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0785 | GTP\n","truncated":false}}
%---
%[output:90e7d3b4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0785 | GTP\n","truncated":false}}
%---
%[output:84a9d554]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0785 | GTP\n","truncated":false}}
%---
%[output:84848b80]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0785 | GTP\n","truncated":false}}
%---
%[output:449f86d6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0785 | GTP\n","truncated":false}}
%---
%[output:05d1cfe4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:74924819]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:19f705cd]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6febcc7e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:65499dda]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:77c20905]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:9b4e4bbe]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:57cb4f7d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0722 | mannose-1-phosphate guanylyltransferase\n","truncated":false}}
%---
%[output:9734b4b6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:297d6391]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:9347105f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4a08f445]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:241d70f7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1616 | TRX1\n","truncated":false}}
%---
%[output:7db09fdc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:13f2c075]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:8c71e29a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:86448e52]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:561d00cf]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0261b8b2]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:32fe4291]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0739 | GDP\n","truncated":false}}
%---
%[output:27467a5e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0543 | cytidine\n","truncated":false}}
%---
%[output:0d10e97c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0972 | ribonucleoside-triphosphate reductase (GTP)\n","truncated":false}}
%---
%[output:9cd34213]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1616 | TRX1\n","truncated":false}}
%---
%[output:413a725d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7bfe67c5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:043644fd]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1616 | TRX1\n","truncated":false}}
%---
%[output:40afbd65]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:32e550c1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:175105f3]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1616 | TRX1\n","truncated":false}}
%---
%[output:8a09c5e5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:5fded4d4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:47d4ec2b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4c92c148]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:305f8159]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:60659c6d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0e14c5c2]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6b430ea7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1175 | GTP\/GDP translocase\n","truncated":false}}
%---
%[output:7d4e00c0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:1f7c0c9f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:2842cedf]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:361e2ea9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:990d89cf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:9b4246af]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:93d7b3da]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:8dbf0ebc]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4a725329]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:537a6777]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:736dca4d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1515de37]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:953373d3]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:65570dfd]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:85c286a4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:8f773133]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0223 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:345b4929]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:049a1148]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:58b9a769]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5c7b6a5b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:935d4911]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:2e871347]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:65922c78]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:4b313b59]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4601cd9e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8a62b961]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1ae2c2c0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:5518494c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:46c12be6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:429b50b5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0999 | L-glutamine\n","truncated":false}}
%---
%[output:4317d4ff]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0224 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:4b498d28]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:0a223522]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:6711d9a0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:7c3dd49b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0180 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:73e00b0c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:8f82df58]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:20b4b6d1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:5afc5320]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:9a3006da]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:69b111c6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:297c7d35]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1322 | phosphate\n","truncated":false}}
%---
%[output:4f4888a9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6c7e01c7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:13888a8b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:5f7623ae]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0722 | mannose-1-phosphate guanylyltransferase\n","truncated":false}}
%---
%[output:228cc991]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:85a4671a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:9cbf7ccd]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0e9ec1ed]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:137345a4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:66efc5a8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:95430e00]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:68502e0a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:47f50b11]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0140 | adenosine deaminase\n","truncated":false}}
%---
%[output:67dfd722]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:6101aaf0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:9900d055]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:04633914]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5a6a5b27]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5132f348]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:8bbaeaa7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0722 | mannose-1-phosphate guanylyltransferase\n","truncated":false}}
%---
%[output:74b09fa8]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:3532d8f4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:6da9730e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:53220182]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:9f2a5703]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:296a6525]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:4b7e3288]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1f760055]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:81b3d865]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:86bdbe83]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0142 | adenosine kinase\n","truncated":false}}
%---
%[output:76147e6c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:53f78135]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:5cd59be2]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:990570f2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:86dec0b2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:5f8ad696]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0972 | ribonucleoside-triphosphate reductase (GTP)\n","truncated":false}}
%---
%[output:8c3ade37]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:830ba580]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:746c2cc8]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0025 | (R)-lactate\n","truncated":false}}
%---
%[output:7932bab8]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:43cdcec2]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:5de07ba7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:79693faf]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:177b0974]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:2788ac7c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:7054bc5f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:486e50df]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:6d7a8c28]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:35420804]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0142 | adenosine kinase\n","truncated":false}}
%---
%[output:14f956e7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:2653a610]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0972 | ribonucleoside-triphosphate reductase (GTP)\n","truncated":false}}
%---
%[output:7b923018]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:3de323c6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6b2fa59c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5d155a74]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:7201ee03]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:669efd3d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:98aeed7b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0ac696e3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:45cda27d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:6d5b8a68]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:888743e0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:7d57a355]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:699e62c0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0142 | adenosine kinase\n","truncated":false}}
%---
%[output:5db1a9c4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:47427404]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0972 | ribonucleoside-triphosphate reductase (GTP)\n","truncated":false}}
%---
%[output:993c6931]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:991be7dc]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:66e6e3b2]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:106ed596]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:0af6b221]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1ff1fb64]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:4449cf59]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:86f92199]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:08679e98]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:5d6a8c74]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:3f793f4c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6341270d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0bfc3e5e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:08f4db37]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1175 | GTP\/GDP translocase\n","truncated":false}}
%---
%[output:969711a9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:30f6ed94]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0142 | adenosine kinase\n","truncated":false}}
%---
%[output:8839ac44]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:682c5025]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:69fcdf80]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:222bf5e9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:9400462c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:89dedec3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:6915dd6e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1043f5ad]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:96039493]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:42943ed3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:3d1b0428]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:21f89dae]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:24d71fb6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:2c595fd3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1175 | GTP\/GDP translocase\n","truncated":false}}
%---
%[output:32239d46]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:1b9e9d14]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1bfff1ee]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:541863d0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:887b999a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:318355ee]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:55b2380c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:7f66ef14]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:3ea9feab]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:9b569f9c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4ddbe073]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:2e5b9d19]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6a7ea151]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:47519b27]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:579ee851]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0224 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:1b4fc589]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:82caaf66]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:41042c5a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:0186f782]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:24b0ee7d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:4658e8b2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:8cab2090]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:7ab46c37]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4717d24d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:61b1f72f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:0b182de3]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:401e9827]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6549e7e4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:46243af1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:9ae1a1c5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:64d222b7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0315 | cytidine kinase (GTP)\n","truncated":false}}
%---
%[output:34b651d7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4fe9691f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0471 | glutamate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:2e634c56]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:2cd8cc14]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8e05d32b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:2dba0354]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:3d338645]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1fe391f9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:5fc8009d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8c5cc371]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:424ca84e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:3ff0463d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:47c2bd0d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:66ea54a2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0972 | ribonucleoside-triphosphate reductase (GTP)\n","truncated":false}}
%---
%[output:7871e7a8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:3000cd9b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5bcf23b0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:98a81625]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:79460cae]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:09552365]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:3282d665]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:757db702]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:81c24b18]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:0fd5ff6a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0d94dcba]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:2af1c3bf]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:98a342c8]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:95bb4093]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:313352d1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0088 | 5PP-IP5 pyrophosphorylation to 4,5-PP2-IP4\n","truncated":false}}
%---
%[output:8eabf50f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0972 | ribonucleoside-triphosphate reductase (GTP)\n","truncated":false}}
%---
%[output:89879ff9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:33678a0f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:033128ec]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:38134e62]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6b2a982b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:7357db6b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:35391af0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0cdce167]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:08d42a6c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:7fe09a6f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5d1821c1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6dadfe07]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0089 | 5PP-IP5 pyrophosphorylation to 5,6-PP2-IP4\n","truncated":false}}
%---
%[output:73701fd9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:31f8eec0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:63f01a3c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:83b542b6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0972 | ribonucleoside-triphosphate reductase (GTP)\n","truncated":false}}
%---
%[output:7a69d12b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:7d663a68]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:177c0085]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:7c6961f3]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8150ca60]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4362622f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1eab916f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:61b7255e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:006b7d49]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:8d6b83de]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:0cb924db]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0358 | diphosphoinositol-1,3,4,6-tetrakisphosphate synthase\n","truncated":false}}
%---
%[output:78a7d36f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:59608aa1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:2268fa9f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:3f51d8ef]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0223 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:66839862]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:256173a4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8c9ad7a0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:36683f27]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:2878aa7b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:06c312fd]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:41ef011f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:0f292aa1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:11510ea1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0404a6c6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_1245 | phosphate transport\n","truncated":false}}
%---
%[output:7778fa4b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:47311a3f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:7cb31358]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1ba6f504]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0223 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:25ad980d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:6a9a85f0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6be39167]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:392f5e13]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0931bdad]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6085c791]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:8099d5d5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6e43ae40]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:9e44e5a2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:8eabf12f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:052f0b5a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:851ebd44]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:6707f046]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:4bba9f25]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4b7c3e6f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:32f0b3d5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0223 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:2949901d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:94750650]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:00deabcb]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:18bdac07]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:469037c9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8f6947db]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:0467a03a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:3ad80065]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:60fc8806]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6b08695d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:58606d62]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:68f68bf5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4561ce6b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:91ed34c2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:00588c55]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0223 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:90e47644]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:99d5bd06]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:622408b6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:2b433293]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:9f09b0ab]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:799dfef4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:68d29d23]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:59deb958]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:1be086de]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6048f989]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:6021bd15]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:0cad9cf5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:6f42dac9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8ff2fb1a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:5c3f3d2e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:9dfa21b8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0224 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:416a46e4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:073e934b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:8ec14757]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:2eca8274]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:932c13d3]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:00443f90]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:940de893]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:528e8955]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:6b622c14]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4f8961bf]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:3e534b3c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:4207ffc4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:228a9b3d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:168b08a2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0224 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:92d3c7d6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:784b19b1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:732589f3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:4e7d4204]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:9fcc6df9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:1eb0b2a5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:6d2328f1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:791d3519]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0e3974be]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0306 | 5,10-methylenetetrahydrofolate\n","truncated":false}}
%---
%[output:6eec7c0b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:4012bb70]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1377 | prephenate\n","truncated":false}}
%---
%[output:6c8a069d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:5d79d1f0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:0a66dd83]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:0007af8a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0224 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:9244e18e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:2ab46e6d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:4cbf7584]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:41a7d503]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:2b5a186a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:76a991f2]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:822f5d23]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:29392050]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:3e14b8bb]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:100e3cf9]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0306 | 5,10-methylenetetrahydrofolate\n","truncated":false}}
%---
%[output:9a25b726]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:35e55fa0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1377 | prephenate\n","truncated":false}}
%---
%[output:3968f335]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:3a6427c1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:87f8bbaf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:26cba808]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:1234bcb4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:7fe28aa9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0224 | ATP adenylyltransferase\n","truncated":false}}
%---
%[output:8a499b09]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:09fd7090]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:64bd7e32]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:15a04c13]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:782d8fe7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:17240c2c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:22a2e1c9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:14400900]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0306 | 5,10-methylenetetrahydrofolate\n","truncated":false}}
%---
%[output:0d452ecd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:6a19ddac]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:2a12cfcf]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:1e37b630]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:6daf0749]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1377 | prephenate\n","truncated":false}}
%---
%[output:63d63ac8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0722 | mannose-1-phosphate guanylyltransferase\n","truncated":false}}
%---
%[output:30153d8b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:56822d5a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:2ae93401]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:072eedfa]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:44734c37]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:8e070853]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:696b0094]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:27371426]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0143 | adenosine monophosphate deaminase\n","truncated":false}}
%---
%[output:00d630db]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:0aaf942e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1399 | pyruvate\n","truncated":false}}
%---
%[output:626ac72d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:798d14ec]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:524928c7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:12c51791]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0722 | mannose-1-phosphate guanylyltransferase\n","truncated":false}}
%---
%[output:189f1913]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:7d97bf50]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:3457c9b7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:360c42be]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:48a09627]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:9695a6a1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:8db21f6e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:7a4f7c5e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0306 | 5,10-methylenetetrahydrofolate\n","truncated":false}}
%---
%[output:81f608f7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1377 | prephenate\n","truncated":false}}
%---
%[output:6f02a1aa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:0a4397a1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:0347b85c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:6cb17c15]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:3584ae6b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:2433509d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:6768c6d0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0722 | mannose-1-phosphate guanylyltransferase\n","truncated":false}}
%---
%[output:566d1878]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:9c3012c0]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:905d46b6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:134c7722]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0799 | H+\n","truncated":false}}
%---
%[output:4ed80847]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:97855979]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0799 | H+\n","truncated":false}}
%---
%[output:9ba08091]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0799 | H+\n","truncated":false}}
%---
%[output:88836c64]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:5740c392]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:968dacab]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0799 | H+\n","truncated":false}}
%---
%[output:7af00410]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:3d7fb2f7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:3b24bcec]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:1fab0bd8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0722 | mannose-1-phosphate guanylyltransferase\n","truncated":false}}
%---
%[output:20cc46a2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:868c9757]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:37d24675]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:2be9bda1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:9b8aa9c0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:34c1bf65]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:4b14734a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:6da99e38]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:3e91ab8f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:392aace5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:7cf5ce99]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:2d771a8f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0799 | H+\n","truncated":false}}
%---
%[output:1852dc51]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0799 | H+\n","truncated":false}}
%---
%[output:00c791a6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:25ad883f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:5cf7cdfd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:1a8aa862]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0972 | ribonucleoside-triphosphate reductase (GTP)\n","truncated":false}}
%---
%[output:035436b8]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0799 | H+\n","truncated":false}}
%---
%[output:50edf1c4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0799 | H+\n","truncated":false}}
%---
%[output:17e50f4e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:6b51a974]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7558289f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:10f0dedd]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:7f298104]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:5127a442]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:8ebfb962]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:98583d82]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:83430433]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:2e884947]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:4aa657e1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0306 | 5,10-methylenetetrahydrofolate\n","truncated":false}}
%---
%[output:8c7454e4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:0277c903]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0972 | ribonucleoside-triphosphate reductase (GTP)\n","truncated":false}}
%---
%[output:871c31af]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:819230de]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1377 | prephenate\n","truncated":false}}
%---
%[output:76352dde]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:5656626e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:44005adb]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:9bd2fcca]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:1cd47f78]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:330726ad]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:82c127b4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:0be398f7]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1205 | NADH\n","truncated":false}}
%---
%[output:69652af1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:1670a534]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:50e2ae76]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:0ab1a0ed]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0972 | ribonucleoside-triphosphate reductase (GTP)\n","truncated":false}}
%---
%[output:9f1d1085]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:1029dbfc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:14af12df]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:6c4b1e3c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:54480bb3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:959ab479]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:6de593d6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0306 | 5,10-methylenetetrahydrofolate\n","truncated":false}}
%---
%[output:64835c8f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:70fb325d]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1377 | prephenate\n","truncated":false}}
%---
%[output:2517942f]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:643f4110]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:880a7d57]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:80707306]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:32684f25]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:48584ea8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:0e8ad17a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0972 | ribonucleoside-triphosphate reductase (GTP)\n","truncated":false}}
%---
%[output:53bbbb36]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:9179d770]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:502c8bcb]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:2155a040]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0373 | acetyl-CoA\n","truncated":false}}
%---
%[output:53e94600]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:3e469d47]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0306 | 5,10-methylenetetrahydrofolate\n","truncated":false}}
%---
%[output:5d2717a6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0850 | indol-3-ylacetaldehyde\n","truncated":false}}
%---
%[output:2dd4e662]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1377 | prephenate\n","truncated":false}}
%---
%[output:4f858b59]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:526df449]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:527afa14]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:3028218a]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0359 | acetaldehyde\n","truncated":false}}
%---
%[output:1df2df2c]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:02b45064]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:50c3e67c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:75c6e2e3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1175 | GTP\/GDP translocase\n","truncated":false}}
%---
%[output:62ddff27]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:573843d4]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:1bfe81f5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7bea5766]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:7e6e3814]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:5447d571]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:567e9ce5]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:2e61e727]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_1203 | NADH\n","truncated":false}}
%---
%[output:02f53d13]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:9b7399e2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:95bf4f63]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8493a6f2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:1c5c27a4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:9e99586e]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:519fda88]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1175 | GTP\/GDP translocase\n","truncated":false}}
%---
%[output:277ff934]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:2af296de]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:3aebab70]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:8011f8e6]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:0ef6f547]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:9a97b80f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:6b485795]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:24730319]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:2faaeec1]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:2edaf821]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:7ed870c2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:7aea7290]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:62e6790b]
%   data: {"dataType":"text","outputData":{"text":"  m_Cs = s_0794 | H+\n","truncated":false}}
%---
%[output:93a03a5a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:4581891b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1175 | GTP\/GDP translocase\n","truncated":false}}
%---
%[output:952247ae]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:1477760e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:0b0e1090]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:7df9862d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:01543104]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:1bb33e93]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:320a0002]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:6721bb78]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:1e84c977]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:1b85f343]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4236 | (R)-lactate hydro-lyase\n","truncated":false}}
%---
%[output:0d64717c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:0641c1b4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:52e989ee]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:5101d99b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:45dc3a0e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1175 | GTP\/GDP translocase\n","truncated":false}}
%---
%[output:3b7f03f2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:50d2133c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:5fd98075]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:45c664b3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:5447bf87]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:9397cf4d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:3815c0ac]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:8432d15c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:31642efb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0959 | pyruvate decarboxylase\n","truncated":false}}
%---
%[output:2e6f265f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:389fa0e7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:2328d1b4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:2548f1b6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:381c8fbb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:17b0fc0c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:79cbb026]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0211 | asparagine synthase (glutamine-hydrolysing)\n","truncated":false}}
%---
%[output:79ca8bff]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:52c9d934]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:107433aa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:154c1ba8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:8b548875]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:879ae7b6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:75c0eeb8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:3597a13b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:7559cdbe]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:890b0ded]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:118e5bad]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:4e7e3558]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:95d70951]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:75d9a855]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:600ae140]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0216 | aspartate transaminase\n","truncated":false}}
%---
%[output:753f5a2e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0153 | adenylosuccinate synthase\n","truncated":false}}
%---
%[output:6f12363c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0154 | adenylyl-sulfate kinase\n","truncated":false}}
%---
%[output:4357f316]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0154 | adenylyl-sulfate kinase\n","truncated":false}}
%---
%[output:89b0a21b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2034 | pyruvate transport\n","truncated":false}}
%---
%[output:727e7828]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0154 | adenylyl-sulfate kinase\n","truncated":false}}
%---
%[output:61ba9f4a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:189e472f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:31eff13c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0154 | adenylyl-sulfate kinase\n","truncated":false}}
%---
%[output:7e49ac2f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:596b991a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:6636682e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:6488ad7a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:4daa28e6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:97569e86]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0216 | aspartate transaminase\n","truncated":false}}
%---
%[output:441bdd92]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:581b3617]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:55253c78]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:194a240d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:039e3cdc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:3bae1260]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:218fc89e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:18a16659]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:70e3cfb1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:229c86d4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:677bcfc6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:57b6daf3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:6e033684]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:4dbf1b8f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:86d8fca5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:486fa5f8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1026 | sulfate adenylyltransferase (ADP)\n","truncated":false}}
%---
%[output:1ba823a0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:59395256]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:05b4ae21]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:7133743f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:77a7a494]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:7bb69dab]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:3c292ff1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:0da47a7e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:508e9a4c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:1c180562]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:8ce17f1d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:9a42a3e1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:9d50b394]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:6a7fb119]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1026 | sulfate adenylyltransferase (ADP)\n","truncated":false}}
%---
%[output:4d42c96e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:65032e78]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:74ed9d9b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7d9e7f10]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:28ff98c7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:8beb7463]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:23bd3a78]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:14ca1706]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:1e45f293]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:77c3a373]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:37f6d7bc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:728f58fe]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:2a0e1bfe]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:3561fc24]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:3ed5f22b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:4e40afdd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:4c054e60]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1026 | sulfate adenylyltransferase (ADP)\n","truncated":false}}
%---
%[output:1e119f98]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:7e5698df]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:9853daa4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:259103e0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:9b05435f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:73c742ef]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:1a62bdd3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:239de070]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:75c55551]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:84737b11]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:25df58dd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:3e5e8f25]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:65e4fb30]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:7994c4da]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:8be7ff00]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_1026 | sulfate adenylyltransferase (ADP)\n","truncated":false}}
%---
%[output:18145e65]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:37e93bf2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:08bcf4c9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:770b3df2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:48f55f49]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:33a4d49e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:00c2472e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:5094a13c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:0fd8ded1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:2262280a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:98e1a147]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:20a94f4b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:3cf17ebb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:1d8af207]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:55a75d8e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:4329d3bc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:24039353]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:29d4f4ee]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:3fd46c54]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:50701ce3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:0aa9dfc0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:5f3227a5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:1a9cd918]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:84e180d3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:155bdb14]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:60bf9001]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:585bacb8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:79366e7b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:6ad7645b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:3dd188b0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:1b8f2037]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:1f4f4185]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:53ef2afa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:8e60c6d6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:1a1439b5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:77447882]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:972acfd4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:863e6489]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:27e6e061]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:2632669c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:2518fe01]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:3e4bc72b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:900138e6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:5dcc9b1e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:8c1ffc83]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:74cb3ea5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:9a2cbfed]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:63e2288f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:27540b19]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:1e9b57f6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:433db78e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:49574648]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:7d6b94c6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:2462d50d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:9b99c78e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:403e15f1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:06d84711]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:4458e01a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:2e3d8f95]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:3adea439]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:4c4f7fcc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:0ef0eed4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:290b8489]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:56671865]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:622c9bcd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:5a23797a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:2c7c8354]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:8470db03]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:1ec704f7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:212b8189]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:5a66f0a1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:99285de2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:05e37099]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:02d1aea6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:79165323]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:30bf1acc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:42733788]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:2ad2d9d7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:03fa98a1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:4c6a361c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:715d7728]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:2acc8372]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:7c87238f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:673b8f00]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:3abb4eb2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:70f71b19]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:112f9527]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:1126bdb0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:3216c2af]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:55ce52e9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:32635763]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:2f053480]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:852c370f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:09f7b1d6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:5502b76e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:968c1a89]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:32059e71]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:41dbb4df]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:82847765]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:3d29f4c3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:997e70a0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7629875f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:3b3381fd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:4a6bb5eb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:77461bde]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:63e060bd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:7f0afe56]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4652 | aldehyde dehydrogenase (1-propanol, NAD)\n","truncated":false}}
%---
%[output:3fd28e93]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:4a28c77e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:36379202]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:3ae9a81b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR204W\n","truncated":false}}
%---
%[output:08857e07]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:52f9e7dd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:2a8cb659]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:6a5c5964]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:6e74992f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:5861f8d1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:60ecbd04]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:157119a4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:2258bc7a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:7809fec5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:60f7fb02]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4661 | aldehyde dehydrogenase (methionol, NAD)\n","truncated":false}}
%---
%[output:292860a7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:935cdf45]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:4698e948]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:74d66d85]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:120fa9c6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:004c6b37]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:52f1b7e0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:0e58077c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:78f00328]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR166C\n","truncated":false}}
%---
%[output:7161fea0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:94e56575]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:2a2ea5f7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:0442126c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4670 | aldehyde dehydrogenase (tyrosol, NAD)\n","truncated":false}}
%---
%[output:408783f0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:829fa097]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:7548aefb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:02e5d532]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:8507bd02]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:4a2ced72]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:3c44993e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:13a4e321]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:43f16d75]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:13e84767]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:88f5a937]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:2a31f5dc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:57298c03]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:29fdabc3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:68a7db27]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:16949b7a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:14dcaa36]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:6ed4ea29]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:09bf4a41]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:21853e51]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:204f654b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:2618260c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:164db8b8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:669fda31]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:3ef2e1eb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:55e314ae]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:5bdea67e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:581cb38d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:36267989]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:28afc39a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:78a227fc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:1b656929]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:22149dcb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:11807403]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:90c72559]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:85cd2d8e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:0f54cee8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:919a356d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0163 | alcohol dehydrogenase (ethanol to acetaldehyde)\n","truncated":false}}
%---
%[output:8a7f913b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0164 | alcohol dehydrogenase (glycerol, NADP)\n","truncated":false}}
%---
%[output:223c18ac]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:10be0c01]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:6b1f8df8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0164 | alcohol dehydrogenase (glycerol, NADP)\n","truncated":false}}
%---
%[output:8297abdc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0164 | alcohol dehydrogenase (glycerol, NADP)\n","truncated":false}}
%---
%[output:6c9bc5e9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:4a395970]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0164 | alcohol dehydrogenase (glycerol, NADP)\n","truncated":false}}
%---
%[output:8f109f19]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:6383fb17]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0164 | alcohol dehydrogenase (glycerol, NADP)\n","truncated":false}}
%---
%[output:888d01d8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:828b775a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:082e92cf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:98ead2ff]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:17b9b2f5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0168 | aldehyde dehydrogenase (2-methylbutanol, NADP)\n","truncated":false}}
%---
%[output:30eedcb2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0168 | aldehyde dehydrogenase (2-methylbutanol, NADP)\n","truncated":false}}
%---
%[output:1c1fa716]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:586740a8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0168 | aldehyde dehydrogenase (2-methylbutanol, NADP)\n","truncated":false}}
%---
%[output:4fb290d9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:09c74c0f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:273c3548]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0168 | aldehyde dehydrogenase (2-methylbutanol, NADP)\n","truncated":false}}
%---
%[output:56ebf098]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0168 | aldehyde dehydrogenase (2-methylbutanol, NADP)\n","truncated":false}}
%---
%[output:0227cfa0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0168 | aldehyde dehydrogenase (2-methylbutanol, NADP)\n","truncated":false}}
%---
%[output:18367dbd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0168 | aldehyde dehydrogenase (2-methylbutanol, NADP)\n","truncated":false}}
%---
%[output:951ec690]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:2439bf58]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:4b2d4ddb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:19458ab2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:6a22151a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR204W\n","truncated":false}}
%---
%[output:15a09977]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0171 | aldehyde dehydrogenase (2-phenylethanol, NADP)\n","truncated":false}}
%---
%[output:98bcffe6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:74d5cca5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:6db91702]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0171 | aldehyde dehydrogenase (2-phenylethanol, NADP)\n","truncated":false}}
%---
%[output:4f1ac81a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0171 | aldehyde dehydrogenase (2-phenylethanol, NADP)\n","truncated":false}}
%---
%[output:465c56b0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0171 | aldehyde dehydrogenase (2-phenylethanol, NADP)\n","truncated":false}}
%---
%[output:917c863d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:855bfac5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0171 | aldehyde dehydrogenase (2-phenylethanol, NADP)\n","truncated":false}}
%---
%[output:1c449c5f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0171 | aldehyde dehydrogenase (2-phenylethanol, NADP)\n","truncated":false}}
%---
%[output:14116a56]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:430b36f2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0171 | aldehyde dehydrogenase (2-phenylethanol, NADP)\n","truncated":false}}
%---
%[output:82773c5e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:9ddd4e05]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:3795ec72]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:7524fa96]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:28e879be]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:3b110b09]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:945efed4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:7e8da264]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:19626f5a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR166C\n","truncated":false}}
%---
%[output:4d9235a3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:3746709a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:10c3b6fe]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:781f3705]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:325fd73f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:57ef70b9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:65ebd627]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:2fcea368]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:3250e04f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:57c0c463]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:8b23d1b7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:59ce63d2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:97514f4c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:6a2a3885]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0174 | aldehyde dehydrogenase (acetylaldehyde, NAD)\n","truncated":false}}
%---
%[output:9a0ba7e5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0174 | aldehyde dehydrogenase (acetylaldehyde, NAD)\n","truncated":false}}
%---
%[output:6ec1226a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0174 | aldehyde dehydrogenase (acetylaldehyde, NAD)\n","truncated":false}}
%---
%[output:7b34277a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0174 | aldehyde dehydrogenase (acetylaldehyde, NAD)\n","truncated":false}}
%---
%[output:1f02be7c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0174 | aldehyde dehydrogenase (acetylaldehyde, NAD)\n","truncated":false}}
%---
%[output:3dee9eb7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:562fe160]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:2ff34beb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:4970f1de]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0174 | aldehyde dehydrogenase (acetylaldehyde, NAD)\n","truncated":false}}
%---
%[output:22b35873]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:4584bb74]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:587a456e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0174 | aldehyde dehydrogenase (acetylaldehyde, NAD)\n","truncated":false}}
%---
%[output:8775cf0d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0174 | aldehyde dehydrogenase (acetylaldehyde, NAD)\n","truncated":false}}
%---
%[output:0af60b05]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0175 | aldehyde dehydrogenase (acetylaldehyde, NADP)\n","truncated":false}}
%---
%[output:9f988f96]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0175 | aldehyde dehydrogenase (acetylaldehyde, NADP)\n","truncated":false}}
%---
%[output:5aa64547]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0175 | aldehyde dehydrogenase (acetylaldehyde, NADP)\n","truncated":false}}
%---
%[output:1844dffe]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:1684026f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0175 | aldehyde dehydrogenase (acetylaldehyde, NADP)\n","truncated":false}}
%---
%[output:2410ef07]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:0642cb04]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:19ebba3d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0176 | aldehyde dehydrogenase (indole-3-acetaldehyde, NAD)\n","truncated":false}}
%---
%[output:23934dd4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0176 | aldehyde dehydrogenase (indole-3-acetaldehyde, NAD)\n","truncated":false}}
%---
%[output:7f835af4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0176 | aldehyde dehydrogenase (indole-3-acetaldehyde, NAD)\n","truncated":false}}
%---
%[output:84fc702d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:310faf1d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:8069be55]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0176 | aldehyde dehydrogenase (indole-3-acetaldehyde, NAD)\n","truncated":false}}
%---
%[output:9a76a14b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:19759201]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0176 | aldehyde dehydrogenase (indole-3-acetaldehyde, NAD)\n","truncated":false}}
%---
%[output:2780284d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0176 | aldehyde dehydrogenase (indole-3-acetaldehyde, NAD)\n","truncated":false}}
%---
%[output:6ea84fe2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:09a52fec]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0176 | aldehyde dehydrogenase (indole-3-acetaldehyde, NAD)\n","truncated":false}}
%---
%[output:10691113]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:2d2c0617]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0176 | aldehyde dehydrogenase (indole-3-acetaldehyde, NAD)\n","truncated":false}}
%---
%[output:52911ca6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:8fc9992f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:3ac70cf6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:9449e096]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:61433de5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:960be5cc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:5740fd4a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:49a8812e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:17334f2a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:7e6f88bd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:23dcab47]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:7aec911d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:15a401e0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:6b67bfda]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:85a69172]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0178 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:41868ece]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0178 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:8c79d1bb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0178 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:8458d661]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:6e4700b6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:88cb7b6b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0178 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:9728b0d8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:7bfceca4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0159 | alcohol acetyltransferase (ethanol)\n","truncated":false}}
%---
%[output:8bb6d2df]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:5f0b414b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR204W\n","truncated":false}}
%---
%[output:7f9ea80b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:378e85db]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0181 | aldehyde dehydrogenase (isoamyl alcohol, NADP)\n","truncated":false}}
%---
%[output:5cab25ad]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0181 | aldehyde dehydrogenase (isoamyl alcohol, NADP)\n","truncated":false}}
%---
%[output:153f9d57]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0181 | aldehyde dehydrogenase (isoamyl alcohol, NADP)\n","truncated":false}}
%---
%[output:01f8cfe2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0181 | aldehyde dehydrogenase (isoamyl alcohol, NADP)\n","truncated":false}}
%---
%[output:9a6daeed]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:376a4837]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0181 | aldehyde dehydrogenase (isoamyl alcohol, NADP)\n","truncated":false}}
%---
%[output:7d6cbc36]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:0d245b0e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0181 | aldehyde dehydrogenase (isoamyl alcohol, NADP)\n","truncated":false}}
%---
%[output:2d001e90]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:258cbbf0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0181 | aldehyde dehydrogenase (isoamyl alcohol, NADP)\n","truncated":false}}
%---
%[output:02ec6fa3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:441f962a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:41ed669b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:1c1766d5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:5c92e797]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0184 | aldehyde dehydrogenase (isobutyl alcohol, NADP)\n","truncated":false}}
%---
%[output:7f2824c9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:028ecf05]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0184 | aldehyde dehydrogenase (isobutyl alcohol, NADP)\n","truncated":false}}
%---
%[output:64b08b37]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0184 | aldehyde dehydrogenase (isobutyl alcohol, NADP)\n","truncated":false}}
%---
%[output:0577262a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR166C\n","truncated":false}}
%---
%[output:4abf8608]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0184 | aldehyde dehydrogenase (isobutyl alcohol, NADP)\n","truncated":false}}
%---
%[output:33254083]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4652 | aldehyde dehydrogenase (1-propanol, NAD)\n","truncated":false}}
%---
%[output:6f834bae]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0184 | aldehyde dehydrogenase (isobutyl alcohol, NADP)\n","truncated":false}}
%---
%[output:95515d2b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:65beeb3a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0184 | aldehyde dehydrogenase (isobutyl alcohol, NADP)\n","truncated":false}}
%---
%[output:3f847acd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0184 | aldehyde dehydrogenase (isobutyl alcohol, NADP)\n","truncated":false}}
%---
%[output:418d10ef]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:2c4b54c6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:24871ea0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:507fe250]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:33ef2632]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:2964fc92]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:032865f6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:4bf1d927]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:015e541f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4661 | aldehyde dehydrogenase (methionol, NAD)\n","truncated":false}}
%---
%[output:54700101]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:4de659ee]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:757be0bb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:3b52104f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:72aa28bb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:2c3a67c0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:4a020190]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:94a5a67a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:61b1cd35]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:4903feb4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:798e14fe]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:43275654]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:7eff6d06]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:8b55b02c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4670 | aldehyde dehydrogenase (tyrosol, NAD)\n","truncated":false}}
%---
%[output:55d8284d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:07444ff8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:00536189]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:8a525ca7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:2cfe2231]
%   data: {"dataType":"text","outputData":{"text":"  r_Cp = r_0185 | aldehyde dehydrogenase (phenylacetaldehyde, NAD)\n","truncated":false}}
%---
%[output:212d5369]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:9dbd0074]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:485e16f8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:053bf4a2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:25036e7c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:8d1e4816]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:93be624e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:0385df11]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:8cec23a0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:53a4f48e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:8fb98f52]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:84b91028]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:718de745]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:6d34af46]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:6687740c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:8be616a4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YGL080W and YGR243W) or (YGL080W and YHR162W)\n","truncated":false}}
%---
%[output:4a51aef5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:695eb2da]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:41988b51]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:2f9890c3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:65cf6bcb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:40c96c1c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:3ea7d2f7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:2eb3b57f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4711 | ethyl-(2S)-lactate esterase, c\n","truncated":false}}
%---
%[output:4e5104f2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:2931141d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOR126C\n","truncated":false}}
%---
%[output:0a989154]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:0f324cf3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:41d3b162]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:8c0b7f1c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:4a498d17]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:0779382b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:57fe7f47]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:23f500d9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:0c8785ae]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:5b00e852]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:8303e6e6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:6d57129a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:893431a7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:25925afd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:319b9f10]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:1f569faa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:35c4906c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:81506b21]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:0eda4363]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:28db04a6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:1906a715]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4712 | ethyl-(2R)-lactate esterase, c\n","truncated":false}}
%---
%[output:23b1e420]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:4748dcab]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:70e25223]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:6c0690ec]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:1002af2a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:59c84627]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:0aaf5ac3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:7b6be62f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:215f5b78]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:25a1fece]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:4066bd48]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:3e9a85b1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:00e4abee]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:0c8514f4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:14d597a3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:165fa37b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:85761337]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:65282246]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:9f85ce8e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:8cab9302]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:56969891]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:03176f71]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:7a5a7375]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4713 | diethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:6030bbda]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:298084ae]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:820ba7c9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:85f424c9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:3666d0d9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:0a01730e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:4f549e01]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:5c2ab548]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:5fefd0d6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:9e1a534b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:4a5d026a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:47c98a4f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:09ff383e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:89fa86fb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:36a3f934]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:814ab8e3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:140975fc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:8543b208]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:4b1ad3e2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:235cc1f2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:46b5436d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:920156cd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4714 | monoethyl-succinate esterase, c\n","truncated":false}}
%---
%[output:325d40d1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:8f4cba31]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:57164c77]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:777f4e51]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:0bfc58fc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:8ce864be]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:4d6ab2bb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:5051f439]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:3db48a68]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:8fdbd353]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:4440b5c3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:00796948]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:72295523]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:2de45dc7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNL316C\n","truncated":false}}
%---
%[output:6b200643]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:04cba7db]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:7debd99d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:7f59b140]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:5143d211]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:763ee520]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:25452f2a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:8895254e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4715 | ethyl-benzoate esterase, c\n","truncated":false}}
%---
%[output:57bdc2e7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:5b353517]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOL059W\n","truncated":false}}
%---
%[output:9c0d1e53]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:34904b6c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:0f74106f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:3d7c3567]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:3f591147]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:86bb41dd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:21e32ea1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:08ee3074]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:56fa33f9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:5ae8aba0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YML120C\n","truncated":false}}
%---
%[output:4c2ecfc1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:184ca560]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:1c34d055]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:87235cb0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:481657da]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:1e43415a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:79158a91]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:9537ee5b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:84ca9418]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:14e0a57e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR204W and YGR255C and YLR201C and YML110C and YOL096C and YOR125C\n","truncated":false}}
%---
%[output:5fedf2da]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:9e92e0ac]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:95ea053f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:8115d309]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:07bf16fc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:8f18c823]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:0a2dd6d6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:839d1460]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:8b08cd40]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:7de4f265]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR051W\n","truncated":false}}
%---
%[output:3f55bbb3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:832fd064]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:051f844b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:2b2c052b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:70e562c9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:3e819d54]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:4301c189]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:1fb1506a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:4d9e75d7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:5021fc8e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOL059W\n","truncated":false}}
%---
%[output:146cdf71]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:8f14604d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:6bcd72cf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:0868bef2]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:80e5921b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4720 | ethyl-isobutyrate esterase, c\n","truncated":false}}
%---
%[output:069a9115]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:7a7853bd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:711d701c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:8980b18f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:37b1271b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YML120C\n","truncated":false}}
%---
%[output:3b79f521]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:86c35573]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:6a418232]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:3151f38a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:2d3209cb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:47182a1c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:74326eae]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4652 | aldehyde dehydrogenase (1-propanol, NAD)\n","truncated":false}}
%---
%[output:0aadbef7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:2dc4940b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:8c671179]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR204W and YGR255C and YLR201C and YML110C and YOL096C and YOR125C\n","truncated":false}}
%---
%[output:44f25d52]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:8450ecbb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:3257f413]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:52161fee]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:03f7ea12]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:796baf51]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:8cd78cde]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_4721 | ethyl-2-methylbutyrate esterase, c\n","truncated":false}}
%---
%[output:24d2f831]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:5e706948]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:18c5326a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR051W\n","truncated":false}}
%---
%[output:87b769b9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:582efaa7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4661 | aldehyde dehydrogenase (methionol, NAD)\n","truncated":false}}
%---
%[output:78984899]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0938 | prephenate dehydratase\n","truncated":false}}
%---
%[output:9d19d009]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:9086769d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:3668fd60]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:849c23e0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:1eae5ec6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:9d7701f4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:1319d506]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOL059W\n","truncated":false}}
%---
%[output:5ebe84f4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YOR126C\n","truncated":false}}
%---
%[output:1170fc24]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:74214652]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:565b37e7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:77bcf499]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4670 | aldehyde dehydrogenase (tyrosol, NAD)\n","truncated":false}}
%---
%[output:5ca92c1b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0938 | prephenate dehydratase\n","truncated":false}}
%---
%[output:37346636]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:47a1bacb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:161deacc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:9790f03d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YML120C\n","truncated":false}}
%---
%[output:98657666]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:75b8978c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:26a36c10]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:0aace77d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:81f95f63]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:1045a877]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:173d3658]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:4b8ba446]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:923ba5a6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0938 | prephenate dehydratase\n","truncated":false}}
%---
%[output:18aca330]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR204W and YGR255C and YLR201C and YML110C and YOL096C and YOR125C\n","truncated":false}}
%---
%[output:3b41965c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:0fd47cf7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:9cc3db8d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:56f30c1a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:37956aad]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:099ccb1e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:439c47fc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:9435a3e3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:979e40b5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:02ec9989]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR051W\n","truncated":false}}
%---
%[output:76b08a33]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:2866f2cd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:9b12a899]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:9180cb67]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNL316C\n","truncated":false}}
%---
%[output:1adecb58]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:539c39d0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:7676bec5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:3baeede7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:88ac9df4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0732 | methylenetetrahydrofolate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:4eebb475]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOL059W\n","truncated":false}}
%---
%[output:6885054a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0939 | prephenate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:5ca6ce01]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:7e0714c4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:405004bc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:86bf3394]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:4e347fec]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:7b97c7eb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:64e087e6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:36834873]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:296d1501]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YML120C\n","truncated":false}}
%---
%[output:68a31698]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:91f2fa98]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:36e48e48]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:85d3febc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:21a84636]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:47f4a7b7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:52b35b6f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:62f90b2a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:56ad8462]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0165 | mitochondrial alcohol dehydrogenase\n","truncated":false}}
%---
%[output:96b1571c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR204W and YGR255C and YLR201C and YML110C and YOL096C and YOR125C\n","truncated":false}}
%---
%[output:11b3f321]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:269dbc6e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:76ad496e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:6a59fba6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:8d70a260]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:67fb38ad]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:60d8a7c5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:537a3514]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:397cea1d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:5e75bd20]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR051W\n","truncated":false}}
%---
%[output:1a27260d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:1b08affb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:2fbc05dc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:6760a959]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:4f62d1b0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:3398b722]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:4f3458c3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:8a12f862]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:67fb9f16]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:542c1d14]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOL059W\n","truncated":false}}
%---
%[output:3f94d695]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:0e28c1f4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0732 | methylenetetrahydrofolate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:4d20d4d2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:0d35b411]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0939 | prephenate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:5d516970]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:2d1f4aba]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:31901581]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:11ea2aa5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:34b551c9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:7755283a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YML120C\n","truncated":false}}
%---
%[output:59aaf96b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0187 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:25725416]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:56f62df7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:73c80221]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:6745ff59]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:126b9bc7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:206fc9ea]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:293723b6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:0a1f48dd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0938 | prephenate dehydratase\n","truncated":false}}
%---
%[output:9f9d03a6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR204W and YGR255C and YLR201C and YML110C and YOL096C and YOR125C\n","truncated":false}}
%---
%[output:84e9f59e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:7030d79f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:519a2448]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:1a2e6670]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:5bbb3fd7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:523f105e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:7afcbba0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:25114548]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2140 | fatty-acyl-CoA synthase (n-C16:0CoA)\n","truncated":false}}
%---
%[output:62270d09]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2141 | fatty-acyl-CoA synthase (n-C18:0CoA)\n","truncated":false}}
%---
%[output:14a4a490]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR051W\n","truncated":false}}
%---
%[output:37372cda]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:27427522]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:0797dd81]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:81282cb9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0938 | prephenate dehydratase\n","truncated":false}}
%---
%[output:13e203b0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:07aebe69]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0960 | pyruvate decarboxylase (acetoin-forming)\n","truncated":false}}
%---
%[output:1dde7715]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:9f7fd7f5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:99d8b6ba]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:38d3fd67]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:53e809aa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:4cf0a0a6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:92dbd16f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:39acd1bc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:86bce8bb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:236d02fe]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:70ac9825]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:85492a15]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:83aba52a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:9aaa57b7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:4440235a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:8ed88977]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:93164d40]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:95207f77]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNL316C\n","truncated":false}}
%---
%[output:2212e609]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:46b65a9f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:1a9ad448]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0171 | aldehyde dehydrogenase (2-phenylethanol, NADP)\n","truncated":false}}
%---
%[output:56bc6827]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0171 | aldehyde dehydrogenase (2-phenylethanol, NADP)\n","truncated":false}}
%---
%[output:27e13ecd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0171 | aldehyde dehydrogenase (2-phenylethanol, NADP)\n","truncated":false}}
%---
%[output:0ba3b505]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:51172f35]
%   data: {"dataType":"text","outputData":{"text":"  r_Cg = r_0171 | aldehyde dehydrogenase (2-phenylethanol, NADP)\n","truncated":false}}
%---
%[output:903da677]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:11de238b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:9cf8c192]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:5d00d2b8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:4248617a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:15a6aabb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4652 | aldehyde dehydrogenase (1-propanol, NAD)\n","truncated":false}}
%---
%[output:922b7bb1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4661 | aldehyde dehydrogenase (methionol, NAD)\n","truncated":false}}
%---
%[output:3c3a274d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4670 | aldehyde dehydrogenase (tyrosol, NAD)\n","truncated":false}}
%---
%[output:59e81f23]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:5fb155a7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:7a06bd5e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:3f557202]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:18e5e103]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:9ff957b8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:9b798597]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:5b43d3ae]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:69fac5ad]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:2f92820f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:0eb7d931]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:62384f37]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:9a71e4b1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:56134f30]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:38c0a4ab]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:12e63ea8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:49e85832]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4652 | aldehyde dehydrogenase (1-propanol, NAD)\n","truncated":false}}
%---
%[output:85204701]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:06567046]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4661 | aldehyde dehydrogenase (methionol, NAD)\n","truncated":false}}
%---
%[output:669248ef]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:29e839f9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:68f0050c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4670 | aldehyde dehydrogenase (tyrosol, NAD)\n","truncated":false}}
%---
%[output:708f4385]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:3788cecf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:4f079df8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:911c203e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:588a2950]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:1e0435d9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:2346fb6d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:2b06f57a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:40ba01b0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:3ed76987]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:0211e6ee]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:08278dd3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:7447593e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:7ab3bdd6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:42569ef8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:85307bdd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4652 | aldehyde dehydrogenase (1-propanol, NAD)\n","truncated":false}}
%---
%[output:30d0fe70]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4661 | aldehyde dehydrogenase (methionol, NAD)\n","truncated":false}}
%---
%[output:2a661c52]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:9bea26c5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YNL316C\n","truncated":false}}
%---
%[output:5cf30e1a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4670 | aldehyde dehydrogenase (tyrosol, NAD)\n","truncated":false}}
%---
%[output:900fb9ea]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:79e43206]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:7824e02a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:3d6625eb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:7052205d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:921fadfa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:4d398825]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:9787e9dc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:37db6da6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YOL059W\n","truncated":false}}
%---
%[output:471fda00]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:63f60293]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:1892b727]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:50c41805]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:9027eeaf]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:41f9ecc5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:684f19e7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:0827b0c0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4652 | aldehyde dehydrogenase (1-propanol, NAD)\n","truncated":false}}
%---
%[output:85ce146a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4661 | aldehyde dehydrogenase (methionol, NAD)\n","truncated":false}}
%---
%[output:9815095a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YML120C\n","truncated":false}}
%---
%[output:242e70e0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4670 | aldehyde dehydrogenase (tyrosol, NAD)\n","truncated":false}}
%---
%[output:3d5b47b9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:7007eef9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:5a9fb54d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:26a401c5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:67352dc9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:4098276c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0186 | aldehyde dehydrogenase (tryptophol, NAD)\n","truncated":false}}
%---
%[output:571dfb46]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:6316c82d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:302d1ff7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDR204W and YGR255C and YLR201C and YML110C and YOL096C and YOR125C\n","truncated":false}}
%---
%[output:0db25b35]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:8d02eb64]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:0d0ed8ba]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:419fb1df]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:473a8cb6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:94582730]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2115 | alcohol dehydrogenase, (acetaldehyde to ethanol)\n","truncated":false}}
%---
%[output:0680536f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:09a7ebfc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4652 | aldehyde dehydrogenase (1-propanol, NAD)\n","truncated":false}}
%---
%[output:0fb20694]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4661 | aldehyde dehydrogenase (methionol, NAD)\n","truncated":false}}
%---
%[output:6727688d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR051W\n","truncated":false}}
%---
%[output:39d9e969]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:9544a24c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4670 | aldehyde dehydrogenase (tyrosol, NAD)\n","truncated":false}}
%---
%[output:71436b10]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0732 | methylenetetrahydrofolate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:687293f5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNL316C\n","truncated":false}}
%---
%[output:88b7751a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:8b769c98]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0939 | prephenate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:351f8354]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:18cb1518]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:31a38038]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0158 | alcohol acetyltransferase (2-methylbutanol)\n","truncated":false}}
%---
%[output:208912a7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:53a6697f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0158 | alcohol acetyltransferase (2-methylbutanol)\n","truncated":false}}
%---
%[output:6f9eddaa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0158 | alcohol acetyltransferase (2-methylbutanol)\n","truncated":false}}
%---
%[output:63246bef]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:8bd430e5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0158 | alcohol acetyltransferase (2-methylbutanol)\n","truncated":false}}
%---
%[output:92bf71cd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:183deec7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0732 | methylenetetrahydrofolate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:22e83f2d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:225d55c1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0939 | prephenate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:10d9813f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:813b4909]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:1c0eda00]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:8e76e35c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0162 | alcohol acetyltransferase (phenylethanol alcohol)\n","truncated":false}}
%---
%[output:0bb19a05]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0162 | alcohol acetyltransferase (phenylethanol alcohol)\n","truncated":false}}
%---
%[output:323066dc]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0162 | alcohol acetyltransferase (phenylethanol alcohol)\n","truncated":false}}
%---
%[output:392ecaa3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:79a40a34]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:89e8d927]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0162 | alcohol acetyltransferase (phenylethanol alcohol)\n","truncated":false}}
%---
%[output:5abb16f8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0732 | methylenetetrahydrofolate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:165fd7a8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:4269db5d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:27cf1add]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0939 | prephenate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:758d7805]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:6ba9b090]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:214a97e9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4717 | benzyl-acetate esterase, c\n","truncated":false}}
%---
%[output:62cbc513]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4717 | benzyl-acetate esterase, c\n","truncated":false}}
%---
%[output:976ca751]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_2034 | pyruvate transport\n","truncated":false}}
%---
%[output:4059712f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:25713139]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:115fb835]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4716 | ethyl-pyruvate esterase, c\n","truncated":false}}
%---
%[output:44484d79]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:78530366]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:6bb7d67f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:2ba914d8]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:41917a46]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:42be268f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:1b864846]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:58d325f3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:69191571]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:67d226f5]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0938 | prephenate dehydratase\n","truncated":false}}
%---
%[output:4f53d0ca]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR204W\n","truncated":false}}
%---
%[output:1e01c3d7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:8f05db7c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0492 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:83c1071b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0773 | NADH:ubiquinone oxidoreductase\n","truncated":false}}
%---
%[output:414617aa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0963 | quinone oxidoreductase\n","truncated":false}}
%---
%[output:93a7c866]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4264 | succinate:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:87734737]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0492 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:47a366aa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0773 | NADH:ubiquinone oxidoreductase\n","truncated":false}}
%---
%[output:937360cd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0963 | quinone oxidoreductase\n","truncated":false}}
%---
%[output:97142c2e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:7c9fbf3d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:85d98ea6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4264 | succinate:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:59849214]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0492 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:9b5dab2c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:2a346737]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0773 | NADH:ubiquinone oxidoreductase\n","truncated":false}}
%---
%[output:41279d4f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0963 | quinone oxidoreductase\n","truncated":false}}
%---
%[output:2bf051b1]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4264 | succinate:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:836ae16e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0492 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:2ae17def]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0773 | NADH:ubiquinone oxidoreductase\n","truncated":false}}
%---
%[output:9bcedccb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0963 | quinone oxidoreductase\n","truncated":false}}
%---
%[output:00391f9b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR166C\n","truncated":false}}
%---
%[output:28e55fdc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:83d6b4aa]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4264 | succinate:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:7e622bad]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0492 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:1805106c]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0773 | NADH:ubiquinone oxidoreductase\n","truncated":false}}
%---
%[output:34f3359e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:200860a3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0963 | quinone oxidoreductase\n","truncated":false}}
%---
%[output:50c5abc0]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4264 | succinate:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:1b25e15a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:62d63bee]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:0efe9db2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:8d713a11]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:55cc439b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:315c5757]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:79644e85]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:6f2a8f53]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:0356a9cb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0080 | 5,10-methylenetetrahydrofolate reductase (NADPH)\n","truncated":false}}
%---
%[output:71b32469]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:92abf472]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0938 | prephenate dehydratase\n","truncated":false}}
%---
%[output:4224c0c4]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0492 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:9e1c051a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:227d68d9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0773 | NADH:ubiquinone oxidoreductase\n","truncated":false}}
%---
%[output:62bfdabe]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0963 | quinone oxidoreductase\n","truncated":false}}
%---
%[output:9c0f3a65]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:115ae1dd]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4264 | succinate:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:652d079f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0160 | alcohol acetyltransferase (isoamyl alcohol)\n","truncated":false}}
%---
%[output:43ba1976]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0160 | alcohol acetyltransferase (isoamyl alcohol)\n","truncated":false}}
%---
%[output:936fd1e6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0160 | alcohol acetyltransferase (isoamyl alcohol)\n","truncated":false}}
%---
%[output:1761c807]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0160 | alcohol acetyltransferase (isoamyl alcohol)\n","truncated":false}}
%---
%[output:02e04187]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:331747a5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:8013802b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0732 | methylenetetrahydrofolate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:3483b6eb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:65011042]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0939 | prephenate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:4d6395c4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:3f93902f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:5a54fc63]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:46ffef91]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0161 | alcohol acetyltransferase (isobutyl alcohol)\n","truncated":false}}
%---
%[output:586c7b34]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0161 | alcohol acetyltransferase (isobutyl alcohol)\n","truncated":false}}
%---
%[output:95b5635e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0161 | alcohol acetyltransferase (isobutyl alcohol)\n","truncated":false}}
%---
%[output:6c2f0fda]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:098ee27a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0161 | alcohol acetyltransferase (isobutyl alcohol)\n","truncated":false}}
%---
%[output:3284150b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:23e1e94b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0732 | methylenetetrahydrofolate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:88baeed7]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0177 | aldehyde dehydrogenase (indole-3-acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:0b38dcf3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:8af7ec3d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0939 | prephenate dehydrogenase (NADP)\n","truncated":false}}
%---
%[output:0c89f0a3]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:2188a004]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0173 | aldehyde dehydrogenase (acetaldehyde, NADP)\n","truncated":false}}
%---
%[output:4a235452]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:49d92f23]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:42e876d9]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:02c59696]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:3cdd45bb]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:9c2f76f5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:60a1cf34]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:5ce40190]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:3401325a]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:965aff62]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0441 | FMN reductase\n","truncated":false}}
%---
%[output:0a1b6ec6]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0472 | glutamate synthase (NADH2)\n","truncated":false}}
%---
%[output:40313976]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR177C or YOR377W\n","truncated":false}}
%---
%[output:0263f737]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0491 | glycerol-3-phosphate dehydrogenase (NAD)\n","truncated":false}}
%---
%[output:56be215f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0546 | homoserine dehydrogenase (NADH)\n","truncated":false}}
%---
%[output:2069d09e]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0770 | NADH dehydrogenase, cytosolic\/mitochondrial\n","truncated":false}}
%---
%[output:0cc76f3f]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_4153 | (R)-Acetoin:NAD+ oxidoreductase\n","truncated":false}}
%---
%[output:3e57e121]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0166 | aldehyde dehydrogenase (2-methylbutanol, NAD)\n","truncated":false}}
%---
%[output:828717eb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:28db91bd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:5769b848]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0169 | aldehyde dehydrogenase (2-phenylethanol, NAD)\n","truncated":false}}
%---
%[output:6287271d]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0179 | aldehyde dehydrogenase (isoamyl alcohol, NAD)\n","truncated":false}}
%---
%[output:1f313628]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YGR204W\n","truncated":false}}
%---
%[output:821ae56b]
%   data: {"dataType":"text","outputData":{"text":"  r_Cr = r_0182 | aldehyde dehydrogenase (isobutyl alcohol, NAD)\n","truncated":false}}
%---
%[output:8c2ff8ac]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:70f8373b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:80af16c7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:2700a183]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:1fd5513c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:7c5feb85]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:044435f6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR204W\n","truncated":false}}
%---
%[output:8af1d860]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:758d4ff5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:60d4e98d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:8385661f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:38d9b9e4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR303C\n","truncated":false}}
%---
%[output:26506441]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR368W or YHR104W\n","truncated":false}}
%---
%[output:7cb4f358]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR368W or YHR104W\n","truncated":false}}
%---
%[output:40eae86a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR368W or YHR104W\n","truncated":false}}
%---
%[output:6904c29f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR368W or YHR104W\n","truncated":false}}
%---
%[output:597c8b53]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YDR368W or YHR104W\n","truncated":false}}
%---
%[output:2ba32851]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR166C\n","truncated":false}}
%---
%[output:2d710021]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR166C\n","truncated":false}}
%---
%[output:00e19fb8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:21bcf00a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:87148d16]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YDR368W or YMR318C\n","truncated":false}}
%---
%[output:62a101f5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YDR368W or YMR318C\n","truncated":false}}
%---
%[output:2568e299]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YDR368W or YMR318C\n","truncated":false}}
%---
%[output:38390145]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YDR368W or YMR318C\n","truncated":false}}
%---
%[output:19c0028f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YDR368W or YMR318C\n","truncated":false}}
%---
%[output:244575a2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YDR368W or YMR318C\n","truncated":false}}
%---
%[output:2d8dd0ab]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YDR368W or YMR318C\n","truncated":false}}
%---
%[output:4d706f09]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:6f9f0bef]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:0b18302a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:532fd258]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:9193f3ab]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:1a572571]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:095f138a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:3427801a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:4d958338]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:315435cc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:21843349]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YPL061W\n","truncated":false}}
%---
%[output:0273481a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:9c40ac4d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:74149295]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:4087f66c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:7c0eeb09]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:9c3c6b35]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:9bfa801e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:099b943f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:8c6512cd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:8b3f9c1b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:3c9b2e6e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:02eb43e1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:710ad041]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:11a18b72]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:92b7cb2d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:6298fdc0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:20c9297c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:44a2937f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:4b2e22b1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:58a76f72]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:4554cfde]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:369dbd47]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:7d67123b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:4ce6d6a3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:565ab4b3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:952a3249]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:8e8aff58]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER073W or YOR374W\n","truncated":false}}
%---
%[output:4732bcc2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER073W or YOR374W\n","truncated":false}}
%---
%[output:67735e6c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER073W or YOR374W\n","truncated":false}}
%---
%[output:3ca11807]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:67786db2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER073W or YOR374W\n","truncated":false}}
%---
%[output:48d7db13]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:41324801]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:519172cb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:28dd2a23]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:0cae8718]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:8c84c6db]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:237afc80]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:920b3d10]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:91ba1cc3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:57608b69]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YOR374W\n","truncated":false}}
%---
%[output:9f4d1525]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:8c7a09b4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:61c0290c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:50e02fcc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:14c9fee3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:5629abae]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:8a5af7da]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:1853aefb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:85466d34]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:13c89cfd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YPL061W\n","truncated":false}}
%---
%[output:410b3828]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER073W or YOR374W\n","truncated":false}}
%---
%[output:52c83729]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER073W or YOR374W\n","truncated":false}}
%---
%[output:9b1ccd13]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER073W or YOR374W\n","truncated":false}}
%---
%[output:0f10e896]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YER073W or YOR374W\n","truncated":false}}
%---
%[output:1970c947]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:3cc911ba]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:03d5607f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7d35eee1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:96019eb4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:1e050781]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:1d78db85]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:2c0e22f2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:80329be6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:0ec1d2e7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:84ddc6e4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:38c0b7bc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:0e0e5c54]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:184af97e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:067abfa9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YLR011W\n","truncated":false}}
%---
%[output:8f276e1e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:711a9150]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:3004c46f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:27f1071c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:13d8528b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:581d631f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:860c0295]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:3ae3d33e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:7661a4a2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:79f68bb1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL171C\n","truncated":false}}
%---
%[output:80f809e8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:972ee56c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:7806015a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:0376d8b4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:834f0ab9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:28592058]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:0d18da98]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:09d3bc8d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:46eb5f70]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:3c15bbbe]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL022W or YOL059W\n","truncated":false}}
%---
%[output:0f3d8361]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:424345d9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:169885ee]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:4882c043]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:75f08f3f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:422fc5c0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:85b96464]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cp) = YMR169C or YMR170C\n","truncated":false}}
%---
%[output:25dc5c5c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:8047e4fd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:3ef7aae4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YJR139C\n","truncated":false}}
%---
%[output:991d30e2]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:826688b5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:7ad6af3d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:93f83329]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:3de47713]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:8964840f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:755bc37c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:01c09aa1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:6a4ef797]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:4a02d1a9]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YDL085W or YMR145C\n","truncated":false}}
%---
%[output:03d278bd]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:26c1a1c3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7e559e7a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:81c21872]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:5bfdc40e]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:4c2a9173]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR204W\n","truncated":false}}
%---
%[output:363a9df5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR166C\n","truncated":false}}
%---
%[output:9f1c778f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:943b0097]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:233869b8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YAL061W\n","truncated":false}}
%---
%[output:1e74e0a1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:6b3f8dbc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGL256W or YMR083W\n","truncated":false}}
%---
%[output:058e3ec4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:1a6e9eb6]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:0f7c2bbf]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:79763e59]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:2c51715a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:9c2549da]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:4c603382]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNL316C\n","truncated":false}}
%---
%[output:39f6d939]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:7531af15]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:91228fd5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:8d55f5b5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:4f4450c5]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:4c1a6b53]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:3bf895e4]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YKL182W and YPL231W\n","truncated":false}}
%---
%[output:66ec09dc]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = (YGL125W and YPL023C) or YGL125W\n","truncated":false}}
%---
%[output:460961bb]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:447ceef3]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YNL316C\n","truncated":false}}
%---
%[output:5e78be45]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:782071b1]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YGR087C or YLR044C or YLR134W\n","truncated":false}}
%---
%[output:134be620]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YOL086C\n","truncated":false}}
%---
%[output:4ea9f130]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:854c1077]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:130a281f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:3fe53062]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:44ea66f7]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:8aa17401]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:033e397f]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:4c27308d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cr) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:9467d873]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:6b70f22b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:71b64e80]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:02b2cb6a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:22516c38]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YBR145W or YDL168W or YOL086C\n","truncated":false}}
%---
%[output:42ed2005]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:4efd6008]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:002fd19b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:9fd412ad]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Cg) = YCR105W or YMR318C\n","truncated":false}}
%---
%[output:8f0ab58b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: For increased performance, remaining outputs are not shown. Consider reducing the number of outputs."}}
%---
