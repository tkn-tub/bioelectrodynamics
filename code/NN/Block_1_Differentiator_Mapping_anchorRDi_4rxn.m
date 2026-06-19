%[text] # Differentiator\_Mapping\_anchorRDi\_4rxn.m
%[text] %\[text\] 4-reaction differentiator topology (with r\\\_Dp re-included): r\\\_DF: m\\\_in --\\\> m\\\_Df (fast tracker production) r\\\_DS: m\\\_in --\\\> m\\\_Ds (slow tracker production) r\\\_Dp: m\\\_Df --\\\> m\\\_d (m\\\_d production from m\\\_Df) r\\\_Di: m\\\_Df + m\\\_Ds + m\\\_d --\\\> sink (3-substrate mass-action sink/integrator)
%[text] Anchored on r\_Di: enumerate every bi-substrate enzyme-mediated reaction. For each substrate pair (s1, s2) take s1 = m\_d candidate, s2 = m\_Ds candidate. Then walk back two hops on the m\_d side (r\_Dp -\> m\_Df -\> r\_DF -\> m\_in), and one hop on the m\_Ds side (r\_DS -\> m\_in), and require the m\_in's to coincide.
% clc; clear; close all;
% 
% % NOTE: model is assumed already loaded (S, model.mets, model.rxns, etc.).
% % If not, uncomment the load block below:
% 
% addingPathParentFolderByName('code');
% initCobraToolbox(false);
% matFileName = 'iMM904.mat';
% cd(fileparts(which(matFileName)));
% model = readCbModel(matFileName);
% S = model.S;
% [nMets, nRxns] = size(S);
%%
%[text] #### Solutions with anchor on reaction r\_Di
results_diff = struct('m_in',{},'m_Df',{},'m_d',{},'m_Ds',{}, ...
                 'r_DF',{},'r_DS',{},'r_Dp',{},'r_Di',{});
nHits = 0;

% Diagnostic counters
nDi_evaluated = 0;
nDi_with_md_singleProd = 0;
nDi_with_rDp_having_singleProdSubstrate = 0;
nDi_full_match = 0;

for r_Di = 1:nRxns
    if REQUIRE_GPR && ~hasGPR(r_Di), continue; end

    subs = find(S(:,r_Di) < 0);
    subs = setdiff(subs, find(isBadMet));
    if numel(subs) < 2, continue; end   % r_Di needs m_d AND m_Ds, both critical
    nDi_evaluated = nDi_evaluated + 1;

    % Try every ordered pair (m_d candidate, m_Ds candidate) among r_Di's substrates
    for i = 1:numel(subs)
        s_md = subs(i);
        if nProducers(s_md) < 1 || nProducers(s_md) > MAX_PRODUCERS_MD, continue; end

        % r_Dp candidates: producers of m_d
        prod_md = find(S(s_md,:) > 0);
        if REQUIRE_GPR, prod_md = prod_md(hasGPR(prod_md)); end
        if isempty(prod_md), continue; end

        % For each r_Dp, look at its substrates -> m_Df candidate(s).
        % m_Df is NOT required to also be a substrate of r_Di -- it only
        % needs to feed r_Dp, which in turn produces m_d.
        for u = 1:numel(prod_md)
            r_Dp = prod_md(u);
            if r_Dp == r_Di, continue; end

            rDp_subs = find(S(:,r_Dp) < 0);
            rDp_subs = setdiff(rDp_subs, find(isBadMet));
            for u2 = 1:numel(rDp_subs)
                s_mDf = rDp_subs(u2);
                if s_mDf == s_md, continue; end
                if nProducers(s_mDf) < 1 || nProducers(s_mDf) > MAX_PRODUCERS_MDF, continue; end

                % r_DF candidates: producers of m_Df
                prod_mDf = find(S(s_mDf,:) > 0);
                if REQUIRE_GPR, prod_mDf = prod_mDf(hasGPR(prod_mDf)); end
                if isempty(prod_mDf), continue; end

                nDi_with_rDp_having_singleProdSubstrate = ...
                    nDi_with_rDp_having_singleProdSubstrate + 1;

                % Now the m_Ds side: r_Di's OTHER substrate
                for j = 1:numel(subs)
                    if j == i, continue; end
                    s_mDs = subs(j);
                    if s_mDs == s_mDf, continue; end
                    if nProducers(s_mDs) < 1 || nProducers(s_mDs) > MAX_PRODUCERS_MDS, continue; end

                    prod_mDs = find(S(s_mDs,:) > 0);
                    if REQUIRE_GPR, prod_mDs = prod_mDs(hasGPR(prod_mDs)); end
                    if isempty(prod_mDs), continue; end

                    % Now require shared input m_in between r_DF and r_DS
                    for u3 = 1:numel(prod_mDf)
                        r_DF = prod_mDf(u3);
                        if ismember(r_DF, [r_Di, r_Dp]), continue; end
                        for v = 1:numel(prod_mDs)
                            r_DS = prod_mDs(v);
                            if ismember(r_DS, [r_Di, r_Dp, r_DF]), continue; end

                            shared = intersect(substrateSet{r_DF}, substrateSet{r_DS});
                            shared = setdiff(shared, find(isBadMet));
                            shared = setdiff(shared, [s_md; s_mDf; s_mDs]);
                            if isempty(shared), continue; end

                            for w = 1:numel(shared)
                                nHits = nHits + 1;
                                results_diff(nHits).m_in = shared(w);
                                results_diff(nHits).m_Df = s_mDf;
                                results_diff(nHits).m_d  = s_md;
                                results_diff(nHits).m_Ds = s_mDs;
                                results_diff(nHits).r_DF = r_DF;
                                results_diff(nHits).r_DS = r_DS;
                                results_diff(nHits).r_Dp = r_Dp;
                                results_diff(nHits).r_Di = r_Di;
                                if nHits >= MAX_HITS, break; end
                            end
                            if nHits >= MAX_HITS, break; end
                        end
                        if nHits >= MAX_HITS, break; end
                    end
                    if nHits >= MAX_HITS, break; end
                end
                if nHits >= MAX_HITS, break; end
            end
            if nHits >= MAX_HITS, break; end
        end
        if nHits >= MAX_HITS, break; end
    end
    if nHits >= MAX_HITS, break; end
end

% fprintf('\n--- Diagnostic ---\n');
% fprintf('  Bi-substrate enzyme rxns evaluated: %d\n', nDi_evaluated);
% fprintf('  ...with viable r_Dp+m_Df chain on one substrate: %d\n', ...
%         nDi_with_rDp_having_singleProdSubstrate);
fprintf('Total candidates (4-reaction topology): %d\n', nHits); %[output:7b520689]
%[text] ## Print first hits
if nHits > 0 %[output:group:8c2592f9]
    Motifs_Diff_4rxn = struct2table(results_diff);
    fprintf('\n=== First candidates (4-reaction differentiator) ===\n'); %[output:630eb9fc]
    for i = 1:min(5, nHits)
        rec = results_diff(i);
        fprintf('\n--- Candidate %d ---\n', i); %[output:44790522] %[output:4d37853d] %[output:2e64ed35] %[output:70cd6412] %[output:4db5bbcf]
        fprintf('  m_in  = %s | %s\n', model.mets{rec.m_in},  model.metNames{rec.m_in}); %[output:08b0cac7] %[output:62726b0d] %[output:34d17fd0] %[output:69463aef] %[output:3f2e141d]
        fprintf('  m_Df  = %s | %s\n', model.mets{rec.m_Df},  model.metNames{rec.m_Df}); %[output:5f776f03] %[output:10790759] %[output:0f75cdbd] %[output:478c6ad5] %[output:4ea7b9bf]
        fprintf('  m_d   = %s | %s\n', model.mets{rec.m_d},   model.metNames{rec.m_d}); %[output:249e19ef] %[output:01609b3e] %[output:840b529a] %[output:6b5e8e75] %[output:6c34a522]
        fprintf('  m_Ds  = %s | %s\n', model.mets{rec.m_Ds},  model.metNames{rec.m_Ds}); %[output:66139710] %[output:693267a2] %[output:243c378b] %[output:2259f7e0] %[output:48c8ed87]
        fprintf('  r_DF  = %s | %s\n', model.rxns{rec.r_DF},  model.rxnNames{rec.r_DF}); %[output:6b45a301] %[output:47240e3c] %[output:45c6f014] %[output:940941b1] %[output:2d316cd9]
        fprintf('  r_DS  = %s | %s\n', model.rxns{rec.r_DS},  model.rxnNames{rec.r_DS}); %[output:2c3c2e7f] %[output:216fd483] %[output:3dc418fa] %[output:40f0eb54] %[output:4b72a925]
        fprintf('  r_Dp  = %s | %s\n', model.rxns{rec.r_Dp},  model.rxnNames{rec.r_Dp}); %[output:7e64c393] %[output:7e099ae3] %[output:3b31ea86] %[output:98030daa] %[output:32f21137]
        fprintf('  r_Di  = %s | %s\n', model.rxns{rec.r_Di},  model.rxnNames{rec.r_Di}); %[output:0fc4f7df] %[output:7fc8cbef] %[output:99af4c72] %[output:2d1b9b74] %[output:9cbd12a6]
        fprintf('  GPR(r_DF) = %s\n', model.grRules{rec.r_DF}); %[output:3f44ba1d] %[output:1f235d00] %[output:50089c91] %[output:9cf917ee] %[output:96e9181b]
        fprintf('  GPR(r_DS) = %s\n', model.grRules{rec.r_DS}); %[output:99b5e676] %[output:9a066707] %[output:4dd33002] %[output:271f4c6a] %[output:45b84307]
        fprintf('  GPR(r_Dp) = %s\n', model.grRules{rec.r_Dp}); %[output:19c34f28] %[output:6ff71f35] %[output:9b620f01] %[output:84a7247b] %[output:6ac1487d]
        fprintf('  GPR(r_Di) = %s\n', model.grRules{rec.r_Di}); %[output:9f36ecd0] %[output:87dd8841] %[output:781ee5d8] %[output:4b47784c] %[output:795ee731]
    end
    save('Block_1_differentiator_candidates_anchorRDi_4rxn.mat', 'Motifs_Diff_4rxn', 'results_diff');
end %[output:group:8c2592f9]



% % Performs four steps:
% %   STEP 1: Print all candidates with full detail (reactions, GPRs, EC, compartment).
% %   STEP 2: Compartment table per candidate; flag same-compartment vs cross-compartment.
% %   STEP 3: Gene-disjointness check across the four reactions.
% %   STEP 4: Deduplicate by (r_DF, r_DS, r_Dp, r_Di) tuple.
% %
% % Output: a ranked table 'Shortlist' with the survivors, sorted by quality.
% %
% % ASSUMES: model and 'results' are already in the workspace from the search script.
% % If not, uncomment the load block below.
% 
% % load('differentiator_candidates_anchorRDi_4rxn.mat');   % loads Motifs_Diff_4rxn, results
% 
% if ~exist('results','var') || isempty(results)
%     error('No ''results'' struct in workspace. Run Differentiator_Mapping_anchorRDi_4rxn.m first.');
% end
% 
% nC = numel(results);
% fprintf('Inspecting %d differentiator candidates.\n\n', nC);
% 
% %% =========================================================
% % STEP 1: Print all candidates in detail
% % =========================================================
% fprintf('==========================================================\n');
% fprintf(' STEP 1: Detailed listing of all %d candidates\n', nC);
% fprintf('==========================================================\n');
% 
% for k = 1:nC
%     rec = results(k);
%     fprintf('\n----------------------------------------\n');
%     fprintf(' Candidate %d / %d\n', k, nC);
%     fprintf('----------------------------------------\n');
% 
%     % Metabolites
%     metsTbl = {'m_in', rec.m_in;
%                'm_Df', rec.m_Df;
%                'm_d',  rec.m_d;
%                'm_Ds', rec.m_Ds};
%     for i = 1:size(metsTbl,1)
%         m = metsTbl{i,2};
%         fprintf('  %s = %-15s | %-40s | formula: %s\n', ...
%             metsTbl{i,1}, ...
%             string(model.mets{m}), ...
%             string(model.metNames{m}), ...
%             string(model.metFormulas{m}));
%     end
%     fprintf('\n');
% 
%     % Reactions
%     rxnsTbl = {'r_DF', rec.r_DF;
%                'r_DS', rec.r_DS;
%                'r_Dp', rec.r_Dp;
%                'r_Di', rec.r_Di};
%     for i = 1:size(rxnsTbl,1)
%         r = rxnsTbl{i,2};
%         fprintf('  %s = %s   |   %s\n', rxnsTbl{i,1}, ...
%             string(model.rxns{r}), string(model.rxnNames{r}));
%         try
%             fprintf('     equation: %s\n', reactionToString_chm(model, r, false));
%         catch
%             % fall back: build the equation from S
%             substrates = find(model.S(:,r) < 0);
%             products   = find(model.S(:,r) > 0);
%             sStr = strjoin(string(model.mets(substrates)), ' + ');
%             pStr = strjoin(string(model.mets(products)),   ' + ');
%             fprintf('     equation: %s --> %s\n', sStr, pStr);
%         end
%         gpr = string(model.grRules{r});
%         if strlength(gpr) == 0, gpr = "(none)"; end
%         fprintf('     GPR: %s\n', gpr);
%         if isfield(model,'eccodes') && numel(model.eccodes) >= r && ~isempty(model.eccodes{r})
%             fprintf('     EC: %s\n', string(model.eccodes{r}));
%         end
%         try
%             loc = getReactionLocationString(model, r);
%             fprintf('     Location: %s\n', loc);
%         catch
%             % helper not on path; fall back to extracting from metabolite IDs
%             mets_r = find(model.S(:,r) ~= 0);
%             comps  = compartmentsOfMets(model, mets_r);
%             fprintf('     Location (inferred): %s\n', strjoin(comps, ','));
%         end
%     end
% end
% 
% %% =========================================================
% % STEP 2: Compartment table per candidate
% % =========================================================
% fprintf('\n\n==========================================================\n');
% fprintf(' STEP 2: Compartment analysis\n');
% fprintf('==========================================================\n');
% 
% compStrings = strings(nC, 4);     % rows: candidates; cols: r_DF, r_DS, r_Dp, r_Di
% sameCompartment = false(nC,1);
% nDistinctComps = zeros(nC,1);
% 
% for k = 1:nC
%     rec = results(k);
%     rxns = [rec.r_DF, rec.r_DS, rec.r_Dp, rec.r_Di];
%     for i = 1:4
%         try
%             cs = string(getReactionLocationString(model, rxns(i)));
%         catch
%             mets_r = find(model.S(:,rxns(i)) ~= 0);
%             comps  = compartmentsOfMets(model, mets_r);
%             cs = strjoin(comps, ',');
%         end
%         compStrings(k,i) = cs;
%     end
%     uc = unique(compStrings(k,:));
%     % Filter out empty strings
%     uc = uc(strlength(uc) > 0);
%     nDistinctComps(k) = numel(uc);
%     sameCompartment(k) = (nDistinctComps(k) == 1);
% end
% 
% fprintf('\n%4s | %-20s | %-20s | %-20s | %-20s | same?\n', ...
%         'idx', 'r_DF compartment', 'r_DS compartment', 'r_Dp compartment', 'r_Di compartment');
% fprintf('%s\n', repmat('-',1,120));
% for k = 1:nC
%     fprintf('%4d | %-20s | %-20s | %-20s | %-20s |  %s\n', ...
%         k, compStrings(k,1), compStrings(k,2), compStrings(k,3), compStrings(k,4), ...
%         ternary(sameCompartment(k),'YES','--'));
% end
% 
% fprintf('\nCompartment-consistent candidates: %d / %d\n', sum(sameCompartment), nC);
% 
% %% =========================================================
% % STEP 3: Gene-disjointness check
% % =========================================================
% fprintf('\n\n==========================================================\n');
% fprintf(' STEP 3: Gene-disjointness check (across r_DF, r_DS, r_Dp, r_Di)\n');
% fprintf('==========================================================\n');
% 
% geneDisjoint = false(nC,1);
% geneOverlapInfo = strings(nC,1);
% 
% for k = 1:nC
%     rec = results(k);
%     rxns = [rec.r_DF, rec.r_DS, rec.r_Dp, rec.r_Di];
%     geneSets = cell(4,1);
%     for i = 1:4
%         gpr = string(model.grRules{rxns(i)});
%         if strlength(gpr) == 0
%             geneSets{i} = strings(0,1);
%         else
%             toks = regexp(gpr, '[A-Za-z0-9_:\-\.]+', 'match');
%             toks = toks(~ismember(lower(toks), {'and','or'}));
%             geneSets{i} = unique(string(toks));
%         end
%     end
% 
%     % Pairwise check
%     isDisj = true;
%     overlapMsg = "";
%     rxnLabels = {'r_DF','r_DS','r_Dp','r_Di'};
%     for i = 1:4
%         for j = i+1:4
%             common = intersect(geneSets{i}, geneSets{j});
%             if ~isempty(common)
%                 isDisj = false;
%                 overlapMsg = overlapMsg + ...
%                     sprintf("%s<>%s share %s; ", rxnLabels{i}, rxnLabels{j}, ...
%                             strjoin(common,','));
%             end
%         end
%     end
%     geneDisjoint(k) = isDisj;
%     geneOverlapInfo(k) = overlapMsg;
% end
% 
% fprintf('\n%4s | %-10s | overlap details\n', 'idx', 'disjoint?');
% fprintf('%s\n', repmat('-',1,80));
% for k = 1:nC
%     if geneDisjoint(k)
%         fprintf('%4d |    YES    |\n', k);
%     else
%         fprintf('%4d |    --     | %s\n', k, geneOverlapInfo(k));
%     end
% end
% fprintf('\nGene-disjoint candidates: %d / %d\n', sum(geneDisjoint), nC);
% 
% %% =========================================================
% % STEP 4: Deduplicate by reaction tuple
% % =========================================================
% fprintf('\n\n==========================================================\n');
% fprintf(' STEP 4: Deduplicate by (r_DF, r_DS, r_Dp, r_Di) tuple\n');
% fprintf('==========================================================\n');
% 
% rxnTuples = zeros(nC,4);
% for k = 1:nC
%     rxnTuples(k,:) = [results(k).r_DF, results(k).r_DS, ...
%                       results(k).r_Dp, results(k).r_Di];
% end
% [~, uniqueIdx] = unique(rxnTuples, 'rows', 'stable');
% isUnique = false(nC,1);
% isUnique(uniqueIdx) = true;
% nUnique = sum(isUnique);
% fprintf('Unique (r_DF, r_DS, r_Dp, r_Di) tuples: %d / %d\n', nUnique, nC);
% fprintf('(Duplicates differ only in m_in choice or other m_Df/m_Ds choices.)\n');
% 
% %% =========================================================
% % Build the shortlist: surviving on all 4 filters
% % =========================================================
% fprintf('\n\n==========================================================\n');
% fprintf(' Final shortlist\n');
% fprintf('==========================================================\n');
% 
% passAll = sameCompartment & geneDisjoint & isUnique;
% nPass   = sum(passAll);
% fprintf('Candidates passing all filters (compartment + gene-disjoint + unique): %d\n', nPass);
% 
% % Also build a "soft" shortlist: just unique + gene-disjoint, allowing cross-compartment.
% % Useful if strict same-compartment yields 0.
% passSoft = isUnique & geneDisjoint;
% nSoft = sum(passSoft);
% fprintf('Soft shortlist (unique + gene-disjoint, any compartment): %d\n', nSoft);
% 
% % Build a table for the strict shortlist
% if nPass > 0
%     sortedIdx = find(passAll);
% elseif nSoft > 0
%     fprintf('\n[Strict empty -> showing soft shortlist instead]\n');
%     sortedIdx = find(passSoft);
% else
%     fprintf('\n[All candidates failed filters; showing all unique candidates]\n');
%     sortedIdx = find(isUnique);
% end
% 
% ShortlistRows = struct('idx',{},'m_in',{},'m_Df',{},'m_d',{},'m_Ds',{}, ...
%                        'r_DF',{},'r_DS',{},'r_Dp',{},'r_Di',{}, ...
%                        'compartments',{},'gene_disjoint',{});
% for kk = 1:numel(sortedIdx)
%     k = sortedIdx(kk);
%     rec = results(k);
%     ShortlistRows(kk).idx  = k;
%     ShortlistRows(kk).m_in = string(model.mets{rec.m_in});
%     ShortlistRows(kk).m_Df = string(model.mets{rec.m_Df});
%     ShortlistRows(kk).m_d  = string(model.mets{rec.m_d});
%     ShortlistRows(kk).m_Ds = string(model.mets{rec.m_Ds});
%     ShortlistRows(kk).r_DF = string(model.rxns{rec.r_DF});
%     ShortlistRows(kk).r_DS = string(model.rxns{rec.r_DS});
%     ShortlistRows(kk).r_Dp = string(model.rxns{rec.r_Dp});
%     ShortlistRows(kk).r_Di = string(model.rxns{rec.r_Di});
%     ShortlistRows(kk).compartments = strjoin(unique(compStrings(k,:)), ',');
%     ShortlistRows(kk).gene_disjoint = geneDisjoint(k);
% end
% Shortlist = struct2table(ShortlistRows);
% disp(Shortlist);
% 
% % Save
% save('differentiator_shortlist.mat', 'Shortlist', 'sameCompartment', ...
%      'geneDisjoint', 'isUnique', 'compStrings', 'geneOverlapInfo');
% 
% fprintf('\nSaved differentiator_shortlist.mat\n');
% 
% %% =========================================================
% % Helper functions
% % =========================================================
% 
% function out = ternary(cond, a, b)
%     if cond, out = a; else, out = b; end
% end
% 
% function comps = compartmentsOfMets(model, metIdxs)
%     comps = strings(numel(metIdxs),1);
%     for ii = 1:numel(metIdxs)
%         id = string(model.mets{metIdxs(ii)});
%         tok = regexp(id, '\[([a-z]+)\]$', 'tokens', 'once');
%         if ~isempty(tok)
%             comps(ii) = tok{1};
%         else
%             comps(ii) = "?";
%         end
%     end
%     comps = unique(comps);
% end

function addingPathParentFolderByName(targetName)
    % Start from the current directory
    currFolder = pwd;
    found = false;
    
    % Continue searching until you reach the root folder
    while true
        % Get the parent folder
        [parentFolder, currentName] = fileparts(currFolder);
        
        % Check if the current folder's name is the target
        if strcmpi(currentName, targetName)
            found = true;
            break;
        end
        
        % If we've reached the root or no change, exit the loop
        if isempty(parentFolder) || strcmp(currFolder, parentFolder)
            break;
        end
        
        % Move one level up
        currFolder = parentFolder;
    end

    if found
        addpath(genpath(currFolder));
        fprintf('Adding matlab path to: %s\n', currFolder);
    else
        error('Folder named "%s" not found in any parent directory.', targetName);
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:7b520689]
%   data: {"dataType":"text","outputData":{"text":"Total candidates (4-reaction topology): 1839\n","truncated":false}}
%---
%[output:630eb9fc]
%   data: {"dataType":"text","outputData":{"text":"\n=== First candidates (4-reaction differentiator) ===\n","truncated":false}}
%---
%[output:44790522]
%   data: {"dataType":"text","outputData":{"text":"\n--- Candidate 1 ---\n","truncated":false}}
%---
%[output:08b0cac7]
%   data: {"dataType":"text","outputData":{"text":"  m_in  = s_1401 | pyruvate\n","truncated":false}}
%---
%[output:5f776f03]
%   data: {"dataType":"text","outputData":{"text":"  m_Df  = s_0182 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:249e19ef]
%   data: {"dataType":"text","outputData":{"text":"  m_d   = s_0233 | 3-methyl-2-oxobutanoate\n","truncated":false}}
%---
%[output:66139710]
%   data: {"dataType":"text","outputData":{"text":"  m_Ds  = s_0376 | acetyl-CoA\n","truncated":false}}
%---
%[output:6b45a301]
%   data: {"dataType":"text","outputData":{"text":"  r_DF  = r_4755 | (R)-2-hydroxyglutarate:pyruvate oxidoreductase\n","truncated":false}}
%---
%[output:2c3c2e7f]
%   data: {"dataType":"text","outputData":{"text":"  r_DS  = r_0961 | pyruvate dehydrogenase\n","truncated":false}}
%---
%[output:7e64c393]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp  = r_1088 | valine transaminase, mitochondiral\n","truncated":false}}
%---
%[output:0fc4f7df]
%   data: {"dataType":"text","outputData":{"text":"  r_Di  = r_0025 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:3f44ba1d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_DF) = YDL178W or YEL071W\n","truncated":false}}
%---
%[output:99b5e676]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_DS) = YBR221C and YER178W and YFL018C and YGR193C and YNL071W\n","truncated":false}}
%---
%[output:19c34f28]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Dp) = YHR208W\n","truncated":false}}
%---
%[output:9f36ecd0]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Di) = YNL104C\n","truncated":false}}
%---
%[output:4d37853d]
%   data: {"dataType":"text","outputData":{"text":"\n--- Candidate 2 ---\n","truncated":false}}
%---
%[output:2e64ed35]
%   data: {"dataType":"text","outputData":{"text":"\n--- Candidate 3 ---\n","truncated":false}}
%---
%[output:70cd6412]
%   data: {"dataType":"text","outputData":{"text":"\n--- Candidate 4 ---\n","truncated":false}}
%---
%[output:4db5bbcf]
%   data: {"dataType":"text","outputData":{"text":"\n--- Candidate 5 ---\n","truncated":false}}
%---
%[output:62726b0d]
%   data: {"dataType":"text","outputData":{"text":"  m_in  = s_0182 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:34d17fd0]
%   data: {"dataType":"text","outputData":{"text":"  m_in  = s_0182 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:69463aef]
%   data: {"dataType":"text","outputData":{"text":"  m_in  = s_0182 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:3f2e141d]
%   data: {"dataType":"text","outputData":{"text":"  m_in  = s_0182 | 2-oxoglutarate\n","truncated":false}}
%---
%[output:10790759]
%   data: {"dataType":"text","outputData":{"text":"  m_Df  = s_0532 | coenzyme A\n","truncated":false}}
%---
%[output:0f75cdbd]
%   data: {"dataType":"text","outputData":{"text":"  m_Df  = s_0532 | coenzyme A\n","truncated":false}}
%---
%[output:478c6ad5]
%   data: {"dataType":"text","outputData":{"text":"  m_Df  = s_0532 | coenzyme A\n","truncated":false}}
%---
%[output:4ea7b9bf]
%   data: {"dataType":"text","outputData":{"text":"  m_Df  = s_1401 | pyruvate\n","truncated":false}}
%---
%[output:87dd8841]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Di) = YNL104C\n","truncated":false}}
%---
%[output:6ff71f35]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Dp) = YAL054C\n","truncated":false}}
%---
%[output:9a066707]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_DS) = YHR208W\n","truncated":false}}
%---
%[output:1f235d00]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_DF) = YDL131W or YDL182W\n","truncated":false}}
%---
%[output:7fc8cbef]
%   data: {"dataType":"text","outputData":{"text":"  r_Di  = r_0025 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:7e099ae3]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp  = r_0113 | acetyl-CoA synthetase\n","truncated":false}}
%---
%[output:216fd483]
%   data: {"dataType":"text","outputData":{"text":"  r_DS  = r_1088 | valine transaminase, mitochondiral\n","truncated":false}}
%---
%[output:47240e3c]
%   data: {"dataType":"text","outputData":{"text":"  r_DF  = r_1838 | homocitrate synthase\n","truncated":false}}
%---
%[output:693267a2]
%   data: {"dataType":"text","outputData":{"text":"  m_Ds  = s_0233 | 3-methyl-2-oxobutanoate\n","truncated":false}}
%---
%[output:01609b3e]
%   data: {"dataType":"text","outputData":{"text":"  m_d   = s_0376 | acetyl-CoA\n","truncated":false}}
%---
%[output:840b529a]
%   data: {"dataType":"text","outputData":{"text":"  m_d   = s_0376 | acetyl-CoA\n","truncated":false}}
%---
%[output:6b5e8e75]
%   data: {"dataType":"text","outputData":{"text":"  m_d   = s_0376 | acetyl-CoA\n","truncated":false}}
%---
%[output:6c34a522]
%   data: {"dataType":"text","outputData":{"text":"  m_d   = s_0376 | acetyl-CoA\n","truncated":false}}
%---
%[output:781ee5d8]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Di) = YNL104C\n","truncated":false}}
%---
%[output:243c378b]
%   data: {"dataType":"text","outputData":{"text":"  m_Ds  = s_0233 | 3-methyl-2-oxobutanoate\n","truncated":false}}
%---
%[output:9b620f01]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Dp) = YML042W\n","truncated":false}}
%---
%[output:2259f7e0]
%   data: {"dataType":"text","outputData":{"text":"  m_Ds  = s_0233 | 3-methyl-2-oxobutanoate\n","truncated":false}}
%---
%[output:4dd33002]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_DS) = YHR208W\n","truncated":false}}
%---
%[output:48c8ed87]
%   data: {"dataType":"text","outputData":{"text":"  m_Ds  = s_0233 | 3-methyl-2-oxobutanoate\n","truncated":false}}
%---
%[output:50089c91]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_DF) = YDL131W or YDL182W\n","truncated":false}}
%---
%[output:45c6f014]
%   data: {"dataType":"text","outputData":{"text":"  r_DF  = r_1838 | homocitrate synthase\n","truncated":false}}
%---
%[output:99af4c72]
%   data: {"dataType":"text","outputData":{"text":"  r_Di  = r_0025 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:940941b1]
%   data: {"dataType":"text","outputData":{"text":"  r_DF  = r_1838 | homocitrate synthase\n","truncated":false}}
%---
%[output:3b31ea86]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp  = r_0254 | carnitine O-acetyltransferase\n","truncated":false}}
%---
%[output:2d316cd9]
%   data: {"dataType":"text","outputData":{"text":"  r_DF  = r_0674 | L-alanine transaminase\n","truncated":false}}
%---
%[output:3dc418fa]
%   data: {"dataType":"text","outputData":{"text":"  r_DS  = r_1088 | valine transaminase, mitochondiral\n","truncated":false}}
%---
%[output:4b47784c]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Di) = YNL104C\n","truncated":false}}
%---
%[output:40f0eb54]
%   data: {"dataType":"text","outputData":{"text":"  r_DS  = r_1088 | valine transaminase, mitochondiral\n","truncated":false}}
%---
%[output:4b72a925]
%   data: {"dataType":"text","outputData":{"text":"  r_DS  = r_1088 | valine transaminase, mitochondiral\n","truncated":false}}
%---
%[output:84a7247b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Dp) = YBR221C and YER178W and YFL018C and YGR193C and YNL071W\n","truncated":false}}
%---
%[output:98030daa]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp  = r_0961 | pyruvate dehydrogenase\n","truncated":false}}
%---
%[output:32f21137]
%   data: {"dataType":"text","outputData":{"text":"  r_Dp  = r_0961 | pyruvate dehydrogenase\n","truncated":false}}
%---
%[output:271f4c6a]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_DS) = YHR208W\n","truncated":false}}
%---
%[output:2d1b9b74]
%   data: {"dataType":"text","outputData":{"text":"  r_Di  = r_0025 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:9cbd12a6]
%   data: {"dataType":"text","outputData":{"text":"  r_Di  = r_0025 | 2-isopropylmalate synthase\n","truncated":false}}
%---
%[output:9cf917ee]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_DF) = YDL131W or YDL182W\n","truncated":false}}
%---
%[output:96e9181b]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_DF) = YLR089C\n","truncated":false}}
%---
%[output:45b84307]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_DS) = YHR208W\n","truncated":false}}
%---
%[output:6ac1487d]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Dp) = YBR221C and YER178W and YFL018C and YGR193C and YNL071W\n","truncated":false}}
%---
%[output:795ee731]
%   data: {"dataType":"text","outputData":{"text":"  GPR(r_Di) = YNL104C\n","truncated":false}}
%---
