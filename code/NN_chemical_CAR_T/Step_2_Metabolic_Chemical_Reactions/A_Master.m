% run_IFFL_master.m  ---------------------------------------------------------
% Single master, three steps:
%   STEP 1 (BRENDA text): reactions whose EC is ACTIVATED by a model metabolite  -> seeds
%   STEP 2 (S matrix):    anchor r_I; grow the two branches BACKWARD to a common
%                         input; fast arm strictly shorter; branches REACTION- and
%                         METABOLITE-disjoint                                     -> motifs
%   STEP 3 (concatenate): pair motifs so D1's output = D2's input; the two
%                         differentiators reaction- and metabolite-disjoint       -> cascade
%
% Needs on path: parse_brenda_flat.m, iffl_coverage_check.m,
% addingPathParentFolderByName.m, COBRA; model + brenda_2026_1.txt under 'code'.
% ---------------------------------------------------------------------------

clc; clearvars; close all;
addingPathParentFolderByName('code');
if exist('readCbModel','file') ~= 2, initCobraToolbox(false); end

modelFile = which('Puniya2021_TNM1055.xml');
if isempty(modelFile), modelFile = which('Puniya2021_TNM1055.mat'); end
assert(~isempty(modelFile),'Model not found under ''code''.');
brendaFile = which('brenda_2026_1.txt');
assert(~isempty(brendaFile),'brenda_2026_1.txt not found under ''code''.');
[~,~,ext] = fileparts(modelFile);
if strcmpi(ext,'.mat'), L=load(modelFile); fn=fieldnames(L); model=L.(fn{1});
else, model = readCbModel(modelFile); end
fprintf('Model %s: %d reactions, %d metabolites\n',modelFile,numel(model.rxns),numel(model.mets));

% ---- options ----
opt = struct();
opt.wantAC            = true;
opt.maxFast           = 3;
opt.maxSlow           = 5;
opt.forwardOnly       = true;
opt.disjointBranches  = true;    % fast & slow share NO reaction
opt.metDisjointBranch = true;    % NEW: fast & slow share NO interior metabolite (besides the input)
opt.cascadeDisjoint   = true;    % D1 & D2 share NO reaction
opt.metDisjointCasc   = true;    % NEW: D1 & D2 share NO metabolite besides the link
opt.topN              = 500;
% currency by metabolite-ID base (compartment stripped):
opt.currency = {'h','h2o','atp','adp','amp','pi','ppi','nad','nadh','nadp','nadph', ...
                'co2','o2','nh4','coa','fad','fadh2','gtp','gdp','udp','utp','ctp', ...
                'so4','na1','k','cl','fe2','fe3','mg2','ca2','zn2'};
% NEW: carrier families excluded by NAME pattern (catches crn + all acyl-carnitines,
% acyl-CoAs, etc., which exact IDs cannot). Remove 'coenzyme a' if too aggressive.
opt.currencyNamePatterns = {'carnitine','coenzyme a'};
% opt.ecMap = <containers.Map rxnID -> 'x.x.x.x'>;

% =============================== STEP 1 ======================================
T = iffl_coverage_check(model, brendaFile, opt);
fprintf('\nSTEP 1: %d seed reactions.\n', height(T));
assert(height(T)>0,'No seeds - nothing to anchor on.');

% =============================== STEP 2 ======================================
M = build_iffl_from_anchor(model, T, opt);
fprintf('STEP 2: %d disjoint IFFL motifs.\n', height(M));
if height(M)>0, writetable(M(:,setdiff(M.Properties.VariableNames,{'metIDs'})),'iffl_motifs.csv'); end

% =============================== STEP 3 ======================================
C = concatenate_diffs(M, model, opt);
fprintf('STEP 3: %d D1->D2 concatenations.\n', height(C));
if height(C)>0
    disp(head(C, min(20,height(C))));
    writetable(C,'iffl_cascade.csv');
    fprintf('Saved iffl_motifs.csv and iffl_cascade.csv\n');
else
    fprintf('No cascade survived the disjointness filters.\n');
end

% =============================== STEP 4 ======================================
% Attach the counter (r_c activated by m_d2 -> m_c) and the neuron
% (r_Nw activated by m_c -> m_out), reusing the BRENDA seed table T as the
% activator->reaction adjacency (same regulatory relay as the IFFL).
F = attach_counter_neuron(C, T, model, opt);
fprintf('STEP 4: %d full feature->counter->neuron pipelines.\n', height(F));
if height(F)>0
    disp(head(F, min(20,height(F))));
    writetable(F,'iffl_full_pipeline.csv');
    fprintf('Saved iffl_full_pipeline.csv\n');
else
    fprintf('No pipeline: no reaction is activated by a D2 output (counter), or none by m_c (neuron).\n');
end

% =============================== STEP 4b =====================================
% Mass-flow housekeeping (NO activation): counter reset r_I1 (consumes m_c),
% input sequestration r_I2 (consumes m_d2), neuron bias r_Nb (produces m_out),
% neuron removal r_Nr (consumes m_out). Pure stoichiometry.
if height(F)>0
    F = augment_housekeeping(F, model, opt);
    writetable(F,'iffl_full_pipeline.csv');     % overwrite with the complete block
    fprintf('STEP 4b: housekeeping attached; complete blocks: %d/%d\n', ...
            sum(F.complete), height(F));
    disp(head(F, min(15,height(F))));
end

% =============================== STEP 5 ======================================
% Rank by a collision-free cAMP->m_in connector and export the best full chain.
if exist('iffl_motifs.csv','file') && exist('iffl_full_pipeline.csv','file')
    try
        export_cAMP_ranked_chain(model, 'iffl_motifs.csv', 'iffl_full_pipeline.csv');
    catch ME
        fprintf('STEP 5 (export) skipped: %s\n', ME.message);
    end
end

% ====================== local function: STEP 2 ==============================
function M = build_iffl_from_anchor(model, T, opt)
    nMet = numel(model.mets);
    curr = currencyMask(model, opt);
    ub = ones(numel(model.rxns),1); if isfield(model,'ub'), ub = model.ub(:); end
    prodByMet = producerIndex(model.S, ub, opt.forwardOnly);
    metNorm = cellfun(@normLigand, model.metNames,'uni',0);
    name2met = containers.Map('KeyType','char','ValueType','any');
    for i=1:nMet, k=metNorm{i};
        if isKey(name2met,k), name2met(k)=[name2met(k) i]; else, name2met(k)=i; end
    end
    rows = {};
    for r = 1:height(T)
        a = find(strcmp(model.rxns, T.rxn{r}),1);  if isempty(a), continue; end
        regs = strsplit(T.regulators{r}, ';');  mFdset = [];
        for g=1:numel(regs), k=strtrim(regs{g});
            if isKey(name2met,k), mFdset=[mFdset name2met(k)]; end %#ok<AGROW>
        end
        mFdset = unique(mFdset(~curr(mFdset)));
        outs = find(model.S(:,a) < 0);  outs = outs(~curr(outs));
        for mFd = mFdset(:)'
            [fD,fP] = backBFS(mFd, opt.maxFast, model.S, prodByMet, curr);
            for md1 = outs(:)'
                if md1==mFd, continue; end
                [sD,sP] = backBFS(md1, opt.maxSlow, model.S, prodByMet, curr);
                common = setdiff(intersect(find(~isnan(fD)),find(~isnan(sD))),[mFd md1]);
                best=[]; bestScore=inf;
                for c = common(:)'
                    if fD(c)<1 || fD(c) >= sD(c), continue; end
                    fc = tracePath(c,fP,mFd,'rxn');  sc = tracePath(c,sP,md1,'rxn');
                    if opt.disjointBranches && ~isempty(intersect(fc,sc)), continue; end
                    if opt.metDisjointBranch
                        fm = tracePath(c,fP,mFd,'met'); sm = tracePath(c,sP,md1,'met');
                        inter = setdiff(intersect(setdiff(fm,c),setdiff(sm,c)), find(curr));
                        if ~isempty(inter), continue; end
                    end
                    s = fD(c)+sD(c);
                    if s<bestScore, bestScore=s; best=c; bFC=fc; bSC=sc; end
                end
                if ~isempty(best)
                    fm = tracePath(best,fP,mFd,'met'); sm = tracePath(best,sP,md1,'met');
                    metset = unique([fm sm]);
                    g = unionGenes(model,[a bFC bSC]);
                    rows(end+1,:) = { T.rxn{r}, T.EC{r}, model.metNames{best}, ...
                        model.metNames{mFd}, model.metNames{md1}, fD(best), sD(best), ...
                        strjoin(model.rxns(bFC),'>'), strjoin(model.rxns(bSC),'>'), g, ...
                        strjoin(model.mets(metset),';') }; %#ok<AGROW>
                end
            end
        end
    end
    M = packMotifs(rows, opt);
end

% ====================== local function: STEP 3 ==============================
function C = concatenate_diffs(M, model, opt)
    if height(M)==0, C = emptyCascade(); return; end
    rxOf  = @(i) [ {M.anchor{i}}, split2(M.fastChain{i}), split2(M.slowChain{i}) ];
    metOf = @(i) split2b(M.metIDs{i});
    rows = {};
    for i = 1:height(M)
        for j = 1:height(M)
            if strcmp(M.anchor{i},M.anchor{j}), continue; end
            if ~strcmp(strtrim(M.m_d1{i}), strtrim(M.input{j})), continue; end
            if opt.cascadeDisjoint && ~isempty(intersect(rxOf(i),rxOf(j))), continue; end
            if opt.metDisjointCasc
                shared = intersect(metOf(i), metOf(j));
                % the link metabolite is allowed to be shared; nothing else
                linkNorm = normLigand(M.m_d1{i});
                shared = shared(~strcmp(cellfun(@normLigand,shared,'uni',0), linkNorm));
                if ~isempty(shared), continue; end
            end
            idx = find(ismember(model.rxns, unique([rxOf(i) rxOf(j)])));
            g = unionGenes(model, idx);
            rows(end+1,:) = { M.anchor{i}, M.input{i}, M.m_Fd{i}, M.m_d1{i}, ...
                              M.anchor{j}, M.m_Fd{j}, M.m_d1{j}, M.m_d1{i}, g }; %#ok<AGROW>
        end
    end
    if isempty(rows), C = emptyCascade();
    else
        C = cell2table(rows,'VariableNames',{'D1_anchor','D1_input','D1_mFd','D1_out', ...
            'D2_anchor','D2_mFd','D2_out','link','nGenesTotal'});
        C = sortrows(C,'nGenesTotal','ascend');
    end
end

% ====================== local function: STEP 4 ==============================
function F = attach_counter_neuron(C, T, model, opt)
    if height(C)==0, F = emptyPipe(); return; end
    curr = currencyMask(model, opt);
    % activator (normalised metabolite name) -> reaction names it activates
    act = containers.Map('KeyType','char','ValueType','any');
    for r = 1:height(T)
        ligs = strsplit(T.regulators{r}, ';');
        for q = 1:numel(ligs)
            k = normLigand(ligs{q}); if isempty(k), continue; end
            if isKey(act,k), act(k) = [act(k) {T.rxn{r}}]; else, act(k) = {T.rxn{r}}; end
        end
    end
    r2i = containers.Map(model.rxns, num2cell(1:numel(model.rxns)));
    rows = {}; CAP = 20000;
    for i = 1:height(C)
        used = {C.D1_anchor{i}, C.D2_anchor{i}};        % don't reuse upstream reactions
        k2 = normLigand(C.D2_out{i});                    % m_d2 activates the counter
        if ~isKey(act,k2), continue; end
        for rc = unique(act(k2))
            rc = rc{1}; if any(strcmp(rc,used)) || ~isKey(r2i,rc), continue; end
            ci = r2i(rc);
            mcs = find(model.S(:,ci) > 0); mcs = mcs(~curr(mcs));   % m_c candidates
            for mc = mcs(:)'
                kc = normLigand(model.metNames{mc});     % m_c activates the neuron
                if ~isKey(act,kc), continue; end
                for rw = unique(act(kc))
                    rw = rw{1}; if any(strcmp(rw,[used {rc}])) || ~isKey(r2i,rw), continue; end
                    wi = r2i(rw);
                    mos = find(model.S(:,wi) > 0); mos = mos(~curr(mos)); % m_out candidates
                    for mo = mos(:)'
                        idx = [r2i(C.D1_anchor{i}) r2i(C.D2_anchor{i}) ci wi];
                        rows(end+1,:) = { C.D1_anchor{i}, C.D1_input{i}, C.D1_mFd{i}, C.D1_out{i}, ...
                            C.D2_anchor{i}, C.D2_mFd{i}, C.D2_out{i}, ...
                            rc, model.metNames{mc}, rw, model.metNames{mo}, unionGenes(model,idx) }; %#ok<AGROW>
                        if size(rows,1) >= CAP, break; end
                    end
                    if size(rows,1) >= CAP, break; end
                end
                if size(rows,1) >= CAP, break; end
            end
            if size(rows,1) >= CAP, break; end
        end
        if size(rows,1) >= CAP, break; end
    end
    if isempty(rows), F = emptyPipe();
    else
        F = cell2table(rows,'VariableNames',{'D1_anchor','D1_input','D1_mFd','link', ...
            'D2_anchor','D2_mFd','m_d2','counter_rc','m_c','neuron_rNw','m_out','nGenesAdded'});
        F = unique(F); F = sortrows(F,'nGenesAdded','ascend');
    end
end
function F = emptyPipe()
    F = cell2table(cell(0,12),'VariableNames',{'D1_anchor','D1_input','D1_mFd','link', ...
        'D2_anchor','D2_mFd','m_d2','counter_rc','m_c','neuron_rNw','m_out','nGenesAdded'});
end

% ====================== local function: STEP 4b =============================
function F = augment_housekeeping(F, model, opt)
    curr = currencyMask(model, opt);
    ub = ones(numel(model.rxns),1); if isfield(model,'ub'), ub = model.ub(:); end
    consByMet = consumerIndex(model.S, ub, opt.forwardOnly);
    prodByMet = producerIndex(model.S, ub, opt.forwardOnly);
    % normalised metabolite name -> all metabolite indices (any compartment)
    name2met = containers.Map('KeyType','char','ValueType','any');
    for i=1:numel(model.mets), k=normLigand(model.metNames{i});
        if isKey(name2met,k), name2met(k)=[name2met(k) i]; else, name2met(k)=i; end
    end
    gcount = @(j) geneCount(model,j);
    n = height(F);
    [reset,seq,bias,rem] = deal(repmat({''},n,1));
    blockGenes = zeros(n,1); complete = false(n,1);
    for r = 1:n
        used = {F.D1_anchor{r},F.D2_anchor{r},F.counter_rc{r},F.neuron_rNw{r}};
        mcIx  = metIx(name2met, F.m_c{r});
        md2Ix = metIx(name2met, F.m_d2{r});
        moIx  = metIx(name2met, F.m_out{r});
        reset{r} = bestCand(consByMet, mcIx,  used, model, gcount);        % r_I1 consumes m_c
        seq{r}   = bestCand(consByMet, md2Ix, used, model, gcount);        % r_I2 consumes m_d2
        rem{r}   = bestCand(consByMet, moIx,  [used reset(r)' seq(r)'], model, gcount); % r_Nr consumes m_out
        bias{r}  = bestCand(prodByMet, moIx,  [used reset(r)' seq(r)' rem(r)'], model, gcount); % r_Nb makes m_out
        picks = [used reset(r)' seq(r)' bias(r)' rem(r)'];
        picks = picks(~cellfun(@isempty,picks));
        idx = find(ismember(model.rxns,picks));
        blockGenes(r) = unionGenes(model, idx);
        complete(r) = ~isempty(reset{r}) && ~isempty(seq{r}) && ~isempty(rem{r}); % bias optional (b=0 ok)
    end
    F.reset_rI1 = reset; F.seq_rI2 = seq; F.bias_rNb = bias; F.rem_rNr = rem;
    F.nGenesBlock = blockGenes; F.complete = complete;
    F = sortrows(F, {'complete','nGenesBlock'}, {'descend','ascend'});
end

function ix = metIx(name2met, nm)
    k = normLigand(nm); if isKey(name2met,k), ix = name2met(k); else, ix = []; end
end
function name = bestCand(byMet, metIxs, excludeNames, model, gcount)
    name = ''; bestG = inf; seen = {};
    for m = metIxs(:)'
        for j = byMet{m}(:)'
            rn = model.rxns{j};
            if any(strcmp(rn, excludeNames)) || any(strcmp(rn,seen)), continue; end
            seen{end+1} = rn; %#ok<AGROW>
            g = gcount(j);
            if g < bestG, bestG = g; name = rn; end
        end
    end
end
function g = geneCount(model,j)
    if isfield(model,'rxnGeneMat')&&~isempty(model.rxnGeneMat), g=nnz(model.rxnGeneMat(j,:)); else, g=1; end
end
function P = consumerIndex(S,ub,forwardOnly)
    nMet=size(S,1); P=cell(nMet,1);
    for m=1:nMet, js=find(S(m,:)<0); if forwardOnly, js=js(ub(js)>0); end, P{m}=js; end
end

% ------------------------------- helpers -----------------------------------
function M = packMotifs(rows, opt)
    vn = {'anchor','EC','input','m_Fd','m_d1','dFast','dSlow','fastChain','slowChain','nGenes','metIDs'};
    if isempty(rows), M = cell2table(cell(0,numel(vn)),'VariableNames',vn); return; end
    M = cell2table(rows,'VariableNames',vn);
    M = sortrows(M,{'nGenes','dSlow'},{'ascend','ascend'});
    if height(M)>opt.topN, M = M(1:opt.topN,:); end
end
function C = emptyCascade()
    C = cell2table(cell(0,9),'VariableNames',{'D1_anchor','D1_input','D1_mFd','D1_out', ...
        'D2_anchor','D2_mFd','D2_out','link','nGenesTotal'});
end
function p = split2(s),  p = strsplit(char(s),'>'); end
function p = split2b(s), p = strsplit(char(s),';'); end

function [depth,pred] = backBFS(target, maxd, S, prodByMet, curr)
    nMet=size(S,1); depth=nan(nMet,1); pred=cell(nMet,1);
    depth(target)=0; frontier=target;
    for d=1:maxd
        nxt=[];
        for mm=frontier(:)'
            for j=prodByMet{mm}(:)'
                for s=find(S(:,j)<0)'
                    if curr(s)||~isnan(depth(s)), continue; end
                    depth(s)=d; pred{s}=[j mm]; nxt(end+1)=s; %#ok<AGROW>
                end
            end
        end
        frontier=nxt; if isempty(frontier), break; end
    end
end
function out = tracePath(c,pred,target,kind)
    % kind='rxn' returns reaction indices along the path; 'met' returns metabolite indices
    rx=[]; mt=c; cur=c;
    while cur~=target
        pc=pred{cur}; if isempty(pc), out=[]; return; end
        rx(end+1)=pc(1); mt(end+1)=pc(2); cur=pc(2); %#ok<AGROW>
    end
    if strcmp(kind,'rxn'), out=rx; else, out=mt; end
end
function P = producerIndex(S,ub,forwardOnly)
    nMet=size(S,1); P=cell(nMet,1);
    for m=1:nMet, js=find(S(m,:)>0); if forwardOnly, js=js(ub(js)>0); end, P{m}=js; end
end
function curr = currencyMask(model, opt)
    base = lower(regexprep(model.mets,'\[[a-z]\]$',''));
    curr = ismember(base, lower(opt.currency));
    if isfield(opt,'currencyNamePatterns') && ~isempty(opt.currencyNamePatterns)
        nm = lower(model.metNames);
        for p = 1:numel(opt.currencyNamePatterns)
            curr = curr | contains(nm, lower(opt.currencyNamePatterns{p}));
        end
    end
end
function n = unionGenes(model,rxnIdx)
    if isfield(model,'rxnGeneMat')&&~isempty(model.rxnGeneMat)
        n=nnz(any(model.rxnGeneMat(rxnIdx,:)~=0,1));
    else, n=numel(rxnIdx); end
end
function s = normLigand(s)
    s=lower(strtrim(char(s)));
    s=regexprep(s,'^((l|d|dl|\(s\)|\(r\)|alpha|beta)-)+','');
    s=regexprep(s,'[^a-z0-9]+',' '); s=strtrim(s);
end
