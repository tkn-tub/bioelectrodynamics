function export_cAMP_ranked_chain(modelIn, motifsFile, pipeFile)
% EXPORT_CAMP_RANKED_CHAIN  Rank full pipelines by a COLLISION-FREE cAMP->m_in
% connector (Block 0 may not touch any block 1-4 signal metabolite) and export
% the best as a human-readable full-chain .txt, in the style of Block_7.
%
%   export_cAMP_ranked_chain()                    % locate model + CSVs via which()
%   export_cAMP_ranked_chain(model,'iffl_motifs.csv','iffl_full_pipeline.csv')
%
% Writes: CAMP_Input_Ranking.csv, Selected_Full_Chain_<n>.txt
%
% -------------------------------------------------------------------------
% NOMENCLATURE (one place: roleTag() below -- edit it to change everything)
% -------------------------------------------------------------------------
% The search builds branches of VARIABLE length, so diagramCircuitSchema's
% fixed slot list cannot be used directly. The rule here is:
%
%   every pool is named after the reaction that produces it:  r_X  ->  m_X
%
%   D1 fast arm :  r_D1_Ff -> m_D1_Ff -> r_D1_Fd -> m_D1_Fd   (gate)
%   D1 slow arm :  r_D1_Sf -> m_D1_Sf -> r_D1_Sd -> m_d1       (output)
%   D1 anchor   :  r_D1_I    saturated degrader of m_d1, activated by m_D1_Fd
%   D1 removal  :  r_D1_Fr   clears m_D1_Fd      (diagram: fast rail removal)
%                  r_D1_lk   leak on m_d1        (diagram: r_lk)
%   D2 likewise with the D2_ tag; its slow arm terminates on m_d2.
%
%   Longer arms get a STEP index on the pool only:  m_D1_Ff1, m_D1_Ff2, ...
%   The differentiator index is always the D<k>_ tag, never a trailing digit,
%   so "m_D1_Ff2" (pool 2 of D1) can never be confused with "m_D2_Ff".
%
%   Interface signals that cross blocks keep their global, untagged names:
%   m_in, m_d1, m_d2, m_c, m_I, m_b, m_out   (same as diagramCircuitSchema).
%
% Mapping to diagramCircuitSchema.m (annihilation topology):
%   m_D1_Ff ~ m_F ,  m_D1_Fd ~ m_P1 ,  m_D1_Sf ~ m_S ,  m_d1 ~ m_d1
%   m_D2_Ff ~ m_F2,  m_D2_Sf ~ m_S2
%   NOTE the polarity difference: in diagramCircuitSchema D2's OUTPUT comes
%   off the FAST rail (r_Fd_2 -> m_d2). In this IFFL master BOTH
%   differentiators output the SLOW rail and gate it with the fast rail, so
%   m_d2 here is the slow-arm endpoint of D2.

    if nargin==1 && (ischar(modelIn)||isstring(modelIn)) && strcmpi(char(modelIn),'selftest')
        rolesSelfTest(); return;                 % nomenclature check, no model needed
    end
    if nargin<1 || isempty(modelIn)
        addingPathParentFolderByName('code');
        if exist('readCbModel','file')~=2, initCobraToolbox(false); end
        mf = which('Puniya2021_TNM1055.xml'); if isempty(mf), mf=which('Puniya2021_TNM1055.mat'); end
        assert(~isempty(mf),'model not found'); [~,~,e]=fileparts(mf);
        if strcmpi(e,'.mat'), L=load(mf); fn=fieldnames(L); model=L.(fn{1}); else, model=readCbModel(mf); end
    else
        model = modelIn;
    end
    if nargin<2 || isempty(motifsFile), motifsFile = which('iffl_motifs.csv'); end
    if nargin<3 || isempty(pipeFile),   pipeFile   = which('iffl_full_pipeline.csv'); end
    M = readtable(motifsFile,'TextType','string');
    F = readtable(pipeFile,'TextType','string');

    % ---- settings ----
    cAMP = 'camp[c]';  AMP = 'amp[c]';   maxSteps = 8;   topN = 3;
    maxEval = 4000;     % cap on complete pipelines evaluated (simplest first)
    inCurrency = {'h','h2o','pi','ppi','co2','o2','nh4','nh3','hco3','so4', ...
                  'na1','k','cl','fe2','fe3','mg2','ca2','zn2'};

    r2i = containers.Map(cellstr(model.rxns), num2cell(1:numel(model.rxns)));
    metByName = containers.Map('KeyType','char','ValueType','double');
    for i=1:numel(model.mets)
        k=char(model.metNames{i}); if ~isKey(metByName,k), metByName(k)=i; end
    end
    ci = find(strcmp(model.mets, cAMP),1);
    assert(~isempty(ci), 'cAMP (%s) not in model.', cAMP);

    % ---- global reachability from cAMP (CHEAP prefilter + diagnostics only) ----
    dist = fwdBFS(model, ci, maxSteps, inCurrency, false, false);

    % ---- per-pipeline COLLISION-FREE Block-0 search ----
    F = F(F.complete==1, :);
    assert(any(strcmp('D1_input', F.Properties.VariableNames)), ...
        ['iffl_full_pipeline.csv has no D1_input column - it is from an OLD master run. ' ...
         'Re-run run_IFFL_master.m and retry.']);
    F = sortrows(F,'nGenesBlock','ascend');          % evaluate simplest first
    rows = {};  nEval = 0;
    for r = 1:height(F)
        if nEval >= maxEval, break; end
        inIx = resolveMet(model, char(F.D1_input(r)));
        if isempty(inIx) || isnan(dist(inIx)), continue; end     % not globally reachable
        d1 = findMotif(M, F.D1_anchor(r), F.D1_input(r), F.link(r));
        d2 = findMotif(M, F.D2_anchor(r), F.link(r),     F.m_d2(r));
        if isempty(d1)||isempty(d2), continue; end
        nEval = nEval + 1;
        roles = buildRoles(model, d1, d2, F(r,:), r2i, false);   % block 1-4 signals
        forbid = false(numel(model.mets),1);
        kk = keys(roles.met);
        if ~isempty(kk), forbid(cell2mat(kk)) = true; end
        forbid(inIx) = false;                      % m_in is the target, allowed
        csteps = constrainedFwdBFS(model, ci, inIx, maxSteps, inCurrency, forbid);
        if isnan(csteps), continue; end            % no clean Block-0 for this circuit
        rows(end+1,:) = { F.D1_anchor(r), F.D2_anchor(r), string(char(F.D1_input(r))), ...
            csteps, F.nGenesBlock(r), r }; %#ok<AGROW>
    end
    fprintf('Evaluated %d pipelines; %d have a collision-free cAMP->m_in connector.\n', nEval, size(rows,1));
    if isempty(rows)
        fprintf('\n*** No complete pipeline input is reachable from cAMP within %d steps. ***\n',maxSteps);
        ai = find(strcmp(model.mets,AMP),1);
        if ~isempty(ai)
            fprintf('AMP (%s) distance from cAMP: %s\n', AMP, num2str(dist(ai)));
        end
        inNames = unique(string(M.input));
        dd = nan(numel(inNames),1);
        for q=1:numel(inNames)
            if isKey(metByName,char(inNames(q))), dd(q)=dist(metByName(char(inNames(q)))); end
        end
        tbl = sortrows(table(inNames,dd,'VariableNames',{'input','cAMP_steps'}),'cAMP_steps');
        fprintf('Nearest motif inputs to cAMP (NaN = unreachable within %d):\n',maxSteps);
        disp(head(tbl, min(15,height(tbl))));
        fprintf(['Diagnosis: EC/GPR is NOT the gate here (connector walks all reactions).\n' ...
                 'If every input is NaN/large, the automated differentiators simply do not\n' ...
                 'anchor near cAMP - the clean adenylate input (AMP) is the curated motif.\n' ...
                 'Raise maxSteps if a few inputs sit just beyond it.\n']);
        writetable(tbl,'CAMP_Input_Distances.csv');
        return;
    end
    R = cell2table(rows,'VariableNames',{'D1_anchor','D2_anchor','input','cAMP_steps', ...
        'nGenesBlock','Frow'});
    R = sortrows(R, {'cAMP_steps','nGenesBlock'}, {'ascend','ascend'});
    writetable(R(:,1:5), 'CAMP_Input_Ranking.csv');
    fprintf('Ranked %d cAMP-reachable pipelines. Best input distance = %d step(s).\n', height(R), R.cAMP_steps(1));

    % ---- export the top-N as full-chain .txt ----
    for n = 1:min(topN,height(R))
        fn = sprintf('Selected_Full_Chain_%d.txt', n);
        writeChain(fn, model, M, F, R(n,:), r2i, ci, AMP, inCurrency, n, height(R));
        fprintf('Wrote %s  (input %s, %d step(s) from cAMP)\n', fn, R.input(n), R.cAMP_steps(n));
    end
end

% =========================================================================
% NOMENCLATURE
% =========================================================================
function s = roleTag(k, base)
% Block-local role name for differentiator k. Change this one function to
% switch conventions, e.g. to the diagramCircuitSchema style:
%    if k==1, s = base; else, s = sprintf('%s_%d',base,k); end
    s = sprintf('D%d_%s', k, base);
end

% =========================================================================
% ROLE TABLE: metabolite roles, reaction roles, per-block print order
% =========================================================================
function R = buildRoles(model, d1, d2, Fr, r2i, withExtras)
% R.met       containers.Map  metIdx -> 'm_...'
% R.metBlock  containers.Map  metIdx -> block number (0 = cross-block interface)
% R.ordIdx/R.ordLab           metabolite legend, in diagram order
% R.blkRxn{b}/R.blkLab{b}     reactions of block b in diagram order
%                             (NaN entry = diagram role with no mapped reaction)
% R.note                      diagnostics (walk failures, role collisions)
    if nargin<6, withExtras = false; end
    curr = legendCurrency(model);
    R.met      = containers.Map('KeyType','double','ValueType','char');
    R.metBlock = containers.Map('KeyType','double','ValueType','double');
    R.ordIdx = []; R.ordLab = {}; R.note = {};
    R.blkRxn = {[],[],[],[]}; R.blkLab = {{},{},{},{}};
    R.blkTitle = { 'BLOCK 1: DIFFERENTIATOR 1 (IFFL, OUTPUT m_d1)', ...
                   'BLOCK 2: DIFFERENTIATOR 2 (INPUT m_d1, OUTPUT m_d2)', ...
                   'BLOCK 3: COUNTER (OUTPUT m_c, SHARED SPECIES m_I)', ...
                   'BLOCK 4: SINGLE NEURON (OUTPUT m_out)' };

    mIn  = resolveMet(model,char(d1.input));   mFd1 = resolveMet(model,char(d1.m_Fd));
    mD1  = resolveMet(model,char(d1.m_d1));    mFd2 = resolveMet(model,char(d2.m_Fd));
    mD2  = resolveMet(model,char(d2.m_d1));    mC   = resolveMet(model,char(Fr.m_c));
    mOut = resolveMet(model,char(Fr.m_out));
    mI   = otherSub(model, strOr(Fr.reset_rI1), r2i, mC,  curr);   % counter shared species
    mB   = otherSub(model, strOr(Fr.bias_rNb),  r2i, mOut, curr);  % bias source

    used = idxOf(r2i, [ {char(d1.anchor)}, splitr(d1.fastChain), splitr(d1.slowChain), ...
                        {char(d2.anchor)}, splitr(d2.fastChain), splitr(d2.slowChain), ...
                        {char(Fr.counter_rc)}, {char(Fr.neuron_rNw)}, ...
                        {strOr(Fr.reset_rI1)}, {strOr(Fr.seq_rI2)}, ...
                        {strOr(Fr.bias_rNb)},  {strOr(Fr.rem_rNr)} ]);

    % ---- m_in is the interface into block 1 ----
    R = addMet(R, model, mIn, 'm_in', 0);

    % ---- BLOCK 1 -------------------------------------------------------
    R = armRoles(R, model, r2i, curr, splitr(d1.fastChain), mIn, mFd1, 1, 'F', 1, ...
                 ['m_' roleTag(1,'Fd')], 1);
    [R, used] = extraRole(R, model, 1, mFd1, used, ['r_' roleTag(1,'Fr')], ...
                 'fast-rail removal (clears the gate)', withExtras);
    R = armRoles(R, model, r2i, curr, splitr(d1.slowChain), mIn, mD1, 1, 'S', 1, 'm_d1', 0);
    R = pushRxn(R, 1, rxIdx(r2i,char(d1.anchor)), ['r_' roleTag(1,'I')]);
    [R, used] = extraRole(R, model, 1, mD1, used, ['r_' roleTag(1,'lk')], ...
                 'output leak on m_d1', withExtras);

    % ---- BLOCK 2 -------------------------------------------------------
    R = armRoles(R, model, r2i, curr, splitr(d2.fastChain), mD1, mFd2, 2, 'F', 2, ...
                 ['m_' roleTag(2,'Fd')], 2);
    [R, used] = extraRole(R, model, 2, mFd2, used, ['r_' roleTag(2,'Fr')], ...
                 'fast-rail removal (clears the gate)', withExtras);
    R = armRoles(R, model, r2i, curr, splitr(d2.slowChain), mD1, mD2, 2, 'S', 2, 'm_d2', 0);
    R = pushRxn(R, 2, rxIdx(r2i,char(d2.anchor)), ['r_' roleTag(2,'I')]);
    [R, used] = extraRole(R, model, 2, mD2, used, ['r_' roleTag(2,'lk')], ...
                 'output leak on m_d2', withExtras);

    % ---- BLOCK 3: counter ---------------------------------------------
    R = addMet(R, model, mC, 'm_c', 0);
    R = addMet(R, model, mI, 'm_I', 0);
    R = pushRxn(R, 3, rxIdx(r2i,char(Fr.counter_rc)),  'r_c');
    R = pushRxn(R, 3, rxIdx(r2i,strOr(Fr.reset_rI1)),  'r_I1');
    R = pushRxn(R, 3, rxIdx(r2i,strOr(Fr.seq_rI2)),    'r_I2');

    % ---- BLOCK 4: neuron ----------------------------------------------
    R = addMet(R, model, mB,   'm_b',   0);
    R = addMet(R, model, mOut, 'm_out', 0);
    R = pushRxn(R, 4, rxIdx(r2i,char(Fr.neuron_rNw)),  'r_Nw');
    R = pushRxn(R, 4, rxIdx(r2i,strOr(Fr.bias_rNb)),   'r_Nb');
    R = pushRxn(R, 4, rxIdx(r2i,strOr(Fr.rem_rNr)),    'r_Nr');

    % ---- diagram checks the search does NOT enforce ---------------------
    if withExtras
        j1 = rxIdx(r2i,strOr(Fr.reset_rI1));  j2 = rxIdx(r2i,strOr(Fr.seq_rI2));
        if ~isnan(j1) && ~isnan(j2)
            s1 = find(model.S(:,j1)<0); s2 = find(model.S(:,j2)<0);
            if isempty(intersect(s1(~curr(s1)), s2(~curr(s2))))
                R.note{end+1} = ['r_I1 and r_I2 share NO non-currency substrate: the diagram ' ...
                    'requires both to draw on the same m_I pool (that shared pool is what ' ...
                    'makes the counter subtractive). m_I above is r_I1''s substrate only.'];
            end
        end
    end
end

% ---- one arm (fast or slow) of differentiator k ------------------------
function R = armRoles(R, model, r2i, curr, chainRx, startIx, endIx, k, kind, blk, endLab, endBlk)
    [mSeq, jSeq] = walkChain(model, chainRx, startIx, endIx, r2i, curr);
    if isempty(jSeq)
        R.note{end+1} = sprintf(['D%d %s arm could not be walked from %s to %s using its own ' ...
            'chain reactions - its pools are left UNNAMED in this report.'], k, armWord(kind), ...
            nameOf(model,startIx), nameOf(model,endIx));
        % still print the chain reactions, unnamed, so nothing is silently dropped
        for q = 1:numel(chainRx)
            jj = rxIdx(r2i, chainRx{q});
            if isnan(jj), continue; end
            R = pushRxn(R, blk, jj, ...
                sprintf('r_%s (unresolved %s-arm step)', roleTag(k,kind), armWord(kind)));
        end
        R = addMet(R, model, endIx, endLab, endBlk);
        return;
    end
    nPool = numel(jSeq) - 1;                    % interior pools
    for t = 1:nPool
        if nPool > 1, base = sprintf('%sf%d', kind, t); else, base = sprintf('%sf', kind); end
        R = addMet(R, model, mSeq(t+1), ['m_' roleTag(k,base)], blk);  % r_X -> m_X
        R = pushRxn(R, blk, jSeq(t),    ['r_' roleTag(k,base)]);
    end
    R = pushRxn(R, blk, jSeq(end), ['r_' roleTag(k, sprintf('%sd',kind))]);
    R = addMet(R, model, endIx, endLab, endBlk);
end
function w = armWord(kind)
    if kind=='F', w = 'fast'; else, w = 'slow'; end
end

% ---- removal/leak reactions the search never selects on ----------------
function [R, used] = extraRole(R, model, blk, metIx, used, lab, what, withExtras)
    if ~withExtras || isempty(metIx), return; end
    j = lowestGeneConsumer(model, metIx, used);
    if isnan(j)
        R = pushRxn(R, blk, NaN, [lab ' -- NOT MAPPED (' what ')']);
        R.note{end+1} = sprintf(['%s has no usable consumer in the model: %s cannot be ' ...
            'realised, so %s never returns to baseline.'], nameOf(model,metIx), lab, nameOf(model,metIx));
    else
        used(end+1) = j; %#ok<AGROW>
        R = pushRxn(R, blk, j, [lab ' (' what '; suggested, NOT part of the selection)']);
    end
end

% ---- small accessors ---------------------------------------------------
function R = addMet(R, model, ix, lab, blk)
    if isempty(ix) || any(isnan(ix)), return; end
    ix = ix(1);
    if isKey(R.met, ix)
        if ~strcmp(R.met(ix), lab)
            R.note{end+1} = sprintf('ROLE COLLISION: %s (%s) is both %s and %s.', ...
                nameOf(model,ix), char(model.mets{ix}), R.met(ix), lab);
        end
        return;
    end
    R.met(ix) = lab; R.metBlock(ix) = blk;
    R.ordIdx(end+1) = ix; R.ordLab{end+1} = lab;
end
function R = pushRxn(R, blk, j, lab)
    if isempty(j), j = NaN; end
    R.blkRxn{blk}(end+1) = j; R.blkLab{blk}{end+1} = lab;
end
function j = rxIdx(r2i, name)
    j = NaN; name = strtrim(char(name));
    if ~isempty(name) && isKey(r2i,name), j = r2i(name); end
end
function v = idxOf(r2i, names)
    v = [];
    for q = 1:numel(names)
        j = rxIdx(r2i, names{q}); if ~isnan(j), v(end+1) = j; end %#ok<AGROW>
    end
    v = unique(v);
end
function s = nameOf(model, ix)
    if isempty(ix), s = '<none>'; else, s = char(model.metNames{ix(1)}); end
end
function j = lowestGeneConsumer(model, metIx, excludeIdx)
    j = NaN; bestG = inf;
    ub = ones(numel(model.rxns),1); if isfield(model,'ub'), ub = model.ub(:); end
    for q = find(model.S(metIx(1),:) < 0)
        if ub(q) <= 0 || any(excludeIdx == q), continue; end
        g = 1;
        if isfield(model,'rxnGeneMat') && ~isempty(model.rxnGeneMat), g = nnz(model.rxnGeneMat(q,:)); end
        if g < bestG, bestG = g; j = q; end
    end
end

% ---- walk a chain: metabolite sequence AND reaction sequence -----------
function [mSeq, jSeq] = walkChain(model, chainRx, startIx, endIx, r2i, curr)
% constrained forward BFS using ONLY the chain's reactions, start -> end
    mSeq = []; jSeq = [];
    if isempty(startIx) || isempty(endIx), return; end
    js = [];
    for q = 1:numel(chainRx)
        j = rxIdx(r2i, chainRx{q}); if ~isnan(j), js(end+1) = j; end %#ok<AGROW>
    end
    if isempty(js), return; end
    predM = containers.Map('KeyType','double','ValueType','double');
    predR = containers.Map('KeyType','double','ValueType','double');
    vis = startIx; frontier = startIx; found = false;
    for d = 1:numel(js)+1
        nx = [];
        for m = frontier(:)'
            for j = js
                if model.S(m,j) < 0
                    for p = find(model.S(:,j) > 0)'
                        if curr(p) || any(vis==p), continue; end
                        predM(p) = m; predR(p) = j; vis(end+1) = p; nx(end+1) = p; %#ok<AGROW>
                        if p == endIx, found = true; end
                    end
                end
            end
        end
        frontier = nx; if found || isempty(frontier), break; end
    end
    if ~isKey(predM, endIx), return; end
    mSeq = endIx; cur = endIx;
    while cur ~= startIx
        jSeq(end+1) = predR(cur); cur = predM(cur); mSeq(end+1) = cur; %#ok<AGROW>
    end
    mSeq = fliplr(mSeq); jSeq = fliplr(jSeq);
end

% =========================================================================
% BFS helpers (unchanged)
% =========================================================================
function [steps, rpath, mpath] = constrainedFwdBFS(model, src, tgt, maxd, inCurrency, forbid)
    steps=NaN; rpath=[]; mpath=[];
    base=lower(regexprep(model.mets,'\[[a-z]\]$','')); curr=ismember(base,lower(inCurrency));
    ub=ones(numel(model.rxns),1); if isfield(model,'ub'), ub=model.ub(:); end
    forbidRxn = full(any(model.S(forbid,:)~=0, 1))';        % reaction touches a forbidden signal
    usable = ~forbidRxn & ub>0;
    nMet=numel(model.mets); dd=nan(nMet,1); pr=nan(nMet,1); pm=nan(nMet,1);
    dd(src)=0; frontier=src;
    for d=1:maxd
        nxt=[];
        for m=frontier(:)'
            js=find(model.S(m,:)<0 & usable(:)');
            for j=js
                for p=find(model.S(:,j)>0)'
                    if curr(p)||~isnan(dd(p)), continue; end
                    dd(p)=d; pr(p)=j; pm(p)=m; nxt(end+1)=p; %#ok<AGROW>
                end
            end
        end
        frontier=nxt; if ~isnan(dd(tgt))||isempty(frontier), break; end
    end
    if isnan(dd(tgt)), return; end
    steps=dd(tgt); cur=tgt; mpath=tgt;
    while cur~=src, rpath(end+1)=pr(cur); cur=pm(cur); mpath(end+1)=cur; end %#ok<AGROW>
    rpath=fliplr(rpath); mpath=fliplr(mpath);
end

function [dist,predRxn,predMet] = fwdBFS(model, src, maxd, inCurrency, reqEC, reqGPR)
    if nargin<5, reqEC=false; end
    if nargin<6, reqGPR=false; end
    okRxn = true(numel(model.rxns),1);
    if reqEC && isfield(model,'rxnECNumbers'), okRxn = okRxn & ~cellfun(@isempty,model.rxnECNumbers(:)); end
    if reqGPR && isfield(model,'grRules'),     okRxn = okRxn & ~cellfun(@isempty,model.grRules(:));     end
    nMet=numel(model.mets); dist=nan(nMet,1); predRxn=nan(nMet,1); predMet=nan(nMet,1);
    base=lower(regexprep(model.mets,'\[[a-z]\]$','')); curr=ismember(base,lower(inCurrency));
    ub=ones(numel(model.rxns),1); if isfield(model,'ub'), ub=model.ub(:); end
    dist(src)=0; frontier=src;
    for d=1:maxd
        nxt=[];
        for m=frontier(:)'
            js=find(model.S(m,:)<0 & ub(:)'>0 & okRxn(:)');
            for j=js
                for p=find(model.S(:,j)>0)'
                    if curr(p)||~isnan(dist(p)), continue; end
                    dist(p)=d; predRxn(p)=j; predMet(p)=m; nxt(end+1)=p; %#ok<AGROW>
                end
            end
        end
        frontier=nxt; if isempty(frontier), break; end
    end
end

% ---------- pair D1 & D2 motifs by the link metabolite and D2 output ----------
function s = findMotif(M, anchor, inName, outName)
% unique match by anchor + input + output (m_d1)
    k = find(strcmp(M.anchor,anchor) & strcmp(M.input,inName) & strcmp(M.m_d1,outName));
    if isempty(k), s=[]; return; end
    s = table2struct(M(k(1),:)); s.row = k(1);
end

function ix = otherSub(model, rxnName, r2i, excludeIx, curr)
    ix=[]; if isempty(rxnName)||~isKey(r2i,rxnName), return; end
    j=r2i(rxnName);
    for m=find(model.S(:,j)<0)'
        if any(m==excludeIx)||curr(m), continue; end
        ix=m; return;
    end
end

function ix = resolveMet(model, name)
% name -> metabolite index, preferring the cytosolic [c] copy
    k=find(strcmp(model.metNames,name));
    if isempty(k), k=find(strcmp(model.mets,name)); end
    if isempty(k), ix=[]; return; end
    kc=k(endsWith(model.mets(k),'[c]'));
    if ~isempty(kc), ix=kc(1); else, ix=k(1); end
end

function curr = legendCurrency(model)
    base = lower(regexprep(model.mets,'\[[a-z]\]$',''));
    list = {'h','h2o','pi','ppi','co2','o2','nh4','nh3','hco3','so4','na1','k','cl', ...
            'fe2','fe3','mg2','ca2','zn2','atp','adp','gtp','gdp','utp','udp','ctp', ...
            'nad','nadh','nadp','nadph','coa','fad','fadh2','q10','q10h2'};
    curr = ismember(base, list);
end

% =========================================================================
% REPORT
% =========================================================================
function writeChain(fn, model, M, F, Rrow, r2i, ci, AMP, inCurrency, n, total)
    fid=fopen(fn,'w'); C=onCleanup(@() fclose(fid)); %#ok<NASGU>
    Fr = F(Rrow.Frow,:);
    d1 = findMotif(M, Fr.D1_anchor, Fr.D1_input, Fr.link);
    d2 = findMotif(M, Fr.D2_anchor, Fr.link,     Fr.m_d2);
    if isempty(d1)||isempty(d2), fprintf('  (chain detail unavailable for row %d)\n',Rrow.Frow); return; end

    R = buildRoles(model, d1, d2, Fr, r2i, true);

    fprintf(fid,'%s\n',repmat('#',1,60));
    fprintf(fid,'FULL CHEMICAL CHAIN (IFFL mapping, cAMP-ranked)\nExported %s\n', datestr(now));
    fprintf(fid,'Selection %d of %d | input %s | cAMP steps %d | block genes %d\n', ...
        n,total,char(Rrow.input),Rrow.cAMP_steps,Rrow.nGenesBlock);
    fprintf(fid,'Nomenclature: r_X produces m_X; differentiator index is the D<k>_ tag,\n');
    fprintf(fid,'trailing digits are STEP indices within one arm (m_D1_Ff2 = D1 fast pool 2).\n');
    fprintf(fid,'%s\n\n',repmat('#',1,60));

    % ---- SIGNAL METABOLITE LEGEND ----
    allRx = [];
    for b = 1:4, allRx = [allRx R.blkRxn{b}]; end %#ok<AGROW>
    allRx = allRx(~isnan(allRx));
    writeLegend(fid, model, R, unique(allRx));

    % ---- REACTION ROLE INDEX ----
    fprintf(fid,'SELECTED MOTIF: REACTION ROLES\n');
    for b = 1:4
        for q = 1:numel(R.blkRxn{b})
            j = R.blkRxn{b}(q); lab = R.blkLab{b}{q};
            if isnan(j)
                fprintf(fid,'  B%d  %-28s | %-12s | %s\n', b, lab, '-', '(no mapped reaction)');
            else
                ec=''; if isfield(model,'rxnECNumbers'), ec=char(model.rxnECNumbers{j}); end
                fprintf(fid,'  B%d  %-28s | %-12s | %s\n', b, lab, model.rxns{j}, ec);
            end
        end
    end
    fprintf(fid,'\n');

    % ---- DIAGNOSTICS ----
    if ~isempty(R.note)
        fprintf(fid,'DIAGNOSTICS\n');
        for q=1:numel(R.note), fprintf(fid,'  ! %s\n', R.note{q}); end
        fprintf(fid,'\n');
    end

    % ---- BLOCK 0 ----
    fprintf(fid,'------------------------------------------------------------\n');
    fprintf(fid,'BLOCK 0: cAMP INPUT CONNECTION (OUTPUT m_in = %s)\n',char(d1.input));
    fprintf(fid,'------------------------------------------------------------\n');
    inIx = resolveMet(model, char(d1.input));
    forbid = false(numel(model.mets),1);
    kk = keys(R.met); if ~isempty(kk), forbid(cell2mat(kk))=true; end
    forbid(inIx)=false;
    [csteps, rpath, mpath] = constrainedFwdBFS(model, ci, inIx, 12, inCurrency, forbid);
    if isnan(csteps)
        fprintf(fid,'  Status: NO collision-free connector within reach (every route consumes a signal).\n\n');
    else
        fprintf(fid,'  Status: reachable, collision-free | conversion steps: %d\n', csteps);
        fprintf(fid,'  Metabolite path: %s\n', strjoin(model.mets(mpath),' -> '));
        fprintf(fid,'  Reaction path  : %s\n', strjoin(model.rxns(rpath),' -> '));
        fprintf(fid,'  Routes through AMP (%s): %d\n\n', AMP, any(strcmp(model.mets(mpath),AMP)));
        for q=1:numel(rpath)
            printRxn(fid, model, rpath(q), sprintf('r_adapter_%d (input connector)',q), R, 0);
        end
    end

    % ---- BLOCKS 1-4 ----
    for b = 1:4, writeBlock(fid, model, R, b); end

    fprintf(fid,'\nASSUMPTIONS AND INTERFACES\n');
    fprintf(fid,'Differentiator subtraction is by activation (IFFL): m_Fd gates degradation of the output.\n');
    fprintf(fid,'Counter and neuron READS are by activation; housekeeping (reset/seq/removal) is mass-flow.\n');
    fprintf(fid,'r_Fr and r_lk are REPORTED here, not selected on: the master does not search for them.\n');
    fprintf(fid,'Connectivity validated: 1. Kinetics validated: 0. Activation strength validated: 0.\n');
    fprintf(fid,'Name-matched BRENDA regulation; verify AC strength and [X] << K_M before building.\n');
end

function writeLegend(fid, model, R, allRx)
    fprintf(fid,'SELECTED MOTIF: SIGNAL METABOLITES\n');
    for q=1:numel(R.ordIdx), metLine(fid, model, R.ordLab{q}, R.ordIdx(q)); end
    seen=containers.Map('KeyType','double','ValueType','logical');
    curr = legendCurrency(model);
    for c=1:numel(allRx)
        j=allRx(c);
        for m=find(model.S(:,j)~=0)'
            if curr(m)||isKey(R.met,m)||(isKey(seen,m)&&seen(m)), continue; end
            seen(m)=true; metLine(fid, model, 'm_aux', m);
        end
    end
    fprintf(fid,'\n');
end
function metLine(fid, model, label, ix)
    id = char(model.mets{ix}); nm = char(model.metNames{ix});
    f = ''; if isfield(model,'metFormulas') && ~isempty(model.metFormulas{ix}), f=char(model.metFormulas{ix}); end
    cm = regexp(id,'\[([a-z])\]$','tokens','once'); if isempty(cm), cm={'?'}; end
    fprintf(fid,'  %-10s | %-20s | %-42s | %-14s | %s\n', label, id, nm, f, cm{1});
end

function writeBlock(fid, model, R, b)
    fprintf(fid,'\n------------------------------------------------------------\n%s\n',R.blkTitle{b});
    fprintf(fid,'------------------------------------------------------------\n\n');
    foreign = {};
    for q = 1:numel(R.blkRxn{b})
        j = R.blkRxn{b}(q); lab = R.blkLab{b}{q};
        if isnan(j)
            fprintf(fid,'  %s\n      (no reaction in the model fills this diagram role)\n\n', lab);
            continue;
        end
        printRxn(fid, model, j, lab, R, b);
        for m = find(model.S(:,j)~=0)'
            if isKey(R.metBlock,m)
                ob = R.metBlock(m);
                if ob~=0 && ob~=b, foreign{end+1} = sprintf('%s (%s, block %d)', R.met(m), model.mets{m}, ob); end %#ok<AGROW>
            end
        end
    end
    if ~isempty(foreign)
        fprintf(fid,'  !! CROSSTALK: this block touches signal species owned by another block:\n');
        u = unique(foreign);
        for q=1:numel(u), fprintf(fid,'       %s\n', u{q}); end
        fprintf(fid,'     (the disjointness filters compare PATH metabolites only, not every\n');
        fprintf(fid,'      substrate/product of the chain reactions, so this is not screened out.)\n\n');
    end
end

function printRxn(fid, model, j, role, R, b)
    ec=''; if isfield(model,'rxnECNumbers'), ec=char(model.rxnECNumbers{j}); end
    nm=''; if isfield(model,'rxnNames'), nm=char(model.rxnNames{j}); end
    lb=-Inf; ub=Inf; if isfield(model,'lb'),lb=model.lb(j);end; if isfield(model,'ub'),ub=model.ub(j);end
    fprintf(fid,'  %s\n      ID  : %s\n      Name: %s\n      EC  : %s\n', role, model.rxns{j}, nm, ec);
    fprintf(fid,'      Bounds: [%g, %g]\n', lb, ub);
    fprintf(fid,'      IDs   : %s\n', eqnString(model, j, 'mets',     R, b));
    fprintf(fid,'      Names : %s\n\n', eqnString(model, j, 'metNames', R, b));
end
function s = eqnString(model, j, field, R, b)
    sub=find(model.S(:,j)<0); pro=find(model.S(:,j)>0);
    f=@(ix) strjoin(arrayfun(@(i) coefName(model,i,abs(model.S(i,j)),field,R,b), ix,'uni',0),' + ');
    s=sprintf('%s  ->  %s', f(sub), f(pro));
end
function t = coefName(model,i,c,field,R,b)
    nm=char(model.(field){i});
    if ~isempty(R) && isKey(R.met,i)
        tag = R.met(i);
        if R.metBlock(i)~=0 && R.metBlock(i)~=b, tag = ['!' tag]; end  % foreign-block signal
        nm=sprintf('%s(%s)', nm, tag);
    end
    if c==1, t=nm; else, t=sprintf('%g %s',c,nm); end
end

function p=splitr(s), p=strsplit(char(s),'>'); end
function s=strOr(x)
% empty bias/removal reads as <missing> (string column) or NaN (cellstr column)
    s = '';
    if isempty(x), return; end
    if iscell(x), x = x{1}; if isempty(x), return; end, end
    if isnumeric(x), return; end                      % NaN placeholder
    if ~ischar(x)
        try, if all(ismissing(x)), return; end, catch, end
    end
    s = strtrim(char(x));
    if strcmpi(s,'NaN'), s = ''; end
end

% =========================================================================
% SELF-TEST: synthetic 2-differentiator circuit, asserts the nomenclature
%   export_cAMP_ranked_chain('selftest')
% =========================================================================
function rolesSelfTest()
    % --- synthetic model: D1 fast arm 3 steps (to exercise pool numbering),
    %     D1 slow arm 2 steps, D2 both arms 2 steps -----------------------
    names = {'m_in','P1','P2','Fd1','Sf1','d1','Ff2','Fd2','Sf2','d2', ...
             'c','I','b','out','waste','H+'};
    model.metNames = names(:);
    model.mets = cellfun(@(s) [lower(regexprep(s,'[^A-Za-z0-9]','')) '[c]'], names(:), 'uni', 0);
    model.mets{16} = 'h[c]';                              % currency
    rx = {'rFf1','rFf2','rFd','rSf','rSd','rI','rFr','rLk', ...
          'rFf_2','rFd_2','rSf_2','rSd_2','rI_2','rc','rI1','rI2','rNw','rNb','rNr'};
    model.rxns = rx(:);
    model.rxnNames = rx(:);
    model.rxnECNumbers = repmat({'1.1.1.1'},numel(rx),1);
    S = zeros(numel(names),numel(rx));
    S = applyRxn(S,1,[2 16],1);      % rFf1 : m_in -> P1 (+H+)
    S = applyRxn(S,2,3,2);           % rFf2 : P1 -> P2
    S = applyRxn(S,3,4,3);           % rFd  : P2 -> Fd1
    S = applyRxn(S,1,5,4);           % rSf  : m_in -> Sf1
    S = applyRxn(S,5,6,5);           % rSd  : Sf1 -> d1
    S = applyRxn(S,6,15,6);          % rI   : d1 -> waste      (anchor D1)
    S = applyRxn(S,4,15,7);          % rFr  : Fd1 -> waste
    S = applyRxn(S,6,15,8);          % rLk  : d1 -> waste
    S = applyRxn(S,6,7,9);           % rFf_2: d1 -> Ff2
    S = applyRxn(S,7,8,10);          % rFd_2: Ff2 -> Fd2
    S = applyRxn(S,6,9,11);          % rSf_2: d1 -> Sf2
    S = applyRxn(S,9,10,12);         % rSd_2: Sf2 -> d2
    S = applyRxn(S,10,15,13);        % rI_2 : d2 -> waste      (anchor D2)
    S = applyRxn(S,10,11,14);        % rc   : d2 -> c
    S = applyRxn(S,[11 12],15,15);   % rI1  : c + I -> waste
    S = applyRxn(S,[10 12],15,16);   % rI2  : d2 + I -> waste
    S = applyRxn(S,12,14,17);        % rNw  : I -> out   (weight, activated by c)
    S = applyRxn(S,13,14,18);        % rNb  : b -> out
    S = applyRxn(S,14,15,19);        % rNr  : out -> waste
    model.S = S;
    model.lb = zeros(numel(rx),1); model.ub = 1000*ones(numel(rx),1);
    model.rxnGeneMat = speye(numel(rx));

    d1 = struct('input','m_in','m_Fd','Fd1','m_d1','d1', ...
                'fastChain','rFf1>rFf2>rFd','slowChain','rSf>rSd','anchor','rI');
    d2 = struct('input','d1','m_Fd','Fd2','m_d1','d2', ...
                'fastChain','rFf_2>rFd_2','slowChain','rSf_2>rSd_2','anchor','rI_2');
    Fr = struct('m_c','c','m_out','out','counter_rc','rc','neuron_rNw','rNw', ...
                'reset_rI1','rI1','seq_rI2','rI2','bias_rNb','rNb','rem_rNr','rNr');
    r2i = containers.Map(model.rxns, num2cell(1:numel(model.rxns)));

    R = buildRoles(model, d1, d2, Fr, r2i, true);

    fprintf('\n--- metabolite roles ---\n');
    for q=1:numel(R.ordIdx)
        fprintf('  %-10s -> %-6s (block %d)\n', R.ordLab{q}, ...
                model.metNames{R.ordIdx(q)}, R.metBlock(R.ordIdx(q)));
    end
    fprintf('--- reaction roles ---\n');
    for b=1:4
        for q=1:numel(R.blkRxn{b})
            j=R.blkRxn{b}(q);
            if isnan(j), nm='<none>'; else, nm=model.rxns{j}; end
            fprintf('  B%d %-34s -> %s\n', b, R.blkLab{b}{q}, nm);
        end
    end
    if ~isempty(R.note)
        fprintf('--- notes ---\n');
        for q=1:numel(R.note), fprintf('  ! %s\n', R.note{q}); end
    end

    lab = @(nm) R.met(find(strcmp(model.metNames,nm),1));
    assert(strcmp(lab('P1'),'m_D1_Ff1'), 'D1 fast pool 1 mislabelled');
    assert(strcmp(lab('P2'),'m_D1_Ff2'), 'D1 fast pool 2 mislabelled');
    assert(strcmp(lab('Fd1'),'m_D1_Fd'), 'D1 gate mislabelled');
    assert(strcmp(lab('Sf1'),'m_D1_Sf'), 'D1 slow pool mislabelled');
    assert(strcmp(lab('d1'),'m_d1') && strcmp(lab('d2'),'m_d2'), 'outputs mislabelled');
    assert(strcmp(lab('Ff2'),'m_D2_Ff'), 'D2 fast pool mislabelled');
    assert(strcmp(lab('Sf2'),'m_D2_Sf'), 'D2 slow pool mislabelled');
    assert(~isKey(R.met, find(strcmp(model.mets,'h[c]'),1)), 'currency got a role');
    want = {'r_D1_Ff1','r_D1_Ff2','r_D1_Fd','r_D1_Fr','r_D1_Sf','r_D1_Sd','r_D1_I','r_D1_lk'};
    got  = regexprep(R.blkLab{1},' .*','');
    assert(isequal(got(:)', want), 'block 1 reaction roles out of order or misnamed');
    assert(R.metBlock(find(strcmp(model.metNames,'Fd2'),1))==2, 'D2 gate assigned to the wrong block');

    % exercise the printing path (legend, tagged equations, crosstalk banner)
    tmp = [tempname '.txt']; fid = fopen(tmp,'w');
    allRx = []; for b=1:4, allRx = [allRx R.blkRxn{b}]; end %#ok<AGROW>
    writeLegend(fid, model, R, unique(allRx(~isnan(allRx))));
    for b = 1:4, writeBlock(fid, model, R, b); end
    fclose(fid);
    txt = fileread(tmp); delete(tmp);
    assert(~isempty(strfind(txt,'(m_D1_Fd)')), 'role tags missing from the printed equations'); %#ok<STREMP>
    fprintf('%s', txt(1:min(1200,numel(txt))));
    fprintf('\nself-test PASSED (nomenclature matches the differentiator diagram)\n\n');
end
function S = applyRxn(S, subs, pros, j)
    S(subs,j) = -1; S(pros,j) = 1;
end
