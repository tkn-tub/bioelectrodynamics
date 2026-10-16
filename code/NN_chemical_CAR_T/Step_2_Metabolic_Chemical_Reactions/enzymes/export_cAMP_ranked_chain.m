function export_cAMP_ranked_chain(modelIn, motifsFile, pipeFile)
% EXPORT_CAMP_RANKED_CHAIN  Rank full pipelines by a COLLISION-FREE cAMP->m_in
% connector (Block 0 may not touch any block 1-4 signal metabolite) and export
% the best as a human-readable full-chain .txt, in the style of Block_7.
%
%   export_cAMP_ranked_chain()                    % locate model + CSVs via which()
%   export_cAMP_ranked_chain(model,'iffl_motifs.csv','iffl_full_pipeline.csv')  % from the master
%
% Writes: CAMP_Input_Ranking.csv, Selected_Full_Chain_<n>.txt

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
    maxEval = 4000;     % cap on how many complete pipelines to evaluate (sorted simplest-first)
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
    [dist, predRxn, predMet] = fwdBFS(model, ci, maxSteps, inCurrency, false, false); %#ok<ASGLU>

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
        if isempty(inIx) || isnan(dist(inIx)), continue; end     % not even globally reachable
        d1 = findMotif(M, F.D1_anchor(r), F.D1_input(r), F.link(r));
        d2 = findMotif(M, F.D2_anchor(r), F.link(r),     F.m_d2(r));
        if isempty(d1)||isempty(d2), continue; end
        nEval = nEval + 1;
        roleMap = buildRoleMap(model, d1, d2, F(r,:), r2i);       % the block 1-4 signals
        forbid = false(numel(model.mets),1);
        forbid(cell2mat(keys(roleMap))) = true;  forbid(inIx) = false;  % m_in is the target, allowed
        csteps = constrainedFwdBFS(model, ci, inIx, maxSteps, inCurrency, forbid);
        if isnan(csteps), continue; end                           % no clean Block-0 for this circuit
        rows(end+1,:) = { F.D1_anchor(r), F.D2_anchor(r), string(char(F.D1_input(r))), csteps, ...
            F.nGenesBlock(r), r }; %#ok<AGROW>
    end
    fprintf('Evaluated %d pipelines; %d have a collision-free cAMP->m_in connector.\n', nEval, size(rows,1));
    if isempty(rows)
        fprintf('\n*** No complete pipeline input is reachable from cAMP within %d steps. ***\n',maxSteps);
        ai = find(strcmp(model.mets,AMP),1);
        if ~isempty(ai)
            fprintf('AMP (%s) distance from cAMP: %s\n', AMP, num2str(dist(ai)));
        end
        % distances of the distinct pipeline inputs (reconstructed) to cAMP
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

% ---------- forward BFS from cAMP: dist(met), predecessor reaction/met ----------
% ---------- collision-free forward BFS: cAMP -> target, forbidding signal mets ----------
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

% ---------- pick a motif row by anchor (and optional required input==linkName) ----------
function s = pickMotif(M, anchor, linkName)
    sel = M(strcmp(M.anchor, anchor), :);
    if ~isempty(linkName), sel = sel(strcmp(sel.input, linkName), :); end
    if isempty(sel), s=[]; return; end
    sel = sortrows(sel,'nGenes','ascend');
    s = table2struct(sel(1,:)); s.row = find(strcmp(M.anchor,anchor) & strcmp(M.input,sel.input(1)) & ...
                                             strcmp(M.m_d1,sel.m_d1(1)),1);
end

% ---------- pair D1 & D2 motifs by the link metabolite and D2 output ----------
function s = findMotif(M, anchor, inName, outName)
% unique match by anchor + input + output (m_d1)
    k = find(strcmp(M.anchor,anchor) & strcmp(M.input,inName) & strcmp(M.m_d1,outName));
    if isempty(k), s=[]; return; end
    s = table2struct(M(k(1),:)); s.row = k(1);
end

% ---------- build role map (metIdx -> label) by walking the IFFL chains ----------
function [roleMap, ordIdx, ordLab] = buildRoleMap(model, d1, d2, Fr, r2i)
    curr = legendCurrency(model);
    mIn  = resolveMet(model,char(d1.input));   mFd  = resolveMet(model,char(d1.m_Fd));
    mD1  = resolveMet(model,char(d1.m_d1));    mFd2 = resolveMet(model,char(d2.m_Fd));
    mD2  = resolveMet(model,char(d2.m_d1));    mC   = resolveMet(model,char(Fr.m_c));
    mOut = resolveMet(model,char(Fr.m_out));
    mI   = otherSub(model, strOr(Fr.reset_rI1), r2i, mC,  curr);   % counter shared species
    mB   = otherSub(model, strOr(Fr.bias_rNb),  r2i, mOut, curr);  % bias source
    % ordered role spec {idx,label}, pools walked from each branch
    specs = [ {mIn,'m_in'} ; ...
              poolSpecs(model,splitr(d1.fastChain),mIn,mFd ,r2i,curr,'m_Ff'  ) ; {mFd ,'m_Fd'  } ; ...
              poolSpecs(model,splitr(d1.slowChain),mIn,mD1 ,r2i,curr,'m_Sf'  ) ; {mD1 ,'m_D1'  } ; ...
              poolSpecs(model,splitr(d2.fastChain),mD1,mFd2,r2i,curr,'m_Ff_2') ; {mFd2,'m_Fd_2'} ; ...
              poolSpecs(model,splitr(d2.slowChain),mD1,mD2 ,r2i,curr,'m_Sf_2') ; {mD2 ,'m_D2'  } ; ...
              {mC,'m_c'} ; {mI,'m_I'} ; {mB,'m_b'} ; {mOut,'m_out'} ];
    roleMap = containers.Map('KeyType','double','ValueType','char');
    ordIdx=[]; ordLab={};
    for q=1:size(specs,1)
        ix=specs{q,1}; lab=specs{q,2};
        if isempty(ix)||isKey(roleMap,ix), continue; end
        roleMap(ix)=lab; ordIdx(end+1)=ix; ordLab{end+1}=lab; %#ok<AGROW>
    end
end

function specs = poolSpecs(model, chainRx, startIx, endIx, r2i, curr, baseLabel)
    specs = cell(0,2);
    seq = walkChainMets(model, chainRx, startIx, endIx, r2i, curr);
    interior = seq(2:end-1);            % exclude start and end
    for t=1:numel(interior)
        lab = baseLabel; if numel(interior)>1, lab=sprintf('%s%d',baseLabel,t); end
        specs = [specs; {interior(t), lab}]; %#ok<AGROW>
    end
end

function seq = walkChainMets(model, chainRx, startIx, endIx, r2i, curr)
% constrained forward BFS using ONLY the chain's reactions, start -> end
    seq=[startIx endIx];
    js=[]; for q=1:numel(chainRx), rn=chainRx{q}; if ~isempty(rn)&&isKey(r2i,rn), js(end+1)=r2i(rn); end, end %#ok<AGROW>
    if isempty(js), return; end
    predM=containers.Map('KeyType','double','ValueType','double'); vis=startIx; frontier=startIx; found=false;
    for d=1:numel(js)+1
        nx=[];
        for m=frontier(:)'
            for j=js
                if model.S(m,j)<0
                    for p=find(model.S(:,j)>0)'
                        if curr(p)||any(vis==p), continue; end
                        predM(p)=m; vis(end+1)=p; nx(end+1)=p; %#ok<AGROW>
                        if p==endIx, found=true; end
                    end
                end
            end
        end
        frontier=nx; if found||isempty(frontier), break; end
    end
    if ~isKey(predM,endIx), return; end
    path=endIx; cur=endIx;
    while cur~=startIx, cur=predM(cur); path(end+1)=cur; end %#ok<AGROW>
    seq=fliplr(path);
end

function ix = otherSub(model, rxnName, r2i, excludeIx, curr)
    ix=[]; if isempty(rxnName)||~isKey(r2i,rxnName), return; end
    j=r2i(rxnName);
    for m=find(model.S(:,j)<0)'
        if m==excludeIx||curr(m), continue; end
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

% ---------- signal-metabolite legend ----------
function writeLegend(fid, model, ordIdx, ordLab, roleMap, allRx, r2i)
    fprintf(fid,'SELECTED MOTIF: SIGNAL METABOLITES\n');
    for q=1:numel(ordIdx), metLine(fid, model, ordLab{q}, ordIdx(q)); end
    % any remaining non-currency chain metabolite not already a role
    curr = legendCurrency(model); seen=containers.Map('KeyType','double','ValueType','logical');
    for c=1:numel(allRx)
        rn=allRx{c}; if isempty(rn)||~isKey(r2i,rn), continue; end
        j=r2i(rn);
        for m=find(model.S(:,j)~=0)'
            if curr(m)||isKey(roleMap,m)||(isKey(seen,m)&&seen(m)), continue; end
            seen(m)=true; metLine(fid, model, 'm_aux', m);
        end
    end
    fprintf(fid,'\n');
end
function metLine(fid, model, label, ix)
    id = char(model.mets{ix}); nm = char(model.metNames{ix});
    f = ''; if isfield(model,'metFormulas') && ~isempty(model.metFormulas{ix}), f=char(model.metFormulas{ix}); end
    cm = regexp(id,'\[([a-z])\]$','tokens','once'); if isempty(cm), cm={'?'}; end
    fprintf(fid,'  %-7s | %-20s | %-42s | %-14s | %s\n', label, id, nm, f, cm{1});
end
function curr = legendCurrency(model)
    base = lower(regexprep(model.mets,'\[[a-z]\]$',''));
    list = {'h','h2o','pi','ppi','co2','o2','nh4','nh3','hco3','so4','na1','k','cl', ...
            'fe2','fe3','mg2','ca2','zn2','atp','adp','gtp','gdp','utp','udp','ctp', ...
            'nad','nadh','nadp','nadph','coa','fad','fadh2','q10','q10h2'};
    curr = ismember(base, list);
end

% ---------- write one full-chain txt ----------
function writeChain(fn, model, M, F, Rrow, r2i, ci, AMP, inCurrency, n, total)
    fid=fopen(fn,'w'); C=onCleanup(@() fclose(fid));
    Fr = F(Rrow.Frow,:);
    % fully-constrained motif lookup (anchor + input + output), from F's own columns
    d1 = findMotif(M, Fr.D1_anchor, Fr.D1_input, Fr.link);
    d2 = findMotif(M, Fr.D2_anchor, Fr.link,     Fr.m_d2);
    if isempty(d1)||isempty(d2), fprintf('  (chain detail unavailable for row %d)\n',Rrow.Frow); return; end
    hr=@() fprintf(fid,'%s\n',repmat('=',1,60));
    fprintf(fid,'%s\n',repmat('#',1,60));
    fprintf(fid,'FULL CHEMICAL CHAIN (IFFL mapping, cAMP-ranked)\nExported %s\n', datestr(now));
    fprintf(fid,'Selection %d of %d | input %s | cAMP steps %d | block genes %d\n', ...
        n,total,char(Rrow.input),Rrow.cAMP_steps,Rrow.nGenesBlock);
    fprintf(fid,'%s\n\n',repmat('#',1,60));

    % ---- SIGNAL METABOLITE LEGEND (role -> id | name | formula | compartment) ----
    allRx = unique([ {char(d1.anchor)}, splitr(d1.fastChain), splitr(d1.slowChain), ...
                     {char(d2.anchor)}, splitr(d2.fastChain), splitr(d2.slowChain), ...
                     {char(Fr.counter_rc), char(Fr.neuron_rNw)}, ...
                     strOr(Fr.reset_rI1), strOr(Fr.seq_rI2), strOr(Fr.bias_rNb), strOr(Fr.rem_rNr) ]);
    [roleMap, ordIdx, ordLab] = buildRoleMap(model, d1, d2, Fr, r2i);
    writeLegend(fid, model, ordIdx, ordLab, roleMap, allRx, r2i);

    % Block 0: COLLISION-FREE cAMP -> m_in path (forbids all block 1-4 signals)
    fprintf(fid,'------------------------------------------------------------\n');
    fprintf(fid,'BLOCK 0: cAMP INPUT CONNECTION (OUTPUT m_in = %s)\n',char(d1.input));
    fprintf(fid,'------------------------------------------------------------\n');
    inIx = resolveMet(model, char(d1.input));
    forbid = false(numel(model.mets),1); forbid(cell2mat(keys(roleMap)))=true; forbid(inIx)=false;
    [csteps, rpath, mpath] = constrainedFwdBFS(model, ci, inIx, 12, inCurrency, forbid);
    if isnan(csteps)
        fprintf(fid,'  Status: NO collision-free connector within reach (every route consumes a signal).\n\n');
    else
        fprintf(fid,'  Status: reachable, collision-free | conversion steps: %d\n', csteps);
        fprintf(fid,'  Metabolite path: %s\n', strjoin(model.mets(mpath),' -> '));
        fprintf(fid,'  Reaction path  : %s\n', strjoin(model.rxns(rpath),' -> '));
        fprintf(fid,'  Routes through AMP (%s): %d\n\n', AMP, any(strcmp(model.mets(mpath),AMP)));
        for j=rpath(:)', printRxn(fid, model, j, 'r_CAMP (input connector)', roleMap); end
    end

    % Blocks 1-4
    blk(fid,'BLOCK 1: DIFFERENTIATOR 1 (IFFL, OUTPUT m_d1)', model, ...
        [splitr(d1.fastChain) splitr(d1.slowChain) {char(d1.anchor)}], ...
        {'fast branch','...','slow branch','r_I (saturated degrader, activated by m_Fd)'}, r2i, roleMap);
    blk(fid,'BLOCK 2: DIFFERENTIATOR 2 (INPUT m_d1, OUTPUT m_d2)', model, ...
        [splitr(d2.fastChain) splitr(d2.slowChain) {char(d2.anchor)}], ...
        {'fast branch','...','slow branch','r_I_2 (activated by m_Fd2)'}, r2i, roleMap);
    blk(fid,'BLOCK 3: COUNTER (OUTPUT m_c)', model, ...
        {char(Fr.counter_rc), strOr(Fr.reset_rI1), strOr(Fr.seq_rI2)}, ...
        {'r_c (activated by m_d2)','r_I1 (reset, consumes m_c)','r_I2 (sequestration, consumes m_d2)'}, r2i, roleMap);
    blk(fid,'BLOCK 4: SINGLE NEURON (OUTPUT m_out)', model, ...
        {char(Fr.neuron_rNw), strOr(Fr.bias_rNb), strOr(Fr.rem_rNr)}, ...
        {'r_Nw (weight, activated by m_c)','r_Nb (bias)','r_Nr (removal, consumes m_out)'}, r2i, roleMap);

    fprintf(fid,'\nASSUMPTIONS AND INTERFACES\n');
    fprintf(fid,'Differentiator subtraction is by activation (IFFL): m_Fd gates degradation of m_d1.\n');
    fprintf(fid,'Counter and neuron READS are by activation; housekeeping (reset/seq/removal) is mass-flow.\n');
    fprintf(fid,'Connectivity validated: 1. Kinetics validated: 0. Activation strength validated: 0.\n');
    fprintf(fid,'Name-matched BRENDA regulation; verify AC strength and [X] << K_M before building.\n');
end

% ---------- helpers ----------
function blk(fid, title, model, rxnNames, roles, r2i, roleMap)
    fprintf(fid,'\n------------------------------------------------------------\n%s\n',title);
    fprintf(fid,'------------------------------------------------------------\n\n');
    for q=1:numel(rxnNames)
        rn=rxnNames{q}; if isempty(rn)||strcmp(rn,'...')||~isKey(r2i,rn), continue; end
        role=''; if q<=numel(roles), role=roles{q}; end
        printRxn(fid, model, r2i(rn), role, roleMap);
    end
end
function printRxn(fid, model, j, role, roleMap)
    ec=''; if isfield(model,'rxnECNumbers'), ec=char(model.rxnECNumbers{j}); end
    nm=''; if isfield(model,'rxnNames'), nm=char(model.rxnNames{j}); end
    lb=-Inf; ub=Inf; if isfield(model,'lb'),lb=model.lb(j);end; if isfield(model,'ub'),ub=model.ub(j);end
    fprintf(fid,'  %s\n      ID  : %s\n      Name: %s\n      EC  : %s\n', role, model.rxns{j}, nm, ec);
    fprintf(fid,'      Bounds: [%g, %g]\n', lb, ub);
    fprintf(fid,'      IDs   : %s\n', eqnString(model, j, 'mets',     roleMap));
    fprintf(fid,'      Names : %s\n\n', eqnString(model, j, 'metNames', roleMap));
end
function s = eqnString(model, j, field, roleMap)
    sub=find(model.S(:,j)<0); pro=find(model.S(:,j)>0);
    f=@(ix) strjoin(arrayfun(@(i) coefName(model,i,abs(model.S(i,j)),field,roleMap), ix,'uni',0),' + ');
    s=sprintf('%s  ->  %s', f(sub), f(pro));
end
function t = coefName(model,i,c,field,roleMap)
    nm=char(model.(field){i});
    if ~isempty(roleMap) && isKey(roleMap,i), nm=sprintf('%s(%s)', nm, roleMap(i)); end
    if c==1, t=nm; else, t=sprintf('%g %s',c,nm); end
end
function ix = metIx(model, name)
    ix=find(strcmp(model.metNames,name),1); if isempty(ix), ix=find(strcmp(model.mets,name),1); end
end
function [rpath,mpath]=tracefwd(model, src, tgt, predRxn, predMet)
    rpath=[]; mpath=tgt; cur=tgt;
    while cur~=src
        rpath(end+1)=predRxn(cur); cur=predMet(cur); mpath(end+1)=cur; %#ok<AGROW>
    end
    rpath=fliplr(rpath); mpath=fliplr(mpath);
end
function p=splitr(s), p=strsplit(char(s),'>'); end
function s=strOr(x)
    if ismissing(x), s=''; return; end         % empty bias/removal reads as <missing>
    s=char(x);
    if strcmpi(s,'NaN')||isempty(s), s=''; end
end
