function out = map_iffl_circuit(model, T, opt)
% MAP_IFFL_CIRCUIT  Map the fixed IFFL circuit of iffl_topology.m onto a COBRA
% model.
%
%   out = map_iffl_circuit(model, T)            % T from iffl_coverage_check
%   out = map_iffl_circuit(model, T, opt)
%   map_iffl_circuit('selftest')                % synthetic model + planted decoys
%
% THE ONLY HARD RULE IS SIGNAL EXCLUSIVITY
%   Each circuit signal may be touched only by the circuit reactions its role
%   declares in iffl_topology.m: m_Ff_1 is produced by r_Ff_1 and consumed by
%   r_Fd_1 and by nothing else in the circuit, and so on for all 13 signals.
%   Catalytic regulators (m_Fd_1, m_Fd_2, m_c) must not appear in the
%   stoichiometry of anything. Everything else below is either structural (it
%   is in the drawn topology) or an ANNOTATION that is reported, not rejected.
%
% STRUCTURAL REQUIREMENTS (part of the diagram, not optional)
%   - two reactions per arm, four per differentiator, plus r_I, r_Fr, r_lk
%   - r_Fr consumes m_Fd_k  (without it the gate never clears)
%   - r_lk consumes m_d_k   (without it the output has no baseline)
%   - the counter's r_I1 and r_I2 draw on the same m_I
%   - m_Fd_k is a BRENDA activator of r_I_k; m_c is an activator of r_Nw
%
% SWITCHABLE FILTERS (opt.<name>, default in brackets)
%   carrierTest            [true]  reject an arm whose net stoichiometry does
%                                  not move mass input -> rail. A conserved
%                                  carrier (glutathione in LALDO+GLYOX) cancels
%                                  and transmits no signal. Keep this on.
%   compartment            ['c']   every circuit reaction has all participants
%                                  here. 'any' disables the restriction but
%                                  NEVER the compartment-qualified matching.
%   distinctEC             [false] reject when two roles share an EC
%   noRegulatoryCrosstalk  [false] reject when a signal gates a circuit
%                                  reaction other than its own target
%   Both of the last two are reported per circuit in out.circuits(k).notes
%   even when they are off, so you can sort instead of filter.
%
% OUTPUT  out.circuits  struct array: met.<role>, rxn.<role> (model indices),
%                       sig, rxnIdx, nGenes, polarity, notes
%         out.blocks    accepted differentiator blocks
%         out.depth     where each D1->D2 join died (see out.depthLabels)

    if nargin==1 && (ischar(model)||isstring(model)) && strcmpi(char(model),'selftest')
        out = selftest(); return;
    end
    if nargin < 3, opt = struct(); end
    opt  = applyDefaults(opt);
    topo = iffl_topology(opt.consumingN, opt.fastArmSteps, opt.miMode, ...
                         opt.catalyticInterstage, opt.gateMode, opt.neuronGate, ...
                         opt.interstageGate);
    % T may be a single activator table, or struct('AC',..,'IN',..) when a stage
    % is built on inhibition.
    T_IN = [];
    if isstruct(T) && isfield(T,'AC')
        if isfield(T,'IN'), T_IN = T.IN; end
        T = T.AC;
    end
    if ~isstruct(T), T = table2struct(T); end
    if ~isempty(T_IN) && ~isstruct(T_IN), T_IN = table2struct(T_IN); end
    needInhib = any(strcmp(topo.gateMode, 'inhibit')) || ...
                ~strcmp(topo.neuronGate, 'activate') || ...
                ~strcmp(topo.interGate, 'activate');
    assert(~needInhib || ~isempty(T_IN), ...
        ['gateMode asks for an inhibition stage but no IN seed table was given.\n' ...
         'Pass T as struct(''AC'', Tac, ''IN'', Tin).']);

    nR = numel(model.rxns); nM = numel(model.mets);
    cmp  = metCompartment(model.mets);
    curr = currencyMask(model, opt);
    sigOK = ~curr(:);
    if ~strcmpi(opt.compartment,'any'), sigOK = sigOK & strcmp(cmp(:), opt.compartment); end
    [usable, rxnWhy] = reactionUsable(model, cmp, opt);
    r2i = containers.Map(cellstr(model.rxns), num2cell(1:nR));

    ctx.prodBy = indexBy(model.S, usable, +1);
    ctx.consBy = indexBy(model.S, usable, -1);
    ctx.sigOK  = sigOK;
    [ctx.actOf, ctx.activates] = activationMaps(model, T, r2i, sigOK, nR, nM);
    if needInhib
        [ctx.inhOf, ctx.inhibits] = activationMaps(model, T_IN, r2i, sigOK, nR, nM);
    else
        ctx.inhOf = cell(nR,1); ctx.inhibits = cell(nM,1);
    end
    % which reactions m_c may drive, per neuronGate
    switch topo.neuronGate
        case 'inhibit', ctx.neuronOf = ctx.inhibits;
        case 'any'
            ctx.neuronOf = cell(nM,1);
            for q = 1:nM, ctx.neuronOf{q} = unique([ctx.activates{q} ctx.inhibits{q}]); end
        otherwise,      ctx.neuronOf = ctx.activates;
    end
    % which reactions m_d1 may drive as the interstage
    switch topo.interGate
        case 'inhibit', ctx.interOf = ctx.inhibits;
        case 'any'
            ctx.interOf = cell(nM,1);
            for q = 1:nM, ctx.interOf{q} = unique([ctx.activates{q} ctx.inhibits{q}]); end
        otherwise,      ctx.interOf = ctx.activates;
    end
    nRegMets = 0;
    for m = 1:nM
        if ~isempty(ctx.activates{m}) && any(usable(ctx.activates{m})), nRegMets = nRegMets + 1; end
    end

    depthLabels = { 'arms of D1 and D2 share a reaction', ...
                    'D1 and D2 share a signal metabolite', ...
                    'no free consumer of m_d2 for r_c', ...
                    'r_c makes no usable m_c', ...
                    'no free consumer of m_c for r_I1', ...
                    'no free second consumer of m_d2 for r_I2', ...
                    'r_I1 / r_I2 have no m_I partner (opt.requireMI)', ...
                    'm_c drives no free reaction (r_Nw; see opt.neuronGate)', ...
                    'r_Nw makes no usable m_out', ...
                    'no second producer of m_out (r_Nb)', ...
                    'no free consumer of m_out (r_Nr)', ...
                    'no clean r_Fr / r_lk left (opt.sinkMode = map)', ...
                    'signal exclusivity violated', ...
                    'optional filter (EC / regulatory)', ...
                    'accepted' };
    depth = zeros(1, numel(depthLabels));
    coll = containers.Map('KeyType','char','ValueType','double');
    viol = containers.Map('KeyType','char','ValueType','double');
    st = struct('seedsSeen',0,'seedsMissing',0,'anchorDrop',zeros(1,5),'blocksInhib',0, ...
                'pairsTotal',0,'pairsLinked',0, ...
                'anchors',0,'arm_carrier',0,'blocks',0,'joined',0, ...
                'gate_pairs',0,'gate_noFr',0,'gate_noArm',0, ...
                'out_pairs',0,'out_noLk',0,'out_noArm',0,'no_common_src',0);

    % ========================= differentiator blocks =========================
    cache2 = containers.Map('KeyType','double','ValueType','any');
    B = struct('mIn',{},'mFf',{},'mFd',{},'mSf',{},'mD',{}, ...
               'rFf',{},'rFd',{},'rSf',{},'rSd',{},'rI',{}, ...
               'frCand',{},'lkCand',{},'nGenes',{},'twoStep',{},'gate',{});
    for t = 1:numel(T)
        rn = strtrim(char(T(t).rxn));
        if ~isKey(r2i, rn)
            st.seedsMissing = st.seedsMissing + 1;   % name not in THIS model
            continue;
        end
        rI = r2i(rn);
        st.seedsSeen = st.seedsSeen + 1;
        if ~usable(rI)
            st.anchorDrop(rxnWhy(rI)) = st.anchorDrop(rxnWhy(rI)) + 1;
            continue;
        end
        st.anchors = st.anchors + 1;
        mDList = rowOf(find(model.S(:,rI) < 0), sigOK);
        for mFd = asRow(ctx.actOf{rI})
            st.gate_pairs = st.gate_pairs + 1;
            frCand = pickConsumers(ctx.consBy, mFd, rI, model, opt.maxHousekeeping);
            if isempty(frCand), st.gate_noFr = st.gate_noFr + 1; continue; end
            if opt.fastArmSteps == 2
                fastA = twoStepInto(mFd, ctx.prodBy, model.S, sigOK, opt.maxArmPaths, cache2);
            else
                fastA = oneStepInto(mFd, ctx.prodBy, model.S, sigOK);
            end
            if isempty(fastA), st.gate_noArm = st.gate_noArm + 1; continue; end
            for mD = asRow(mDList)
                if mD == mFd, continue; end
                st.out_pairs = st.out_pairs + 1;
                lkCand = pickConsumers(ctx.consBy, mD, rI, model, opt.maxHousekeeping);
                if isempty(lkCand), st.out_noLk = st.out_noLk + 1; continue; end
                slowA = twoStepInto(mD, ctx.prodBy, model.S, sigOK, opt.maxArmPaths, cache2);
                if isempty(slowA), st.out_noArm = st.out_noArm + 1; continue; end
                shared = asRow(intersect(fastA(:,4), slowA(:,4)));
                if isempty(shared), st.no_common_src = st.no_common_src + 1; end
                for si = 1:numel(shared)
                    src = shared(si);
                    if src==mFd || src==mD, continue; end
                    fr = fastA(fastA(:,4)==src, :);
                    sr = slowA(slowA(:,4)==src, :);
                    for a = 1:size(fr,1)
                        for b = 1:size(sr,1)
                            [blk, why] = tryBlock(model, src, fr(a,:), sr(b,:), rI, ...
                                                  mFd, mD, frCand, lkCand, opt);
                            if isempty(blk)
                                if strcmp(why,'carrier'), st.arm_carrier = st.arm_carrier+1; end
                                continue;
                            end
                            B(end+1) = blk; %#ok<AGROW>
                            if numel(B) >= opt.maxBlocks, break; end
                        end
                        if numel(B) >= opt.maxBlocks, break; end
                    end
                    if numel(B) >= opt.maxBlocks, break; end
                end
                if numel(B) >= opt.maxBlocks, break; end
            end
            if numel(B) >= opt.maxBlocks, break; end
        end
        if numel(B) >= opt.maxBlocks, break; end
    end
    B = dedupeBlocks(B);
    if ~isempty(B), [~,ord] = sort([B.nGenes]); B = B(ord); end
    if numel(B) > opt.topBlocks, B = B(1:opt.topBlocks); end
    BI = emptyBlock();
    if needInhib
        BI = buildInhibitionBlocks(model, T_IN, ctx, usable, r2i, sigOK, opt);
        BI = dedupeBlocks(BI);
        if ~isempty(BI), [~,o2] = sort([BI.nGenes]); BI = BI(o2); end
        if numel(BI) > opt.topBlocks, BI = BI(1:opt.topBlocks); end
    end
    st.blocks = numel(B); st.blocksInhib = numel(BI);
    vprint(opt,'differentiator blocks (activation): %d  (%d arms rejected as carriers)\n', ...
        numel(B), st.arm_carrier);
    if needInhib
        vprint(opt,'differentiator blocks (inhibition): %d\n', numel(BI));
    end
    if strcmp(topo.gateMode{1},'inhibit'), B1 = BI; else, B1 = B; end
    if strcmp(topo.gateMode{2},'inhibit'), B2 = BI; else, B2 = B; end
    samePool = strcmp(topo.gateMode{1}, topo.gateMode{2});
    if opt.verbose && isKey(cache2, 0) && cache2(0) > 0
        fprintf(['  WARNING: the 2-step arm enumeration hit opt.maxArmPaths (%d) on %d\n' ...
                 '  target(s). Those arms are TRUNCATED, so the block pool is arbitrarily\n' ...
                 '  thinned. Raise opt.maxArmPaths before concluding anything.\n'], ...
                 opt.maxArmPaths, cache2(0));
    end
    if opt.verbose && ~isempty(B), reportPool(model, B); end

    % ============================ assembly ===================================
    C = emptyCircuit();
    for i = 1:numel(B1)
        for j = 1:numel(B2)
            if samePool && i==j, continue; end
            st.pairsTotal = st.pairsTotal + 1;
            if topo.catalytic
                rA2s = catalyticLinks(model, ctx, B1(i), B2(j));
                if isempty(rA2s), continue; end           % m_d1 gates nothing making m_in_2
                st.pairsLinked = st.pairsLinked + 1;
            else
                if B2(j).mIn ~= B1(i).mD, continue; end   % m_d1 is block 2's input
                st.pairsLinked = st.pairsLinked + 1;
                rA2s = NaN;
            end
            st.joined = st.joined + 1;
            c = []; d = 1; why = '';
            for rA2 = asRow(rA2s)
                [c, d, why] = assemble(model, B1(i), B2(j), topo, opt, ctx, rA2);
                if ~isempty(c), break; end
            end
            depth(d) = depth(d) + 1;
            if ~isempty(why) && isempty(c)
                if strncmp(why,'EXCL: ',6)
                    k = why(7:end);
                    if isKey(viol, k), viol(k) = viol(k) + 1; else, viol(k) = 1; end
                else
                    if isKey(coll, why), coll(why) = coll(why) + 1; else, coll(why) = 1; end
                end
            end
            if ~isempty(c), C(end+1) = c; end %#ok<AGROW>
            if numel(C) >= opt.maxCircuits, break; end
        end
        if numel(C) >= opt.maxCircuits, break; end
    end
    if ~isempty(C), [~,ord] = sort([C.nGenes]); C = C(ord); end

    out.nRegMets = nRegMets;
    out.collKeys = keys(coll);
    if coll.Count > 0, out.collVals = asRow(cell2mat(values(coll))); else, out.collVals = []; end
    out.violKeys = keys(viol);
    if viol.Count > 0, out.violVals = asRow(cell2mat(values(viol))); else, out.violVals = []; end
    out.armCapHits = 0; if isKey(cache2,0), out.armCapHits = cache2(0); end
    out.poolA = poolSpecies(B); out.poolI = poolSpecies(BI);
    out.circuits = C; out.blocks = B; out.blocksInhib = BI; out.stats = st;
    out.depth = depth; out.depthLabels = depthLabels;
    out.topology = topo; out.opt = opt;
    vprint(opt,'circuits with exclusive signals: %d\n', numel(C));
    if opt.verbose, printDepth(depth, depthLabels, st); printColl(coll); end
end

% ======================== one differentiator block ==========================
function [blk, why] = tryBlock(model, src, fastRow, slowRow, rI, mFd, mD, frCand, lkCand, opt)
    blk = []; why = '';
    rFf = fastRow(1); mFf = fastRow(2); rFd = fastRow(3);
    rSf = slowRow(1); mSf = slowRow(2); rSd = slowRow(3);
    twoStep = ~isnan(rFf);
    if twoStep
        mets = [src mFf mFd mSf mD];
        rxns = [rFf rFd rSf rSd rI];
        localAllowed = { [1 2], [2 3], [1 4], [4 5], 5 };
    else
        mets = [src mFd mSf mD];
        rxns = [rFd rSf rSd rI];
        localAllowed = { [1 2], [1 3], [3 4], 4 };
    end
    if numel(unique(mets)) ~= numel(mets), why='mets'; return; end
    if numel(unique(rxns)) ~= numel(rxns), why='rxns'; return; end
    if opt.carrierTest
        if twoStep && ~armCarriesMass(model.S, rFf, rFd, src, mFd), why='carrier'; return; end
        if ~armCarriesMass(model.S, rSf, rSd, src, mD),  why='carrier'; return; end
    end
    % block-local signal exclusivity (cheap prefilter for the global test)
    for q = 1:numel(rxns)
        p = find(model.S(:,rxns(q)) ~= 0);
        if ~isempty(setdiff(intersect(p, mets), mets(localAllowed{q}))), why='local'; return; end
    end
    blk.mIn=src; blk.mFf=mFf; blk.mFd=mFd; blk.mSf=mSf; blk.mD=mD;
    blk.rFf=rFf; blk.rFd=rFd; blk.rSf=rSf; blk.rSd=rSd; blk.rI=rI;
    blk.frCand = frCand;        % r_Fr candidates, resolved at assembly time
    blk.lkCand = lkCand;        % r_lk candidates, resolved at assembly time
    blk.nGenes = unionGenes(model, rxns);
    blk.twoStep = twoStep;
    blk.gate = 'activate';
end

function ok = armCarriesMass(S, r1, r2, srcIx, railIx)
    net = S(:,r1) + S(:,r2);
    ok = net(srcIx) < 0 && net(railIx) > 0;
end

% ============================== assembly ====================================
function [circ, depth, why] = assemble(model, b1, b2, topo, opt, ctx, rA2)
% Housekeeping (r_Fr, r_lk) is chosen LAST, once every other reaction is fixed:
% r_lk_1 and r_Ff_2 / r_Sf_2 are all consumers of m_d1, so choosing r_lk_1
% inside block 1 made it collide with block 2 and silently killed the join.
    circ = []; depth = 1; why = '';
    % Role vectors are derived from the topology, never hard-coded: with
    % fastArmSteps, catalyticInterstage, miMode and gateMode all variable, an
    % enumerated vector per combination is how ordering bugs get in.
    blocks = {b1, b2};
    arms = []; armNames = {};
    for k = 1:2
        b = blocks{k}; sfx = sprintf('_%d', k);
        if topo.fastSteps == 2
            arms(end+1) = b.rFf; armNames{end+1} = ['r_Ff' sfx]; %#ok<AGROW>
        end
        arms(end+1) = b.rFd; armNames{end+1} = ['r_Fd' sfx]; %#ok<AGROW>
        arms(end+1) = b.rSf; armNames{end+1} = ['r_Sf' sfx]; %#ok<AGROW>
        arms(end+1) = b.rSd; armNames{end+1} = ['r_Sd' sfx]; %#ok<AGROW>
        if ~isnan(b.rI)                      % absent on an inhibition stage
            arms(end+1) = b.rI; armNames{end+1} = ['r_I' sfx]; %#ok<AGROW>
        end
    end
    if topo.catalytic
        arms(end+1) = rA2; armNames{end+1} = 'r_A_2';
    end
    base = b1.mIn; baseNames = {'m_in'};
    for k = 1:2
        b = blocks{k}; sfx = sprintf('_%d', k);
        if k == 2 && topo.catalytic
            base(end+1) = b.mIn; baseNames{end+1} = 'm_in_2'; %#ok<AGROW>
        end
        if topo.fastSteps == 2
            base(end+1) = b.mFf; baseNames{end+1} = ['m_Ff' sfx]; %#ok<AGROW>
        end
        base(end+1) = b.mFd; baseNames{end+1} = ['m_Fd' sfx]; %#ok<AGROW>
        base(end+1) = b.mSf; baseNames{end+1} = ['m_Sf' sfx]; %#ok<AGROW>
        base(end+1) = b.mD;  baseNames{end+1} = sprintf('m_d%d', k); %#ok<AGROW>
    end
    why = firstDup(arms, armNames);
    if ~isempty(why), return; end
    depth = 2;
    why = firstDup(base, baseNames);
    if ~isempty(why), return; end

    consBy = ctx.consBy; prodBy = ctx.prodBy; sigOK = ctx.sigOK;
    for rc = asRow(freeOf(consBy{b2.mD}, arms))
        depth = max(depth,4);
        for mC = asRow(rowOf(find(model.S(:,rc) > 0), sigOK))
            if any(mC == base), continue; end
            % m_c must not be touched by any reaction already fixed. Without
            % this the join is built and then fails exclusivity ("r_I_2
            % touches m_c"), wasting every remaining combination.
            if any(model.S(mC, arms) ~= 0), continue; end
            depth = max(depth,5);
            for rI1 = asRow(freeOf(consBy{mC}, [arms rc]))
                depth = max(depth,6);
                for rI2 = asRow(freeOf(consBy{b2.mD}, [arms rc rI1]))
                    depth = max(depth,7);
                    switch topo.miMode
                        case 'none'
                            mIs = NaN;                    % one pass, no m_I role
                        case 'shared'
                            mIs = intersect(find(model.S(:,rI1)<0), find(model.S(:,rI2)<0));
                            mIs = rowOf(mIs, sigOK); mIs = setdiff(mIs, [base mC]);
                        otherwise
                            mIs = find(model.S(:,rI1)<0); % r_I2 keeps its own partner
                            mIs = rowOf(mIs, sigOK); mIs = setdiff(mIs, [base mC]);
                    end
                    if isempty(mIs), continue; end
                    depth = max(depth,8);
                    for mI = asRow(mIs)
                        for rNw = asRow(freeOf(ctx.neuronOf{mC}, [arms rc rI1 rI2]))
                            if ~opt.consumingN && model.S(mC, rNw) ~= 0, continue; end
                            depth = max(depth,9);
                            for mOut = asRow(rowOf(find(model.S(:,rNw) > 0), sigOK))
                                if any(mOut == [base mC mI]), continue; end
                                depth = max(depth,10);
                                for rNb = asRow(freeOf(prodBy{mOut}, [arms rc rI1 rI2 rNw]))
                                    mbs = rowOf(find(model.S(:,rNb) < 0), sigOK);
                                    mbs = setdiff(mbs, [base mC mI mOut]);
                                    if isempty(mbs), continue; end
                                    depth = max(depth,11);
                                    for rNr = asRow(freeOf(consBy{mOut}, [arms rc rI1 rI2 rNw rNb]))
                                        depth = max(depth,12);
                                        for mB = asRow(mbs)
                                            sig = buildSig(topo, b1, b2, mC, mI, mB, mOut);
                                            if numel(unique(sig)) ~= numel(sig), continue; end
                                            fixed = [arms rc rI1 rI2 rNw rNb rNr];
                                            hk = pickHousekeeping(model, b1, b2, sig, fixed, opt.sinkMode);
                                            if isempty(hk), continue; end
                                            depth = max(depth,13);
                                            rxn = buildRxn(topo, b1, b2, hk, rA2, rc, rI1, rI2, rNw, rNb, rNr);
                                            [okx, vx] = exclusiveSignals(model, sig, rxn, topo);
                                            if ~okx
                                                if isempty(why), why = ['EXCL: ' vx]; end
                                                continue;
                                            end
                                            depth = max(depth,14);
                                            notes = annotate(model, sig, rxn, topo, ctx, opt);
                                            if opt.noRegulatoryCrosstalk && notes.nRegCollisions > 0, continue; end
                                            if opt.distinctEC && notes.nSharedEC > 0, continue; end
                                            depth = 15; why = '';
                                            if topo.catalytic
                                                isens = senseOf(ctx, b1.mD, rA2);
                                            else
                                                isens = 'n/a (mass-flow link)';
                                            end
                                            circ = packCircuit(model, topo, sig, rxn, notes, ...
                                                senseOf(ctx, mC, rNw), isens);
                                            return;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function sig = buildSig(topo, b1, b2, mC, mI, mB, mOut)
% Signal vector IN TOPOLOGY ORDER, from a role map rather than an enumerated
% list: with fastArmSteps, catalyticInterstage, miMode and gateMode all
% variable, a hand-written vector per combination is how ordering bugs get in.
    M = containers.Map('KeyType','char','ValueType','double');
    M('m_in') = b1.mIn;
    blocks = {b1, b2};
    for k = 1:2
        b = blocks{k}; sfx = sprintf('_%d', k);
        if k == 2 && topo.catalytic, M('m_in_2') = b.mIn; end
        if topo.fastSteps == 2, M(['m_Ff' sfx]) = b.mFf; end
        M(['m_Fd' sfx]) = b.mFd;
        M(['m_Sf' sfx]) = b.mSf;
        M(sprintf('m_d%d', k)) = b.mD;
    end
    M('m_c') = mC; M('m_b') = mB; M('m_out') = mOut;
    if ~strcmp(topo.miMode,'none') && ~isnan(mI), M('m_I') = mI; end
    sig = zeros(1, numel(topo.metRoles));
    for q = 1:numel(topo.metRoles), sig(q) = M(topo.metRoles{q}); end
end

function rxn = buildRxn(topo, b1, b2, hk, rA2, rc, rI1, rI2, rNw, rNb, rNr)
% Reaction vector IN TOPOLOGY ORDER. hk = [r_Fr_1 r_lk_1 r_Fr_2 r_lk_2]; an
% entry may be NaN, meaning an added sink rather than a mapped reaction.
    R = containers.Map('KeyType','char','ValueType','double');
    blocks = {b1, b2};
    for k = 1:2
        b = blocks{k}; sfx = sprintf('_%d', k);
        if topo.fastSteps == 2, R(['r_Ff' sfx]) = b.rFf; end
        R(['r_Fd' sfx]) = b.rFd;
        R(['r_Sf' sfx]) = b.rSf;
        R(['r_Sd' sfx]) = b.rSd;
        if ~isnan(b.rI), R(['r_I' sfx]) = b.rI; end
        R(['r_Fr' sfx]) = hk(2*k-1);
        R(['r_lk' sfx]) = hk(2*k);
    end
    if topo.catalytic, R('r_A_2') = rA2; end
    R('r_c') = rc; R('r_I1') = rI1; R('r_I2') = rI2;
    R('r_Nw') = rNw; R('r_Nb') = rNb; R('r_Nr') = rNr;
    rxn = zeros(1, numel(topo.rxnRoles));
    for q = 1:numel(topo.rxnRoles), rxn(q) = R(topo.rxnRoles{q}); end
end

function s = firstDup(v, names)
% which two roles got the same model object - the histogram of this is what
% tells you whether the search space or the policy is the problem
    s = '';
    for a = 1:numel(v)-1
        for b = a+1:numel(v)
            if v(a) == v(b), s = sprintf('%s = %s', names{a}, names{b}); return; end
        end
    end
end
function B = buildInhibitionBlocks(model, T, ctx, usable, r2i, sigOK, opt)
% An inhibition stage anchors on the reaction that PRODUCES m_d and is
% INHIBITED by the gate - the mirror image of the activation anchor, which
% consumes m_d and is activated by the gate. There is no r_I: m_d falls
% because its production stops, and r_lk drains it. m_Fd is still cleared by
% r_Fr, so the gate returns to baseline and the response is a pulse.
    B = emptyBlock();
    cache2 = containers.Map('KeyType','double','ValueType','any');
    for t = 1:numel(T)
        rn = strtrim(char(T(t).rxn));
        if ~isKey(r2i, rn), continue; end
        rSd = r2i(rn);
        if ~usable(rSd), continue; end
        mDList  = asRow(rowOf(find(model.S(:,rSd) > 0), sigOK));   % products
        mSfList = asRow(rowOf(find(model.S(:,rSd) < 0), sigOK));   % substrates
        if isempty(mDList) || isempty(mSfList), continue; end
        for mFd = asRow(ctx.inhOf{rSd})
            frCand = pickConsumers(ctx.consBy, mFd, rSd, model, opt.maxHousekeeping);
            if isempty(frCand), continue; end
            if opt.fastArmSteps == 2
                fastA = twoStepInto(mFd, ctx.prodBy, model.S, sigOK, opt.maxArmPaths, cache2);
            else
                fastA = oneStepInto(mFd, ctx.prodBy, model.S, sigOK);
            end
            if isempty(fastA), continue; end
            for mD = mDList
                if mD == mFd, continue; end
                lkCand = pickConsumers(ctx.consBy, mD, rSd, model, opt.maxHousekeeping);
                if isempty(lkCand), continue; end
                for mSf = mSfList
                    if mSf == mD || mSf == mFd, continue; end
                    slowA = oneStepInto(mSf, ctx.prodBy, model.S, sigOK);
                    if isempty(slowA), continue; end
                    for src = asRow(intersect(fastA(:,4), slowA(:,4)))
                        if src==mFd || src==mD || src==mSf, continue; end
                        fr = fastA(fastA(:,4)==src, :);
                        sr = slowA(slowA(:,4)==src, :);
                        for a = 1:size(fr,1)
                            for b = 1:size(sr,1)
                                blk = tryInhibBlock(model, src, fr(a,:), sr(b,1:3), ...
                                                    rSd, mFd, mSf, mD, frCand, lkCand, opt);
                                if isempty(blk), continue; end
                                B(end+1) = blk; %#ok<AGROW>
                                if numel(B) >= opt.maxBlocks, return; end
                            end
                        end
                    end
                end
            end
        end
    end
end

function blk = tryInhibBlock(model, src, fastRow, slowRow, rSd, mFd, mSf, mD, frCand, lkCand, opt)
    blk = [];
    rFf = fastRow(1); mFf = fastRow(2); rFd = fastRow(3);
    rSf = slowRow(3);
    twoStep = ~isnan(rFf);
    if twoStep, mets = [src mFf mFd mSf mD]; else, mets = [src mFd mSf mD]; end
    if numel(unique(mets)) ~= numel(mets), return; end
    if twoStep, rxns = [rFf rFd rSf rSd]; else, rxns = [rFd rSf rSd]; end
    if numel(unique(rxns)) ~= numel(rxns), return; end
    if opt.carrierTest
        if twoStep && ~armCarriesMass(model.S, rFf, rFd, src, mFd), return; end
        if ~armCarriesMass(model.S, rSf, rSd, src, mD), return; end
    end
    % block-local signal exclusivity, same rule as the activation builder
    if twoStep
        localAllowed = { [1 2], [2 3], [1 4], [4 5] };   % rFf rFd rSf rSd
    else
        localAllowed = { [1 2], [1 3], [3 4] };          % rFd rSf rSd
    end
    for q = 1:numel(rxns)
        p = find(model.S(:,rxns(q)) ~= 0);
        if ~isempty(setdiff(intersect(p, mets), mets(localAllowed{q}))), return; end
    end
    blk.mIn=src; blk.mFf=mFf; blk.mFd=mFd; blk.mSf=mSf; blk.mD=mD;
    blk.rFf=rFf; blk.rFd=rFd; blk.rSf=rSf; blk.rSd=rSd; blk.rI=NaN;
    blk.frCand = frCand; blk.lkCand = lkCand;
    blk.nGenes = unionGenes(model, rxns);
    blk.twoStep = twoStep; blk.gate = 'inhibit';
end

function p = poolSpecies(B)
% the distinct species a block pool can offer in each role - the number that
% actually limits how many disjoint pairs exist
    if isempty(B)
        p = struct('gates',[],'outs',[],'ins',[],'n',0); return;
    end
    p.gates = unique([B.mFd]);
    p.outs  = unique([B.mD]);
    p.ins   = unique([B.mIn]);
    p.n     = numel(B);
end

function B = emptyBlock()
    B = struct('mIn',{},'mFf',{},'mFd',{},'mSf',{},'mD',{}, ...
               'rFf',{},'rFd',{},'rSf',{},'rSd',{},'rI',{}, ...
               'frCand',{},'lkCand',{},'nGenes',{},'twoStep',{},'gate',{});
end

function js = catalyticLinks(model, ctx, b1, b2)
% Reactions that m_d1 ACTIVATES and that PRODUCE m_in_2.
%
% Two conditions, both load-bearing:
%   1. the link must not touch m_d1 - the interstage READS D1's output and
%      never draws on it, which is what makes it high-impedance;
%   2. among all the signals of both stages it may touch ONLY m_in_2, because
%      that is r_A_2's entire allowance in iffl_topology. Checking this here
%      rather than at the exclusivity test lets the search try another link
%      instead of discarding the whole join.
    js = [];
    sigSet = blockSignals(b1);
    s2 = blockSignals(b2);
    sigSet = [sigSet s2];
    for j = asRow(ctx.interOf{b1.mD})
        if model.S(b2.mIn, j) <= 0, continue; end
        if model.S(b1.mD, j) ~= 0, continue; end
        touched = intersect(find(model.S(:,j) ~= 0), sigSet);
        if ~isempty(setdiff(touched, b2.mIn)), continue; end
        js(end+1) = j; %#ok<AGROW>
    end
end

function v = blockSignals(b)
    v = [b.mIn b.mFd b.mSf b.mD];
    if ~isnan(b.mFf), v = [v b.mFf]; end
    v = unique(v);
end
function hk = pickHousekeeping(model, b1, b2, sig, fixed, sinkMode)
% r_Fr_1, r_lk_1, r_Fr_2, r_lk_2 - each must touch, among the 13 signals, only
% the one it is supposed to consume. Candidates are individually clean, so any
% four distinct survivors work.
% The four lists are pairwise DISJOINT: cleanCands keeps only reactions that
% touch exactly one signal (their own target), and the four targets m_Fd_1,
% m_d1, m_Fd_2, m_d2 are distinct. So no assignment conflict is possible and a
% greedy pick is already optimal - a failure here means some list is EMPTY.
%
% An empty list is not necessarily a dead end. In the diagram r_Fr and r_lk
% both discharge to the empty set: they are REMOVAL steps, not conversions to
% something the cell needs. In a cell-free reactor you supply them yourself -
% outflow in a CSTR, a protease, a phosphatase - so they need not correspond to
% an existing metabolic reaction at all. sinkMode controls that:
%   'map'    every removal must be a clean model reaction (strictest)
%   'prefer' map one when a clean reaction exists, otherwise declare a sink
%   'sink'   always declare added sinks, never map them
    hk = [];
    lists = { cleanCands(model, b1.frCand, sig, b1.mFd), ...
              cleanCands(model, b1.lkCand, sig, b1.mD ), ...
              cleanCands(model, b2.frCand, sig, b2.mFd), ...
              cleanCands(model, b2.lkCand, sig, b2.mD ) };
    used = fixed; pick = nan(1,4);
    for q = 1:4
        if strcmp(sinkMode,'sink'), continue; end        % leave as NaN = added sink
        p = freeOf(lists{q}, used);
        if isempty(p)
            if strcmp(sinkMode,'map'), return; end       % strict: no clean reaction
            continue;                                    % 'prefer': fall back to a sink
        end
        pick(q) = p(1); used = [used pick(q)]; %#ok<AGROW>
    end
    hk = pick;
end
function c = cleanCands(model, cand, sig, ownMet)
    keep = false(1, numel(cand));
    for q = 1:numel(cand)
        p = find(model.S(:, cand(q)) ~= 0);
        keep(q) = isempty(setdiff(intersect(p, sig), ownMet));
    end
    c = cand(keep);
end

function [ok, viol] = exclusiveSignals(model, sig, rxn, topo)
% THE rule: a circuit signal is touched only by the reactions of its own roles,
% with the declared sign. Everything the reactions do to NON-signal species is
% left alone.
    ok = false; viol = '';
    r = rxn(~isnan(rxn));
    if numel(unique(r)) ~= numel(r), viol = 'duplicate reaction'; return; end
    for j = 1:numel(rxn)
        if isnan(rxn(j)), continue; end              % added sink: nothing to check
        col = model.S(:, rxn(j));
        bad = setdiff(intersect(find(col ~= 0), sig), sig(topo.allowed{j}));
        if ~isempty(bad)
            viol = sprintf('%s touches %s', topo.rxnRoles{j}, ...
                topo.metRoles{find(sig == bad(1), 1)});
            return;
        end
        req = topo.subs{j};
        if isfield(topo,'reqSubs'), req = topo.reqSubs{j}; end
        for q = 1:numel(req)
            if col(sig(strcmp(topo.metRoles, req{q}))) >= 0
                viol = sprintf('%s does not consume %s', topo.rxnRoles{j}, req{q}); return;
            end
        end
        for q = 1:numel(topo.prods{j})
            if col(sig(strcmp(topo.metRoles, topo.prods{j}{q}))) <= 0
                viol = sprintf('%s does not produce %s', topo.rxnRoles{j}, topo.prods{j}{q}); return;
            end
        end
    end
    ok = true;
end

function notes = annotate(model, sig, rxn, topo, ctx, opt)
% Reported, never rejected unless the matching opt flag is on.
    notes.regCollisions = {};
    for j = 1:numel(rxn)
        if isnan(rxn(j)), continue; end
        acts = intersect(ctx.actOf{rxn(j)}, sig);
        k = find(strcmp(topo.regTarget, topo.rxnRoles{j}), 1);
        if isempty(k), allowedAct = [];
        else, allowedAct = sig(strcmp(topo.metRoles, topo.regulator{k})); end
        for m = asRow(setdiff(acts, allowedAct))
            notes.regCollisions{end+1} = sprintf('%s also gates %s', ...
                topo.metRoles{sig==m}, topo.rxnRoles{j}); %#ok<AGROW>
        end
    end
    notes.nRegCollisions = numel(notes.regCollisions);
    e = cell(1, numel(rxn));
    for q = 1:numel(rxn)
        if ~isnan(rxn(q)) && isfield(model,'rxnECNumbers')
            e{q} = strtrim(char(model.rxnECNumbers{rxn(q)}));
        else, e{q} = ''; end
    end
    notes.sharedEC = {}; present = e(~cellfun(@isempty, e));
    u = unique(present);
    for q = 1:numel(u)
        k = find(strcmp(e, u{q}));
        if numel(k) > 1
            notes.sharedEC{end+1} = sprintf('%s shared by %s', u{q}, ...
                strjoin(topo.rxnRoles(k), ', ')); %#ok<AGROW>
        end
    end
    notes.nSharedEC = numel(notes.sharedEC);
    base = lower(regexprep(model.mets, '\[[a-z0-9]+\]$', ''));
    isCurr = ismember(base, lower(opt.currency));
    seen = []; shared = [];
    for q = 1:numel(rxn)
        if isnan(rxn(q)), continue; end
        p = find(model.S(:, rxn(q)) ~= 0)';
        p = p(~isCurr(p)); p = setdiff(p, sig);
        shared = [shared asRow(intersect(p, seen))]; %#ok<AGROW>
        seen = asRow(union(seen, p));
    end
    notes.sharedAux = unique(shared);
    % how many reactions OUTSIDE the circuit also touch each signal (only
    % matters in a whole cell; in a cell-free pot those enzymes are absent)
    notes.competitors = zeros(1, numel(sig));
    for q = 1:numel(sig)
        notes.competitors(q) = nnz(model.S(sig(q), :) ~= 0) - ...
            nnz(ismember(rxn(~isnan(rxn)), find(model.S(sig(q),:) ~= 0)));
    end
end

function c = packCircuit(model, topo, sig, rxn, notes, neuronSense, interSense)
    c.met = struct(); c.rxn = struct();
    for q = 1:numel(topo.metRoles), c.met.(topo.metRoles{q}) = sig(q); end
    for q = 1:numel(topo.rxnRoles), c.rxn.(topo.rxnRoles{q}) = rxn(q); end
    c.sig = sig; c.rxnIdx = rxn;
    c.nGenes = unionGenes(model, rxn);
    c.nSinks = sum(isnan(rxn));
    c.neuronSense = neuronSense;
    c.interSense  = interSense;
    c.polarity = topo.stagePolarity;
    % an inhibiting interstage is an extra sign flip between the stages
    % strncmp, NOT strfind: the ambiguous label
    % 'activates (also annotated as an inhibitor)' contains the substring
    % 'inhib' but describes an ACTIVATING edge, so a substring test flips the
    % sign of a perfectly good circuit.
    c.interGain = 1;
    if strncmp(interSense, 'inhibit', 7), c.interGain = -1; end
    c.netPolarity = topo.netPolarity * c.interGain;
    c.notes = notes;
end
function C = emptyCircuit()
    C = struct('met',{},'rxn',{},'sig',{},'rxnIdx',{},'nGenes',{},'nSinks',{}, ...
               'neuronSense',{},'interSense',{},'polarity',{},'interGain',{}, ...
               'netPolarity',{},'notes',{});
end

function s = senseOf(ctx, mC, rNw)
% Which regulatory relation actually links m_c to r_Nw in THIS circuit. With
% neuronGate = 'any' the policy allows either, so the report must state what
% was found rather than what was permitted.
    isAct = any(ctx.actOf{rNw} == mC);
    isInh = ~isempty(ctx.inhOf{rNw}) && any(ctx.inhOf{rNw} == mC);
    if isAct && ~isInh,     s = 'activates';
    elseif isInh && ~isAct, s = 'inhibits';
    elseif isAct && isInh,  s = 'activates (also annotated as an inhibitor)';
    else,                   s = 'regulates';
    end
end

% ============================== helpers =====================================
function v = asRow(x)
% MATLAB and Octave disagree about the orientation set-functions return for
% vector inputs, and `for x = columnVector` binds the WHOLE column to x, which
% then blows up at the next `||`. full() because a sparse scalar is also not a
% valid short-circuit operand. Every loop variable in this file comes through
% here, so neither failure can recur.
    v = reshape(full(x), 1, []);
end
function v = rowOf(ix, mask)
    ix = asRow(ix);
    if nargin > 1 && ~isempty(ix)
        ix = ix(logical(asRow(mask(ix))));
    end
    v = asRow(ix);
end
function v = freeOf(list, used)
    v = asRow(setdiff(asRow(list), asRow(used)));
end

function A = oneStepInto(target, prodBy, S, sigOK)
% rows [NaN NaN rDown src] : src --rDown--> target  (one-reaction fast arm)
    A = zeros(0,4);
    for rD = asRow(prodBy{target})
        for s = asRow(rowOf(find(S(:,rD) < 0), sigOK))
            if s == target, continue; end
            A(end+1,:) = [NaN NaN rD s]; %#ok<AGROW>
        end
    end
end
function A = twoStepInto(target, prodBy, S, sigOK, cap, cache)
% rows [rUp mMid rDown src] : src --rUp--> mMid --rDown--> target
    if isKey(cache, target), A = cache(target); return; end
    A = zeros(0,4);
    for rD = asRow(prodBy{target})
        for mMid = asRow(rowOf(find(S(:,rD) < 0), sigOK))
            if mMid == target, continue; end
            for rU = asRow(prodBy{mMid})
                if rU == rD, continue; end
                for s = asRow(rowOf(find(S(:,rU) < 0), sigOK))
                    if s==mMid || s==target, continue; end
                    A(end+1,:) = [rU mMid rD s]; %#ok<AGROW>
                    if size(A,1) >= cap
                        if isKey(cache,0), cache(0) = cache(0)+1; else, cache(0) = 1; end
                        cache(target) = A; return;
                    end
                end
            end
        end
    end
    cache(target) = A;
end

function cand = pickConsumers(consBy, metIx, excludeRxn, model, maxN)
    js = setdiff(consBy{metIx}, excludeRxn);
    if isempty(js), cand = []; return; end
    g = zeros(1,numel(js));
    for q = 1:numel(js), g(q) = geneCount(model, js(q)); end
    [~,ord] = sort(g); js = js(ord);
    cand = js(1:min(maxN, numel(js)));
    cand = cand(:)';
end

function P = indexBy(S, usable, sgn)
    nM = size(S,1); P = cell(nM,1);
    for m = 1:nM
        if sgn > 0, js = find(S(m,:) > 0); else, js = find(S(m,:) < 0); end
        P{m} = js(usable(js));
    end
end

function [actOf, activates] = activationMaps(model, T, r2i, sigOK, nR, nM)
    actOf = cell(nR,1); activates = cell(nM,1);
    for q = 1:nR, actOf{q} = []; end
    for q = 1:nM, activates{q} = []; end
    name2met = containers.Map('KeyType','char','ValueType','any');
    for i = 1:nM
        if ~sigOK(i), continue; end
        k = normLigand(model.metNames{i});
        if isKey(name2met,k), name2met(k) = [name2met(k) i]; else, name2met(k) = i; end
    end
    for t = 1:numel(T)
        rn = strtrim(char(T(t).rxn));
        if ~isKey(r2i, rn), continue; end
        j = r2i(rn);
        regs = strsplit(char(T(t).regulators), ';');
        for g = 1:numel(regs)
            k = normLigand(regs{g});
            if isempty(k) || ~isKey(name2met,k), continue; end
            ix = name2met(k);
            actOf{j} = asRow(unique([actOf{j} asRow(ix)]));
            for m = asRow(ix), activates{m} = asRow(unique([activates{m} j])); end
        end
    end
end

function [usable, why] = reactionUsable(model, cmp, opt)
% why(j) = 0 usable, 1 not forward, 2 boundary (exchange/demand/sink),
%          3 not entirely in opt.compartment, 4 transport (species unchanged),
%          5 no EC number
    nR = numel(model.rxns);
    usable = false(nR,1); why = zeros(nR,1);
    ub = ones(nR,1); if isfield(model,'ub'), ub = model.ub(:); end
    base = regexprep(model.mets, '\[[a-z0-9]+\]$', '');
    anyCmp = strcmpi(opt.compartment, 'any');
    for j = 1:nR
        if ub(j) <= 0, why(j) = 1; continue; end
        if opt.requireEC && isfield(model,'rxnECNumbers') && ...
           isempty(strtrim(char(model.rxnECNumbers{j}))), why(j) = 5; continue; end
        p = find(model.S(:,j) ~= 0);
        if numel(p) < 2, why(j) = 2; continue; end            % exchange / demand / sink
        if ~anyCmp && ~all(strcmp(cmp(p), opt.compartment)), why(j) = 3; continue; end
        subs = base(model.S(:,j) < 0); pros = base(model.S(:,j) > 0);
        if ~isempty(intersect(subs, pros)), why(j) = 4; continue; end  % species unchanged
        usable(j) = true;
    end
end

function curr = currencyMask(model, opt)
    base = lower(regexprep(model.mets, '\[[a-z0-9]+\]$', ''));
    curr = ismember(base, lower(opt.currency));
    nm = lower(model.metNames);
    for p = 1:numel(opt.currencyNamePatterns)
        pat = lower(opt.currencyNamePatterns{p});
        for i = 1:numel(nm)
            if ~curr(i) && ~isempty(strfind(nm{i}, pat)), curr(i) = true; end %#ok<STREMP>
        end
    end
end

function cmp = metCompartment(mets)
    cmp = regexp(mets, '\[([a-z0-9]+)\]$', 'tokens', 'once');
    cmp = cellfun(@(c) char([c{:}]), cmp, 'uni', 0);
end

function B = dedupeBlocks(B)
    if isempty(B), return; end
    key = zeros(numel(B), 8);
    for q = 1:numel(B)
        key(q,:) = [fz(B(q).rFf) B(q).rFd B(q).rSf B(q).rSd B(q).rI B(q).mFd B(q).mD B(q).mIn];
    end
    [~, ia] = unique(key, 'rows', 'first');
    B = B(sort(ia));
end

function v = fz(x), if isnan(x), v = 0; else, v = x; end, end
function n = unionGenes(model, rxnIdx)
    rxnIdx = rxnIdx(~isnan(rxnIdx));
    if isfield(model,'rxnGeneMat') && ~isempty(model.rxnGeneMat)
        n = nnz(any(model.rxnGeneMat(rxnIdx,:) ~= 0, 1));
    else
        n = numel(unique(rxnIdx));
    end
end
function g = geneCount(model, j)
    if isfield(model,'rxnGeneMat') && ~isempty(model.rxnGeneMat)
        g = nnz(model.rxnGeneMat(j,:));
    else, g = 1; end
end
function s = normLigand(s)
    s = lower(strtrim(char(s)));
    s = regexprep(s, '^((l|d|dl|\(s\)|\(r\)|alpha|beta)-)+', '');
    s = regexprep(s, '[^a-z0-9]+', ' ');
    s = strtrim(s);
end
function vprint(opt, varargin)
    if opt.verbose, fprintf(varargin{:}); end
end
function printDepth(depth, labels, st)
    fprintf('\nD1->D2 joins tried: %d. Where each one stopped:\n', st.joined);
    for q = 1:numel(labels)
        if depth(q) > 0
            fprintf('  %6d  %s\n', depth(q), labels{q});
        end
    end
    if depth(end) == 0
        fprintf(['The deepest stage reached is the one to attack. Signal exclusivity is\n' ...
                 'the only hard rule; if the block is upstream of that, it is the search\n' ...
                 'space, not the crosstalk policy.\n']);
    end
end

function reportPool(model, B)
    gates = unique([B.mFd]); outs = unique([B.mD]); ins = unique([B.mIn]);
    fprintf('  distinct m_Fd (gates) %d | distinct m_d (outputs) %d | distinct m_in %d\n', ...
        numel(gates), numel(outs), numel(ins));
    k = zeros(1,numel(gates));
    for q = 1:numel(gates), k(q) = sum([B.mFd] == gates(q)); end
    [k,ord] = sort(k,'descend'); gates = gates(ord);
    fprintf('  blocks per gate (top 5):');
    for q = 1:min(5,numel(gates))
        fprintf(' %s x%d;', model.mets{gates(q)}, k(q));
    end
    fprintf('\n');
    % how many blocks END at each metabolite vs START from it: a hinge needs both
    hinge = intersect(outs, ins);
    fprintf('  usable hinges (a metabolite that is both some block''s m_d and another''s m_in): %d\n', ...
        numel(hinge));
end
function printColl(coll)
    if coll.Count == 0, return; end
    ks = keys(coll); vs = cell2mat(values(coll));
    [vs,ord] = sort(vs,'descend'); ks = ks(ord);
    fprintf('\nWhich two roles claimed the same object (top %d):\n', min(10,numel(ks)));
    for q = 1:min(10,numel(ks)), fprintf('  %6d  %s\n', vs(q), ks{q}); end
    fprintf('See the VERDICT in IFFL_Diagnostics.txt for what this implies.\n');
end
function opt = applyDefaults(opt)
    d.compartment = 'c';
    d.consumingN = false;
    d.fastArmSteps = 2;      % 1 collapses the fast arm to r_Fd only
    d.requireSharedMI = true;% r_I1 and r_I2 must draw on ONE m_I pool
    d.requireMI = true;      % false drops m_I entirely (see iffl_topology)
    d.catalyticInterstage = false;  % true: m_d1 GATES r_A_2 -> m_in_2 (not consumed)
    d.gateMode = {'activate','activate'};  % per stage; {'activate','inhibit'} = hybrid
    d.neuronGate = 'activate';  % 'inhibit' = negative weight, 'any' = either
    d.interstageGate = 'activate';  % how m_d1 drives r_A_2; 'any' flips sign per circuit
    d.carrierTest = true;
    d.distinctEC = false;             % reported, not rejected
    d.noRegulatoryCrosstalk = false;  % reported, not rejected
    d.requireEC = false;
    d.maxArmPaths = 2000;
    d.maxHousekeeping = 80;  % candidates per r_Fr / r_lk slot
    d.sinkMode = 'prefer';   % map / prefer / sink  (see pickHousekeeping)
    d.maxBlocks = 20000;
    d.topBlocks = 400;
    d.maxCircuits = 200;
    d.verbose = true;
    d.currency = {'h','h2o','atp','adp','amp','pi','ppi','nad','nadh','nadp','nadph', ...
                  'co2','o2','nh4','nh3','hco3','coa','fad','fadh2','gtp','gdp','udp', ...
                  'utp','ctp','so4','na1','k','cl','fe2','fe3','mg2','ca2','zn2', ...
                  'q10','q10h2','gthrd','gthox','lgt_S','trdrd','trdox','thf','mlthf'};
    d.currencyNamePatterns = {'carnitine','coenzyme a','glutathione','thioredoxin'};
    f = fieldnames(d);
    for q = 1:numel(f)
        if ~isfield(opt, f{q}), opt.(f{q}) = d.(f{q}); end
    end
    if ~opt.requireMI,           opt.miMode = 'none';
    elseif opt.requireSharedMI,  opt.miMode = 'shared';
    else,                        opt.miMode = 'separate';
    end
end

% ================================ selftest ==================================
function out = selftest()
    names = {'m_in','Ff1','Fd1','Sf1','d1','Ff2','Fd2','Sf2','d2','c','I','b','out', ...
             'w','H','X','carrier','P','camp','junk','junk2', ...
             'pre2','in2','Ff2b','Fd2b','Sf2b','d2b', ...
             'Fd2c','Sf2c','d2c'};
    model.metNames = names(:);
    model.mets = cellfun(@(s) [lower(s) '[c]'], names(:), 'uni', 0);
    model.mets{15} = 'h[c]';
    rx = {'rFf1','rFd1','rFr1','rSf1','rSd1','rIa','rlk1', ...
          'rFf2','rFd2','rFr2','rSf2','rSd2','rIb','rlk2', ...
          'rc','rI1','rI2','rNw','rNb','rNr', ...
          'dCarrier1','dCarrier2','dIcarrier','dIcross','spareFr','rAdapt', ...
          'dFast1','dFast2','rA2','rFf2b','rFd2b','rSf2b','rSd2b','rIb2','rlk2b','rFr2b', ...
          'rc2','rI2b','iFd2','iSf2','iSd2','iFr2','ilk2','irc','irI2'};
    model.rxns = rx(:); model.rxnNames = rx(:);
    ecs = arrayfun(@(k) sprintf('9.9.9.%d',k), 1:numel(rx), 'uni', 0);
    model.rxnECNumbers = ecs(:);
    S = zeros(numel(names), numel(rx));
    S = rxn(S,1 ,1,[2 15]);       S = rxn(S,2 ,2,3);        S = rxn(S,3 ,3,14);
    S = rxn(S,4 ,1,4);            S = rxn(S,5 ,4,5);        S = rxn(S,6 ,5,14);
    S = rxn(S,7 ,5,20);           S = rxn(S,8 ,5,6);        S = rxn(S,9 ,6,7);
    S = rxn(S,10,7,14);           S = rxn(S,11,5,8);        S = rxn(S,12,8,9);
    S = rxn(S,13,9,14);           S = rxn(S,14,9,21);       S = rxn(S,15,9,10);
    S = rxn(S,16,[10 11],14);     S = rxn(S,17,[9 11],14);  S = rxn(S,18,16,13);
    S = rxn(S,19,12,13);          S = rxn(S,20,13,14);
    S = rxn(S,21,[1 17],18);      S = rxn(S,22,18,[17 14]); % conserved carrier arm
    S = rxn(S,23,5,14);                                     % anchor gated by carrier
    S = rxn(S,24,5,[14 10]);                                % anchor that also makes m_c
    S = rxn(S,25,3,14);           S = rxn(S,26,19,1);
    S = rxn(S,27,1,3);            S = rxn(S,28,5,7);   % one-step fast arms
    % an independent D2 branch, reachable ONLY through the catalytic link
    S = rxn(S,29,22,23);          % rA2  : pre2 -> in2      (gated by d1)
    S = rxn(S,30,23,24);          % rFf2b: in2 -> Ff2b
    S = rxn(S,31,24,25);          % rFd2b: Ff2b -> Fd2b
    S = rxn(S,32,23,26);          % rSf2b: in2 -> Sf2b
    S = rxn(S,33,26,27);          % rSd2b: Sf2b -> d2b
    S = rxn(S,34,27,14);          % rIb2 : d2b -> w         (gated by Fd2b)
    S = rxn(S,35,27,21);          % rlk2b: d2b -> junk2
    S = rxn(S,36,25,14);          % rFr2b: Fd2b -> w
    S = rxn(S,37,27,10);          % rc2  : d2b -> c          (counter production)
    S = rxn(S,38,[27 11],14);     % rI2b : d2b + I -> w      (counter sequestration)
    % an INHIBITION-anchored D2: iSd2 produces d2c and is inhibited by Fd2c
    S = rxn(S,39,23,28);          % iFd2 : in2  -> Fd2c      (fast arm, makes the gate)
    S = rxn(S,40,23,29);          % iSf2 : in2  -> Sf2c      (slow arm step 1)
    S = rxn(S,41,29,30);          % iSd2 : Sf2c -> d2c       (ANCHOR: inhibited by Fd2c)
    S = rxn(S,42,28,14);          % iFr2 : Fd2c -> w         (clears the gate)
    S = rxn(S,43,30,21);          % ilk2 : d2c  -> junk2     (output leak)
    S = rxn(S,44,30,10);          % irc  : d2c  -> c         (counter production)
    S = rxn(S,45,[30 11],14);     % irI2 : d2c + I -> w      (counter sequestration)
    model.S = S;
    model.lb = zeros(numel(rx),1); model.ub = 1000*ones(numel(rx),1);
    model.rxnGeneMat = speye(numel(rx));

    T = struct('rxn',{'rIa','rIb','rNw','dIcarrier','dIcross','rA2','rIb2'}, ...
               'EC',{'9.9.9.6','9.9.9.13','9.9.9.18','9.9.9.23','9.9.9.24','9.9.9.29','9.9.9.34'}, ...
               'regulators',{'fd1','fd2','c','carrier','fd1','d1','fd2b; ff2b'});

    out = map_iffl_circuit(model, T, struct('verbose',true));
    one = map_iffl_circuit(model, T, struct('verbose',false,'fastArmSteps',1));
    assert(~isempty(one.circuits), 'selftest FAILED: one-step fast arm found nothing');
    o = one.circuits(1);
    assert(numel(o.sig)==11 && numel(o.rxnIdx)==18, 'one-step role vectors wrong size');
    assert(strcmp(model.rxns{o.rxn.r_Fd_1},'dFast1'), 'one-step fast arm 1 wrong');
    assert(strcmp(model.rxns{o.rxn.r_Fd_2},'dFast2'), 'one-step fast arm 2 wrong');
    assert(~isfield(o.met,'m_Ff_1'), 'm_Ff_1 should not exist with fastArmSteps=1');
    fprintf('one-step variant PASSED (%d circuit(s), %d roles)\n', ...
        numel(one.circuits), numel(o.rxnIdx));
    assert(~isempty(out.circuits), 'selftest FAILED: no circuit found');
    c = out.circuits(1); topo = out.topology;
    fprintf('\n--- selftest circuit ---\n');
    for q = 1:numel(topo.rxnRoles)
        fprintf('  %-8s -> %-9s   %-8s -> %s\n', topo.rxnRoles{q}, model.rxns{c.rxnIdx(q)}, ...
            vget(topo.metRoles,q), vget2(model.mets, c.sig, q));
    end
    nm = @(f) model.rxns{c.rxn.(f)};
    assert(strcmp(nm('r_Ff_1'),'rFf1') && strcmp(nm('r_Fd_1'),'rFd1'), 'fast arm 1 wrong');
    assert(strcmp(nm('r_I_1'),'rIa'), 'anchor 1 wrong');
    assert(strcmp(nm('r_lk_1'),'rlk1'), 'r_lk_1 collided with block 2');
    assert(strcmp(nm('r_lk_2'),'rlk2') && strcmp(nm('r_c'),'rc'), 'block 2 / counter wrong');
    assert(model.mets{c.met.m_I}(1)=='i', 'counter shared species wrong');
    flat = {};
    for k = 1:numel(out.circuits), flat = [flat model.rxns(out.circuits(k).rxnIdx)']; end %#ok<AGROW>
    assert(~any(strcmp(flat,'dCarrier1')), 'carrier decoy accepted: net-mass test broken');
    assert(~any(strcmp(flat,'dIcross')),  'crosstalk decoy accepted: exclusivity test broken');
    % --- hybrid: D1 on activation, D2 on INHIBITION -------------------------
    Tin = struct('rxn',{'iSd2'}, 'EC',{'9.9.9.41'}, 'regulators',{'fd2c'});
    hyb = map_iffl_circuit(model, struct('AC',T,'IN',Tin), ...
        struct('verbose',false,'catalyticInterstage',true,'fastArmSteps',1, ...
               'requireMI',false,'gateMode',{{'activate','inhibit'}}));
    assert(~isempty(hyb.circuits), 'selftest FAILED: hybrid found no circuit');
    hc = hyb.circuits(1); ht = hyb.topology;
    assert(~any(strcmp(ht.rxnRoles,'r_I_2')), 'r_I_2 must be absent on an inhibition stage');
    assert(any(strcmp(ht.rxnRoles,'r_Fr_2')), 'r_Fr_2 must remain: the gate has to be cleared');
    assert(strcmp(model.rxns{hc.rxn.r_Sd_2},'iSd2'), 'inhibition anchor wrong');
    assert(strcmp(model.rxns{hc.rxn.r_Fr_2},'iFr2'), 'gate removal not mapped');
    k = find(strcmp(ht.regTarget,'r_Sd_2'),1);
    assert(~isempty(k) && strcmp(ht.regSense{k},'inhibits'), 'r_Sd_2 should be inhibited');
    % the gate must be CONSUMED by r_Fr_2, and NOT by the reaction it regulates
    mFd2 = hc.met.m_Fd_2;
    assert(model.S(mFd2, hc.rxn.r_Fr_2) < 0, 'm_Fd_2 is not consumed by r_Fr_2');
    assert(model.S(mFd2, hc.rxn.r_Sd_2) == 0, 'm_Fd_2 must not appear in r_Sd_2 stoichiometry');
    fprintf('hybrid gate modes PASSED (%d circuit(s); m_Fd_2 cleared by r_Fr_2, not consumed by r_Sd_2)\n', ...
        numel(hyb.circuits));

    % --- catalytic interstage: m_d1 must be READ, never consumed -------------
    cat = map_iffl_circuit(model, T, struct('verbose',false,'catalyticInterstage',true));
    assert(~isempty(cat.circuits), 'selftest FAILED: catalytic interstage found nothing');
    cc = cat.circuits(1);
    assert(strcmp(model.rxns{cc.rxn.r_A_2},'rA2'), 'r_A_2 not the planted link');
    assert(strcmp(model.mets{cc.met.m_in_2},'in2[c]'), 'm_in_2 wrong');
    mD1 = cc.met.m_d1;
    for q = 1:numel(cc.rxnIdx)
        nmq = cat.topology.rxnRoles{q};
        if any(strcmp(nmq,{'r_I_1','r_lk_1','r_Sd_1'})), continue; end
        assert(model.S(mD1, cc.rxnIdx(q)) == 0, ...
            sprintf('%s touches m_d1 - interstage is not high-impedance', nmq));
    end
    fprintf('catalytic interstage PASSED (%d circuit(s); m_d1 read, not consumed)\n', ...
        numel(cat.circuits));

    if exist('export_iffl_chain','file')
        out.model = model;
        f = export_iffl_chain(model, out, 1, 'camp[c]');
        txt = fileread(f); delete(f);
        assert(~isempty(strfind(txt,'r_Ff_1')), 'export: role notation missing'); %#ok<STREMP>
        assert(~isempty(strfind(txt,'(m_Fd_1)')), 'export: role tags missing'); %#ok<STREMP>
        fprintf('export check PASSED (%d bytes)\n', numel(txt));
    end
    fprintf('\nselftest PASSED: %d circuit(s); carrier arm and m_c-producing anchor rejected\n\n', ...
        numel(out.circuits));
end
function S = rxn(S, j, subs, pros)
    S(subs,j) = -1; S(pros,j) = 1;
end
function s = vget(c, q), if q <= numel(c), s = c{q}; else, s = ''; end, end
function s = vget2(mets, sig, q)
    if q <= numel(sig), s = mets{sig(q)}; else, s = ''; end
end
