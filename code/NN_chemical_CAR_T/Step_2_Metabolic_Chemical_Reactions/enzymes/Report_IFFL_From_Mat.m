function Report_IFFL_From_Mat(matFile, outFile, topN, prefix)
% REPORT_IFFL_FROM_MAT  Rebuild the whole IFFL diagnostic from the checkpoint
% written by A_Master_IFFL, with no model load, no COBRA and no BRENDA parse.
%
%   Report_IFFL_From_Mat                                  % iffl_stage2.mat
%   Report_IFFL_From_Mat('iffl_stage2.mat')
%   Report_IFFL_From_Mat('iffl_stage2.mat','IFFL_Diagnostics.txt', 3)
%
% Everything it prints also goes to the .txt, so a session that closes on you
% loses nothing. It ends with a VERDICT: the one change worth making next,
% derived from the funnel rather than guessed.

    if nargin < 1 || isempty(matFile), matFile = 'iffl_stage2.mat'; end
    if nargin < 2 || isempty(outFile), outFile = 'IFFL_Diagnostics.txt'; end
    if nargin < 3 || isempty(topN),    topN    = 3; end
    if nargin < 4 || isempty(prefix),  prefix  = 'Selected_IFFL'; end

    assert(exist(matFile,'file') == 2, ...
        ['Checkpoint %s not found. Run A_Master_IFFL first - it writes it even\n' ...
         'when zero circuits survive, which is exactly the case worth inspecting.'], matFile);
    L = load(matFile);
    assert(isfield(L,'out') && isfield(L,'mdl'), ...
        '%s is not an A_Master_IFFL checkpoint (needs ''out'' and ''mdl'').', matFile);
    out = L.out; mdl = L.mdl;

    fid = fopen(outFile,'w'); cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    say = @(varargin) both(fid, varargin{:});

    bar = repmat('=',1,70);
    say('%s\nIFFL MAPPING DIAGNOSTIC\n%s\n', bar, bar);
    say('checkpoint : %s\n', matFile);
    say('written    : %s\n', dstr(L));
    say('topology   : %s\n', out.topology.version);
    o = out.opt;
    say('settings   : fastArmSteps %d | compartment %s | carrierTest %d | maxArmPaths %d\n', ...
        o.fastArmSteps, o.compartment, o.carrierTest, o.maxArmPaths);
    say('             distinctEC %d | noRegulatoryCrosstalk %d | consumingN %d\n\n', ...
        o.distinctEC, o.noRegulatoryCrosstalk, o.consumingN);

    % ---------------- funnel ----------------
    st = out.stats; nJoin = sum(out.depth);
    say('%s\nFUNNEL\n%s\n', bar, bar);
    say('  seed reactions read ............. %d\n', fld(st,'seedsSeen',-1));
    if fld(st,'seedsMissing',0) > 0
        say(['    !! %d seed reactions are NOT IN THIS MODEL. The seed table was built\n' ...
             '       against a different model - almost certainly a cached\n' ...
             '       iffl_stage1.mat from before compartments were collapsed.\n' ...
             '       Delete iffl_stage1.mat and re-run; THIS RESULT IS NOT VALID.\n'], ...
             st.seedsMissing);
    end
    if isfield(st,'anchorDrop') && any(st.anchorDrop)
        lbl = {'not usable in the forward direction (ub <= 0)', ...
               'boundary reaction (exchange / demand / sink)', ...
               sprintf('not entirely in compartment [%s]', fld(o,'compartment','c')), ...
               'transport: a species crosses unchanged', ...
               'no EC number (opt.requireEC)'};
        say('    dropped as unusable:\n');
        for q = 1:5
            if st.anchorDrop(q) > 0
                say('      %6d  %s\n', st.anchorDrop(q), lbl{q});
            end
        end
    end
    say('  anchors examined ................ %d\n', st.anchors);
    say('  arms rejected as carriers ....... %d\n', st.arm_carrier);
    say('  differentiator blocks accepted .. %d\n', st.blocks);
    if fld(st,'pairsTotal',0) > 0
        say('  block pairs enumerated .......... %d\n', st.pairsTotal);
        say('    ...with a usable interstage ... %d  (%.1f%%)\n', ...
            st.pairsLinked, 100*st.pairsLinked/max(st.pairsTotal,1));
    end
    say('  D1->D2 joins tried .............. %d\n', nJoin);
    say('  circuits with exclusive signals . %d\n', numel(out.circuits));
    if isfield(out,'armCapHits') && out.armCapHits > 0
        say(['\n  WARNING: the arm enumeration hit opt.maxArmPaths (%d) on %d target(s).\n' ...
             '  Those arms are TRUNCATED, so the block pool is arbitrarily thinned and\n' ...
             '  every number below is a lower bound. Raise maxArmPaths and re-run before\n' ...
             '  drawing any conclusion.\n'], o.maxArmPaths, out.armCapHits);
    end
    if isfield(st,'gate_pairs')
        say('\n  GATE FUNNEL (why so few distinct m_Fd survive)\n');
        say('    metabolites that activate a usable reaction ..... %d\n', fld(out,'nRegMets',-1));
        say('    (anchor, gate) pairs examined ................... %d\n', st.gate_pairs);
        say('      dropped: gate has no consumer for r_Fr ........ %d\n', st.gate_noFr);
        say('      dropped: no %d-step fast arm into the gate .... %d\n', ...
            fld(o,'fastArmSteps',2), st.gate_noArm);
        say('    (anchor, output) pairs examined ................. %d\n', st.out_pairs);
        say('      dropped: output has no consumer for r_lk ...... %d\n', st.out_noLk);
        say('      dropped: no 2-step slow arm into the output ... %d\n', st.out_noArm);
        say('      dropped: arms share no common source m_in ..... %d  <- largest\n', st.no_common_src);
    end
    say('\n');

    % ---------------- block pool ----------------
    B = out.blocks;
    say('%s\nBLOCK POOL\n%s\n', bar, bar);
    if isempty(B)
        say('  No blocks at all. The problem is upstream of the cascade: either no\n');
        say('  anchor has a usable 2-step arm into both its gate and its output, or the\n');
        say('  compartment gate is too tight. Try opt.compartment = ''any''.\n\n');
    else
        gates = unique([B.mFd]); outs = unique([B.mD]); ins = unique([B.mIn]);
        say('  distinct m_Fd (gates) %d | distinct m_d (outputs) %d | distinct m_in %d\n\n', ...
            numel(gates), numel(outs), numel(ins));
        k = zeros(1,numel(gates));
        for q = 1:numel(gates), k(q) = sum([B.mFd] == gates(q)); end
        [k,ord] = sort(k,'descend'); gates = gates(ord);
        say('  blocks per gate (top %d of %d):\n', min(20,numel(gates)), numel(gates));
        for q = 1:min(20,numel(gates))
            say('    %6d  %-18s %s\n', k(q), mdl.mets{gates(q)}, mdl.metNames{gates(q)});
        end
        hinge = intersect(outs, ins);
        pi_ = fld(out,'poolI',struct('n',0));
        if isfield(pi_,'n') && pi_.n > 0
            pa = out.poolA;
            say('\n  INHIBITION pool: %d blocks | %d gates | %d outputs | %d inputs\n', ...
                pi_.n, numel(pi_.gates), numel(pi_.outs), numel(pi_.ins));
            say('    gates: ');
            for q = 1:min(8,numel(pi_.gates)), say('%s ', mdl.mets{pi_.gates(q)}); end
            if numel(pi_.gates) > 8, say('...'); end
            say('\n');
            shG = numel(intersect(pa.gates, pi_.gates));
            shO = numel(intersect(pa.outs,  pi_.outs));
            shI = numel(intersect(pa.ins,   pi_.ins));
            say('    shared with the activation pool: %d/%d gates, %d outputs, %d inputs\n', ...
                shG, numel(pi_.gates), shO, shI);
            if shG == 0
                say(['    The two pools use DISJOINT gates, so the hybrid did diversify.\n' ...
                     '    Any remaining collision is on arms or intermediates, not gates.\n']);
            elseif shG >= 0.5*numel(pi_.gates)
                say(['    The pools largely SHARE their gates, so the hybrid did not\n' ...
                     '    diversify: inhibition reaches the same chemistry by another route.\n']);
            end
        end
        say('\n  hinges (a species that is some block''s m_d AND another''s m_in): %d\n', numel(hinge));
        if ~isempty(hinge)
            hk = zeros(numel(hinge),2);
            for q = 1:numel(hinge)
                hk(q,1) = sum([B.mD]  == hinge(q));
                hk(q,2) = sum([B.mIn] == hinge(q));
            end
            [~,ord] = sort(min(hk,[],2),'descend'); hinge = hinge(ord); hk = hk(ord,:);
            say('    %-18s %-30s  ending  starting\n', 'id', 'name');
            for q = 1:min(15,numel(hinge))
                say('    %-18s %-30s  %6d  %8d\n', mdl.mets{hinge(q)}, ...
                    trunc(mdl.metNames{hinge(q)},30), hk(q,1), hk(q,2));
            end
        end
        say('\n');
    end

    % ---------------- where the joins stopped ----------------
    say('%s\nWHERE EACH D1->D2 JOIN STOPPED\n%s\n', bar, bar);
    tot = sum(out.depth);
    say('  (%d joins in total)\n', tot);
    for q = 1:numel(out.depthLabels)
        if out.depth(q) > 0
            say('  %7d / %d  (%5.1f%%)  %s\n', out.depth(q), tot, ...
                100*out.depth(q)/max(tot,1), out.depthLabels{q});
        end
    end
    say('\n');

    % ---------------- collisions ----------------
    say('%s\nWHICH TWO ROLES CLAIMED THE SAME OBJECT\n%s\n', bar, bar);
    ck = out.collKeys; cv = out.collVals;
    if isempty(ck)
        say('  none\n\n');
        topColl = '';
    else
        [cv,ord] = sort(cv,'descend'); ck = ck(ord);
        den = sum(cv);
        say('  (%d joins stopped on a role collision)\n', den);
        for q = 1:min(25,numel(ck))
            say('  %7d / %d  (%5.1f%%)  %s\n', cv(q), den, 100*cv(q)/max(den,1), ck{q});
        end
        topColl = ck{1};
        say('\n');
    end

    % ---------------- exclusivity violations ----------------
    vk = fld(out,'violKeys',{}); vv = fld(out,'violVals',[]);
    if ~isempty(vk)
        say('%s\nWHICH ROLE TOUCHED A SIGNAL IT MUST NOT\n%s\n', bar, bar);
        [vv,ord] = sort(vv,'descend'); vk = vk(ord);
        den = sum(vv);
        say('  (%d joins stopped here)\n', den);
        for q = 1:min(25,numel(vk))
            say('  %7d / %d  (%5.1f%%)  %s\n', vv(q), den, 100*vv(q)/max(den,1), vk{q});
        end
        say('\n');
    end

    % ---------------- verdict ----------------
    say('%s\nVERDICT\n%s\n', bar, bar);
    verdict(say, out, mdl, topColl);

    % ---------------- circuits ----------------
    if ~isempty(out.circuits)
        say('\n%s\nCIRCUITS (best %d of %d)\n%s\n', bar, min(topN,numel(out.circuits)), ...
            numel(out.circuits), bar);
        topo = out.topology;
        for n = 1:min(topN, numel(out.circuits))
            c = out.circuits(n);
            say('\n  circuit %d | genes %d | regulatory collisions %d | shared EC %d\n', ...
                n, c.nGenes, c.notes.nRegCollisions, c.notes.nSharedEC);
            for q = 1:numel(topo.rxnRoles)
                say('    %-8s %-14s %s\n', topo.rxnRoles{q}, mdl.rxns{c.rxnIdx(q)}, ...
                    trunc(ecOf(mdl, c.rxnIdx(q)), 14));
            end
        end
        if exist('export_iffl_chain','file') == 2
            say('\n');
            for n = 1:min(topN, numel(out.circuits))
                try
                    f = export_iffl_chain(mdl, out, n, 'camp[c]', prefix);
                    say('  wrote %s\n', f);
                catch ME
                    say('  export of circuit %d failed: %s\n', n, ME.message);
                end
            end
        end
    end
    say('\nWritten to %s\n', outFile);
end

% --------------------------------------------------------------------------
function verdict(say, out, mdl, topColl)
    if ~isempty(out.circuits)
        say('  %d circuit(s) survived. Nothing to unblock - go to the kinetics.\n', numel(out.circuits));
        return;
    end
    if isfield(out,'armCapHits') && out.armCapHits > 0
        say(['  maxArmPaths was hit. Fix that FIRST (opt.maxArmPaths = 20000) and\n' ...
             '  re-run: every count above is a truncated lower bound.\n']);
        return;
    end
    d = out.depth;
    if sum(d) == 0
        say(['  No D1->D2 join was even attempted: no block''s output is another\n' ...
             '  block''s input, so there is no hinge at all. Widen the block pool\n' ...
             '  (opt.fastArmSteps = 1, then opt.compartment = ''any'').\n']);
        return;
    end
    % The DEEPEST stage reached matters more than the biggest bucket: those joins
    % got furthest, so whatever stops them is the last thing between you and a
    % circuit. Report it first, then the largest bucket when it is a different one.
    nz = find(d(1:end-1) > 0);
    if isempty(nz), say('  Nothing to report.\n'); return; end
    deepest = nz(end);
    [~, biggest] = max(d(1:end-1));
    stages = deepest;
    if biggest ~= deepest
        tot = sum(d);
        deepFrac = d(deepest) / max(tot,1);
        if deepFrac < 0.02
            % The deepest bucket is statistically irrelevant. Leading with it
            % points at a rounding error while the real blocker sits upstream.
            stages = [biggest deepest];
            say(['  The deepest stage reached holds only %d of %d joins (%.1f%%), so it is\n' ...
                 '  NOT the story - a handful of pairs got far by luck. The dominant\n' ...
                 '  blocker is reported first; the deep one second, for completeness.\n'], ...
                 d(deepest), tot, 100*deepFrac);
        else
            stages = [deepest biggest];
            say(['  Reported DEEPEST-first. The %d joins at the deepest stage got furthest,\n' ...
                 '  so that is the last thing standing between you and a circuit. The larger\n' ...
                 '  bucket below is a shallower failure affecting a different set of joins.\n'], ...
                 d(deepest));
        end
    end
    for k = stages
        say('\n  [%d joins] %s\n', d(k), out.depthLabels{k});
        stageAdvice(say, out, mdl, topColl, k);
    end
end

function stageAdvice(say, out, mdl, topColl, k)
    o = out.opt;
    if k <= 2
        if ~isempty(strfind(topColl, 'm_Fd_1 = m_Fd_2')) %#ok<STREMP>
            say(['    Both stages want the same gate: stage 2''s fast arm would drive\n' ...
                 '    stage 1''s degrader. Real conflict, but look at WHY.\n']);
            [conc, nm] = gateConc(out, mdl);
            if conc > 0.5
                say(['    %.0f%% of blocks use the single gate %s - the pool has collapsed\n' ...
                     '    onto one subsystem, so both stages are forced into it. The\n' ...
                     '    activator annotation is not the binding constraint; the gate\n' ...
                     '    funnel above shows what is. Widen the structure.\n'], 100*conc, nm);
            end
        elseif ~isempty(strfind(topColl, 'm_Sf')) || ~isempty(strfind(topColl, 'm_Ff')) %#ok<STREMP>
            say(['    The stages collide on an intermediate pool, not on a gate. This is\n' ...
                 '    the hinge forcing both stages into the same neighbourhood of the\n' ...
                 '    network: D2 must CONSUME D1''s output, so it starts where D1 ended.\n']);
            if ~fld(o,'catalyticInterstage',false)
                say(['    opt.catalyticInterstage = true removes this failure mode by\n' ...
                     '    construction: m_d1 gates r_A_2 -> m_in_2 instead of being consumed,\n' ...
                     '    so D2 no longer has to start where D1 ended.\n']);
            elseif o.fastArmSteps == 2
                say('    Try opt.fastArmSteps = 1 first (cheap, already implemented).\n');
            else
                say(['    fastArmSteps is already 1, so the cheap lever is spent. Next:\n' ...
                     '      opt.maxArmPaths = 20000 ; opt.compartment = ''any''\n' ...
                     '    If those do not clear it, the mass-flow hinge is the problem and\n' ...
                     '    the fix is structural: make the interstage CATALYTIC so m_d1\n' ...
                     '    gates a reaction producing D2''s own input instead of being\n' ...
                     '    consumed by D2. D2 is then free to sit anywhere in metabolism.\n']);
            end
        elseif ~isempty(strfind(topColl, 'r_')) %#ok<STREMP>
            say(['    The collision is on a REACTION. r_I_1, r_lk_1 and both of D2''s arm\n' ...
                 '    entries all consume m_d1, so the hinge is crowded and a thin pool\n' ...
                 '    drafts one enzyme twice. Widen: opt.maxArmPaths = 20000, then\n' ...
                 '    opt.compartment = ''any''.\n']);
        else
            say('    The stages overlap on %s - widen the pool and re-run.\n', topColl);
        end
    elseif k >= 3 && k <= 6
        say(['    The counter plumbing is missing: m_d2 needs a consumer that makes a\n' ...
             '    usable m_c, and m_c needs a consumer of its own. Check that m_d2 is not\n' ...
             '    a dead-end species in this model.\n']);
    elseif k == 7
        say('    Current mode: opt.miMode = %s.\n', fld(o,'miMode','shared'));
        if strcmp(fld(o,'miMode','shared'), 'shared')
            say(['    r_I1 (m_c + m_I) and r_I2 (m_d2 + m_I) must draw on ONE shared m_I and\n' ...
                 '    no pair does. Try opt.requireSharedMI = false first.\n']);
        else
            say(['    r_I1 has NO second non-currency substrate at all, so there is no m_I\n' ...
                 '    to name - sharing was never the issue. In this model the consumers of\n' ...
                 '    a lipid m_c are hydrolases whose other substrate is water, which is\n' ...
                 '    currency: they are unimolecular in signal terms.\n' ...
                 '    That is not a defect. With the partner clamped - water, or an m_I pool\n' ...
                 '    replenished externally to m_I_base as your own counter spec does - the\n' ...
                 '    reset is first order in m_c either way, and you program its rate by\n' ...
                 '    enzyme loading, not by [m_I]. So drop the role:\n' ...
                 '      opt.requireMI = false\n' ...
                 '    Keep it true only if m_I must be a finite resource the two removals\n' ...
                 '    genuinely compete for - then it belongs in the circuit and this model\n' ...
                 '    cannot supply it.\n']);
        end
    elseif k >= 8 && k <= 11
        say(['    The neuron is the blocker: m_c must gate a reaction producing an m_out\n' ...
             '    that has both a second producer (r_Nb, the bias) and a consumer (r_Nr).\n' ...
             '    The bias is usually what is missing, and a circuit can run with b = 0 -\n' ...
             '    say the word and I will make r_Nb optional.\n']);
    elseif k == 12
        say(['    These joins built a COMPLETE circuit - cascade, counter and neuron all\n' ...
             '    mapped - and failed only on the removal steps. The four slots draw on\n' ...
             '    DISJOINT candidate pools, so this is not a bad assignment: at least one\n' ...
             '    pool is simply EMPTY, i.e. every consumer of that species also touches\n' ...
             '    another circuit signal, which is unavoidable in a dense lipid cycle.\n' ...
             '    But r_Fr and r_lk discharge to the empty set in your diagram: they are\n' ...
             '    removals, not conversions, and in a cell-free reactor you supply them\n' ...
             '    yourself - outflow, a protease, a phosphatase. They need not be model\n' ...
             '    reactions at all:\n' ...
             '      opt.sinkMode = ''prefer''  (map one when a clean reaction exists,\n' ...
             '                                otherwise declare an added sink)\n' ...
             '    Their rate constants then become free parameters - which is what the\n' ...
             '    baseline and gate-clearing conditions need anyway.\n']);
    elseif k == 13
        say(['    The crosstalk rule is finally the thing rejecting: whole circuits exist\n' ...
             '    but one reaction touches a signal it must not. Never relax it.\n']);
        vk = fld(out,'violKeys',{}); vv = fld(out,'violVals',[]);
        if isempty(vk), return; end
        [~,ix] = max(vv); top = vk{ix};
        say('    Dominant violation: %s\n', top);
        if ~isempty(strfind(top,'r_c ')) || ~isempty(strfind(top,'r_I1')) || ...
           ~isempty(strfind(top,'r_I2')) %#ok<STREMP>
            say(['    A COUNTER reaction is the offender. The counter was attached last and\n' ...
                 '    had the fewest candidates to choose from, so it got whatever consumer\n' ...
                 '    of m_d2 existed - and in the lipid cycle that consumer feeds straight\n' ...
                 '    back into a differentiator pool. Raising opt.maxHousekeeping will not\n' ...
                 '    help: r_c is not a housekeeping pick. The counter needs to sit away\n' ...
                 '    from the lipid cycle, which is what the catalytic interstage buys you.\n']);
        elseif ~isempty(strfind(top,'r_Nw')) || ~isempty(strfind(top,'r_Nb')) || ...
               ~isempty(strfind(top,'r_Nr')) %#ok<STREMP>
            say(['    A NEURON reaction is the offender - the last block attached, same\n' ...
                 '    story. m_out or m_b is colliding with an upstream signal.\n']);
        elseif ~isempty(strfind(top,'does not consume')) || ...
               ~isempty(strfind(top,'does not produce')) %#ok<STREMP>
            say(['    This is a SIGN failure, not a collision: the reaction runs the wrong\n' ...
                 '    way round for the role. Check whether that reaction is reversible in\n' ...
                 '    the model and whether forcing the forward direction is legitimate for\n' ...
                 '    it; if so, the arm enumeration should be allowed to use it reversed.\n']);
        else
            if ~fld(o,'catalyticInterstage',false)
                say(['    A differentiator reaction re-touches another stage''s signal - the\n' ...
                     '    hinge again, one layer deeper: the stages are network neighbours, so\n' ...
                     '    their reactions keep brushing each other''s pools. Set\n' ...
                     '    opt.catalyticInterstage = true to decouple them.\n']);
            else
                say(['    A differentiator reaction re-touches another stage''s signal even\n' ...
                     '    with the catalytic interstage, so the two stages are still landing\n' ...
                     '    in the same subsystem by coincidence rather than by necessity.\n' ...
                     '    Widen: opt.maxArmPaths = 20000, then opt.compartment = ''any''.\n']);
            end
        end
    else
        say(['    Only the optional filters reject now. Set opt.distinctEC = false and\n' ...
             '    opt.noRegulatoryCrosstalk = false; both are reported per circuit anyway.\n']);
    end
end

function [conc, nm] = gateConc(out, mdl)
    conc = 0; nm = '';
    B = out.blocks; if isempty(B), return; end
    g = unique([B.mFd]); k = zeros(1,numel(g));
    for q = 1:numel(g), k(q) = sum([B.mFd] == g(q)); end
    [k,ix] = max(k); conc = k / numel(B); nm = mdl.mets{g(ix)};
end

function v = fld(s, f, dflt)
    if isfield(s, f), v = s.(f); else, v = dflt; end
end
function both(fid, fmt, varargin)
    fprintf(fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end
function s = trunc(s, n)
    s = char(s);
    if numel(s) > n, s = [s(1:n-1) '~']; end
end
function s = ecOf(mdl, j)
    s = '';
    if isfield(mdl,'rxnECNumbers'), s = strtrim(char(mdl.rxnECNumbers{j})); end
end
function s = dstr(L)
    if isfield(L,'stamp'), s = char(L.stamp); else, s = '(no timestamp)'; end
end
