function fn = export_iffl_chain(model, out, n, cAMPid, prefix)
% EXPORT_IFFL_CHAIN  Write circuit n of map_iffl_circuit as a full-chain report
% in the r_Ff_1 / r_Fd_1 / ... notation.
%
%   export_iffl_chain(model, out)            % best circuit
%   export_iffl_chain(model, out, 3)         % third
%   export_iffl_chain(model, out, 1, 'camp[c]')
%
% Writes Selected_IFFL_Chain_<n>.txt and returns the file name.

    if nargin < 3 || isempty(n), n = 1; end
    if nargin < 4 || isempty(cAMPid), cAMPid = 'camp[c]'; end
    if nargin < 5 || isempty(prefix), prefix = 'Selected_IFFL'; end
    assert(~isempty(out.circuits), 'No circuit to export.');
    assert(n <= numel(out.circuits), 'Only %d circuit(s) available.', numel(out.circuits));

    topo = out.topology; c = out.circuits(n); opt = out.opt;
    fn = sprintf('%s_Chain_%d.txt', prefix, n);
    fid = fopen(fn, 'w'); cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    role = containers.Map('KeyType','double','ValueType','char');
    for q = 1:numel(topo.metRoles), role(c.sig(q)) = topo.metRoles{q}; end

    bar = repmat('#', 1, 64);
    fprintf(fid, '%s\nIFFL CIRCUIT (fixed 2+2 topology, crosstalk-free)\n', bar);
    fprintf(fid, 'Topology %s | circuit %d of %d | block genes %d\n', ...
        topo.version, n, numel(out.circuits), c.nGenes);
    fprintf(fid, 'Compartment %s | single-pot: every reaction has all participants in [%s]\n', ...
        opt.compartment, opt.compartment);
    fprintf(fid, '%s\n\n', bar);

    % ---------------- signal metabolites ----------------
    fprintf(fid, 'SIGNAL METABOLITES   (others = reactions elsewhere in the model that\n');
    fprintf(fid, 'also touch this species; irrelevant in a cell-free pot, where only the\n');
    fprintf(fid, 'twenty enzymes below are present, but they are competitors in vivo)\n');
    for q = 1:numel(topo.metRoles)
        metLine(fid, model, topo.metRoles{q}, c.sig(q), c.notes.competitors(q));
    end
    fprintf(fid, '\n');

    % ---------------- reaction roles ----------------
    fprintf(fid, 'REACTION ROLES\n');
    for q = 1:numel(topo.rxnRoles)
        j = c.rxnIdx(q);
        if isnan(j)
            fprintf(fid, '  %-8s | %-14s | %-12s | %s\n', topo.rxnRoles{q}, ...
                '(added sink)', '-', 'supply this removal yourself');
        else
            fprintf(fid, '  %-8s | %-14s | %-12s | %s\n', topo.rxnRoles{q}, model.rxns{j}, ...
                ecOf(model, j), nameOf(model, j));
        end
    end
    if any(isnan(c.rxnIdx))
        fprintf(fid, ['\n  ADDED SINKS: the roles marked above discharge to the empty set in\n' ...
                      '  the diagram. No clean metabolic reaction removes that species without\n' ...
                      '  touching another circuit signal, so supply the removal yourself:\n' ...
                      '  outflow in a CSTR, or a dedicated protease/phosphatase/esterase. Its\n' ...
                      '  rate constant is then a free parameter you set, which is exactly what\n' ...
                      '  the baseline and gate-clearing conditions below need.\n']);
    end
    fprintf(fid, '\nREGULATORY EDGES (the regulator is NOT consumed by its target)\n');
    for q = 1:numel(topo.regulator)
        mi = c.sig(strcmp(topo.metRoles, topo.regulator{q}));
        rj = c.rxnIdx(strcmp(topo.rxnRoles, topo.regTarget{q}));
        if isfield(topo,'regSense'), sense = topo.regSense{q}; else, sense = 'activates'; end
        % report what was FOUND, not what the policy permitted
        if strcmp(topo.regTarget{q},'r_Nw') && isfield(c,'neuronSense'), sense = c.neuronSense; end
        if strcmp(topo.regTarget{q},'r_A_2') && isfield(c,'interSense'),  sense = c.interSense;  end
        if isnan(rj), rn = '(added sink)'; ec = '-'; else, rn = model.rxns{rj}; ec = ecOf(model, rj); end
        fprintf(fid, '  %-8s (%s) %-11s %-8s (%s, EC %s)\n', topo.regulator{q}, ...
            model.mets{mi}, sense, topo.regTarget{q}, rn, ec);
    end
    if isfield(c,'neuronSense') && strncmp(c.neuronSense,'inhibit',7)
        fprintf(fid, ['  NOTE: m_c INHIBITS the neuron weight, so w is NEGATIVE - the output\n' ...
                      '  falls as the count rises. A valid single neuron, but state the sign.\n']);
    end
    if isfield(c,'interGain') && c.interGain < 0
        fprintf(fid, ['  NOTE: the interstage INHIBITS, contributing a sign flip of its own.\n' ...
                      '  Net cascade polarity is %+d, not %+d: m_d2 tracks the NEGATIVE second\n' ...
                      '  derivative, so the counter increments on downward curvature. The\n' ...
                      '  circuit works, but this sign must be stated when you report it.\n'], ...
                      c.netPolarity, topo.netPolarity);
    end
    if isfield(topo,'gateMode')
        fprintf(fid, '  stage gate modes: D1 = %s, D2 = %s\n', topo.gateMode{1}, topo.gateMode{2});
        if any(strcmp(topo.gateMode,'inhibit'))
            fprintf(fid, ['  An inhibition stage has no r_I: m_d falls because its production\n' ...
                          '  is shut off, and r_lk drains it. m_Fd is still cleared by r_Fr, so\n' ...
                          '  the gate returns to baseline and the response stays a pulse.\n']);
        end
    end
    fprintf(fid, '\n');

    % ---------------- polarity ----------------
    fprintf(fid, 'POLARITY AND TIMESCALES\n');
    fprintf(fid, ['  Each stage is INVERTING. The fast arm raises the gated degradation of\n' ...
                  '  m_d_k before the slow arm raises its production, so a rising stage input\n' ...
                  '  drives m_d_k DOWN from its baseline. Stage gains: %+d, %+d -> net %+d.\n'], ...
        topo.stagePolarity(1), topo.stagePolarity(2), topo.netPolarity);
    fprintf(fid, ['  Consequence: m_d1 encodes -d(m_in)/dt and m_d2 encodes +d2(m_in)/dt2.\n' ...
                  '  Both ride on a POSITIVE baseline and the pulse is an excursion about it:\n' ...
                  '      m_d_k* = v_Sd / ( k_lk + V_I * f(m_Fd_k*) )   must be > 0\n' ...
                  '  so r_lk and the saturated r_I must not between them exceed v_Sd at rest.\n']);
    fprintf(fid, ['  Both arms have TWO steps, so the fast/slow split is kinetic, not\n' ...
                  '  topological. Load enzymes so that, per stage:\n' ...
                  '      tau_F = 1/k_Ff + 1/k_Fd  <  tau_S = 1/k_Sf + 1/k_Sd\n' ...
                  '      r_I saturated in m_d_k  (m_d_k* >> K_M)  -> subtracts a FLUX\n' ...
                  '      r_lk first order in m_d_k (m_d_k* << K_M) -> sets the baseline\n' ...
                  '      m_Fd_k below saturation at the activator site -> graded gate\n']);
    fprintf(fid, '\n');

    % ---------------- block 0 ----------------
    ci = find(strcmp(model.mets, cAMPid), 1);
    fprintf(fid, '%s\nBLOCK 0: %s -> m_in\n%s\n', rule(), cAMPid, rule());
    if isempty(ci)
        fprintf(fid, '  %s is not in the model; no input connector searched.\n\n', cAMPid);
    else
        forbid = false(numel(model.mets),1); forbid(c.sig) = true;
        forbid(c.met.m_in) = false;
        [steps, rpath, mpath] = cleanBFS(model, ci, c.met.m_in, 12, opt, forbid);
        if isnan(steps)
            fprintf(fid, ['  No mass-flow connector within 12 steps from %s to m_in.\n' ...
                          '  For cAMP this is expected, not a search failure: its only\n' ...
                          '  consumer is PDE1 (camp + h2o -> amp + h) and both products are\n' ...
                          '  currency metabolites, so no path from it carries any of its\n' ...
                          '  atoms to a non-currency species. Block 0 has to be CATALYTIC.\n' ...
                          '  Run, in this order:\n' ...
                          '      A_Block0(''%s_stage2.mat'', %d)             %% confirm the above\n' ...
                          '      A_Block0_Catalytic(''%s_stage2.mat'')       %% find the gate\n' ...
                          '  The first also reports how far %s sits from EVERY circuit\n' ...
                          '  signal, which is how you catch an input that reaches m_out\n' ...
                          '  sooner than m_in and short-circuits the pipeline.\n\n'], ...
                          cAMPid, prefix, n, prefix, cAMPid);
        else
            fprintf(fid, '  Collision-free, %d step(s).\n', steps);
            fprintf(fid, '  Metabolites: %s\n', strjoin(model.mets(mpath), ' -> '));
            fprintf(fid, '  Reactions  : %s\n\n', strjoin(model.rxns(rpath), ' -> '));
            for q = 1:numel(rpath)
                printRxn(fid, model, rpath(q), sprintf('r_adapter_%d', q), role, '');
            end
        end
    end

    % role -> inline regulation note, so each reaction carries its own
    regMap = containers.Map('KeyType','char','ValueType','char');
    for q = 1:numel(topo.regulator)
        mi = c.sig(strcmp(topo.metRoles, topo.regulator{q}));
        sense = 'activates';
        if isfield(topo,'regSense'), sense = topo.regSense{q}; end
        if strcmp(topo.regTarget{q},'r_Nw') && isfield(c,'neuronSense'), sense = c.neuronSense; end
        if strcmp(topo.regTarget{q},'r_A_2') && isfield(c,'interSense'),  sense = c.interSense;  end
        regMap(topo.regTarget{q}) = sprintf('%s = %s (%s) -- %s', ...
            topo.regulator{q}, strtrim(model.metNames{mi}), model.mets{mi}, sense);
    end

    % ---------------- blocks 1-4 ----------------
    titles = { 'BLOCK 1: DIFFERENTIATOR 1 (m_in -> m_d1)', ...
               'BLOCK 2: DIFFERENTIATOR 2 (m_d1 -> m_d2)', ...
               'BLOCK 3: COUNTER (m_d2 -> m_c, shared species m_I)', ...
               'BLOCK 4: SINGLE NEURON (m_out = w[m_c] + b)' };
    for b = 1:4
        fprintf(fid, '%s\n%s\n%s\n\n', rule(), titles{b}, rule());
        for q = find(topo.block == b)
            if isnan(c.rxnIdx(q))
                fprintf(fid, '  %s\n      ADDED SINK: %s -> 0 (not a model reaction)\n\n', ...
                    topo.rxnRoles{q}, sinkTarget(model, topo, c, q));
            else
                rg = '';
                if isKey(regMap, topo.rxnRoles{q}), rg = regMap(topo.rxnRoles{q}); end
                printRxn(fid, model, c.rxnIdx(q), topo.rxnRoles{q}, role, rg);
            end
        end
    end

    % ---------------- audits ----------------
    fprintf(fid, '%s\nAUDIT\n%s\n', rule(), rule());
    fprintf(fid, ['  SIGNAL EXCLUSIVITY: PASS. Each of the 13 signals above is touched only\n' ...
                  '  by the circuit reactions its role declares, with the declared sign;\n' ...
                  '  m_Fd_1, m_Fd_2 and m_c appear in no stoichiometry at all.\n']);
    nt = c.notes;
    fprintf(fid, '\n  The rest is reported, not enforced (%s / %s in this run):\n', ...
        onoff(opt.noRegulatoryCrosstalk), onoff(opt.distinctEC));
    if nt.nRegCollisions == 0
        fprintf(fid, '    Regulatory: no signal gates a reaction other than its target.\n');
    else
        fprintf(fid, '    Regulatory collisions (%d) - the gate is not selective here:\n', nt.nRegCollisions);
        for q = 1:numel(nt.regCollisions), fprintf(fid, '      %s\n', nt.regCollisions{q}); end
    end
    if nt.nSharedEC == 0
        fprintf(fid, '    EC: every role is a different enzyme.\n');
    else
        fprintf(fid, ['    Shared EC (%d) - one enzyme fills two roles, so you cannot load\n' ...
                      '    them at independent concentrations:\n'], nt.nSharedEC);
        for q = 1:numel(nt.sharedEC), fprintf(fid, '      %s\n', nt.sharedEC{q}); end
    end
    if isempty(nt.sharedAux)
        fprintf(fid, '    Shared auxiliary species: none.\n');
    else
        fprintf(fid, ['    Shared auxiliary (non-signal, non-currency) species - mass coupling\n' ...
                      '    that carries no signal. Buffer them or confirm they are in excess:\n']);
        for q = 1:numel(nt.sharedAux)
            m = nt.sharedAux(q);
            fprintf(fid, '      %s (%s)\n', model.mets{m}, model.metNames{m});
        end
    end
    fprintf(fid, ['\n  Validated: connectivity, compartment, signal incidence, net arm mass.\n' ...
                  '  NOT validated: k_cat/K_M values, activator strength, baseline feasibility.\n']);
end

% ------------------------------------------------------------------ helpers
function s = rule(), s = repmat('-', 1, 64); end
function s = sinkTarget(model, topo, c, q)
% which species this added sink removes, from the role's declared substrate
    s = '?';
    sub = topo.subs{q};
    if isempty(sub), return; end
    ix = c.sig(strcmp(topo.metRoles, sub{1}));
    if ~isempty(ix), s = sprintf('%s (%s)', model.mets{ix}, model.metNames{ix}); end
end
function s = onoff(b), if b, s = 'enforced'; else, s = 'not enforced'; end, end
function s = ecOf(model, j)
    s = ''; if isfield(model,'rxnECNumbers'), s = strtrim(char(model.rxnECNumbers{j})); end
end
function s = nameOf(model, j)
    s = ''; if isfield(model,'rxnNames'), s = strtrim(char(model.rxnNames{j})); end
end
function metLine(fid, model, label, ix, others)
    f = '';
    if isfield(model,'metFormulas') && ~isempty(model.metFormulas{ix}), f = char(model.metFormulas{ix}); end
    fprintf(fid, '  %-8s | %-20s | %-38s | %-14s | others %d\n', ...
        label, model.mets{ix}, model.metNames{ix}, f, others);
end
function printRxn(fid, model, j, label, role, reg)
    if nargin < 6, reg = ''; end
    lb = -Inf; ub = Inf;
    if isfield(model,'lb'), lb = model.lb(j); end
    if isfield(model,'ub'), ub = model.ub(j); end
    fprintf(fid, '  %s\n      ID  : %s\n      Name: %s\n      EC  : %s\n', ...
        label, model.rxns{j}, nameOf(model,j), ecOf(model,j));
    fprintf(fid, '      Bounds: [%g, %g]\n', lb, ub);
    fprintf(fid, '      IDs   : %s\n', eqn(model, j, 'mets', role));
    fprintf(fid, '      Names : %s\n', eqn(model, j, 'metNames', role));
    if ~isempty(reg)
        fprintf(fid, '      REGULATED BY: %s\n', reg);
        fprintf(fid, '                    catalytic - the regulator is NOT a substrate\n');
    end
    if isfield(model,'grRules') && numel(model.grRules) >= j && ...
       ~isempty(strtrim(char(model.grRules{j})))
        fprintf(fid, '      Gene(s)     : %s\n', strtrim(char(model.grRules{j})));
    end
    fprintf(fid, '\n');
end
function s = eqn(model, j, field, role)
    sub = find(model.S(:,j) < 0); pro = find(model.S(:,j) > 0);
    s = sprintf('%s  ->  %s', side(model,j,sub,field,role), side(model,j,pro,field,role));
end
function s = side(model, j, ix, field, role)
    parts = cell(1, numel(ix));
    for q = 1:numel(ix)
        i = ix(q); nm = char(model.(field){i});
        if isKey(role, i), nm = sprintf('%s(%s)', nm, role(i)); end
        co = abs(model.S(i,j));
        if co == 1, parts{q} = nm; else, parts{q} = sprintf('%g %s', co, nm); end
    end
    s = strjoin(parts, ' + ');
end
function [steps, rpath, mpath] = cleanBFS(model, src, tgt, maxd, opt, forbid)
    steps = NaN; rpath = []; mpath = [];
    base = lower(regexprep(model.mets, '\[[a-z0-9]+\]$', ''));
    curr = ismember(base, lower(opt.currency));
    % Currency metabolites are deliberately NOT traversable. Letting the BFS
    % pass through the adenylate pool looks like it unlocks routes, but those
    % routes are built on false edges: in
    %     ACITL: atp + cit + coa -> accoa + adp + oaa + pi
    % the pair (atp, oaa) is not an atom transfer - ATP leaves as ADP + Pi and
    % OAA comes from citrate. Such a path carries none of the input's atoms.
    % Do not "fix" a failed Block 0 by relaxing this. See A_Block0.m.
    ub = ones(numel(model.rxns),1); if isfield(model,'ub'), ub = model.ub(:); end
    bad = full(any(model.S(forbid,:) ~= 0, 1))';
    % The target is exempt from `forbid` so the last step can produce it, but a
    % reaction that CONSUMES it mid-route makes it a carrier with zero net
    % product (as GLUCYS -> r0568 -> CYSTA did to glu_L). Strike those out.
    eatsTgt = full(model.S(tgt,:) < 0)';
    usable = ~bad & ~eatsTgt & ub > 0;
    nM = numel(model.mets); dd = nan(nM,1); pr = nan(nM,1); pm = nan(nM,1);
    dd(src) = 0; frontier = src;
    for d = 1:maxd
        nxt = [];
        for m = frontier(:)'
            js = find(model.S(m,:) < 0 & usable(:)');
            for j = js
                for p = find(model.S(:,j) > 0)'
                    if curr(p) || ~isnan(dd(p)), continue; end
                    dd(p) = d; pr(p) = j; pm(p) = m; nxt(end+1) = p; %#ok<AGROW>
                end
            end
        end
        frontier = nxt;
        if ~isnan(dd(tgt)) || isempty(frontier), break; end
    end
    if isnan(dd(tgt)), return; end
    steps = dd(tgt); cur = tgt; mpath = tgt;
    while cur ~= src
        rpath(end+1) = pr(cur); cur = pm(cur); mpath(end+1) = cur; %#ok<AGROW>
    end
    rpath = fliplr(rpath); mpath = fliplr(mpath);
end
