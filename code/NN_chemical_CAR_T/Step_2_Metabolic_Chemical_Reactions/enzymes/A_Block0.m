function R = A_Block0(matFile, n, opts)
% A_BLOCK0  Find the best cAMP -> m_in connector for a mapped circuit, and when
% no clean one exists, say exactly which circuit signals stand in the way.
%
%   A_Block0                                   % IFFL_catalytic_stage2.mat, circuit 1
%   A_Block0('IFFL_catalytic_stage2.mat', 1)
%   A_Block0(matFile, 1, struct('maxTolerate',3))
%   A_Block0('selftest')
%
% Runs off the checkpoint, so it needs no model reload and no re-search.
%
% WHY THIS EXISTS
%   export_iffl_chain refuses any Block 0 route that touches a circuit signal,
%   and then reports only "no collision-free connector". That is true but
%   useless: a route colliding at ONE known point may be perfectly acceptable
%   (if it brushes m_b, which you supply externally anyway) or fatal (if it
%   touches m_out and short-circuits the computation). You cannot tell which
%   from a refusal.
%
% WHAT IT DOES
%   1. Unconstrained reachability: is m_in reachable from cAMP at all, and how
%      far? If not, the problem is connectivity, not crosstalk.
%   2. Strict search: the current rule, forbidding every signal.
%   3. Signal proximity: how far cAMP is from EVERY circuit signal, not just
%      m_in. If the input reaches a downstream signal sooner than it reaches
%      m_in, Block 0 is not merely blocked - it is a SHORT CIRCUIT, and the
%      input would bypass the computation entirely. This matters here: cAMP
%      catabolism runs camp -> AMP -> adenosine -> adenine, and the neuron
%      uses adn (m_b) and ade (m_out).
%   4. Tolerant search: allow 1, then 2, ... signals to be touched, and report
%      the smallest set that admits a route, with the severity of each.
%
% WRITES <tag>_Block0.txt alongside the checkpoint's own outputs.

    if nargin >= 1 && (ischar(matFile)||isstring(matFile)) && strcmpi(char(matFile),'selftest')
        R = selftest(); return;
    end
    if nargin < 1 || isempty(matFile), matFile = 'IFFL_catalytic_stage2.mat'; end
    if nargin < 2 || isempty(n), n = 1; end
    if nargin < 3, opts = struct(); end
    if ~isfield(opts,'cAMP'),        opts.cAMP = 'camp[c]'; end
    if ~isfield(opts,'maxSteps'),    opts.maxSteps = 14; end
    if ~isfield(opts,'maxTolerate'), opts.maxTolerate = 2; end
    % allowThrough un-marks currency metabolites so the BFS may route through
    % them. DEFAULT EMPTY, AND IT SHOULD STAY THAT WAY. Routing through the
    % adenylate pool creates FALSE EDGES: in
    %     ACITL:  atp + cit + coa -> accoa + adp + oaa + pi
    % the BFS reads (atp -> oaa) as a step, but ATP's atoms leave as ADP + Pi
    % while OAA's come from citrate. No atom crosses. Paths built on such
    % edges transmit nothing, however clean they look. Set this only to
    % explore, never to accept a route.
    % opts.natural = true when Block 0 is DESCRIPTIVE - you want the pathway a
    % real cell would use from cAMP to m_in, and will not be loading these
    % enzymes into a reactor. Currency hops are then legitimate (in a cell the
    % adenylate pool really does carry this flux), and signal collisions are
    % reported for information rather than used to disqualify a route.
    if ~isfield(opts,'natural'), opts.natural = false; end
    if ~isfield(opts,'allowThrough')
        if opts.natural
            opts.allowThrough = {'amp','adp','atp','adn','ppi'};
        else
            opts.allowThrough = {};
        end
    end
    if ~isfield(opts,'allowCurrency'), opts.allowCurrency = false; end
    if ~isfield(opts,'outFile'),     opts.outFile = ''; end
    if ~isfield(opts,'verbose'),     opts.verbose = true; end

    assert(exist(matFile,'file')==2, 'Checkpoint %s not found.', matFile);
    L = load(matFile);
    assert(isfield(L,'out') && isfield(L,'mdl'), '%s is not a mapping checkpoint.', matFile);
    out = L.out; mdl = L.mdl;
    assert(~isempty(out.circuits), 'That checkpoint holds no circuits.');
    assert(n <= numel(out.circuits), 'Only %d circuit(s) in %s.', numel(out.circuits), matFile);

    if isempty(opts.outFile)
        [~,base] = fileparts(matFile);
        base = regexprep(base, '_stage2$', '');
        opts.outFile = [base '_Block0.txt'];
    end
    fid = fopen(opts.outFile,'w'); cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    say = @(varargin) both(fid, opts, varargin{:});

    R = runBlock0(say, mdl, out, n, opts);
    R.outFile = opts.outFile;
    say('\nWritten to %s\n', opts.outFile);
end

% ===========================================================================
function R = runBlock0(say, mdl, out, n, opts)
    c = out.circuits(n); topo = out.topology;
    sig = c.sig; roles = topo.metRoles;
    mIn = c.met.m_in;
    src = find(strcmp(mdl.mets, opts.cAMP), 1);

    bar = repmat('=',1,70);
    say('%s\nBLOCK 0 ANALYSIS  (%s -> m_in)\n%s\n', bar, opts.cAMP, bar);
    say('circuit %d of %d | m_in = %s (%s)\n\n', n, numel(out.circuits), ...
        mdl.mets{mIn}, strtrim(mdl.metNames{mIn}));

    if isempty(src)   % fall back to a base-name match in any compartment
        base = lower(regexprep(mdl.mets, '\[[a-z0-9]+\]$', ''));
        want = lower(regexprep(opts.cAMP, '\[[a-z0-9]+\]$', ''));
        alt = find(strcmp(base, want));
        if isempty(alt)
            say('  %s is not in the model, under any compartment. Nothing to search.\n', ...
                opts.cAMP);
            R = struct('reachable',false); return;
        end
        src = alt(1);
        say('  NOTE: %s not found; using %s instead.\n\n', opts.cAMP, mdl.mets{src});
        opts.cAMP = mdl.mets{src};
    end

    curr = currencyMask(mdl, out.opt, opts);
    currStrict = currencyMask(mdl, out.opt, ...
                     struct('allowThrough',{{}},'allowCurrency',false));

    % ---- 0. what can cAMP even become? ------------------------------------
    say('0. DIRECT FATE OF %s  (every usable reaction consuming it)\n', opts.cAMP);
    cons = find(mdl.S(src,:) < 0 & reshape(usableMask(mdl),1,[]));
    if isempty(cons)
        say('   None. %s has no usable consuming reaction, so Block 0 is\n', opts.cAMP);
        say('   impossible by mass flow whatever else we relax.\n\n');
    end
    nBlockedByCurr = 0;
    for j = reshape(cons,1,[])
        ps = find(mdl.S(:,j) > 0);
        tags = cell(1,numel(ps));
        for q = 1:numel(ps)
            t = 'free';
            if ismember(ps(q), sig)
                t = sprintf('SIGNAL %s', roles{find(sig==ps(q),1)});
            elseif currStrict(ps(q))
                t = 'currency';
                if ~curr(ps(q)), t = 'currency (allowed through)'; end
            end
            tags{q} = sprintf('%s [%s]', mdl.mets{ps(q)}, t);
        end
        say('   %-14s %s\n', mdl.rxns{j}, strjoin(tags, ', '));
        if ~isempty(ps) && all(currStrict(ps)) && ~any(ismember(ps, sig))
            nBlockedByCurr = nBlockedByCurr + 1;
        end
    end
    if nBlockedByCurr > 0 && nBlockedByCurr == numel(cons) && opts.natural
        say(['\n   NOTE. Every reaction consuming %s produces only currency\n' ...
             '   metabolites, so the route below necessarily passes through the\n' ...
             '   adenylate pool. In a cell that is exactly what happens, so it is\n' ...
             '   fine for a descriptive pathway. It would NOT be a valid design for\n' ...
             '   an input stage you actually load, because a currency pool is\n' ...
             '   buffered and carries none of cAMP''s atoms onward - run without\n' ...
             '   natural mode if you ever need that stricter answer.\n'], opts.cAMP);
    elseif nBlockedByCurr > 0 && nBlockedByCurr == numel(cons)
        say(['\n   STRUCTURAL RESULT. Every reaction consuming %s produces nothing\n' ...
             '   but currency metabolites. Currency pools are buffered and shared, and\n' ...
             '   a path through one carries no atoms of the input, so there is NO\n' ...
             '   atom-transferring route from %s to any non-currency species. That\n' ...
             '   holds for every candidate m_in, not just this one: a mass-flow\n' ...
             '   Block 0 from %s does not exist in this model.\n' ...
             '\n' ...
             '   Block 0 must therefore be a CATALYTIC coupling - %s gating a\n' ...
             '   reaction that produces m_in, the way r_A_2 gates m_in_2 without\n' ...
             '   being consumed. Run A_Block0_Catalytic to search for one.\n'], ...
             opts.cAMP, opts.cAMP, opts.cAMP, opts.cAMP);
        if ~isempty(opts.allowThrough)
            say(['\n   WARNING: opts.allowThrough = {%s} is set, so the searches below\n' ...
                 '   may route through those pools. Any route they return is likely to\n' ...
                 '   be built on FALSE EDGES - a (substrate, product) pair in a reaction\n' ...
                 '   that transfers no atoms, e.g. (atp, oaa) in\n' ...
                 '       ACITL: atp + cit + coa -> accoa + adp + oaa + pi\n' ...
                 '   Check the net stoichiometry in section 2: if the input appears on\n' ...
                 '   the product side (as adp, say), it contributed nothing to m_in and\n' ...
                 '   the route is an artifact. Treat what follows as exploration only.\n'], ...
                 strjoin(opts.allowThrough, ', '));
        end
    end
    say('\n');

    % ---- 1. unconstrained reachability ------------------------------------
    [d0, rp0, mp0] = fwdBFS(mdl, src, mIn, opts.maxSteps, curr, false(numel(mdl.mets),1));
    if opts.natural
        say('1. NATURAL ROUTE  <-- THE ANSWER (currency hops allowed, as in a cell)\n');
    else
        say('1. UNCONSTRAINED REACHABILITY (ignoring every circuit signal)\n');
    end
    if isnan(d0)
        say(['   m_in is NOT reachable from %s within %d steps, even ignoring every\n' ...
             '   circuit signal. This is a CONNECTIVITY problem, not a crosstalk one,\n' ...
             '   so tolerating signals cannot help. Next levers, in order:\n' ...
             '     opts.allowThrough  = {''amp'',''adp'',''atp'',''adn''}   (currently: %s)\n' ...
             '     opts.allowCurrency = true                     (drop the filter)\n' ...
             '     opts.maxSteps      = %d                       (currently %d)\n' ...
             '   If none of those help, this m_in simply cannot be driven from cAMP\n' ...
             '   by mass flow, and Block 0 must be a regulatory coupling instead.\n\n'], ...
             opts.cAMP, opts.maxSteps, ...
             iff(isempty(opts.allowThrough), 'none', strjoin(opts.allowThrough, ',')), ...
             opts.maxSteps + 6, opts.maxSteps);
        R = struct(); R.reachable = false; return;
    end
    say('   reachable in %d step(s).\n', d0);
    say('   metabolites: %s\n', strjoin(mdl.mets(mp0), ' -> '));
    say('   reactions  : %s\n', strjoin(mdl.rxns(rp0), ' -> '));
    if opts.natural
        say('\n');
        for q = 1:numel(rp0)
            showRxn(say, mdl, rp0(q), sprintf('step %d', q));
            hits = setdiff(intersect(find(mdl.S(:,rp0(q)) ~= 0), sig), mIn);
            for m = reshape(hits,1,[])
                rl = roles{find(sig==m,1)};
                if mdl.S(m,rp0(q)) < 0, verb = 'consumes'; else, verb = 'produces'; end
                say('        note: also %s %s (%s) - fine for a descriptive\n', ...
                    verb, rl, strtrim(mdl.metNames{m}));
                say('              pathway, but see section 2b before you build it\n');
            end
        end
        net0 = sum(mdl.S(:, rp0), 2);
        say('\n   NET (unit flux per step): %s\n', netEqn(mdl, net0, curr));
        say('   net %s: %+g\n\n', mdl.mets{mIn}, full(net0(mIn)));
    else
        say(['   This route IGNORES circuit signals. It is here only to prove that\n' ...
             '   m_in can be reached at all - it is NOT a Block 0 candidate. Section 2\n' ...
             '   says whether it is usable, and 2b says what disqualifies it if not.\n' ...
             '   If Block 0 is DESCRIPTIVE rather than something you will build, pass\n' ...
             '   struct(''natural'',true) and this route becomes the answer.\n\n']);
    end

    % ---- 2. strict: the rule the exporter uses ----------------------------
    forbidAll = false(numel(mdl.mets),1);
    forbidAll(sig) = true; forbidAll(mIn) = false;
    [d1, rp1, mp1] = fwdBFS(mdl, src, mIn, opts.maxSteps, curr, forbidAll);
    if opts.natural
        say(['2. SIGNAL-FREE ALTERNATIVE (only matters if you ever build Block 0;\n' ...
             '   ignore it for a descriptive pathway)\n']);
    else
        say('2. STRICT (no circuit signal may be touched)  <-- THE USABLE ROUTE\n');
    end
    clean = ~isnan(d1);
    if clean
        say('   clean route found in %d step(s).\n', d1);
        say('   metabolites: %s\n', strjoin(mdl.mets(mp1), ' -> '));
        say('   reactions  : %s\n\n', strjoin(mdl.rxns(rp1), ' -> '));
        for q = 1:numel(rp1)
            showRxn(say, mdl, rp1(q), sprintf('r_adapter_%d', q));
        end
        net = sum(mdl.S(:, rp1), 2);
        say('\n   NET STOICHIOMETRY (unit flux through each step):\n');
        say('        %s\n', netEqn(mdl, net, curr));
        say('   net %s: %+g   <-- must be positive, or the route is a carrier\n', ...
            mdl.mets{mIn}, full(net(mIn)));
        say('\n');
    else
        say('   no route. Every path consumes or produces a circuit signal.\n\n');
    end

    % ---- 2b. why was the shorter unconstrained route rejected? ------------
    if (opts.natural && ~isnan(d0)) || (clean && d0 < d1)
        if opts.natural
            say(['2b. WHERE THE NATURAL ROUTE MEETS THE CIRCUIT (informational)\n' ...
                 '   The %d-step route in section 1 touches these circuit signals.\n' ...
                 '   Harmless if Block 0 stays on paper; worth a sentence in the\n' ...
                 '   write-up, since it says the circuit sits in the same metabolic\n' ...
                 '   neighbourhood as its own input:\n'], d0);
        else
            say(['2b. WHY NOT THE SHORTER %d-STEP ROUTE FROM SECTION 1?\n' ...
                 '   Because these reactions in it touch a circuit signal:\n'], d0);
        end
        anyHit = false;
        for j = reshape(rp0,1,[])
            hits = intersect(find(mdl.S(:,j) ~= 0), sig);
            hits = setdiff(hits, mIn);
            if isempty(hits), continue; end
            anyHit = true;
            for m = reshape(hits,1,[])
                rl = roles{find(sig==m,1)};
                if mdl.S(m,j) < 0, verb = 'CONSUMES'; else, verb = 'PRODUCES'; end
                say('     %-14s %s %s (%s)\n', mdl.rxns{j}, verb, rl, ...
                    strtrim(mdl.metNames{m}));
                say('                    severity: %s\n', severity(rl));
            end
            showRxn(say, mdl, j, 'rejected');
        end
        if ~anyHit
            say(['     ...none, which should not happen. If you see this, the two\n' ...
                 '     searches disagree and there is a bug - send me this report.\n']);
        elseif opts.natural
            say(['   In a cell these fluxes coexist, so none of this invalidates the\n' ...
                 '   section 1 pathway. It only means the natural route could not be\n' ...
                 '   reused as a built input stage without rework.\n']);
        else
            say(['   Those are hard collisions, not preferences: the signal would be\n' ...
                 '   drained or spiked by the input stage. Use the section 2 route.\n']);
        end
        say('\n');
    end

    % ---- 3. how close is cAMP to every OTHER signal? ----------------------
    % A signal nearer to cAMP than m_in is a short circuit: the input would
    % reach the middle or the end of the pipeline before reaching its start.
    say('3. SIGNAL PROXIMITY (unconstrained distance from %s to each signal)\n', opts.cAMP);
    prox = nan(numel(sig),1);
    for i = 1:numel(sig)
        prox(i) = fwdBFS(mdl, src, sig(i), opts.maxSteps, curr, false(numel(mdl.mets),1));
    end
    [~,ord] = sort(prox);
    short = {};
    for i = reshape(ord,1,[])
        if isnan(prox(i))
            tag = sprintf('%8s', 'unreach');
        else
            tag = sprintf('%6d  ', prox(i));
        end
        flag = '';
        if sig(i) == mIn
            flag = '   <- m_in, the intended entry point';
        elseif ~isnan(prox(i)) && prox(i) < d0
            flag = '   <<< SHORT CIRCUIT: nearer than m_in';
            short{end+1} = roles{i}; %#ok<AGROW>
        end
        say('   %s %-10s %-18s %s%s\n', tag, roles{i}, mdl.mets{sig(i)}, ...
            strtrim(mdl.metNames{sig(i)}), flag);
    end
    say('\n');
    if ~isempty(short) && opts.natural
        say(['   %d signal(s) sit closer to %s than m_in does: %s. For a descriptive\n' ...
             '   pathway this is an observation, not a fault - you are not adding cAMP\n' ...
             '   to a reactor, so nothing leaks anywhere. It is still worth a line in\n' ...
             '   the write-up: the input lies nearer the circuit''s own output species\n' ...
             '   than its input, which is why cAMP could not serve as a built input.\n\n'], ...
            numel(short), opts.cAMP, strjoin(short, ', '));
    elseif ~isempty(short)
        say(['   READ THIS FIRST. %d signal(s) sit CLOSER to %s than m_in does:\n' ...
             '   %s.\n' ...
             '   That is not a Block 0 problem, it is a circuit-level one. Whatever\n' ...
             '   connector you build, cAMP will also leak into those species by a\n' ...
             '   shorter route, and the shorter route wins. If one of them is m_out\n' ...
             '   or m_b, the input is wired straight to the answer and the four\n' ...
             '   blocks in between compute nothing observable.\n\n'], ...
            numel(short), opts.cAMP, strjoin(short, ', '));
    else
        say(['   No signal is nearer to %s than m_in. The pipeline is correctly\n' ...
             '   ordered with respect to the input - good.\n\n'], opts.cAMP);
    end

    if clean
        R = struct(); R.reachable = true; R.clean = true; R.steps = d1;
        R.prox = prox; R.shortCircuit = {short};
        R.route = {mdl.rxns(rp1)}; R.routeMets = {mdl.mets(mp1)};
        say('   Strict route exists, so no tolerance is needed. Done.\n');
        return;
    end

    % ---- 4. tolerant: smallest set of signals that must be touched --------
    cand = setdiff(1:numel(sig), find(sig == mIn));   % m_in is the target
    say('4. TOLERANT (allow a small number of signals to be touched)\n');
    found = {};
    for k = 1:min(opts.maxTolerate, numel(cand))
        combos = nchoosek(cand, k);
        for r = 1:size(combos,1)
            tol = combos(r,:);
            fb = forbidAll; fb(sig(tol)) = false;
            [dk, rpk, mpk] = fwdBFS(mdl, src, mIn, opts.maxSteps, curr, fb);
            if ~isnan(dk)
                found(end+1,:) = {k, tol, dk, rpk, mpk}; %#ok<AGROW>
            end
        end
        if ~isempty(found), break; end       % smallest k wins; stop there
    end

    if isempty(found)
        say(['   No route even tolerating up to %d signal(s). The circuit occupies\n' ...
             '   the corridor between %s and m_in too thoroughly for any Block 0.\n'], ...
             opts.maxTolerate, opts.cAMP);
        R = struct(); R.reachable = true; R.clean = false;
        R.prox = prox; R.shortCircuit = {short}; R.best = {}; R.all = {};
        return;
    end

    % rank by steps, then report
    steps = cell2mat(found(:,3));
    [~,ord] = sort(steps); found = found(ord,:);
    say('   %d route(s) exist if %d signal(s) may be touched. Best first:\n\n', ...
        size(found,1), found{1,1});

    show = min(5, size(found,1));
    for q = 1:show
        tol = found{q,2}; dk = found{q,3}; rpk = found{q,4}; mpk = found{q,5};
        nm = cell(1,numel(tol));
        for t = 1:numel(tol)
            nm{t} = sprintf('%s (%s)', roles{tol(t)}, mdl.mets{sig(tol(t))});
        end
        say('   [%d] %d step(s), touches: %s\n', q, dk, strjoin(nm, ', '));
        say('       metabolites: %s\n', strjoin(mdl.mets(mpk), ' -> '));
        say('       reactions  : %s\n', strjoin(mdl.rxns(rpk), ' -> '));
        % where exactly is the collision, and in which direction?
        for t = 1:numel(tol)
            m = sig(tol(t));
            for j = rpk
                v = mdl.S(m, j);
                if v < 0
                    say('       !! %s CONSUMES %s (%s)\n', mdl.rxns{j}, ...
                        roles{tol(t)}, strtrim(mdl.metNames{m}));
                elseif v > 0
                    say('       !! %s PRODUCES %s (%s)\n', mdl.rxns{j}, ...
                        roles{tol(t)}, strtrim(mdl.metNames{m}));
                end
            end
            say('       severity: %s\n', severity(roles{tol(t)}));
        end
        say('\n');
    end

    say(['   HOW TO READ THIS. A Block 0 that brushes m_b is usually tolerable:\n' ...
         '   the bias is supplied externally, so buffer it and the perturbation is\n' ...
         '   absorbed. One that touches a GATE (m_Fd_1, m_Fd_2, m_c) or an OUTPUT\n' ...
         '   (m_d1, m_d2, m_out) is not: those carry the computation, nothing else\n' ...
         '   buffers them, and the input would be writing directly into the result.\n']);

    R = struct(); R.reachable = true; R.clean = false;
    R.prox = prox; R.shortCircuit = {short};
    R.best = found(1,:); R.all = found;
end

% ---------------------------------------------------------------- helpers
function s = severity(role)
    if any(strcmp(role, {'m_out','m_d1','m_d2'}))
        s = 'SEVERE - a computed signal; Block 0 would write into the result';
    elseif any(strcmp(role, {'m_Fd_1','m_Fd_2','m_c'}))
        s = 'SEVERE - a gate; nothing buffers it, so the regulation shifts';
    elseif strcmp(role, 'm_b')
        s = 'MILD - the bias is supplied externally; buffer it and carry on';
    elseif strcmp(role, 'm_in_2')
        s = 'MODERATE - D2 input; shifts stage 2 baseline';
    else
        s = 'MODERATE - an arm intermediate; shifts one stage baseline';
    end
end

function showRxn(say, mdl, j, label)
% one reaction, with its equation and whatever annotation the model carries
    ec = ''; gn = '';
    if isfield(mdl,'rxnECNumbers') && numel(mdl.rxnECNumbers) >= j
        ec = strtrim(char(mdl.rxnECNumbers{j}));
    end
    if isfield(mdl,'grRules') && numel(mdl.grRules) >= j
        gn = strtrim(char(mdl.grRules{j}));
    end
    nm = '';
    if isfield(mdl,'rxnNames') && numel(mdl.rxnNames) >= j
        nm = strtrim(char(mdl.rxnNames{j}));
    end
    say('     [%s] %s', label, mdl.rxns{j});
    if ~isempty(nm), say('  %s', nm); end
    say('\n');
    say('        %s\n', eqn(mdl, j));
    if ~isempty(ec), say('        EC: %s\n', ec); end
    if ~isempty(gn), say('        Gene(s): %s\n', gn); end
end

function s = netEqn(mdl, net, curr)
% net equation, currency folded away so the real transformation is visible
    keep = find(net ~= 0 & ~curr(:));
    lhs = keep(net(keep) < 0); rhs = keep(net(keep) > 0);
    s = sprintf('%s  ->  %s', netSide(mdl, net, lhs), netSide(mdl, net, rhs));
    hid = find(net ~= 0 & curr(:));
    if ~isempty(hid)
        s = sprintf('%s      [+ currency: %s]', s, strjoin(mdl.mets(hid), ', '));
    end
end

function s = netSide(mdl, net, ix)
    if isempty(ix), s = '0'; return; end
    p = cell(1,numel(ix));
    for q = 1:numel(ix)
        co = abs(full(net(ix(q))));
        if co == 1, p{q} = mdl.mets{ix(q)};
        else, p{q} = sprintf('%g %s', co, mdl.mets{ix(q)}); end
    end
    s = strjoin(p, ' + ');
end

function s = eqn(mdl, j)
    subs = find(mdl.S(:,j) < 0); prods = find(mdl.S(:,j) > 0);
    s = sprintf('%s  ->  %s', side(mdl, j, subs), side(mdl, j, prods));
end

function s = side(mdl, j, ix)
    if isempty(ix), s = '0'; return; end
    p = cell(1,numel(ix));
    for q = 1:numel(ix)
        co = abs(mdl.S(ix(q), j));
        if co == 1, p{q} = mdl.mets{ix(q)};
        else, p{q} = sprintf('%g %s', co, mdl.mets{ix(q)}); end
    end
    s = strjoin(p, ' + ');
end

function s = iff(c, a, b)
    if c, s = a; else, s = b; end
end

function u = usableMask(mdl)
    nR = size(mdl.S,2);
    if isfield(mdl,'ub'), u = mdl.ub(:) > 0; else, u = true(nR,1); end
end

function [d, rpath, mpath] = fwdBFS(mdl, src, tgt, maxd, curr, forbid)
% Delegates to block0_route, which owns the search: the forbidden-metabolite
% rule, the carrier rule (no reaction may CONSUME the target), and currency
% skipping all live there. This file and export_iffl_chain.m used to run
% separate BFS copies, which is how the chain file came to report "no
% connector" while this report found a 3-step route. One implementation now.
    R = block0_route(mdl, struct('src', src, 'tgt', tgt, 'maxSteps', maxd, ...
                                 'curr', curr, 'forbid', forbid));
    d = R.steps; rpath = R.rpath; mpath = R.mpath;
end

function curr = currencyMask(mdl, opt, opts)
% Asks block0_route for the mask it would build, so the currency list has a
% single definition. A self-targeted call returns immediately and does no
% search; only R.curr is wanted.
    at = {};
    if isstruct(opts) && isfield(opts,'allowThrough'), at = opts.allowThrough; end
    if isstruct(opts) && isfield(opts,'allowCurrency') && opts.allowCurrency
        curr = false(size(mdl.S,1),1); return;
    end
    R = block0_route(mdl, struct('src', 1, 'tgt', 1, 'opt', opt, ...
                                 'allowThrough', {at}));
    curr = R.curr;
end

function both(fid, opts, fmt, varargin)
    if opts.verbose, fprintf(fmt, varargin{:}); end
    fprintf(fid, fmt, varargin{:});
end

% ================================ selftest ==================================
function R = selftest()
% camp -> a -> b -> target, but the only route passes through a circuit signal
    mdl.mets     = {'camp[c]';'a[c]';'b[c]';'tgt[c]';'sigX[c]';'h[c]'};
    mdl.metNames = {'cAMP';'A';'B';'Target';'SignalX';'H+'};
    mdl.rxns     = {'R1';'R2';'R3';'R4'};
    mdl.rxnNames = mdl.rxns;
    S = zeros(6,4);
    S(1,1) = -1; S(2,1) = 1;             % R1 : camp -> a
    S(2,2) = -1; S(5,2) = 1;             % R2 : a -> sigX        (touches the signal)
    S(5,3) = -1; S(3,3) = 1;             % R3 : sigX -> b        (touches the signal)
    S(3,4) = -1; S(4,4) = 1;             % R4 : b -> tgt
    mdl.S = S; mdl.lb = zeros(4,1); mdl.ub = 1000*ones(4,1);

    topo.metRoles = {'m_in','m_out'};
    c.sig = [4 5];                       % m_in = tgt, m_out = sigX
    c.met.m_in = 4;
    out.circuits = c; out.topology = topo; out.opt = struct('currency',{{'h'}});

    tmp = [tempname '.mat'];
    save(tmp, 'out', 'mdl', '-v7');
    R = A_Block0(tmp, 1, struct('cAMP','camp[c]','verbose',false, ...
                                'outFile',[tempname '.txt']));
    delete(tmp);

    assert(R.reachable, 'target should be reachable when signals are ignored');
    assert(~R.clean, 'the strict search must fail: the only route uses sigX');
    assert(~isempty(R.best), 'a tolerant route should exist at k = 1');
    tol = R.best{2};
    assert(numel(tol) == 1 && tol == 2, 'the tolerated signal should be m_out');
    assert(isequal(R.prox(:)', [4 2]), 'proximities should be tgt 4, sigX 2');
    assert(isequal(R.shortCircuit{1}, {'m_out'}), ...
        'm_out sits 2 steps from camp but m_in sits 4: that is a short circuit');
    % --- scenario 2: the currency bottleneck (camp -> amp -> tgt) ----------
    m2.mets     = {'camp[c]';'amp[c]';'tgt[c]';'h[c]'};
    m2.metNames = {'cAMP';'AMP';'Target';'H+'};
    m2.rxns     = {'PDE';'AMPASE'};
    m2.rxnNames = m2.rxns;
    T = zeros(4,2);
    T(1,1) = -1; T(2,1) = 1;             % PDE    : camp -> amp
    T(2,2) = -1; T(3,2) = 1;             % AMPASE : amp  -> tgt
    m2.S = T; m2.lb = zeros(2,1); m2.ub = 1000*ones(2,1);

    t2.metRoles = {'m_in'};
    c2.sig = 3; c2.met.m_in = 3;
    out2.circuits = c2; out2.topology = t2;
    out2.opt = struct('currency',{{'h','amp'}});   % amp IS currency, as in the real opt

    tmp2 = [tempname '.mat'];
    out = out2; mdl = m2; save(tmp2, 'out', 'mdl', '-v7');

    blocked = A_Block0(tmp2, 1, struct('cAMP','camp[c]','verbose',false, ...
                  'allowThrough',{{}}, 'outFile',[tempname '.txt']));
    opened  = A_Block0(tmp2, 1, struct('cAMP','camp[c]','verbose',false, ...
                  'allowThrough',{{'amp'}}, 'outFile',[tempname '.txt']));
    delete(tmp2);

    assert(~blocked.reachable, ...
        'with amp as currency the BFS must die at step 1, not find a route');
    assert(opened.reachable && opened.clean && opened.steps == 2, ...
        'allowing amp through must expose the 2-step camp -> amp -> tgt route');

    % --- scenario 3: a short dirty route vs a longer clean one (section 2b) -
    m3.mets     = {'camp[c]';'amp[c]';'tgt[c]';'sigG[c]';'x[c]';'y[c]'};
    m3.metNames = {'cAMP';'AMP';'Target';'GateG';'X';'Y'};
    m3.rxns     = {'PDE';'FAST';'S1';'S2';'S3'};
    m3.rxnNames = {'phosphodiesterase';'fast but dirty';'s1';'s2';'s3'};
    U = zeros(6,5);
    U(1,1) = -1; U(2,1) = 1;                           % PDE  : camp -> amp
    U(2,2) = -1; U(4,2) = -1; U(3,2) = 1;              % FAST : amp + sigG -> tgt
    U(2,3) = -1; U(5,3) = 1;                           % S1   : amp -> x
    U(5,4) = -1; U(6,4) = 1;                           % S2   : x   -> y
    U(6,5) = -1; U(3,5) = 1;                           % S3   : y   -> tgt
    m3.S = U; m3.lb = zeros(5,1); m3.ub = 1000*ones(5,1);

    t3.metRoles = {'m_in','m_Fd_1'};
    c3.sig = [3 4]; c3.met.m_in = 3;
    out3.circuits = c3; out3.topology = t3;
    out3.opt = struct('currency',{{'amp'}});

    tmp3 = [tempname '.mat']; txt3 = [tempname '.txt'];
    out = out3; mdl = m3; save(tmp3, 'out', 'mdl', '-v7');
    mixed = A_Block0(tmp3, 1, struct('cAMP','camp[c]','verbose',false, ...
                  'allowThrough',{{'amp'}}, 'outFile',txt3));
    rep3 = fileread(txt3);
    delete(tmp3); delete(txt3);

    assert(mixed.clean && mixed.steps == 4, ...
        'the clean route should be the 4-step one, not the 2-step dirty one');
    assert(isequal(mixed.route{1}(:)', {'PDE','S1','S2','S3'}), ...
        'the printed route should be the clean one');
    assert(~isempty(strfind(rep3, '2b. WHY NOT THE SHORTER 2-STEP ROUTE')), ...
        'section 2b must explain the rejection of the shorter route');
    assert(~isempty(strfind(rep3, 'FAST')) && ~isempty(strfind(rep3, 'CONSUMES m_Fd_1')), ...
        'section 2b must name FAST and the gate it consumes');

    % --- scenario 4: the carrier trap (GLUCYS -> r0568 -> CYSTA in miniature) -
    % A route exists that consumes the target mid-way and regenerates it at the
    % end: net production of tgt is zero. It must be rejected, and the honest
    % longer route taken instead.
    m4.mets     = {'camp[c]';'amp[c]';'tgt[c]';'w[c]';'v[c]'};
    m4.metNames = {'cAMP';'AMP';'Target';'W';'V'};
    m4.rxns     = {'PDE';'EAT';'GIVEBACK';'LONG1';'LONG2'};
    m4.rxnNames = {'phosphodiesterase';'consumes target';'returns target'; ...
                   'honest step 1';'honest step 2'};
    V = zeros(5,5);
    V(1,1) = -1; V(2,1) = 1;                  % PDE      : camp -> amp
    V(2,2) = -1; V(3,2) = -1; V(4,2) = 1;     % EAT      : amp + tgt -> w
    V(4,3) = -1; V(3,3) = 1;                  % GIVEBACK : w -> tgt
    V(2,4) = -1; V(5,4) = 1;                  % LONG1    : amp -> v
    V(5,5) = -1; V(3,5) = 1;                  % LONG2    : v -> tgt
    m4.S = V; m4.lb = zeros(5,1); m4.ub = 1000*ones(5,1);

    t4.metRoles = {'m_in'};
    c4.sig = 3; c4.met.m_in = 3;
    out4.circuits = c4; out4.topology = t4;
    out4.opt = struct('currency',{{'amp'}});

    tmp4 = [tempname '.mat']; txt4 = [tempname '.txt'];
    out = out4; mdl = m4; save(tmp4, 'out', 'mdl', '-v7');
    carr = A_Block0(tmp4, 1, struct('cAMP','camp[c]','verbose',false, ...
                  'allowThrough',{{'amp'}}, 'outFile',txt4));
    rep4 = fileread(txt4);
    delete(tmp4); delete(txt4);

    assert(carr.clean, 'the honest route should still be found');
    assert(isequal(carr.route{1}(:)', {'PDE','LONG1','LONG2'}), ...
        'the carrier route (PDE->EAT->GIVEBACK) must be rejected, not returned');
    assert(isempty(strfind(rep4, 'EAT')), ...
        'no reaction consuming the target may appear in the chosen route');
    assert(~isempty(strfind(rep4, 'net tgt[c]: +1')), ...
        'the net stoichiometry line must show positive target production');

    % --- scenario 5: natural mode (descriptive pathway, currency allowed) ---
    % Mirrors the real case: camp's only fate is amp (currency), and the route
    % onward to tgt also brushes a circuit signal. In natural mode the route
    % must be RETURNED as the answer, with the collision reported but not
    % treated as disqualifying.
    m5.mets     = {'camp[c]';'amp[c]';'prpp[c]';'tgt[c]';'outsig[c]'};
    m5.metNames = {'cAMP';'AMP';'PRPP';'Target';'OutSignal'};
    m5.rxns     = {'PDE1';'APRTrev';'PPAT'};
    m5.rxnNames = {'phosphodiesterase';'aprt reverse';'amidotransferase'};
    m5.rxnECNumbers = {'3.1.4.17';'2.4.2.7';'2.4.2.14'};
    W = zeros(5,3);
    W(1,1) = -1; W(2,1) = 1;                 % PDE1    : camp -> amp
    W(2,2) = -1; W(3,2) = 1; W(5,2) = 1;     % APRTrev : amp -> prpp + outsig
    W(3,3) = -1; W(4,3) = 1;                 % PPAT    : prpp -> tgt
    m5.S = W; m5.lb = zeros(3,1); m5.ub = 1000*ones(3,1);

    t5.metRoles = {'m_in','m_out'};
    c5.sig = [4 5]; c5.met.m_in = 4;
    out5.circuits = c5; out5.topology = t5;
    out5.opt = struct('currency',{{'amp'}});

    tmp5 = [tempname '.mat']; txt5 = [tempname '.txt'];
    out = out5; mdl = m5; save(tmp5, 'out', 'mdl', '-v7');

    strictRun = A_Block0(tmp5, 1, struct('verbose',false,'outFile',[tempname '.txt']));
    natRun = A_Block0(tmp5, 1, struct('verbose',false,'natural',true,'outFile',txt5));
    rep5 = fileread(txt5);
    delete(tmp5); delete(txt5);

    assert(~strictRun.reachable, ...
        'without natural mode, amp is currency and camp is a dead end');
    assert(natRun.reachable, 'natural mode must find the 3-step in-cell route');
    assert(~isempty(strfind(rep5, 'NATURAL ROUTE')) && ...
           ~isempty(strfind(rep5, '<-- THE ANSWER')), ...
        'natural mode must present the route as the answer, not as a reject');
    assert(~isempty(strfind(rep5, 'PPAT')) && ~isempty(strfind(rep5, '2.4.2.14')), ...
        'each step of the natural route needs its enzyme detail printed');
    assert(~isempty(strfind(rep5, 'WHERE THE NATURAL ROUTE MEETS THE CIRCUIT')), ...
        'the signal collision must be reported as informational');

    fprintf(['\nselftest PASSED (5 scenarios).\n' ...
             '  1) strict search fails; one tolerated signal admits a route; the\n' ...
             '     offending signal is named; and m_out is flagged as a short\n' ...
             '     circuit for sitting nearer to cAMP than m_in.\n' ...
             '  2) with amp treated as currency, cAMP is a dead end at step 1;\n' ...
             '     allowing it through reveals the 2-step route. This is the bug\n' ...
             '     behind the "no collision-free connector" verdict.\n' ...
             '  3) a short signal-touching route is rejected in favour of a longer\n' ...
             '     clean one, the clean route is printed, and section 2b names the\n' ...
             '     reaction and the gate that disqualified the short one.\n' ...
             '  4) a route that consumes the target mid-way and hands it back at the\n' ...
             '     end (net zero m_in - the carrier trap) is rejected, and the net\n' ...
             '     stoichiometry of the chosen route is reported.\n' ...
             '  5) natural mode returns the in-cell route through the currency pool\n' ...
             '     as THE ANSWER, with enzyme detail per step, and reports the\n' ...
             '     signal collision as informational rather than disqualifying.\n\n']);
end
