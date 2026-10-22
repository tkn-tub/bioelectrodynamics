function R = A_Block0_Catalytic(stage2, stage1, n, opts)
% A_BLOCK0_CATALYTIC  Find a catalytic Block 0: cAMP GATES a reaction that
% makes m_in, instead of being consumed into it.
%
%   A_Block0_Catalytic                                    % catalytic study, circuit 1
%   A_Block0_Catalytic('IFFL_catalytic_stage2.mat')
%   A_Block0_Catalytic('IFFL_catalytic_stage2.mat','iffl_stage1.mat',1)
%   A_Block0_Catalytic('selftest')
%
% WHY A_Block0 CANNOT ANSWER THIS
%   cAMP's only consumer in the model is PDE1 (camp + h2o -> amp + h), and
%   both products are currency metabolites. A path through a currency pool
%   carries none of the input's atoms, so no mass-flow route from cAMP
%   transmits anything - for ANY choice of m_in. Block 0 has to be
%   regulatory, exactly like r_A_2 in the catalytic interstage: cAMP sits at
%   an allosteric site, is not consumed, and modulates the rate of a reaction
%   that produces m_in. cAMP is a textbook allosteric effector, so this is
%   the natural reading of the block anyway.
%
% WHAT IT DOES
%   1. Lists every reaction cAMP regulates in this model, per BRENDA, with
%      the sense (activator / inhibitor) and whether the reaction is usable.
%   2. DIRECT candidates: a regulated reaction that PRODUCES m_in and touches
%      no other circuit signal. This is the ideal Block 0 - one stage, no
%      adapter chain.
%   3. INDIRECT candidates: a regulated reaction whose product reaches m_in
%      by a signal-free, carrier-free path through NON-currency species only.
%   4. If cAMP regulates nothing usable, says so plainly - then Block 0 is
%      impossible and m_in must be fed directly.
%
%   Currency metabolites are NEVER traversed here. That restriction is what
%   keeps the paths real; see the note above.
%
% WRITES <tag>_Block0_Catalytic.txt

    if nargin >= 1 && (ischar(stage2)||isstring(stage2)) && strcmpi(char(stage2),'selftest')
        R = selftest(); return;
    end
    if nargin < 1 || isempty(stage2), stage2 = 'IFFL_catalytic_stage2.mat'; end
    if nargin < 2 || isempty(stage1), stage1 = 'iffl_stage1.mat'; end
    if nargin < 3 || isempty(n), n = 1; end
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'cAMP'),     opts.cAMP = 'camp[c]'; end
    if ~isfield(opts,'maxSteps'), opts.maxSteps = 6; end
    if ~isfield(opts,'maxShow'),  opts.maxShow = 12; end
    if ~isfield(opts,'outFile'),  opts.outFile = ''; end
    if ~isfield(opts,'verbose'),  opts.verbose = true; end

    assert(exist(stage2,'file')==2, 'Checkpoint %s not found.', stage2);
    L = load(stage2);
    assert(isfield(L,'out') && isfield(L,'mdl'), '%s is not a mapping checkpoint.', stage2);
    out = L.out; mdl = L.mdl;
    assert(~isempty(out.circuits), 'That checkpoint holds no circuits.');
    assert(n <= numel(out.circuits), 'Only %d circuit(s) in %s.', numel(out.circuits), stage2);

    assert(exist(stage1,'file')==2, ...
        ['Seed table %s not found. It holds the BRENDA regulator annotations\n' ...
         'this search needs; run A_Master once to produce it.'], stage1);
    S1 = load(stage1);
    assert(isfield(S1,'T'), '%s holds no activator table T.', stage1);
    Tac = S1.T;
    Tin = []; if isfield(S1,'Tin'), Tin = S1.Tin; end

    if isempty(opts.outFile)
        [~,base] = fileparts(stage2);
        base = regexprep(base, '_stage2$', '');
        opts.outFile = [base '_Block0_Catalytic.txt'];
    end
    fid = fopen(opts.outFile,'w'); cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    say = @(varargin) both(fid, opts, varargin{:});

    R = run(say, mdl, out, Tac, Tin, n, opts);
    R.outFile = opts.outFile;
    say('\nWritten to %s\n', opts.outFile);
end

% ===========================================================================
function R = run(say, mdl, out, Tac, Tin, n, opts)
    c = out.circuits(n); roles = out.topology.metRoles;
    sig = c.sig; mIn = c.met.m_in;

    bar = repmat('=',1,70);
    say('%s\nCATALYTIC BLOCK 0  (%s gates a reaction making m_in)\n%s\n', ...
        bar, opts.cAMP, bar);
    say('circuit %d | m_in = %s (%s)\n\n', n, mdl.mets{mIn}, strtrim(mdl.metNames{mIn}));

    src = find(strcmp(mdl.mets, opts.cAMP), 1);
    if isempty(src)
        say('  %s is not in the model.\n', opts.cAMP);
        R = struct(); R.found = false; return;
    end

    % ---- 1. what does cAMP regulate? --------------------------------------
    key  = normLigand(mdl.metNames{src});
    r2i  = containers.Map(mdl.rxns, num2cell(1:numel(mdl.rxns)));
    [jAC, nAC] = regulatedBy(Tac, key, r2i);
    [jIN, nIN] = regulatedBy(Tin, key, r2i);

    say('1. REACTIONS REGULATED BY %s (BRENDA ligand "%s")\n', opts.cAMP, key);
    say('   activator annotations: %d row(s) -> %d reaction(s) in this model\n', ...
        nAC, numel(jAC));
    say('   inhibitor annotations: %d row(s) -> %d reaction(s) in this model\n\n', ...
        nIN, numel(jIN));

    if isempty(jAC) && isempty(jIN)
        say(['   %s regulates NOTHING in this model, so a catalytic Block 0 is\n' ...
             '   impossible too. Taken with the mass-flow result from A_Block0,\n' ...
             '   cAMP cannot drive this network at all. The remaining options are:\n' ...
             '     - feed m_in (%s) directly as the input; in a cell-free reactor\n' ...
             '       the input is whatever you pipette, so Block 0 is not needed;\n' ...
             '     - pick a different input species that DOES regulate something\n' ...
             '       producing m_in (I can rank candidates if you want);\n' ...
             '     - accept that the BRENDA coverage is the limit here, which the\n' ...
             '       ligand-matching diagnostic already showed is substantial.\n'], ...
             opts.cAMP, mdl.mets{mIn});
        R = struct(); R.found = false; R.nAC = nAC; R.nIN = nIN;
        R.direct = []; R.indirect = {}; return;
    end

    usable = usableMask(mdl);
    curr = currencyMask(mdl, out.opt);     % STRICT: never traversed
    forbid = false(numel(mdl.mets),1); forbid(sig) = true; forbid(mIn) = false;

    cands = [reshape(jAC,1,[]) reshape(jIN,1,[])];
    sense = [repmat({'activates'},1,numel(jAC)) repmat({'inhibits'},1,numel(jIN))];

    % ---- 2. direct: the regulated reaction itself makes m_in --------------
    say('2. DIRECT CANDIDATES (regulated reaction produces m_in itself)\n');
    direct = []; dsense = {};
    for q = 1:numel(cands)
        j = cands(q);
        if ~usable(j), continue; end
        if mdl.S(mIn, j) <= 0, continue; end
        other = intersect(find(mdl.S(:,j) ~= 0), sig);
        other = setdiff(other, mIn);
        if ~isempty(other)
            continue;   % would touch another circuit signal
        end
        direct(end+1) = j; dsense{end+1} = sense{q}; %#ok<AGROW>
    end
    if isempty(direct)
        say('   none.\n\n');
    else
        say('   %d found. These are ideal: one stage, no adapter chain.\n\n', numel(direct));
        for q = 1:min(numel(direct), opts.maxShow)
            showRxn(say, mdl, direct(q), sprintf('r_A_0 (%s)', dsense{q}));
            say('        polarity: %s\n', polarityNote(dsense{q}));
            say('\n');
        end
    end

    % ---- 3. indirect: regulated reaction -> clean path -> m_in -----------
    say('3. INDIRECT CANDIDATES (regulated reaction, then a clean path to m_in)\n');
    indirect = {};
    for q = 1:numel(cands)
        j = cands(q);
        if ~usable(j), continue; end
        if any(direct == j), continue; end
        touches = intersect(find(mdl.S(:,j) ~= 0), sig);
        if ~isempty(setdiff(touches, mIn)), continue; end
        ps = find(mdl.S(:,j) > 0);
        ps = ps(~curr(ps));            % a currency product carries no signal
        for p = reshape(ps,1,[])
            [d, rp, mp] = fwdBFS(mdl, p, mIn, opts.maxSteps, curr, forbid);
            if isnan(d), continue; end
            indirect{end+1,1} = struct('rxn',j,'sense',sense{q},'via',p, ...
                                       'steps',d,'rpath',rp,'mpath',mp); %#ok<AGROW>
            break;
        end
    end
    if isempty(indirect)
        say('   none within %d step(s).\n\n', opts.maxSteps);
    else
        st = cellfun(@(s) s.steps, indirect);
        [~,ord] = sort(st); indirect = indirect(ord);
        say('   %d found, shortest first.\n\n', numel(indirect));
        for q = 1:min(numel(indirect), opts.maxShow)
            s = indirect{q};
            say('   [%d] %s %s, then %d step(s) via %s\n', q, ...
                mdl.rxns{s.rxn}, s.sense, s.steps, mdl.mets{s.via});
            showRxn(say, mdl, s.rxn, sprintf('r_A_0 (%s)', s.sense));
            say('        then: %s\n', strjoin(mdl.rxns(s.rpath), ' -> '));
            net = sum(mdl.S(:, s.rpath), 2);
            say('        net m_in over the chain: %+g\n', full(net(mIn)));
            say('        polarity: %s\n\n', polarityNote(s.sense));
        end
    end

    % ---- verdict ----------------------------------------------------------
    say('%s\nVERDICT\n%s\n', bar, bar);
    if ~isempty(direct)
        say(['  Use a DIRECT candidate. Block 0 becomes a single gated reaction:\n' ...
             '  %s loads at the allosteric site of %s, is not consumed, and sets\n' ...
             '  the rate at which m_in is made. Structurally identical to r_A_2, so\n' ...
             '  the topology table needs one new role (r_A_0) and no new signal.\n'], ...
             opts.cAMP, mdl.rxns{direct(1)});
    elseif ~isempty(indirect)
        say(['  No single regulated reaction makes m_in, but %d route(s) reach it\n' ...
             '  through non-currency intermediates. Each added step is another\n' ...
             '  enzyme in the pot and another species to check for crosstalk, so\n' ...
             '  prefer the shortest.\n'], numel(indirect));
    else
        say(['  %s regulates %d reaction(s) here, but none of them makes m_in or\n' ...
             '  reaches it cleanly within %d steps. Either raise opts.maxSteps, or\n' ...
             '  accept that this m_in cannot be driven from %s and feed it directly.\n'], ...
             opts.cAMP, numel(cands), opts.maxSteps, opts.cAMP);
    end

    R = struct(); R.found = ~isempty(direct) || ~isempty(indirect);
    R.nAC = nAC; R.nIN = nIN; R.direct = direct; R.indirect = {indirect};
end

% ---------------------------------------------------------------- helpers
function [js, nrow] = regulatedBy(T, key, r2i)
% reactions whose BRENDA regulator list contains `key`
    js = []; nrow = 0;
    if isempty(T), return; end
    for t = 1:numel(T)
        regs = strsplit(char(T(t).regulators), ';');
        hit = false;
        for g = 1:numel(regs)
            if strcmp(normLigand(regs{g}), key), hit = true; break; end
        end
        if ~hit, continue; end
        nrow = nrow + 1;
        rn = strtrim(char(T(t).rxn));
        if isKey(r2i, rn), js(end+1) = r2i(rn); end %#ok<AGROW>
    end
    js = unique(js);
end

function s = polarityNote(sense)
    if strncmp(sense,'inhibit',7)
        s = ['INVERTING - m_in falls as cAMP rises, so this flips the sign of the ' ...
             'whole chain; re-check the net polarity of the circuit'];
    else
        s = 'non-inverting - m_in rises with cAMP, polarity of the chain unchanged';
    end
end

function s = normLigand(s)
    s = lower(strtrim(char(s)));
    s = regexprep(s, '^((l|d|dl|\(s\)|\(r\)|alpha|beta)-)+', '');
    s = regexprep(s, '[^a-z0-9]+', ' ');
    s = strtrim(s);
end

function curr = currencyMask(mdl, opt)
    base = lower(regexprep(mdl.mets, '\[[a-z0-9]+\]$', ''));
    list = {'h','h2o','pi','ppi','co2','o2','nh4','nh3','hco3','so4','na1','k','cl', ...
            'fe2','fe3','mg2','ca2','zn2','atp','adp','amp','gtp','gdp','utp','udp', ...
            'ctp','nad','nadh','nadp','nadph','coa','fad','fadh2','q10','q10h2'};
    if isstruct(opt) && isfield(opt,'currency')
        list = unique([list, lower(opt.currency)]);
    end
    curr = ismember(base, list);
    if isstruct(opt) && isfield(opt,'currencyNamePatterns')
        nm = lower(mdl.metNames);
        for p = 1:numel(opt.currencyNamePatterns)
            pat = lower(opt.currencyNamePatterns{p});
            for i = 1:numel(nm)
                if ~curr(i) && ~isempty(strfind(nm{i}, pat)), curr(i) = true; end %#ok<STREMP>
            end
        end
    end
end

function u = usableMask(mdl)
    nR = size(mdl.S,2);
    if isfield(mdl,'ub'), u = mdl.ub(:) > 0; else, u = true(nR,1); end
end

function [d, rpath, mpath] = fwdBFS(mdl, src, tgt, maxd, curr, forbid)
% Same contract as A_Block0's: forbidden metabolites block a reaction
% outright, and any reaction CONSUMING the target is struck out so the
% target cannot act as a carrier.
    d = NaN; rpath = []; mpath = [];
    if src == tgt, d = 0; mpath = tgt; return; end
    nR = size(mdl.S,2);
    ub = ones(nR,1); if isfield(mdl,'ub'), ub = mdl.ub(:); end
    bad = false(nR,1);
    if any(forbid), bad = full(any(mdl.S(forbid,:) ~= 0, 1))'; end
    eatsTgt = full(mdl.S(tgt,:) < 0)';
    usable = ~bad & ~eatsTgt & ub > 0;
    nM = size(mdl.S,1);
    dd = nan(nM,1); pr = nan(nM,1); pm = nan(nM,1);
    dd(src) = 0; frontier = src;
    for k = 1:maxd
        nxt = [];
        for m = reshape(frontier,1,[])
            js = find(mdl.S(m,:) < 0 & usable(:)');
            for j = js
                for p = find(mdl.S(:,j) > 0)'
                    if curr(p) || ~isnan(dd(p)), continue; end
                    dd(p) = k; pr(p) = j; pm(p) = m; nxt(end+1) = p; %#ok<AGROW>
                end
            end
        end
        frontier = nxt;
        if ~isnan(dd(tgt)) || isempty(frontier), break; end
    end
    if isnan(dd(tgt)), return; end
    d = dd(tgt); cur = tgt; mpath = tgt;
    while cur ~= src
        rpath(end+1) = pr(cur); cur = pm(cur); mpath(end+1) = cur; %#ok<AGROW>
    end
    rpath = fliplr(rpath); mpath = fliplr(mpath);
end

function showRxn(say, mdl, j, label)
    nm = ''; ec = ''; gn = '';
    if isfield(mdl,'rxnNames') && numel(mdl.rxnNames) >= j
        nm = strtrim(char(mdl.rxnNames{j}));
    end
    if isfield(mdl,'rxnECNumbers') && numel(mdl.rxnECNumbers) >= j
        ec = strtrim(char(mdl.rxnECNumbers{j}));
    end
    if isfield(mdl,'grRules') && numel(mdl.grRules) >= j
        gn = strtrim(char(mdl.grRules{j}));
    end
    say('     [%s] %s', label, mdl.rxns{j});
    if ~isempty(nm), say('  %s', nm); end
    say('\n        %s\n', eqn(mdl, j));
    if ~isempty(ec), say('        EC: %s\n', ec); end
    if ~isempty(gn), say('        Gene(s): %s\n', gn); end
end

function s = eqn(mdl, j)
    s = sprintf('%s  ->  %s', side(mdl, j, find(mdl.S(:,j) < 0)), ...
                              side(mdl, j, find(mdl.S(:,j) > 0)));
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

function both(fid, opts, fmt, varargin)
    if opts.verbose, fprintf(fmt, varargin{:}); end
    fprintf(fid, fmt, varargin{:});
end

% ================================ selftest ==================================
function R = selftest()
% MAKEIN is activated by cAMP and makes m_in directly  -> direct candidate
% FARAWAY is activated by cAMP and makes w, w -> m_in  -> indirect candidate
% DIRTY   is activated by cAMP but also eats a gate    -> must be rejected
    mdl.mets     = {'camp[c]';'tgt[c]';'w[c]';'gate[c]';'p[c]';'amp[c]'};
    mdl.metNames = {'cAMP';'Target';'W';'GateG';'P';'AMP'};
    mdl.rxns     = {'MAKEIN';'FARAWAY';'DIRTY';'WTOTGT';'PDE1'};
    mdl.rxnNames = {'makes the input';'makes w';'dirty maker';'w to target';'pde'};
    mdl.rxnECNumbers = {'1.1.1.1';'2.2.2.2';'3.3.3.3';'4.4.4.4';'3.1.4.17'};
    S = zeros(6,5);
    S(5,1) = -1; S(2,1) = 1;                 % MAKEIN  : p -> tgt
    S(5,2) = -1; S(3,2) = 1;                 % FARAWAY : p -> w
    S(5,3) = -1; S(4,3) = -1; S(2,3) = 1;    % DIRTY   : p + gate -> tgt
    S(3,4) = -1; S(2,4) = 1;                 % WTOTGT  : w -> tgt
    S(1,5) = -1; S(6,5) = 1;                 % PDE1    : camp -> amp
    mdl.S = S; mdl.lb = zeros(5,1); mdl.ub = 1000*ones(5,1);

    topo.metRoles = {'m_in','m_Fd_1'};
    c.sig = [2 4]; c.met.m_in = 2;
    out.circuits = c; out.topology = topo;
    out.opt = struct('currency',{{'amp'}});

    T = struct('rxn', {'MAKEIN','FARAWAY','DIRTY'}, ...
               'regulators', {'cAMP','cAMP','cAMP'});
    Tin = struct('rxn', {}, 'regulators', {});

    t2 = [tempname '.mat']; t1 = [tempname '.mat']; txt = [tempname '.txt'];
    save(t2, 'out', 'mdl', '-v7');
    save(t1, 'T', 'Tin', '-v7');
    R = A_Block0_Catalytic(t2, t1, 1, struct('verbose',false,'outFile',txt));
    rep = fileread(txt);
    delete(t2); delete(t1); delete(txt);

    assert(R.found, 'a catalytic Block 0 should be found');
    assert(isequal(R.direct, 1), 'MAKEIN alone should be the direct candidate');
    assert(isempty(strfind(rep, '[r_A_0 (activates)] DIRTY')), ...
        'DIRTY eats a gate and must not be offered');
    ind = R.indirect{1};
    assert(numel(ind) == 1 && ind{1}.rxn == 2, 'FARAWAY should be the indirect one');
    assert(ind{1}.steps == 1, 'FARAWAY reaches tgt in one further step');
    assert(R.nAC == 3, 'three activator rows mention cAMP');

    fprintf(['\nselftest PASSED.\n' ...
             '  - a cAMP-activated reaction producing m_in is found as DIRECT\n' ...
             '  - one that reaches m_in via a non-currency intermediate is INDIRECT\n' ...
             '  - one that also consumes a gate is rejected\n' ...
             '  - activator rows are counted even when several map to one reaction\n\n']);
end
