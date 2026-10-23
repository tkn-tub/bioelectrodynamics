function R = block0_route(mdl, o)
% BLOCK0_ROUTE  The one Block 0 path search, shared by A_Block0.m (the detailed
% analysis report) and export_iffl_chain.m (the BLOCK 0 section of the chain
% file). Having a single implementation is the point: the two used to run
% separate searches and could therefore disagree, which is exactly what
% happened - the chain file reported "no connector" while the Block 0 report
% found a 3-step route.
%
%   R = block0_route(mdl, o)
%   block0_route('selftest')
%
% o fields (all optional unless marked):
%   src       REQUIRED  source metabolite index (cAMP)
%   tgt       REQUIRED  target metabolite index (m_in)
%   opt                 the study's opt struct, for opt.currency and
%                       opt.currencyNamePatterns. Used only to build `curr`.
%   curr                logical mask of untraversable currency metabolites.
%                       Supply it to override what `opt` would produce.
%   natural   false     true un-marks the adenylate pool so the route may hop
%                       through it, as flux in a real cell does. See WHY below.
%   allowThrough        base names to un-mark; defaults to the adenylate pool
%                       when natural is true, {} otherwise.
%   forbid              logical mask of metabolites no reaction may touch.
%                       Defaults to none. The target is force-allowed.
%   sig                 circuit signal indices, for collision reporting.
%   roles               role names parallel to sig.
%   maxSteps  14        search depth.
%
% R fields:
%   steps     NaN if unreachable, else path length
%   rpath     reaction indices in order
%   mpath     metabolite indices in order, starting at src
%   curr      the currency mask actually used
%   net       net stoichiometry over rpath at unit flux (sparse column)
%   netTgt    net(tgt), which MUST be > 0 for a usable route
%   collisions struct array: .rxn .met .role .verb, one per signal touch
%
% WHY `natural` EXISTS
%   cAMP's only consumer is PDE1 (camp + h2o -> amp + h) and both products are
%   currency metabolites. With currency untraversable the search dies at step 1.
%   That is the correct answer for an input stage you intend to BUILD: a path
%   through a buffered pool carries none of the input's atoms, and worse, pairs
%   like (atp, oaa) in
%       ACITL: atp + cit + coa -> accoa + adp + oaa + pi
%   are false edges - ATP leaves as ADP + Pi while OAA comes from citrate.
%   But when Block 0 is DESCRIPTIVE - the pathway a cell would use, not
%   something loaded into a pot - those hops are legitimate, because in a cell
%   that flux genuinely occurs. Hence a flag rather than a fixed rule.
%
% THE CARRIER RULE
%   The target is deliberately exempt from `forbid`, so the last reaction may
%   produce it. That exemption must NOT extend to consuming it: a route where
%   an intermediate eats the target and a later step hands it back makes the
%   target a CARRIER with zero net product. Every reaction consuming the target
%   is therefore struck out, which also guarantees netTgt >= 1.

    if nargin >= 1 && (ischar(mdl)||isstring(mdl)) && strcmpi(char(mdl),'selftest')
        R = selftest(); return;
    end

    assert(isstruct(o) && isfield(o,'src') && isfield(o,'tgt'), ...
        'block0_route needs o.src and o.tgt.');
    if ~isfield(o,'natural'),  o.natural = false; end
    if ~isfield(o,'maxSteps'), o.maxSteps = 14; end
    if ~isfield(o,'sig'),      o.sig = []; end
    if ~isfield(o,'roles'),    o.roles = {}; end
    if ~isfield(o,'allowThrough')
        if o.natural, o.allowThrough = {'amp','adp','atp','adn','ppi'};
        else,         o.allowThrough = {}; end
    end
    if ~isfield(o,'opt'), o.opt = struct(); end

    nM = size(mdl.S,1);
    if isfield(o,'curr') && ~isempty(o.curr)
        curr = logical(o.curr(:));
    else
        curr = currencyMask(mdl, o.opt, o.allowThrough);
    end
    if isfield(o,'forbid') && ~isempty(o.forbid)
        forbid = logical(o.forbid(:));
    else
        forbid = false(nM,1);
    end
    forbid(o.tgt) = false;          % the last step has to be allowed to make it

    [d, rpath, mpath] = bfs(mdl, o.src, o.tgt, o.maxSteps, curr, forbid);

    R = struct();
    R.steps = d; R.rpath = rpath; R.mpath = mpath; R.curr = curr;
    R.collisions = struct('rxn',{},'met',{},'role',{},'verb',{});
    if isnan(d)
        R.net = sparse(nM,1); R.netTgt = 0; return;
    end
    R.net = sum(mdl.S(:, rpath), 2);
    R.netTgt = full(R.net(o.tgt));

    % which circuit signals does the route touch, and how?
    if ~isempty(o.sig)
        for j = reshape(rpath,1,[])
            hits = setdiff(intersect(find(mdl.S(:,j) ~= 0), o.sig), o.tgt);
            for m = reshape(hits,1,[])
                k = find(o.sig == m, 1);
                if ~isempty(o.roles) && numel(o.roles) >= k
                    rl = o.roles{k};
                else
                    rl = sprintf('signal %d', k);
                end
                if mdl.S(m,j) < 0, vb = 'consumes'; else, vb = 'produces'; end
                R.collisions(end+1) = struct('rxn',j,'met',m,'role',rl,'verb',vb);
            end
        end
    end
end

% ---------------------------------------------------------------- internals
function [d, rpath, mpath] = bfs(mdl, src, tgt, maxd, curr, forbid)
    d = NaN; rpath = []; mpath = [];
    if src == tgt, d = 0; mpath = tgt; return; end
    nR = size(mdl.S,2);
    ub = ones(nR,1); if isfield(mdl,'ub'), ub = mdl.ub(:); end
    bad = false(nR,1);
    if any(forbid), bad = full(any(mdl.S(forbid,:) ~= 0, 1))'; end
    eatsTgt = full(mdl.S(tgt,:) < 0)';      % the carrier rule; see the header
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

function curr = currencyMask(mdl, opt, allowThrough)
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
    if ~isempty(allowThrough)
        curr(ismember(base, lower(allowThrough))) = false;
    end
end

% ================================ selftest ==================================
function R = selftest()
% The real case in miniature: camp -> amp (currency) -> prpp -> tgt, where the
% middle step also makes m_out and the last also eats a gate. Plus a carrier
% trap that must be refused.
    mdl.mets     = {'camp[c]';'amp[c]';'prpp[c]';'tgt[c]';'outsig[c]';'gate[c]';'w[c]'};
    mdl.metNames = {'cAMP';'AMP';'PRPP';'Target';'Adenine';'Gln';'W'};
    mdl.rxns     = {'PDE1';'APRTrev';'PPAT';'EAT';'GIVEBACK'};
    mdl.rxnNames = mdl.rxns;
    S = zeros(7,5);
    S(1,1) = -1; S(2,1) = 1;                            % PDE1     : camp -> amp
    S(2,2) = -1; S(3,2) = 1; S(5,2) = 1;                % APRTrev  : amp -> prpp + outsig
    S(3,3) = -1; S(6,3) = -1; S(4,3) = 1;               % PPAT     : prpp + gate -> tgt
    S(2,4) = -1; S(4,4) = -1; S(7,4) = 1;               % EAT      : amp + tgt -> w
    S(7,5) = -1; S(4,5) = 1;                            % GIVEBACK : w -> tgt
    mdl.S = S; mdl.lb = zeros(5,1); mdl.ub = 1000*ones(5,1);

    opt = struct('currency', {{'amp'}});
    sig = [4 5 6]; roles = {'m_in','m_out','m_Fd_1'};

    % strict: amp is currency, so camp is a dead end
    a = block0_route(mdl, struct('src',1,'tgt',4,'opt',opt,'sig',sig,'roles',{roles}));
    assert(isnan(a.steps), 'strict mode must not traverse the currency pool');

    % natural: the 3-step route, with both collisions reported
    b = block0_route(mdl, struct('src',1,'tgt',4,'opt',opt,'natural',true, ...
                                 'sig',sig,'roles',{roles}));
    assert(b.steps == 3, 'natural mode should find camp -> amp -> prpp -> tgt');
    assert(isequal(mdl.rxns(b.rpath)', {'PDE1','APRTrev','PPAT'}), 'wrong route');
    assert(b.netTgt == 1, 'net target production must be +1');
    assert(numel(b.collisions) == 2, 'two signal collisions expected');
    vs = sort({b.collisions.verb});
    assert(isequal(vs, {'consumes','produces'}), 'one consume and one produce');
    rl = sort({b.collisions.role});
    assert(isequal(rl, {'m_Fd_1','m_out'}), 'collisions should name m_out and m_Fd_1');

    % the carrier trap (EAT eats tgt, GIVEBACK returns it) must never be used
    assert(~any(strcmp(mdl.rxns(b.rpath), 'EAT')), ...
        'no reaction consuming the target may appear in a route');

    % forbidding the gate kills the only natural route
    fb = false(7,1); fb(6) = true;
    c = block0_route(mdl, struct('src',1,'tgt',4,'opt',opt,'natural',true, ...
                                 'forbid',fb,'sig',sig,'roles',{roles}));
    assert(isnan(c.steps), 'forbidding the gate should block the route');

    R = b;
    fprintf(['\nblock0_route selftest PASSED.\n' ...
             '  - strict mode treats the currency pool as untraversable\n' ...
             '  - natural mode finds the 3-step in-cell route, netTgt = +1\n' ...
             '  - both signal collisions reported with the right verbs\n' ...
             '  - the carrier trap is refused\n' ...
             '  - forbid is honoured\n\n']);
end
