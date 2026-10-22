function topo = iffl_topology(consumingN, fastSteps, miMode, catalytic, gateMode, neuronGate, interGate)
% IFFL_TOPOLOGY  Canonical role names, wiring and anti-crosstalk table for the
% two-differentiator -> counter -> neuron circuit.
%
%   topo = iffl_topology()        % fast arm of two reactions (the drawing)
%   topo = iffl_topology(false,1) % fast arm of ONE reaction
%
% NOTATION  <role>_<block>, block 1 and block 2:
%
%   m_in --r_Ff_1--> m_Ff_1 --r_Fd_1--> m_Fd_1 --r_Fr_1--> 0
%                                        |
%                                        | catalytic (gate, NOT consumed)
%                                        v
%   m_in --r_Sf_1--> m_Sf_1 --r_Sd_1--> m_d1 --r_I_1--> 0
%                                        |  \--r_lk_1--> 0
%                                        \--> block 2 (m_d1 is its m_in)
%
%   block 2 repeats the same roles with the _2 suffix and outputs m_d2.
%   counter:  r_c : m_d2 -> m_c ;  r_I1 : m_c + m_I -> 0 ;  r_I2 : m_d2 + m_I -> 0
%   neuron :  r_Nw : -> m_out (gated by m_c) ; r_Nb : m_b -> m_out ; r_Nr : m_out -> 0
%
% fastSteps = 1 collapses the fast arm to r_Fd_k : m_in -> m_Fd_k, dropping
% m_Ff_k and r_Ff_k. The loop is still incoherent and feed-forward; the fast
% arm is simply one conversion shorter, which is what "fast" wanted anyway. It
% also removes two species per stage from the distinctness requirement, so use
% it when the cascade keeps failing because D1 and D2 claim the same pool.
%
% THE ANTI-CROSSTALK RULE lives here and nowhere else. `subs` and `prods` list
% the ONLY signal species each reaction may touch stoichiometrically. Any other
% contact between a circuit reaction and a circuit signal is crosstalk and the
% search rejects the candidate: it is what forbids r_I_1 producing m_c, r_c
% consuming m_d1, or r_Fr_1 consuming m_c.
%
% Catalytic edges (regulator -> target) are deliberately ABSENT from subs/prods:
% a gate that also appears in the stoichiometry is being consumed, which is the
% one thing the IFFL primitive must not do.

    if nargin < 1 || isempty(consumingN), consumingN = false; end
    if nargin < 2 || isempty(fastSteps),  fastSteps  = 2; end
    if nargin < 3 || isempty(miMode),     miMode     = 'shared'; end
    if islogical(miMode)                      % backward compatible
        if miMode, miMode = 'shared'; else, miMode = 'separate'; end
    end
    assert(any(strcmp(miMode,{'shared','separate','none'})), ...
        'miMode must be shared, separate or none.');
    assert(any(fastSteps == [1 2]), 'fastSteps must be 1 or 2.');
    if nargin < 4 || isempty(catalytic), catalytic = false; end
    % gateMode{k} says how stage k's gate acts:
    %   'activate' - m_Fd_k ACTIVATES r_I_k, which degrades m_d_k (original)
    %   'inhibit'  - m_Fd_k INHIBITS r_Sd_k, the reaction that PRODUCES m_d_k.
    %                No gated degrader is needed; m_d_k falls because its
    %                production stops and r_lk_k drains it. Still an incoherent
    %                feed-forward loop, and m_Fd_k is still cleared by r_Fr_k.
    % A hybrid {'activate','inhibit'} draws the two stages' gates from DIFFERENT
    % BRENDA fields, which is what stops both stages claiming the same gate.
    if nargin < 5 || isempty(gateMode), gateMode = {'activate','activate'}; end
    if ischar(gateMode) || isstring(gateMode), gateMode = {char(gateMode), char(gateMode)}; end
    assert(numel(gateMode) == 2, 'gateMode needs one entry per stage.');
    for q = 1:2
        assert(any(strcmp(gateMode{q}, {'activate','inhibit'})), ...
            'gateMode entries must be activate or inhibit.');
    end

    % neuronGate: how m_c drives the neuron weight.
    %   'activate' - m_c activates r_Nw, giving a POSITIVE weight w
    %   'inhibit'  - m_c inhibits r_Nw, giving a NEGATIVE weight. A single
    %                neuron with a negative weight is perfectly ordinary; it
    %                just means the output falls as the count rises.
    %   'any'      - accept either, and report which one was found
    if nargin < 6 || isempty(neuronGate), neuronGate = 'activate'; end
    assert(any(strcmp(neuronGate, {'activate','inhibit','any'})), ...
        'neuronGate must be activate, inhibit or any.');
    topo.neuronGate = neuronGate;
    % interGate: how m_d1 drives the interstage r_A_2.
    %   'activate' - m_d1 activates r_A_2, so m_in_2 RISES with m_d1 (gain +1)
    %   'inhibit'  - m_d1 inhibits r_A_2, so m_in_2 FALLS with m_d1 (gain -1),
    %                which flips the sign of the whole cascade
    %   'any'      - accept either; the sign is computed per circuit and
    %                reported, never assumed
    if nargin < 7 || isempty(interGate), interGate = 'activate'; end
    assert(any(strcmp(interGate, {'activate','inhibit','any'})), ...
        'interGate must be activate, inhibit or any.');
    topo.interGate = interGate;
    topo.gateMode  = gateMode;
    topo.catalytic = catalytic;
    topo.fastSteps = fastSteps;
    topo.miMode    = miMode;
    topo.sharedMI  = strcmp(miMode,'shared');
    topo.metRoles = {'m_in'};
    topo.rxnRoles = {};
    topo.subs = {}; topo.prods = {};
    % CATALYTIC INTERSTAGE. With catalytic = false the cascade is mass-flow:
    % D2's arms CONSUME m_d1, which forces D2 to start where D1 ended and pins
    % both stages into the same neighbourhood of the network.
    % With catalytic = true, m_d1 instead ACTIVATES r_A_2, which makes D2's own
    % input m_in_2 from an independent precursor. m_d1 is read, never drawn on -
    % the high-impedance interstage - and D2 is free to sit anywhere.
    if catalytic
        inRole = {'m_in','m_in_2'};
    else
        inRole = {'m_in','m_d1'};
    end

    for k = 1:2
        if k == 2 && catalytic
            topo.metRoles = [topo.metRoles, {'m_in_2'}];
            topo.rxnRoles = [topo.rxnRoles, {'r_A_2'}];
            topo.subs  = [topo.subs,  {{}}];        % precursor is not a signal
            topo.prods = [topo.prods, {{'m_in_2'}}];
        end
        sfx = sprintf('_%d', k);
        if k == 1, dRole = 'm_d1'; else, dRole = 'm_d2'; end
        FdRole = ['m_Fd' sfx]; SfRole = ['m_Sf' sfx]; FfRole = ['m_Ff' sfx];
        gated = strcmp(gateMode{k}, 'activate');   % false => no r_I on this stage
        if fastSteps == 2
            topo.metRoles = [topo.metRoles, {FfRole, FdRole, SfRole, dRole}];
            rr = {['r_Ff' sfx], ['r_Fd' sfx], ['r_Fr' sfx], ['r_Sf' sfx], ['r_Sd' sfx]};
            ss = {{inRole{k}}, {FfRole}, {FdRole}, {inRole{k}}, {SfRole}};
            pp = {{FfRole},    {FdRole}, {},       {SfRole},    {dRole}};
        else
            topo.metRoles = [topo.metRoles, {FdRole, SfRole, dRole}];
            rr = {['r_Fd' sfx], ['r_Fr' sfx], ['r_Sf' sfx], ['r_Sd' sfx]};
            ss = {{inRole{k}}, {FdRole}, {inRole{k}}, {SfRole}};
            pp = {{FdRole},    {},       {SfRole},    {dRole}};
        end
        if gated
            rr = [rr, {['r_I' sfx]}]; ss = [ss, {{dRole}}]; pp = [pp, {{}}];
        end
        rr = [rr, {['r_lk' sfx]}]; ss = [ss, {{dRole}}]; pp = [pp, {{}}];
        topo.rxnRoles = [topo.rxnRoles, rr];
        topo.subs  = [topo.subs,  ss];
        topo.prods = [topo.prods, pp];
    end

    % m_I is the bimolecular partner of the counter's two removal reactions. In
    % this model most consumers of a lipid m_c are hydrolases whose only other
    % substrate is water, i.e. they are unimolecular in non-currency terms, so
    % miMode = 'none' drops the role: with the partner clamped (water, or an
    % externally replenished pool) the reset is first order in m_c either way,
    % and the rate is programmed by enzyme loading rather than by [m_I].
    if strcmp(miMode,'none')
        topo.metRoles = [topo.metRoles, {'m_c','m_b','m_out'}];
        subI1 = {'m_c'}; subI2 = {'m_d2'};
    else
        topo.metRoles = [topo.metRoles, {'m_c','m_I','m_b','m_out'}];
        subI1 = {'m_c','m_I'}; subI2 = {'m_d2','m_I'};
    end
    topo.rxnRoles = [topo.rxnRoles, {'r_c','r_I1','r_I2','r_Nw','r_Nb','r_Nr'}];
    topo.subs  = [topo.subs,  {{'m_d2'}, subI1, subI2, {}, {'m_b'}, {'m_out'}}];
    topo.prods = [topo.prods, {{'m_c'},  {},    {},    {'m_out'}, {'m_out'}, {}}];

    % consuming_N = true makes the neuron weight reaction draw on m_c (then it is
    % no longer a programmable, non-depleting read).
    if consumingN, topo.subs{ strcmp(topo.rxnRoles,'r_Nw') } = {'m_c'}; end

    % reqSubs = the substrates a role MUST have (checked for sign). subs is the
    % wider set it is ALLOWED to touch. They differ only for r_I2 when the
    % counter's m_I is not required to be shared: r_I2 may still use m_I, it
    % just no longer has to.
    topo.reqSubs = topo.subs;
    if strcmp(miMode,'separate')
        topo.reqSubs{ strcmp(topo.rxnRoles,'r_I2') } = {'m_d2'};
    end

    % the gate regulates r_I when it activates, r_Sd when it inhibits
    topo.regulator = {}; topo.regTarget = {}; topo.regSense = {};
    for k = 1:2
        sfx = sprintf('_%d', k);
        topo.regulator{end+1} = ['m_Fd' sfx];
        if strcmp(gateMode{k}, 'activate')
            topo.regTarget{end+1} = ['r_I' sfx];  topo.regSense{end+1} = 'activates';
        else
            topo.regTarget{end+1} = ['r_Sd' sfx]; topo.regSense{end+1} = 'inhibits';
        end
    end
    % With consumingN the neuron READS m_c by consuming it, so there is no
    % catalytic edge to declare - m_c is r_Nw's substrate instead, and listing
    % it as a regulator would demand both at once.
    if ~consumingN
        topo.regulator{end+1} = 'm_c'; topo.regTarget{end+1} = 'r_Nw';
        switch neuronGate
            case 'inhibit', topo.regSense{end+1} = 'inhibits';
            case 'any',     topo.regSense{end+1} = 'regulates';
            otherwise,      topo.regSense{end+1} = 'activates';
        end
    end
    if catalytic
        topo.regulator = [topo.regulator, {'m_d1'}];
        topo.regTarget = [topo.regTarget, {'r_A_2'}];
        switch interGate
            case 'inhibit', topo.regSense = [topo.regSense, {'inhibits'}];
            case 'any',     topo.regSense = [topo.regSense, {'regulates'}];
            otherwise,      topo.regSense = [topo.regSense, {'activates'}];
        end
    end

    % allowed{j} = indices into metRoles that reaction j may touch
    topo.allowed = cell(size(topo.rxnRoles));
    for j = 1:numel(topo.rxnRoles)
        names = [topo.subs{j}, topo.prods{j}];
        ix = zeros(1, numel(names));
        for q = 1:numel(names), ix(q) = find(strcmp(topo.metRoles, names{q})); end
        topo.allowed{j} = ix;
    end

    topo.block = zeros(size(topo.rxnRoles));
    for j = 1:numel(topo.rxnRoles)
        r = topo.rxnRoles{j};
        if ~isempty(regexp(r,'_1$','once')), topo.block(j) = 1;
        elseif ~isempty(regexp(r,'_2$','once')), topo.block(j) = 2;
        elseif any(strcmp(r,{'r_c','r_I1','r_I2'})), topo.block(j) = 3;
        else, topo.block(j) = 4; end
    end

    % POLARITY. Each IFFL stage is inverting: the fast arm raises the gated
    % degradation of m_d before the slow arm raises its production, so a step up
    % in the stage input drives a NEGATIVE-going excursion of the stage output
    % about its baseline. Two stages in series therefore give a net +1.
    topo.stagePolarity = [-1 -1];
    topo.netPolarity   = prod(topo.stagePolarity);
    if catalytic, link = 'cat'; else, link = 'mass'; end
    gm = [gateMode{1}(1) gateMode{2}(1)];      % 'aa', 'ai', ...
    topo.version = sprintf('iffl_fixed_%dplus2_%sI_%slink_%sgate_N%s_A%s_v7', ...
        fastSteps, miMode, link, gm, ngTag(neuronGate), ngTag(interGate));
end

function tag = ngTag(ng)
    switch ng
        case 'inhibit', tag = 'i';
        case 'any',     tag = 'x';
        otherwise,      tag = 'a';
    end
end
