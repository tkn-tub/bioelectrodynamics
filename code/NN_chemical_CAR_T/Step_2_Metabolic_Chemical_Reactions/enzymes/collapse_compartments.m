function [mdl, rep] = collapse_compartments(model, opts)
% COLLAPSE_COMPARTMENTS  Turn a compartmentalised GEM into a genuine single-pot
% model, for a cell-free reactor where compartments do not exist.
%
%   [mdl, rep] = collapse_compartments(model)
%   collapse_compartments('selftest')
%
% WHY NOT JUST SET opt.compartment = 'any'
%   The compartment tag is on the METABOLITES, not only the reactions. akg[c]
%   and akg[m] are the same chemical. Allowing any compartment without merging
%   lets the search treat them as two distinct signals, so it can build a
%   "circuit" whose fast arm uses akg[c] and whose slow arm uses akg[m] - one
%   molecule, in one pot, counted twice. Signal exclusivity cannot catch that,
%   because as far as the index arithmetic is concerned they are two species.
%
% WHAT THIS DOES
%   1. Merges every metabolite sharing a base id into one species, and gives it
%      the suffix [c] so downstream code needing a compartment still works.
%   2. Re-sums each reaction's stoichiometry over the merged species.
%   3. Drops reactions that become empty or null (A -> A). Those are exactly
%      the transporters: moving a species between compartments is a no-op once
%      the compartments are gone.
%   4. Collapses reactions that become stoichiometrically identical, keeping
%      one and recording how many merged.
%
% WHAT IT COSTS YOU
%   Genuine biology is lost along with the artefacts. A proton gradient across
%   a membrane, or the same reaction running with different kinetics in two
%   compartments, cannot be represented. For a single-pot reactor that is the
%   correct simplification; for anything claiming to describe an intact cell it
%   is not. Say which you are doing in the paper.

    if nargin == 1 && (ischar(model) || isstring(model)) && strcmpi(char(model),'selftest')
        [mdl, rep] = selftest(); return;
    end
    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'verbose'),   opts.verbose = true; end
    if ~isfield(opts,'keepFirst'), opts.keepFirst = true; end   % name/EC from first copy

    nM = numel(model.mets); nR = numel(model.rxns);

    % ---- 1. merge metabolites by base id -----------------------------------
    base = regexprep(model.mets, '\[[a-z0-9]+\]$', '');
    % built by hand: unique(...,'stable') does not return the third output in
    % every MATLAB/Octave build, and this is the index arithmetic everything
    % else depends on.
    seenMet = containers.Map('KeyType','char','ValueType','double');
    map = zeros(nM,1); uBase = {}; firstIx = [];
    for i = 1:nM
        b = base{i};
        if isKey(seenMet, b)
            map(i) = seenMet(b);
        else
            uBase{end+1} = b; firstIx(end+1) = i; %#ok<AGROW>
            seenMet(b) = numel(uBase);
            map(i) = numel(uBase);
        end
    end
    nM2 = numel(uBase);

    mdl = struct();
    mdl.mets = cellfun(@(b) [b '[c]'], uBase(:), 'uni', 0);
    mdl.metNames = model.metNames(firstIx(:));
    for f = {'metFormulas','metCHEBIID','metKEGGID','metHMDBID','metPubChemID','metInChIString'}
        if isfield(model, f{1}), mdl.(f{1}) = model.(f{1})(firstIx); end
    end

    % ---- 2. re-sum stoichiometry over merged species -----------------------
    S = model.S;
    A = sparse(map(:), (1:nM)', 1, nM2, nM);           % merge operator
    S2 = A * S;                                        % nM2 x nR

    % ---- 3. drop null reactions (pure transport) and empties ---------------
    nnzPerRxn = sum(S2 ~= 0, 1);
    keep = nnzPerRxn > 1;                 % a null or single-species column is not a conversion
    nNull = sum(~keep);

    % ---- 4. collapse reactions that became identical -----------------------
    idx = find(keep);
    sig = cell(1, numel(idx));
    for q = 1:numel(idx)
        j = idx(q);
        r = find(S2(:,j) ~= 0);
        v = full(S2(r,j));
        sig{q} = sprintf('%d:%g;', [r(:)'; v(:)']);
    end
    seenRxn = containers.Map('KeyType','char','ValueType','double');
    keepIdx = []; mergedCount = [];
    for q = 1:numel(idx)
        k = sig{q};
        if isKey(seenRxn, k)
            w = seenRxn(k); mergedCount(w) = mergedCount(w) + 1;
        else
            keepIdx(end+1) = idx(q); mergedCount(end+1) = 1; %#ok<AGROW>
            seenRxn(k) = numel(keepIdx);
        end
    end
    nDup = numel(idx) - numel(keepIdx);

    mdl.S = S2(:, keepIdx);
    mdl.rxns = model.rxns(keepIdx);
    for f = {'rxnNames','rxnECNumbers','subSystems','grRules'}
        if isfield(model, f{1}), mdl.(f{1}) = model.(f{1})(keepIdx); end
    end
    if isfield(model,'lb'), mdl.lb = model.lb(keepIdx); else, mdl.lb = zeros(numel(keepIdx),1); end
    if isfield(model,'ub'), mdl.ub = model.ub(keepIdx); else, mdl.ub = 1000*ones(numel(keepIdx),1); end
    if isfield(model,'rxnGeneMat') && ~isempty(model.rxnGeneMat)
        mdl.rxnGeneMat = model.rxnGeneMat(keepIdx, :);
    end
    if isfield(model,'genes'), mdl.genes = model.genes; end
    mdl.rxnMergedCount = mergedCount(:);

    % ---- 5. report ---------------------------------------------------------
    rep.nMetBefore = nM;  rep.nMetAfter = nM2;
    rep.nRxnBefore = nR;  rep.nRxnAfter = numel(keepIdx);
    rep.nNullDropped = nNull; rep.nDuplicateMerged = nDup;
    cmp = regexp(model.mets, '\[([a-z0-9]+)\]$', 'tokens', 'once');
    cmp = cellfun(@(c) char([c{:}]), cmp, 'uni', 0);
    rep.compartments = unique(cmp(~cellfun(@isempty, cmp)));

    if opts.verbose
        fprintf('\n--- collapse_compartments ------------------------------------\n');
        fprintf('  compartments merged ......... %s\n', strjoin(rep.compartments, ', '));
        fprintf('  metabolites ................. %d -> %d  (%d duplicates merged)\n', ...
            nM, nM2, nM - nM2);
        fprintf('  reactions ................... %d -> %d\n', nR, rep.nRxnAfter);
        fprintf('    dropped as null/transport . %d\n', nNull);
        fprintf('    merged as now-identical ... %d\n', nDup);
        fprintf(['  Every metabolite now carries [c]; the model looks like a one-pot\n' ...
                 '  cytosolic system, so opt.compartment = ''c'' remains correct.\n']);
    end
end

% ================================ selftest ==================================
function [mdl, rep] = selftest()
    % two compartments, one shared species, one transporter, one duplicated rxn
    model.mets     = {'a[c]';'b[c]';'akg[c]';'akg[m]';'d[m]';'h[c]'};
    model.metNames = {'A';'B';'2-Oxoglutarate';'2-Oxoglutarate';'D';'H+'};
    model.rxns     = {'R1';'R2';'T1';'R3';'R4'};
    model.rxnNames = model.rxns;
    model.rxnECNumbers = {'1.1.1.1';'2.2.2.2';'';'3.3.3.3';'1.1.1.1'};
    S = zeros(6,5);
    S(1,1) = -1; S(3,1) = 1;        % R1 : a[c] -> akg[c]
    S(2,2) = -1; S(4,2) = 1;        % R2 : b[c] -> akg[m]
    S(3,3) = -1; S(4,3) = 1;        % T1 : akg[c] -> akg[m]   (pure transport)
    S(4,4) = -1; S(5,4) = 1;        % R3 : akg[m] -> d[m]
    S(1,5) = -1; S(3,5) = 1;        % R4 : a[c] -> akg[c]     (duplicate of R1)
    model.S = S;
    model.lb = zeros(5,1); model.ub = 1000*ones(5,1);

    [mdl, rep] = collapse_compartments(model, struct('verbose', true));

    assert(rep.nMetAfter == 5, 'akg[c] and akg[m] should merge to one species');
    assert(any(strcmp(mdl.mets, 'akg[c]')), 'merged species should carry [c]');
    assert(~any(strcmp(mdl.rxns, 'T1')), 'the transporter should vanish as a null reaction');
    assert(rep.nNullDropped == 1, 'exactly one null reaction expected');
    assert(rep.nDuplicateMerged == 1, 'R4 should merge into R1');
    assert(rep.nRxnAfter == 3, 'R1(+R4), R2, R3 should remain');
    % R2 now also makes the merged akg, so it is a genuine second producer
    j = find(strcmp(mdl.rxns,'R2'));
    k = find(strcmp(mdl.mets,'akg[c]'));
    assert(mdl.S(k, j) > 0, 'R2 should produce the merged 2-oxoglutarate');
    fprintf('\nselftest PASSED: compartment duplicates merged, transport removed,\n');
    fprintf('identical reactions collapsed, and the merged species is shared.\n\n');
end
