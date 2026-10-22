function R = diagnose_ligand_matching(model, flatFile, opts)
% DIAGNOSE_LIGAND_MATCHING  How much of the BRENDA activator annotation is lost
% to name matching, and how much of that loss is recoverable?
%
%   R = diagnose_ligand_matching(model, 'brenda_2026_1.txt')
%   diagnose_ligand_matching('selftest')
%
% The mapping pipeline matches BRENDA activator LIGAND NAMES against
% model.metNames after normalisation. That silently discards every activator
% whose name happens to be written differently in the two sources. This script
% quantifies the damage and, crucially, separates it into:
%
%   VARIANT     - a spelling / formatting variant of a compound that IS in the
%                 model. Detectable from the strings alone, so certainly
%                 recoverable by an identifier cross-reference.
%   UNRESOLVED  - no string-detectable counterpart. This is EITHER a compound
%                 genuinely absent from the model (a plant metabolite, a drug)
%                 OR a true synonym with no shared words, such as
%                 'adenosine triphosphate' against a model that writes 'ATP'.
%                 No string method can tell those two apart - only an
%                 identifier cross-reference can.
%
% So VARIANT is a FLOOR on what ChEBI mapping would recover, never a ceiling.
% The actionable output is the ranked UNRESOLVED list: a handful of ligand
% names usually account for most of the lost reactions, and resolving those by
% hand is far cheaper than building the full pipeline.
%
% Writes Ligand_Match_Report.csv (every distinct ligand, matched or not, with
% its best model candidate) for manual review.

    if nargin == 1 && (ischar(model) || isstring(model)) && strcmpi(char(model),'selftest')
        R = selftest(); return;
    end
    if nargin < 3, opts = struct(); end
    if ~isfield(opts,'field'),     opts.field = 'AC'; end
    if ~isfield(opts,'minScore'),  opts.minScore = 0.50; end   % token-overlap threshold
    if ~isfield(opts,'topN'),      opts.topN = 15; end        % manual-review shortlist
    if ~isfield(opts,'csvFile'),   opts.csvFile = 'Ligand_Match_Report.csv'; end
    if ~isfield(opts,'verbose'),   opts.verbose = true; end

    nR = numel(model.rxns); nM = numel(model.mets);

    % ---- 0. what identifier annotation does the model already carry? -------
    idFields = {'metCHEBIID','metKEGGID','metInChIString','metInChIKey', ...
                'metSmiles','metPubChemID','metMetaNetXID','metHMDBID'};
    R.annot = struct();
    vprint(opts,'\n--- model identifier annotation -------------------------------\n');
    anyID = false;
    for q = 1:numel(idFields)
        f = idFields{q};
        if isfield(model, f)
            v = model.(f);
            n = sum(~cellfun(@(x) isempty(strtrim(char(x))), v(:)));
            R.annot.(f) = n;
            vprint(opts,'  %-16s present, %d/%d metabolites populated\n', f, n, nM);
            if n > 0, anyID = true; end
        else
            R.annot.(f) = -1;
            vprint(opts,'  %-16s absent\n', f);
        end
    end
    R.hasIdentifiers = anyID;
    if anyID
        vprint(opts,['  -> the model ALREADY carries identifiers. The model half of an\n' ...
                     '     ID-based cross-reference is free; only the BRENDA half needs work.\n']);
    else
        vprint(opts,['  -> no identifiers on the model. BiGG ids would have to be resolved\n' ...
                     '     externally (BiGG metabolites file or MetaNetX chem_xref).\n']);
    end

    % ---- 1. EC per reaction, same rule the pipeline uses -------------------
    EC = repmat({''}, nR, 1);
    if isfield(model,'rxnECNumbers'), EC = cellfun(@char, model.rxnECNumbers, 'uni', 0); end
    reg = parse_brenda_flat(flatFile, {opts.field});

    hasEC  = find(~cellfun(@isempty, EC));
    ecInBr = hasEC(cellfun(@(e) isKey(reg,e) && ~isempty(reg(e).(opts.field)), EC(hasEC)));

    % ---- 2. model metabolite name index ------------------------------------
    metNorm = cellfun(@normLigand, model.metNames, 'uni', 0);
    metSet  = containers.Map('KeyType','char','ValueType','double');
    for i = 1:nM
        if ~isempty(metNorm{i}) && ~isKey(metSet, metNorm{i}), metSet(metNorm{i}) = i; end
    end
    tokIx = tokenIndex(metNorm);          % token -> metabolite indices

    % ---- 3. walk the reactions exactly as iffl_coverage_check does ---------
    ligCount = containers.Map('KeyType','char','ValueType','double');
    ligHit   = containers.Map('KeyType','char','ValueType','double');
    nRxnAny = 0; nRxnMatched = 0;
    for jj = 1:numel(ecInBr)
        j = ecInBr(jj);
        ligs = reg(EC{j}).(opts.field);
        if isempty(ligs), continue; end
        nRxnAny = nRxnAny + 1;
        anyHit = false;
        for q = 1:numel(ligs)
            L = ligs{q};
            if isempty(L), continue; end
            if isKey(ligCount, L), ligCount(L) = ligCount(L) + 1; else, ligCount(L) = 1; end
            hit = isKey(metSet, L);
            if hit
                anyHit = true;
                if ~isKey(ligHit, L), ligHit(L) = 1; end
            end
        end
        if anyHit, nRxnMatched = nRxnMatched + 1; end
    end

    % ---- 4. for every unmatched ligand, is there a near counterpart? -------
    ligs = keys(ligCount);
    nLig = numel(ligs);
    rows = cell(nLig, 5);
    nMatched = 0; nVariant = 0; nUnresolved = 0;
    variantWeight = 0; unresolvedWeight = 0;
    for q = 1:nLig
        L = ligs{q}; cnt = ligCount(L);
        if isKey(ligHit, L)
            nMatched = nMatched + 1;
            rows(q,:) = {L, cnt, 'matched', model.metNames{metSet(L)}, 1.00};
            continue;
        end
        [bestIx, score] = bestCandidate(L, tokIx, metNorm);
        if ~isempty(bestIx) && score >= opts.minScore
            nVariant = nVariant + 1; variantWeight = variantWeight + cnt;
            rows(q,:) = {L, cnt, 'variant', model.metNames{bestIx}, score};
        else
            nUnresolved = nUnresolved + 1; unresolvedWeight = unresolvedWeight + cnt;
            if isempty(bestIx), nm = ''; score = 0; else, nm = model.metNames{bestIx}; end
            rows(q,:) = {L, cnt, 'unresolved', nm, score};
        end
    end

    % ---- 5. report ---------------------------------------------------------
    vprint(opts,'\n--- reaction level (field %s) ---------------------------------\n', opts.field);
    vprint(opts,'  reactions whose EC has activators in BRENDA .. %d\n', nRxnAny);
    vprint(opts,'  ...with at least one name matched ............ %d\n', nRxnMatched);
    vprint(opts,'  ...lost entirely to name matching ............ %d\n', nRxnAny - nRxnMatched);

    vprint(opts,'\n--- ligand level ---------------------------------------------\n');
    vprint(opts,'  distinct activator ligand names ............. %d\n', nLig);
    vprint(opts,'  matched by name ............................. %d\n', nMatched);
    vprint(opts,'  spelling variant of a model metabolite ...... %d   (%d mentions)\n', nVariant, variantWeight);
    vprint(opts,'  unresolved: absent OR a true synonym ........ %d   (%d mentions)\n', nUnresolved, unresolvedWeight);
    vprint(opts,['\n  The variant count is a FLOOR on what an identifier cross-reference\n' ...
                 '  would recover, not a ceiling: a synonym sharing no words with the\n' ...
                 '  model name (ATP vs adenosine triphosphate) lands in ''unresolved''\n' ...
                 '  and no string method can distinguish it from a genuine absence.\n']);

    % ---- Pareto shortlist: which names cost the most reactions? ------------
    unIx = find(strcmp(rows(:,3), 'unresolved'));
    if ~isempty(unIx)
        cnts = cell2mat(rows(unIx,2));
        [cnts, ord] = sort(cnts, 'descend'); unIx = unIx(ord);
        k = min(opts.topN, numel(unIx));
        vprint(opts,'\n--- highest-value unresolved ligands (resolve these by hand first) ---\n');
        vprint(opts,'  %-42s %8s  %s\n', 'ligand', 'rxns', 'nearest model name (score)');
        for q = 1:k
            r = unIx(q);
            vprint(opts,'  %-42s %8d  %s (%.2f)\n', trunc(rows{r,1},42), rows{r,2}, ...
                trunc(rows{r,4},30), rows{r,5});
        end
        top = sum(cnts(1:k));
        vprint(opts,'\n  These %d names account for %d of the %d unresolved mentions (%.0f%%).\n', ...
            k, top, unresolvedWeight, 100*top/max(unresolvedWeight,1));
    end

    writeCsv(opts.csvFile, rows, opts);

    R.nRxnWithAC = nRxnAny; R.nRxnMatched = nRxnMatched;
    R.nLigands = nLig; R.nMatched = nMatched;
    R.nVariant = nVariant; R.nUnresolved = nUnresolved;
    R.rows = rows;
    R.status = rows(:,3); R.ligand = rows(:,1);
end

function writeCsv(fname, rows, opts)
% plain fprintf, so this runs under both MATLAB and Octave
    fid = fopen(fname, 'w');
    if fid < 0, vprint(opts,'\n  Could not open %s for writing.\n', fname); return; end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'ligand,nReactions,status,bestModelMatch,score\n');
    for q = 1:size(rows,1)
        fprintf(fid, '"%s",%d,%s,"%s",%.3f\n', esc(rows{q,1}), rows{q,2}, ...
            rows{q,3}, esc(rows{q,4}), rows{q,5});
    end
    vprint(opts,'\n  Wrote %s. A high token score is a hint, not proof of identity -\n', fname);
    vprint(opts,'  check the ''variant'' rows before trusting them.\n');
end
function s = esc(s), s = strrep(char(s), '"', '""'); end
function s = trunc(s, n)
    s = char(s);
    if numel(s) > n, s = [s(1:n-1) '~']; end
end

% ------------------------------------------------------------------ helpers
function ix = tokenIndex(metNorm)
    ix = containers.Map('KeyType','char','ValueType','any');
    for i = 1:numel(metNorm)
        for t = strsplit(metNorm{i}, ' ')
            k = t{1};
            if isempty(k) || numel(k) < 3, continue; end     % skip 1-2 char noise
            if isKey(ix,k), ix(k) = [ix(k) i]; else, ix(k) = i; end
        end
    end
end

function [bestIx, bestScore] = bestCandidate(L, tokIx, metNorm)
% Jaccard over token sets, restricted to metabolites sharing a token.
    bestIx = []; bestScore = 0;
    lt = strsplit(L, ' ');
    lt = lt(cellfun(@(t) numel(t) >= 3, lt));
    if isempty(lt), return; end
    cand = [];
    for q = 1:numel(lt)
        if isKey(tokIx, lt{q}), cand = [cand tokIx(lt{q})]; end %#ok<AGROW>
    end
    if isempty(cand), return; end
    cand = unique(cand);
    for i = cand(:)'
        mt = strsplit(metNorm{i}, ' ');
        mt = mt(cellfun(@(t) numel(t) >= 3, mt));
        if isempty(mt), continue; end
        inter = numel(intersect(lt, mt));
        uni   = numel(union(lt, mt));
        s = inter / uni;
        if s > bestScore, bestScore = s; bestIx = i; end
    end
end

function s = normLigand(s)
    s = lower(strtrim(char(s)));
    s = regexprep(s, '^((l|d|dl|\(s\)|\(r\)|alpha|beta)-)+', '');
    s = regexprep(s, '[^a-z0-9]+', ' ');
    s = strtrim(s);
end

function vprint(opts, varargin)
    if opts.verbose, fprintf(varargin{:}); end
end

% ================================ selftest ==================================
function R = selftest()
    model.metNames = {'L-Cysteine';'D-Glucose 6-Phosphate';'Pyruvate'; ...
                      'Phosphatidylinositol 4,5-Bisphosphate';'ATP'};
    model.mets = {'cys_L[c]';'g6p[c]';'pyr[c]';'pail45p_hs[c]';'atp[c]'};
    model.rxns = {'R1';'R2';'R3'};
    model.rxnECNumbers = {'1.1.1.1';'2.2.2.2';'3.3.3.3'};
    model.metCHEBIID = {'CHEBI:17561';'CHEBI:4170';'';'';'CHEBI:15422'};

    tmp = [tempname '.txt']; fid = fopen(tmp,'w');
    % R1: 'cysteine' matches L-Cysteine exactly after normalisation.
    % R2: 'glucose 6-phosphate dipotassium salt' shares tokens with
    %     'D-Glucose 6-Phosphate' -> VARIANT (string-detectable).
    % R2: 'adenosine triphosphate' shares NO token with the model's 'ATP',
    %     so it lands in UNRESOLVED even though the compound is present.
    %     That is the documented blind spot of the string heuristic.
    % R3: 'colchicine' is genuinely nowhere in the model -> UNRESOLVED.
    fprintf(fid, ['ID\t1.1.1.1\nAC\t#1# cysteine <1>\n///\n' ...
                  'ID\t2.2.2.2\nAC\t#1# glucose 6-phosphate dipotassium salt <1>\n' ...
                  'AC\t#1# adenosine triphosphate <1>\n///\n' ...
                  'ID\t3.3.3.3\nAC\t#1# colchicine <1>\n///\n']);
    fclose(fid);

    R = diagnose_ligand_matching(model, tmp, struct('csvFile', [tempname '.csv']));
    delete(tmp);

    assert(R.nRxnWithAC == 3, 'reaction count wrong');
    assert(R.nRxnMatched == 1, 'matched-reaction count wrong');
    g = @(name) char(R.status{strcmp(R.ligand, name)});
    assert(strcmp(g('cysteine'),'matched'), 'cysteine should match by name');
    assert(strcmp(g('glucose 6 phosphate dipotassium salt'),'variant'), ...
        'the dipotassium-salt form should be flagged as a spelling variant');
    assert(strcmp(g('adenosine triphosphate'),'unresolved'), ...
        'a no-shared-token synonym must land in unresolved - that is the blind spot');
    assert(strcmp(g('colchicine'),'unresolved'), 'colchicine is absent from the model');
    fprintf(['\nselftest PASSED: matched / variant / unresolved all classified,\n' ...
             'including the synonym blind spot (adenosine triphosphate vs ATP).\n\n']);
end
