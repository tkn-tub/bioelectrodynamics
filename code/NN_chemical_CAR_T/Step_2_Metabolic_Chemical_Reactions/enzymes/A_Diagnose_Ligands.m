function R = A_Diagnose_Ligands(arg1, arg2)
% A_DIAGNOSE_LIGANDS  Self-contained: how much of the BRENDA activator
% annotation is lost to name matching, and how much of it is recoverable?
%
%   A_Diagnose_Ligands                       % find model + BRENDA automatically
%   A_Diagnose_Ligands('selftest')           % synthetic check, no data needed
%   A_Diagnose_Ligands(modelFile, brendaFile)% explicit paths
%
% Nothing else is required on the path: the BRENDA parser and every helper are
% local functions in this file. COBRA is used only to read an .xml model.
%
% WRITES (so a closing session loses nothing):
%   Ligand_Diagnostic_log.txt   full console transcript
%   Ligand_Match_Report.csv     every distinct ligand, its status, best match
%   ligand_diagnostic.mat       the result struct R
%
% ---------------------------------------------------------------------------
% WHAT IT MEASURES
%
% iffl_coverage_check keeps a reaction only when at least one BRENDA activator
% NAME matches a model metabolite NAME after normalisation. Everything else is
% discarded silently. This script splits that loss into:
%
%   matched     the name matched outright.
%   variant     a spelling / formatting variant of a compound that IS in the
%               model ('glucose 6-phosphate dipotassium salt' against
%               'D-Glucose 6-Phosphate'). Certainly recoverable.
%   unresolved  no string-detectable counterpart. EITHER genuinely absent from
%               the model (a plant metabolite, a drug) OR a true synonym with
%               no shared words, such as 'adenosine triphosphate' against a
%               model that writes 'ATP'. No string method separates those two;
%               only an identifier cross-reference can.
%
% So `variant` is a FLOOR on what a ChEBI/KEGG cross-reference would recover,
% never a ceiling. The actionable output is the ranked `unresolved` list: a
% handful of names usually account for most of the lost reactions, and
% resolving those by hand is far cheaper than building the full pipeline.
% ---------------------------------------------------------------------------

    if nargin >= 1 && (ischar(arg1) || isstring(arg1)) && strcmpi(char(arg1), 'selftest')
        R = selftest(); return;
    end

    diary off;
    logFile = 'Ligand_Diagnostic_log.txt';
    if exist(logFile, 'file'), delete(logFile); end
    diary(logFile); diary on;
    closeDiary = onCleanup(@() diary('off'));     % closes the log even on error
    fprintf('Ligand matching diagnostic, %s\n', datestr(now));

    % ---- locate the inputs -------------------------------------------------
    modelFile = ''; brendaFile = '';
    if nargin >= 1 && ~isempty(arg1), modelFile  = char(arg1); end
    if nargin >= 2 && ~isempty(arg2), brendaFile = char(arg2); end

    if isempty(modelFile)
        addPathByName('code');
        modelFile = which('Puniya2021_TNM1055.xml');
        if isempty(modelFile), modelFile = which('Puniya2021_TNM1055.mat'); end
    end
    assert(~isempty(modelFile) && exist(modelFile,'file') == 2, ...
        'Model not found. Pass it explicitly: A_Diagnose_Ligands(''path/to/model.xml'', ''path/to/brenda.txt'')');

    if isempty(brendaFile), brendaFile = which('brenda_2026_1.txt'); end
    assert(~isempty(brendaFile) && exist(brendaFile,'file') == 2, ...
        'brenda_2026_1.txt not found. Pass it explicitly as the second argument.');

    % ---- load the model ----------------------------------------------------
    [~,~,ext] = fileparts(modelFile);
    if strcmpi(ext, '.mat')
        L = load(modelFile); fn = fieldnames(L); model = L.(fn{1});
    else
        if exist('readCbModel','file') ~= 2, initCobraToolbox(false); end
        model = readCbModel(modelFile);
    end
    fprintf('Model  : %s\n', modelFile);
    fprintf('BRENDA : %s\n', brendaFile);
    fprintf('         %d reactions, %d metabolites\n', numel(model.rxns), numel(model.mets));

    % ---- run ---------------------------------------------------------------
    R = runDiagnostic(model, brendaFile, struct());
    R.modelFile = modelFile; R.brendaFile = brendaFile; R.stamp = datestr(now);

    try
        save('ligand_diagnostic.mat', 'R', '-v7.3');
        fprintf('\nSaved ligand_diagnostic.mat\n');
    catch ME
        fprintf('\nCould not save ligand_diagnostic.mat: %s\n', ME.message);
    end
    fprintf('Done %s. Files: %s, %s, ligand_diagnostic.mat\n', ...
        datestr(now), logFile, R.csvFile);
    diary off;
end

% ===========================================================================
% THE DIAGNOSTIC
% ===========================================================================
function R = runDiagnostic(model, flatFile, opts)
    if ~isfield(opts,'field'),    opts.field = 'AC';  end   % 'IN' for inhibitors
    if ~isfield(opts,'minScore'), opts.minScore = 0.50; end
    if ~isfield(opts,'topN'),     opts.topN = 20; end
    if ~isfield(opts,'csvFile'),  opts.csvFile = 'Ligand_Match_Report.csv'; end
    if ~isfield(opts,'verbose'),  opts.verbose = true; end

    nR = numel(model.rxns); nM = numel(model.mets);

    % ---- 0. identifier annotation already on the model ---------------------
    idFields = {'metCHEBIID','metKEGGID','metInChIString','metInChIKey', ...
                'metSmiles','metPubChemID','metMetaNetXID','metHMDBID'};
    R.annot = struct(); anyID = false;
    vprint(opts,'\n--- model identifier annotation -------------------------------\n');
    for q = 1:numel(idFields)
        f = idFields{q};
        if isfield(model, f)
            v = model.(f);
            n = 0;
            for i = 1:numel(v)
                if ~isempty(strtrim(char(v{i}))), n = n + 1; end
            end
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
        vprint(opts,['  -> the model ALREADY carries identifiers, so the model half of an\n' ...
                     '     ID cross-reference is free. Only the BRENDA half needs work.\n']);
    else
        vprint(opts,['  -> no identifiers on the model. BiGG ids would have to be resolved\n' ...
                     '     externally (BiGG metabolites file, or MetaNetX chem_xref).\n']);
    end

    % ---- 1. EC per reaction, exactly as the pipeline does it ----------------
    EC = repmat({''}, nR, 1);
    if isfield(model,'rxnECNumbers')
        for j = 1:nR, EC{j} = char(model.rxnECNumbers{j}); end
    end
    reg = parseBrendaFlat(flatFile, {opts.field});

    hasEC = []; ecInBr = [];
    for j = 1:nR
        if isempty(strtrim(EC{j})), continue; end
        hasEC(end+1) = j; %#ok<AGROW>
        e = strtrim(EC{j});
        if isKey(reg, e) && ~isempty(reg(e).(opts.field)), ecInBr(end+1) = j; end %#ok<AGROW>
    end
    vprint(opts,'\n--- reaction level (field %s) ---------------------------------\n', opts.field);
    vprint(opts,'  reactions total ............................. %d\n', nR);
    vprint(opts,'  reactions with an EC number ................. %d\n', numel(hasEC));
    vprint(opts,'  ...whose EC has %s entries in BRENDA ......... %d\n', opts.field, numel(ecInBr));

    % ---- 2. model metabolite name index ------------------------------------
    metNorm = cell(nM,1);
    for i = 1:nM, metNorm{i} = normLigand(model.metNames{i}); end
    metSet = containers.Map('KeyType','char','ValueType','double');
    for i = 1:nM
        k = metNorm{i};
        if ~isempty(k) && ~isKey(metSet, k), metSet(k) = i; end
    end
    tokIx = tokenIndex(metNorm);

    % ---- 3. walk the reactions ---------------------------------------------
    ligCount = containers.Map('KeyType','char','ValueType','double');
    ligHit   = containers.Map('KeyType','char','ValueType','double');
    nRxnAny = 0; nRxnMatched = 0;
    for jj = 1:numel(ecInBr)
        j = ecInBr(jj);
        ligs = reg(strtrim(EC{j})).(opts.field);
        if isempty(ligs), continue; end
        nRxnAny = nRxnAny + 1; anyHit = false;
        for q = 1:numel(ligs)
            L = ligs{q};
            if isempty(L), continue; end
            if isKey(ligCount, L), ligCount(L) = ligCount(L) + 1; else, ligCount(L) = 1; end
            if isKey(metSet, L)
                anyHit = true;
                if ~isKey(ligHit, L), ligHit(L) = 1; end
            end
        end
        if anyHit, nRxnMatched = nRxnMatched + 1; end
    end
    vprint(opts,'  ...with at least one ligand name matched .... %d\n', nRxnMatched);
    vprint(opts,'  ...lost entirely to name matching ........... %d\n', nRxnAny - nRxnMatched);

    % ---- 4. classify every distinct ligand ---------------------------------
    ligs = keys(ligCount); nLig = numel(ligs);
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

    vprint(opts,'\n--- ligand level ---------------------------------------------\n');
    vprint(opts,'  distinct activator ligand names ............. %d\n', nLig);
    vprint(opts,'  matched by name ............................. %d\n', nMatched);
    vprint(opts,'  spelling variant of a model metabolite ...... %d   (%d mentions)\n', nVariant, variantWeight);
    vprint(opts,'  unresolved: absent OR a true synonym ........ %d   (%d mentions)\n', nUnresolved, unresolvedWeight);
    vprint(opts,['\n  The variant count is a FLOOR on what an identifier cross-reference\n' ...
                 '  would recover, not a ceiling: a synonym sharing no words with the\n' ...
                 '  model name (adenosine triphosphate vs ATP) lands in ''unresolved'',\n' ...
                 '  and no string method can tell it from a genuine absence.\n']);

    % ---- 5. Pareto shortlist ------------------------------------------------
    unIx = [];
    for q = 1:nLig
        if strcmp(rows{q,3}, 'unresolved'), unIx(end+1) = q; end %#ok<AGROW>
    end
    if ~isempty(unIx)
        cnts = zeros(1, numel(unIx));
        for q = 1:numel(unIx), cnts(q) = rows{unIx(q), 2}; end
        [cnts, ord] = sort(cnts, 'descend'); unIx = unIx(ord);
        k = min(opts.topN, numel(unIx));
        vprint(opts,'\n--- highest-value unresolved ligands (resolve these first) ----\n');
        vprint(opts,'  %-40s %6s  %s\n', 'ligand', 'rxns', 'nearest model name (score)');
        for q = 1:k
            r = unIx(q);
            vprint(opts,'  %-40s %6d  %s (%.2f)\n', trunc(rows{r,1},40), rows{r,2}, ...
                trunc(rows{r,4},32), rows{r,5});
        end
        top = sum(cnts(1:k));
        vprint(opts,'\n  These %d names account for %d of the %d unresolved mentions (%.0f%%).\n', ...
            k, top, unresolvedWeight, 100*top/max(unresolvedWeight,1));
        if 100*top/max(unresolvedWeight,1) > 50
            vprint(opts,['  Over half the loss sits in that short list - resolve those by hand\n' ...
                         '  and you recover most of the annotation with no infrastructure.\n']);
        end
    end

    writeCsv(opts.csvFile, rows, opts);

    R.nRxnTotal = nR; R.nRxnWithEC = numel(hasEC);
    R.nRxnWithAC = nRxnAny; R.nRxnMatched = nRxnMatched;
    R.nLigands = nLig; R.nMatched = nMatched;
    R.nVariant = nVariant; R.nUnresolved = nUnresolved;
    R.variantMentions = variantWeight; R.unresolvedMentions = unresolvedWeight;
    R.rows = rows; R.ligand = rows(:,1); R.status = rows(:,3);
    R.csvFile = opts.csvFile; R.field = opts.field;
end

% ===========================================================================
% BRENDA FLAT FILE PARSER  (inlined so this script stands alone)
% ===========================================================================
function reg = parseBrendaFlat(flatFile, fields)
    if nargin < 2, fields = {'AC','IN'}; end
    skipLig = {'more','additional information',''};
    reg = containers.Map('KeyType','char','ValueType','any');
    fid = fopen(flatFile, 'r', 'n', 'UTF-8');
    assert(fid >= 0, 'Cannot open BRENDA flat file: %s', flatFile);
    clean = onCleanup(@() fclose(fid)); %#ok<NASGU>

    ec = ''; buf = ''; TAB = char(9);
    while true
        raw = fgetl(fid);
        atEnd = ~ischar(raw);
        isCont = ~atEnd && ~isempty(raw) && (raw(1) == ' ' || raw(1) == TAB);
        if isCont
            buf = [buf ' ' strtrim(raw)]; %#ok<AGROW>
            continue;
        end
        [ec, reg] = flushLine(buf, ec, reg, fields, skipLig, TAB);
        if atEnd, break; end
        if strcmp(strtrim(raw), '///'), buf = ''; continue; end
        buf = raw;
    end
    ks = keys(reg);
    for i = 1:numel(ks)
        s = reg(ks{i});
        for f = 1:numel(fields), s.(fields{f}) = unique(s.(fields{f})); end
        reg(ks{i}) = s;
    end
end

function [ec, reg] = flushLine(line, ec, reg, fields, skipLig, TAB)
    if isempty(line), return; end
    t = strfind(line, TAB);
    if isempty(t), return; end
    tag  = line(1:t(1)-1);
    body = line(t(1)+1:end);
    if strcmp(tag, 'ID')
        ec = strtrim(regexprep(body, '\s.*', ''));
        if ~isKey(reg, ec)
            s = struct();
            for f = 1:numel(fields), s.(fields{f}) = {}; end
            reg(ec) = s;
        end
        return;
    end
    if any(strcmp(tag, fields)) && ~isempty(ec)
        body = regexprep(body, '^#[^#]*#\s*', '');
        lig  = regexp(body, '^[^(\{<]*', 'match', 'once');
        lig  = normLigand(lig);
        if ~any(strcmp(lig, skipLig)) && isRealLigand(lig)
            s = reg(ec); s.(tag){end+1} = lig; reg(ec) = s;
        end
    end
end

% ===========================================================================
% HELPERS
% ===========================================================================
function ix = tokenIndex(metNorm)
    ix = containers.Map('KeyType','char','ValueType','any');
    for i = 1:numel(metNorm)
        parts = strsplit(metNorm{i}, ' ');
        for q = 1:numel(parts)
            k = parts{q};
            if numel(k) < 3, continue; end
            if isKey(ix, k), ix(k) = [ix(k) i]; else, ix(k) = i; end
        end
    end
end

function [bestIx, bestScore] = bestCandidate(L, tokIx, metNorm)
% Jaccard over token sets, restricted to metabolites sharing at least one token
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
        s = numel(intersect(lt, mt)) / numel(union(lt, mt));
        if s > bestScore, bestScore = s; bestIx = i; end
    end
end

function writeCsv(fname, rows, opts)
    fid = fopen(fname, 'w');
    if fid < 0, vprint(opts, '\n  Could not open %s for writing.\n', fname); return; end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'ligand,nReactions,status,bestModelMatch,score\n');
    for q = 1:size(rows,1)
        fprintf(fid, '"%s",%d,%s,"%s",%.3f\n', esc(rows{q,1}), rows{q,2}, ...
            rows{q,3}, esc(rows{q,4}), rows{q,5});
    end
    vprint(opts, '\n  Wrote %s. A high token score is a hint, not proof of identity -\n', fname);
    vprint(opts, '  check the ''variant'' rows before trusting them.\n');
end

function ok = isRealLigand(s)
% Reject parse fragments. A BRENDA AC/IN line that wraps, or whose ligand field
% is malformed, can yield a bare token like 'n', '3' or '2'. Those appeared in
% the IN run with 400+ mentions each. A one-character or purely numeric name is
% never a real compound, and if it happened to match a model metabolite it
% would create a phantom gate.
    s = char(s);
    ok = numel(s) >= 3 && ~isempty(regexp(s, '[a-z]{2}', 'once'));
end
function s = esc(s), s = strrep(char(s), '"', '""'); end

function s = trunc(s, n)
    s = char(s);
    if numel(s) > n, s = [s(1:n-1) '~']; end
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

function addPathByName(targetName)
% walk up from pwd until a folder with this name is found, then add it
    currFolder = pwd; found = false;
    while true
        [parentFolder, currentName] = fileparts(currFolder);
        if strcmpi(currentName, targetName), found = true; break; end
        if isempty(parentFolder) || strcmp(currFolder, parentFolder), break; end
        currFolder = parentFolder;
    end
    if found
        addpath(genpath(currFolder));
        fprintf('Added to path: %s\n', currFolder);
    else
        fprintf('Note: no parent folder named "%s" found; relying on the current path.\n', targetName);
    end
end

% ===========================================================================
% SELFTEST
% ===========================================================================
function R = selftest()
    model.metNames = {'L-Cysteine';'D-Glucose 6-Phosphate';'Pyruvate'; ...
                      'Phosphatidylinositol 4,5-Bisphosphate';'ATP'};
    model.mets = {'cys_L[c]';'g6p[c]';'pyr[c]';'pail45p_hs[c]';'atp[c]'};
    model.rxns = {'R1';'R2';'R3'};
    model.rxnECNumbers = {'1.1.1.1';'2.2.2.2';'3.3.3.3'};
    model.metCHEBIID = {'CHEBI:17561';'CHEBI:4170';'';'';'CHEBI:15422'};

    tmp = [tempname '.txt']; fid = fopen(tmp, 'w');
    % R1: 'cysteine' matches L-Cysteine exactly after normalisation.
    % R2: 'glucose 6-phosphate dipotassium salt' shares tokens with
    %     'D-Glucose 6-Phosphate' -> VARIANT (string-detectable).
    % R2: 'adenosine triphosphate' shares NO token with the model's 'ATP', so
    %     it lands in UNRESOLVED even though the compound is present. That is
    %     the documented blind spot of the string heuristic.
    % R3: 'colchicine' is genuinely nowhere in the model -> UNRESOLVED.
    fprintf(fid, ['ID\t1.1.1.1\nAC\t#1# cysteine <1>\n///\n' ...
                  'ID\t2.2.2.2\nAC\t#1# glucose 6-phosphate dipotassium salt <1>\n' ...
                  'AC\t#1# adenosine triphosphate <1>\n///\n' ...
                  'ID\t3.3.3.3\nAC\t#1# colchicine <1>\n///\n']);
    fclose(fid);

    R = runDiagnostic(model, tmp, struct('csvFile', [tempname '.csv']));
    delete(tmp);

    assert(R.nRxnWithAC == 3, 'reaction count wrong');
    assert(R.nRxnMatched == 1, 'matched-reaction count wrong');
    g = @(name) char(R.status{strcmp(R.ligand, name)});
    assert(strcmp(g('cysteine'), 'matched'), 'cysteine should match by name');
    assert(strcmp(g('glucose 6 phosphate dipotassium salt'), 'variant'), ...
        'the dipotassium-salt form should be flagged as a spelling variant');
    assert(strcmp(g('adenosine triphosphate'), 'unresolved'), ...
        'a no-shared-token synonym must land in unresolved - that is the blind spot');
    assert(strcmp(g('colchicine'), 'unresolved'), 'colchicine is absent from the model');
    % --- degenerate inputs must not crash on the real machine -------------
    m2.metNames = {'Pyruvate';'ATP'}; m2.mets = {'pyr[c]';'atp[c]'};
    m2.rxns = {'R1'}; m2.rxnECNumbers = {''};          % no EC anywhere
    tmp2 = [tempname '.txt']; f2 = fopen(tmp2,'w');
    fprintf(f2, 'ID\t1.1.1.1\nAC\t#1# cysteine <1>\n///\n'); fclose(f2);
    R2 = runDiagnostic(m2, tmp2, struct('verbose', false, 'csvFile', [tempname '.csv']));
    delete(tmp2);
    assert(R2.nRxnWithAC == 0 && R2.nLigands == 0, 'empty-EC path mishandled');

    % --- a model with no identifier fields at all -------------------------
    m3 = rmfield(model, 'metCHEBIID');
    R3 = runDiagnostic(m3, tmp3Make(), struct('verbose', false, 'csvFile', [tempname '.csv']));
    assert(R3.hasIdentifiers == false, 'should report no identifiers');

    fprintf(['\nselftest PASSED: matched / variant / unresolved all classified,\n' ...
             'including the synonym blind spot (adenosine triphosphate vs ATP),\n' ...
             'plus the no-EC and no-identifier degenerate paths.\n\n']);
end

function p = tmp3Make()
    p = [tempname '.txt']; f = fopen(p,'w');
    fprintf(f, 'ID\t1.1.1.1\nAC\t#1# cysteine <1>\n///\n'); fclose(f);
end
