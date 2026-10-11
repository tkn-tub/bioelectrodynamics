function R = brendaKineticsTable(txtFile, xlsxFile, email, password, organism)
%BRENDAKINETICSTABLE  Motif report (.txt)  ->  BRENDA  ->  kinetics workbook (.xlsx)
%
%   brendaKineticsTable('Block_7_Selected_Motifs.txt', 'kinetics_BRENDA.xlsx', email, password)
%   brendaKineticsTable(..., 'Homo sapiens')     % organism to prefer (this is the default)
%   R = brendaKineticsTable(...)                 % also returns the rows as a struct array
%
%   Needs brendaQuery.m on the MATLAB path.
%
%   The workbook has three sheets:
%     values           one row per reaction, laid out like the kinetics.xlsx template
%     KM details       every substrate of every reaction, with n / min / median / max
%     SA MW details    the same spread for specific activity and molecular weight
%
%   Which BRENDA value is used (change the settings below to do it differently):
%     K_M  median of the wild-type values for the reaction's own substrate
%     SA   maximum of the wild-type values (the purest preparation reported)
%     MW   median of the subunit masses (SDS-PAGE / sequence), so that
%          k_cat = SA * MW / 60 is per subunit
%   Rate constant, gain and time constant (formulas in the workbook; [E] is
%   written as 0, so G and tau stay empty until real concentrations are typed in), following Eqs. (23)-(27) of the paper:
%     k   = k_cat [E] / K_M                  linear regime, substrate << K_M
%     tau = 1 / sum(k of the removal reactions) = K_M,r / (k_cat,r [E]r)
%     G   = k / sum(k of the removal reactions) = k_f / k_r
%   The removal reactions of a row are the reactions that consume its product
%   (column "Removed by").
%
%   Each one prefers the chosen organism and falls back to all organisms; K_M
%   falls back further to the enzyme's other substrates (rough placeholder,
%   written as OTHER SUBSTRATES). The Notes column says which case applied.

    if nargin < 5 || isempty(organism), organism = 'Homo sapiens'; end
    if exist('brendaQuery', 'file') ~= 2
        error('brendaKineticsTable:noBrendaQuery', 'brendaQuery.m must be on the MATLAB path.');
    end

    % ------------------------------ settings ------------------------------
    cfg.organism  = strtrim(char(organism));
    cfg.kmStat    = 'median';              % 'median' | 'geomean' | 'min' | 'max'
    cfg.kcatStat  = 'max';
    cfg.saStat    = 'max';
    cfg.mwStat    = 'median';
    cfg.kmOtherSubstrates = true;          % false: leave K_M empty when BRENDA has no value for the substrate
    cfg.cacheFile = 'brenda_cache.mat';    % raw BRENDA replies; delete the file to query again
    cfg.useCache  = true;

    % EC numbers in the model annotation (Puniya2021_TNM1055.xml) that are wrong
    % for the reaction.
    %   reaction ID, correct EC ('' = none exists), why the model's EC is wrong, use the correct EC?
    % true : BRENDA is queried with the correct EC and the row is annotated as corrected
    % false: the model's EC is queried and the row is only annotated as wrong
    cfg.ecFix = { ...
        'RNDR3',   '1.17.4.1', 'ribonucleoside-diphosphate reductase is 1.17.4.1',  true; ...
        'PUNP1',   '2.4.2.1',  'purine-nucleoside phosphorylase is 2.4.2.1',        true; ...
        'FTHFL',   '6.3.4.3',  'formate-tetrahydrofolate ligase is 6.3.4.3 (6.3.2.17 is folylpolyglutamate synthase)', false; ...
        'ASPGLUm', '',         'this is a transporter with no EC number (4.1.1.15 is glutamate decarboxylase)', false};

    % Group names in the report -> titles in column A
    cfg.groupTitle = { ...
        'adapter',          'cAMP Input Adapter'; ...
        'differentiator 1', 'Chemical Differentiator 1'; ...
        'differentiator 2', 'Chemical Differentiator 2'; ...
        'counter',          'Counter'; ...
        'NN',               'Neural Network'};

    % Never used as the K_M substrate / used only when nothing else is there.
    % Names exactly as they appear in the report.
    cfg.neverSubstrate = {'Proton', 'Water'};
    cfg.currency = {'Adenosine Triphosphate', 'Adenosine Diphosphate', 'Orthophosphate', 'Diphosphate', ...
        'Nicotinamide Adenine Dinucleotide', 'Nicotinamide Adenine Dinucleotide - Reduced', ...
        'Nicotinamide Adenine Dinucleotide Phosphate', 'Nicotinamide Adenine Dinucleotide Phosphate - Reduced'};

    % Name in the report -> names BRENDA uses for the same compound.
    % Compared ignoring case and punctuation. Add a line when the Notes column
    % says "other substrates" for a compound that BRENDA does list.
    cfg.synonyms = { ...
        '3,5-Cyclic Adenosine Monophosphate', {'cAMP', '3'',5''-cAMP', '3'',5''-cyclic AMP', 'cyclic AMP', 'adenosine 3'',5''-cyclic monophosphate'}; ...
        'Adenosine Monophosphate',            {'AMP', '5''-AMP', 'adenosine 5''-monophosphate'}; ...
        'Adenosine',                          {'adenosine'}; ...
        'Orthophosphate',                     {'phosphate', 'orthophosphate', 'Pi'}; ...
        'Alpha-D-Ribose 1-Phosphate',         {'alpha-D-ribose 1-phosphate', 'D-ribose 1-phosphate', 'ribose 1-phosphate'}; ...
        'Alpha-D-Ribose 5-Phosphate',         {'D-ribose 5-phosphate', 'ribose 5-phosphate', 'alpha-D-ribose 5-phosphate'}; ...
        'D-Xylulose 5-Phosphate',             {'D-xylulose 5-phosphate', 'xylulose 5-phosphate'}; ...
        'Glyceraldehyde 3-Phosphate',         {'D-glyceraldehyde 3-phosphate', 'glyceraldehyde 3-phosphate', 'DL-glyceraldehyde 3-phosphate'}; ...
        'Sedoheptulose 7-Phosphate',          {'sedoheptulose 7-phosphate', 'D-sedoheptulose 7-phosphate'}; ...
        'Dihydroxyacetone Phosphate',         {'dihydroxyacetone phosphate', 'glycerone phosphate', 'DHAP'}; ...
        'D-Erythrose 4-Phosphate',            {'D-erythrose 4-phosphate', 'erythrose 4-phosphate'}; ...
        'Sedoheptulose 1,7-Bisphosphate',     {'sedoheptulose 1,7-bisphosphate', 'D-sedoheptulose 1,7-bisphosphate'}; ...
        'Uridine Diphosphate',                {'UDP'}; ...
        'D-Fructose 6-Phosphate',             {'D-fructose 6-phosphate', 'fructose 6-phosphate', 'beta-D-fructose 6-phosphate'}; ...
        'Uridine-5''-TrIphosphate',           {'UTP'}; ...
        'Adenosine Triphosphate',             {'ATP', 'MgATP2-', 'MgATP', 'Mg-ATP'}; ...
        'Cytidine 5''-Monophosphate',         {'CMP'}; ...
        'Cytidine-5''-DIphosphate',           {'CDP'}; ...
        'Reduced Thioredoxin',                {'thioredoxin', 'reduced thioredoxin'}; ...
        'Deoxycytidylic Acid',                {'dCMP'}; ...
        'Deoxycytidine-5''-Diphosphate',      {'dCDP'}; ...
        'L-Aspartate',                        {'L-aspartate', 'aspartate', 'L-Asp', 'L-aspartic acid'}; ...
        'L-Glutamate',                        {'L-glutamate', 'glutamate', 'L-Glu'}; ...
        'L-Glutamine',                        {'L-glutamine', 'glutamine', 'L-Gln'}; ...
        '5-Phospho-Alpha-D-Ribose 1-Diphosphate', {'5-phospho-alpha-D-ribose 1-diphosphate', 'PRPP', '5-phosphoribosyl 1-diphosphate', '5-phospho-D-ribose 1-diphosphate', 'phosphoribosyl diphosphate', '5-phosphoribosyl-1-pyrophosphate'}; ...
        'Formate',                            {'formate'}; ...
        '5,6,7,8-Tetrahydrofolate',           {'tetrahydrofolate', '5,6,7,8-tetrahydrofolate', '(6S)-tetrahydrofolate', 'THF'}; ...
        '10-Formyltetrahydrofolate',          {'10-formyltetrahydrofolate', 'N10-formyltetrahydrofolate', '(6R)-10-formyltetrahydrofolate'}; ...
        '5-Amino-1-(5-Phospho-D-Ribosyl)Imidazole-4-Carboxamide', {'AICAR', '5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxamide', '5-aminoimidazole-4-carboxamide ribonucleotide', '5-amino-4-imidazolecarboxamide ribotide'}; ...
        '(S)-Lactaldehyde',                   {'(S)-lactaldehyde', 'L-lactaldehyde', 'lactaldehyde'}; ...
        'Nicotinamide Adenine Dinucleotide Phosphate - Reduced', {'NADPH'}; ...
        'Nicotinamide Adenine Dinucleotide',  {'NAD+', 'NAD'}; ...
        '(R)-Propane-1,2-Diol',               {'(R)-propane-1,2-diol', '(R)-1,2-propanediol', '1,2-propanediol', 'propane-1,2-diol'}; ...
        '(S)-Propane-1,2-Diol',               {'(S)-propane-1,2-diol', '(S)-1,2-propanediol', '1,2-propanediol', 'propane-1,2-diol'}; ...
        'Hydroxyacetone',                     {'hydroxyacetone', 'acetol'}; ...
        '3-Phospho-D-Glyceroyl Phosphate',    {'3-phospho-D-glyceroyl phosphate', '1,3-bisphosphoglycerate', '1,3-diphosphoglycerate', '1,3-bisphospho-D-glycerate'}};
    % ----------------------------------------------------------------------

    % 1) reactions from the report
    rx = parseReport(txtFile);
    if isempty(rx)
        error('brendaKineticsTable:noReactions', 'No reactions found in %s.', txtFile);
    end
    rx = applyEcFixes(rx, cfg);
    fprintf('%d reactions read from %s\n', numel(rx), txtFile);

    % 2) BRENDA, once per EC number (all organisms; the filtering is done here)
    data = fetchAll(unique({rx(~cellfun('isempty', {rx.ec})).ec}), email, password, cfg);

    % 3) choose one value of each kind per reaction
    [rx.km] = deal([]);  [rx.kcat] = deal([]);  [rx.sa] = deal([]);  [rx.mw] = deal([]);  [rx.primary] = deal(0);
    for i = 1:numel(rx)
        rx(i) = chooseValues(rx(i), data, cfg);
    end

    rx = findRemovers(rx, cfg);

    % 4) workbook
    writeWorkbook(xlsxFile, rx, cfg);
    fprintf('Wrote %s\n\n', xlsxFile);
    printSummary(rx);
    if nargout > 0, R = rx; end
end

% =========================================================================
%  Report parsing
% =========================================================================
function rx = parseReport(txtFile)
    txt   = fileread(txtFile);
    lines = regexp(txt, '\r?\n', 'split');
    rx = struct('motif', {}, 'label', {}, 'group', {}, 'id', {}, 'name', {}, ...
                'ecFile', {}, 'names', {});
    k = 0;  inRx = false;  motif = 1;
    for i = 1:numel(lines)
        L = lines{i};
        mo = regexp(L, 'SELECTED (?:COMPLETE )?MOTIF\s+(\d+)', 'tokens', 'once');
        if ~isempty(mo), motif = str2double(mo{1}); inRx = false; continue; end
        if ~isempty(regexp(L, '^\s*(---|###)', 'once')), inRx = false; continue; end

        a = regexp(L, '^\s*Adapter reaction\s+(\d+)\s*\([^)]*\)\s*:\s*(\S+)\s*\|\s*(.*?)\s*$', 'tokens', 'once');
        m = regexp(L, '^\s*(r_\w+)\s*\(([^:()]+):.*\)\s*$', 'tokens', 'once');
        if ~isempty(a)
            k = k + 1;  inRx = true;
            rx(k).motif = motif;  rx(k).label = ['Adapter ' a{1}];  rx(k).group = 'adapter';
            rx(k).id = a{2};  rx(k).name = a{3};  rx(k).ecFile = '';  rx(k).names = '';
        elseif ~isempty(m)
            k = k + 1;  inRx = true;
            rx(k).motif = motif;  rx(k).label = m{1};  rx(k).group = strtrim(m{2});
            rx(k).id = '';  rx(k).name = '';  rx(k).ecFile = '';  rx(k).names = '';
        elseif inRx
            f = regexp(L, '^\s*(ID|Names|Name|EC)\s*:\s*(.*?)\s*$', 'tokens', 'once');
            if isempty(f), continue; end
            switch f{1}
                case 'ID',    rx(k).id     = f{2};
                case 'Name',  rx(k).name   = f{2};
                case 'EC',    rx(k).ecFile = f{2};
                case 'Names', rx(k).names  = f{2};
            end
        end
    end

    for i = 1:numel(rx)
        [rx(i).subs, rx(i).tagged] = parseSubstrates(rx(i).names, 1);
        [rx(i).prods, rx(i).prodTagged] = parseSubstrates(rx(i).names, 2);
    end
end

function [subs, tagged] = parseSubstrates(namesLine, side)
% One side (1 = left, 2 = right) of "A(m_in) + B  ->  C + D"; role tags such as (m_in) are removed.
    subs = {};  tagged = false(1, 0);
    if isempty(namesLine), return; end
    sides = regexp(namesLine, '\s*(?:<=>|<->|->|=>)\s*', 'split');
    if numel(sides) < side, return; end
    terms = regexp(strtrim(sides{side}), '\s\+\s', 'split');
    for i = 1:numel(terms)
        s = strtrim(terms{i});
        if isempty(s), continue; end
        t = regexp(s, '^(.*?)\((m_\w+|adapter_\d+|[^()]*input)\)$', 'tokens', 'once');
        if ~isempty(t), s = strtrim(t{1}); end
        s = regexprep(s, '^\d+(\.\d+)?\s+', '');          % stoichiometric coefficient
        subs{end+1}   = s;                                %#ok<AGROW>
        tagged(end+1) = ~isempty(t);                      %#ok<AGROW>
    end
end

function rx = applyEcFixes(rx, cfg)
    for i = 1:numel(rx)
        ecs = regexp(rx(i).ecFile, '\d+\.\d+\.\d+\.\d+', 'match');
        rx(i).ec = '';  rx(i).ecNote = '';
        if ~isempty(ecs)
            rx(i).ec = ecs{1};
            if numel(ecs) > 1
                rx(i).ecNote = ['report lists ' strjoin(ecs, ', ') '; first one used'];
            end
        else
            rx(i).ecNote = 'no EC number in report';
        end
        j = find(strcmp(cfg.ecFix(:, 1), rx(i).id), 1);
        if ~isempty(j)
            if cfg.ecFix{j, 4}
                rx(i).ecNote = sprintf(['EC CORRECTED: the model annotation says %s, which is wrong for this reaction; ' ...
                    '%s. Values in this row are for the corrected EC'], rx(i).ec, cfg.ecFix{j, 3});
                rx(i).ec = cfg.ecFix{j, 2};
            else
                rx(i).ecNote = sprintf(['EC %s FROM THE MODEL ANNOTATION IS WRONG for this reaction: %s. ' ...
                    'Values in this row are not for this enzyme'], rx(i).ec, cfg.ecFix{j, 3});
            end
        end
    end
end

% =========================================================================
%  BRENDA
% =========================================================================
function data = fetchAll(ecs, email, password, cfg)
    cache = struct();
    if cfg.useCache && exist(cfg.cacheFile, 'file') == 2
        S = load(cfg.cacheFile);
        if isfield(S, 'cache'), cache = S.cache; end
    end
    data = struct();
    for i = 1:numel(ecs)
        key = ecKey(ecs{i});
        if isfield(cache, key)
            fprintf('  EC %-12s from %s\n', ecs{i}, cfg.cacheFile);
            data.(key) = tidy(cache.(key));
            continue
        end
        fprintf('  EC %-12s querying BRENDA (%d of %d)\n', ecs{i}, i, numel(ecs));
        try
            out = brendaQuery(ecs{i}, email, password);        % no organism = all organisms
            data.(key) = tidy(out);
            D = data.(key);
            if cfg.useCache && (numel(D.km.val) + numel(D.kc.val) + numel(D.sa.val) + numel(D.mw.val)) > 0
                cache.(key) = out;                             % empty replies are not cached
                save(cfg.cacheFile, 'cache');
            end
        catch err
            warning('brendaKineticsTable:queryFailed', 'EC %s: %s', ecs{i}, err.message);
            data.(key) = tidy(struct());
            data.(key).error = err.message;
        end
        pause(1.1);                                            % BRENDA: at most one request per second
    end
end

function key = ecKey(ec)
    key = ['ec_' strrep(ec, '.', '_')];
end

function D = tidy(out)
% brendaQuery output (tables) -> plain columns, with the per-row flags used below.
    D.error = '';
    KM = getPart(out, 'KM');  SA = getPart(out, 'SpecificActivity');  MW = getPart(out, 'MolecularWeight');
    TN = getPart(out, 'TurnoverNumber');

    n = nRows(TN);
    D.kc.val    = rangeValue(numCol(TN, 'turnoverNumber', n), numCol(TN, 'turnoverNumberMaximum', n));
    D.kc.subKey = normName(txtCol(TN, 'substrate', n));
    D.kc.org    = txtCol(TN, 'organism', n);
    D.kc.ok     = D.kc.val > 0 & isfinite(D.kc.val) & ~isMutant(txtCol(TN, 'commentary', n));

    n = nRows(KM);
    D.km.val    = rangeValue(numCol(KM, 'kmValue', n), numCol(KM, 'kmValueMaximum', n));
    D.km.sub    = txtCol(KM, 'substrate', n);
    D.km.subKey = normName(D.km.sub);
    D.km.org    = txtCol(KM, 'organism', n);
    D.km.ok     = D.km.val > 0 & isfinite(D.km.val) & ~isMutant(txtCol(KM, 'commentary', n));

    n = nRows(SA);
    D.sa.val = rangeValue(numCol(SA, 'specificActivity', n), numCol(SA, 'specificActivityMaximum', n));
    D.sa.org = txtCol(SA, 'organism', n);
    D.sa.ok  = D.sa.val > 0 & isfinite(D.sa.val) & ~isMutant(txtCol(SA, 'commentary', n));

    n = nRows(MW);
    com = txtCol(MW, 'commentary', n);
    D.mw.val     = rangeValue(numCol(MW, 'molecularWeight', n), numCol(MW, 'molecularWeightMaximum', n));
    D.mw.org     = txtCol(MW, 'organism', n);
    D.mw.ok      = D.mw.val > 0 & isfinite(D.mw.val);
    D.mw.subunit = ~cellfun('isempty', regexpi(com, ...
                     'sds|sequence|subunit|monomer|calculat|deduced|\d\s*\*\s*\d', 'once'));
end

function T = getPart(out, name)
    T = [];
    if isstruct(out) && isfield(out, name), T = out.(name); end
end

function n = nRows(T)
    n = 0;
    if isa(T, 'table')
        n = size(T, 1);
    elseif isstruct(T)
        f = fieldnames(T);
        if ~isempty(f), n = numel(T.(f{1})); end
    end
end

function tf = hasCol(T, name)
    if isa(T, 'table')
        tf = any(strcmp(T.Properties.VariableNames, name));
    else
        tf = isstruct(T) && isfield(T, name);
    end
end

function v = numCol(T, name, n)
    v = nan(n, 1);
    if n == 0 || ~hasCol(T, name), return; end
    x = T.(name);
    if ~isnumeric(x), x = str2double(cellstr(x)); end
    v = double(x(:));
end

function v = txtCol(T, name, n)
    v = repmat({''}, n, 1);
    if n == 0 || ~hasCol(T, name), return; end
    v = cellstr(T.(name));
    v = v(:);
end

function v = rangeValue(lo, hi)
% BRENDA gives either one value or a range lo..hi; a range is replaced by its midpoint.
    v = lo;
    isRange = hi > lo & lo > 0 & isfinite(hi);
    v(isRange) = (lo(isRange) + hi(isRange)) / 2;
end

function tf = isMutant(commentary)
% Engineered or mutated enzymes: the word "mutant"/"mutated", or a substitution such as K45A.
    kw = ~cellfun('isempty', regexpi(commentary, 'mutant|mutat', 'once'));
    aa = '[ACDEFGHIKLMNPQRSTVWY]';
    pm = ~cellfun('isempty', regexp(commentary, ...
            ['(?<![A-Za-z0-9])' aa '\d{2,4}' aa '(?![A-Za-z0-9])'], 'once'));
    tf = kw | pm;
end

function s = normName(s)
% Lower case, letters and digits only: "D-Fructose 6-Phosphate" -> "dfructose6phosphate".
    s = lower(regexprep(s, '[^A-Za-z0-9]', ''));
end

% =========================================================================
%  Choosing the values
% =========================================================================
function r = chooseValues(r, data, cfg)
    nSub = numel(r.subs);
    r.km = repmat(emptyPick(''), 1, nSub);
    r.kcat = r.km;
    r.sa = emptyPick('');  r.mw = emptyPick('');
    r.primary = primarySubstrate(r, cfg);

    if isempty(r.ec)
        for j = 1:nSub, r.km(j) = emptyPick('not queried'); end
        r.kcat = r.km;
        r.sa = emptyPick('not queried');  r.mw = emptyPick('not queried');
        return
    end
    D = data.(ecKey(r.ec));
    if ~isempty(D.error)
        msg = ['BRENDA error: ' D.error];
        for j = 1:nSub, r.km(j) = emptyPick(msg); end
        r.kcat = r.km;
        r.sa = emptyPick(msg);  r.mw = emptyPick(msg);
        return
    end

    org   = cfg.organism;
    never = normName(cfg.neverSubstrate);

    % K_M and k_cat: substrate first, organism second
    isOrg = strcmpi(D.km.org, org);
    isOrgK = strcmpi(D.kc.org, org);
    for j = 1:nSub
        if any(strcmp(normName(r.subs{j}), never))
            r.km(j) = emptyPick('not a K_M substrate');
            r.kcat(j) = r.km(j);
            continue
        end
        isSub = ismember(D.km.subKey, synonymKeys(r.subs{j}, cfg));
        levels = { ...
            D.km.ok & isSub & isOrg,  org; ...
            D.km.ok & isSub,          'other organisms'; ...
            D.km.ok & isOrg,          [org ', OTHER SUBSTRATES']; ...
            D.km.ok,                  'other organisms, OTHER SUBSTRATES'};
        if ~cfg.kmOtherSubstrates, levels = levels(1:2, :); end
        r.km(j) = pickFirst(D.km.val, cfg.kmStat, 'no K_M in BRENDA for this substrate', levels);

        isSub = ismember(D.kc.subKey, synonymKeys(r.subs{j}, cfg));
        levels = { ...
            D.kc.ok & isSub & isOrgK,  org; ...
            D.kc.ok & isSub,           'other organisms'; ...
            D.kc.ok & isOrgK,          [org ', OTHER SUBSTRATES']; ...
            D.kc.ok,                   'other organisms, OTHER SUBSTRATES'};
        if ~cfg.kmOtherSubstrates, levels = levels(1:2, :); end
        r.kcat(j) = pickFirst(D.kc.val, cfg.kcatStat, 'no k_cat in BRENDA for this substrate', levels);
    end

    % Specific activity
    isOrg = strcmpi(D.sa.org, org);
    r.sa = pickFirst(D.sa.val, cfg.saStat, 'no SA in BRENDA', { ...
        D.sa.ok & isOrg,  org; ...
        D.sa.ok,          'other organisms'});

    % Molecular weight [Da]: subunit mass first, organism second
    isOrg = strcmpi(D.mw.org, org);
    r.mw = pickFirst(D.mw.val, cfg.mwStat, 'no MW in BRENDA', { ...
        D.mw.ok & D.mw.subunit & isOrg,  [org ', subunit']; ...
        D.mw.ok & D.mw.subunit,          'other organisms, subunit'; ...
        D.mw.ok & isOrg,                 [org ', native or unspecified']; ...
        D.mw.ok,                         'other organisms, native or unspecified'});
end

function rx = findRemovers(rx, cfg)
% For every reaction: its product (the motif species the report tags on the
% right-hand side) and the reactions of the same motif that consume it.
    skip = [normName(cfg.neverSubstrate), normName(cfg.currency)];
    for i = 1:numel(rx)
        rx(i).product = '';  rx(i).removers = [];  rx(i).removerOtherKm = [];
        key = normName(rx(i).prods);
        k = find(rx(i).prodTagged & ~ismember(key, skip), 1);
        if isempty(k), continue; end
        rx(i).product = rx(i).prods{k};
        for j = 1:numel(rx)
            if j == i || rx(j).motif ~= rx(i).motif, continue; end
            if any(rx(j).tagged & strcmp(normName(rx(j).subs), key{k}))
                rx(i).removers(end+1) = j;
                % true when row j's K_M is for another of its substrates
                other = rx(j).primary == 0 || ~strcmp(normName(rx(j).subs{rx(j).primary}), key{k});
                rx(i).removerOtherKm(end+1) = other;
            end
        end
    end
end

function p = emptyPick(basis)
    p = struct('value', NaN, 'n', 0, 'min', NaN, 'median', NaN, 'max', NaN, 'basis', basis);
end

function p = pickFirst(val, how, noneText, levels)
% levels: {row mask, description}; the first mask with any rows wins.
    p = emptyPick(noneText);
    for i = 1:size(levels, 1)
        v = val(levels{i, 1});
        if isempty(v), continue; end
        p.n = numel(v);  p.min = min(v);  p.median = median(v);  p.max = max(v);
        switch how
            case 'median',  p.value = p.median;
            case 'max',     p.value = p.max;
            case 'min',     p.value = p.min;
            case 'geomean', p.value = exp(mean(log(v)));
            otherwise, error('brendaKineticsTable:stat', 'Unknown statistic "%s".', how);
        end
        p.basis = levels{i, 2};
        return
    end
end

function keys = synonymKeys(name, cfg)
    key  = normName(name);
    keys = {key};
    for i = 1:size(cfg.synonyms, 1)
        if strcmp(normName(cfg.synonyms{i, 1}), key)
            more = normName(cfg.synonyms{i, 2});
            keys = [keys, more(:).'];                     %#ok<AGROW>
        end
    end
end

function j = primarySubstrate(r, cfg)
% The substrate whose K_M goes on the "values" sheet: the motif species the
% report tags on the left-hand side, cofactors last, protons and water never.
    j = 0;
    if isempty(r.subs), return; end
    key   = normName(r.subs);
    never = ismember(key, normName(cfg.neverSubstrate));
    cur   = ismember(key, normName(cfg.currency));
    order = { r.tagged & ~cur & ~never, r.tagged & cur & ~never, ~cur & ~never, ~never };
    for i = 1:numel(order)
        k = find(order{i}, 1);
        if ~isempty(k), j = k; return; end
    end
end

% =========================================================================
%  Workbook
% =========================================================================
function writeWorkbook(xlsxFile, rx, cfg)
    manyMotifs = numel(unique([rx.motif])) > 1;

    % ---- sheet 1: values (same layout as the template) -------------------
    % styles: 0 plain, 1 centred + thin line below, 2 header (medium line below),
    %         3 header centred, 4 group label, 5 thin line below, 6 wrapped text,
    %         7 wrapped text + thin line below
    c = {};                                             % rows of {ref, kind, value, style}
    c(end+1, :) = {'E1', 's', 'Constants', 1};
    for col = 'FGHIJKLM', c(end+1, :) = {[col '1'], 's', '', 1}; end       %#ok<AGROW>
    hdr = {'Reaction pair', '', 'Enzyme', 'EC code', 'K_M [mM]', 'SA [{mu}mol/min/mg]', ...
            'MW [mg/{mu}mol]', 'k_cat [s^-1]', '[E] [mM]', 'k = k_cat [E] / K_M [s^-1]', ...
            'Removed by (* = K_M there is for another substrate)', 'G', '{tau} [s]', 'Notes'};
    for j = 1:numel(hdr)
        c(end+1, :) = {[colName(j) '2'], 's', hdr{j}, 2 + (j <= 2)}; %#ok<AGROW>
    end
    merges = {'E1:M1', 'A2:B2'};

    row = 2;  i = 1;
    while i <= numel(rx)
        % rows i..last form one group
        last = i;
        while last < numel(rx) && rx(last+1).motif == rx(i).motif && strcmp(rx(last+1).group, rx(i).group)
            last = last + 1;
        end
        first = row + 1;
        for q = i:last
            row = row + 1;  r = rx(q);  n = sprintf('%d', row);
            under = (q == last);                        % line under the last row of a group
            if q == i
                c(end+1, :) = {['A' n], 's', groupTitle(r, cfg, manyMotifs), 4}; %#ok<AGROW>
            else
                c(end+1, :) = {['A' n], 's', '', 4};    %#ok<AGROW>
            end
            c(end+1, :) = {['B' n], 's', r.label, 5 * under};         %#ok<AGROW>
            c(end+1, :) = {['C' n], 's', r.name,  5 * under};         %#ok<AGROW>
            c(end+1, :) = {['D' n], 's', r.ec,    5 * under};         %#ok<AGROW>
            km = NaN;
            if r.primary > 0, km = r.km(r.primary).value; end
            mw = r.mw.value / 1000;                     % Da -> mg/umol
            c(end+1, :) = {['E' n], 'n', km,         5 * under};      %#ok<AGROW>
            c(end+1, :) = {['F' n], 'n', r.sa.value, 5 * under};      %#ok<AGROW>
            c(end+1, :) = {['G' n], 'n', mw,         5 * under};      %#ok<AGROW>
            % k_cat [1/s] = SA [umol/min/mg] * MW [mg/umol] / 60
            c(end+1, :) = {['H' n], 'f', {sprintf('IF(AND(ISNUMBER(F%s),ISNUMBER(G%s)),F%s*G%s/60,"")', n, n, n, n), ...
                                          r.sa.value * mw / 60}, 5 * under}; %#ok<AGROW>
            c(end+1, :) = {['I' n], 'n', 0, 5 * under};               %#ok<AGROW>   % [E]: 0 for now, type the real value here
            % k = k_cat [E] / K_M
            kcN = kcatRef(cfg, row);
            f = sprintf('IF(AND(ISNUMBER(%s),ISNUMBER(I%s),ISNUMBER(E%s)),IF(E%s>0,%s*I%s/E%s,""),"")', kcN, n, n, n, kcN, n, n);
            k0 = NaN;
            if km > 0 && isfinite(kcatValue(cfg, r, mw)), k0 = 0; end   % value shown while [E] = 0
            c(end+1, :) = {['J' n], 'f', {f, k0}, 5 * under};         %#ok<AGROW>
            % removal reactions: reaction q is on row q + 2
            names = cell(1, numel(r.removers));  refs = names;
            for t = 1:numel(r.removers)
                names{t} = rx(r.removers(t)).label;
                if r.removerOtherKm(t), names{t} = [names{t} '*']; end
                refs{t} = sprintf('J%d', r.removers(t) + 2);
            end
            c(end+1, :) = {['K' n], 's', strjoin(names, ', '), 5 * under}; %#ok<AGROW>
            if isempty(refs)
                c(end+1, :) = {['L' n], 's', '', 5 * under};          %#ok<AGROW>
                c(end+1, :) = {['M' n], 's', '', 5 * under};          %#ok<AGROW>
            elseif numel(refs) == 1
                % one removal reaction on row rr: Eqs. (26) and (27) written out
                %   G   = k_cat,f [E]f K_M,r / (k_cat,r [E]r K_M,f)
                %   tau = K_M,r / (k_cat,r [E]r)
                rr = r.removers(1) + 2;  kcR = kcatRef(cfg, rr);
                fG = sprintf('IFERROR((%s*I%s*E%d)/(%s*I%d*E%s),"")', kcN, n, rr, kcR, rr, n);
                fT = sprintf('IFERROR(E%d/(%s*I%d),"")', rr, kcR, rr);
                c(end+1, :) = {['L' n], 'f', {fG, NaN}, 5 * under};   %#ok<AGROW>
                c(end+1, :) = {['M' n], 'f', {fT, NaN}, 5 * under};   %#ok<AGROW>
            else
                % several removal reactions: k_r is the sum of their rate constants
                allK = strjoin(refs, ',');  m = numel(refs);
                fG = sprintf('IFERROR(IF(COUNT(%s)=%d,J%s/SUM(%s),""),"")', allK, m, n, allK);
                fT = sprintf('IFERROR(IF(COUNT(%s)=%d,1/SUM(%s),""),"")', allK, m, allK);
                c(end+1, :) = {['L' n], 'f', {fG, NaN}, 5 * under};   %#ok<AGROW>
                c(end+1, :) = {['M' n], 'f', {fT, NaN}, 5 * under};   %#ok<AGROW>
            end
            c(end+1, :) = {['N' n], 's', noteText(r), 6 + under};     %#ok<AGROW>
        end
        if row > first, merges{end+1} = sprintf('A%d:A%d', first, row); end %#ok<AGROW>
        i = last + 1;
    end
    sheet1 = sheetXml(c, merges, [26 12 44 10 11 18 15 12 10 24 30 10 10 70]);

    % ---- sheet 2: K_M and k_cat for every substrate ----------------------
    c = {};  row = 1;
    hdr = {'Reaction', 'Reaction ID', 'EC used', 'Substrate', 'On values sheet', 'K_M [mM]', ...
            'n', 'min', 'median', 'max', 'K_M taken from'};
    for j = 1:numel(hdr), c(end+1, :) = {[colName(j) '1'], 's', hdr{j}, 2}; end %#ok<AGROW>
    for q = 1:numel(rx)
        r = rx(q);
        for j = 1:numel(r.subs)
            if any(strcmp(normName(r.subs{j}), normName(cfg.neverSubstrate))), continue; end
            row = row + 1;  p = r.km(j);  mark = '';
            if j == r.primary, mark = 'yes'; end
            vals = {r.label, r.id, r.ec, r.subs{j}, mark, p.value, p.n, p.min, p.median, p.max, p.basis};
            c = [c; rowCells(row, vals)];               %#ok<AGROW>
        end
    end
    sheet2 = sheetXml(c, {}, [14 12 10 48 14 11 6 11 11 11 44]);

    % ---- sheet 3: SA and MW ----------------------------------------------
    c = {};  row = 1;
    hdr = {'Reaction', 'Reaction ID', 'EC in model', 'EC used', ...
            'SA [{mu}mol/min/mg]', 'n', 'min', 'median', 'max', 'SA taken from', ...
            'MW [Da]', 'n', 'min', 'median', 'max', 'MW taken from', 'EC note'};
    for j = 1:numel(hdr), c(end+1, :) = {[colName(j) '1'], 's', hdr{j}, 2}; end %#ok<AGROW>
    for q = 1:numel(rx)
        r = rx(q);  row = row + 1;
        vals = {r.label, r.id, r.ecFile, r.ec, ...
                r.sa.value, r.sa.n, r.sa.min, r.sa.median, r.sa.max, r.sa.basis, ...
                r.mw.value, r.mw.n, r.mw.min, r.mw.median, r.mw.max, r.mw.basis, r.ecNote};
        c = [c; rowCells(row, vals)];                   %#ok<AGROW>
    end
    sheet3 = sheetXml(c, {}, [14 12 12 10 18 6 11 11 11 30 11 6 11 11 11 40 70]);

    % ---- package ----------------------------------------------------------
    ns   = 'http://schemas.openxmlformats.org/';
    decl = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>';
    ws   = 'application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml';
    types = [decl '<Types xmlns="' ns 'package/2006/content-types">' ...
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>' ...
        '<Default Extension="xml" ContentType="application/xml"/>' ...
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>' ...
        '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="' ws '"/>' ...
        '<Override PartName="/xl/worksheets/sheet2.xml" ContentType="' ws '"/>' ...
        '<Override PartName="/xl/worksheets/sheet3.xml" ContentType="' ws '"/>' ...
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>' ...
        '</Types>'];
    rels = [decl '<Relationships xmlns="' ns 'package/2006/relationships">' ...
        '<Relationship Id="rId1" Type="' ns 'officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>' ...
        '</Relationships>'];
    workbook = [decl '<workbook xmlns="' ns 'spreadsheetml/2006/main" xmlns:r="' ns 'officeDocument/2006/relationships">' ...
        '<sheets><sheet name="values" sheetId="1" r:id="rId1"/>' ...
        '<sheet name="KM details" sheetId="2" r:id="rId2"/>' ...
        '<sheet name="SA MW details" sheetId="3" r:id="rId3"/></sheets>' ...
        '<calcPr fullCalcOnLoad="1"/></workbook>'];
    rt = [ns 'officeDocument/2006/relationships/'];
    wbRels = [decl '<Relationships xmlns="' ns 'package/2006/relationships">' ...
        '<Relationship Id="rId1" Type="' rt 'worksheet" Target="worksheets/sheet1.xml"/>' ...
        '<Relationship Id="rId2" Type="' rt 'worksheet" Target="worksheets/sheet2.xml"/>' ...
        '<Relationship Id="rId3" Type="' rt 'worksheet" Target="worksheets/sheet3.xml"/>' ...
        '<Relationship Id="rId4" Type="' rt 'styles" Target="styles.xml"/>' ...
        '</Relationships>'];
    names = {'[Content_Types].xml', '_rels/.rels', 'xl/workbook.xml', 'xl/_rels/workbook.xml.rels', ...
             'xl/styles.xml', 'xl/worksheets/sheet1.xml', 'xl/worksheets/sheet2.xml', 'xl/worksheets/sheet3.xml'};
    parts = {types, rels, workbook, wbRels, stylesXml(), sheet1, sheet2, sheet3};

    % Written to a temporary file first: Java resolves relative paths against
    % its own start-up folder, not against MATLAB's current folder.
    tmp = [tempname '.xlsx'];
    writeZip(tmp, names, parts);
    [ok, msg] = movefile(tmp, xlsxFile, 'f');
    if ~ok
        error('brendaKineticsTable:write', 'Could not write %s (is it open in Excel?): %s', xlsxFile, msg);
    end
end

function ref = kcatRef(cfg, row) %#ok<INUSL>
% Cell holding k_cat (= SA x MW / 60) for a sheet row.
    ref = sprintf('H%d', row);
end

function v = kcatValue(cfg, r, mw) %#ok<INUSL>
    v = r.sa.value * mw / 60;
end

function t = groupTitle(r, cfg, manyMotifs)
    t = r.group;
    j = find(strcmpi(cfg.groupTitle(:, 1), r.group), 1);
    if ~isempty(j), t = cfg.groupTitle{j, 2}; end
    if manyMotifs, t = sprintf('Motif %d: %s', r.motif, t); end
end

function s = noteText(r)
    if isempty(r.ec), s = ['Not queried. ' r.ecNote]; return; end
    p = {};
    if ~isempty(r.ecNote), p{end+1} = r.ecNote; end
    if r.primary > 0
        k = r.km(r.primary);
        p{end+1} = sprintf('K_M for %s: %s', r.subs{r.primary}, withN(k));
    end
    p{end+1} = ['SA: ' withN(r.sa)];
    p{end+1} = ['MW: ' withN(r.mw)];
    s = strjoin(p, ' | ');
end

function s = withN(p)
    s = p.basis;
    if p.n > 0, s = sprintf('%s (n=%d)', p.basis, p.n); end
end

function c = rowCells(row, vals)
    c = cell(numel(vals), 4);
    for j = 1:numel(vals)
        kind = 's';
        if isnumeric(vals{j}), kind = 'n'; end
        c(j, :) = {sprintf('%s%d', colName(j), row), kind, vals{j}, 0};
    end
end

function s = colName(j)
    s = char('A' + j - 1);                               % A..Z is enough here
end

function xml = sheetXml(cells, merges, widths)
% cells: rows of {ref, kind, value, style}; kind is 's' text, 'n' number,
% 'f' formula with value = {formula, cached result}. Cells must be in row order.
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main';
    out = cell(1, size(cells, 1) + 8);  k = 0;
    k = k + 1;  out{k} = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="' ns '"><cols>'];
    w = '';
    for j = 1:numel(widths)
        w = [w sprintf('<col min="%d" max="%d" width="%g" customWidth="1"/>', j, j, widths(j))]; %#ok<AGROW>
    end
    k = k + 1;  out{k} = [w '</cols><sheetData>'];

    current = 0;
    for i = 1:size(cells, 1)
        ref = cells{i, 1};  kind = cells{i, 2};  v = cells{i, 3};  st = cells{i, 4};
        rowNum = str2double(regexp(ref, '\d+$', 'match', 'once'));
        s = '';
        if rowNum ~= current
            if current > 0, s = '</row>'; end
            s = [s sprintf('<row r="%d">', rowNum)];     %#ok<AGROW>
            current = rowNum;
        end
        a = sprintf(' r="%s"', ref);
        if st > 0, a = [a sprintf(' s="%d"', st)]; end   %#ok<AGROW>
        switch kind
            case 'n'
                if isempty(v) || ~isfinite(v)
                    if st > 0, s = [s '<c' a '/>']; end  %#ok<AGROW>   % empty cell, keeps its border
                else
                    s = [s '<c' a '><v>' sprintf('%.15g', v) '</v></c>']; %#ok<AGROW>
                end
            case 'f'
                if isfinite(v{2})
                    s = [s '<c' a '><f>' xmlText(v{1}) '</f><v>' sprintf('%.15g', v{2}) '</v></c>']; %#ok<AGROW>
                else
                    s = [s '<c' a ' t="str"><f>' xmlText(v{1}) '</f><v></v></c>']; %#ok<AGROW>
                end
            otherwise
                if isempty(v)
                    if st > 0, s = [s '<c' a '/>']; end  %#ok<AGROW>
                else
                    s = [s '<c' a ' t="inlineStr"><is><t xml:space="preserve">' xmlText(v) '</t></is></c>']; %#ok<AGROW>
                end
        end
        k = k + 1;  out{k} = s;
    end
    tail = '';
    if current > 0, tail = '</row>'; end
    tail = [tail '</sheetData>'];
    if ~isempty(merges)
        tail = [tail sprintf('<mergeCells count="%d">', numel(merges))];
        for j = 1:numel(merges), tail = [tail '<mergeCell ref="' merges{j} '"/>']; end %#ok<AGROW>
        tail = [tail '</mergeCells>'];
    end
    k = k + 1;  out{k} = [tail '</worksheet>'];
    xml = [out{1:k}];
end

function xml = stylesXml()
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main';
    bd   = @(style) ['<border><left/><right/><top/><bottom style="' style '"><color auto="1"/></bottom><diagonal/></border>'];
    xf   = @(border, align) ['<xf numFmtId="0" fontId="0" fillId="0" borderId="' border '" xfId="0"' ...
                             ' applyBorder="1" applyAlignment="1">' align '</xf>'];
    xml = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?><styleSheet xmlns="' ns '">' ...
        '<fonts count="1"><font><sz val="11"/><name val="Calibri"/><family val="2"/></font></fonts>' ...
        '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>' ...
        '<borders count="3"><border><left/><right/><top/><bottom/><diagonal/></border>' bd('thin') bd('medium') '</borders>' ...
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>' ...
        '<cellXfs count="8">' ...
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>' ...
        xf('1', '<alignment horizontal="center"/>') ...
        xf('2', '<alignment/>') ...
        xf('2', '<alignment horizontal="center"/>') ...
        xf('1', '<alignment horizontal="center" vertical="center" wrapText="1"/>') ...
        xf('1', '<alignment/>') ...
        xf('0', '<alignment vertical="top" wrapText="1"/>') ...
        xf('1', '<alignment vertical="top" wrapText="1"/>') ...
        '</cellXfs>' ...
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>' ...
        '</styleSheet>'];
end

function s = xmlText(s)
% Escape for XML and keep the result pure ASCII: characters above 126 become
% numeric references, so the file does not depend on any text encoding.
% {mu} and {tau} stand for the Greek letters in the headers.
    s = char(s);
    s = s(:).';
    code = double(s);
    s(code < 32 & code ~= 9 & code ~= 10) = ' ';
    s = strrep(s, '&', '&amp;');
    s = strrep(s, '<', '&lt;');
    s = strrep(s, '>', '&gt;');
    code = double(s);
    if any(code > 126)
        piece = cell(1, numel(s));
        for i = 1:numel(s)
            if code(i) <= 126
                piece{i} = s(i);
            elseif (code(i) >= 55296 && code(i) <= 57343) || code(i) >= 65534
                piece{i} = '?';
            else
                piece{i} = sprintf('&#%d;', code(i));
            end
        end
        s = [piece{:}];
    end
    s = strrep(s, '{mu}', '&#956;');
    s = strrep(s, '{tau}', '&#964;');
end

function writeZip(zipFile, names, contents)
% An .xlsx file is a zip archive of XML parts.
    fos = javaObject('java.io.FileOutputStream', zipFile);
    zos = javaObject('java.util.zip.ZipOutputStream', fos);
    for i = 1:numel(names)
        zos.putNextEntry(javaObject('java.util.zip.ZipEntry', names{i}));
        bytes = int8(contents{i});                       % ASCII only, see xmlText
        zos.write(bytes, int32(0), int32(numel(bytes)));
        zos.closeEntry();
    end
    zos.close();
end

% =========================================================================
function printSummary(rx)
    fprintf('%-11s %-10s %-12s %-12s %-12s %-12s\n', 'Reaction', 'EC', 'K_M [mM]', 'SA', 'MW [Da]', 'k_cat [1/s]');
    for i = 1:numel(rx)
        r = rx(i);  km = NaN;
        if r.primary > 0, km = r.km(r.primary).value; end
        fprintf('%-11s %-10s %-12s %-12s %-12s %-12s\n', r.label, r.ec, num(km), num(r.sa.value), ...
                num(r.mw.value), num(r.sa.value * r.mw.value / 1000 / 60));
    end
    fprintf('\n');
end

function s = num(x)
    if isfinite(x), s = sprintf('%.4g', x); else, s = '-'; end
end
