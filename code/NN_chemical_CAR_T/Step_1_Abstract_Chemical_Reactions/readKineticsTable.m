function T = readKineticsTable(xlsxFile, E, kcatSource)
%READKINETICSTABLE  Import the "values" sheet written by brendaKineticsTable
%   and compute the rate constant, gain and time constant of every reaction.
%
%   T = readKineticsTable('kinetics_BRENDA.xlsx')
%       uses the enzyme concentrations typed into the [E] column of the workbook
%   T = readKineticsTable('kinetics_BRENDA.xlsx', 1e-3)
%       uses [E] = 1e-3 mM for every reaction
%   T = readKineticsTable('kinetics_BRENDA.xlsx', Evec)
%       one [E] per reaction, in the row order of the sheet (NaN = take the sheet value)
%   T = readKineticsTable(..., 'BRENDA')
%       which k_cat to use: 'SAxMW' (default, as in the paper) | 'BRENDA' | 'BRENDA else SAxMW'
%
%   Everything is computed here in MATLAB from K_M, SA, MW, k_cat and [E];
%   the formula cells of the workbook are not relied upon.
%
%     k_cat = SA * MW / 60                              [1/s]   (or the BRENDA value)
%     k     = k_cat [E] / K_M                           [1/s]   linear regime, Eq. (21)
%     tau   = 1 / k_r  = K_M,r / (k_cat,r [E]r)         [s]     Eq. (27)
%     G     = k_f / k_r
%           = k_cat,f [E]f K_M,r / (k_cat,r [E]r K_M,f)         Eq. (26)
%
%   f is the reaction on the row; r is the reaction named in the column
%   "Removed by". When several reactions are named there, k_r is the sum of
%   their rate constants.
%
%   Recover single constants by reaction name:
%     T{'r_S1', 'KM'}      T{'r_No', 'kcat'}      T('r_Ni', :)
%     KM = T.KM;  kcat = T.kcat;          % whole columns, in sheet order
%
%   T has one row per reaction: Group, Reaction, Enzyme, EC, KM, SA, MW,
%   kcat_SAxMW, kcat_BRENDA, kcat, E, k, RemovedBy, G, tau.

    if nargin < 2, E = []; end
    if nargin < 3 || isempty(kcatSource), kcatSource = 'SAxMW'; end

    raw = readcell(xlsxFile, 'Sheet', 'values');
    S = kineticsFromCells(raw, E, kcatSource);

    T = struct2table(S, 'AsArray', true);
    if numel(unique(T.Reaction)) == height(T)
        T.Properties.RowNames = T.Reaction;        % allows T{'r_S1', 'KM'} and T('r_No', :)
    end
    disp(T(:, {'Reaction', 'KM', 'kcat', 'E', 'k', 'RemovedBy', 'G', 'tau'}));
end

% =========================================================================
function S = kineticsFromCells(raw, E, kcatSource)
% raw: the sheet as a cell array (numbers, text, anything else = empty).

    % ---- locate the header row and the columns by their titles -----------
    hdrRow = 0;
    for i = 1:min(size(raw, 1), 20)
        if any(cellfun(@(x) startsWithText(x, 'K_M'), raw(i, :)))
            hdrRow = i;  break
        end
    end
    if hdrRow == 0
        error('readKineticsTable:header', 'Could not find the header row (a cell starting with "K_M").');
    end
    hdr = cellfun(@txt, raw(hdrRow, :), 'UniformOutput', false);
    col.group  = 1;
    col.label  = 2;
    col.enzyme = findCol(hdr, 'Enzyme');
    col.ec     = findCol(hdr, 'EC code');
    col.km     = findCol(hdr, 'K_M');
    col.sa     = findCol(hdr, 'SA');
    col.mw     = findCol(hdr, 'MW');
    col.kcatB  = find(strncmp(hdr, 'k_cat BRENDA', 12), 1);    % absent in newer workbooks
    col.E      = findCol(hdr, '[E]');
    col.rem    = findCol(hdr, 'Removed by');

    % ---- rows -------------------------------------------------------------
    rows = (hdrRow + 1):size(raw, 1);
    rows = rows(~cellfun('isempty', cellfun(@txt, raw(rows, col.label), 'UniformOutput', false)));
    n = numel(rows);
    if n == 0, error('readKineticsTable:empty', 'No reaction rows found.'); end

    Ein = nan(n, 1);
    if isscalar(E)
        Ein(:) = E;
    elseif ~isempty(E)
        if numel(E) ~= n
            error('readKineticsTable:E', 'E has %d values but the sheet has %d reactions.', numel(E), n);
        end
        Ein = double(E(:));
    end

    S = struct('Group', cell(n, 1));
    group = '';
    for i = 1:n
        r = rows(i);
        g = txt(raw{r, col.group});
        if ~isempty(g), group = g; end                 % merged cell: the title is on its first row only
        S(i).Group       = group;
        S(i).Reaction    = txt(raw{r, col.label});
        S(i).Enzyme      = txt(raw{r, col.enzyme});
        S(i).EC          = txt(raw{r, col.ec});
        S(i).KM          = num(raw{r, col.km});         % mM
        S(i).SA          = num(raw{r, col.sa});         % umol/min/mg
        S(i).MW          = num(raw{r, col.mw});         % mg/umol
        S(i).kcat_SAxMW  = S(i).SA * S(i).MW / 60;      % 1/s
        S(i).kcat_BRENDA = NaN;                          % 1/s, only if the workbook has that column
        if ~isempty(col.kcatB), S(i).kcat_BRENDA = num(raw{r, col.kcatB}); end
        switch kcatSource
            case 'SAxMW',  S(i).kcat = S(i).kcat_SAxMW;
            case 'BRENDA', S(i).kcat = S(i).kcat_BRENDA;
            case 'BRENDA else SAxMW'
                S(i).kcat = S(i).kcat_BRENDA;
                if isnan(S(i).kcat), S(i).kcat = S(i).kcat_SAxMW; end
            otherwise
                error('readKineticsTable:kcat', 'Unknown k_cat source "%s".', kcatSource);
        end
        S(i).E = num(raw{r, col.E});                    % mM
        if ~isnan(Ein(i)), S(i).E = Ein(i); end
        S(i).k = NaN;                                   % 1/s
        if S(i).KM > 0, S(i).k = S(i).kcat * S(i).E / S(i).KM; end
        S(i).RemovedBy = txt(raw{r, col.rem});
        S(i).G = NaN;  S(i).tau = NaN;
    end

    % ---- G and tau, Eqs. (26) and (27) -------------------------------------
    labels = {S.Reaction};
    for i = 1:n
        if isempty(S(i).RemovedBy), continue; end
        names = strtrim(regexp(S(i).RemovedBy, ',', 'split'));
        names = regexprep(names, '\*$', '');            % "*" only flags the K_M substrate
        % the removal reaction must belong to the same motif: take the nearest row with that label
        kr = 0;
        for j = 1:numel(names)
            idx = find(strcmp(labels, names{j}));
            if isempty(idx)
                warning('readKineticsTable:label', '%s: removal reaction "%s" not found.', S(i).Reaction, names{j});
                kr = NaN;  break
            end
            [~, nearest] = min(abs(idx - i));
            kr = kr + S(idx(nearest)).k;
        end
        if kr > 0
            S(i).tau = 1 / kr;                          % Eq. (27)
            S(i).G   = S(i).k / kr;                     % Eq. (26)
        end
    end
end

% =========================================================================
function j = findCol(hdr, prefix)
    j = find(strncmp(hdr, prefix, numel(prefix)), 1);
    if isempty(j)
        error('readKineticsTable:column', 'Column "%s..." not found in the header row.', prefix);
    end
end

function tf = startsWithText(x, prefix)
    s = txt(x);
    tf = strncmp(s, prefix, numel(prefix));
end

function s = txt(x)
% Cell content as text ('' for numbers, empty or missing cells).
    s = '';
    if ischar(x)
        s = strtrim(x);
    elseif isa(x, 'string') && isscalar(x) && ~ismissing(x)
        s = strtrim(char(x));
    end
end

function v = num(x)
% Cell content as a number (NaN for text, empty or missing cells).
    v = NaN;
    if isnumeric(x) && isscalar(x)
        v = double(x);
    elseif ischar(x) || (isa(x, 'string') && isscalar(x) && ~ismissing(x))
        v = str2double(x);
    end
end
