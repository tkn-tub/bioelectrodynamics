function reg = parse_brenda_flat(flatFile, fields)
% PARSE_BRENDA_FLAT  Parse the BRENDA flat text file (release 2026.1 format) into
% a map  EC -> struct with one cellstr of NORMALISED ligand names per field.
%
%   reg = parse_brenda_flat('brenda_2026_1.txt')         % AC + IN
%   reg = parse_brenda_flat(file, {'AC'})                % activators only
%   reg('1.1.1.27').AC   -> {'fructose 1 6 bisphosphate', 'oxamate', ...}
%
% Format (per the 2026.1 README):
%   * EC block ends with a line '///'.
%   * Each line: 2/3-letter acronym, then a TAB, then the content.
%   * A line beginning with whitespace is a CONTINUATION of the previous line.
%   * In the content: protein refs in #...#, commentary in (...),
%     field-special in {...}, citations in <...>.
%   For AC/IN the ligand name is the text after the #...# tag and before the
%   first '(', '{' or '<'. BRENDA placeholders ('more','additional information')
%   are skipped.

    if nargin < 2, fields = {'AC','IN'}; end
    skipLig = {'more','additional information',''};
    reg = containers.Map('KeyType','char','ValueType','any');
    fid = fopen(flatFile,'r','n','UTF-8');
    assert(fid>=0, 'Cannot open BRENDA flat file: %s', flatFile);
    clean = onCleanup(@() fclose(fid));

    ec = ''; buf = '';                       % buf = current logical (joined) line
    TAB = char(9);
    while true
        raw = fgetl(fid);
        atEnd = ~ischar(raw);
        isCont = ~atEnd && ~isempty(raw) && (raw(1)==' ' || raw(1)==TAB);
        if isCont
            buf = [buf ' ' strtrim(raw)];    %#ok<AGROW>  join continuation
            continue;
        end
        % a new logical line (or EOF / '///') starts -> flush the previous buffer
        [ec, reg] = flushLine(buf, ec, reg, fields, skipLig, TAB);
        if atEnd, break; end
        if strcmp(strtrim(raw),'///'), buf = ''; ec_done = true; continue; end %#ok<NASGU>
        buf = raw;
    end

    % de-duplicate
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
    tag = line(1:t(1)-1);
    body = line(t(1)+1:end);
    if strcmp(tag,'ID')
        ec = strtrim(regexprep(body,'\s.*',''));          % EC up to first space
        if ~isKey(reg,ec)
            s = struct(); for f=1:numel(fields), s.(fields{f})={}; end
            reg(ec) = s;
        end
        return;
    end
    if any(strcmp(tag,fields)) && ~isempty(ec)
        body = regexprep(body,'^#[^#]*#\s*','');           % strip #protein refs#
        lig  = regexp(body,'^[^(\{<]*','match','once');     % up to first ( { <
        lig  = normLigand(lig);
        if ~any(strcmp(lig, skipLig))
            s = reg(ec); s.(tag){end+1} = lig; reg(ec) = s; %#ok<AGROW>
        end
    end
end

function s = normLigand(s)
    s = lower(strtrim(char(s)));
    s = regexprep(s,'^((l|d|dl|\(s\)|\(r\)|alpha|beta)-)+','');
    s = regexprep(s,'[^a-z0-9]+',' ');
    s = strtrim(s);
end
