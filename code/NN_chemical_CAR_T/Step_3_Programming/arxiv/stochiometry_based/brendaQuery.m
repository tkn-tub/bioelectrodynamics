function out = brendaQuery(ecNumber, email, password, organism)
%BRENDAQUERY  Query BRENDA (SOAP API) for KM, specific activity and molecular weight.
%
%   out = brendaQuery('1.1.1.1', 'you@uni.edu', 'yourPassword')
%   out = brendaQuery('1.1.1.1', 'you@uni.edu', 'yourPassword', 'Saccharomyces cerevisiae')
%
%   out.KM, out.SpecificActivity, out.MolecularWeight are tables (empty if no entries).
%   Units as in BRENDA: KM in mM, specific activity in umol/min/mg, MW in Da.
%   BRENDA uses -999 for "no numeric value" (see the commentary column).
%
%   BRENDA asks for at most one request per second; this function pauses
%   between its three requests. Add your own pause when looping over EC numbers.

    if nargin < 4, organism = ''; end
    email    = strtrim(char(email));
    password = strtrim(char(password));
    organism = strtrim(char(organism));
    ecNumber = strtrim(char(ecNumber));

    pw = sha256hex(password);      % BRENDA expects the SHA-256 hex digest, not the plain password

    % Request format: ONE string "email,hash,field*value#field*value#"
    query = ['ecNumber*' ecNumber '#'];
    if ~isempty(organism)
        query = [query 'organism*' organism '#'];
    end
    param = [email ',' pw ',' query];

    out.KM = callBrenda('getKmValue', param);
    pause(1.1);
    out.SpecificActivity = callBrenda('getSpecificActivity', param);
    pause(1.1);
    out.MolecularWeight = callBrenda('getMolecularWeight', param);
end

% -------------------------------------------------------------------------
function T = callBrenda(method, param)
    endpoint = 'https://www.brenda-enzymes.org/soap/brenda_server.php';

    envelope = ['<?xml version="1.0" encoding="UTF-8"?>' ...
        '<SOAP-ENV:Envelope' ...
        ' xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"' ...
        ' xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/"' ...
        ' xmlns:xsd="http://www.w3.org/2001/XMLSchema"' ...
        ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"' ...
        ' xmlns:ns1="https://www.brenda-enzymes.org/soap"' ...
        ' SOAP-ENV:encodingStyle="http://schemas.xmlsoap.org/soap/encoding/">' ...
        '<SOAP-ENV:Body><ns1:' method '>' ...
        '<param0 xsi:type="xsd:string">' xmlEscape(param) '</param0>' ...
        '</ns1:' method '></SOAP-ENV:Body></SOAP-ENV:Envelope>'];

    import matlab.net.http.*
    % Set Payload (raw bytes) instead of Data, so MATLAB does not try to
    % convert the body itself (for text/xml it would insist on a DOM object).
    body = MessageBody();
    body.Payload = unicode2native(envelope, 'UTF-8');
    req = RequestMessage('POST', ...
        [HeaderField('Content-Type', 'text/xml; charset=utf-8'), ...
         HeaderField('SOAPAction', '""')], ...
        body);
    opts = HTTPOptions('ConnectTimeout', 60, 'ConvertResponse', false);
    resp = req.send(endpoint, opts);

    d = resp.Body.Data;                 % text or raw bytes, depending on MATLAB version
    if isstring(d) || ischar(d)
        xml = char(strjoin(string(d), ''));
    else
        xml = native2unicode(uint8(d(:))', 'UTF-8');
    end

    if resp.StatusCode ~= StatusCode.OK || contains(xml, 'Fault>')
        msg = regexp(xml, '<faultstring[^>]*>(.*?)</faultstring>', 'tokens', 'once');
        if isempty(msg), msg = {xml(1:min(end, 500))}; end
        error('brendaQuery:soapFault', '%s failed (HTTP %d): %s', ...
              method, double(resp.StatusCode), msg{1});
    end

    T = parseResult(resultString(xml));
end

% -------------------------------------------------------------------------
function s = resultString(xml)
% Pull the returned string out of the SOAP envelope (DOM handles the XML escapes).
    doc = xmlread(org.xml.sax.InputSource(java.io.StringReader(xml)));
    nodes = doc.getElementsByTagName('return');
    if nodes.getLength == 0
        nodes = doc.getElementsByTagName('SOAP-ENV:Body');
    end
    if nodes.getLength == 0
        s = '';
    else
        s = strtrim(char(nodes.item(0).getTextContent));
    end
end

% -------------------------------------------------------------------------
function T = parseResult(s)
% BRENDA returns:  field*value#field*value#...!field*value#...!...
% Entries are separated by "!" and fields by "#". Commentaries can contain
% those characters themselves, so only split where a known field name follows.
    if isempty(s), T = table(); return; end

    keys = ['ecNumber|kmValueMaximum|kmValue|specificActivityMaximum|specificActivity|' ...
            'molecularWeightMaximum|molecularWeight|substrate|commentary|organism|' ...
            'ligandStructureId|literature'];

    entries = regexp(s, '!(?=ecNumber\*)', 'split');
    recs = {};
    for i = 1:numel(entries)
        e = regexprep(entries{i}, '[#!]+$', '');
        if isempty(e), continue; end
        fields = regexp(e, ['#(?=(?:' keys ')\*)'], 'split');
        r = struct();
        for j = 1:numel(fields)
            tok = regexp(fields{j}, '^(\w+)\*(.*)$', 'tokens', 'once');
            if isempty(tok), continue; end
            r.(matlab.lang.makeValidName(tok{1})) = string(tok{2});
        end
        if ~isempty(fieldnames(r)), recs{end+1} = r; end %#ok<AGROW>
    end
    if isempty(recs), T = table(); return; end

    % Union of fields, so records with missing fields don't break the table
    f = {};
    for i = 1:numel(recs), f = union(f, fieldnames(recs{i}), 'stable'); end
    for i = 1:numel(recs)
        for j = 1:numel(f)
            if ~isfield(recs{i}, f{j}), recs{i}.(f{j}) = ""; end
        end
        recs{i} = orderfields(recs{i}, f);
    end
    T = struct2table([recs{:}], 'AsArray', true);

    numCols = intersect(T.Properties.VariableNames, ...
        {'kmValue', 'kmValueMaximum', 'specificActivity', 'specificActivityMaximum', ...
         'molecularWeight', 'molecularWeightMaximum'});
    for j = 1:numel(numCols)
        T.(numCols{j}) = str2double(T.(numCols{j}));
    end
end

% -------------------------------------------------------------------------
function h = sha256hex(str)
    md = java.security.MessageDigest.getInstance('SHA-256');
    d  = typecast(md.digest(unicode2native(str, 'UTF-8')), 'uint8');
    h  = lower(reshape(dec2hex(d, 2)', 1, []));
end

function s = xmlEscape(s)
    s = strrep(s, '&', '&amp;');
    s = strrep(s, '<', '&lt;');
    s = strrep(s, '>', '&gt;');
end
