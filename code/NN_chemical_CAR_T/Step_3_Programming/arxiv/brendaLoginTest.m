function brendaLoginTest(email, password)
%BRENDALOGINTEST  Try BRENDA's SOAP login in every documented request style
%   and print the server's raw answer for each, to find out whether the
%   problem is the request format or the account.
%
%   brendaLoginTest('you@uni.edu', 'yourPassword')

    email    = strtrim(char(email));
    password = strtrim(char(password));

    md = java.security.MessageDigest.getInstance('SHA-256');
    d  = typecast(md.digest(unicode2native(password, 'UTF-8')), 'uint8');
    pw = lower(reshape(dec2hex(d, 2)', 1, []));

    % Checks on the inputs (nothing secret is printed)
    fprintf('Email: "%s" (%d characters)\n', email, numel(email));
    fprintf('Password: %d characters, non-ASCII: %d, characters other than letters/digits: %d\n\n', ...
        numel(password), any(password > 127), sum(~isstrprop(password, 'alphanum')));

    endpoint = 'https://www.brenda-enzymes.org/soap/brenda_server.php';
    query    = 'ecNumber*1.1.1.1#organism*Homo sapiens#';

    % Style A: one element per field (Python 3 / zeep example in the BRENDA docs)
    names  = {'email','password','ecNumber','organism','kmValue','kmValueMaximum', ...
              'substrate','commentary','ligandStructureId','literature'};
    values = {email, pw, 'ecNumber*1.1.1.1', 'organism*Homo sapiens', 'kmValue*', ...
              'kmValueMaximum*', 'substrate*', 'commentary*', 'ligandStructureId*', 'literature*'};
    bodyA = '';
    for i = 1:numel(names)
        bodyA = [bodyA sprintf('<%s xsi:type="xsd:string">%s</%s>', ...
                 names{i}, esc(values{i}), names{i})]; %#ok<AGROW>
    end

    % Style B: a single string "email,hash,query" (Perl / PHP / Python 2 examples)
    bodyB = sprintf('<param0 xsi:type="xsd:string">%s</param0>', ...
                    esc([email ',' pw ',' query]));

    % Style B with a lower-cased email, in case the address is case-sensitive
    bodyC = sprintf('<param0 xsi:type="xsd:string">%s</param0>', ...
                    esc([lower(email) ',' pw ',' query]));

    tests = { ...
      'A  named fields, namespace ws.brenda',        'ws.brenda',                            bodyA; ...
      'B  single string, namespace .../soap',        'https://www.brenda-enzymes.org/soap',  bodyB; ...
      'B2 single string, namespace ws.brenda',       'ws.brenda',                            bodyB; ...
      'C  single string, lower-case email',          'https://www.brenda-enzymes.org/soap',  bodyC};

    for t = 1:size(tests, 1)
        fprintf('--- %s ---\n', tests{t, 1});
        try
            [code, txt] = post(endpoint, tests{t, 2}, tests{t, 3});
            fault = regexp(txt, '<faultstring[^>]*>(.*?)</faultstring>', 'tokens', 'once');
            if ~isempty(fault)
                fprintf('HTTP %d, server says: %s\n\n', code, fault{1});
            else
                fprintf('HTTP %d, NO ERROR. First 300 characters of the answer:\n%s\n\n', ...
                        code, txt(1:min(end, 300)));
            end
        catch err
            fprintf('Request failed in MATLAB: %s\n\n', err.message);
        end
        pause(1.5);     % BRENDA asks for at most one request per second
    end
end

function [code, txt] = post(endpoint, ns, body)
    import matlab.net.http.*
    envelope = ['<?xml version="1.0" encoding="UTF-8"?>' ...
        '<SOAP-ENV:Envelope' ...
        ' xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"' ...
        ' xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/"' ...
        ' xmlns:xsd="http://www.w3.org/2001/XMLSchema"' ...
        ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"' ...
        ' xmlns:ns1="' ns '"' ...
        ' SOAP-ENV:encodingStyle="http://schemas.xmlsoap.org/soap/encoding/">' ...
        '<SOAP-ENV:Body><ns1:getKmValue>' body '</ns1:getKmValue>' ...
        '</SOAP-ENV:Body></SOAP-ENV:Envelope>'];
    mb = MessageBody();
    mb.Payload = unicode2native(envelope, 'UTF-8');
    req = RequestMessage('POST', ...
        [HeaderField('Content-Type', 'text/xml; charset=utf-8'), ...
         HeaderField('SOAPAction', '""')], mb);
    resp = req.send(endpoint, HTTPOptions('ConnectTimeout', 60, 'ConvertResponse', false));
    code = double(resp.StatusCode);
    d = resp.Body.Data;
    if isstring(d) || ischar(d)
        txt = char(strjoin(string(d), ''));
    else
        txt = native2unicode(uint8(d(:))', 'UTF-8');
    end
end

function s = esc(s)
    s = strrep(strrep(strrep(s, '&', '&amp;'), '<', '&lt;'), '>', '&gt;');
end
