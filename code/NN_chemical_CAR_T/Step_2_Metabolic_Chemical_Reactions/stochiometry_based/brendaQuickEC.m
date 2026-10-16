function T = brendaQuickEC(ec,organism,email,password,parameters)
%BRENDAQUICKEC  Look up one EC number in BRENDA: Km, kcat, SA and MW.
%
%   T = brendaQuickEC('4.3.2.9','',email,password)
%   T = brendaQuickEC('4.3.2.9','Homo sapiens',email,password)
%   T = brendaQuickEC(ec,organism,email,password,{'Km','kcat'})
%
%   organism = '' queries every organism.
%
%   Why this exists rather than Brenda_Kinetics.m: zeep builds its SOAP
%   operations dynamically through __getattr__, and MATLAB's dot notation
%   resolves only real attributes on a Python object. So
%
%       client.service.getKmValue(...)
%
%   raises "Unrecognized method, property, or field 'getKmValue' for class
%   py.zeep.proxy.ServiceProxy". The call has to be made INSIDE Python, with
%   only a built-in string crossing back into MATLAB. That is what
%   brendaKineticsTable does, and this is the same bridge for a single EC.
%
%   Returns one table with a Parameter column, so Km, kcat, SA and MW rows
%   sit together and are easy to filter.

if nargin < 2, organism = ''; end
if nargin < 5 || isempty(parameters), parameters = {'Km','kcat','SA','MW'}; end
ec = char(ec); organism = char(organism);
assert(~isempty(regexp(ec,'^\d+\.\d+\.\d+\.\d+$','once')), ...
    'Use a complete EC number, for example 4.3.2.9.');
assert(~isempty(email) && ~isempty(password), ...
    'Supply your registered BRENDA email and password.');

T = table();
for k = 1:numel(parameters)
    [rows,status] = localFetch(ec,organism,email,password,parameters{k});
    fprintf('  %-5s : %s\n',parameters{k},status);
    T = [T;rows]; %#ok<AGROW>
end
clear email password

if isempty(T)
    fprintf(['\nBRENDA holds no %s record for EC %s%s.\n', ...
        'Check that the code is current: a transferred entry (EC 2.3.2.4 -> 4.3.2.9, ', ...
        'say) carries no data. https://enzyme.expasy.org/EC/%s.txt says so directly.\n'], ...
        strjoin(parameters,'/'),ec,localOrganismText(organism),ec);
    return
end
fprintf('\nEC %s%s: %d record(s).\n',ec,localOrganismText(organism),height(T));
summary = groupsummary(T,{'parameter','substrate'});
disp(summary);
end

%% =======================================================================
function [T,status] = localFetch(ec,organism,email,password,parameter)
hasSubstrate = true;
switch parameter
    case 'Km'
        operation = 'getKmValue';        field = 'kmValue';
        maximum   = 'kmValueMaximum';    unit  = 'mM';
    case 'kcat'
        operation = 'getTurnoverNumber'; field = 'turnoverNumber';
        maximum   = 'turnoverNumberMaximum'; unit = '1/s';
    case 'SA'
        operation = 'getSpecificActivity'; field = 'specificActivity';
        maximum   = 'specificActivityMaximum'; unit = 'umol/min/mg';
        hasSubstrate = false;
    case 'MW'
        operation = 'getMolecularWeight'; field = 'molecularWeight';
        maximum   = 'molecularWeightMaximum'; unit = 'Da';
        hasSubstrate = false;
    otherwise
        error('Unknown parameter %s.',parameter);
end
T = localEmpty(parameter,unit);
try
    % The zeep client and its proxies never leave Python. MATLAB receives a
    % plain str holding JSON, which is why the dot-notation failure above
    % cannot happen here.
    persistent bridge
    if isempty(bridge), bridge = py.dict; end
    bridge.update(py.dict(pyargs('email',email,'password',password, ...
        'ec',ec,'organism',organism,'operation',operation, ...
        'field',field,'maximum',maximum,'has_substrate',hasSubstrate)));
    code = { ...
        'import json'; ...
        'phase = "import dependencies"'; ...
        'try:'; ...
        '    import hashlib, time, zeep'; ...
        '    from zeep.transports import Transport'; ...
        '    from zeep.helpers import serialize_object'; ...
        '    phase = "load WSDL"'; ...
        '    if "client" not in globals():'; ...
        '        client = zeep.Client("https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl", settings=zeep.Settings(strict=False), transport=Transport(timeout=60, operation_timeout=120))'; ...
        '    phase = "resolve operation"'; ...
        '    soap_call = getattr(client.service, operation)'; ...
        '    password_hash = hashlib.sha256(password.encode("utf-8")).hexdigest()'; ...
        '    phase = "execute request"'; ...
        '    time.sleep(1.1)'; ...
        '    if has_substrate:'; ...
        '        args = [email, password_hash, "ecNumber*" + ec, field + "*", maximum + "*", "substrate*", "commentary*", "organism*" + organism, "ligandStructureId*", "literature*"]'; ...
        '    else:'; ...
        '        args = [email, password_hash, "ecNumber*" + ec, field + "*", maximum + "*", "commentary*", "organism*" + organism, "literature*"]'; ...
        '    raw = soap_call(*args)'; ...
        '    phase = "serialize response"'; ...
        '    response_json = json.dumps({"ok": True, "records": serialize_object(raw)})'; ...
        'except Exception as exc:'; ...
        '    response_json = json.dumps({"ok": False, "phase": phase, "error_type": type(exc).__module__ + "." + type(exc).__name__})'; ...
        'finally:'; ...
        '    for key in ("email", "password", "password_hash", "raw", "soap_call", "args"):'; ...
        '        globals().pop(key, None)'};
    py.exec(strjoin(code,newline),bridge);
    response = jsondecode(char(bridge.get('response_json')));
    if ~response.ok
        status = sprintf('query failed during %s (%s)',response.phase,response.error_type);
        return
    end
    records = response.records;
    if isempty(records), status = 'no records'; return; end
    for i = 1:numel(records)
        if iscell(records), r = records{i}; else, r = records(i); end
        if ~isstruct(r), continue; end
        value = localText(r,field);
        T(end+1,:) = {string(ec),localText(r,'organism'),localText(r,'substrate'), ...
            str2double(value),value,localText(r,'commentary'), ...
            localText(r,'literature'),string(parameter),string(unit)}; %#ok<AGROW>
    end
    T.value(~isfinite(T.value) | T.value <= 0) = NaN;
    status = sprintf('OK: %d record(s)',height(T));
catch ME
    T = localEmpty(parameter,unit);
    status = sprintf('failed (%s); check the local Python/zeep setup',ME.identifier);
end
end

function T = localEmpty(parameter,unit)
T = table(strings(0,1),strings(0,1),strings(0,1),nan(0,1),strings(0,1), ...
    strings(0,1),strings(0,1),repmat(string(parameter),0,1),repmat(string(unit),0,1), ...
    'VariableNames',{'EC','organism','substrate','value','rawValue', ...
                     'commentary','literature','parameter','unit'});
end

function value = localText(record,field)
value = "";
if ~isfield(record,field) || isempty(record.(field)), return; end
v = record.(field);
if ischar(v) || (isstring(v) && isscalar(v)), value = string(v);
else, value = string(jsonencode(v)); end
end

function text = localOrganismText(organism)
if isempty(organism), text = " (all organisms)";
else, text = " in " + string(organism); end
end
