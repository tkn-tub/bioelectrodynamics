function [Km, Kcat] = Brenda_Kinetics(EC, organism, email, password)
% BRENDA_KINETICS Retrieve measured kinetic records by EC number.
%
% Requires a registered BRENDA account and a supported Python installation
% with zeep available to MATLAB's pyenv interpreter. No generated WSDL class.
%
% Example, with credentials supplied locally:
%   [Km,Kcat] = Brenda_Kinetics('2.7.4.14','Homo sapiens',email,password);
%   disp(Km); disp(Kcat);
%
% Or set BRENDA_EMAIL and BRENDA_PASSWORD in your local environment:
%   [Km,Kcat] = Brenda_Kinetics('2.7.4.14','Homo sapiens');
%
% Use organism='' to query all organisms explicitly. Empty human results
% are retained; the function never silently substitutes another organism.
% Results include substrate, commentary and references. No value is selected
% automatically. value/valueMaximum are in mM for Km and 1/s for Kcat.
% A missing maximum stays NaN; it does not imply a second measured value.
%
% To query the attached motif once per unique EC (after auditing annotations):
%   txt = fileread('Full_Motif_1(1).txt');
%   ec = regexp(txt,'(?m)^EC:\s*(\d+\.\d+\.\d+\.\d+)\s*$','tokens');
%   ec = unique(cellfun(@(x)x{1},ec,'UniformOutput',false),'stable');
%   Km_all = table(); Kcat_all = table();
%   for j = 1:numel(ec)
%       [a,b] = Brenda_Kinetics(ec{j},'Homo sapiens',email,password);
%       Km_all = [Km_all;a]; Kcat_all = [Kcat_all;b];
%   end
%   writetable(Km_all,'Brenda_Km_records.csv');
%   writetable(Kcat_all,'Brenda_kcat_records.csv');
%
% The attached motif is the older three-reaction/Hill-NN topology. Querying
% its EC codes does not validate the newer five-reaction/affine topology.
% Audit FTHFL: its named reaction corresponds to EC 6.3.4.3, whereas the file
% lists 6.3.2.17. ASPGLUm is a transport reaction; its listed EC 4.1.1.15
% denotes glutamate decarboxylase, not that shuttle. No silent corrections.
%
% Sources (checked 2026-10-08):
% https://www.brenda-enzymes.org/soap.php
% https://www.brenda-enzymes.org/datafields.php
% https://www.mathworks.com/help/matlab/call-python-libraries.html
% No authenticated BRENDA query or native MATLAB execution was available
% during preparation. Review assay conditions before assigning parameters.

if nargin < 2
    organism = 'Homo sapiens';
end
if nargin < 3
    email = getenv('BRENDA_EMAIL');
end
if nargin < 4
    password = getenv('BRENDA_PASSWORD');
end
EC = char(EC);
organism = char(organism);
email = char(email);
password = char(password);
assert(isrow(EC) && ~isempty(regexp(EC,'^\d+\.\d+\.\d+\.\d+$','once')), ...
    'Use a complete EC number, for example 2.7.4.14.');
assert(~isempty(email) && ~isempty(password), ...
    'Supply your registered BRENDA email/password locally or through environment variables.');

% Reuse the service description for successive EC queries.
persistent client
if isempty(client)
    zeep = py.importlib.import_module('zeep');
    transportModule = py.importlib.import_module('zeep.transports');
    transport = transportModule.Transport(pyargs('timeout',60,'operation_timeout',120));
    settings = zeep.Settings(pyargs('strict',false));
    client = zeep.Client('https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl', ...
        pyargs('settings',settings,'transport',transport));
end

% BRENDA requires the SHA-256 digest of the UTF-8 password.
passwordText = py.str(password);
hashObject = py.hashlib.sha256(passwordText.encode('utf-8'));
passwordHash = char(hashObject.hexdigest());
clear password passwordText hashObject

% Each field is a separate positional argument, as in BRENDA's Zeep example.
% BRENDA permits no more than one API request per second. Wait before every
% call, including the first, to separate calls across successive function uses.
pause(1.1);
rawKm = client.service.getKmValue(email,passwordHash, ...
    ['ecNumber*' EC],'kmValue*','kmValueMaximum*','substrate*', ...
    'commentary*',['organism*' organism],'ligandStructureId*','literature*');
pause(1.1);
rawKcat = client.service.getTurnoverNumber(email,passwordHash, ...
    ['ecNumber*' EC],'turnoverNumber*','turnoverNumberMaximum*','substrate*', ...
    'commentary*',['organism*' organism],'ligandStructureId*','literature*');
clear passwordHash

Km = recordsToTable(rawKm,'kmValue','kmValueMaximum','Km','mM');
Kcat = recordsToTable(rawKcat,'turnoverNumber','turnoverNumberMaximum','kcat','1/s');
end

function T = recordsToTable(raw, valueField, maximumField, parameter, unit)
% Serialize SOAP objects without losing commentary or literature references.
helper = py.importlib.import_module('zeep.helpers');
native = helper.serialize_object(raw);
records = jsondecode(char(py.json.dumps(native)));
n = numel(records);
T = table(strings(n,1),strings(n,1),strings(n,1),nan(n,1),nan(n,1), ...
    strings(n,1),strings(n,1),strings(n,1),strings(n,1),strings(n,1), ...
    repmat(string(parameter),n,1),repmat(string(unit),n,1), ...
    'VariableNames',{'EC','organism','substrate','value','valueMaximum', ...
    'rawValue','rawValueMaximum','commentary','literature','ligandStructureId', ...
    'parameter','unit'});
for i = 1:n
    if iscell(records)
        r = records{i};
    else
        r = records(i);
    end
    assert(isstruct(r),'Unexpected SOAP result: expected kinetic record objects.');
    T.EC(i) = textField(r,'ecNumber');
    T.organism(i) = textField(r,'organism');
    T.substrate(i) = textField(r,'substrate');
    T.commentary(i) = textField(r,'commentary');
    T.literature(i) = textField(r,'literature');
    T.ligandStructureId(i) = textField(r,'ligandStructureId');
    T.rawValue(i) = textField(r,valueField);
    T.rawValueMaximum(i) = textField(r,maximumField);
    T.value(i) = str2double(T.rawValue(i));
    T.valueMaximum(i) = str2double(T.rawValueMaximum(i));
end
% Preserve original strings but exclude nonnumeric/nonpositive sentinel values
% from the numerical columns. Ranges and missing values are never averaged.
T.value(~isfinite(T.value) | T.value <= 0) = NaN;
T.valueMaximum(~isfinite(T.valueMaximum) | T.valueMaximum <= 0) = NaN;
end

function text = textField(record, field)
text = "";
if ~isfield(record,field) || isempty(record.(field))
    return
end
value = record.(field);
if ischar(value) || (isstring(value) && isscalar(value))
    text = string(value);
else
    text = string(jsonencode(value)); % reference arrays remain identifiable
end
end
