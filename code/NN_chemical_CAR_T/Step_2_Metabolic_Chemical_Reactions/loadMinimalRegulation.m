function evidence = loadMinimalRegulation(model,filename)
% Exact model IDs, signed reaction direction, effect, and nonempty evidence.
% No interactions are inferred from GPRs, EC codes, names, or model.S.
evidence = table(zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1), ...
    'VariableNames',{'MetIndex','RxnIndex','Direction','Effect','Evidence'});
if isempty(filename) || strlength(string(filename)) == 0, return; end
if ~isfile(filename), error('Regulation:MissingFile','Cannot find %s.',filename); end
raw = readtable(filename,'TextType','string');
required = {'RegulatorMetID','TargetRxnID','TargetDirection','Effect','Evidence'};
if ~all(ismember(required,raw.Properties.VariableNames))
    error('Regulation:Columns','Regulation file requires columns: %s.',strjoin(required,', '));
end
if isempty(raw), return; end
[hasMet,m] = ismember(string(raw.RegulatorMetID),string(model.mets(:)));
[hasRxn,r] = ismember(string(raw.TargetRxnID),string(model.rxns(:)));
d = str2double(string(raw.TargetDirection));
effectText = lower(strtrim(string(raw.Effect)));
effect = zeros(height(raw),1);
effect(effectText == "activation") = 1; effect(effectText == "inhibition") = -1;
source = strtrim(string(raw.Evidence));
valid = hasMet & hasRxn & ismember(d,[-1 1]) & effect ~= 0 & ...
    ~ismissing(source) & strlength(source) > 0;
if ~all(valid)
    error('Regulation:InvalidRows', ...
        'Invalid IDs/direction/effect/evidence in regulation rows: %s.',strjoin(string(find(~valid)),', '));
end
evidence = table(m,r,d,effect,source,'VariableNames',evidence.Properties.VariableNames);
end
