function saved = loadMinimalCheckpoint(filename,model,tableName)
% Validate current polarity schemas and the original model's ID ordering.
saved = load(filename);
if strcmp(tableName,'Motifs_Diff_Block_2')
    expected = "five_reaction_differentiators_v2"; stage = 2;
elseif strcmp(tableName,'Motifs_Counter_Block_3')
    expected = "five_reaction_regulated_full_v2"; stage = 3;
elseif ismember(string(tableName),["Motifs_NN_Block_4","Motifs_NN_Block_5","Motifs_NN_Block_6"])
    expected = "five_reaction_regulated_full_v2"; stage = 4;
else
    error('MinimalLoad:TableName','Unsupported checkpoint table %s.',tableName);
end
if ~isfield(saved,'TopologyVersion') || ~isequal(string(saved.TopologyVersion),expected)
    error('MinimalLoad:OldTopology', ...
        '%s has an incompatible topology. Expected %s; regenerate old downstream checkpoints.',filename,expected);
end
if ~isfield(saved,tableName) || ~istable(saved.(tableName))
    error('MinimalLoad:MissingTable','%s lacks %s.',filename,tableName);
end
if ~all(isfield(saved,{'EC_ModelMetIDs','EC_ModelRxnIDs','InputAdapterPaths'})) || ...
        ~isequal(string(saved.EC_ModelMetIDs(:)),string(model.mets(:))) || ...
        ~isequal(string(saved.EC_ModelRxnIDs(:)),string(model.rxns(:)))
    error('MinimalLoad:ModelMismatch','%s uses another model/index ordering.',filename);
end
rows = saved.(tableName);
mapLeak = ismember('r_CsLeak',rows.Properties.VariableNames);
schema = minimalCircuitSchema(stage,mapLeak);
required = [schema.mets,schema.rxns,schema.dirs,{'adapterPathIndex','m_d','m_d1','m_d2','m_X2'}];
if ~all(ismember(required,rows.Properties.VariableNames)) || ismember('r_CsSupply',rows.Properties.VariableNames)
    error('MinimalLoad:Schema','%s does not contain the current reaction roles.',filename);
end
if any(rows.m_d ~= rows.m_Y1 | rows.m_d1 ~= rows.m_Y1 | rows.m_d2 ~= rows.m_Y2 | rows.m_X2 ~= rows.m_Y1)
    error('MinimalLoad:Aliases','Saved aliases must identify m_d=m_d1=Y1 and m_d2=Y2.');
end
if stage >= 3 && (~isfield(saved,'CounterCsLeakMapped') || ~isscalar(saved.CounterCsLeakMapped) || ...
        ~isequal(logical(saved.CounterCsLeakMapped),mapLeak))
    error('MinimalLoad:BufferMode','Saved Cs-support metadata disagree with the reaction columns.');
end
indices = table2array(rows(:,schema.mets)); reactions = table2array(rows(:,schema.rxns));
directions = table2array(rows(:,schema.dirs));
if ~isnumeric(indices) || ~isnumeric(reactions) || ~isnumeric(directions) || ...
        any(~isfinite(indices) | indices < 1 | indices > size(model.S,1) | indices ~= floor(indices),'all') || ...
        any(~isfinite(reactions) | reactions < 1 | reactions > size(model.S,2) | reactions ~= floor(reactions),'all') || ...
        any(~ismember(directions,[-1 1]),'all')
    error('MinimalLoad:Indices','Saved indices or directions are invalid for this model.');
end
end
