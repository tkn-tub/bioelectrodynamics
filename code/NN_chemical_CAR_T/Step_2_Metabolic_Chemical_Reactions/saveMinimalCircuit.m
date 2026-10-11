function saveMinimalCircuit(filename,tableName,motifs,model,paths,settings,info,diagnostics,evidence)
% One numeric/logical table, one set of annotations, no model or caches.
checkpointSettings = settings;
checkpointSettings.stage = 4;
if strcmp(tableName,'Motifs_Counter_Block_3'), checkpointSettings.stage = 3; end
payload = fiveReactionDifferentiatorCheckpoint(model,motifs,diagnostics,paths,checkpointSettings,info);
payload.TopologyVersion = "five_reaction_regulated_full_v2";
payload.TopologyDescription = "Ten differentiator reactions; a regulated counter with external Cs support; four Hill-NN production/removal pairs.";
payload.CounterCsLeakMapped = ismember('r_CsLeak',motifs.Properties.VariableNames);
schema = minimalCircuitSchema(checkpointSettings.stage,payload.CounterCsLeakMapped);
payload.FullReactionRoles = string(schema.rxns);
payload.ProtectedMetaboliteRoles = string(schema.mets);
if payload.CounterCsLeakMapped
    payload.CounterSupportAssumption = "External positive Cs supply; basal Cs removal is mapped.";
else
    payload.CounterSupportAssumption = "The complete net Cs buffer term k_s*(Cs_baseline-Cs) is supplied externally.";
end
payload.(tableName) = motifs;
payload.Diagnostics = diagnostics;
payload.RegulatoryEvidenceInput = evidence;
if isfield(settings,'regulationMode'), payload.RegulationMode = string(settings.regulationMode); end
saveCompactCheckpoint(filename,payload);
end
