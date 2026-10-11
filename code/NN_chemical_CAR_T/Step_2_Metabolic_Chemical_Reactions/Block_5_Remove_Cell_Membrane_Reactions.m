%% Current polarity/Hill schema; include adapter reactions by default.
if ~exist('FILTER_ADAPTER_MEMBRANE','var'), FILTER_ADAPTER_MEMBRANE = true; end
if ~exist('Motifs_NN_Block_4','var')
    b5Saved = loadMinimalCheckpoint('Block_4_Hill_NN.mat',model,'Motifs_NN_Block_4');
    Motifs_NN_Block_4 = b5Saved.Motifs_NN_Block_4;
    InputAdapterPaths = b5Saved.InputAdapterPaths;
    RegulatoryEvidenceInput = b5Saved.RegulatoryEvidenceInput;
    clear b5Saved;
end
[b5Membrane,b5CompartmentLabels] = classifyMinimalMembrane(model);
b5Schema = minimalCircuitSchema(4,ismember('r_CsLeak',Motifs_NN_Block_4.Properties.VariableNames));
b5Indices = table2array(Motifs_NN_Block_4(:,b5Schema.rxns));
if any(b5Indices < 1 | b5Indices > numel(b5Membrane) | b5Indices ~= floor(b5Indices),'all')
    error('Block5:ModelIndex','Invalid selected reaction index.');
end
b5Membership = reshape(b5Membrane(b5Indices(:)),size(b5Indices));
b5Keep = ~any(b5Membership,2);
if FILTER_ADAPTER_MEMBRANE
    for p = 1:numel(InputAdapterPaths)
        if any(b5Membrane(InputAdapterPaths(p).reactionIndices))
            b5Keep(Motifs_NN_Block_4.adapterPathIndex == InputAdapterPaths(p).pathIndex) = false;
        end
    end
end
Motifs_NN_Block_5 = Motifs_NN_Block_4(b5Keep,:);
Motifs_NN_Block_5.motifIndexBeforeMembraneFilter = find(b5Keep);
results_nn_no_membrane = table2struct(Motifs_NN_Block_5);
fullJoinRows = results_nn_no_membrane;
Membrane_Filter_Report_Block_5 = table((1:height(Motifs_NN_Block_4)).',b5Keep, ...
    'VariableNames',{'OriginalIndex','Retained'});
b5Info = struct('rowsBefore',height(Motifs_NN_Block_4),'rowsAfter',height(Motifs_NN_Block_5));
b5Settings = struct('filterAdapterMembrane',FILTER_ADAPTER_MEMBRANE);
saveMinimalCircuit('Block_5_no_cell_membrane.mat','Motifs_NN_Block_5', ...
    Motifs_NN_Block_5,model,InputAdapterPaths,b5Settings,b5Info, ...
    Membrane_Filter_Report_Block_5,RegulatoryEvidenceInput);
fprintf('Membrane filter: %d -> %d motifs. Internal organelle membranes retained.\n', ...
    b5Info.rowsBefore,b5Info.rowsAfter);
clear b5Membrane b5Membership b5CompartmentLabels b5Schema b5Indices b5Keep b5Info b5Settings;
