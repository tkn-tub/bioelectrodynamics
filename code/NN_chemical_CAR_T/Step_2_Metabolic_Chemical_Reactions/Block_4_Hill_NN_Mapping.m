%% Hill-production NN: C -> Ni -> {Na,Nb} -> output through regulation.
% Four production/removal pairs. C,Ni,Na,Nb are NOT consumed by the regulated
% downstream production reactions. No fictitious material-flow arrows.
minimalDownstreamOptions;
downstreamOptions.maxHits = MAX_HITS_NN;
downstreamOptions.poolCap = MAX_NN_POOL;
if ~exist('Motifs_Counter_Block_3','var')
    b4Saved = loadMinimalCheckpoint('Block_3_counter.mat',model,'Motifs_Counter_Block_3');
    Motifs_Counter_Block_3 = b4Saved.Motifs_Counter_Block_3;
    InputAdapterPaths = b4Saved.InputAdapterPaths;
    clear b4Saved;
end
[Motifs_NN_Block_4,NN_Search_Diagnostics_Block_4,InputAdapterPaths, ...
    b4Info,RegulatoryEvidenceInput] = extendMinimalCircuit( ...
    model,isBadMet,InputAdapterPaths,Motifs_Counter_Block_3,4,downstreamOptions);
results_nn = table2struct(Motifs_NN_Block_4);
saveMinimalCircuit('Block_4_Hill_NN.mat','Motifs_NN_Block_4', ...
    Motifs_NN_Block_4,model,InputAdapterPaths,downstreamOptions,b4Info, ...
    NN_Search_Diagnostics_Block_4,RegulatoryEvidenceInput);
clear b4Info downstreamOptions;
