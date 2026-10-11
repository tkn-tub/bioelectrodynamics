%% Regulated counter: Cp, Cg, Cr. Cs support is external.
% Optional MAP_COUNTER_CS_LEAK=true maps an additional basal Cs removal.
% Only m_d2=Y2 inhibits Cp and Cr. These are NONCONSUMING links.
minimalDownstreamOptions;
if ~exist('Motifs_Diff_Block_2','var')
    b3Saved = loadMinimalCheckpoint('Block_2_differentiator.mat',model,'Motifs_Diff_Block_2');
    Motifs_Diff_Block_2 = b3Saved.Motifs_Diff_Block_2;
    InputAdapterPaths = b3Saved.InputAdapterPaths;
    clear b3Saved;
end
[Motifs_Counter_Block_3,Counter_Search_Diagnostics_Block_3,InputAdapterPaths, ...
    b3Info,RegulatoryEvidenceInput] = extendMinimalCircuit( ...
    model,isBadMet,InputAdapterPaths,Motifs_Diff_Block_2,3,downstreamOptions);
results_counter = table2struct(Motifs_Counter_Block_3);
saveMinimalCircuit('Block_3_counter.mat','Motifs_Counter_Block_3', ...
    Motifs_Counter_Block_3,model,InputAdapterPaths,downstreamOptions,b3Info, ...
    Counter_Search_Diagnostics_Block_3,RegulatoryEvidenceInput);
clear b3Info downstreamOptions;
