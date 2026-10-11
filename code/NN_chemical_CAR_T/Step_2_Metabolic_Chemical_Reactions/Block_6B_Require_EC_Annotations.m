%% Optional compatibility gate for new full motifs; the searches enforce ECs.
if ~exist('BLOCK0_REQUIRE_EC','var'), BLOCK0_REQUIRE_EC = false; end
b6bMask = getReactionECAnnotations(model);
InputAdapterPaths = requireECAnnotatedAdapters(model,InputAdapterPaths,b6bMask,BLOCK0_REQUIRE_EC);
b6bSchema = minimalCircuitSchema(4,ismember('r_CsLeak',Motifs_NN_Block_6.Properties.VariableNames));
[Motifs_NN_Block_6,b6bIndex] = requireECAnnotatedRows(Motifs_NN_Block_6,b6bMask,b6bSchema.rxns,InputAdapterPaths);
Motifs_NN_Block_6.motifIndexBeforeECFilter = b6bIndex;
results_nn_no_crosstalk = table2struct(Motifs_NN_Block_6);
fullJoinRows = results_nn_no_crosstalk;
nF = height(Motifs_NN_Block_6);
saveMinimalCircuit('Block_6_EC_annotated.mat','Motifs_NN_Block_6', ...
    Motifs_NN_Block_6,model,InputAdapterPaths,struct(),struct('rowsAfter',nF), ...
    table(b6bIndex),RegulatoryEvidenceInput);
clear b6bMask b6bSchema b6bIndex;
