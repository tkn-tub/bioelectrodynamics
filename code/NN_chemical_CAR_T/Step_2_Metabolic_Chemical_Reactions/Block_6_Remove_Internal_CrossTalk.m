%% Protect algorithm species within selected reactions; external edges allowed.
minimalDownstreamOptions;
RegulatoryEvidenceInput = loadMinimalRegulation(model,REGULATION_FILE);
if ~exist('Motifs_NN_Block_5','var')
    b6Saved = loadMinimalCheckpoint('Block_5_no_cell_membrane.mat',model,'Motifs_NN_Block_5');
    Motifs_NN_Block_5 = b6Saved.Motifs_NN_Block_5;
    InputAdapterPaths = b6Saved.InputAdapterPaths;
    clear b6Saved;
end
b6Keep = false(height(Motifs_NN_Block_5),1);
b6Reasons = strings(height(Motifs_NN_Block_5),1);
b6Schema = minimalCircuitSchema(4,ismember('r_CsLeak',Motifs_NN_Block_5.Properties.VariableNames));
b6EC = getReactionECAnnotations(model);
InputAdapterPaths = requireECAnnotatedAdapters(model,InputAdapterPaths,b6EC,BLOCK0_REQUIRE_EC);
[~,b6EligibleIndex] = requireECAnnotatedRows(Motifs_NN_Block_5,b6EC,b6Schema.rxns,InputAdapterPaths);
b6Eligible = false(height(Motifs_NN_Block_5),1); b6Eligible(b6EligibleIndex) = true;
b6BoundedGraph = minimalDirectedGraph(model,isBadMet,false);
for r = 1:numel(b6Schema.rxns)
    b6Eligible = b6Eligible & ismember([Motifs_NN_Block_5.(b6Schema.rxns{r}), ...
        Motifs_NN_Block_5.(b6Schema.dirs{r})],[b6BoundedGraph.rxn,b6BoundedGraph.dir],'rows');
end
b6Supported = false(height(Motifs_NN_Block_5),1);
b6Pulse = tic;
for k = 1:height(Motifs_NN_Block_5)
    if ~b6Eligible(k)
        b6Reasons(k) = "Missing EC, infeasible direction or ineligible adapter"; continue;
    end
    [b6Keep(k),b6Reasons(k)] = checkMinimalCircuit(model,table2struct(Motifs_NN_Block_5(k,:)),4);
    b6Supported(k) = minimalRegulatoryEdges(table2struct(Motifs_NN_Block_5(k,:)),4,RegulatoryEvidenceInput);
    if string(REGULATION_MODE) == "evidence" && ~b6Supported(k)
        b6Keep(k) = false; b6Reasons(k) = "Missing supplied regulatory evidence";
    end
    if toc(b6Pulse) >= 5
        fprintf('Cross-talk filter: %.1f%% | %d/%d\n',100*k/height(Motifs_NN_Block_5),k,height(Motifs_NN_Block_5));
        b6Pulse = tic;
    end
end
Motifs_NN_Block_6 = Motifs_NN_Block_5(b6Keep,:);
Motifs_NN_Block_6.RegulatoryEvidenceComplete = b6Supported(b6Keep);
Motifs_NN_Block_6.KineticsValidated = false(height(Motifs_NN_Block_6),1);
Motifs_NN_Block_6.motifIndexBeforeCrossTalkFilter = find(b6Keep);
results_nn_no_crosstalk = table2struct(Motifs_NN_Block_6);
fullJoinRows = results_nn_no_crosstalk;
nF = height(Motifs_NN_Block_6);
CrossTalk_Filter_Report_Block_6 = table((1:height(Motifs_NN_Block_5)).',b6Keep,b6Reasons, ...
    'VariableNames',{'OriginalIndex','Retained','Reason'});
b6Info = struct('rowsBefore',height(Motifs_NN_Block_5),'rowsAfter',nF);
saveMinimalCircuit('Block_6_no_internal_crosstalk.mat','Motifs_NN_Block_6', ...
    Motifs_NN_Block_6,model,InputAdapterPaths,downstreamOptions,b6Info, ...
    CrossTalk_Filter_Report_Block_6,RegulatoryEvidenceInput);
fprintf('Internal cross-talk filter: %d -> %d motifs.\n',b6Info.rowsBefore,nF);
clear b6Keep b6Reasons b6Pulse b6Info b6Schema b6EC b6EligibleIndex b6Eligible b6BoundedGraph b6Supported downstreamOptions;
