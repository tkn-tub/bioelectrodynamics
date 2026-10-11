%% Resume your successful Block 2 run in the current MATLAB workspace.
% Current 5+5 differentiators; no COBRA reinitialization or repeat of Blocks 0-2.
addpath(fileparts(mfilename('fullpath')));
if ~exist('model','var') || ~isfield(model,'S')
    error('MinimalContinue:Model', ...
        'Run after your differentiator master with the same model still loaded.');
end
if ~exist('isBadMet','var') || numel(isBadMet) ~= size(model.S,1)
    error('MinimalContinue:CurrencyMask','Prepare the currency mask for the same model.');
end
continueBlock2 = loadMinimalCheckpoint('Block_2_differentiator.mat',model,'Motifs_Diff_Block_2');
Motifs_Diff_Block_2 = continueBlock2.Motifs_Diff_Block_2;
InputAdapterPaths = continueBlock2.InputAdapterPaths;
clear continueBlock2;
minimalDownstreamOptions;
if ~exist('RUN_BLOCK3','var'), RUN_BLOCK3 = true; end
if ~exist('RUN_BLOCK4','var'), RUN_BLOCK4 = true; end
if ~exist('RUN_BLOCK5','var'), RUN_BLOCK5 = true; end
if ~exist('RUN_BLOCK6','var'), RUN_BLOCK6 = true; end
if ~exist('FILTER_ADAPTER_MEMBRANE','var'), FILTER_ADAPTER_MEMBRANE = true; end
if ~exist('MOTIF_TO_PRINT','var'), MOTIF_TO_PRINT = 1; end
if ~exist('PRINT_SELECTED_MOTIF','var'), PRINT_SELECTED_MOTIF = true; end

fprintf('\nContinuing from %d Block 2 cascades.\n',height(Motifs_Diff_Block_2));
if RUN_BLOCK3
    Block_3_Counter_Mapping;
else
    continueSaved = loadMinimalCheckpoint('Block_3_counter.mat',model,'Motifs_Counter_Block_3');
    Motifs_Counter_Block_3 = continueSaved.Motifs_Counter_Block_3;
    InputAdapterPaths = continueSaved.InputAdapterPaths;
    localCheckBufferMode(continueSaved,MAP_COUNTER_CS_LEAK);
    localCheckAncestry(Motifs_Counter_Block_3,Motifs_Diff_Block_2,'block2_idx',2);
end
if RUN_BLOCK4
    Block_4_Hill_NN_Mapping;
else
    continueSaved = loadMinimalCheckpoint('Block_4_Hill_NN.mat',model,'Motifs_NN_Block_4');
    Motifs_NN_Block_4 = continueSaved.Motifs_NN_Block_4;
    InputAdapterPaths = continueSaved.InputAdapterPaths;
    RegulatoryEvidenceInput = continueSaved.RegulatoryEvidenceInput;
    localCheckBufferMode(continueSaved,MAP_COUNTER_CS_LEAK);
    localCheckAncestry(Motifs_NN_Block_4,Motifs_Counter_Block_3,'block3_idx',3);
end
if RUN_BLOCK5
    Block_5_Remove_Cell_Membrane_Reactions;
else
    continueSaved = loadMinimalCheckpoint('Block_5_no_cell_membrane.mat',model,'Motifs_NN_Block_5');
    Motifs_NN_Block_5 = continueSaved.Motifs_NN_Block_5;
    InputAdapterPaths = continueSaved.InputAdapterPaths;
    RegulatoryEvidenceInput = continueSaved.RegulatoryEvidenceInput;
    localCheckBufferMode(continueSaved,MAP_COUNTER_CS_LEAK);
    localCheckFilterParent(Motifs_NN_Block_5,Motifs_NN_Block_4,'motifIndexBeforeMembraneFilter');
end
if RUN_BLOCK6
    Block_6_Remove_Internal_CrossTalk;
else
    continueSaved = loadMinimalCheckpoint('Block_6_no_internal_crosstalk.mat',model,'Motifs_NN_Block_6');
    Motifs_NN_Block_6 = continueSaved.Motifs_NN_Block_6;
    InputAdapterPaths = continueSaved.InputAdapterPaths;
    RegulatoryEvidenceInput = continueSaved.RegulatoryEvidenceInput;
    localCheckBufferMode(continueSaved,MAP_COUNTER_CS_LEAK);
    localCheckFilterParent(Motifs_NN_Block_6,Motifs_NN_Block_5,'motifIndexBeforeCrossTalkFilter');
    if isfield(continueSaved,'RegulationMode') && string(continueSaved.RegulationMode) ~= string(REGULATION_MODE)
        error('MinimalContinue:CachedPolicy','Set RUN_BLOCK6=true to apply the new regulation mode.');
    end
end
results_nn_no_crosstalk = table2struct(Motifs_NN_Block_6);
fullJoinRows = results_nn_no_crosstalk;
nF = height(Motifs_NN_Block_6);
Regulation_Requirements = writeMinimalRegulationRequirements(Motifs_NN_Block_6,model,RegulatoryEvidenceInput);
fprintf('\nCURRENT FIVE-REACTION PIPELINE COMPLETE\n');
fprintf('Block 2 cascades             : %d\n',height(Motifs_Diff_Block_2));
fprintf('Block 3 counter joins        : %d\n',height(Motifs_Counter_Block_3));
fprintf('Block 4 Hill NN joins        : %d\n',height(Motifs_NN_Block_4));
fprintf('After membrane filter       : %d\n',height(Motifs_NN_Block_5));
fprintf('After internal cross-talk   : %d\n',nF);
fprintf('Supplied individual regulatory links complete: %d\n',sum(Motifs_NN_Block_6.RegulatoryEvidenceComplete));
fprintf('Kinetically validated       : 0 (not tested by this search)\n');
fprintf('Unverified unique regulatory requirements: %d\n',sum(~Regulation_Requirements.Supported));
if PRINT_SELECTED_MOTIF
    PRINT_RESULTS_FILE = 'Block_6_no_internal_crosstalk.mat';
    Block_7_Print_Selected_Motif;
end
clear continueSaved downstreamOptions;

function localCheckBufferMode(saved,mapLeak)
if logical(saved.CounterCsLeakMapped) ~= logical(mapLeak)
    error('MinimalContinue:CachedBuffer', ...
        'Cached Cs support differs from MAP_COUNTER_CS_LEAK. Regenerate Blocks 3-6 or restore that setting.');
end
end

function localCheckAncestry(rows,parent,indexField,stage)
schema = minimalCircuitSchema(stage,ismember('r_CsLeak',parent.Properties.VariableNames));
fields = [schema.mets,schema.rxns,schema.dirs,{'adapterPathIndex'}];
localCompareParent(rows,parent,indexField,fields);
end

function localCheckFilterParent(rows,parent,indexField)
schema = minimalCircuitSchema(4,ismember('r_CsLeak',parent.Properties.VariableNames));
localCompareParent(rows,parent,indexField,[schema.mets,schema.rxns,schema.dirs,{'adapterPathIndex'}]);
end

function localCompareParent(rows,parent,indexField,fields)
if ~ismember(indexField,rows.Properties.VariableNames)
    error('MinimalContinue:CachedParent','Cached rows lack their saved parent index %s.',indexField);
end
idx = rows.(indexField);
if any(~isfinite(idx) | idx < 1 | idx > height(parent) | idx ~= floor(idx)) || ...
        ~isequal(table2array(rows(:,fields)),table2array(parent(idx,fields)))
    error('MinimalContinue:CachedParent', ...
        'Cached downstream rows do not match the loaded upstream results. Regenerate the downstream blocks.');
end
end
