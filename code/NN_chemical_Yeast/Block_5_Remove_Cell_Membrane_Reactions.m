%% Block 5: remove complete motifs containing cell-membrane reactions
%
% This filter is applied to the full Block 4 differentiator + counter + NN
% table. A complete motif is rejected when ANY canonical reaction role:
%
%   1) explicitly involves the "cell envelope" compartment, OR
%   2) contains both a cytoplasmic and an extracellular metabolite.
%
% Internal membranes (mitochondrial, ER, Golgi, vacuolar, and so on) are
% intentionally retained. Thus, this is an outer-cell-boundary filter; it is
% not a filter for every multi-compartment reaction in the model.
%
% The canonical reaction roles checked are:
%
%   First differentiator:
%       r_S1, r_Sg1, r_Ss1, r_F1, r_Fg1, r_Fs1, r_I1
%   Second differentiator:
%       r_S2, r_Sg2, r_Ss2, r_F2, r_Fg2, r_Fs2, r_I2
%   Counter:
%       r_C, r_CI1, r_CI2
%   Single-neuron NN:
%       r_omega, r_b
%
% Compatibility aliases such as r_Cp/r_Cr/r_Cg and r_Nomega/r_Nb are not
% checked separately because they point to the same canonical reactions.
%
% Expected master-file sequence:
%   load Block_4_single_neuron_NN.mat
%   Block_5_Remove_Cell_Membrane_Reactions
%
% Required workspace variables:
%   model, including model.S, model.compNames, model.metComps and model.rxns
% Block 4 must provide either:
%   Motifs_NN_Block_4      full table (preferred), or
%   results_nn             full struct array
%
% Saved output:
%   Block_5_no_cell_membrane.mat
%       Motifs_NN_Block_5                 filtered full table
%       results_nn_no_membrane            filtered full struct array
%       Membrane_Filter_Report_Block_5    rejected rows and reasons
%       Cell_Membrane_Reactions_Block_5   classified model reactions
%       indexBeforeMembraneFilter         retained original row indices
%       keepMotif                         logical mask over Block 4 rows
%       isCellMembraneRxn                 logical mask over model reactions
%       membraneReactionReason            classification details by reaction

%% Display controls
if ~exist('SHOW_MEMBRANE_FILTER_PROGRESS','var')
    SHOW_MEMBRANE_FILTER_PROGRESS = true;
end

if ~exist('MEMBRANE_PROGRESS_STEP_PERCENT','var')
    MEMBRANE_PROGRESS_STEP_PERCENT = 10;
end
MEMBRANE_PROGRESS_STEP_PERCENT = max(1, ...
    min(100,round(MEMBRANE_PROGRESS_STEP_PERCENT)));

% Rejected motifs can be numerous. Detailed printing is therefore disabled
% by default; the complete report is saved in the MAT file regardless.
if ~exist('PRINT_REJECTED_MOTIFS','var')
    PRINT_REJECTED_MOTIFS = false;
end

if ~exist('MAX_REJECTED_TO_PRINT','var')
    MAX_REJECTED_TO_PRINT = 20;
end
MAX_REJECTED_TO_PRINT = max(0,floor(MAX_REJECTED_TO_PRINT));

if ~exist('PRINT_RETAINED_INDICES','var')
    PRINT_RETAINED_INDICES = false;
end

%% Validate the model
if ~exist('model','var')
    error('Block5:MissingModel', ...
        'Load the COBRA model before running Block 5.');
end

requiredModelFields = {'S','compNames','metComps','rxns'};
missingModelFields = requiredModelFields( ...
    ~isfield(model,requiredModelFields));

if ~isempty(missingModelFields)
    error('Block5:IncompleteModel', ...
        'model is missing required field(s): %s.', ...
        strjoin(missingModelFields,', '));
end

[nMets,nRxns] = size(model.S);

rawCompNames = model.compNames;
if ischar(rawCompNames)
    rawCompNames = cellstr(rawCompNames);
end
compNameLabels = string(rawCompNames(:));
nCompartments = numel(compNameLabels);

if numel(model.metComps) ~= nMets
    error('Block5:InvalidMetComps', ...
        ['model.metComps has %d entries, but model.S has %d ', ...
         'metabolite rows.'],numel(model.metComps),nMets);
end

metComps = double(model.metComps(:));
if any(~isfinite(metComps) | metComps < 1 | ...
       metComps > nCompartments | metComps ~= floor(metComps))
    error('Block5:InvalidMetCompIndex', ...
        ['model.metComps contains an invalid compartment index. ', ...
         'Every entry must be an integer from 1 to %d.'],nCompartments);
end

rawRxnIds = model.rxns;
if ischar(rawRxnIds)
    rawRxnIds = cellstr(rawRxnIds);
end
rxnIds = string(rawRxnIds(:));

if numel(rxnIds) ~= nRxns
    error('Block5:InvalidReactionIds', ...
        ['model.rxns has %d entries, but model.S has %d ', ...
         'reaction columns.'],numel(rxnIds),nRxns);
end

%% Read the full Block 4 table without changing it
if exist('Motifs_NN_Block_4','var') && istable(Motifs_NN_Block_4)
    sourceMotifs = Motifs_NN_Block_4;
elseif exist('results_nn','var') && isstruct(results_nn)
    sourceMotifs = struct2table(results_nn(:));
else
    error('Block5:MissingBlock4', ...
        ['Load Block_4_single_neuron_NN.mat before running Block 5. ', ...
         'Expected Motifs_NN_Block_4 or results_nn.']);
end

nMotifs = height(sourceMotifs);
if nMotifs == 0
    error('Block5:EmptyBlock4', ...
        'The full Block 4 table contains no motif rows to filter.');
end

rxnFields = { ...
    'r_S1','r_Sg1','r_Ss1','r_F1','r_Fg1','r_Fs1','r_I1', ...
    'r_S2','r_Sg2','r_Ss2','r_F2','r_Fg2','r_Fs2','r_I2', ...
    'r_C','r_CI1','r_CI2', ...
    'r_omega','r_b'};

sourceVariables = sourceMotifs.Properties.VariableNames;
missingRxnFields = rxnFields(~ismember(rxnFields,sourceVariables));
if ~isempty(missingRxnFields)
    error('Block5:IncompleteBlock4', ...
        ['The full Block 4 table is missing canonical reaction ', ...
         'field(s): %s.'],strjoin(missingRxnFields,', '));
end

%% Locate the outer-cell compartments
compNames = lower(strtrim(compNameLabels));

% "cytosol" and "extracellular space" are accepted as common synonyms.
idxCytoplasm = find(ismember(compNames,["cytoplasm","cytosol"]));
idxEnvelope = find(ismember(compNames, ...
    ["cell envelope","cell membrane","plasma membrane"]));
idxExtra = find(ismember(compNames, ...
    ["extracellular","extracellular space"]));

fprintf('\n============================================================\n');
fprintf('BLOCK 5: OUTER CELL-MEMBRANE FILTER\n');
fprintf('============================================================\n');
fprintf('  cytoplasm/cytosol indices : ');
fprintf('%d ',idxCytoplasm);
fprintf('\n');
fprintf('  cell-envelope indices     : ');
fprintf('%d ',idxEnvelope);
fprintf('\n');
fprintf('  extracellular indices     : ');
fprintf('%d ',idxExtra);
fprintf('\n');

if isempty(idxCytoplasm)
    warning('Block5:MissingCytoplasm', ...
        'Could not identify a cytoplasm or cytosol compartment.');
end
if isempty(idxEnvelope)
    warning('Block5:MissingCellEnvelope', ...
        ['Could not identify a cell envelope compartment. The filter ', ...
         'can still detect direct cytoplasm-to-extracellular reactions.']);
end
if isempty(idxExtra)
    warning('Block5:MissingExtracellular', ...
        'Could not identify an extracellular compartment.');
end

%% Classify every model reaction once
% The previous motif-by-motif implementation repeatedly called find() on the
% same reaction. Here each model reaction is inspected only once. Motif rows
% are then filtered using their already-stored reaction indices.
isCellMembraneRxn = false(nRxns,1);
reactionCompartmentText = strings(nRxns,1);
membraneReactionReason = strings(nRxns,1);

if isfield(model,'rxnNames')
    rawRxnNames = model.rxnNames;
    if ischar(rawRxnNames)
        rawRxnNames = cellstr(rawRxnNames);
    end
    if numel(rawRxnNames) == nRxns
        rxnNames = string(rawRxnNames(:));
    else
        rxnNames = strings(nRxns,1);
    end
else
    rxnNames = strings(nRxns,1);
end

if SHOW_MEMBRANE_FILTER_PROGRESS
    fprintf('\nClassifying model reactions:   0%% (0/%d)',nRxns);
    nextProgress = MEMBRANE_PROGRESS_STEP_PERCENT;
end

for rIdx = 1:nRxns
    metIdx = find(model.S(:,rIdx) ~= 0);

    if ~isempty(metIdx)
        reactionCompIdx = unique(metComps(metIdx),'stable');
        reactionCompartmentText(rIdx) = strjoin( ...
            compNameLabels(reactionCompIdx)," <-> ");

        usesCellEnvelope = any(ismember(reactionCompIdx,idxEnvelope));
        hasCytoplasm = any(ismember(reactionCompIdx,idxCytoplasm));
        hasExtracellular = any(ismember(reactionCompIdx,idxExtra));
        crossesCellBoundary = hasCytoplasm && hasExtracellular;

        isCellMembraneRxn(rIdx) = ...
            usesCellEnvelope || crossesCellBoundary;

        if isCellMembraneRxn(rIdx)
            reasons = strings(0,1);
            if usesCellEnvelope
                reasons(end+1,1) = "uses cell envelope"; %#ok<SAGROW>
            end
            if crossesCellBoundary
                reasons(end+1,1) = ...
                    "connects cytoplasm and extracellular"; %#ok<SAGROW>
            end
            membraneReactionReason(rIdx) = strjoin(reasons," + ");
        end
    end

    if SHOW_MEMBRANE_FILTER_PROGRESS
        percent = floor(100*rIdx/max(nRxns,1));
        if percent >= nextProgress || rIdx == nRxns
            fprintf('\rClassifying model reactions: %3d%% (%d/%d)', ...
                percent,rIdx,nRxns);
            while nextProgress <= percent
                nextProgress = nextProgress + ...
                    MEMBRANE_PROGRESS_STEP_PERCENT;
            end
            drawnow;
        end
    end
end

if SHOW_MEMBRANE_FILTER_PROGRESS
    fprintf('\n');
end

cellMembraneRxnIdx = find(isCellMembraneRxn);
Cell_Membrane_Reactions_Block_5 = table( ...
    cellMembraneRxnIdx,rxnIds(cellMembraneRxnIdx), ...
    rxnNames(cellMembraneRxnIdx), ...
    reactionCompartmentText(cellMembraneRxnIdx), ...
    membraneReactionReason(cellMembraneRxnIdx), ...
    'VariableNames',{ ...
        'reactionIndex','reactionId','reactionName', ...
        'compartments','classificationReason'});

%% Assemble the canonical reaction-index matrix
nRoles = numel(rxnFields);
motifRxnIdx = zeros(nMotifs,nRoles);

for f = 1:nRoles
    values = sourceMotifs.(rxnFields{f});
    if ~isnumeric(values) || ~isvector(values) || numel(values) ~= nMotifs
        error('Block5:InvalidReactionColumn', ...
            ['Column %s must contain one numeric reaction index per ', ...
             'motif row.'],rxnFields{f});
    end
    motifRxnIdx(:,f) = double(values(:));
end

invalidEntries = ~isfinite(motifRxnIdx) | motifRxnIdx < 1 | ...
    motifRxnIdx > nRxns | motifRxnIdx ~= floor(motifRxnIdx);

if any(invalidEntries(:))
    [badRow,badField] = find(invalidEntries,1,'first');
    error('Block5:InvalidReactionIndex', ...
        ['Block 4 row %d has invalid reaction index %g in %s. ', ...
         'Filtering stopped so an incomplete motif is not retained.'], ...
        badRow,motifRxnIdx(badRow,badField),rxnFields{badField});
end

%% Reject a complete motif when any of its 19 reactions is classified
badRoleMask = reshape(isCellMembraneRxn(motifRxnIdx(:)), ...
    nMotifs,nRoles);
keepMotif = ~any(badRoleMask,2);

indexBeforeMembraneFilter = find(keepMotif);
removedIdx = find(~keepMotif);

%% Build a complete rejection report
membraneReason = strings(nMotifs,1);

for k = 1:numel(removedIdx)
    motifIdx = removedIdx(k);
    badFields = find(badRoleMask(motifIdx,:));
    details = strings(numel(badFields),1);

    for q = 1:numel(badFields)
        fieldIdx = badFields(q);
        rIdx = motifRxnIdx(motifIdx,fieldIdx);
        details(q) = sprintf( ...
            '%s: %s | %s | compartments: %s | test: %s', ...
            rxnFields{fieldIdx},rxnIds(rIdx),rxnNames(rIdx), ...
            reactionCompartmentText(rIdx), ...
            membraneReactionReason(rIdx));
    end

    membraneReason(motifIdx) = strjoin(details,"; ");
end

nRejectedReactions = sum(badRoleMask(removedIdx,:),2);
Membrane_Filter_Report_Block_5 = table( ...
    removedIdx,nRejectedReactions,membraneReason(removedIdx), ...
    'VariableNames',{ ...
        'motifIndexBeforeMembraneFilter', ...
        'nCellMembraneReactionRoles', ...
        'reason'});

%% Retain the full rows, their order, and their original indices
Motifs_NN_Block_5 = sourceMotifs(keepMotif,:);
Motifs_NN_Block_5.motifIndexBeforeMembraneFilter = ...
    indexBeforeMembraneFilter;

results_nn_no_membrane = table2struct(Motifs_NN_Block_5);
results_nn_no_membrane = results_nn_no_membrane(:);

% Compatibility aliases for older downstream printing/filtering sections.
% These aliases exist in the workspace but are not duplicated in the MAT file.
fullJoinRows = results_nn_no_membrane;
nF = numel(fullJoinRows);

%% Save the filtered full connection table
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end
outputFile = fullfile(scriptFolder,'Block_5_no_cell_membrane.mat');

save(outputFile, ...
    'Motifs_NN_Block_5','results_nn_no_membrane', ...
    'Membrane_Filter_Report_Block_5', ...
    'Cell_Membrane_Reactions_Block_5', ...
    'indexBeforeMembraneFilter','keepMotif', ...
    'isCellMembraneRxn','membraneReactionReason','rxnFields', ...
    '-v7.3');

%% Summary
fprintf('\n############################################################\n');
fprintf('CELL-MEMBRANE FILTER SUMMARY\n');
fprintf('############################################################\n');
fprintf('Model reactions classified as cell-membrane : %d / %d\n', ...
    sum(isCellMembraneRxn),nRxns);
fprintf('Canonical reaction roles checked per motif  : %d\n',nRoles);
fprintf('Motifs before filter                        : %d\n',nMotifs);
fprintf('Motifs retained                             : %d\n', ...
    sum(keepMotif));
fprintf('Motifs removed                              : %d\n', ...
    sum(~keepMotif));
fprintf('Filtered full-table columns                 : %d\n', ...
    width(Motifs_NN_Block_5));
fprintf('############################################################\n');

if PRINT_RETAINED_INDICES
    fprintf('\nOriginal Block 4 row indices retained:\n');
    fprintf('%d ',indexBeforeMembraneFilter);
    fprintf('\n');
end

if PRINT_REJECTED_MOTIFS
    fprintf('\n============================================================\n');
    fprintf('REJECTED MOTIFS AND CELL-MEMBRANE REACTIONS\n');
    fprintf('============================================================\n');

    if isempty(removedIdx)
        fprintf('No motif contained an outer cell-membrane reaction.\n');
    else
        nToPrint = min(numel(removedIdx),MAX_REJECTED_TO_PRINT);
        for k = 1:nToPrint
            motifIdx = removedIdx(k);
            fprintf('\nOriginal Block 4 motif %d rejected because of:\n', ...
                motifIdx);
            fprintf('  %s\n',membraneReason(motifIdx));
        end

        if nToPrint < numel(removedIdx)
            fprintf(['\nPrinted %d of %d rejected motifs. See ', ...
                     'Membrane_Filter_Report_Block_5 for all reasons.\n'], ...
                    nToPrint,numel(removedIdx));
        end
    end
end

fprintf('\n%d motifs remain for downstream filtering.\n',nF);
fprintf('Saved Block 5 results to:\n%s\n',outputFile);
