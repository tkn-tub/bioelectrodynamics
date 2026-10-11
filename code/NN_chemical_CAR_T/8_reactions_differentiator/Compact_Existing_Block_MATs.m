function conversionReport = Compact_Existing_Block_MATs(sourceFolder,destinationFolder)
%COMPACT_EXISTING_BLOCK_MATS Convert Blocks 1, 2, 5, 6 without rerunning searches.
% Run from the project folder:
%   Compact_Existing_Block_MATs(pwd)
% Originals are NEVER overwritten. Copies with the SAME filenames are saved
% in sourceFolder/compact_results. Move/copy them into the project folder only
% after backing up originals and installing the updated master and scripts.
% Only result tables and small diagnostics are loaded from old checkpoints.

if nargin < 1 || isempty(sourceFolder), sourceFolder = pwd; end
if nargin < 2 || isempty(destinationFolder)
    destinationFolder = fullfile(sourceFolder,'compact_results');
end
sourceFolder = char(string(sourceFolder));
destinationFolder = char(string(destinationFolder));
if ~isfolder(destinationFolder), mkdir(destinationFolder); end

fileNames = { ...
    'Block_1_differentiator.mat', ...
    'Block_2_differentiator.mat', ...
    'Block_5_no_cell_membrane.mat', ...
    'Block_6_no_internal_crosstalk.mat'};
tableNames = { ...
    'Motifs_Diff_Block_1','Motifs_Diff_Block_2', ...
    'Motifs_NN_Block_5','Motifs_NN_Block_6'};
legacyNames = { ...
    'results_diff','results_diff2', ...
    'results_nn_no_membrane','results_nn_no_crosstalk'};
metadata = { ...
    {'Differentiator_Input_Anchor_Block_1','BLOCK0_RESULTS_FILE_USED', ...
     'inputAnchorIdx','inputAnchorIDs','REQUIRE_DISTINCT_ADAPTER_RXNS', ...
     'DIFF_INPUT_MET_ID','inputMetIdx','inputConsumerRxns'}, ...
    {'uniqueInputPairs','pairOfBlock1','pairSolutionCounts', ...
     'potentiallyCappedPairs','runInParallel', ...
     'FILTER_INTERNAL_CROSSTALK_DURING_SEARCH','nRejectedCrossTalk2', ...
     'MAX_HITS_BLOCK2','MAX_STAGE2_PER_PAIR', ...
     'MAX_BRANCH_OPTIONS_PER_INPUT','MAX_STAGE2_OPTION_PAIRS_TESTED'}, ...
    {'Membrane_Filter_Report_Block_5','Cell_Membrane_Reactions_Block_5', ...
     'indexBeforeMembraneFilter','keepMotif','isCellMembraneRxn','rxnFields'}, ...
    {'CrossTalk_Filter_Report_Block_6','indexBeforeCrossTalkFilter', ...
     'keepCrossTalkFree','rxnFields','metFields','allowedConsumers'}};

reportRows = struct('File',{},'OriginalMB',{},'CompactMB',{},'OutputFile',{});
for b = 1:numel(fileNames)
    sourceFile = fullfile(sourceFolder,fileNames{b});
    if b == 2 && ~isfile(sourceFile)
        sourceFile = fullfile(sourceFolder,'Block_2_second_differentiator.mat');
    end
    if ~isfile(sourceFile)
        warning('Checkpoint:MissingFile','Skipping missing file: %s',sourceFile);
        continue;
    end

    destinationFile = fullfile(destinationFolder,fileNames{b});
    % Reject equal folders, symbolic-link aliases and already-existing
    % destination files. The converter never modifies any existing MAT file.
    if isfile(destinationFile)
        error('Checkpoint:DestinationExists', ...
            ['Destination already exists: %s. Choose a new destination ', ...
             'folder; existing files will not be overwritten.'],destinationFile);
    end

    info = whos('-file',sourceFile);
    names = {info.name};
    keepMetadata = metadata{b}(ismember(metadata{b},names));
    if ismember(tableNames{b},names)
        payload = load(sourceFile,tableNames{b},keepMetadata{:});
    elseif ismember(legacyNames{b},names)
        payload = load(sourceFile,legacyNames{b},keepMetadata{:});
        payload.(tableNames{b}) = struct2table(payload.(legacyNames{b})(:));
        payload = rmfield(payload,legacyNames{b});
    else
        error('Checkpoint:MissingMotifs', ...
            '%s lacks both %s and %s.',sourceFile,tableNames{b},legacyNames{b});
    end
    if ~istable(payload.(tableNames{b}))
        error('Checkpoint:InvalidMotifs','%s must contain a table.',tableNames{b});
    end
    payload.CheckpointStorageVersion = 1;
    originalInfo = dir(sourceFile);
    compactInfo = saveCompactCheckpoint(destinationFile,payload);
    entry = struct;
    entry.File = string(fileNames{b});
    entry.OriginalMB = originalInfo.bytes/1e6;
    entry.CompactMB = compactInfo.bytes/1e6;
    entry.OutputFile = string(destinationFile);
    reportRows(end+1) = entry; %#ok<AGROW>
    clear payload;
end
conversionReport = struct2table(reportRows);
disp(conversionReport);
fprintf('Original files were left unchanged.\n');
end
