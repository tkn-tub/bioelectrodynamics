function test_compact_checkpoints
% Run in MATLAB after adding this folder to the path. No COBRA is needed.
testFolder = tempname;
mkdir(testFolder);
cleanup = onCleanup(@() rmdir(testFolder,'s'));

sampleTable = table((1:20).',(101:120).', ...
    'VariableNames',{'m_in','r_S'});
checkpoint = struct;
checkpoint.Motifs_Diff_Block_1 = sampleTable;
testFile = fullfile(testFolder,'table_only.mat');
saveCompactCheckpoint(testFile,checkpoint);
saved = load(testFile);
assert(isequal(saved.Motifs_Diff_Block_1,sampleTable));
assert(numel(fieldnames(saved)) == 1);

% Empty filter results must retain column names and types for later loading.
checkpoint.Motifs_Diff_Block_1 = sampleTable([],:);
saveCompactCheckpoint(testFile,checkpoint);
saved = load(testFile);
assert(isequal(saved.Motifs_Diff_Block_1,sampleTable([],:)));

fileNames = {'Block_1_differentiator.mat','Block_2_differentiator.mat', ...
    'Block_5_no_cell_membrane.mat','Block_6_no_internal_crosstalk.mat'};
tableNames = {'Motifs_Diff_Block_1','Motifs_Diff_Block_2', ...
    'Motifs_NN_Block_5','Motifs_NN_Block_6'};
structNames = {'results_diff','results_diff2', ...
    'results_nn_no_membrane','results_nn_no_crosstalk'};
for b = 1:4
    legacy = struct;
    legacy.(tableNames{b}) = sampleTable;
    legacy.(structNames{b}) = table2struct(sampleTable);
    legacy.rawResults = legacy.(structNames{b});
    legacy.stage2ResultsByPair = {legacy.(structNames{b})};
    legacy.fullJoinRows_beforeCrossTalkFilter = legacy.(structNames{b});
    sourceFile = fullfile(testFolder,fileNames{b});
    save(sourceFile,'-struct','legacy','-v7.3');
end

report = Compact_Existing_Block_MATs(testFolder);
assert(height(report) == 4);
for b = 1:4
    compactFile = fullfile(testFolder,'compact_results',fileNames{b});
    saved = load(compactFile);
    assert(isequal(saved.(tableNames{b}),sampleTable));
    assert(~isfield(saved,structNames{b}));
    assert(~isfield(saved,'rawResults'));
    assert(~isfield(saved,'stage2ResultsByPair'));
    assert(~isfield(saved,'fullJoinRows_beforeCrossTalkFilter'));
    original = load(fullfile(testFolder,fileNames{b}));
    assert(isfield(original,structNames{b}));
end

% A second conversion must not replace an existing destination checkpoint.
try
    Compact_Existing_Block_MATs(testFolder);
    error('Test:ExpectedError','Expected DestinationExists error.');
catch err
    assert(strcmp(err.identifier,'Checkpoint:DestinationExists'));
end
fprintf('All compact checkpoint tests passed.\n');
clear cleanup;
end
