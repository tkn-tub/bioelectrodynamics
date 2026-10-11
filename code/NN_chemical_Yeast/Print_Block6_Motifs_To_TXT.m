%% Print selected Block 6 motifs and save them in a text file
%
% Run this script after the master file, while `model` is still available
% in the MATLAB workspace.  COBRA is not initialized again.

%% Make all project scripts and MAT files visible
addingPathParentFolderByName('code');

%% Select motifs
% Examples:
%   MOTIFS_TO_PRINT = 1;
%   MOTIFS_TO_PRINT = [1 4 10];
%   MOTIFS_TO_PRINT = 1:5;
if ~exist('MOTIFS_TO_PRINT','var') || isempty(MOTIFS_TO_PRINT)
    MOTIFS_TO_PRINT = 1;
end

% When several Block 0 paths reach the same m_in, select which one to print.
if ~exist('ADAPTER_MATCH_TO_PRINT','var') || ...
        isempty(ADAPTER_MATCH_TO_PRINT)
    ADAPTER_MATCH_TO_PRINT = 1;
end

% Default text-file name. A complete path can also be supplied beforehand.
if ~exist('OUTPUT_TXT_FILE','var') || isempty(OUTPUT_TXT_FILE)
    OUTPUT_TXT_FILE = 'Block_7_Selected_Motifs.txt';
end
OUTPUT_TXT_FILE = char(string(OUTPUT_TXT_FILE));

%% Load the yeast model through COBRA when it is not already in memory
if ~exist('model','var') || ~isstruct(model) || ...
        ~isfield(model,'S') || ~isfield(model,'mets') || ...
        ~isfield(model,'rxns')

    fprintf('Initializing COBRA and loading the yeast-GEM model...\n');
    initCobraToolbox(false);
    model = readCbModel('model_yeastGEM_v9_0_2.mat');
end

%% Validate the model required by Block 7
if ~exist('model','var') || ~isstruct(model) || ...
        ~isfield(model,'S') || ~isfield(model,'mets') || ...
        ~isfield(model,'rxns')
    error('PrintMotifs:MissingModel', ...
        'The yeast-GEM model could not be loaded through COBRA.');
end

%% Load saved results using the project path
load Block_6_no_internal_crosstalk.mat
load Block_0_Input_Adapter.mat

%% Determine how many strict motifs are available
if exist('results_nn_no_crosstalk','var') && ...
        isstruct(results_nn_no_crosstalk)
    nStrictMotifs = numel(results_nn_no_crosstalk);
elseif exist('Motifs_NN_Block_6','var') && istable(Motifs_NN_Block_6)
    nStrictMotifs = height(Motifs_NN_Block_6);
else
    error('PrintMotifs:InvalidBlock6File', ...
        ['Block_6_no_internal_crosstalk.mat does not contain ', ...
         'results_nn_no_crosstalk or Motifs_NN_Block_6.']);
end

if nStrictMotifs == 0
    error('PrintMotifs:NoMotifs', ...
        'Block 6 contains no motifs to print.');
end

%% Validate the requested motif indices
MOTIFS_TO_PRINT = unique(double(MOTIFS_TO_PRINT(:).'),'stable');

validSelection = isfinite(MOTIFS_TO_PRINT) & ...
    MOTIFS_TO_PRINT == floor(MOTIFS_TO_PRINT) & ...
    MOTIFS_TO_PRINT >= 1 & ...
    MOTIFS_TO_PRINT <= nStrictMotifs;

if ~all(validSelection)
    error('PrintMotifs:InvalidSelection', ...
        'MOTIFS_TO_PRINT must contain integers from 1 through %d.', ...
        nStrictMotifs);
end

%% Start a new text report
fileID = fopen(OUTPUT_TXT_FILE,'w');
if fileID < 0
    error('PrintMotifs:CannotWriteFile', ...
        'Could not create the text file: %s',OUTPUT_TXT_FILE);
end
fclose(fileID);

diary(OUTPUT_TXT_FILE);
diaryCleanup = onCleanup(@() diary('off'));

fprintf('\n############################################################\n');
fprintf('STRICT BLOCK 6 MOTIF REPORT\n');
fprintf('Available strict motifs : %d\n',nStrictMotifs);
fprintf('Selected motif indices  : %s\n',num2str(MOTIFS_TO_PRINT));
fprintf('############################################################\n');

%% Invoke Block 7 once for each requested motif
MOTIF_SOURCE = "block6";

for motifPosition = 1:numel(MOTIFS_TO_PRINT)
    MOTIF_TO_PRINT = MOTIFS_TO_PRINT(motifPosition);
    Block_7_Print_Selected_Motif
end

%% Close the report
diary off;
clear diaryCleanup;

fprintf('\nSaved %d selected motif report(s) to:\n%s\n', ...
    numel(MOTIFS_TO_PRINT),OUTPUT_TXT_FILE);
