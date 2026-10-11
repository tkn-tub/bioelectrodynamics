%% Print selected strict Block 6 motifs in the detailed A/B report format.
% Put this file and Block_7_Print_Selected_Motif.m together in your code folder.
% The first export uses the SAME model already loaded by the master, or a
% local model cache. If neither exists, the local XML can be read via COBRA.
addingPathParentFolderByName('code');
addpath(fileparts(mfilename('fullpath')));

% Run from the folder containing the NEW polarity/Hill Block 6 checkpoint.
% The printer rejects older checkpoints and displays the actual resolved path.
MOTIF_TO_PRINT = 1; % one row; or e.g. [1 4 8] to print several selected rows
PRINT_RESULTS_FILE = 'Block_6_no_internal_crosstalk.mat';
PRINT_TXT_FILE = 'Block_7_Selected_Motifs.txt';
PRINT_MODEL_FILE = 'Puniya2021_TNM1055.xml';
PRINT_MODEL_CACHE_FILE = 'Block_7_Print_Model.mat';
PRINT_ALLOW_XML_LOAD = true;
PRINT_INITIALIZE_COBRA_WHEN_NEEDED = false; % true only when your session needs initialization
PRINT_SAVE_MODEL_CACHE = true;
ADAPTER_MATCH_TO_PRINT = 1; % used only when the motif has no stored adapterPathIndex

Block_7_Print_Selected_Motif;

function addingPathParentFolderByName(targetName)
currFolder = pwd;
while true
    [parentFolder,currentName] = fileparts(currFolder);
    if strcmpi(currentName,targetName)
        addpath(genpath(currFolder));
        fprintf('Adding matlab path to: %s\n',currFolder);
        return;
    end
    if isempty(parentFolder) || strcmp(currFolder,parentFolder)
        error('PrintPath:MissingParent','Folder named "%s" not found in any parent directory.',targetName);
    end
    currFolder = parentFolder;
end
end
