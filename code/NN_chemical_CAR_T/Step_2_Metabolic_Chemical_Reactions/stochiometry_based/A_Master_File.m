%% Fresh metabolic circuit mapping: edit this section, then run this file.
% HPC: matlab -batch "A_Master_File"
% Keep all companion .m files in this folder (inside your code tree).
clc;
clearvars;
close all;

masterFolder = fileparts(mfilename('fullpath'));
addpath(masterFolder);
addingPathParentFolderByName('code');

%% Model discovery: which() finds files in any folder added above.
% HPC_MODEL_FILE = 'model_yeastGEM_v9_0_2.mat';
HPC_MODEL_FILE = 'Puniya2021_TNM1055.xml';
MODEL_VARIABLE = ''; % MAT only: leave empty if there is just one model.

%% Search settings
cfg = circuitSearchDefaults();
cfg.FORWARD_ONLY = true;      % Both orientations, only where bounds permit.
cfg.USE_PARALLEL = true;     % Requires Parallel Computing Toolbox.
cfg.NUM_WORKERS = [];         % []: use pool/profile default; or e.g. 8.
cfg.PARALLEL_PROFILE = '';    % '': current default; or HPC profile name.
cfg.ONE_MOTIF_PER_INPUT = true; % One circuit for each distinct m_in.
cfg.MAX_MOTIFS = 10^6;          % Up to 20 different inputs; Inf: all feasible inputs.
cfg.MAX_SECONDS = 3600;       % Inf: no search time limit.
cfg.VERBOSE = true;

% Complete EC annotations are the default evidence for enzymatic mediation.
% EC_OR_GPR is a weaker alternative when EC annotations are incomplete.
cfg.ENZYME_EVIDENCE = 'EC_OR_GPR';    % 'EC' or 'EC_OR_GPR'.
% By default, any common EC number OR gene excludes two reactions.
% This is conservative: EC numbers identify activities, not unique proteins.
% Optional curated enzyme IDs: name an nRxn-by-1 model field containing IDs.
% With this field, its IDs replace EC grouping; shared genes still conflict.
cfg.ENZYME_ID_FIELD = '';

% Signal metabolites must remain distinct and participate only in the
% selected reactions specified by circuitPattern. Non-signal substrates,
% cofactors and products may be shared freely between selected reactions.
% This checks the selected circuit, not other reactions elsewhere in the model.
% Their assumed saturation is a kinetic assumption, not verified by this search.
cfg.METABOLITE_SCOPE = 'roles'; % 'all' additionally restricts shared auxiliaries.
cfg.METABOLITE_IDENTITY = 'chemical'; % 'chemical' or 'model_id'.
% Auxiliary exemptions below are needed only if scope is changed to 'all'.
cfg.AUXILIARY_REUSE_ALLOWED = {}; % Example: {'h2o','h'}.
cfg.EXCLUDED_ROLE_METABOLITES = {}; % Example: {'h2o','h'}.
cfg.COMPARTMENT = '';         % '': any; 'c': all reaction species must be c.

% The arrow m_c -> r_Nw is interpreted as a consumed substrate by default.
cfg.NEURON_INPUT_MODE = 'substrate'; % 'substrate' or 'activation'.
% Activation mode REQUIRES supplied {reactionID, metaboliteID} pairs.
% Optional third column: an evidence/source description. Never inferred from S.
cfg.ACTIVATION_PAIRS = cell(0,2);

% Optional anchors, using existing names from your example report:
cfg.FIXED_METABOLITES = struct();
cfg.FIXED_REACTIONS = struct();
% cfg.FIXED_METABOLITES.m_in = 'camp[c]';
% cfg.FIXED_REACTIONS.r_Ff = 'HMR_4672';

%% cAMP -> m_in connection ranking
cfg.CAMP_SOURCE_METABOLITES = {}; % Auto-detect; or exact model IDs {'camp[c]'}.
cfg.CAMP_SOURCE_COMPARTMENT = 'c'; % Set to the model's cytosol code; '' = any.
cfg.CAMP_ALLOW_TRANSPORT = false;
cfg.CAMP_MAX_DEPTH = 10;         % Conversion steps; 0 = only exact cAMP inputs.
cfg.CAMP_EXCLUDED_PATH_METABOLITES = {}; % Extra path-node IDs/base IDs to omit.
PRINT_CAMP_CHAIN_TO_CONSOLE = true; % Print the closest returned connection.

%% Load, search, and print; no Blocks 0-7 or bridge runners are needed.
[model, modelFile] = loadCircuitModel(HPC_MODEL_FILE, MODEL_VARIABLE);
cfg.MODEL_FILE = modelFile;
[motifs, searchInfo, pattern] = findMetabolicCircuit(model, cfg);

outputFolder = fullfile(masterFolder, 'Circuit_Search_Results', ...
    datestr(now, 'yyyymmdd_HHMMSS'));
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end
outputFile = fullfile(outputFolder, 'Selected_Circuit_Motifs.txt');
writeMetabolicCircuitReport(outputFile, model, motifs, searchInfo, pattern, cfg);
resultsFile = fullfile(outputFolder,'Circuit_Search_Results.mat');
% Preserve the search before ranking, so a cAMP setting error cannot lose it.
save(resultsFile, ...
    'model', 'motifs', 'searchInfo', 'pattern', 'cfg', '-v7');
fprintf('\nSearch status: %s; motifs returned: %d\n', ...
    searchInfo.status, numel(motifs));
fprintf('Text report: %s\n', outputFile);

%% Rank the returned solutions and print the full connector reactions.
[campRanking,campPaths,campRankInfo] = rankCircuitInputsByCAMP(model,motifs,cfg,pattern);
campReportFile = fullfile(outputFolder,'CAMP_Input_Ranking.txt');
writeCAMPInputRankingReport(campReportFile,model,motifs,campRanking,campPaths,campRankInfo,pattern,cfg);
% Include each motif's connector in the same report as its circuit reactions.
writeMetabolicCircuitReport(outputFile,model,motifs,searchInfo,pattern,cfg, ...
    campRanking,campPaths,campRankInfo);
writetable(campRanking,fullfile(outputFolder,'CAMP_Input_Ranking.csv'));
save(resultsFile,'campRanking','campPaths','campRankInfo','-v7','-append');
disp(campRanking(:,{'Rank','MotifIndex','InputMetabolite','Distance','Status'}));
fprintf('cAMP ranking and reaction chains: %s\n',campReportFile);
if PRINT_CAMP_CHAIN_TO_CONSOLE && ~isempty(campPaths) && isfinite(campPaths(1).distance)
    bestMotif = campRanking.MotifIndex(1);
    fprintf('\nClosest returned input: motif %d\n',bestMotif);
    printCAMPConnection(1,model,motifs(bestMotif),campPaths(1),campRankInfo,pattern,cfg);
elseif PRINT_CAMP_CHAIN_TO_CONSOLE && ~isempty(campPaths)
    fprintf('No finite cAMP-to-input connection was found; see the ranking report.\n');
end

