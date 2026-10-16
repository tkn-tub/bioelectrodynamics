function [campRanking,campPaths,campRankInfo,outputFolder] = Rank_Existing_Circuit_Results(SAVED_RESULTS_FILE,options)
%RANK_EXISTING_CIRCUIT_RESULTS Rank saved motifs without repeating the search.
% Run with no arguments to select the newest results beside this function.
% Or: Rank_Existing_Circuit_Results('/full/path/Circuit_Search_Results.mat')
% Optional options struct overrides CAMP_* fields, e.g. CAMP_MAX_DEPTH = 8.
if nargin < 1, SAVED_RESULTS_FILE = ''; end
if nargin < 2, options = struct(); end
if isstring(SAVED_RESULTS_FILE) && isscalar(SAVED_RESULTS_FILE)
    SAVED_RESULTS_FILE = char(SAVED_RESULTS_FILE);
end
if ~ischar(SAVED_RESULTS_FILE) || size(SAVED_RESULTS_FILE,1) > 1
    error('Circuit:SavedResults','Use a text MAT file path or an empty path.');
end
if ~isstruct(options) || ~isscalar(options)
    error('Circuit:CAMPOptions','Options must be a scalar struct of CAMP_* settings.');
end
masterFolder = fileparts(mfilename('fullpath'));
addpath(masterFolder);

if isempty(SAVED_RESULTS_FILE)
    candidates = dir(fullfile(masterFolder,'Circuit_Search_Results'));
    latestTime = -Inf; chosenFile = '';
    for k = 1:numel(candidates)
        if ~candidates(k).isdir || any(strcmp(candidates(k).name,{'.','..'})), continue; end
        f = fullfile(masterFolder,'Circuit_Search_Results',candidates(k).name,'Circuit_Search_Results.mat');
        details = dir(f);
        if ~isempty(details) && details(1).datenum > latestTime
            chosenFile = f; latestTime = details(1).datenum;
        end
    end
    if isempty(chosenFile), chosenFile = which('Circuit_Search_Results.mat'); end
    if isempty(chosenFile)
        error('Circuit:SavedResults','Set SAVED_RESULTS_FILE to your saved Circuit_Search_Results.mat.');
    end
    SAVED_RESULTS_FILE = chosenFile;
end
saved = load(SAVED_RESULTS_FILE);
if ~all(isfield(saved,{'model','motifs','cfg','pattern','searchInfo'}))
    error('Circuit:SavedResults','Expected model, motifs, cfg, pattern and searchInfo in saved results.');
end
model = saved.model; motifs = saved.motifs; searchInfo = saved.searchInfo;
cfg = circuitSearchDefaults();
settings = fieldnames(saved.cfg);
for k = 1:numel(settings), cfg.(settings{k}) = saved.cfg.(settings{k}); end
overrideNames = fieldnames(options);
for k = 1:numel(overrideNames)
    field = overrideNames{k};
    if ~strncmp(field,'CAMP_',5) || ~isfield(cfg,field)
        error('Circuit:CAMPOptions','Unknown cAMP ranking option: %s',field);
    end
    cfg.(field) = options.(field);
end
pattern = circuitPattern(cfg.NEURON_INPUT_MODE);
if ~isequal(pattern.T,saved.pattern.T) || ~isequal(pattern.rxnNames,saved.pattern.rxnNames)
    error('Circuit:SavedPattern','Saved results have a different topology; they cannot be relabelled by index.');
end
[campRanking,campPaths,campRankInfo] = rankCircuitInputsByCAMP(model,motifs,cfg,pattern);
outputFolder = fullfile(masterFolder,'CAMP_Ranking_Results',datestr(now,'yyyymmdd_HHMMSS'));
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
reportFile = fullfile(outputFolder,'CAMP_Input_Ranking.txt');
writeCAMPInputRankingReport(reportFile,model,motifs,campRanking,campPaths,campRankInfo,pattern,cfg);
writetable(campRanking,fullfile(outputFolder,'CAMP_Input_Ranking.csv'));
writeMetabolicCircuitReport(fullfile(outputFolder,'Selected_Circuit_Motifs_Relabelled.txt'), ...
    model,motifs,searchInfo,pattern,cfg,campRanking,campPaths,campRankInfo);
save(fullfile(outputFolder,'Circuit_CAMP_Ranking.mat'),'model','motifs','searchInfo', ...
    'cfg','pattern','campRanking','campPaths','campRankInfo','SAVED_RESULTS_FILE','-v7');
fprintf('Ranked saved results: %s\n',SAVED_RESULTS_FILE);
disp(campRanking(:,{'Rank','MotifIndex','InputMetabolite','Distance','Status'}));
fprintf('cAMP ranking report: %s\n',reportFile);
if ~isempty(campPaths) && isfinite(campPaths(1).distance)
    bestMotif = campRanking.MotifIndex(1);
    fprintf('\nClosest returned input: motif %d\n',bestMotif);
    printCAMPConnection(1,model,motifs(bestMotif),campPaths(1),campRankInfo,pattern,cfg);
elseif ~isempty(campPaths)
    fprintf('No finite cAMP-to-input connection was found; see the ranking report.\n');
end
end
