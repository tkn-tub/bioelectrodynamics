%% Print strict Block 6 motifs using the matching metabolic model
%
% Run with the MATLAB current folder containing the CAR-T results and
% Block_7_Print_Selected_Motif.m. The TXT is saved in that same folder.
% Run it independently, or after the master. CAR-T is the default model.
% Optional overrides before running:
%   BLOCK0_REQUIRE_EC = false; % passive adapter; Blocks 1-4 remain strict
%   PRINT_BLOCK6_FILE = 'Block_6_EC_annotated.mat'; % optional EC-filtered file
%   MOTIFS_TO_PRINT = [1 4 10];
%   ADAPTER_MATCH_TO_PRINT = 1;
%   OUTPUT_TXT_FILE = 'Block_7_Selected_Motifs.txt';
%
% The configured model is loaded afresh to avoid a stale yeast model in the
% workspace. COBRA initialization is invoked only if readCbModel is missing.
% Old and compact checkpoints are supported. Result files are resolved only
% in the current folder, never by searching all project folders on the path.
% The model XML is found by filename on the path added below.
% EC filtering runs automatically AFTER loading the model and result files.
% No model variable or previous master-file run is required in the workspace.

addingPathParentFolderByName('code');

PRINT_RESULTS_DIR = pwd;
addpath(PRINT_RESULTS_DIR);
PRINT_MODEL_FILE = 'Puniya2021_TNM1055.xml';
if ~exist('BLOCK0_REQUIRE_EC','var'), BLOCK0_REQUIRE_EC = false; end
if ~exist('PRINT_BLOCK6_FILE','var') || isempty(PRINT_BLOCK6_FILE)
    PRINT_BLOCK6_FILE = 'Block_6_no_internal_crosstalk.mat';
end
PRINT_BLOCK6_FILE = char(string(PRINT_BLOCK6_FILE));
if ~isempty(fileparts(PRINT_BLOCK6_FILE))
    error('PrintMotifs:InvalidCheckpointName', ...
        'PRINT_BLOCK6_FILE must be a filename in the current folder.');
end
if ~exist('MOTIFS_TO_PRINT','var') || isempty(MOTIFS_TO_PRINT)
    MOTIFS_TO_PRINT = 1;
end
if ~exist('ADAPTER_MATCH_TO_PRINT','var') || isempty(ADAPTER_MATCH_TO_PRINT)
    ADAPTER_MATCH_TO_PRINT = 1;
end
if ~exist('OUTPUT_TXT_FILE','var') || isempty(OUTPUT_TXT_FILE)
    OUTPUT_TXT_FILE = 'Block_7_Selected_Motifs.txt';
end
% Keep the selected filename, always writing to the current folder.
[~,reportOutputName,reportOutputExt] = fileparts(char(string(OUTPUT_TXT_FILE)));
if isempty(reportOutputName)
    error('PrintMotifs:InvalidOutputName','Provide a TXT output filename.');
end
OUTPUT_TXT_FILE = fullfile(PRINT_RESULTS_DIR,[reportOutputName,reportOutputExt]);

reportBlock6File = fullfile(PRINT_RESULTS_DIR,PRINT_BLOCK6_FILE);
if strcmpi(PRINT_BLOCK6_FILE,'Block_6_EC_annotated.mat') && ~isfile(reportBlock6File)
    % Build the EC-filtered checkpoint on the first independent export.
    reportBlock6File = fullfile(PRINT_RESULTS_DIR,'Block_6_no_internal_crosstalk.mat');
end
reportBlock0File = fullfile(PRINT_RESULTS_DIR,'Block_0_Input_Adapter.mat');
reportPrinterFile = fullfile(PRINT_RESULTS_DIR,'Block_7_Print_Selected_Motif.m');
reportFilterFile = fullfile(PRINT_RESULTS_DIR,'Block_6B_Require_EC_Annotations.m');
reportECHelperFile = fullfile(PRINT_RESULTS_DIR,'getReactionECAnnotations.m');
for reportRequiredFile = {reportBlock6File,reportPrinterFile,reportFilterFile,reportECHelperFile}
    if ~isfile(reportRequiredFile{1})
        error('PrintMotifs:MissingProjectFile', ...
            ['Required file not found: %s. Run this script from the ', ...
             'folder containing the matching Block 0/6 files and Block 7.'], ...
            reportRequiredFile{1});
    end
end

reportResolvedModel = reportResolveModelFile(PRINT_MODEL_FILE,PRINT_RESULTS_DIR);
fprintf('\nModel file : %s\nBlock 6    : %s\n', ...
    reportResolvedModel,reportBlock6File);
if exist('readCbModel','file') ~= 2
    if exist('initCobraToolbox','file') ~= 2
        error('PrintMotifs:MissingCOBRA', ...
            'COBRA is not on the MATLAB path. Add it before running this script.');
    end
    initCobraToolbox(false);
end
reportLoadedModel = readCbModel(reportResolvedModel);
reportValidateModel(reportLoadedModel);

% Load only the named results. No model or stale compatibility structs are
% imported from the general workspace or from other MAT variables.
reportVariables = whos('-file',reportBlock6File);
reportVariableNames = {reportVariables.name};
reportECVariables = {'EC_ModelMetIDs','EC_ModelRxnIDs','EC_Reaction_Annotations'};
if all(ismember(reportECVariables,reportVariableNames))
    reportECMetadata = load(reportBlock6File,reportECVariables{:});
    if ~isequal(string(reportLoadedModel.mets(:)),reportECMetadata.EC_ModelMetIDs(:)) || ...
            ~isequal(string(reportLoadedModel.rxns(:)),reportECMetadata.EC_ModelRxnIDs(:))
        error('PrintMotifs:ECModelMismatch', ...
            'The loaded model IDs or their ordering disagree with the EC-filtered checkpoint.');
    end
    reportECAnnotations = reportECMetadata.EC_Reaction_Annotations;
    if ~istable(reportECAnnotations) || ...
            ~all(ismember({'reactionID','completeECNumbers'},reportECAnnotations.Properties.VariableNames)) || ...
            ~isequal(string(reportECAnnotations.reactionID(:)),string(reportLoadedModel.rxns(:)))
        error('PrintMotifs:InvalidECMetadata','EC annotation metadata does not match this model.');
    end
    % Recheck EC annotations against the freshly loaded model below.
    clear reportECMetadata reportECAnnotations;
end
if ismember('Motifs_NN_Block_6',reportVariableNames)
    reportSaved = load(reportBlock6File,'Motifs_NN_Block_6');
    reportMotifs = reportSaved.Motifs_NN_Block_6;
elseif ismember('results_nn_no_crosstalk',reportVariableNames)
    reportSaved = load(reportBlock6File,'results_nn_no_crosstalk');
    reportMotifs = struct2table(reportSaved.results_nn_no_crosstalk(:));
else
    error('PrintMotifs:MissingMotifs', ...
        'The Block 6 file contains neither its motif table nor its result structs.');
end
clear reportSaved;
if ~istable(reportMotifs) || height(reportMotifs) == 0
    error('PrintMotifs:NoStrictMotifs', ...
        ['The selected checkpoint contains no motifs to print. ', ...
         'Inspect EC_Filter_Report if the EC filter was applied.']);
end
% The EC-filtered file includes adapter paths under the configured passive
% adapter policy. Computational reaction roles remain fully EC annotated.
if ismember('InputAdapterPaths',reportVariableNames)
    reportSavedAdapter = load(reportBlock6File,'InputAdapterPaths');
    reportAdapterSource = reportBlock6File;
else
    if ~isfile(reportBlock0File)
        error('PrintMotifs:MissingProjectFile', ...
            'Block_0_Input_Adapter.mat is missing from the current folder.');
    end
    reportSavedAdapter = load(reportBlock0File,'InputAdapterPaths');
    reportAdapterSource = reportBlock0File;
end
if ~isfield(reportSavedAdapter,'InputAdapterPaths') || ...
        ~isstruct(reportSavedAdapter.InputAdapterPaths)
    error('PrintMotifs:MissingAdapterPaths', ...
        'The Block 0 file must contain the InputAdapterPaths struct array.');
end
reportAdapterPaths = reportSavedAdapter.InputAdapterPaths(:);
fprintf('Adapter source: %s\n',reportAdapterSource);
clear reportSavedAdapter;

% Verify the source topology before writing a filtered checkpoint. Run the
% filter with an explicit model and results in a private function workspace,
% rather than calling it at the top of this script before model loading.
reportValidateMotifs(reportLoadedModel,reportMotifs);
[reportMotifs,reportAdapterPaths,reportECAnnotations] = reportApplyECFilter( ...
    reportLoadedModel,reportMotifs,reportAdapterPaths,reportFilterFile,BLOCK0_REQUIRE_EC);
reportBlock6File = fullfile(PRINT_RESULTS_DIR,'Block_6_EC_annotated.mat');
if height(reportMotifs) == 0
    error('PrintMotifs:NoECAnnotatedMotifs', ...
        ['No motif with EC-annotated computational reactions survived. Diagnostics were saved ', ...
         'in Block_6_EC_annotated.mat. The TXT file was not changed.']);
end
% Display the complete codes resolved from the current model's annotations.
reportLoadedModel.eccodes = cellstr(reportECAnnotations.completeECNumbers);
clear reportECAnnotations;
if ~isnumeric(MOTIFS_TO_PRINT) || ~isreal(MOTIFS_TO_PRINT) || ...
        ~isvector(MOTIFS_TO_PRINT)
    error('PrintMotifs:InvalidSelection','MOTIFS_TO_PRINT must be a numeric vector.');
end
MOTIFS_TO_PRINT = unique(double(MOTIFS_TO_PRINT(:).'),'stable');
if any(~isfinite(MOTIFS_TO_PRINT) | MOTIFS_TO_PRINT < 1 | ...
        MOTIFS_TO_PRINT ~= floor(MOTIFS_TO_PRINT) | ...
        MOTIFS_TO_PRINT > height(reportMotifs))
    error('PrintMotifs:InvalidSelection', ...
        'Select EC-filtered motif rows from 1 through %d.',height(reportMotifs));
end
reportValidateAdapters(reportLoadedModel,reportMotifs(MOTIFS_TO_PRINT,:), ...
    reportAdapterPaths,ADAPTER_MATCH_TO_PRINT);
fprintf('Validated model: %d metabolites, %d reactions.\n',size(reportLoadedModel.S));
reportWriteTXT(reportLoadedModel,reportMotifs,reportAdapterPaths, ...
    MOTIFS_TO_PRINT,ADAPTER_MATCH_TO_PRINT,reportPrinterFile, ...
    char(string(OUTPUT_TXT_FILE)),reportResolvedModel,reportBlock6File,BLOCK0_REQUIRE_EC);

%% Local functions
function [filteredMotifs,filteredPaths,annotations] = reportApplyECFilter( ...
        model,Motifs_NN_Block_6,InputAdapterPaths,filterFile,BLOCK0_REQUIRE_EC)
% Block 6B sees only the model and checkpoint passed by this export run.
run(filterFile);
filteredMotifs = Motifs_NN_Block_6;
filteredPaths = InputAdapterPaths;
annotations = EC_Reaction_Annotations;
end

function filePath = reportResolveModelFile(requested,resultsFolder)
requested = char(string(requested));
if isfile(requested), filePath = requested; return; end
candidate = fullfile(resultsFolder,requested);
if isfile(candidate), filePath = candidate; return; end
matches = which(requested,'-all');
if ischar(matches), matches = cellstr(matches); end
matches = unique(matches,'stable');
if numel(matches) == 1 && isfile(matches{1})
    filePath = matches{1};
else
    error('PrintMotifs:ModelFileResolution', ...
        ['Could not resolve a unique model file %s. Put the matching XML ', ...
         'in the current folder or make one copy visible through ', ...
         'addingPathParentFolderByName(''code'').'],requested);
end
end

function reportValidateModel(model)
if ~isstruct(model) || ~isscalar(model) || ...
        ~all(isfield(model,{'S','mets','rxns'}))
    error('PrintMotifs:InvalidModel','The loaded model must contain S, mets and rxns.');
end
if numel(model.mets) ~= size(model.S,1) || numel(model.rxns) ~= size(model.S,2)
    error('PrintMotifs:ModelSizeMismatch','The model ID lists disagree with model.S.');
end
end

function reportValidateMotifs(model,motifs)
metFields = {'m_in','m_S1','m_F1','m_Sg1','m_Fg1', ...
    'm_S2','m_F2','m_Sg2','m_d','m_c','m_Cs','m_out'};
rxnFields = {'r_S1','r_Sg1','r_Ss1','r_F1','r_Fg1','r_Fs1','r_I1', ...
    'r_S2','r_Sg2','r_Ss2','r_F2','r_Fg2','r_Fs2','r_I2', ...
    'r_C','r_CI1','r_CI2','r_omega','r_b'};
for group = 1:2
    if group == 1, fields = metFields; limit = size(model.S,1);
    else, fields = rxnFields; limit = size(model.S,2); end
    for f = 1:numel(fields)
        field = fields{f};
        if ~ismember(field,motifs.Properties.VariableNames)
            error('PrintMotifs:IncompleteMotif','Missing motif column: %s',field);
        end
        values = motifs.(field);
        if ~isnumeric(values) || ~isreal(values) || size(values,2) ~= 1
            error('PrintMotifs:InvalidColumn','%s must be a numeric index column.',field);
        end
        bad = find(~isfinite(values) | values < 1 | values > limit | ...
            values ~= floor(values),1);
        if ~isempty(bad)
            error('PrintMotifs:WrongModelOrResults', ...
                ['%s=%g is outside the loaded model range 1:%d. Model has ', ...
                 '%d metabolites and %d reactions. Use the exact model and ', ...
                 'Block 0/6 files from the same run; indices cannot be clipped ', ...
                 'or reused between yeast and CAR-T models.'], ...
                 field,values(bad),limit,size(model.S,1),size(model.S,2));
        end
    end
end

% Range checks alone cannot detect a different model with enough rows.
% Check the required signed edges of the abstract chemical topology too.
edges = { ...
    'm_in','r_S1',-1; 'm_S1','r_S1',1; ...
    'm_S1','r_Sg1',-1; 'm_Sg1','r_Sg1',1; 'm_S1','r_Ss1',-1; ...
    'm_in','r_F1',-1; 'm_F1','r_F1',1; ...
    'm_F1','r_Fg1',-1; 'm_Fg1','r_Fg1',1; 'm_F1','r_Fs1',-1; ...
    'm_Sg1','r_I1',-1; 'm_Fg1','r_I1',-1; ...
    'm_Sg1','r_S2',-1; 'm_S2','r_S2',1; ...
    'm_S2','r_Sg2',-1; 'm_Sg2','r_Sg2',1; 'm_S2','r_Ss2',-1; ...
    'm_Fg1','r_F2',-1; 'm_F2','r_F2',1; ...
    'm_F2','r_Fg2',-1; 'm_d','r_Fg2',1; 'm_F2','r_Fs2',-1; ...
    'm_Sg2','r_I2',-1; 'm_d','r_I2',-1; ...
    'm_d','r_C',-1; 'm_c','r_C',1; ...
    'm_c','r_CI1',-1; 'm_Cs','r_CI1',-1; ...
    'm_d','r_CI2',-1; 'm_Cs','r_CI2',-1; ...
    'm_c','r_omega',-1; 'm_out','r_omega',1; 'm_out','r_b',1};
for e = 1:size(edges,1)
    indices = sub2ind(size(model.S),double(motifs.(edges{e,1})), ...
        double(motifs.(edges{e,2})));
    if any(full(model.S(indices))*edges{e,3} <= 0)
        error('PrintMotifs:TopologyMismatch', ...
            ['Saved indices do not reproduce the required %s / %s edge ', ...
             'in this model. Check the exact model version, ordering and ', ...
             'result folder before printing.'],edges{e,1},edges{e,2});
    end
end
end

function reportValidateAdapters(model,motifs,paths,choice)
if ~isnumeric(choice) || ~isscalar(choice) || ~isfinite(choice) || ...
        ~isreal(choice) || choice < 1 || choice ~= floor(choice)
    error('PrintMotifs:InvalidAdapterSelection', ...
        'ADAPTER_MATCH_TO_PRINT must be a positive integer.');
end
if isempty(paths) || ~all(isfield(paths, ...
        {'anchorMetIdx','inputMetIdx','metaboliteIndices','reactionIndices'}))
    error('PrintMotifs:InvalidAdapters','Block 0 has no complete adapter paths.');
end
anchors = double([paths.anchorMetIdx]);
for row = 1:height(motifs)
    matches = find(anchors == motifs.m_in(row));
    if numel(matches) < choice
        error('PrintMotifs:AdapterMismatch', ...
            ['No requested matching Block 0 path for m_in=%g. Block 0 and ', ...
             'Block 6 must come from the same run.'],motifs.m_in(row));
    end
    path = paths(matches(choice));
    if ~isscalar(path.inputMetIdx) || ~isscalar(path.anchorMetIdx)
        error('PrintMotifs:InvalidAdapterEndpoint','Block 0 endpoints must be scalar indices.');
    end
    mets = [double(path.inputMetIdx),double(path.anchorMetIdx), ...
        double(path.metaboliteIndices(:).')];
    rxns = double(path.reactionIndices(:).');
    if any(~isfinite(mets) | mets < 1 | mets > size(model.S,1) | mets ~= floor(mets)) || ...
            any(~isfinite(rxns) | rxns < 1 | rxns > size(model.S,2) | rxns ~= floor(rxns))
        error('PrintMotifs:AdapterIndexMismatch','Block 0 path indices do not match this model.');
    end
    for check = {'input','anchor'}
        idField = [check{1},'MetID'];
        idxField = [check{1},'MetIdx'];
        if isfield(path,idField) && ...
                string(model.mets(path.(idxField))) ~= string(path.(idField))
            error('PrintMotifs:AdapterIDMismatch', ...
                'Block 0 %s does not match the model ID at its saved index.',idField);
        end
    end

    % Verify the complete adapter path against saved IDs when available.
    pathMetIdx = double(path.metaboliteIndices(:));
    if isfield(path,'metabolitePath') && strlength(string(path.metabolitePath)) > 0
        savedMetIDs = strtrim(split(string(path.metabolitePath)," -> "));
        actualMetIDs = string(model.mets(pathMetIdx));
        if numel(savedMetIDs) ~= numel(actualMetIDs) || ...
                any(savedMetIDs(:) ~= actualMetIDs(:))
            error('PrintMotifs:AdapterPathIDMismatch', ...
                'Saved Block 0 metabolite IDs do not match this model ordering.');
        end
    end
    if isfield(path,'reactionPath') && strlength(string(path.reactionPath)) > 0
        savedRxnIDs = strtrim(split(string(path.reactionPath)," -> "));
        savedRxnIDs = regexprep(savedRxnIDs,'\s*\(reverse\)$','');
        actualRxnIDs = string(model.rxns(rxns));
        if numel(savedRxnIDs) ~= numel(actualRxnIDs) || ...
                any(savedRxnIDs(:) ~= actualRxnIDs(:))
            error('PrintMotifs:AdapterPathIDMismatch', ...
                'Saved Block 0 reaction IDs do not match this model ordering.');
        end
    end
end
end

function reportWriteTXT(model,motifTable,InputAdapterPaths,selection, ...
        ADAPTER_MATCH_TO_PRINT,printerFile,outputFile,modelFile,block6File,BLOCK0_REQUIRE_EC)
% Private function workspace prevents stale structs or models from other runs
% reaching Block 7. The standard printer still sees its original row indices.
results_nn_no_crosstalk = table2struct(motifTable); %#ok<NASGU>
MOTIF_SOURCE = "block6"; %#ok<NASGU>
outputFolder = fileparts(outputFile);
if isempty(outputFolder), outputFolder = pwd; end
if ~isfolder(outputFolder)
    error('PrintMotifs:MissingOutputFolder','Output folder not found: %s',outputFolder);
end
temporaryReport = [tempname(outputFolder),'.txt'];
cleanup = onCleanup(@() reportCleanupTXT(temporaryReport));
try
    diary(temporaryReport);
    fprintf('\nSTRICT BLOCK 6 MOTIF REPORT\n');
    fprintf('Model: %s\nResults: %s\n',modelFile,block6File);
    fprintf('EC codes required in Blocks 1-4. Block 0 EC requirement: %d\n',BLOCK0_REQUIRE_EC);
    fprintf('Available motifs: %d\nSelected checkpoint rows: %s\n', ...
        height(motifTable),num2str(selection));
    for reportPosition = 1:numel(selection)
        MOTIF_TO_PRINT = selection(reportPosition); %#ok<NASGU>
        if ismember('motifIndexBeforeECFilter',motifTable.Properties.VariableNames)
            fprintf('\nBlock 6 row before EC filter: %d\n', ...
                motifTable.motifIndexBeforeECFilter(MOTIF_TO_PRINT));
        end
        run(printerFile);
    end
    diary off;
    [ok,message] = movefile(temporaryReport,outputFile,'f');
    if ~ok, error('PrintMotifs:CannotFinalize','%s',message); end
catch err
    diary off;
    rethrow(err);
end
clear cleanup;
fprintf('\nSaved %d selected motif report(s) to:\n%s\n',numel(selection),outputFile);
end

function reportCleanupTXT(filePath)
diary off;
if isfile(filePath), delete(filePath); end
end
