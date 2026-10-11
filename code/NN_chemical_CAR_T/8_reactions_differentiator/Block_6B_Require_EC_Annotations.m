%% Require complete EC annotations for the computational motif
%
% Run AFTER Block 6 and BEFORE Block 7, with the matching model loaded.
% If result variables are absent, this script loads the Block 0/6 MAT files
% from MATLAB's current folder. Existing result files are not modified.
%
% Output in the current folder: Block_6_EC_annotated.mat
% The checkpoint stores one motif table, eligible adapter paths and reports.
% Compatibility structs are created in memory for Block_7_Print_Selected_Motif.
% Blocks 1-4 require EC codes with four numeric components, e.g. 2.3.3.8.
% Block 0 is a passive adapter: BLOCK0_REQUIRE_EC defaults to false.
% Empty annotations, placeholders and partial codes fail in controlled roles.
% No EC numbers are inferred from reaction names or gene associations.
% For independent TXT export, run Print_Block6_Motifs_To_TXT; it loads the
% model and passes it to this filter after loading the result files.

if ~exist('BLOCK0_REQUIRE_EC','var'), BLOCK0_REQUIRE_EC = false; end
if ~exist('model','var') || ~isstruct(model) || ~isscalar(model) || ...
        ~all(isfield(model,{'S','mets','rxns'}))
    error('ECFilter:MissingModel', ...
        ['Load the model before calling this filter in the master, or run ', ...
         'Print_Block6_Motifs_To_TXT for independent model loading and export.']);
end

if exist('Motifs_NN_Block_6','var') && istable(Motifs_NN_Block_6)
    ec6Source = Motifs_NN_Block_6;
elseif exist('results_nn_no_crosstalk','var') && ...
        isstruct(results_nn_no_crosstalk)
    ec6Source = struct2table(results_nn_no_crosstalk(:));
else
    ec6File = fullfile(pwd,'Block_6_no_internal_crosstalk.mat');
    if ~isfile(ec6File)
        error('ECFilter:MissingBlock6','Block 6 MAT file is missing from the current folder.');
    end
    ec6Vars = whos('-file',ec6File);
    if ismember('Motifs_NN_Block_6',{ec6Vars.name})
        ec6Loaded = load(ec6File,'Motifs_NN_Block_6');
        ec6Source = ec6Loaded.Motifs_NN_Block_6;
    elseif ismember('results_nn_no_crosstalk',{ec6Vars.name})
        ec6Loaded = load(ec6File,'results_nn_no_crosstalk');
        ec6Source = struct2table(ec6Loaded.results_nn_no_crosstalk(:));
    else
        error('ECFilter:MissingMotifTable','Block 6 contains no recognized motif table.');
    end
    clear ec6Loaded;
end

if ~exist('InputAdapterPaths','var')
    ec6File = fullfile(pwd,'Block_0_Input_Adapter.mat');
    if ~isfile(ec6File)
        error('ECFilter:MissingBlock0','Block 0 MAT file is missing from the current folder.');
    end
    ec6Loaded = load(ec6File,'InputAdapterPaths');
    if ~isfield(ec6Loaded,'InputAdapterPaths')
        error('ECFilter:MissingAdapterPaths','Block 0 contains no InputAdapterPaths.');
    end
    InputAdapterPaths = ec6Loaded.InputAdapterPaths;
    clear ec6Loaded;
end

[Motifs_NN_Block_6,InputAdapterPaths,EC_Filter_Report, ...
    EC_Adapter_Filter_Report,EC_Reaction_Annotations] = ...
    ec6Filter(model,ec6Source,InputAdapterPaths,BLOCK0_REQUIRE_EC);

% Keep the original ordering and row provenance; prevent Block 7's automatic
% fallback to motifs that have not passed the EC filter.
results_nn_no_crosstalk = table2struct(Motifs_NN_Block_6);
results_nn_no_crosstalk = results_nn_no_crosstalk(:);
fullJoinRows = results_nn_no_crosstalk;
nF = height(Motifs_NN_Block_6);
MOTIF_SOURCE = "block6";

ec6Checkpoint = struct;
ec6Checkpoint.Motifs_NN_Block_6 = Motifs_NN_Block_6;
ec6Checkpoint.InputAdapterPaths = InputAdapterPaths;
ec6Checkpoint.EC_Filter_Report = EC_Filter_Report;
ec6Checkpoint.EC_Adapter_Filter_Report = EC_Adapter_Filter_Report;
ec6Checkpoint.EC_Reaction_Annotations = EC_Reaction_Annotations;
ec6Checkpoint.EC_ModelMetIDs = string(model.mets(:));
ec6Checkpoint.EC_ModelRxnIDs = string(model.rxns(:));
ec6Checkpoint.EC_RequireCompleteCodes = true;
ec6Checkpoint.EC_RequireAdapterCodes = logical(BLOCK0_REQUIRE_EC);
ec6Checkpoint.CheckpointStorageVersion = 1;
ec6Save(fullfile(pwd,'Block_6_EC_annotated.mat'),ec6Checkpoint);

fprintf('\nEC ANNOTATION FILTER\n');
fprintf('Reactions with complete EC codes: %d / %d\n', ...
    sum(EC_Reaction_Annotations.hasCompleteEC),height(EC_Reaction_Annotations));
fprintf('Adapter paths before / retained : %d / %d\n', ...
    height(EC_Adapter_Filter_Report),numel(InputAdapterPaths));
fprintf('Adapter EC annotation required  : %d\n',BLOCK0_REQUIRE_EC);
fprintf('Motifs before / retained        : %d / %d\n',height(ec6Source),nF);
fprintf('Rejected for motif reaction ECs : %d\n', ...
    sum(EC_Filter_Report.nRolesWithoutEC > 0));
fprintf('Rejected for no eligible adapter: %d\n', ...
    sum(EC_Filter_Report.nRolesWithoutEC == 0 & ...
        ~EC_Filter_Report.hasEligibleAdapter));
if nF == 0
    fprintf(['No motif with EC-annotated computational reactions survived. The saved reports ', ...
        'identify the missing annotations; skip Block 7 for this run.\n']);
end
clear ec6Source ec6Checkpoint ec6File ec6Loaded ec6Vars;

%% Local functions
function [motifs,paths,motifReport,adapterReport,rxnReport] = ec6Filter(model,motifs,paths,requireAdapterEC)
[nMets,nRxns] = size(model.S);
if numel(model.mets) ~= nMets || numel(model.rxns) ~= nRxns
    error('ECFilter:ModelSizeMismatch','model.S, mets and rxns have inconsistent sizes.');
end
if ~istable(motifs) || ~isstruct(paths)
    error('ECFilter:InvalidResults','Expected a motif table and an adapter struct array.');
end
roles = {'r_S1','r_Sg1','r_Ss1','r_F1','r_Fg1','r_Fs1','r_I1', ...
    'r_S2','r_Sg2','r_Ss2','r_F2','r_Fg2','r_Fs2','r_I2', ...
    'r_C','r_CI1','r_CI2','r_omega','r_b'};
if height(motifs) > 0
    required = [roles,{'m_in'}];
    missing = setdiff(required,motifs.Properties.VariableNames);
    if ~isempty(missing)
        error('ECFilter:MissingRoles','Missing motif fields: %s.',strjoin(missing,', '));
    end
end

[hasEC,codes,annotationFields] = getReactionECAnnotations(model);
rxnReport = table((1:nRxns).',string(model.rxns(:)),hasEC,codes,annotationFields, ...
    'VariableNames',{'reactionIndex','reactionID','hasCompleteEC', ...
    'completeECNumbers','annotationFields'});

n = height(motifs);
missingRoles = zeros(n,1,'uint16');
for k = 1:numel(roles)
    if n == 0, break; end
    idx = motifs.(roles{k});
    ec6CheckIndices(idx,nRxns,n,roles{k});
    missingRoles = missingRoles + uint16(~hasEC(double(idx(:))));
end
metRoles = {'m_in','m_S1','m_F1','m_Sg1','m_Fg1','m_S2','m_F2', ...
    'm_Sg2','m_Fg2','m_d','m_c','m_Cs','m_out'};
for k = 1:numel(metRoles)
    if n > 0 && ismember(metRoles{k},motifs.Properties.VariableNames)
        ec6CheckIndices(motifs.(metRoles{k}),nMets,n,metRoles{k});
    end
end

paths = paths(:);
np = numel(paths);
pathAnchor = zeros(np,1);
pathNumber = (1:np).';
missingPathRxns = zeros(np,1);
pathHasEC = true(np,1);
for p = 1:np
    path = paths(p);
    if ~all(isfield(path,{'anchorMetIdx','reactionIndices'}))
        error('ECFilter:IncompleteAdapter','Adapter %d lacks its anchor or reaction indices.',p);
    end
    ec6CheckIndices(path.anchorMetIdx,nMets,1,'adapter.anchorMetIdx');
    pathAnchor(p) = double(path.anchorMetIdx);
    idx = path.reactionIndices;
    if ~isempty(idx)
        ec6CheckIndices(idx(:),nRxns,numel(idx),'adapter.reactionIndices');
    end
    % These saved IDs guard against using the yeast model for CAR-T paths,
    % or a differently ordered model with the same dimensions.
    if isfield(path,'anchorMetID') && ...
            string(path.anchorMetID) ~= string(model.mets(pathAnchor(p)))
        error('ECFilter:ModelMismatch','Adapter %d anchor ID disagrees with this model.',p);
    end
    if isfield(path,'reactionPath') && strlength(string(path.reactionPath)) > 0
        savedIDs = strtrim(split(string(path.reactionPath),"->"));
        savedIDs = regexprep(savedIDs,'\s*\(reverse\)$','');
        if ~isequal(savedIDs(:),string(model.rxns(double(idx(:)))))
            error('ECFilter:ModelMismatch','Adapter %d reaction IDs disagree with this model.',p);
        end
    end
    if isfield(path,'pathIndex') && ~isempty(path.pathIndex)
        pathNumber(p) = double(path.pathIndex);
    end
    missingPathRxns(p) = sum(~hasEC(double(idx(:))));
    pathHasEC(p) = missingPathRxns(p) == 0;
end
pathEligible = pathHasEC | ~requireAdapterEC;

% Match exactly as Block 7 does: Block 0 endpoint equals motif m_in.
% A motif needs at least one eligible path; a different rejected path to the
% same endpoint does not invalidate that eligible connection.
hasAdapter = false(n,1);
if n > 0
    hasAdapter = ismember(double(motifs.m_in(:)),pathAnchor(pathEligible));
    if ismember('adapterPathIndex',motifs.Properties.VariableNames)
        ec6CheckIndices(motifs.adapterPathIndex,Inf,n,'adapterPathIndex');
        hasAdapter = ismember( ...
            [double(motifs.m_in(:)),double(motifs.adapterPathIndex(:))], ...
            [pathAnchor(pathEligible),pathNumber(pathEligible)],'rows');
    end
end
keep = missingRoles == 0 & hasAdapter;
motifReport = table((1:n).',keep,missingRoles,hasAdapter, ...
    'VariableNames',{'motifIndexBeforeECFilter','pass', ...
    'nRolesWithoutEC','hasEligibleAdapter'});

if ~ismember('motifIndexBeforeECFilter',motifs.Properties.VariableNames)
    motifs.motifIndexBeforeECFilter = (1:n).';
end
motifs = motifs(keep,:);
keepPath = false(np,1);
if height(motifs) > 0
    keepPath = pathEligible & ismember(pathAnchor,double(motifs.m_in(:)));
    if ismember('adapterPathIndex',motifs.Properties.VariableNames)
        keepPath = pathEligible & ismember([pathAnchor,pathNumber], ...
            [double(motifs.m_in(:)),double(motifs.adapterPathIndex(:))],'rows');
    end
end
adapterReport = table((1:np).',pathNumber,pathAnchor,pathHasEC, ...
    missingPathRxns,pathEligible,keepPath, ...
    'VariableNames',{'adapterArrayIndexBeforeECFilter','pathIndex', ...
    'anchorMetIdx','hasCompleteEC','nReactionsWithoutEC','eligible','retained'});
paths = paths(keepPath);
end

function ec6CheckIndices(idx,limit,count,role)
if ~isnumeric(idx) || ~isreal(idx) || numel(idx) ~= count || ...
        (count > 0 && ~isvector(idx)) || ...
        any(~isfinite(idx(:)) | idx(:) < 1 | idx(:) > limit | ...
        idx(:) ~= floor(idx(:)))
    error('ECFilter:InvalidIndex', ...
        '%s contains invalid indices. Check that the model matches the results.',role);
end
end

function ec6Save(filePath,payload)
temporaryFile = [tempname(pwd),'.mat'];
cleanup = onCleanup(@() ec6DeleteTemporary(temporaryFile));
save(temporaryFile,'-struct','payload','-v7');
[ok,message] = movefile(temporaryFile,filePath,'f');
if ~ok, error('ECFilter:SaveFailed','%s',message); end
clear cleanup;
fprintf('\nSaved EC-filtered checkpoint: %s\n',filePath);
end

function ec6DeleteTemporary(filePath)
if isfile(filePath), delete(filePath); end
end
