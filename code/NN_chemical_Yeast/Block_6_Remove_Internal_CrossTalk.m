%% Block 6: remove motifs containing internal pathway cross-talk
%
% This filter starts from the motifs retained by Block 5 and checks only
% consumption edges AMONG the 19 selected reactions of each motif. It does
% not reject a motif merely because another, unselected T-cell reaction can
% consume one of its metabolites.
%
% A protected signal metabolite may be consumed only by its intended
% reaction role(s). Any additional consumption edge from another selected
% reaction is classified as internal cross-talk.
%
% Expected master-file sequence:
%   load Block_5_no_cell_membrane.mat
%   Block_6_Remove_Internal_CrossTalk
%
% Saved output:
%   Block_6_no_internal_crosstalk.mat
%       Motifs_NN_Block_6
%       results_nn_no_crosstalk
%       CrossTalk_Filter_Report_Block_6
%       indexBeforeCrossTalkFilter
%       keepCrossTalkFree

%% Display controls
if ~exist('SHOW_CROSSTALK_FILTER_PROGRESS','var')
    SHOW_CROSSTALK_FILTER_PROGRESS = true;
end

if ~exist('CROSSTALK_PROGRESS_STEP_PERCENT','var')
    CROSSTALK_PROGRESS_STEP_PERCENT = 10;
end
CROSSTALK_PROGRESS_STEP_PERCENT = max(1, ...
    min(100,round(CROSSTALK_PROGRESS_STEP_PERCENT)));

if ~exist('PRINT_REJECTED_CROSSTALK_MOTIFS','var')
    PRINT_REJECTED_CROSSTALK_MOTIFS = false;
end

if ~exist('MAX_CROSSTALK_REJECTED_TO_PRINT','var')
    MAX_CROSSTALK_REJECTED_TO_PRINT = 20;
end
MAX_CROSSTALK_REJECTED_TO_PRINT = max(0, ...
    floor(MAX_CROSSTALK_REJECTED_TO_PRINT));

%% Validate the model
if ~exist('model','var') || ~isfield(model,'S') || ...
        ~isfield(model,'mets') || ~isfield(model,'rxns')
    error('Block6:MissingModel', ...
        'model.S, model.mets and model.rxns are required.');
end

[nMets,nRxns] = size(model.S);
if numel(model.mets) ~= nMets || numel(model.rxns) ~= nRxns
    error('Block6:ModelSizeMismatch', ...
        'model.mets/model.rxns do not agree with model.S.');
end

metIds = string(model.mets(:));
rxnIds = string(model.rxns(:));

if isfield(model,'metNames') && numel(model.metNames) == nMets
    metNames = string(model.metNames(:));
else
    metNames = metIds;
end

if isfield(model,'rxnNames') && numel(model.rxnNames) == nRxns
    rxnNames = string(model.rxnNames(:));
else
    rxnNames = rxnIds;
end

%% Read the Block 5 motifs without changing their order
if exist('Motifs_NN_Block_5','var') && istable(Motifs_NN_Block_5)
    sourceMotifs = Motifs_NN_Block_5;
elseif exist('results_nn_no_membrane','var') && ...
        isstruct(results_nn_no_membrane)
    sourceMotifs = struct2table(results_nn_no_membrane(:));
elseif exist('fullJoinRows','var') && isstruct(fullJoinRows)
    sourceMotifs = struct2table(fullJoinRows(:));
else
    error('Block6:MissingBlock5', ...
        ['Load Block_5_no_cell_membrane.mat before running Block 6. ', ...
         'Expected Motifs_NN_Block_5 or results_nn_no_membrane.']);
end

nMotifs = height(sourceMotifs);
if nMotifs == 0
    error('Block6:EmptyBlock5', ...
        'Block 5 contains no motifs to test for cross-talk.');
end

%% Canonical reaction roles in the current complete motif
rxnFields = { ...
    'r_S1','r_Sg1','r_Ss1','r_F1','r_Fg1','r_Fs1','r_I1', ...
    'r_S2','r_Sg2','r_Ss2','r_F2','r_Fg2','r_Fs2','r_I2', ...
    'r_C','r_CI1','r_CI2', ...
    'r_omega','r_b'};

%% Protected signal metabolites and their intended consumers
%
% First differentiator:
%   m_in --r_S1--> m_S1 --r_Sg1--> m_Sg1
%                         \--r_Ss1--> ignored
%        --r_F1--> m_F1 --r_Fg1--> m_Fg1
%                         \--r_Fs1--> ignored
%   m_Sg1 + m_Fg1 --r_I1--> ignored
%
% Second differentiator:
%   m_Sg1 --r_S2--> m_S2 --r_Sg2--> m_Sg2
%                           \--r_Ss2--> ignored
%   m_Fg1 --r_F2--> m_F2 --r_Fg2--> m_d
%                           \--r_Fs2--> ignored
%   m_Sg2 + m_d --r_I2--> ignored
%
% Counter and neuron:
%   m_d --r_C--> m_c
%   m_c + m_Cs --r_CI1--> ignored
%   m_d + m_Cs --r_CI2--> ignored
%   m_c --r_omega--> m_out
%   buffered substrates --r_b--> m_out

metFields = { ...
    'm_in', ...
    'm_S1','m_F1','m_Sg1','m_Fg1', ...
    'm_S2','m_F2','m_Sg2','m_d', ...
    'm_c','m_Cs','m_out'};

allowedConsumers = { ...
    {'r_S1','r_F1'}, ...                 % m_in
    {'r_Sg1','r_Ss1'}, ...              % m_S1
    {'r_Fg1','r_Fs1'}, ...              % m_F1
    {'r_I1','r_S2'}, ...                % m_Sg1
    {'r_I1','r_F2'}, ...                % m_Fg1
    {'r_Sg2','r_Ss2'}, ...              % m_S2
    {'r_Fg2','r_Fs2'}, ...              % m_F2
    {'r_I2'}, ...                        % m_Sg2
    {'r_I2','r_C','r_CI2'}, ...         % m_d (= m_Fg2)
    {'r_CI1','r_omega'}, ...            % m_c (= m_1)
    {'r_CI1','r_CI2'}, ...              % m_Cs (= m_s)
    cell(1,0)};                          % m_out: terminal output

sourceVariables = sourceMotifs.Properties.VariableNames;
missingRxnFields = rxnFields(~ismember(rxnFields,sourceVariables));
missingMetFields = metFields(~ismember(metFields,sourceVariables));

if ~isempty(missingRxnFields)
    error('Block6:MissingReactionRoles', ...
        'Block 5 motifs are missing reaction field(s): %s.', ...
        strjoin(missingRxnFields,', '));
end
if ~isempty(missingMetFields)
    error('Block6:MissingMetaboliteRoles', ...
        'Block 5 motifs are missing metabolite field(s): %s.', ...
        strjoin(missingMetFields,', '));
end

%% Assemble and validate the reaction/metabolite index matrices
nRxnRoles = numel(rxnFields);
nMetRoles = numel(metFields);
motifRxnIdx = zeros(nMotifs,nRxnRoles);
motifMetIdx = zeros(nMotifs,nMetRoles);

for f = 1:nRxnRoles
    values = sourceMotifs.(rxnFields{f});
    if ~isnumeric(values) || ~isvector(values) || numel(values) ~= nMotifs
        error('Block6:InvalidReactionColumn', ...
            'Column %s must contain one numeric index per motif.', ...
            rxnFields{f});
    end
    motifRxnIdx(:,f) = double(values(:));
end

for f = 1:nMetRoles
    values = sourceMotifs.(metFields{f});
    if ~isnumeric(values) || ~isvector(values) || numel(values) ~= nMotifs
        error('Block6:InvalidMetaboliteColumn', ...
            'Column %s must contain one numeric index per motif.', ...
            metFields{f});
    end
    motifMetIdx(:,f) = double(values(:));
end

badRxnEntry = ~isfinite(motifRxnIdx) | motifRxnIdx < 1 | ...
    motifRxnIdx > nRxns | motifRxnIdx ~= floor(motifRxnIdx);
badMetEntry = ~isfinite(motifMetIdx) | motifMetIdx < 1 | ...
    motifMetIdx > nMets | motifMetIdx ~= floor(motifMetIdx);

if any(badRxnEntry(:))
    [row,field] = find(badRxnEntry,1,'first');
    error('Block6:InvalidReactionIndex', ...
        'Motif %d contains invalid index %g in %s.', ...
        row,motifRxnIdx(row,field),rxnFields{field});
end
if any(badMetEntry(:))
    [row,field] = find(badMetEntry,1,'first');
    error('Block6:InvalidMetaboliteIndex', ...
        'Motif %d contains invalid index %g in %s.', ...
        row,motifMetIdx(row,field),metFields{field});
end

allowedRoleIdx = cell(nMetRoles,1);
for m = 1:nMetRoles
    [found,roleIdx] = ismember(allowedConsumers{m},rxnFields);
    if any(~found)
        error('Block6:InvalidAllowedConsumer', ...
            'Unknown allowed reaction role for %s: %s.', ...
            metFields{m},strjoin(allowedConsumers{m}(~found),', '));
    end
    allowedRoleIdx{m} = roleIdx;
end

%% Test every possible protected-metabolite/selected-reaction edge
nCrossTalkEdges = zeros(nMotifs,1);
badMetRoleMask = false(nMotifs,nMetRoles);
crossTalkReason = strings(nMotifs,1);

nRolePairs = nMetRoles*nRxnRoles;
rolePairsDone = 0;
nextProgress = CROSSTALK_PROGRESS_STEP_PERCENT;

if SHOW_CROSSTALK_FILTER_PROGRESS
    fprintf('\nBlock 6 internal cross-talk scan:   0%% (0/%d role pairs)', ...
        nRolePairs);
end

for m = 1:nMetRoles
    metIdx = motifMetIdx(:,m);
    allowedIdx = allowedRoleIdx{m};

    if isempty(allowedIdx)
        allowedRxnIdx = zeros(nMotifs,0);
    else
        allowedRxnIdx = motifRxnIdx(:,allowedIdx);
    end

    for r = 1:nRxnRoles
        rxnIdx = motifRxnIdx(:,r);
        linearIdx = sub2ind([nMets,nRxns],metIdx,rxnIdx);
        consumesMet = full(model.S(linearIdx)) < 0;

        if isempty(allowedRxnIdx)
            isAllowed = false(nMotifs,1);
        else
            % Match actual reaction IDs, retaining the behavior of the
            % original filter when two role names happen to alias one ID.
            isAllowed = any(rxnIdx == allowedRxnIdx,2);
        end

        badRows = find(consumesMet & ~isAllowed);
        if ~isempty(badRows)
            nCrossTalkEdges(badRows) = nCrossTalkEdges(badRows) + 1;
            badMetRoleMask(badRows,m) = true;

            for q = 1:numel(badRows)
                row = badRows(q);
                mIdx = metIdx(row);
                rIdx = rxnIdx(row);
                detail = sprintf( ...
                    '%s [%s | %s] is also consumed by %s [%s | %s]', ...
                    metFields{m},metIds(mIdx),metNames(mIdx), ...
                    rxnFields{r},rxnIds(rIdx),rxnNames(rIdx));

                if strlength(crossTalkReason(row)) == 0
                    crossTalkReason(row) = detail;
                else
                    crossTalkReason(row) = ...
                        crossTalkReason(row) + "; " + detail;
                end
            end
        end

        rolePairsDone = rolePairsDone + 1;
        if SHOW_CROSSTALK_FILTER_PROGRESS
            percent = floor(100*rolePairsDone/nRolePairs);
            if percent >= nextProgress || rolePairsDone == nRolePairs
                fprintf('\rBlock 6 internal cross-talk scan: %3d%% (%d/%d role pairs)', ...
                    percent,rolePairsDone,nRolePairs);
                while nextProgress <= percent
                    nextProgress = nextProgress + ...
                        CROSSTALK_PROGRESS_STEP_PERCENT;
                end
                drawnow;
            end
        end
    end
end

if SHOW_CROSSTALK_FILTER_PROGRESS
    fprintf('\n');
end

keepCrossTalkFree = nCrossTalkEdges == 0;
nCrossTalkMets = zeros(nMotifs,1);
rejectedIdx = find(~keepCrossTalkFree);

% Count unique model metabolites, rather than merely distinct role labels.
for q = 1:numel(rejectedIdx)
    row = rejectedIdx(q);
    nCrossTalkMets(row) = numel(unique( ...
        motifMetIdx(row,badMetRoleMask(row,:)),'stable'));
end

crossTalkReason(keepCrossTalkFree) = ...
    "(no forbidden internal consumption edges)";

CrossTalk_Filter_Report_Block_6 = table( ...
    (1:nMotifs).',keepCrossTalkFree,nCrossTalkMets, ...
    nCrossTalkEdges,crossTalkReason, ...
    'VariableNames',{ ...
        'motifIndexBeforeCrossTalkFilter','pass', ...
        'nCrossTalkMetabolites','nCrossTalkEdges','reason'});

%% Retain cross-talk-free motifs in their existing order
indexBeforeCrossTalkFilter = find(keepCrossTalkFree);
Motifs_NN_Block_6 = sourceMotifs(keepCrossTalkFree,:);
Motifs_NN_Block_6.motifIndexBeforeCrossTalkFilter = ...
    indexBeforeCrossTalkFilter;

results_nn_no_crosstalk = table2struct(Motifs_NN_Block_6);
results_nn_no_crosstalk = results_nn_no_crosstalk(:);

% Compatibility variables for later printing and filtering code.
fullJoinRows_beforeCrossTalkFilter = table2struct(sourceMotifs);
fullJoinRows_beforeCrossTalkFilter = ...
    fullJoinRows_beforeCrossTalkFilter(:);
fullJoinRows = results_nn_no_crosstalk;
nF = numel(fullJoinRows);

%% Save
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end
outputFile = fullfile(scriptFolder, ...
    'Block_6_no_internal_crosstalk.mat');

save(outputFile, ...
    'Motifs_NN_Block_6','results_nn_no_crosstalk', ...
    'CrossTalk_Filter_Report_Block_6', ...
    'fullJoinRows_beforeCrossTalkFilter', ...
    'indexBeforeCrossTalkFilter','keepCrossTalkFree', ...
    'nCrossTalkMets','nCrossTalkEdges','crossTalkReason', ...
    'rxnFields','metFields','allowedConsumers','-v7.3');

%% Summary
fprintf('\n############################################################\n');
fprintf('BLOCK 6 INTERNAL CROSS-TALK FILTER\n');
fprintf('############################################################\n');
fprintf('Motifs tested   : %d\n',nMotifs);
fprintf('Motifs retained : %d\n',sum(keepCrossTalkFree));
fprintf('Motifs removed  : %d\n',sum(~keepCrossTalkFree));
fprintf('Reaction roles  : %d\n',nRxnRoles);
fprintf('Protected mets  : %d\n',nMetRoles);
fprintf('############################################################\n');

if PRINT_REJECTED_CROSSTALK_MOTIFS && ~isempty(rejectedIdx)
    nToPrint = min(numel(rejectedIdx),MAX_CROSSTALK_REJECTED_TO_PRINT);
    fprintf('\nFIRST %d CROSS-TALK REJECTIONS\n',nToPrint);
    for q = 1:nToPrint
        row = rejectedIdx(q);
        fprintf('\nMotif %d | metabolites: %d | edges: %d\n', ...
            row,nCrossTalkMets(row),nCrossTalkEdges(row));
        reasons = split(crossTalkReason(row),"; ");
        for k = 1:numel(reasons)
            fprintf('  %d. %s\n',k,reasons(k));
        end
    end
end

fprintf('\n%d motifs remain after internal cross-talk filtering.\n',nF);
fprintf('Saved Block 6 results to:\n%s\n',outputFile);
