%% Block 7: print one selected cAMP-to-NN motif
%
% Set MOTIF_TO_PRINT to one row of the selected result source. The script prints:
%   A) the first Block 0 path connecting cAMP to this motif's m_in;
%   B) the selected 19-reaction Block 1-4 motif.
%
% MOTIF_SOURCE controls which filtered result is printed:
%   "auto"   : use Block 6 when it has survivors; otherwise use Block 5
%   "block6" : require a strict cross-talk-free Block 6 survivor
%   "block5" : use a membrane-filtered Block 5 motif (cross-talk unchecked)
%
% The model must already be loaded. If result variables are not in the
% workspace, this script loads the corresponding MAT files from its own
% directory or the current directory.

%% Selection
if ~exist('MOTIF_TO_PRINT','var') || isempty(MOTIF_TO_PRINT)
    MOTIF_TO_PRINT = 1;
end

% If several Block 0 paths reach the selected m_in, print the first one.
if ~exist('ADAPTER_MATCH_TO_PRINT','var') || isempty(ADAPTER_MATCH_TO_PRINT)
    ADAPTER_MATCH_TO_PRINT = 1;
end

% Recommended: prefer strict Block 6 results, but still allow inspection of
% Block 5 motifs when the strict cross-talk filter retained no rows.
if ~exist('MOTIF_SOURCE','var') || isempty(MOTIF_SOURCE)
    MOTIF_SOURCE = "auto";
end

motifSourceRequest = lower(strtrim(string(MOTIF_SOURCE)));
if ~isscalar(motifSourceRequest) || ...
        ~ismember(motifSourceRequest,["auto","block6","block5"])
    error('Printer:InvalidMotifSource', ...
        'MOTIF_SOURCE must be "auto", "block6", or "block5".');
end

%% Validate model
if ~exist('model','var') || ~isfield(model,'S') || ...
        ~isfield(model,'mets') || ~isfield(model,'rxns')
    error('Printer:MissingModel', ...
        'Load the same COBRA model before printing a motif.');
end

[nMets,nRxns] = size(model.S);
[~,printerECValues] = getReactionECAnnotations(model);

%% Locate result files
printerFolder = fileparts(mfilename('fullpath'));
if isempty(printerFolder)
    printerFolder = pwd;
end

%% Obtain rows from the requested filter stage
printableRows = struct([]);
printableSource = "";
fellBackFromBlock6 = false;

% First try the strict Block 6 output unless Block 5 was explicitly chosen.
if motifSourceRequest ~= "block5"
    if exist('results_nn_no_crosstalk','var') && ...
            isstruct(results_nn_no_crosstalk)
        block6Rows = results_nn_no_crosstalk(:);
    elseif exist('Motifs_NN_Block_6','var') && istable(Motifs_NN_Block_6)
        block6Rows = table2struct(Motifs_NN_Block_6);
        block6Rows = block6Rows(:);
    else
        block6Rows = struct([]);
        resultFile = localFindFile( ...
            'Block_6_no_internal_crosstalk.mat',printerFolder);
        if ~isempty(resultFile)
            B6 = load(resultFile);
            if isfield(B6,'results_nn_no_crosstalk') && ...
                    isstruct(B6.results_nn_no_crosstalk)
                block6Rows = B6.results_nn_no_crosstalk(:);
            elseif isfield(B6,'Motifs_NN_Block_6') && ...
                    istable(B6.Motifs_NN_Block_6)
                block6Rows = table2struct(B6.Motifs_NN_Block_6);
                block6Rows = block6Rows(:);
            end
        end
    end

    if ~isempty(block6Rows)
        printableRows = block6Rows;
        printableSource = "Block 6: strict cross-talk-filtered motifs";
    elseif motifSourceRequest == "block6"
        error('Printer:NoBlock6Motifs', ...
            ['Block 6 retained zero motifs. Choose MOTIF_SOURCE = ', ...
             '"block5" to inspect membrane-filtered motifs, or relax ', ...
             'the Block 6 cross-talk rules and rerun Block 6.']);
    else
        fellBackFromBlock6 = true;
    end
end

% Explicit Block 5 choice, or automatic fallback when Block 6 is empty.
if isempty(printableRows) && ...
        (motifSourceRequest == "block5" || fellBackFromBlock6)
    if exist('results_nn_no_membrane','var') && ...
            isstruct(results_nn_no_membrane)
        block5Rows = results_nn_no_membrane(:);
    elseif exist('Motifs_NN_Block_5','var') && istable(Motifs_NN_Block_5)
        block5Rows = table2struct(Motifs_NN_Block_5);
        block5Rows = block5Rows(:);
    else
        block5Rows = struct([]);
        resultFile = localFindFile('Block_5_no_cell_membrane.mat', ...
            printerFolder);
        if ~isempty(resultFile)
            B5 = load(resultFile);
            if isfield(B5,'results_nn_no_membrane') && ...
                    isstruct(B5.results_nn_no_membrane)
                block5Rows = B5.results_nn_no_membrane(:);
            elseif isfield(B5,'Motifs_NN_Block_5') && ...
                    istable(B5.Motifs_NN_Block_5)
                block5Rows = table2struct(B5.Motifs_NN_Block_5);
                block5Rows = block5Rows(:);
            end
        end
    end

    if ~isempty(block5Rows)
        printableRows = block5Rows;
        printableSource = ...
            "Block 5: membrane-filtered motifs (cross-talk not removed)";
        if fellBackFromBlock6
            warning('Printer:Block5Fallback', ...
                ['Block 6 retained zero motifs. Printing a Block 5 motif ', ...
                 'instead. This motif is NOT guaranteed to be free of ', ...
                 'internal cross-talk.']);
        end
    end
end

nAvailableMotifs = numel(printableRows);
if nAvailableMotifs == 0
    error('Printer:NoMotifs', ...
        ['No printable motif was found in the selected source. Confirm ', ...
         'that Block 5 produced survivors and saved ', ...
         'Block_5_no_cell_membrane.mat.']);
end

if ~isscalar(MOTIF_TO_PRINT) || ~isfinite(MOTIF_TO_PRINT) || ...
        MOTIF_TO_PRINT < 1 || MOTIF_TO_PRINT ~= floor(MOTIF_TO_PRINT) || ...
        MOTIF_TO_PRINT > nAvailableMotifs
    error('Printer:InvalidSelection', ...
        'MOTIF_TO_PRINT must be one integer from 1 to %d.', ...
        nAvailableMotifs);
end

rec = printableRows(MOTIF_TO_PRINT);

%% Obtain the Block 0 adapter paths
if ~exist('InputAdapterPaths','var') || ~isstruct(InputAdapterPaths)
    block0File = localFindFile('Block_0_Input_Adapter.mat',printerFolder);
    if isempty(block0File)
        InputAdapterPaths = struct([]);
        warning('Printer:MissingBlock0Results', ...
            ['Block_0_Input_Adapter.mat was not found. The selected ', ...
             'Block 1-4 motif will be printed without its cAMP adapter.']);
    else
        B0 = load(block0File,'InputAdapterPaths');
        if isfield(B0,'InputAdapterPaths')
            InputAdapterPaths = B0.InputAdapterPaths(:);
        else
            InputAdapterPaths = struct([]);
        end
    end
end

%% Current metabolite and reaction roles
metFields = { ...
    'm_in', ...
    'm_S1','m_F1','m_Sg1','m_Fg1', ...
    'm_S2','m_F2','m_Sg2','m_d', ...
    'm_c','m_Cs','m_out'};

metLabels = { ...
    'm_in   (Block 0 endpoint / differentiator 1 input)', ...
    'm_S1   (differentiator 1 slow tracker)', ...
    'm_F1   (differentiator 1 fast tracker)', ...
    'm_Sg1  (differentiator 1 upper reporter / Block 2 input)', ...
    'm_Fg1  (differentiator 1 lower reporter / Block 2 input)', ...
    'm_S2   (differentiator 2 upper tracker)', ...
    'm_F2   (differentiator 2 lower tracker)', ...
    'm_Sg2  (differentiator 2 upper reporter)', ...
    'm_d    (m_Fg2: differentiator 2 bottom output / counter input)', ...
    'm_c    (counter output / NN input)', ...
    'm_Cs   (shared counter species)', ...
    'm_out  (single-neuron output)'};

rxnFields = { ...
    'r_S1','r_Sg1','r_Ss1','r_F1','r_Fg1','r_Fs1','r_I1', ...
    'r_S2','r_Sg2','r_Ss2','r_F2','r_Fg2','r_Fs2','r_I2', ...
    'r_C','r_CI1','r_CI2','r_omega','r_b'};

rxnLabels = { ...
    'r_S1     (differentiator 1: upper production)', ...
    'r_Sg1    (differentiator 1: upper reporter)', ...
    'r_Ss1    (differentiator 1: upper removal)', ...
    'r_F1     (differentiator 1: lower production)', ...
    'r_Fg1    (differentiator 1: lower reporter)', ...
    'r_Fs1    (differentiator 1: lower removal)', ...
    'r_I1     (differentiator 1: joint consumption)', ...
    'r_S2     (differentiator 2: upper production)', ...
    'r_Sg2    (differentiator 2: upper reporter)', ...
    'r_Ss2    (differentiator 2: upper removal)', ...
    'r_F2     (differentiator 2: lower production)', ...
    'r_Fg2    (differentiator 2: bottom reporter)', ...
    'r_Fs2    (differentiator 2: lower removal)', ...
    'r_I2     (differentiator 2: joint consumption)', ...
    'r_C      (counter: m_d to m_c)', ...
    'r_CI1    (counter: consumes m_c and m_Cs)', ...
    'r_CI2    (counter: consumes m_d and m_Cs)', ...
    'r_omega  (single neuron: weighted input)', ...
    'r_b      (single neuron: constant bias)'};

requiredFields = [metFields,rxnFields];
missingFields = requiredFields(~isfield(rec,requiredFields));
if ~isempty(missingFields)
    error('Printer:IncompleteMotif', ...
        'The selected motif is missing field(s): %s.', ...
        strjoin(missingFields,', '));
end

%% ================================================================
%% A) Print the first matching cAMP -> m_in adapter path
%% ================================================================

fprintf('\n\n############################################################\n');
fprintf('A) cAMP INPUT ADAPTER FOR SELECTED MOTIF %d\n',MOTIF_TO_PRINT);
fprintf('############################################################\n');

adapterPrinted = false;
if ~isempty(InputAdapterPaths)
    adapterAnchors = double([InputAdapterPaths.anchorMetIdx]).';
    matchingAdapterPaths = find(adapterAnchors == double(rec.m_in));

    if isempty(matchingAdapterPaths)
        warning('Printer:NoMatchingAdapter', ...
            ['No Block 0 path ends at the selected motif input %s. ', ...
             'The main motif will still be printed.'], ...
            localMetId(model,rec.m_in));
    else
        if ~isscalar(ADAPTER_MATCH_TO_PRINT) || ...
                ADAPTER_MATCH_TO_PRINT < 1 || ...
                ADAPTER_MATCH_TO_PRINT ~= floor(ADAPTER_MATCH_TO_PRINT) || ...
                ADAPTER_MATCH_TO_PRINT > numel(matchingAdapterPaths)
            error('Printer:InvalidAdapterSelection', ...
                ['ADAPTER_MATCH_TO_PRINT must be one integer from 1 to ', ...
                 '%d for the selected motif input.'], ...
                numel(matchingAdapterPaths));
        end

        adapterArrayIdx = matchingAdapterPaths(ADAPTER_MATCH_TO_PRINT);
        adapter = InputAdapterPaths(adapterArrayIdx);
        adapterPrinted = true;

        fprintf('Selected final motif row       : %d / %d\n', ...
            MOTIF_TO_PRINT,nAvailableMotifs);
        fprintf('Matching Block 0 paths to m_in : %d\n', ...
            numel(matchingAdapterPaths));
        fprintf('Printed matching path          : %d\n', ...
            ADAPTER_MATCH_TO_PRINT);
        if isfield(adapter,'pathIndex') && ~isempty(adapter.pathIndex)
            fprintf('Stored Block 0 path index      : %d\n', ...
                adapter.pathIndex);
        end

        fprintf('\nInput        : %s | %s\n', ...
            localMetId(model,adapter.inputMetIdx), ...
            localMetName(model,adapter.inputMetIdx));
        fprintf('Fixed m_in   : %s | %s\n', ...
            localMetId(model,rec.m_in),localMetName(model,rec.m_in));

        if isfield(adapter,'metabolitePath')
            fprintf('Metabolites  : %s\n',string(adapter.metabolitePath));
        end
        if isfield(adapter,'reactionPath')
            fprintf('Reactions    : %s\n',string(adapter.reactionPath));
        end

        if isfield(adapter,'metaboliteIndices') && ...
                isfield(adapter,'reactionIndices')
            adapterMets = double(adapter.metaboliteIndices(:)).';
            adapterRxns = double(adapter.reactionIndices(:)).';
            if isfield(adapter,'reactionDirections') && ...
                    numel(adapter.reactionDirections) == numel(adapterRxns)
                adapterDirections = ...
                    double(adapter.reactionDirections(:)).';
            else
                adapterDirections = ones(size(adapterRxns));
            end

            adapterLabelMap = containers.Map( ...
                'KeyType','double','ValueType','char');
            for q = 1:numel(adapterMets)
                if q == 1
                    tag = 'cAMP input';
                elseif q == numel(adapterMets)
                    tag = 'm_in';
                else
                    tag = sprintf('adapter_%d',q-1);
                end
                adapterLabelMap(adapterMets(q)) = tag;
            end

            fprintf('\n--- Block 0 reactions ---\n');
            for q = 1:numel(adapterRxns)
                rIdx = adapterRxns(q);
                direction = adapterDirections(q);
                if direction < 0
                    directionText = 'reverse';
                else
                    directionText = 'forward';
                end

                fprintf('\n  Adapter reaction %d (%s): %s | %s\n', ...
                    q,directionText,localRxnId(model,rIdx), ...
                    localRxnName(model,rIdx));
                fprintf('      GPR: %s\n',localGPR(model,rIdx));
                fprintf('      EC : %s\n',localEC(rIdx,printerECValues));
                fprintf('      Compartment(s): %s\n', ...
                    localReactionCompartments(model,rIdx));
                localPrintReaction(model,rIdx,adapterLabelMap,direction);
            end
        end
    end
end

if ~adapterPrinted
    fprintf('No Block 0 adapter path was printed.\n');
end

%% ================================================================
%% B) Print the selected Block 1-4 motif
%% ================================================================

fprintf('\n\n############################################################\n');
fprintf('B) SELECTED COMPLETE MOTIF %d OF %d\n', ...
    MOTIF_TO_PRINT,nAvailableMotifs);
fprintf('Source: %s\n',printableSource);
if fellBackFromBlock6
    fprintf(['WARNING: Block 6 had no survivors; this Block 5 motif ', ...
             'may contain internal cross-talk.\n']);
end
fprintf('############################################################\n');

provenanceFields = { ...
    'motifIndexBeforeCrossTalkFilter', ...
    'motifIndexBeforeMembraneFilter', ...
    'block3_idx','block2_idx','block1_idx', ...
    'nn_candidate_idx','counter_input_idx','input_pair_idx'};

for q = 1:numel(provenanceFields)
    field = provenanceFields{q};
    if isfield(rec,field)
        value = rec.(field);
        if isnumeric(value) && isscalar(value) && isfinite(value)
            fprintf('%-38s : %g\n',field,value);
        end
    end
end

fprintf('\n--- Metabolites ---\n');
for f = 1:numel(metFields)
    field = metFields{f};
    mIdx = double(rec.(field));
    localValidateMetIndex(mIdx,nMets,field);
    fprintf('  %-64s | %-18s | %-40s | %s\n', ...
        metLabels{f},localMetId(model,mIdx), ...
        localMetName(model,mIdx),localMetFormula(model,mIdx));
end

% Map every motif metabolite index to its role label. Merge labels when two
% roles are aliases of the same model metabolite, for example m_d=m_Fg2.
labelMap = containers.Map('KeyType','double','ValueType','char');
for f = 1:numel(metFields)
    mIdx = double(rec.(metFields{f}));
    role = metFields{f};
    if isKey(labelMap,mIdx)
        labelMap(mIdx) = [labelMap(mIdx),'=',role];
    else
        labelMap(mIdx) = role;
    end
end

fprintf('\n--- Reactions ---\n');
for f = 1:numel(rxnFields)
    field = rxnFields{f};
    rIdx = double(rec.(field));
    localValidateRxnIndex(rIdx,nRxns,field);

    fprintf('\n  %-54s\n',rxnLabels{f});
    fprintf('      ID  : %s\n',localRxnId(model,rIdx));
    fprintf('      Name: %s\n',localRxnName(model,rIdx));
    fprintf('      GPR : %s\n',localGPR(model,rIdx));
    fprintf('      EC  : %s\n',localEC(rIdx,printerECValues));
    fprintf('      Compartment(s): %s\n', ...
        localReactionCompartments(model,rIdx));
    localPrintReaction(model,rIdx,labelMap,1);
end

fprintf('\n############################################################\n');
fprintf('END OF SELECTED MOTIF %d\n',MOTIF_TO_PRINT);
fprintf('############################################################\n');

%% Local functions
function filePath = localFindFile(fileName,scriptFolder)
candidates = {fullfile(scriptFolder,fileName),fullfile(pwd,fileName)};
filePath = '';
for k = 1:numel(candidates)
    if isfile(candidates{k})
        filePath = candidates{k};
        return;
    end
end
end

function localValidateMetIndex(idx,nMets,field)
if ~isscalar(idx) || ~isfinite(idx) || idx < 1 || ...
        idx > nMets || idx ~= floor(idx)
    error('Printer:InvalidMetaboliteIndex', ...
        'Field %s contains invalid metabolite index %g.',field,idx);
end
end

function localValidateRxnIndex(idx,nRxns,field)
if ~isscalar(idx) || ~isfinite(idx) || idx < 1 || ...
        idx > nRxns || idx ~= floor(idx)
    error('Printer:InvalidReactionIndex', ...
        'Field %s contains invalid reaction index %g.',field,idx);
end
end

function value = localMetId(model,idx)
value = string(model.mets(idx));
end

function value = localMetName(model,idx)
value = localTextField(model,'metNames',idx,localMetId(model,idx));
end

function value = localMetFormula(model,idx)
value = localTextField(model,'metFormulas',idx,localMetId(model,idx));
end

function value = localRxnId(model,idx)
value = string(model.rxns(idx));
end

function value = localRxnName(model,idx)
value = localTextField(model,'rxnNames',idx,localRxnId(model,idx));
end

function value = localGPR(model,idx)
value = "(none annotated)";
if isfield(model,'grRules') && numel(model.grRules) >= idx
    rawValues = model.grRules;
    candidate = strtrim(string(rawValues(idx)));
    if strlength(candidate) > 0
        value = candidate;
        return;
    end
end
if isfield(model,'rules') && numel(model.rules) >= idx
    rawValues = model.rules;
    candidate = strtrim(string(rawValues(idx)));
    if strlength(candidate) > 0
        value = candidate;
    end
end
end

function value = localEC(idx,codes)
value = "(none annotated)";
if strlength(codes(idx)) > 0
    value = codes(idx);
end
end

function value = localTextField(model,field,idx,fallback)
value = string(fallback);
if isfield(model,field) && numel(model.(field)) >= idx
    rawValues = model.(field);
    candidate = strtrim(string(rawValues(idx)));
    if strlength(candidate) > 0
        value = candidate;
    end
end
end

function text = localReactionCompartments(model,rIdx)
metIdx = find(model.S(:,rIdx) ~= 0);
if isempty(metIdx)
    text = "(none)";
    return;
end

if isfield(model,'metComps') && isfield(model,'compNames') && ...
        numel(model.metComps) == size(model.S,1)
    compIdx = unique(double(model.metComps(metIdx)),'stable');
    valid = compIdx >= 1 & compIdx <= numel(model.compNames);
    compIdx = compIdx(valid);
    if ~isempty(compIdx)
        text = strjoin(string(model.compNames(compIdx))," <-> ");
        return;
    end
end

% Fallback when the SBML importer did not create model.metComps.
symbols = strings(numel(metIdx),1);
for k = 1:numel(metIdx)
    metID = char(localMetId(model,metIdx(k)));
    token = regexp(metID,'\[([^\]]+)\]$','tokens','once');
    if isempty(token)
        token = regexp(metID,'_([^_]+)$','tokens','once');
    end
    if isempty(token)
        symbols(k) = "unknown";
    else
        symbols(k) = string(token{1});
    end
end
text = strjoin(unique(symbols,'stable')," <-> ");
end

function localPrintReaction(model,rIdx,labelMap,direction)
col = full(model.S(:,rIdx))*direction;
subIdx = find(col < 0);
prodIdx = find(col > 0);

fprintf('      Eqn   : ');
localPrintSide(model,subIdx,col,true,labelMap,'formula');
fprintf('  ->  ');
localPrintSide(model,prodIdx,col,false,labelMap,'formula');
fprintf('\n');

fprintf('      Names : ');
localPrintSide(model,subIdx,col,true,labelMap,'name');
fprintf('  ->  ');
localPrintSide(model,prodIdx,col,false,labelMap,'name');
fprintf('\n');
end

function localPrintSide(model,idxList,col,isSubstrate,labelMap,kind)
if isempty(idxList)
    fprintf('(none)');
    return;
end

for k = 1:numel(idxList)
    mIdx = idxList(k);
    if isSubstrate
        coefficient = -col(mIdx);
    else
        coefficient = col(mIdx);
    end

    if strcmp(kind,'formula')
        component = localMetFormula(model,mIdx);
    else
        component = localMetName(model,mIdx);
    end

    if abs(coefficient-1) < 1e-12
        coefficientText = '';
    else
        coefficientText = sprintf('%.3g ',coefficient);
    end

    if isKey(labelMap,double(mIdx))
        token = sprintf('%s%s(%s)',coefficientText,component, ...
            labelMap(double(mIdx)));
    else
        token = sprintf('%s%s',coefficientText,component);
    end
    fprintf('%s',token);

    if k < numel(idxList)
        fprintf(' + ');
    end
end
end
