function test_five_reaction_differentiators(runParallel)
% Native MATLAB tests using independently specified chemical equations.
% No COBRA or real model required. Optional allocated-pool test with true.
if nargin < 1, runParallel = false; end
addpath(fileparts(mfilename('fullpath')));
[model,bad,paths,o] = fixture();
[one,~,~,info] = mapFiveReactionDifferentiators(model,bad,paths,table(),o);
assert(height(one) == 1 && one.m_in == 2 && one.m_M1 == 3 && one.m_Y1 == 4 && one.m_F1_delay == 5);
assert(isequal(one{:,{'r_S1','r_S2','r_F','r_I','r_R'}},12:16));
assert(info.globalHitLimitReached && one.adapterPathIndex == 11);
o.stage = 2;
[two,diagnostics,~,info] = mapFiveReactionDifferentiators(model,bad,paths,one,o);
assert(height(two) == 1 && two.m_X2 == 4 && two.m_d1 == 4 && two.m_d == 4 && two.m_d2 == 7);
assert(isequal(two{:,{'r_S1_2','r_S2_2','r_F2','r_I2','r_R2'}},17:21));
assert(info.globalHitLimitReached && diagnostics.InputIndex == 4);
schema = minimalCircuitSchema(2);
assert(checkMinimalCircuit(model,table2struct(two),2));
assert(numel(schema.rxns) == 10 && numel(schema.mets) == 6);

% Native searches, not only the schema, must reject the old reaction graphs.
oldSlow = model; oldSlow.S(7,18) = 0;
[none,~,~,~] = mapFiveReactionDifferentiators(oldSlow,bad,paths,one,o);
assert(isempty(none));
oldInhibition = model; oldInhibition.S(6,20) = 0; oldInhibition.S(4,20) = -1;
[none,~,~,~] = mapFiveReactionDifferentiators(oldInhibition,bad,paths,one,o);
assert(isempty(none));
oldFirst = model; oldFirst.S(5,15) = 0; oldFirst.S(2,15) = -1;
firstOptions = o; firstOptions.stage = 1;
[none,~,~,~] = mapFiveReactionDifferentiators(oldFirst,bad,paths,table(),firstOptions);
assert(isempty(none));
oldRows = one; oldRows(:,{'m_F1_delay','r_I','r_R','dir_I','dir_R'}) = [];
assertError(@() mapFiveReactionDifferentiators(model,bad,paths,oldRows,o),'MinimalDiff:OldBlock1');

% All ten controlled roles require EC and feasible bounds before hit counting.
for reaction = 12:21
    missing = model; missing.eccodes{reaction} = '';
    blocked = model; blocked.ub(reaction) = 0;
    reverse = model; reverse.S(:,reaction) = -reverse.S(:,reaction);
    reverse.lb(reaction) = -1000; reverse.ub(reaction) = 0;
    if reaction <= 16
        [emptyEC,~,~,~] = mapFiveReactionDifferentiators(missing,bad,paths,table(),firstOptions);
        [emptyBounds,~,~,~] = mapFiveReactionDifferentiators(blocked,bad,paths,table(),firstOptions);
        [signedOne,~,~,~] = mapFiveReactionDifferentiators(reverse,bad,paths,table(),firstOptions);
        [signedTwo,~,~,~] = mapFiveReactionDifferentiators(reverse,bad,paths,signedOne,o);
    else
        [emptyEC,~,~,~] = mapFiveReactionDifferentiators(missing,bad,paths,one,o);
        [emptyBounds,~,~,~] = mapFiveReactionDifferentiators(blocked,bad,paths,one,o);
        [signedTwo,~,~,~] = mapFiveReactionDifferentiators(reverse,bad,paths,one,o);
    end
    assert(isempty(emptyEC) && isempty(emptyBounds));
    assert(height(signedTwo) == 1 && checkMinimalCircuit(reverse,table2struct(signedTwo),2));
end
strict = o; strict.requireAdapterEC = true;
assertError(@() mapFiveReactionDifferentiators(model,bad,paths,one,strict),'MinimalDiff:NoAdapters');

% Other network producers are allowed, but selected cross-talk is rejected.
external = model; external.S(:,22) = sparse([11 7],[1 1],[-1 1],11,1);
external.rxns{22} = 'external_Y2'; external.rxnNames{22} = 'external Y2 production';
external.eccodes{22} = ''; external.grRules{22} = ''; external.lb(22) = 0; external.ub(22) = 1000;
[allowed,~,~,~] = mapFiveReactionDifferentiators(external,bad,paths,one,o);
assert(height(allowed) == 1);
cross = model; cross.S(5,19) = -1;
[none,~,~,crossInfo] = mapFiveReactionDifferentiators(cross,bad,paths,one,o);
assert(isempty(none) && crossInfo.nRejectedJoins > 0);

% Structural and join budgets stop independently of successful matches.
limited = o; limited.maxChecks = 1;
[~,diag,~,limitedInfo] = mapFiveReactionDifferentiators(model,bad,paths,one,limited);
assert(diag.CheckLimitReached && limitedInfo.perInputSearchTruncated);
limited = o; limited.maxHits = 100; limited.maxJoinChecks = 1;
[limitedRows,~,~,limitedInfo] = mapFiveReactionDifferentiators(model,bad,paths,[one;one],limited);
assert(height(limitedRows) == 1 && limitedInfo.joinChecks == 1 && limitedInfo.joinCheckLimitReached);

folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder,'s'));
payload = fiveReactionDifferentiatorCheckpoint(model,two,diagnostics,paths,o,info);
payload.Motifs_Diff_Block_2 = two;
saveCompactCheckpoint(fullfile(folder,'roundtrip.mat'),payload);
loaded = loadMinimalCheckpoint(fullfile(folder,'roundtrip.mat'),model,'Motifs_Diff_Block_2');
assert(isequal(loaded.Motifs_Diff_Block_2,two));
assert(~isfield(loaded,'model') && ~isfield(loaded,'results_diff2'));
report = renderFiveReactionDifferentiators(loaded,1);
assert(contains(report,'r_I') && contains(report,'m_F1_delay') && contains(report,'TEN-REACTION'));
localScripts(model,bad,paths,folder);
if runParallel
    if exist('gcp','file') ~= 2 || isempty(gcp('nocreate'))
        error('FiveReactionTest:Pool','Start your allocated MATLAB pool first.');
    end
    p = firstOptions; p.useParallel = true;
    [parallelOne,~,~,parallelInfo] = mapFiveReactionDifferentiators(model,bad,paths,table(),p);
    assert(isequal(one,parallelOne) && parallelInfo.parallelUsed);
    p.stage = 2;
    [parallelTwo,~,~,parallelInfo] = mapFiveReactionDifferentiators(model,bad,paths,one,p);
    assert(isequal(two,parallelTwo) && parallelInfo.parallelUsed);
end
fprintf('Five-reaction differentiator regression tests passed.\n');
clear cleanup;
end

function localScripts(model,isBadMet,InputAdapterPaths,folder)
REQUIRE_GPR = false; SHOW_PROGRESS = false; BLOCK0_REQUIRE_EC = false;
USE_PARALLEL = false; USE_PARALLEL_BLOCK2 = false;
MAX_HITS = 1; MAX_HITS_BLOCK2 = 1; BLOCK1_OUTPUT_DIR = folder; BLOCK2_OUTPUT_DIR = folder;
previous = pwd; cleanup = onCleanup(@() cd(previous)); cd(folder);
Block_1_Differentiator_Mapping;
one = load('Block_1_differentiator.mat');
assert(one.Motifs_Diff_Block_1.m_d1 == 4 && one.Motifs_Diff_Block_1.m_F1_delay == 5);
clear Motifs_Diff_Block_1 results_diff;
Block_2_Differentiator_Mapping;
two = load('Block_2_differentiator.mat'); assert(two.Motifs_Diff_Block_2.m_d2 == 7);
Print_Five_Reaction_Differentiators;
assert(isfile('Differentiator_Motif_1.txt'));
clear cleanup;
end

function [model,bad,paths,o] = fixture()
ids = {'camp[c]';'X[c]';'M1[c]';'Y1[c]';'F1[c]';'M2[c]';'Y2[c]';'W1[c]';'W2[c]';'W3[c]';'currency[c]'};
% Independently declared equations from the October 6 kinetics discussion.
left = {2,3,2,[5 4],5,4,6,4,[6 7],7};
right = {3,4,5,8,8,6,7,7,9,10};
core = sparse(11,10);
for r = 1:10, core(left{r},r) = -1; core(right{r},r) = 1; end
adapter = sparse([1 2],[1 1],[-1 1],11,1);
model = struct('S',[adapter core core],'mets',{ids},'metNames',{ids});
model.rxns = cellstr(["adapter";"no_ec_"+string((1:10).');"annotated_"+string((1:10).')]);
model.rxnNames = model.rxns; model.grRules = repmat({'gene'},21,1);
model.lb = zeros(21,1); model.ub = repmat(1000,21,1);
model.eccodes = [repmat({''},11,1);repmat({'1.1.1.1'},10,1)];
bad = false(11,1); bad(8:11) = true;
paths = struct('pathIndex',11,'inputMetIdx',1,'inputMetID',"camp[c]", ...
    'anchorMetIdx',2,'anchorMetID',"X[c]",'reactionIndices',1, ...
    'reactionPath',"adapter",'metabolitePath',"camp[c] -> X[c]");
o = struct('stage',1,'maxHits',1,'maxInputs',inf,'perInputCap',100,'maxChecks',inf, ...
    'maxJoinChecks',inf,'requireGPR',false,'requireAdapterEC',false,'showProgress',false, ...
    'useParallel',false,'autoPool',false,'workers',[],'distinctAdapter',true, ...
    'distinctStageRxns',true,'distinctStageMets',true,'crossTalk',true);
end

function assertError(action,identifier)
try
    action(); error('FiveReactionTest:ExpectedError','Expected %s.',identifier);
catch failure
    assert(strcmp(failure.identifier,identifier),failure.message);
end
end
