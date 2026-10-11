function test_polarity_differentiators(runParallel)
% MATLAB regression tests. No COBRA, real model or external files required.
% Run test_polarity_differentiators(true) with an already allocated pool to
% include a parallel hit-limit test. All generated files use a temp folder.
if nargin < 1, runParallel = false; end
addpath(fileparts(mfilename('fullpath')));
[model,bad,paths,settings] = fixture();
[one,~,~,info] = mapPolarityDifferentiators(model,bad,paths,table(),settings);
assert(height(one) == 1 && one.m_in == 2 && one.m_M1 == 3 && one.m_Y1 == 4);
assert(isequal(one{:,{'r_S1','r_S2','r_F'}},[10 11 12]));
assert(info.globalHitLimitReached && one.adapterPathIndex == 11);
settings.stage = 2;
[two,diagnostics,~,info] = mapPolarityDifferentiators(model,bad,paths,one,settings);
assert(height(two) == 1 && two.m_X2 == two.m_Y1 && two.m_d == 4 && two.m_d2 == 6);
assert(isequal(two{:,{'r_S1_2','r_S2_2','r_F2','r_I2','r_R2'}},[13 14 15 16 17]));
assert(info.globalHitLimitReached && diagnostics.InputIndex == 4);
% The unparameterized chemical graph permits both slow/fast role assignments.
% The search must keep these distinct, because M2 and Y2 are different outputs.
unlimited = settings; unlimited.maxHits = inf;
[both,~,~,~] = mapPolarityDifferentiators(model,bad,paths,one,unlimited);
assert(height(both) == 2 && numel(unique(both.m_Y2)) == 2);

% Direction signs must cover all three removal roles, including shared removal.
for reaction = [14 16 17]
    reverse = model; reverse.S(:,reaction) = -reverse.S(:,reaction);
    reverse.lb(reaction) = -1000; reverse.ub(reaction) = 0;
    [reversed,~,~,~] = mapPolarityDifferentiators(reverse,bad,paths,one,settings);
    assert(height(reversed) == 1);
    if reaction == 14, assert(reversed.dir_S2_2 == -1); end
    if reaction == 16, assert(reversed.dir_I2 == -1); end
    if reaction == 17, assert(reversed.dir_R2 == -1); end
end

% The old second-stage chain and its old Y1+Y2 sink are rejected separately.
oldChain = model; oldChain.S(4,15) = 0; oldChain.S(5,15) = -1;
[none,~,~,~] = mapPolarityDifferentiators(oldChain,bad,paths,one,settings);
assert(isempty(none));
oldSink = model; oldSink.S(5,16) = 0; oldSink.S(4,16) = -1;
[none,~,~,~] = mapPolarityDifferentiators(oldSink,bad,paths,one,settings);
assert(isempty(none));
oldSlowOutput = model; oldSlowOutput.S(6,14) = 1;
[none,~,~,~] = mapPolarityDifferentiators(oldSlowOutput,bad,paths,one,settings);
assert(isempty(none));

% Y2 removal r_R2 is mandatory. Missing EC or blocked flux cannot fill MAX_HITS.
for reaction = 13:17
    noEC = model; noEC.eccodes{reaction} = '';
    [none,~,~,~] = mapPolarityDifferentiators(noEC,bad,paths,one,settings);
    assert(isempty(none));
    blocked = model; blocked.ub(reaction) = 0;
    [none,~,~,~] = mapPolarityDifferentiators(blocked,bad,paths,one,settings);
    assert(isempty(none));
end
strict = settings; strict.requireAdapterEC = true;
assertError(@() mapPolarityDifferentiators(model,bad,paths,one,strict),'MinimalDiff:NoAdapters');

% Existing three-reaction Block 1 checkpoints can be renamed and extended.
legacy = one;
old = {'r_11','r_12','r_13','dir_11','dir_12','dir_13'};
current = {'r_S1','r_S2','r_F','dir_S1','dir_S2','dir_F'};
for k = 1:numel(old), legacy.(old{k}) = legacy.(current{k}); end
legacy(:,current) = [];
[compatible,~,~,~] = mapPolarityDifferentiators(model,bad,paths,legacy,settings);
assert(isequal(two,compatible));

% Producers outside the selected motif remain allowed; selected cross-talk fails.
external = model;
external.S(:,18) = sparse([10 6],[1 1],[-1 1],10,1);
external.rxns{18} = 'external_Y2_production'; external.rxnNames{18} = 'external production';
external.eccodes{18} = ''; external.grRules{18} = '';
external.lb(18) = 0; external.ub(18) = 1000;
[outside,~,~,~] = mapPolarityDifferentiators(external,bad,paths,one,settings);
assert(height(outside) == 1);
cross = model; cross.S(3,15) = -1;
[none,~,~,crossInfo] = mapPolarityDifferentiators(cross,bad,paths,one,settings);
assert(isempty(none) && crossInfo.nRejectedJoins > 0);

% Both structural-search and join budgets terminate independently of hit counts.
limited = settings; limited.maxChecks = 1;
[~,diag,limitedInfo] = localRun(model,bad,paths,one,limited);
assert(diag.CheckLimitReached && limitedInfo.perInputSearchTruncated);
limited = settings; limited.maxHits = 100; limited.maxJoinChecks = 1;
[limitedRows,~,~,limitedInfo] = mapPolarityDifferentiators(model,bad,paths,[one;one],limited);
assert(height(limitedRows) == 1 && limitedInfo.joinChecks == 1 && limitedInfo.joinCheckLimitReached);

% MAT round-trip contains the table once, with selected annotation dictionaries.
folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder,'s'));
checkpoint = polarityDifferentiatorCheckpoint(model,two,diagnostics,paths,settings,info);
checkpoint.Motifs_Diff_Block_2 = two;
saveCompactCheckpoint(fullfile(folder,'roundtrip.mat'),checkpoint);
loaded = load(fullfile(folder,'roundtrip.mat'));
assert(isequal(loaded.Motifs_Diff_Block_2,two));
assert(~isfield(loaded,'model') && ~isfield(loaded,'results_diff2'));
report = renderPolarityDifferentiators(loaded,1);
assert(contains(report,'r_R2') && contains(report,'EIGHT-REACTION'));
localScriptRoundTrip(model,bad,paths,folder);
if runParallel
    if exist('gcp','file') ~= 2 || isempty(gcp('nocreate'))
        error('PolarityTest:Pool','Start your allocated MATLAB pool first.');
    end
    [~,~,~,parallelSettings] = fixture(); parallelSettings.useParallel = true;
    [parallelOne,~,~,parallelInfo] = mapPolarityDifferentiators(model,bad,paths,table(),parallelSettings);
    assert(isequal(one,parallelOne) && parallelInfo.parallelUsed);
    parallelSettings.stage = 2;
    [parallelTwo,~,~,parallelInfo] = mapPolarityDifferentiators(model,bad,paths,one,parallelSettings);
    assert(isequal(two,parallelTwo) && parallelInfo.parallelUsed);
end
fprintf('Polarity differentiator regression tests passed.\n');
clear cleanup;
end

function [rows,diag,info] = localRun(model,bad,paths,one,settings)
[rows,diag,~,info] = mapPolarityDifferentiators(model,bad,paths,one,settings);
end

function localScriptRoundTrip(model,isBadMet,InputAdapterPaths,folder)
REQUIRE_GPR = false; SHOW_PROGRESS = false;
USE_PARALLEL = false; USE_PARALLEL_BLOCK2 = false; BLOCK0_REQUIRE_EC = false;
MAX_HITS = 1; MAX_HITS_BLOCK2 = 1;
BLOCK1_OUTPUT_DIR = folder; BLOCK2_OUTPUT_DIR = folder;
previous = pwd; cleanup = onCleanup(@() cd(previous)); cd(folder);
Block_1_Differentiator_Mapping;
first = load('Block_1_differentiator.mat');
assert(height(first.Motifs_Diff_Block_1) == 1 && first.Motifs_Diff_Block_1.m_d == 4);
clear Motifs_Diff_Block_1 results_diff;
Block_2_Differentiator_Mapping;
second = load('Block_2_differentiator.mat');
assert(height(second.Motifs_Diff_Block_2) == 1 && second.Motifs_Diff_Block_2.m_d2 == 6);
Print_Polarity_Differentiators;
assert(isfile('Differentiator_Motif_1.txt'));
clear cleanup;
end

function [model,bad,paths,s] = fixture()
ids = {'camp[c]';'X[c]';'M1[c]';'Y1[c]';'M2[c]';'Y2[c]';'W1[c]';'W2[c]';'W3[c]';'currency[c]'};
% These eight equations are independently specified from the attached diagram.
substrates = {2,3,[2 4],4,5,4,[5 6],6};
products = {3,4,7,5,8,6,8,9};
core = sparse(10,8);
for r = 1:8, core(substrates{r},r) = -1; core(products{r},r) = 1; end
adapter = sparse([1 2],[1 1],[-1 1],10,1);
model.S = [adapter core core];
model.mets = ids; model.metNames = ids;
model.rxns = cellstr(["adapter_missing_ec";"no_ec_"+string((1:8).');"annotated_"+string((1:8).')]);
model.rxnNames = model.rxns;
model.lb = zeros(17,1); model.ub = repmat(1000,17,1);
model.eccodes = [repmat({''},9,1);repmat({'1.1.1.1'},8,1)];
model.grRules = repmat({'gene'},17,1);
bad = false(10,1); bad(7:10) = true;
paths = struct('pathIndex',11,'inputMetIdx',1,'inputMetID',"camp[c]", ...
    'anchorMetIdx',2,'anchorMetID',"X[c]",'reactionIndices',1, ...
    'reactionPath',"adapter_missing_ec",'metabolitePath',"camp[c] -> X[c]");
s = struct('stage',1,'maxHits',1,'maxInputs',inf,'perInputCap',100, ...
    'maxChecks',inf,'maxJoinChecks',inf,'requireGPR',false,'requireAdapterEC',false, ...
    'showProgress',false,'useParallel',false,'autoPool',false, ...
    'workers',[],'distinctAdapter',true,'distinctStageRxns',true, ...
    'distinctStageMets',true,'crossTalk',true);
end

function assertError(action,identifier)
try
    action(); error('PolarityTest:ExpectedError','Expected %s.',identifier);
catch failure
    assert(strcmp(failure.identifier,identifier),failure.message);
end
end
