function test_minimal_differentiators(runParallel)
% No COBRA or real model required. Use true with an allocated MATLAB pool
% for the additional asynchronous search test. Writes only to a temp folder.
if nargin < 1, runParallel = false; end
source = fileparts(mfilename('fullpath')); addpath(source);
[model,bad,paths,settings] = fixture();
[one,~,~,info] = mapMinimalDifferentiator(model,bad,paths,table(),settings);
assert(height(one) == 1 && one.m_in == 2 && one.m_M1 == 3 && one.m_Y1 == 4);
assert(all(one{:,{'r_11','r_12','r_13'}} == [8 9 10]));
assert(info.globalHitLimitReached && one.adapterPathIndex == 11);
settings.stage = 2;
[two,diagnostics,~,info] = mapMinimalDifferentiator(model,bad,paths,one,settings);
assert(height(two) == 1 && two.m_X2 == two.m_Y1 && two.m_Y2 == 6 && two.m_d == 6);
assert(all(two{:,{'r_21','r_22','r_23'}} == [11 12 13]));
assert(info.globalHitLimitReached && diagnostics.InputIndex == 4);

% Reversible reactions retain canonical IDs with an explicit feasible sign.
reverse = model; reverse.S(:,10) = -reverse.S(:,10);
reverse.lb(10) = -1000; reverse.ub(10) = 0;
settings.stage = 1;
[rev,~,~,~] = mapMinimalDifferentiator(reverse,bad,paths,table(),settings);
assert(height(rev) == 1 && rev.dir_13 == -1);

% Missing EC and blocked flux are excluded before MAX_HITS is counted.
noEC = model; noEC.eccodes{10} = '';
[none,~,~,~] = mapMinimalDifferentiator(noEC,bad,paths,table(),settings);
assert(isempty(none));
blocked = model; blocked.ub(9) = 0;
[none,~,~,~] = mapMinimalDifferentiator(blocked,bad,paths,table(),settings);
assert(isempty(none));
strict = settings; strict.requireAdapterEC = true;
assertError(@() mapMinimalDifferentiator(model,bad,paths,table(),strict),'MinimalDiff:NoAdapters');

% The joint sink must consume both X and Y. Two independent sinks fail.
separate = model; separate.S(2,10) = 0;
[none,~,~,~] = mapMinimalDifferentiator(separate,bad,paths,table(),settings);
assert(isempty(none));
extraProducer = model; extraProducer.S(3,10) = 1;
[none,~,~,~] = mapMinimalDifferentiator(extraProducer,bad,paths,table(),settings);
assert(isempty(none));

% External producers are unrestricted; selected unintended consumption fails.
external = model;
external.S(:,14) = sparse([9 4],[1 1],[-1 1],size(model.S,1),1);
external.rxns{14} = 'external_producer'; external.rxnNames{14} = 'external producer';
external.eccodes{14} = ''; external.grRules{14} = ''; external.lb(14) = 0; external.ub(14) = 1000;
[outside,~,~,~] = mapMinimalDifferentiator(external,bad,paths,table(),settings);
assert(height(outside) == 1);
cross = model; cross.S(3,11) = -1; settings.stage = 2;
[none,~,~,info] = mapMinimalDifferentiator(cross,bad,paths,one,settings);
assert(isempty(none) && info.nRejectedJoins > 0);
settings.maxChecks = 1; settings.stage = 1;
[~,diag,~,info] = mapMinimalDifferentiator(model,bad,paths,table(),settings);
assert(diag.CheckLimitReached && info.perInputSearchTruncated);

% Round-trip compact checkpoint, then render without model or COBRA.
settings.maxChecks = inf;
folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder,'s'));
checkpoint = minimalDifferentiatorCheckpoint(model,two,diagnostics,paths,settings,info);
checkpoint.Motifs_Diff_Block_2 = two;
saved = fullfile(folder,'Block_2_differentiator.mat');
saveCompactCheckpoint(saved,checkpoint);
loaded = load(saved);
assert(isequal(loaded.Motifs_Diff_Block_2,two));
assert(~isfield(loaded,'model') && ~isfield(loaded,'results_diff2'));
report = renderMinimalDifferentiator(loaded,1);
assert(contains(report,'r_23') && contains(report,'adapter_missing_ec'));
localScriptRoundTrip(model,bad,paths,folder);
if runParallel
    if exist('gcp','file') ~= 2 || isempty(gcp('nocreate'))
        error('MinimalTest:Pool','Start your allocated MATLAB pool first.');
    end
    [~,~,~,settings] = fixture(); settings.useParallel = true;
    [parallelOne,~,~,info] = mapMinimalDifferentiator(model,bad,paths,table(),settings);
    assert(isequal(parallelOne,one) && info.parallelUsed);
    secondPath = paths;
    secondPath.pathIndex = 12; secondPath.anchorMetIdx = 4;
    secondPath.anchorMetID = "Y1[c]"; secondPath.reactionIndices = [1 8 9];
    secondPath.reactionPath = "adapter_missing_ec -> good_1 -> good_2";
    secondPath.metabolitePath = "camp[c] -> X1[c] -> M1[c] -> Y1[c]";
    [capped,~,~,info] = mapMinimalDifferentiator(model,bad,[paths;secondPath],table(),settings);
    assert(height(capped) == 1 && info.nInputs == 2 && info.nInputsCompleted == 1);
end
fprintf('Minimal differentiator regression tests passed.\n');
clear cleanup;
end

function localScriptRoundTrip(model,isBadMet,InputAdapterPaths,folder)
REQUIRE_GPR = false; SHOW_PROGRESS = false;
USE_PARALLEL = false; USE_PARALLEL_BLOCK2 = false;
BLOCK0_REQUIRE_EC = false;
MAX_HITS = 1; MAX_HITS_BLOCK2 = 1;
BLOCK1_OUTPUT_DIR = folder; BLOCK2_OUTPUT_DIR = folder;
previous = pwd; cleanup = onCleanup(@() cd(previous)); cd(folder);
Block_1_Differentiator_Mapping;
first = load('Block_1_differentiator.mat');
assert(height(first.Motifs_Diff_Block_1) == 1);
% Exercise the automatic Block 1 checkpoint load, not the in-memory struct.
clear Motifs_Diff_Block_1 results_diff;
Block_2_Differentiator_Mapping;
second = load('Block_2_differentiator.mat');
assert(height(second.Motifs_Diff_Block_2) == 1 && second.Motifs_Diff_Block_2.m_d == 6);
Print_Minimal_Differentiators;
assert(isfile('Differentiator_Motif_1.txt'));
clear cleanup;
end

function [model,bad,paths,s] = fixture()
ids = {'camp[c]';'X1[c]';'M1[c]';'Y1[c]';'M2[c]';'Y2[c]';'W1[c]';'W2[c]';'currency[c]'};
substrates = {2,3,[2 4],4,5,[4 6]};
products = {3,4,7,5,6,8};
core = sparse(numel(ids),6);
for r = 1:6, core(substrates{r},r) = -1; core(products{r},r) = 1; end
adapter = sparse([1 2],[1 1],[-1 1],numel(ids),1);
model.S = [adapter core core];
model.mets = ids; model.metNames = ids;
model.rxns = cellstr(["adapter_missing_ec";"bad_"+string((1:6).');"good_"+string((1:6).')]);
model.rxnNames = model.rxns;
model.lb = zeros(13,1); model.ub = repmat(1000,13,1);
model.eccodes = [repmat({''},7,1);repmat({'1.1.1.1'},6,1)];
model.grRules = repmat({'gene'},13,1);
bad = false(9,1); bad(7:9) = true;
paths = struct('pathIndex',11,'inputMetIdx',1,'inputMetID',"camp[c]", ...
    'anchorMetIdx',2,'anchorMetID',"X1[c]",'reactionIndices',1, ...
    'reactionPath',"adapter_missing_ec",'metabolitePath',"camp[c] -> X1[c]");
s = struct('stage',1,'maxHits',1,'maxInputs',inf,'perInputCap',100, ...
    'maxChecks',inf,'requireGPR',false,'requireAdapterEC',false, ...
    'showProgress',false,'useParallel',false,'autoPool',false, ...
    'workers',[],'distinctAdapter',true,'distinctStageRxns',true, ...
    'distinctStageMets',true,'crossTalk',true);
end

function assertError(action,identifier)
try
    action(); error('MinimalTest:ExpectedError','Expected %s.',identifier);
catch failure
    assert(strcmp(failure.identifier,identifier),failure.message);
end
end
