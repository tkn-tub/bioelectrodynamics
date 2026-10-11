function test_ec_search_update(runParallel)
%TEST_EC_SEARCH_UPDATE Tests require no COBRA installation or real model.
% Run test_ec_search_update for serial tests. To exercise parfor too, first
% start your allocated HPC pool and then run test_ec_search_update(true).
if nargin < 1, runParallel = false; end
sourceFolder = fileparts(mfilename('fullpath'));
addpath(sourceFolder);

% Missing annotations in the first source must not hide a complete code in
% another source. Nested values, class 7, EC prefixes and partial codes occur
% in imported models and need different decisions.
model = struct('S',sparse(1,8),'rxns',{{'a';'b';'c';'d';'e';'f';'g';'h'}});
model.eccodes = {'';'(none annotated)';'2.3.3.-';'7.1.1.2'; ...
    'EC1.2.3.4';{'1.2.3.4','2.7.1.1'};'1.2.3.4.5';string(missing)};
model.rxnECNumbers = {'';'3.1.4.17';'';'';'';'';'';''};
[hasEC,codes] = getReactionECAnnotations(model);
assert(isequal(hasEC,[false;true;false;true;true;true;false;false]));
assert(codes(2) == "3.1.4.17" && codes(6) == "1.2.3.4; 2.7.1.1");
model = rmfield(model,{'eccodes','rxnECNumbers'});
assert(~any(getReactionECAnnotations(model)));

% Provenance survives removing an earlier upstream row.
rows = table([1;2;4],[2;4;5],'VariableNames',{'r_A','r_B'});
[rows,index] = requireECAnnotatedRows(rows,hasEC,{'r_A','r_B'});
assert(isequal(index,[2;3]) && height(rows) == 2);

% Reverse labels are annotations on a direction, not part of reaction IDs.
model = struct('S',sparse([-1;1]),'mets',{{'a[c]';'b[c]'}}, ...
    'rxns',{{'rev'}},'eccodes',{{'1.2.3.4'}});
path = struct('anchorMetIdx',1,'anchorMetID',"a[c]", ...
    'reactionIndices',1,'reactionPath',"rev (reverse)",'pathIndex',9);
[paths,index] = requireECAnnotatedAdapters(model,path);
assert(numel(paths) == 1 && index == 1 && paths.pathIndex == 9);
model.eccodes = {''};
assert(numel(requireECAnnotatedAdapters(model,path)) == 1);
assert(isempty(requireECAnnotatedAdapters(model,path,false,true)));
path.anchorMetID = "wrong[c]";
try
    requireECAnnotatedAdapters(model,path);
    error('Test:ExpectedFailure','A mismatched model should be rejected.');
catch err
    assert(strcmp(err.identifier,'EC:AdapterModelMismatch'));
end

% The pipeline fixture puts unannotated duplicates BEFORE every eligible
% reaction and uses hit limits of one with REQUIRE_GPR=false. Post-filtering
% capped raw hits would produce no result; filtering before caps produces
% a complete, EC-annotated connection through Blocks 0-4.
testFolder = tempname;
mkdir(testFolder);
cleanup = onCleanup(@() rmdir(testFolder,'s'));
files = dir(fullfile(sourceFolder,'*.m'));
for k = 1:numel(files)
    copyfile(fullfile(sourceFolder,files(k).name),testFolder);
end
serialTable = localPipeline(testFolder,false,false);
strictTable = localPipeline(testFolder,false,true);
assert(isequal(serialTable,strictTable));
if runParallel
    if exist('gcp','file') ~= 2 || isempty(gcp('nocreate'))
        error('Test:MissingPool','Start your allocated parallel pool before the parallel test.');
    end
    parallelTable = localPipeline(testFolder,true,false);
    assert(isequal(serialTable,parallelTable));
end
fprintf('EC annotation and hit-limit regression tests passed.\n');
clear cleanup;
end

function result = localPipeline(folder,parallelMode,requireAdapterEC)
ids = {'camp[c]';'amp[c]';'a[c]';'s1[c]';'f1[c]';'sg1[c]';'fg1[c]'; ...
    's2[c]';'f2[c]';'sg2[c]';'d[c]';'c[c]';'cs[c]';'out[c]';'q[c]';'junk[c]'};
substrates = {1,2,3,4,4,3,5,5,[6 7],6,8,8,7,9,9,[10 11],11,[12 13],[11 13],12,15,15};
products = {2,3,4,6,16,5,7,16,16,8,10,16,9,11,16,16,12,16,16,14,14,13};
net = sparse(numel(ids),numel(substrates));
for r = 1:numel(substrates)
    net(substrates{r},r) = -1;
    net(products{r},r) = 1;
end
nr = size(net,2);
model = struct;
model.S = [net,net];
model.mets = ids;
model.metNames = ids;
model.rxns = cellstr(["unannotated_"+string((1:nr).');"annotated_"+string((1:nr).')]);
model.rxnNames = model.rxns;
model.grRules = repmat({'gene1'},2*nr,1);
model.genes = {'gene1'};
model.lb = zeros(2*nr,1);
model.ub = ones(2*nr,1)*1000;
model.eccodes = [repmat({''},nr,1);repmat({'1.1.1.1'},nr,1)];
S = model.S;
isBadMet = false(numel(ids),1);
isBadMet([2,16]) = true;
REQUIRE_GPR = false;
BLOCK0_REQUIRE_EC = requireAdapterEC;
SHOW_PROGRESS = false;
FILTER_INTERNAL_CROSSTALK_DURING_SEARCH = false;
AUTO_START_PARPOOL = false;
USE_PARALLEL = parallelMode;
USE_PARALLEL_BLOCK2 = parallelMode;
USE_PARALLEL_COUNTER = parallelMode;
USE_PARALLEL_NN = parallelMode;
BLOCK0_MAX_DEPTH = 2;
BLOCK0_MAX_PATHS = 1;
BLOCK0_MAX_PRINT_PATHS = 0;
MAX_HITS = 1;
MAX_HITS_BLOCK2 = 1;
MAX_STAGE2_PER_PAIR = 1;
MAX_BRANCH_OPTIONS_PER_INPUT = 100;
MAX_STAGE2_OPTION_PAIRS_TESTED = 100;
MAX_HITS_COUNTER = 1;
MAX_COUNTER_CANDIDATES_PER_INPUT = 1;
MAX_HITS_NN = 1;
MAX_NN_CANDIDATES_PER_INPUT = 1;
MAX_FULL_MOTIFS_NN = 1;
FULL_TABLE_SAFETY_CAP = 1;
run(fullfile(folder,'Block_0_Input_Adapter.m'));
assert(numel(InputAdapterPaths) == 1);
if requireAdapterEC
    assert(all(InputAdapterPaths(1).reactionIndices > nr));
else
    assert(all(InputAdapterPaths(1).reactionIndices <= nr));
end
run(fullfile(folder,'Block_1_Differentiator_Mapping.m'));
assert(height(Motifs_Diff_Block_1) == 1);
run(fullfile(folder,'Block_2_Differentiator_Mapping.m'));
assert(height(Motifs_Diff_Block_2) == 1);
run(fullfile(folder,'Block_3_Counter_Mapping.m'));
assert(height(Motifs_Counter_Block_3) == 1);
run(fullfile(folder,'Block_4_Single_Neuron_NN_Mapping.m'));
assert(height(Motifs_NN_Block_4) == 1);
roles = Motifs_NN_Block_4.Properties.VariableNames;
roles = roles(startsWith(roles,'r_'));
assert(all(Motifs_NN_Block_4{:,roles} > nr));
assert(isequal(model.S,[net,net]) && isequal(S,model.S));
result = Motifs_NN_Block_4;
end
