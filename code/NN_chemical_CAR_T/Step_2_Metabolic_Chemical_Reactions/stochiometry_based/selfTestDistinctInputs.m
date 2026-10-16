function selfTestDistinctInputs(runParallelTest)
%SELFTESTDISTINCTINPUTS Representative search and compact result round-trip.
if nargin < 1, runParallelTest = false; end
[model,cfg,p] = fixture();
anchored = cfg;
for i = 1:numel(p.metNames)
    anchored.FIXED_METABOLITES.(p.metNames{i}) = model.mets{i};
end
anchored.ONE_MOTIF_PER_INPUT = false;
allMatches = findMetabolicCircuit(model,anchored);
assert(numel(allMatches) == 2,'Expected two circuits using the same anchored input.');
anchored.ONE_MOTIF_PER_INPUT = true;
[one,info] = findMetabolicCircuit(model,anchored);
assert(numel(one) == 1 && info.distinctInputsReturned == 1 && info.exhaustive);
fprintf('PASS: abandon an input after its first valid circuit.\n');

twin = duplicateModel(model);
[two,info] = findMetabolicCircuit(twin,cfg);
twinInfo = info;
assert(info.exhaustive && numel(two) == 2);
assert(numel(unique(arrayfun(@(m) m.metIndex(1),two))) == 2);
limited = cfg; limited.MAX_MOTIFS = 1;
one = findMetabolicCircuit(twin,limited);
assert(numel(one) == 1);
fprintf('PASS: distinct inputs and a limit on their representatives.\n');

alias = twin; n = numel(model.mets);
alias.mets{n+1} = 'm_in[m]'; alias.metComps(n+1) = 2; alias.comps{2} = 'm';
[one,info] = findMetabolicCircuit(alias,cfg);
assert(info.exhaustive && numel(one) == 1 && info.distinctInputsReturned == 1);
byID = cfg; byID.METABOLITE_IDENTITY = 'model_id';
twoByID = findMetabolicCircuit(alias,byID);
assert(numel(twoByID) == 2);
fprintf('PASS: compartment aliases remain in a single chemical input group.\n');

if runParallelTest
    parallelCfg = cfg; parallelCfg.USE_PARALLEL = true; parallelCfg.NUM_WORKERS = 2;
    [parallel,info] = findMetabolicCircuit(twin,parallelCfg);
    assert(info.parallelUsed && numel(parallel) == 2);
    assert(isequal(sort(arrayfun(@(m) m.metIndex(1),two)), ...
        sort(arrayfun(@(m) m.metIndex(1),parallel))));
    [parallel,info] = findMetabolicCircuit(alias,parallelCfg);
    assert(info.parallelUsed && numel(parallel) == 1);
    fprintf('PASS: disjoint parallel input groups and chemical aliases.\n');
end

testFolder = tempname; mkdir(testFolder);
cleanup = onCleanup(@() rmdir(testFolder,'s')); %#ok<NASGU>
sourceFile = fullfile(testFolder,'Circuit_Search_Results.mat');
model = twin; pattern = p; motifs = two; searchInfo = twinInfo;
% Exercise model structures, sparse stoichiometry and saved circuit variables.
save(sourceFile,'model','motifs','searchInfo','pattern','cfg','-v7.3');
[compactFile,stats] = Compact_Circuit_Search_Results(sourceFile);
assert(exist(sourceFile,'file') == 2 && exist(compactFile,'file') == 2);
assert(stats.verifiedVariables == 5);
saved = load(compactFile);
assert(isequaln(saved.model,model) && isequaln(saved.motifs,motifs));
assert(isequaln(saved.searchInfo,searchInfo) && isequaln(saved.cfg,cfg));
fprintf('PASS: compact results retain the model and all saved variables.\n');
fprintf('All distinct-input/compact-result MATLAB tests passed.\n');
end

function [model,cfg,p] = fixture()
p = circuitPattern('substrate'); S = p.T;
mets = cellfun(@(x) [x,'[c]'],p.metNames,'UniformOutput',false);
for j = 1:18
    if ~any(S(:,j) > 0)
        S(end+1,j) = 1; mets{end+1,1} = sprintf('waste_%02d[c]',j); %#ok<AGROW>
    end
end
model.S = sparse([S,S(:,1)]); model.mets = mets; model.metNames = mets;
model.metFormulas = repmat({''},numel(mets),1);
model.rxns = arrayfun(@(j) sprintf('synthetic_%02d',j),(1:19)','UniformOutput',false);
model.rxnNames = model.rxns;
model.rxnECNumbers = arrayfun(@(j) sprintf('1.1.1.%d',j),(1:19)','UniformOutput',false);
model.grRules = arrayfun(@(j) sprintf('G%d',j),(1:19)','UniformOutput',false);
model.lb = zeros(19,1); model.ub = 1000*ones(19,1);
model.comps = {'c'}; model.metComps = ones(numel(mets),1);
cfg = circuitSearchDefaults(); cfg.VERBOSE = false;
cfg.MAX_MOTIFS = Inf; cfg.MAX_SECONDS = 120; cfg.USE_PARALLEL = false;
end

function twin = duplicateModel(model)
twin = model; twin.S = blkdiag(model.S,model.S);
for field = {'mets','metNames','rxns','rxnNames','grRules'}
    f = field{1};
    twin.(f) = [model.(f);cellfun(@(x) ['copy_',x],model.(f),'UniformOutput',false)];
end
twin.metFormulas = [model.metFormulas;model.metFormulas];
twin.metComps = [model.metComps;model.metComps];
twin.lb = [model.lb;model.lb]; twin.ub = [model.ub;model.ub];
twin.rxnECNumbers = [model.rxnECNumbers;arrayfun(@(j) ...
    sprintf('2.1.1.%d',j),(1:numel(model.rxns))','UniformOutput',false)];
end
