function selfTestCircuitSearch(runParallelTest)
%SELFTESTCIRCUITSEARCH Regression checks using synthetic (not cell) chemistry.
% Run: selfTestCircuitSearch
% Optional pool check: selfTestCircuitSearch(true)
if nargin < 1, runParallelTest = false; end
[model,cfg,p] = fixture();
cfg.MAX_MOTIFS = 1;
[found,~,~] = findMetabolicCircuit(model,cfg);
assert(~isempty(found),'Unanchored synthetic topology was not found.');
assert(isequal(sign(full(model.S(found(1).metIndex,found(1).reactionIndex))),p.T));
fprintf('PASS: unanchored full topology.\n');

cfg = anchorMetabolites(cfg,p);
expect(model,cfg,true,'anchored full topology');

bad = model; bad.rxnECNumbers{2} = bad.rxnECNumbers{1};
expect(bad,cfg,false,'shared EC classification rejected');
bad = model; bad.grRules{2} = bad.grRules{1};
expect(bad,cfg,false,'shared gene rejected even with different ECs');
bad = model; bad.grRules{1} = '(G1 or G2)';
expect(bad,cfg,false,'overlapping alternative GPR genes rejected');

bad = model; bad.rxnECNumbers{1} = '';
expect(bad,cfg,false,'missing enzymatic annotation rejected');
relaxed = cfg; relaxed.ENZYME_EVIDENCE = 'EC_OR_GPR';
expect(bad,relaxed,true,'explicit GPR evidence mode');
bad = model; bad.rxnECNumbers{1} = '1.1.1.-';
expect(bad,cfg,false,'incomplete EC annotation rejected');
bad = model; bad.rxnNames{1} = 'spontaneous conversion';
expect(bad,cfg,false,'spontaneous annotation rejected');

bad = model; bad.ub(1) = 0;
expect(bad,cfg,false,'closed reaction rejected');
bad = model; bad.S(:,1) = -bad.S(:,1); bad.lb(1) = -1000; bad.ub(1) = 0;
expect(bad,cfg,false,'reverse-only reaction rejected in forward mode');
reverse = cfg; reverse.FORWARD_ONLY = false;
[found,~,~] = findMetabolicCircuit(bad,reverse);
assert(~isempty(found) && found(1).direction(1) == -1);
fprintf('PASS: allowed reverse orientation accepted.\n');
bad.lb(1) = 0; bad.ub(1) = 1000;
expect(bad,reverse,false,'reverse orientation forbidden by bounds rejected');

bad = model; bad.S(13,1) = 1;
expect(bad,cfg,false,'extra participation by a named role rejected');
badCfg = cfg; badCfg.FIXED_METABOLITES.m_Sf = 'm_Ff[c]';
expect(model,badCfg,false,'same metabolite assigned to two roles rejected');
badCfg = cfg;
badCfg.FIXED_REACTIONS.r_Ff = 'synthetic_01';
badCfg.FIXED_REACTIONS.r_Sf = 'synthetic_01';
expect(model,badCfg,false,'same original reaction assigned twice rejected');

cofactor = model;
cofactor.S(end+1,:) = -1;
cofactor.mets{end+1} = 'cofactor[c]';
cofactor.metNames{end+1} = 'Synthetic shared cofactor';
cofactor.metFormulas{end+1} = '';
cofactor.metComps(end+1) = 1;
expect(cofactor,cfg,false,'shared auxiliary rejected in strict mode');
rolesOnly = cfg; rolesOnly.METABOLITE_SCOPE = 'roles';
expect(cofactor,rolesOnly,true,'shared auxiliary allowed in roles mode');
exempt = cfg; exempt.AUXILIARY_REUSE_ALLOWED = {'cofactor'};
expect(cofactor,exempt,true,'explicit auxiliary exemption');

curated = model;
curated.rxnECNumbers(:) = {'1.1.1.1'};
curated.rxnEnzymeIDs = arrayfun(@(x) sprintf('protein_%02d',x),(1:18)', ...
    'UniformOutput',false);
customCfg = cfg; customCfg.ENZYME_ID_FIELD = 'rxnEnzymeIDs';
expect(curated,customCfg,true,'curated proteins distinguish identical EC activities');
curated.rxnEnzymeIDs{2} = curated.rxnEnzymeIDs{1};
expect(curated,customCfg,false,'shared curated enzyme rejected');

% Distinct compartments cannot disguise the same labelled chemical.
alias = model;
alias.S(end+1,:) = alias.S(3,:); alias.S(3,:) = 0;
alias.mets{end+1} = 'm_Ff[m]';
alias.metNames{end+1} = 'Alias of m_Ff in another compartment';
alias.metFormulas{end+1} = '';
alias.metComps(end+1) = 2; alias.comps{2} = 'm';
aliasCfg = cfg; aliasCfg.FIXED_METABOLITES.m_Sf = 'm_Ff[m]';
expect(alias,aliasCfg,false,'chemical alias cannot fill two named roles');
aliasCfg.METABOLITE_IDENTITY = 'model_id';
expect(alias,aliasCfg,true,'explicit model-ID identity mode');

% Activation uses supplied evidence and other substrates for r_Nw.
activated = model; activated.S(10,16) = 0;
activated.S(end+1,16) = -1;
activated.mets{end+1} = 'neuron_feed[c]';
activated.metNames{end+1} = 'Synthetic neuron substrate';
activated.metFormulas{end+1} = '';
activated.metComps(end+1) = 1;
activationCfg = cfg; activationCfg.NEURON_INPUT_MODE = 'activation';
mustError(@() findMetabolicCircuit(activated,activationCfg),'Circuit:ActivationEvidence');
activationCfg.ACTIVATION_PAIRS = {'synthetic_16','m_c[c]','synthetic test only'};
expect(activated,activationCfg,true,'supplied activation evidence accepted');

% Exercise text output, coefficients, no-hit output, and MAT discovery.
testFolder = tempname; mkdir(testFolder);
cleanup = onCleanup(@() rmdir(testFolder,'s')); %#ok<NASGU>
scaled = model; scaled.S(1,1) = -2.5;
[found,info,pattern] = findMetabolicCircuit(scaled,cfg);
reportFile = fullfile(testFolder,'report.txt');
writeMetabolicCircuitReport(reportFile,scaled,found,info,pattern,cfg);
report = fileread(reportFile);
assert(contains(report,'2.5 m_in[c](m_in)'));
assert(contains(report,'r_Ff_2') && contains(report,'r_Nr') && contains(report,'m_out'));
bad = model; bad.rxnECNumbers{1} = '';
[found,info,pattern] = findMetabolicCircuit(bad,cfg);
writeMetabolicCircuitReport(reportFile,bad,found,info,pattern,cfg);
assert(contains(fileread(reportFile),'NO COMPLETE MOTIF RETURNED'));
modelFile = fullfile(testFolder,'synthetic.mat');
save(modelFile,'model');
[loaded,~] = loadCircuitModel(modelFile);
assert(isequal(loaded.S,model.S));
secondModel = model; %#ok<NASGU>
save(modelFile,'model','secondModel');
mustError(@() loadCircuitModel(modelFile),'Circuit:MATVariable');
[loaded,~] = loadCircuitModel(modelFile,'model');
assert(isequal(loaded.S,model.S));
fprintf('PASS: text reports, stoichiometric coefficients, and MAT loader.\n');

if runParallelTest
    % Two disconnected copies guarantee a nontrivial root partition even
    % if propagation completely solves a single anchored copy.
    twin = model;
    twin.S = blkdiag(model.S,model.S);
    for field = {'mets','metNames','rxns','rxnNames','grRules'}
        f = field{1};
        twin.(f) = [model.(f);cellfun(@(x) ['copy_',x],model.(f),'UniformOutput',false)];
    end
    twin.metFormulas = [model.metFormulas;model.metFormulas];
    twin.metComps = [model.metComps;model.metComps];
    twin.lb = [model.lb;model.lb]; twin.ub = [model.ub;model.ub];
    twin.rxnECNumbers = [model.rxnECNumbers;arrayfun(@(j) ...
        sprintf('2.1.1.%d',j),(1:18)','UniformOutput',false)];
    parallelCfg = cfg; parallelCfg.USE_PARALLEL = true;
    parallelCfg.FIXED_METABOLITES = struct();
    parallelCfg.NUM_WORKERS = 2;
    parallelCfg.MAX_MOTIFS = Inf;
    serialCfg = parallelCfg; serialCfg.USE_PARALLEL = false;
    serial = findMetabolicCircuit(twin,serialCfg);
    [parallel,info,~] = findMetabolicCircuit(twin,parallelCfg);
    assert(info.parallelUsed,'No parallel pool was used.');
    assert(isequal(sortrows(signatures(serial)),sortrows(signatures(parallel))));
    fprintf('PASS: serial/parallel exhaustive result equality.\n');
end
fprintf('All requested MATLAB self-tests passed.\n');
end

function [model,cfg,p] = fixture()
p = circuitPattern('substrate');
S = p.T;
mets = cellfun(@(x) [x,'[c]'],p.metNames,'UniformOutput',false);
for j = 1:18
    if ~any(S(:,j) > 0)
        S(end+1,j) = 1; %#ok<AGROW>
        mets{end+1,1} = sprintf('waste_%02d[c]',j); %#ok<AGROW>
    end
end
model.S = sparse(S); model.mets = mets;
model.metNames = mets; model.metFormulas = repmat({''},numel(mets),1);
model.rxns = arrayfun(@(j) sprintf('synthetic_%02d',j),(1:18)','UniformOutput',false);
model.rxnNames = model.rxns;
model.rxnECNumbers = arrayfun(@(j) sprintf('1.1.1.%d',j),(1:18)','UniformOutput',false);
model.grRules = arrayfun(@(j) sprintf('G%d',j),(1:18)','UniformOutput',false);
model.lb = zeros(18,1); model.ub = 1000*ones(18,1);
model.comps = {'c'}; model.metComps = ones(numel(mets),1);
cfg = circuitSearchDefaults(); cfg.VERBOSE = false;
cfg.ONE_MOTIF_PER_INPUT = false; % Preserve checks for complete circuit enumeration.
cfg.METABOLITE_SCOPE = 'all'; % Explicitly exercise the stricter optional mode.
cfg.MAX_MOTIFS = 1; cfg.MAX_SECONDS = 120;
cfg.MODEL_FILE = 'SYNTHETIC TEST NETWORK (not cell chemistry)';
end
function cfg = anchorMetabolites(cfg,p)
for i = 1:13, cfg.FIXED_METABOLITES.(p.metNames{i}) = [p.metNames{i},'[c]']; end
end
function expect(model,cfg,wanted,label)
[found,info,~] = findMetabolicCircuit(model,cfg);
assert(strcmp(info.status,'exhausted') || ~isempty(found), ...
    'Search stopped before this test could establish a result: %s',label);
assert(~isempty(found) == wanted,'Unexpected result: %s',label);
fprintf('PASS: %s.\n',label);
end
function mustError(f,identifier)
try
    f();
catch e
    assert(strcmp(e.identifier,identifier),'Unexpected error: %s',e.message); return;
end
error('Circuit:SelfTest','Expected error %s was not raised.',identifier);
end
function rows = signatures(motifs)
rows = zeros(numel(motifs),49);
for k = 1:numel(motifs)
    rows(k,:) = [motifs(k).metIndex;motifs(k).reactionIndex;motifs(k).direction]';
end
end

