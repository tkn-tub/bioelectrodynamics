function test_polarity_pipeline()
% Native MATLAB integration/regression tests, without COBRA or a real model.
addpath(fileparts(mfilename('fullpath')));
temporary = tempname; mkdir(temporary); previous = pwd;
cleanup = onCleanup(@() localCleanup(previous,temporary)); cd(temporary);
[model,isBadMet,InputAdapterPaths,up,options] = fixture();
[counter,~,~,info,none] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,up,3,options);
assert(height(counter) == 1 && counter.m_c == 10 && counter.m_Cs == 11);
assert(counter.r_Cp == 10 && counter.r_Cg == 11 && counter.r_Cr == 12);
assert(~ismember('r_CsSupply',counter.Properties.VariableNames));
assert(~ismember('r_CsLeak',counter.Properties.VariableNames));
[full,diag,~,~,none] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,counter,4,options);
assert(height(full) == 1 && full.m_Ni == 14 && full.m_Na == 16 && full.m_Nb == 18 && full.m_out == 20);
row = table2struct(full); assert(checkMinimalCircuit(model,row,4));
schema = minimalCircuitSchema(4);
assert(numel(schema.rxns) == 19 && numel(schema.mets) == 11);
assert(all(strcmp(schema.regulators(1:2),'m_d2')) && ~any(strcmp(schema.regulators,'m_Y1')));
assert(row.m_d == row.m_Y1 && row.m_d2 == row.m_Y2);
assert(~full.RegulatoryEvidenceComplete && ~full.KineticsValidated);

% Optional basal Cs removal is a real fourth reaction, with no mapped supply.
withLeak = options; withLeak.mapCsLeak = true;
[counterLeak,~,~,~,~] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,up,3,withLeak);
assert(height(counterLeak) == 1 && counterLeak.r_CsLeak == 21);
[fullLeak,~,~,~,~] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,counterLeak,4,withLeak);
assert(height(fullLeak) == 1 && checkMinimalCircuit(model,table2struct(fullLeak),4));
leakSchema = minimalCircuitSchema(4,true);
assert(numel(leakSchema.rxns) == 20);

% Missing individual interactions are not invented by evidence mode.
strict = options; strict.regulationMode = "evidence";
[emptyCounter,~,~,~,~] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,up,3,strict);
assert(isempty(emptyCounter));
met = [6;6;10;14;14;16;18]; rxn = [10;12;13;15;17;19;19];
effect = ["inhibition";"inhibition";repmat("activation",5,1)];
raw = table(string(model.mets(met)),string(model.rxns(rxn)),ones(7,1),effect, ...
    repmat("synthetic fixture, not biological evidence",7,1), ...
    'VariableNames',{'RegulatorMetID','TargetRxnID','TargetDirection','Effect','Evidence'});
writetable(raw,'test_regulation.csv'); strict.regulationFile = 'test_regulation.csv';
[verifiedCounter,~,~,~,evidence] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,up,3,strict);
[verifiedFull,~,~,~,~] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,verifiedCounter,4,strict);
assert(height(verifiedFull) == 1 && verifiedFull.RegulatoryEvidenceComplete && ~verifiedFull.KineticsValidated);
missing = model; missing.eccodes{13} = '';
[noNN,~,~,~,~] = extendMinimalCircuit(missing,isBadMet,InputAdapterPaths,verifiedCounter,4,strict);
assert(isempty(noNN));

% Sign handling and cross-talk account for the updated second differentiator.
reverse = model; reverse.S(:,11) = -reverse.S(:,11); reverse.lb(11) = -1000; reverse.ub(11) = 0;
[reversed,~,~,~,~] = extendMinimalCircuit(reverse,isBadMet,InputAdapterPaths,up,3,options);
assert(height(reversed) == 1 && reversed.dir_Cg == -1);
cross = model; cross.S(row.m_Y1,row.r_Cp) = -1;
assert(~checkMinimalCircuit(cross,row,4));
external = model; external.S(:,22) = sparse([9 4],[1 1],[-1 1],size(model.S,1),1);
assert(checkMinimalCircuit(external,row,4));
[membrane,~] = classifyMinimalMembrane(model); assert(~any(membrane));
transport = model; transport.S(22,row.r_No) = -1;
[membrane,~] = classifyMinimalMembrane(transport); assert(membrane(row.r_No));

% Current Block 2 loads unchanged; legacy version tags and fake schema fail.
payload = polarityDifferentiatorCheckpoint(model,up,table(),InputAdapterPaths,struct('stage',2),struct());
payload.Motifs_Diff_Block_2 = up;
saveCompactCheckpoint('Block_2_differentiator.mat',payload);
loaded = loadMinimalCheckpoint('Block_2_differentiator.mat',model,'Motifs_Diff_Block_2');
assert(isequal(loaded.Motifs_Diff_Block_2,up));
old = payload; old.TopologyVersion = "minimal_three_reaction_v1";
saveCompactCheckpoint('old.mat',old);
assertError(@() loadMinimalCheckpoint('old.mat',model,'Motifs_Diff_Block_2'),'MinimalLoad:OldTopology');
fake = payload; fake.Motifs_Diff_Block_2(:,{'r_R2'}) = [];
saveCompactCheckpoint('fake.mat',fake);
assertError(@() loadMinimalCheckpoint('fake.mat',model,'Motifs_Diff_Block_2'),'MinimalLoad:Schema');
wrong = model; wrong.mets{2} = 'another_input[c]';
assertError(@() loadMinimalCheckpoint('Block_2_differentiator.mat',wrong,'Motifs_Diff_Block_2'),'MinimalLoad:ModelMismatch');

% Compact full checkpoint and the detailed report retain current reaction roles.
saveMinimalCircuit('saved.mat','Motifs_NN_Block_6',full,model,InputAdapterPaths,options,info,diag,none);
saved = loadMinimalCheckpoint('saved.mat',model,'Motifs_NN_Block_6');
assert(~isfield(saved,'model') && ~isfield(saved,'results_nn'));
assert(~saved.CounterCsLeakMapped && contains(saved.CounterSupportAssumption,'complete net Cs'));
saveMinimalCircuit('leak.mat','Motifs_NN_Block_6',fullLeak,model,InputAdapterPaths,withLeak,info,diag,none);
savedLeak = loadMinimalCheckpoint('leak.mat',model,'Motifs_NN_Block_6');
assert(savedLeak.CounterCsLeakMapped);

% Continuation runs Blocks 3-7 from the existing Block 2; no solver initialization.
REQUIRE_GPR = false; REGULATION_MODE = "design"; REGULATION_FILE = '';
MAX_HITS_COUNTER = 1; MAX_HITS_NN = 1; MAX_COUNTER_POOL = 100; MAX_NN_POOL = 100;
MAP_COUNTER_CS_LEAK = false; SHOW_PROGRESS = false;
FILTER_INTERNAL_CROSSTALK_DURING_SEARCH = true;
PRINT_ALLOW_XML_LOAD = false; PRINT_SAVE_MODEL_CACHE = false;
MOTIF_TO_PRINT = 1;
Resume_Polarity_Pipeline;
assert(height(Motifs_NN_Block_6) == 1 && isfile('Block_7_Selected_Motifs.txt'));
report = fileread('Block_7_Selected_Motifs.txt');
assert(contains(report,'r_S1_2') && contains(report,'r_R2') && contains(report,'joint M2 and Y2 removal'));
assert(~contains(report,'r_CsSupply') && ~contains(report,'r_21'));

% Reusing current caches is allowed only if they refer to the loaded parents.
RUN_BLOCK3 = false; RUN_BLOCK4 = false; RUN_BLOCK5 = false; RUN_BLOCK6 = false;
Continue_Minimal_Pipeline;
changed = load('Block_2_differentiator.mat');
changed.Motifs_Diff_Block_2.adapterPathIndex = 2;
secondPath = changed.InputAdapterPaths; secondPath.pathIndex = 2;
changed.InputAdapterPaths = [changed.InputAdapterPaths;secondPath];
saveCompactCheckpoint('Block_2_differentiator.mat',changed);
assertError(@() localResume(model,isBadMet),'MinimalContinue:CachedParent');

% A zero-survivor checkpoint remains printable, without requesting a model.
Motifs_NN_Block_4 = Motifs_NN_Block_4([],:);
Block_5_Remove_Cell_Membrane_Reactions;
Block_6_Remove_Internal_CrossTalk;
Block_7_Print_Selected_Motif;
assert(isempty(Motifs_NN_Block_6));
fprintf('Polarity continuation regression tests passed.\n');
clear cleanup;
end

function localResume(model,isBadMet) %#ok<INUSD>
RUN_BLOCK3 = false; RUN_BLOCK4 = false; RUN_BLOCK5 = false; RUN_BLOCK6 = false;
REQUIRE_GPR = false; SHOW_PROGRESS = false; PRINT_SELECTED_MOTIF = false;
Continue_Minimal_Pipeline;
end

function [model,bad,paths,up,o] = fixture()
ids = {'camp[c]';'X[c]';'M1[c]';'Y1[c]';'M2[c]';'Y2[c]';'W1[c]';'W2[c]'; ...
    'PC[c]';'C[c]';'Cs[c]';'WC[c]';'PNi[c]';'Ni[c]';'PNa[c]';'Na[c]'; ...
    'PNb[c]';'Nb[c]';'PNo[c]';'Nout[c]';'WN[c]';'extra[e]'};
sub = {1,2,3,[2 4],4,5,4,[5 6],6,9,[10 11],11,13,14,15,16,17,18,19,20,11};
prod = {2,3,4,7,5,8,6,8,8,10,12,12,14,21,16,21,18,21,20,21,21};
S = sparse(numel(ids),numel(sub));
for k = 1:numel(sub), S(sub{k},k) = -1; S(prod{k},k) = 1; end
model = struct('S',S,'mets',{ids},'metNames',{ids});
model.rxns = cellstr("r"+string((1:21).')); model.rxnNames = model.rxns;
model.eccodes = [{''};repmat({'1.1.1.1'},20,1)]; model.grRules = repmat({'gene'},21,1);
model.lb = zeros(21,1); model.ub = repmat(1000,21,1);
bad = true(numel(ids),1); bad([2:6,10,11,14,16,18,20]) = false;
paths = struct('pathIndex',1,'anchorMetIdx',2,'anchorMetID',"X[c]", ...
    'inputMetIdx',1,'inputMetID',"camp[c]",'reactionIndices',1, ...
    'reactionPath',"r1",'metabolitePath',"camp[c] -> X[c]");
s = struct('stage',1,'maxHits',1,'maxInputs',inf,'perInputCap',100, ...
    'maxChecks',10000,'requireGPR',false,'requireAdapterEC',false, ...
    'showProgress',false,'useParallel',false,'autoPool',false,'workers',[], ...
    'distinctAdapter',true,'distinctStageRxns',true,'distinctStageMets',true,'crossTalk',true);
[one,~,~,~] = mapPolarityDifferentiators(model,bad,paths,table(),s);
s.stage = 2;
[up,~,~,~] = mapPolarityDifferentiators(model,bad,paths,one,s);
assert(height(up) == 1);
o = struct('maxHits',1,'poolCap',100,'maxChecks',100000,'maxJoinChecks',100000, ...
    'requireGPR',false,'requireAdapterEC',false,'regulationMode',"design",'regulationFile','', ...
    'optionsCap',50,'nodeCap',100,'nodesPerRole',40,'allowCurrencyCs',false, ...
    'distinctAdapter',true,'crossTalk',true,'showProgress',false,'mapCsLeak',false);
end

function assertError(action,identifier)
try
    action(); error('PolarityTest:ExpectedError','Expected %s.',identifier);
catch failure
    assert(strcmp(failure.identifier,identifier),failure.message);
end
end

function localCleanup(previous,temporary)
cd(previous);
if isfolder(temporary), rmdir(temporary,'s'); end
end
