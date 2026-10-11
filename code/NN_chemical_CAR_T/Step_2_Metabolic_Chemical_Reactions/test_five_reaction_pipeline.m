function test_five_reaction_pipeline()
% Synthetic native MATLAB integration test: mapping -> counter/NN -> TXT.
addpath(fileparts(mfilename('fullpath')));
temporary = tempname; mkdir(temporary); previous = pwd;
cleanup = onCleanup(@() localCleanup(previous,temporary)); cd(temporary);
ids = {'camp[c]';'X[c]';'M1[c]';'Y1[c]';'F1[c]';'M2[c]';'Y2[c]';'W1[c]';'W2[c]'; ...
    'PC[c]';'C[c]';'Cs[c]';'WC[c]';'PNi[c]';'Ni[c]';'PNa[c]';'Na[c]'; ...
    'PNb[c]';'Nb[c]';'PNo[c]';'Nout[c]';'WN[c]';'extra[e]'};
left = {1,2,3,2,[5 4],5,4,6,4,[6 7],7,10,[11 12],12,14,15,16,17,18,19,20,21,12};
right = {2,3,4,5,8,8,6,7,7,9,9,11,13,13,15,22,17,22,19,22,21,22,22};
S = sparse(numel(ids),numel(left));
for k = 1:numel(left), S(left{k},k) = -1; S(right{k},k) = 1; end
model = struct('S',S,'mets',{ids},'metNames',{ids});
model.rxns = cellstr("r"+string((1:23).')); model.rxnNames = model.rxns;
model.eccodes = [{''};repmat({'1.1.1.1'},22,1)]; model.grRules = repmat({'gene'},23,1);
model.lb = zeros(23,1); model.ub = repmat(1000,23,1);
isBadMet = true(numel(ids),1); isBadMet([2:7,11,12,15,17,19,21]) = false;
InputAdapterPaths = struct('pathIndex',1,'anchorMetIdx',2,'anchorMetID',"X[c]", ...
    'inputMetIdx',1,'inputMetID',"camp[c]",'reactionIndices',1, ...
    'reactionPath',"r1",'metabolitePath',"camp[c] -> X[c]");
REQUIRE_GPR = false; SHOW_PROGRESS = false; USE_PARALLEL = false; USE_PARALLEL_BLOCK2 = false;
MAX_HITS = 1; MAX_HITS_BLOCK2 = 1; MAX_HITS_COUNTER = 1; MAX_HITS_NN = 1;
REGULATION_MODE = "design"; REGULATION_FILE = ''; MAP_COUNTER_CS_LEAK = false;
MAX_COUNTER_POOL = 100; MAX_NN_POOL = 100;
PRINT_ALLOW_XML_LOAD = false; PRINT_SAVE_MODEL_CACHE = false;
BLOCK1_OUTPUT_DIR = pwd; BLOCK2_OUTPUT_DIR = pwd; MOTIF_TO_PRINT = 1;
Block_1_Differentiator_Mapping;
Block_2_Differentiator_Mapping;
Resume_Five_Reaction_Pipeline;
assert(height(Motifs_NN_Block_6) == 1);
row = table2struct(Motifs_NN_Block_6);
assert(checkMinimalCircuit(model,row,4));
assert(row.r_Cp == 12 && row.r_Cg == 13 && row.r_Cr == 14);
schema = minimalCircuitSchema(4);
assert(numel(schema.rxns) == 21 && numel(schema.mets) == 12);
assert(all(strcmp(schema.regulators(1:2),'m_d2')));
assert(~ismember('r_CsLeak',Motifs_NN_Block_6.Properties.VariableNames));
assert(~Motifs_NN_Block_6.KineticsValidated && ~Motifs_NN_Block_6.RegulatoryEvidenceComplete);
saved = loadMinimalCheckpoint('Block_6_no_internal_crosstalk.mat',model,'Motifs_NN_Block_6');
assert(~isfield(saved,'model') && ~isfield(saved,'fullJoinRows'));
report = fileread('Block_7_Selected_Motifs.txt');
assert(contains(report,'m_F1_delay') && contains(report,'joint F1 and Y1 removal'));
assert(contains(report,'M2 to Y2') && ~contains(report,'r_CsSupply'));
old = saved; old.TopologyVersion = "polarity_regulated_full_v1";
saveCompactCheckpoint('old.mat',old);
assertError(@() loadMinimalCheckpoint('old.mat',model,'Motifs_NN_Block_6'),'MinimalLoad:OldTopology');

% Cache ancestry is verified, and optional CsLeak keeps all ten diff roles.
RUN_BLOCK3 = false; RUN_BLOCK4 = false; RUN_BLOCK5 = false; RUN_BLOCK6 = false;
Continue_Minimal_Pipeline;
MAP_COUNTER_CS_LEAK = true;
Resume_Five_Reaction_Pipeline;
assert(height(Motifs_NN_Block_6) == 1 && Motifs_NN_Block_6.r_CsLeak == 23);
leakSchema = minimalCircuitSchema(4,true); assert(numel(leakSchema.rxns) == 22);
assert(checkMinimalCircuit(model,table2struct(Motifs_NN_Block_6),4));

% Cross-talk within the algorithm is rejected; outside producers are allowed.
cross = model; cross.S(row.m_F1_delay,row.r_Cp) = 1;
assert(~checkMinimalCircuit(cross,row,4));
external = model; external.S(:,24) = sparse([10 4],[1 1],[-1 1],size(S,1),1);
assert(checkMinimalCircuit(external,row,4));
transport = model; transport.S(23,row.r_No) = -1;
[membrane,~] = classifyMinimalMembrane(transport); assert(membrane(row.r_No));
fprintf('Five-reaction full pipeline regression tests passed.\n');
clear cleanup;
end

function assertError(action,identifier)
try
    action(); error('FiveReactionTest:ExpectedError','Expected %s.',identifier);
catch failure
    assert(strcmp(failure.identifier,identifier),failure.message);
end
end

function localCleanup(previous,temporary)
cd(previous);
if isfolder(temporary), rmdir(temporary,'s'); end
end
