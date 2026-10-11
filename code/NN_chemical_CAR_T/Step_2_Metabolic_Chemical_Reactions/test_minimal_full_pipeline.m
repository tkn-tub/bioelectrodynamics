function test_minimal_full_pipeline()
% Synthetic regression tests only; no COBRA/model download is required.
folder = fileparts(mfilename('fullpath')); addpath(folder);
temporary = tempname; mkdir(temporary);
previous = pwd; cleanup = onCleanup(@() localCleanup(previous,temporary)); cd(temporary);
[model,isBadMet,InputAdapterPaths,up,options] = fixture();
[counter,~,~,info,none] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,up,3,options);
assert(height(counter) == 1 && counter.m_c == 10 && counter.m_Cs == 11);
assert(~counter.RegulatoryEvidenceComplete && ~counter.KineticsValidated && info.hitLimitReached);
assert(counter.r_Cp == 8 && counter.r_Cg == 9 && counter.r_Cr == 10);
assert(counter.r_CsSupply == 11 && counter.r_CsLeak == 12);
[full,diag,~,~,none] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,counter,4,options);
assert(height(full) == 1 && full.m_Ni == 15 && full.m_Na == 17 && full.m_Nb == 19 && full.m_out == 21);
row = table2struct(full);
assert(checkMinimalCircuit(model,row,4));
schema = minimalCircuitSchema(4);
assert(numel(schema.regulators) == 7 && all(strcmp(schema.regulators(1:2),'m_d2')));
assert(~any(strcmp(schema.regulators,'m_Y1')) && row.m_d2 == row.m_Y2);
assert(~full.RegulatoryEvidenceComplete && ~full.KineticsValidated);
assert(model.S(row.m_c,row.r_Ni) == 0 && model.S(row.m_Ni,row.r_Na) == 0);

% No regulatory interactions are fabricated in evidence mode.
strict = options; strict.regulationMode = "evidence";
[emptyCounter,~,~,~,~] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,up,3,strict);
assert(isempty(emptyCounter));

% Supply seven individual synthetic links; Y1 does not regulate the counter.
met = [6;6;10;15;15;17;19]; rxn = [8;10;13;15;17;19;19];
effect = ["inhibition";"inhibition";repmat("activation",5,1)];
raw = table(string(model.mets(met)),string(model.rxns(rxn)),ones(7,1),effect, ...
    repmat("synthetic fixture; not biological evidence",7,1),'VariableNames', ...
    {'RegulatorMetID','TargetRxnID','TargetDirection','Effect','Evidence'});
writetable(raw,'test_regulation.csv'); strict.regulationFile = 'test_regulation.csv';
[counterVerified,~,~,~,evidence] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,up,3,strict);
[fullVerified,~,~,~,~] = extendMinimalCircuit(model,isBadMet,InputAdapterPaths,counterVerified,4,strict);
assert(height(fullVerified) == 1 && fullVerified.RegulatoryEvidenceComplete);
assert(~fullVerified.KineticsValidated); % Individual links do not prove gate/Hill kinetics.

% A missing EC cannot fill a hit budget. The passive adapter remains exempt.
missing = model; missing.eccodes{13} = '';
[noNN,~,~,~,~] = extendMinimalCircuit(missing,isBadMet,InputAdapterPaths,counterVerified,4,strict);
assert(isempty(noNN));

% Signed reverse gate, bounded flux, and selected cross-talk.
reverse = model; reverse.S(:,9) = -reverse.S(:,9); reverse.lb(9) = -1000; reverse.ub(9) = 0;
[counterReverse,~,~,~,~] = extendMinimalCircuit(reverse,isBadMet,InputAdapterPaths,up,3,options);
assert(height(counterReverse) == 1 && counterReverse.dir_Cg == -1);
cross = model; cross.S(row.m_Y1,row.r_Cp) = -1;
[ok,~] = checkMinimalCircuit(cross,row,4); assert(~ok);
external = model; external.S(:,21) = sparse([9 4],[1 1],[-1 1],size(model.S,1),1);
external.rxns{21} = 'outside'; external.rxnNames{21} = 'outside';
external.eccodes{21} = ''; external.grRules{21} = '';
external.lb(21) = 0; external.ub(21) = 1000;
assert(checkMinimalCircuit(external,row,4)); % Other network producers are allowed.

% One-row membrane classification must retain a ROW, not 19 column entries.
[membrane,~] = classifyMinimalMembrane(model); assert(~any(membrane));
transport = model; transport.S(23,row.r_No) = -1;
[membrane,~] = classifyMinimalMembrane(transport); assert(membrane(row.r_No));
internal = model; internal.S(24,row.r_No) = -1;
internal.mets{24} = 'internal[m]'; internal.metNames{24} = 'internal';
[membrane,~] = classifyMinimalMembrane(internal); assert(~membrane(row.r_No));

% Compact round-trip and a printer requiring neither model nor COBRA.
saveMinimalCircuit('saved.mat','Motifs_NN_Block_6',full,model,InputAdapterPaths,options,info,diag,none);
saved = load('saved.mat'); assert(~isfield(saved,'model') && ~isfield(saved,'results_nn'));
report = renderMinimalFullMotif(saved,1);
assert(contains(report,'PENDING') && contains(report,'r_CsSupply') && contains(report,'r_NoDecay'));
saveMinimalCircuit('verified.mat','Motifs_NN_Block_6',fullVerified,model,InputAdapterPaths,strict,info,diag,evidence);
verified = load('verified.mat'); report = renderMinimalFullMotif(verified,1);
assert(~contains(report,'PENDING') && contains(report,'Kinetics validated: 0'));

% Script-level smoke test, including the new filters and selected TXT export.
Motifs_Diff_Block_2 = up;
REQUIRE_GPR = false; REGULATION_MODE = "design"; REGULATION_FILE = '';
MAX_HITS_COUNTER = 1; MAX_HITS_NN = 1;
MAX_COUNTER_POOL = 100; MAX_NN_POOL = 100;
SHOW_PROGRESS = false; FILTER_INTERNAL_CROSSTALK_DURING_SEARCH = true;
Block_3_Counter_Mapping;
Block_4_Hill_NN_Mapping;
Block_5_Remove_Cell_Membrane_Reactions;
assert(height(Motifs_NN_Block_5) == 1);
Block_6_Remove_Internal_CrossTalk;
assert(height(Motifs_NN_Block_6) == 1 && ~Motifs_NN_Block_6.RegulatoryEvidenceComplete);
MOTIF_TO_PRINT = 1; Block_7_Print_Selected_Motif;
assert(isfile('Full_Motif_1.txt'));

% No survivor is a valid empty checkpoint, not a Block 7 exception.
Motifs_NN_Block_4 = Motifs_NN_Block_4([],:);
Block_5_Remove_Cell_Membrane_Reactions;
Block_6_Remove_Internal_CrossTalk;
Block_7_Print_Selected_Motif;
assert(isempty(Motifs_NN_Block_6));
fprintf('Full minimal-pipeline regression tests passed.\n');
clear cleanup;
end

function [model,bad,paths,up,o] = fixture()
ids = {'camp[c]';'X1[c]';'M1[c]';'Y1[c]';'M2[c]';'Y2[c]';'W1[c]';'W2[c]'; ...
    'PC[c]';'C[c]';'Cs[c]';'PCs[c]';'WC[c]';'PNi[c]';'Ni[c]';'PNa[c]';'Na[c]'; ...
    'PNb[c]';'Nb[c]';'PNo[c]';'Nout[c]';'WN[c]';'extra[e]'};
sub = {1,2,3,[2 4],4,5,[4 6],9,[10 11],11,12,11,14,15,16,17,18,19,20,21};
prod = {2,3,4,7,5,6,8,10,13,13,11,13,15,22,17,22,19,22,21,22};
net = sparse(numel(ids),numel(sub));
for k = 1:numel(sub), net(sub{k},k) = -1; net(prod{k},k) = 1; end
model = struct;
model.S = net; model.mets = ids; model.metNames = ids;
model.rxns = cellstr("r"+string((1:20).')); model.rxnNames = model.rxns;
model.eccodes = [{''};repmat({'1.1.1.1'},19,1)];
model.grRules = repmat({'gene'},20,1);
model.lb = zeros(20,1); model.ub = ones(20,1)*1000;
bad = true(numel(ids),1); bad([2:6,10,11,15,17,19,21]) = false;
paths = struct('pathIndex',1,'anchorMetIdx',2,'anchorMetID',"X1[c]", ...
    'inputMetIdx',1,'inputMetID',"camp[c]",'reactionIndices',1, ...
    'reactionPath',"r1",'metabolitePath',"camp[c] -> X1[c]");
settings = struct('stage',1,'maxHits',1,'maxInputs',inf,'perInputCap',100, ...
    'maxChecks',10000,'requireGPR',false,'requireAdapterEC',false, ...
    'showProgress',false,'useParallel',false,'autoPool',false,'workers',[], ...
    'distinctAdapter',true,'distinctStageRxns',true,'distinctStageMets',true,'crossTalk',true);
[one,~,~,~] = mapMinimalDifferentiator(model,bad,paths,table(),settings);
settings.stage = 2;
[up,~,~,~] = mapMinimalDifferentiator(model,bad,paths,one,settings);
assert(height(up) == 1);
o = struct('maxHits',1,'poolCap',100,'maxChecks',100000,'maxJoinChecks',100000, ...
    'requireGPR',false,'requireAdapterEC',false,'regulationMode',"design", ...
    'regulationFile','','optionsCap',50,'nodeCap',100,'nodesPerRole',40, ...
    'allowCurrencyCs',false,'distinctAdapter',true,'crossTalk',true,'showProgress',false);
end

function localCleanup(previous,temporary)
cd(previous);
if isfolder(temporary), rmdir(temporary,'s'); end
end
