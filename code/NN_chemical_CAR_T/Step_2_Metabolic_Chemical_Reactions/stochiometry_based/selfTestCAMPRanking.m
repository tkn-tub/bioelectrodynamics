function selfTestCAMPRanking()
%SELFTESTCAMPRANKING Synthetic checks for directed cAMP input ranking.
% Run from this folder: selfTestCAMPRanking
% These tests establish graph-search behavior, not biological feasibility.
[model,motif,cfg,p] = fixture();

exact = model;
exact.mets{1} = 'camp[c]'; exact.metNames{1} = '3'',5''-cyclic AMP';
exact.mets{14} = 'other[c]'; exact.metNames{14} = 'Other compound';
[ranking,paths] = rankCircuitInputsByCAMP(exact,motif,cfg,p);
assert(ranking.Distance(1) == 0 && strcmp(paths(1).status,'exact_input'));
assert(isequal(paths(1).metIndex(:),1) && isempty(paths(1).reactionIndex));
fprintf('PASS: exact cAMP input has distance zero.\n');

one = addBridge(model,14,1);
[ranking,paths] = rankCircuitInputsByCAMP(one,motif,cfg,p);
assert(ranking.Distance(1) == 1);
assert(isequal(paths(1).metIndex(:),[14;1]));
assert(isequal(paths(1).direction(:),1));
namedSource = one; namedSource.mets{14} = 'opaque_source[c]';
[ranking,~] = rankCircuitInputsByCAMP(namedSource,motif,cfg,p);
assert(ranking.Distance(1) == 1,'A precise 3'',5''-cAMP name must identify an opaque ID.');
fprintf('PASS: one directed enzymatic conversion.\n');

two = addBridge(addBridge(model,14,15),15,1);
[ranking,paths] = rankCircuitInputsByCAMP(two,motif,cfg,p);
assert(ranking.Distance(1) == 2 && numel(paths(1).reactionIndex) == 2);
assert(isequal(paths(1).metIndex(:),[14;15;1]));
assert(isempty(intersect(paths(1).reactionIndex,motif.reactionIndex)), ...
    'Selected circuit reactions cannot be reused as bridge steps.');
fprintf('PASS: two conversions with reconstructed path.\n');

wrong = addBridge(model,1,14);
expectNoFinitePath(wrong,motif,cfg,p,'wrong irreversible direction');
reverse = wrong; reverse.lb(end) = -1000; reverse.ub(end) = 0;
[ranking,paths] = rankCircuitInputsByCAMP(reverse,motif,cfg,p);
assert(ranking.Distance(1) == 1 && paths(1).direction(1) == -1);
forwardOnly = cfg; forwardOnly.FORWARD_ONLY = true;
expectNoFinitePath(reverse,motif,forwardOnly,p,'reverse disallowed by forward-only setting');
closed = one; closed.lb(end) = 0; closed.ub(end) = 0;
expectNoFinitePath(closed,motif,cfg,p,'closed reaction bounds');
unannotated = one; unannotated.rxnECNumbers{end} = '';
expectNoFinitePath(unannotated,motif,cfg,p,'missing required EC evidence');
spontaneous = one; spontaneous.rxnNames{end} = 'Spontaneous conversion';
expectNoFinitePath(spontaneous,motif,cfg,p,'spontaneous reaction');
fprintf('PASS: reaction direction and bounds govern reachability.\n');

currency = addBridge(addBridge(model,14,16),16,1);
expectNoFinitePath(currency,motif,cfg,p,'ATP cannot provide an intermediate shortcut');
amp = two; amp.mets{15} = 'amp[c]'; amp.metNames{15} = 'AMP';
[ranking,~] = rankCircuitInputsByCAMP(amp,motif,cfg,p);
assert(ranking.Distance(1) == 2,'AMP must remain an allowed intermediate.');
adenosine = model;
adenosine.mets{15} = 'amp[c]'; adenosine.metNames{15} = 'AMP';
adenosine.mets{19} = 'ado[c]'; adenosine.metNames{19} = 'Adenosine';
adenosine = addBridge(addBridge(addBridge(adenosine,14,15),15,19),19,1);
[ranking,~] = rankCircuitInputsByCAMP(adenosine,motif,cfg,p);
assert(ranking.Distance(1) == 3,'Adenosine must remain an allowed intermediate.');
fprintf('PASS: ATP shortcuts excluded while AMP and adenosine are retained.\n');

isomer = one;
isomer.mets{14} = '23camp[c]'; isomer.metNames{14} = '2'',3''-cyclic AMP';
[ranking,paths] = rankCircuitInputsByCAMP(isomer,motif,cfg,p);
assert(isinf(ranking.Distance(1)) && strcmp(paths(1).status,'source_missing'));
isomer.mets{14} = 'camp[c]'; % A generic ID cannot override explicit isomer evidence.
[ranking,~,info] = rankCircuitInputsByCAMP(isomer,motif,cfg,p);
assert(isinf(ranking.Distance(1)) && isempty(info.sourceMetaboliteIndices));
fprintf('PASS: 2'',3''-cyclic AMP is not mistaken for the cAMP source.\n');

sideSignal = two; sideSignal.S(13,19) = 1; % m_out appears only as a side product.
expectNoFinitePath(sideSignal,motif,cfg,p,'selected neuron output participation');
sharedEC = one; sharedEC.rxnECNumbers{end} = sharedEC.rxnECNumbers{1};
expectNoFinitePath(sharedEC,motif,cfg,p,'bridge enzyme overlaps circuit EC');
sharedGene = one; sharedGene.grRules{end} = sharedGene.grRules{1};
expectNoFinitePath(sharedGene,motif,cfg,p,'bridge enzyme overlaps circuit gene');
fprintf('PASS: other signal participation and circuit enzyme conflicts excluded.\n');

internalConflict = two;
internalConflict.rxnECNumbers{20} = internalConflict.rxnECNumbers{19};
[ranking,paths] = rankCircuitInputsByCAMP(internalConflict,motif,cfg,p);
assert(ranking.Distance(1) == 2 && paths(1).enzymeConflict);
assert(ranking.BridgeEnzymeConflict(1));
fprintf('PASS: enzyme reuse within the bridge is reported.\n');

transport = one;
transport.mets{1} = 'input[m]'; transport.metComps(1) = 2;
expectNoFinitePath(transport,motif,cfg,p,'transport gate');
transportCfg = cfg; transportCfg.CAMP_ALLOW_TRANSPORT = true;
[ranking,paths] = rankCircuitInputsByCAMP(transport,motif,transportCfg,p);
assert(ranking.Distance(1) == 1 && isequal(paths(1).metIndex(:),[14;1]));
fprintf('PASS: transport requires explicit permission.\n');

% Copies of the same chemical must still support transport when requested.
sameChemical = addBridge(addBridge(model,14,18),18,1);
expectNoFinitePath(sameChemical,motif,cfg,p,'same-chemical compartment transport');
[ranking,paths] = rankCircuitInputsByCAMP(sameChemical,motif,transportCfg,p);
assert(ranking.Distance(1) == 2 && isequal(paths(1).metIndex(:),[14;18;1]));
assert(strcmp(ranking.SourceMetabolite{1},'camp[c]'), ...
    'The compartment filter must keep camp[m] from becoming a direct source.');
fprintf('PASS: same-chemical transports retain model species and source compartment.\n');

zeroDepth = cfg; zeroDepth.CAMP_MAX_DEPTH = 0;
[ranking,paths] = rankCircuitInputsByCAMP(exact,motif,zeroDepth,p);
assert(ranking.Distance(1) == 0 && strcmp(paths(1).status,'exact_input'));
expectNoFinitePath(one,motif,zeroDepth,p,'depth zero permits only an exact input');
shallow = cfg; shallow.CAMP_MAX_DEPTH = 1;
[ranking,paths] = rankCircuitInputsByCAMP(two,motif,shallow,p);
assert(isinf(ranking.Distance(1)) && strcmp(paths(1).status,'not_reached_within_depth'));
atLimit = cfg; atLimit.CAMP_MAX_DEPTH = 2;
[ranking,~] = rankCircuitInputsByCAMP(two,motif,atLimit,p);
assert(ranking.Distance(1) == 2,'A path exactly at the depth limit must be retained.');
deep = cfg; deep.CAMP_MAX_DEPTH = Inf;
[ranking,paths] = rankCircuitInputsByCAMP(model,motif,deep,p);
assert(isinf(ranking.Distance(1)) && strcmp(paths(1).status,'unreachable'));
fprintf('PASS: depth truncation is distinguished from graph disconnection.\n');

missing = model;
missing.mets{14} = 'missing_source[c]'; missing.metNames{14} = 'No cAMP here';
[ranking,paths] = rankCircuitInputsByCAMP(missing,motif,cfg,p);
assert(isinf(ranking.Distance(1)) && strcmp(paths(1).status,'source_missing'));
fprintf('PASS: missing source is reported without a fabricated path.\n');

[twin,second] = addSecondCircuit(two,p);
[ranking,paths] = rankCircuitInputsByCAMP(twin,[motif;second],cfg,p);
assert(isequal(ranking.MotifIndex(:),[2;1]));
assert(isequal(ranking.Distance(:),[1;2]));
assert(paths(1).metIndex(end) == second.metIndex(1));
assert(paths(2).metIndex(end) == motif.metIndex(1));
[ranking,~] = rankCircuitInputsByCAMP(two,[motif;motif],cfg,p);
assert(isequal(ranking.MotifIndex(:),[1;2]),'Equal distances must preserve discovery order.');
fprintf('PASS: ranking preserves original motif indices and stable ties.\n');
[ranking,paths,info] = rankCircuitInputsByCAMP(model,motif([]),cfg,p);
assert(height(ranking) == 0 && isempty(paths) && strcmp(info.status,'no_motifs'));
fprintf('PASS: no motifs produces an empty ranking.\n');
% Report the complete reactions, including co-substrates not on the path.
testFolder = tempname; mkdir(testFolder);
cleanup = onCleanup(@() rmdir(testFolder,'s')); %#ok<NASGU>
scaled = one; scaled.S(14,19) = -2.5; scaled.S(19,19) = -1;
[ranking,paths,rankInfo] = rankCircuitInputsByCAMP(scaled,motif,cfg,p);
reportFile = fullfile(testFolder,'camp_report.txt');
writeCAMPInputRankingReport(reportFile,scaled,motif,ranking,paths,rankInfo,p,cfg);
report = fileread(reportFile);
assert(contains(report,'r_CAMP_1') && contains(report,'2.5 camp[c](m_cAMP)'));
assert(contains(report,'other_aux[c]') && contains(report,'input[c](m_in)'));
assert(contains(report,'Original bounds:') && contains(report,'Eqn   :') && ...
    contains(report,'Names :') && contains(report,'IDs   :'));
info.metadata = rankInfo.metadata; info.status = 'exhausted'; info.exhaustive = true;
info.elapsedSeconds = 0; info.nodes = 0; info.parallelUsed = false;
circuitReport = fullfile(testFolder,'circuit_report.txt');
writeMetabolicCircuitReport(circuitReport,scaled,motif,info,p,cfg,ranking,paths,rankInfo);
assert(contains(fileread(circuitReport),'cAMP INPUT CONNECTION (OUTPUT m_in)'));
[ranking,paths,rankInfo] = rankCircuitInputsByCAMP(two,motif,cfg,p);
writeCAMPInputRankingReport(reportFile,two,motif,ranking,paths,rankInfo,p,cfg);
report = fileread(reportFile);
assert(contains(report,'r_CAMP_2') && contains(report,'bridge[c](m_CAMP_1)'));
fprintf('PASS: full connector equations, path labels and embedded circuit report.\n');
fprintf('All cAMP-ranking MATLAB self-tests passed.\n');
end

function [model,motif,cfg,p] = fixture()
p = circuitPattern('substrate');
S = sparse(19,18); S(1:13,:) = p.T;
model.mets = arrayfun(@(i) sprintf('signal_%02d[c]',i),(1:13)', ...
    'UniformOutput',false);
model.mets{1} = 'input[c]';
model.mets = [model.mets;{'camp[c]';'bridge[c]';'atp[c]'; ...
    'camp23[c]';'camp[m]';'other_aux[c]'}];
model.metNames = model.mets;
model.metNames{14} = '3'',5''-cyclic AMP';
model.metNames{16} = 'ATP'; model.metNames{17} = '2'',3''-cyclic AMP';
model.metNames{18} = '3'',5''-cyclic AMP';
for j = 1:18
    if ~any(S(:,j) > 0)
        S(end+1,j) = 1;
        model.mets{end+1,1} = sprintf('waste_%02d[c]',j);
        model.metNames{end+1,1} = sprintf('Private waste %02d',j);
    end
end
model.S = S;
model.metFormulas = repmat({''},numel(model.mets),1);
model.comps = {'c';'m'}; model.metComps = ones(numel(model.mets),1);
model.metComps(18) = 2;
model.rxns = arrayfun(@(j) sprintf('circuit_%02d',j),(1:18)', ...
    'UniformOutput',false);
model.rxnNames = model.rxns;
model.rxnECNumbers = arrayfun(@(j) sprintf('7.1.1.%d',j),(1:18)', ...
    'UniformOutput',false);
model.grRules = arrayfun(@(j) sprintf('CIRCUIT_G%d',j),(1:18)', ...
    'UniformOutput',false);
model.lb = zeros(18,1); model.ub = 1000*ones(18,1);
motif.metIndex = (1:13)'; motif.orientedIndex = (1:18)';
motif.reactionIndex = (1:18)'; motif.direction = ones(18,1);
cfg = circuitSearchDefaults(); cfg.FORWARD_ONLY = false;
cfg.ENZYME_EVIDENCE = 'EC'; cfg.VERBOSE = false;
cfg.CAMP_SOURCE_METABOLITES = {}; cfg.CAMP_SOURCE_COMPARTMENT = 'c';
cfg.CAMP_ALLOW_TRANSPORT = false; cfg.CAMP_EXCLUDED_PATH_METABOLITES = {};
cfg.CAMP_MAX_DEPTH = 6;
end

function model = addBridge(model,substrate,product)
j = size(model.S,2)+1;
model.S(:,j) = sparse(size(model.S,1),1);
model.S(substrate,j) = -1; model.S(product,j) = 1;
model.rxns{j,1} = sprintf('bridge_%02d',j);
model.rxnNames{j,1} = sprintf('Synthetic bridge %02d',j);
model.rxnECNumbers{j,1} = sprintf('1.1.1.%d',100+j);
model.grRules{j,1} = sprintf('BRIDGE_G%d',j);
model.lb(j,1) = 0; model.ub(j,1) = 1000;
end

function [model,motif] = addSecondCircuit(model,p)
nmet = numel(model.mets); nrxn = numel(model.rxns);
metIndex = [15;(nmet+1:nmet+12)'];
for i = 2:13
    row = metIndex(i);
    model.mets{row,1} = sprintf('second_signal_%02d[c]',i);
    model.metNames{row,1} = sprintf('Second signal %02d',i);
    model.metFormulas{row,1} = ''; model.metComps(row,1) = 1;
end
model.S(numel(model.mets),nrxn+18) = 0;
model.S(metIndex,nrxn+(1:18)) = p.T;
for j = 1:18
    col = nrxn+j;
    if ~any(model.S(:,col) > 0)
        row = size(model.S,1)+1; model.S(row,col) = 1;
        model.mets{row,1} = sprintf('second_waste_%02d[c]',j);
        model.metNames{row,1} = sprintf('Second private waste %02d',j);
        model.metFormulas{row,1} = ''; model.metComps(row,1) = 1;
    end
    model.rxns{col,1} = sprintf('second_circuit_%02d',j);
    model.rxnNames{col,1} = model.rxns{col};
    model.rxnECNumbers{col,1} = sprintf('6.1.1.%d',j);
    model.grRules{col,1} = sprintf('SECOND_G%d',j);
    model.lb(col,1) = 0; model.ub(col,1) = 1000;
end
motif.metIndex = metIndex; motif.orientedIndex = (nrxn+(1:18))';
motif.reactionIndex = motif.orientedIndex; motif.direction = ones(18,1);
end

function expectNoFinitePath(model,motif,cfg,p,label)
[ranking,paths] = rankCircuitInputsByCAMP(model,motif,cfg,p);
assert(isinf(ranking.Distance(1)) && isempty(paths(1).reactionIndex), ...
    'Unexpected path for %s.',label);
end
