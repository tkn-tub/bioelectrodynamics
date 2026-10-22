% A_Master_Consumed.m -----------------------------------------------------------
% Fixed-topology, crosstalk-free IFFL mapping.
%
%   STEP 1  BRENDA activators present in the model   -> seed table T  (cached)
%   STEP 2  map the WHOLE circuit at once            -> <tag>_stage2.mat
%   STEP 3  diagnostic + export, read back from disk -> <tag>_Diagnostics.txt
%
% NOTHING IS LOST IF THE SESSION DIES. Every run writes:
%   <tag>_run_log.txt     the full console transcript (diary)
%   iffl_stage1.mat       the BRENDA seed table, reused on the next run
%   <tag>_stage2.mat      the search result + a trimmed model copy
%   <tag>_Diagnostics.txt the funnel, histograms and verdict
%   <tag>_Chain_<n>.txt   the exported circuits
%   where <tag> = IFFL_consumed, so this study never overwrites the other
% To re-read any of it later without MATLAB redoing the work:
%   Report_IFFL_From_Mat('IFFL_consumed_stage2.mat', '', 3, 'IFFL_consumed')
%
% The search is not staged (motifs -> cascade -> counter -> neuron) and
% deliberately so: every crosstalk failure came from a later stage picking a
% reaction that collided with an earlier stage's signals. The whole circuit is
% accepted or rejected as a unit against one incidence table.
%
% Needs on path: iffl_topology.m, map_iffl_circuit.m, export_iffl_chain.m,
% Report_IFFL_From_Mat.m, parse_brenda_flat.m, iffl_coverage_check.m,
% addingPathParentFolderByName.m, COBRA; model + brenda_2026_1.txt under 'code'.
% ---------------------------------------------------------------------------

% Companion study to A_Master.m. Same model, same seeds, same filters - the one
% difference is that the signals BETWEEN blocks are consumed rather than read
% catalytically. Outputs are tagged so both studies live in the same folder.
clc; clearvars; close all;
diary off;
studyTag = 'IFFL_consumed';   % every output of this study carries this prefix
logFile  = [studyTag '_run_log.txt'];
if exist(logFile,'file'), delete(logFile); end
diary(logFile); diary on;
fprintf('IFFL run started %s\n', datestr(now));
cleanupDiary = onCleanup(@() diary('off'));    % closes the log even on an error

addingPathParentFolderByName('code');
if exist('readCbModel','file') ~= 2, initCobraToolbox(false); end

modelFile = which('Puniya2021_TNM1055.xml');
if isempty(modelFile), modelFile = which('Puniya2021_TNM1055.mat'); end
assert(~isempty(modelFile), 'Model not found under ''code''.');
[~,~,ext] = fileparts(modelFile);
if strcmpi(ext,'.mat'), L = load(modelFile); fn = fieldnames(L); model = L.(fn{1});
else, model = readCbModel(modelFile); end
fprintf('Model %s: %d reactions, %d metabolites\n', modelFile, numel(model.rxns), numel(model.mets));

% ---- options --------------------------------------------------------------
% The ONLY hard rule is signal exclusivity: each circuit signal is touched only
% by the circuit reactions its role declares (iffl_topology.m). distinctEC and
% noRegulatoryCrosstalk are REPORTED per circuit instead of rejecting it.
opt = struct();
opt.collapseCompartments  = true;    % merge all compartments into one pot FIRST.
                                     % Do this rather than opt.compartment='any':
                                     % merging makes akg[c] and akg[m] one species,
                                     % ignoring compartments makes them two.
opt.compartment           = 'c';     % after collapsing, everything is [c]
opt.consumingN            = false;   % m_c gates r_Nw without being consumed
opt.gateMode              = {'activate','inhibit'};
                                     % per stage. {'activate','activate'} is the
                                     % original. 'inhibit' makes that stage's gate
                                     % a BRENDA INHIBITOR of the reaction producing
                                     % m_d, which drops r_I for that stage and draws
                                     % its gate from a different species pool.
opt.interstageGate        = 'any';   % m_d1 may ACTIVATE or INHIBIT r_A_2. Only
                                     % 1.3%% of block pairs had an activating link,
                                     % so this is the pairing bottleneck. An
                                     % inhibiting link flips the cascade sign; the
                                     % net polarity is computed and reported per
                                     % circuit rather than assumed.
opt.neuronGate            = 'any';   % m_c may ACTIVATE r_Nw (positive weight w)
                                     % or INHIBIT it (negative w). A single neuron
                                     % with a negative weight is ordinary, so there
                                     % is no reason to search only half the pool.
opt.catalyticInterstage   = true;    % m_d1 GATES r_A_2 -> m_in_2 instead of being
                                     % consumed by D2. Breaks the mass-flow hinge that
                                     % pins both stages into one subsystem, and is the
                                     % high-impedance interstage by construction.
opt.fastArmSteps          = 1;       % 1 = fast arm is r_Fd only (no m_Ff)
opt.requireMI             = false;   % false drops m_I: r_I1/r_I2 are first order in
                                     % m_c/m_d2, rate programmed by enzyme loading.
opt.requireSharedMI       = false;   % only read when requireMI is true: forces both
                                     % removals onto ONE finite m_I pool.
opt.carrierTest           = true;    % reject arms that regenerate their own rail
opt.distinctEC            = false;   % report, do not reject, enzyme reuse
opt.noRegulatoryCrosstalk = false;   % report, do not reject, non-selective gates
opt.requireEC             = false;
opt.maxArmPaths           = 2000;    % per rail; a WARNING prints if it binds
opt.sinkMode              = 'prefer';% r_Fr / r_lk go to the empty set in the diagram:
                                     % map a clean model reaction when one exists,
                                     % otherwise declare an ADDED sink (CSTR outflow,
                                     % protease, phosphatase) with a free rate constant.
                                     % 'map' = strict, 'sink' = always added.
opt.maxHousekeeping       = 80;      % candidates kept per r_Fr / r_lk slot
opt.topBlocks             = 400;
opt.maxCircuits           = 200;
opt.verbose               = true;

% ---- optional: collapse to a single pot -----------------------------------
if isfield(opt,'collapseCompartments') && opt.collapseCompartments
    assert(exist('collapse_compartments','file')==2, ...
        'collapse_compartments.m not found on the path.');
    [model, collapseRep] = collapse_compartments(model);  %#ok<ASGLU>
end


% ---- THIS STUDY: the inter-block signals are CONSUMED, not read ------------
% m_d1 is drawn on by D2's arms, m_d2 by r_c, m_c by r_Nw. The differentiator
% GATES (m_Fd_1, m_Fd_2) stay catalytic - that is the IFFL primitive itself,
% and each is still cleared by its own r_Fr.
% Set last, so these win over the block above whatever it says.
opt.catalyticInterstage = false;      % no r_A_2: D2 consumes m_d1 directly
opt.consumingN          = true;       % r_Nw consumes m_c instead of being gated
opt.interstageGate      = 'activate'; % unused when catalyticInterstage is false
opt.neuronGate          = 'activate'; % unused when consumingN is true

% =============================== STEP 1 ====================================
% The BRENDA parse is slow and does not depend on opt, so it is cached.
% Delete iffl_stage1.mat to force a re-parse.
% The seed table is matched against THIS model's reaction and metabolite
% names, so it is only valid for the model that produced it. Collapsing
% compartments changes both. The signature below makes a stale cache
% impossible to use by accident.
modelSig = sprintf('%d_%d_%d', numel(model.rxns), numel(model.mets), ...
                   double(opt.collapseCompartments));
useCache = false;
if exist('iffl_stage1.mat','file')
    S1 = load('iffl_stage1.mat');
    if isfield(S1,'modelSig') && strcmp(S1.modelSig, modelSig)
        useCache = true;
    else
        if isfield(S1,'modelSig'), was = S1.modelSig; else, was = '(unstamped)'; end
        fprintf(['\nSTEP 1: cached seeds are STALE (built for model %s, this model is %s).\n' ...
                 '        Recomputing against the current model.\n'], was, modelSig);
    end
end
needIN = any(strcmp(opt.gateMode, 'inhibit')) || ~strcmp(opt.neuronGate, 'activate') ...
         || ~strcmp(opt.interstageGate, 'activate');
if useCache && (~needIN || isfield(S1,'Tin'))
    T = S1.T; if needIN, Tin = S1.Tin; end
    fprintf('\nSTEP 1: %d activator seeds (cache valid for this model).\n', height(T));
    if needIN, fprintf('        %d inhibitor seeds.\n', height(Tin)); end
else
    brendaFile = which('brenda_2026_1.txt');
    assert(~isempty(brendaFile), 'brenda_2026_1.txt not found under ''code''.');
    T = iffl_coverage_check(model, brendaFile, struct('wantAC', true));
    stamp = datestr(now); %#ok<NASGU>
    if needIN
        Tin = iffl_coverage_check(model, brendaFile, struct('wantAC', false));
        save('iffl_stage1.mat', 'T', 'Tin', 'stamp', 'modelSig', '-v7.3');
        fprintf('\nSTEP 1: %d activator seeds, %d inhibitor seeds (saved).\n', ...
            height(T), height(Tin));
    else
        save('iffl_stage1.mat', 'T', 'stamp', 'modelSig', '-v7.3');
        fprintf('\nSTEP 1: %d activator seeds (saved to iffl_stage1.mat).\n', height(T));
    end
end
assert(height(T) > 0, 'No seeds - nothing to anchor on.');
if needIN
    assert(height(Tin) > 0, 'No inhibitor seeds - gateMode asks for an inhibition stage.');
    seedArg = struct('AC', T, 'IN', Tin);
else
    seedArg = T;
end

% =============================== STEP 2 ====================================
out = map_iffl_circuit(model, seedArg, opt);
fprintf('STEP 2: %d circuits.\n', numel(out.circuits));

mdl   = trim_model(model);           % only the fields the report/export need
stamp = datestr(now);                %#ok<NASGU>
save([studyTag '_stage2.mat'], 'out', 'mdl', 'opt', 'stamp', '-v7.3');
fprintf('Saved %s_stage2.mat\n', studyTag);

% =============================== STEP 3 ====================================
% Everything below reads the checkpoint, so you can re-run this line alone in a
% fresh session: Report_IFFL_From_Mat([studyTag '_stage2.mat'], '', 3, studyTag)
Report_IFFL_From_Mat([studyTag '_stage2.mat'], [studyTag '_Diagnostics.txt'], 3, studyTag);
if ~isempty(out.circuits), write_circuit_csv(mdl, out, [studyTag '_circuits.csv']); end

fprintf('\nDone %s. Files: %s_run_log.txt, %s_Diagnostics.txt, %s_stage2.mat\n', ...
    datestr(now), studyTag, studyTag, studyTag);
diary off;

% ---------------------------------------------------------------------------
function m = trim_model(model)
% keep the checkpoint small but complete enough for the report and the export
    keep = {'mets','metNames','metFormulas','rxns','rxnNames','rxnECNumbers', ...
            'lb','ub','S','rxnGeneMat'};
    m = struct();
    for q = 1:numel(keep)
        if isfield(model, keep{q}), m.(keep{q}) = model.(keep{q}); end
    end
end

function write_circuit_csv(mdl, out, fname)
    topo = out.topology;
    fid = fopen(fname, 'w'); c0 = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'rank,nGenes,addedSinks,regCollisions,sharedEC');
    for q = 1:numel(topo.rxnRoles), fprintf(fid, ',%s', topo.rxnRoles{q}); end
    for q = 1:numel(topo.metRoles), fprintf(fid, ',%s', topo.metRoles{q}); end
    fprintf(fid, '\n');
    for n = 1:numel(out.circuits)
        c = out.circuits(n);
        fprintf(fid, '%d,%d,%d,%d,%d', n, c.nGenes, c.nSinks, c.notes.nRegCollisions, c.notes.nSharedEC);
        for q = 1:numel(c.rxnIdx)
            if isnan(c.rxnIdx(q)), fprintf(fid, ',SINK');
            else, fprintf(fid, ',%s', mdl.rxns{c.rxnIdx(q)}); end
        end
        for q = 1:numel(c.sig),    fprintf(fid, ',%s', mdl.mets{c.sig(q)}); end
        fprintf(fid, '\n');
    end
end
