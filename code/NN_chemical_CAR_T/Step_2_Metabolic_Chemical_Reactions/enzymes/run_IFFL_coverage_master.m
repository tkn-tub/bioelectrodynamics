% run_IFFL_coverage_master.m
% Master file: locate the Puniya2021 model and the BRENDA flat file (both kept
% under the 'code' folder), then run the IFFL coverage check.
%
% Needs on the path (added below): iffl_coverage_check.m, parse_brenda_flat.m,
% addingPathParentFolderByName.m, the COBRA Toolbox.

clc; clearvars; close all;

addingPathParentFolderByName('code');          % puts 'code' and its subfolders on the path

% ---- make sure the COBRA Toolbox is ready (skip if already initialised) ----
if exist('readCbModel','file') ~= 2
    initCobraToolbox(false);                   % false = no auto-update
end

% ---- locate the metabolic model under 'code' via which() ----
modelFile = which('Puniya2021_TNM1055.xml');   % change the name if yours differs
if isempty(modelFile)
    modelFile = which('Puniya2021_TNM1055.mat');   % fallback: a saved COBRA model
end
assert(~isempty(modelFile), ...
    'Model file not found on the path. Put Puniya2021_TNM1055.(xml|mat) under ''code''.');

[~,~,ext] = fileparts(modelFile);
if strcmpi(ext,'.mat')
    S = load(modelFile);                       % expect one struct field holding the model
    fn = fieldnames(S);  model = S.(fn{1});
else
    model = readCbModel(modelFile);
end
fprintf('Loaded model %s  (%d reactions, %d metabolites)\n', ...
        modelFile, numel(model.rxns), numel(model.mets));

% ---- locate the BRENDA flat file under 'code' via which() ----
brendaFile = which('brenda_2026_1.txt');       % same pattern as BRENDA_credentials.txt
assert(~isempty(brendaFile), ...
    'brenda_2026_1.txt not found on the path. Put it under ''code''.');
fprintf('Using BRENDA flat file: %s\n', brendaFile);

% ---- (optional) EC override from your kinetics workbook, if the SBML lacks ECs ----
opts = struct('wantAC', true);                 % activators (set false for inhibitors)
% opts.ecMap = <containers.Map rxnID -> 'x.x.x.x'>;   % fill if model.rxnECNumbers is empty

% ---- run the coverage check ----
T = iffl_coverage_check(model, brendaFile, opts);

fprintf('\nDone. %d seed reactions; details in iffl_coverage_seeds.csv\n', height(T));
