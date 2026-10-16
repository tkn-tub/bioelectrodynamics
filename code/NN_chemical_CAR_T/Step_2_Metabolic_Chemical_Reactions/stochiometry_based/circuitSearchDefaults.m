function cfg = circuitSearchDefaults()
%CIRCUITSEARCHDEFAULTS Options shared by the master and direct function calls.
cfg.FORWARD_ONLY = true;
cfg.USE_PARALLEL = false;
cfg.NUM_WORKERS = [];
cfg.PARALLEL_PROFILE = '';
cfg.MAX_MOTIFS = 20;
cfg.ONE_MOTIF_PER_INPUT = true; % Stop an input's search after its first circuit.
cfg.MAX_SECONDS = 3600;
cfg.VERBOSE = true;
cfg.ENZYME_EVIDENCE = 'EC';
cfg.ENZYME_ID_FIELD = '';
cfg.METABOLITE_SCOPE = 'roles';
cfg.METABOLITE_IDENTITY = 'chemical';
cfg.AUXILIARY_REUSE_ALLOWED = {};
cfg.EXCLUDED_ROLE_METABOLITES = {};
cfg.COMPARTMENT = '';
cfg.NEURON_INPUT_MODE = 'substrate';
cfg.ACTIVATION_PAIRS = cell(0,2);
cfg.FIXED_METABOLITES = struct();
cfg.FIXED_REACTIONS = struct();
cfg.MODEL_FILE = '';
cfg.TOLERANCE = 1e-9;
% Rank returned circuits by directed conversion distance from cAMP to m_in.
cfg.CAMP_SOURCE_METABOLITES = {}; % Empty: identify 3',5'-cAMP by ID/name.
cfg.CAMP_SOURCE_COMPARTMENT = 'c'; % Cytosol; '' permits any compartment.
cfg.CAMP_ALLOW_TRANSPORT = false;
cfg.CAMP_EXCLUDED_PATH_METABOLITES = {}; % Extra IDs/base IDs; hubs excluded too.
cfg.CAMP_MAX_DEPTH = 6; % 0: only exact m_in=cAMP; Inf: unbounded graph distance.
end

