function T = iffl_coverage_check(model, flatFile, opts)
% IFFL_COVERAGE_CHECK  Is the IFFL (activator-gated degradation) mappable in a model?
%
%   T = iffl_coverage_check(model, 'brenda_flat.txt')
%   T = iffl_coverage_check()            % runs a built-in synthetic SELF-TEST
%
% It answers the cheap gate before any delay-chain search:
%   "For each reaction that could serve as r_I (degrades a candidate output m_d1),
%    does BRENDA list, for that reaction's EC, an ACTIVATOR that is also a
%    metabolite present in the model (a candidate m_Fd)?"
%
% INPUTS
%   model     COBRA model (readCbModel of Puniya2021_TNM1055.xml). Needs:
%             .S .rxns .mets .metNames ; EC via .rxnECNumbers (preferred) or set
%             opts.ecFromWorkbook to supply EC per rxn from your kinetics table.
%   flatFile  path to the BRENDA flat text file.
%
% OPTIONS (opts struct, all optional)
%   .wantAC   true  (use activators; set false to use inhibitors instead)
%   .ecMap    containers.Map rxnID->EC to override/fill EC numbers
%   .xref     containers.Map modelMetName(normalised)->ChEBI/KEGG id, and a
%             parallel ligand xref, to match by ID instead of name (recommended).
%
% OUTPUT  table T of seeds: {rxn, EC, role, matched_regulators, nConsumed}
%         plus a printed coverage report. Zero rows with good coverage = this
%         path is a dead end for Puniya2021 (fall back to annihilation).

    if nargin==0, T = selftest(); return; end
    if nargin<3, opts = struct(); end
    if ~isfield(opts,'wantAC'), opts.wantAC = true; end
    field = 'AC'; if ~opts.wantAC, field = 'IN'; end

    nR = numel(model.rxns);
    % ---- EC per reaction ----
    EC = repmat({''},nR,1);
    if isfield(model,'rxnECNumbers'), EC = cellfun(@char,model.rxnECNumbers,'uni',0); end
    if isfield(opts,'ecMap')
        for j=1:nR, if isKey(opts.ecMap,model.rxns{j}), EC{j}=opts.ecMap(model.rxns{j}); end, end
    end
    % normalise model metabolite names for matching
    metN = cellfun(@normLigand, model.metNames, 'uni', 0);
    modelMetSet = containers.Map(metN, num2cell(1:numel(metN)));

    reg = parse_brenda_flat(flatFile, {field});

    % ---- coverage diagnostics FIRST (so a zero is interpretable) ----
    hasEC   = find(~cellfun(@isempty,EC));
    ecInBr  = hasEC(cellfun(@(e) isKey(reg,e) && ~isempty(reg(e).(field)), EC(hasEC)));
    fprintf('\n--- IFFL coverage (field %s) ---\n', field);
    fprintf('reactions total ............ %d\n', nR);
    fprintf('reactions with an EC ....... %d\n', numel(hasEC));
    fprintf('  ...whose EC has %s in BRENDA %d\n', field, numel(ecInBr));

    % ---- seeds: regulator of r_I's EC is also a model metabolite ----
    rows = {};
    for jj = 1:numel(ecInBr)
        j = ecInBr(jj); e = EC{j};
        ligs = reg(e).(field);
        hit = ligs(cellfun(@(L) isKey(modelMetSet,L), ligs));
        if ~isempty(hit)
            consumed = sum(model.S(:,j) < 0);      % candidate m_d1 count
            rows(end+1,:) = {model.rxns{j}, e, field, strjoin(hit,'; '), consumed}; %#ok<AGROW>
        end
    end
    if isempty(rows)
        T = cell2table(cell(0,5),'VariableNames',{'rxn','EC','role','regulators','nConsumed'});
    else
        T = cell2table(rows,'VariableNames',{'rxn','EC','role','regulators','nConsumed'});
        T = sortrows(T,'nConsumed','descend');
    end
    fprintf('SEEDS (candidate r_I with a %s that is a model metabolite = m_Fd): %d\n', field, height(T));
    if height(T)>0, disp(head(T, min(20,height(T)))); end
    fprintf(['NOTE: name-matched (approximate). For the definitive run pass opts.xref\n' ...
             '      to match BRENDA ligands<->model metabolites by ChEBI/KEGG id\n' ...
             '      (BiGG ids -> ChEBI via MetaNetX chem_xref or the BiGG metabolites file).\n']);
    try, writetable(T,'iffl_coverage_seeds.csv'); catch, end
end

function s = normLigand(s)
    s = lower(strtrim(char(s)));
    s = regexprep(s,'^((l|d|dl|\(s\)|\(r\)|alpha|beta)-)+','');
    s = regexprep(s,'[^a-z0-9]+',' ');
    s = strtrim(s);
end

function T = selftest()
% tiny synthetic model + fake BRENDA to prove the counter before real data
    model.rxns = {'r_I_cand';'r_pyk';'r_x'};
    model.mets = {'nadh';'pyr';'adp';'pep';'fbp';'xmp'};
    model.metNames = {'NADH';'pyruvate';'ADP';'PEP';'fructose 1,6-bisphosphate';'XMP'};
    model.rxnECNumbers = {'1.1.1.27';'2.7.1.40';'9.9.9.9'};
    model.S = zeros(6,3); model.S([1 2],1)=-1; model.S([3 4],2)=-1; % consumed sets
    tmp = [tempname '.txt']; fid=fopen(tmp,'w');
    fprintf(fid,['ID\t1.1.1.27\nAC\t#1# fructose 1,6-bisphosphate {act} <1>\n///\n' ...
                 'ID\t2.7.1.40\nAC\t#1# fructose 1,6-bisphosphate {act} <1>\n///\n' ...
                 'ID\t9.9.9.9\n///\n']);
    fclose(fid);
    T = iffl_coverage_check(model, tmp, struct('wantAC',true));
    assert(height(T)==2 && all(ismember({'r_I_cand','r_pyk'}, T.rxn)), 'self-test FAILED');
    fprintf('\nself-test PASSED (found the 2 planted seeds, skipped the unregulated one)\n');
    delete(tmp);
end
