function [Values, Km, Kcat, Reactions] = brendaKineticsTable(txtFile, xlsxFile, varargin)
% BRENDA kinetics workbook for FULL CHEMICAL CHAIN and legacy motif TXT files.
% Updated 2026-10-16. Keeps every printed reaction occurrence, including transport.
%
% brendaKineticsTable('Selected_Full_Chain_1.txt','kinetics_Full_Chain_1.xlsx');
% Offline layout check (no credentials or live measurements):
% brendaKineticsTable('Selected_Full_Chain_1.txt','kinetics_Full_Chain_1.xlsx', ...
%     'QueryBRENDA',false,'UseReferenceValues',false);
% Unnamed steps receive r_B<block>_<step>; repeated r_CAMP labels receive
% r_CAMP_1, r_CAMP_2, ... . Reaction mapping retains original headings and IDs.
% Bounds do not establish orientation relative to model.S: when Direction is
% absent, +1 means the printed left-to-right equation only (see Notes).
% G/tau are generated only for explicit legacy producer/consumer labels;
% abbreviated branch headings do not specify these pairs.
% Local EC recovery/annotation review reads ID/RN/SN/SY/RE fields from the
% BRENDA flat file. Unique equation matches take priority over exact enzyme
% names/synonyms. Original ECs, candidates and decisions remain in Reaction
% mapping and Notes. Generic model reaction IDs are not BRENDA EC identifiers.
% Local lookup recovers ECs only; the existing API queries supply kinetics.
% Input TXT and CredentialsFile are resolved with which, so files inside
% folders added with addpath(genpath(...)) need not be in the current folder.
% Explicit existing paths also work. Output XLSX location is caller-supplied
% (a bare output filename is written in the current working folder).
% The embedded reference snapshot belongs to the older network; it applies
% only on an exact ID/equation/direction match. Offline mode therefore creates
% an unpopulated template for this new chain. Live queries need credentials.txt.
% Self-contained: no separate Brenda_Kinetics.m is required.
%
% brendaKineticsTable('Block_7_Selected_Motifs.txt','kinetics_BRENDA.xlsx',email,password);
% brendaKineticsTable('Block_7_Selected_Motifs(7).txt','kinetics_BRENDA_updated.xlsx');
% Optional name/value arguments follow the filenames or the email/password pair:
%   'Organism'         'Homo sapiens' (use '' explicitly for all organisms)
%   'CredentialsFile'  'credentials.txt' (two lines: email, then password)
%   'QueryBRENDA'      true (false skips live queries; references can still apply)
%   'UseLocalEC'       true: search the local BRENDA flat file before API calls
%   'LocalBrendaFile'  'brenda_2026_1.txt' (resolved with which)
%   'ECOverrides'      N-by-2 strings: model ID or row label, complete EC(s)
%                     e.g. ["PYK","2.7.1.40"]. Explicit overrides take priority.
% Local equation matching ignores water/protons and recognises documented
% aliases/nucleotide abbreviations, but retains all other stoichiometry.
% A unique name-only match is labelled UNCONFIRMED BY EQUATION. Multiple
% candidates are not assigned automatically. Missing local file: retain TXT
% annotations and warn; inspect Reaction mapping before using the kinetics.
%   'PrimarySubstrate' struct, e.g. struct('r_I','L-glutamate')
%   'SubstrateAliases' N-by-2 string/cell array: model name, BRENDA name
%   'PreferredEC'      struct, e.g. struct('r_adapter_1','3.1.4.17')
%   'UseReferenceValues' true: apply the documented 2026-10-12 selections for
%                        these 23 reaction IDs/equations. False restores
%                        strict live-record matching without surrogates.
%   'ResolveTransferredEC' true: when an EC code returns NO Km, kcat, SA and
%                        MW record at all, read https://enzyme.expasy.org/EC/
%                        <code>.txt and follow a "Transferred entry" to the
%                        live code, then re-query BRENDA under it. An EC that
%                        has been split or renumbered carries no BRENDA data
%                        under its old number: EC 2.3.2.4 returns nothing,
%                        while its successor EC 4.3.2.9 returns 110 records.
%                        Needs outbound HTTPS to ExPASy; set false to skip.
%                        The old -> new mapping is written to the Query status
%                        and EC resolution sheets and into the Notes column.
%   'ECReplacements'   N-by-2 string array, e.g. ["2.3.2.4","4.3.2.9"]: state
%                        the live code yourself when ExPASy is unreachable
%                        (firewalled HPC). Checked before the ExPASy lookup and
%                        applied under the same "all four parameters empty"
%                        condition.
%   'SurrogateSubstrate' struct, e.g. struct('r_I','5-L-Glutamyl-L-alanine'):
%                        the ligand to score the K_M substitution against when
%                        the reaction's own substrate has no record.
%
% ESCALATION ORDER for each row. Every step is recorded in the Notes column and
% summarised in the machine-readable "Evidence tier" column, so no number is
% ever reported as a clean match when it is not:
%   1. exact substrate name, requested organism
%   2. exact substrate name, ANY organism            (OrganismFallback)
%   3. relaxed name / nucleotide abbreviation, requested organism, then ANY
%                                                    (RelaxedSubstrateMatch)
%   4. k_cat from any ligand of the same EC          (KcatAnySubstrate)
%   5. k_cat = SA*MW/60                              (KcatSource = SAxMW)
%   6. K_M of a surrogate ligand of the same EC      (SubstituteKm)
%   7. nothing at all for the code: resolve a transferred/deleted EC at ExPASy
%      and restart at step 1 under the live code     (ResolveTransferredEC)
% Step 7 fires only when Km, kcat, SA and MW are ALL empty, which is the
% signature of an obsolete code; a code that returns some data is current by
% definition and is never re-looked-up.
% With QueryBRENDA=false and UseReferenceValues=true, the embedded snapshot
% builds reference-matched rows without Python, credentials or API access.
% Reference choices may differ in organism, substrate, assay or direction;
% two kcat estimates use paired SA/MW records. All choices are documented.
% Known cAMP/cGMP names are matched as synonyms. A single reactant remaining
% after excluding water/protons is used when no primary substrate is tagged.
% For r_adapter_1 consuming cAMP, EC 3.1.4.17 is preferred when annotated.
% PreferredEC overrides this default; all annotated EC records are retained.
%
% Explicit email/password arguments take precedence. If both are omitted or
% empty, credentials are read from CredentialsFile. No default credentials,
% environment lookup, recovery, printing, or credential files are generated.
% Each unique complete EC is queried once for Km and kcat. All measured
% records, commentary, literature and organism labels are retained.
% Numeric candidates are medians of positive scalar values, NOT validated
% reaction parameters; maximum/range values are not pooled or averaged.
% Other substrate names use exact matching after motif tags/coefficients
% are removed. Supply SubstrateAliases for additional known synonyms.
% With UseReferenceValues=false, unresolved multiple EC matches and reverse
% directions leave values blank.
% Live queries do not infer activation kinetics. Reference mode supplements
% the summary with an explicit sourced snapshot, separate from live records.
% The output workbook is rebuilt and replaced on each
% run: old values and manual edits at xlsxFile are not retained.
% Requires MATLAB with Python/zeep for live queries; offline mode needs neither.
% SOAP reference: https://www.brenda-enzymes.org/soap.php

if nargin < 1 || isempty(txtFile), txtFile = 'Selected_Full_Chain_1.txt'; end
if nargin < 2 || isempty(xlsxFile), xlsxFile = 'kinetics_Full_Chain_1.xlsx'; end
% Accept the requested positional credentials AND the existing option syntax.
email = ''; password = ''; positionalCredentials = false;
optionNames = {'Organism','CredentialsFile','QueryBRENDA', ...
    'PrimarySubstrate','SubstrateAliases','PreferredEC','UseReferenceValues', ...
    'OrganismFallback','RelaxedSubstrateMatch','KcatAnySubstrate','SubstituteKm', ...
    'AllowReverse','ResolveTransferredEC','SurrogateSubstrate','ECReplacements','KcatSource', ...
    'UseLocalEC','LocalBrendaFile','ECOverrides'};
options = varargin;
isTextScalar = @(v)(ischar(v)&&(isrow(v)||isempty(v))) || ...
    (isstring(v)&&isscalar(v)&&~ismissing(v));
if ~isempty(options) && isTextScalar(options{1}) && ...
        ~any(strcmpi(char(options{1}),optionNames))
    assert(numel(options)>=2 && isTextScalar(options{2}), ...
        'Supply email and password together, followed by optional name/value arguments.');
    email = char(options{1}); password = char(options{2});
    assert(isempty(email)==isempty(password), ...
        'Supply both email and password, or leave both empty to read credentials.txt.');
    positionalCredentials = ~isempty(email);
    options = options(3:end);
end
clear varargin % Avoid keeping an extra copy of positional credentials.
p = inputParser;
p.PartialMatching = false;
addParameter(p,'Organism','Homo sapiens',@(v)ischar(v)||isstring(v));
addParameter(p,'CredentialsFile','credentials.txt',@(v)ischar(v)||isstring(v));
addParameter(p,'QueryBRENDA',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'UseLocalEC',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'LocalBrendaFile','brenda_2026_1.txt',@(v)ischar(v)||isstring(v));
addParameter(p,'ECOverrides',strings(0,2),@(v)isstring(v)||iscell(v));
addParameter(p,'PrimarySubstrate',struct(),@isstruct);
addParameter(p,'PreferredEC',struct(),@isstruct);
addParameter(p,'UseReferenceValues',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'SubstrateAliases',strings(0,2),@(v)isstring(v)||iscell(v));
addParameter(p,'OrganismFallback',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'RelaxedSubstrateMatch',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'KcatAnySubstrate',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'SubstituteKm',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'AllowReverse',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'ResolveTransferredEC',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'SurrogateSubstrate',struct(),@isstruct);
addParameter(p,'ECReplacements',strings(0,2),@(v)isstring(v)||iscell(v));
addParameter(p,'KcatSource','SAxMW', ...
    @(v)any(strcmpi(char(v),{'SAxMW','BRENDA','BRENDA else SAxMW'})));
parse(p,options{:}); opt = p.Results;
clear options
aliases = string(opt.SubstrateAliases);
assert(size(aliases,2)==2,'SubstrateAliases must have two columns.');
ecReplacements = string(opt.ECReplacements);
assert(size(ecReplacements,2)==2, ...
    'ECReplacements must have two columns: old EC, replacement EC.');
% Caller aliases take precedence over the shipped defaults.
defaults = defaultAliases();
keep = ~ismember(normalize(defaults(:,1)),normalize(aliases(:,1)));
aliases = [aliases; defaults(keep,:)];
txtFile = resolveInputFile(txtFile,'Reaction TXT');
assert(~strcmp(char(txtFile),char(xlsxFile)),'Input and output must differ.');
Reactions = parseMotifs(fileread(txtFile));
assert(~isempty(Reactions),'No reaction records found in the TXT file.');
fprintf('Parsed %d reaction occurrences across %d motif(s).\n', ...
    height(Reactions),numel(unique(Reactions.Motif)));
Reactions = recoverLocalEC(Reactions,opt,aliases);

% Validate credentials before touching the output workbook.
if opt.QueryBRENDA && ~positionalCredentials
    [email,password] = readCredentials(opt.CredentialsFile);
end
Km = emptyRecords('Km','mM'); Kcat = emptyRecords('kcat','1/s');
SArec = emptyRecords('SA','umol/min/mg'); MWrec = emptyRecords('MW','Da');
ecs = strings(0,1);
for i=1:height(Reactions)
    ecs = [ecs; string(regexp(char(Reactions.EC(i)), ...
        '(?<![\d.])\d+\.\d+\.\d+\.\d+(?![\d.])','match'))']; %#ok<AGROW>
end
ecs = unique(ecs,'stable');
status = cell(numel(ecs)+1,8);
status(1,:) = {'EC before transfer','EC actually queried','EC resolution', ...
    'Organism requested','Km query','kcat query','SA query','MW query'};
% Old EC in the motif file -> the code BRENDA was finally queried under.
% Identity for every code that already returned data.
ecMap = containers.Map('KeyType','char','ValueType','char');
ecResolution = {'EC before transfer','EC actually queried','Verdict', ...
    'Trail','BRENDA records before EC transfer', ...
    'BRENDA records under the live code','Accepted name of the live entry'};
for j=1:numel(ecs)
    ec = char(ecs(j)); queried = ec; verdict = ''; trail = ''; liveName = '';
    before = 0; after = 0;
    fprintf('EC %d/%d: %s\n',j,numel(ecs),ec);
    if opt.QueryBRENDA
        % Query each parameter independently: one failure must not discard
        % the others. Many EC codes have no record at all for the requested
        % organism (EC 6.3.2.12 has no human entry: humans lack the enzyme),
        % so an empty result is retried without the organism restriction.
        [a,sa] = fetchBothPools(ec,opt,email,password,'Km');
        [b,sb] = fetchBothPools(ec,opt,email,password,'kcat');
        [c,sc] = fetchBothPools(ec,opt,email,password,'SA');
        [d,sd] = fetchBothPools(ec,opt,email,password,'MW');
        before = height(a)+height(b)+height(c)+height(d);
        % LAST RESORT, and only when the code is completely silent: an EC
        % that has been transferred or deleted carries no BRENDA data under
        % its old number, and no amount of organism/substrate relaxation can
        % recover it. ExPASy is the authority on the renumbering.
        manual = '';
        if ~isempty(ecReplacements)
            hit = find(strtrim(ecReplacements(:,1))==string(ec),1);
            if ~isempty(hit), manual = char(strtrim(ecReplacements(hit,2))); end
        end
        if before==0 && (opt.ResolveTransferredEC || ~isempty(manual))
            if ~isempty(manual)
                % A locked-down machine may not reach ExPASy at all; the
                % caller can then state the replacement directly.
                newEC = manual;
                verdict = sprintf('caller-supplied ECReplacements: %s -> %s',ec,manual);
                trail = [ec ' -> ' manual];
            else
                [newEC,verdict,trail,liveName] = resolveTransferredEC(ec);
            end
            if ~isempty(newEC) && ~strcmp(newEC,ec)
                fprintf('  EC %s is not current: %s. Re-querying BRENDA as EC %s.\n', ...
                    ec,verdict,newEC);
                [a,sa] = fetchBothPools(newEC,opt,email,password,'Km');
                [b,sb] = fetchBothPools(newEC,opt,email,password,'kcat');
                [c,sc] = fetchBothPools(newEC,opt,email,password,'SA');
                [d,sd] = fetchBothPools(newEC,opt,email,password,'MW');
                after = height(a)+height(b)+height(c)+height(d);
                queried = newEC;
                fprintf('  EC %s: %d record(s) recovered after the transfer.\n',newEC,after);
            else
                fprintf('  EC %s returned nothing. Verdict: %s\n',ec,verdict);
            end
        elseif before==0
            verdict = 'not checked (ResolveTransferredEC = false)';
        else
            verdict = 'current entry (BRENDA returned records)';
        end
        Km=[Km;a]; Kcat=[Kcat;b]; SArec=[SArec;c]; MWrec=[MWrec;d]; %#ok<AGROW>
    else
        sa = 'Not queried (offline mode; see Reference choices)';
        sb = sa; sc = sa; sd = sa;
        verdict = 'not checked (offline mode)';
    end
    ecMap(ec) = queried;
    if strcmp(queried,ec), resolutionText = verdict;
    else, resolutionText = ['TRANSFERRED: ' ec ' -> ' queried '. ' verdict]; end
    status(j+1,:) = {ec,queried,resolutionText,char(opt.Organism),sa,sb,sc,sd};
    ecResolution(end+1,:) = {ec,queried,verdict,trail,before,after,liveName}; %#ok<AGROW>
end
clear email password

headers = {'Reaction pair','Reaction','Enzyme','EC code','K_M [mM]', ...
    'SA [μmol/min/mg]','MW [mg/μmol]','k_cat [s^-1]','[E] [mM]', ...
    'k = k_cat [E] / K_M [s^-1]', ...
    'Removed by (* = K_M there is for another substrate)','G','τ [s]','Notes', ...
    'EC queried','EC transferred from','Substrate used','Organism used', ...
    'Evidence tier'};
% Columns O..S are provenance only. They sit AFTER the Notes column so that
% every Excel formula below keeps referring to E..J by the same letters.
values = cell(height(Reactions)+1,19); values(1,:) = headers;
details = {'Reaction','Reaction ID','EC used','Substrate','On values sheet', ...
    'K_M [mM]','n','min','median','max','K_M taken from'};
saMW = {'Reaction','Reaction ID','EC in model','EC used','SA [μmol/min/mg]', ...
    'n','min','median','max','SA taken from','MW [Da]','n','min','median', ...
    'max','MW taken from','EC note'};
referenceChoices={'Reaction','Reaction ID','EC used','K_M [mM]','k_cat [s^-1]', ...
    'Km substrate','kcat substrate','Km organism','kcat organism', ...
    'Km literature','kcat literature','Km method','kcat method','Km n','kcat n', ...
    'Km assay','kcat assay','Km source','kcat source','Assumptions','Selected date'};
snapshot=referenceSnapshot();
rowReason=repmat("",height(Reactions),1);   % why a row ended up without values

for i=1:height(Reactions)
    row=i+1; R=Reactions(i,:);
    [substrates,tags] = splitSide(R.NamesEquation,'left');
    primary = ""; primarySource = "";
    if isfield(opt.PrimarySubstrate,char(R.Reaction))
        primary = string(opt.PrimarySubstrate.(char(R.Reaction)));
        assert(any(normalize(substrates)==normalize(primary)), ...
            'PrimarySubstrate for %s must be a printed reactant.',char(R.Reaction));
        primarySource = "caller";
    elseif sum(tags~="")==1
        primary = substrates(find(tags~="",1));
        primarySource = "single tagged reactant";
    elseif numel(substrates)==1
        primary = substrates(1); primarySource = "only reactant";
    elseif R.PathSubstrate~=""
        % Adapter chains print no tag on the consumed intermediate: the
        % Metabolite path / Reaction path lines identify it unambiguously.
        primary = R.PathSubstrate; primarySource = "adapter metabolite path";
    elseif sum(tags~="")>1
        % Cancellation and sequestration steps are genuinely bi-substrate.
        % Take the first tagged reactant and say so, rather than withholding.
        primary = substrates(find(tags~="",1));
        primarySource = "first of "+string(sum(tags~=""))+" tagged reactants";
    elseif all(tags=="")
        % Do not let water/protons or cofactors obscure the carrier.
        remaining=~ismember(substrateKey(substrates),currencyKeys());
        if sum(remaining)==1
            primary=substrates(remaining); primarySource="only non-currency reactant";
        end
    end
    codes=string(regexp(char(R.EC),'(?<![\d.])\d+\.\d+\.\d+\.\d+(?![\d.])','match'));
    % Translate every annotated code into the code BRENDA was queried under.
    % The motif file keeps the model's annotation; the records are filed under
    % the live code, so selection must use the live one.
    liveCodes=codes; transferNote=""; transferredFrom=strings(0,1);
    for t=1:numel(codes)
        key=char(codes(t));
        if isKey(ecMap,key)
            lc=string(ecMap(key));
            if lc~=codes(t)
                liveCodes(t)=lc; transferredFrom(end+1,1)=codes(t)+" -> "+lc; %#ok<AGROW>
                transferNote=transferNote+"EC "+codes(t)+" is NOT a current entry: " + ...
                    "ExPASy transfers it to EC "+lc+", and BRENDA holds the data only " + ...
                    "under the new number. All constants below were queried as EC "+lc+". ";
            end
        end
    end
    tiers=strings(0,1);
    if ~isempty(transferredFrom), tiers(end+1,1)="TRANSFERRED EC ("+strjoin(transferredFrom,", ")+")"; end
    substrateUsed=""; organismUsed="";
    preferredEC=""; preferenceNote="";
    if isfield(opt.PreferredEC,char(R.Reaction))
        preferredEC=string(opt.PreferredEC.(char(R.Reaction)));
        assert(isscalar(preferredEC) && any(codes==preferredEC), ...
            'PreferredEC for %s must be one of its annotated EC codes.',char(R.Reaction));
        preferenceNote="Explicit EC preference: "+preferredEC+". ";
        preferredEC=liveCodes(find(codes==preferredEC,1));  % query the live code
    elseif R.Reaction=="r_adapter_1" && ...
            substrateKey(aliasName(primary,aliases))=="camp" && any(codes=="3.1.4.17")
        preferenceNote="Default cAMP adapter EC preference: 3.1.4.17 (override with PreferredEC). ";
        preferredEC=liveCodes(find(codes=="3.1.4.17",1));
    end
    candidates = {}; chosenKm=NaN; chosenKcat=NaN; chosenEC=""; ref=[];
    kmProvenance=""; kcatProvenance="";
    for ec=liveCodes
        for s=substrates'
            % Try the printed model name AND its alias: an alias that misses
            % must not hide a record the original name would have matched.
            names=unique([s;aliasName(s,aliases)],'stable');
            ka=Km([],:); kb=Kcat([],:); pa=""; pb="";
            for q=names'
                if isempty(ka) || ~any(isfinite(ka.value))
                    [ka,pa]=selectRecords(Km,ec,q,string(opt.Organism), ...
                        opt.RelaxedSubstrateMatch,opt.OrganismFallback);
                end
                if isempty(kb) || ~any(isfinite(kb.value))
                    [kb,pb]=selectRecords(Kcat,ec,q,string(opt.Organism), ...
                        opt.RelaxedSubstrateMatch,opt.OrganismFallback);
                end
            end
            v=ka.value(isfinite(ka.value));
            selected = normalize(s)==normalize(primary) && ...
                (preferredEC=="" || ec==preferredEC);
            onSheet='no';
            if selected && ~isempty(v) && (R.Direction==1 || opt.AllowReverse), onSheet='candidate'; end
            if isempty(v), stats={[],0,[],[],[]};
            else, stats={median(v),numel(v),min(v),median(v),max(v)}; end
            details(end+1,:) = [{char(R.Reaction),char(R.ID),char(ec),char(s),onSheet}, ...
                stats,{char(strjoin(unique(ka.organism,'stable'),'; '))}]; %#ok<AGROW>
            c=kb.value(isfinite(kb.value));
            if selected && (~isempty(v)||~isempty(c))
                vm=NaN; cm=NaN; subUsed=s; orgUsed="";
                if ~isempty(v)
                    vm=median(v); kmProvenance=pa;
                    subUsed=strjoin(unique(ka.substrate,'stable'),'; ');
                    orgUsed=strjoin(unique(ka.organism,'stable'),'; ');
                elseif ~isempty(c)
                    orgUsed=strjoin(unique(kb.organism,'stable'),'; ');
                end
                if ~isempty(c), cm=median(c); kcatProvenance=pb; end
                candidates(end+1,:)={ec,vm,cm,numel(v),numel(c),subUsed,orgUsed}; %#ok<AGROW>
            end
        end
    end
    % Keep only the preferred EC when one is set, so a second annotated EC
    % does not force the summary to be withheld.
    if preferredEC~="" && ~isempty(candidates)
        candidates = candidates(string(candidates(:,1))==preferredEC,:);
    end
    note="Candidate medians require assay, enzyme and direction review. ";
    if ~opt.QueryBRENDA, note="Offline template: no BRENDA measurements queried. "; end
    note=note+preferenceNote;
    note=note+R.ECLookupNote+" ";
    if R.ECSource~="TXT", tiers(end+1,1)=R.ECSource; end
    if R.DirectionBasis=="printed left-to-right only"
        note=note+"Direction not supplied: +1 denotes the printed equation only; " + ...
            "orientation relative to model.S and BRENDA was not established. ";
    end
    if R.PrintedHeading~=""
        note=note+"TXT heading: "+R.PrintedHeading+". ";
    end
    if R.Direction~=1 && opt.AllowReverse
        % Filled on request. BRENDA measures the direction the enzyme is
        % written in, so these constants describe the forward reaction and
        % transfer to the reverse only through K_eq (Haldane relation).
        note=note+"REVERSE in the model (Direction = "+string(R.Direction)+ ...
            "): constants below are BRENDA's FORWARD values, reused without a " + ...
            "K_eq correction. Verify before quoting. ";
    end
    if R.Direction~=1 && ~opt.AllowReverse
        note=note+"NOT forward relative to model.S (Direction = "+string(R.Direction)+ ...
            "): Km/kcat withheld. BRENDA measures the forward direction, and forward " + ...
            "constants do not transfer to the reverse without K_eq (Haldane). " + ...
            "Set AllowReverse = true to fill the row anyway. ";
        rowReason(i)="reverse direction";
    elseif ~isempty(candidates)
        % Several annotated EC codes can match. Prefer the one with the most
        % Km records rather than withholding the row, and record the choice.
        if size(candidates,1)>1
            [~,best]=max(cell2mat(candidates(:,4)));
            note=note+sprintf('%d EC candidates; kept the one with most Km records (set PreferredEC to override). ', ...
                size(candidates,1));
        else
            best=1;
        end
        chosenEC=candidates{best,1}; chosenKm=candidates{best,2}; chosenKcat=candidates{best,3};
        substrateUsed=candidates{best,6}; organismUsed=candidates{best,7};
        note=note+"Primary substrate: "+primary+" ("+primarySource+"); EC used: "+chosenEC+". ";
        note=note+sprintf('Km n=%d; kcat n=%d. ',candidates{best,4},candidates{best,5});
        if kmProvenance~=""
            note=note+"Km match: "+kmProvenance+". "; tiers(end+1,1)="K_M by "+kmProvenance;
        end
        if kcatProvenance~=""
            note=note+"kcat match: "+kcatProvenance+". "; tiers(end+1,1)="k_cat by "+kcatProvenance;
        end
        note=note+"Medians may pool different isoforms/assays; not a paired enzyme-specific measurement. ";
        % k_cat is often absent for the exact substrate even where K_M exists.
        % Fall back to the same EC across all substrates, clearly flagged.
        if ~isfinite(chosenKcat) && opt.KcatAnySubstrate
            pool=Kcat(Kcat.EC==string(chosenEC),:);
            if ~opt.OrganismFallback || string(opt.Organism)~=""
                sameOrg=pool(normalize(pool.organism)==normalize(string(opt.Organism)),:);
                if any(isfinite(sameOrg.value)), pool=sameOrg; end
            end
            cc=pool.value(isfinite(pool.value));
            if ~isempty(cc)
                chosenKcat=median(cc);
                note=note+sprintf(['k_cat UNMATCHED to the primary substrate: median of %d records ' ...
                    'for EC %s across all substrates. Treat as an order-of-magnitude surrogate. '], ...
                    numel(cc),chosenEC);
                tiers(end+1,1)="k_cat UNMATCHED substrate (any ligand of EC "+chosenEC+")";
            end
        end
    elseif isempty(codes)
        note=note+"NO EC CODE annotated for this reaction in the motif file: nothing can be queried. ";
        rowReason(i)="no EC code in the motif file";
        tiers(end+1,1)="no EC code";
    elseif opt.SubstituteKm && primary~="" && (R.Direction==1 || opt.AllowReverse)
        % BRENDA often has no record for the ligand this reaction consumes
        % (EC 2.7.4.6 has ADP/ATP/GTP/dTDP, never dGDP). Rather than leaving
        % the row empty, take the best-sampled ligand of the same EC and say
        % so: this is the case the "*" in the "Removed by" header refers to.
        wanted=primary; surrogateSource="nearest native ligand by name similarity";
        if isfield(opt.SurrogateSubstrate,char(R.Reaction))
            wanted=string(opt.SurrogateSubstrate.(char(R.Reaction)));
            surrogateSource="caller-specified surrogate";
        end
        [subKm,subLigand,subN,subEC,subOrg] = substituteKm(Km,liveCodes, ...
            string(opt.Organism),opt.OrganismFallback,wanted);
        if isfinite(subKm)
            chosenKm=subKm; chosenEC=subEC;
            substrateUsed=subLigand; organismUsed=subOrg;
            note=note+sprintf(['K_M SUBSTITUTED: BRENDA holds no record for %s under EC %s. ' ...
                'Using the ligand %s instead (%d records; %s), a DIFFERENT substrate of the ' ...
                'SAME enzyme. Chromogenic and fluorogenic assay substrates were excluded. ' ...
                'This K_M describes that surrogate, NOT %s: declare it in the methods and ' ...
                'mark this row with "*". '], ...
                primary,subEC,subLigand,subN,surrogateSource,primary);
            tiers(end+1,1)="K_M SUBSTITUTED ligand ("+subLigand+")";
        else
            note=note+"No record for "+primary+" nor any other ligand of EC "+strjoin(cellstr(liveCodes),'/')+". ";
            rowReason(i)="EC has no usable K_M record at all";
            tiers(end+1,1)="no K_M record under any ligand";
        end
    elseif primary==""
        note=note+"Primary substrate unresolved: set PrimarySubstrate or review KM details. ";
        rowReason(i)="primary substrate unresolved";
        tiers(end+1,1)="primary substrate unresolved";
    else
        note=note+"No matched records for "+primary+" under EC "+strjoin(cellstr(liveCodes),'/')+ ...
            "; review Km records, Query status and SubstrateAliases. ";
        tiers(end+1,1)="no matched K_M record";
    end
    if opt.UseReferenceValues
        hit=find(string({snapshot.reaction})==R.Reaction & string({snapshot.id})==R.ID);
        if numel(hit)==1
            candidate=snapshot(hit);
            compatible=any([codes liveCodes]==string(candidate.ec)) && ...
                normalize(candidate.equation)==normalize(R.NamesEquation) && ...
                candidate.direction==R.Direction;
            % Explicit caller choices take precedence over the fixed snapshot.
            if isfield(opt.PrimarySubstrate,char(R.Reaction)), compatible=false; end
            if preferredEC~="" && preferredEC~=string(candidate.ec), compatible=false; end
            if compatible
                ref=candidate;
                chosenKm=ref.km.value; chosenKcat=ref.kcat.value; chosenEC=string(ref.ec);
                substrateUsed=string(ref.km.substrate);
                organismUsed="Km: "+string(ref.km.organism)+"; kcat: "+string(ref.kcat.organism);
                note=string(ref.note)+" "+R.ECLookupNote;
                referenceChoices(end+1,:)={ref.reaction,ref.id,ref.ec,chosenKm,chosenKcat, ...
                    ref.km.substrate,ref.kcat.substrate,ref.km.organism,ref.kcat.organism, ...
                    ref.km.literature,ref.kcat.literature,ref.km.kind,ref.kcat.kind, ...
                    ref.km.n,ref.kcat.n,ref.km.conditions,ref.kcat.conditions, ...
                    ref.km.source,ref.kcat.source,ref.assumptions,'2026-10-12'}; %#ok<AGROW>
                % Mark the original matching view as superseded, then list
                % the exact reference used (including product-side references).
                for d=2:size(details,1)
                    if strcmp(details{d,1},char(R.Reaction)) && strcmp(details{d,2},char(R.ID))
                        details{d,5}='not selected (reference mode)';
                    end
                end
                details(end+1,:)={char(R.Reaction),char(R.ID),ref.ec,ref.km.substrate, ...
                    'reference selected',chosenKm,ref.km.n,[],chosenKm,[], ...
                    [ref.km.organism '; ' ref.km.kind '; literature ' ref.km.literature]}; %#ok<AGROW>
            end
        end
    end
    if R.Reaction=="r_Nw"
        note=note+"m_c activation is a design assumption; enzyme Km does not specify its activation law. ";
    end
    if numel(codes)>1, note=note+"Multiple model EC annotations retained without correction. "; end
    values(row,1:4)={sprintf('Motif %d / Block %d',R.Motif,R.Block), ...
        char(R.Reaction),char(R.Name),char(R.EC)};
    if isfinite(chosenKm), values{row,5}=chosenKm; end
    % Specific activity and molecular weight give k_cat = SA*MW/60 whenever a
    % direct turnover number is missing. BRENDA reports both per EC, so this
    % recovers rows the turnover query alone leaves blank.
    % Specific activity and molecular weight are reported per enzyme, not per
    % substrate, so they can be filled from the EC code even when no K_M
    % record matched the printed substrate name.
    % 'SAxMW' (the default) reproduces the original workbook: k_cat is always
    % computed as SA*MW/60 and a BRENDA turnover number is never written to
    % column H. It is still reported in the note so the two can be compared.
    useSAxMW = strcmpi(char(opt.KcatSource),'SAxMW');
    if useSAxMW && isfinite(chosenKcat)
        note=note+sprintf('BRENDA turnover number %.4g 1/s was found but NOT used (KcatSource = SAxMW). ',chosenKcat);
        chosenKcat=NaN;
    end
    saValue=NaN; mwValue=NaN; ecSAMW=chosenEC;
    if ecSAMW==""
        if preferredEC~=""
            ecSAMW=preferredEC;
        elseif numel(liveCodes)==1
            ecSAMW=liveCodes(1);
        elseif ~isempty(liveCodes)
            have=arrayfun(@(e) sum(SArec.EC==e & isfinite(SArec.value)),liveCodes);
            [most,which]=max(have);
            if most>0, ecSAMW=liveCodes(which); end
        end
    end
    if ecSAMW~=""
        [saValue,saOrg,saN]=medianFor(SArec,ecSAMW,string(opt.Organism),opt.OrganismFallback);
        [mwValue,mwOrg,mwN]=medianFor(MWrec,ecSAMW,string(opt.Organism),opt.OrganismFallback);
        if ~isfinite(chosenKcat) && isfinite(saValue) && isfinite(mwValue)
            note=note+sprintf(['k_cat ESTIMATED as SA*MW/60 from %g umol/min/mg (n=%d; %s) and ' ...
                '%g Da (n=%d; %s), medians for EC %s. '], ...
                saValue,saN,saOrg,mwValue,mwN,mwOrg,ecSAMW);
            note=note+mwSubunitHint(MWrec,ecSAMW);
            tiers(end+1,1)="k_cat ESTIMATED as SA*MW/60 (EC "+ecSAMW+")";
            if organismUsed=="", organismUsed="SA: "+saOrg+"; MW: "+mwOrg; end
        end
    end
    if isfinite(saValue), values{row,6}=saValue; end
    if isfinite(mwValue), values{row,7}=mwValue/1000; end % Da -> mg/umol
    if ~isempty(ref) && ~isempty(ref.sa) && ~useSAxMW
        values{row,6}=ref.sa; values{row,7}=ref.mw;
        values{row,8}=sprintf('=F%d*G%d/60',row,row);
    elseif isfinite(chosenKcat), values{row,8}=chosenKcat;
    elseif isfinite(saValue) && isfinite(mwValue)
        values{row,6}=saValue; values{row,7}=mwValue/1000;   % mg/umol
        values{row,8}=sprintf('=F%d*G%d/60',row,row);
    else
        % Retain template conversion for optional manual SA/MW entry only.
        values{row,8}=sprintf('=IF(AND(ISNUMBER(F%d),ISNUMBER(G%d)),F%d*G%d/60,"")',row,row,row,row);
    end
    values{row,9}=0; % User supplies enzyme concentration, as in original workbook.
    values{row,10}=sprintf(['=IF(AND(ISNUMBER(H%d),ISNUMBER(I%d),ISNUMBER(E%d)),' ...
        'IF(E%d>0,H%d*I%d/E%d,""),"")'],row,row,row,row,row,row,row);
    if ~isempty(ref), tiers(end+1,1)="REFERENCE SNAPSHOT (2026-10-12), not a live query"; end
    if ~opt.QueryBRENDA, tiers(end+1,1)="offline: nothing queried"; end
    if isempty(tiers), tiers="no value"; end
    values{row,14}=char(transferNote+note);
    % Provenance of the numbers in E..H, so a reader can filter the workbook
    % instead of parsing the Notes prose.
    values{row,15}=char(chosenEC);
    if isempty(transferredFrom), values{row,16}='';
    else, values{row,16}=char(strjoin(transferredFrom,'; ')); end
    if substrateUsed=="" && primary~="", substrateUsed=primary; end
    values{row,17}=char(substrateUsed);
    values{row,18}=char(organismUsed);
    values{row,19}=char(strjoin(unique(tiers,'stable'),' | '));
    saMW(end+1,:) = [{char(R.Reaction),char(R.ID),char(R.OriginalEC),char(chosenEC)}, ...
        cell(1,12),{'No selected SA/MW records; see Query status and raw records.'}]; %#ok<AGROW>
    if ecSAMW~=""
        saMW{end,4}=char(ecSAMW);
        saMW{end,17}='SA and MW selected independently; medians may represent different assays.';
        for parameterIndex=1:2
            if parameterIndex==1, pool=SArec; firstColumn=5;
            else, pool=MWrec; firstColumn=11; end
            [v,org,n]=medianFor(pool,ecSAMW,string(opt.Organism),opt.OrganismFallback);
            subset=pool(pool.EC==ecSAMW & isfinite(pool.value),:);
            same=subset(normalize(subset.organism)==normalize(string(opt.Organism)),:);
            if ~isempty(same), subset=same;
            elseif ~opt.OrganismFallback && string(opt.Organism)~="", subset=subset([],:); end
            if isfinite(v)
                saMW(end,firstColumn:firstColumn+5)={v,n,min(subset.value),v,max(subset.value),char(org)};
            end
        end
    end
    if ~isempty(ref)
        saMW{end,17}=['Reference choice: ' ref.kcat.kind '; see Reference choices.'];
        if ~isempty(ref.sa)
            saMW(end,5:16)={ref.sa,1,ref.sa,ref.sa,ref.sa,ref.kcat.source, ...
                ref.mw*1000,1,ref.mw*1000,ref.mw*1000,ref.mw*1000,ref.kcat.source};
        end
    end
end

% Only first-order pool producer/consumer pairs receive G/tau formulas.
% Cancellation and sequestration are multi-substrate steps: summing their
% kcat[E]/Km as first-order losses would require co-substrate assumptions.
% October 2026 nomenclature: each delay pool is produced by its feed
% reaction and consumed by its delay reaction (m_F by r_Ff then r_Fd, etc.).
pairs={'r_Ff','r_Fd';'r_Sf','r_Sd';'r_Ff_2','r_Fd_2';'r_Sf_2','r_Sd_2'};
for motif=unique(Reactions.Motif)'
    for j=1:size(pairs,1)
        a=find(Reactions.Motif==motif & Reactions.Reaction==string(pairs{j,1}));
        b=find(Reactions.Motif==motif & Reactions.Reaction==string(pairs{j,2}));
        if numel(a)==1 && numel(b)==1
            ra=a+1; rb=b+1;
            values{ra,11}=pairs{j,2};
            values{ra,12}=sprintf('=IFERROR(IF(COUNT(J%d,J%d)=2,J%d/J%d,""),"")',ra,rb,ra,rb);
            values{ra,13}=sprintf('=IFERROR(IF(ISNUMBER(J%d),1/J%d,""),"")',rb,rb);
            values{ra,14}=[values{ra,14} 'G/tau assume buffered co-substrates and first-order pool kinetics.'];
        end
    end
end
Values=cell2table(values(2:end,:),'VariableNames', ...
    {'ReactionPair','Reaction','Enzyme','EC','Km','SA','MW','Kcat','E','k', ...
     'RemovedBy','G','tau','Notes','ECqueried','ECtransferredFrom', ...
     'SubstrateUsed','OrganismUsed','EvidenceTier'});

% Coverage report: name every reaction that is still missing a constant, so
% an empty workbook cell is never silent.
isNum=@(v)isnumeric(v)&&isscalar(v)&&isfinite(v);
haveKm=cellfun(isNum,values(2:end,5));
haveKcat=cellfun(isNum,values(2:end,8)) | ...
    (cellfun(isNum,values(2:end,6)) & cellfun(isNum,values(2:end,7)));
fprintf('\nCoverage: K_M %d/%d, k_cat %d/%d, both %d/%d.\n', ...
    sum(haveKm),numel(haveKm),sum(haveKcat),numel(haveKcat), ...
    sum(haveKm&haveKcat),numel(haveKm));
gaps=find(~(haveKm&haveKcat));
for g=gaps(:)'
    missing="K_M and k_cat";
    if haveKm(g), missing="k_cat"; elseif haveKcat(g), missing="K_M"; end
    why=rowReason(g);
    if why=="", why="see the Notes column"; end
    fprintf('  %-14s EC %-12s missing %-14s (%s)\n', ...
        values{g+1,2},values{g+1,4},missing,why);
end
if opt.UseReferenceValues
    fprintf('Reference snapshot applied to %d of %d rows (matched on reaction ID, equation and direction).\n', ...
        size(referenceChoices,1)-1,height(Reactions));
end
transferred=find(cellfun(@(a,b)~strcmp(a,b),ecResolution(2:end,1),ecResolution(2:end,2)));
if ~isempty(transferred)
    fprintf('\nObsolete EC codes resolved through ExPASy (%d):\n',numel(transferred));
    for t=transferred(:)'
        fprintf('  %-10s -> %-10s %s (%d record(s) recovered)\n', ...
            ecResolution{t+1,1},ecResolution{t+1,2},ecResolution{t+1,7},ecResolution{t+1,6});
    end
    fprintf('Correct these codes in the model annotation: the motif file still carries the old ones.\n');
end

% Always write a fresh workbook: stale old-motif rows/sheets must not survive.
folder=fileparts(char(xlsxFile)); if isempty(folder), folder=pwd; end
assert(isfolder(folder),'Output folder does not exist.');
tempFile=[tempname(folder) '.xlsx'];
cleanup=onCleanup(@()deleteIfExists(tempFile));
writecell(values,tempFile,'Sheet','values');
writecell(details,tempFile,'Sheet','KM details');
writecell(saMW,tempFile,'Sheet','SA MW details');
writetable(Km,tempFile,'Sheet','Km records');
writetable(Kcat,tempFile,'Sheet','kcat records');
writetable(SArec,tempFile,'Sheet','SA records');
writetable(MWrec,tempFile,'Sheet','MW records');
writetable(Reactions,tempFile,'Sheet','Reaction mapping');
writecell(status,tempFile,'Sheet','Query status');
writecell(ecResolution,tempFile,'Sheet','EC resolution');
writecell(referenceChoices,tempFile,'Sheet','Reference choices');
installExcelFormulas(tempFile,values);
[ok,msg]=movefile(tempFile,char(xlsxFile),'f'); assert(ok,'%s',msg);
fprintf('Workbook written: %s\n',char(xlsxFile));
end

function T=parseMotifs(txt)
% ID fields, rather than role headings, delimit reaction occurrences. This
% handles repeated r_CAMP headings, abbreviated branches and legacy exports.
template=struct('Motif',NaN,'Block',NaN,'Reaction',"",'Role',"", ...
    'ID',"",'Name',"",'GPR',"",'EC',"",'Compartment',"", ...
    'Direction',NaN,'FormulaEquation',"",'NamesEquation',"",'IDsEquation',"", ...
    'PathSubstrate',"",'PrintedHeading',"",'Bounds',"",'DirectionBasis',"");
records=repmat(template,0,1);
lines=regexp(txt,'\r\n|\n|\r','split'); motif=1; block=NaN;
metPath=containers.Map('KeyType','double','ValueType','any');
rxnPath=containers.Map('KeyType','double','ValueType','any');
current=struct(); pending="";
for i=1:numel(lines)
    line=lines{i}; trimmed=strtrim(line);
    m=regexp(line,'^\s*MOTIF\s+(\d+)\s*$','tokens','once');
    b=regexp(line,'^\s*BLOCK\s+(\d+)\s*:','tokens','once');
    isEnd=strcmp(trimmed,'ASSUMPTIONS AND INTERFACES');
    if ~isempty(m)||~isempty(b)||isEnd
        if ~isempty(fieldnames(current)), records(end+1)=current; current=struct(); end %#ok<AGROW>
        pending="";
    end
    if isEnd, break; end
    if ~isempty(m), motif=str2double(m{1}); block=NaN; end
    if ~isempty(b), block=str2double(b{1}); end
    pm=regexp(line,'^\s*Metabolite path\s*:\s*(.*)$','tokens','once');
    pr=regexp(line,'^\s*Reaction path\s*:\s*(.*)$','tokens','once');
    if ~isempty(pm), metPath(motif)=string(regexp(strtrim(pm{1}),'\s*->\s*','split')); end
    if ~isempty(pr), rxnPath(motif)=string(regexp(strtrim(pr{1}),'\s*->\s*','split')); end
    f=regexp(line,'^\s*(ID|Name|GPR|EC|Compartment\(s\)|Direction|Eqn|Names|IDs|Bounds)\s*:\s*(.*)$','tokens','once');
    if ~isempty(f) && strcmp(f{1},'ID')
        if ~isempty(fieldnames(current)), records(end+1)=current; end %#ok<AGROW>
        current=template; current.Motif=motif; current.Block=block;
        current.PrintedHeading=pending;
        role=regexp(char(pending),'^(r_[A-Za-z0-9_]+)\s*(?:\((.*)\))?$','tokens','once');
        if ~isempty(role)
            current.Reaction=string(role{1});
            if numel(role)>1, current.Role=string(role{2}); end
        else
            current.Role=pending;
        end
        pending="";
    end
    if ~isempty(f) && ~isempty(fieldnames(current))
        switch f{1}
            case 'Compartment(s)', field='Compartment';
            case 'Eqn', field='FormulaEquation';
            case 'Names', field='NamesEquation';
            case 'IDs', field='IDsEquation';
            otherwise, field=f{1};
        end
        if strcmp(field,'Direction')
            d=regexp(f{2},'^[+-]?\d+','match','once');
            current.Direction=str2double(d); current.DirectionBasis="explicit TXT Direction";
        else, current.(field)=string(strtrim(f{2})); end
    elseif isempty(f) && isfinite(block) && ~isempty(trimmed) && ...
            isempty(regexp(trimmed,'^[-#]+$','once')) && ...
            isempty(regexp(trimmed,':','once')) && isempty(b) && isempty(m)
        % The last standalone heading before an ID belongs to that record.
        pending=string(trimmed);
    end
end
if ~isempty(fieldnames(current)), records(end+1)=current; end
if isempty(records), T=table(); return; end
T=struct2table(records);
assert(all(T.ID~="" & T.NamesEquation~="" & isfinite(T.Block)), ...
    'Incomplete reaction record: every ID needs a block and Names equation.');
% Allocate deterministic labels within each motif/block; do not infer missing
% feed/delay roles from "fast branch", "slow branch" or "...".
for i=1:height(T)
    if T.Reaction(i)==""
        step=sum(T.Motif(1:i)==T.Motif(i) & T.Block(1:i)==T.Block(i));
        T.Reaction(i)="r_B"+string(T.Block(i))+"_"+string(step);
    end
    if isnan(T.Direction(i))
        T.Direction(i)=1; T.DirectionBasis(i)="printed left-to-right only";
    end
    if T.Compartment(i)==""
        tokens=regexp(char(T.IDsEquation(i)),'\[([^\]]+)\]','tokens');
        comps=string(cellfun(@(v)v{1},tokens,'UniformOutput',false));
        T.Compartment(i)=strjoin(unique(comps,'stable'),'; ');
    end
end
base=T.Reaction;
for i=1:height(T)
    same=T.Motif==T.Motif(i) & base==base(i);
    if sum(same)>1
        T.Reaction(i)=base(i)+"_"+string(sum(same(1:i)));
    end
end
for k=unique(T.Motif)'
    labels=T.Reaction(T.Motif==k);
    assert(numel(unique(labels))==numel(labels), ...
        'Generated reaction labels collide; give the TXT reactions unique r_ labels.');
end
for i=1:height(T)
    k=T.Motif(i);
    if T.Block(i)~=0 || ~isKey(rxnPath,k) || ~isKey(metPath,k), continue; end
    rp=rxnPath(k); mp=metPath(k); pos=find(rp==T.ID(i));
    if numel(pos)~=1 || pos>=numel(mp), continue; end
    T.PathSubstrate(i)=nameForId(T.NamesEquation(i),T.IDsEquation(i),mp(pos));
end
end

function name=nameForId(namesEquation,idsEquation,metaboliteId)
% The Names and IDs equations are printed with the same term order, so the
% position of an ID gives the printed name of the same species.
name="";
try
    [ids,~]  = splitSide(idsEquation,'left');
    [nms,~]  = splitSide(namesEquation,'left');
catch
    return
end
if numel(ids)~=numel(nms), return; end
hit=find(ids==metaboliteId,1);
if ~isempty(hit), name=nms(hit); end
end

function [names,tags]=splitSide(equation,side)
halves=regexp(char(equation),'\s*->\s*','split');
assert(numel(halves)==2,'Expected a printed -> equation.');
if strcmp(side,'left'), part=halves{1}; else, part=halves{2}; end
% Split spaced plus delimiters; preserve H+ and NADP(+).
parts=regexp(part,'\s+\+\s+','split'); names=strings(numel(parts),1); tags=names;
for i=1:numel(parts)
    s=strtrim(parts{i}); tag=regexp(s,'\((m_[A-Za-z0-9_]+)\)\s*$','tokens','once');
    if ~isempty(tag), tags(i)=string(tag{1}); s=regexprep(s,'\(m_[A-Za-z0-9_]+\)\s*$',''); end
    names(i)=string(regexprep(strtrim(s),'^\d+(?:\.\d+)?\s+',''));
end
end

function name=aliasName(name,aliases)
hit=find(normalize(aliases(:,1))==normalize(name));
assert(numel(hit)<=1,'Duplicate substrate alias.');
if ~isempty(hit), name=aliases(hit,2); end
end
function s=normalize(s)
s=lower(strtrim(string(s))); s=regexprep(s,'\s+',' ');
end
function s=substrateKey(s)
% Explicit synonym groups, grounded in the model and BRENDA record names.
% Keep cAMP and cGMP separate; do not use fuzzy/punctuation-free matching.
s=normalize(s);
s(ismember(s,["camp","3',5'-camp","3',5'-cyclic amp"]))="camp";
s(ismember(s,["cgmp","3',5'-cgmp","3',5'-cyclic gmp"]))="cgmp";
end
function [T,provenance]=selectRecords(T,ec,substrate,organism,relaxed,organismFallback)
% Ranked selection. A single hard filter on EC + exact ligand name + exact
% organism leaves most rows empty, because many EC codes have no record for
% the requested organism at all and model metabolite names rarely equal
% BRENDA ligand names. Try the strictest combination first and relax in a
% fixed order, reporting which step produced the records.
provenance="no records";
isEC = T.EC==ec;
exactSub = isEC & substrateKey(T.substrate)==substrateKey(substrate);
if relaxed, looseSub = isEC & looseMatch(T.substrate,substrate);
else,       looseSub = false(height(T),1); end
sameOrg = normalize(T.organism)==normalize(organism);
attempts = {exactSub & sameOrg,'exact substrate, requested organism'; ...
            exactSub,          'exact substrate, ANY organism'; ...
            looseSub & sameOrg,'relaxed substrate name, requested organism'; ...
            looseSub,          'relaxed substrate name, ANY organism'};
if organism=="",       attempts=attempts([2 4],:); end
if ~organismFallback,  attempts=attempts(~contains(string(attempts(:,2)),"ANY organism"),:); end
for a=1:size(attempts,1)
    mask=attempts{a,1};
    if any(mask & isfinite(T.value))
        provenance=string(attempts{a,2}); T=T(mask,:); return
    end
end
T=T([],:);
end

function [m,organisms,n]=medianFor(T,ec,organism,organismFallback)
% Median of a substrate-free parameter (SA, MW) for one EC, preferring the
% requested organism and falling back to all organisms.
m=NaN; organisms=""; n=0;
rows=T(T.EC==string(ec),:);
if isempty(rows), return; end
if organism~=""
    same=rows(normalize(rows.organism)==normalize(organism),:);
    if any(isfinite(same.value)), rows=same;
    elseif ~organismFallback, return; end
end
rows=rows(isfinite(rows.value),:);
if isempty(rows), return; end
m=median(rows.value); n=height(rows);
organisms=strjoin(unique(rows.organism,'stable'),'; ');
end

function text=mwSubunitHint(MWrec,ec)
% k_cat = SA*MW/60 needs the mass that carries ONE active site. BRENDA mixes
% native-oligomer and subunit masses in the same field, and the commentary is
% the only place that says which. A native dimer mass inflates k_cat twofold,
% so the count is reported rather than silently averaged.
rows=MWrec(MWrec.EC==string(ec) & isfinite(MWrec.value),:);
if isempty(rows), text=""; return; end
c=lower(rows.commentary);
subunit=sum(contains(c,"subunit")|contains(c,"monomer")|contains(c,"sds"));
native=sum(contains(c,"native")|contains(c,"dimer")|contains(c,"tetramer")| ...
           contains(c,"oligomer")|contains(c,"gel filtration"));
text=sprintf(['CHECK THE MASS BASIS: of %d MW records, %d mention a subunit/monomer ' ...
    'and %d a native/oligomeric mass. SA*MW/60 is only a per-site turnover if MW is ' ...
    'the mass bearing ONE active site; a native oligomer mass overestimates k_cat by ' ...
    'the number of subunits (human GGCT, EC 4.3.2.9, is a homodimer of ~21 kDa ' ...
    'subunits, so the native mass doubles it). See the MW records sheet commentary. '], ...
    height(rows),subunit,native);
end

function tf=looseMatch(candidates,wanted)
% Compare the printed name, and also the abbreviation built from it, because
% genome-scale models print "Deoxycytidine-5'-Triphosphate" where BRENDA
% lists "dCTP". Both spellings of both sides are tried.
a=compactName(candidates);                b=compactName(wanted);
aAbbrev=compactName(nucleotideAbbrev(candidates));
bAbbrev=compactName(nucleotideAbbrev(wanted));
tf=false(size(a));
for i=1:numel(a)
    tf(i) = pairMatch(a(i),b)       || pairMatch(aAbbrev(i),b) || ...
            pairMatch(a(i),bAbbrev) || pairMatch(aAbbrev(i),bAbbrev);
end
end

function tf=pairMatch(x,y)
% Equality at any length; containment only when the shorter name is long
% enough to be specific, so that "ATP" never matches "dATP".
tf=false;
if strlength(x)<2 || strlength(y)<2, return; end
if x==y, tf=true; return; end
short=x; long=y;
if strlength(x)>strlength(y), short=y; long=x; end
tf = strlength(short)>=6 && contains(long,short);
end

function out=nucleotideAbbrev(names)
% Build the BRENDA-style abbreviation of a nucleotide name. Returns "" for
% anything that is not a recognisable nucleoside phosphate.
names=string(names); out=strings(size(names));
bases=["adenosine","A";"guanosine","G";"cytidine","C";"uridine","U"; ...
       "thymidine","T";"inosine","I";"xanthosine","X"];
for i=1:numel(names)
    s=regexprep(lower(string(names(i))),'[^a-z0-9]',' ');
    deoxy=contains(s,"deoxy");  s=erase(s,"deoxy");
    b="";
    for j=1:size(bases,1)
        if contains(s,bases(j,1)), b=bases(j,2); break; end
    end
    if b=="", continue; end
    if     contains(s,"triphosphate"),  p="TP";
    elseif contains(s,"diphosphate"),   p="DP";
    elseif contains(s,"monophosphate"), p="MP";
    elseif contains(s,"phosphate"),     p="MP";
    else,  continue; end
    if deoxy, out(i)="d"+b+p; else, out(i)=b+p; end
end
end

function s=compactName(s)
s=lower(string(s));
s=regexprep(s,'^(l|d|dl|r|s)[-\s]','');
s=regexprep(s,'ic\s+acid\s*$','ate');     % "succinic acid" and "succinate" are one ligand
s=regexprep(s,'[^a-z0-9]','');
end

function [value,ligand,n,ecUsed,organisms]=substituteKm(Km,codes,organism,organismFallback,wanted)
% Surrogate ligand of the SAME enzyme, for the case where BRENDA has no record
% at all for the substrate this reaction consumes.
%
% Picking the best-sampled ligand alone is not good enough. EC 4.3.2.9 (the
% successor of 2.3.2.4) has 40 distinct Km ligands and no gamma-glutamyl-L-
% cysteine at all; the list is dominated by synthetic amines and chromogenic
% assay substrates (4-nitroanilides, methylumbelliferyl and phenoxazinyl
% esters, benzylamine, neohexylamine). Those report the kinetics of the assay,
% not of the native step, so they are excluded outright.
%
% Among what is left, a ligand is scored by how much of the wanted name it
% shares (word overlap on the chemical name), and record count breaks ties.
% For gamma-glutamyl-L-cysteine that selects 5-L-glutamyl-L-alpha-
% aminobutyrate / gamma-glutamyl-L-alanine: same gamma-glutamyl donor, a
% near-isosteric acceptor side chain.
if nargin<5, wanted=""; end
value=NaN; ligand=""; n=0; ecUsed=""; organisms=""; bestScore=-Inf;
for ec=codes
    pool=Km(Km.EC==ec & isfinite(Km.value),:);
    if isempty(pool), continue; end
    if organism~="" && ~organismFallback
        pool=pool(normalize(pool.organism)==normalize(organism),:);
    elseif organism~=""
        same=pool(normalize(pool.organism)==normalize(organism),:);
        if ~isempty(same), pool=same; end
    end
    if isempty(pool), continue; end
    names=unique(pool.substrate,'stable');
    names=names(~isArtificialLigand(names));
    if isempty(names), continue; end
    for k=1:numel(names)
        count=sum(pool.substrate==names(k));
        % Similarity dominates; the record count only orders equal matches.
        score=100*nameOverlap(names(k),wanted)+min(count,20);
        if score>bestScore
            bestScore=score; n=count; ligand=names(k); ecUsed=ec;
            rows=pool(pool.substrate==ligand,:);
            value=median(rows.value);
            organisms=strjoin(unique(rows.organism,'stable'),'; ');
        end
    end
end
end

function tf=isArtificialLigand(names)
% Chromogenic/fluorogenic reporters and BRENDA's own placeholder. These
% measure the assay, not the metabolic step, and must never set a K_M used in
% a simulation of the native reaction.
markers=["nitroanilide","nitroanilin","umbellifer","phenoxazin","coumarin", ...
         "fluoresc","naphthyl","nitrophenyl","anilide","benzyloxycarbonyl", ...
         "carbobenzoxy","dansyl","resorufin"];
s=lower(string(names));
tf = s=="more" | s=="" | s=="?";
for m=markers, tf = tf | contains(s,m); end
end

function f=nameOverlap(candidate,wanted)
% Fraction of the wanted name's chemical words that the candidate also uses.
% Crude on purpose: it is a ranking aid for a choice that is flagged in the
% workbook anyway, not a structural similarity measure.
f=0;
if wanted=="" , return; end
tok=@(s) unique(string(regexp(lower(char(s)),'[a-z]{3,}','match')));
stop=["acid","the","and","alpha","beta","gamma","delta"];
a=setdiff(tok(candidate),stop); b=setdiff(tok(wanted),stop);
if isempty(b), return; end
f=numel(intersect(a,b))/numel(b);
end

function k=currencyKeys()
% Cofactors and carriers that never act as the signal-carrying substrate.
k=["h2o","water","h+","proton","h","sodium","na+","na1","atp","adp","amp","pi","phosphate", ...
   "ppi","diphosphate","nad","nadh","nadp","nadph","nadp(+)","nad(+)", ...
   "coa","coenzymea","co2","carbondioxide","o2","nh4","ammonium"];
end

function A=defaultAliases()
% Model metabolite name -> the BRENDA ligand name most often reported.
% Extend as needed; caller-supplied SubstrateAliases override these.
A=[ "THF","tetrahydrofolate"; ...
    "dihydrofolic acid","7,8-dihydrofolate"; ...
    "7,8-dihydropteroate","dihydropteroate"; ...
    "3-phospho-serine","O-phospho-L-serine"; ...
    "3-phospho-hydroxypyruvate","3-phosphohydroxypyruvate"; ...
    "glyceraldehyde 3-phosphate","D-glyceraldehyde 3-phosphate"; ...
    "keto-phenylpyruvate","phenylpyruvate"; ...
    "3-(4-hydroxyphenyl)pyruvate","4-hydroxyphenylpyruvate"; ...
    "xanthosine-5-phosphate","XMP"; ...
    "1-C-(indol-3-yl)glycerol 3-phosphate","indole-3-glycerol phosphate"; ...
    "acetyl-CoA","acetyl-CoA"; ...
    "coenzyme A","CoA"; ...
    "2-oxoglutarate","2-oxoglutarate"; ...
    "2-oxobutanoate","2-oxobutyrate"; ...
    "bicarbonate","hydrogencarbonate"; ...
    "L-cystathionine","cystathionine"; ...
    "tRNA(Phe)","tRNAPhe"; ...
    "tRNA(Ala)","tRNAAla"; ...
    "tRNA(Cys)","tRNACys"; ...
    "L-Arginosuccinic Acid","argininosuccinate"; ...
    "L-Argininosuccinic Acid","argininosuccinate"; ...
    "D-Glucose 1-Phosphate","alpha-D-glucose 1-phosphate"; ...
    "D-Fructose 6-Phosphate","D-fructose 6-phosphate"; ...
    "D-Glucosamine-6-Phosphate","D-glucosamine 6-phosphate"; ...
    "Uridine-5'-Diphosphate-D-Galactose","UDP-alpha-D-galactose"; ...
    "L-Aspartate","L-aspartate"];
end
function T=recoverLocalEC(T,opt,aliases)
% Stream the flat file once: retain only metadata and candidates for this
% selected network, rather than loading the full kinetics database in memory.
n=height(T); T.OriginalEC=T.EC;
T.ECSource=repmat("TXT",n,1); T.LocalECCandidates=strings(n,1);
T.LocalBrendaFile=strings(n,1);
T.ECLookupNote=repmat("EC annotation retained from TXT; not verified locally.",n,1);
file="";
if opt.UseLocalEC
    file=string(which(char(opt.LocalBrendaFile)));
    if file=="" && isfile(opt.LocalBrendaFile), file=string(opt.LocalBrendaFile); end
    if file=="" || ~isfile(file)
        warning('brendaKineticsTable:LocalBrendaMissing', ...
            'Local BRENDA file "%s" not found. Add its folder to the MATLAB path or set LocalBrendaFile. TXT EC annotations retained.', ...
            char(opt.LocalBrendaFile));
    else
        T.LocalBrendaFile(:)=file;
        fprintf('Local EC lookup: scanning %s\n',char(file));
        nameKeys=localEnzymeKey(T.Name);
        organismKey=localEnzymeKey(string(opt.Organism));
        % A model name may explicitly append the requested organism name.
        if organismKey~=""
            for i=1:n
                suffix=" "+organismKey;
                if endsWith(nameKeys(i),suffix)
                    nameKeys(i)=extractBefore(nameKeys(i),strlength(nameKeys(i))-strlength(suffix)+1);
                end
            end
        end
        equationKeys=strings(n,2);
        for i=1:n
            [equationKeys(i,1),equationKeys(i,2)]=localEquationKeys(T.NamesEquation(i),aliases);
        end
        hits=repmat(struct('EC',strings(0,1),'Kind',strings(0,1)),n,1);
        fid=fopen(char(file),'r','n','UTF-8');
        assert(fid>=0,'Cannot open local BRENDA file.'); guard=onCleanup(@()fclose(fid));
        entry=localEmptyEntry(); field=""; payload=""; scanned=0;
        while true
            line=fgetl(fid);
            if ~ischar(line), break; end
            line=erase(line,char(65279)); % optional UTF-8 BOM
            if startsWith(line,sprintf('\t'))
                if field~="", payload=payload+" "+strtrim(string(line)); end
                continue
            end
            entry=localStoreField(entry,field,payload); field=""; payload="";
            id=regexp(line,'^ID\s+([^\s]+)','tokens','once');
            if strcmp(strtrim(line),'///') || ~isempty(id)
                if entry.EC~=""
                    hits=localMatchEntry(hits,entry,nameKeys,equationKeys,aliases);
                    scanned=scanned+1;
                    if mod(scanned,1000)==0, fprintf('  Scanned %d EC entries.\n',scanned); end
                end
                entry=localEmptyEntry();
                if ~isempty(id) && ~isempty(regexp(id{1},'^\d+\.\d+\.\d+\.\d+$','once'))
                    entry.EC=string(id{1});
                end
                continue
            end
            if entry.EC=="", continue; end % e.g. ID spontaneous
            f=regexp(line,'^(RN|SN|SY|RE)\s+(.*)$','tokens','once');
            if ~isempty(f), field=string(f{1}); payload=string(f{2}); end
        end
        entry=localStoreField(entry,field,payload);
        if entry.EC~=""
            hits=localMatchEntry(hits,entry,nameKeys,equationKeys,aliases);
            scanned=scanned+1;
        end
        clear guard
        for i=1:n
            ec=hits(i).EC; kind=hits(i).Kind;
            T.LocalECCandidates(i)=strjoin(unique(ec,'stable'),'; ');
            strong=unique(ec(startsWith(kind,"equation")),'stable');
            named=unique(ec(kind=="exact enzyme name/synonym"),'stable');
            chosen=""; basis="";
            if numel(strong)==1
                chosen=strong; basis="local BRENDA equation match (water/protons omitted; either direction)";
            elseif isempty(strong) && numel(named)==1
                chosen=named; basis="local BRENDA exact name/synonym; UNCONFIRMED BY EQUATION";
            end
            if isfield(opt.PreferredEC,char(T.Reaction(i)))
                T.ECLookupNote(i)="Caller PreferredEC takes priority; local candidates: "+T.LocalECCandidates(i)+".";
            elseif chosen~=""
                T.EC(i)=chosen; T.ECSource(i)=basis;
                old=T.OriginalEC(i); if old=="", old="(missing)"; end
                T.ECLookupNote(i)="EC lookup: "+old+" -> "+chosen+"; "+basis+". ";
                if T.OriginalEC(i)~="" && T.OriginalEC(i)~=chosen
                    T.ECLookupNote(i)=T.ECLookupNote(i)+"CONFLICT with original TXT annotation; inspect Reaction mapping. ";
                end
            elseif ~isempty(ec)
                T.ECLookupNote(i)="Ambiguous local EC candidates: "+T.LocalECCandidates(i)+ ...
                    ". No automatic assignment; original TXT EC retained. Use ECOverrides after review.";
            else
                T.ECLookupNote(i)="No exact name/synonym or equation match in local BRENDA file; TXT EC retained. " + ...
                    "An absent match does not establish that the reaction lacks an EC classification.";
            end
        end
        fprintf('  Scanned %d EC entries; local evidence selected for %d/%d reactions.\n', ...
            scanned,sum(T.ECSource~="TXT"),n);
    end
end
% Model IDs become useful here as explicit mapping keys, never as BRENDA IDs.
overrides=string(opt.ECOverrides);
assert(size(overrides,2)==2,'ECOverrides must have two columns: model ID/row label, EC.');
assert(numel(unique(overrides(:,1)))==size(overrides,1),'Duplicate ECOverrides keys.');
for j=1:size(overrides,1)
    match=T.ID==overrides(j,1) | T.Reaction==overrides(j,1);
    assert(any(match),'ECOverrides key %s matches no model ID or row label.',char(overrides(j,1)));
    codes=regexp(char(overrides(j,2)),'\d+\.\d+\.\d+\.\d+','match');
    remainder=regexprep(char(overrides(j,2)),'\d+\.\d+\.\d+\.\d+','');
    assert(~isempty(codes) && isempty(regexprep(remainder,'[\s;,/]','')), ...
        'ECOverrides must contain complete EC numbers separated by semicolons, commas or slashes.');
    T.EC(match)=strjoin(string(codes),'; '); T.ECSource(match)="caller ECOverrides";
    T.ECLookupNote(match)="Explicit ECOverrides for "+overrides(j,1)+"; original TXT EC remains in Reaction mapping.";
end
end

function E=localEmptyEntry()
E=struct('EC',"",'Names',strings(0,1),'Equations',strings(0,1));
end

function E=localStoreField(E,field,payload)
if E.EC=="" || field=="" || payload=="", return; end
% Protein-specific commentary and literature tags also occur after RE fields.
payload=regexprep(payload,'^\s*#[^#]*#\s*','');
payload=regexprep(payload,'\s+\(#[\s\S]*$','');
payload=regexprep(payload,'\s*<[^>]*>','');
payload=regexprep(payload,'\s*\{[^}]*\}','');
payload=strtrim(payload);
if field=="RE"
    E.Equations(end+1,1)=payload;
else
    % Keep parentheses belonging to names, e.g. dehydrogenase (NAD).
    if payload~="", E.Names(end+1,1)=payload; end
end
end

function hits=localMatchEntry(hits,E,nameKeys,equationKeys,aliases)
matchedName=ismember(nameKeys,localEnzymeKey(E.Names)) & nameKeys~="";
matchedEquation=false(size(matchedName));
for j=1:numel(E.Equations)
    [forward,reverse]=localEquationKeys(E.Equations(j),aliases);
    if forward=="", continue; end
    matchedEquation=matchedEquation | (equationKeys(:,1)~="" & ...
        (equationKeys(:,1)==forward | equationKeys(:,1)==reverse));
end
for i=find(matchedName | matchedEquation)'
    hits(i).EC(end+1,1)=E.EC;
    if matchedEquation(i), hits(i).Kind(end+1,1)="equation";
    else, hits(i).Kind(end+1,1)="exact enzyme name/synonym"; end
end
end

function key=localEnzymeKey(name)
% Punctuation/case differences only; no substring or fuzzy-name assignment.
key=lower(strtrim(string(name))); key=replace(key,["γ","Γ"],"gamma");
key=regexprep(key,'[^a-z0-9]+',' '); key=strtrim(regexprep(key,'\s+',' '));
end

function [forward,reverse]=localEquationKeys(equation,aliases)
forward=""; reverse="";
equation=regexprep(char(equation),'\(m_[A-Za-z0-9_]+\)','');
parts=regexp(equation,'\s*(?:<=>|<->|->|=)\s*','split');
if numel(parts)~=2, return; end
left=localSideKey(parts{1},aliases); right=localSideKey(parts{2},aliases);
if left=="" || right=="" || left==right, return; end
forward=left+" => "+right; reverse=right+" => "+left;
end

function key=localSideKey(side,aliases)
key=""; parts=regexp(char(side),'\s+\+\s+','split');
counts=containers.Map('KeyType','char','ValueType','double');
for j=1:numel(parts)
    term=strtrim(parts{j}); stoich=1;
    coefficient=regexp(term,'^(\d+(?:\.\d+)?)\s+(.+)$','tokens','once');
    if ~isempty(coefficient), stoich=str2double(coefficient{1}); term=coefficient{2}; end
    name=normalize(aliasName(string(term),aliases));
    pairs=["water","h2o";"proton","h+";"sodium","na+"; ...
        "orthophosphate","phosphate";"pi","phosphate"; ...
        "carbon dioxide","co2";"diphosphate","ppi"; ...
        "3,5-cyclic adenosine monophosphate","camp"; ...
        "3',5'-camp","camp";"3',5'-cyclic amp","camp"; ...
        "3,5-cyclic guanosine monophosphate","cgmp"; ...
        "3',5'-cgmp","cgmp";"3',5'-cyclic gmp","cgmp"];
    hit=find(pairs(:,1)==name,1); if ~isempty(hit), name=pairs(hit,2); end
    if ismember(name,["h2o","h+"]), continue; end
    if ~isempty(regexp(char(name), ...
            '^(?:deoxy)?(?:adenosine|guanosine|cytidine|uridine|thymidine|inosine|xanthosine) (?:mono|di|tri)phosphate$','once'))
        abbreviation=nucleotideAbbrev(name);
        if abbreviation~="", name=lower(abbreviation); end
    end
    if name=="", return; end
    k=char(name);
    if isKey(counts,k), counts(k)=counts(k)+stoich; else, counts(k)=stoich; end
end
names=sort(string(keys(counts))); pieces=strings(size(names));
for j=1:numel(names), pieces(j)=sprintf('%.12g*%s',counts(char(names(j))),char(names(j))); end
key=strjoin(pieces,' + ');
end

function [email,password]=readCredentials(file)
file=resolveInputFile(file,'Credentials');
lines=regexp(fileread(file),'\r\n|\n|\r','split');
% Only strip email whitespace; spaces inside a password must be preserved.
assert(numel(lines)>=2 && ~isempty(strtrim(lines{1})) && ~isempty(lines{2}), ...
    'Credentials file must contain email on line 1 and password on line 2.');
email=strtrim(lines{1}); password=lines{2};
end
function file=resolveInputFile(requested,label)
% Respect MATLAB's current-folder/path search order after addpath(genpath).
% which accepts non-MATLAB filenames when the extension is supplied.
assert((ischar(requested) && isrow(requested)) || ...
    (isstring(requested) && isscalar(requested) && ~ismissing(requested)), ...
    'brendaKineticsTable:InvalidInputFile','%s filename must be a text scalar.',label);
requested=char(requested);
file=which(requested);
if isempty(file) && isfile(requested)
    % Fallback for an existing explicit path that which did not resolve.
    file=requested;
end
if isempty(file) || ~isfile(file)
    error('brendaKineticsTable:InputFileNotFound', ...
        ['%s file "%s" was not found. Add its folder to the MATLAB path ' ...
         '(for example, addingPathParentFolderByName(''code'')), ' ...
         'or supply its full path.'],label,requested);
end
end
function T=emptyRecords(parameter,unit)
T=table(strings(0,1),strings(0,1),strings(0,1),nan(0,1),nan(0,1), ...
    strings(0,1),strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    repmat(string(parameter),0,1),repmat(string(unit),0,1), ...
    'VariableNames',{'EC','organism','substrate','value','valueMaximum', ...
    'rawValue','rawValueMaximum','commentary','literature','ligandStructureId','parameter','unit'});
end
function [T,status]=fetchBothPools(ec,opt,email,password,parameter)
% Fetch the requested organism AND all organisms, and keep both. A query
% restricted to one organism can return plenty of records yet none for the
% ligand this reaction needs (EC 2.7.4.6 in human has ADP, ATP, GTP and dTDP
% but no dGDP), so retrying only when the first query is EMPTY is not enough.
% selectRecords still ranks the requested organism first, so nothing changes
% for reactions that already matched.
[T,status]=fetchRecords(ec,char(opt.Organism),email,password,parameter);
if ~opt.OrganismFallback || isempty(char(opt.Organism)), return; end
if ~startsWith(status,'OK') && ~startsWith(status,'No records'), return; end  % real failure
[T2,status2]=fetchRecords(ec,'',email,password,parameter);
if height(T2)>0
    % Drop duplicates of the organism-specific rows before pooling.
    key  = @(X) string(X.organism)+"|"+string(X.substrate)+"|"+string(X.rawValue);
    T2   = T2(~ismember(key(T2),key(T)),:);
    T    = [T;T2];
    status=[status '; ALL ORGANISMS: ' status2];
end
end

function [T,status]=fetchRecords(ec,organism,email,password,parameter)
hasSubstrate=true;
switch parameter
    case 'Km'
        operation='getKmValue'; field='kmValue'; maximum='kmValueMaximum'; unit='mM';
    case 'kcat'
        operation='getTurnoverNumber'; field='turnoverNumber'; maximum='turnoverNumberMaximum'; unit='1/s';
    case 'SA'
        operation='getSpecificActivity'; field='specificActivity';
        maximum='specificActivityMaximum'; unit='umol/min/mg'; hasSubstrate=false;
    case 'MW'
        operation='getMolecularWeight'; field='molecularWeight';
        maximum='molecularWeightMaximum'; unit='Da'; hasSubstrate=false;
    otherwise
        error('Unknown BRENDA parameter %s.',parameter);
end
T=emptyRecords(parameter,unit);
try
    % Zeep's dynamic classes/proxies must stay inside Python. MATLAB only
    % receives a built-in str containing JSON, never a SOAP object or callable.
    persistent bridge
    if isempty(bridge), bridge=py.dict; end
    bridge.update(py.dict(pyargs('email',email,'password',password, ...
        'ec',ec,'organism',organism,'operation',operation, ...
        'field',field,'maximum',maximum,'has_substrate',hasSubstrate)));
    clear email password
    code={ ...
        'import json'; ...
        'phase = "import Python dependencies"'; ...
        'try:'; ...
        '    import hashlib, time, zeep'; ...
        '    from zeep.transports import Transport'; ...
        '    from zeep.helpers import serialize_object'; ...
        '    phase = "load BRENDA WSDL"'; ...
        '    if "client" not in globals():'; ...
        '        client = zeep.Client("https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl", settings=zeep.Settings(strict=False), transport=Transport(timeout=60, operation_timeout=120))'; ...
        '    phase = "resolve SOAP operation"'; ...
        '    soap_call = getattr(client.service, operation)'; ...
        '    password_hash = hashlib.sha256(password.encode("utf-8")).hexdigest()'; ...
        '    phase = "execute SOAP request"'; ...
        '    time.sleep(1.1)'; ...
        '    if has_substrate:'; ...
        '        args = [email, password_hash, "ecNumber*" + ec, field + "*", maximum + "*", "substrate*", "commentary*", "organism*" + organism, "ligandStructureId*", "literature*"]'; ...
        '    else:'; ...
        '        args = [email, password_hash, "ecNumber*" + ec, field + "*", maximum + "*", "commentary*", "organism*" + organism, "literature*"]'; ...
        '    raw = soap_call(*args)'; ...
        '    phase = "serialize SOAP response"'; ...
        '    response_json = json.dumps({"ok": True, "records": serialize_object(raw)})'; ...
        'except Exception as exc:'; ...
        '    response_json = json.dumps({"ok": False, "phase": phase, "error_type": type(exc).__module__ + "." + type(exc).__name__})'; ...
        'finally:'; ...
        '    for key in ("email", "password", "password_hash", "raw", "soap_call", "args"):'; ...
        '        globals().pop(key, None)'};
    py.exec(strjoin(code,newline),bridge);
    response=jsondecode(char(bridge.get('response_json')));
    if ~response.ok
        status=sprintf('Query failed during %s (%s).',response.phase,response.error_type);
        warning('brendaKineticsTable:QueryFailed','%s / %s: %s',ec,parameter,status);
        return
    end
    records=response.records;
    if isempty(records), status='No records'; return; end
    assert(isstruct(records)||iscell(records),'Unexpected SOAP result type.');
    for i=1:numel(records)
        if iscell(records), r=records{i}; else, r=records(i); end
        assert(isstruct(r),'Unexpected SOAP record type.');
        rv=textField(r,field); rm=textField(r,maximum);
        T(end+1,:)={textField(r,'ecNumber'),textField(r,'organism'), ...
            textField(r,'substrate'),str2double(rv),str2double(rm),rv,rm, ...
            textField(r,'commentary'),textField(r,'literature'), ...
            textField(r,'ligandStructureId'),string(parameter),string(unit)}; %#ok<AGROW>
        if T.EC(end)=="", T.EC(end)=string(ec); end
    end
    T.value(~isfinite(T.value)|T.value<=0)=NaN;
    T.valueMaximum(~isfinite(T.valueMaximum)|T.valueMaximum<=0)=NaN;
    status=sprintf('OK: %d records',height(T));
catch ME
    T=emptyRecords(parameter,unit);
    % Never expose remote exception messages that might include request secrets.
    status=['Query failed (' ME.identifier '); check local Python/zeep and BRENDA access.'];
    warning('brendaKineticsTable:QueryFailed','%s / %s: %s',ec,parameter,status);
end
end
function [newEC,verdict,trail,acceptedName]=resolveTransferredEC(ec)
% Ask ExPASy whether an EC number is still current, and follow it if not.
%
% Why this is needed. BRENDA files its data under the CURRENT EC number only.
% When an entry is renumbered the old number keeps appearing in genome-scale
% models (BiGG, Recon, the CAR-T model used here) for years, and every query
% under it returns nothing: no Km, no kcat, no SA, no MW, with no error to
% distinguish "obsolete code" from "enzyme never characterised". EC 2.3.2.4 is
% exactly this case -- it was transferred to EC 4.3.2.9 (gamma-glutamyl-
% cyclotransferase), which holds 110 records.
%
% The ExPASy flat file states the verdict in one line:
%       ID   2.3.2.4
%       DE   Transferred entry: 4.3.2.9
% A transfer can chain, so it is followed transitively with a hop cap, and a
% transfer that fans out to several codes is NOT followed: picking one of them
% is a scientific decision, not a lookup. Results are cached per session.
persistent cache
newEC=''; verdict=''; trail=''; acceptedName='';
if isempty(cache), cache=containers.Map('KeyType','char','ValueType','any'); end
if isKey(cache,ec)
    c=cache(ec); newEC=c{1}; verdict=c{2}; trail=c{3}; acceptedName=c{4}; return
end
chain={ec}; current=ec; options=weboptions('Timeout',25,'ContentType','text');
for hop=1:5
    try
        page=webread(sprintf('https://enzyme.expasy.org/EC/%s.txt',current),options);
    catch ME
        verdict=sprintf('ExPASy unreachable for EC %s (%s); code left unchanged', ...
            current,ME.identifier);
        trail=strjoin(chain,' -> ');
        cache(ec)={'',verdict,trail,''};
        return
    end
    page=char(page);
    % DE lines hold either the accepted name or the transfer/deletion notice.
    % lineanchors makes ^/$ per line; dotexceptnewline stops '.' from running
    % past the DE line (MATLAB's regexp is dot-all by default).
    de=regexp(page,'^DE\s+(.*)$','tokens','lineanchors','dotexceptnewline');
    de=strtrim(strjoin(cellfun(@(t)strtrim(t{1}),de,'UniformOutput',false),' '));
    % The notice ends in a full stop ("Transferred entry: 4.3.2.9."), so the
    % lookahead must reject a following DIGIT only, not a following period.
    targets=unique(regexp(de,'(?<![\d.])\d+\.\d+\.\d+\.\d+(?!\d)','match'),'stable');
    if contains(lower(de),'deleted entry')
        verdict=sprintf('ExPASy: EC %s is a DELETED entry; there is no successor and no data exists under it',current);
        break
    elseif contains(lower(de),'transferred entry') && numel(targets)==1
        chain{end+1}=targets{1}; %#ok<AGROW>
        current=targets{1};
        continue
    elseif contains(lower(de),'transferred entry') && numel(targets)>1
        verdict=sprintf(['ExPASy: EC %s was SPLIT into %s. Choosing between successors is a ' ...
            'biochemical decision, not a lookup: set the right one in the motif file.'], ...
            current,strjoin(targets,', '));
        break
    else
        acceptedName=de;
        if strcmp(current,ec)
            verdict=sprintf('ExPASy: EC %s is current ("%s"), so BRENDA simply holds no measurement for it',current,de);
        else
            verdict=sprintf('ExPASy: transferred to EC %s ("%s")',current,de);
            newEC=current;
        end
        break
    end
end
if isempty(verdict)
    verdict=sprintf('ExPASy: transfer chain from EC %s did not terminate within 5 hops',ec);
    newEC='';
end
trail=strjoin(chain,' -> ');
cache(ec)={newEC,verdict,trail,acceptedName};
end

function value=textField(record,field)
value="";
if ~isfield(record,field)||isempty(record.(field)), return; end
v=record.(field);
if ischar(v)||(isstring(v)&&isscalar(v)), value=string(v);
else, value=string(jsonencode(v)); end
end
function deleteIfExists(file)
if isfile(file), delete(file); end
end

function installExcelFormulas(file,values)
% Explicit OpenXML formula nodes: writecell can export formula strings as
% text when Excel is unavailable (e.g. on Linux/HPC). No Excel COM required.
root=tempname; mkdir(root); guard=onCleanup(@()rmdir(root,'s'));
unzip(file,root);
sheetFile=fullfile(root,'xl','worksheets','sheet1.xml');
doc=xmlread(sheetFile); cells=doc.getElementsByTagName('c');
ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main';
for i=0:cells.getLength-1
    cellNode=cells.item(i); address=char(cellNode.getAttribute('r'));
    token=regexp(address,'^([A-Z]+)(\d+)$','tokens','once');
    if isempty(token), continue; end
    col=0;
    for ch=token{1}, col=26*col+double(ch)-double('A')+1; end
    row=str2double(token{2});
    if row>size(values,1)||col>size(values,2), continue; end
    v=values{row,col};
    if ~(ischar(v)&&~isempty(v)&&v(1)=='='), continue; end
    % Remove shared/inline text and any stale cache, then install the formula.
    while cellNode.hasChildNodes, cellNode.removeChild(cellNode.getFirstChild); end
    cellNode.removeAttribute('t');
    f=doc.createElementNS(ns,'f'); f.appendChild(doc.createTextNode(v(2:end)));
    cellNode.appendChild(f);
    % Cache the two explicit SA/MW estimates for readcell on Linux/HPC.
    if col==8 && isnumeric(values{row,6}) && isscalar(values{row,6}) && ...
            isnumeric(values{row,7}) && isscalar(values{row,7})
        cached=values{row,6}*values{row,7}/60;
        if isfinite(cached)
            node=doc.createElementNS(ns,'v');
            node.appendChild(doc.createTextNode(sprintf('%.17g',cached)));
            cellNode.appendChild(node);
        end
    end
end
xmlwrite(sheetFile,doc);
workbookFile=fullfile(root,'xl','workbook.xml');
doc=xmlread(workbookFile); nodes=doc.getElementsByTagName('calcPr');
if nodes.getLength==0
    calc=doc.createElementNS(ns,'calcPr'); doc.getDocumentElement.appendChild(calc);
else, calc=nodes.item(0); end
calc.setAttribute('calcMode','auto'); calc.setAttribute('fullCalcOnLoad','1');
calc.setAttribute('forceFullCalc','1'); xmlwrite(workbookFile,doc);
archive=[tempname '.zip']; archiveGuard=onCleanup(@()deleteIfExists(archive));
entries=dir(root); entries=entries(~ismember({entries.name},{'.','..'}));
zip(archive,{entries.name},root);
[ok,msg]=movefile(archive,file,'f'); assert(ok,'%s',msg);
end

function records=referenceSnapshot()
% Frozen, cited BRENDA selections approved for exploratory simulation.
% Numeric values are never represented as fresh query results.
payload={ ...
    '[{"reaction":"r_adapter_1","id":"r_0033","ec":"3.1.4.17","model_ec":"3.1.4.17; 3.1.4.53","equation":"3'',5''-cyclic AMP + H2O  ->  AMP + H+","direction":1,"km":{"value":0.00285,"substrate":"3'',5''-cAMP","organism":"Homo sapiens","literature":"135219; 653640; 667684; 668133; 668990; 669403; 678221; 679784","conditions":"Median of 14 positive records; multiple assays/variants may be included.","kind":"median","n":14,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.4.17"},"kcat":{"value":0.49,"substrate":"3'',5''-cAMP","organism":"Homo sapiens","literature":"667684; 669403","conditions":"Median of 4 positive records; multiple assays/variants may be included.","kind":"median","n":4,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.4.17"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 3.1.4.17. K_M=0.00285 mM for 3'',5''-cAMP; Homo sapiens; BRENDA literature 135219; 653640; 667684; 668133; 668990; 669403; 678221; 679784; Median of 14 positive records; multiple assays/variants may be included. k_cat=0.49 s^-1 for 3'',5''-cAMP; Homo sapiens; BRENDA literature 667684; 669403; Median of 4 positive records; multiple assays/variants may be included. Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized."},'; ...
    '{"reaction":"r_adapter_2","id":"r_0208","ec":"6.3.4.5","model_ec":"6.3.4.5","equation":"(N(omega)-L-arginino)succinic acid + AMP + diphosphate + H+  ->  ATP + L-aspartate + L-citrulline","direction":-1,"km":{"value":0.051,"substrate":"ATP","organism":"Homo sapiens","literature":"1598","conditions":"wild-type enzyme","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=6.3.4.5"},"kcat":{"value":3.01,"substrate":"ATP + citrulline + aspartate","organism":"Homo sapiens","literature":"1581","conditions":"Estimated from purified-enzyme SA=4.2 umol/min/mg and subunit MW=43 mg/umol: SA*MW/60; assumes one active site per subunit.","kind":"derived SA*MW/60","n":1,"source":"https://www.brenda-enzymes.org/literature.php?e=6.3.4.5&r=1581"},"sa":4.2,"mw":43,"note":"PROVISIONAL REFERENCE; EC selected: 6.3.4.5. K_M=0.051 mM for ATP; Homo sapiens; BRENDA literature 1598; wild-type enzyme k_cat=3.01 s^-1 for ATP + citrulline + aspartate; Homo sapiens; BRENDA literature 1581; Estimated from purified-enzyme SA=4.2 umol/min/mg and subunit MW=43 mg/umol: SA*MW/60; assumes one active site per subunit. BRENDA''s reference reaction runs in the opposite direction, consuming ATP, aspartate and citrulline. Model consumes argininosuccinate, AMP and diphosphate. Forward constants used as reverse-model assumptions; not measured AMP/reverse kinetics. Km and kcat come from different studies. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"BRENDA''s reference reaction runs in the opposite direction, consuming ATP, aspartate and citrulline. Model consumes argininosuccinate, AMP and diphosphate. Forward constants used as reverse-model assumptions; not measured AMP/reverse kinetics. Km and kcat come from different studies."},'; ...
    '{"reaction":"r_adapter_3","id":"r_0211","ec":"6.3.5.4","model_ec":"6.3.5.4","equation":"ATP + H2O + L-aspartate + L-glutamine  ->  AMP + diphosphate + H+ + L-asparagine + L-glutamate(m_d1)","direction":1,"km":{"value":1.9,"substrate":"L-glutamine","organism":"Homo sapiens","literature":"660816","conditions":"pH 8, C-terminally tagged recombinant enzyme","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=6.3.5.4"},"kcat":{"value":1.7,"substrate":"glutamine","organism":"Homo sapiens","literature":"660816","conditions":"pH 8, C-terminally tagged recombinant enzyme","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=6.3.5.4"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 6.3.5.4. K_M=1.9 mM for L-glutamine; Homo sapiens; BRENDA literature 660816; pH 8, C-terminally tagged recombinant enzyme k_cat=1.7 s^-1 for glutamine; Homo sapiens; BRENDA literature 660816; pH 8, C-terminally tagged recombinant enzyme Glutamine selected from multiple reactants. ATP and aspartate are not separately parameterized. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Glutamine selected from multiple reactants. ATP and aspartate are not separately parameterized."},'; ...
    '{"reaction":"r_adapter_4","id":"r_0918","ec":"2.6.1.52","model_ec":"2.6.1.52","equation":"3-phospho-hydroxypyruvate + L-glutamate(m_d1)  ->  2-oxoglutarate + 3-phospho-serine","direction":1,"km":{"value":1.03,"substrate":"L-glutamate","organism":"Synechocystis sp.","literature":"739017","conditions":"pH 9.0; 30 C.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/literature.php?e=2.6.1.52&r=739017"},"kcat":{"value":2.49,"substrate":"L-glutamate","organism":"Synechocystis sp.","literature":"739017","conditions":"pH 9.0; 30 C.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/literature.php?e=2.6.1.52&r=739017"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 2.6.1.52. K_M=1.03 mM for L-glutamate; Synechocystis sp.; BRENDA literature 739017; pH 9.0; 30 C. k_cat=2.49 s^-1 for L-glutamate; Synechocystis sp.; BRENDA literature 739017; pH 9.0; 30 C. Non-human reference pair; assay conditions may differ from the simulation. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Non-human reference pair; assay conditions may differ from the simulation."},'; ...
    '{"reaction":"r_adapter_5","id":"r_0917","ec":"3.1.3.3","model_ec":"3.1.3.3","equation":"3-phospho-serine + H2O  ->  L-serine(m_in) + phosphate","direction":1,"km":{"value":0.036,"substrate":"L-phosphoserine","organism":"Homo sapiens","literature":"81054","conditions":"Assay conditions not specified in this record.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.3.3"},"kcat":{"value":1.508,"substrate":"L-phosphoserine","organism":"Porphyromonas gingivalis","literature":"682506","conditions":"pH 7.5; 30 C; 5 mM MgCl2; one reported scalar selected.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/literature.php?e=3.1.3.3&r=682506"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 3.1.3.3. K_M=0.036 mM for L-phosphoserine; Homo sapiens; BRENDA literature 81054; Assay conditions not specified in this record. k_cat=1.508 s^-1 for L-phosphoserine; Porphyromonas gingivalis; BRENDA literature 682506; pH 7.5; 30 C; 5 mM MgCl2; one reported scalar selected. Model 3-phospho-serine interpreted as L-phosphoserine. Km and kcat are from DIFFERENT organisms/studies; unpaired modelling reference. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Model 3-phospho-serine interpreted as L-phosphoserine. Km and kcat are from DIFFERENT organisms/studies; unpaired modelling reference."},'; ...
    '{"reaction":"r_Ff","id":"r_0309","ec":"4.2.1.22","model_ec":"4.2.1.22","equation":"L-homocysteine + L-serine(m_in)  ->  H2O + L-cystathionine(m_F)","direction":1,"km":{"value":2.13,"substrate":"L-serine","organism":"Homo sapiens","literature":"649869; 650017; 663623; 665593; 677837; 680084; 698182; 701064; 746768; 747922; 770220","conditions":"Median of 34 positive records; multiple assays/variants may be included.","kind":"median","n":34,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=4.2.1.22"},"kcat":{"value":5.4,"substrate":"L-serine","organism":"Homo sapiens","literature":"649869; 650017; 663623; 665593; 677837; 680084; 701064; 746768; 770220","conditions":"Median of 27 positive records; multiple assays/variants may be included.","kind":"median","n":27,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=4.2.1.22"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 4.2.1.22. K_M=2.13 mM for L-serine; Homo sapiens; BRENDA literature 649869; 650017; 663623; 665593; 677837; 680084; 698182; 701064; 746768; 747922; 770220; Median of 34 positive records; multiple assays/variants may be included. k_cat=5.4 s^-1 for L-serine; Homo sapiens; BRENDA literature 649869; 650017; 663623; 665593; 677837; 680084; 701064; 746768; 770220; Median of 27 positive records; multiple assays/variants may be included. Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized."},'; ...
    '{"reaction":"r_Fd","id":"r_0310","ec":"4.4.1.1","model_ec":"4.4.1.1","equation":"H2O + L-cystathionine(m_F)  ->  2-oxobutanoate + ammonium + L-cysteine(m_P1)","direction":1,"km":{"value":0.55,"substrate":"L-cystathionine","organism":"Homo sapiens","literature":"34508; 690967; 693241","conditions":"Median of 6 positive records; multiple assays/variants may be included.","kind":"median","n":6,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=4.4.1.1"},"kcat":{"value":1.7,"substrate":"L-cystathionine","organism":"Homo sapiens","literature":"690967; 693241","conditions":"Median of 5 positive records; multiple assays/variants may be included.","kind":"median","n":5,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=4.4.1.1"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 4.4.1.1. K_M=0.55 mM for L-cystathionine; Homo sapiens; BRENDA literature 34508; 690967; 693241; Median of 6 positive records; multiple assays/variants may be included. k_cat=1.7 s^-1 for L-cystathionine; Homo sapiens; BRENDA literature 690967; 693241; Median of 5 positive records; multiple assays/variants may be included. Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized."},'; ...
    '{"reaction":"r_Sf","id":"r_1055","ec":"4.2.1.20","model_ec":"4.2.1.20","equation":"1-C-(indol-3-yl)glycerol 3-phosphate + L-serine(m_in)  ->  glyceraldehyde 3-phosphate(m_S) + H2O + L-tryptophan","direction":1,"km":{"value":0.49,"substrate":"L-serine","organism":"Escherichia coli","literature":"771571","conditions":"Wild-type TrpA/TrpB complex; pH 7.5; temperature unspecified.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/literature.php?e=4.2.1.20&r=771571"},"kcat":{"value":2.83,"substrate":"L-serine","organism":"Escherichia coli","literature":"771571","conditions":"Wild-type TrpA/TrpB complex; pH 7.5; temperature unspecified.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/literature.php?e=4.2.1.20&r=771571"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 4.2.1.20. K_M=0.49 mM for L-serine; Escherichia coli; BRENDA literature 771571; Wild-type TrpA/TrpB complex; pH 7.5; temperature unspecified. k_cat=2.83 s^-1 for L-serine; Escherichia coli; BRENDA literature 771571; Wild-type TrpA/TrpB complex; pH 7.5; temperature unspecified. Non-human reference pair; indoleglycerol phosphate cosubstrate not separately parameterized. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Non-human reference pair; indoleglycerol phosphate cosubstrate not separately parameterized."},'; ...
    '{"reaction":"r_Sd","id":"r_4211","ec":"4.3.3.6","model_ec":"4.3.3.6; 3.5.1.2","equation":"glyceraldehyde 3-phosphate(m_S) + L-glutamine + aldehydo-D-ribose 5-phosphate  ->  H+ + 3 H2O + L-glutamate(m_d1) + phosphate + pyridoxal 5''-phosphate","direction":1,"km":{"value":1.06,"substrate":"D-glyceraldehyde 3-phosphate","organism":"Geobacillus kaustophilus","literature":"748203","conditions":"pH 8.5; 37 C.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=4.3.3.6"},"kcat":{"value":0.00017,"substrate":"D-glyceraldehyde 3-phosphate","organism":"Geobacillus kaustophilus","literature":"748203","conditions":"pH 8.5; 37 C.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=4.3.3.6"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 4.3.3.6. K_M=1.06 mM for D-glyceraldehyde 3-phosphate; Geobacillus kaustophilus; BRENDA literature 748203; pH 8.5; 37 C. k_cat=0.00017 s^-1 for D-glyceraldehyde 3-phosphate; Geobacillus kaustophilus; BRENDA literature 748203; pH 8.5; 37 C. EC 4.3.3.6 selected for PLP synthesis; alternative annotation 3.5.1.2 not used. Non-human pair; model glyceraldehyde phosphate interpreted as D isomer. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"EC 4.3.3.6 selected for PLP synthesis; alternative annotation 3.5.1.2 not used. Non-human pair; model glyceraldehyde phosphate interpreted as D isomer."},'; ...
    '{"reaction":"r_I","id":"r_0460","ec":"6.3.2.2","model_ec":"6.3.2.2","equation":"ATP + L-cysteine(m_P1) + L-glutamate(m_d1)  ->  ADP + H+ + L-gamma-glutamyl-L-cysteine + phosphate","direction":1,"km":{"value":0.07,"substrate":"L-Cys","organism":"Homo sapiens","literature":"714265","conditions":"wild-type, CGL holoenzyme, pH 8.0, 37?C","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=6.3.2.2"},"kcat":{"value":8.2,"substrate":"L-cysteine","organism":"Mus musculus","literature":"660892","conditions":"Holoenzyme; pH 7.2; 37 C.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/literature.php?e=6.3.2.2&r=660892"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 6.3.2.2. K_M=0.07 mM for L-Cys; Homo sapiens; BRENDA literature 714265; wild-type, CGL holoenzyme, pH 8.0, 37?C k_cat=8.2 s^-1 for L-cysteine; Mus musculus; BRENDA literature 660892; Holoenzyme; pH 7.2; 37 C. Cysteine reference chosen for cancellation. Km is human, kcat is mouse; unpaired modelling values. ATP/glutamate dependence requires additional parameters. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Cysteine reference chosen for cancellation. Km is human, kcat is mouse; unpaired modelling values. ATP/glutamate dependence requires additional parameters."},'; ...
    '{"reaction":"r_Fr","id":"r_0313","ec":"6.1.1.16","model_ec":"6.1.1.16","equation":"ATP + L-cysteine(m_P1) + tRNA(Cys)  ->  AMP + Cys-tRNA(Cys) + diphosphate","direction":1,"km":{"value":0.0153,"substrate":"L-cysteine","organism":"Homo sapiens","literature":"675435","conditions":"Median of 2 positive records; multiple assays/variants may be included.","kind":"median","n":2,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=6.1.1.16"},"kcat":{"value":14.200000000000001,"substrate":"L-cysteine","organism":"Homo sapiens","literature":"675435","conditions":"Median of 2 positive records; multiple assays/variants may be included.","kind":"median","n":2,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=6.1.1.16"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 6.1.1.16. K_M=0.0153 mM for L-cysteine; Homo sapiens; BRENDA literature 675435; Median of 2 positive records; multiple assays/variants may be included. k_cat=14.2 s^-1 for L-cysteine; Homo sapiens; BRENDA literature 675435; Median of 2 positive records; multiple assays/variants may be included. Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized."},'; ...
    '{"reaction":"r_Ff_2","id":"r_0851","ec":"2.6.1.39","model_ec":"2.6.1.39; 2.6.1.57","equation":"keto-phenylpyruvate + L-glutamate(m_d1)  ->  2-oxoglutarate + L-phenylalanine(m_F2)","direction":-1,"km":{"value":12.5,"substrate":"L-glutamate","organism":"Homo sapiens","literature":"639983","conditions":"pH 8.0, 37?C, isoenzyme AadAT-II","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=2.6.1.39"},"kcat":{"value":12.8,"substrate":"L-glutamate","organism":"Saccharomyces cerevisiae","literature":"721430","conditions":"pH 7.0; temperature unspecified.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=2.6.1.39"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 2.6.1.39. K_M=12.5 mM for L-glutamate; Homo sapiens; BRENDA literature 639983; pH 8.0, 37?C, isoenzyme AadAT-II k_cat=12.8 s^-1 for L-glutamate; Saccharomyces cerevisiae; BRENDA literature 721430; pH 7.0; temperature unspecified. Model Direction=-1. Glutamate records are not verified for the model phenylpyruvate reaction. Human Km and yeast kcat are an unpaired reference; EC 2.6.1.57 not used. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Model Direction=-1. Glutamate records are not verified for the model phenylpyruvate reaction. Human Km and yeast kcat are an unpaired reference; EC 2.6.1.57 not used."},'; ...
    '{"reaction":"r_Fd_2","id":"r_2117","ec":"2.6.1.7","model_ec":"2.6.1.58; 2.6.1.7","equation":"L-phenylalanine(m_F2) + pyruvate  ->  keto-phenylpyruvate + L-alanine(m_d2)","direction":1,"km":{"value":1.7,"substrate":"L-phenylalanine","organism":"Homo sapiens","literature":"658649","conditions":"Median of 1 positive records; multiple assays/variants may be included.","kind":"median","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=2.6.1.7"},"kcat":{"value":91,"substrate":"L-phenylalanine","organism":"Homo sapiens","literature":"658649","conditions":"Median of 1 positive records; multiple assays/variants may be included.","kind":"median","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=2.6.1.7"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 2.6.1.7. K_M=1.7 mM for L-phenylalanine; Homo sapiens; BRENDA literature 658649; Median of 1 positive records; multiple assays/variants may be included. k_cat=91 s^-1 for L-phenylalanine; Homo sapiens; BRENDA literature 658649; Median of 1 positive records; multiple assays/variants may be included. Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized."},'; ...
    '{"reaction":"r_Sf_2","id":"r_0989","ec":"1.5.1.10","model_ec":"1.5.1.10","equation":"H+ + L-allysine + L-glutamate(m_d1) + NADPH  ->  H2O + L-saccharopine + NADP(+)(m_S2)","direction":1,"km":{"value":2.32,"substrate":"L-saccharopine","organism":"Saccharomyces cerevisiae","literature":"437757","conditions":"Scalar saccharopine value.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=1.5.1.10"},"kcat":{"value":224.49999999999997,"substrate":"L-saccharopine","organism":"Saccharomyces cerevisiae","literature":"437757","conditions":"Estimated from purified-enzyme SA=269.4 umol/min/mg and monomer MW=50 mg/umol: SA*MW/60; assumes one active site per monomer.","kind":"derived SA*MW/60","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=1.5.1.10"},"sa":269.4,"mw":50,"note":"PROVISIONAL REFERENCE; EC selected: 1.5.1.10. K_M=2.32 mM for L-saccharopine; Saccharomyces cerevisiae; BRENDA literature 437757; Scalar saccharopine value. k_cat=224.5 s^-1 for L-saccharopine; Saccharomyces cerevisiae; BRENDA literature 437757; Estimated from purified-enzyme SA=269.4 umol/min/mg and monomer MW=50 mg/umol: SA*MW/60; assumes one active site per monomer. Model makes saccharopine; reference consumes it. Reverse/substrate surrogate despite model Direction=1. Reference does not provide glutamate or NADPH Km. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Model makes saccharopine; reference consumes it. Reverse/substrate surrogate despite model Direction=1. Reference does not provide glutamate or NADPH Km."},'; ...
    '{"reaction":"r_Sd_2","id":"r_0939","ec":"1.3.1.13","model_ec":"1.3.1.13","equation":"NADP(+)(m_S2) + prephenate  ->  3-(4-hydroxyphenyl)pyruvate(m_N2) + carbon dioxide + NADPH","direction":1,"km":{"value":0.38,"substrate":"prephenate","organism":"Methanothrix harundinacea","literature":"763088","conditions":"Wild-type; pH 7.6; temperature unspecified.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=1.3.1.13"},"kcat":{"value":2.44,"substrate":"prephenate","organism":"Methanothrix harundinacea","literature":"763088","conditions":"Wild-type; pH 7.6; temperature unspecified.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=1.3.1.13"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 1.3.1.13. K_M=0.38 mM for prephenate; Methanothrix harundinacea; BRENDA literature 763088; Wild-type; pH 7.6; temperature unspecified. k_cat=2.44 s^-1 for prephenate; Methanothrix harundinacea; BRENDA literature 763088; Wild-type; pH 7.6; temperature unspecified. Non-human prephenate reference; model tagged input is NADP(+), whose Km is not represented by this value. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Non-human prephenate reference; model tagged input is NADP(+), whose Km is not represented by this value."},'; ...
    '{"reaction":"r_I_2","id":"r_2119","ec":"2.6.1.7","model_ec":"2.6.1.57; 2.6.1.7","equation":"3-(4-hydroxyphenyl)pyruvate(m_N2) + L-alanine(m_d2)  ->  L-tyrosine + pyruvate","direction":1,"km":{"value":1.5,"substrate":"hydroxyphenylpyruvate","organism":"Homo sapiens","literature":"702779","conditions":"in 100 mM potassium phosphate buffer (pH 7.4), at 37?C","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=2.6.1.7"},"kcat":{"value":0.39,"substrate":"hydroxyphenylpyruvate","organism":"Homo sapiens","literature":"702779","conditions":"in 100 mM potassium phosphate buffer (pH 7.4), at 37?C","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=2.6.1.7"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 2.6.1.7. K_M=1.5 mM for hydroxyphenylpyruvate; Homo sapiens; BRENDA literature 702779; in 100 mM potassium phosphate buffer (pH 7.4), at 37?C k_cat=0.39 s^-1 for hydroxyphenylpyruvate; Homo sapiens; BRENDA literature 702779; in 100 mM potassium phosphate buffer (pH 7.4), at 37?C EC 2.6.1.7 selected. Hydroxyphenylpyruvate reference used for the model hydroxyphenylpyruvate; donor/cosubstrate equivalence to L-alanine not established. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"EC 2.6.1.7 selected. Hydroxyphenylpyruvate reference used for the model hydroxyphenylpyruvate; donor/cosubstrate equivalence to L-alanine not established."},'; ...
    '{"reaction":"r_Sr_2","id":"r_4669","ec":"4.1.1.1","model_ec":"4.1.1.1","equation":"3-(4-hydroxyphenyl)pyruvate(m_N2) + H+  ->  carbon dioxide + (4-hydroxyphenyl)acetaldehyde","direction":1,"km":{"value":0.68,"substrate":"pyruvate","organism":"Zymomonas mobilis","literature":"650009","conditions":"Wild-type PDC; 30 C.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/literature.php?e=4.1.1.1&r=650009"},"kcat":{"value":113,"substrate":"pyruvate","organism":"Zymomonas mobilis","literature":"650009","conditions":"Wild-type PDC; 30 C.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/literature.php?e=4.1.1.1&r=650009"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 4.1.1.1. K_M=0.68 mM for pyruvate; Zymomonas mobilis; BRENDA literature 650009; Wild-type PDC; 30 C. k_cat=113 s^-1 for pyruvate; Zymomonas mobilis; BRENDA literature 650009; Wild-type PDC; 30 C. Non-human, DIFFERENT substrate: pyruvate constants used as surrogates for 3-(4-hydroxyphenyl)pyruvate. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Non-human, DIFFERENT substrate: pyruvate constants used as surrogates for 3-(4-hydroxyphenyl)pyruvate."},'; ...
    '{"reaction":"r_c","id":"r_0156","ec":"2.6.1.44","model_ec":"2.6.1.44","equation":"glyoxylate + L-alanine(m_d2)  ->  L-glycine(m_c) + pyruvate","direction":1,"km":{"value":31.0,"substrate":"L-alanine","organism":"Homo sapiens","literature":"640079; 640081; 640087; 684989; 704634; 721718; 739610; 758860; 759779","conditions":"Median of 32 positive records; multiple assays/variants may be included.","kind":"median","n":32,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=2.6.1.44"},"kcat":{"value":25.5,"substrate":"L-alanine","organism":"Homo sapiens","literature":"684989; 704634; 721718; 739610; 758860; 759779","conditions":"Median of 28 positive records; multiple assays/variants may be included.","kind":"median","n":28,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=2.6.1.44"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 2.6.1.44. K_M=31 mM for L-alanine; Homo sapiens; BRENDA literature 640079; 640081; 640087; 684989; 704634; 721718; 739610; 758860; 759779; Median of 32 positive records; multiple assays/variants may be included. k_cat=25.5 s^-1 for L-alanine; Homo sapiens; BRENDA literature 684989; 704634; 721718; 739610; 758860; 759779; Median of 28 positive records; multiple assays/variants may be included. Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Existing candidate medians retained; not a matched enzyme/assay pair. Other substrates are not separately parameterized."},'; ...
    '{"reaction":"r_I1","id":"r_4218","ec":"3.1.1.96","model_ec":"3.1.1.96","equation":"H+ + L-glycine(m_c) + tRNA(Ala)(m_I)  ->  H2O + Glycyl-tRNA(Ala)","direction":-1,"km":{"value":0.00112,"substrate":"D-Tyr-tRNATyr","organism":"Thermus thermophilus","literature":"769120","conditions":"Wild-type; pH 7.5; 37 C.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.1.96"},"kcat":{"value":0.37,"substrate":"D-Tyr-tRNATyr","organism":"Thermus thermophilus","literature":"769120","conditions":"Wild-type; pH 7.5; 37 C.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.1.96"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 3.1.1.96. K_M=0.00112 mM for D-Tyr-tRNATyr; Thermus thermophilus; BRENDA literature 769120; Wild-type; pH 7.5; 37 C. k_cat=0.37 s^-1 for D-Tyr-tRNATyr; Thermus thermophilus; BRENDA literature 769120; Wild-type; pH 7.5; 37 C. Non-human, DIFFERENT substrate and opposite reaction: measured deacylation reference applied to model Gly-tRNA formation. Does not establish reverse synthesis or glycine-binding kinetics. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Non-human, DIFFERENT substrate and opposite reaction: measured deacylation reference applied to model Gly-tRNA formation. Does not establish reverse synthesis or glycine-binding kinetics."},'; ...
    '{"reaction":"r_I2","id":"r_0157","ec":"6.1.1.7","model_ec":"6.1.1.7","equation":"ATP + L-alanine(m_d2) + tRNA(Ala)(m_I)  ->  Ala-tRNA(Ala) + AMP + diphosphate","direction":1,"km":{"value":3.43,"substrate":"L-alanine","organism":"Homo sapiens","literature":"767641","conditions":"pH 8, 37?C","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=6.1.1.7"},"kcat":{"value":0.067,"substrate":"L-alanine","organism":"Homo sapiens","literature":"767641","conditions":"pH 8, 37?C","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=6.1.1.7"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 6.1.1.7. K_M=3.43 mM for L-alanine; Homo sapiens; BRENDA literature 767641; pH 8, 37?C k_cat=0.067 s^-1 for L-alanine; Homo sapiens; BRENDA literature 767641; pH 8, 37?C Alanine chosen from multiple tagged reactants. ATP and tRNA kinetics are not separately parameterized. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Alanine chosen from multiple tagged reactants. ATP and tRNA kinetics are not separately parameterized."},'; ...
    '{"reaction":"r_Nw","id":"r_0002","ec":"1.1.2.4","model_ec":"1.1.2.4","equation":"(R)-lactate + 2 ferricytochrome c  ->  2 ferrocytochrome c + pyruvate(m_out)","direction":1,"km":{"value":0.0682,"substrate":"(R)-lactate","organism":"Oryza sativa Indica Group","literature":"763102","conditions":"pH 8.35; 25 C; phenazine methosulfate and dichlorophenolindophenol acceptors.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=1.1.2.4"},"kcat":{"value":2.6,"substrate":"(R)-lactate","organism":"Oryza sativa Indica Group","literature":"763102","conditions":"pH 8.35; 25 C; phenazine methosulfate and dichlorophenolindophenol acceptors.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=1.1.2.4"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 1.1.2.4. K_M=0.0682 mM for (R)-lactate; Oryza sativa Indica Group; BRENDA literature 763102; pH 8.35; 25 C; phenazine methosulfate and dichlorophenolindophenol acceptors. k_cat=2.6 s^-1 for (R)-lactate; Oryza sativa Indica Group; BRENDA literature 763102; pH 8.35; 25 C; phenazine methosulfate and dichlorophenolindophenol acceptors. Non-human pair uses artificial electron acceptors, not model ferricytochrome c. Does not determine the designed m_c activation law. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Non-human pair uses artificial electron acceptors, not model ferricytochrome c. Does not determine the designed m_c activation law."},'; ...
    '{"reaction":"r_Nb","id":"r_0674","ec":"2.6.1.2","model_ec":"2.6.1.2","equation":"2-oxoglutarate(m_b) + L-alanine  ->  L-glutamate + pyruvate(m_out)","direction":1,"km":{"value":0.20500000000000002,"substrate":"2-oxoglutarate","organism":"Homo sapiens","literature":"636722","conditions":"Median of mitochondrial/cytosolic measurements; pH 7.8; 30 C.","kind":"median","n":2,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=2.6.1.2"},"kcat":{"value":6.09,"substrate":"2-oxoglutarate","organism":"Saccharomyces cerevisiae","literature":"721437","conditions":"L-alanine cosubstrate; pH 8.0; 25 C.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/enzyme.php?ecno=2.6.1.2"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 2.6.1.2. K_M=0.205 mM for 2-oxoglutarate; Homo sapiens; BRENDA literature 636722; Median of mitochondrial/cytosolic measurements; pH 7.8; 30 C. k_cat=6.09 s^-1 for 2-oxoglutarate; Saccharomyces cerevisiae; BRENDA literature 721437; L-alanine cosubstrate; pH 8.0; 25 C. Human Km median and yeast kcat; unpaired organism/assay reference. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Human Km median and yeast kcat; unpaired organism/assay reference."},'; ...
    '{"reaction":"r_Nr","id":"r_0016","ec":"2.2.1.6","model_ec":"2.2.1.6","equation":"2-oxobutanoate + H+ + pyruvate(m_out)  ->  (S)-2-acetyl-2-hydroxybutanoate + carbon dioxide","direction":1,"km":{"value":2.6,"substrate":"pyruvate","organism":"Escherichia coli","literature":"395904","conditions":"Native isoenzyme II.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/literature.php?e=2.2.1.6&r=395904"},"kcat":{"value":66.7,"substrate":"pyruvate","organism":"Escherichia coli","literature":"395904","conditions":"Native enzyme.","kind":"measured","n":1,"source":"https://www.brenda-enzymes.org/literature.php?e=2.2.1.6&r=395904"},"sa":null,"mw":null,"note":"PROVISIONAL REFERENCE; EC selected: 2.2.1.6. K_M=2.6 mM for pyruvate; Escherichia coli; BRENDA literature 395904; Native isoenzyme II. k_cat=66.7 s^-1 for pyruvate; Escherichia coli; BRENDA literature 395904; Native enzyme. Non-human pyruvate/acetolactate assay used for the model pyruvate + 2-oxobutanoate reaction; cosubstrate/product specificity remains an assumption. Sources and record counts: Reference choices sheet. Selected 2026-10-12.","assumptions":"Non-human pyruvate/acetolactate assay used for the model pyruvate + 2-oxobutanoate reaction; cosubstrate/product specificity remains an assumption."}]'; ...
    };
records=jsondecode(strjoin(payload,newline));
end
