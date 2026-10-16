function [outFile,Report] = completeMotifEC(motifFile,varargin)
%COMPLETEMOTIFEC  Fill the blank EC codes in a selected-motif TXT.
%
%   [outFile,Report] = completeMotifEC('Block_7_Selected_Motifs.txt')
%   [outFile,Report] = completeMotifEC(file,'Overwrite',true)
%
%   The file is located with which(), so a bare name works anywhere under
%   the folder added by addingPathParentFolderByName('code').
%
%   For every reaction whose "EC  :" line is empty, the BiGG universal
%   reaction endpoint is queried by the reaction ID and the returned EC code
%   is written into the file. A manual table takes precedence over BiGG, for
%   reactions BiGG does not know (Recon/HMR/RE identifiers often have no
%   universal entry).
%
%   EVERY reaction is checked, not only the blank ones, and the report says
%   whether the code in the file AGREES with BiGG. That matters: in the
%   motif this was written for, two of the six Block 1 reactions carry an EC
%   from the wrong class. PYK (pyruvate kinase) is annotated 3.4.2.1, which
%   is a carboxypeptidase, and METTRANSc (methionine transaminase) is
%   annotated 3.1.1.3, a lipase. A wrong code is worse than a blank one,
%   because BRENDA answers it with kinetics for a different enzyme.
%   Disagreements are REPORTED but not changed unless Overwrite is true.
%
%   OPTIONS
%     'Manual'      struct or N-by-2 cell/string: reaction ID -> EC code.
%                   Defaults to the three resolved for this motif; pass
%                   struct() to disable.
%     'UseBiGG'     true. false works offline, from Manual only.
%     'Overwrite'   false. true also replaces codes that disagree with BiGG.
%     'OutFile'     default: <name>_ECcompleted.txt beside the input
%     'SupplementCsv' default: EC_Supplement.csv beside the output, in the
%                   ReactionID,EC form that applyECSupplement reads
%     'UpdateConflictKeys' true. Keeps the "Enzyme conflict keys" line
%                   consistent with the EC that was written.
%     'Verbose'     true
%
%   Report has one row per reaction: Role, ID, Name, ECInFile, ECFromBiGG,
%   ECWritten, Source, Status.

p = inputParser;
addParameter(p,'Manual',localDefaultManual());
addParameter(p,'UseBiGG',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'Overwrite',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'OutFile','',@(v)ischar(v)||isstring(v));
addParameter(p,'SupplementCsv','',@(v)ischar(v)||isstring(v));
addParameter(p,'UpdateConflictKeys',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'Verbose',true,@(v)islogical(v)&&isscalar(v));
parse(p,varargin{:}); opt = p.Results;
manual = localManualMap(opt.Manual);

%% ---- locate the file ----------------------------------------------------
motifFile = localResolve(motifFile);
say(opt,'Motif file: %s\n',motifFile);
lines = splitlines(string(fileread(motifFile)));

if isempty(char(opt.OutFile))
    [folder,base,ext] = fileparts(motifFile);
    opt.OutFile = fullfile(folder,[base '_ECcompleted' ext]);
end
if isempty(char(opt.SupplementCsv))
    opt.SupplementCsv = fullfile(fileparts(char(opt.OutFile)),'EC_Supplement.csv');
end

%% ---- find every reaction record ----------------------------------------
isRole = ~cellfun(@isempty,regexp(lines,'^\s+(r_[A-Za-z0-9_]+)\s+\(','once'));
roleRows = find(isRole);
assert(~isempty(roleRows),'No "  r_<role> (...)" reaction blocks found in %s.',motifFile);

n = numel(roleRows);
Role = strings(n,1); ID = strings(n,1); RxName = strings(n,1);
ECInFile = strings(n,1); ecRow = zeros(n,1); keyRow = zeros(n,1);
for k = 1:n
    tok = regexp(lines(roleRows(k)),'^\s+(r_[A-Za-z0-9_]+)\s+\(','tokens','once');
    Role(k) = tok{1};
    last = numel(lines);
    if k < n, last = roleRows(k+1) - 1; end
    span = roleRows(k):last;
    ID(k)     = localField(lines,span,'ID');
    RxName(k) = localField(lines,span,'Name');
    [ECInFile(k),ecRow(k)] = localFieldRow(lines,span,'EC');
    keyRow(k) = localRowOf(lines,span,'Enzyme conflict keys');
end
say(opt,'Reactions found: %d (%d with a blank EC).\n',n,sum(strlength(ECInFile)==0));

%% ---- resolve an EC for every reaction ----------------------------------
ECFromBiGG = strings(n,1);
ECWritten  = strings(n,1);
Source     = strings(n,1);
Status     = strings(n,1);
cache = containers.Map('KeyType','char','ValueType','any');

for k = 1:n
    if opt.UseBiGG
        ECFromBiGG(k) = localBiGGEC(ID(k),cache,opt);
    end
    manualEC = "";
    if isKey(manual,char(ID(k))), manualEC = string(manual(char(ID(k)))); end

    blank = strlength(ECInFile(k)) == 0;
    candidate = manualEC;
    source = "manual";
    if strlength(candidate) == 0
        candidate = ECFromBiGG(k); source = "BiGG";
    end

    if blank && strlength(candidate) > 0
        ECWritten(k) = candidate;  Source(k) = source;  Status(k) = "filled";
    elseif blank
        ECWritten(k) = "";  Source(k) = "";  Status(k) = "STILL MISSING";
    elseif strlength(ECFromBiGG(k)) > 0 && ECFromBiGG(k) ~= ECInFile(k)
        if opt.Overwrite
            ECWritten(k) = candidate;  Source(k) = source;  Status(k) = "REPLACED (disagreed)";
        else
            ECWritten(k) = ECInFile(k); Source(k) = "file"; Status(k) = "DISAGREES with BiGG";
        end
    elseif strlength(manualEC) > 0 && manualEC ~= ECInFile(k) && opt.Overwrite
        ECWritten(k) = manualEC; Source(k) = "manual"; Status(k) = "REPLACED (manual)";
    else
        ECWritten(k) = ECInFile(k); Source(k) = "file"; Status(k) = "unchanged";
    end
end

%% ---- rewrite the EC lines ----------------------------------------------
changed = 0;
for k = 1:n
    if ECWritten(k) == ECInFile(k) || strlength(ECWritten(k)) == 0, continue; end
    lines(ecRow(k)) = regexprep(lines(ecRow(k)),'(^\s*EC\s*:\s*).*$', ...
        "$1" + ECWritten(k));
    if opt.UpdateConflictKeys && keyRow(k) > 0
        lines(keyRow(k)) = localMergeKeys(lines(keyRow(k)),ECWritten(k),ECInFile(k));
    end
    changed = changed + 1;
end

Report = table(Role,ID,RxName,ECInFile,ECFromBiGG,ECWritten,Source,Status, ...
    'VariableNames',{'Role','ID','Name','ECInFile','ECFromBiGG','ECWritten','Source','Status'});
disp(Report);

bad = Report.Status == "DISAGREES with BiGG";
if any(bad)
    warning('completeMotifEC:Disagreement', ...
        ['%d reaction(s) carry an EC that BiGG contradicts. They were NOT changed. ' ...
         'A wrong EC sends the BRENDA query to a different enzyme, so resolve these ' ...
         'before querying: %s'],sum(bad),strjoin(Report.ID(bad),', '));
end
stillMissing = Report.Status == "STILL MISSING";
if any(stillMissing)
    warning('completeMotifEC:Missing', ...
        ['No EC found for %s. BiGG has no universal entry for Recon/HMR/RE style ' ...
         'identifiers; resolve these by hand (UniProt via the GPR genes, or ' ...
         'ExPASy by reaction) and pass them through ''Manual''.'], ...
        strjoin(Report.ID(stillMissing),', '));
end

%% ---- write ---------------------------------------------------------------
fid = fopen(char(opt.OutFile),'w','n','UTF-8');
assert(fid ~= -1,'Cannot write %s',char(opt.OutFile));
fprintf(fid,'%s\n',lines{:});
fclose(fid);
outFile = char(opt.OutFile);
say(opt,'\n%d EC line(s) rewritten.\nWrote %s\n',changed,outFile);

% Supplement for applyECSupplement: every reaction with a resolved code.
hasEC = strlength(Report.ECWritten) > 0;
if any(hasEC)
    writetable(table(Report.ID(hasEC),Report.ECWritten(hasEC), ...
        'VariableNames',{'ReactionID','EC'}),char(opt.SupplementCsv));
    say(opt,'Wrote %s (%d reactions) for MASTER_EC_SUPPLEMENT.\n', ...
        char(opt.SupplementCsv),sum(hasEC));
end
end

%% =======================================================================
function ec = localBiGGEC(reactionID,cache,opt)
% Universal BiGG reaction record -> first complete EC code, "" if unknown.
ec = "";
key = char(reactionID);
if isempty(key), return; end
if isKey(cache,key), ec = cache(key); return; end
url = ['http://bigg.ucsd.edu/api/v2/universal/reactions/' key];
try
    record = webread(url,weboptions('Timeout',20,'ContentType','json'));
catch
    cache(key) = ec;    % absent or unreachable; the report says STILL MISSING
    say(opt,'  BiGG: no record for %s\n',key);
    return
end
if isstruct(record) && isfield(record,'database_links')
    links = record.database_links;
    names = fieldnames(links);
    hit = names(~cellfun(@isempty,regexpi(names,'^EC','once')));
    for h = 1:numel(hit)
        entry = links.(hit{h});
        if iscell(entry), entry = entry{1}; end
        if isstruct(entry) && isfield(entry,'id')
            code = string(entry(1).id);
            if ~isempty(regexp(code,'^\d+\.\d+\.\d+\.\d+$','once'))
                ec = code; break
            end
        end
    end
end
cache(key) = ec;
if strlength(ec) > 0, say(opt,'  BiGG: %-12s -> EC %s\n',key,ec); end
end

function map = localDefaultManual()
% Resolved for the motif this function was written for. BiGG gives 2.7.4.3
% for ADK1 and 2.7.1.40 for PYK; METTRANSc has no BiGG EC, and ExPASy
% EC 2.6.1.88 "methionine transaminase" catalyses exactly its equation
% (L-methionine + 2-oxocarboxylate -> 4-methylsulfanyl-2-oxobutanoate +
% L-amino acid). Pass struct() to disable these.
map = struct('ADK1','2.7.4.3','PYK','2.7.1.40','METTRANSc','2.6.1.88');
end

function map = localManualMap(value)
map = containers.Map('KeyType','char','ValueType','any');
if isempty(value), return; end
if isstruct(value)
    f = fieldnames(value);
    for k = 1:numel(f), map(f{k}) = string(value.(f{k})); end
    return
end
value = string(value);
assert(size(value,2) == 2,'Manual must be a struct or an N-by-2 list of ID,EC.');
for k = 1:size(value,1), map(char(value(k,1))) = value(k,2); end
end

function line = localMergeKeys(line,newEC,oldEC)
% Keep "Enzyme conflict keys" in step with the EC that was written.
body = regexprep(line,'^\s*Enzyme conflict keys\s*:\s*','');
parts = strtrim(split(body,";"));
parts = parts(strlength(parts) > 0);
if strlength(oldEC) > 0
    parts = parts(parts ~= "EC:" + oldEC);
end
parts = unique(["EC:" + newEC; parts(:)],'stable');
line = "      Enzyme conflict keys: " + strjoin(parts,"; ");
end

function value = localField(lines,span,name)
[value,~] = localFieldRow(lines,span,name);
end

function [value,row] = localFieldRow(lines,span,name)
value = ""; row = 0;
hit = find(~cellfun(@isempty, ...
    regexp(lines(span),['^\s*' name '\s*:'],'once')),1);
if isempty(hit), return; end
row = span(hit);
value = strtrim(regexprep(lines(row),['^\s*' name '\s*:\s*'],''));
end

function row = localRowOf(lines,span,name)
row = 0;
hit = find(~cellfun(@isempty, ...
    regexp(lines(span),['^\s*' name '\s*:'],'once')),1);
if ~isempty(hit), row = span(hit); end
end

function full = localResolve(name)
name = char(name);
if isfile(name)
    info = dir(name); full = fullfile(info(1).folder,info(1).name); return
end
found = which(name,'-all');
if ~isempty(found)
    full = found{1};
    if numel(found) > 1
        warning('completeMotifEC:Shadowed', ...
            '%d copies of %s are on the path; using\n  %s',numel(found),name,full);
    end
    return
end
error('completeMotifEC:NotFound', ...
    ['Motif file not found: %s\nIt is neither relative to %s nor on the MATLAB ' ...
     'search path. Give the extension, or pass the full path.'],name,pwd);
end

function say(opt,varargin)
if opt.Verbose, fprintf(varargin{:}); end
end
