function outFile = writeMotifsForBrenda(rankingFile,motifsFile,varargin)
%WRITEMOTIFSFORBRENDA  Merge the ranking and the motif listing into one TXT
%   that brendaKineticsTable can parse.
%
%   outFile = writeMotifsForBrenda('CAMP_Input_Ranking.txt', ...
%                                  'Selected_Circuit_Motifs.txt')
%   outFile = writeMotifsForBrenda(ranking,motifs,'Select','reachable')
%
%   Writes into the subfolder Circuit_Search_Results beside this file, so the
%   result sits with the rest of the search output.
%
%   The motif listing already uses the MOTIF / BLOCK / r_<role> layout that
%   brendaKineticsTable expects, but three things stop it being usable:
%
%   1. The cAMP connector sits under the heading "cAMP INPUT CONNECTION",
%      not "BLOCK 0:", so parseMotifs gives those reactions no block number
%      and the record is rejected.
%   2. There is no "Metabolite path:" / "Reaction path:" pair, which is how
%      brendaKineticsTable works out WHICH substrate each adapter step
%      consumes. Without it those rows come back "primary substrate
%      unresolved". Both lines are derived here from the printed Path.
%   3. All 59 motifs are listed, including the 56 whose distance is Inf.
%      Querying BRENDA for every one of them is pointless.
%
%   OPTIONS
%     'Select'      "best" (default) the single best-ranked motif
%                   "reachable"      every motif with a finite distance
%                   "all"            every motif in the listing
%                   numeric vector   ORIGINAL motif numbers, in that order
%     'OutFolder'   default: Circuit_Search_Results beside this file
%     'OutName'     default: Block_7_Selected_Motifs.txt
%     'RenameAdapters'  true (default) renames r_CAMP_k to r_adapter_k
%     'Verbose'     true
%
%   The motifs are renumbered 1..N in the output, in ranking order, while the
%   original motif number is kept in a comment line so nothing is lost.

p = inputParser;
addParameter(p,'Select',"best");
addParameter(p,'OutFolder','',@(v)ischar(v)||isstring(v));
addParameter(p,'OutName','Block_7_Selected_Motifs.txt',@(v)ischar(v)||isstring(v));
addParameter(p,'RenameAdapters',true,@(v)islogical(v)&&isscalar(v));
addParameter(p,'Verbose',true,@(v)islogical(v)&&isscalar(v));
parse(p,varargin{:}); opt = p.Results;

% Resolve through the MATLAB path, so a bare file name works anywhere under
% the folder added by addingPathParentFolderByName('code').
rankingFile = localResolve(rankingFile,'Ranking file');
motifsFile  = localResolve(motifsFile,'Motif listing');
say(opt,'Ranking file : %s\n',rankingFile);
say(opt,'Motif listing: %s\n',motifsFile);

if isempty(char(opt.OutFolder))
    here = fileparts(mfilename('fullpath'));
    if isempty(here), here = pwd; end
    opt.OutFolder = fullfile(here,'Circuit_Search_Results');
end
if ~isfolder(opt.OutFolder), mkdir(opt.OutFolder); end
outFile = fullfile(char(opt.OutFolder),char(opt.OutName));

%% ---- ranking: original motif number, rank, distance, status -------------
rankLines = splitlines(string(fileread(rankingFile)));
hit = regexp(rankLines, ...
    '^Rank\s+(\d+)\s*\|\s*original motif\s+(\d+)\s*\|\s*distance\s+(\S+)\s*\|\s*(\S+)', ...
    'tokens','once');
keep = ~cellfun(@isempty,hit);
hit = hit(keep);
Ranking = table( ...
    cellfun(@(t) str2double(t{1}),hit), ...
    cellfun(@(t) str2double(t{2}),hit), ...
    cellfun(@(t) str2double(t{3}),hit), ...
    cellfun(@(t) string(t{4}),hit), ...
    'VariableNames',{'Rank','OriginalMotif','Distance','Status'});
assert(~isempty(Ranking),'No "Rank N | original motif M" lines found in %s.',char(rankingFile));
Ranking = sortrows(Ranking,{'Distance','Rank','OriginalMotif'});
say(opt,'Ranking: %d motifs, %d reachable (finite distance).\n', ...
    height(Ranking),sum(isfinite(Ranking.Distance)));

%% ---- motif listing: one text block per MOTIF ----------------------------
motifLines = splitlines(string(fileread(motifsFile)));
starts = find(~cellfun(@isempty,regexp(motifLines,'^MOTIF\s+\d+\s*$','once')));
assert(~isempty(starts),'No "MOTIF n" headings found in %s.',char(motifsFile));
motifNumber = zeros(numel(starts),1);
for k = 1:numel(starts)
    tok = regexp(motifLines(starts(k)),'^MOTIF\s+(\d+)','tokens','once');
    motifNumber(k) = str2double(tok{1});
end
stops = [starts(2:end)-2; numel(motifLines)];   % drop the trailing "====" rule
stops = max(stops,starts);

%% ---- which motifs to write ---------------------------------------------
if isnumeric(opt.Select)
    wanted = opt.Select(:);
else
    switch lower(string(opt.Select))
        case "best"
            wanted = Ranking.OriginalMotif(1);
        case "reachable"
            wanted = Ranking.OriginalMotif(isfinite(Ranking.Distance));
        case "all"
            wanted = Ranking.OriginalMotif;
        otherwise
            error('Select must be "best", "reachable", "all" or a vector of motif numbers.');
    end
end
assert(~isempty(wanted),'No motif selected. With "reachable", none has a finite distance.');

%% ---- write --------------------------------------------------------------
fid = fopen(outFile,'w','n','UTF-8');
assert(fid ~= -1,'Cannot write %s',outFile);
closer = onCleanup(@() fclose(fid));

fprintf(fid,'############################################################\n');
fprintf(fid,'SELECTED CHEMICAL CIRCUIT MOTIFS FOR THE BRENDA QUERY\n');
fprintf(fid,'Built by writeMotifsForBrenda on %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf(fid,'Ranking source : %s\n',char(rankingFile));
fprintf(fid,'Motif source   : %s\n',char(motifsFile));
fprintf(fid,'Selection      : %s\n',localSelectText(opt.Select));
fprintf(fid,'Motifs written : %d of %d\n',numel(wanted),numel(starts));
fprintf(fid,'Block 0 is the cAMP-to-m_in connector, printed in full.\n');
fprintf(fid,'############################################################\n');

written = 0;
for w = 1:numel(wanted)
    source = find(motifNumber == wanted(w),1);
    if isempty(source)
        warning('writeMotifsForBrenda:Missing', ...
            'Motif %d is ranked but absent from the listing; skipped.',wanted(w));
        continue
    end
    rankRow = Ranking(Ranking.OriginalMotif == wanted(w),:);
    written = written + 1;
    body = motifLines(starts(source):stops(source));
    body = localConnectorToBlock0(body,opt.RenameAdapters);

    fprintf(fid,'\n============================================================\n');
    fprintf(fid,'MOTIF %d\n',written);
    fprintf(fid,'============================================================\n');
    fprintf(fid,'original motif: %d\n',wanted(w));
    if ~isempty(rankRow)
        fprintf(fid,'input connection rank: %d\n',rankRow.Rank(1));
        fprintf(fid,'cAMP conversion steps: %g\n',rankRow.Distance(1));
        fprintf(fid,'input connection status: %s\n',rankRow.Status(1));
    end
    % body(1) is the "MOTIF n" line of the source file, replaced above.
    for L = 2:numel(body)
        fprintf(fid,'%s\n',body(L));
    end
    say(opt,'  motif %d (original %d, distance %g) written.\n', ...
        written,wanted(w),rankRow.Distance(1));
end
assert(written > 0,'Nothing was written: no selected motif was found in the listing.');

fprintf(fid,'\nASSUMPTIONS AND INTERFACES\n');
fprintf(fid,'Block 0 reactions are the route to the circuit; no rate constant is read for them.\n');
fprintf(fid,'Kinetics validated: 0. Enzyme parameters and feasibility require separate evaluation.\n');
clear closer

say(opt,'\nWrote %d motif(s) to\n  %s\n',written,outFile);
say(opt,'Next: brendaKineticsTable(''%s'',''kinetics_BRENDA.xlsx'',email,password);\n', ...
    char(opt.OutName));
end

%% =======================================================================
function body = localConnectorToBlock0(body,renameAdapters)
% Give the cAMP connector a BLOCK 0 heading and the two path lines that
% brendaKineticsTable uses to identify each adapter step's own substrate.
head = find(~cellfun(@isempty,regexp(body,'^\s*cAMP INPUT CONNECTION','once')),1);
if isempty(head), return; end
body(head) = "BLOCK 0: cAMP INPUT CONNECTION (OUTPUT m_in)";

% Metabolite path from the printed "Path: a -> b -> c" line.
pathRow = find(~cellfun(@isempty,regexp(body,'^\s*Path\s*:','once')),1);
if isempty(pathRow), return; end
metPath = strtrim(regexprep(body(pathRow),'^\s*Path\s*:\s*',''));
if strlength(metPath) == 0, return; end

% Reaction path from the IDs of the connector steps, in printed order.
rxnRows = find(~cellfun(@isempty,regexp(body,'^\s+r_(CAMP|adapter)_\d+\s+\(','once')));
rxnIDs = strings(0,1);
for k = 1:numel(rxnRows)
    idRow = rxnRows(k) + find(~cellfun(@isempty, ...
        regexp(body(rxnRows(k)+1:min(rxnRows(k)+6,numel(body))),'^\s*ID\s*:','once')),1);
    if isempty(idRow), continue; end
    rxnIDs(end+1,1) = strtrim(regexprep(body(idRow),'^\s*ID\s*:\s*','')); %#ok<AGROW>
end
if isempty(rxnIDs), return; end

insert = ["  Metabolite path: " + metPath; ...
          "  Reaction path: " + strjoin(rxnIDs," -> ")];
body = [body(1:pathRow); insert; body(pathRow+1:end)];

if renameAdapters
    body = regexprep(body,'\br_CAMP_(\d+)\b','r_adapter_$1');
end
end

function full = localResolve(name,what)
% A direct path first, then the MATLAB search path. which() finds a data
% file only when it is given with its extension, which these always are.
name = char(name);
if isfile(name), full = localFullPath(name); return; end
found = which(name);
if ~isempty(found) && isfile(found), full = found; return; end
% which() can return the first of several shadowing copies; say so.
allFound = which(name,'-all');
if ~isempty(allFound)
    full = allFound{1};
    if numel(allFound) > 1
        warning('writeMotifsForBrenda:Shadowed', ...
            '%d copies of %s are on the path; using\n  %s', ...
            numel(allFound),name,full);
    end
    return
end
error('writeMotifsForBrenda:NotFound', ...
    ['%s not found: %s\nIt is neither a path relative to %s nor on the MATLAB ', ...
     'search path. Check the name and extension, or pass the full path. ', ...
     'addingPathParentFolderByName(''code'') only adds folders that exist ', ...
     'under that tree when it runs.'],what,name,pwd);
end

function full = localFullPath(name)
info = dir(name);
full = fullfile(info(1).folder,info(1).name);
end

function text = localSelectText(value)
if isnumeric(value)
    text = "explicit motif list: " + strjoin(string(value(:)).',", ");
else
    text = string(value);
end
end

function say(opt,varargin)
if opt.Verbose, fprintf(varargin{:}); end
end
