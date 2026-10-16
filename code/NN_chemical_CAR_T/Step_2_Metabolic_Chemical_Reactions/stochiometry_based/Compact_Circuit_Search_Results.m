function [compactFile,stats] = Compact_Circuit_Search_Results(sourceFile,compactFile)
%COMPACT_CIRCUIT_SEARCH_RESULTS Lossless conversion to compressed MAT v7.
% Compact_Circuit_Search_Results('/full/path/Circuit_Search_Results.mat')
% Creates Circuit_Search_Results_compact.mat, retaining every variable.
if nargin < 1 || isempty(sourceFile)
    error('Circuit:CompactSource','Supply the path of your existing results MAT file.');
end
sourceFile = char(sourceFile);
[ok,sourceAttributes] = fileattrib(sourceFile);
if ~ok || sourceAttributes.directory
    error('Circuit:CompactSource','Cannot find the results file: %s',sourceFile);
end
sourceFile = sourceAttributes.Name;
[folder,base] = fileparts(sourceFile);
if nargin < 2 || isempty(compactFile)
    compactFile = fullfile(folder,[base,'_compact.mat']);
else
    compactFile = char(compactFile);
end
if exist(compactFile,'file') ~= 0
    error('Circuit:CompactDestination','Choose a new output file; it already exists: %s',compactFile);
end
entries = whos('-file',sourceFile);
if any([entries.bytes] >= 2^31)
    error('Circuit:CompactSize','A variable exceeds MAT v7 capacity; retain v7.3 for that dataset.');
end
data = load(sourceFile);
[outputFolder,~,~] = fileparts(compactFile);
if isempty(outputFolder), outputFolder = pwd; end
temporaryFile = [tempname(outputFolder),'.mat'];
cleanup = onCleanup(@() removeTemporary(temporaryFile)); %#ok<NASGU>
save(temporaryFile,'-struct','data','-v7');
% Verify values as well as variable names before making the output available.
restored = load(temporaryFile);
names = fieldnames(data);
if ~isequal(sort(names),sort(fieldnames(restored)))
    error('Circuit:CompactVerification','Variable names changed during conversion.');
end
for k = 1:numel(names)
    if ~isequaln(data.(names{k}),restored.(names{k}))
        error('Circuit:CompactVerification','Variable changed during conversion: %s',names{k});
    end
end
if exist(compactFile,'file') ~= 0
    error('Circuit:CompactDestination','Output appeared during conversion; choose a new destination.');
end
[ok,message] = movefile(temporaryFile,compactFile);
if ~ok, error('Circuit:CompactMove','Could not finalize compact file: %s',message); end
original = dir(sourceFile); compact = dir(compactFile);
stats.originalBytes = original.bytes; stats.compactBytes = compact.bytes;
stats.percentSaved = 100*(1-compact.bytes/original.bytes);
stats.verifiedVariables = numel(names);
fprintf('Original: %.2f MiB | compact: %.2f MiB | reduction: %.1f%%\n', ...
    original.bytes/2^20,compact.bytes/2^20,stats.percentSaved);
fprintf('Verified %d unchanged variables. Compact results: %s\n',numel(names),compactFile);
end

function removeTemporary(file)
if exist(file,'file') == 2, delete(file); end
end
