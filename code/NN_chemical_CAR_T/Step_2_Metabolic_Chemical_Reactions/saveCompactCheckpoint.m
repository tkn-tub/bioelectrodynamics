function fileInfo = saveCompactCheckpoint(filePath,payload)
%SAVECOMPACTCHECKPOINT Save a scalar payload struct as named MAT variables.
% Store column-based tables once, not their row-struct duplicates or caches.
% Version 7 compression is used unless a variable exceeds 1.8 GB in memory.
% A temporary file protects an existing checkpoint if serialization fails.

if ~isstruct(payload) || ~isscalar(payload)
    error('Checkpoint:InvalidPayload','payload must be a scalar struct.');
end
filePath = char(string(filePath));
[folder,~,~] = fileparts(filePath);
if isempty(folder), folder = pwd; end
if ~isfolder(folder), mkdir(folder); end

fields = fieldnames(payload);
versionFlag = '-v7';
for k = 1:numel(fields)
    value = payload.(fields{k}); %#ok<NASGU>
    variableInfo = whos('value');
    if variableInfo.bytes >= 1.8e9
        versionFlag = '-v7.3';
        break;
    end
end

temporaryFile = [tempname(folder),'.mat'];
temporaryCleanup = onCleanup(@() localRemoveTemporary(temporaryFile));
save(temporaryFile,'-struct','payload',versionFlag);
[ok,message] = movefile(temporaryFile,filePath,'f');
if ~ok
    error('Checkpoint:MoveFailed','Could not finalize %s: %s', ...
        filePath,message);
end
clear temporaryCleanup;
fileInfo = dir(filePath);
fprintf('Compact checkpoint: %.2f MB (%s)\n%s\n', ...
    fileInfo.bytes/1e6,versionFlag,filePath);
end

function localRemoveTemporary(filePath)
if isfile(filePath), delete(filePath); end
end
