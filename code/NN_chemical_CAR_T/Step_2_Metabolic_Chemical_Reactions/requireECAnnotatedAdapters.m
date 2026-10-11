function [paths,originalIndices] = requireECAnnotatedAdapters(model,paths,hasEC,requireEC)
%REQUIREECANNOTATEDADAPTERS Validate passive paths; optionally require ECs.
% Keep stored pathIndex values and validate IDs against the loaded model.
if nargin < 3, hasEC = getReactionECAnnotations(model); end
if nargin < 4, requireEC = false; end
if ~isscalar(requireEC) || ~(islogical(requireEC) || isnumeric(requireEC)) || ...
        ~isreal(requireEC) || ~ismember(requireEC,[0,1])
    error('EC:InvalidAdapterPolicy','BLOCK0_REQUIRE_EC must be true or false.');
end
if ~isstruct(paths)
    error('EC:InvalidAdapters','InputAdapterPaths must be a struct array.');
end
paths = paths(:);
keep = true(numel(paths),1);
for p = 1:numel(paths)
    path = paths(p);
    if ~all(isfield(path,{'anchorMetIdx','reactionIndices'}))
        error('EC:IncompleteAdapter','Adapter %d lacks its anchor or reactions.',p);
    end
    anchor = path.anchorMetIdx;
    idx = path.reactionIndices;
    if ~isnumeric(anchor) || ~isreal(anchor) || ~isscalar(anchor) || ...
            ~isfinite(anchor) || anchor < 1 || anchor > size(model.S,1) || ...
            anchor ~= floor(anchor)
        error('EC:InvalidAdapterIndex','Adapter %d has an invalid endpoint.',p);
    end
    if ~isnumeric(idx) || ~isreal(idx) || ...
            any(~isfinite(idx(:)) | idx(:) < 1 | idx(:) > numel(hasEC) | ...
            idx(:) ~= floor(idx(:)))
        error('EC:InvalidAdapterIndex','Adapter %d has invalid reaction indices.',p);
    end
    if isfield(path,'anchorMetID') && ...
            string(path.anchorMetID) ~= string(model.mets(anchor))
        error('EC:AdapterModelMismatch','Adapter %d endpoint ID disagrees with this model.',p);
    end
    if isfield(path,'reactionPath') && strlength(string(path.reactionPath)) > 0
        savedIDs = strtrim(split(string(path.reactionPath),"->"));
        savedIDs = regexprep(savedIDs,'\s*\(reverse\)$','');
        if ~isequal(savedIDs(:),string(model.rxns(double(idx(:)))))
            error('EC:AdapterModelMismatch','Adapter %d reaction IDs disagree with this model.',p);
        end
    end
    keep(p) = ~requireEC || all(hasEC(double(idx(:))));
end
originalIndices = find(keep);
paths = paths(keep);
end
