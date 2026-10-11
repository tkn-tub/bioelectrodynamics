function [rows,originalIndices] = requireECAnnotatedRows(rows,hasEC,roles,paths)
%REQUIREECANNOTATEDROWS Check upstream rows before new hit budgets are used.
% Original row numbers are returned for stable block*_idx provenance.
% When paths are supplied, require an eligible Block 0 path to m_in too.
if ~istable(rows) && ~isstruct(rows)
    error('EC:InvalidRows','Expected a motif table or struct array.');
end
hasEC = logical(hasEC(:));
if istable(rows), n = height(rows); else, rows = rows(:); n = numel(rows); end
keep = true(n,1);
for k = 1:numel(roles)
    if n == 0, break; end
    role = roles{k};
    if istable(rows)
        if ~ismember(role,rows.Properties.VariableNames)
            error('EC:MissingRole','Missing upstream reaction role %s.',role);
        end
        idx = rows.(role);
    else
        if ~isfield(rows,role)
            error('EC:MissingRole','Missing upstream reaction role %s.',role);
        end
        idx = vertcat(rows.(role));
    end
    if ~isnumeric(idx) || ~isreal(idx) || numel(idx) ~= n || ...
            any(~isfinite(idx(:)) | idx(:) < 1 | idx(:) > numel(hasEC) | ...
            idx(:) ~= floor(idx(:)))
        error('EC:InvalidIndex','%s contains invalid indices for this model.',role);
    end
    keep = keep & hasEC(double(idx(:)));
end
if nargin >= 4 && n > 0
    if istable(rows), inputs = rows.m_in; else, inputs = vertcat(rows.m_in); end
    anchors = double([paths.anchorMetIdx]).';
    keep = keep & ismember(double(inputs(:)),anchors);
    if istable(rows)
        hasPathIndex = ismember('adapterPathIndex',rows.Properties.VariableNames);
    else
        hasPathIndex = isfield(rows,'adapterPathIndex');
    end
    if hasPathIndex && ~isempty(paths)
        if ~isfield(paths,'pathIndex')
            error('EC:MissingPathIndex','Saved adapter path indices are required.');
        end
        if istable(rows), idx = rows.adapterPathIndex; else, idx = vertcat(rows.adapterPathIndex); end
        keep = keep & ismember([double(inputs(:)),double(idx(:))], ...
            [anchors,double([paths.pathIndex]).'],'rows');
    end
end
originalIndices = find(keep);
if istable(rows), rows = rows(keep,:); else, rows = rows(keep); end
end
