function [passes,nForbiddenEdges,firstReason] = ...
    passesInternalCrossTalk(S,rec)
%PASSESINTERNALCROSSTALK Test one partial or complete algorithm motif.
%
% Only the canonical reactions already present in REC are inspected. A
% protected algorithm metabolite may be consumed only by its intended
% reaction role(s). Reactions outside REC are deliberately ignored, and
% production by additional network reactions is not restricted.
%
% This is the same internal-consumption rule used by Block 6, but it can be
% called while Blocks 2-4 are assembling candidates so rejected candidates
% do not consume the final MAX_HITS limits.

persistent rxnNames2 metNames2 allowed2
persistent rxnNames3 metNames3 allowed3
persistent rxnNames4 metNames4 allowed4

if isempty(rxnNames2)
    rxnNames2 = { ...
        'r_S1','r_Sg1','r_Ss1','r_F1','r_Fg1','r_Fs1','r_I1', ...
        'r_S2','r_Sg2','r_Ss2','r_F2','r_Fg2','r_Fs2','r_I2'};
    metNames2 = { ...
        'm_in','m_S1','m_F1','m_Sg1','m_Fg1', ...
        'm_S2','m_F2','m_Sg2','m_d'};
    allowed2 = { ...
        [1 4],[2 3],[5 6],[7 8],[7 11], ...
        [9 10],[12 13],14,14};

    rxnNames3 = [rxnNames2,{'r_C','r_CI1','r_CI2'}];
    metNames3 = [metNames2,{'m_c','m_Cs'}];
    allowed3 = { ...
        [1 4],[2 3],[5 6],[7 8],[7 11], ...
        [9 10],[12 13],14,[14 15 17],16,[16 17]};

    rxnNames4 = [rxnNames3,{'r_omega','r_b'}];
    metNames4 = [metNames3,{'m_out'}];
    allowed4 = { ...
        [1 4],[2 3],[5 6],[7 8],[7 11], ...
        [9 10],[12 13],14,[14 15 17],[16 18],[16 17],[]};
end

% Select the fixed layout for the partial motif currently being assembled.
if isfield(rec,'r_b')
    rxnFields = rxnNames4;
    metFields = metNames4;
    allowedRoles = allowed4;
elseif isfield(rec,'r_CI2')
    rxnFields = rxnNames3;
    metFields = metNames3;
    allowedRoles = allowed3;
elseif isfield(rec,'r_I2')
    rxnFields = rxnNames2;
    metFields = metNames2;
    allowedRoles = allowed2;
else
    error('CrossTalk:IncompleteMotif', ...
        'The candidate does not contain a complete Block 2 motif.');
end

rxnIdx = zeros(1,numel(rxnFields));
for r = 1:numel(rxnFields)
    rxnIdx(r) = double(rec.(rxnFields{r}));
end
metIdx = zeros(1,numel(metFields));
for m = 1:numel(metFields)
    metIdx(m) = double(rec.(metFields{m}));
end

[nMets,nRxns] = size(S);
if any(~isfinite(rxnIdx)) || any(rxnIdx < 1) || ...
        any(rxnIdx > nRxns) || any(rxnIdx ~= floor(rxnIdx))
    error('CrossTalk:InvalidReactionIndex', ...
        'The candidate contains an invalid canonical reaction index.');
end
if any(~isfinite(metIdx)) || any(metIdx < 1) || ...
        any(metIdx > nMets) || any(metIdx ~= floor(metIdx))
    error('CrossTalk:InvalidMetaboliteIndex', ...
        'The candidate contains an invalid protected-metabolite index.');
end

consumes = full(S(metIdx,rxnIdx)) < 0;
passes = true;
nForbiddenEdges = 0;
firstReason = "";

for m = 1:numel(metFields)
    allowedIdx = rxnIdx(allowedRoles{m});
    forbidden = consumes(m,:) & ~ismember(rxnIdx,allowedIdx);
    if ~any(forbidden), continue; end

    passes = false;
    nForbiddenEdges = nForbiddenEdges + sum(forbidden);
    if strlength(firstReason) == 0
        firstBad = find(forbidden,1,'first');
        firstReason = sprintf('%s is consumed by unintended role %s', ...
            metFields{m},rxnFields{firstBad});
    end

    % The search blocks request only the logical output. Return at the first
    % violation in that common performance-sensitive case.
    if nargout < 2
        return;
    end
end
end
