%% Independent printer: no model load and no COBRA initialization required.
% Edit this single row number; both files and the TXT stay in the local folder.
minimalMotifSelection = 1;
minimalData = load('Block_2_differentiator.mat');
if ~isfield(minimalData,'TopologyVersion') || ...
        minimalData.TopologyVersion ~= "minimal_three_reaction_v1"
    error('MinimalPrint:OldFile','This printer requires the new Block 2 checkpoint.');
end
if isempty(minimalData.Motifs_Diff_Block_2)
    fprintf('No connected two-differentiator motif is available.\n');
else
    minimalText = renderMinimalDifferentiator(minimalData,minimalMotifSelection);
    fprintf('%s',minimalText);
    minimalFile = sprintf('Differentiator_Motif_%d.txt',minimalMotifSelection);
    minimalHandle = fopen(minimalFile,'w','n','UTF-8');
    if minimalHandle < 0, error('MinimalPrint:Write','Cannot write %s.',minimalFile); end
    minimalCleanup = onCleanup(@() fclose(minimalHandle));
    fprintf(minimalHandle,'%s',minimalText);
    clear minimalCleanup;
    fprintf('Saved %s\n',minimalFile);
end
clear minimalData minimalText minimalHandle;
