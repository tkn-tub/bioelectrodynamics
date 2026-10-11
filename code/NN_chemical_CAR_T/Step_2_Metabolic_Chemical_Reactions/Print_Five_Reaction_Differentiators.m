%% Independent local TXT export for the updated Blocks 1 and 2.
% Does not load a model or initialize COBRA.
addpath(fileparts(mfilename('fullpath')));
if ~exist('POLARITY_MOTIF_TO_PRINT','var'), POLARITY_MOTIF_TO_PRINT = 1; end
polarityPrintData = load('Block_2_differentiator.mat');
if ~isfield(polarityPrintData,'TopologyVersion') || ...
        ~isequal(string(polarityPrintData.TopologyVersion),"five_reaction_differentiators_v2")
    error('PolarityPrint:Topology','Load the new five-plus-five-reaction Block 2 checkpoint.');
end
if isempty(polarityPrintData.Motifs_Diff_Block_2)
    fprintf('No connected two-differentiator motif is available.\n');
else
    polarityPrintText = renderFiveReactionDifferentiators(polarityPrintData,POLARITY_MOTIF_TO_PRINT);
    polarityPrintFile = sprintf('Differentiator_Motif_%d.txt',POLARITY_MOTIF_TO_PRINT);
    polarityPrintHandle = fopen(polarityPrintFile,'w','n','UTF-8');
    if polarityPrintHandle < 0, error('PolarityPrint:Write','Cannot write %s.',polarityPrintFile); end
    polarityPrintCleanup = onCleanup(@() fclose(polarityPrintHandle));
    fprintf(polarityPrintHandle,'%s',polarityPrintText); clear polarityPrintCleanup;
    fprintf('%s\nSaved %s in the current folder.\n',polarityPrintText,polarityPrintFile);
end
clear polarityPrintData polarityPrintText polarityPrintFile polarityPrintHandle;
