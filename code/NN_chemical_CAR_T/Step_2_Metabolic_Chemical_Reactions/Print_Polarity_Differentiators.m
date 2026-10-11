%% Independent local TXT export for the updated Blocks 1 and 2.
% Does not load a model or initialize COBRA.
addpath(fileparts(mfilename('fullpath')));
POLARITY_MOTIF_TO_PRINT = 1;
polarityPrintData = load('Block_2_differentiator.mat');
if ~isfield(polarityPrintData,'TopologyVersion') || ...
        ~isequal(string(polarityPrintData.TopologyVersion),"polarity_3plus5_differentiators_v1")
    error('PolarityPrint:Topology','Load the new three-plus-five-reaction Block 2 checkpoint.');
end
if isempty(polarityPrintData.Motifs_Diff_Block_2)
    fprintf('No connected two-differentiator motif is available.\n');
else
    polarityPrintText = renderPolarityDifferentiators(polarityPrintData,POLARITY_MOTIF_TO_PRINT);
    polarityPrintFile = sprintf('Differentiator_Motif_%d.txt',POLARITY_MOTIF_TO_PRINT);
    polarityPrintHandle = fopen(polarityPrintFile,'w','n','UTF-8');
    if polarityPrintHandle < 0, error('PolarityPrint:Write','Cannot write %s.',polarityPrintFile); end
    polarityPrintCleanup = onCleanup(@() fclose(polarityPrintHandle));
    fprintf(polarityPrintHandle,'%s',polarityPrintText); clear polarityPrintCleanup;
    fprintf('%s\nSaved %s in the current folder.\n',polarityPrintText,polarityPrintFile);
end
clear polarityPrintData polarityPrintText polarityPrintFile polarityPrintHandle;
