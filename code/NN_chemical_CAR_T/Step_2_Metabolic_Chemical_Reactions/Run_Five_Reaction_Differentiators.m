%% Rerun only Blocks 1-2 after loading your model, currency mask and Block 0.
% Keep your master settings; this script does not initialize COBRA.
addpath(fileparts(mfilename('fullpath')));
Block_1_Differentiator_Mapping;
if isempty(Motifs_Diff_Block_1)
    fprintf('Block 1 saved zero motifs. Inspect its search diagnostics before running Block 2.\n');
else
    Block_2_Differentiator_Mapping;
    Print_Five_Reaction_Differentiators;
end
