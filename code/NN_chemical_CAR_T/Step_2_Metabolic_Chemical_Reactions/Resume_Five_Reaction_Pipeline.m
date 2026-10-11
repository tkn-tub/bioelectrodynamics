%% Run after your master stops with model/isBadMet still in the workspace.
% Reuses CURRENT Block_2_differentiator.mat. It never runs Blocks 0, 1 or 2.
addpath(fileparts(mfilename('fullpath')));
RUN_BLOCK3 = true;
RUN_BLOCK4 = true;
RUN_BLOCK5 = true;
RUN_BLOCK6 = true;
if ~exist('MAP_COUNTER_CS_LEAK','var'), MAP_COUNTER_CS_LEAK = false; end
% false: three mapped counter reactions, complete external Cs buffer support.
% true : positive external Cs supply, with additional mapped basal Cs removal.
Continue_Minimal_Pipeline;
