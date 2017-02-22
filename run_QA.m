%%% Run quality analysis %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script *must* be run from within the 'hackathon' directory %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% Add 'fun' directory to user's search path
addpath(fullfile(pwd,'fun'));

%%%%%%%%%%%%%%%%%%%% SET SOME INPUT ARGUMENTS %%%%%%%%%%%%%%%%%%%%%
show_plots = true; % true = show pre-/post-alignment Cho subtraction artifact plots
indSim = 41; % file number in GannetLoad output structure that corresponds to the simulated data (will usually be the last file [#41])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load output data
[filename, pathname] = uigetfile({'*.mat','MAT-files (*.mat)'},'Select hackathon output data');
load(fullfile(pathname, filename));

% Run QA
QA.in_vivo = QA_InVivo(MRS, show_plots, indSim);
QA.sim = QA_Sim(MRS, indSim);

% Save QA data with filename corresponding to loaded output data file (*_run#_QA.mat)
[pathname2, filename2, ext] = fileparts(fullfile(pathname, filename));
newFilename = fullfile(pathname2, [filename2 '_QA' ext]);
save(newFilename,'QA');
