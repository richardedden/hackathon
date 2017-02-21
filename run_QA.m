%%% Run quality analysis functions %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script *must* be run from within the "hackathon" directory %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% Add directories to user's search path
addpath(fullfile(pwd,'fun'));

%%%%%%%%%%%% Set some initial variables %%%%%%%%%%%%
show_plots = true; % true = show Cho subtraction artifact figures
indSim = 41; % file number in MRS_struct that corresponds to the simulated data (usually the last file run [#41])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load output data
[filename,pathname] = uigetfile({'*.mat','MAT-files (*.mat)'},'Select hackathon output data');
load(fullfile(pathname,filename));

% Run QA
QA.in_vivo = QA_InVivo(MRS, show_plots, indSim);
QA.sim = QA_Sim(MRS, indSim);

% Save QA data with filename corresponding to loaded output data file (*_run#_QA.mat)
[pathname2,filename2,ext] = fileparts(fullfile(pathname,filename));
newFilename = fullfile(pathname2,[filename2 '_QA' ext]);
save(newFilename,'QA');
