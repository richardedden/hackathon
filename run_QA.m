%%%%%%%%%%%%%%%%%%%%%% Run quality analysis %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script *must* be run from within the 'hackathon' directory %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%%%%%%%%%%%%%%%%%% Set some input arguments %%%%%%%%%%%%%%%%%%%%%
showPlots = true; % true = show pre-/post-alignment Cho subtraction artifact plots (in vivo data)
                   %        and freq/phase offset estimates (simulated data)
simInd = 41; % file number in GannetLoad output structure that corresponds to the simulated data (will usually be the last file [#41])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add 'fun' directory to user's search path
addpath(fullfile(pwd,'fun'));

% Load GannetLoad output data
[filename, pathname] = uigetfile({'*.mat','MAT-files (*.mat)'},'Select GannetLoad output data');
load(fullfile(pathname, filename));

% Run QA
QA.in_vivo = QA_InVivo(MRS, showPlots, simInd);
QA.sim = QA_Sim(MRS, showPlots, simInd);

% Save QA data with filename corresponding to GannetLoad output data file (*_run#_QA.mat)
[pathname2, filename2, ext] = fileparts(fullfile(pathname, filename));
newFilename = fullfile(pathname2, [filename2 '_QA' ext]);
save(newFilename,'QA');
