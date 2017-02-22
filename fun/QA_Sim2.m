function out = QA_Sim(MRS_struct, indSim)
% f_results is vector of Navg freq offset estimates or minus the correction factors (in Hz)
% ph_results is vector of Navg phase offset estimates (in deg)

% Load freq and phase offsets added to simulated data (ground truths)
load(fullfile(pwd, 'data', 'simulated', 'simHermesWithDrift_Feb20.mat'));
fs = fs.';
fs = fs(:);
phs = phs.';
phs = phs(:);

% Calculate QA metrics
out.error.mean_f = mean(MRS_struct.out.f_results(indSim,:)'-fs);
out.error.rms_f = sqrt(mean((MRS_struct.out.f_results(indSim,:)'-fs).^2));
out.error.std_f = std(abs(MRS_struct.out.f_results(indSim,:)'-fs));

out.error.mean_ph = mean(MRS_struct.out.ph_results(indSim,:)'-phs);
out.error.rms_ph = sqrt(mean((MRS_struct.out.ph_results(indSim,:)'-phs).^2));
out.error.std_ph = std(abs(MRS_struct.out.ph_results(indSim,:)'-phs));

% Normalize outcomes relative to doing nothing
% 1 = perfect; 0 = did nothing; < 0 = worse than doing nothing
rms_f0 = sqrt(mean((zeros(320,1)-fs).^2));
rms_ph0 = sqrt(mean((zeros(320,1)-phs).^2));
out.quality.f = 1 - out.error.rms_f/rms_f0;
out.quality.ph = 1 - out.error.rms_ph/rms_ph0;

out.quality.overall = (out.quality.f + out.quality.ph)/2;