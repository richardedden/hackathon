function out = QA_sim(f_results, ph_results)
% f_results is vector of 320 freq offset estimates (in Hz) or minus the correction factors
% ph_results is vector of 320 phase offset estimates (in deg)

load(fullfile(pwd, 'Data', 'Simulated', 'simHermesWithDrift_Feb16.mat'));
fs = fs.';
fs = fs(:);
phs = phs.';
phs = phs(:);

out.error.mean_f = mean(f_results(:)-fs);
out.error.rms_f = sqrt(mean((f_results(:)-fs).^2));
out.error.std_f = std(abs(f_results(:)-fs));

out.error.mean_ph = mean(ph_results(:)-phs);
out.error.rms_ph = sqrt(mean((ph_results(:)-phs).^2));
out.error.std_ph = std(abs(ph_results(:)-phs));

% Normalise outcomes relative to doing nothing
% Factor of 1 is perfect; 0 is did nothing; < 0 is worse than doing nothing
rms_F_error0 = sqrt(mean((zeros(320,1)-fs).^2));
rms_PH_error0 = sqrt(mean((zeros(320,1)-phs).^2));
out.quality.f = 1 - out.error.rms_f/rms_F_error0;
out.quality.ph = 1 - out.error.rms_ph/rms_PH_error0;

out.quality.overall = (out.quality.f + out.quality.ph)/2;

end




