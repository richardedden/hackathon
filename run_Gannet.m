%%% Batch-run Gannet %%%

% Add directories to user's search path
clear;
addpath(fullfile(pwd,'Gannet3.0'), fullfile(pwd,'Data'));

% Create list of files to process
exp = {'GO','IW','KC','MS'};
subj = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10'};

metab = cell(numel(exp)*numel(subj),1);
water = metab;
c = 1;
for ii = 1:length(exp)
    for jj = 1:length(subj)        
        metab{c} = fullfile(pwd, 'Data', 'In Vivo', exp{ii}, subj{jj}, [exp{ii} '_' subj{jj} '_HERMES_act.sdat']);
        water{c} = fullfile(pwd, 'Data', 'In Vivo', exp{ii}, subj{jj}, [exp{ii} '_' subj{jj} '_HERMES_ref.sdat']);
        c = c + 1;
    end
end
metab{end+1} = fullfile(pwd, 'Data', 'Simulated', 'simHermesWithDrift_Feb16.mat');
water{end+1} = fullfile(pwd, 'Data', 'Simulated', 'simHermesWithDrift_Feb16.mat');

% Run Gannet
MRS = GannetLoad(metab, water);
% filename = ['hackathon_data_' datestr(clock,'yymmdd-HHMM')];
filename = 'hackathon_data';
save(filename,'MRS');