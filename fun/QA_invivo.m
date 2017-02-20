function out = QA_invivo(MRS_struct)

% Set some initial variables
n = numel(MRS_struct.gabafile);
freq = MRS_struct.spec.freq;
GABA_diff = MRS_struct.spec.vox1.GABAGlx.diff;
GABA_diff_noalign = MRS_struct.spec.vox1.GABAGlx.diff_noalign;
GSH_diff = MRS_struct.spec.vox1.GSH.diff;
GSH_diff_noalign = MRS_struct.spec.vox1.GSH.diff_noalign;

% Frequency ranges
lb = find(freq >= 3.18);
ub = find(freq <= 3.20);
ChoRange = intersect(lb,ub);
lb = find(freq >= 10);
ub = find(freq <= 11);
noiseRange = intersect(lb,ub);

for ii = 1:n
    
    subplot(2,1,1);
    plot(freq(:), GABA_diff_noalign(ii,:), 'r', freq(:), GABA_diff(ii,:), 'b');
    set(gca,'xdir', 'reverse', 'xlim', [3.2-0.5 3.2+0.5]);
    legend('PRE','POST');
    
    out.SA.GABA.prealign(ii) = std(GABA_diff_noalign(ii,ChoRange));
    out.SA.GABA.postalign(ii) = std(GABA_diff(ii,ChoRange));
    out.GABA_noise(ii,:) = std(GABA_diff(ii,noiseRange));
    
    subplot(2,1,2);
    plot(freq(:), GSH_diff_noalign(ii,:), 'r', freq(:), GSH_diff(ii,:), 'b');
    set(gca,'xdir', 'reverse', 'xlim', [3.2-0.5 3.2+0.5]);
    
    out.SA.GSH.prealign(ii) = std(GSH_diff_noalign(ii,ChoRange));
    out.SA.GSH.postalign(ii) = std(GSH_diff(ii,ChoRange));
    out.GSH_noise(ii,:) = std(GSH_diff(ii,noiseRange));
    
    pause;
    
end

