function out = QA_InVivo(MRS_struct, showPlots, simInd)

% Set some initial variables
n = 1:numel(MRS_struct.gabafile); % spectra to run QA on
n(simInd) = []; % exclude simulated data
freq = MRS_struct.spec.freq;
GABA_diff_noalign = real(MRS_struct.spec.vox1.GABAGlx.diff_noalign);
GABA_diff = real(MRS_struct.spec.vox1.GABAGlx.diff);
GSH_diff_noalign = real(MRS_struct.spec.vox1.GSH.diff_noalign);
GSH_diff = real(MRS_struct.spec.vox1.GSH.diff);

grey = [0.6 0.6 0.6];
shading = 0.3;
w = 0.3;
h = 0.65;
l = (1-w)/2;
b = 0.3;

for ii = n
    
    % Set frequency ranges
    lb = find(freq(ii,:) >= 3.185-0.05+0.02);
    ub = find(freq(ii,:) <= 3.185+0.05+0.02);
    ChoRange = intersect(lb,ub);
    lb = find(freq(ii,:) >= 10);
    ub = find(freq(ii,:) <= 11);
    noiseRange = intersect(lb,ub);
    lb = find(freq(ii,:) >= 3.2-0.75);
    ub = find(freq(ii,:) <= 3.2+0.75);
    plotRange = intersect(lb,ub);
    
    % Calculate QA metrics
    out.SA.GABA.prealign_std(ii) = std(GABA_diff_noalign(ii,ChoRange));
    out.SA.GABA.prealign_max(ii) = max(GABA_diff_noalign(ii,ChoRange));
    out.SA.GABA.postalign_std(ii) = std(GABA_diff(ii,ChoRange));
    out.SA.GABA.postalign_max(ii) = max(GABA_diff(ii,ChoRange));
    out.noise.GABA.prealign(ii) = std(detrend(GABA_diff_noalign(ii,noiseRange)));
    out.noise.GABA.postalign(ii) = std(detrend(GABA_diff(ii,noiseRange)));
    out.SA.GABA.prealign_SNR(ii) = out.SA.GABA.prealign_max(ii)/out.noise.GABA.prealign(ii);
    out.SA.GABA.postalign_SNR(ii) = out.SA.GABA.postalign_max(ii)/out.noise.GABA.postalign(ii);
    
    out.SA.GSH.prealign_std(ii) = std(GSH_diff_noalign(ii,ChoRange));
    out.SA.GSH.prealign_max(ii) = max(GSH_diff_noalign(ii,ChoRange));
    out.SA.GSH.postalign_std(ii) = std(GSH_diff(ii,ChoRange));
    out.SA.GSH.postalign_max(ii) = max(GSH_diff(ii,ChoRange));
    out.noise.GSH.prealign(ii) = std(detrend(GSH_diff_noalign(ii,noiseRange)));
    out.noise.GSH.postalign(ii) = std(detrend(GSH_diff(ii,noiseRange)));
    out.SA.GSH.prealign_SNR(ii) = out.SA.GSH.prealign_max(ii)/out.noise.GSH.prealign(ii);
    out.SA.GSH.postalign_SNR(ii) = out.SA.GSH.postalign_max(ii)/out.noise.GSH.postalign(ii);
    
    % Normalize outcomes relative to pre-aligned data
    % 1 = perfect; 0 = did nothing; < 0 = worse than no alignment
    out.SA.GABA.quality(ii) = 1 - out.SA.GABA.postalign_std(ii)/out.SA.GABA.prealign_std(ii);
    out.SA.GSH.quality(ii) = 1 - out.SA.GSH.postalign_std(ii)/out.SA.GSH.prealign_std(ii);
    
    % Plots
    if showPlots
        
        h1 = figure(111);
        set(h1, 'Units', 'normalized', 'OuterPosition', [l b w h], ...
            'Name', 'QA of In Vivo Data', 'NumberTitle', 'off');
        % GABA DIFF spectra
        subplot(2,1,1);
        hold on;
        patch([freq(ii,ChoRange), fliplr(freq(ii,ChoRange))], ...
            [repmat(1e6, size(freq(ii,ChoRange))), fliplr(repmat(-1e6, size(freq(ii,ChoRange))))], 1, ...
            'facecolor', grey+(1-grey)*(1-shading), 'edgecolor', 'none');
        h2 = plot(freq(ii,:), GABA_diff_noalign(ii,:), 'r', freq(ii,:), GABA_diff(ii,:), 'b');
        hold off;
        set(gca, 'xdir', 'reverse', 'tickdir', 'out', 'xlim', [3.2-0.75 3.2+0.75], ...
            'ylim', [min([GABA_diff_noalign(ii,plotRange), GABA_diff(ii,plotRange)]), max([GABA_diff_noalign(ii,plotRange), GABA_diff(ii,plotRange)])]);
        text(0.5, 0.85, 'Cho SA', 'Rotation', 90, 'HorizontalAlignment', 'right', 'Units', 'normalized');
        title(sprintf('PRESS SPACE TO CONTINUE\n\nGABA'));
        legend(h2, {'PRE','POST'}, 'location', 'northeast');
        legend boxoff;
        
        % GSH DIFF spectra
        subplot(2,1,2);
        hold on;
        patch([freq(ii,ChoRange), fliplr(freq(ii,ChoRange))], ...
            [repmat(1e6, size(freq(ii,ChoRange))), fliplr(repmat(-1e6, size(freq(ii,ChoRange))))], 1, ...
            'facecolor', grey+(1-grey)*(1-shading), 'edgecolor', 'none');
        plot(freq(ii,:), GSH_diff_noalign(ii,:), 'r', freq(ii,:), GSH_diff(ii,:), 'b');
        hold off;
        title('GSH');
        set(gca, 'xdir', 'reverse', 'tickdir', 'out', 'xlim', [3.2-0.75 3.2+0.75], ...
            'ylim', [min([GSH_diff_noalign(ii,plotRange), GSH_diff(ii,plotRange)]), max([GSH_diff_noalign(ii,plotRange), GSH_diff(ii,plotRange)])]);
        
        pause;
        clf;
        
    end
    
end

out.SA.GABA.overall_quality = mean(out.SA.GABA.quality);
out.SA.GSH.overall_quality = mean(out.SA.GSH.quality);
out.SA.overall_quality = (out.SA.GABA.overall_quality + out.SA.GSH.overall_quality)/2;

close all;

