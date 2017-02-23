function Gannetplotprepostalign(MRS_struct, reg, specno)
% Plots pre and post alignment spectra in MRSLoadPfiles
% 110214:  Scale spectra by the peak _height_ of water
%          Plot multiple spectra as a stack - baselines offset
%            by mean height of GABA

ii=MRS_struct.ii; % MM

for kk = 1:length(reg)
    
    if MRS_struct.p.HERMES
        numspec = 4;
        
        SpectraToPlot = [eval(['MRS_struct.spec.', reg{kk}, sprintf('.%s',MRS_struct.p.target),'.diff(specno,:)']); ...
            eval(['MRS_struct.spec.', reg{kk}, sprintf('.%s',MRS_struct.p.target),'.diff_noalign(specno,:)']); ...
            eval(['MRS_struct.spec.', reg{kk}, sprintf('.%s',MRS_struct.p.target2),'.diff(specno,:)']); ...
            eval(['MRS_struct.spec.', reg{kk}, sprintf('.%s',MRS_struct.p.target2),'.diff_noalign(specno,:)']);];
        
        % Estimate baseline from between GABAGlx or Lac and GSH. The values might be changed depending on the future choice of metabolites -- MGSaleh
        if strcmp(MRS_struct.p.target2, 'Lac') % Defining different limits for diferent target -- MGSaleh 2016
            z=abs(MRS_struct.spec.freq(ii,:)-1.5);
            Glx_right=find(min(z)==z);
            z=abs(MRS_struct.spec.freq(ii,:)-1.0);
            GABA_left=find(min(z)==z);
            z=abs(MRS_struct.spec.freq(ii,:)-0.5);
            GABA_right=find(min(z)==z);
        else
            z=abs(MRS_struct.spec.freq(ii,:)-3.1);
            Glx_right=find(min(z)==z);
            z=abs(MRS_struct.spec.freq(ii,:)-2.9);
            GABA_left=find(min(z)==z);
            z=abs(MRS_struct.spec.freq(ii,:)-2.8);
            GABA_right=find(min(z)==z);
        end
        
        specbaseline = (mean(real(SpectraToPlot(1,Glx_right:GABA_left)),2));
        if isnan(specbaseline) % MM (for simulated data)
            [GABA_left, Glx_right] = deal(Glx_right, GABA_left);
            specbaseline = (mean(real(SpectraToPlot(1,Glx_right:GABA_left)),2));
        end
    else
        numspec = 2;
        
        %To determine the output depending on the type of acquistion used -- MGSaleh 2016
        SpectraToPlot = [eval(['MRS_struct.spec.', reg{kk}, sprintf('.%s',MRS_struct.p.target),'.diff(specno,:)']); ...
            eval(['MRS_struct.spec.', reg{kk}, sprintf('.%s',MRS_struct.p.target),'.diff_noalign(specno,:)']);];
        
        % Estimate baseline from between Glx and GABA
        z=abs(MRS_struct.spec.freq(ii,:)-3.6);
        Glx_right=find(min(z)==z);
        z=abs(MRS_struct.spec.freq(ii,:)-3.3);
        GABA_left=find(min(z)==z);
        z=abs(MRS_struct.spec.freq(ii,:)-2.8);
        GABA_right=find(min(z)==z);
        
        specbaseline = (mean(real(SpectraToPlot(1,Glx_right:GABA_left)),2));
    end
    
    
    % SpectraToPlot = [MRS_struct.spec.diff(specno,:); MRS_struct.spec.diff_noalign(specno,:)];
    
    % % Estimate baseline from between Glx and GABA
    % z=abs(MRS_struct.spec.freq(ii,:)-3.6);
    % Glx_right=find(min(z)==z);
    % z=abs(MRS_struct.spec.freq(ii,:)-3.3);
    % GABA_left=find(min(z)==z);
    % z=abs(MRS_struct.spec.freq(ii,:)-2.8);
    % GABA_right=find(min(z)==z);;
    % specbaseline = (mean(real(SpectraToPlot(1,Glx_right:GABA_left)),2));
    
    
    % Added by MGSaleh 2016
    if MRS_struct.p.HERMES
        
        % averaged gaba height across all scans - to estimate stack spacing
        gabaheight = abs(max(SpectraToPlot(1,Glx_right:GABA_right),[],2));
        if isempty(gabaheight) % MM (for simulated data)
            [GABA_right, Glx_right] = deal(Glx_right, GABA_right);
            gabaheight = abs(max(SpectraToPlot(1,Glx_right:GABA_right),[],2));
        end
        gabaheight = mean(gabaheight);
        plotstackoffset = [ 0 : (numspec-1) ]';
        
        if strcmp(MRS_struct.p.target2, 'Lac')  % Defining different limits for diferent target -- MGSaleh 2016
            plotstackoffset = plotstackoffset * 0.5 * gabaheight;
        else
            if ~strcmp(MRS_struct.gabafile{ii}((end-3):end),'.mat') % MM (for simulated data)
                plotstackoffset = plotstackoffset * 1.75 * gabaheight; %RE % 3.0
            else
                plotstackoffset = plotstackoffset * 1.3 * gabaheight; %RE % 3.0
            end
        end
        plotstackoffset = plotstackoffset - specbaseline;
        
        aa=1.2;
        plot(MRS_struct.spec.freq(ii,:), aa*real(SpectraToPlot((1),:)),'b',MRS_struct.spec.freq(ii,:), aa*real(SpectraToPlot((2),:)),'r');
        hold on
        shft=repmat(plotstackoffset, [1 length(SpectraToPlot(1,:))]);
        SpectraToPlot(3:4,:) = SpectraToPlot(3:4,:) + [max(shft,[],1); max(shft,[],1)] ;
        plot(MRS_struct.spec.freq(ii,:), aa*real(SpectraToPlot((3),:)),'b',MRS_struct.spec.freq(ii,:), aa*real(SpectraToPlot((4),:)),'r');
        hold off
        
        % yaxismax
        % yaxismin
        if strcmp(MRS_struct.p.target2, 'Lac')
            yaxismax = (numspec + 1.0) * 0.5 * gabaheight; % top spec + 0.5 * height of gaba %Changed slightly by MGSaleh to accomodate both GSH and GSH/Lac -- 2016
        else
            %yaxismax = (numspec + 1.0) * 3.0 * gabaheight; % top spec + 1.0 * height of gaba %Changed slightly by MGSaleh to accomodate both GSH and GABAGlx/GSH -- 2016
            if ~strcmp(MRS_struct.gabafile{ii}((end-3):end),'.mat') % MM (for simulated data)
                yaxismax = (numspec + 1.0) * 1.95 * gabaheight;
            else
                yaxismax = (numspec + 1.0) * 1.35 * gabaheight;
            end
        end
        yaxismin =  -2.0* gabaheight; % extend 2* gaba heights below zero %Changed slightly by MGSaleh to accomodate both GSH and GABAGlx/Lac -- 2016
        
        if (yaxismax<yaxismin)
            dummy=yaxismin;
            yaxismin=yaxismax;
            yaxismax=dummy;
        end
    else
        % averaged gaba height across all scans - to estimate stack spacing
        %gabaheight = abs(max(SpectraToPlot(1,Glx_right:GABA_right),[],2));
        gabaheight = abs(max(SpectraToPlot([1 2],Glx_right:GABA_right),[],2)); % Changed for better view -- MGSaleh
        %gabaheight = mean(gabaheight);
        gabaheight = max(gabaheight); % Changed for better view -- MGSaleh
        plotstackoffset = [ 0 : (numspec-1) ]';
        plotstackoffset = plotstackoffset * gabaheight;
        plotstackoffset = plotstackoffset - specbaseline;
        
        SpectraToPlot = SpectraToPlot + ...
            repmat(plotstackoffset, [ 1  length(SpectraToPlot(1,:))]);
        
        plot(MRS_struct.spec.freq(ii,:), real(SpectraToPlot(1,:)),'b',MRS_struct.spec.freq(ii,:), real(SpectraToPlot(2,:)),'r');
        % yaxismax
        % yaxismin
        
        yaxismax = 1.5*abs(max(max(real(SpectraToPlot([1 2],Glx_right:GABA_right)),[],2))); % I removed the axis code below to zoom into the data -- MGSaleh 2016
        %yaxismax = (numspec + 1.0) *gabaheight; % top spec + 2* height of gaba %Changed slightly by MGSaleh to accomodate both GSH and GABAGlx/Lac -- 2016
        yaxismin = -10.0*abs(min(min(real(SpectraToPlot([1 2],Glx_right:GABA_right)),[],2))); % I removed the axis code below to zoom into the data -- MGSaleh 2016
        %yaxismin = - 2.0* gabaheight; % extend 2* gaba heights below zero %Changed slightly by MGSaleh to accomodate both GSH and GABAGlx/Lac -- 2016
        if (yaxismax<yaxismin)
            dummy=yaxismin;
            yaxismin=yaxismax;
            yaxismax=dummy;
        end
    end
    
    legendtxt = {'post', 'pre'};
    hl=legend(legendtxt);
    set(hl,'EdgeColor',[1 1 1]);
    set(gca,'XDir','reverse');
    oldaxis = axis;
    
    axis([0 5  yaxismin yaxismax])
    
end
