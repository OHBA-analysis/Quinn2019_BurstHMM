%% Overview
%
% This script will create a simple simulation of 1/f noise with additive bursts
% at two different frequencies. The simulated time-courses are filtered to one
% of two frequencies bands and the amplitude envelope of each band is computed.
% A two state HMM is inferred on each of these to identify periods of high
% average amplitude envelope within the signal.
%
% A summary figure is created and saved out into the config.figpath directory

%% Set up and check paths.

config = hmm_0_initialise;

%% Create simulated signal

[data,x,time_vect,sample_rate] = hmm_util_get_simulation;

%% HMM estimation

% The simulated data is first filtered to 20-30Hz and 20-45Hz
[B,A] = butter( 5, [20 30]./(sample_rate/2));
data(2,:) = filtfilt(B,A,data(1,:)')';

[B,A] = butter( 5, [20 45]./(sample_rate/2));
data(3,:) = filtfilt(B,A,data(1,:)')';

% The amplitude envelope is computed using the hilbert transform
envelopes = abs(hilbert(data')');

% Here we generate an amplitude threshold using the Shin 2018 method
threshold = 2*median(envelopes,2);

% HMM - Main loop
Gamma_emb = cell(3,1);
for ii = 2:3

    % HMM options
    % more details can be found here:
    % https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#estimation

    options = struct();
    options.initrep = 10; % Set to be large as this is a short simulation
    options.K = 2;
    options.standardise = 0;
    options.verbose = 1;
    options.Fs = sample_rate;
    options.useMEX = 1;
    options.zeromean = 0;
    options.dropstates = 1;
    options.order = 0;
    options.DirichletDiag = 1e10; % set very large as we don't have much data
    
    % HMM inference - we only store the Gamma time-course of posterior probabilities
    T = length(data);
    [hmm, Gamma_emb{ii},~,vpath] = hmmmar(envelopes(ii,:)',T,options);

end


%% Create a summary figure

% Colours - http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
cols = [228,26,28; 128,128,128] ./ 255;

timelims = [7 13];
time_vect2 = time_vect - timelims(1);
font_size = 12;

figure('Position',[100 100 1024 768])
ax1 = axes('Position',[.075 .85 .725 .15],'Visible',false);hold on; grid on
plot(time_vect,data(1,:).*2,'k')
plot(time_vect,x(:,2)-4,'Color',[27,158,119]./255)
plot(time_vect,x(:,3)-6,'Color',[ 117,112,179]./255)
%legend({'Data','High Burst','Low Burst'})
title('Time Series','FontSize',font_size)
axis('tight');
xlim(timelims);
text(timelims(1)-.5,0,'Data','FontSize',font_size,'FontWeight','bold')
text(timelims(1)-.6,-3.8,'Burst Type 1','FontSize',font_size,'FontWeight','bold')
text(timelims(1)-.6,-5.9,'Burst Type 2','FontSize',font_size,'FontWeight','bold')


ybase = [.475 .05];
idx = [2 4];
titles = {'20-30Hz Amplitude Envelope',...
          '20-45Hz Amplitude Envelope'};

for ii = 1:2
    
    ax2 = axes('Position',[.075 ybase(ii) .725 .14]);hold on
    h = area(ax2,time_vect,Gamma_emb{ii+1});
    xlim(ax2,timelims);
    %xlabel('Time (seconds','FontSize',font_size);
    ylabel({'HMM-State','Probability'},'FontSize',font_size);
    ax2.YAxis.FontSize = font_size;
    ax2.XAxis.FontSize = font_size;
    xticklabels(ax2,1:7)
    xlabel('Time (seconds)')

    ax3 = axes('Position',[.075 ybase(ii)+.225 .725 .1],'Visible',false);hold on
    ax3.YAxis.Visible = 'on';
    plot(ax3,time_vect,data(ii+1,:),'k')
    plot(ax3,time_vect,envelopes(ii+1,:),'r')
    xlim(ax3,timelims);
    plot(ax3,[time_vect(1) time_vect(end)],[threshold(ii+1) threshold(ii+1)],'b--')
    title(titles{ii},'Visible',true,'FontSize',font_size+2)
    ylabel({'Filtered','Data'},'FontSize',font_size);

    ax4 = axes('Position',[.075 ybase(ii)+.15 .725 .075],'Visible',false);hold on
    area(ax4, time_vect,envelopes(ii+1,:)>threshold(ii+1))
    xlim(ax4,timelims);
    ylim([-.25 1.25])
    ylabel(ax4,{'Amplitude','Threshold',''},'FontSize',font_size);
    ax4.YAxis.Visible = 'on';
    yticks([]);

    ax6 = axes('Position',[.845 ybase(ii)+.165 .15 .125],'defaultAxesFontSize',font_size);hold on
    histogram(ax6, envelopes(ii+1,:), linspace(0,1),...
        'DisplayStyle','stairs','EdgeColor','k',...
                'linewidth',2)
    inds = envelopes(ii+1,:)>threshold(ii+1);
    histogram(ax6, envelopes(ii+1,inds), linspace(0,1),'EdgeColor','none')
    yl = ylim;
    plot([threshold(ii+1) threshold(ii+1)],yl,'b--')
    ylabel('Num Occurrences','FontSize',font_size)

    
    ax7 = axes('Position',[.845 ybase(ii) .15 .125],'defaultAxesFontSize',font_size);hold on
    for jj = 1:2
        inds = Gamma_emb{ii+1}(:,jj) > threshold(ii+1);
        env = envelopes(ii+1,:)';
        if mean(env(inds)) > .4
            histogram(ax7, env(inds), linspace(0,1),...
                'DisplayStyle','stairs','EdgeColor',cols(1,:),...
                'linewidth',2)
            h(jj).FaceColor = cols(1,:);
        else
            histogram(ax7, env(inds), linspace(0,1),...
                'DisplayStyle','stairs','EdgeColor',cols(2,:),...
                'linewidth',2)
            h(jj).FaceColor = cols(2,:);
        end
    end
    xlabel('Amplitude','FontSize',font_size)
    ylabel('Num Occurrences','FontSize',font_size)
    
end

% save figure
figpath = fullfile(config.figpath,'hmm_2_envelope_figure');
saveas(gcf,figpath,'png') 
saveas(gcf,figpath,'tiff') 