%% Overview
%
% This script loads some data from the Wakeman & Henson (2015) dataset
% preprocessed with the methods from Quinn et al 2018.
%
% A single parcel time-course from motor cortex is downsampled to 100Hz and
% bad segments are removed from the data before the time-course is modelled
% by a single time-delay embedded HMM with 6 states.
%
% Finally, the task-evoked state time-courses are computed and a
% HMM-regularised image of the task response created.
%
% summary figures are created and saved out into the config.figpath directory
%
% References
%
% Quinn, A. J., Vidaurre, D., Abeysuriya, R., Becker, R., Nobre, A. C., &
% Woolrich, M. W. (2018). Task-Evoked Dynamic Network Analysis Through Hidden
% Markov Modeling. Frontiers in Neuroscience, 12.
% https://doi.org/10.3389/fnins.2018.00603
%
% Wakeman, D. G., and Henson, R. N. (2015). A multi-subject, multi-modal human
% neuroimaging dataset. Sci. Data 2:150001. doi: 10.1038/sdata.2015.1

%% Set up and check paths.

config = hmm_0_initialise;

%% Load in data
% The data file contains...
%
% data_trl - [ntrials x nsamples] matrix of data
% data_hmm - [ntrials*nsamples x 1] vector of data (reshaped from data_trl)
% T        - [ntrials x t] vector indicating the length of each trial
% wavelet_evoked - [nfreqs x nsamples] matrix of wavelet transformed data
% wavelet_freq - [nfreqs x 1] vector of frequencies

load( fullfile(config.basepath,'hmm_WakemanHenson_data'))


%% HMM Inference and post-stats

% HMM Inference - this takes ~15 minutes on a fast computer

tic % set up a timer 

% HMM options
% more details can be found here:
% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#estimationoptions = struct();
options = [];
options.initrep = 3;
options.initcyc = 1e3;
options.K = 6;
options.standardise = 0;
options.verbose = 1;
options.Fs = sample_rate;
options.useMEX = 1;
options.zeromean = 0;
options.dropstates = 1;
options.useParallel = 0; % set to 1 if you have access to the parallel processing toolbox

options.order = 0;
options.embeddedlags = -11:11;
options.covtype = 'full';

[hmm, gamma_emb,~,vpath,GammaInit] = hmmmar(data_hmm,T,options);

toc % Stop timer and print elapsed time

% Compute state-wise multi-taper power spectrum
options.win = 128;
options.verbose = 0; % Don't need the output here
spec_emb = hmmspectramt(data_hmm,T,gamma_emb,options);

% Compute state life-times and covert from samples to milliseconds
lifetimes = getStateLifeTimes(vpath,T,options);

% Compute state interval-times and covert from samples to milliseconds
intervaltimes = getStateIntervalTimes(vpath,T,options);

for ii = 2:length(T)
    for jj = 1:options.K
        lifetimes{1,jj} = cat(2,lifetimes{1,jj},lifetimes{ii,jj});
        intervaltimes{1,jj} = cat(2,intervaltimes{1,jj},intervaltimes{ii,jj});
    end
end
lifetimes = {lifetimes{1,:}};
intervaltimes = {intervaltimes{1,:}};

lifetimes_ms = cellfun(@(x) x*(1/sample_rate),lifetimes,'UniformOutput',false);
intervaltimes_ms = cellfun(@(x) x*(1/sample_rate),intervaltimes,'UniformOutput',false);

% Compute fractional occupancy
fractional_occupancy = getFractionalOccupancy(vpath,size(vpath,1),options);


% % correct Gamma time-course for lags in  model
gamma_emb = padGamma(gamma_emb,T,options);

gamma_emb_trl = reshape(gamma_emb,size(data_trl,2),size(data_trl,1),options.K);
evoked_time_vect = linspace(-1,2.33,351)';

%% Summary figure

% Colours from:
% http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9
cols = {[228,26,28],...
        [55,126,184],...
        [77,175,74],...
        [152,78,163],...
        [255,127,1],...
        [166,86,40],...
        [247,129,191]};
    
font_size = 14;

psd = cell2mat({spec_emb.state(:).psd});
[~,state_inds] = sort(sum(psd,1),'descend');

cols2 = cell(5,1);
for ii = 1:options.K
    cols2{state_inds(ii)} = cols{ii} / 255;
end

%
figure('Position',[100 100 1280 1024])

trial_inds = [25,330,200];

for ii = 1:2
    ax = subplot(6,3,ii);
    ax.Position(2) = ax.Position(2) + .05;
    plot( evoked_time_vect,data_trl(trial_inds(ii),:),'k','linewidth',.5 )
    axis('tight')
    %plot(time_vect(inds),data_pad(inds),'k','linewidth',.5);axis('tight')
    %xlim([time_vect(inds(1)) time_vect(inds(end))])
    set(ax,'FontSize',font_size)
    
    ax = subplot(6,3,[4,7]+(ii-1));
    ax.Position(2) = ax.Position(2) + .1;
    ax.Position(4) = ax.Position(4) -.05;
    %contourf(time_vect(inds),wf,abs(wt(:,inds)).^.75,'linestyle','none')
    contourf( evoked_time_vect,wavelet_freq,wavelet_trl(:,:,trial_inds(ii)),'linestyle','none')
    ylabel('Frequency (Hz)')
    title('Wavelet Transform')
    set(ax,'FontSize',font_size)
    
    if ii == 1
        cb = colorbar;
        cb.Position = cb.Position + 1e-10;
        cb.Position(1) = .63;
        cb.Label.String = 'Amplitude';
        cb.Label.FontSize = 14;
    end
    
    ax = subplot(6,3,10+(ii-1));
    ax.Position(2) = ax.Position(2) + .1;
    h = area(evoked_time_vect,squeeze(gamma_emb_trl(:,trial_inds(ii),:)));
    xlim([evoked_time_vect(1) evoked_time_vect(end)])
    axis('tight');grid on
    for jj = 1:options.K
        text(81.5,1-(.25*(jj-1)),['State ' num2str(jj)],'FontSize',14,...
            'Color',cols2{jj},'FontWeight','bold')
    end
    xlabel('Time (seconds)')
    ylabel('Probability');
    title('State Posterior Probability')
    set(ax,'FontSize',font_size)
    for jj = 1:options.K
        h(jj).FaceColor = cols2{jj};
    end
end

ax = subplot(6,3,[3,6,9]);
ax.Position(1) = ax.Position(1) + .06;
hold on;grid on
for ii = 1:options.K
    plot(spec_emb.state(ii).f,abs(spec_emb.state(ii).psd),...
            'linewidth',3,'Color',cols2{ii})
    %h(ii).FaceColor = cols2{ii};
end
%set(ax,'YScale','log');
title('State Spectra','FontSize',font_size)
xlabel('Frequency (Hz)')
ylabel('Power')
set(ax,'FontSize',font_size)

ax = subplot(6,3,13);
ax.Position(4) = ax.Position(4) + .08;
ax.Position(2) = ax.Position(2) - .02;
distributionPlot(lifetimes_ms,'color',cols2)
yticks(0:.25:1)
grid on;ylim([0 1])
title('State Lifetimes','FontSize',font_size)
xlabel('State')
ylabel('Lifetime (ms)')
set(ax,'FontSize',font_size)

ax = subplot(6,3,14);
ax.Position(4) = ax.Position(4) + .08;
ax.Position(2) = ax.Position(2) - .02;
distributionPlot(intervaltimes_ms,'color',cols2)
yticks(0:1:6);ylim([0 3.5])
grid on
title('State Interval-times','FontSize',font_size)
xlabel('State')
ylabel('Interval Time (s)')
set(ax,'FontSize',font_size)

ax = subplot(6,3,15);hold on
ax.Position(4) = ax.Position(4) + .08;
ax.Position(2) = ax.Position(2) - .02;
for ii = 1:options.K
    bar(ii,fractional_occupancy(ii),'FaceColor',cols2{ii});
end
xticks(1:options.K);grid on
title('State Fractional Occupancy','FontSize',font_size)
xlabel('State')
ylabel('Fractional Occupancy')
set(ax,'FontSize',font_size)

for ii = 1:options.K
    ax = subplot(6,6,30+ii);
    imagesc(hmm.state(ii).W.S_W)
    ax.Position(2) = ax.Position(2) - .05;
    set(ax,'XTick',[1,13,24],'XTickLabel',[-12 0 12],'FontSize',font_size)
    set(ax,'YTick',[1,13,24],'YTickLabel',[-12 0 12])
    title(['State ' num2str(ii)])
end

% save figures
figpath = fullfile(config.figpath,'hmm_4_realdata_figure');
saveas(gcf,figpath,'png') 
saveas(gcf,figpath,'tiff') 


%% Task-evoked state analysis


gamma_evoked = squeeze(nanmean(gamma_emb_trl,2));
gamma_baseline = squeeze(mean(gamma_evoked(5:95,:),1));

hmm_reg_evoked = (gamma_evoked-gamma_baseline) * (psd.*fractional_occupancy)';

wavelet_evoked = nanmean(wavelet_trl,3);
wavelet_baseline = mean(wavelet_evoked(:,25:75),2);

% Create a summary figure

figure('Position',[100 100 1024 341])
ax = subplot(131);hold on; grid on
ax.Position(1) = ax.Position(1) - .05;
ax.Position(3) = ax.Position(3) + .05;
ax.Position(4) = ax.Position(4) - .1;
contourf(evoked_time_vect,wavelet_freq, wavelet_evoked-wavelet_baseline,48,'linestyle','none' );
cl = [min(min(wavelet_evoked-wavelet_baseline)),max(max(wavelet_evoked-wavelet_baseline))];
hmm_util_redblue_colourmap(gca,cl(1),cl(2));
ylim([0 35])
colorbar
xlabel('Time (seconds)')
ylabel('Frequency (Hz')
title({'Wavelet','TF Transform','',''})
plot([0 0],ylim,'k--')
text(0,37,{'Stimulus','onset'},'HorizontalAlignment','center','FontSize',9)
plot([.932 .932],ylim,'k--')
text(.932,37,{'Average','Reaction Time'},'HorizontalAlignment','center','FontSize',9)

ax = subplot(132);grid on; hold on;
ax.Position(3) = ax.Position(3) + .05;
ax.Position(4) = ax.Position(4) - .1;
contourf(evoked_time_vect,spec_emb.state(1).f, hmm_reg_evoked',48,'linestyle','none' );
cl = [min(min(hmm_reg_evoked)),max(max(hmm_reg_evoked))];
hmm_util_redblue_colourmap(gca,cl(1),cl(2));
ylim([0 35])
colorbar
xlabel('Time (seconds)')
ylabel('Frequency (Hz')
title({'HMM-Regularised','TF Transform','',''})
plot([0 0],ylim,'k--')
text(0,37,{'Stimulus','onset'},'HorizontalAlignment','center')
plot([.932 .932],ylim,'k--')
text(.932,37,{'Average','Reaction Time'},'HorizontalAlignment','center')

ax = subplot(133);hold on;grid on
ax.Position(1) = ax.Position(1) + .05;
ax.Position(4) = ax.Position(4) - .1;
for ii = 1:options.K
    plot(evoked_time_vect,gamma_evoked(:,ii) - gamma_baseline(ii),...
            'Color',cols2{ii},'linewidth',2)
end
l = legend({'State 1','State 2','State 3','State 4','State 5','State 6'})
l.Position(1) = l.Position(1) +.035;

xlabel('Time (seconds')
ylabel({'Relativel Occupancy','Baseline-corrected'})
title({'Task-evoked','Fractional-Occupancy','',''});
plot([0 0],ylim,'k--','HandleVisibility','off')
text(0,.165,{'Stimulus','onset'},'HorizontalAlignment','center')
plot([.932 .932],ylim,'k--','HandleVisibility','off')
text(.932,.165,{'Average','Reaction Time'},'HorizontalAlignment','center')

% save figures
figpath = fullfile(config.figpath,'hmm_4_realdata_spectrum_figure');
saveas(gcf,figpath,'png') 
saveas(gcf,figpath,'tiff') 

%
figure('Position',[100 100 1024 512])
for k = 1:options.K
    ax = subplot(3,6,k);
    x = gamma_evoked(:,k) * psd(:,k)';
    contourf(ax,evoked_time_vect,spec_emb.state(k).f,x', 48, 'linestyle','none' )
    title(['State ' num2str(k)]);
    if k == 1
        ylabel('Frequency (Hz)')
    end
    
    ax = subplot(3,6,k+6);
    ax.Position(2) = ax.Position(2) + .15;
    ax.Position(4) = .1;
    plot(ax,evoked_time_vect,gamma_evoked(:,k),'Color',cols2{k},'linewidth',2);
    grid on
    xlim(ax,[evoked_time_vect(1) evoked_time_vect(end)]);
    if k == 1
        ylabel({'Fractional','Occupancy'})
    end
    
    ax = subplot(2,6,k+6);
    ax.Position(4) = ax.Position(4)+.1;
    x = squeeze(gamma_emb_trl(:,:,k) < .5);
    imagesc(ax,evoked_time_vect,1:size(gamma_emb_trl,2),x');
    colormap(ax,'bone');
    
    if k == 1
        ylabel('Trial')
    end
    xlabel('Time (Seconds)')
end

% save figures
figpath = fullfile(config.figpath,'hmm_4_realdata_statespectrum_figure');
saveas(gcf,figpath,'png') 
saveas(gcf,figpath,'tiff') 

%% Compute evoked temporal statistics from state lifetimes

embedding_lags = 11;

lifetimes_trl = getStateLifeTimes(vpath,T,options);
lifetimes_trl_ms = cellfun(@(x) x*(1/sample_rate),lifetimes_trl,'UniformOutput',false);

lt_trl = zeros(351,6,878);
gam = gamma_emb_trl(embedding_lags:end-embedding_lags,:,:);

vp = padVP(vpath,T,options);
vp = reshape(vp,351,878);

lives = zeros( 878,351,6 );

for ii = 1:878
    for jj = 1:6
        
        state_starts = find(diff(vp(:,ii)==jj) == 1);
        state_stops = find(diff(vp(:,ii)==jj) == -1);
        
        for kk = 1:length(state_starts)
            if state_starts(kk)<=embedding_lags*2 || state_stops(kk) >= 351-embedding_lags*2
                continue
            end
            lives( ii,state_starts(kk):state_stops(kk),jj ) = lifetimes_trl_ms{ii,jj}(kk);
        end
    end
end
lives(lives==0) = nan;

%% Summary figure

ind = 2; % This needs to be specified to match the rebound state in your HMM
figure('Position',[100 100 768 1024]);
ax = subplot(321);hold on
ax.Position(1) = ax.Position(1) - .025;
plot( evoked_time_vect, gamma_evoked(:,ind),'k','linewidth',2,'color',cols2{ind})
xlabel('Time (seconds)')
ylabel('Fractional Occupancy')
title({'Task-evoked','Fractional-Occupancy','',''})

plot([0 0],ylim,'k--')
text(0,.275,{'Stimulus','onset'},'HorizontalAlignment','center')
plot([.932 .932],ylim,'k--')
text(.932,.275,{'Average','Reaction Time'},'HorizontalAlignment','center')
set(ax,'FontSize',font_size)

xlim([-1,2.33])
grid on
set(gca,'FontSize',font_size)

ax = subplot(3,2,[3,5]);hold on
ax.Position(1) = ax.Position(1) - .025;
ax.Position(2) = .175;
ax.Position(4) = .375;
x = gamma_evoked(:,ind) * psd(:,ind)';
ctf = contourf(evoked_time_vect,spec_emb.state(ind).f,x', 48, 'linestyle','none' );
title({'HMM Regularised','Spectrum','',''},'FontSize',font_size)
xlabel('Time (seconds)')
ylabel('Frequency (Hz)')
cb = colorbar('Position',[.1 .1 .33 .015],'Location','South');
cb.Label.String = 'Amplitude';

set(gca,'FontSize',font_size)

plot([0 0],ylim,'k--')
text(0,52,{'Stimulus','onset'},'HorizontalAlignment','center')
plot([.932 .932],ylim,'k--')
text(.932,52,{'Average','Reaction Time'},'HorizontalAlignment','center')
set(ax,'FontSize',font_size)

ax = subplot(122);hold on
ax.Position(1) = ax.Position(1) - .0325;
plot(evoked_time_vect,lives(:,:,ind)','linewidth',1)
plot(evoked_time_vect,nanmedian(lives(:,:,ind),1)','k','linewidth',2)
xlabel('Time (seconds)')
ylabel('Lifetime (seconds)')
xlim([-1,2.33])

title({'Task-Evoked','State Lifetimes','',''},'FontSize',font_size)

plot([0 0],ylim,'k--')
text(0,1.85,{'Stimulus','onset'},'HorizontalAlignment','center')
plot([.932 .932],ylim,'k--')
text(.932,1.85,{'Average','Reaction Time'},'HorizontalAlignment','center')
set(ax,'FontSize',font_size)

axes('Position',[.91 ax.Position(2) .075, ax.Position(4)])
distributionPlot({lifetimes_ms{ind}},'color',cols2{ind},'histOpt',1)
grid on;
ylim([0 1.8])
xlabel('State')
xticklabels(ind)
set(gca,'FontSize',font_size)

% save figure
figpath = fullfile(config.figpath,'hmm_4_realdata_evokedlifetime_figure');
saveas(gcf,figpath,'png') 
saveas(gcf,figpath,'tiff') 
    
