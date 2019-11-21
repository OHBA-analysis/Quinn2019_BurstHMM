%% Overview
%
% This script will create a simple simulation of 1/f noise with additive bursts
% at two different frequencies. The simulated time-course is modelled with a
% single time-delay embedded HMM using three states.

% A summary figure is created and saved out into the config.figpath directory

%% Set up and check paths.

config = hmm_0_initialise;

%% Create simulated signal

% the TDE-HMM will infer more parameters than the AE-HMM, so this example
% uses a longer version of the simulation.
[data,x,time_vect,sample_rate]  = hmm_util_get_simulation('long');

%% Wavelet transform

[wt,wf] = cwt(data,'amor',sample_rate);

%% HMM estimation
% The HMM is esimated on the time-domain data so no further preprocessing is
% required.

% HMM options
% more details can be found here:
% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#estimation

options = struct();
options.repetitions = 4;
options.initrep = 5;
options.K = 3;
options.standardise = 0;
options.verbose = 1;
options.Fs = sample_rate;
options.useMEX = 1;
options.zeromean = 0;
options.dropstates = 1;
options.order = 0;
options.embeddedlags = -11:11;
options.covtype = 'full';

% The inference itself
T = length(data);
[hmm, Gamma_emb,~,vpath] = hmmmar(data,T,options);

% correct Gamma time-course for lags in  model
Gamma_emb = padGamma(Gamma_emb,T,options);
options.win = 1024;

% Compute post-hoc spectra for each state
fit_emb = hmmspectramt(data,T,Gamma_emb,options);

%% Create a summary figure

% Colours - http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
cols = [27,158,119;117,112,179;200,200,200] ./ 255;

timelims = [7 13];
time_inds = find((time_vect>timelims(1)) .* (time_vect<timelims(2)));
time_vect2 = time_vect - timelims(1);
font_size = 12;

figure('Position',[100 100 1024 768])
ax1 = axes('Position',[.075 .84 .825 .15],'Visible',false);hold on; grid on
plot(time_vect,data.*2,'k')
plot(time_vect,x(:,2)-4,'Color',cols(1,:))
plot(time_vect,x(:,3)-6,'Color',cols(2,:))
title('Time Series','FontSize',font_size)
axis('tight');
xlim(timelims);
text(timelims(1)-.5,0,'Data','FontSize',font_size,'FontWeight','bold')
text(timelims(1)-.54,-3.8,'Burst Type 1','FontSize',font_size,'FontWeight','bold')
text(timelims(1)-.54,-5.9,'Burst Type 2','FontSize',font_size,'FontWeight','bold')

ax2 = axes('Position',[.075 .55 .825 .25]);
contourf(time_vect(time_inds),flipud(wf),flipud(abs(wt(:,time_inds))),36,'linestyle','none')
xlim(timelims);
cb = colorbar;
cb.Position = cb.Position + 1e-10;
cb.Position(1) = .925;
cb.Label.String = 'Amplitude';
cb.Label.FontSize = 14;

ax2.YAxis.FontSize = font_size;
ax2.XAxis.FontSize = font_size;
set(ax2,'XMinorTick','on','YTick',0:10:90)

title('Wavelet Transform','FontSize',font_size+2)
ylabel('Frequency (Hz)','FontSize',font_size)

% match cols
psd = cell2mat({fit_emb.state(:).psd});
state_ind = [0,0,0];
[~,state_ind(1)] = max(psd(95,:));
[~,state_ind(2)] = max(psd(155,:));
state_ind(3) = setdiff(1:3,state_ind);

ax3 = axes('Position',[0.0750 0.3650 0.8250 0.1243]);
h = area(time_vect,Gamma_emb);
for ii = 1:3
    h(ii).FaceColor = cols(state_ind(ii),:);
end
axis('tight')
title('State-wise Posterior-Probabilities','FontSize',font_size+2)
xlim(timelims);

ax3.YAxis.FontSize = font_size;
ax3.XAxis.FontSize = font_size;
xlabel('Time (seconds)','FontSize',font_size)
ylabel('Probability','FontSize',font_size)
set(ax3,'XMinorTick','on')
text(timelims(2)+.15,.75,'State 1','FontSize',font_size,...
        'Color',cols(state_ind(1),:),'FontWeight','bold')
text(timelims(2)+.15,.5,'State 2','FontSize',font_size,...
        'Color',cols(state_ind(2),:)-.2,'FontWeight','bold')
text(timelims(2)+.15,.25,'State 3','FontSize',font_size,...
        'Color',cols(state_ind(3),:),'FontWeight','bold')
%
ax4 = axes('Position',[.075 .075 .175 .2088]);%subplot(5,4,17);
hold on;
for ii = 1:options.K
    plot(fit_emb.state(ii).f,100*abs(fit_emb.state(ii).psd),...
            'linewidth',3,'Color',cols(state_ind(ii),:))
end
title('State Spectra (MT)','FontSize',font_size+2)
ax4.YAxis.FontSize = font_size;
ax4.XAxis.FontSize = font_size;
xlabel('Frequency (Hz)');
ylabel({'Power Spectral';'Density (au)'},'FontSize',font_size);
grid on
xlim([0 75])

ax = subplot(5,4,18);
imagesc(options.embeddedlags,options.embeddedlags,flipud(hmm.state(1).W.S_W))
ax.Position(2) = .075;
ax.Position(4) = .2088;
set(ax,'XTick',[-11  0 11],'XTickLabel',[-11  0 11],'YTick',[-11  0 11],'YTickLabel',[11  0 -11])
xlabel('Lag');ylabel('Lag')
title({'State 1 Autocovariance'},'FontSize',font_size+2)
ax.YAxis.FontSize = font_size;
ax.XAxis.FontSize = font_size;

ax = subplot(5,4,19);
imagesc(options.embeddedlags,options.embeddedlags,flipud(hmm.state(2).W.S_W))
ax.Position(2) = .075;
ax.Position(4) = .2088;
set(ax,'XTick',[-11  0 11],'XTickLabel',[-11  0 11],'YTick',[-11  0 11],'YTickLabel',[11  0 -11])
xlabel('Lag');ylabel('Lag')
title({'State 2 Autocovariance'},'FontSize',font_size+2)
ax.YAxis.FontSize = font_size;
ax.XAxis.FontSize = font_size;

ax = subplot(5,4,20);
imagesc(options.embeddedlags,options.embeddedlags,flipud(hmm.state(3).W.S_W))
ax.Position(2) = .075;
ax.Position(4) = .2088;
set(ax,'XTick',[-11  0 11],'XTickLabel',[-11  0 11],'YTick',[-11  0 11],'YTickLabel',[11  0 -11])
xlabel('Lag');ylabel('Lag');
title({'State 3 Autocovariance'},'FontSize',font_size+2)
ax.YAxis.FontSize = font_size;
ax.XAxis.FontSize = font_size;
axis square

% save figure
figpath = fullfile(config.figpath,'hmm_fig3_embedded_simu');
saveas(gcf,figpath,'png')
saveas(gcf,figpath,'tiff')
