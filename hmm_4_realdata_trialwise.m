config = hmm_0_initialise;

%% Load dataset

load( fullfile(config.basepath,'hmm_Nowak2016_data.mat') );

%% Prepare data for HMM and compute wavelet transform.

sample_rate = 500;
ds_factor = 6;
dsample_rate = sample_rate/ds_factor;

data = zeros(505008,1);
T = [];
v = [];
wt = zeros(55,334,1512);
subj = [];
mt = [];

trl_ind = 1;
time_ind = 1;
for ii = 1:33

    % Find indices for good/bad trials - including some special cases
    % missed by the automatic algorithm
    [~,I] = rmoutliers(sum(motor_data{ii}.^2),'gesd');
    if ii == 19; I(47) = 1; end % special case
    if ii == 29; I(48) = 1; end % special case
    if ii == 31; I(48) = 1; end % special case
    goods = find(~I);

    % Downsamaple data
    dat = resample(motor_data{ii}(:,goods),1,ds_factor);

    % Concatenate z-transformed data into full dataset
    new_time_ind = time_ind-1+size(dat,1)*size(goods,2);
    data(time_ind:new_time_ind,1) = zscore(dat(:));

    % Concatenate tracked variables
    T = [ T; ones(size(dat,2),1)*size(dat,1)];
    subj = [ subj; ones(size(dat,2),1)*ii ];
    mt = [mt; MT_all{ii}(1,goods)'];
    v = [ v; sum(motor_data{ii}.^2)'];

    % Compute wavelet transform
    new_trl_ind = trl_ind + length(goods);
    for jj = 1:length(goods)
        wt(:,:,trl_ind+jj-1) = abs(cwt(dat(:,jj)));
    end

    % Update indices
    trl_ind = new_trl_ind;
    time_ind = new_time_ind;
end

% Get frequency vector for wavelet transform
[~,freq_vect_cwt] = cwt(dat(:,jj),dsample_rate);
f_inds_cwt = find((freq_vect_cwt > 15) .* (freq_vect_cwt < 30));

% Create time vectors and baselines
time_vect = linspace(-2,2,T(1));
bl_inds = find((time_vect > -1.5) .* (time_vect < -1));
bl = nanmean(wt(:,bl_inds,:),[2,3]);


%% HMM Inference
%
% This loop runs the HMM inference and saves the output.
% Details for the useage of hmmmar and its options can be found on the
% hmmmar wiki page (https://github.com/OHBA-analysis/HMM-MAR/wiki)

for K = (6)
    for emb = (11)

        tic % Start timer

        % Set options.
        options = struct();
        options.initrep = 6;
        options.K = K;
        options.standardise = 1;
        options.verbose = 1;
        options.Fs = dsample_rate;
        options.useMEX = 1;
        options.zeromean = 0;
        options.dropstates = 1;
        options.useParallel = 0;

        options.order = 0;
        options.embeddedlags = -emb:emb;
        options.covtype = 'full';

        % Infer HMM
        [hmm, Gamma_emb,~,vpath] = hmmmar(data,T,options);

        % Save inferred time-courses and observation models
        nm = ['beta_burst_ers_K' num2str(K) '_emb' num2str(emb)];
        save(nm,'hmm','Gamma_emb','vpath','options');

        toc % Print time elapased since tic

    end
end

%% Compute HMM Metrics

% Load previously run HMM into workspace
K = 6;
emb = 11;
nm = ['beta_burst_ers_K' num2str(K) '_emb' num2str(emb)];
load(nm)

%-------------------------
% Compute temporal metrics

% Compute state life-times
lifetimes = getStateLifeTimes(vpath,T,options);

% Compute state interval-times
intervaltimes = getStateIntervalTimes(vpath,T,options);

% Unstack lifetimes and interval times (these are returned in
for ii = 2:length(T)
    for jj = 1:options.K
        lifetimes{1,jj} = cat(2,lifetimes{1,jj},lifetimes{ii,jj});
        intervaltimes{1,jj} = cat(2,intervaltimes{1,jj},intervaltimes{ii,jj});
    end
end
lifetimes = {lifetimes{1,:}};
intervaltimes = {intervaltimes{1,:}};

% Convert lifetimes and interval times to milliseconds.
lifetimes_ms = cellfun(@(x) x*(1/dsample_rate),lifetimes,'UniformOutput',false);
intervaltimes_ms = cellfun(@(x) x*(1/dsample_rate),intervaltimes,'UniformOutput',false);

% Compute fractional occupancy
fractional_occupancy = getFractionalOccupancy(vpath,size(vpath,1),options);

% correct Gamma time-course for lags in  model
gamma_emb = padGamma(Gamma_emb,T,options);

% Compute task evoked gamma time-course
gamma_evoked = reshape(gamma_emb,T(1),[],options.K);


%-------------------------
% Compute spectral metrics

% Compute multi-taper spectral estimate
options.win = 128;
options.verbose = 0;

options_mt = struct('Fs',dsample_rate); % Sampling rate - for the 25subj it is 300
options_mt.fpass = [1 40];  % band of frequency you're interested in
options_mt.tapers = [4 7]; % taper specification - leave it with default values
options_mt.p = 0; %0.01; % interval of confidence
options_mt.win = round(2 * dsample_rate); % multitaper window
options_mt.to_do = [1 0]; % turn off pdc
options_mt.order = 0;
options_mt.embeddedlags = options.embeddedlags;

fit_emb = hmmspectramt(data,T,gamma_emb,options_mt);

% Compute spectrum directly from lags in observation model
f = dsample_rate*(0:(length(options.embeddedlags)/2))/length(options.embeddedlags);
freq = (0:length(options.embeddedlags)-1).*(dsample_rate/length(options.embeddedlags));

states = 1:6;
obs_spec = zeros( length(options.embeddedlags),K);
for ii = 1:K
    obs_spec(:,ii) = abs(fft(hmm.state(ii).W.S_W((length(options.embeddedlags)+1)/2,:)));
end
task_obs_spec = obs_spec(1:options.embeddedlags(end)+1,states)*squeeze(mean(gamma_evoked(:,:,states),2))';


%-------------------------
% Compute event metrics


% Compute continuous beta power from the wavelet transform (this is used to
% compute the beta amplitude for each individual visit to a state)
thresh = 2/3;
G = gamma_evoked(emb+1:end-emb-1,:,:); % remove padding
d = reshape(data,T(1),size(T,1));
d = d(emb+1:end-emb-1,:); % remove padding
bp = wt(f_inds_cwt,emb+1:end-emb-1,:);


% Compute the lifetime and amplitude of each event and the
mask = G > thresh; % binary mask of state visits

% Some preallocation
lifetimes_vect = zeros(size(mask))*nan;
lt_vect = zeros(size(mask))*nan;
numevents_vect = zeros(size(mask))*nan;
amps = zeros(size(mask))*nan;

% Define connectivity for finding contiguous state visits in masked gamma
conn = [ 0 1 0; 0 1 0; 0 1 0];

% Main loop
im_pix = numel(mask(:,:,1));
for k = 1:K
    % Use bwconncomp to find contiguous state visits.
    CC = bwconncomp(mask(:,:,k),conn);
    for ii = 1:CC.NumObjects
        % Replace each state visit with its length in samples
        lifetimes_vect(CC.PixelIdxList{ii}+im_pix*(k-1)) = numel(CC.PixelIdxList{ii});
        % Replace first sample of each visit with a one
        numevents_vect(CC.PixelIdxList{ii}(1)+im_pix*(k-1)) = numel(CC.PixelIdxList{ii});
        % Replace each visit with the mean beta power of that time-window
        amps(CC.PixelIdxList{ii}+im_pix*(k-1)) = mean(mean(bp(:,CC.PixelIdxList{ii})));
    end
end
%
lifetimes_vect_ms = 1000*(lifetimes_vect/dsample_rate);

numevents = squeeze(nansum(numevents_vect>1,1));
ltevents = squeeze(nanmean(numevents_vect,1));


%-------------------------
% Colours

% Some colours from Set1 in colorbrewer.org
cols = {[228,26,28]./255,...
    [55,126,184]./255,...
    [77,175,74]./255,...
    [152,78,163]./255,...
    [255,127,1]./255,...
    [166,86,40]./255,...
    [247,129,191]./255};


%% FIGURE 4

figure('Position',[100 100 1512 768]) ;

% Hardcoded trials to plot
trl_inds = {378*334:379*334, 1337*334:1338*334};
trls = [378,1337];

% Axes widths and heights for state visits.
bw = [.05 .3 .05 .3 .05 .3];
bh = [.708 .708 .37 .37 .032 .032];

% Main loop to plot state observation models
for k = 1:options.K

    % Add group of three axes
    [ax1,ax2,ax3] = hmm_util_add_state_axes([bw(k) bh(k)],1.4);

    % Plot observation model power spectrum

    ylim(ax1,[0 8.5e-5])
    area(ax1,f,obs_spec(1:options.embeddedlags(end)+1,k),...
        'FaceColor',cols{k},'EdgeColor',cols{k},'FaceAlpha',.25);
    plot(ax1,f,obs_spec(1:options.embeddedlags(end)+1,k),'Color',cols{k},'linewidth',3)
    xlabel(ax1,'Frequency (Hz)');ylabel(ax1,'Power')
    xticks(ax1,[0,10,20,30,40]);

    % Plot line from middle of autocovariance
    xlim(ax2,[-emb,emb]);xticks(ax2,[-11,-5,0,5,11])
    plot(ax2,options.embeddedlags,hmm.state(k).W.S_W((length(options.embeddedlags)+1)/2,:),'color',cols{k},'linewidth',2)
    ylabel(ax2,'Covariance');

    % Plot whole autocovariance matrix
    xlim(ax3,[-emb,emb]);ylim(ax3,[-emb emb])
    imagesc(ax3,options.embeddedlags,options.embeddedlags,hmm.state(k).W.S_W)
    xticks(ax3,[-11,-5,0,5,11]);yticks([-11,-5,0,5,11])
    xlabel(ax3,'Lag');ylabel(ax3,'Lag')
    axis(ax3,'square')
end


xshift = .05;
time_vect = linspace(-2,2,335);

% Loop for plotting example trials
for ii = 2:-1:1

    % Plot MEG signal time-course
    inds = trl_inds{ii};
    ax1 = subplot(3,4,ii+2);hold on;grid on
    ax1.Position(1) = ax1.Position(1)+xshift;
    ax1.Position(2) = ax1.Position(2)-.08;
    plot(time_vect,data(inds),'k')
    ylim([-6 6])
    ylabel('MEG sensor signal (a.u)');

    % Plot movement time box
    pos = ax1.Position;
    pos(2) = pos(2) + .25;
    pos(4) = .04;
    ax1b = axes('Position',pos);
    mov = zeros(size(time_vect));
    mov_inds = (time_vect>-mt(trls(ii))/1000) .* (time_vect<0);
    mov(find(mov_inds)) = 1;
    h = area(ax1b,time_vect,mov);
    h(1).FaceColor = [.7 .7 .7];
    ylim(ax1b,[0,1.25])
    yticks(ax1b,[]);
    set(ax1b, 'visible', 'off')
    ax1b.XAxis.Visible = 'on';
    ylabel(ax1b,{'Movement','time'},'Visible','on')

    % Plot wavelet transform of trial
    ax2 = subplot(3,4,ii+6);
    ax2.Position(2) = ax2.Position(2)-.05;
    ax2.Position(1) = ax2.Position(1)+xshift;
    [wt2,wf2] = cwt(data(inds),dsample_rate);
    contourf(time_vect,wf2,abs(wt2),24,'linestyle','none')
    hold on;grid on;
    ylabel('Frequency (Hz)')

    % Plot state time course
    ax3 = subplot(3,4,ii+10);
    ax3.Position(2) = ax3.Position(2);
    ax3.Position(1) = ax3.Position(1)+xshift;
    gg = gamma_emb(inds,:);
    gg(1:emb,:) = nan;
    gg(end-emb:end,:) = nan;
    h = area(time_vect,gg);
    ylim([0,1])
    for k = 1:options.K
        h(k).FaceColor = cols{k};
    end
    xlabel('Time (seconds)')
    ylabel('Probability')
    grid on;hold(ax3,'on')

    % Add movement time dashed lines to all plots
    plot(ax1, [-mt(trls(ii))/1000,-mt(trls(ii))/1000],[-6,6],':k','linewidth',3)
    plot(ax1, [0,0],[-6,6],':k','linewidth',3)
    plot(ax2, [-mt(trls(ii))/1000,-mt(trls(ii))/1000],[1,35],':k','linewidth',3)
    plot(ax2, [0,0],[1,35],':k','linewidth',3)
    plot(ax3, [-mt(trls(ii))/1000,-mt(trls(ii))/1000],[0,1],':k','linewidth',3)
    plot(ax3, [0,0],[0,1],':k','linewidth',3)

end

% Save
figpath = fullfile(config.figpath,'hmm_fig4_statesummary');
saveas(gcf,figpath,'png')
saveas(gcf,figpath,'tiff')


%% FIGURE 5

% State indices of the alpha/beta and beta states
states_to_plot = [5 3];

% Find indeces for time and frequency windows.
f_inds_cwt = find((freq_vect_cwt > 15) .* (freq_vect_cwt < 30));
t_inds = find( (time_vect(emb+1:end-emb-1)>1.25) .* (time_vect(emb+1:end-emb-1)<2) );
t_inds2 = find( (time_vect(emb+1:end-emb-1)>0) .* (time_vect(emb+1:end-emb-1)<1) );

% Preallocate arrays
subj_sp = zeros(size(wt,1),33);
subj_fo = zeros(K,33);
subj_po = zeros(1,33);
subj_am = zeros(K,33);

% Compute metrics for each subject
for ii = 1:33
    trls = subj == ii;
    subj_sp(:,ii) = squeeze(mean(wt(:,:,trls),[2,3]));
    subj_fo(:,ii) = squeeze(mean(G(:,trls,:),[1,2]));
    subj_po(:,ii) = squeeze(mean(wt(f_inds_cwt,:,trls),[1,2,3]));
    subj_am(:,ii) = squeeze(nanmedian(amps(:,trls,:),[1,2]));
end

%------------
% Main Figure

figure('Position',[100 100 1152 768]);

% State power spectrum
subplot(241);grid on;hold on
for ii = 1:options.K
    plot(f,obs_spec(1:options.embeddedlags(end)+1,ii),'linewidth',2,'color',cols{ii} )
end
xlabel('Frequency (Hz)')
ylabel('Power (a.u.)')
title({'Observed state spectra',''})
yticklabels(0:9)

% State fractional occupancy
subplot(242);grid on;hold on
for ii = 1:options.K
    bar(ii,fractional_occupancy(ii),'FaceColor',cols{ii});
end
xticks(1:options.K);
grid on
xlabel('State')
ylabel('Fractional Occupancy')
title({'State Fractional Occupancy',''})

% State lifetimes
subplot(243)
distributionPlot(lifetimes_ms,'color',cols)
grid on
ylim([0,1])
xlabel('State')
ylabel('Lifetime (ms)')
title({'State Lifetime',''})

% State interval times
subplot(244)
distributionPlot(intervaltimes_ms,'color',cols)
ylim([0,3.5])
xlabel('State')
ylabel('Interval Time (s)')
title({'State Intervaltime',''})

% State wavelet power spectrum
subplot(234);grid on;hold on
p = patch([15 15 30 30],[0 2.2e-3 2.2e-3 0],[ .8 .8 .8],'FaceAlpha',.5,'LineStyle','none');
plot(freq_vect_cwt,subj_sp,'linewidth',1.5)
xlim([2,36]);grid on;
ylim([0,2.2e-3])
xlabel('Frequency (Hz)')
ylabel('Power (a.u.)')
title({'Per-Subject Spectrum',''})

% State to plot - might be different for your inference.
k = states_to_plot(1);

% Compute linear fit and correlation between state FO and wavelet power
md = fitlm(subj_fo(k,:)',subj_po','linear','VarNames',{'FO','Power'});
betas = table2array(md.Coefficients(:,1));
[r,p] = corrcoef( subj_fo(k,:)',subj_po' );

% Plot relationship between state FO and wavelet power.
subplot(2,3,5);grid on;hold on
plot(subj_fo(k,:),subj_po,'o','color',cols{k},'MarkerFaceColor',cols{k})
x = linspace(min(subj_fo(k,:)),max(subj_fo(k,:)));
plot( x, betas(1) + betas(2).*x,'k')
xlabel('Fractional Occupancy')
ylabel('Power (a.u.)')
xlim([0,.4]);ylim([0,2e-3])
text(.15,.5e-3,['r(31)=' num2str(round(r(1,2),3)) ' p<.001'],'FontSize',12)
title({'State 4 FO and Power',''})

% State to plot - might be different for your inference.
k = states_to_plot(2);

% Compute linear fit and correlation between state FO and wavelet power
md = fitlm(subj_fo(k,:)',subj_po','linear','VarNames',{'FO','Power'});
betas = table2array(md.Coefficients(:,1));
[r,p] = corrcoef( subj_fo(k,:)',subj_po' );

% Plot relationship between state FO and wavelet power.
subplot(2,3,6);grid on;hold on
plot(subj_fo(k,:),subj_po,'o','color',cols{k},'MarkerFaceColor',cols{k})
x = linspace(min(subj_fo(k,:)),max(subj_fo(k,:)));
plot( x, betas(1) + betas(2).*x,'k')
xlabel('Fractional Occupancy')
ylabel('Power (a.u.)')
xlim([0,.4]);ylim([0,2e-3])
text(.15,.5e-3,['r(31)=' num2str(round(r(1,2),3)) ' p=' num2str(round(p(1,2),4))],'FontSize',12)
title({'State 5 FO and Power',''})

% Adjust fontsize
set(findall(gcf,'type','text'),'FontSize',15)

% Save figures
figpath = fullfile(config.figpath,'hmm_fig5_temporalstats');
saveas(gcf,figpath,'png')
saveas(gcf,figpath,'tiff')


%% FIGURE 6

figure('Position',[100 100 1512 768])
time_vect = linspace(-2,2,334);

% Some factors
step = .165; % Space between left-hand-side of columns
scale = [.1 .275]; % Size of HMM tf plot at top of column
for k = 1:6

    % Add axes for HMM TF plot
    [ax31,ax32,ax33] = hmm_util_add_tf_axes( [.025+(k-1)*step,.65],scale );

    % Plot HMM regularised TF
    xx = squeeze(mean(gamma_evoked(:,:,k),2)) * obs_spec(1:options.embeddedlags(end)+1,k)';
    xx(1:emb,:) = nan;
    xx(end-emb:end,:) = nan;
    contourf(ax31,time_vect,f,xx', 48, 'linestyle','none' )
    title(ax31,{['State ' num2str(k)],' TF Transform'});
    ylabel(ax31,'Frequency (Hz)')
    xlabel(ax31,'Time (secs)')
    xlim(ax31,[-2 2])
    xticks(ax31,-2:.5:2);
    cb = colorbar(ax31);
    cb.Position(1) = cb.Position(1) + .025;
    cb.Position(2) = cb.Position(2) + .025;
    cb.Position(4) = cb.Position(4) - .05;

    % Plot evoked state time-course
    plot(ax33,time_vect,squeeze(mean(gamma_evoked(:,:,k),2)),'color',cols{k},'linewidth',2);
    % Plot state spectrum
    plot(ax32,obs_spec(1:options.embeddedlags(end)+1,k),f,'color',cols{k},'linewidth',2);

    ylabel(ax32,'Frequency (Hz)')
    xlabel(ax32,{'Power',''})
    ylabel(ax33,'FO')
    xlabel(ax33,'Time (secs)')
    xticks(ax33,-2:1:2)

    % Plot state raster
    ax = axes( 'Position', [.025+(k-1)*step + scale(1)*.2, .075, .8*scale(1), .5]);
    x = squeeze(gamma_evoked(:,1:4:length(T),k) < .5);
    imagesc(ax,time_vect,1:size(gamma_evoked,2),x');
    colormap(ax,'bone');
    xlabel('Time (secs)')
    ylabel('Trial')
    title('State Visits')
    xticks(ax,-2:1:2)

end
set(findall(gcf,'type','text'),'FontSize',14)

% Save
figpath = fullfile(config.figpath,'hmm_fig6_evokedstates');
saveas(gcf,figpath,'png')
saveas(gcf,figpath,'tiff')


%% Figure 7
states_to_plot = [5 3];

time_vect = linspace(-2,2,T(1));

lt = lifetimes_vect_ms;
time_vect_emb = time_vect(emb+1:end-emb-1);
lt(time_vect_emb<-1.5,:,:) = nan;
lt(time_vect_emb>1.5,:,:) = nan;

figure('Position',[100 100 1512 768])

% Plot wavelet transform
ax = axes('Position',[0.4300 0.3200 0.1800 0.4800]);
contourf( time_vect,freq_vect_cwt,nanmean(wt,3)-bl,48,'linestyle','none')
colorbar;
hmm_util_redblue_colourmap(ax,min(min(nanmean(wt,3)-bl)),max(max(nanmean(wt,3)-bl)));
grid on;ylim([0 35])
title('Wavelet Transform','FontSize',24)
ylabel('Frequency (Hz)')
xlabel('Time (secs)')
xticks(-2:2)
grid on;hold on
plot([0,0],[0,35],'k:','linewidth',2)
xticks(-2:.5:2);


% Plot HMM regularised TF for all states.
[ax31,ax32,ax33] = hmm_util_add_tf_axes( [.085,.2],[.225 .6] );

% HMM TF
states = 1:K;
task_bl = mean(task_obs_spec(:,bl_inds,:),[2,3]);
xx = task_obs_spec-task_bl;
xx(:,1:emb) = nan;
xx(:,end-emb:end) = nan;
contourf(ax31,time_vect,f,xx,48,'linestyle','none')
cb = colorbar(ax31);
cb.Position(1) = .32;
cb.Position(2) = cb.Position(2) + .025;
cb.Position(4) = cb.Position(4) - .05;
hmm_util_redblue_colourmap(ax31,-1.7e-5,1.e-5);grid on;ylim([0 35])
ylim(ax31,[0 35])
ylabel(ax31,'Frequency (Hz)')
xlabel(ax31,'Time (secs)')
xticks(ax31,-2:2)
caxis(ax31,[-1.7e-5 1e-5])
title(ax31,'HMM Transform','FontSize',24)
plot(ax31,[0,0],[0,35],'k:','linewidth',2)
xticks(ax31,-2:.5:2);
grid(ax31,'on')

% HMM evoked time-courses and state spectra
for ii = 1:6
    plot(ax33,time_vect,squeeze(mean(gamma_evoked(:,:,ii),2)),'color',cols{ii},'linewidth',2);
    plot(ax32,obs_spec(1:options.embeddedlags(end)+1,ii),f,'color',cols{ii},'linewidth',2);
end
ylim(ax33,[0,.45])
ylim(ax32,[0 35])
ylabel(ax32,'Frequency (Hz)')
xlabel(ax32,{'Power',''})

ylabel(ax33,'FO')
xlabel(ax33,'Time (secs)')

% Compute time-window contrasts
bl_inds = find((time_vect > -1.5) .* (time_vect < -.75));
rb_inds = find((time_vect > .5) .* (time_vect < 1.5));

subj_gamma = zeros(334,6,33);
subj_gamma_effect = zeros(6,33);
subj_lt = zeros(311,6,33);
subj_am = zeros(311,6,33);

for ss = 1:33
    subj_gamma(:,:,ss) = squeeze(mean(gamma_evoked(:,subj==ss,:),2));
    subj_gamma_effect(:,ss) = mean(subj_gamma(rb_inds,:,ss)) - mean(subj_gamma(bl_inds,:,ss));
    subj_lt(:,:,ss) = squeeze(nanmean(lt(:,subj==ss,:),2));
    subj_am(:,:,ss) = squeeze(nanmean(amps(:,subj==ss,:),2));
end

[g_h,g_p,g_ci,g_stats] = ttest2( squeeze(mean(subj_gamma(rb_inds,:,:)))',squeeze(mean(subj_gamma(bl_inds,:,:)))');
edge = 25;

bl_inds = find((time_vect_emb(edge:end-edge) > -1.5) .* (time_vect_emb(edge:end-edge) < -.75));
rb_inds = find((time_vect_emb(edge:end-edge) > .5) .* (time_vect_emb(edge:end-edge) < 1.5));
[l_h,l_p,l_ci,l_stats] = ttest2( squeeze(nanmean(subj_lt(rb_inds,:,:)))',squeeze(nanmean(subj_lt(bl_inds,:,:)))');
[a_h,a_p,a_ci,a_stats] = ttest2( squeeze(nanmean(subj_am(rb_inds,:,:)))',squeeze(nanmean(subj_am(bl_inds,:,:)))');


% Plot evoked FO time-series
ax = subplot(333);grid on; hold on
p = patch([-1.5 -1.5 -.75 -.75],[0 .3 .3 0],[ .8 .8 .8],'FaceAlpha',.5,'LineStyle','none');
p = patch([.5 .5 1.5 1.5],[0 .3 .3 0],[ .8 .8 .8],'FaceAlpha',.5,'LineStyle','none');

ax.Position(3) = ax.Position(3) - .05;
for k = 1:2
    l = mean(gamma_evoked(emb+1:end-emb-1,:,states_to_plot(k)),2);
    e = std(gamma_evoked(emb+1:end-emb-1,:,states_to_plot(k)),[],2) ./ sqrt(size(l,1));
    fill( [time_vect_emb fliplr(time_vect_emb)], [l+e;flipud(l-e)]',cols{states_to_plot(k)},'FaceAlpha',.5,'LineStyle','none')
    plot(ax,time_vect_emb,l,'linewidth',2,'color',cols{states_to_plot(k)});
    xlabel('Time (secs)')
    ylabel('FO')
    title('Fractional Occupancy')

end

% Plot evoked FO contrasts
pos = ax.Position;
pos(1) = .885;
pos(3) = .075;
axes('Position',pos);grid on;hold on;ylim([-3.5 3.5]);
for k = 1:2
    bar(k,g_stats.tstat(states_to_plot(k)),'FaceColor',cols{states_to_plot(k)});
    if g_p(states_to_plot(k)) < .01
        text(k,g_stats.tstat(states_to_plot(k))+.2,'**','horizontalalign','center');
    elseif g_p(states_to_plot(k)) < .05
        text(k,g_stats.tstat(states_to_plot(k))+.2,'*','horizontalalign','center');
    end
end
title('Rebound-Baseline')
ylabel('t-value')
xticks(1:2);
xticklabels(4:5)
xlabel('State')

% Plot evoked lifetime time-series
ax = subplot(336);grid on;hold on
ax.Position(3) = ax.Position(3) - .05;
p = patch([-1.5 -1.5 -.75 -.75],[200 600 600 200],[ .8 .8 .8],'FaceAlpha',.5,'LineStyle','none');
p = patch([.5 .5 1.5 1.5],[200 600 600 200],[ .8 .8 .8],'FaceAlpha',.5,'LineStyle','none');

for k = 1:2
    l = squeeze(nanmean(lt(edge:end-edge,:,states_to_plot(k)),2));
    e = squeeze(nanstd(lt(edge:end-edge,:,states_to_plot(k)),[],2)) ./ sqrt(size(l,1));
    fill(ax, [time_vect_emb(edge:end-edge) fliplr(time_vect_emb(edge:end-edge))], [l+e;flipud(l-e)]',cols{states_to_plot(k)},'FaceAlpha',.5,'LineStyle','none')
    plot(ax,time_vect_emb(edge:end-edge),l,'color',cols{states_to_plot(k)},'linewidth',2)
    grid on
    xlim([-2 2])
    xlabel('Time (secs)')
    ylabel('Lifetime (ms)')
    title('Lifetimes')
end
ylim([200 550])

% plot evoked lifetime contrasts
pos = ax.Position;
pos(1) = .885;
pos(3) = .075;
axes('Position',pos);grid on;hold on;ylim([-3.5 3.5]);
for k = 1:2
    bar(k,l_stats.tstat(states_to_plot(k)),'FaceColor',cols{states_to_plot(k)});
    if l_p(states_to_plot(k)) < .01
        text(k,l_stats.tstat(states_to_plot(k))+.2,'**','horizontalalign','center');
    elseif l_p(states_to_plot(k)) < .05
        text(k,l_stats.tstat(states_to_plot(k))+.2,'*','horizontalalign','center');
    end
end
title('Rebound-Baseline')
ylabel('t-value')
xticks(1:2);
xticklabels(4:5)
xlabel('State')

% Plot evoked amplitude time-series
ax = subplot(339);grid on; hold on
ax.Position(3) = ax.Position(3) - .05;
p = patch([-1.5 -1.5 -.75 -.75],[1e-3 2e-3 2e-3 1e-3],[ .8 .8 .8],'FaceAlpha',.5,'LineStyle','none');
p = patch([.5 .5 1.5 1.5],[1e-3 2e-3 2e-3 1e-3],[ .8 .8 .8],'FaceAlpha',.5,'LineStyle','none');

for k = 1:2
    l = squeeze(nanmean(amps(edge:end-edge,:,states_to_plot(k)),2));
    e = squeeze(nanstd(amps(edge:end-edge,:,states_to_plot(k)),[],2)) ./ sqrt(size(l,1));
    fill(ax, [time_vect_emb(edge:end-edge) fliplr(time_vect_emb(edge:end-edge))], [l+e;flipud(l-e)]',cols{states_to_plot(k)},'FaceAlpha',.5,'LineStyle','none');
    plot(ax,time_vect_emb(edge:end-edge),l,'color',cols{states_to_plot(k)},'linewidth',2)
    grid on
    xlim([-2 2])
    xlabel('Time (secs)')
    ylabel('Amplitude')
    title('Amplitude')
end

% Plot evoked amplitude contrasts.
pos = ax.Position;
pos(1) = .885;
pos(3) = .075;
axes('Position',pos);grid on;hold on;ylim([-3.5 3.5]);
for k = 1:2
    bar(k,a_stats.tstat(states_to_plot(k)),'FaceColor',cols{states_to_plot(k)});
    if a_p(states_to_plot(k)) < .01
        text(k,a_stats.tstat(states_to_plot(k))+.2,'**','horizontalalign','center');
    elseif a_p(states_to_plot(k)) < .05
        text(k,a_stats.tstat(states_to_plot(k))+.2,'*','horizontalalign','center');
    end
end
title('Rebound-Baseline')
ylabel('t-value')
xticks(1:2);
xticklabels(4:5)
xlabel('State')

% Adjust font size
set(findall(gcf,'type','text'),'FontSize',14)

% Save figures
figpath = fullfile(config.figpath,'hmm_fig7_taskstats');
saveas(gcf,figpath,'png')
saveas(gcf,figpath,'tiff')
