%% Overview
%
% This script will create a set of simple oscillatory time-series with
% different amplitudes and temporal dynamics. The power spectrum of each
% example is computed and a summary figure is created and saved out into the config.figpath directory

%% Set up and check paths.

config = hmm_0_initialise;

%% Create simulated signals & estimate Power spectra

seconds = 4;
sample_rate = 500;
time_vect = linspace(0,seconds,seconds*sample_rate)';
freq = 4;

base = sin(2*pi*freq*time_vect);
x = zeros(length(base),3,4);

% Static Amplitude Increase
amplitudes = [.5 1 2];
for ii = 1:3
    x(:,ii,1) = base.*amplitudes(ii);
end
[freq_vect,PSD,ASD] = hmm_util_get_spectra(x(:,:,1),sample_rate);

% Event Amplitude Increase
durations = [500];
amplitudes = [.5 1 2];
for ii = 1:3
    inds = (1:durations(1)) + 1000-(durations(1)/2);
    x(inds,ii,2) = base(1:durations(1)).*amplitudes(ii);
end
[freq_vect,PSD2,ASD2] = hmm_util_get_spectra(x(:,:,2),sample_rate);

% Event Duration Increase
durations = [500 1000 1500];
for ii = 1:3
    inds = (1:durations(ii)) + 1000-(durations(ii)/2);
    x(inds,ii,3) = base(1:durations(ii));
end
[freq_vect,PSD3,ASD3] = hmm_util_get_spectra(x(:,:,3),sample_rate);

% Number of Events Increase
starts = [1 500 1000];
for ii = 1:3
    if ii > 0
        x(750:750+499,ii,4) = base(1:500);
    end
    if ii > 1
        x(1:1+499,ii,4) = base(1:500);
    end
    if ii > 2
        x(1500:1500+499,ii,4) = base(1:500);
    end        
end
[freq_vect,PSD4,ASD4] = hmm_util_get_spectra(x(:,:,4),sample_rate);

PSD = cat(3,PSD,PSD2,PSD3,PSD4);
ASD = cat(3,ASD,ASD2,ASD3,ASD4);

%% Make summary figure

freq_lims = [0,2250;0,325;0,325;0,325];
ts_lims = [2.1,2.1,2.1,2.1];
X = [freq_vect,fliplr(freq_vect)];

figure('Position',[100 100 768 768])
set(gcf, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 16)
for ii = 1:3
    for jj = 1:4
        ax = subplot(4,3,ii+((jj-1)*3));
        ax.Position(4) = .11;
        Y = [PSD(:,ii,jj);zeros(2000,1)];
        patch(X,Y,'blue',...
        'EdgeColor','blue','FaceAlpha',.25,'LineWidth',2)
        xlim([1 7]);grid on;ylim(freq_lims(jj,:))
        if ii == 1
            ylabel('Power fft(x)^2/N')
        end
        
        if jj == 1
            ax.Position(2) = ax.Position(2)+.015;
        elseif jj == 4
            ax.Position(2) = ax.Position(2)-.015;
            xlabel('Frequency (Hz)')
        end
        
        pos = ax.Position;
        ax2 = axes('Position',[pos(1) pos(2)+.1 pos(3) .075]);hold on
        plot( ones(2000,1),'Color',[.7 .7 .7])
        plot( -ones(2000,1),'Color',[.7 .7 .7])
        plot( x(:,ii,jj),'k','linewidth',2)
        ylim([-ts_lims(jj) ts_lims(jj)])
        set(gca,'xaxisLocation','top')

        set(ax2,'Visible',false)
        ax2.YAxis.Visible='on';
        ax2.YAxis.TickValues = [-1 0 1];
        ax2.YAxis.FontSize=16;
        
        if jj == 1
            plot( 251:750,ones(500,1)*2,'k','linewidth',3)
            text( 500,2.5,'1 Second','HorizontalAlignment','center')

        end
    end
end

% save figure
figpath = fullfile(config.figpath,'hmm_1_dynamics_figure');
saveas(gcf,figpath,'png');
saveas(gcf,figpath,'tiff');
