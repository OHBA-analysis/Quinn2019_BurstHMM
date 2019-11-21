function [ax1,ax2,ax3] = add_hmm_tf( bl, scale )


pos1 = [bl(1)+.2*scale(1) bl(2)+.2*scale(2) .8*scale(1) .8*scale(2)];
pos2 = [bl(1)             bl(2)+.2*scale(2) .2*scale(1) .8*scale(2)];
pos3 = [bl(1)+.2*scale(1) bl(2)             .8*scale(1) .2*scale(2)];

ax1 = axes('Position',pos1); % TF
grid on;hold on
set(ax1,'XTickLabels',[],'YTickLabels',[]);

ax2 = axes('Position',pos2); % F
grid on;hold on;
set(ax2, 'XDir', 'reverse','xaxisLocation','top');

ax3 = axes('Position',pos3); % T
grid on;hold on;
set(ax3,'yaxisLocation','right')



end