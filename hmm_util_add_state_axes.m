function [ax1,ax2,ax3] = hmm_util_add_state_axes( bl, scale )

pos1 = [ bl(1)+0.0794*scale  bl(2)+.01              0.0603*scale    0.1577*scale];
pos2 = [ bl(1)               bl(2)+.05+0.075*scale  0.0603*scale    0.0760*scale];
pos3 = [ bl(1)               bl(2)+.01              0.0603*scale    0.0760*scale];

%figure

ax1 = axes('Position',pos1); % Spectrum
grid on;hold on
ax2 = axes('Position',pos2); % Chirp
grid on;hold on
ax3 = axes('Position',pos3); % autocov
grid on;hold on

end