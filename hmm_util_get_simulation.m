function [data,x,time_vect,sample_rate] = hmm_util_get_simulation
% hmm_util_get_simulation - Generate simulated time-series

rng(42)

% A noise time-course is created by direct pole-placement. This creates a filte
% with an approximately 1/f frequency profile.
seconds = 20;
sample_rate = 256;
time_vect = linspace(0,seconds,seconds*sample_rate)';

% the value in poly is the filter root. it can vary between 0 < x < 1 where 0
% is white noise and 1 is an extremely sloped spectrum
% a is then the denominator polynomial of a digital filter which is used to
% create the signal
a = -poly(.92);

noise = .2*randn(size(time_vect));
x = filter(1,a,noise);


% Define burst occurances and durations
starts = [100 500 900 1500 1800 2400 3000 3500 4000 4300 4800] / 2;
starts2 = starts + 2560;
duration = [150 100 50 100 50 200 50 100 100 200 50] + 1;

% Create burst time-courses
x2 = zeros(size(x));
x3 = zeros(size(x));    

for ii = 1:length(starts)
    % Add slow burst
    tmp = sin( 2*pi*25*time_vect(starts(ii):starts(ii)+duration(ii)));
    tmp = tmp.*tukeywin(length(tmp),.25);
    x2(starts(ii):starts(ii)+duration(ii)) = x2(starts(ii):starts(ii)+duration(ii)) + tmp;
    
    % Add fast burst
    tmp = sin( 2*pi*40*time_vect(starts2(ii):starts2(ii)+duration(ii)));
    tmp = tmp.*tukeywin(length(tmp),.25);
    x3(starts2(ii):starts2(ii)+duration(ii)) = x3(starts2(ii):starts2(ii)+duration(ii)) + tmp;
end

% sum the noise and burst time-courses to create the final simulation
amp2 = .63;
amp3 = .551;
data = 1*x' + amp2.*x2' + amp3.*x3';

% Concatenate different parts of simulation
x = cat(2,x,x2,x3);

end
