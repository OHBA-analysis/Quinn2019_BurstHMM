function [freq_vect,PSD,ASD] = hmm_util_get_spectra(x,sample_rate,nfft)

if nargin < 3 || isempty(nfft)
    nfft = [];
end
if nargin < 2 || isempty(sample_rate)
    sample_rate = 1;
end
fx = fftshift(fft(x,nfft,1),1);
freq_vect = sample_rate*linspace(-size(x,1)/2,size(x,1)/2,size(fx,1))./size(x,1);

PSD = (abs(fx).^2) ./ size(x,1);
ASD = abs(fx) ./ size(x,1);


end

