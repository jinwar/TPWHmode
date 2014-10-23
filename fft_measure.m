function [amp phase] = fft_measure(data,delta,periods)

fftdata =fft(data);
N = length(data);
T = N*delta;

if mod(N,2)
     faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
else
     faxis = [0:N/2,-N/2+1:-1]*(1/T);
end

freqs = 1./periods;

amp = interp1(faxis,abs(fftdata),freqs,'nearest');
phase = interp1(faxis,angle(fftdata),freqs,'nearest');
