function plot_spectrum(data,delta,frange,fignum)

N = length(data);
T = N*delta;
if ~exist('frange','var')
	frange = [0 1/2/delta];
end
if ~exist('fignum','var')
	fignum = 38;
end


if mod(N,2)
     faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
else
     faxis = [0:N/2,-N/2+1:-1]*(1/T);
end

fftdata = fft(data);

figure(fignum)
clf
subplot(1,2,1)
plot(faxis,abs(fftdata),'o');
xlim(frange);
subplot(1,2,2)
plot(faxis,angle(fftdata),'o');
xlim(frange);
