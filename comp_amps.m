clear;

fund_path = 'pa5_0_CSmeasure';
first_path = 'pa5_1_CSmeasure';
setup_parameters;
periods = parameters.periods;

CSfiles = dir(fullfile(fund_path,'*.mat'));

for ie=1:length(CSfiles)
	fund(ie) = load(fullfile(fund_path,CSfiles(ie).name));
	first(ie) = load(fullfile(first_path,CSfiles(ie).name));
end

for ip=1:length(fund(ie).eventcs.autocor(1).amp)
	for ie=1:length(fund)
		clear amps
		for ista = 1:length(fund(ie).eventcs.autocor)
			amps(ista) = fund(ie).eventcs.autocor(ista).amp(ip); 
		end
		fund_mean(ip,ie) = mean(amps);
		clear amps
		for ista = 1:length(fund(ie).eventcs.autocor)
			amps(ista) = first(ie).eventcs.autocor(ista).amp(ip); 
		end
		first_mean(ip,ie) = mean(amps);
	end
end

figure(58)
clf
for ie = 1:length(fund)
	for ip=1:size(first_mean,1)
		subplot(3,3,ip)
		hold on
		semilogy(ie,first_mean(ip,ie)./fund_mean(ip,ie),'x')
		title(['First/fund, ',num2str(periods(ip)),' s'],'fontsize',18);
	end
end
