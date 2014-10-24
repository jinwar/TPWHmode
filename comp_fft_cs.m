clear

figure(48)
clf
hold on

load pa5mod

is_cs_dtp = 1;
Wph = 30;
iprange = [2:4];

for ip=iprange

setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
sample_range = parameters.maxstadist;
center_stla = 9.3887;
center_stlo = -144.88;
periods = parameters.periods;
min_cs_num = 5;

v1_0 = interp1(pa5mod.T_first,pa5mod.phv_first,periods(ip));
v2_0 = interp1(pa5mod.T_fund,pa5mod.phv_fund,periods(ip));


csmatfiles = dir(fullfile('CSmeasure',['*_cs_',parameters.component,'.mat']));

% Gather event information 
clear event_data event_parastr
event_num = 0;
for ie = 1:length(csmatfiles)
%for ie = 1:5
	filename = fullfile('CSmeasure',csmatfiles(ie).name)
	load(filename);
	evla = eventcs.evla;
	evlo = eventcs.evlo;
	stlas = eventcs.stlas;
	stlos = eventcs.stlos;

	stadists = distance(stlas,stlos,center_stla,center_stlo);
	[temp center_sta] = min(stadists);

	stadists = distance(stlas,stlos,stlas(center_sta),stlos(center_sta));
	stadists = deg2km(stadists);

	nbstaids = find(stadists < sample_range);

	% reconstruct the data format to local array index
	CSnum = 0;
	localcs = [];
	localcs.stlas = stlas(nbstaids);
	localcs.stlos = stlos(nbstaids);
	for ista = 1:length(localcs.stlas)
		localcs.dists(ista) = vdist(localcs.stlas(ista),localcs.stlos(ista),evla,evlo)/1e3;
	end
	localcs.evla = evla;
	localcs.evlo = evlo;
	localcs.center_sta = find(nbstaids == center_sta);
	localcs.period = parameters.periods(ip);
	localcs.id = eventcs.eventmatfile;
	for ista = 1:length(nbstaids)
		localcs.amps(ista) = eventcs.autocor(nbstaids(ista)).fft_amp(ip);
		localcs.phases(ista) = eventcs.autocor(nbstaids(ista)).fft_phase(ip);
		localcs.cs_amps(ista) = eventcs.autocor(nbstaids(ista)).amp(ip).^0.5;
		localcs.cs_dtps(ista) = 0;
		localcs.isgood(ista) = -1;
	end
	localcs.isgood(localcs.center_sta) = 1;
	% convert phase to dtp
	dphases = localcs.phases - localcs.phases(localcs.center_sta);
	localcs.dtps = -dphases/2/pi*periods(ip);
	ddists = localcs.dists - localcs.dists(center_sta);
	localcs.dtps = corr_cycle_skip(ddists,localcs.dtps,periods(ip),parameters.refphv(ip));
	localcs.fft_dtps = localcs.dtps;
	% normalize the amplitude
%	normamp = localcs.amps(localcs.center_sta);
	normamp = median(localcs.amps);
	localcs.true_amps = localcs.amps;
	for ista = 1:length(nbstaids)
		localcs.amps(ista) = localcs.amps(ista)./normamp;
	end

	% find the CS measurements that contains the center station
	for ics = 1:length(eventcs.CS)
		if eventcs.CS(ics).sta1==center_sta 
			staid = find(nbstaids==eventcs.CS(ics).sta2);
			localcs.cs_dtps(staid) = -eventcs.CS(ics).dtp(ip);
			localcs.isgood(staid) = eventcs.CS(ics).isgood(ip);
			localcs.cohere(staid) = eventcs.CS(ics).cohere(ip);
		elseif eventcs.CS(ics).sta2==center_sta 
			staid = find(nbstaids==eventcs.CS(ics).sta1);
			localcs.cs_dtps(staid) = eventcs.CS(ics).dtp(ip);
			localcs.isgood(staid) = eventcs.CS(ics).isgood(ip);
			localcs.cohere(staid) = eventcs.CS(ics).cohere(ip);
		end
	end % ics
	if is_cs_dtp
		localcs.dtps = localcs.cs_dtps;
	end
	ind = find(localcs.isgood<0);
	localcs.stlas(ind) = [];
	localcs.stlos(ind) = [];
	localcs.dists(ind) = [];
	localcs.amps(ind) = [];
	localcs.dtps(ind) = [];
	localcs.cs_dtps(ind) = [];
	localcs.cs_amps(ind) = [];
	localcs.true_amps(ind) = [];
	localcs.fft_dtps(ind) = [];
	localcs.isgood(ind) = [];
	localcs.cohere(ind) = [];
	localcs.id = eventcs.id;
	stadists = distance(localcs.stlas,localcs.stlos,center_stla,center_stlo);
	[temp center_sta] = min(stadists);
	localcs.center_sta = center_sta;
    localcs.azi = azimuth(center_stla,center_stlo,evla,evlo)+180;
	para = polyfit(localcs.dtps,localcs.dists,1);
	localcs.OPW_v = para(1);
	localcs.amps = localcs.amps(:);
	localcs.dtps = localcs.dtps(:);
	localcs.dists = localcs.dists(:);
	localcs.dtps = corr_cycle_skip(ddists,localcs.dtps,periods(ip),localcs.OPW_v);
	localcs.amp_diff = max(localcs.amps)./min(localcs.amps);
	localcs.Wph = Wph;
	event_data(ie) = localcs;
end

cmap = colormap('lines');
cmapx = linspace(1,length(periods),size(cmap,1));
color = interp1(cmapx,cmap,ip);
for ie=1:length(event_data)
	ddists = event_data(ie).dists-event_data(ie).dists(event_data(ie).center_sta);
	amps = event_data(ie).amps;
	dtp_diff = event_data(ie).cs_dtps-event_data(ie).fft_dtps;
	subplot(1,2,1)
	hold on
	ax(ip) = plot(ddists,dtp_diff,'x','color',color);
	subplot(1,2,2)
	hold on
	ax(ip) = plot(amps,dtp_diff,'x','color',color);
	label(ip) = {num2str(periods(ip))};
end

end
legend(ax(iprange),label(iprange));
