
% this program try to follow the two plane wave method, except for 
% allowing the phase velocity of the second wave different from the
% first one

clear;

load pa5mod

isfakedata = 0;
test_r = 0.03;
test_N = 5;
Wph = 30;

for ip=[1:8]

setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
sample_range = parameters.maxstadist;
center_stla = 9.3887;
center_stlo = -144.88;
periods = parameters.periods;
min_cs_num = 7;

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
		localcs.amps(ista) = eventcs.autocor(nbstaids(ista)).amp(ip)^0.5;
		localcs.dtps(ista) = 0;
		localcs.isgood(ista) = 0;
	end
	localcs.isgood(localcs.center_sta) = 1;
	% normalize the amplitude
%	normamp = localcs.amps(localcs.center_sta);
	normamp = median(localcs.amps);
	for ista = 1:length(nbstaids)
		localcs.amps(ista) = localcs.amps(ista)./normamp;
	end

	% find the CS measurements that contains the center station
	for ics = 1:length(eventcs.CS)
		if eventcs.CS(ics).sta1==center_sta 
			staid = find(nbstaids==eventcs.CS(ics).sta2);
			localcs.dtps(staid) = -eventcs.CS(ics).dtp(ip);
			localcs.isgood(staid) = eventcs.CS(ics).isgood(ip);
			localcs.cohere(staid) = eventcs.CS(ics).cohere(ip);
		elseif eventcs.CS(ics).sta2==center_sta 
			staid = find(nbstaids==eventcs.CS(ics).sta1);
			localcs.dtps(staid) = eventcs.CS(ics).dtp(ip);
			localcs.isgood(staid) = eventcs.CS(ics).isgood(ip);
			localcs.cohere(staid) = eventcs.CS(ics).cohere(ip);
		end
	end % ics
	ind = find(localcs.isgood<0);
	localcs.stlas(ind) = [];
	localcs.stlos(ind) = [];
	localcs.dists(ind) = [];
	localcs.amps(ind) = [];
	localcs.dtps(ind) = [];
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
	ddists = localcs.dists - localcs.dists(center_sta);
	localcs.dtps = corr_cycle_skip(ddists,localcs.dtps,periods(ip),localcs.OPW_v);
	localcs.Wph = Wph;
	event_data(ie) = localcs;
end

%%% check bad events
% remove events with too few cs measurements
for ie=1:length(event_data)
	csnum(ie) = length(event_data(ie).isgood);
end
ind = find(csnum < min_cs_num);
for ie=1:length(ind)
	disp(['delete event: ',event_data(ind(ie)).id,' for too few cs']);
end
event_data(ind) = [];

% correct for azimuthal variation of amplitude due to radiation pattern
event_data = correct_amp_azi(event_data);

% remove events with large errors
v1 = v1_0;
v2 = v2_0;
event_parastr = fit_event_para(v1,v2,event_data);
for ie=1:length(event_data)
	errs = TPW_err_array(event_parastr(ie),event_data(ie),1); 
	event_data(ie).err = sum(errs.^2)./length(event_data(ie).stlas);
%	pause
end
event_errs = [event_data.err];
ind = find(event_errs > median(event_errs)+std(event_errs));
disp(['Initial model error:',num2str(TPW_vel_err(v1,v2,event_data))]);
for ie=1:length(ind)
	disp(['delete large err event: ',event_data(ind(ie)).id]);
end
%event_data(ind) = [];


% fake Two plane wave measurement
if isfakedata
	event_data = fake_TPW_measure(v1_0,v2_0,event_data);
end

% grid search through velocity space

v1 = v1_0;
v2 = v2_0;

v1_array = linspace(v1*(1-test_r),v1*(1+test_r),test_N);
v2_array = linspace(v2*(1-test_r),v2*(1+test_r),test_N);

for iv1 = 1:length(v1_array)
	for iv2 = 1:length(v2_array)
		v1 = v1_array(iv1);
		v2 = v2_array(iv2);
		if v1<v2
			errmat(iv1,iv2) = NaN;
			continue;
		end
		errmat(iv1,iv2) = TPW_vel_err(v1,v2,event_data);
		disp(sprintf('v1: %f, v2: %f: %f',v1,v2,errmat(iv1,iv2)));
	end
end
			
figure(39)
clf
[xi yi] = ndgrid(v1_array,v2_array);
surface(xi,yi,errmat);
xlim([v1_array(1) v1_array(end)]);
ylim([v2_array(1) v2_array(end)]);
colorbar
caxis([nanmin(errmat(:)) nanmedian(errmat(find(errmat<nanmedian(errmat(:)))))]);

ind = find(errmat==min(errmat(:)));
disp(sprintf('%f %f',xi(ind),yi(ind)));
stemp = sprintf('v1: %f v2: %f, v10: %f v2_0 %f',xi(ind),yi(ind),v1_0,v2_0);
title(stemp)
save2pdf(['errmat_',num2str(ip),'.pdf'],gcf,100);

v1 = xi(ind); v2=yi(ind);

v1=v1_0;
v2=v2_0;

TPW_vel_err(v1,v2,event_data)

event_parastr = fit_event_para(v1,v2,event_data);
for ie=1:length(event_data)
	errs = TPW_err_array(event_parastr(ie),event_data(ie),1); 
	event_data(ie).err = sum(errs.^2);
%%	pause
end

filename = ['workspace_ip_',num2str(ip)];
save(filename);

end   %loop of ip
