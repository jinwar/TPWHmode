function [sta_phv sta_theta] = OPW_inv(eventcs,ip)

clear;
load ./CSmeasure/201112140504_cs_LHT.mat;
ip = 4;

setup_parameters
sample_range = parameters.maxstadist;

evla = eventcs.evla;
evlo = eventcs.evlo;
stlas = eventcs.stlas;
stlos = eventcs.stlos;

for ista = 1:length(stlas)

	center_sta = ista;
	center_stla = stlas(ista);
	center_stlo = stlos(ista);

	stadists = distance(stlas,stlos,center_stla,center_stlo);
	stadists = deg2km(stadists);

	nbstaids = find(stadists < sample_range);

	% reconstruct the data format to local array index
	CSnum = 0;
	localcs = [];
	localcs.stlas = stlas(nbstaids);
	localcs.stlos = stlos(nbstaids);
	localcs.dists = deg2km(distance(localcs.stlas,localcs.stlos,evla,evlo));
	localcs.evla = evla;
	localcs.evlo = evlo;
	localcs.center_sta = find(nbstaids == center_sta);
	localcs.period = parameters.periods(ip);
	for jsta = 1:length(nbstaids)
		localcs.amps(jsta) = eventcs.autocor(nbstaids(jsta)).amp(ip)^0.5;
		localcs.dtps(jsta) = 0;
		localcs.isgood(jsta) = 0;
	end
	localcs.isgood(localcs.center_sta) = 1;
	% normalize the amplitude
	normamp = localcs.amps(localcs.center_sta);
	for jsta = 1:length(nbstaids)
		localcs.amps(jsta) = localcs.amps(jsta)./normamp;
	end

	% find the CS measurements that contains the center station
	for ics = 1:length(eventcs.CS)
		if eventcs.CS(ics).sta1==center_sta 
			staid = find(nbstaids==eventcs.CS(ics).sta2);
			localcs.dtps(staid) = -eventcs.CS(ics).dtp(ip);
			localcs.isgood(staid) = eventcs.CS(ics).isgood(ip);
		elseif eventcs.CS(ics).sta2==center_sta 
			staid = find(nbstaids==eventcs.CS(ics).sta1);
			localcs.dtps(staid) = eventcs.CS(ics).dtp(ip);
			localcs.isgood(staid) = eventcs.CS(ics).isgood(ip);
		end
	end % ics

	% define the initial model
	center_la = localcs.stlas(localcs.center_sta);
	center_lo = localcs.stlos(localcs.center_sta);
	[epi_dist baz] = distance(center_la,center_lo,evla,evlo);

	if length(find(localcs.isgood)) > 5
		para0(1) = parameters.refphv(ip);
		para0(2) = baz+180;
		[para,resnorm,residual,exitflag] = lsqnonlin(@(para) OPW_vel_err_array(para,localcs),para0);
		sta_phv(ista) = para(1);
		sta_theta(ista) = para(2);
	else
		sta_phv(ista) = NaN;
		sta_theta(ista) = NaN;
	end
end
