
% this program try to follow the two plane wave method, except for 
% allowing the phase velocity of the second wave different from the
% first one

clear;

setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
ip = 4;
sample_range = parameters.maxstadist;

csmatfiles = dir(fullfile(CSmeasure,['*_cs_',parameters.component,'.mat']));
for ie = 1:length(csmatfiles)
	load CSmeasure/201202021334_cs_BHT.mat
	center_sta = 12;
	evla = eventcs.evla;
	evlo = eventcs.evlo;
	stlas = eventcs.stlas;
	stlos = eventcs.stlos;

	stadists = distance(stlas,stlos,stlas(center_sta),stlos(center_sta));
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
	for ista = 1:length(nbstaids)
		localcs.amps(ista) = eventcs.autocor(nbstaids(ista)).amp(ip)^0.5;
		localcs.dtps(ista) = 0;
		localcs.isgood(ista) = 0;
	end
	localcs.isgood(localcs.center_sta) = 1;
	% normalize the amplitude
	normamp = localcs.amps(localcs.center_sta);
	for ista = 1:length(nbstaids)
		localcs.amps(ista) = localcs.amps(ista)./normamp;
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

	xnode = lalim(1):0.1:lalim(2);
	ynode = lolim(1):0.1:lolim(2);
	[t_surf xi yi] = gridfit(localcs.stlas,localcs.stlos,localcs.dtps,xnode,ynode);
	[amp_surf xi yi] = gridfit(localcs.stlas,localcs.stlos,localcs.amps,xnode,ynode);
	%figure(37)
	%clf
	%subplot(1,2,1)
	%worldmap(lalim,lolim);
	%surfacem(xi,yi,t_surf);
	%colorbar('southoutside');
	%plotm(localcs.stlas,localcs.stlos,'kv');
	%title('Travel time');
	%subplot(1,2,2)
	%worldmap(lalim,lolim);
	%surfacem(xi,yi,amp_surf);
	%colorbar('southoutside');
	%plotm(localcs.stlas,localcs.stlos,'kv');
	%title('amplitude');

	%figure(38)
	%clf
	%subplot(1,2,1)
	%plot(localcs.dists,localcs.dtps,'x')
	%title('Travel time');
	%subplot(1,2,2)
	%plot(localcs.dists,localcs.amps,'x')
	%title('amplitude');

	v1 = 4.65;
	phi1 = 0;
	A1 = (max(localcs.amps)+min(localcs.amps))/2;
	v2 = 4.95;
	phi2 = 0;
	A2 = (max(localcs.amps)-min(localcs.amps))/2;
	theta1 = baz+180;
	theta2 = baz+180;
	para0 = [v1 A1 phi1 v2 phi2 A2 theta1 theta2];

	localcs.v1 = v1;
	localcs.v2 = v2;
	%N = 100;
	%phi_array = linspace(0,2*pi,N);
	%for iy = 1:N
	%	phi2 = phi_array(iy);
	%	para0 = [v1 A1 v2 phi2 A2 theta1 theta2];
	%	errs = TPW_err(para0,localcs);
	%	errmat(iy) = sum(errs.^2);
	%end
	%%figure(58)
	%%clf
	%%plot(phi_array,errmat)
	%%pause
	%
	%phi2 = phi_array(find(errmat == min(errmat)));

	parastr0 = para2str(para0,localcs);
	r = 0.05;
	parastrL.v1 = v1*(1-r);
	parastrL.phi1 = 0;
	parastrL.A1 = 0;
	parastrL.v2 = v2*(1-r);
	parastrL.phi2 = 0;
	parastrL.A2 = 0;
	parastrL.theta1 = theta1-10;
	parastrL.theta2 = theta2-10;

	parastrU.v1 = v1*(1+r);
	parastrU.phi1 = 2*pi;
	parastrU.A1 = max(localcs.amps);
	parastrU.v2 = v2*(1+r);
	parastrU.phi2 = 2*pi;
	parastrU.A2 = max(localcs.amps);
	parastrU.theta1 = theta1+10;
	parastrU.theta2 = theta2+10;

	paraN = 3;
	para0 = str2para(parastr0,paraN);
	paraL = str2para(parastrL,paraN);
	paraU = str2para(parastrU,paraN);
	errs0 = TPW_err(para0,localcs);
	%options = optimset('lsqnonlin');
	%options.MaxFunEvals = 1000;
	%options.PlotFcns = {'@optimplotfval'};
	%[para,resnorm,residual]= lsqnonlin(@(para) TPW_err(para,localcs),para0,paraL,paraU,options);
	[para,residual,exitflag,output] = simulannealbnd(@(para) TPW_err(para,localcs),para0,paraL,paraU);
	TPW_comp(para,localcs);
	parastr = para2str(para,localcs)
	residual
end
