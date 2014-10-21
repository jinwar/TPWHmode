
% this program try to follow the two plane wave method, except for 
% allowing the phase velocity of the second wave different from the
% first one

clear;

setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
ip = 3;
sample_range = parameters.maxstadist;
center_stla = 9.3887;
center_stlo = -144.88;
v1_0 = 4.7;
v2_0 = 5.57;
r = 0.05;
v1_0 = 5.0;
v2_0 = 4.65;

init_num = 0;
for A2 = [0.25 0.5 0.75 1]
	for ang = [0 pi/2 pi 1.5*pi]
		init_num = init_num + 1;
		w2_para0(init_num,:) = [A2 ang];
	end
end
%w2_para0 = [0.25 pi];
w2_paraL = [0 0];
w2_paraU = [1 2*pi];
test_N = 10;

csmatfiles = dir(fullfile('CSmeasure',['*_cs_',parameters.component,'.mat']));

% Gather event information 
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
	localcs.dists = deg2km(distance(localcs.stlas,localcs.stlos,evla,evlo));
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
	normamp = max(localcs.amps);
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
	ind = find(localcs.isgood<0);
	localcs.stlas(ind) = [];
	localcs.stlos(ind) = [];
	localcs.dists(ind) = [];
	localcs.amps(ind) = [];
	localcs.dtps(ind) = [];
	localcs.isgood(ind) = [];
	stadists = distance(localcs.stlas,localcs.stlos,center_stla,center_stlo);
	[temp center_sta] = min(stadists);
	localcs.center_sta = center_sta;

	event_data(ie) = localcs;
end

% fake Two plane wave measurement
%event_data = fake_TPW_measure(event_data);

% Get a good average velocity model
[vel_para,resnorm,residual] = lsqnonlin(@(para) OPW_vel_sumerr(para,event_data),v1_0);
vel_para

% inverse again by allowing off great-circle arrival angle
para0(1) = vel_para(1);
for ie=1:length(event_data)
	center_la = event_data(ie).stlas(event_data(ie).center_sta);
	center_lo = event_data(ie).stlos(event_data(ie).center_sta);
	theta1 = azimuth(center_la,center_lo,event_data(ie).evla,event_data(ie).evlo)+180;
	para0(ie+1) = theta1;
	event_data(ie).gcazi = theta1;
	vel_para = lsqnonlin(@(para) OPW_vel_err_array(para,event_data(ie)),[v1_0,theta1]);
	event_data(ie).OPW_v1 = vel_para(1);
	event_data(ie).OPW_phi = vel_para(2);
end

[vel_para,resnorm,residual] = lsqnonlin(@(para) OPW_vel_sumerr(para,event_data),para0);

for ie=1:length(event_data)
	event_data(ie).azi = vel_para(ie+1);
	event_data(ie).v1 = vel_para(1);
	event_data(ie).v2 = v2_0;
end
OPW_vel_sumerr(vel_para,event_data);
%keyboard

%A2_array = 0:0.05:1.5;
%phi2_array = 0:0.1:2*pi;
%[xi yi] = ndgrid(A2_array,phi2_array);
%for i = 1:length(A2_array)
%	for j = 1:length(phi2_array)
%	para = [A2_array(i) phi2_array(j)];
%	errmat(i,j) = sum(TPW_err(para,event_data(1)).^2);
%	end
%end
%figure(38)
%clf
%surface(xi,yi,errmat);
%colormap('grey');
%shading flat
%colorbar

% invert the parameter for the second wave
for ie=1:length(event_data)
	opts1=  optimset('display','off');
	minerr = 1e9;
	disp(event_data(ie).id);
	for itest = 1:init_num
		[para,resnorm,residual] = lsqnonlin(@(para) TPW_err(para,event_data(ie)),w2_para0(itest,:),w2_paraL,w2_paraU,opts1);
		if minerr > resnorm
			minerr = resnorm;
			min_para = para;
			disp(sprintf('Find smaller error: %f %f: %f',para(1),para(2),resnorm));
		end
	end
	para = min_para;
%	disp([num2str(para),' ',num2str(residual)]);
	event_data(ie).A2 = para(1);
	event_data(ie).phi2 = para(2);
end
% setup the parameters
for ie=1:length(event_data)
	event_parastr(ie) = para2str([event_data(ie).A2 event_data(ie).phi2],event_data(ie));
end
para0 = [event_data(1).v1 event_data(1).v2];
[vel_para,resnorm,residual] = lsqnonlin(@(para) TPW_vel_err(para,event_parastr,event_data),para0);
vel_para

v1 = v1_0;
v2 = v2_0;

v1_array = linspace(v1*0.9,v1*1.1,test_N);
v2_array = linspace(v2*0.9,v2*1.1,test_N);

%h = waitbar(0,'grid searching')
for iv1 = 1:length(v1_array)
%	waitbar(iv1/length(v1_array),h);
	for iv2 = 1:length(v2_array)
		for ie=1:length(event_data)
			event_data(ie).v1 = v1_array(iv1);
			event_data(ie).v2 = v2_array(iv2);
			opts1=  optimset('display','off');
			minerr = 1e9;
			for itest = 1:init_num
				[para,resnorm,residual] = lsqnonlin(@(para) TPW_err(para,event_data(ie)),w2_para0(itest,:),w2_paraL,w2_paraU,opts1);
				if minerr > resnorm
					minerr = resnorm;
					min_para = para;
				end
			end
			para = min_para;
		%	disp([num2str(para),' ',num2str(residual)]);
			event_data(ie).A2 = para(1);
			event_data(ie).phi2 = para(2);
			event_parastr(ie) = para2str(para,event_data(ie));
		end
		errarray = TPW_vel_err([v1_array(iv1) v2_array(iv2)],event_parastr,event_data);
		errmat(iv1,iv2) = sum(errarray.^2);
		disp(sprintf('v1: %f, v2: %f: %f',v1_array(iv1),v2_array(iv2),errmat(iv1,iv2)));
	end
end
			
figure(39)
clf
[xi yi] = ndgrid(v1_array,v2_array);
surface(xi,yi,errmat);
colorbar
caxis([0 median(errmat(find(errmat<median(errmat(:)))))]);

ind = find(errmat==min(errmat(:)));
disp(sprintf('%f %f',xi(ind),yi(ind)));

v1 = xi(ind); v2=yi(ind);

%v1=v2_0;
%v2=v1_0;

for ie=1:length(event_data)
	event_data(ie).v1 = v1;
	event_data(ie).v2 = v2;
	opts1=  optimset('display','off');
	minerr = 1e9;
	for itest = 1:init_num
		[para,resnorm,residual] = lsqnonlin(@(para) TPW_err(para,event_data(ie)),w2_para0(itest,:),w2_paraL,w2_paraU,opts1);
		if minerr > resnorm
			minerr = resnorm;
			min_para = para;
		end
	end
	para = min_para;
%	disp([num2str(para),' ',num2str(residual)]);
	event_data(ie).A2 = para(1);
	event_data(ie).phi2 = para(2);
	event_parastr(ie) = para2str(para,event_data(ie));
end

for ie=1:length(event_data)
	errs = TPW_err_array(event_parastr(ie),event_data(ie),1); 
	sum(errs.^2);
	pause;
end
