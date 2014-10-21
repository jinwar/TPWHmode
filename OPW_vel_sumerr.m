function sumerr = OPW_vel_sumerr(para,event_data,isfigure)

if ~exist('isfigure','var')
	isfigure = 0;
end

v1 = para(1);

% in the situation that invert the velocity only
sumerr = [];
for ie=1:length(event_data)
	% gather information
	T = event_data(ie).period;
	stlas = event_data(ie).stlas;
	stlos = event_data(ie).stlos;
	dists = event_data(ie).dists;
	center_sta = event_data(ie).center_sta;
	center_la = stlas(center_sta);
	center_lo = stlos(center_sta);
	if length(para)==1
		theta1 = azimuth(center_la,center_lo,event_data(ie).evla,event_data(ie).evlo)+180;
	else
		theta1 = para(ie+1);
	end

	% transfer to xy coor
	kx1 = 2*pi/v1/T*cosd(theta1);
	ky1 = 2*pi/v1/T*sind(theta1);
	[delta_stlas delta_stlos] = distance_InArray(center_la,center_lo,stlas,stlos);
	stxs = deg2km(delta_stlas);
	stys = deg2km(delta_stlos);
	
	%calculate the phase misfit
	i = sqrt(-1);
%	wave1 = exp(i*(stxs.*kx1+stys.*ky1));
	k1 = 2*pi/v1/T;
	wave1 = exp(i*k1.*dists);
	final_wave = wave1;
	phi_pre = angle(final_wave);
	dphi_pre = phi_pre - phi_pre(center_sta);
	dphi_pre = wrapTo2Pi(dphi_pre);
	dphi_obs = event_data(ie).dtps/T*2*pi;
	dphi_obs = wrapTo2Pi(dphi_obs);
	dph = (dphi_pre - dphi_obs);
	dph = wrapTo2Pi(dph);
	ind = find(dph>pi);
	if ~isempty(ind)
		dph(ind) = dph(ind)-2*pi;
	end
	sumerr = [sumerr; dph(:)];
	if isfigure
		figure(84)
		clf
		subplot(1,2,1)
		hold on
		plot(event_data(ie).dists,dphi_pre,'rx');
		plot(event_data(ie).dists,dphi_obs,'o');
		subplot(1,2,2)
		hold on
		plot(event_data(ie).dists,dph,'rx');
%		keyboard;
	end
end

