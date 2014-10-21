function errs = OPW_vel_err_array(para,localcs,isfigure)

if ~exist('isfigure','var')
	isfigure = 0;
end

v1 = para(1);

% in the situation that invert the velocity only
		% gather information
	T = localcs.period;
	stlas = localcs.stlas;
	stlos = localcs.stlos;
	center_sta = localcs.center_sta;
	center_la = stlas(center_sta);
	center_lo = stlos(center_sta);
	if length(para)==1
		theta1 = azimuth(center_la,center_lo,localcs.evla,localcs.evlo)+180;
	else
		theta1 = para(2);
	end
	% transfer to xy coor
	kx1 = 2*pi/v1/T*cosd(theta1);
	ky1 = 2*pi/v1/T*sind(theta1);
	stxs = deg2km(stlas-center_la);
	stys = deg2km((stlos-center_lo)*cosd(center_la));
	%calculate the phase misfit
	i = sqrt(-1);
	wave1 = exp(i*(stxs.*kx1+stys.*ky1));
	final_wave = wave1;
	phi_pre = angle(final_wave);
	dphi_pre = phi_pre - phi_pre(center_sta);
	dphi_pre = wrapTo2Pi(dphi_pre);
	dphi_obs = localcs.dtps/T*2*pi;
	dphi_obs = wrapTo2Pi(dphi_obs);
	dph = (dphi_pre - dphi_obs);
	dph = wrapTo2Pi(dph);
	ind = find(dph>pi);
	if ~isempty(ind)
		dph(ind) = dph(ind)-2*pi;
	end
	errs = dph;

if isfigure
	figure(57)
	clf
	subplot(1,2,1)
	hold on
	plot(localcs.dists,dphi_pre,'x')
	plot(localcs.dists,dphi_obs,'o')
	subplot(1,2,2)
	hold on
	plot(localcs.dists,dph,'x')
end
	
