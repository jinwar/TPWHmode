function event_data_out = correct_amp_azi(event_data)
% by assuming the the amplitude follows a form as A = ax^2+bx+c+dy
% where x is epicentral distance and y is azimuth of 
% the event to the station. And correct the amplitude to 
% the central azimuth of the array.
% written by Ge Jin, jinwar@gmail.com

isfigure = 0;
is_corr_phase = 0;


event_data_out = event_data;
for ie=1:length(event_data)
	stlas = event_data(ie).stlas;
	stlos = event_data(ie).stlos;
	evla = event_data(ie).evla;
	evlo = event_data(ie).evlo;
	center_la = stlas(event_data(ie).center_sta);
	center_lo = stlos(event_data(ie).center_sta);
	amps = event_data(ie).amps;
	dtps = event_data(ie).dtps;
	[dists azis] = distance(evla,evlo,stlas,stlos);
	ft = fittype('a*dists.^2+b*dists+c+d*azis','independent',{'dists','azis'});
	para_amp = fit([dists(:),azis(:)],amps(:),ft);
	cent_azis = azis;
	cent_azis(:) = azimuth(evla,evlo,center_la,center_lo);
	post_amps = feval(para_amp,[dists(:),cent_azis(:)]);
	event_data_out(ie).amps = post_amps;
	event_data_out(ie).ori_amps = event_data(ie).amps;

	if is_corr_phase
		para_dtp = fit([dists(:),azis(:)],dtps(:),ft);
		cent_azis = azis;
		cent_azis(:) = azimuth(evla,evlo,center_la,center_lo);
		post_dtps = feval(para_dtp,[dists(:),cent_azis(:)]);
		event_data_out(ie).dtps = post_dtps;
		event_data_out(ie).ori_dtps = event_data(ie).dtps;
	end


	if isfigure
		figure(58)
		clf
		subplot(1,2,1)
		hold on
		pre_amps = feval(para_amp,[dists(:),azis(:)]);
		cent_azis = azis;
		cent_azis(:) = mean(azis(:));
		post_amps = feval(para_amp,[dists(:),cent_azis(:)]);
		plot(dists,pre_amps,'x');
		plot(dists,amps,'o');
		plot(dists,post_amps,'rx');
		if is_corr_phase
		subplot(1,2,2)
		hold on
		pre_dtps = feval(para_dtp,[dists(:),azis(:)]);
		cent_azis = azis;
		cent_azis(:) = mean(azis(:));
		post_dtps = feval(para_dtp,[dists(:),cent_azis(:)]);
		plot(dists,pre_dtps,'x');
		plot(dists,dtps,'o');
		plot(dists,post_dtps,'rx');
		end
		pause
	end
end  % loop of event 

end % end of function correct_amp_azi

function errs = ampfun(para,dists,azis,amps)
	x = dists;
	y = azis;
	a = para(1); b=para(2); c=para(3); d=para(4);
	pre_amps = a*x.^2+b*x+c+d*y;
	errs = amps-pre_amps;
end
