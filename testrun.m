
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exam the spectrum of modes
%
rehash;
eventfile = '201208270437_LHT.mat';
eventfile = '201209051442_LHT.mat';
fund = load(fullfile('../fake_sta_00/eventmat',eventfile));
first = load(fullfile('../fake_sta_11/eventmat',eventfile));
mix = load(fullfile('../fake_sta_01/eventmat',eventfile));

frange = [0.018 0.022];
staid = 9;
plot_spectrum(fund.event.stadata(staid).data,1,frange,1);
title('fund')
plot_spectrum(first.event.stadata(staid).data,1,frange,2);
title('first')
plot_spectrum(mix.event.stadata(staid).data,1,frange,3);
title('mix')
%plot_spectrum(first.event.stadata(staid).data+fund.event.stadata(staid).data,1,frange,4);
%title('fund+first')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exam the relation between A and OPW_v
%
if 0
for ie=1:length(event_data)
	para = polyfit(event_data(ie).amps,event_data(ie).dists,1);
	A_grad(ie) = para(1);
	disp([num2str(A_grad(ie)),' ',num2str(event_data(ie).OPW_v)]);
end
figure(48)
clf
plot(A_grad,[event_data.OPW_v],'o');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the comparison between initial model and final model
if 0

v1 = v1_0; v2=v2_0;
event_parastr_0 = fit_event_para(v1,v2,event_data);
v1 = xi(ind); v2=yi(ind);
event_parastr = fit_event_para(v1,v2,event_data);
for ie = 1:length(event_data)
	errs = TPW_err_array(event_parastr_0(ie),event_data(ie)); 
	event_err_0(ie) = sum(errs.^2);
	errs = TPW_err_array(event_parastr(ie),event_data(ie)); 
	event_err(ie) = sum(errs.^2);
	para = polyfit(event_data(ie).dists,event_data(ie).dtps,1);
	polyerrs = polyval(para,event_data(ie).dists)-event_data(ie).dtps;
	disp(sprintf('ie: %d, stanum: %d, init: %f final: %f, polyerr: %f',ie,length(event_data(ie).dtps),event_err_0(ie),event_err(ie),sum(polyerrs.^2)));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to plot all the event's data with a linear fit

%figure(1)
%clf
%hold on
%for ie = 1:length(event_data)
%	color = rand(1,3);
%	dists = event_data(ie).dists;
%	dists = dists-mean(dists);
%	plot(dists,event_data(ie).dtps,'x','color',color);
%	para = polyfit(dists,event_data(ie).dtps,1);
%	plot(dists,polyval(para,dists),'color',color);
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This section is for grid search for the error surface of the event-based parameter inversion of A2 and phi
%
%

%A2_array = linspace(0,1,10);
%phi2_array = linspace(0,2*pi,20);
%
%for iA2 = 1:length(A2_array)
%	for iphi2 = 1:length(phi2_array)
%		event_parastr(ie).A2 = A2_array(iA2);
%		event_parastr(ie).phi2 = phi2_array(iphi2);
%		errs = TPW_err_array(event_parastr(ie),event_data(ie),1);
%		errmat_ie(iA2,iphi2) = sum(errs.^2);
%		pause
%	end
%end
%
%figure(58)
%clf
%[xi yi] = ndgrid(A2_array,phi2_array);
%surface(xi,yi,errmat_ie);
%colorbar
%
%ind = find(errmat_ie == min(errmat_ie(:)));
%event_parastr(ie).A2 = xi(ind);
%event_parastr(ie).phi2 = yi(ind);
%errs = TPW_err_array(event_parastr(ie),event_data(ie),1);
