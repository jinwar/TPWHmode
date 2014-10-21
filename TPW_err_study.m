
clear;

ip = 6;

filename = ['workspace_ip_',num2str(ip),'.mat'];
load(filename);

figure(39)
clf
[xi yi] = ndgrid(v1_array,v2_array);
surface(xi,yi,errmat);
xlim([v1_array(1) v1_array(end)]);
ylim([v2_array(1) v2_array(end)]);
colorbar
caxis([nanmin(errmat(:)) nanmedian(errmat(find(errmat<nanmedian(errmat(:)))))]);

v1 = v1_0; v2=v2_0;
event_parastr_0 = fit_event_para(v1,v2,event_data);
v1 = xi(ind); v2=yi(ind);
event_parastr = fit_event_para(v1,v2,event_data);
for ie = 1:length(event_data)
	errs = TPW_err_array(event_parastr_0(ie),event_data(ie)); 
	event_err_0(ie) = sum(errs.^2);
	event_amperr_0(ie) = sum(errs(1:end/2).^2);
	event_pherr_0(ie) = sum(errs(end/2+1:end).^2);
	errs = TPW_err_array(event_parastr(ie),event_data(ie)); 
	event_err(ie) = sum(errs.^2);
	event_amperr(ie) = sum(errs(1:end/2).^2);
	event_pherr(ie) = sum(errs(end/2+1:end).^2);
	para = polyfit(event_data(ie).dists,event_data(ie).dtps,1);
	polyerrs = polyval(para,event_data(ie).dists)-event_data(ie).dtps;
	disp(sprintf('ie: %d, stanum: %d, init: %f final: %f , polyerr: %f',ie,length(event_data(ie).dtps),event_err_0(ie),event_err(ie),sum(polyerrs.^2)));
end

disp(['init: ',num2str(sum(event_err_0)),' final: ',num2str(sum(event_err))]);
ie = 5;

errs = TPW_err_array(event_parastr(ie),event_data(ie),48); 
errs = TPW_err_array(event_parastr_0(ie),event_data(ie),49); 
