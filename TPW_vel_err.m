function errs_all = TPW_vel_err(para,event_parastr,event_data)
% error function for the velocity inversion.

errs_all = [];

for ie = 1:length(event_data)
	event_parastr(ie).v1 = para(1);
	event_parastr(ie).v2 = para(2);
	errs = TPW_err_array(event_parastr(ie),event_data(ie));
	errs_all = [errs_all(:); errs(:)];
end

%residual = sum(errs_all.^2);
