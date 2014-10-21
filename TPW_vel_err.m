function sum_errs = TPW_vel_err(v1,v2,event_data)
% error function for the velocity inversion.

errs_all = [];

event_parastr = fit_event_para(v1,v2,event_data);

for ie = 1:length(event_data)
	errs = TPW_err_array(event_parastr(ie),event_data(ie));
	errs_all = [errs_all(:); errs(:)];
end

sum_errs = sum(errs_all.^2);

%residual = sum(errs_all.^2);
