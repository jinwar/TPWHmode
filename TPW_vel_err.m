function sum_errs = TPW_vel_err(v1,v2,event_data)
% error function for the velocity inversion.

event_parastr = fit_event_para(v1,v2,event_data);

for ie = 1:length(event_data)
	errs = TPW_err_array(event_parastr(ie),event_data(ie));
	errs_all(ie) = sum(errs.^2);
end

sum_errs = errs_all(:);

%residual = sum(errs_all.^2);
