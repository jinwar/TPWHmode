function sumerr = TPW_err(para,localcs);

parastr = para2str(para,localcs);
errs = TPW_err_array(parastr,localcs);

sumerr = sum(errs.^2);
