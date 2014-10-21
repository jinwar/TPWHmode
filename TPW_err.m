function sumerr = TPW_err(para,localcs,isfigure);

if ~exist('isfigure','var')
	isfigure = 0;
end

parastr = para2str(para,localcs);
errs = TPW_err_array(parastr,localcs,isfigure);

%sumerr = sum(errs.^2);
sumerr = errs;
