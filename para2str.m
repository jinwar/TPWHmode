function parastr = para2str(para,localcs)

if length(para) == 2
	parastr.A1 = 1;
	parastr.phi1 = 0;
	parastr.A2 = para(1);
	parastr.phi2 = para(2);
	parastr.v1 = localcs.v1;
	parastr.v2 = localcs.v2;
end
