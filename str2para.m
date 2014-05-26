function para = str2para(parastr,pnum)

if pnum == 3
	para(1) = parastr.A1;
	para(2) = parastr.A2;
	para(3) = parastr.phi2;
end
	
if pnum == 4
	para(1) = parastr.phi1;
 	para(2) = parastr.A1;
	para(3) = parastr.phi2;
	para(4) = parastr.A2;
end

if pnum == 5
	para(1) = parastr.A1;
	para(2) = parastr.A2;
	para(3) = parastr.phi2;
	para(4) = parastr.v1;
	para(5) = parastr.v2;
end
