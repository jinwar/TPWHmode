function parastr = para2str(para,localcs)

	evla = localcs.evla;
	evlo = localcs.evlo;

if length(para) == 8
	parastr.v1 = para(1);
	parastr.phi1 = para(2);
	parastr.A1 = para(3);
	parastr.v2 = para(4);
	parastr.phi2 = para(5);
	parastr.A2 = para(6);
	parastr.theta1 = para(7);
	parastr.theta2 = para(8);
end

if length(para) == 7
	parastr.v1 = para(1);
	parastr.A1 = para(2);
	parastr.v2 = para(3);
	parastr.phi2 = para(4);
	parastr.A2 = para(5);
	parastr.theta1 = para(6);
	parastr.theta2 = para(7);
	parastr.phi1 = 0;
	return;
end
	
if length(para) == 6
	parastr.v1 = para(1);
	parastr.phi1 = para(2);
	parastr.A1 = para(3);
	parastr.v2 = para(4);
	parastr.phi2 = para(5);
	parastr.A2 = para(6);
	center_la = localcs.stlas(localcs.center_sta);
	center_lo = localcs.stlos(localcs.center_sta);
	[epi_dist baz] = distance(center_la,center_lo,evla,evlo);
	parastr.theta1 = baz+180;
	parastr.theta2 = baz+180;
	return;
end

if length(para) == 5
	parastr.A1 = para(1);
	parastr.A2 = para(2);
	parastr.phi2 = para(3);
	parastr.v1 = para(4);
	parastr.v2 = para(5);
	center_la = localcs.stlas(localcs.center_sta);
	center_lo = localcs.stlos(localcs.center_sta);
	[epi_dist baz] = distance(center_la,center_lo,evla,evlo);
	parastr.theta1 = baz+180;
	parastr.theta2 = baz+180;
	parastr.phi1 = 0;
	return;
end

%if length(para) == 4
%	parastr.phi1 = para(1);
%	parastr.A1 = para(2);
%	parastr.phi2 = para(3);
%	parastr.A2 = para(4);
%	center_la = localcs.stlas(localcs.center_sta);
%	center_lo = localcs.stlos(localcs.center_sta);
%	[epi_dist baz] = distance(center_la,center_lo,evla,evlo);
%	parastr.theta1 = baz+180;
%	parastr.theta2 = baz+180;
%	parastr.v1 = localcs.v1;
%	parastr.v2 = localcs.v2;
%	return;
%end

if length(para) == 3
	parastr.A1 = para(1);
	parastr.A2 = para(2);
	parastr.phi2 = para(3);
	parastr.phi1 = 0;
	center_la = localcs.stlas(localcs.center_sta);
	center_lo = localcs.stlos(localcs.center_sta);
	[epi_dist baz] = distance(center_la,center_lo,evla,evlo);
	parastr.theta1 = baz+180;
	parastr.theta2 = baz+180;
	parastr.v1 = localcs.v1;
	parastr.v2 = localcs.v2;
	return;
end
