function errs = TPW_err_array(parastr,localcs,isfigure);

if ~exist('isfigure','var')
	isfigure = 0;
end

setup_parameters;

v1 = parastr.v1;
phi1 = parastr.phi1;
A1 = parastr.A1;
v2 = parastr.v2;
phi2 = parastr.phi2;
A2 = parastr.A2;
Wph = localcs.Wph;

center_sta = localcs.center_sta;
stlas = localcs.stlas;
stlos = localcs.stlos;
dists = localcs.dists;
center_la = stlas(center_sta);
center_lo = stlos(center_sta);
baz = azimuth(center_la,center_lo,localcs.evla,localcs.evlo);
Taz = baz - 90;
T = localcs.period;

k1 = 2*pi/v1/T;
k2 = 2*pi/v2/T;

i = sqrt(-1);
wave1 = A1.*exp(i*(dists.*k1+phi1));
wave2 = A2.*exp(i*(dists.*k2+phi2));

%phi1 phi2
final_wave = wave1 + wave2;

% amplitude misfit
A_pre = abs(final_wave);
A_obs = localcs.amps;
A_pre = A_pre./sum(A_pre).*sum(A_obs);
dA = A_pre(:) - A_obs(:);
%dA = dA./sum(abs(dA));

% phase misfit
phi_pre = angle(final_wave);
dphi_pre = phi_pre - phi_pre(center_sta);
%disp(phi_pre(end))
%disp(phi_pre(center_sta))
%disp(dphi_pre(end))
dphi_pre = wrapTo2Pi(dphi_pre);
dphi_obs = localcs.dtps/T*2*pi;
dphi_obs = wrapTo2Pi(dphi_obs);
dph = (dphi_pre(:) - dphi_obs(:));
dph = wrapTo2Pi(dph);
ind = find(dph>pi);
if ~isempty(ind)
	dph(ind) = dph(ind) - 2*pi;
end
dph = dph./pi;
%dph = dph./sum(abs(dph));
errs = [dA(:); dph(:)*Wph];

if isfigure
	figure(isfigure)
	clf
	subplot(1,3,1)
	hold on
	plot(localcs.dists,A_pre,'x')
	plot(localcs.dists,A_obs,'o')
	if isfield(localcs,'ori_amps')
		plot(localcs.dists,localcs.ori_amps,'ro')
	end
	subplot(1,3,2)
	hold on
	ddists = localcs.dists-localcs.dists(center_sta);
	phase_var = localcs.dtps-ddists/v1;
	plot(localcs.dists,phase_var,'o')
	para = polyfit(localcs.dists(:),localcs.dtps(:),1);
%	plot(localcs.dists,polyval(para,localcs.dists));
	plot(localcs.dists,phase_var+dph/2*T,'x')
	if isfield(localcs,'ori_dtps')
		ori_phase_var = localcs.ori_dtps-ddists/v1;
		plot(localcs.dists,ori_phase_var,'ro')
	end
	title(['OPW vel: ',num2str(1./para(1))]);
	subplot(1,3,3)
	hold on
	polyerrs = (polyval(para,localcs.dists)-localcs.dtps)/T*2*Wph;
	plot(errs,'x')
	plot([length(dA)+1:length(errs)],polyerrs,'rx');
	title(num2str(sum(errs.^2)));
end

