function errs = TPW_comp(para,localcs);

setup_parameters;
parastr = para2str(para,localcs);
v1 = parastr.v1;
phi1 = parastr.phi1;
A1 = parastr.A1;
v2 = parastr.v2;
phi2 = parastr.phi2;
A2 = parastr.A2;
theta1 = parastr.theta1;
theta2 = parastr.theta2;

center_sta = localcs.center_sta;
stlas = localcs.stlas;
stlos = localcs.stlos;
center_la = stlas(center_sta);
center_lo = stlos(center_sta);
baz = azimuth(center_la,center_lo,localcs.evla,localcs.evlo);
Taz = baz - 90;
T = localcs.period;

% transfer the coordinator
stxs = deg2km(stlas-center_la);
stys = deg2km((stlos-center_lo)*cosd(center_la));

kx1 = 2*pi/v1/T*cosd(theta1);
ky1 = 2*pi/v1/T*sind(theta1);
kx2 = 2*pi/v2/T*cosd(theta2);
ky2 = 2*pi/v2/T*sind(theta2);

i = sqrt(-1);
wave1 = A1.*exp(i*(stxs.*kx1+stys.*ky1+phi1))*cosd(Taz-theta1-90);
wave2 = A2.*exp(i*(stxs.*kx2+stys.*ky2+phi2))*cosd(Taz-theta2-90);

final_wave = wave1 + wave2;

% amplitude misfit
A_pre = abs(final_wave);
dA = (A_pre - localcs.amps)./localcs.amps;
%dA = dA./sum(abs(dA));

% phase misfit
phi_pre = angle(final_wave);
dphi_pre = phi_pre - phi_pre(center_sta);
dphi_pre = wrapTo2Pi(dphi_pre);
dphi_obs = localcs.dtps/T*2*pi;
dphi_obs = wrapTo2Pi(dphi_obs);
dph = (dphi_pre - dphi_obs)/2/pi;
%dph = dph./sum(abs(dph));
Wph = 1;
errs = [dA(:); dph(:)*Wph];

figure(28)
clf
subplot(1,2,1)
hold on
plot(localcs.dists,dphi_obs,'x')
plot(localcs.dists,dphi_pre,'ro')
title('phase')
subplot(1,2,2)
hold on
plot(localcs.dists,localcs.amps,'x')
plot(localcs.dists,A_pre,'ro')
title('amplitude')
figure(29)
clf
subplot(1,2,2)
hold on
plot(localcs.dists,dA,'rx')
title('amplitude')
subplot(1,2,1)
hold on
plot(localcs.dists,dph,'rx')
title('phase')
