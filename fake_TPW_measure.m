function outevents = fake_TPW_measure(v1,v2,event_data)


for ie=1:length(event_data)
	localcs = event_data(ie);
	center_sta = localcs.center_sta;
	stlas = localcs.stlas;
	stlos = localcs.stlos;
	dists = localcs.dists;
	center_la = stlas(center_sta);
	center_lo = stlos(center_sta);
	baz = azimuth(center_la,center_lo,localcs.evla,localcs.evlo);
	Taz = baz - 90;
	T = localcs.period;
	theta1 = baz+180;
	theta2 = baz+180;
	A1 = 1;
	A2 = 2*(rand);
	phi1 = 0;
	phi2 = 2*pi*(rand);
	localcs.synA2 = A2;
	localcs.synph2 = phi2;

% transfer the coordinator
	
k1 = 2*pi/v1/T;
k2 = 2*pi/v2/T;

i = sqrt(-1);
wave1 = A1.*exp(i*(dists.*k1+phi1));
wave2 = A2.*exp(i*(dists.*k2+phi2));

	final_wave = wave1 + wave2;

	amps = abs(final_wave);
	dtps = angle(final_wave)/2/pi*T;
	dtps = dtps - dtps(center_sta);
	localcs.amps = amps;
	localcs.dtps = dtps;
	
	outevents(ie) = localcs;

end  % event loop:1

