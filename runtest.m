
A_test = 0:0.02:0.3; 
dphi_test = 0:0.1:2*pi;
for iA = 1:length(A_test)
	for iph = 1:length(dphi_test)
		A = A_test(iA);
		dphi = dphi_test(iph);
		errs = TPW_err([A dphi],event_data(ie)); 
		errmat(iA,iph) = sum(errs.^2);
	end
end

[xi yi] = ndgrid(A_test,dphi_test);
		
figure(48)
clf
surface(xi,yi,errmat)
