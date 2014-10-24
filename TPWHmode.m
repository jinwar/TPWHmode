
% this program try to follow the two plane wave method, except for 
% allowing the phase velocity of the second wave different from the
% first one

clear;

load pa5mod

isfakedata = 0;
is_azi_corr = 0;
is_cs_dtp = 1;
is_rm_badevent = 0;
max_err_tol = 1e6;

test_r = 0.03;
test_N = 5;
Wph = 30;

for ip=[2 3 5 6]

setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
sample_range = parameters.maxstadist;
center_stla = 9.3887;
center_stlo = -144.88;
periods = parameters.periods;
min_cs_num = 5;

v1_0 = interp1(pa5mod.T_first,pa5mod.phv_first,periods(ip));
v2_0 = interp1(pa5mod.T_fund,pa5mod.phv_fund,periods(ip));

event_data = gather_eventdata(center_stla,center_stlo,ip,Wph);
if ~is_cs_dtp
	for ie=1:length(event_data);
		event_data(ie).dtps = event_data(ie).fft_dtps;
		normamp = median(event_data(ie).amps);
		event_data(ie).amps = event_data(ie).fft_amps./normamp;
	end
end

%%% check bad events
% remove events with too few cs measurements
for ie=1:length(event_data)
	csnum(ie) = length(event_data(ie).isgood);
end
ind = find(csnum < min_cs_num);
for ie=1:length(ind)
	disp(['delete event: ',event_data(ind(ie)).id,' for too few cs']);
end
event_data(ind) = [];

% correct for azimuthal variation of amplitude due to radiation pattern
if is_azi_corr
	event_data = correct_amp_azi(event_data);
end

% remove events with large errors
v1 = v1_0;
v2 = v2_0;
event_parastr = fit_event_para(v1,v2,event_data);
for ie=1:length(event_data)
	errs = TPW_err_array(event_parastr(ie),event_data(ie),1); 
	event_data(ie).err = sum(errs.^2);
	disp(sprintf('id: %s, err: %f, amp diff: %f',event_data(ie).id,event_data(ie).err,event_data(ie).amp_diff));
%	pause
end
event_errs = [event_data.err];
ind = find(event_errs > median(event_errs)+std(event_errs));
for ie=1:length(ind)
	disp(['Found large err event: ',event_data(ind(ie)).id]);
end
%keyboard
if is_rm_badevent
	event_data(ind) = []; disp('removed');
end

ind = find([event_data.err]>max_err_tol);
if ~isempty(ind)
	event_data(ind) = []; disp('large error event removed');
end

disp(['sum error after removal: ',num2str(sum([event_data.err]))]);


% fake Two plane wave measurement
if isfakedata
	event_data = fake_TPW_measure(v1_0,v2_0,event_data);
end

% grid search through velocity space

v1 = v1_0;
v2 = v2_0;

v1_array = linspace(v1*(1-test_r),v1*(1+test_r),test_N);
v2_array = linspace(v2*(1-test_r),v2*(1+test_r),test_N);

errmat_all = zeros(length(event_data),length(v1_array),length(v2_array));

for iv1 = 1:length(v1_array)
	for iv2 = 1:length(v2_array)
		v1 = v1_array(iv1);
		v2 = v2_array(iv2);
		if v1<v2
			errmat_all(:,iv1,iv2) = NaN;
			errmat(iv1,iv2)=NaN;
			continue;
		end
		errmat_all(:,iv1,iv2) = TPW_vel_err(v1,v2,event_data);
		errmat(iv1,iv2) = sum(errmat_all(:,iv1,iv2));
		disp(sprintf('v1: %f, v2: %f: %f',v1,v2,sum(errmat_all(:,iv1,iv2))));
	end
end

figure(39)
clf
[xi yi] = ndgrid(v1_array,v2_array);
surface(xi,yi,errmat);
xlim([v1_array(1) v1_array(end)]);
ylim([v2_array(1) v2_array(end)]);
colorbar
caxis([nanmin(errmat(:)) nanmedian(errmat(find(errmat<nanmedian(errmat(:)))))]);

ind = find(errmat==min(errmat(:)));
disp(sprintf('%f %f',xi(ind),yi(ind)));
stemp = sprintf('v1: %f v2: %f, v10: %f v2_0 %f',xi(ind),yi(ind),v1_0,v2_0);
title(stemp)
save2pdf(['errmat_',num2str(ip),'.pdf'],gcf,100);
%
%v1 = xi(ind); v2=yi(ind);
%
%v1=v1_0;
%v2=v2_0;
%
%TPW_vel_err(v1,v2,event_data)
%
%event_parastr = fit_event_para(v1,v2,event_data);
%for ie=1:length(event_data)
%	errs = TPW_err_array(event_parastr(ie),event_data(ie),1); 
%	event_data(ie).err = sum(errs.^2);
%%%	pause
%end

filename = ['workspace_ip_',num2str(ip)];
save(filename);

end   %loop of ip
