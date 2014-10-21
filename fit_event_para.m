function event_parastr = fit_event_para(v1,v2,event_data)

is_allow_reverse = 1;

	init_num = 0;
for A2 = [0.5]
	for ang = [pi/2 pi 1.5*pi]
		init_num = init_num + 1;
		w2_para0(init_num,:) = [A2 ang];
	end
end
%w2_para0 = [0.25 pi];
w2_paraL = [0 0];
w2_paraU = [0.95 2*pi];

for ie=1:length(event_data)
	event_data(ie).v1 = v1;
	event_data(ie).v2 = v2;
	opts1=  optimset('display','off');
	minerr = 1e9;
	for itest = 1:init_num
		[para,resnorm,residual] = lsqnonlin(@(para) TPW_err(para,event_data(ie)),w2_para0(itest,:),w2_paraL,w2_paraU,opts1);
		if minerr > resnorm
			minerr = resnorm;
			min_para = para;
		end
	end
	is_reverse = 0;
	if is_allow_reverse
		event_data(ie).v2 = v1;
		event_data(ie).v1 = v2;
		for itest = 1:init_num
			[para,resnorm,residual] = lsqnonlin(@(para) TPW_err(para,event_data(ie)),w2_para0(itest,:),w2_paraL,w2_paraU,opts1);
			if minerr > resnorm
				minerr = resnorm;
				min_para = para;
				is_reverse = 1;
			end
		end
	end
	para = min_para;
%	disp([num2str(para),' ',num2str(residual)]);
	event_data(ie).A2 = para(1);
	event_data(ie).phi2 = para(2);
	if ~is_reverse
		event_data(ie).v1 = v1;
		event_data(ie).v2 = v2;
	else
		disp('reversed')
		event_data(ie).v2 = v1;
		event_data(ie).v1 = v2;
	end
	event_parastr(ie) = para2str(para,event_data(ie));
end

