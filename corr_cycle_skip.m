function cor_tp = corr_cycle_skip(ddist,tp,period,refv)

predict_tp = ddist./refv;

for i=1:length(tp)
	while tp(i)-predict_tp(i) > period
		tp(i) = tp(i)-period;
	end
	while predict_tp(i) - tp(i) > period
		tp(i) = tp(i)+period;
	end
end

cor_tp = tp;
