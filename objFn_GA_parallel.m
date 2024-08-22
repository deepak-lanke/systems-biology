function error=objFn_GA_parallel(P, odeFcn, Time, a0, Mean)
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[~,a]=ode15s(@(t,a) odeFcn(t, a, P),Time,a0,options);

multiple=a(:,6);
% error = calculate_kldiv(multiple,Mean);
rmse=rms(multiple-Mean(:,1));
error=mean(rmse);
end