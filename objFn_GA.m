function error=objFn_GA(P)
global  Time Mean errMat Pars_all

a0=[0.0870 0.0870 2.7427 2.7427 30.3761 0.0258 1.9425];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[~,a]=ode15s(@(t,a)Single_Cell_Model_GA(t,a,P),Time,a0,options);

multiple=a(:,6);
% PDFMean = Mean / sum(Mean);
% PDFmultiple = multiple / sum(multiple);
% error = calculate_kldiv(PDFMean,PDFmultiple);
rmse=rms(multiple-Mean(:,1));
error=mean(rmse);
% error = sum(Mean .* log(Mean ./ multiple), 'omitnan');
% errMat = [errMat;error];
% Pars_all = [Pars_all;P];
% disp(error)
end