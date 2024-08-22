clc;clear;close all;
tic
% availableGPUs = gpuDeviceCount("available");
% parpool('Processes',availableGPUs);
% spmd
%     gpuDevice
% end
%% Importing Data
global R0 ka R ktf kts m1 m2 k1 k2 k3 k4 k5 k6 k7 k8 km1 km2 km3 km4 km5 km6 num b br Time Mean
f = 1;
% data = readtable("13_Deepak.csv");
load("project2data.mat")
X = []; Y = []; Area = []; Means = [];
for i = 5:4:86
    X = [X table2array(data(2,i))];
end
num = numel(X);
for i = 6:4:86
    Y = [Y table2array(data(2,i))];
end
Vertices = [X;Y];
for i = 3:4:86
    Area = [Area table2array(data(2,i))];
end
for i = 4:4:86
    Means = [Means table2array(data(:,i))];
end
% Normalized means
NMeans = rescale(Means,'InputMin',min(Means),'InputMax',max(Means));
Time = table2array(data(:,1));

P=[0.15	3 0.465	0.5	0.02 0.125	0.875	0.0001	0.6894 4.525	65.352...
	0.05 0.0048951 183.6	4.85	0.788	0.0450	14.545	0.670...
    0.02511 0.05];
a0=[0.0870    0.0870    2.7427    2.7427   30.3761    0.0258    1.9425];
%% Multiple cell model beta & beta' calculations
dist = pdist(Vertices','euclidean');
dist = squareform(dist);
figure(f)
heatmap(dist,"Colormap",jet)
title('Distances')
min_distance=min(min(dist));
max_distance=max(max(dist));
dist=(dist-min_distance)/(max_distance-min_distance);
kdiffin = zeros(num,num);
kdiffout_cells = zeros(num,num);
for i = 1:1:num
    for j = 1:1:num
        if(0<dist(i,j))&&(dist(i,j)<0.25)
            kdiffin(i,j) =0.3;
            kdiffout_cells(i,j)=0.3*0.7;
        elseif (0.25<dist(i,j))&&(dist(i,j)<0.5)
            kdiffin(i,j) =.25;
            kdiffout_cells(i,j) =.25*0.7;
        elseif(0.5<dist(i,j))&&(dist(i,j)<0.75)
            kdiffin(i,j) =.2;
            kdiffout_cells(i,j)=.2*0.7;
        elseif(0.75<dist(i,j))&&(dist(i,j)<=1)
            kdiffin(i,j) =0.15;
            kdiffout_cells(i,j)=0.15*0.7;
        end
    end
end
b = kdiffin; br = kdiffout_cells;
figure(f+1)
heatmap(b,"Colormap",jet);
title('Beta')
figure(f+2)
heatmap(br,"Colormap",jet);
title('Beta dash')
figure(f+3)
scatter(X,Y,500,'filled')
% colorbar;
for i = 1:length(X)
    text(X(i), Y(i), num2str(i), 'Color', 'red', 'FontSize', 8);
end
%% Before GA multiple cell ode solving and plotting
H0 = [];
for i = 1:1:num
    H0 = [H0 a0];
    % H0(7*i-1) = NMeans(1,i);
end
[t, H] = ode15s(@(t,a)HighPackingModel(t,a,P), Time, H0);

figure('Name',"Before GA");
for i = 1:1:num
    subplot(3,7,i)
    plot(t,H(:,7*i-1),'LineWidth',2)
    hold on
    plot(Time,NMeans(:,i),'LineWidth',2)
    hold off
    xlabel('time')
    ylabel('Ca^{+2}')
    % legend("Simulated","Experimental")
    tit = sprintf('Cell %d',i);
    title(tit)
end
legend("Simulated","Experimental")
%% GA
rng default % For reproducibility
X = 0.05;
Y = 3;
LL=P*X;
UL=P*Y;
% OP = [];
% for i = 1:1:num
%     Mean = NMeans(:,i);
%     options = optimoptions('ga','PopulationSize',4,'MaxGenerations',10, 'TolFun', 1e-6, 'TolCon', 1e-6);%, 'FunctionTolerance', 1e-12);
%     [fpar, fval] = ga(@objFn_GA,21,[],[],[],[],LL,UL,[],options);
%     i
%     col = sprintf("Cell %1d",i);
%     OP = [OP,col,fpar];
% end
%% Parallel computing
% Define the ODE function handle outside of the loop
odeFcn = @(t,a,P) Single_Cell_Model_GA(t, a, P);

% Initialize OP matrix
OP = [];
% Iterate over the indices using a parfor loop
parfor i = 1:num
    Mean = NMeans(:,i);
    options = optimoptions('ga','PopulationSize',20,'MaxGenerations',100, 'TolFun', 1e-6, 'TolCon', 1e-6);
    % Store the ODE function handle and its parameters in temporary variables
    odeFcn_temp = odeFcn;
    Time_temp = Time;
    a0_temp = a0;
    % Define the objective function handle
    objFn_handle = @(P) objFn_GA_parallel(P, odeFcn_temp, Time_temp, a0_temp, Mean);
    % Call GA optimization
    [fpar, fval] = ga(objFn_handle, 21, [], [], [], [], LL, UL, [], options);
    % sprintf('Cell %d fitness = %f',i,fval)
    % Store the result
    col = sprintf("Cell %1d",i)
    OP = [OP; fpar];
end
%% After GA Multiple cell ode solving and plotting
[t, H_GA] = ode15s(@(t,a)HighPackingModel_GA(t,a,OP), Time, H0);

figure('Name','After GA');
for i = 1:1:num
    subplot(3,7,i)
    plot(t,H_GA(:,7*i-1),'LineWidth',2)
    hold on
    plot(Time,NMeans(:,i),'LineWidth',2)
    hold off
    xlabel('time')
    ylabel('Ca^{+2}')
    % legend("Simulated","Experimental")
    tit = sprintf('Cell %d',i);
    title(tit)
end
legend("Simulated","Experimental")
%%
delete(gcp('nocreate'));
toc