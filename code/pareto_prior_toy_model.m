% PARETO PRIOR TOY MODEL

% The goal of this script is to generate data from our model with a nominal
% set of parameters where there are two sources of measurements: N(t) and
% phi(t)

% We will fit the model using the artificially generated data from the
% model with some additive noise.

% We will then find the pareto front and use this to obtain a prior on the
% model parameters (say using parameters w/in some distance of the front).

% Eventually validate on newly generated data...
close all; clear all; clc;
%% Load in the current best fitting parameters just for reference
ptest = load('../out/ptest.mat');
ptest = struct2cell(ptest);
ptest = cell2mat(ptest);
P = num2cell(ptest);

[phi0, carcapN, carcapphi, rs, alpha, zr, ds, zd, k, kdrug, gtot] = deal(P{:});

% Reset the "nominal" parameters that we will be fitting for.
%6 fit parameters
phi0= 0.85;
rs = 0.02;
alpha = 0.1;
zr = 0.2;
ds = 0.03;
zd = 0.3;
p = [ phi0, carcapN,rs,alpha, zr, ds, zd];

% Define the inputs needed to generate data
tlong = 0:4:672;
dt = 1;
tdrug = 1;
Cdoxmax = 1000;
Cdox = 200;
tgen = 0:1:tlong(end);
U1=k*Cdox*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
N0 = 3000;
% Generate data from the nominal parameters
[Nsri, tcrit, Ncrit] = fwd_Greene_model2(p, tlong, N0, U1, dt, tdrug);
pnom = [phi0, rs, alpha, zr, ds, zd];
figure;
plot(tlong, Nsri(:,1), 'b-', 'LineWidth', 2)
hold on
plot(tlong, Nsri(:,2), 'g-', 'LineWidth',2)
plot(tlong, Nsri(:,3), 'r-', 'LineWidth',2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('Number of cells')
xlim([0 tlong(end)])
legend('N(t)', 'S(t)', 'R(t)', 'Location', 'NorthWest')
legend boxoff
title('Toy model from nominal parameters')
%% Add noise and look at your generated measurement data for fitting
eta = 200;
Nsrdata = [];
Nsrdata(:,2:3) = Nsri(:,2:3) + normrnd(0, eta,[length(tlong) 2]);
Nsrdata(Nsrdata(:,2)<0,2)=0;
Nsrdata(:,1) = sum(Nsrdata(:,2:3),2);

phi_t = Nsrdata(:,2)./Nsrdata(:,1);

figure;
subplot(1,2,1)
plot(tlong, Nsrdata(:,1), 'b.')
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
legend('N(t) simulated data', 'Location', 'NorthWest')
legend boxoff
ylabel('Number of cells')
xlim([0 tlong(end)])
title('N(t) simulated data for fitting')
subplot(1,2,2)
plot(tlong, phi_t, 'g.')
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
legend('\phi(t) simulated data', 'Location', 'NorthWest')
legend boxoff
ylabel('Proportion of sensitive cells (\phi(t))')
xlim([0 tlong(end)])
title('\phi(t) simulated data for fitting')

iphi = [1, 126, 169];
phi_data = phi_t(iphi);
tbot = tlong(iphi);

figure;
subplot(1,2,1)
plot(tlong, Nsrdata(:,1), 'b.')
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
legend('N(t) simulated data', 'Location', 'NorthWest')
legend boxoff
ylabel('Number of cells')
xlim([0 tlong(end)])
title('N(t) simulated data for fitting')
subplot(1,2,2)
plot(tbot, phi_data, 'g*', 'LineWidth', 3)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
legend('\phi(t) simulated data', 'Location', 'NorthWest')
legend boxoff
ylabel('Proportion of sensitive cells (\phi(t))')
xlim([0 tlong(end)])
title('\phi(t) simulated data for fitting')
%% Fit the data
pset = [carcapN, carcapN]; % here phi and N from the same carcap
psetID = [ 3, 4]; % carcapN, carcapphi(=carcapN)
pfitID = [ 1, 2, 5, 6, 7, 8]; % corresponds to phi0, rs, alpha, zr, ds, zd

phi0guess = 0.7;
rsguess = 0.01;
zrguess = 0.3;
alphaguess = 0.05;
dsguess = 0.02;
zdguess = 0.2;
theta = [phi0guess, rsguess, zrguess, alphaguess, dsguess, zdguess];

%phi0, rs, alpha, zr, ds, and zd
pbounds =  [ 0,1; 0,1; 0,1; 0,1; 0,1; 0,1]; 
% For now, set an arbitrary lambda
lambda = 0.9;
sigmafit = eta*ones(length(tlong),1);
phisigfit = eta./Nsrdata(iphi,1);
lengthvec = horzcat(length(tlong), length(U1));
lengthvecphi = horzcat(length(tbot), length(U1));
%% Test out using your fitting function for a single value of lambda
[pbest,N_model, phi_model, negLL, err_N, err_phi] = fit_fxn_Greenephi_Nprof(Nsrdata(:,1),sigmafit,phi_data, phisigfit, pfitID, psetID, theta, pset, tlong,tbot, U1, U1, lengthvec,lengthvecphi, N0,N0,lambda, pbounds);



% Simulate the phi(t) trajectory from all time points
phi_model_long = simmodelgreenephi2(pbest, tlong, N0, pset, U1, [length(tlong) length(tgen)], pfitID, psetID);


figure;
subplot(1,2,1)
plot(tlong, Nsrdata(:,1), 'b.')
hold on
plot(tlong, N_model(:,1), 'b-', 'LineWidth', 2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
legend('N(t) simulated data', 'N(t) fit', 'Location', 'NorthWest')
legend boxoff
ylabel('Number of cells')
xlim([0 tlong(end)])
title('N(t) simulated data for fitting')
subplot(1,2,2)
plot(tbot, phi_data, 'g*', 'LineWidth', 3)
hold on
plot(tlong, phi_model_long, 'g-', 'LineWidth', 2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
legend('\phi(t) simulated data','\phi(t) fit data', 'Location', 'NorthWest')
legend boxoff
ylabel('Proportion of sensitive cells (\phi(t))')
xlim([0 tlong(end)])
title('\phi(t) simulated data for fitting')
%% Loop through lambdas and do this for a bunch and record the parameters
lambdavec = linspace(0.8, 1, 500);
for i = 1:length(lambdavec)
    lambdai = lambdavec(i);
[pbesti,N_model, phi_model, negLL, err_Ni, err_phii] = fit_fxn_Greenephi_Nprof(Nsrdata(:,1),sigmafit,phi_data, phisigfit, pfitID, psetID, theta, pset, tlong,tbot, U1, U1, lengthvec,lengthvecphi, N0,N0,lambdai, pbounds);
pmat(i,:) = pbesti;
err_Nvec(i) = err_Ni;
err_phivec(i) = err_phii;
end
%% Just plot these for now

figure;
plot(err_Nvec, err_phivec, '.')
xlabel('weighted error in N(t)')
ylabel('weighted error in \phi(t)')
title ('Pareto front search using \lambda')
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
paramnames = {'phi0', 'rs', 'zr', 'alpha', 'ds', 'zd'};
xlim([2e2 1e3])
ylim([1e-2 1])

figure;
for j = 1:length(pbest)
    subplot(2, 3, j)
    hist(pmat(:,j))
    hold on
    plot(pnom(j), 1, 'r*', 'LineWidth',5)
    xlabel('Parameter value')
    ylabel('Frequency')
    legend('Parameter distribution', 'Nominal parameter value')
    legend boxoff
    title(paramnames(j))
    set(gca,'FontSize',20,'LineWidth',1.5)
end




