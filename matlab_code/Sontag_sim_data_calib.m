% Calibrate simulated data from a new simple model.

% The goal of this script is to verify the data integration framework- 
% i.e. the objective function that combines lower time resolution estimates
% of population composition (i.e. phi(t) data) with longitudinal bulk 
% population data (i.e. N(t)).

% In this simulated case, we want to generate data from the Sontag model,
% and compare the accuracy of estimating the model parameters from N(t) 
% alone and from N(t) and phi(t)

% We can do this for a number of different "scenarios" i.e data acqusition
% and treatments

% The model:
% dS/dt = gsS(1-(S+R)/K) -alphaS-dsu(t)S
% dR/dt = grR(1-(S+R)/K)+alphaS-dru(t)R

% All parameters: 
% phi0, rs, carcap, alpha, rr, ds, dr

% Scenario 1: N(t) only
% Scenario 2: N(t) and phi(t)

% Let's compare the general accuracy of parameter estimation for both of
% these scenarios 
close all; clear all; clc
%% Set and store parameters
nsamps = 100;

pfit = load('../out/pfitN');
pfit = struct2cell(pfit);
pfit = cell2mat(pfit);
P = num2cell(pfit);
% These are our parameter estimates, which we want to make sure are in our
% ranges
[phi0f, rsf, carcapN, alphaf, rrf, dsf, drf, k, kdrug] = deal(P{:});
pfitf = [phi0f, rsf, alphaf, rrf, dsf, drf];


% Set these arbitrarily 
phidom =linspace(0.6,1,100);
rsdom = linspace(0.15, 0., 100);
carcapdom = linspace(4.8e4, 5.5e4, 100);
alphadom = linspace(0, 0.2, 100);
rrdom = linspace(0, 0.1, 100);
dsdom = linspace(0,0.3, 100);
drdom = linspace(0,0.2, 100);
carcapphi = 20e6;

%Set them around the best fitting parameter values
phidom = linspace(0.5*phi0f, 1, 100);
rsdom = linspace(0.5*rsf, 1.5*rsf, 100);
carcapdom = linspace(4.8e4, 5.5e4, 100);
alphadom = linspace(0.5*alphaf, 1.5*alphaf, 100);
rrdom = linspace(0.5*rrf, 1.5*rrf, 100);
dsdom = linspace(0.5*dsf, 1.5*dsf, 100);
drdom = linspace(0.5*drf, 1.5*drf, 100);

carcapphi = 20e6;

for i = 1:nsamps
    % Draw independent random samples from parameter space 
    phi0=randsample(phidom,1);
    rs = randsample(rsdom,1);
    carcapN = randsample(carcapdom,1);
    alpha = randsample(alphadom,1);
    rr = randsample(rrdom,1);
    % don't let resistant grow faster than sensitive population
    if rs< rr
        rr=rs;
    end
    ds = randsample(dsdom,1);
    dr = randsample(drdom,1);
    % don't let resistant death rate be greater than sensitive death rate
     if dr> ds
        dr=ds;
    end
    pallstore(i,:) = [phi0, rs, carcapN, carcapphi, alpha, rr, ds ,dr];
end

psetID = [ 3, 4]; % corresponds to carcap
pfitID = [1, 2, 5, 6, 7, 8]; % call other parameters
psetstore = pallstore(:,psetID);
pfitstore = pallstore(:,pfitID);
paramnames = {'\phi_{0}', 'r_{S}', '\alpha', 'r_{R}', 'd_{S}', 'd_{R}'};
%% Get what we need from the real data
% In this experiment, we have longitudinal data from treatments at 8 doses
% (including 0)and scRNAseq data for just 1 treatment scenario.

% load in the longitudinal data to get the correct time and effective dose
% vectors
S = load('../out/trajsumfit231.mat');
trajsum= S.trajsum;
pfit = load('../out/pfitN');
pfit = struct2cell(pfit);
pfit = cell2mat(pfit);
P = num2cell(pfit);
% These are our parameter estimates, which we want to make sure are in our
% ranges
[phi0f, rsf, carcapN, alphaf, rrf, dsf, drf, k, kdrug] = deal(P{:});
pfitf = [phi0f, rsf, alphaf, rrf, dsf, drf];

phi_est_filename = '../data/phi_t_est_pyth.csv';
phi_est = readtable(phi_est_filename);
tphiexp = phi_est.t; % time points of phi collection
ntrt = phi_est.ncells;
sigtech = 1e-2;
phitrt = phi_est.phi_t; % don't need
N0phi = 0.8*0.24e6; % set this because this is what we think we seeded for this experiment
Cdoxphi = 550;
Cdoxmax = 1000;
tgen = [0:1:tphiexp(end)];
Uphi=k*Cdoxphi*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
carcapphi =20e6; % set this based on final outgrowth 



% Get what we need from real data
sigmafit = []; % keep this but don't need it. 
tN = [];
ydatafit = [];
N0s = [];
lengtht = [];
lengthU = [];
Uvec = [];
dt = 1;
dosevec = [ 1 3 5];
for j=1:length(dosevec)
    i = dosevec(j);
    % Run to vary t time sampling
 %in = 1:1:length(trajsum(i).Nstd(1:end));
iN = 1:10:length(trajsum(i).Nstd(1:end));
sigfit = trajsum(i).Nstd(iN);
sigmafit = vertcat(sigmafit,sigfit);
tsim = trajsum(i).tvec(iN);
tN = vertcat(tN, tsim);
N0s = vertcat(N0s,trajsum(i).Nmean(1));
lengtht = vertcat(lengtht, length(trajsum(i).Nmean(iN)));
lengthU = vertcat(lengthU, length(trajsum(i).U));


Uvec = vertcat(Uvec, trajsum(i).U');
end
lengthvec = horzcat(lengtht, lengthU);
Cdoxvec = [];
for i = 1:length(trajsum)
        Cdox = trajsum(i).Cdox;
        Cdoxvec = vertcat(Cdoxvec,Cdox);
end
Cdoxdata = Cdoxvec(dosevec);

figure;
for i = 1:length(pfitID)
    subplot(2, 3, i)
    histogram(pfitstore(:,i) ,10)
    hold on
    plot(pfitf(i), 1, '*', 'LineWidth', 10)
    xlabel([paramnames(i)])
    ylabel('frequency')
    if i ==6
    legend('samples', '\theta')
    legend boxoff
    end
    set(gca,'FontSize',20,'LineWidth',1.5)
end

%% Generate simulated data 
% set an initial value for eta
eta = 2.*mean(sigmafit);

% SET the time resolution of the phisimualted data
% start with time sampling every 4 hours
ikeep = [1:400:tphiexp(end)+1];
ikeep = [1, 1201, 1601];
tphi = tgen(ikeep);
lengthvecphi = [length(tphi), length(tgen)];
phisigfit = ones(lengthvecphi(1),1).*sigtech; % don't need





for i = 1:nsamps
% Set up functions to generate the data
piter = pfitstore(i,:);
pset = psetstore(i,:);
modelfunN = @(p)simmodelgreene2(p, tN, N0s, pset, Uvec, lengthvec, pfitID, psetID);
% generates phi data with specified resolution
modelfunphi = @ (p)simmodelgreenephi2(p, tphi, N0phi, pset, Uphi, lengthvecphi, pfitID, psetID);
Nmod = modelfunN(piter);
Nsim = Nmod + 4.*normrnd(0,1).*sigmafit + normrnd(0, eta,[length(tN) 1]); % make data systematically off
%Nsim = Nmod  + normrnd(0, eta,[length(tN) 1]); % make them systematically off
Nsim(Nsim<0)=0;
phimod = modelfunphi(piter);
phisim = phimod + normrnd(0, 1e-2, [length(tphi) 1]);
phisim(phisim<0)=0;


% store the generated data:
Nsimstore(i,:) = Nsim;
phisimstore(i,:) = phisim';

end
%% Fit the simulated data 
phi0guess = 0.9;
alphaguess =  0.1;
rrguess = 1e-3;
rsguess = 0.7.*0.03; % expression that relates gs, gr, and gtot
dsguess = 0.2;
drguess = 0.1*dsguess;
theta = [phi0guess, rsguess, alphaguess, rrguess, dsguess, drguess];
theta = pfitf;
pbounds =  [ 0,1; 0,1; 0,1; 0,1; 0,1; 0,1]; 
theta = pfitf;
%% Using lsqnonlin
% for i = 1:nsamps
%     Nsim = Nsimstore(i,:);
%     phisim = phisimstore(i,:);
%     fhNphi = @(p)objfun(p, Nsim', phisim', tN, tphi,N0s, N0phi, Uvec, Uphi, lengthvec, lengthvecphi, pfitID, psetID, pset);
%     fhN = @(p)objfunN(p, Nsim', tN,N0s, Uvec, lengthvec, pfitID, psetID, pset);
% 
%     [pbest] = lsqnonlin(fhNphi, theta, pbounds(:,1), pbounds(:,2));
%     [pbestN]=lsqnonlin(fhN, theta, pbounds(:,1), pbounds(:,2));
% 
%     modelfunN = @(p)simmodelgreene2(p, tN, N0s, pset, Uvec, lengthvec, pfitID, psetID); 
%     modelfunphi = @ (p)simmodelgreenephi2(p, tphi, N0phi, pset, Uphi, lengthvecphi, pfitID, psetID);
% 
%     model_N = modelfunN(pbest);
%     model_phi = modelfunphi(pbest);
%     model_N_N = modelfunN(pbestN);
%     model_phi_N = modelfunphi(pbestN);
% 
%     peststore(i,:) = pbest;
%     peststoreN(i,:) = pbestN;
% 
%     Neststore(i,:) = model_N;
%     phieststore(i,:) = model_phi;
%     NeststoreN(i,:) = model_N_N;
%     phieststoreN(i,:) = model_phi_N;
% 
% end
%% Using original objective function
ct_int = 0;
for i = 1:nsamps
    Nsim = Nsimstore(i,:);
    phisim = phisimstore(i,:);
    
[pbest,model_N, model_phi, negLL, pbestN, model_N_N, model_phi_N] = fit_fxn_phi_N(Nsim', sigmafit,phisim', phisigfit, pfitID, psetID, theta, pset, tN, tphi, Uvec, Uphi, lengthvec,lengthvecphi, N0s, N0phi, pbounds);

peststore(i,:) = pbest;
peststoreN(i,:) = pbestN;

Neststore(i,:) = model_N;
phieststore(i,:) = model_phi;
NeststoreN(i,:) = model_N_N;
phieststoreN(i,:) = model_phi_N;
pfset = pfitstore(i,:);

    CCC_vec(i,1) = f_CCC([model_N, Nsim'], 0.05);
    CCC_vec(i,2) = f_CCC([model_phi, phisim'], 0.05);
    CCC_vec(i,3)=f_CCC([pfset', pbest'], 0.05);
    
    CCC_vecN(i,1) = f_CCC([model_N_N, Nsim'], 0.05);
    CCC_vecN(i,2) = f_CCC([model_phi_N, phisim'], 0.05);
    CCC_vecN(i,3)=f_CCC([pfset', pbestN'], 0.05);
    
    if CCC_vec(i,3)>CCC_vecN(i,3)
        ct_int = ct_int +1;
    end

end
%% Plot some comparisons
Nsimlong = reshape(Nsimstore,[size(Nsimstore,1)*size(Nsimstore,2), 1]);
Nestlong = reshape(Neststore, [size(Neststore,1)*size(Neststore,2),1]);
NestlongN = reshape(NeststoreN, [size(NeststoreN,1)*size(NeststoreN,2),1]);
phisimlong = reshape(phisimstore, [size(phisimstore,1)*size(phisimstore,2),1]);
phiestlong = reshape(phieststore, [size(phieststore,1)*size(phieststore,2),1]);
phiestlongN = reshape(phieststoreN, [size(phieststoreN,1)*size(phieststoreN,2),1]);
psimlong = reshape(pfitstore, [size(pfitstore,1)*size(pfitstore,2),1]);
pestlong = reshape(peststore, [size(peststore,1)*size(peststore,2),1]);
pestlongN = reshape(peststoreN, [size(peststoreN,1)*size(peststoreN,2),1]);


CCCNphi(1) = f_CCC([Nestlong, Nsimlong], 0.05);
CCCNphi(2) = f_CCC([phiestlong, phisimlong], 0.05);
CCCN(1) = f_CCC([NestlongN, Nsimlong], 0.05);
CCCN(2) = f_CCC([phiestlongN, phisimlong], 0.05);
CCCpar(1) = f_CCC([pestlong, psimlong], 0.05);
CCCpar(2) = f_CCC([pestlongN, psimlong], 0.05);

% Concordance Plots for N(t) & phi(t)
figure;
plot(Nsimlong, Nestlong, 'k.')
hold on
plot(Nsimlong, NestlongN, 'r.')
plot([min(Nestlong), max(Nestlong)],[min(Nestlong), max(Nestlong)], 'b-', 'LineWidth',2)
xlim([min(Nestlong), max(Nestlong)])
ylim ([min(Nestlong), max(Nestlong)])
xlabel('Simulated N(t)')
ylabel('Estimated N(t)')
title(['CCC_{integrated fit}= ',num2str(round(CCCNphi(1),4)), ', CCC_{N(t) only}= ',num2str(round(CCCN(1),4))])

set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(phisimlong, phiestlong, 'k.')
hold on
plot(phisimlong, phiestlongN, 'r.')
plot([min(phiestlong), max(phiestlong)],[min(phiestlong), max(phiestlong)], 'b-', 'LineWidth',2)
xlim( [min(phiestlong), max(phiestlong)])
ylim( [min(phiestlong), max(phiestlong)])
xlabel('Simulated \phi(t)')
ylabel('Estimated \phi(t)')
title(['CCC_{integrated fit}= ',num2str(round(CCCNphi(2),4)), ', CCC_{N(t) only}= ',num2str(round(CCCN(2),4))])
set(gca,'FontSize',20,'LineWidth',1.5)

%% Concordance plots for individual parameters
figure;
for i = 1:length(pfitID)
    CCCipar(i) = f_CCC([peststore(:,i), pfitstore(:,i)], 0.05);
    CCCiparN(i) = f_CCC([peststoreN(:,i), pfitstore(:,i)], 0.05);
    
    
    subplot(2,3,i)
    plot(pfitstore(:,i), peststore(:,i), 'k.')
    hold on
    plot(pfitstore(:,i), peststoreN(:,i), 'r.')
    plot([min(pfitstore(:,i)) max(pfitstore(:,i))],[min(pfitstore(:,i)) max(pfitstore(:,i))], 'b-', 'LineWidth', 2)
    xlabel(['set ', paramnames(i)])
    ylabel(['estimated ', paramnames(i)])
    xlim([min(pfitstore(:,i)) max(pfitstore(:,i))])
    ylim( [min(pfitstore(:,i)) max(pfitstore(:,i))])
    legend('integrated fit', 'N(t) only', 'Location', 'NorthWest')
    legend boxoff
    title(['CCC in' paramnames(i)])
    set(gca,'FontSize',20,'LineWidth',1.5)
end

figure;
for i = 1:length(pfitID)
    subplot(2,3,i)
    bar([CCCipar(i) CCCiparN(i)])
    xlabel([paramnames(i)])
    ylabel('CCC')
    set(gca,'XTickLabel',{'integrated fit', 'N(t) fit'});
    set(gca,'FontSize',20,'LineWidth',1.5)
    
end
%% An example 
figure;
id=22;
subplot(1,2,1)
plot(tN, Nsimstore(id,:), '*')
hold on
plot(tN,Neststore(id,:), 'k.')
plot(tN, NeststoreN(id,:), 'r.')
legend('simulated data', 'integrated fit', 'N(t) only fit', 'Location', 'NorthWest')
legend boxoff
xlabel('time(hours)')
ylabel('N(t)')
set(gca,'FontSize',20,'LineWidth',1.5)
subplot(1,2,2)
plot(tphi, phisimstore(id,:), 'g*', 'LineWidth',5)
hold on
plot(tphi, phieststore(id,:), 'k*')
plot(tphi, phieststoreN(id,:), 'r*')
ylim([0 1.2])
legend('simulated data', 'integrated fit', 'N(t) only fit', 'Location', 'NorthWest')
legend boxoff
xlabel('time(hours)')
ylabel('\phi(t)')
set(gca,'FontSize',20,'LineWidth',1.5)

%% Next up- which fit is better able to predict new treatments?
% 1. Generate new treatment data
% 2. Predict resposne to new treatment using params from integrated fit and
% params from N(t) only fit and compare results. 

% Get what we need from real data
sigmafit = []; % keep this but don't need it. 
tN = [];
ydatafit = [];
N0s = [];
lengtht = [];
lengthU = [];
Uvec = [];
dt = 1;
dosevec = [2 4 6 7];
for j=1:length(dosevec)
    i = dosevec(j);
sigmafit = vertcat(sigmafit,trajsum(i).Nstd(1:end));
tN = vertcat(tN, trajsum(i).tvec(1:end));
N0s = vertcat(N0s,trajsum(i).Nmean(1));
lengtht = vertcat(lengtht, length(trajsum(i).Nmean));
lengthU = vertcat(lengthU, length(trajsum(i).U));
Uvec = vertcat(Uvec, trajsum(i).U');
end
lengthvec = horzcat(lengtht, lengthU);
Cdoxvec = [];
for i = 1:length(trajsum)
        Cdox = trajsum(i).Cdox;
        Cdoxvec = vertcat(Cdoxvec,Cdox);
end
Cdoxdata = Cdoxvec(dosevec);



% Generate simulated data 
% set an initial value for eta
eta = 0.25*mean(sigmafit);


%where p = [phi0, rs, carcapN, carcapphi, alpha, rr, ds, dr];
for i = 1:nsamps
% Set up functions to generate the data
piter = pfitstore(i,:);
pset = psetstore(i,:);
modelfunN = @(p)simmodelgreene2(p, tN, N0s, pset, Uvec, lengthvec, pfitID, psetID); 
modelfunphi = @ (p)simmodelgreenephi2(p, tphi, N0phi, pset, Uphi, lengthvecphi, pfitID, psetID);
Nmod = modelfunN(piter);
Nsim = Nmod + normrnd(0, eta,[length(tN) 1]);
Nsim(Nsim<0)=0;
phimod = modelfunphi(piter);
phisim = phimod + normrnd(0, 1e-3, [length(tphi) 1]);
phisim(phisim<0)=0;

% store the generated data:
Nsimstorepred(i,:) = Nsim;
phisimstorepred(i,:) = phisim';

end
%% Compare data to predicted treatments

for i = 1:nsamps
    % pest predictions
    
    pest = peststore(i,:);
    pestN = peststoreN(i,:);
    pset = psetstore(i,:);
    modelfunN = @(p)simmodelgreene2(p, tN, N0s, pset, Uvec, lengthvec, pfitID, psetID); 
    modelfunphi = @ (p)simmodelgreenephi2(p, tphi, N0phi, pset, Uphi, lengthvecphi, pfitID, psetID);
    
    Npred = modelfunN(pest);
    phipred = modelfunphi(pest);
    NpredN = modelfunN(pestN);
    phipredN = modelfunphi(pestN);
    
    % store the predicted data
    Npredstore(i,:) = Npred;
    phipredstore(i,:) = phipred';
    NpredstoreN(i,:) = NpredN;
    phipredstoreN(i,:) = phipredN';
end
%% An example prediction
figure;
id=13;
subplot(1,2,1)
plot(tN, Nsimstorepred(id,:), '*')
hold on
plot(tN,Npredstore(id,:), 'k.')
plot(tN, NpredstoreN(id,:), 'r.')
legend('simulated data', 'integrated fit prediction', 'N(t) only fit prediction', 'Location', 'NorthWest')
legend boxoff
xlabel('time(hours)')
ylabel('N(t)')
set(gca,'FontSize',20,'LineWidth',1.5)
subplot(1,2,2)
plot(tphi, phisimstore(id,:), 'g*', 'LineWidth',5)
hold on
plot(tphi, phieststore(id,:), 'k*')
plot(tphi, phieststoreN(id,:), 'r*')
ylim([0 1.2])
legend('simulated data', 'integrated fit', 'N(t) only fit', 'Location', 'NorthWest')
legend boxoff
xlabel('time(hours)')
ylabel('\phi(t)')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Plot some comparisons
Nsimpredlong = reshape(Nsimstorepred,[size(Nsimstorepred,1)*size(Nsimstorepred,2), 1]);
Npredlong = reshape(Npredstore, [size(Npredstore,1)*size(Npredstore,2),1]);
NpredlongN = reshape(NpredstoreN, [size(NpredstoreN,1)*size(NpredstoreN,2),1]);
phisimpredlong = reshape(phisimstorepred, [size(phisimstorepred,1)*size(phisimstorepred,2),1]);
phipredlong = reshape(phipredstore, [size(phipredstore,1)*size(phipredstore,2),1]);
phipredlongN = reshape(phipredstoreN, [size(phipredstoreN,1)*size(phipredstoreN,2),1]);



CCCNphipred(1) = f_CCC([Npredlong, Nsimpredlong], 0.05);
CCCNphipred(2) = f_CCC([phipredlong, phisimpredlong], 0.05);
CCCNpred(1) = f_CCC([NpredlongN, Nsimpredlong], 0.05);
CCCNpred(2) = f_CCC([phipredlongN, phisimpredlong], 0.05);


% Concordance Plots for N(t) & phi(t)
figure;
plot(Nsimpredlong, Npredlong, 'k.')
hold on
plot(Nsimpredlong, NpredlongN, 'r.')
plot([min(Npredlong), max(Npredlong)],[min(Npredlong), max(Npredlong)], 'b-', 'LineWidth',2)
xlim( [min(Npredlong), max(Npredlong)])
ylim( [min(Npredlong), max(Npredlong)])
xlabel('Simulated N(t)')
ylabel('Predicted N(t)')
title(['CCC_{pred integrated fit}= ',num2str(round(CCCNphipred(1),4)), ', CCC_{pred N(t) only}= ',num2str(round(CCCNpred(1),4))])

set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(phisimpredlong, phipredlong, 'k.')
hold on
plot(phisimpredlong, phipredlongN, 'r.')
plot([min(phisimpredlong), max(phisimpredlong)],[min(phisimpredlong), max(phisimpredlong)], 'b-', 'LineWidth',2)
xlim([min(phisimpredlong), max(phisimpredlong)])
ylim([min(phisimpredlong), max(phisimpredlong)])
xlabel('Simulated \phi(t)')
ylabel('Predicted \phi(t)')
title(['CCC_{pred integrated fit}= ',num2str(round(CCCNphipred(2),4)), ', CCC_{pred N(t) only}= ',num2str(round(CCCNpred(2),4))])
set(gca,'FontSize',20,'LineWidth',1.5)
