% Try to integrate CLL data into this framework!

%First pass at fitting the N(t) data and the phi(t) estimates (just
% transcribed in from the latest python results from scRNAseq)



% It uses the N0s, time vectors, and error in data for weighting from the
% actual data.
 close all; clear all; clc
 
 %% Load in data from 231 dose response and from fit
S = load('../out/CLLdataresp.mat');
CLL= S.CLLdata;
tdata = CLL(1).time;
Ndata = CLL(1).rawN;

KuntT25 = load('../out/KuntCLLT25.mat');
KuntT25 = struct2cell(KuntT25);
KuntT25 = cell2mat(KuntT25);
carcapphi = KuntT25;
carcapN =2e6; % set this based on final outgrowth 


% We will use this carcapNf only, and the k and kdrug to be consistent
kdrug = 0.0175*0.75; % make the drug last longer
k = 0.5;
% load in the scRNAseq estimates of phi(t) and the corresponding time
% vectors from python

phi_est_filename = '../data/phi_t_CLL.csv';
phi_est = readtable(phi_est_filename);
tbot = phi_est.t; % this is 29 days 
phitrt = phi_est.phi_t;
ntrt = phi_est.ncells;
sigtech = 1e-2;
phisigfit = [phitrt.*(1-phitrt)./ntrt] + sigtech;
N0phi = 1e7; % set this because this is what we think we seeded for this experiment
phi0= phitrt(1);
S0phi = phi0*N0phi;
R0phi = (1-phi0)*N0phi;
Cfluphi = 5; % I have no idea what this trea
Cflumax = 50;
tgen = [0:1:tbot(end)];
tvec = [0:1:tdata(end)];
Uphi=k*Cfluphi*exp(-kdrug*(tgen))/(0.1*Cflumax);
Unt =k*Cfluphi*exp(-kdrug*(tvec))/(0.1*Cflumax);
lengthvecphi = [length(tbot), length(tgen)];
lengthvecN = [length(tdata), length(tvec)];

%% Plot the phi(t) data from single cell sequencing output
figure;
subplot(2,1,1)
errorbar(tbot, phitrt, phisigfit/2,  'go', 'LineWidth', 4)
hold on
errorbar(tbot, 1-phitrt, phisigfit/2, 'ro', 'LineWidth', 4)
legend('\phi_{sens}(t)', '\phi_{res}(t)', 'Location', 'NorthWest')
legend boxoff
xlabel('time(hours)')
ylabel(' \phi(t) data')
xlim([ 0 tbot(end)])
ylim([ 0 1.2])
title('CLL \phi(t) for dosing for scRNAseq expt')
set(gca,'FontSize',20,'LineWidth',1.5)
subplot(2,1,2)
plot(tgen, Uphi, 'b-', 'LineWidth', 3)
xlabel('time (hours)')
ylabel('Effective dose U(t)')
title('CLL: U(t) for dose for scRNAseq expt')
xlim([ 0 tbot(end)])
set(gca,'FontSize',20,'LineWidth',1.5)
%% Plot N(t) curves
dt = 1;
% input time vectors for each different dose response
Cdoxvec = [];
figure;
subplot(2,1,1)
plot(tdata, Ndata, 'k*-','LineWidth', 3 )
legend('TP0')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([ 0 tvec(end)])
xlabel('time(hours)')
ylabel('N(t)')
title('CLL N(t) for longitudinal expt')
subplot(2,1,2)
plot(tvec, Unt, 'b-', 'LineWidth', 3)
xlabel('time (hours)')
ylabel('Effective dose U(t)')
title('CLL: U(t) for dose for longitudinal')
xlim([ 0 tvec(end)])
set(gca,'FontSize',20,'LineWidth',1.5)

%% Compile dosed data
% Here we're going to generate dosed data and output N(t) and phi(t)
% Change pset to only phi0 and carcaps

psetID = [ 3, 4]; % phi0, carcapN, carcapphi
pfitID = [ 1, 2, 5, 6, 7, 8]; % corresponds to phi0, rs, alpha, rr, ds, dr

% Get what we need from real data
% For right now- we're just going to use TP0 and have one Unt treatment,
% but in the future we could see about trying to use the CXCR4+ and CD18+
% IF we make the very big assumption that those represent the sensitive and
% resistant populations, and we set those as an initial condition of all
% one phenotype 
sigmafit = 1e4*ones(length(tdata),1);
ytimefit = tdata;
ydatafit = Ndata;
N0s = Ndata(1);
Uvec = Unt;
%% Set the parameters for setting and fitting
phi0guess = phi0;
gr = CLL(1).params(2);
pset = [carcapN, carcapphi];
rstar = phi0/(1-phi0);
zrguess = 0.4;
rsguess =  1.8*gr;
alphaguess = 0.1;
dsguess = 0.04;
zdguess = 0.1;
theta = [phi0guess, rsguess, zrguess, alphaguess, dsguess, zdguess];

pbounds =  [ 0,1; 0, 1; 0,1; 0,1; 0,1; 0,1]; 
%Give this function both Ntrt and phitrt
% can toggle the amount that we weigh each portion..
% lambda =0-- fit on N(t) only. if lambda =1, fit on phi(t) only
lambda = 0.5;
% This function internally instead of actually fitting rr and dr, fits the ratio 
[pbest,N_model, phi_model, negLL, err_N, err_phi] = fit_fxn_Greenephi_Nprof(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Uphi, lengthvecN,lengthvecphi, N0s,N0phi,lambda, pbounds);
% Adjust transformed parameter estimates so we capture estimated value
% of rr and dr (since these are what we saved).
phi0 = pbest(1);
rs =pbest(2);
alpha = pbest(3);
zr = pbest(4);
rr = rs*zr;
ds = pbest(5);
zd = pbest(6);
dr = ds*zd;
N_model_long = simmodelgreene2(pbest,tvec, N0s, pset, Unt, [length(tvec) length(tvec)], pfitID, psetID);
phi_model_long = simmodelgreenephi2(pbest, tgen, N0phi, pset, Uphi, [length(tgen) length(tgen)], pfitID, psetID);

CCC_vec(1) = f_CCC([N_model, ydatafit], 0.05);
CCC_vec(2) = f_CCC([phi_model, phitrt], 0.05);

psave = [phi0, carcapN, carcapphi, rs, alpha, zr, ds, zd, k, kdrug];
save('../out/pCLL.mat', 'psave')
%%
figure;
subplot(1,2,1)
errorbar(ytimefit, ydatafit, 1.96*sigmafit/2, 'b*')
hold on
plot(tvec, N_model_long, 'k-', 'LineWidth',2)
xlabel ('time (hours)')
ylabel(' N(t)')
legend ( 'N(t) data & 95% CI','model fit', 'Location', 'NorthWest')
legend box off
title (['N(t), CCC_{N}=', num2str(CCC_vec(1))])
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,2,2)
%plot(tbot, phitrt, 'g*', 'LineWidth',5)
hold on
errorbar(tbot, phitrt, phisigfit/2,  'go', 'LineWidth', 4)
plot(tgen, phi_model_long,'k-', 'LineWidth',2)
xlabel ('time (hours)')
ylabel(' \phi_{sens}(t)')
legend ('\phi_{sens}(t) data','model fit', 'Location', 'NorthWest')
legend box off
title (['\phi(t), CCC_{\phi}=', num2str(CCC_vec(2))])
set(gca,'FontSize',20,'LineWidth',1.5)
ylim([0 1.2])
xlim([0 tbot(end)])


