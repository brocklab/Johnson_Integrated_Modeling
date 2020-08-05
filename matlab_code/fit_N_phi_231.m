% This script loads in the combined datasets for each treatment from the
% 231 cells and the phi(t) estimates from the scRNAseq classifier output.
% Goal here is to perform a joint calibration to the data and compare the
% results from including the phi(t) estimates and not including them.

close all; clear all; clc;
%% Load in the combined N(t) data set

S = load('../out/trajsumfit231.mat');
trajsum= S.trajsum;
% Use the fit parameters from N(t) only to make the parameter domains
pfitN = load('../out/pfitN.mat');
pfitN = struct2cell(pfitN);
pfitN = cell2mat(pfitN);
P = num2cell(pfitN);
% Can use some of these as first guesses/ballpark ranges of what to expect

[phi0f, rsf, carcapN, alphaf, rrf, dsf, drf, k, kdrug] = deal(P{:});
% We will use this carcapNf only, and the k and kdrug to be consistent

%% Load in the scRNAseq estimates of phi(t)
% Also add in the corresponding scenario parameters

% Run this line for the PCA estimates
%phi_est_filename = '../data/phi_t_est_pythPCA.csv';
phi_est_filename = '../data/phi_t_est_pythSVM.csv';
% Run this line for the SVM estimates
phi_est = readtable(phi_est_filename);
tphi = phi_est.t; % time points of phi collection
phitrt = phi_est.phi_t;
ntrt = phi_est.ncells;
sigtech = 1e-2;
phisigfit = [phitrt.*(1-phitrt)./ntrt] + sigtech;
N0phi = 0.8*0.24e6; % set this because this is what we think we seeded for this experiment
Cdoxphi = 550;
Cdoxmax = 1000;
tgen = [0:1:tphi(end)];
Uphi=k*Cdoxphi*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
lengthvecphi = [length(tphi), length(tgen)];
carcapphi =20e6; % set this based on final outgrowth 

figure;
errorbar(tphi, phitrt, 1.96.*phisigfit,  'go', 'LineWidth', 4)
hold on
errorbar(tphi, 1-phitrt, 1.96.*phisigfit, 'ro', 'LineWidth', 4)
legend('\phi(t)=\phi_{S}(t)', '1-\phi(t)=\phi_{R}(t)', 'Location', 'NorthWest')
legend boxoff
xlabel('time (hours)')
ylabel(' \phi(t)')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([ 0 tphi(end)+10])
ylim([ 0 1.2])
%title('\phi(t) for dosing for scRNAseq expt')
figure;
plot(tgen, Uphi, 'b-', 'LineWidth', 3)
xlabel('time (hours)')
ylabel('Effective dose u(t)')
title('U(t) for dose for scRNAseq expt')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([ 0 tphi(end)])

%% Plot the mean N(t) data 
  figure;
  dosevec = [ 1, 3, 5];
  dosevecall = 1:7;
for j=1:length(dosevecall)
    i = dosevecall(j);
 
     subplot(2,1,2)
         plot(trajsum(i).tvec, trajsum(i).Nmean, 'color', trajsum(i).color, 'LineWidth', 2)
         hold on
         %text(trajsum(i).tvec(end-4), trajsum(i).Nmean(end)+1.96*trajsum(i).Nstd(end), [num2str(trajsum(i).Cdox),' nM'], 'FontSize', 14)
         %for m = 1:6
             %plot(trajsum(i).tvec, trajsum(i).Nfit(:,m), 'color', trajsum(i).color, 'LineWidth',1)
         %end
          errorbar(trajsum(i).tvec, trajsum(i).Nmean, 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
         %plot(trajsum(i).tvec, trajsum(i).Nmean - 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
         
         xlabel('time (hours)')
        ylabel('N(t) experimental data')
        %title('N(t) data')
        set(gca,'FontSize',20,'LineWidth',1.5)
        xlim([0 trajsum(7).tvec(end)])
        ylim([0 6e4])
        dt = 1;
       subplot(2,1,1)
       ttest = [];
       ttest = 0:dt:trajsum(i).tvec(end);
       plot(ttest, trajsum(i).U,'*', 'color',trajsum(i).color, 'LineWidth',1)
        hold on
        xlabel('time (hours)')
        ylabel('Effective dose (u(t))')
        %title('Effective dose (u(t)) ')
        set(gca,'FontSize',20,'LineWidth',1.5)
        legend('0 nM', '25 nM', '50 nM', '75 nM','100 nM', '150 nM', '200 nM', '300 nM', '500 nM','1000 nM', 'Location', 'NorthEast')
        legend boxoff
        xlim([0 trajsum(7).tvec(end)])
end
%% Make plot of trit versus dox
figure;
for i = 1:length(trajsum)
    errorbar(trajsum(i).Cdox, trajsum(i).tcrit, 1.96*trajsum(i).tcritstd, '*', 'color', trajsum(i).color', 'LineWidth',2)
hold on
 set(gca,'FontSize',20,'LineWidth',1.5)
 %legend('0 nM', '25 nM', '50 nM', '75 nM','100 nM', '150 nM', '200 nM', '300 nM', '500 nM','1000 nM', 'Location', 'NorthEast')
 xlabel('[Dox]')
 ylabel('t_{crit}')
 xlim([0 350])
end
%% Compile dosed data and fit it using puntfit
% Here we're going to generate dosed data and output N(t) and phi(t)
% Change pset to only phi0 and carcap
psetID = [ 3, 4]; % carcapN, carcapphi
pfitID = [ 1, 2, 5, 6, 7, 8]; % corresponds to phi0, rs, alpha, rr, ds, dr

% Get what we need from real data
sigmafit = [];
tN = [];
ydatafit = [];
N0s = [];
lengtht = [];
lengthU = [];
Uvec = [];

for j=1:length(dosevec)
    i = dosevec(j);
sigmafit = vertcat(sigmafit,trajsum(i).Nstd(1:end));
tN = vertcat(tN, trajsum(i).tvec(1:end));
ydatafit = vertcat(ydatafit, trajsum(i).Nmean(1:end));
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
%% Now fit your data using both Ntrt and phitrt

pset = [ carcapN, carcapphi];
% Set your guess to the best fitting parameter values from N

phi0guess = 0.8;
alphaguess =  0.01;
%rrguess = 1e-3;
zrguess = 0.5;
rsguess = 0.7.*0.03; % expression that relates gs, gr, and gtot
dsguess = 0.2;
%drguess = 0.1*dsguess;
zdguess = 0.5;
theta = [phi0guess, rsguess, alphaguess, zrguess, dsguess, zdguess];
pbounds =  [ 0,1; 0,1; 0,1; 0,1; 0,1; 0,1]; 

% objfun loads in all of the necessary scenario parameters and then outputs
% an error which is stacked with N(t) error on top and phi(t) error below
fhNphi = @(p)objfun(p, ydatafit, phitrt, tN, tphi,N0s, N0phi, Uvec, Uphi, lengthvec, lengthvecphi, pfitID, psetID, pset, sigmafit, phisigfit);
fhN = @(p)objfunN(p, ydatafit, tN,N0s, Uvec, lengthvec, pfitID, psetID, pset, sigmafit);


%% Try using lsqnonlin instead
[pbest] = lsqnonlin(fhNphi, theta, pbounds(:,1), pbounds(:,2));
[pbestN]=lsqnonlin(fhN, theta, pbounds(:,1), pbounds(:,2));

modelfunN = @(p)simmodelgreene2(p, tN, N0s, pset, Uvec, lengthvec, pfitID, psetID); 
modelfunphi = @ (p)simmodelgreenephi2(p, tphi, N0phi, pset, Uphi, lengthvecphi, pfitID, psetID);

model_N = modelfunN(pbest);
model_phi = modelfunphi(pbest);
model_N_N = modelfunN(pbestN);
model_phi_N = modelfunphi(pbestN);
wres= fhNphi(pbest);
wresN = fhN(pbestN);
JNphi = sum(wres);
JN = sum(wresN);
figure;
plot(1:1:length(wres), wres)

%% Original objective function
%Give this function both Ntrt and phitrt
%[pbest,model_N, model_phi, negLL, pbestN, model_N_N, model_phi_N] = fit_fxn_phi_N(ydatafit, sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, tN, tphi, Uvec, Uphi, lengthvec,lengthvecphi, N0s, N0phi, pbounds);

%% Plot model fit versus data 
phi0 = pbest(1);
rs =pbest(2);
alpha = pbest(3);
zr = pbest(4);
rr = zr*rs;
ds = pbest(5);
zd = pbest(6);
dr = zd*ds;
for i = 1:length(trajsum)
 N_modeli = [];  
N_modelit = simmodelgreene2(pbest, trajsum(i).tvec, trajsum(i).Nmean(1), pset, trajsum(i).U, [length(trajsum(i).Nmean) length(trajsum(i).U)], pfitID, psetID);
trajsum(i).Nmodel = N_modelit;
N_modelN = simmodelgreene2(pbestN, trajsum(i).tvec, trajsum(i).Nmean(1), pset, trajsum(i).U, [length(trajsum(i).Nmean) length(trajsum(i).U)], pfitID, psetID);
trajsum(i).NmodelN = N_modelN;
end
phi_model_long = simmodelgreenephi2(pbest, tgen, N0phi, pset, Uphi, [length(tgen) length(tgen)], pfitID, psetID);
phi_model_long_N = simmodelgreenephi2(pbestN, tgen, N0phi, pset, Uphi, [length(tgen) length(tgen)], pfitID, psetID);

CCC_Nphi(1) = f_CCC([model_N, ydatafit], 0.05);
CCC_Nphi(2) = f_CCC([model_phi, phitrt], 0.05);
CCC_N(1) = f_CCC([model_N_N, ydatafit], 0.05);
CCC_N(2) = f_CCC([model_phi_N, phitrt], 0.05);


pNphi = [phi0, carcapN, carcapphi, rs, alpha, zr, ds, zd, k, kdrug];
pN = [pbestN(1), carcapN, carcapphi, pbestN(2), pbestN(3), pbestN(4), pbestN(5), pbestN(6), k, kdrug];
save('../out/pNphi.mat', 'pNphi')
save('../out/pN.mat', 'pN')
save('../out/trajsumfit231.mat', 'trajsum');

% Plot this fit, then search within a reasonable lambda and plot those fits.
figure;
hold on
%plot(ytimefit, N_model, 'k*', 'LineWidth',3)
for j = 1:length(dosevec)
    i = dosevec(j);
    errorbar(trajsum(i).tvec, trajsum(i).Nmean, 1.96*trajsum(i).Nstd, '*','color', trajsum(i).color')
    
end
for j = 1:length(dosevec)
    i = dosevec(j);
plot(trajsum(i).tvec, trajsum(i).Nmodel, 'k-', 'LineWidth', 6)
%plot(trajsum(i).tvec, trajsum(i).NmodelN, 'r-', 'LineWidth',6)
    %text(trajsum(i).tvec(end), trajsum(i).Nmean(end), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'], 'FontSize', 14)
end
%plot(ytimefit, N_model, 'k*', 'LineWidth',3)
 %plot(ytimefit, ydatafit-1.96*sigmafit, 'b.')
%text(ytimeunt(20), Nfitunt(20,ind), ['CCC_{Ntrt}=', num2str(CCC_Ntrt(ind)),', CCC_{pfit}=', num2str(CCC_ptrt(ind))])
xlabel ('time (hours)')
ylabel(' N(t)')
ylim([ 0 5e4])
xlim([0 trajsum(i).tvec(end)])
legend(' 0 nM', '50 nM', '100 nM', 'model fit', 'Location', 'NorthWest')
%legend ( 'model fit N(t), \lambda*', 'N(t) data', 'Location', 'NorthWest')
legend box off
%title (['N(t), CCC_{N}=', num2str(CCC_vec(1))])
set(gca,'FontSize',26,'LineWidth',1.5)

figure;
%plot(tbot, phitrt, 'g*', 'LineWidth',5)
hold on
plot(tgen, phi_model_long,'k-', 'LineWidth',3)
errorbar(tphi, phitrt, 1.96*phisigfit,  'go', 'LineWidth', 3)

%plot(tbot, phitrt + 1.96*phisigfit, 'k.')
%plot(tbot, phitrt - 1.96*phisigfit, 'k.')
%plot(tbot, modelfunphi(pbestf), 'r', 'LineWidth', 2)
xlabel ('time (hours)')
ylabel(' \phi(t)')
legend ('integrated model fit', '\phi(t) data', 'Location', 'NorthEast')
legend box off
%title (['\phi(t), CCC_{integrated fit}=', num2str(CCC_Nphi(2))])
set(gca,'FontSize',26,'LineWidth',1.5)
ylim([0 1.2])
xlim([0 1656])

figure;
%plot(tbot, phitrt, 'g*', 'LineWidth',5)
hold on
plot(tgen, phi_model_long_N, 'r', 'LineWidth', 3)
errorbar(tphi, phitrt, 1.96*phisigfit,  'go', 'LineWidth', 3)

%plot(tbot, phitrt + 1.96*phisigfit, 'k.')
%plot(tbot, phitrt - 1.96*phisigfit, 'k.')
%plot(tbot, modelfunphi(pbestf), 'r', 'LineWidth', 2)
xlabel ('time (hours)')
ylabel(' \phi(t)')
legend ('N(t) model fit', '\phi(t) data', 'Location', 'NorthEast')
legend box off
%title (['\phi(t), CCC_{N(t) fit}=', num2str(CCC_N(2))])
set(gca,'FontSize',26,'LineWidth',1.5)
ylim([0 1.2])
xlim([0 1656])
%% Bootstrapping for CI
Ndatares = model_N - ydatafit;
phidatares = model_phi-phitrt;
nN = length(Ndatares);
nphi = length(phidatares);
% Run a loop to creat and fit synthetic data
nruns = 100;
rand_index_matrix_N  = randi([1 nN], [nN,nruns]);% store
rand_index_matrix_phi = randi([1 nphi], [nphi, nruns]);
pbigboot = zeros(nruns, length(pbest));
pbigbootN = zeros(nruns, length(pbest));
for i = 1:nruns
    % create synthetic data by:
    % randomly sample with replacement from residual list and add to the
    % model fit from original pbest
    rand_int_N  = rand_index_matrix_N(:,i);
    rand_int_phi = rand_index_matrix_phi(:,i);
    YNsim = model_N + Ndatares(rand_int_N);
    YNsim(YNsim<=0)=0;
    Yphisim = model_phi + phidatares(rand_int_phi);
    sigmaN = sigmafit(rand_int_N);
    sigmaphi =phisigfit(rand_int_phi);
   
    fhNphi = @(p)objfun(p, YNsim, Yphisim, tN, tphi,N0s, N0phi, Uvec, Uphi, lengthvec, lengthvecphi, pfitID, psetID, pset, sigmaN, sigmaphi);
    fhN = @(p)objfunN(p, YNsim, tN,N0s, Uvec, lengthvec, pfitID, psetID, pset, sigmaN);
   
     [pbesti] = lsqnonlin(fhNphi, pbest, pbounds(:,1), pbounds(:,2));
     [pbestNi]=lsqnonlin(fhN, pbestN, pbounds(:,1), pbounds(:,2));
    testmod= modelfunN(pbesti);
    testmodphi = modelfunphi(pbesti);
    
    pbigboot(i,:) = pbesti; % gather full model inputs as vector
    pbigbootN(i,:) = pbestNi;
    % For pbest, generate model  

    
    

end

CI =prctile(pbigboot, [2.5, 97.5], 1)
CIN = prctile(pbigbootN, [2.5, 97.5], 1)
% quick test 
figure;
plot(tN, YNsim, '*')
hold on
plot(tN, testmod, '.')

figure;
plot(tphi, Yphisim, '*')
hold on
plot(tphi, testmodphi, '.')
%% 
paramnames = {'\phi_{0}','r_{s}', '\alpha','r_{r}/r_{s}', 'd_{s}', 'd_{r}/d_{s}'};
figure;
[S,AX,BigAx,H,HAx] = plotmatrix(pbigboot);
title(BigAx,'A Comparison of Bootstrap Parameter Estimates from Integrated Fit')
for i = 1:length(paramnames)
ylabel(AX(i,1),paramnames{i})
xlabel(AX(6,i),paramnames{i})
end

figure;
[S,AX,BigAx,H,HAx] = plotmatrix(pbigbootN);
title(BigAx,'A Comparison of Bootstrap Parameter Estimates from N(t) Fit')
for i = 1:length(paramnames)
ylabel(AX(i,1),paramnames{i})
xlabel(AX(6,i),paramnames{i})
end
