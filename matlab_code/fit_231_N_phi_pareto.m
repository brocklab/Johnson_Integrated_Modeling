% First pass at fitting the N(t) data and the phi(t) estimates (just
% transcribed in from the latest python results from scRNAseq)



% It uses the N0s, time vectors, and error in data for weighting from the
% actual data.
 close all; clear all; clc
 
 %% Load in data from 231 dose response and from fit
S = load('../out/trajfit231.mat');
traj= S.traj;
% Only look at the first dose data from one experiment
ntot = length(traj);
%
for i = 1:ntot
matches= ismember(traj(i).celltype, '231');
matchesdate = ismember(traj(i).date, '5-6-19');
ind(i) = all(matches, 'all');
ind2(i) = all(matchesdate, 'all');

end
indall = and(ind,ind2);
trajM = traj(indall);
traj = trajM;


% Use the fit parameters from N(t) only to make the parameter domains
p4fit = load('../out/p4fit231.mat');
p4fit = struct2cell(p4fit);
p4fit = cell2mat(p4fit);
P = num2cell(p4fit);
% Can use some of these as first guesses/ballpark ranges of what to expect
[rsf, carcapNf, alphaf, rrf, dsf, drf, k, kdrugi, gtot] = deal(P{:});
% We will use this carcapNf only, and the k and kdrug to be consistent
kdrug = kdrugi*0.75; % make the drug last longer
% load in the scRNAseq estimates of phi(t) and the corresponding time
% vectors from python

phi_est_filename = '../data/phi_t_est_pyth.csv';
phi_est = readtable(phi_est_filename);
tbot = phi_est.t;
phitrt = phi_est.phi_t;
ntrt = phi_est.ncells;
sigtech = 1e-2;
phisigfit = [phitrt.*(1-phitrt)./ntrt] + sigtech;
N0phi = 0.8*0.24e6; % set this because this is what we think we seeded for this experiment
phi0= phitrt(1);
S0phi = phi0*N0phi;
R0phi = (1-phi0)*N0phi;
Cdoxphi = 550;
Cdoxmax = 1000;
tgen = [0:1:tbot(end)];
Ub=k*Cdoxphi*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
lengthvecphi = [length(tbot), length(tgen)];
carcapphi =20e6; % set this based on final outgrowth 

%% Plot the phi(t) data from single cell sequencing output
figure;
errorbar(tbot, phitrt, 5*phisigfit,  'go', 'LineWidth', 4)
hold on
errorbar(tbot, 1-phitrt, 5*phisigfit, 'ro', 'LineWidth', 4)
legend('\phi_{sens}(t)', '\phi_{res}(t)', 'Location', 'NorthWest')
legend boxoff
xlabel('time (hours)')
ylabel(' \phi(t) data')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([ 0 1666])
ylim([ 0 1.2])
%title('\phi(t) for dosing for scRNAseq expt')
figure;
plot(tgen, Ub, 'b-', 'LineWidth', 3)
xlabel('time (hours)')
ylabel('Effective dose U(t)')
title('U(t) for dose for scRNAseq expt')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([ 0 1656])


%% Make cleaned N(t) data from range of dox concentrations
uniqdose = [];
doselist = [];
for i = 1:length(traj)
    if ~isempty(traj(i).dosenum)
    if traj(i).dosenum==1 || traj(i).dosenum == 0
    doselist =vertcat(doselist, traj(i).dose);
    end
    end
end
uniqdose= unique(doselist);
% Make a new structure which combines each dox concentration
 for i = 1:length(uniqdose)
    trajsum(i).Cdox = [];
    trajsum(i).Nmat = [];
    trajsum(i).nreps = 0;
    trajsum(i).tmat = [];
    trajsum(i).tcritvec= [];
 end

 for i = 1:length(uniqdose) % number of unique seed numbers
    for j = 1:length(traj)
        date = {'5-6-19'}; % pull from the same experiment: first treat
        if contains(traj(j).date, date)  % only want data from this run
                if traj(j).dose == uniqdose(i)
                    trajsum(i).nreps = trajsum(i).nreps +1;
                    trajsum(i).Cdox = traj(j).dose;
                    trajsum(i).color = traj(j).color;
                    trajsum(i).tmat = horzcat(trajsum(i).tmat,traj(j).time);
                    trajsum(i).Nmat = horzcat(trajsum(i).Nmat, traj(j).rawN);
                    trajsum(i).tdose = traj(j).tdose;
                    % Need to think on how to deal with sensoring!
                    trajsum(i).tcritvec = horzcat(trajsum(i).tcritvec, traj(j).tcrit);
                end
        end
    end
 end
 %% Again, clean data for fitting... (cut off data before carcap)
 for i = 1:length(trajsum)
%     if i ==1
%     Nfin = 5e4;
%     else
    Nfin = 3.5e4;
    %end
     N = trajsum(i).Nmat;
    t = trajsum(i).tmat;
    i0 = find(t(:,1)>trajsum(i).tdose,1,'first'); % I arbitrarily search for a maximum in the first 200 hours
    iend = find(N(:,1)>=Nfin,1, 'first');
    if ~isempty(iend)
    tfit = t(i0:iend,:)-t(i0, :); 
    Nfit = N(i0:iend, :);
    end
    if isempty(iend)
        tfit = t(i0:end, :)-t(i0, :); 
        Nfit = N(i0:end, :);
    end
    
    trajsum(i).tfit =round(tfit,0);
    trajsum(i).Nfit =Nfit;
    trajsum(i).tvec = trajsum(i).tfit(:,1);
end
% Test and set U(t) curves
dt = 1;
% input time vectors for each different dose response
Cdoxvec = [];
figure;
for i = 1:length(trajsum)
    ttest = [];
    ttest = 0:dt:trajsum(i).tvec(end);
    Cdox = trajsum(i).Cdox;
    Cdoxvec = vertcat(Cdoxvec,Cdox);
    trajsum(i).U = k*Cdox*exp(-kdrug*(ttest))/(0.1*Cdoxmax);  
    subplot(1, length(trajsum),i)
    plot(ttest, trajsum(i).U, 'b-', 'LineWidth',2)
    ylim([0 5])
    xlim([ 0 ttest(end)])
    xlabel('time (hours)')
    ylabel('Effective dose (U(t))')
    title([num2str(trajsum(i).Cdox),' nM Dox'])
end
% Now find mean and standard deviation vectors
for i = 1:length(trajsum)
    trajsum(i).Nmean = mean(trajsum(i).Nfit,2);
    trajsum(i).tvec = round(trajsum(i).tfit(:,1),0);
    trajsum(i).Nstd = std(trajsum(i).Nfit,0,2);
    trajsum(i).tcritmean = mean(trajsum(i).tcritvec);
    trajsum(i).tcritstd = std(trajsum(i).tcritvec);
end
% Plot the average data 
  figure;
  dosevec = [ 1, 4, 7, 9];
  dosevecall = 1:10;
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
        xlim([0 trajsum(10).tvec(end)])
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
        xlim([0 trajsum(10).tvec(end)])
save('../out/trajsumfit231.mat', 'trajsum')
end

%% Make plot of trit verus dox
figure;
for i = 1:length(trajsum)
    errorbar(trajsum(i).Cdox, trajsum(i).tcritmean, 1.96*trajsum(i).tcritstd, '*', 'color', trajsum(i).color', 'LineWidth',2)
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
psetID = [ 3, 4]; % phi0, carcapN, carcapphi
pfitID = [ 1, 2, 5, 6, 7, 8]; % corresponds to phi0, rs, alpha, rr, ds, dr

% Get what we need from real data
sigmafit = [];
ytimefit = [];
ydatafit = [];
N0s = [];
lengtht = [];
lengthU = [];
Uvec = [];
dosevec = [ 1, 4, 7, 9];
for j=1:length(dosevec)
    i = dosevec(j);
sigmafit = vertcat(sigmafit,trajsum(i).Nstd(1:end));
ytimefit = vertcat(ytimefit, trajsum(i).tvec(1:end));
ydatafit = vertcat(ydatafit, trajsum(i).Nmean(1:end));
N0s = vertcat(N0s,trajsum(i).Nmean(1));
lengtht = vertcat(lengtht, length(trajsum(i).Nmean));
lengthU = vertcat(lengthU, length(trajsum(i).U));
Uvec = vertcat(Uvec, trajsum(i).U');
end
lengthvec = horzcat(lengtht, lengthU);
Cdoxdata = Cdoxvec(dosevec);
%% Now fit your data using both Ntrt and phitrt


rr_to_pop_ratio= 0.03/0.0392;
zrdata = rr_to_pop_ratio;
zddata = 0.1;
%pset = [phi0, carcapNf, carcapphi, zrdata];
pset = [ carcapNf, carcapphi];
rstar = phi0/(1-phi0);
zrguess = 0.4;
rsguess =  0.026;
alphaguess = 0.1;
dsguess = 0.04;
zdguess = 0.1;
phi0guess = phi0;
theta = [phi0guess, rsguess, zrguess, alphaguess, dsguess, zdguess];

%phi0, rs, alpha, zr, ds, and zd
% Unfortunately it appears that the 
%pbounds =  [ 0.5*gtot, 2; 0,1; 0,2; 0.0,0.1; 0,1]; 
pbounds =  [ 0,1; 0,1; 0,1; 0,1; 0,1; 0,1]; 
%Give this function both Ntrt and phitrt
% can toggle the amount that we weigh each portion..
% lambda =0-- fit on N(t) only. if lambda =1, fit on phi(t) only
lambda = 0.5;
%% RUN TO FIT ON N(t) only
lambda = 0;

%% Try multistart function with a range of initial guesses 
% Set your range of start values
thetarange = [0.7, 1; 0.02, 0.03; 0,1; 0.01 0.5; 0 0.1; 0,1];
[pbestvec,N_model, phi_model, negLLvec, err_N, err_phi] = fit_fxn_Greenephi_Nms(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta,thetarange, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);
% Want to output the different parameter sets from the multistart
% optimization 
[negLL, ind] = min(negLLvec);
pbest = pbestvec(ind, :);

% Plot parameter distributions from multistart
paramnames = {'\phi_{0}','r_{S}', 'r_{R}/r_{S} ratio', '\alpha', 'd_{S}', 'd_{R}/d_{S} ratio'};
figure;
for i = 1:length(pbest)
    subplot(1,length(pbest), i)
    plot(pbestvec(:,i), negLLvec,'o', 'LineWidth',2)
    hold on
    %xlim([min(lambdavec), max(lambdavec)])
    xlabel('parameter value')
    ylabel('J(\theta)')
    title(paramnames(i));
    set(gca,'FontSize',20,'LineWidth',1.5)
    ylim([ 0 100])
end


%% Old single start function- DON'T RUN
% This function internally instead of actually fitting rr and dr, fits the ratio 
[pbest,N_model, phi_model, negLL, err_N, err_phi] = fit_fxn_Greenephi_Nprof(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);
% Adjust transformed parameter estimates so we capture estimated value
% of rr and dr (since these are what we saved).
%% 
phi0 = pbest(1);
rs =pbest(2);
alpha = pbest(3);
zr = pbest(4);
rr = rs*zr;
ds = pbest(5);
zd = pbest(6);
dr = ds*zd;
for i = 1:length(trajsum)
 N_modeli = [];  
N_modelit = simmodelgreene2(pbest, trajsum(i).tvec, trajsum(i).Nmean(1), pset, trajsum(i).U, [length(trajsum(i).Nmean) length(trajsum(i).U)], pfitID, psetID);
trajsum(i).Nmodel = N_modelit;
end
phi_model_long = simmodelgreenephi2(pbest, tgen, N0phi, pset, Ub, [length(tgen) length(tgen)], pfitID, psetID);

CCC_vec(1) = f_CCC([N_model, ydatafit], 0.05);
CCC_vec(2) = f_CCC([phi_model, phitrt], 0.05);
fvalbest = negLL


psave = [phi0, carcapNf, carcapphi, rs, alpha, zr, ds, zd, k, kdrug, gtot];
save('../out/ptest.mat', 'psave')
%% Plot this fit, then search within a reasonable lambda and plot those fits.
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
    %text(trajsum(i).tvec(end), trajsum(i).Nmean(end), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'], 'FontSize', 14)
end
%plot(ytimefit, N_model, 'k*', 'LineWidth',3)
 %plot(ytimefit, ydatafit-1.96*sigmafit, 'b.')
%text(ytimeunt(20), Nfitunt(20,ind), ['CCC_{Ntrt}=', num2str(CCC_Ntrt(ind)),', CCC_{pfit}=', num2str(CCC_ptrt(ind))])
xlabel ('time (hours)')
ylabel(' N(t)')
ylim([ 0 5e4])
xlim([0 trajsum(i).tvec(end)])
legend(' 0 nM', '75 nM', '200 nM', '500 nM', 'model fit(\lambda*)', 'Location', 'NorthEast')
%legend ( 'model fit N(t), \lambda*', 'N(t) data', 'Location', 'NorthWest')
legend box off
%title (['N(t), CCC_{N}=', num2str(CCC_vec(1))])
set(gca,'FontSize',26,'LineWidth',1.5)

figure;
%plot(tbot, phitrt, 'g*', 'LineWidth',5)
hold on
plot(tgen, phi_model_long,'k-', 'LineWidth',3)
errorbar(tbot, phitrt, 1.96*phisigfit,  'go', 'LineWidth', 3)

%plot(tbot, phitrt + 1.96*phisigfit, 'k.')
%plot(tbot, phitrt - 1.96*phisigfit, 'k.')
%plot(tbot, modelfunphi(pbestf), 'r', 'LineWidth', 2)
xlabel ('time (hours)')
ylabel(' \phi(t)')
legend ('model fit \phi(t), \lambda*', '\phi(t) data', 'Location', 'NorthEast')
legend box off
%title (['\phi(t), CCC_{\phi}=', num2str(CCC_vec(2))])
set(gca,'FontSize',26,'LineWidth',1.5)
ylim([0 1.2])
xlim([0 1656])
%% Model vs measurement N(t) colored by dose
figure;
hold on
%plot(ytimefit, N_model, 'k*', 'LineWidth',3)
for j = 1:length(dosevec)
    i = dosevec(j);
    plot(trajsum(i).Nmean, trajsum(i).Nmodel, '*','color', trajsum(i).color)

end
plot([0 max(ydatafit)], [ 0 max(ydatafit)], 'k-', 'LineWidth', 2)
legend(' 0 nM', '75 nM', '200 nM', '500 nM', 'Location', 'NorthWest')
legend boxoff
xlim([0 max(ydatafit)])
ylim([0 max(ydatafit)])
set(gca,'FontSize',26,'LineWidth',1.5)
xlabel ('Measured N(t)')
ylabel('Calibrated N(t)')
%title(['CCC_{N(t)} calibrated=', num2str(CCC_vec(1))])
%% Make a plot of the theoretical critical time vs dox concentration
tbig = [0:1:2000];
for i = 1:length(Cdoxvec)
    Ui = k*Cdoxvec(i)*exp(-kdrug*(tbig))/(0.1*Cdoxmax); 
    pi = [ phi0, carcapNf,rs,alpha, zr, ds, zd];
    N0i = trajsum(i).Nmean(1);
    tdrug =1;
    [Nsri, tcrit(i), Ncrit] = fwd_Greene_model2(pi, tbig, N0i, Ui, dt, tdrug);
end
% Dox vs tcrit model and measured
figure;
hold on
dosevecpred = [2 3 5 6 8 10]
tcritcalib = [];
for j = 1:length(dosevec)
    i = dosevec(j);
    errorbar(trajsum(i).Cdox, trajsum(i).tcritmean, 1.96*trajsum(i).tcritstd, '*', 'color', trajsum(i).color', 'LineWidth', 3)
    %text(trajsum(i).Cdox+5, trajsum(i).tcritmean-10, ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'], 'FontSize', 14)
    tcritcalib = vertcat(tcritcalib, trajsum(i).tcritmean);
end
plot(Cdoxvec, tcrit, 'k-', 'LineWidth', 3)
CCC_calib = f_CCC([tcrit(dosevec)', tcritcalib], 0.05)
legend(' 0 nM', '75 nM', '200 nM', '500 nM', 'model t_{crit}(\lambda*)', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',26,'LineWidth',1.5)
xlabel('[Dox]')
xlim([ 0 550])
ylabel('t_{crit}')
%title(['CCC_{tcrit} calibrated=', num2str(CCC_calib)])
%% Critical time predicted versus data
tcritlist = [];
tcritpred = [];
figure;
for j = 1:length(dosevecpred)
    i = dosevecpred(j);
plot(trajsum(i).tcritvec, tcrit(i)*ones(length(trajsum(i).tcritvec),1), '*','color', trajsum(i).color, 'LineWidth', 3)
tcritlist = horzcat(tcritlist, trajsum(i).tcritvec);
tcritpred = vertcat(tcritpred, tcrit(i)*ones(length(trajsum(i).tcritvec),1));
hold on
end
plot([1:1:(max(tcrit))], [1:1:(max(tcrit))], 'k-', 'LineWidth',2)
legend('25 nM', '50 nM', '100 nM', '150 nM', '300 nM', '1000 nM', 'Location', 'NorthWest')
legend boxoff
xlabel('measured t_{crit} (hours)')
ylabel('predicted t_{crit} (hours)')
xlim([ 0 max(tcrit)])
ylim([0 max(tcrit)])
set(gca,'FontSize',26,'LineWidth',1.5)
CCC_pred = f_CCC([tcritpred, tcritlist'], 0.05)
%title(['CCC_{tcrit} predicted=', num2str(CCC_pred)])

%% Model calibrated tcrit vs measured t crit
figure;
for j = 1:length(dosevec)
    i = dosevec(j);
plot(trajsum(i).tcritvec, tcrit(i)*ones(length(trajsum(i).tcritvec),1), '*','color', trajsum(i).color, 'LineWidth', 3)
tcritlist = horzcat(tcritlist, trajsum(i).tcritvec);
tcritpred = vertcat(tcritpred, tcrit(i)*ones(length(trajsum(i).tcritvec),1));
hold on
end
plot([1:1:(max(trajsum(i).tcritvec))], [1:1:(max(trajsum(i).tcritvec))], 'k-', 'LineWidth',2)
xlabel('Measured t_{crit} (hours)')
ylabel('Calibrated t_{crit} (hours)')
xlim([ 0 max(trajsum(i).tcritvec)])
ylim([0 max(trajsum(i).tcritvec)])
set(gca,'FontSize',20,'LineWidth',1.5)
title(['CCC_{tcrit} calibrated=', num2str(CCC_calib)])


%% Now flip through, varying lambda each fit
% Eventually, we want to get more of these
lambdaN = 0;
lambdaphi = 1;
lambdavec = linspace(lambdaN, lambdaphi, 1000);
lambdavec = lambdavec';
[vals, ord] = sort(abs(lambdavec-lambda));
lambdavec =lambdavec(ord); % order the lambdas to search by decreasing distance from the center
%% Use homotypic continuation to find new pareto front parameter sets
%since we already performed randomization of the initial guess over a rnage
%of reasonable thetas, we will start there 
for i = 1:length(lambdavec)
    lambdai = lambdavec(i);
    % Set the initial guess of the new lambda as the best fitting value
    % from previous lambda
    if i ==1
    pguess = pbest;
    else
    pguess = pbesti(i-1,:);
    end
    
    [pbesti(i,:),N_modeli(:,i), phi_modeli(:,i), negLLi(:,i), weightederr_Ni, weightederr_phii] = fit_fxn_Greenephi_Nprof(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, pguess, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambdai, pbounds);
    sumsqerrN(i,1) =weightederr_Ni;
    sumsqerrphi(i,1) =weightederr_phii;
    % Also calculate CCC in N and CCC in phi
    CCC_phi(i,1) = f_CCC([phi_modeli(:,i), phitrt], 0.05);
    CCC_N(i,1) = f_CCC([N_modeli(:,i), ydatafit], 0.05);
end

%% Save everything into one table
pareto = table(sumsqerrphi, sumsqerrN, pbesti, CCC_phi, CCC_N, lambdavec);

%% Plot error in N vs. error in phi from lambda search
figure;
plot(err_phi, err_N, 'r*', 'LineWidth', 4)
hold on
scatter(sumsqerrphi, sumsqerrN,[],lambdavec, 'filled')
plot(err_phi, err_N, 'r*', 'LineWidth', 8)
hcb=colorbar;
title(hcb, '\lambda')
xlabel('E_{\phi} (weighted sum-of-squares error in \phi(t))')
ylabel('E_{N} (weighted sum-of-squares error in N(t))')
set(gca,'FontSize',26,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
%xlim([ .8 1e2])
%ylim([ 1e4 2e4])
legend('\theta*', 'Location', 'NorthWest')
legend boxoff


%% Filter by CCC
iphi = CCC_phi>0.8;
iN = CCC_N>0.8;
ikeep = and(iphi, iN);
pareto_filtered = pareto(ikeep,:);

figure;
plot(err_phi, err_N, 'r*', 'LineWidth',8)
hold on
plot( sumsqerrphi(ikeep),sumsqerrN(ikeep), 'c*')
plot( sumsqerrphi(~ikeep),sumsqerrN(~ikeep), 'm*')
plot(err_phi, err_N, 'r*', 'LineWidth',3)
%plot(x,y, 'b-','LineWidth', 3)
ylabel ('E_{N}')
xlabel('E_{\phi}')
%ylim([min(sumsqerrphi), max(sumsqerrphi)])
%xlim([min(sumsqerrN), max(sumsqerrN)])
legend('\theta*','CCC>0.8','CCC<0.8', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
%title('Parameter sets filtered by CCC')
% xlim([ 1 1e2])
% ylim([ 1.3e4 2e4])
%xlim([ .8 1e2])
%ylim([ 1e4 2e4])

%% Try something different, instead, find pareto solutions iteratively
indkeep = [];
for i = 1:height(pareto_filtered)
    iN= [];
    iphi = [];
    % 1. Find the theta that have lower EN than current theta
        iN = pareto_filtered.sumsqerrN <pareto_filtered.sumsqerrN(i);
     % 2. Find the theta that have lower Ephi than current theta
        iphi = pareto_filtered.sumsqerrphi<pareto_filtered.sumsqerrphi(i);
        ibetter = and(iphi, iN);
        if nnz(ibetter) ==0
            indkeep = vertcat(indkeep, i);
        end
end
%% Save the pareto solutions
pareto_in = pareto_filtered(indkeep, :);
writetable(pareto_in,'../out/pbest_table.csv')
%% Plot the parameter sets that we keep
figure;
scatter(pareto_in.sumsqerrphi, pareto_in.sumsqerrN,[], pareto_in.lambdavec, 'filled')
hold on
plot(err_phi, err_N, 'r*', 'LineWidth',3)
plot( sumsqerrphi(ikeep),sumsqerrN(ikeep), 'c*')
%plot(pareto_in.sumsqerrphi, pareto_in.sumsqerrN, 'b*')
scatter(pareto_in.sumsqerrphi, pareto_in.sumsqerrN,[], pareto_in.lambdavec, 'filled')
plot(err_phi, err_N, 'r*', 'LineWidth',8)
hcb=colorbar;
title(hcb,'\lambda')
ylabel ('E_{N}')
xlabel('E_{\phi}')
%ylim([ min(sumsqerrphi), max(sumsqerrphi)])
%xlim([min(sumsqerrN), max(sumsqerrN)])
%xlim([ 1 0.5e2])
%ylim([ 1.3e4 2e4])
%xlim([ .8 1e2])
ylim([ 1e4 2e4])

%set(gca,'FontSize',20,'LineWidth',1.5)
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
legend('\theta_{pareto front}', '\theta*', 'non-pareto \theta', 'Location', 'NorthEast')
legend boxoff
%title('Parameter sets within pareto boundary')
%%
figure;
scatter(pareto_in.sumsqerrphi, pareto_in.sumsqerrN,[], pareto_in.lambdavec, 'filled')
hold on
plot(err_phi, err_N, 'r*', 'LineWidth',8)
hcb=colorbar;
title(hcb,'\lambda')
ylabel ('E_{N} (error in N(t))')
xlabel('E_{\phi} (error in \phi(t))')
%ylim([ min(sumsqerrphi), max(sumsqerrphi)])
%xlim([min(sumsqerrN), max(sumsqerrN)])
xlim([ min(pareto_in.sumsqerrphi) max(pareto_in.sumsqerrphi)])
ylim([ min(pareto_in.sumsqerrN) max(pareto_in.sumsqerrN)])
%set(gca,'FontSize',20,'LineWidth',1.5)
set(gca,'FontSize',24,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
legend('\theta_{pareto front}', '\theta*', 'non-pareto \theta', 'Location', 'NorthEast')
legend boxoff
%title('Parameter sets within pareto boundary')

%% Plot parameter distributions
%paramnames = {'\phi_{0}','r_{s}', 'r_{r}/r_{s} ratio', '\alpha', 'd_{s}', 'd_{r}/d_{s} ratio'};
figure;
for i = 1:length(pbest)
    subplot(2,length(pbest), i)
    scatter(pareto_in.lambdavec, pareto_in.pbesti(:, i),[],pareto_in.CCC_N, 'filled')
    hcb=colorbar;
    title(hcb, 'CCC_{N}')
    set(hcb, 'Ticks', [])
    hold on
    %xlim([min(lambdavec), max(lambdavec)])
    xlabel('\lambda')
    %ylabel('parameter value')
    ylabel(paramnames(i));
    set(gca,'FontSize',20,'LineWidth',1.5, 'XTick', [], 'YTick', [])
    subplot(2,length(pbest), i+length(pbest))
    scatter(pareto_in.lambdavec, pareto_in.pbesti(:, i),[],pareto_in.CCC_phi, 'filled')
    hcb = colorbar;
    title(hcb, 'CCC_{\phi}')
    set(hcb, 'Ticks', [])
    hold on
    %xlim([min(lambdavec), max(lambdavec)])
    xlabel('\lambda')
    %ylabel('parameter value')
    ylabel(paramnames(i));
    set(gca,'FontSize',20,'LineWidth',1.5, 'XTick', [], 'Ytick', [])
end

%% Make histograms of parameter distributions
paramlabs = {'\phi_{0}*','r_{s}*', 'r_{r}/r_{s}*', '\alpha*', 'd_{s}*', 'd_{r}/d_{s}*'};
%paramnames = {' \phi_{0}','r_{s}', 'r_{r}/r_{s}', '\alpha', 'd_{s}', 'd_{r}/d_{s}'};
CI = load('../out/CIpbest.mat')
CI = struct2cell(CI);
CI = cell2mat(CI);

psets = pareto_in.pbesti

figure;
for i = 1:length(pbest)
    subplot(2,length(pbest)/2, i)
    plot(pbest(i), 0, 'r*', 'LineWidth', 5)
    hold on
    histogram(psets(:, i), 25, 'Normalization', 'probability')
    hold on
    plot([pbest(i)], [0], 'r*', 'LineWidth', 10)
    %plot([ CI(i,1) CI(i,1)], [0 0.4], 'g-', 'LineWidth', 3)
    %plot([ CI(i,2) CI(i,2)], [0 0.4], 'g-', 'LineWidth', 3)
    %plot(mean(param_table{:,i}),0, 'g*', 'LineWidth', 3)
    %plot(median(param_table{:,i}),0, 'm*', 'LineWidth', 3)
    hold on
    %xlim([0.99*CI(i,1), 1.1*CI(i,2)])
    %ylim([ 0 0.4])
    xlabel(paramnames(i))
    ylabel('probability')
    %title(paramnames(i));
    %legend([paramlabs(i)], 'FontSize', 14)
    %legend boxoff
    set(gca,'FontSize',24,'LineWidth',1.5)
    
end
%% Plot some examples of fits along the pareto front
% find the highest weighting N data set
[minlam, indmin] = min(pareto_in.lambdavec);
[maxlam, indmax] = max(pareto_in.lambdavec);

for i = 1:length(trajsum)
 N_modeliNmin = [];  
N_modeliNmin = simmodelgreene2(pareto_in.pbesti(indmin,:), trajsum(i).tvec, trajsum(i).Nmean(1), pset, trajsum(i).U, [length(trajsum(i).Nmean) length(trajsum(i).U)], pfitID, psetID);
NmodelN = simmodelgreene2(pbest, trajsum(i).tvec, trajsum(i).Nmean(1), pset, trajsum(i).U, [length(trajsum(i).Nmean) length(trajsum(i).U)], pfitID, psetID);
N_modeliNmax = simmodelgreene2(pareto_in.pbesti(indmax,:), trajsum(i).tvec, trajsum(i).Nmean(1), pset, trajsum(i).U, [length(trajsum(i).Nmean) length(trajsum(i).U)], pfitID, psetID);
trajsum(i).NmodelN = NmodelN;
trajsum(i).NmodelNmin = N_modeliNmin;
trajsum(i).NmodelNmax = N_modeliNmax;
end
phi_model_longNmin= simmodelgreenephi2(pareto_in.pbesti(indmin,:), tgen, N0phi, pset, Ub, [length(tgen) length(tgen)], pfitID, psetID);
phi_model_longN= simmodelgreenephi2(pbest, tgen, N0phi, pset, Ub, [length(tgen) length(tgen)], pfitID, psetID);
phi_model_longNmax= simmodelgreenephi2(pareto_in.pbesti(indmax,:), tgen, N0phi, pset, Ub, [length(tgen) length(tgen)], pfitID, psetID);
figure;
hold on
for j = 1:length(dosevec)
    i = dosevec(j);
plot(trajsum(i).tvec, trajsum(i).NmodelN, 'k-', 'LineWidth', 2)
hold on
plot(trajsum(i).tvec, trajsum(i).NmodelNmin, 'b-', 'LineWidth', 2)
plot(trajsum(i).tvec, trajsum(i).NmodelNmax, 'r-', 'LineWidth', 2)
   
end
for j = 1:length(dosevec)
    i = dosevec(j);
    errorbar(trajsum(i).tvec, trajsum(i).Nmean, 1.96*trajsum(i).Nstd, '*','color', trajsum(i).color')
    
end

xlabel ('time (hours)')
ylabel(' N(t)')
ylim([ 0 5e4])
xlim([0 trajsum(i).tvec(end)])
legend('model fit(\lambda=\lambda*=0.5)', 'model fit(\lambda=0.3)','model fit(\lambda=0.9)' , 'Location', 'NorthEast')
%legend ( 'model fit N(t), \lambda*', 'N(t) data', 'Location', 'NorthWest')
legend box off
%title (['N(t), CCC_{N}=', num2str(CCC_vec(1))])
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
%plot(tbot, phitrt, 'g*', 'LineWidth',5)
hold on
plot(tgen, phi_model_longN,'k-', 'LineWidth',3)
plot(tgen, phi_model_longNmin,'b-', 'LineWidth',3)
plot(tgen, phi_model_longNmax,'r-', 'LineWidth',3)
errorbar(tbot, phitrt, 1.96*phisigfit,  'go', 'LineWidth', 3)

%plot(tbot, phitrt + 1.96*phisigfit, 'k.')
%plot(tbot, phitrt - 1.96*phisigfit, 'k.')
%plot(tbot, modelfunphi(pbestf), 'r', 'LineWidth', 2)
xlabel ('time (hours)')
ylabel(' \phi_{sens}(t)')
legend ('\lambda=\lambda*=0.5', '\lambda=0.3','\lambda = 0.9','\phi_{sens}(t) data', 'Location', 'NorthWest')
legend box off
%title (['\phi(t), CCC_{\phi}=', num2str(CCC_vec(2))])
set(gca,'FontSize',20,'LineWidth',1.5)
ylim([0 1.2])
xlim([ 0 1656])



%% REPEAT BUT TRY Only using N(t) data 
lambdatry = 0;
% This function internally instead of actually fitting rr and dr, fits the ratio 

[pbestN,N_modelN, phi_modelN, negLLN, err_NN, err_phiN] = fit_fxn_Greenephi_Nprof(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambdatry, pboundsprof);
% Adjust transformed parameter estimates so we capture estimated value
% of rr and dr (since these are what we saved).

for i = 1:length(trajsum)
 N_modeliN = [];  
N_modeliN = simmodelgreene2(pbestN, trajsum(i).tvec, trajsum(i).Nmean(1), pset, trajsum(i).U, [length(trajsum(i).Nmean) length(trajsum(i).U)], pfitID, psetID);
trajsum(i).NmodelN = N_modeliN;
end
phi_model_longN= simmodelgreenephi2(pbestN, tgen, N0phi, pset, Ub, [length(tgen) length(tgen)], pfitID, psetID);

CCC_vecN(1) = f_CCC([N_modelN, ydatafit], 0.05);
CCC_vecN(2) = f_CCC([phi_modelN, phitrt], 0.05);
fvalbestN = negLLN

%% Plot this fit, then search within a reasonable lambda and plot those fits.
figure;
hold on

for j = 1:length(dosevec)
    i = dosevec(j);
    errorbar(trajsum(i).tvec, trajsum(i).Nmean, 1.96*trajsum(i).Nstd, '*','color', trajsum(i).color')
    
end
for j = 1:length(dosevec)
    i = dosevec(j);
plot(trajsum(i).tvec, trajsum(i).NmodelN, 'k-', 'LineWidth', 6)
    %text(trajsum(i).tvec(end), trajsum(i).Nmean(end), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'], 'FontSize', 14)
end
%plot(ytimefit, N_model, 'k*', 'LineWidth',3)
 %plot(ytimefit, ydatafit-1.96*sigmafit, 'b.')
%text(ytimeunt(20), Nfitunt(20,ind), ['CCC_{Ntrt}=', num2str(CCC_Ntrt(ind)),', CCC_{pfit}=', num2str(CCC_ptrt(ind))])
xlabel ('time (hours)')
ylabel(' N(t)')
ylim([ 0 5e4])
xlim([0 trajsum(i).tvec(end)])
legend(' 0 nM', '75 nM', '200 nM', '500 nM', 'model fit(\lambda=0)', 'Location', 'NorthEast')
%legend ( 'model fit N(t), \lambda*', 'N(t) data', 'Location', 'NorthWest')
legend box off
%title (['N(t), CCC_{N}=', num2str(CCC_vec(1))])
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
%plot(tbot, phitrt, 'g*', 'LineWidth',5)
hold on
plot(tgen, phi_model_longN,'k-', 'LineWidth',3)
errorbar(tbot, phitrt, 1.96*phisigfit,  'go', 'LineWidth', 3)

%plot(tbot, phitrt + 1.96*phisigfit, 'k.')
%plot(tbot, phitrt - 1.96*phisigfit, 'k.')
%plot(tbot, modelfunphi(pbestf), 'r', 'LineWidth', 2)
xlabel ('time (hours)')
ylabel(' \phi_{sens}(t)')
legend ('model fit \phi_{sens}(t), \lambda=0', '\phi_{sens}(t) data', 'Location', 'NorthWest')
legend box off
%title (['\phi(t), CCC_{\phi}=', num2str(CCC_vec(2))])
set(gca,'FontSize',20,'LineWidth',1.5)
ylim([0 1.2])
xlim([ 0 1656])
%% Plot all pareto front solutions fit to the data...

