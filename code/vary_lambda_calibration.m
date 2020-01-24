% This script furthers tests the calibration method by randomly sampling parameter
% space, generating in silico data, and running the calibration using
% phi(t) & N(t), returning the parameter CCCs, the CCCs on the outputs, the
% CCCs on the parameters, and the mean parameter error as we downsample the
% phi(t) time points. 

%THIS SCRIPT CALIBRATES THE DATA USING ONLY 3 TIME POINTS (CONSISTENT WITH
%OUR CURRENT DATA, AND VARIES THE REGULARIZATION TERM, LAMBDA, THAT
%DETERMINES THE RELATIVE WEIGHT IN FITTING OF N(T) VS PHI(T) DATA
 


% It uses the N0s, time vectors, and error in data for weighting from the
% actual data.
 close all; clear all; clc
 
 %% Load in data
S = load('../out/trajfit.mat');
traj= S.traj;
% Separate by dose
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
 end

 for i = 1:length(uniqdose) % number of unique seed numbers
    for j = 1:length(traj)
        date = {'8-16-18'}; % pull from the same experiment: first treat
        if contains(traj(j).date, date)  % only want data from this run
                if traj(j).dose == uniqdose(i)
                    trajsum(i).nreps = trajsum(i).nreps +1;
                    trajsum(i).Cdox = traj(j).dose;
                    trajsum(i).color = traj(j).color;
                    trajsum(i).tmat = horzcat(trajsum(i).tmat,traj(j).time);
                    trajsum(i).Nmat = horzcat(trajsum(i).Nmat, traj(j).rawN);
                    trajsum(i).tdose = traj(j).tdose;
                end
        end
    end
 end
 % Again, clean data for fitting...
 for i = 1:length(trajsum)
    if i ==1
    Nfin = 5.5e4;
    else
    Nfin = 4e4;
    end
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
kdrug = 0.0175;
k = 0.5;
dt = 1;
% input time vectors for each different dose response
figure;
for i = 1:length(trajsum)
    ttest = [];
    ttest = 0:dt:trajsum(i).tvec(end);
    Cdox = trajsum(i).Cdox;
    trajsum(i).U = k*Cdox*exp(-kdrug*(ttest)); 
    subplot(1, length(trajsum),i)
    plot(ttest, trajsum(i).U, 'b-', 'LineWidth',2)
    ylim([0 300])
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
end
% Plot the average data 
  figure;
 for i = 1:6%length(trajsum)
     subplot(2,1,1)
         plot(trajsum(i).tvec, trajsum(i).Nmean, 'color', trajsum(i).color, 'LineWidth', 2)
         hold on
         text(trajsum(i).tvec(end-10), trajsum(i).Nmean(end-10), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'])
         plot(trajsum(i).tvec, trajsum(i).Nmean + trajsum(i).Nstd, 'color', trajsum(i).color)
         plot(trajsum(i).tvec, trajsum(i).Nmean - trajsum(i).Nstd, 'color', trajsum(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title('N(t) for different single pulse treatments')
        dt = 1;
       subplot(2,1,2)
       ttest = [];
       ttest = 0:dt:trajsum(i).tvec(end);
       plot(ttest, trajsum(i).U,'.', 'color',trajsum(i).color, 'LineWidth',1)
        hold on
        xlabel('time (hours)')
        ylabel('Effective dose U(t)')
        title('U(t) for different single pulse treatments')
 end
%% Set & store known parameters
nsamps = 500;

pbounds = [0,1; 0, 1; 0, 1; 0,1]; 
phidom =linspace(0.7,1,100);
rsdom = linspace(0.001, 0.1, 100);
carcapdom = linspace(4.8e4, 5.5e4, 100);
alphadom = linspace(1e-10, 0.01, 100);
rrdom = linspace(1e-10, 0.001, 100);
dsdom = linspace(1e-10,0.01, 100);

for i = 1:nsamps
    phi0=randsample(phidom,1);
    rs = randsample(rsdom,1);
    carcap = randsample(carcapdom,1);
    alpha = randsample(alphadom,1);
    rr = randsample(rrdom,1);
    % don't let resistant grow faster than sensitive population
    if rs< rr
        rr=rs;
    end
    ds = randsample(dsdom,1);
    dr = 0;
    pallstore(i,:) = [phi0, rs,carcap, alpha, rr, ds ,dr];
    rstar = phi0/(1-phi0);
    puntfstore(i,:) = [carcap];
    pfitstore(i,:) = [ rs, alpha, rr, ds];
end
%% Generate untreated control data
% Fit this only to N(t) 
sigmaunt = trajsum(1).Nstd(1:end);
ytimeunt = trajsum(1).tvec;
Uunt = trajsum(1).U;
N0unt = trajsum(1).Nmean(1);
tdrug = 1;
eta=500;

for i = 1:nsamps
    P= num2cell(pallstore(i,:));
    [phi0, rs, carcap, alpha, rr, ds ,dr]= deal(P{:});
    S0 = N0unt*phi0;
    R0 = N0unt*(1-phi0);
    P= num2cell(pallstore(i,:));
    pit = [S0,R0, rs, carcap, alpha, rr, ds, dr];
    % Generate in silico data and store it
    [Nsrunt, ~,~] = fwd_Greene_model(pit, ytimeunt, Uunt, dt, tdrug);
    Nunt = Nsrunt(:,1) + normrnd(0, eta,[length(ytimeunt) 1]);
    phiunt = Nsrunt(:,2)./Nunt;
    Nstoreunt(:,i)= Nunt;
    phistoreunt(:,i) = phiunt;
    
    % Fit to insilico data and store fit parameters
    [punt, Nmodunt] = fit_untreated(Nunt,ytimeunt, sigmaunt);
    % punt is gtot, carcap
    % gtot is not a parameter we set, so dont include this when comparing
    % parameter error but fine to keep track of it.
    puntfit(i,:) = punt;
    Nfitunt(:,i) = Nmodunt;

end   
%%
figure;
ind = 10;
plot(ytimeunt, Nstoreunt(:,ind), '*')
hold on
plot(ytimeunt, Nfitunt(:,ind), '-')
plot(ytimeunt, Nstoreunt(:,ind) + 1.96*sigmaunt, 'k-')
plot(ytimeunt, Nstoreunt(:,ind) - 1.96*sigmaunt, 'k-')
%text(ytimeunt(20), Nfitunt(15,ind), ['CCC_{N}=', num2str(round(CCC_N(ind),3)),', % error_{carcap}=', num2str(round(pct_err_carcap(ind),2)), '%'],'FontSize',14)
xlabel ('time (hours)')
ylabel(' N(t)')
legend ('in silico data', 'model fit', '95% CI on data', 'Location', 'NorthWest')
legend box off
title ('Example fit to untreated control simulated data')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Generate dosed data and fit it using puntfit as you vary lambda
% Here we're going to generate dosed data and output N(t) and phi(t)
psetID = [1, 3, 7];
pfitID = [2, 4, 5, 6];
% Get what we need from real data
sigmafit = [];
ytimefit = [];
N0s = [];
lengtht = [];
lengthU = [];
Uvec = [];
for i = 2:6%length(trajsum)
sigmafit = vertcat(sigmafit,trajsum(i).Nstd(1:end));
ytimefit = vertcat(ytimefit, trajsum(i).tvec(1:end));
N0s = vertcat(N0s,trajsum(i).Nmean(1));
lengtht = vertcat(lengtht, length(trajsum(i).Nmean));
lengthU = vertcat(lengthU, length(trajsum(i).U));
Uvec = vertcat(Uvec, trajsum(i).U');
end
lengthvec = horzcat(lengtht, lengthU);

% Generate the dosed data 
% Simulate the effect of a very strong pulse treatment on phi(t) and save
% that output for the longest time vector
% We are going to pretend we can capture phi(t) for 8 weeks every 4 hours
% to start
Cdox = 50;
tgen = [0:1:1344];
% Vary lambda (linearly or log space, try both)
lambdavec = [ 0.1;0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9];
%lambdavec = logspace(-1.5,0,21)
%lambdavec = linspace(0,1,11)

tres = 672;
ntps = 1;

%% Rewrite this so you generate all the data, and then subsample for fitting:
Ub=k*Cdox*exp(-kdrug*(tgen));
lam = 50;

% First generate data at tbot = tgen
for i = 1:nsamps
    % Generate N(t) data-- keep this constant
        P=num2cell(pallstore(i,:));
        [phi0, rs, carcap, alpha, rr, ds, dr]= deal(P{:});
        pset = [phi0, carcap, dr];
        pfset = [rs, alpha, rr, ds];
        Ntrt = [];
        phitrt = [];
        for j = 2:6
            Nsr = [];
            U = trajsum(j).U;
            tvec = trajsum(j).tvec;
            N0 = trajsum(j).Nmean(1);
            S0 = phi0*N0;
            R0 = (1-phi0)*N0;
            pit = [S0, R0, rs, carcap, alpha, rr, ds, dr]; 
            % Generate data for N(t)
            [Nsr, ~, ~] = fwd_Greene_model(pit, tvec, U, dt, tdrug);
            Nsrdat = Nsr(:,2:3) + normrnd(0, eta,[length(tvec) 2]);
            ind0=Nsrdat<0;
            Nsrdat(ind0)=0;
            Nsr(:,1) = Nsrdat(:,1) + Nsrdat(:,2);
            Ntrt = vertcat(Ntrt, Nsr(:,1));
        end
         % Generate bottlenecked data for phi(t)
            N0phi = 2e3; % set this because this is what we think we seeded
            S0 = phi0*N0phi;
            R0 = (1-phi0)*N0phi;
            pit2 = [S0, R0, rs, carcap, alpha, rr, ds, dr];
            [Nb, ~,~]=fwd_Greene_model(pit2, tgen, Ub, dt, tdrug);
            Nbnoise=Nb; % set it as this then replace it
            Nbnoise(:,2:3) = Nb(:,2:3) + normrnd(0, lam,[length(tgen) 2]);
            % remove negative numbers
            for j =1:length(Nbnoise)
                if Nbnoise(j,2) <0
                    Nbnoise(j,2)=0;
                end
                if Nbnoise(j,3) <0
                    Nbnoise(j,3) =0;
                end
            end
            Nbnoise(:,1) = Nbnoise(:,2) + Nbnoise(:,3);
            phitrt = Nbnoise(:,2)./Nbnoise(:,1);
            phisig =0.1*ones(length(tgen),1);

            Ntrtstore(:,i) = Ntrt; % store N(t) data for that parameter set
            phitrtstore(:,i) = phitrt; % store phi(t) data for that parameter set
end



%% Take stored data and fit with each lambda value and the same time res.
for m = 1:length(lambdavec)
     %m is indicator for the lambda
    % Need to clear the variables that change size every iteration with
    % changing tbot size
    phitrt = [];
    phimodel = [];
    phismooth = [];
    phisigfit = [];
    phifittrt = [];
    phigiven = [];
    phifit = [];
    tbot = [];
    itime = [];
    % Set your itime based on the time resolution of the mth iteration
     tbot = [0:tres:tgen(end)];
    itime = find(ismember(tgen,tbot));
    % use that to grab from phtrtstore
        for i = 1:nsamps 
        Ntrt = Ntrtstore(:,i); 
        phitrt = phitrtstore(itime,i);
        phisigfit = phisig(itime);
    
                % Now fit your in silico data using both Ntrt and phitrt
                % first fit Ntrt
                gtot = puntfit(i,1);
                carcapfit = puntfit(i,2);
                pset = [phi0, carcapfit, dr];
                rrguess = 1e-2*gtot;
                rstar = phi0/(1-phi0);
                rsguess =  ((rstar+1)*gtot - rrguess)./rstar;
                theta = [rsguess, 0.0035, rrguess, 0.001];

                %Give this function both Ntrt and phitrt
                %
                % can toggle the amount that we weigh each portion..
                % lambda =0-- fit on N(t) only. if lambda =1, fit on phi(t) only
                lambda = lambdavec(m);
                lengthvecphi = [length(tbot), length(tgen)];
               [pbestf,N_model, phi_model, negLL] = fit_fxn_Greenephi_N(Ntrt,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);
    
                pfittrt(i,:) = pbestf;
                Nfittrt(:,i) = N_model;
                % this changes each m iteration
                phifittrt(:,i) = phi_model;
                CCC_vec(i,1) = f_CCC([N_model, Ntrt], 0.05);
                CCC_vec(i,2) = f_CCC([phi_model, phitrt], 0.05);
                CCC_vec(i,3)=f_CCC([pfset', pbestf'], 0.05);
        end
            
            % make the vectors that keep track of the accuracy for each tbot length
            nfit = size(pfitstore,2);
            ndata = size(Ntrtstore,1);
            nphi = length(phi_model);
            pgiven = reshape(pfitstore, [nsamps*nfit, 1]);
            pfit = reshape(pfittrt, [nsamps*nfit, 1]);
            Ngiven = reshape(Ntrtstore, [nsamps*ndata, 1]);
            Nfit = reshape(Nfittrt, [nsamps*ndata,1]);
            phigiven = reshape(phitrtstore(itime,:), [nsamps*nphi,1]);
            phifit = reshape(phifittrt, [nsamps*nphi,1]);
            index = isnan(pfit);
    CCC_p(m) = f_CCC([pgiven(index==0), pfit(index==0)], 0.05);
    ind = isnan(Nfit);
    CCC_N(m)= f_CCC([Ngiven(ind==0), Nfit(ind==0)], 0.05);
    ind = isnan(phifit);
    ind2 = isnan(phigiven);
    indall = or(ind,ind2);
    CCC_phi(m) = f_CCC([phigiven(indall==0), phifit(indall==0)], 0.05);

    param_err = [];
    % Calculate average percent parameter error
    index = isnan(pfit);
    param_err = 100*abs(pgiven(index==0)-pfit(index==0))./pgiven(index==0);

    median_param_err(m) = median(param_err);
    % Find the CCC and the param error for each parameters
    perr = [];
        for k = 1:nfit
            index = [];
            ind = [];
            pg = pfitstore(:,k); 
            pf = pfittrt(:,k);
            ind = isnan(pf);
            CCC_pmat(k,m) = f_CCC([pg(ind==0), pf(ind==0)], 0.05);
            perr = 100*abs(pg(ind==0)-pf(ind==0))./pf(ind==0);
            index = perr == Inf;
            perr= perr(index==0);
            avg_perr(k) =median(perr);


        end
        
       mean_param_err(m) = mean(avg_perr);
  
        
end
%% Plot some things

figure;
subplot(1,3,1)
plot(lambdavec, CCC_p, 'ro-', 'LineWidth',3)
xlabel('\lambda')
ylabel('CCC_{parameters}')
xlim([min(lambdavec) max(lambdavec)])
ylim([0.9 1])
title('Parameter concordance')
set(gca,'FontSize',20,'LineWidth',1.5)
subplot(1,3,2)
plot(lambdavec, CCC_N, 'bo-', 'LineWidth',3)
xlim([min(lambdavec) max(lambdavec)])
xlabel('\lambda')
ylabel('CCC_{N}')
ylim([0.9 1])
title('N(t) concordance')
set(gca,'FontSize',20,'LineWidth',1.5)
subplot(1,3,3)
plot(lambdavec, CCC_phi, 'go-', 'LineWidth',3)
xlabel('\lambda')
xlim([min(lambdavec) max(lambdavec)])
ylim([0.9 1])
ylabel('CCC_{\phi}')
title('\phi(t) concordance')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(lambdavec, mean_param_err,'ko-', 'LineWidth',3)
xlabel('\lambda')
ylabel('Average parameter error')
xlim([min(lambdavec) max(lambdavec)])
ylim([14 17])
title('Average parameter error vs. \lambda')
set(gca,'FontSize',20,'LineWidth',1.5)




figure;
for i = 1:nfit
    plot(lambdavec, CCC_pmat(i,:), '-o', 'LineWidth',2)
    hold on
end
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('\lambda')
xlim([min(lambdavec) max(lambdavec)])
ylim([0 1.4])
ylabel('CCC for each parameter')
title('CCC vs. \lambda')
legend('rs', '\alpha', 'rr', 'ds') 
legend boxoff
