% This script loads in the traj structure which contains the dose response
% for all conditions for the 231 data. It creates the trajsum structure which will contain a
% mean vector and U(t) vector for each different concentration tested. This
% will be used as a first pass to calibrate to the model.

% Adjusted this to include phi0 from scRNAseq estimates, but still only use
% N(t)
 close all; clear all; clc
 
 %% Load in data structure 
S = load('../out/trajfit231.mat');
traj= S.traj;
date = {'5-6-19'}; % date of 231 single pulse experiment
for j = 1:length(traj)
    if traj(j).numdoses == 1 && contains(traj(j).date, date)
        doses(j) = traj(j).dose;
    end
end
uniqdose= unique(doses);
% Make a new structure which combines each dox concentration
    % percapita growth rate 
    % variance of per capita growth rate
 %%
 % find groups by dose
 for i = 1:length(uniqdose)
    trajsum(i).Cdox = [];
    trajsum(i).Nmat = [];
    trajsum(i).nreps = 0;
    trajsum(i).tmat = [];
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
                end
        end
    end
 end
% Plot the raw data
 for i = 1:length(trajsum)
     for j = 1: trajsum(i).nreps
         plot(trajsum(i).tmat(:,j), trajsum(i).Nmat(:,j), 'color', trajsum(i).color)
         hold on
     end
 end
 xlabel('time (hours)')
 ylabel('N(t)')
 title('N(t) for different single pulse treatments')
 % Again, clean data for fitting...
 for i = 1:length(trajsum)
    if i ==1
    Nfin = 5e4;
    else
    Nfin = 3.5e4;
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
  %% Plot the cleaned data
  figure;
 for i = 1:length(trajsum)
     for j = 1: trajsum(i).nreps
         plot(trajsum(i).tfit(:,j), trajsum(i).Nfit(:,j), 'color', trajsum(i).color)
         hold on
     end
 end
 xlabel('time (hours)')
 ylabel('N(t)')
 title('N(t) for different single pulse treatments')
 %% Test and set U(t) curves
% Set some arbitrary things as inputs

kdrug = 0.0175;
k = 0.5;
dt = 1;
Cdoxmax = 1000;
% input time vectors for each different dose response
figure;
for i = 1:length(trajsum)
    ttest = [];
    ttest = 0:dt:trajsum(i).tvec(end);
    Cdox = trajsum(i).Cdox;
    trajsum(i).U = k*Cdox*exp(-kdrug*(ttest))/(0.1*Cdoxmax); 
    subplot(1, length(trajsum),i)
    plot(ttest, trajsum(i).U, 'b-', 'LineWidth',2)
    ylim([0 1])
    xlim([ 0 ttest(end)])
    xlabel('time (hours)')
    ylabel('Effective dose (U(t))')
    title([num2str(trajsum(i).Cdox),' nM Dox'])
end
%% Now find mean and standard deviation vectors

for i = 1:length(trajsum)
    trajsum(i).Nmean = mean(trajsum(i).Nfit,2);
    trajsum(i).tvec = round(trajsum(i).tfit(:,1),0);
    trajsum(i).Nstd = std(trajsum(i).Nfit,0,2);
end
% Plot the average data 
  figure;
 for i = 1:length(trajsum)
     subplot(2,1,1)
         plot(trajsum(i).tvec, trajsum(i).Nmean, 'color', trajsum(i).color, 'LineWidth', 2)
         hold on
         text(trajsum(i).tvec(end-10), trajsum(i).Nmean(end-10), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'])
         plot(trajsum(i).tvec, trajsum(i).Nmean + 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
         plot(trajsum(i).tvec, trajsum(i).Nmean - 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title('N(t) for different single pulse treatments')
         set(gca,'FontSize',20,'LineWidth',1.5)
        dt = 1;
       subplot(2,1,2)
       ttest = [];
       ttest = 0:dt:trajsum(i).tvec(end);
       plot(ttest, trajsum(i).U,'.', 'color',trajsum(i).color, 'LineWidth',1)
        hold on
        xlabel('time (hours)')
        ylabel('Effective dose U(t)')
        title('U(t) for different single pulse treatments')
         set(gca,'FontSize',20,'LineWidth',1.5)
 end


%% Test forward model
% set parameters of forward model
carcap = 5.2409e4;
S0=0.9*trajsum(1).Nmean(1); % set initial conditions (assume all cells sensitive to start)
R0=0.1*trajsum(1).Nmean(1); 
rs = 0.0287;
rr = 0.001;
ds = 0.2;
dr = 0;
alpha =0.1; %0.0001;


p = [ S0, R0, rs, carcap, alpha, rr, ds, dr];


Nsr = [];
tsum = [];
ttestsum = [];
tdrug = 1; % since we are monitoring first treatment only
ttest = [];
tvec = [];

for i = 1:6
    tvec = round(trajsum(i).tvec,0);
    U = trajsum(i).U;
    [Nsri, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U, dt, tdrug);
    tsum = vertcat(tsum, tvec);
    Nsr = vertcat(Nsr, Nsri);
end



figure;
plot(tsum, Nsr(:,1), 'b*')
hold on
plot(trajsum(1).tvec, trajsum(1).Nmean, 'r*')
xlabel('time (hours)')
ylabel('Model predicted response to pulse treat')
title('Test ability to generate model trajectories')
%% Fit untreated control data by calling separate function 
% grab untreated data
sigma = trajsum(1).Nstd(1:end);
ydataf = trajsum(1).Nmean;
ytimef = trajsum(1).tvec;


% Run through fitting function, get parameters (gs & carcap) and simulated
% forward model
[punt, singexpmodel] = fit_untreated(ydataf,ytimef, sigma);
gtot = punt(1);
carcap = punt(2);
    figure;
    plot(ytimef, ydataf, '*', 'LineWidth', 3)
    hold on
    plot(ytimef, singexpmodel,'-', 'LineWidth', 3')
    plot(ytimef, ydataf + 1.96*sigma, 'b-')
    plot(ytimef, ydataf - 1.96*sigma, 'b-')
    xlabel ('time (hours)')
    ylabel(' N(t)')
    legend ('Untreated control data', 'model fit', '95% CI on data', 'Location', 'NorthWest')
    legend box off
    title ('Fit for g_{tot} & K using untreated control')
    set(gca,'FontSize',20,'LineWidth',1.5)

 %% Now use this and fit for gr, ds, and alpha from single pulse treatments

%dr=0;
%phi0 = 0.92; % Set this from the scRNAseq from pre-treatment
rstar= phi0/(1-phi0);
pset = [carcap];
%params = [phi, rs, carcap, alpha, rr, ds, dr];
% Now we want to allow dr to be non-zero
psetID = [3];

%% Bayesian fit for alpha, rr, and ds using all treatments

% The simplified  model looks like this
%
% THE MODEL:
% dS/dt = rs(1-(S+R)/K)*S - alpha*u(t)*S - ds*u(t)*S
% dR/dt = rr(1-(S+R)/K)*R + alpha*u(t)*S- dr*u(t)*R ;
% 
% We will fit N= S +R trajectories for doses 10, 20, 35, 50, & 75
% Make vectors of inputs from 2:6 of traj sum corresponding to doses desired to
% fit
sigmafit = [];
ydatafit = [];
ytimefit = [];
N0s = [];
lengtht = [];
lengthU = [];
Uvec = [];

% fit on doses 10, 20, 35, 50 & 75 nM dox
% Grab the stdev, ydata, initial cell number, and U from the traj data
% structure 
dosefit = [1,4,7,9];
for j = 1:4%length(trajsum)
    i = dosefit(j);
sigmafit = vertcat(sigmafit,trajsum(i).Nstd(1:end));
ydatafit = vertcat(ydatafit, trajsum(i).Nmean(1:end));
ytimefit = vertcat(ytimefit, trajsum(i).tvec(1:end));
N0s = vertcat(N0s,trajsum(i).Nmean(1));
lengtht = vertcat(lengtht, length(trajsum(i).Nmean));
lengthU = vertcat(lengthU, length(trajsum(i).U));
Uvec = vertcat(Uvec, trajsum(i).U');
end
lengthvec = horzcat(lengtht, lengthU);

% Declare the parameters to be fit and give them an initial guess

pfID = [ 1, 2, 4, 5, 6, 7];
phi0guess = 0.9;
alphaguess =  0.1;
rrguess = 1e-3;
rsguess = ((rstar+1)*gtot - rrguess)./rstar; % expression that relates gs, gr, and gtot
dsguess = 0.2;
drguess = 0.1*dsguess;
pfitguess = [phi0guess, rsguess, alphaguess, rrguess, dsguess, drguess]; 

%test out the initial guess parameters 
modelfun = @(p)simmodelgreene(p, ytimefit, N0s, pset, Uvec, lengthvec, pfID, psetID); % single exponential model with death  


% Test out your forward function
figure;
plot(ytimefit, modelfun(pfitguess), 'r.','LineWidth', 2)
hold on
plot(ytimefit, ydatafit, 'b*', 'LineWidth', 3)
plot(ytimefit, ydatafit + 1.96*sigmafit, 'g.')
plot(ytimefit, ydatafit - 1.96*sigmafit, 'g.')
xlabel ('time (hours)')
ylabel(' N(t)')
legend ('first guess', 'datat')
legend box off
title ('Test guess for fit for rs, \alpha, rr and ds')


%% Try to run your fitting function
% Set some guesses for your phitrt and phisigfit
%phitrt = [0.1, 0.7];
%phisigfit = [0.1, 0.1];
Ntrt = ydatafit;
% set bounds for rs, alpha, rr, ds
pbounds = [ 0,1; 0,0.08; 0, 1; 0, 1; 0,1; 0,1];
%lambda = 0.9;
%N0phi = [phi0*N0, (1-phi0)*N0];
%tbot = [30 1008];
% Need Ub and lengthvecphi-- which are the effective dose and the size of
% the data for the phi data aspects


%[pbestf,N_model, negLL, pbestGD, N_modelGD, negLLGD]
% Here is where we will replace fit_fix_Greene with fit_fxn_Greenephi_N
%[pbestf,N_model, phi_model, negLL] = fit_fxn_Greenephi_N(Ntrt,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);
%(ydatafit, sigmafit, pfID, psetID, pfit, pset, time, Uvec, lengthvec, N0s, pbounds)
[pbestf,N_model, negLL] = fit_fxn_Greene(ydatafit,sigmafit, pfID, psetID, pfitguess, pset, ytimefit, Uvec, lengthvec, N0s, pbounds);
CCC = corrcoef(N_model,ydatafit);
phi0 = pbestf(1);
rs = pbestf(2);
alpha = pbestf(3);  
rr = pbestf(4);
ds = pbestf(5);
dr = pbestf(6);
%%
chi_sq = sum(((N_model-ydatafit)./sigmafit).^2)
%chi_sqGD = sum(((N_modelGD-ydatafit)./sigmafit).^2)
%%    
    figure;
    plot(ytimefit, ydatafit, 'b*', 'LineWidth', 3)
    hold on
    plot(ytimefit, N_model,'r*', 'LineWidth', 3')
    plot(ytimefit, ydatafit + 1.96*sigmafit, 'g.')
    plot(ytimefit, ydatafit - 1.96*sigmafit, 'g.')
    xlabel ('time (hours)')
    ylabel(' N(t)')
    legend ('data from mult treatments', 'model fit')
    legend box off
    title ('Fit for \alpha, rr, ds and dr')
    set(gca,'FontSize',20,'LineWidth',1.5)
%% Plot calibrated data
N0 = 2e3;
tdrug = 1;
p= [phi0*N0, (1-phi0)*N0, rs, carcap, alpha, rr, ds, dr];
 figure;
 for j = 1:4%length(trajsum)
     i = dosefit(j);
     subplot(2,1,1)
         plot(trajsum(i).tvec, trajsum(i).Nmean, 'color', trajsum(i).color, 'LineWidth', 2)
         hold on
         plot(trajsum(i).tvec, trajsum(i).Nmean + trajsum(i).Nstd, 'color', trajsum(i).color)
         plot(trajsum(i).tvec, trajsum(i).Nmean - trajsum(i).Nstd, 'color', trajsum(i).color)
       
         xlabel('time (hours)')
        ylabel('N(t)')
        title('Model calibration Trial Run N(t) only')
        dt = 1;
        tvec = [];
        Nsri = [];
        tvec = trajsum(i).tvec;
        U = trajsum(i).U;
         pi = p;
         pi(1) = phi0*trajsum(i).Nmean(1);
         pi(2) = (1-phi0)*trajsum(i).Nmean(1);
        [Nsri, tcrit, Ncrit] = fwd_Greene_model(pi, tvec, U, dt, tdrug);
        plot(tvec, Nsri(:,1), 'color','r', 'LineWidth',2)
        %plot(tcrit, Ncrit, '*', 'LineWidth',2)
        %text(tcrit +5, Ncrit + 5, ['t_{crit}=', num2str(tcrit),' hours'])
        text(trajsum(i).tvec(end-10), trajsum(i).Nmean(end-10), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'], 'FontSize', 14)
        set(gca,'FontSize',20,'LineWidth',1.5)
        trajsum(i).Nmod1pulse = Nsri;

        subplot(2,1,2)
       ttest = [];
       ttest = 0:dt:trajsum(i).tvec(end);
       plot(ttest, trajsum(i).U,'.', 'color',trajsum(i).color, 'LineWidth',1)
        hold on
        xlabel('time (hours)')
        ylabel('Effective dose U(t)')
        title('Effective dose of each pulse treatment')
        set(gca,'FontSize',20,'LineWidth',1.5)
 end
 %% Save calibrated parameters from 231s with 4 doses 25-150 nM
 
p4fit= [rs, carcap, alpha, rr, ds, dr, k ,kdrug, gtot];
save('../out/p4fit231', 'p4fit')
% Save trajsum
%p= [N0, 0, rs, carcap, alpha, rr, ds, dr];
for i = 1:length(trajsum)
    pi = p;
    pi(1) = trajsum(i).Nmean(1);
    tvec = [];
    Nsri = [];
    tvec = trajsum(i).tvec;
    U = trajsum(i).U;
    [Nsri, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U, dt, tdrug);
   trajsum(i).Nmod1pulse = Nsri;
    trajsum(i).params = pi;
    trajsum(i).N0 = trajsum(i).Nmean(1);
    trajsum(i).rs= rs;
    trajsum(i).carcap = carcap;
    trajsum(i).alpha = alpha;
    trajsum(i).rr = rr;
    trajsum(i).ds = ds;
    trajsum(i).dr = 0;
    trajsum(i).kdrug =kdrug;
end

%% Save the parameter estimates
p4fit= [rs, carcap, alpha, rr, ds, dr, k, kdrug];
save('../out/p4fit231', 'p4fit')
save('../out/trajsumfit231.mat', 'trajsum')