% This script loads in Grant's raw dose response data and begins performing
% some preliminary visualization/fitting

close all; clear all; clc;
%% Load in data structure 
S = load('../out/trajraw.mat');
traj= S.traj;


%% Flip through raw data trajectories
% ctrl C to end this loop without iterating through all!
figure;
for j = 400:length(traj)
    plot(traj(j).time, traj(j).rawN, 'LineWidth', 2, 'color', num2str(traj(j).color))
    xlabel ('time (hours)')
    ylabel('N(t)')
    title(['Raw data, dose =' num2str(traj(j).dose), 'nM, WPT=', num2str(traj(j).WPT)])
    pause
end
%% Plot trajectories together by color

figure;
for j = 1:60%length(traj)
    plot(traj(j).time, traj(j).rawN, 'LineWidth', 1, 'color', num2str(traj(j).color))
    hold on
    xlabel ('time (hours)')
    ylabel('N(t)')
    title('Raw data colored by WPT')
    xlim([0,traj(j).time(end)])
end
% Look at variation among replicates
for j = 1:length(traj)
    if ~isempty(traj(j).WPT)
    WPT(j) = traj(j).WPT;
    end
    if traj(j).numdoses == 1
        doses(j) = traj(j).dose;
    end
end
colorsets = varycolor(length(unique(WPT))+1);
uniqWPT= unique(WPT);
uniqdose= unique(doses);

figure;
for i = 1:length(traj)
    for j = 1:length(uniqWPT)
    if traj(i).WPT == uniqWPT(j)
        if traj(i).dosenum == 2
        subplot(1,length(uniqWPT), j)
        plot(traj(i).time, traj(i).rawN, 'color', num2str(traj(i).color))
        xlabel ('time (hours)')
        ylabel('N(t)')
        title([num2str(uniqWPT(j)), ' WPT'])
        hold on
        end
    end
    end
end


%% Now want to "clean" data for fitting
% What does this mean? This is where we can make decisions, sort of
% arbitrarily/ might be deata set dependent. 

% For the first data set, let's set N0 at t0 as the time after drug
% We'll set Nend as either the last data point or when the cell number
% exceeds 2e4
Nfin=2e4;
for i = 1:length(traj)
    N = traj(i).rawN;
    t = traj(i).time;
    i0 = find(t>traj(i).tdose,1,'first'); % I arbitrarily search for a maximum in the first 200 hours
    N0 = N(i0); % set N0 as maximum over first 200 hours
    iend = find(N>=Nfin,1, 'first');
    if ~isempty(iend)
    tfit = t(i0:iend)-t(i0); 
    Nfit = N(i0:iend);
    end
    if isempty(iend)
        tfit = t(i0:end)-t(i0); 
        Nfit = N(i0:end);
    end
    traj(i).tfit =tfit;
    traj(i).Nfit =Nfit;
end
%% Plot data to be fit
figure;
for i = 1:length(traj)
plot(traj(i).tfit, traj(i).Nfit, 'color', num2str(traj(i).color))
hold on
end
xlabel('time (hours)')
ylabel('N(t)')
title('Data for fitting')

figure;
for i = 1:length(traj)
    for j = 1:length(uniqWPT)
    if traj(i).WPT == uniqWPT(j)
        if traj(i).dosenum == 2
        subplot(1,length(uniqWPT), j)
        plot(traj(i).tfit, traj(i).Nfit, 'color', num2str(traj(i).color))
        xlabel ('time (hours)')
        ylabel('N(t)')
        title([num2str(uniqWPT(j)), ' WPT'])
        hold on
        end
    end
    end
end

%% Find tcrit if present
% Set NCRIT Here
for i = 1:length(traj)
    icrit = [];
    N0 = traj(i).Nfit(1);
    N = traj(i).Nfit;
    tfit = traj(i).tfit;
    Ncrit= 1.2*N0;
    traj(i).Ncrit = Ncrit;
    icrit = find(N>Ncrit,1, 'first');
    if ~isempty(icrit)
    tcrit= tfit(icrit);
    traj(i).tcrit = tcrit;
    end
    if isempty(icrit)
        traj(i).tcrit = [];
    end
end

%% Plot data for single dose response
figure;
for i = 1:length(traj)
    for j = 1:length(uniqdose)
    if traj(i).dose == uniqdose(j)
        if traj(i).numdoses ==1 || traj(i).numdoses == 0
        subplot(1,length(uniqdose), j)
        plot(traj(i).tfit, traj(i).Nfit, 'color', num2str(traj(i).color))
        hold on
        if ~isempty(traj(i).tcrit)
        plot(traj(i).tcrit, traj(i).Ncrit, 'k*', 'LineWidth', 3)
        end
        xlabel ('time (hours)')
        ylabel('N(t)')
        title([num2str(uniqdose(j)), ' nM'])
        hold on
        end
    end
    end
end

figure;
for i = 1:length(traj)
    for j = 1:length(uniqdose)
    if traj(i).dose == uniqdose(j)
        if traj(i).numdoses ==1 || traj(i).numdoses == 0
        if ~isempty(traj(i).tcrit)   
        plot(traj(i).dose, traj(i).tcrit, 'k*')
        hold on
        end
        xlabel('Concentration Dox','FontSize',20)
        ylabel('t_{crit}','FontSize',20)
        set(gca,'FontSize',20,'LineWidth',1.5)
        title('Dox Concentration vs. Critical time')
        hold on
        end
    end
    end
end

%% Use the new N and T vectors to perform fitting
% The model looks like this
%
% $$ N(t) = N_0 [ \phi e^{gt} + (1-\phi) e^{-kt} ] $$
% 
% where:
%
% * $N_0$ is the initial cell number
% * $\phi$ is the initial resistant fraction
% * g>0 is the resistant growth rate
% * k>0 is the kill rate (actually -net "growth rate" on treatment) on sensitive cells
%

% Define transforms 
% single exponential
pfxform1 = @(pval)[1].*log(pval); %'forward' parameter transform into Reals
pbxform1 = @(phat)[1].*exp(phat);  %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% double exponential
pfxform2 = @(pval)[0 1 1].*log(pval)+[1 0 0].*log(pval./(1-pval)); %'forward' parameter transform into Reals
pbxform2 = @(phat)[0 1 1].*exp(phat)+[1 0 0].*(exp(phat)./(1+exp(phat)));  %'backward' parameter transform into model space

sigma = 5; % this should be data driven! i.e. what is our confidence in samples 
for j = 1:length(traj)
    ydata = traj(j).Nfit;
    N0 = ydata(1);
    ydata = ydata(2:end);
    ytime = traj(j).tfit;
    ytime = ytime(2:end);
    
    % Set up forward models, fit all three nested versions of model
    modelfungoodd = @(p)simmodeld(p, ytime, N0); % single exponential model with death  
    modelfungoodg = @(p)simmodelg(p, ytime, N0); % single exponential model with growth
    modelfungood2 = @(p)simmodel2(p,ytime, N0); % double exponential model function with ytime
     
    % INITIAL GUESSES BASED ON DATA

        kguess = -(yfxform(ydata(5)) - yfxform(ydata(1)))/(ytime(5));
   if kguess <= 0
        kguess = 1e-5;
    end
    gguess = (yfxform(ydata(end))-yfxform(ydata(end-5)))/(ytime(end)-ytime(end-5)); 
    % alter initial guesses to prevent NaNs and zeros
    if isnan(gguess)
        gguess = 1e-5;
    end
    if isinf(gguess)
        gguess = 0.8;
    end
    if gguess <= 0
        gguess = 1e-5;
    end
    phiguess = (ydata(end)/ydata(1)).*exp(-gguess*ytime(end));
    if phiguess <= 0
       phiguess = 1e-5;
    end
    if phiguess>=0
        phiguess = 0.99;
    end
    % Initial guess matrices
    theta1 = [kguess]; % and k
    thetag = [gguess];
    theta2 = [phiguess, gguess, kguess];
    
    % Write log likelihood function based on assumption of normally
    % distributed sampling error
    
    % Goal: maximize the probability of the data given the model. Normpdf
    % will output a probability of the data (x- 1st argument), given the
    % mean(expectation, here the model), and the variance, at each time point. 
    % take the log and minimize the NLL to maximize likeihood of data given
    % the model
    
    % single exponential with death
    loglikelihoodd = @(phat)sum(log(normpdf(yfxform(ydata),yfxform(modelfungoodd(pbxform1(phat))), sigma)));
    % single exponential with growth
    loglikelihoodg = @(phat)sum(log(normpdf(yfxform(ydata),yfxform(modelfungoodg(pbxform1(phat))), sigma)));
    % double exponential
    loglikelihood2 = @(phat)sum(log(normpdf(yfxform(ydata),yfxform(modelfungood2(pbxform2(phat))), sigma)));
  
    % Write objective functions for each model
    objfund = @(phat)-loglikelihoodd(phat);
    objfung = @(phat)-loglikelihoodg(phat);
    objfun2 = @(phat)-loglikelihood2(phat);
    phatbestd = fminsearch(objfund, pfxform1(theta1)); 
    phatbestg = fminsearch(objfung, pfxform1(thetag));
    phatbest2 = fminsearch(objfun2, pfxform2(theta2));
    
    % Save some stuff
    % save biexp stuff (most relevant so will put first)
    traj(j).params2 = pbxform2(phatbest2);
    traj(j).biexpmodel = simmodel2(pbxform2(phatbest2), ytime, N0); % mode
    residuals2 = traj(j).biexpmodel-ydata;
    ybar = mean(ydata);
    Rsq2 = 1-(sum((residuals2).^2)./(sum((ybar-ydata).^2)));
    traj(j).Rsq2 = Rsq2;
    num_params2 = 3;
    n = length(ytime);
    AIC2 = -2*loglikelihood2(phatbest2) + 2*num_params2;
    traj(j).AIC2 = AIC2;
    if isempty(traj(j).tcrit)
        ext = 1000;
        text = vertcat(ytime, (round(ytime(end),0):4:round(ytime(end),0)+ext)');
        traj(j).biexpmodel = simmodel2(pbxform2(phatbest2), text, N0); % model
        traj(j).tmod = text;
        icrit = find(traj(j).biexpmodel>traj(j).Ncrit, 1, 'first');
        traj(j).tcrit= text(icrit);
    end
    
    % save single exponential death
    traj(j).paramsd = pbxform1(phatbestd);
    traj(j).expmodeld = simmodeld(pbxform1(phatbestd), ytime, N0); % mode
    residualsd = traj(j).expmodeld-ydata;
    ybar = mean(ydata);
    Rsqd = 1-(sum((residualsd).^2)./(sum((ybar-ydata).^2)));
    traj(j).Rsqd = Rsqd;
    num_paramsd = 1;
    n = length(ytime);
    AICd = -2*loglikelihoodd(phatbestd) + 2*num_paramsd;
    traj(j).AICd = AICd;
    
    % save single exponential growth
    traj(j).paramsg = pbxform1(phatbestg);
    traj(j).expmodelg = simmodelg(pbxform1(phatbestg), ytime, N0); % mode
    residualsg = traj(j).expmodelg-ydata;
    ybar = mean(ydata);
    Rsqg = 1-(sum((residualsg).^2)./(sum((ybar-ydata).^2)));
    traj(j).Rsqg = Rsqg;
    num_paramsg = 1;
    n = length(ytime);
    AICg = -2*loglikelihoodg(phatbestg) + 2*num_paramsg;
    traj(j).AICg= AICg;
    
    % Code for best fitting model
    AICs = [AIC2, AICd, AICg];
    RSQs = [Rsq2, Rsqd, Rsqg];
    % 1 is biexp model
    if min(AICs)==AIC2
        traj(j).bfmod=1;   
    end
    % 2 is single exponential death model
    if min(AICs)==AICd
        traj(j).bfmod =2;
    end
    % 3 is single exponential growth model
    if min(AICs)==AICg
        traj(j).bfmod =3;
    end
    % 0 is bad fit
    if max(RSQs) < 0.6
        traj(j).bfmod = 0;
    end
    
end


%% Flip through model fit compared to data
% Show best fitting "chosen" model

figure;
for j = 120:length(traj)
    hold off
    plot(traj(j).tfit, traj(j).Nfit, 'LineWidth', 2, 'color', num2str(traj(j).color))
    hold on
    
        if isempty(traj(j).tmod)
    plot(traj(j).tfit(2:end), traj(j).biexpmodel, 'LineWidth', 2, 'color', 'k')
        end
        if ~isempty(traj(j).tmod)
     plot(traj(j).tmod, traj(j).biexpmodel, 'LineWidth', 2, 'color', 'k')   
        end
    title(['Data fit to biexp model, \phi=', num2str(traj(j).params2(1)),', g=', num2str(traj(j).params2(2)),', k=', num2str(traj(j).params2(3))])
   
    
%     if traj(j).bfmod ==2
%     plot(traj(j).tfit(2:end), traj(j).expmodeld, 'LineWidth', 2, 'color', 'k')
%     title(['Data fit to exp death model, k=', num2str(traj(j).paramsd(1))])
%     end
%     
%     if traj(j).bfmod ==3
%     plot(traj(j).tfit(2:end), traj(j).expmodelg, 'LineWidth', 2, 'color', 'k')
%     title(['Data fit to exp growth model, g=', num2str(traj(j).paramsg(1))])
%     end
    
    if ~isempty(traj(j).tcrit)
        plot(traj(j).tcrit, traj(j).Ncrit, 'k*')
    end
    
%     xlim( [ 0 traj(j).tfit(end)])
    xlabel ('time (hours)')
    ylabel('N(t)')
    pause
end
%% Critical time vs dox concetration
figure;
for i = 1:length(traj)
    for j = 1:length(uniqdose)
    if traj(i).dose == uniqdose(j)
        if traj(i).numdoses ==1 || traj(i).numdoses == 0
        if ~isempty(traj(i).tcrit)   
        plot(traj(i).dose, traj(i).tcrit, 'k*')
        hold on
        end
        xlabel('Concentration Dox','FontSize',20)
        ylabel('t_{crit}','FontSize',20)
        set(gca,'FontSize',20,'LineWidth',1.5)
        title('Dox Concentration vs. Critical time')
        hold on
        end
    end
    end
end


%% Record model parameters of interest
for j = 1:length(traj)
    if traj(j).bfmod ==1
        params = traj(j).params2;
        traj(j).phi = params(1);
        traj(j).g = params(2);
        traj(j).k = params(3);
    end
    
      if traj(j).bfmod ==2
        params = traj(j).paramsd;
        traj(j).phi = 1;
        traj(j).g = 0;
        traj(j).k = params(1);
      end
     if traj(j).bfmod ==3
        params = traj(j).paramsg;
        traj(j).phi = 0;
        traj(j).g = params(1);
        traj(j).k = 0;
     end
     if traj(j).bfmod ==4
        traj(j).N0fit = [];
        traj(j).phi = [];
        traj(j).g = [];
        traj(j).k = [];
     end
    

end

%% Save the fitted data structure
% this saves the fitted data structure, obtained from the raw data
% structure (run load_raw_data.m)
save('../out/trajfit.mat', 'traj')