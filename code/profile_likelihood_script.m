%function [profiles] = profile_likelihood(params, tsamp, N0, N,mudatavec, vardatavec, modelcode, factor, numpoints)
%fit_fxn_Greenephi_N2(ydatafit,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);
%% For CLL
Ub = Uphi;
lengthvec = lengthvecN;
%% RUN this to perform profile likelihood on joint calibration (N(t) and phi(t))
nfitparams = length(pfitID);
        profile = [];
factor0 = 1.3;
numpoints = 30;
params = pbest; % rs, alpha, zr, ds, zd];
threshold = chi2inv(0.95,length(pfitID))/2 + negLL;
profiles = [];


pboundsprof = [0, Inf; 0, Inf; 0, 1; 0, Inf; 0, Inf; 0, 1];
profiles = [];
%%
ikeep = [];
for k = 1:nfitparams
profindex = pfitID(k); % PROFILE the kth fit parameter 
    factor = factor0;
    if k ==1
        factor = 0.2*factor0;
    end
    if k == 5 % increase factor for ds and dr/ds ratio 
        factor = 2*factor0;
    end
    if k==6
        factor = 8*factor0;
    end
    profrangeDown = linspace((params(k)), (params(k)*(1-factor)),numpoints)'; 
    profrangeUp = linspace((params(k)), (params(k)*(1+factor)),numpoints)';
    % split into up and down so we can use last fitted value as starting value for next run
    profrange = [profrangeDown;profrangeUp];
    profrange = sort(profrange);
    currfvals = [];
    currparams = [];
    currflags = [];
    paramstemp = [];
    profile = [];

    for m = 1:length(profrange)
        [m] %track progress
        currp = profrange(m);
        % Change pfitID and psetID so that pfit is all params except
        % currently profiled param and psetID includes currently profiled
        % param
        [psetIDcurr, ordp] = sort([psetID, profindex]);
        ikeep = pfitID~=profindex; % identify all indices corresponding to
        % parameters not being profiled
        pfitIDcurr = pfitID(ikeep);
        thetacurr = pbest(ikeep);
        pboundscurr = pboundsprof(ikeep, :);
        pcomb = [pset, currp];
        psetcurr = pcomb(ordp);
        [paramstemp,~, ~, fvaltemp, ~, ~]=fit_fxn_Greenephi_Nprof(ydatafit,sigmafit,phitrt, phisigfit, pfitIDcurr, psetIDcurr, thetacurr, psetcurr, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pboundscurr);
        
        %[fvaltemp, paramstemp] = ML_fitnegLL(params, tsamp, N0, N, mudatavec, vardatavec, modelcode, profindex, currp);
        % fminsearch will out put the values of dguess that give the lowest
        % objfun_donly
      
        currfvals = [currfvals; fvaltemp];
        %currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
        currparams = [currparams; [paramstemp(1:k-1),currp,paramstemp(k:end)]];
    end
    
    profile = horzcat(profrange, real(currfvals), currparams);
 


    profiles(:,:,k) = profile;
% each profile has columns: profiled parameter value, resulting
        % cost-function (e.g. RSS, negLL) value, and
        % then columns for each of the other parameter estimates.
    

ilo = [];
 ilo=find(profiles(:,2,k)<=threshold);
    if isempty(ilo)
        display('No CI reached')
    else
        if ilo(end) ==numpoints*2 ||ilo(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CI(k,:) = [NaN, NaN];
        else
             CI(k,:) = [profiles(ilo(1)-1,1,k), profiles(ilo(end),1,k)];
        end
    end
end



%% Plot the profiles joint fitting 
plist = {'r_{s}', '\alpha', 'rr/rs', 'd_{s}', 'dr/ds'};
plist = {'\phi_{0}','r_{S}', 'r_{S}/r_{S} ratio', '\alpha', 'd_{S}', 'd_{R}/d_{S} ratio'};
threshold = chi2inv(0.95,length(pfitID))/2 + negLL;
figure;
for i = 1:length(pbest)
subplot(2,3,i)
plot([CI(i,1) CI(i,1)], [0 50], 'g-', 'LineWidth', 3)
hold on
plot([CI(i,2) CI(i,2)], [0 50], 'g-', 'LineWidth',3)
plot(profiles(:,1, i), profiles(:,2, i),'b-', 'LineWidth', 3)
hold on
plot(pbest(i), negLL, 'r*', 'LineWidth', 4)
plot([profiles(1,1,i) profiles(end,1,i)],[threshold threshold],'r--')
xlabel(plist(i))
ylabel('J(\theta)')
set(gca,'FontSize',24,'LineWidth',1.5)
legend('95% CI')
legend boxoff

%title(plist(i))
ylim([0 1.3*threshold])
%xlim([0.5*pbest(i) 2*pbest(i)])
% if i == 1
%      xlim([0.8 0.9])
% end
% else
xlim([0.7*CI(i,1),1.1*CI(i,2)])
% end
% xlim([CI(i,1)-0.1*pbest(i),CI(i,2)+0.1*pbest(i)])
% ylim([ 0 50])
end
%% Save the CI from the profile likelihood of N(t) and phi(t)
% uncertainty to model uncertainty
save('../out/CIpbest.mat', 'CI')


%% RUN THIS TO PLOT PROFILE LIKELIHOOD FROM N fitting only

nfitparams = length(pfitID);
        profile = [];
factor0 = 0.6;
numpoints = 30;
params = pbestN; % rs, alpha, zr, ds, zd];
threshold = chi2inv(0.95,length(pfitID))/2 + negLLN;
profiles = [];

pboundsprof = [-Inf, Inf; -Inf, Inf; 0, 1; -Inf, Inf; -Inf, Inf; 0, 1];
profiles = [];
ikeep = [];
for k = 1:nfitparams
profindex = pfitID(k); % PROFILE the kth fit parameter 
    factor = factor0;
    if k ==1
        factor = 0.1*factor0;
    end
    if k == 5 % increase factor for ds and dr/ds ratio 
        factor = 2*factor0;
    end
    if k==6
        factor = 3*factor0;
    end
    profrangeDown = linspace((params(k)), (params(k)*(1-factor)),numpoints)'; 
    profrangeUp = linspace((params(k)), (params(k)*(1+factor)),numpoints)';
    % split into up and down so we can use last fitted value as starting value for next run
    profrange = [profrangeDown;profrangeUp];
    profrange = sort(profrange);
    currfvals = [];
    currparams = [];
    currflags = [];
    paramstemp = [];
    profile = [];

    for m = 1:length(profrange)
        [m] %track progress
        currp = profrange(m);
        % Change pfitID and psetID so that pfit is all params except
        % currently profiled param and psetID includes currently profiled
        % param
        [psetIDcurr, ordp] = sort([psetID, profindex]);
        ikeep = pfitID~=profindex; % identify all indices corresponding to
        % parameters not being profiled
        pfitIDcurr = pfitID(ikeep);
        thetacurr = theta(ikeep);
        pboundscurr = pboundsprof(ikeep, :);
        pcomb = [pset, currp];
        psetcurr = pcomb(ordp);
        [paramstemp,~, ~, fvaltemp, ~, ~]=fit_fxn_Greenephi_Nprof(ydatafit,sigmafit,phitrt, phisigfit, pfitIDcurr, psetIDcurr, thetacurr, psetcurr, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambdatry, pboundscurr);
        
        %[fvaltemp, paramstemp] = ML_fitnegLL(params, tsamp, N0, N, mudatavec, vardatavec, modelcode, profindex, currp);
        % fminsearch will out put the values of dguess that give the lowest
        % objfun_donly
      
        currfvals = [currfvals; fvaltemp];
        %currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
        currparams = [currparams; [paramstemp(1:k-1),currp,paramstemp(k:end)]];
    end
    
    profile = horzcat(profrange, real(currfvals), currparams);
 


    profiles(:,:,k) = profile;
% each profile has columns: profiled parameter value, resulting
        % cost-function (e.g. RSS, negLL) value, and
        % then columns for each of the other parameter estimates.
    

ilo = [];
 ilo=find(profiles(:,2,k)<=threshold);
    if isempty(ilo)
        display('No CI reached')
    else
        if ilo(end) ==numpoints*2 ||ilo(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CI(k,:) = [NaN, NaN];
        else
             CI(k,:) = [profiles(ilo(1)-1,1,k), profiles(ilo(end),1,k)];
        end
    end
end



%% Plot the profiles for the N fitting only 
plist = {'r_{s}', '\alpha', 'rr/rs', 'd_{s}', 'dr/ds'};
plist = {'\phi_{0}','r_{s}', 'r_{r}/r_{s} ratio', '\alpha', 'd_{s}', 'd_{r}/d_{s} ratio'};
threshold = chi2inv(0.95,length(pfitID))/2 + negLLN;
figure;
for i = 1:length(pbest)
subplot(2,3,i)
%plot([CI(i,1) CI(i,1)], [90 threshold], 'g-', 'LineWidth', 3)
hold on
%plot([CI(i,2) CI(i,2)], [90 threshold], 'g-', 'LineWidth',3)
plot(profiles(:,1, i), profiles(:,2, i),'b-', 'LineWidth', 3)
hold on
plot(pbestN(i), negLLN, 'r*', 'LineWidth', 4)
plot([profiles(1,1,i) profiles(end,1,i)],[threshold threshold],'r--')
xlabel(plist(i))
ylabel('J(\theta)')
set(gca,'FontSize',16,'LineWidth',1.5)
%legend('95% CI')
%legend boxoff
if i ==2 || i ==3 || i ==5
    xlim([0.8*pbestN(i) 1.2*pbestN(i)]) 
else
xlim([profiles(1,1,i) profiles(end,1,i)])
end
ylim([negLLN-50, threshold+1e2])
%title(plist(i))
%ylim([75 90])
% if i == 1
%     xlim([CI(i,1)-.01,CI(i,2)+.01])
% else
% xlim([0.8*CI(i,1),1.1*CI(i,2)])
% end
%xlim([CI(i,1)-0.1*pbest(i),CI(i,2)+0.1*pbest(i)])
%ylim([ 90 120])
end

%% Plot parameter relationships from joint N(t) and phi(t) fitting
figure;
for i = 1:length(pbest)
    subplot(2,3,i)
    set(gca,'FontSize',20,'LineWidth',1.5)
    hold on
    plot(profiles(:,1,i),profiles(:,3:end,i),'LineWidth',2)
    plot(pbest(i),pbest,'r*')
    xlabel(plist{i})
    ylabel('Estimated Parameter Value')
    legend(plist, 'FontSize', 8)
    legend boxoff
    title(plist(i))
    ylim([0 1])
    xlim([CI(i,1) CI(i,2)])
end
%% Is there another good way to plot the parameters and their CI?
