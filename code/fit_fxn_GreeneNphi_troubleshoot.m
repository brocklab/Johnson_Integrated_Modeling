%function [pbest,model_N, negLL] = fit_fxn_GreeneNphi(ydatafit, sigmafit,phitrt, phisigfit, pfID, psetID, pfit, pset, time, tbot, Uvec, Ub, lengthvec, N0s,N0phi, pbounds) 
%[pbest,model_N, negLL, pbestGD, model_NGD, negLLGD]
% This function

% Delete this section once remake into a function
ydatafit = Ntrt;
time = ytimefit;
pfit = theta;
pfID = pfitID;
N0phi = 2e3;

%Write a fucntion that fits any combination of the paramters in the
%simplest Greene model with the following parameters:
%p = [ S0, R0, rs, carcap, alpha, rr, ds, dr];
% To determine which parameters are set and which ones are fit, pID tells
% us this with a 1 for a parameter that must be fit and a 0 for a parameter
% that is set.

% TEST

% Find pfit
% P = num2cell(pID); 
% Initialize parameters
prop = 0;
rs = 0;
carcap = 0;
alpha = 0;
rr = 0;
ds = 0;
dr = 0;
params = [prop, rs, carcap, alpha, rr, ds, dr];
for i = 1:length(params)
    indset= find(ismember(psetID, i));
    if ~isempty(indset)
    params(i) = pset(indset);
    end
    indfit = find(ismember(pfID,i));
    if ~isempty(indfit)
    params(i) = pfit(indfit);
    end
end
P = num2cell(params);
[prop, rs, carcap, alpha, rr, ds, dr] = deal(P{:});
 %%
% Define transforms 
% for number of variables in pfit
nfit = length(pfit);
pfxform = @(pval)ones(1,nfit).*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)ones(1,nfit).*exp(phat);  %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
yfxformp = @(y)log(y./(1-y));
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output
ybxformp = @(yhat) exp(yhat)./(1+exp(yhat));


theta = pfit; 

% write a function that takes set parameters and fit parameters, combines
% them, and runs the model forward for the Uvec and 0s provided
modelfunN = @(p)simmodelgreene(p, time, N0s, pset, Uvec, lengthvec, pfID, psetID); 
modelfunphi = @ (p)simmodelgreenephi(p, tbot, N0phi, pset, Ub, lengthvecphi, pfID, psetID);

%% Test your modelfunphi
% this is what your forward function should generate

model_phi = modelfunphi(pfset);
model_phi_theta = modelfunphi(theta);
model_N = modelfunN(pfset);
modelfit = modelfunN(pbest2);
modelfitphi = modelfunphi(pbest2);
% phitrt is the data generated from those same parameters
figure;
plot(tbot, phitrt, '*','color', 'g')
hold on
plot(tbot, model_phi, '-','color', 'b', 'LineWidth',2)
plot(tbot, modelfitphi, 'm-')
plot(tbot, model_phi_theta, '-','color', 'r', 'LineWidth',2)
plot(tbot, model_phi + 1.96*phisigfit, '-', 'color', 'k')
plot(tbot, model_phi - 1.96*phisigfit, '-', 'color', 'k')
xlabel('time(hours)')
ylabel('\phi(t)')
title('Model fxn versus model generated data')
legend('in silico data', 'model fxn', 'model fxn fit params','95% CI on data')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(ytimefit, ydatafit, '*','color', 'g')
hold on
plot(ytimefit, Ntrtstore(:,ind), 'k.')
plot(ytimefit, model_N, 'o','color', 'b', 'LineWidth',2)
plot(ytimefit, modelfit, 'ro')
xlabel('time(hours)')
ylabel('N(t)')
title('Model fxn versus model generated data')
legend('in silico data', 'in silico data', 'model fxn known params','model fxn fit params')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
%% Plot likelihood in time of data given model
% reset theta
theta2 = [pfset(1)+.03, pfset(2)-0.0002, pfset(3)+0.0003, pfset(4)+0.001]
%%
sigmafit2 = sigmafit;
loglikelihoodphi_t2= @(phat)(log(2*pi.*(phisigfit).^2)+((modelfunphi(pbxform(phat))-phitrt)./phisigfit).^2);
loglikelihoodN_t2= @(phat)(log(2*pi.*(sigmafit2).^2)+((modelfunN(pbxform(phat))-ydatafit)./sigmafit2).^2);
probs_phi = loglikelihoodphi_t2(pfxform(pfset))
probs_N = loglikelihoodN_t2(pfxform(pfset))

figure;
plot(tbot, loglikelihoodphi_t2(pfxform(pfset)), '*')
hold on
plot(ytimefit, loglikelihoodN_t2(pfxform(pfset)),'g.')
legend ('NLL phi', 'NLL N' )
xlabel('time(hours)')
ylabel('NLL(mod;data)')
title('NLL for \phi & N')
set(gca,'FontSize',20,'LineWidth',1.5)
phisum = sum(loglikelihoodphi_t2(pfxform(pfset)))
Nsum = sum(loglikelihoodN_t2(pfxform(pfset)))
%%
lambda = 0.5;

loglikelihoodphi2= @(phat)sum(log(2*pi.*(phisigfit).^2)+((modelfunphi(pbxform(phat))-phitrt)./phisigfit).^2);
loglikelihoodN2= @(phat)sum(log(2*pi.*(sigmafit2).^2)+((modelfunN(pbxform(phat))-ydatafit)./sigmafit2).^2);

objfun2 = @(phat)((1-lambda)*loglikelihoodN2(phat) + lambda*loglikelihoodphi2(phat));

phatbest2 = fminsearch(objfun2, pfxform(theta));
pbest2 = pbxform(phatbest2)
negLL2 = objfun(phatbest2)

figure;
plot(tbot, loglikelihoodphi_t2(pfxform(pfset)), '*')
hold on
plot(tbot, loglikelihoodphi_t2(pfxform(pbest2)),'*')
legend ('true params', 'fit params')
xlabel('time(hours)')
ylabel('NLL(mod|\phi)')
title('Probability at each time point')
set(gca,'FontSize',20,'LineWidth',1.5)
fitsum = sum(loglikelihoodphi_t(pfxform(pbest)))
truesum = sum(loglikelihoodphi_t(pfxform(pfset)))
figure;
plot(ytimefit, loglikelihoodN_t2(pfxform(pfset)), '*')
hold on
plot(ytimefit, loglikelihoodN_t2(pfxform(pbest2)),'*')
legend ('true params', 'fit params')
xlabel('time(hours)')
ylabel('NLL(mod|N)')
title('NLL at each time point')
set(gca,'FontSize',20,'LineWidth',1.5)
phicontrib = sum(-loglikelihoodphi_t(pfxform(pbest)))
Ncontrib = sum(-loglikelihoodN_t(pfxform(pbest)))
negLL = objfun(phatbest)
%%
loglikelihoodN_t = @(phat)((log(normpdf(yfxform(ydatafit),yfxform(modelfunN(pbxform(phat))), sigmafit))));
%loglikelihoodN= @(phat)sum(((ydatafit-modelfunN(pbxform(phat))).^2)./(2.*sigmafit.^2) + log(sigmafit.^2) + 0.5*log(2*pi));
%loglikelihoodN = @(phat)0.5.*(sum(log(normpdf((ydatafit),(modelfunN(pbxform(phat))), sigmafit))));
NLLN = @(phat)-loglikelihoodN_t(phat)
negLLNtrue = NLLN(pfxform(pfset))
negLLNguess = NLLN(pfxform(theta))

loglikelihoodphi =@(phat)0.5.*(sum(log(normpdf(yfxformp(phitrt),yfxformp(modelfunphi(pbxform(phat))), phisigfit))));
loglikelihood = @(phat)loglikelihoodphi(phat);

%% 
%J = @(phat) (sum(((mudatavec-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    %sum(((vardatavec-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));

   
    % Write objective functions for each model
    objfun = @(phat)-loglikelihood(phat); 
    phatbest = fminsearch(objfun, pfxform(theta));
   
    pbest = pbxform(phatbest);
%%
    model_N = modelfunN(pbest2);
    model_phi = modelfunphi(pbest2);
   negLL = objfun(phatbest2);
   %%
   figure;
   subplot(1,2,1)
   plot(ytimefit, model_N, '.')
   hold on
   plot(ytimefit, ydatafit, '.')
   subplot(1,2,2)
   plot(tbot, model_phi, 'LineWidth', 2)
   hold on
   plot(tbot, phitrt,'.', 'LineWidth', 2)
   %% Try gradient descent method for finding phi
   phisigfit= 0.1*ones(length(tbot),1);
   
    tol = 1e-5;
    lambda =2;
    it_count = 0;
    negLL_c = 0;
%     % want to minimize your objective function, which in some cases is the
%     % sum of your errors, and here is the negative loglikelihood.
    max_iters = 2000;
    deltap = 1e-8;
    spant = length(tbot);
%     %Initialize best guess
    negLLbest = objfun(pfxform(theta));
    negLLinit = negLLbest;
    phiinit = modelfunphi(theta);
    pe = theta;
    pevec = [];
    delvec = [];
    weightvec = 1./(phisigfit.^2);
    W =diag(weightvec);
    eta = 1e-2; 
    
    while negLLbest>tol && it_count<max_iters
        it_count = it_count +1;
        
        % Calculate Gradient of your cost function (negLL)
        for n = 1:nfit
            ptemp = pe;
            ptemp(n) = ptemp(n) + deltap;
            negLLt = objfun(pfxform(ptemp)); % find temporary model output
            dobjdp= (negLLt-negLLbest)./deltap;
            % Calculate the change in your parameters
            del(n) = eta.*(negLLt-negLLbest)./deltap;
  
        end
            delvec = vertcat(delvec, del);
            
            % Update parameters by adding del
            % reset phattemp to pe (back to your current parameters that
            % were probed while calculating the Jacobian
            ptemp = pe;
            pevec= vertcat(pevec,pe);
            % transform into reals to impose real bounds 
            for n = 1:nfit
                ptemp(n)= ptemp(n)-del(n);
                if ptemp(n)< pbounds(n,1)
                    ptemp(n) = pbounds(n,1);
                elseif ptemp(n)>pbounds(n,2)
                    ptemp(n) = pbounds(n,2);
                end
            end
            
            %Now have an updated phattemp
            % Evaluate the model
            negLLt = objfun(pfxform(ptemp));
            if negLLt < negLLbest
                negLL_c(it_count) = negLLt;
                negLLbest = negLLt;
                pe = ptemp;
                lambda = lambda/2;
                eta = 0.99*eta;
            else
                negLL_c(it_count) = negLLt;
                lambda = lambda*4;
                eta = 1.01*eta;
            end
    end
    
    pbestGD = pe;
    model_hiGD = modelfunphi(pbestGD);
    negLLGD = negLLbest;
%end