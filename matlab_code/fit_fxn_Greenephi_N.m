function [pbest,model_N, model_phi, negLL] = fit_fxn_Greenephi_N(ydatafit, sigmafit,phitrt, phisigfit, pfID, psetID, pfit, pset, time, tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s, N0phi, lambda, pbounds) 
%[pbest,model_N, negLL, pbestGD, model_NGD,                      (Ntrt,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);
% This function

% Delete this section once remake into a function


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
%

loglikelihoodphi= @(phat)sum(log(2*pi.*(phisigfit).^2)+((modelfunphi(pbxform(phat))-phitrt)./phisigfit).^2);
loglikelihoodN= @(phat)sum(log(2*pi.*(sigmafit).^2)+((modelfunN(pbxform(phat))-ydatafit)./sigmafit).^2);
objfun=@(phat)((1-lambda)*loglikelihoodN(phat) + lambda*loglikelihoodphi(phat));
LB = pfxform(pbounds(:,1)');
UB = pfxform(pbounds(:,2)');
phatbest = fminsearchbnd(objfun, pfxform(theta), LB, UB);
pbest = pbxform(phatbest);



    model_N = modelfunN(pbest);
    model_phi = modelfunphi(pbest);
   negLL = objfun(phatbest);
   
   % Try gradient descent method for finding phi
   
%     tol = 1e-5;
%     lambda =2;
%     it_count = 0;
%     negLL_c = 0;
%     % want to minimize your objective function, which in some cases is the
%     % sum of your errors, and here is the negative loglikelihood.
%     max_iters = 2000;
%     deltap = 1e-8;
%     spant = length(tbot);
%      %Initialize best guess
%     negLLbest = objfun(pfxform(theta));
%     negLLinit = negLLbest;
%     phiinit = modelfunphi(theta);
%     pe = theta;
%     pevec = [];
%     delvec = [];
%     weightvec = 1./(phisigfit.^2);
%     W =diag(weightvec);
%     eta = 1e-2; 
    
%     while negLLbest>tol && it_count<max_iters
%         it_count = it_count +1;
%         
%         % Calculate Gradient of your cost function (negLL)
%         for n = 1:nfit
%             ptemp = pe;
%             ptemp(n) = ptemp(n) + deltap;
%             negLLt = objfun(pfxform(ptemp)); % find temporary model output
%             dobjdp= (negLLt-negLLbest)./deltap;
%             % Calculate the change in your parameters
%             del(n) = eta.*(negLLt-negLLbest)./deltap;
%   
%         end
%             delvec = vertcat(delvec, del);
%             
%             % Update parameters by adding del
%             % reset phattemp to pe (back to your current parameters that
%             % were probed while calculating the Jacobian
%             ptemp = pe;
%             pevec= vertcat(pevec,pe);
%             % transform into reals to impose real bounds 
%             for n = 1:nfit
%                 ptemp(n)= ptemp(n)-del(n);
%                 if ptemp(n)< pbounds(n,1)
%                     ptemp(n) = pbounds(n,1);
%                 elseif ptemp(n)>pbounds(n,2)
%                     ptemp(n) = pbounds(n,2);
%                 end
%             end
%             
%             %Now have an updated phattemp
%             % Evaluate the model
%             negLLt = objfun(pfxform(ptemp));
%             if negLLt < negLLbest
%                 negLL_c(it_count) = negLLt;
%                 negLLbest = negLLt;
%                 pe = ptemp;
%                 lambda = lambda/2;
%                 eta = 0.99*eta;
%             else
%                 negLL_c(it_count) = negLLt;
%                 lambda = lambda*4;
%                 eta = 1.01*eta;
%             end
%     end
%     
%     pbestGD = pe;
%     model_hiGD = modelfunphi(pbestGD);
%     negLLGD = negLLbest;
end