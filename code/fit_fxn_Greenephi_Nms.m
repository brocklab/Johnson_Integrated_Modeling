function [pbestvec,model_N, model_phi, negLLvec, werr_N,werr_phi] = fit_fxn_Greenephi_Nms(ydatafit, sigmafit,phitrt, phisigfit, pfID, psetID, pfit, thetarange, pset, time, tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s, N0phi, lambda, pbounds) 
%[pbest,model_N, negLL, pbestGD, model_NGD,                      (Ntrt,sigmafit,phitrt, phisigfit, pfitID, psetID, theta, pset, ytimefit,tbot, Uvec, Ub, lengthvec,lengthvecphi, N0s,N0phi,lambda, pbounds);
% This function needs to redoes the fit by searching for the optimal zr and
% zd which are the ratios of rr/rs and dr/ds respectively
% Therefore it will output a pbest in this form, and then we will convert
% it to the parameter estimates in a next step.

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
carcap1 = 0;
carcap2 = 0;
alpha = 0;
zr = 0;
ds = 0;
zd = 0;
params = [prop, rs, carcap1, carcap2, alpha, zr, ds, zd];
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
[prop, rs, carcap1, carcap2, alpha, zr, ds, zd] = deal(P{:});
 
% Define transforms: Now these must be logs and logits since zr and zd are
% bound between 0 and 1 (3&5)
% for number of variables in pfit
nfit = length(pfit);

% This needs to be generalizable based on the psetID and pfitID
% or we could just add an additional input to the function that is ikeep,
% and ikeep will act to reduce these vectors of 1s and 0s
pfxform = @(pval)log(pval);%+vec2curr.*log(pval./(1-pval)); %'forward' parameter transform into Reals
pbxform = @(phat)exp(phat);%+vec2curr.*(exp(phat)./(1+exp(phat)));  %'backward' parameter transform into model space

yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% Write a loop that will randomly select initial guesses from the bounds of
% thetarange.
Nstarts = 100;
Xmat = rand(length(pfit), Nstarts);
for i = 1:Nstarts
% Make a matrix where each column is a initial guess
thetamat(:,i) = thetarange(:,1) + Xmat(:,i).*(thetarange(:,2)-thetarange(:,1));
end

nsampsN= length(ydatafit);
nsampsphi = length(phitrt);

    % write a function that takes set parameters and fit parameters, combines
    % them, and runs the model forward for the Uvec and 0s provided
    modelfunN = @(p)simmodelgreene2(p, time, N0s, pset, Uvec, lengthvec, pfID, psetID);
    modelfunphi = @ (p)simmodelgreenephi2(p, tbot, N0phi, pset, Ub, lengthvecphi, pfID, psetID);
    %
    
    loglikelihoodphi= @(phat)sum(log(2*pi.*(phisigfit).^2)+((modelfunphi(pbxform(phat))-phitrt)./phisigfit).^2);
    loglikelihoodN= @(phat)sum(log(2*pi.*(sigmafit).^2)+((modelfunN(pbxform(phat))-ydatafit)./sigmafit).^2);
    %What if the objective function was just minimized error weighted by uncertainty??
    weightederrN =@(phat)sum(((modelfunN(pbxform(phat))-ydatafit).^2)./(sigmafit.^2));
    weightederrphi =@(phat)sum(((modelfunphi(pbxform(phat))-phitrt).^2)./(phisigfit.^2));
    % What if the objective function was just the minimized error, with no
    % weighting?
    %weightederrN =@(phat)sum(((modelfunN(pbxform(phat))-ydatafit).^2));
    %weightederrphi =@(phat)sum(((modelfunphi(pbxform(phat))-phitrt).^2));
    % Make your objective function weighted by the number of samples in
    % each data set
    objfun= @(phat)(1-lambda).*(1./nsampsN).*weightederrN(phat) + lambda.*(1./nsampsphi).*weightederrphi(phat);
    
    %objfunLL=@(phat)((1-lambda)*loglikelihoodN(phat) + lambda*loglikelihoodphi(phat));
    LB = pfxform(pbounds(:,1)');
    UB = pfxform(pbounds(:,2)');
    % This should make the search more exhaustive...
    options = optimset('MaxIter', 1e4, 'MaxFunEvals', 1e6, 'TolFun', 1e-6, 'TolX', 1e-6);
%INITALIZE with current theta
    phatbest = fminsearchbnd(objfun, pfxform(pfit), LB, UB, options);
    negLLvec(1) = objfun(phatbest);
    pbestvec(1,:) = pbxform(phatbest);   
    
for i = 2:Nstarts
    theta = thetamat(:,i); 
    
    phatbest = fminsearchbnd(objfun, pfxform(theta), LB, UB, options);
    negLLvec(i) = objfun(phatbest);
    pbestvec(i,:) = pbxform(phatbest);
end
    [negLLbest, ind] = min(negLLvec);
    pbest = pbestvec(ind, :);



   % model_N = modelfunN(pbest);
   % model_phi = modelfunphi(pbest);
   Nfwd = modelfunN(pbest);
    phifwd = modelfunphi(pbest);
    model_N = Nfwd;
    model_phi = phifwd;
    
    werr_N = sum(((Nfwd-ydatafit).^2)./(sigmafit.^2));
    werr_phi=sum(((phifwd-phitrt).^2)./(phisigfit.^2));
    
   
  %werr_N = sum(((model_N-ydatafit).^2)./sigmafit);
  %werr_phi = sum(((model_phi-phitrt).^2)./phisigfit);
   
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