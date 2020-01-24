%[pbest,model_N] = fit_fxn_Greene(ydatafit, sigmafit, pfID, psetID, pfit, pset, time, Uvec, lengthvec, N0s)
%ydatafit,sigmafit, pfID, psetID, pfitguess, pset, ytimefit, Uvec, lengthvec, N0s
pfit = pfitguess;
time = ytimefit;

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
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output


theta = pfitguess; 

% write a function that takes set parameters and fit parameters, combines
% them, and runs the model forward for the Uvec and 0s provided
modelfun = @(ptest)simmodelgreene(ptest, time, N0s, pset, Uvec, lengthvec, pfID, psetID); 

loglikelihood = @(phat)sum(log(normpdf(yfxform(ydatafit),yfxform(modelfun(pbxform(phat))), sigmafit)));
  
   
    % Write objective functions for each model
    objfun = @(phat)-loglikelihood(phat); 
    phatbest = fminsearch(objfun, pfxform(theta));
    
    pbest = pbxform(phatbest);

    model_N = modelfun(pbest);
    negLLbest = objfun(phatbest)
    negLLinit = objfun(pfxform(theta))
    dif = negLLinit-negLLbest
    
    
    %% Consider instead using Levenberg Marquadt algorithm to maximize
    % likelihood
    
    tol = 1e-5;
    lambda =2;
    it_count = 0;
    negLL_c = 0;
    % want to minimize your objective function, which in some cases is the
    % sum of your errors, and here is the negative loglikelihood.
    max_iters = 2000;
    deltap = 1e-8;
    spant = length(time);
    %Initialize best guess
    negLLbest = objfun(pfxform(theta));
    negLLinit = negLLbest;
    Ntinit = modelfun(theta);
    pe = theta;
    pbounds = [0, 1; 0, 1; 0,1]; % set your parameter bounds pretty widely
    pevec = [];
    delvec = [];
    weightvec = 1./(sigmafit.^2);
    W =diag(weightvec);
    eta = 1e-2;
    %% 
    
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
 %% Plot some things to demonstrate algorithm convergence
    
pbestLM = pe;
Nbest = modelfun(pbestLM);
chi_sq = sum(((Nbest-ydatafit)./sigmafit).^2)
 figure;
 plot(time, Nbest, 'b*', 'LineWidth',1)
 hold on
 plot(time, ydatafit, 'k*', 'LineWidth', 2)
 plot(time, ydatafit+1.96*sigmafit, 'k.', 'LineWidth',1)
  plot(time, ydatafit-1.96*sigmafit, 'k.', 'LineWidth',1)
 xlabel('time (hours)')
 ylabel('N(t)')
 legend('best N(t)', 'real data', '95% CI in data')
 legend boxoff
 title(['Data versus model for gradient descent, \chi^{2}=', num2str(chi_sq)])
 
 figure;
 plot(1:1:max_iters, negLL_c, 'b*')
 hold on
 plot(1, negLLinit, 'r*')
 legend('current negLL', 'initial negLL')
 xlabel('iteration')
 ylabel('objective function (negative log likelihood)')
 title('Convergence of objective function')
 


 
 
 