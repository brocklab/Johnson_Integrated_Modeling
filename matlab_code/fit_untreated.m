function [punt, singexpmodel] = fit_untreated(ydataf,ytimef, sigma)
% FIT UNTREATED CONTROL

%The simplified  model looks like this
%
% $$ N(t) = N0*K*(1/(N0 + (K-N0)exp(-gst));
% 
% where:
%
% * $N_0$ is the initial cell number
% * gs is the sensitive cell growth rate 
% * K is the carrying capacity


% Define transforms 
% single exponential
pfxform1 = @(pval)[1 1].*log(pval); %'forward' parameter transform into Reals
pbxform1 = @(phat)[1 1].*exp(phat);  %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output


N0 = ydataf(1);
ydata = ydataf(1:end);

ytime = ytimef(1:end);

% Set up forward models, fit all three nested versions of model
modelfun1 = @(p)simmodel1(p, ytime, N0); % single exponential model with carrying capacity  


    % INITIAL GUESSES BASED ON DATA

    gguess = (yfxform(ydata(end))-yfxform(ydata(end-5)))/(ytime(end)-ytime(end-5)); 
    % alter initial guesses to prevent NaNs and zeros
    Kguess = ydata(end);
    if isnan(gguess)
        gguess = 1e-5;
    end
    if isinf(gguess)
        gguess = 0.8;
    end
    if gguess <= 0
        gguess = 1e-5;
    end
    
    % Initial guess matrices
    theta1 = [gguess, Kguess]; % and k

    
    % Write log likelihood function based on assumption of normally
    % distributed sampling error
    
    % Goal: maximize the probability of the data given the model. Normpdf
    % will output a probability of the data (x- 1st argument), given the
    % mean(expectation, here the model), and the variance, at each time point. 
    % take the log and minimize the NLL to maximize likeihood of data given
    % the model
    
    % single exponential with carrying capacity
    loglikelihood1 = @(phat)sum(log(normpdf(yfxform(ydata),yfxform(modelfun1(pbxform1(phat))), sigma)));
    % single exponential with growth
   
    % Write objective functions for each model
    objfun1 = @(phat)-loglikelihood1(phat); 
    phatbest1 = fminsearch(objfun1, pfxform1(theta1));
    
    pi = pbxform1(phatbest1);
    gs = pi(1);
    carcap = pi(2); 
    punt = [gs, carcap];
    singexpmodel = simmodel1(pbxform1(phatbest1), ytime, N0);
end

