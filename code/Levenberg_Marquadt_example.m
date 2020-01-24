%% Problem 3a (1D optimization)
clear all; clc;
load N_1d-1
ktrue = k;
Dtrue = 0.01;
carcap = 1;
dx = 0.1; 
dt = 0.02;
sx = 100;
x = 0:.1:10;
num_ps = length(ktrue) + length(Dtrue);
pbounds = zeros(num_ps, 2);
pbounds(1,1) = 0;
pbounds(1,2) = 5;
pbounds (2:num_ps, 1) = -2;
pbounds(2:num_ps,2) = 8;
ptrue = vertcat(Dtrue, ktrue);    
% Noiseless
N_s_inf = N_s;
Nm = zeros(100,3); % looking at time points 1, 2, and 3 for now 
Nm(:,1) = N_s(:,1);
% SNR 16
N_s_snr_16 = N_s.*random('Normal',1,1/16,[100 10]);

% SNR 4
N_s_snr_4 = N_s.*random('Normal',1,1/4,[100 10]);
%% Check FDM

% Make sure FDM model works by comparing ktrue and Dtrue approximation to
% measurement with no noise
Ntrue=N_s;
Ntest = Rxn_Diff_CD1D( N_s, 2, ktrue, Dtrue );

figure(1)
hold off
plot(1:100,Ntest(:,2./dt),'LineWidth',1.5);
hold on
plot(1:100,Ntrue(:,2),'LineWidth',1.5);

set(gca,'LineWidth',1.5,'FontSize',12);
xlabel('x','FontSize',20)
ylabel('Number of Cells','FontSize',20)
title('Problem 3a, Compare guess to measured no noise', 'FontSize',14)
legend('FDM Model N',' Noiseless Measured N')
  
%% Initialize Parameters and guesses

 % Parameters in model are k and D0, k is a 100 x 1, Do is single value
 % g stands for initial guess
 kg = ones(sx,1);
 Dg = .05;
 lambda = 2;
 tol = 1e-5;
 time = 4;
 % looking at time points 1, 2, and 3 for now 
 % Initialize guess for N with first time point from data
 %% First let Ntrue = N_s (noiseless data)
Ntrue =N_s;
Ng  = Rxn_Diff_CD1D( N_s, time, kg, Dg ); % this outputs 200 time steps of N, where the 100th and 200th are the 
                                          % measurements at time points 2
                                          % and 3 ( t = 2 days and t = 4
                                          % days)
Nbest = Ng; % initial best guess
% difference between model at t = 2 days and measurement at t = 2days
err = sum((Ng(:,2./dt)-Ntrue(:,2)).^2) + sum((Ng(:,4./dt)-Ntrue(:,4)).^2);% sum squared error between second and third time point
ke = kg; % best current values of k, let ke = that parameter guess
De = Dg; % best current value of D, let De = that parameter guess
pe = vertcat(De, ke);

% This is initializing step, now do LM method

%%
it_count = 0;
err_c = 0;
while err>tol && it_count<1000
    it_count = it_count +1;
    
% Calculate Jacobian
    for n = 1:num_ps
        pt = pe;
        pt(n) = pt(n) + 1e-8;
        kt = pt(2:num_ps);
        Dt = pt(1);
        Nt= Rxn_Diff_CD1D( N_s, time, kt, Dt );  
        J(1:sx,n) = (Nt(:,2/dt)-Nbest(:,2/dt))/(1e-8); % TP1
        J(sx+1:2*sx, n) = (Nt(:,(4/dt))-Nbest(:,(4/dt)))/(1e-8); % TP2
    end
 
  % Calculate change in parameters
    err_v(1:sx,1) = Ntrue(:,2)-Nbest(:,2/dt); % yhat - ymeasure % evaluate error
    err_v(sx+1:2*sx,1) = Ntrue(:,4)-Nbest(:,4/dt);
    del=(J'*J+lambda*diag(diag(J'*J)))\((J'*err_v));


    
    % Update parameters by adding del
    % reset pt to pe
     pt = pe; % reset pt to pe (current parameter values, probed during Jacobian calculation)
    for n = 1:num_ps
        pt(n) = pe(n) + del(n);
        if pt(n) < pbounds(n,1)
            pt(n) = pbounds(n,1);
        elseif pt(n) > pbounds(n,2)
            pt(n) = pbounds(n,2);
        end
    end
% Now have an updated pt
     
    % Evaluate the model
    Dt = pt(1);
    kt = pt(2:num_ps);
    Nt = Rxn_Diff_CD1D( N_s, time, kt, Dt ); 
    err_t = sum((Nt(:,2./dt)-Ntrue(:,2)).^2) + sum((Nt(:,4./dt)-Ntrue(:,4)).^2);% new error from change in parameters
    if err_t < err
        err_c(it_count) = err_t;
        err = err_t;
        pe = pt;
        lambda = lambda/2;
        Nbest = Nt;
    else
        err_c(it_count) = err_t;
        lambda = lambda*4;
    end
    
end
%%
% (3) Calculate parameter error
err_in_params_noiseless = pt-ptrue;
err_noiseless = sum(err_in_params_noiseless);
% (2) Predict future tumor growth (guess based on initial paramters)
% Use pt to predict tumor growth out to day 10
time_pred = 10;
N_pred  = Rxn_Diff_CD1D( N_s, time_pred, kt, Dt );
%%
subplot(1,2,1)
hold off
plot(1:100,N_pred(:,3./dt),'LineWidth',1.5);
hold on
plot(1:100,Ntrue(:,3),'LineWidth',1.5);
set(gca,'LineWidth',1.5,'FontSize',12);
xlabel('x','FontSize',20)
ylabel('Number of Cells','FontSize',20)
title('3 day prediction on noiseless data', 'FontSize',14)
legend('Best Model N',' Noiseless Measured N')

subplot(1,2,2)
hold off
plot(1:100,N_pred(:,6./dt),'LineWidth',1.5);
hold on
plot(1:100,Ntrue(:,6),'LineWidth',1.5);
set(gca,'LineWidth',1.5,'FontSize',12);
xlabel('x','FontSize',20)
ylabel('Number of Cells','FontSize',20)
title('6 day prediction on noiseless data', 'FontSize',14)
legend('Best Model N',' Noiseless Measured N')
%% %% Next let Ntrue = N_s_snr_16 
Ntrue =N_s_snr_16;
Ng  = Rxn_Diff_CD1D( N_s, time, kg, Dg ); % this outputs 200 time steps of N, where the 100th and 200th are the 
                                          % measurements at time points 2
                                          % and 3 ( t = 2 days and t = 4
                                          % days)
Nbest = Ng; % initial best guess
% difference between model at t = 2 days and measurement at t = 2days
err = sum((Ng(:,2./dt)-Ntrue(:,2)).^2) + sum((Ng(:,4./dt)-Ntrue(:,4)).^2);% sum squared error between second and third time point
ke = kg; % best current values of k, let ke = that parameter guess
De = Dg; % best current value of D, let De = that parameter guess
pe = vertcat(De, ke);

% This is initializing step, now do LM method

%%
it_count = 0;
err_c = 0;
while err>tol && it_count<1000
    it_count = it_count +1;
    
% Calculate Jacobian
    for n = 1:num_ps
        pt = pe;
        pt(n) = pt(n) + 1e-8;
        kt = pt(2:num_ps);
        Dt = pt(1);
        Nt= Rxn_Diff_CD1D( N_s, time, kt, Dt );  
        J(1:sx,n) = (Nt(:,2/dt)-Nbest(:,2/dt))/(1e-8);
        J(sx+1:2*sx, n) = (Nt(:,(4/dt))-Nbest(:,(4/dt)))/(1e-8);
    end
 
  % Calculate change in parameters
    err_v(1:sx,1) = Ntrue(:,2)-Nbest(:,2/dt); % yhat - ymeasure % evaluate error
    err_v(sx+1:2*sx,1) = Ntrue(:,4)-Nbest(:,4/dt);
    del=(J'*J+lambda*diag(diag(J'*J)))\((J'*err_v));


    
    % Update parameters by adding del
    % reset pt to pe
     pt = pe; % reset pt to pe (current parameter values, probed during Jacobian calculation)
    for n = 1:num_ps
        pt(n) = pe(n) + del(n);
        if pt(n) < pbounds(n,1)
            pt(n) = pbounds(n,1);
        elseif pt(n) > pbounds(n,2)
            pt(n) = pbounds(n,2);
        end
    end
% Now have an updated pt
     
    % Evaluate the model
    Dt = pt(1);
    kt = pt(2:num_ps);
    Nt = Rxn_Diff_CD1D( N_s, time, kt, Dt ); 
    err_t = sum((Nt(:,2./dt)-Ntrue(:,2)).^2) + sum((Nt(:,4./dt)-Ntrue(:,4)).^2);% new error from change in parameters
    if err_t < err
        err_c(it_count) = err_t;
        err = err_t;
        pe = pt;
        lambda = lambda/2;
        Nbest = Nt;
    else
        err_c(it_count) = err_t;
        lambda = lambda*4;
    end
    
end
%%
% (3) Calculate parameter error
err_in_params_snr_16 = pt-ptrue;
err_snr_16= sum(err_in_params_snr_16);
% (2) Predict future tumor growth (guess based on initial paramters)
% Use pt to predict tumor growth out to day 10
time_pred = 10;
N_pred_16  = Rxn_Diff_CD1D( N_s, time_pred, kt, Dt );
%%
subplot(1,2,1)
hold off
plot(1:100,N_pred_16(:,3./dt),'LineWidth',1.5);
hold on
plot(1:100,Ntrue(:,3),'LineWidth',1.5);
set(gca,'LineWidth',1.5,'FontSize',12);
xlabel('x','FontSize',20)
ylabel('Number of Cells','FontSize',20)
title('3 day prediction on SNR 16 data', 'FontSize',14)
legend('Best Model N',' Noiseless Measured N')

subplot(1,2,2)
hold off
plot(1:100,N_pred_16(:,6./dt),'LineWidth',1.5);
hold on
plot(1:100,Ntrue(:,6),'LineWidth',1.5);
set(gca,'LineWidth',1.5,'FontSize',12);
xlabel('x','FontSize',20)
ylabel('Number of Cells','FontSize',20)
title('6 day prediction on SNR 16 data', 'FontSize',14)
legend('Best Model N',' Noiseless Measured N')
%% %% Next let Ntrue = N_s_snr_4 
Ntrue =N_s_snr_4;
Ng  = Rxn_Diff_CD1D( N_s, time, kg, Dg ); % this outputs 200 time steps of N, where the 100th and 200th are the 
                                          % measurements at time points 2
                                          % and 3 ( t = 2 days and t = 4
                                          % days)
Nbest = Ng; % initial best guess
% difference between model at t = 2 days and measurement at t = 2days
err = sum((Ng(:,2./dt)-Ntrue(:,2)).^2) + sum((Ng(:,4./dt)-Ntrue(:,4)).^2);% sum squared error between second and third time point
ke = kg; % best current values of k, let ke = that parameter guess
De = Dg; % best current value of D, let De = that parameter guess
pe = vertcat(De, ke);

% This is initializing step, now do LM method

%%
it_count = 0;
err_c = 0;
while err>tol && it_count<1000
    it_count = it_count +1;
    
% Calculate Jacobian
    for n = 1:num_ps
        pt = pe;
        pt(n) = pt(n) + 1e-8;
        kt = pt(2:num_ps);
        Dt = pt(1);
        Nt= Rxn_Diff_CD1D( N_s, time, kt, Dt );  
        J(1:sx,n) = (Nt(:,2/dt)-Nbest(:,2/dt))/(1e-8);
        J(sx+1:2*sx, n) = (Nt(:,(4/dt))-Nbest(:,(4/dt)))/(1e-8);
    end
 
  % Calculate change in parameters
    err_v(1:sx,1) = Ntrue(:,2)-Nbest(:,2/dt); % yhat - ymeasure % evaluate error
    err_v(sx+1:2*sx,1) = Ntrue(:,4)-Nbest(:,4/dt);
    del=(J'*J+lambda*diag(diag(J'*J)))\((J'*err_v));


    
    % Update parameters by adding del
    % reset pt to pe
     pt = pe; % reset pt to pe (current parameter values, probed during Jacobian calculation)
    for n = 1:num_ps
        pt(n) = pe(n) + del(n);
        if pt(n) < pbounds(n,1)
            pt(n) = pbounds(n,1);
        elseif pt(n) > pbounds(n,2)
            pt(n) = pbounds(n,2);
        end
    end
% Now have an updated pt
     
    % Evaluate the model
    Dt = pt(1);
    kt = pt(2:num_ps);
    Nt = Rxn_Diff_CD1D( N_s, time, kt, Dt ); 
    err_t = sum((Nt(:,2./dt)-Ntrue(:,2)).^2) + sum((Nt(:,4./dt)-Ntrue(:,4)).^2);% new error from change in parameters
    if err_t < err
        err_c(it_count) = err_t;
        err = err_t;
        pe = pt;
        lambda = lambda/2;
        Nbest = Nt;
    else
        err_c(it_count) = err_t;
        lambda = lambda*4;
    end
    
end
%%
% (3) Calculate parameter error
err_in_params_snr_4 = pt-ptrue;
err_snr_4= sum(err_in_params_snr_4);

% (2) Predict future tumor growth (guess based on initial paramters)
% Use pt to predict tumor growth out to day 10
time_pred = 10;
N_pred_4  = Rxn_Diff_CD1D( N_s, time_pred, kt, Dt );

subplot(1,2,1)
hold off
plot(1:100,N_pred_4(:,3./dt),'LineWidth',1.5);
hold on
plot(1:100,Ntrue(:,3),'LineWidth',1.5);
set(gca,'LineWidth',1.5,'FontSize',12);
xlabel('x','FontSize',20)
ylabel('Number of Cells','FontSize',20)
title('3 day prediction on SNR 4 data', 'FontSize',14)
legend('Best Model N',' Noiseless Measured N')

subplot(1,2,2)
hold off
plot(1:100,N_pred_4(:,6./dt),'LineWidth',1.5);
hold on
plot(1:100,Ntrue(:,6),'LineWidth',1.5);
set(gca,'LineWidth',1.5,'FontSize',12);
xlabel('x','FontSize',20)
ylabel('Number of Cells','FontSize',20)
title('6 day prediction on SNR 4 data', 'FontSize',14)
legend('Best Model N',' Noiseless Measured N')

%% Problem 3b (Parameter Optimization Noiseless)
clear all; clc
load N_2d
warning off;
ktrue = k;
Dtrue = 0.05;
carcap = 100;
    dx = 0.1; % in mm
    dy = 0.1; % in mm
    sx = 10/dx;  %no. grid points
    sy = 10/dy; %no. grid points
    dt = 0.01; %in days
    dims(1) = dx;
    dims(2) = dy;
    
   LB = zeros(26,1);  % Lower Bounds
   UB = 10*ones(26,1); % Upper Bounds
   UB(end) = (1/4)/(1/dx^2 + 1/dy^2)/dt;
   params = ones(26,1); params(end) = 0.5*UB(end); % Initial Guess...
   options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
   % Start with Nmeas = noiseless data
   
   Nmeas = N_snr_Inf;

   
   % Will be using time points 1, 2, and 3 to compare to 2D reaction
   % diffusion model at time points 1, 2, and 3
   Ntest(:,:,1) = Nmeas(:,:,1);
   
   % Run FDM assuming we know the initial conditions to be true, just as a
   % test
   Dtest = Dtrue*ones(100,100);
   ktest = ktrue;
   for t = 1:4/dt-1
   Ntest(:,:,t+1) = rd_fdm_center_v1(Ntest(:,:,t),Dtest,ktest,carcap,dims,dt);
   end
   
   % Compares FDM to measured data with correct parameters
   % compare second time point of N and 2nd day of Ntest
   figure(3)
   subplot(1,2,1)
   imagesc(Nmeas(:,:,2))
   title('Nmeasured at 48 hours')
   subplot(1,2,2)
   imagesc(Ntest(:,:,2/dt))
   title('FDM prediction at 48 hours')
 % compare third t
   figure(4)
   subplot(1,2,1)
   imagesc(Nmeas(:,:,3))
   title('Nmeasured at 4 days')
   subplot(1,2,2)
   imagesc(Ntest(:,:,4/dt))
   title('FDM prediction at 4 days')
 %% Now use LSQnonlin to create vector of model-data from measured data and FDM model

% note that params is an initial guess of 25 k values and 1 D value. Need
% to change the 25 k values to be 100 x 100 proliferation points 

% Make a cell containing all of the Nmeasured data sets for each SNR ratio
 Nmeas_cell = {N_snr_Inf, N_snr_64, N_snr_32, N_snr_16, N_snr_8, N_snr_4, N_snr_1};

for i = 1:7
    Nmeas_i = Nmeas_cell{i};
    opt_params(:,i) = lsqnonlin(@fit_data, params, LB, UB, options, Nmeas_i, carcap, dims, dt);
    
kbest(:,i) = opt_params(1:25,i);
Dbest(i) = opt_params(26,i);
params_best(:,i) = vertcat(kbest(:,i),Dbest(:,i));
Dbestbig = Dbest(i)*ones(100,100);
kbestbig= prolif_func(kbest(:,i));

kbestmat(:,:,i) = kbestbig;
% Run model forward with optimized parameters out to 10 time points (20
% days)
N_best = zeros(100,100,2000);
N_best(:,:,1) = Nmeas_i(:,:,1);

for t = 1:20/dt-1
   N_best(:,:,t+1) = rd_fdm_center_v1(N_best(:,:,t),Dbestbig,kbestbig,carcap,dims,dt);
end

N_best_final(:,:,i)= N_best(:,:,2000);

%Calculate sum-squared error in prediction
N_for_err = N_best(:,:,1:200:2000);
err_tot_all_time = N_for_err - Nmeas;
sum_sq_err_all_time(i) =  sum(sum(sum((err_tot_all_time).^2)));

% Calculate sum squared error in params

sum_sq_err_params(i) = sum(sum((kbestmat(:,:,i)-ktrue).^2))+ ((Dbest(i)-Dtrue).^2);

end

%%
   figure(3)
   hold off
   subplot(1,2,1)
   imagesc(N_snr_Inf(:,:,10))
   title('Noiseless N at 10th time')
   subplot(1,2,2)
   imagesc(N_best_final(:,:,1))
   title('Noiseless Model Fit N at 10th time')
   
   figure(4)
   hold off
   subplot(1,2,1)
   imagesc(N_snr_4(:,:,10))
   title('SNR = 4 N at 10th time')
   
   subplot(1,2,2)
   imagesc(N_best_final(:,:,6))
   title('SNR = 4 N Model at 10th time')  