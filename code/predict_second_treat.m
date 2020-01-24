% This script loads in the MCF-7 single pulse treatment parameter estimates
% and makes a prediction for a second dose treatment that we actually have
% data on...

% We start by simulating completely the effect of a second treatment using
% the forward model and the fitted parameters. We choose a dosing schedule
% that has been observed in our data.

close all; clear all; clc
%% Load in data structure and parameters
S = load('../out/trajfit.mat');
traj= S.traj;

S = load('../out/trajsumfit.mat');
trajsum = S.trajsum;

p4fit = load('../out/p4fit.mat');
p4fit = struct2cell(p4fit);
p4fit = cell2mat(p4fit);
P = num2cell(p4fit);
[rs, carcap, alpha, rr, ds, dr] = deal(P{:});
%% Run this forward for a single pulse treatment at 75 nM 
% Again, assume R0=0; dr = 0 and
kdrug = 0.0175;
k = 0.5;
tdrug =1; % stands for first treatment
dt = 1;


% Pull from trajsum dose = 75
for i = 1:length(trajsum)
    if trajsum(i).Cdox == 75
        tvec = trajsum(i).tvec;
        Nmean = trajsum(i).Nmean;
        Nstd = trajsum(i).Nstd;
        U1 = trajsum(i).U;
        S0 = trajsum(i).Nmean(1);
    end

end
tvecdat = tvec;
R0 = 0;
p = [ S0, R0, rs, carcap, alpha, rr, ds, dr];
tlong = 0:1:1344;
tvec = tlong; % simulate long term treatment dynamics
Cdox = 75;
U1 = k*Cdox*exp(-kdrug*(tlong));
[Nsri, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U1, dt, tdrug);
fracSens = Nsri(end,2)./Nsri(end,1);
sensfrac_t = Nsri(:,2)./ Nsri(:,1);
resfrac_t = Nsri(:,3)./Nsri(:,1);
figure;
plot(tvecdat, Nmean,'k*', 'LineWidth', 2)
hold on
plot(tvec, Nsri(:,1), 'LineWidth',2)
plot(tvecdat, Nmean + Nstd, 'color', 'k')
plot(tvecdat, Nmean - Nstd, 'color', 'k')
xlabel ('time (hours)')
ylabel ('N(t)')
title('N(t) for single 75 nM pulse treatment')
legend ('data mean', 'model', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(tvec, Nsri(:,1),'b', 'LineWidth',2)
hold on
plot(tvec, Nsri(:,2), 'g','LineWidth',2)
plot(tvec, Nsri(:,3),'r', 'LineWidth',2)
text(tvec(2), Nsri(2,2),['\phi_{sens_i}=', num2str(1)])
text(tvec(end-20), Nsri(end-20,2),['\phi_{sens_f}=', num2str(fracSens)])
xlabel ('time (hours)')
ylabel ('N(t)')
title('S(t) & R(t)for single 75 nM pulse treatment')
legend ('N(t)', 'S(t)','R(t)', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(tvec, sensfrac_t,'g', 'LineWidth',3)
hold on
plot(tvec, resfrac_t, 'r','LineWidth',3)
text(tvec(end-100), sensfrac_t(end-100), ['\phi_{sens_t=8WPT}=', num2str(fracSens)], 'FontSize', 14)
text(tvec(end-100), resfrac_t(end-100), ['\phi_{res_t=8WPT}=', num2str(1-fracSens)], 'FontSize', 14)
plot(tvec(end-100), sensfrac_t(end-100), 'k*', 'LineWidth',5)
plot(tvec(end-100), resfrac_t(end-100), 'k*', 'LineWidth',5)
text(tvec(30), sensfrac_t(30), ['\phi_{sens_t=30hrs}=', num2str(sensfrac_t(30))], 'FontSize', 14)
text(tvec(30), resfrac_t(30), ['\phi_{res_t=30hrs}=', num2str(resfrac_t(30))], 'FontSize', 14)
plot(tvec(30), sensfrac_t(30), 'k*', 'LineWidth',5)
plot(tvec(504), resfrac_t(504), 'k*', 'LineWidth',5)
text(tvec(504), sensfrac_t(504), ['\phi_{sens_t=3WPT}=', num2str(sensfrac_t(504))], 'FontSize', 14)
text(tvec(504), resfrac_t(504), ['\phi_{res_t=3WPT}=', num2str(resfrac_t(504))], 'FontSize', 14)
plot(tvec(504), sensfrac_t(504), 'k*', 'LineWidth',5)
plot(tvec(504), resfrac_t(504), 'k*', 'LineWidth',5)
xlim([ 0 tlong(end)])
ylabel ('Proportion of cells')
title('\phi_{S} & \phi_{R} following single 75 nM pulse treatment')
% legend ( '\phi_{S}','\phi_{R}', 'Location', 'NorthWest')
% legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)

%% Test function to create an easy to work with combined data structure
filter_criteria = 'treatmentnum';
ntreat = 2;
dose = 75;
date = '12-19-18';
[traj2] = comb_data_fxn(traj, filter_criteria, ntreat, dose, date)

%% Write a function to find phi at a certain time
% If repeat dose is within the time frame measured in the first pusle
% treatment

  Cdox = 75;
for i = 1:length(traj2)
   % Find the Number of sensitive and resisitant cells predicted at start
   % of second treatment
    tdose = traj2(i).WPT*24*7;
    tfirst = 0:1:tdose;
    U1 = k*Cdox*exp(-kdrug*(tfirst));
    p1 = [ S0, R0, rs, carcap, alpha, rr, ds, dr];
    [Nsri, tcrit, Ncrit] = fwd_Greene_model(p1, tfirst, U1, dt, tdrug);
    ind = find(ismember(tfirst, tdose), 1, 'first');
    fracS= Nsri(ind,2)/Nsri(ind,1);
    N01(1,i) = traj2(i).Nmean(1);
    S01(1,i) = fracS*N01(1,i);
    R01(1,i) = (1-fracS)*N01(1,i);

% Simulate second dose 
    tfin = 880;
    tvec2 = 0:dt:tfin;
    U2 = k*Cdox*exp(-kdrug*(tvec2));
    tdrug = 1; % only simulating one dose)
    p2 = [S01(1,i), R01(1,i), p(3:end)];
    [Nsr2, tcrit2, Ncrit2] = fwd_Greene_model(p2, tvec2, U2, dt, tdrug);
    traj2(i).Nmodpred = Nsr2;
    traj2(i).tvecmod = tvec2;




subplot(2, length(traj2)/2, i)

    plot(traj2(i).tvec, traj2(i).Nmean, 'k*')
    hold on
    plot (traj2(i).tvecmod, traj2(i).Nmodpred(:,1),'b-' ,'LineWidth',2)
    plot(traj2(i).tvec, traj2(i).Nmean + traj2(i).Nstd, 'color', 'k')
    plot(traj2(i).tvec, traj2(i).Nmean -traj2(i).Nstd, 'color', 'k')
    title (['WPT=', num2str(traj2(i).WPT)])
    xlabel ('time (hours)')
    xlim([0 traj2(i).tvec(end)])
    ylabel ('N(t)')
    title(['WPT=', num2str(traj2(i).WPT)])
    legend ('data mean', 'model', 'Location', 'NorthWest')
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
    
end
        
