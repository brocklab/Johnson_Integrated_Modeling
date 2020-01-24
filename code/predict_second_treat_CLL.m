% This script loads in the CLL single pulse treatment parameter estimates
% and makes a prediction for a second dose treatment that we actually have
% data on...

% We start by simulating completely the effect of a second treatment using
% the forward model and the fitted parameters. We choose a dosing schedule
% that has been observed in our data.

close all; clear all; clc
%% Load in data structure and parameters
S = load('../out/CLLdataresp.mat');
CLL= S.CLLdata;
tdata = CLL(1).time;
Ndata = CLL(1).rawN;
N0 = Ndata(1);

NdataFM1 = CLL(4).rawN;
NdataFM7 = CLL(7).rawN;

% load in the scRNAseq estimates of phi(t) and the corresponding time
% vectors from python

phi_est_filename = '../data/phi_t_CLL.csv';
phi_est = readtable(phi_est_filename);
tbot = phi_est.t; % this is 29 days 
phitrt = phi_est.phi_t;
ntrt = phi_est.ncells;
sigtech = 1e-2;
phisigfit = [phitrt.*(1-phitrt)./ntrt] + sigtech;
N0phi = 1e7; % set this because this is what we think we seeded for this experiment
phi0= phitrt(1);
S0phi = phi0*N0phi;
R0phi = (1-phi0)*N0phi;
Cfluphi = 5; % I have no idea what this trea
Cflumax = 50;

toutgrowth = 29*24;
tlong = 0:4:toutgrowth;
tvec = [0:1:toutgrowth];


ptest = load('../out/pCLL.mat');
ptest = struct2cell(ptest);
ptest = cell2mat(ptest);
P = num2cell(ptest);

[phi0, carcapN, carcapphi, rs, alpha, zr, ds, zd, k, kdrug] = deal(P{:});
Unt =k*Cfluphi*exp(-kdrug*(tvec))/(0.1*Cflumax);
%% Run this forward for single pulse treatment of 5uM flu in 12 well

p = [ phi0, carcapN,rs,alpha, zr, ds, zd];
toutgrowth = 29*24;
tlong = 0:4:toutgrowth;
dt = 1;
tdrug = 1;
Cfluphi = 5; % I have no idea what this trea
Cflumax = 50;
tgen = 0:1:tlong(end);

 % simulate long term treatment dynamics
[Nsri, tcrit, Ncrit] = fwd_Greene_model2(p, tlong, N0, Unt, dt, tdrug);
fracSens = Nsri(end,2)./Nsri(end,1);
sensfrac_t = Nsri(:,2)./ Nsri(:,1);
resfrac_t = Nsri(:,3)./Nsri(:,1);
figure;
plot(tdata, Ndata,'k*', 'LineWidth', 2)
hold on
plot(tlong, Nsri(:,1), 'b-', 'LineWidth',3)
xlabel ('time (hours)')
ylabel ('N(t)')
title('N(t) for single 5 uM pulse treatment')
legend ('data', 'model', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0 toutgrowth])
ylim([ 0 Ndata(end)+0.5e6])

figure;
plot(tdata, Ndata,'k*', 'LineWidth', 2)
hold on
plot(tlong, Nsri(:,1),'b', 'LineWidth',2)
hold on
plot(tlong, Nsri(:,2), 'g','LineWidth',2)
plot(tlong, Nsri(:,3),'r', 'LineWidth',2)
text(tlong(2), Nsri(2,2),['\phi_{sens_i}=', num2str(1)], 'FontSize', 14)
text(tlong(end-20), Nsri(end-20,2),['\phi_{sens_f}=', num2str(fracSens)], 'FontSize', 14)
xlim([ 0 toutgrowth])
xlabel ('time (hours)')
ylabel ('N(t)')
title('S(t) & R(t)for single 5uM pulse treatment')
legend ('N(t) data', 'N(t)', 'S(t)','R(t)', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
ylim([ 0 Ndata(end)+0.5e6])

figure;
plot(tlong, sensfrac_t,'g', 'LineWidth',3)
hold on
plot(tlong, resfrac_t, 'r','LineWidth',3)
text(tlong(end-20), sensfrac_t(end-10)-.1, ['\phi_{sens_t=29d}=', num2str(fracSens)], 'FontSize', 14)
text(tlong(end-20), resfrac_t(end-10)+.1, ['\phi_{res_t=29d}=', num2str(1-fracSens)], 'FontSize', 14)
plot(tlong(end-10), sensfrac_t(end-10), 'k*', 'LineWidth',5)
plot(tlong(end-10), resfrac_t(end-10), 'k*', 'LineWidth',5)
xlim([ 0 toutgrowth])
ylabel ('Proportion of cells')
title('\phi_{S} & \phi_{R} following single 5uM pulse treatment')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Now run forward for the scRNAseq experiment
 % simulate long term treatment dynamics
 p(2) = carcapphi;

[Nsri, tcrit, Ncrit] = fwd_Greene_model2(p, tlong, N0phi, Unt, dt, tdrug);
fracSens = Nsri(end,2)./Nsri(end,1);
sensfrac_t = Nsri(:,2)./ Nsri(:,1);
resfrac_t = Nsri(:,3)./Nsri(:,1);
figure;
plot(tdata, Ndata,'k*', 'LineWidth', 2)
hold on
plot(tlong, Nsri(:,1), 'b-', 'LineWidth',3)
xlabel ('time (hours)')
ylabel ('N(t)')
title('N(t) for single 5 uM pulse treatment')
legend ('data', 'model', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0 toutgrowth])
ylim([ 0 Ndata(end)+0.5e6])

figure;
hold on
plot(tlong, Nsri(:,1),'b', 'LineWidth',2)
hold on
plot(tlong, Nsri(:,2), 'g','LineWidth',2)
plot(tlong, Nsri(:,3),'r', 'LineWidth',2)
text(tlong(2), Nsri(2,2),['\phi_{sens_i}=', num2str(1)], 'FontSize', 14)
text(tlong(end-20), Nsri(end-20,2),['\phi_{sens_f}=', num2str(fracSens)], 'FontSize', 14)
xlim([ 0 toutgrowth])
xlabel ('time (hours)')
ylabel ('N(t)')
title('S(t) & R(t)for single 5uM pulse treatment from scRNAseq')
legend ( 'N(t)', 'S(t)','R(t)', 'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)


figure;
plot(tlong, sensfrac_t,'g', 'LineWidth',3)
hold on
plot(tlong, resfrac_t, 'r','LineWidth',3)
text(tlong(end-20), sensfrac_t(end-10)-.1, ['\phi_{sens_t=29d}=', num2str(fracSens)], 'FontSize', 14)
text(tlong(end-20), resfrac_t(end-10)+.1, ['\phi_{res_t=29d}=', num2str(1-fracSens)], 'FontSize', 14)
plot(tlong(end-10), sensfrac_t(end-10), 'k*', 'LineWidth',5)
plot(tlong(end-10), resfrac_t(end-10), 'k*', 'LineWidth',5)
xlim([ 0 toutgrowth])
ylabel ('Proportion of cells')
title('\phi_{S} & \phi_{R} following single 5uM pulse treatment from scRNAseq')
set(gca,'FontSize',20,'LineWidth',1.5)


%% Simulate repeat treatment and then compare it to data
tdose = toutgrowth;
tfirst = 0:1:tdose;
tsecond = 0:1:480; % monitor for two weeks after
tdrug = [0; tdose];
ttot = horzcat(tfirst, tsecond(2:end) + tfirst(end));
U1 =k*Cfluphi*exp(-kdrug*(tfirst))/(0.1*Cflumax);
U2 =k*Cfluphi*exp(-kdrug*(tsecond))/(0.1*Cflumax);
U = horzcat(U1(1:end-1), U2);
    % Run the model forward for the two different U(t)s and times of drug
    % for each variest treatment interval
p2 = p;
p2(1) = sensfrac_t(end);
N02 = 1e6;
p2(2) = carcapN; % these experiments now done in a 12 well
[Nsri, tcrit, Ncrit] = fwd_Greene_model2(p2, tsecond, N02, U2, dt, tdrug(1));

figure;
 hold on
plot (tsecond, Nsri(:,1),'b-' ,'LineWidth',2)
plot(tsecond, Nsri(:,2), 'g-', 'LineWidth',2)
plot(tsecond, Nsri(:,3), 'r-', 'LineWidth', 2)
plot(tdata, NdataFM1,'*', 'color', CLL(4).color, 'LineWidth', 2)
plot(tdata, NdataFM7, '*', 'color', CLL(7).color, 'LineWidth',2)
plot(tdata, Ndata, '*', 'color', CLL(1).color, 'LineWidth',2)
title (['Prediction of repeat treatment at 29 days'])
xlabel ('time (hours)')
ylabel ('N(t)')
legend ('N(t)','S(t)', 'R(t)','FM1', 'FM7','TP0',  'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0 tsecond(end)])

