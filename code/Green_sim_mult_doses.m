% Run simulations for different dosing regimens using the Greene model
% Start by assuming:
% 1. all cells are sensitive at the start
% 2. resistant cells are invincible (dr = 0)

close all; clear all; clc
%%
ptest = load('../out/ptest.mat');
ptest = struct2cell(ptest);
ptest = cell2mat(ptest);
P = num2cell(ptest);
[phi0, carcapNf, carcapphi, rs, alpha, zr, ds, zd, k, kdrug, gtot] = deal(P{:});
% reset the parametes to be near these values but not them exactly
phi0 = 0.8;
carcapN = 5e4;
rs = 0.023;
alpha = 0.0;
zr = 0.2;
ds = 0.05;
zd = 0.1;
Cdoxmax = 1000;
% Set up the conditions to simulate with pulsed treatment 

dt = 1; % this corresponds to hours

t1 = 168; % 2 weeks

dt = 1;
tdrug = 1;
Cdox = 200;
tgen = 0:1:t1;
U1=k*Cdox*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
N0 = 3000;
 pi = [phi0, carcapN,rs, alpha, zr, ds,zd];
 %
tvec = horzcat(tgen(1:end-1), tgen(1:end-1) + t1, tgen (1:end-1)+ 2*t1, tgen + 3*t1);
tdrug = [tgen(1); (tgen(1) + tgen(end)); tgen(end)+tgen(end)];
% set original U as U1 to indicate first dose
U2=k*Cdox*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
U3=k*Cdox*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
U4=k*Cdox*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
U = horzcat(U1(1:end-1), U2(1:end-1), U3(1:end-1), U4);

[Nsrpulse, tcrit, Ncrit] = fwd_Greene_model2(pi, tvec, N0, U, dt, tdrug);

Ucons = mean(U)*ones(length(U));
[Nsrcons, tcrit, Ncrit] = fwd_Greene_model2(pi, tvec, N0, Ucons, dt, tdrug);
% Add to the trajsum data structure the model fit for that Cdox (whether it



figure;
subplot(2,1,2)
plot(tvec,Nsrpulse(:,1),'LineWidth',3, 'color','b');
hold on
plot(tvec,Nsrcons(:,1), 'LineWidth', 3, 'color', 'k')
%plot(tvec, Nsrulse(:,2), 'LineWidth', 3,'color', 'g')
%plot(tvec, Nsrpulse(:,3), 'LineWidth', 3, 'color', 'r')
%plot(tcrit, Ncrit, 'k*', 'LineWidth', 3)
plot([0 tvec(end)], [2*N0, 2*N0], 'r--', 'LineWidth', 2)
%legend('total cell number', 'sensitive', 'resistant', 'critical N', 'Location', 'NorthWest')
legend('pulsed', 'constant', 'Location', 'NorthWest')
legend boxoff
xlim([ 0 tvec(end)])
ylim([ 0 3*N0])
xlabel('time (hours)','FontSize',20)
ylabel('N(t)','FontSize',20)
title(['\alpha= ', num2str(alpha)])
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(2,1,1)
plot(tvec, U,'LineWidth',3, 'color', 'b')
hold on
plot(tvec, Ucons, 'LineWidth', 3, 'color', 'k')
xlim([ 0 tvec(end)])
ylim([ 0 1.3])
xlabel('time (hours)','FontSize',20)
ylabel('Effective dose u(t)','FontSize',20)
legend('pulsed', 'constant')
legend boxoff 
set(gca,'FontSize',20,'LineWidth',1.5)
%title('Pulse of 200 nM dox every week')

%% Maximum achievable critical time versus alpha 
pstar = pi;
alphavec = linspace(0,0.01, 2);
dosevec = linspace(0, 1000, 10);
tcritstar = [];
for i = 1:length(alphavec)
    pstar(4) = alphavec(i);
    for j = 1:length(dosevec)
        tgen = 0:1:5*t1;
        Uj=k*dosevec(j)*exp(-kdrug*(tgen))/(0.1*Cdoxmax);
        [Nsri, tcritstar(i,j), Ncrit] = fwd_Greene_model2(pstar,tgen, N0, Uj, dt, 0);
    end
end 

figure;
for i = 1:length(alphavec)
    plot(dosevec, tcritstar(i,:), '*-', 'LineWidth', 2)
    hold on
end
xlabel('[Dox]')
ylabel('Critical time')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
for j = 1:length(dosevec)
    plot(alphavec, tcritstar(:,j), '*-', 'LineWidth', 2)
    hold on
end
xlabel('\alpha')
ylabel('Critical time')
set(gca,'FontSize',20,'LineWidth',1.5)

%% Plot just the first two doses
figure;
subplot(1,2,1)
plot(tvec,Nsrpulse(:,1),'LineWidth',3, 'color','b');
hold on
%plot(tvec,Nsrcons(:,1), 'LineWidth', 3, 'color', 'k')
plot(tvec, Nsrpulse(:,2), 'LineWidth', 3,'color', 'g')
plot(tvec, Nsrpulse(:,3), 'LineWidth', 3, 'color', 'r')
%plot(tcrit, Ncrit, 'k*', 'LineWidth', 3)
%plot([0 tvec(end)], [2*N0, 2*N0], 'r--', 'LineWidth', 2)
%legend('total cell number', 'sensitive', 'resistant', 'critical N', 'Location', 'NorthWest')
legend('N(t)', 'S(t)', 'R(t)', 'Location', 'NorthWest')
legend boxoff
xlim([ 0 2*t1])
%ylim([ 0 3*N0])
xlabel('time (hours)','FontSize',20)
ylabel('Number of cells','FontSize',20)
%title(['\alpha= ', num2str(alpha)])
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,2,2)
plot(tvec, U,'LineWidth',3)
hold on
%plot(tvec, Ucons, 'LineWidth', 3, 'color', 'k')
xlim([ 0 2*t1-1])
ylim([ 0 1.3])
xlabel('time (hours)','FontSize',20)
ylabel('Effective dose (u(t))','FontSize',20)
%legend('pulsed', 'constant')
%legend boxoff 
set(gca,'FontSize',20,'LineWidth',1.5)
%title('Pulse of 200 nM dox every week')
%% Simulate second dose at two weeks

t2vec = (0:4:336)'; %
tvec = vertcat(t1vec, t2vec(2:end) + t1vec(end));
tdrug = [t1vec(1); (t2vec(1) + t1vec(end))];
% set original U as U1 to indicate first dose

ttest2 = (0:dt:t2vec(end))';
ttestvec = vertcat(ttest', (ttest2(2:end) + 336));
U2 = k*Cdox*exp(-kdrug*(ttest2-t)); % minus one because t starts at 1
U = vertcat(U1(1:end-1), U2);
[Nsr, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U, dt, tdrug);
N = Nsr(:,1);
S = Nsr(:,2);
R = Nsr(:,3);

figure;
subplot(1,2,1)
plot(tvec,N,'LineWidth',3, 'color','b');
hold on
plot(tvec, S, 'LineWidth', 3,'color', 'g')
plot(tvec, R, 'LineWidth', 3, 'color', 'r')
for i = 1:length(tdrug)
plot(tcrit(i) + tdrug(i)-1, Ncrit(i), 'k*', 'LineWidth',3)
text(tcrit(i)+tdrug(i)+2, Ncrit(i), ['t_{crit}=', num2str(tcrit(i)), ' hrs post treat'])
end
legend('total cell number', 'sensitive', 'resistant', 'critical N', 'Location', 'NorthWest')
legend boxoff
xlim([ 0 tvec(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Total Cell Number','FontSize',20)
title('Multiple treatment response')
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,2,2)
plot(ttestvec, U,'LineWidth',3)
xlim([ 0 ttestvec(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Effective dose','FontSize',20)
set(gca,'FontSize',20,'LineWidth',1.5)
title('Dosing regimen')

%% Simulate multiple doses
% Set up easy way to make time and dose vector
int_treat = [2; 2; 2]; % intervals between treatments
cum_treat = cumsum(int_treat);
totweeks = sum(int_treat) + 4; % monitor for two weeks after last treatment
tdrug = [1; cum_treat*24*7];
tvec = [1:1: totweeks*7*24]';
Cdox = [ 75; 75; 75; 75];
acc_dose = cumsum(Cdox)
Uvec = [];
Utot = zeros([length(tvec),1]);
for i = 1:length(tdrug)
% start with zeros 1:tdrug(1)
Uinit = zeros([tdrug(i)-1,1]);
tin = tvec(tdrug(i):end)-tdrug(i);
Udrug = k*Cdox(i)*exp(-kdrug*(tin-1)); % minus one because t starts at 1
% add on each additional treatment
Upulse = vertcat(Uinit, Udrug);
Utot = Utot + Upulse;

end
alpha = 1e-4;
pset = [S0, R0, rs, carcap];
pfit = [ alpha, rr, ds, dr];
p = [ pset, pfit];


[Nsr, tcrit, Ncrit] = fwd_Greene_model(p, tvec, Utot, dt, tdrug);
N = Nsr(:,1);
S = Nsr(:,2);
R = Nsr(:,3);

figure;
subplot(1,3,1)
plot(tvec,N,'LineWidth',3, 'color','b');
hold on
plot(tvec, S, 'LineWidth', 3,'color', 'g')
plot(tvec, R, 'LineWidth', 3, 'color', 'r')
for i = 1:length(tdrug)
plot(tcrit(i) + tdrug(i)-1, Ncrit(i), 'k*', 'LineWidth',3)
text(tcrit(i)+tdrug(i)+2, Ncrit(i), ['t_{crit}=', num2str(tcrit(i)), ' hrs post treat'])
end
legend('total cell number', 'sensitive', 'resistant', 'critical N', 'Location', 'NorthWest')
legend boxoff
xlim([ 0 tvec(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Total Cell Number','FontSize',20)
title('Multiple treatment response')
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,3,2)
plot(tvec, Utot,'LineWidth',3)
xlim([ 0 tvec(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Effective dose (U(t))','FontSize',20)
set(gca,'FontSize',20,'LineWidth',1.5)
title('Dosing regimen')

subplot(1,3,3)
plot(acc_dose, tcrit, 'k*', 'LineWidth',3)
xlabel('Accumulated Dose (nM)','FontSize',20)
ylabel('T_{crit} post last treatment','FontSize',20)
title(['T_{crit} vs. Accumulated Dose. alpha= ',num2str(alpha)])
set(gca,'FontSize',20,'LineWidth',1.5)