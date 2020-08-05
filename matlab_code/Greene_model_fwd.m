% James Greene model run forward
% This script is to code up the forward model presented in the James Greene
% paper which consists of sensitive and resistant cell compartments. Both
% subpopulations have a different growth rate and drug induced death rate.
% Both the drug induced death rate and the drug-induced resistance term
% alpha are functions of the effective dose, u(t). We will assume an
% exponential decay of effective dose after drug is administered

% We want to run the forward model in response to pulse treatment at
% different doses, and capture the output 

% Assumptions:
% logistic growth + log-kill hypothesis
% no spontaneous resistance (eta =0)
% no resensitization (gamma = 0), resistance is irreversible)
% effective dose decays exponentially (at time of treatment set u(t) =
% C0exp(-kt)
% rr <rs (resistant cells grow more slowly)

% THE MODEL:
% dS/dt = rs(1-(S+R)/K)*S - alpha*u(t)*S - ds*u(t)*S
% dR/dt = rr(1-(S+R)/K)*R + alpha*u(t)*S- dr*u(t)*R

close all; clear all; clc
dt = 1;
carcap = 2e4; % K
tend = 500; % monitor for 500 hours
ttot = 1:dt:tend;
S = zeros([tend,1]);
R = zeros([tend,1]);
N = zeros([tend,1]);
U = zeros([tend,1]);
S(1)=2e3; % set initial conditions (assume all cells sensitive to start)
R(1) = 0; 
Cdox = 75; % nm doxorubicin
rs = 0.01;
rr = 0.005;
ds = 0.002;
dr = 0;
alpha = 0.0001;
kdrug = 0.025;
N(1) = S(1)+R(1);
t = 1; % time in hours
k = 0.5; 
U=k*Cdox*exp(-kdrug*(ttot-1));
for t = 2:tend
    growth_s = rs*S(t-1)*(1-((S(t-1)+R(t-1))/carcap));
    growth_r = rr*R(t-1)*(1-((S(t-1)+R(t-1))/carcap));
    death_s = ds*U(t-1)*S(t-1);
    death_r = dr*U(t-1)*R(t-1);
    trans_s = alpha*U(t-1)*S(t-1);
   
    S(t) = S(t-1) + dt*(growth_s - death_s - trans_s);
    R(t) = R(t-1) + dt *(growth_r - death_r + trans_s);
    N(t) = S(t) + R(t);
end
Ncrit = N(1)*2;

icrit = find(N>Ncrit,1, 'first');
tcrit= ttot(icrit);

subplot(1,2,1)
plot(1:t,N,'LineWidth',3, 'color','b');
hold on
plot(1:t, S, 'LineWidth', 3,'color', 'g')
plot(1:t, R, 'LineWidth', 3, 'color', 'r')
plot(tcrit, Ncrit, 'k*')
legend('total cell number', 'sensitive', 'resistant', 'critical N', 'Location', 'NorthWest')
legend boxoff
xlim([ 0 ttot(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Total Cell Number','FontSize',20)
set(gca,'FontSize',20,'LineWidth',1.5)
drawnow
subplot(1,2,2)
plot(1:t, U,'LineWidth',3)
xlim([ 0 ttot(end)])
xlabel('Time (hours)','FontSize',20)
ylabel('Effective dose','FontSize',20)
set(gca,'FontSize',20,'LineWidth',1.5)


%% Write this as a function to run model forward
%P = num2cell(p); 
%[V_0, k] = deal(P{:}); % our parameters
% for now, write a function that takes all parameters in
S0 = 2e3;
R0 = 0;
pset = [S0, R0, rs, carcap];
pfit = [ alpha, rr, ds, dr];
p = [ pset, pfit];

% Other inputs that will be known
tend = 1500; % monitor for 500 hours
ttot = 1:2:tend;
tvec = ttot; % time vector of data measurements (Can vary from dt)
ttest = 1:dt:tend;
U=k*Cdox*exp(-kdrug*(ttest'-1)); % U needs to be in terms of dt 
dt = 1; % how frequently to update model
Ncrit = (S0+R0)*2;

% Look at tcrit as a function of dose
Cdoxvec = [10:10:250];
alphavec = [ 0; 1e-4; 5e-4; 1e-3; 5e-3; 1e-2];
Nsr = [];
tcrit = [];
for j = 1:length(alphavec)
for i = 1:length(Cdoxvec)
U(:,i)=k*Cdoxvec(i)*exp(-kdrug*(ttest'-1));
pj = p;
pj(5) = alphavec(j);
tdrug = 1;
[Nsr(:,:,i), tcrit(i,j), Ncrit(i,j)] = fwd_Greene_model(pj, tvec, U(:,i), dt, tdrug);
end
end


figure;
for j = 1:length(alphavec)
plot(Cdoxvec, tcrit(:,j),'*-', 'LineWidth',3)
hold on
xlabel('Concentration Dox','FontSize',20)
ylabel('t_{crit}','FontSize',20)
set(gca,'FontSize',20,'LineWidth',1.5)
end
legend('\alpha = 0', '\alpha=1e-4','\alpha=5e-4', '\alpha = 1e-3','\alpha = 5e-3','\alpha =1e-2', 'Location', 'NorthWest')
legend boxoff
title('Dox concentration vs. tritical time for different drug induction rates (\alpha)')


%% Plot example for a single alpha

figure;

for i = 1:length(Cdoxvec)
subplot(1,2,1)
plot(tvec,Nsr(:,1,i),'LineWidth',3);
hold on

plot(tcrit(i,end), Ncrit(i,end), 'k*')
legend(['Cdox =', num2str(Cdoxvec(i)),' nM'])
legend boxoff
xlabel('Time (hours)','FontSize',20)
ylabel('Total Cell Number','FontSize',20)
set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,2,2)
plot(tvec, U(:,i),'LineWidth',3)
hold on 
legend(['Cdox =', num2str(Cdoxvec(i)),' nM'])
legend boxoff
xlabel('Time (hours)','FontSize',20)
ylabel('Effective dose','FontSize',20)
set(gca,'FontSize',20,'LineWidth',1.5)
end

