% This script loads in the traj structure which contains the dose response
% for all conditions. It creates the trajsum structure which will contain a
% mean vector and U(t) vector for each different concentration tested. This
% will be used to calibrate to the model.

% Start by writing a model that generates a prediction for a given time
% vector and U(t) for multiple concentrations and concatenates them. (Will
% turn this into a function).
close all; clear all; clc
%%
% Set some arbitrary things as inputs
Cdox = [10 100];
kdrug = 0.15;
k = 0.5;
dt = 1;
% input time vectors for each different dose response
timevec1 = 0:4:336; 
timevec2 = 0:4:332;
trajsum(1).tvec = timevec1;
trajsum(2).tvec = timevec2;

for i = 1:2
    ttest = [];
    ttest = 1:dt:trajsum(i).tvec(end);
    trajsum(i).U = k*Cdox(i)*exp(-kdrug*(ttest)); 
end

% set parameters of forward model
carcap = 5e4;
S0=2e3; % set initial conditions (assume all cells sensitive to start)
R0 = 0; 
Cdox = 75; % nm doxorubicin
rs = 0.0287;
rr = 0.005;
ds = 0.002;
dr = 0;
alpha = 0.0001;

pset = [S0, R0, rs, carcap];
pfit = [ alpha, rr, ds, dr];
p = [ pset, pfit];


Nsr = [];
tsum = [];
tdrug = 1; % since we are monitoring first treatment only

for i = 1:2
    tvec = trajsum(i).tvec(2:end)';
    U = trajsum(i).U;
    [Nsri, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U, dt, tdrug);
    tsum = vertcat(tsum, tvec);
    
    Nsr = vertcat(Nsr, Nsri);
end

figure;
plot(tsum, Nsr(:,1), 'b*')
xlabel('time (hours)')
ylabel('Model predicted response to pulse treat')
title('Test ability to generate model trajectories')

