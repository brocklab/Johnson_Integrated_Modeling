% Calibrate simulated data from a new simple model.

% The goal of this script is to demonstrate the generalizability of the
% data integration framework- i.e. the objective function that combines
% lower time resolution estimates of population composition (i.e. phi(t)
% data) with longitudinal bulk population data (i.e. N(t)).

% In this simulated case, we want to generate data from a simpler model,
% i.e. the biexponential evolutionary model, and compare the accuracy of
% estimating the model parameters from N(t) alone and from N(t) and phi(t)

% The model:
% S(t) = N0(1-fr0)*exp(-kt)
% R(t) = N0*fr0*exp(gt)
% N(t) = S(t) + R(t)
% phi(t) = S(t)./(S(t)+R(t))

% Unknown parameters: 
% g, k, and fr0,

% Input parameter:
% N0

% Scenario 1: N(t) only
% Scenario 2: N(t) and phi(t)

% Let's compare the general accuracy of parameter estimation for both of
% these scenarios 
close all; clear all; clc
%% Set and store the parameters
nsamps = 100;

gdom = linspace(0,0.1, nsamps);
kdom = linspace(0,0.1,nsamps);
fr0dom = linspace(0,1,nsamps);
N0dom = linspace(1e3,1e6,nsamps);

for i = 1:nsamps
    g = randsample(gdom,1);
    k = randsample(kdom,1);
    fr0 = randsample(fr0dom,1);
    N0 = randsample(N0dom,1);
    
    pallstore(i,:) = [g, k, fr0, N0];
end

psetid = 4;
pfitid = [1, 2, 3];
psetstore = pallstore(:,psetid);
pfitstore = pallstore(:,pfitid);
%% Generate simulated data
eta = 50;
time = 0:1:100;
ind= [1, 75, 101]; % time of phi(t) data
tphi = time(ind)
for i = 1:nsamps
    P = num2cell(pallstore(i,:));
    [g, k, fr0, N0] = deal(P{:});
 
    prun = [g, k, fr0, N0]
    % Generate in silico data and store it
    [NSR] = bi_exp_model(prun, time);
    N_tdata = NSR(:,1) + normrd(0, eta,[length(time) 1])
    phi_tdata = NSR(:,2) + normrd(0,era, [length(time) 1])./N_tdata;
    phi_tobs = phi_tdata(ind)
    % Fit to in silico data using N(t)only
    lambda = 0;
    [pest, Nest] = fit_bi_exp(N_tdata, time, phi_tobs, tphi, lambda, sigma);
end
