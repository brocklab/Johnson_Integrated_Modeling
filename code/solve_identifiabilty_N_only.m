% Script to attempt to solve non-linear system of ODEs needed that describe
% the injectivity of the parameter mapping of the model of drug-induced
% resistance

%parameters to solve for
syms S0 real;
syms R0 real;
syms pr real;

% observable variables
syms hx0 real;
syms Lfh real;
syms Lf2h real;

% set up system of 3 equations with 3 unknowns
eq1 = S0 + R0 ==hx0;
eq2 = -(S0).^2 + S0 - (1+pr)*S0*R0 - pr*(R0).^2 + pr*R0==Lfh;
eq3 = (1-S0-R0)*S0*(-2*S0+1-(1+pr)*R0) + pr*(1-S0-R0)*R0*(-(1+pr)*S0-2*pr*R0+pr)== Lf2h;

% isolate the three unknowns
[solS0, solR0, solpr] = solve(eq1, eq2, eq3, S0, R0, pr)

%% Now introduce the other parameters and the other Lie derivatives
syms alpha real;
syms ds real;
syms dr real;

% Lie derivatives continued
syms Lgh real;
syms Lg2h real;
syms LfLgh real;
syms LgLfh real;

eq4 = -ds*S0 - dr*R0==Lgh;
eq5 = (alpha +ds)*ds*S0 - dr*alpha*S0 + R0*(dr.^2)== Lg2h;
eq6 = -ds*(1-S0-R0)*S0 - dr*pr*(1-S0-R0)*R0==LfLgh;
eq7 = (alpha + ds)*S0*(2*S0-1+(1+pr)*R0) + (alpha*S0 - dr*R0)*(-(1+pr)*S0-2*pr*R0 + pr)== LgLfh;

[solalpha, solds, soldr] = solve(eq4, eq5, eq6, alpha, ds, dr)