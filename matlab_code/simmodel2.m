function [Y] = simmodel2(p,T, V_0)
% SIMMODEL returns Y values corresponding to times T for parameter vector p
% for double exponential model

% $$ V(t) = V_0 [ \phi e^{gt} + (1-\phi) e^{-kt} ] $$
% where:
% * $V_0$ is the initial volume
% * $\phi$ is the initial resistant fraction
% * g>0 is the resistant growth rate
% * k>0 is the kill rate (actually -net "growth rate" on treatment) on sensitive cells

P = num2cell(p); 
[phi, g, k] = deal(P{:}); % our parameters

TT = T;
Ntimes = length(TT);

% initialize solutions

%Y(1,1) = (V_0.*((1-phi)*exp(-k*TT(1)) + phi*exp(g*TT(1))));

for j = 1:length(TT)
    if TT(j)<=0
        Y(j,1) = V_0;
    end
    if TT(j)>0
    Y(j,1)= (V_0.*((1-phi)*exp(-k*TT(j)) + phi*exp(g*TT(j))));
    end
end

% return only predictions at time points originally passed in T
ikeep = find(ismember(TT,T));
Y = Y(ikeep,1);
end