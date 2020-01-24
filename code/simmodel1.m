function [Y] = simmodel1(p,T, V_0)
% SIMMODEL returns Y values corresponding to times T for parameter vector p
% for single exponential model

% Single Exponential
% $$ N(t) = N0*K*(1/(N0 + (K-N0)exp(-gst));
% where:
% * $N_0$ is the initial cell numer
% * K is the carrying capacity
% * gs>0 is the growth rate on all cells

P = num2cell(p); 
[ gs, K] = deal(P{:}); % our parameters

TT = T';
Ntimes = length(TT);

% initialize solutions
Y(1,1) = V_0.*exp(TT(1)*gs);

for j = 1:length(TT)
    if TT(j)<=0
        Y(j,1) = V_0;
    end
    if TT(j) > 0 
    Y(j,1) = V_0*K*(1./(V_0 + ((K-V_0)*exp(-gs*TT(j)))));
    % need to write something to account for when Y(j,1) <= 0...
    end
    
  
end

% return only predictions at time points originally passed in T
ikeep = find(ismember(TT,T));
Y = Y(ikeep,1);
end