function [Nsrdat, tcrit, Ncrit] = fwd_Greene_model(p, tvec, U, dt,tdrug)
% This function runs the Greene model of S & R cells with drug-induced
% resistance forward

P = num2cell(p); 
[S0, R0, rs, carcap, alpha, rr, ds, dr] = deal(P{:}); % our parameters
tend = tvec(end);
ttot = (0:dt:tend)';

S(1)=S0; 
R(1) = R0; 
N(1) = S(1) + R(1);

for t = 2:length(ttot)
    growth_s = rs*S(t-1)*(1-((S(t-1)+R(t-1))/carcap));
    growth_r = rr*R(t-1)*(1-((S(t-1)+R(t-1))/carcap));
    death_s = ds*U(t-1)*S(t-1);
    death_r = dr*U(t-1)*R(t-1);
    trans_s = alpha*U(t-1)*S(t-1);
   
    S(t) = S(t-1) + dt*(growth_s - death_s - trans_s);
    if S(t) <=0
        S(t) = 0;
    end
 
    R(t) = R(t-1) + dt *(growth_r - death_r + trans_s);
    
    if R(t) <=0
        R(t) = 0;
    end
    N(t) = S(t) + R(t);

end

Nsr = horzcat(N', S', R');
 ikeep = find(ismember(ttot,tvec));
 Nsrdat= Nsr(ikeep,:);

 % find multiple critical times
for i = 1:length(tdrug)
    % for each tdrug, search interval after it
    if i == length(tdrug)
        ilow = find(ttot<ttot(end));
    else
    ilow = find(ttot<tdrug(i+1)); % all indices before next drug
    end
    ihigh = find(ttot>=tdrug(i));% all indices after this dose
    ind = find(ismember(ilow, ihigh));
    Npd = N(ind);
    Ncrit(i) = 1.2*Npd(1);
    icrit = find(Npd>Ncrit(i),1, 'first');

    if isempty(icrit)
        tcrit(i) = 0;
    else
        tcrit(i)= ttot(icrit);
    end
end


end