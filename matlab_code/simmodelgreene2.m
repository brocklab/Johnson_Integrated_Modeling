function [Y] = simmodelgreene2(pfit,ytimefit, N0s, pset, Uvec, lengthvec, pfID, psetID)
% SIMMODELGreene returns Y values corresponding to times T for parameter vector p
% for James Greene model
% THE MODEL:
% dS/dt = rs(1-(S+R)/K)*S - alpha*u(t)*S - ds*u(t)*S
% dR/dt = rr(1-(S+R)/K)*R + alpha*u(t)*S- dr*u(t)*R ;
% 
% Define inputs
phi = 0;
rs = 0;
carcap1 = 0;
carcap2 = 0;
alpha = 0;
zr = 0;
ds = 0;
zd = 0;
params = [phi, rs, carcap1, carcap2, alpha, zr, ds, zd];
for i = 1:length(params)
    indset= find(ismember(psetID, i));
    if ~isempty(indset)
    params(i) = pset(indset);
    end
    indfit = find(ismember(pfID,i));
    if ~isempty(indfit)
    params(i) = pfit(indfit);
    end
end
P = num2cell(params);
[phi, rs, carcap1,carcap2,alpha, zr, ds, zr] = deal(P{:});


Nsr = [];
tsum = [];
tdrug = 1; % since we are monitoring first treatment only
istart = vertcat(1,cumsum(lengthvec(1:end-1,1))+1);
iend = cumsum(lengthvec(:,1));

ist = vertcat(1,cumsum(lengthvec(1:end-1,2))+1);
ie = cumsum(lengthvec(:,2));

for i = 1:size(lengthvec,1)
    S0 = phi*N0s(i);
    R0 = (1-phi)*N0s(i);
    % vary these based on what we're fitting
    rr = rs*zr;
    dr = ds*zd;
    p = [ S0, R0, rs, carcap1, alpha, rr, ds, dr];
    
    tvec = round(ytimefit(istart(i):iend(i)),0);
    U = Uvec(ist(i):ie(i));
    dt = 1;
    [Nsri, ~, ~] = fwd_Greene_model(p, tvec, U, dt, tdrug);
    tsum = vertcat(tsum, tvec);
    
    Nsr = vertcat(Nsr, Nsri);
end



Y = Nsr(:,1);
end