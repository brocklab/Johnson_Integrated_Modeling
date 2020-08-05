function [err] = objfun(theta,yN, yphi, timeN, timephi,N0s, N0phi, UN, Uphi, lengthvecN, lengthvecphi, pfitID, psetID, pset, sigmafit, phisigfit)
% This function returns the difference between the model (N and phi)
% together and the measurement 
prop = 0;
rs = 0;
carcapN = 0;
carcapphi = 0;
alpha = 0;
zr = 0;
ds = 0;
zd = 0;
params = [prop, rs, carcapN, carcapphi, alpha, zr, ds, zd];
for i = 1:length(params)
    indset= find(ismember(psetID, i));
    if ~isempty(indset)
    params(i) = pset(indset);
    end
    indfit = find(ismember(pfitID,i));
    if ~isempty(indfit)
    params(i) = theta(indfit);
    end
end


modelfunN = @(p)simmodelgreene2(p, timeN, N0s, pset, UN, lengthvecN, pfitID, psetID); 
modelfunphi = @ (p)simmodelgreenephi2(p, timephi, N0phi, pset, Uphi, lengthvecphi, pfitID, psetID);
nNdata = length(yN);
nphidata = length(yphi);
err = vertcat((1/nNdata).*(modelfunN(theta)-yN).^2./(sigmafit.^2), ((1/nphidata).*(modelfunphi(theta)-yphi).^2./(phisigfit.^2)));
%loglikelihoodN = @(pval)(log(normpdf((yN),modelfunN(pval), sigmafit)));
%loglikelihoodphi = @(pval)(log(normpdf(yphi,modelfunphi(pval), phisigfit)));
%err = vertcat(loglikelihoodN(theta)./nNdata, loglikelihoodphi(theta)./nphidata);

end

