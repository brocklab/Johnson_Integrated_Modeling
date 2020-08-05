function [err] = objfunN(theta,yN, timeN,N0s, UN, lengthvecN, pfitID, psetID, pset,sigmafit)
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
params = [prop, rs, carcapN,carcapphi, alpha, zr, ds, zd];
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

nNdata = length(yN);

err = ((modelfunN(theta)-yN).^2)./(sigmafit.^2);

end

