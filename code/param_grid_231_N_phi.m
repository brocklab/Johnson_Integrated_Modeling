% param_grid_231_N_phi

% This script is intended to run in parallel to fit_231_N_phi
% We want to find the pareto boundary but not by fitting- simply by
% evaluating cost functions in a grid of parameter space and outputting the
% error in N(t) and error in phi(t)

% Then, we will identify the parameters that are within some optimal range
% of error in N(t) & phi(t) just by looking at this boundary and its
% corresponding parameters

% Start by making the grid of fit parameters:
% Your fit parameters are:
% rs, zr, alpha, ds, zd
theta = [rsguess,alphaguess, zrguess, dsguess, zdguess];

% Vary them within some reasonable range:
npts = 10; % start with a course grid
rsvec = linspace(0.5*gtot, 2*gtot, npts);
zrvec = linspace(1e-3, 0.8, npts);
alphavec = linspace(1e-3, 0.5, npts);
dsvec = linspace(1e-3, 0.1, npts);
zdvec = linspace(1e-3, 0.8, npts);

pmat = vertcat(rsvec, alphavec, zrvec, dsvec, zdvec)
%%
% write a function that takes set parameters and fit parameters, combines
% them, and runs the model forward for the Uvec and 0s provided
modelfunN = @(p)simmodelgreene2(p, ytimefit, N0s, pset, Uvec, lengthvec, pfitID, psetID); 
modelfunphi = @ (p)simmodelgreenephi2(p, tbot, N0phi, pset, Ub, lengthvecphi, pfitID, psetID);
%%
paretomat = [];

for i = 1:npts
    for j = 1:npts
        for k = 1:npts
            for l = 1:npts
                for m = 1:npts

    ptemp = [rsvec(i), alphavec(k), zrvec(j), dsvec(l), zdvec(m)];
    
    Nfwd = modelfunN(ptemp);
    phifwd = modelfunphi(ptemp);
    
    weighted_errN = sum(((Nfwd-ydatafit).^2)./sigmafit);
    weighted_errphi=sum(((phifwd-phitrt).^2)./phisigfit);
    
    % make a matrix that stacks parametrs, then weighted errors in N and
    % phi
    paretocol(1:length(ptemp),1) = ptemp'; % vertical parameters
    paretocol(length(ptemp)+1,1) = weighted_errN;
    paretocol(length(ptemp)+2,1) = weighted_errphi;
    
    paretomat = horzcat(paretomat, paretocol);
                end
            end
        end
    end
end
%% Randomly sample space logarithmically
paretomat = [];
% Vary them within some reasonable range:
npts = 1e7; % start with a course grid
rsvec = linspace((min(param_table.rsvals)), (max(param_table.rsvals)), npts);
zrvec = linspace((min(param_table.rr_rs_ratio)), (max(param_table.rr_rs_ratio)), npts);
alphavec = linspace((min(param_table.alphavals)), (max(param_table.alphavals)), npts);
dsvec = linspace((min(param_table.dsvals)), (max(param_table.dsvals)), npts);
zdvec = linspace(min(param_table.dr_ds_ratio), max(param_table.dr_ds_ratio), npts);
%%
for i = 1:1e6
    % generate 5 random intergers between 1 and npts
    inds = randi([1, npts],5,1);
    ptemp = [rsvec(inds(1)), alphavec(inds(2)), zrvec(inds(3)), dsvec(inds(4)), zdvec(inds(5))];
    
    Nfwd = modelfunN(ptemp);
    phifwd = modelfunphi(ptemp);
    
    weighted_errN = sum(((Nfwd-ydatafit).^2)./(sigmafit.^2));
    weighted_errphi=sum(((phifwd-phitrt).^2)./(phisigfit.^2));
    
    % make a matrix that stacks parametrs, then weighted errors in N and
    % phi
    paretocol(1:length(ptemp),1) = ptemp'; % vertical parameters
    paretocol(length(ptemp)+1,1) = weighted_errN;
    paretocol(length(ptemp)+2,1) = weighted_errphi;
    
    paretomat = horzcat(paretomat, paretocol);
end
    
%% Make into a table and export
rsvals = paretomat(1,:)';
alphavals = paretomat(2,:)';
rr_rs_ratio= paretomat(3,:)';
dsvals = paretomat(4,:)';
dr_ds_ratio = paretomat(5,:)';
err_N = paretomat(6,:)';
err_phi = paretomat(7,:)';
%%
pareto_table_big = table(rsvals, alphavals, rr_rs_ratio, dsvals, dr_ds_ratio, err_N, err_phi);
save('../out/paramgridbig.mat', 'pareto_table_big');

%% Plot the error in N and phi
figure
plot(paretomat(length(ptemp)+1, :), paretomat(length(ptemp)+2, :), '.')
hold on
%plot(sum(((err_N.^2)./sigmafit)), sum(((err_phi.^2)./phisigfit)), 'r*', 'LineWidth', 3)
plot(sumsqerrN,sumsqerrphi,'r.')
xlabel('weighted error in N(t)')
ylabel('weighted error in \phi(t)')
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
%legend('manual search', 'optimizer')
title('Manual search for pareto boundary')

figure
plot(paretomat(length(ptemp)+1, :), paretomat(length(ptemp)+2, :), '.')
hold on
plot(sumsqerrN,sumsqerrphi,'r.')
xlabel('weighted error in N(t)')
ylabel('weighted error in \phi(t)')
legend('manual search', 'optimizer')
xlim([ 0 0.6e7])
ylim([0 1])
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log', 'Yscale', 'log')
title('Low error parameter sets')

%% Make a smaller pareto table with low error values 
paretosel = [];

for i = 1:length(paretomat)
    if paretomat(length(ptemp)+1,i) < 1e10 && paretomat(length(ptemp)+2,i) <1e-2
        paretocol = paretomat(:,i);
        paretosel = horzcat(paretosel, paretocol);
    end
        
end
%% Plot some things
pointsize = 30;
figure;
scatter3(paretosel(1,:), paretosel(2,:), paretosel(3,:), pointsize, paretosel(length(ptemp)+1,:), 'filled')
xlabel('r_{s}')
ylabel('\alpha')
zlabel('rr/rs ratio')
colorbar
title('params plotted by error in N(t)')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
scatter3(paretosel(1,:), paretosel(4,:), paretosel(5,:), pointsize, paretosel(length(ptemp)+1,:), 'filled')
xlabel('r_{s}')
ylabel('d_{s}')
zlabel('ds/dr ratio')
colorbar
title('params plotted by error in N(t)')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
scatter3(paretosel(1,:), paretosel(2,:), paretosel(3,:), pointsize, paretosel(length(ptemp)+2,:), 'filled')
xlabel('r_{s}')
ylabel('\alpha')
zlabel('rr/rs ratio')
colorbar
title('params plotted by error in \phi(t)')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
scatter3(paretosel(1,:), paretosel(4,:), paretosel(5,:), pointsize, paretosel(length(ptemp)+2,:), 'filled')
xlabel('r_{s}')
ylabel('d_{s}')
zlabel('ds/dr ratio')
colorbar
title('params plotted by error in \phi(t)')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Look at all parameters
pointsize = 30;
figure;
scatter3(paretomat(1,:), paretomat(2,:), paretomat(3,:), pointsize, paretomat(length(ptemp)+1,:), 'filled')
xlabel('r_{s}')
ylabel('\alpha')
zlabel('rr/rs ratio')
colorbar
title('params plotted by error in N(t)')
set(gca,'FontSize',20,'LineWidth',1.5)

figure
scatter3(paretomat(1,:), paretomat(2,:), paretomat(3,:), pointsize, paretomat(length(ptemp)+2,:), 'filled')
xlabel('r_{s}')
ylabel('\alpha')
zlabel('rr/rs ratio')
colorbar
title('params plotted by error in \phi(t)')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
scatter3(paretomat(1,:),paretomat(4,:), paretomat(5,:), pointsize, paretomat(length(ptemp)+2,:), 'filled')
xlabel('r_{s}')
ylabel('d_{s}')
zlabel('ds/dr ratio')
colorbar
title('params plotted by error in \phi(t)')
set(gca,'FontSize',20,'LineWidth',1.5)
%%
pnames = {'rs', '\alpha', 'rr/rs', 'ds', 'dr/ds'};
    
figure;
for i = 1:length(ptemp)
    subplot(2, length(ptemp), i)
    plot(paretosel(i,:), paretosel(length(ptemp)+1,:), 'r.')
    xlim([pmat(i,1) pmat(i,end)])
    xlabel([pnames(i), 'value'])
    ylabel('error in N(t)')
    set(gca,'FontSize',20,'LineWidth',1.5)
    subplot(2, length(ptemp), i+length(ptemp))
    plot(paretosel(i,:), paretosel(length(ptemp)+2,:), 'g.')
    xlim([pmat(i,1) pmat(i,end)])
    xlabel([pnames(i),'value'])
    ylabel('error in \phi(t)')
    set(gca,'FontSize',20,'LineWidth',1.5)
end

figure;
for i = 1:length(ptemp)
    subplot(2, length(ptemp), i)
    plot(paretomat(i,:), paretomat(length(ptemp)+1,:), 'r.')
    xlabel([pnames(i), 'value'])
    ylabel('error in N(t)')
    set(gca,'FontSize',20,'LineWidth',1.5)
    subplot(2, length(ptemp), i+length(ptemp))
    plot(paretomat(i,:), paretomat(length(ptemp)+2,:), 'g.')
    xlabel([pnames(i),'value'])
    ylabel('error in \phi(t)')
    set(gca,'FontSize',20,'LineWidth',1.5)
end
%% How variable are the parameters within this error range?
figure;
for i = 1:length(ptemp)
    subplot(2, length(ptemp), i)
    scatter(1:1:length(paretosel), paretosel(i,:), pointsize, paretosel(6,:), 'filled')
    colorbar
    xlabel('iteration')
    ylabel([pnames(i), 'value']) 
    ylim([pmat(i,1) pmat(i,end)])
    set(gca,'FontSize',20,'LineWidth',1.5)
    subplot(2, length(ptemp), i+length(ptemp))
    scatter(1:1:length(paretosel), paretosel(i,:), pointsize, paretosel(7,:), 'filled')
    colorbar
    xlabel('iteration')
    ylabel([pnames(i), 'value']) 
    ylim([pmat(i,1) pmat(i,end)])
    set(gca,'FontSize',20,'LineWidth',1.5)
    
end
%% Make into a table and export
rsvals = paretosel(1,:)';
alphavals = paretosel(2,:)';
rr_rs_ratio= paretosel(3,:)';
dsvals = paretosel(4,:)';
dr_ds_ratio = paretosel(5,:)';
err_N = paretosel(6,:)';
err_phi = paretosel(7,:)';

pareto_table = table(rsvals, alphavals, rr_rs_ratio, dsvals, dr_ds_ratio, err_N, err_phi)
save('../out/paramgrid.mat', 'pareto_table')
%% 
np = length(ptemp);
itvec = 1:1:np*np;
ct_iter = 0;
figure;
for i = 1: length(pbest)
    for j = 1:length(ptemp)
        ct_iter = ct_iter+1;
        subplot(length(ptemp), length(ptemp),itvec(ct_iter))
        scatter(paretosel(i,:), paretosel(j,:),10, paretosel(6,:), 'filled')
        colorbar
        set(gca,'FontSize',10,'LineWidth',1.5)
        xlabel([pbestnames(i)])
        ylabel([pbestnames(j)])
        if i ==j
            hist(paretosel(i,:))
            xlabel([pnames(i)])
        end
    end
end
