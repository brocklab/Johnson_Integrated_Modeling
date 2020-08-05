% This script loads in the 231 single pulse treatment parameter estimates
% and makes a prediction for a new dose that wasn't used for calibration.

% We then can compare the prediction of the new dose to the experimental
% data.

% We start by simulating completely the effect of a second treatment using
% the forward model and the fitted parameters. We choose a dosing schedule
% that has been observed in our data.

close all; clear all; clc
%% Load in data structure and parameters
S = load('../out/trajfit231.mat');
traj= S.traj;

S = load('../out/trajsumfit231.mat');
trajsum = S.trajsum;

ptest = load('../out/pNphi.mat');
ptest = struct2cell(ptest);
ptest = cell2mat(ptest);
P = num2cell(ptest);

ptestN = load('../out/pN.mat');
ptestN = struct2cell(ptestN);
ptestN = cell2mat(ptestN);
PN = num2cell(ptestN);

% % Also don't really need the confidence intervals
% CI = load('../out/CIpbest.mat');
% CI = struct2cell(CI);
% CI = cell2mat(CI);
% CIphi0 = CI(1,:);
% CIrs = CI(2,:);
% CIalpha = CI(3,:);
% CIzr = CI(4,:);
% CIds = CI(5,:);
% CIzd=CI(6,:);

[phi0, carcapN, carcapphi, rs, alpha, zr, ds, zd, k, kdrug] = deal(P{:});
[phi0N, carcapN, carcapphi, rsN, alphaN, zrN, dsN, zdN, k, kdrug,] = deal(PN{:});

dosevec = [ 1 3 5];
dosevecnew = [2,4,6,7];
%% Now compare predictive power of new doses 
figure;
for i = 1:length(dosevecnew)
    j = dosevecnew(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    %plot(trajsum(j).tvec, trajsum(j).Nmean, '*', 'color', trajsum(j).color, 'LineWidth', 2)
    hold on
end
    %text(trajsum(j).tvec(end), trajsum(j).Nmean(end), ['Cdox= ', num2str(trajsum(j).Cdox), ' nM'], 'FontSize', 12)
for i = 1:length(dosevecnew)
    j = dosevecnew(i);
    plot(trajsum(j).tvec, trajsum(j).Nmodel, 'k-', 'LineWidth',5)
    %plot(trajsum(j).tvec, trajsum(j).NmodelN, 'r-', 'LineWidth',2,'color', trajsum(j).color)
    hold on
end
legend('25 nM', '75 nM','150 nM', '200 nM', 'Location', 'NorthWest')
legend boxoff
xlabel ('time (hours)')
ylabel ('N(t)')
%title('Prediction of N(t) with integrated dit')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
for i = 1:length(dosevecnew)
    j = dosevecnew(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    plot(trajsum(j).tvec, trajsum(j).Nmean, '*', 'color', trajsum(j).color, 'LineWidth', 2)
    hold on
    text(trajsum(j).tvec(end), trajsum(j).Nmean(end), ['Cdox= ', num2str(trajsum(j).Cdox), ' nM'], 'FontSize', 12)
  
    %plot(trajsum(j).tvec, trajsum(j).Nmodel, 'k-', 'LineWidth',5)
    plot(trajsum(j).tvec, trajsum(j).NmodelN, 'r-', 'LineWidth',5)
    
end
xlabel ('time (hours)')
ylabel ('N(t)')
title('Prediction of N(t) using just N(t)')
set(gca,'FontSize',20,'LineWidth',1.5)

%% subplots of individual doses
NpredNphi = [];
NpredN= [];
Ndat = [];
figure;
for i = 1:length(dosevecnew)
    subplot(2, 2, i)
    j = dosevecnew(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    hold on
    
    plot(trajsum(j).tvec, trajsum(j).Nmodel, 'k-', 'LineWidth',4)
    %plot(trajsum(j).tvec, trajsum(j).NmodelN, 'r-', 'LineWidth',4)

    CCCNphi(i) = f_CCC([trajsum(j).Nmodel(:,1), trajsum(j).Nmean], 0.05);
    NpredNphi = vertcat(NpredNphi, trajsum(j).Nmodel(:,1));
    CCCN(i) = f_CCC([trajsum(j).NmodelN(:,1), trajsum(j).Nmean], 0.05);
    NpredN = vertcat(NpredN, trajsum(j).NmodelN(:,1));
    Ndat = vertcat(Ndat, trajsum(j).Nmean);

    xlim([0 trajsum(j).tvec(end)])
    xlabel ('time (hours)')
    ylabel ('N(t)')
    %title(['Cdox= ', num2str(trajsum(j).Cdox), ' nM, CCC= ', num2str(round(CCCNphi(i),3))])
    title(['[Dox]= ', num2str(trajsum(j).Cdox), ' nM'])
    legend('N(t) data', 'integrated model' ,'FontSize', 12,'Location', 'NorthWest')
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
end

figure;
for i = 1:length(dosevecnew)
    subplot(2, 2, i)
    j = dosevecnew(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    hold on
    
    %plot(trajsum(j).tvec, trajsum(j).Nmodel, 'k-', 'LineWidth',4)
    plot(trajsum(j).tvec, trajsum(j).NmodelN, 'r-', 'LineWidth',4)

    CCCNphi(i) = f_CCC([trajsum(j).Nmodel(:,1), trajsum(j).Nmean], 0.05);
    NpredNphi = vertcat(NpredNphi, trajsum(j).Nmodel(:,1));
    CCCN(i) = f_CCC([trajsum(j).NmodelN(:,1), trajsum(j).Nmean], 0.05);
    NpredN = vertcat(NpredN, trajsum(j).NmodelN(:,1));
    Ndat = vertcat(Ndat, trajsum(j).Nmean);

    xlim([0 trajsum(j).tvec(end)])
    xlabel ('time (hours)')
    ylabel ('N(t)')
    title(['Cdox= ', num2str(trajsum(j).Cdox), ' nM, CCC= ', num2str(round(CCCN(i),3))])
    %title(['[Dox]= ', num2str(trajsum(j).Cdox), ' nM'])
    legend('N(t) data', 'N(t) model' ,'FontSize', 12,'Location', 'NorthWest')
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
end
CCC_Nphiall = f_CCC([NpredNphi, Ndat], 0.05);
CCC_Nall = f_CCC([NpredN, Ndat], 0.05);
%% Concordance plots
figure;
hold on
for i = 1:length(dosevecnew)
    j = dosevecnew(i);
   plot(trajsum(j).Nmean, trajsum(j).Nmodel, '*', 'LineWidth',4, 'color', trajsum(j).color) 

hold on
end
plot([0 max(NpredNphi)], [0 max(NpredNphi)], 'k-', 'LineWidth', 3)
xlim([0 max(NpredNphi)])
ylim([0 max(NpredNphi)])
xlabel('N(t) data')
ylabel('Model Predicted N(t)')
%title(['CCC_{integrated fit} predicted doses=', num2str(round(CCC_Nphiall,3))])
legend( '25 nM', '75nM', '150nM', '200 nM',  'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
hold on
for i = 1:length(dosevecnew)
    j = dosevecnew(i);
   plot(trajsum(j).Nmean, trajsum(j).NmodelN, '*', 'LineWidth',4, 'color', trajsum(j).color) 
hold on
end
plot([0 max(NpredN)], [0 max(NpredN)], 'r-', 'LineWidth', 3)
xlim([0 max(NpredN)])
ylim([0 max(NpredN)])
xlabel('N(t) data')
ylabel('Model Predicted N(t)')
title(['CCC_{N(t) fit} predicted doses=', num2str(round(CCC_Nall,3))])
legend( '25 nM', '75nM', '150nM', '200 nM',  'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
%%
%% Start with calibration data comparison 
figure;
for i = 1:length(dosevec)
    j = dosevec(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    plot(trajsum(j).tvec, trajsum(j).Nmean, '*', 'color', trajsum(j).color, 'LineWidth', 2)
    hold on
    text(trajsum(j).tvec(end), trajsum(j).Nmean(end), ['Cdox= ', num2str(trajsum(j).Cdox), ' nM'], 'FontSize', 12)
  
    plot(trajsum(j).tvec, trajsum(j).Nmodel, 'k-', 'LineWidth',5)
    %plot(trajsum(j).tvec, trajsum(j).NmodelN, 'r-', 'LineWidth',2,'color', trajsum(j).color)
    
end
xlabel ('time (hours)')
ylabel ('N(t)')
title('Calibration to N(t) with both data sets')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
for i = 1:length(dosevec)
    j = dosevec(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    plot(trajsum(j).tvec, trajsum(j).Nmean, '*', 'color', trajsum(j).color, 'LineWidth', 2)
    hold on
    text(trajsum(j).tvec(end), trajsum(j).Nmean(end), ['Cdox= ', num2str(trajsum(j).Cdox), ' nM'], 'FontSize', 12)
  
    %plot(trajsum(j).tvec, trajsum(j).Nmodel, 'k-', 'LineWidth',5)
    plot(trajsum(j).tvec, trajsum(j).NmodelN, 'r-', 'LineWidth',5)
    
end
xlabel ('time (hours)')
ylabel ('N(t)')
title('Calibration to N(t) with just N(t)')
set(gca,'FontSize',20,'LineWidth',1.5)

%% subplots of individual doses
NpredNphi = [];
NpredN= [];
Ndat = [];
figure;
for i = 1:length(dosevec)
    subplot(1,3, i)
    j = dosevec(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    hold on
    
    plot(trajsum(j).tvec, trajsum(j).Nmodel, 'k-', 'LineWidth',4)
    %plot(trajsum(j).tvec, trajsum(j).NmodelN, 'r-', 'LineWidth',4)

    CCCNphi(i) = f_CCC([trajsum(j).Nmodel(:,1), trajsum(j).Nmean], 0.05);
    NpredNphi = vertcat(NpredNphi, trajsum(j).Nmodel(:,1));
    CCCN(i) = f_CCC([trajsum(j).NmodelN(:,1), trajsum(j).Nmean], 0.05);
    NpredN = vertcat(NpredN, trajsum(j).NmodelN(:,1));
    Ndat = vertcat(Ndat, trajsum(j).Nmean);

    xlim([0 trajsum(j).tvec(end)])
    xlabel ('time (hours)')
    ylabel ('N(t)')
    %title(['Cdox= ', num2str(trajsum(j).Cdox), ' nM, CCC= ', num2str(round(CCC(i),3)), ', CCC_{pareto}=', num2str(round(CCCpareto(i),3))])
    title(['[Dox]= ', num2str(trajsum(j).Cdox), ' nM'])
    legend('N(t) data', 'integrated model' ,'FontSize', 12,'Location', 'NorthWest')
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
end

figure;
for i = 1:length(dosevec)
    subplot(1, 3, i)
    j = dosevec(i);
    errorbar(trajsum(j).tvec, trajsum(j).Nmean,  1.96*trajsum(j).Nstd/2, '*', 'color', trajsum(j).color)
    hold on
    
    %plot(trajsum(j).tvec, trajsum(j).Nmodel, 'k-', 'LineWidth',4)
    plot(trajsum(j).tvec, trajsum(j).NmodelN, 'r-', 'LineWidth',4)

    CCCNphi(i) = f_CCC([trajsum(j).Nmodel(:,1), trajsum(j).Nmean], 0.05);
    NpredNphi = vertcat(NpredNphi, trajsum(j).Nmodel(:,1));
    CCCN(i) = f_CCC([trajsum(j).NmodelN(:,1), trajsum(j).Nmean], 0.05);
    NpredN = vertcat(NpredN, trajsum(j).NmodelN(:,1));
    Ndat = vertcat(Ndat, trajsum(j).Nmean);

    xlim([0 trajsum(j).tvec(end)])
    xlabel ('time (hours)')
    ylabel ('N(t)')
    %title(['Cdox= ', num2str(trajsum(j).Cdox), ' nM, CCC= ', num2str(round(CCC(i),3)), ', CCC_{pareto}=', num2str(round(CCCpareto(i),3))])
    title(['[Dox]= ', num2str(trajsum(j).Cdox), ' nM'])
    legend('N(t) data', 'N(t) model' ,'FontSize', 12,'Location', 'NorthWest')
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
end
CCC_Nphiall = f_CCC([NpredNphi, Ndat], 0.05);
CCC_Nall = f_CCC([NpredN, Ndat], 0.05);
%% Concordance plots
figure;
hold on
for i = 1:length(dosevec)
    j = dosevec(i);
   plot(trajsum(j).Nmean, trajsum(j).Nmodel, '*', 'LineWidth',4, 'color', trajsum(j).color) 

hold on
end
plot([0 max(NpredNphi)], [0 max(NpredNphi)], 'k-', 'LineWidth', 3)
xlim([0 max(NpredNphi)])
ylim([0 max(NpredNphi)])
xlabel('N(t) data')
ylabel('Model Fit N(t)')
%title(['CCC_{integrated fit} fit doses=', num2str(round(CCC_Nphiall,4))])
legend( '0 nM', '50nM', '100 nM',  'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
hold on
for i = 1:length(dosevec)
    j = dosevec(i);
   plot(trajsum(j).Nmean, trajsum(j).NmodelN, '*', 'LineWidth',4, 'color', trajsum(j).color) 
hold on
end
plot([0 max(NpredN)], [0 max(NpredN)], 'r-', 'LineWidth', 3)
xlim([0 max(NpredN)])
ylim([0 max(NpredN)])
xlabel('N(t) data')
ylabel('Model Fit N(t)')
title(['CCC_{N(t) fit} fit doses=', num2str(round(CCC_Nall,4))])
legend( '0 nM', '50nM', '100nM',  'Location', 'NorthWest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
