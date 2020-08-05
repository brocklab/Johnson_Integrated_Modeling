% This script is meant to be used to perform a global sensitivity analysis
% on the parameter ranges derived from the calibrated parameters.  The
% results of this sensitivity analysis will reveal the most important
% parameters of the system, causing the greatest variation in outputs, for
% the area of parameter space for which the model is able to replicate the
% experimental data and the area of uncertainty.

% The important parameters have been shown to be drivers in changing the
% stability of steady states in mathematical models and may drive future
% experimental investigations and/or provide support to theories about the
% biology. Could also be used to reveal potential simplified versions of
% the model. 

% We start by doing this on the simplest model of:
% THE MODEL:
% dS/dt = rs(1-(S+R)/K)*S - alpha*u(t)*S - ds*u(t)*S
% dR/dt = rr(1-(S+R)/K)*R + alpha*u(t)*S- dr*u(t)*R

% For now start by keeping the u(t) parameters out of the analysis, but can
% put these in if we need to later

close all; clear all; clc
%% Load in data structure and parameters
S = load('../out/trajfit231.mat');
traj= S.traj;
Ssum = load('../out/trajsumfit231.mat');
trajsum = Ssum.trajsum;

% Use the fit parameters from N(t) only to make the parameter domains
% We shouldn't really need these... but maybe keep just as a test
ptest = load('../out/ptest.mat');
ptest = struct2cell(ptest);
ptest = cell2mat(ptest);
P = num2cell(ptest);
% Also don't really need the confidence intervals
CI = load('../out/CIpbest.mat');
CI = struct2cell(CI);
CI = cell2mat(CI);
CIrs = CI(1,:);
CIalpha = CI(2,:);
CIzr = CI(3,:);
CIds = CI(4,:);
CIzd=CI(5,:);
% Can use some of these as first guesses/ballpark ranges of what to expect
[phi0f, carcapNf, carcapphif, rsf, alphaf, zrdata, dsf, zdf, k, kdrug, gtot]= deal(P{:});
%[rsf, carcapNf, alphaf, rrf, dsf, drf, k, kdrug, gtot] = deal(P{:});
% We will use this carcapNf only, and the k and kdrug to be consistent

phi_est_filename = '../data/phi_t_est.csv';
phi_est = readtable(phi_est_filename);
tbot = phi_est.t;
phitrt = phi_est.phi_t;
tscRNAseq = tbot(end)+4*24; % add a 4 day buffer since we need to see reach 2*N0
ntrt = phi_est.ncells;

%% Plot the average data for the set of controls (u(t)s) we want to use
  figure;
 for i = 1:length(trajsum)
     subplot(2,1,2)
         plot(trajsum(i).tvec, trajsum(i).Nmean, 'color', trajsum(i).color, 'LineWidth', 2)
         hold on
         text(trajsum(i).tvec(end-10), trajsum(i).Nmean(end-10), ['C_{dox}= ', num2str(trajsum(i).Cdox),' nM'])
         plot(trajsum(i).tvec, trajsum(i).Nmean + 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
         plot(trajsum(i).tvec, trajsum(i).Nmean - 1.96*trajsum(i).Nstd, 'color', trajsum(i).color)
        xlabel('time (hours)')
        ylabel('N(t)')
        title('N(t)')
        set(gca,'FontSize',20,'LineWidth',1.5)
        dt = 1;
       subplot(2,1,1)
       ttest = [];
       ttest = 0:dt:trajsum(i).tvec(end);
       plot(ttest, trajsum(i).U,'.', 'color',trajsum(i).color, 'LineWidth',1)
        hold on
        xlabel('time (hours)')
        ylabel('Effective dose U(t)')
        title('U(t)')
        set(gca,'FontSize',20,'LineWidth',1.5)
 end
 %% Plot experimental data of phi(t)
 sigtech = 0.5*1e-1;
phisigfit = [phitrt.*(1-phitrt)./ntrt] + sigtech;
figure;
errorbar(tbot, phitrt, phisigfit, 'g*', 'LineWidth', 3)
hold on
errorbar(tbot, 1-phitrt, phisigfit, 'r*', 'LineWidth', 3)
legend('\phi(t)=\phi_{S}(t)', '1-phi(t)=\phi_{R}(t)', 'Location', 'Northwest')
legend('boxoff')
ylim([-.1, 1.1])
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('\phi(t)')
xlim([-50,tscRNAseq]) 
%title('Phenotypic composition estimates from scRNAseq')
%% Run loop to generate output needed from each control u(t)
% If we want to output the effect of individual parameter changes on the
% "model behavior" and we quantify model behavior as some observable
% feature about the treatment response-- then we could say let's have this
% output tcrit versus dox concentration & phi(tseq) versus treatment (i.e.
% phenotypic composition at the time of treatment
Cdoxmax = 1000;
Cdoxvec = [0, 75, 200, 500]; % find tcrit and phi for the relevant experimental dox concentrations
%Cdoxvec =[0:10:500];
%Cdoxvec = 200; % make this just one output

tvec1=[0:3:600];
tgen1 =[0:1:tvec1(end)];
tvec2 =[0:1:1.5*tscRNAseq];
tgen2 = [0:1:1.5*tscRNAseq]; % make this real long so that we always observe tcrit even for high doses
tdrug = 1;
for i = 1:length(Cdoxvec)
U1(:,i)=k*Cdoxvec(i)*exp(-kdrug*(tgen1))/(0.1*Cdoxmax);
dt = 1;
tdrug = 1;
N0 = 2e3;
p1 = [ptest(1:2), ptest(4:8)];
% Use this to find the critical times as a function of dose in 96 well
[Nsrdat1(:,:,i), tcrit1(i), Ncrit1(i)]=fwd_Greene_model2(p1, tvec1, N0, U1(:,i), dt,tdrug);
% Use this to find the phi(@tscRNAseq) as a function of dose in expansion
U2(:,i)=k*Cdoxvec(i)*exp(-kdrug*(tgen2))/(0.1*Cdoxmax);
p2 =p1;
% replace carcaph
p2(2) = carcapphif;
[Nsrdat2(:,:,i), tcrit2(i), Ncrit2(i)]=fwd_Greene_model2(p2, tgen2, N0, U2(:,i), dt,tdrug);

% Use this to find phi@tSCRNAseq
% This could be found at any time, we just pick last bc last is the time of
% scNAseq
phi_sc(i) = Nsrdat2(tcrit2(i)-1,2,i)/Nsrdat2(tcrit2(i)-1,1,i);
end
%% Plot an example and then the resulting tcrit and phi(scRNAseq) vs. Cdox.
figure;
subplot(1,3,1)
plot(tgen1, U1(:,2), 'k', 'LineWidth', 2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('Effective dose (u(t))')
title(['Dosing input, Cdox= ', num2str(Cdoxvec(2)), ' nM'])
xlim([0, tgen1(end)])

subplot(1,3,2)
plot(tvec1, Nsrdat1(:,1,2), 'b', 'LineWidth', 2)
hold on
plot(tvec1,Nsrdat1(:,2,2), 'g', 'LineWidth', 2)
plot(tvec1, Nsrdat1(:,3,2), 'r', 'LineWidth',2)
plot(tcrit1(2), Ncrit1(2), 'k*', 'LineWidth', 2)
text(tcrit1(2), Ncrit1(2), ['t_{crit}= ',num2str(tcrit1(2)),' hours'], 'FontSize',14,'HorizontalAlignment', 'left')
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('N(t)')
title('N(t) from 96 well')
xlim([0, tgen1(end)])

subplot(1,3,3)
plot(tvec2, Nsrdat2(:,1,2), 'b', 'LineWidth', 2)
hold on
plot(tvec2,Nsrdat2(:,2,2), 'g', 'LineWidth', 2)
plot(tvec2, Nsrdat2(:,3,2), 'r', 'LineWidth',2)
plot(tcrit2(2), Ncrit2(2), 'k*', 'LineWidth', 2)
text(tcrit2(2), Ncrit2(2), ['\phi(@tcrit)=',num2str(phi_sc(2))], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 14)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('Number of cells')
title('N(t) from scRNAseq')
set(gca,'FontSize',20,'LineWidth',1.5)
xlim([0, tgen2(end)])

figure;
subplot(1,2,1)
plot(Cdoxvec, tcrit1, 'b-', 'LineWidth', 2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('[Dox]')
ylabel('t_{crit}')
title('t_{crit}')
xlim([0 Cdoxvec(end)])
subplot(1,2,2)
plot(Cdoxvec, phi_sc, 'm-', 'LineWidth', 2)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('[Dox]')
ylabel('\phi_{s}')
title('\phi_{s}(t=t_{crit})')
xlim([0 Cdoxvec(end)])
%% Plot some example trajectories

rs = 0.01;
rr = 0.005;
ds = 0.04;
dr = 0;
alpha = 0.001;
pex = [0.8, carcapNf, rs, alpha, 0.5, ds, 0.1];
[Nsrdat, tcrit, Ncrit]=fwd_Greene_model2(pex, tvec1, N0, U1(:,2), dt,tdrug);

figure;
plot(tvec1, Nsrdat(:,1), 'b', 'LineWidth', 4)
hold on
%plot(tvec1,Nsrdat(:,2), 'g-', 'LineWidth', 2)
%plot(tvec1, Nsrdat(:,3), 'r-', 'LineWidth',2)
legend('N(t)', 'S(t)', 'R(t)', 'Location', 'Northwest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('N(t)')
xlim([0, tvec1(end)])

figure;
plot(tvec1, (Nsrdat(:,2)./Nsrdat(:,1)), 'g', 'LineWidth', 4)
hold on
plot(tvec1, (Nsrdat(:,3)./Nsrdat(:,1)), 'r', 'LineWidth',4)
legend('\phi(t)=\phi_{S}(t)', '1-\phi(t)=\phi_{R}(t)', 'Location', 'Northwest')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('time (hours)')
ylabel('\phi(t)')
xlim([0, tvec1(end)])

%% Perform a local sensitivity analysis
% tornado_plot of sensitivities of model at the default values for other parameters
 

par = ptest(1:8);
parlow = par/2;
parhigh = par*2;

paramnames = {'\phi_{s}(0)', 'K_{N}', 'K_{\phi}', 'r_{s}', '\alpha', 'r_{r}/r_{s} ratio', 'd_{s}', 'd_{r}/d_{s} ratio'};

Npar = length(par);

% Write in line functions that output tcrit and phi_sc 
modelfuntcrit=@(p)drug_ind_wrapper(p,Cdoxvec);
modelfunphisc=@(p)drug_ind_wrapper2(p,Cdoxvec);

defvaltcrit = modelfuntcrit(par)';
defvalphisc = modelfunphisc(par)';

for j = 1:Npar
	plow = par; phigh = par;
    % Replace each parameter one by one with high and low value
	plow(j)=parlow(j);
    % record the result of lowering the jth parameter
	lovaltcrit(:,j) = modelfuntcrit(plow)';
    lovalphisc(:,j) = modelfunphisc(plow)';
    % Replace jth parameter with high value
	phigh(j) = parhigh(j);
    % record the result of raising the jth parameter
	hivaltcrit(:,j) = modelfuntcrit(phigh)';
    hivalphisc(:,j) = modelfunphisc(phigh)';
end

%% now plot the results for tcrit
defval = defvaltcrit;
loval = lovaltcrit;
hival = hivaltcrit;

figure;
sq_diffs_low = (defval-loval).^2;
sq_diffs_high = (hival - defval).^2;
diffs_lows = sum(sq_diffs_low,1);
diffs_high = sum(sq_diffs_high,1);
diffs = diffs_high + diffs_lows ; 
[~,jsort] = sort(diffs);
def_vals_zeros = zeros(1,Npar);
for j = 1:Npar
	%plot([loval(jsort(j)) hival(jsort(j))],j*[1 1],'b-o','linewidth',4,'MarkerFaceColor','b');
	plot([-diffs_lows(jsort(j)) diffs_high(jsort(j))],j*[1 1],'b-o','linewidth',4,'MarkerFaceColor','b');
    plot(-diffs_lows(jsort(j)),j,'bo','MarkerFaceColor','g');
	plot(diffs_high(jsort(j)),j,'bo','MarkerFaceColor','r');
	plot(def_vals_zeros,j,'ko','MarkerFaceColor','k','MarkerSize',5);
	hold on;
	if diffs_lows(jsort(j)) < diffs_high(jsort(j))
		loalign = 'right'; highalign = 'left';
	else
		loalign = 'left';highalign = 'right';
	end
	%text(defval,j,num2str(par{paramstochange(jsort(j)),'value'}),'horizontalalignment','center','verticalalignment','bottom');
	%text(diffs_lows(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j))),'low') '  '],'HorizontalAlignment',loalign,'color','g', 'FontSize',12);
	%text(diffs_high(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j))),'high') '  '],'HorizontalAlignment',highalign,'color','r', 'FontSize',12);
    %text(diffs_high(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j)))) '  '],'HorizontalAlignment','right','color','k', 'FontSize',12);
end

hold on
plot(0*[1 1],[0 Npar+0.5],'k-');
%set(gca,'Ytick',[1:Npar],'Yticklabel',par(jsort), 'FontSize', 15);
set(gca,'Ytick',[1:Npar],'Yticklabel',char(paramnames(jsort)), 'FontSize', 15);
%set(gca,'Xlim',[-20 120],'Ylim',[0.5 Npar+0.5]);
xlabel('t_{crit}');
grid on
title('Local Sensitivity in t_{crit}')
set(gca,'FontSize',16)
%% Now for phisc
defval = defvalphisc;
loval = lovalphisc;
hival = hivalphisc;
figure;
sq_diffs_low = (defval-loval).^2;
sq_diffs_high = (hival - defval).^2;
diffs_lows = sum(sq_diffs_low,1);
diffs_high = sum(sq_diffs_high,1);
diffs = diffs_high + diffs_lows ; 
[~,jsort] = sort(diffs);
def_vals_zeros = zeros(1,Npar);
for j = 1:Npar
	%plot([loval(jsort(j)) hival(jsort(j))],j*[1 1],'b-o','linewidth',4,'MarkerFaceColor','b');
	plot([-diffs_lows(jsort(j)) diffs_high(jsort(j))],j*[1 1],'b-o','linewidth',4,'MarkerFaceColor','b');
    plot(-diffs_lows(jsort(j)),j,'bo','MarkerFaceColor','g');
	plot(diffs_high(jsort(j)),j,'bo','MarkerFaceColor','r');
	plot(def_vals_zeros,j,'ko','MarkerFaceColor','k','MarkerSize',5);
	hold on;
	if diffs_lows(jsort(j)) < diffs_high(jsort(j))
		loalign = 'right'; highalign = 'left';
	else
		loalign = 'left';highalign = 'right';
	end
	%text(defval,j,num2str(par{paramstochange(jsort(j)),'value'}),'horizontalalignment','center','verticalalignment','bottom');
	%text(diffs_lows(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j))),'low') '  '],'HorizontalAlignment',loalign,'color','g', 'FontSize',12);
	%text(diffs_high(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j))),'high') '  '],'HorizontalAlignment',highalign,'color','r', 'FontSize',12);
    %text(diffs_high(jsort(j)),j,['  ' num2str(char(paramnames(jsort(j)))) '  '],'HorizontalAlignment','right','color','k', 'FontSize',12);
end

hold on
plot(0*[1 1],[0 Npar+0.5],'k-');
%set(gca,'Ytick',[1:Npar],'Yticklabel',par(jsort), 'FontSize', 15);
set(gca,'Ytick',[1:Npar],'Yticklabel',char(paramnames(jsort)), 'FontSize', 15);
%set(gca,'Xlim',[-20 120],'Ylim',[0.5 Npar+0.5]);
xlabel('\phi_{s} @ 10 weeks');
grid on
title('Local Sensitivity in \phi_{s}')
set(gca,'FontSize',16)



%% Next step is to use the modelfuns to perform a global sensitivity analysis
% Using the procedure outlined in Angela's slides
% We will separately use modelfuntcrit and modelfunphisc to generate our
% model output trajectoriies

N=5e3;
% 1. Generate two sets of random numbers on (0,1)
Xmat = rand([Npar N]);
Zmat = rand([Npar N]);
% 2. Convert each row to match the distributions for each parameter...
% I am a little confused what this means, but I think I can likely assume
% that the distribution on parameters is uniformly distributed between soem
% reasonable upper and lower bounds (that I set)
%paramnames = {'\phi_{s}(0)', 'K1', 'K2', 'r_{s}', '\alpha', '(rr/rs)', 'd_{s}', '(dr/ds)'};
% Set upper and lower bounds of these parameters and then sample from it
pbounds(1,:)= [0,1]; % phis0
pbounds(2,:) = [3e4, 7e4]; % carrying capacity of cells in a 96 well
pbounds(3,:) = [15e6, 25e6]; %carrying capacity of cells in a 10 cm dish
pbounds(4,:) = [0.5*gtot, 1.5*gtot]; % sensitive cell growth rate per hour
pbounds(5,:) = [0, 0.5]; % alpha (degree of drug induced resistance)
pbounds(6,:) = [0,1.5]; % resistant to sensitive cell growth rate ratio
pbounds(7,:) = [0.01, 0.1]; % ds
pbounds(8,:)= [0, 1]; % %resistant to sensitive cell death rate ratio 
% How to sample from uniform distribution between upper and lower values
for i = 1:N
Xmatp(:,i) = pbounds(:,1) + Xmat(:,i).*(pbounds(:,2)-pbounds(:,1));
Zmatp(:,i) = pbounds(:,1) + Zmat(:,i).*(pbounds(:,2)-pbounds(:,1));
end
% Run N simulations to calculate the outputs' mean and variance with first
% set of samples (Xmatp)
for i=1:N
Omattcrit(:,i) = modelfuntcrit(Xmatp(:,i)');
Omatphi(:,i) = modelfunphisc(Xmatp(:,i)');
end
%% Calculate mean across the columns (i.e. all parameter sets)
mutcrit = mean(Omattcrit,2); 
sigmasqtcrit = std(Omattcrit,1, 2).^2;
vartcrit = var(Omattcrit,1,2);
sigmasq = (sum(Omattcrit.^2,2)./N) - mutcrit.^2; % check that Matlab's variance calculation is same

muphi = mean(Omatphi,2); 
sigmasqphi = std(Omatphi,1, 2).^2;
sigmasqphitest = (sum(Omatphi.^2,2)./N) - muphi.^2;
varphi = var(Omatphi,1,2);
%% Check that the three ways you calculate variance are the same
figure;
plot(Cdoxvec, mutcrit, 'b','LineWidth', 3)
hold on
plot(Cdoxvec, mutcrit-sigmasqtcrit, 'b--','LineWidth', 2)
plot(Cdoxvec, mutcrit+sigmasqtcrit, 'b--','LineWidth', 2)
plot(Cdoxvec, mutcrit-sigmasq, 'r--','LineWidth', 2)
plot(Cdoxvec, mutcrit+sigmasq, 'r--','LineWidth', 2)
plot(Cdoxvec, mutcrit-vartcrit, 'm--','LineWidth', 3)
plot(Cdoxvec, mutcrit+vartcrit, 'm--','LineWidth', 3)
ylabel('\mu & \sigma^{2}')
title('t_{crit} over Cdox range')
set(gca,'FontSize',20,'LineWidth',1.5)

figure;
plot(Cdoxvec, muphi, 'b','LineWidth', 3)
hold on
plot(Cdoxvec, muphi-sigmasqphi, 'b--','LineWidth', 2)
plot(Cdoxvec, muphi+sigmasqphi, 'b--','LineWidth', 2)
plot(Cdoxvec, muphi-sigmasqphitest, 'r--','LineWidth', 2)
plot(Cdoxvec, muphi+sigmasqphitest, 'r--','LineWidth', 2)
plot(Cdoxvec, muphi-varphi, 'm--','LineWidth', 3)
plot(Cdoxvec, muphi+varphi, 'm--','LineWidth', 3)
ylabel('\mu & \sigma^{2}')
title('\phi over Cdox range')
set(gca,'FontSize',20,'LineWidth',1.5)
%% Make appropriate rearrangements of the two random sets to calculate the
% for the each parameter
% The main effects matrix will have xs from the row corresponding to itself
% and z in all others, and the total effects matrix will have zs from the
% row corresponding to itself and x in all others.
for i = 1:Npar
    MEmat(:,:)=Zmatp; % start by filling the main effects mat with zs
    MEmat(i,:) = Xmatp(i,:); % replace the kth row with a row of xs
    TEmat(:,:) = Xmatp; % same thing: start by filling with Xs
    TEmat(i,:) = Zmatp(i,:); % repalce kth row with a row of zs
    
    % Run 2N simulations and calculate main & total effects for each
    % parameter by running forward model for each column in the MEmat and
    % TEmat
    for j=1:N
    OMtcrit(:,j) = modelfuntcrit(MEmat(:,j)');
    OTtcrit(:,j) = modelfuntcrit(TEmat(:,j)');
    OMphi(:,j) = modelfunphisc(MEmat(:,j)');
    OTphi(:,j) = modelfunphisc(TEmat(:,j)');
    end
    
    % Now you have your Main and Total effects outputs, find the SI
    for j = 1:N
        mat_to_sumMt(:,j) = Omattcrit(:,j).*OMtcrit(:,j);
        mat_to_sumTt(:,j) = (Omattcrit(:,j) - OTtcrit(:,j)).^2;
        mat_to_sumMp(:,j) = Omatphi(:,j).*OMphi(:,j);
        mat_to_sumTp(:,j) = (Omatphi(:,j)-OTphi(:,j)).^2;
    end
    
    % Using slide 153 of definition of main effects
        MEtcrit(:,i)=((sum(mat_to_sumMt,2)./N)-(mutcrit.^2))./vartcrit;
        TEtcrit(:,i) = sum(mat_to_sumTt,2)./(2*N.*vartcrit);
         MEphi(:,i)=((sum(mat_to_sumMp,2)./N)-(muphi.^2))./varphi;
        TEphi(:,i) = sum(mat_to_sumTp,2)./(2*N.*varphi);
    
% These should be the same!

% Compare this to the easy version (i.e. just find the variance of the
% Total and Main effect matrices (p x N)that you generated where p =
% length(Cdoxvec)
% This might need to be sample variance
    MEcheckt(:,i) = (var(OMtcrit,1,2));
    TEcheckt(:,i) = (var(OTtcrit,1,2));
    MEcheckp(:,i) = (var(OMphi,1,2))./varphi;
    TEcheckp(:,i) = (var(OTphi,1,2))./varphi;
end 
%% Make bar graphs of the total effects on for each dose

figure;
bar(TEtcrit')
hold on
plot([0 9], [0.05 0.05], 'r--')
%ylim([0 1])
%xlim([0 9])
set(gca,'FontSize',20,'LineWidth',1.5)
set(gca, 'Xtick', [1:1:8])
set(gca,'XTickLabel', paramnames)
legend('0 nM', '75 nM', '200 nM', '500 nM','threshold','FontSize', 12, 'Location', 'NorthWest')
legend boxoff
title('Total Effects in t_{crit}')
%%
figure;
bar(TEphi')
hold on
plot([0 9], [0.05 0.05], 'r--')
ylim([0 1])
xlim([0 9])
set(gca,'FontSize',20,'LineWidth',1.5)
set(gca, 'Xtick', [1:1:8])
set(gca,'XTickLabel', paramnames)
legend('0 nM', '75 nM', '200 nM', '500 nM','threshold','FontSize', 12, 'Location', 'NorthWest')
legend boxoff
title('Total Effects in \phi_{s} at t_{crit}')
%% Main effects bar graph
figure;
bar(MEtcrit')
hold on
plot([0 9], [0.05 0.05], 'r--')
ylim([-0.3 1])
xlim([0 9])
set(gca,'FontSize',20,'LineWidth',1.5)
set(gca, 'Xtick', [1:1:8])
set(gca,'XTickLabel', paramnames)
legend('0 nM', '75 nM', '200 nM', '500 nM', 'threshold', 'FontSize', 12, 'Location', 'NorthWest')
legend boxoff
title('Main Effects in t_{crit}')

figure;
bar(MEphi')
hold on
plot([0 9], [0.05 0.05], 'r--')
ylim([-0.3 1])
xlim([0 9])
set(gca,'FontSize',20,'LineWidth',1.5)
set(gca, 'Xtick', [1:1:8])
set(gca,'XTickLabel', paramnames)
legend('0 nM', '75 nM', '200 nM', '500 nM','threshold', 'FontSize', 12, 'Location', 'NorthWest')
legend boxoff
title('Main Effects in \phi_{s} at t_{crit}')
%% Additional effects bar graph
AEtcrit = TEtcrit-MEtcrit;
AEphi = TEphi-MEphi;

figure;
bar(AEtcrit')
hold on
plot([0 9], [0.05 0.05], 'r--')
ylim([0 1])
xlim([0 9])
set(gca,'FontSize',20,'LineWidth',1.5)
set(gca, 'Xtick', [1:1:8])
set(gca,'XTickLabel', paramnames)
legend('0 nM', '75 nM', '200 nM', '500 nM','threshold', 'FontSize', 12, 'Location', 'NorthWest')
legend boxoff
title('Additional Effects in t_{crit}')

figure;
bar(AEphi')
hold on
plot([0 9], [0.05 0.05], 'r--')
ylim([0 1])
xlim([0 9])
set(gca,'FontSize',20,'LineWidth',1.5)
set(gca, 'Xtick', [1:1:8])
set(gca,'XTickLabel', paramnames)
legend('0 nM', '75 nM', '200 nM', '500 nM','threshold', 'FontSize', 12, 'Location', 'NorthWest')
legend boxoff
title('Additional Effects in \phi_{s} at t_{crit}')

 %% Plot the SI of the main effects and the total effects for each parameter

 
 figure;
 for i = 1
 plot(Cdoxvec, MEtcrit(:,i), 'b-', 'LineWidth', 2)
     hold on
plot(Cdoxvec, MEcheckt(:,i), 'r--', 'LineWidth', 2)
plot(Cdoxvec, MEcheckt(:,i)-MEtcrit(:,i), 'm--')
title('Main Effects in tcrit')
 end
 
  figure;
 for i = 1
 plot(Cdoxvec, TEtcrit(:,i), 'b-', 'LineWidth', 2)
     hold on
plot(Cdoxvec, TEcheckt(:,i), 'r--', 'LineWidth', 2)
 end
 title('Total effects in tcrit')
  
 figure;
 for i = 1
 plot(Cdoxvec, MEphi(:,i), 'b-', 'LineWidth', 2)
     hold on
plot(Cdoxvec, MEcheckp(:,i), 'r--', 'LineWidth', 2)
title('Main Effects in phi')
 end
 
  figure;
 for i = 1
 plot(Cdoxvec, TEphi(:,i), 'b-', 'LineWidth', 2)
     hold on
plot(Cdoxvec, TEcheckp(:,i), 'r--', 'LineWidth', 2)
title('Main Effects in phi')
 end
 %%
 
 figure;
 subplot(1,2,1)
 for i = 1:Npar
 hold on
 plot(Cdoxvec, MEtcrit(:,i), '-', 'LineWidth', 2)
 end
xlabel ('C_{dox}')
ylabel('Main Effects on t_{crit}')
title('T_{crit} Main Effects')
legend(paramnames)
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
subplot(1,2,2)
for i = 1:Npar
 plot(Cdoxvec, TEtcrit(:,i), '-', 'LineWidth', 2)
 hold on
end
xlabel ('C_{dox}')
ylabel('Total Effects on t_{crit}')
title('T_{crit} Total Effects')
legend(paramnames)
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
 %%
  figure;
     subplot(1,2,1)
     for i = 1:Npar
     plot(Cdoxvec, MEphi(:,i), '-', 'LineWidth', 2)
     hold on
     end
    xlabel ('C_{dox}')
    ylabel('Main Effects on \phi_{s} at t_{crit}')
    title('\phi_{s} Main Effects')
    legend(paramnames)
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
    subplot(1,2,2)
    for i = 1:Npar
     plot(Cdoxvec, TEphi(:,i), '-', 'LineWidth', 2)
     hold on
    end
    xlabel ('C_{dox}')
    ylabel('Total Effects on \phi_{s} at t_{crit}')
    title('\phi_{s} Total Effects')
    legend(paramnames)
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
 %%   
    figure;
 subplot(1,2,1)
 for i = 1:Npar
 hold on
 plot(Cdoxvec, MEcheckt(:,i), '-', 'LineWidth', 2)
 end
xlabel ('C_{dox}')
ylabel('Main Effects on t_{crit}')
title('T_{crit} Main Effects')
legend(paramnames)
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
subplot(1,2,2)
for i = 1:Npar
 plot(Cdoxvec, TEcheckt(:,i), '-', 'LineWidth', 2)
 hold on
end
xlabel ('C_{dox}')
ylabel('Total Effects on t_{crit}')
title('T_{crit} Total Effects')
legend(paramnames)
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
 %%
  figure;
     subplot(1,2,1)
     for i = 1:Npar
     plot(Cdoxvec, MEcheckp(:,i), '-', 'LineWidth', 2)
     hold on
     end
    xlabel ('C_{dox}')
    ylabel('Main Effects on \phi_{s} at t_{crit}')
    title('\phi_{s} Main Effects')
    legend(paramnames)
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
    subplot(1,2,2)
    for i = 1:Npar
     plot(Cdoxvec, TEcheckp(:,i), '-', 'LineWidth', 2)
     hold on
    end
    xlabel ('C_{dox}')
    ylabel('Total Effects on \phi_{s} at t_{crit}')
    title('\phi_{s} Total Effects')
    legend(paramnames)
    legend boxoff
    set(gca,'FontSize',20,'LineWidth',1.5)
