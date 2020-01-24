% Functional growth rate and drug sensitivity analysis

close all; clear all; clc
%%  Load in the 231 AA161 and AA170 growth dynamics (N(t)) data from Didi

[N, T] =xlsread('../data/Recall_growth_analysis.xlsx');
[N1, T1] = xlsread('../data/recall_dose_response.xlsx');
K = 5e4;

for i = 1:size(N,2)-1 %each data set runs
        recall(i).time = N(:,1);
    recall(i).conf = N(:, 1+i);
    recall(i).Nest = K.*(N(:,1+i)./100);
end
    for i = 1:12
        recall(i).lin = 'AA161';
        recall(i).lintype = 'res';
    end
    for i = 13:24
        recall(i).lin = 'AA170';
        recall(i).lintype = 'sens';
    end


%% Plot the raw data
figure;
hold on
plot(recall(1).time, recall(1). Nest, 'm-')
plot(recall(13).time, recall(13).Nest, 'g-')
for i = 1:12
    plot(recall(i).time, recall(i).Nest, 'm-')
    hold on
end
for i = 13:24
    plot(recall(i).time, recall(i).Nest, 'g-')
end
legend('AA161 resistant lineage', 'AA170 sensitive lineage', 'Location', 'Northwest')
legend boxoff
xlabel ('time (hours)')
ylabel('N(t)/N_{0}')
xlim([0 300])
set(gca,'FontSize',20,'LineWidth',1.5)
%% Fit each trajectory to logistic growth model

for i = 1:length(recall)
    sigma =1;
    K = 5e4;
     [g, singexpmodel] = fit_logistic(recall(i).Nest,recall(i).time, sigma, K)
     recall(i).g = g;
     recall(i).Nmod = singexpmodel;
end
figure;
hold on
plot(recall(1).time, recall(1). Nest, 'r.')
plot(recall(13).time, recall(13).Nest, 'g.')
for i = 1:12
    plot(recall(i).time, recall(i).Nest, 'r.')
    plot(recall(i).time, recall(i).Nmod, 'r-')
    
    hold on
end
for i = 13:24
    plot(recall(i).time, recall(i).Nest, 'g.')
    plot(recall(i).time, recall(i).Nmod, 'g-')
end
legend('AA161 resistant lineage', 'AA170 sensitive lineage', 'Location', 'Northwest')
legend boxoff
xlabel ('time (hours)')
ylabel('N(t)')

set(gca,'FontSize',20,'LineWidth',1.5)
%% Compile into summary structure
for i = 1:2
    recallsum(i).time = N(:,1);
end
recallsum(1).Nmat = [];
recallsum(1).gvec = [];
for i = 1
    recallsum(i).lin = 'AA1161';
    recallsum(i).lintype = 'res';
    for j = 1:12
    recallsum(i).Nmat =horzcat(recallsum(1).Nmat,recall(j).Nest);
    recallsum(i).gvec = horzcat(recallsum(1).gvec, recall(j).g);
    end
end
recallsum(2).Nmat = [];
recallsum(2).gvec = [];
for i = 2
    recallsum(i).lin = 'AA170';
    recallsum(i).lintype = 'sens';
    for j = 13:24
    recallsum(i).Nmat =horzcat(recallsum(2).Nmat,recall(j).Nest);
    recallsum(i).gvec = horzcat(recallsum(2).gvec, recall(j).g);
    end
end
for i = 1:2
    recallsum(i).Nmean = mean(recallsum(i).Nmat,2);
    recallsum(i).Nmean = recallsum(i).Nmean(1:50)
    recallsum(i).Nstd = std(recallsum(i).Nmat,0,2);
    recallsum(i).Nstd = recallsum(i).Nstd(1:50)
    sigma =recallsum(i).Nstd;
    K = 5e4;
     [g, singexpmodel] = fit_logistic(recallsum(i).Nmean,recallsum(i).time(1:50), sigma, K)
     recallsum(i).g = g;
     recallsum(i).Nmod = singexpmodel;
     recallsum(i).gmean = mean(recallsum(i).gvec);
     recallsum(i).gstd = std(recallsum(i).gvec);
end

figure;
errorbar(recallsum(1).time(1:50), recallsum(1).Nmean, 1.96*recallsum(1).Nstd, 'm')
hold on
errorbar(recallsum(2).time(1:50), recallsum(2).Nmean, 1.96*recallsum(2).Nstd, 'g')
plot(recallsum(1).time(1:50), recallsum(1).Nmod, 'm-', 'LineWidth',2)
text(recallsum(1).time(30), recallsum(1).Nmod(end-20), ['g_{res}= ',num2str(round(recallsum(1).g,3)), 'per hour'], 'FontSize', 16)
text(recallsum(2).time(30), recallsum(2).Nmod(end-20), ['g_{sens}= ',num2str(round(recallsum(2).g,3)), 'per hour'], 'FontSize', 16)
hold on
plot(recallsum(2).time(1:50), recallsum(2).Nmod, 'g-', 'LineWidth',2)
legend('AA161 resistant lineage', 'AA170 sensitive lineage', 'Location', 'Northwest')
legend boxoff
xlabel('time (hours)')
ylabel('N(t)')
xlim([0 190])
set(gca,'FontSize',20,'LineWidth',1.5)
%% Add in the dosed data
for i = 1:12
    recallsum(1).dosevec = N1(1:12, 1);
    recallsum(1).viabvec = N1(1:12,3);
end
for i = 13:24
    recallsum(2).dosevec = N1(13:24, 1);
    recallsum(2).viabvec = N1(13:24,3);
end

for i = 1:2
    recallsum(i).viab = [mean(recallsum(i).viabvec(1:6)), mean(recallsum(i).viabvec(7:12))];
    recallsum(i).stdviab = [std(recallsum(i).viabvec(1:6)), std(recallsum(i).viabvec(7:12))];
    
end
%% Run a t-test to compare groups
[h1,pgrowth,ci1,stats] = ttest2(recallsum(1).gvec, recallsum(2).gvec);
[h2,pdrug1,ci2,stats] = ttest2(recallsum(1).viabvec(1:6), recallsum(2).viabvec(1:6));
[h3,pdrug2,ci3,stats] = ttest2(recallsum(1).viabvec(7:12), recallsum(2).viabvec(7:12));

%% Make a bar graph of growth rate
linnames = {'AA161', 'AA170'}
gvec = [recallsum(1).gmean, recallsum(2).gmean];
gerr = (1.96./2).*[recallsum(1).gstd, recallsum(2).gstd];

figure;
% subplot(2, 2, 1:2)
% errorbar(recallsum(1).time(1:50), recallsum(1).Nmean, 1.96*recallsum(1).Nstd, 'm')
% hold on
% errorbar(recallsum(2).time(1:50), recallsum(2).Nmean, 1.96*recallsum(2).Nstd, 'g')
% plot(recallsum(1).time(1:50), recallsum(1).Nmod, 'm-', 'LineWidth',2)
% text(recallsum(1).time(30), recallsum(1).Nmod(end-20), ['g_{res}= ',num2str(round(recallsum(1).g,3)), ' cells/hr'], 'FontSize', 16)
% text(recallsum(2).time(30), recallsum(2).Nmod(end-20), ['g_{sens}= ',num2str(round(recallsum(2).g,3)), ' cells/hr'], 'FontSize', 16)
% hold on
% plot(recallsum(2).time(1:50), recallsum(2).Nmod, 'g-', 'LineWidth',2)
% legend('AA161 resistant lineage', 'AA170 sensitive lineage', 'Location', 'Northwest')
% legend boxoff
% xlabel('time (hours)')
% ylabel('N(t)')
% xlim([0 190])
% set(gca,'FontSize',20,'LineWidth',1.5)

subplot(1,2,1)
for i = 1:length(gvec)
    h=bar(i, gvec(i));
    hold on
    er = errorbar(i,gvec(i),gerr(i), gerr(i), 'LineWidth', 2);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    hold on
    if i == 1
        set(h,'FaceColor','m');
    elseif i==2 
        set(h,'FaceColor','g');
    end
end
plot([1,2], [1 1]*0.016, '-k', 'LineWidth',2)
plot([1.4, 1.5, 1.6], [1 1 1]*0.016*1.15, '*k')
set(gca, 'Xtick', [1:1:2])
%ylim([ 0 0.018])
set(gca,'XTickLabel',linnames )
set(gca,'FontSize',20,'LineWidth',1.5)
ylabel ('growth rate')

viabvec = [recallsum(1).viab(1), recallsum(2).viab(1) recallsum(1).viab(2), recallsum(2).viab(2)];
viaberr = (1.96./2).*[ recallsum(1).stdviab(1), recallsum(2).stdviab(1) recallsum(1).stdviab(2), recallsum(2).stdviab(2)];

dosename = ['LD25 (400 nM)','', '', 'LD75 (2.5uM)', ''];
subplot(1,2,2)
for i = 1:4
    if i ==3 || i ==4
    h=bar(i+1, viabvec(i));
    hold on
    er = errorbar(i+1,viabvec(i),viaberr(i), viaberr(i), 'LineWidth', 2);
    else
     h=bar(i, viabvec(i));
    hold on
    er = errorbar(i,viabvec(i),viaberr(i), viaberr(i), 'LineWidth', 2);   
    end
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    hold on
    if i == 1 || i == 3
        set(h,'FaceColor','m');
    elseif i==2  || i == 4
        set(h,'FaceColor','g');
    end
end
plot([4,5], [1 1]*75, '-k', 'LineWidth',2)
plot(4.5, [1 1]*75*1.15, '*k')

set(gca,'XTickLabel',[] )
set(gca,'FontSize',20,'LineWidth',1.5)
ylabel ('% viability @ 48h')
ylim ([0 100])
%% Run a t-test to compare groups

