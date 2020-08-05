% Plot the dose-response of the resistant and sensitive lineages to the dox
% treatment

close all; clear all; clc

%import data
[N, T] = xlsread('../data/recall_dose_response.xlsx');
% ic25 = 400 nM
% ic75 = 2.5 uM
figure;
plot(N(1:12, 1), N(1:12, 3), 'm*', 'LineWidth',4)
hold on
plot(N(13:end, 1), N(13:end, 3), 'g*', 'LineWidth', 4)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel(' [Dox] for 48h')
ylabel('% Viability after 48h')
legend('resistant lineage AA161', 'sensitive lineage AA170')
legend boxoff

figure;
errorbar(N(1, 1), mean(N(1:6, 3)), 1.96*std(N(1:6,3))./2, 'm*', 'LineWidth',4)
hold on
errorbar(N(1, 1),mean(N(13:18, 3)), 1.96*std(N(13:18,3))./2, 'g*', 'LineWidth', 4)
errorbar(N(7, 1), mean(N(7:12, 3)), 1.96*std(N(7:12,3))./2, 'm*', 'LineWidth',4)
errorbar(N(7, 1), mean(N(19:24, 3)), 1.96*std(N(19:24,3))./2, 'g*', 'LineWidth',4)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel(' [Dox] for 48h')
ylabel('% Viability after 48h')
legend('resistant lineage AA161', 'sensitive lineage AA170')
legend boxoff
xlim([ 0 100])