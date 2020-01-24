% This script analyzes the MCF-7 data-- lot that can be done here
% So far just wanted to plot some things

% THIS IS ANOTHER TEST

close all; clear all; clc;
%% Load in data structure 
S = load('../out/trajfit.mat');
traj= S.traj;
%% Plot some things that might be interesting
% Just playing around here

figure;
subplot(1,3,1)
for i = 1:length(traj)
    if ~isempty(traj(i).WPT) && traj(i).bfmod == 1
    plot(traj(i).WPT, traj(i).phi, 'r*')
    hold on
    xlabel('weeks post prior treatment')
    ylabel ('resistant fraction (\phi)')
    title('Resistant fraction versus WPT')
    end
end
subplot(1,3,2)
for i =1:length(traj)
    if ~isempty(traj(i).WPT) && traj(i).bfmod == 1
    plot(traj(i).WPT, traj(i).k, 'b*')
    hold on
    xlabel('weeks post prior treatment')
    ylabel ('death rate (k)')
    title('Death rate versus WPT')
    end
end

subplot(1,3,3)
for i =1:length(traj)
    if ~isempty(traj(i).WPT) && traj(i).bfmod == 1
    plot(traj(i).WPT, traj(i).g, 'g*')
    hold on
    xlabel('weeks post prior treatment')
    ylabel ('resistant regrowth rate (g)')
    title('Resistant regrowth rate versus WPT')
    end
end

%% Plot comparing number of doses
% Just playing around here

figure;
subplot(1,3,1)
for i = 1:length(traj)
    if ~isempty(traj(i).numdoses) && traj(i).bfmod == 1
    plot(traj(i).numdoses, traj(i).phi, 'r*')
    hold on
    xlabel('number of doses received')
    ylabel ('resistant fraction (\phi)')
    title('Resistant fraction versus number of doses')
    end
end
subplot(1,3,2)
for i =1:length(traj)
    if ~isempty(traj(i).numdoses) && traj(i).bfmod == 1
    plot(traj(i).numdoses, traj(i).k, 'b*')
    hold on
    xlabel('number of doses recevied')
    ylabel ('death rate (k)')
    title('Death rate versus number of doses')
    end
end

subplot(1,3,3)
for i =1:length(traj)
    if ~isempty(traj(i).numdoses) && traj(i).bfmod == 1
    plot(traj(i).numdoses, traj(i).g, 'g*')
    hold on
    xlabel('weeks post prior treatment')
    ylabel ('resistant regrowth rate (g)')
    title('Resistant growth rate versus number of doses')
    end
end

%% Plot tcrit versus number of doses
figure;
for i = 1:length(traj)
    hold on
    if ~isempty(traj(i).tcrit)
    if traj(i).dose == 75 && traj(i).dosenum ==1
        plot(traj(i).accdose, traj(i).tcrit, 'k*', 'LineWidth',3)
    end
    if traj(i).dose == 75 && traj(i).dosenum == 2 && traj(i).WPT ==2
        plot(traj(i).accdose, traj(i).tcrit, 'k*', 'LineWidth',3)
    end
    if traj(i).dose == 75 && traj(i).dosenum == 3 && traj(i).WPT == 2
        if traj(i).doseints(1)==2 
        plot(traj(i).accdose, traj(i).tcrit, 'k*', 'LineWidth', 3)
        end
    end
    end
end
title('T_{crit} for 75 nM pulse every two weeks', 'FontSize', 20)
xlabel('Accumulated Dose (nM)', 'FontSize', 20)
ylabel('T_{crit} (N>1.2*N_{treat})', 'FontSize',20)