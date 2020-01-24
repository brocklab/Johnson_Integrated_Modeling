% Make dose history figure

% This script will be used to make visualizations of the dosing history for
% each dose response curve in the data set, using the information in traj.
% mat.

close all; clear all; clc;
%% Load in data structure 
S = load('../out/trajraw.mat');
traj= S.traj;
%% pick wells (rows in traj) to graph
% keeps track of the "index" in the subplot
count = 0;
% active represents the list of i values for the wells being plotted
active = [121, 127, 133, 139, 145, 151, 157, 163, 169, 175, 1081, 1093, 1105,1117,1129, 1141, 1153, 1165, 1177, 1189, ...
        1249, 1237, 1225, 1213, 1201];
% checks every i in traj
for i = 1:length(traj)
    % checks if i is one to be graphed
    if ismember(i, active)
        
        % advances index in subplot
        count = count + 1;
        
        % mess with count manually here
        if i == 1081
             count = count + 1;
        elseif i == 1093
            count = count + 5;
        elseif i == 1129
            count = count - 10;
        elseif i == 1141
            count = count + 1;
        elseif i == 1201
            count = count + 11;
        end
        if i > 1201
            count = count - 3;
        end
        
        %% time and dosage vectors to plot
        n = traj(i).numdoses;

        % two ways to set time scale (x-axis): auto scale or manual set
        % comment out one or the other
        
            % Below: tvec truncated after last dose
            % 2 days before first dose, 2 days after last dose
            % tvec = 0:1:(sum(traj(i).doseintdays) + 4 + (3 * n));

            % Below: tvec equal to constant
            if i < 1200
                tvec = 0:25;
            else
                tvec = 0:75;
            end
        % equalize vector length of y values
        uvec = 0.*tvec;

        % see "dosegraphindex.pdf" in documentation to understand this index jumping/modifying
        ind = 4;
        for j = 1:n
            % modify tvec
            for k = ind:length(tvec)
                tvec(k) = tvec(k) - 1;
            end
            for k = (ind + 2):length(tvec)
                tvec(k) = tvec(k) - 1;
            end
            % modify uvec
            % final dose
            if j == n  
                uvec(ind) = traj(i).dose;
                uvec(ind + 1) = traj(i).dose;
            % not final dose
            else
                uvec(ind) = traj(i).prevdose(j);
                uvec(ind + 1) = traj(i).prevdose(j);
                ind = ind + 3 + traj(i).doseintdays(j);
            end
        end
        %% plot and format
        % first two numbers are dimensions of overall plot: count =
        % ordering index
        if i < 1200
            subplot(3,10,count);
        else
            subplot(3,10, [count (count + 1)]);
        end
        % generates plots
        area(tvec,uvec);
        
        %sets y axis
        ylim([0 300]);
        
        % individual labeling 
        %{
        xlabel('time (days)')
        ylabel('Dose')
        title('Treatment History Example')
        %}
        
        % formatting options
        set(gca,'FontSize',10,'LineWidth',1.5,'XTickLabel',[],'YTickLabel',[])
    end
end

%% following block is Kaitlyn's old time and dosage vectors
%{
% Start my making a time vector that spans the pre-treatment history + two
% weeks into the current treatment. So here that will be
tvec = 0:1:(7*(sum(traj(i).doseints)+2));

% Next make your dose vec. Start with a vector of zeros the length of your
% time vector
Uvec = 0.*tvec;
% Add in pulse treatments at the corresponding days
% First pulse treatment starts at 0
ind1 = find(ismember(tvec, 0:traj(i).doseduration/24))
Uvec(ind1) = traj(i).prevdose(1)
%% Second pulse treatment start at variable places given by the 

for j = 1:length(traj(i).doseints)
    ind = find(ismember(tvec, (traj(i).doseints(j)*7):((traj(i).doseints(j)*7)+traj(i).doseduration/24)));
    Uvec(ind)= traj(i).prevdose(j);
end
%}