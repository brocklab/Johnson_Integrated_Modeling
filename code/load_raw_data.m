% This script loads in data sets from Grant's various experimental
% conditions of drug response. The goal of this script is just to capture
% the input variables and data for each well in the data set.
close all; clear all;
% load in adjusted data set (remove top part and record relevant info
% below)
[raw(1).N, raw(1).T] =xlsread('../data/GH1831_MCF7_PredModel_2nd_Dose_Varying_Interval_v2.xls');
[raw(2).N, raw(2).T] =xlsread('../data/GH1830_MCF7_PredModel_3rd_Dose_Vary_Interval_v2.xls');
[raw(3).N, raw(3).T] =xlsread('../data/GH1818_MCF7_PredModel_1st_Dose_Analysis_v6.xls');
[raw(4).N, raw(4).T] =xlsread('../data/GH1832_MCF7_PredModel_3rd_Dose_4-4_Intervals_v1.xls');
[raw(5).N, raw(5).T] =xlsread('../data/GH1820_MCF7_PredModel_2nd_Dose_Varying_Interval_v3.xls');
[raw(6).N, raw(6).T] =xlsread('../data/GH1829_MCF7_PredModel_2nd_Dose_4_Week_Interval_v1.xls');
[raw(7).N, raw(7).T] =xlsread('../data/GH1819_MCF7_PredModel_3rd_Treatment_Analysis_v2.xls');
[raw(8).N, raw(8).T] =xlsread('../data/GH1817_MCF7_PredModel_2nd_Dose_Analysis_v3.xls');
[raw(9).N, raw(9).T] =xlsread('../data/GH1815_MCF7_PredModel_Effect_of_1st_Dose_Analysis_v2.xls');
[raw(10).N, raw(10).T] =xlsread('../data/GH1835_BT474_PredModel_2nd_Dose_Vary_Interval_v1.xls');
[raw(11).N, raw(11).T] =xlsread('../data/GH1902_MCF7_PredModel_4th_Dose_4_Week_Intervals_v1.xls');
[raw(12).N, raw(12).T] =xlsread('../data/GH1903_MCF7_PredModel_5th_Dose_4_Week_Intervals_v1.xls');
[raw(13).N, raw(13).T] =xlsread('../data/GH1827_BT474_PredModel_3rd_Dose_2_Week_Int_v1.xls');
[raw(14).N, raw(14).T] =xlsread('../data/GH1828_BT474_PredModel_2nd_Dose_4_WPT1_v1.xls');
[raw(15).N, raw(15).T] =xlsread('../data/GH1825_BT474_PredModel_1st_Dose_0-150_nM_v2.xls');
[raw(16).N, raw(16).T] =xlsread('../data/GH1826_BT474_PredModel_2nd_Dose_2_Week_Int_v1.xls');
[raw(17).N, raw(17).T] =xlsread('../data/GH1904_231_PredModel_1st_Dose_Ex_1_v1.xls');
[raw(18).N, raw(18).T] =xlsread('../data/GH1909_231_PredModel_1st_Dose_Ex_2_v2.xls');
[raw(19).N, raw(19).T] =xlsread('../data/GH1910_MCF7_PredModel_2nd_Dose_Vary_Int_Short_Ex_1_v2.xls');
[raw(20).N, raw(20).T] =xlsread('../data/GH1911_MCF7_PredModel_2nd_Dose_Vary_Int_Short_Ex_2_v2.xls');
[raw(21).N, raw(21).T] =xlsread('../data/GH1912_MCF7_PredModel_Serial_Doses_2_wk_Int_v2.xls');
[raw(22).N, raw(22).T] =xlsread('../data/GH1919_BT474_PredModel_Vary_Int_Ex1_v2.xls');
[raw(23).N, raw(23).T] =xlsread('../data/GH1920_BT474_PredModel_Vary_Int_Ex2_v1.xls');
[raw(24).N, raw(24).T] =xlsread('../data/GH1907_231_PredModel_Vary_Int_Ex_1_v2.xls');
[raw(25).N, raw(25).T] =xlsread('../data/GH1908_231_PredModel_Vary_Int_Ex2_v1.xls');
[raw(26).N, raw(26).T] =xlsread('../data/GH1914_231_PredModel_Vary_Int_Ex3_v1.xls');
[raw(27).N, raw(27).T] =xlsread('../data/GH1918_231_PredModel_Serial_Dosing_2_wk_Int_v2.xls');

% this reads in the numeric values in an excel sheet (N1), and the text
% values (T1)
% Make a "size" matrix that keeps track of the number of wells and time
% points of each data set loaded. This will be used to add to the structure
% once more data sets are loaded in
% 1st column in number of time points, 2nd column number of wells in that
% data set, third column is cumulative number of wells

sz = zeros(length(raw),3);
for i = 1:length(raw)
    sz(i,1:2) = size(raw(i).N);
    sz(i,2) = sz(i,2) - 1;
end
sz(:,3) = cumsum(sz(:,2));
%% Structure using fors and ifs
size = sz(end,3);
traj = struct(  'time', zeros(size,1), 'rawN', zeros(size,1), 'N0true', zeros(size,1), ...
                'date', strings(size,1), 'tdose', zeros(size,1), 'doseduration', zeros(size,1), ...
                'celltype', strings(size,1), 'drug', strings(size,1), 'seed', zeros(size,1), ...
                'welllabel', strings(size,1), 'well', strings(size,1), 'column', strings(size,1), ...
                'dose', zeros(size,1), 'prevdose', zeros(size,1), 'accdose', zeros(size,1), ...
                'dosenum', zeros(size,1), 'numdoses', zeros(size,1), ...
                'doseints', zeros(size,1), 'doseintdays', zeros(size,1), 'WPT', zeros(size,1), ...
                'code', zeros(size,1));
for i = 1:length(raw) %each data set runs
    if i == 1
        x = 1;
    else
        x = sz(i-1,3)+1;
    end
    for j = x:sz(i,3) % set bound in traj for data set
        % find well #, specific to the data set, and load data
        k = j + 1 - x;
        traj(j).time = raw(i).N(1:end,1);
        traj(j).rawN = raw(i).N(1:end,k+1);
        traj(j).N0true = traj(j).rawN(1);
        % load and mess with raw well labels to get actual well
        traj(j).welllabel =  raw(i).T(1,k+1);
        wellhalf = extractAfter(traj(j).welllabel,"(");
        if i < 19 || i == 24 % these have column headers with two parentheses then well #
            wellact = extractAfter(wellhalf,"(");
        else % these only have one set of brackets
            wellact = wellhalf;
        end
        traj(j).well = strtok(wellact, {')'});
        traj(j).column = strtok(traj(j).well, {'B','C','D','E','F','G'});
        % subject to change; currently all data follows these
        traj(j).celltype = 'MCF-7';
        traj(j).drug = 'dox';
        traj(j).doseduration = 24;
        % individual data set parameters
            if i == 1       % GH1831: 2nd dose varying interval
                traj(j).code = 1831;
                traj(j).dose = 75;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).prevdose = 75;
                traj(j).tdose = 60;
                traj(j).date = '12-19-18';
                traj(j).seed = 2000;
                if ismember(traj(j).column, '2')
                    traj(j).WPT = 8;
                    traj(j).doseints = 8;
                elseif ismember(traj(j).column, '3')
                    traj(j).WPT = 7;
                    traj(j).doseints = 7;
                elseif ismember(traj(j).column, '4')
                    traj(j).WPT = 6;
                    traj(j).doseints = 6;
                elseif ismember(traj(j).column, '5')
                    traj(j).WPT = 5;
                    traj(j).doseints = 5;
                elseif ismember(traj(j).column, '6')
                    traj(j).WPT = 4;
                    traj(j).doseints = 4;
                elseif ismember(traj(j).column, '7')
                    traj(j).WPT = 3;
                    traj(j).doseints = 3;
                elseif ismember(traj(j).column, '8')
                    traj(j).WPT = 2;
                    traj(j).doseints = 1;
                elseif ismember(traj(j).column, '9')
                    traj(j).WPT = 1;
                    traj(j).doseints = 1;
                elseif ismember(traj(j).column, '10')
                    traj(j).prevdose = [75, 100, 75, 75, 75]; 
                    traj(j).doseints = [2, 4, 3, 4, 8];
                    traj(j).dosenum = 6;
                    traj(j).numdoses = 6;
                    traj(j).WPT = 8;
                elseif ismember(traj(j).column, '11')
                    traj(j).dosenum = 1;
                    traj(j).numdoses = 1;
                    traj(j).WPT = [];
                    traj(j).prevdose = [];
                    traj(j).doseints = [];
                end
            elseif i == 2   %GH1830; 3rd dose varying interval
                traj(j).code = 1830;
                traj(j).date = '12-19-18';
                traj(j).tdose = 69;
                traj(j).seed = 2000;
                traj(j).dosenum = 3;
                traj(j).numdoses = 3;
                traj(j).prevdose = [75,75];
                traj(j).dose = 75;
                if ismember(traj(j).column, '2')
                    traj(j).doseints = [4,4];
                    traj(j).WPT = 4;
                elseif ismember(traj(j).column, '3')
                    traj(j).doseints = [4,3];
                    traj(j).WPT = 3;
                elseif ismember(traj(j).column, '4')
                    traj(j).doseints = [3,4];
                    traj(j).WPT = 4;
                elseif ismember(traj(j).column, '5')
                    traj(j).doseints = [4,2];
                    traj(j).WPT = 2;
                elseif ismember(traj(j).column, '6')
                    traj(j).doseints = [3,3];
                    traj(j).WPT = 3;
                elseif ismember(traj(j).column, '7')
                    traj(j).doseints = [3,2];
                    traj(j).WPT = 2;
                elseif ismember(traj(j).column, '8')
                    traj(j).doseints = [2,3];
                    traj(j).WPT = 3;
                elseif ismember(traj(j).column, '9')
                    traj(j).doseints = [2,2];
                    traj(j).WPT = 2;
                elseif ismember(traj(j).column, '10')
                    traj(j).doseints = [2,1];
                    traj(j).WPT = 4;
                elseif ismember(traj(j).column, '11')
                    traj(j).dosenum = 1;
                    traj(j).numdoses = 1;
                    traj(j).WPT = [];
                    traj(j).prevdose = [];
                    traj(j).doseints = [];
                end
            elseif i == 3 %GH1818 single dose
                traj(j).code = 1818;
                traj(j).date = '8-16-18';
                traj(j).tdose = 73;
                traj(j).seed = 2000;
                traj(j).dosenum = 1;
                traj(j).prevdose = [];
                traj(j).WPT = [];
                traj(j).numdoses = 1;
                traj(j).doseints = [];
                if ismember(traj(j).column, '2')
                    traj(j).dose = 300;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 150;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 125;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 100;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 75;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 35;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 20;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 10;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end
            elseif i == 4 %GH1832 3rd dose varying dosage
                traj(j).code = 1832;
                traj(j).date = '12-18-18';
                traj(j).tdose = 69;
                traj(j).seed = 2000;
                traj(j).dosenum = 3;
                traj(j).numdoses = 3;
                traj(j).WPT = 4;
                traj(j).prevdose = [75,75];
                traj(j).doseints = [4,4];
                if ismember(traj(j).column, '2')
                    traj(j).dose = 300;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 150;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 125;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 100;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 75;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 35;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 20;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 10;
                elseif ismember(traj(j).column, '11')
                    traj(j).dose = 0;
                end
            elseif i == 5 %GH1820: 2nd/3rd dose varying interval
                traj(j).code = 1820;
                traj(j).date = '9-15-18';
                traj(j).tdose = 63;
                traj(j).seed = 3000;
                traj(j).dose = 75;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).prevdose = 75;
                if ismember(traj(j).column, '2')
                    traj(j).dosenum = 3;
                    traj(j).numdoses = 3;
                    traj(j).WPT = 4;
                    traj(j).prevdose = [75,100];
                    traj(j).doseints = [6,4];
                elseif ismember(traj(j).column, '3')
                    traj(j).dosenum = 3;
                    traj(j).numdoses = 3;
                    traj(j).WPT = 1;
                    traj(j).prevdose = [75,50];
                    traj(j).doseints = [5,1];
                elseif ismember(traj(j).column, '4')
                    traj(j).dosenum = 3;
                    traj(j).numdoses = 3;
                    traj(j).WPT = 1;
                    traj(j).prevdose = [75,35];
                    traj(j).doseints = [4,1];
                elseif ismember(traj(j).column, '5')
                    traj(j).WPT = 6;
                    traj(j).doseints = 6;
                elseif ismember(traj(j).column, '6')
                    traj(j).WPT = 5;
                    traj(j).doseints = 5;
                elseif ismember(traj(j).column, '7')
                    traj(j).WPT = 4;
                    traj(j).doseints = 4;
                elseif ismember(traj(j).column, '8')
                    traj(j).WPT = 3;
                    traj(j).doseints = 3;
                elseif ismember(traj(j).column, '9')
                    traj(j).WPT = 2;
                    traj(j).doseints = 2;
                elseif ismember(traj(j).column, '10')
                    traj(j).WPT = 1;
                    traj(j).doseints = 1;
                elseif ismember(traj(j).column, '11') 
                    traj(j).dosenum = 1;
                    traj(j).numdoses = 1;
                    traj(j).WPT = [];
                    traj(j).prevdose = [];
                    traj(j).doseints = [];
                end    
            elseif i == 6 %GH1829: 2nd dose varying dosage
                traj(j).code = 1829;
                traj(j).date = '12-11-18';
                traj(j).tdose = 68.38;
                traj(j).seed = 2000;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).WPT = 4;
                traj(j).prevdose = 75;
                traj(j).doseints = 4;
                if ismember(traj(j).column, '2')
                    traj(j).dose = 300;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 150;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 125;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 100;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 75;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 35;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 20;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 10;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end    
            elseif i == 7 %GH1839 3rd dose varying dosage
                traj(j).code = 1839;
                traj(j).date = '8-23-18';
                traj(j).tdose = 72;
                traj(j).seed = 2000;
                traj(j).dosenum = 3;
                traj(j).numdoses = 3;
                traj(j).WPT = 2;
                traj(j).prevdose = [75,100];
                traj(j).doseints = [2,2];
                if ismember(traj(j).column, '2')
                    traj(j).dose = 300;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 150;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 125;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 100;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 75;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 35;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 20;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 10;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end    
            elseif i == 8 %GH1817: 2nd dose varying dosage
                traj(j).code = 1817;
                traj(j).date = '8-19-18';
                traj(j).tdose = 68;
                traj(j).seed = 2000;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).WPT = 2;
                traj(j).prevdose = 75;
                traj(j).doseints = 2;
                if ismember(traj(j).column, '2')
                    traj(j).dose = 300;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 150;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 125;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 100;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 75;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 35;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 20;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 10;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end
            elseif i == 9 %GH1815: 2nd dose varying initial dosage
                traj(j).code = 1815;
                traj(j).date = '7-12-18';
                traj(j).tdose = 76.55;
                traj(j).seed = 2000;
                traj(j).numdoses = 2;
                traj(j).dosenum = 2;
                traj(j).doseints = 2;
                traj(j).WPT = 2;
                traj(j).dose = 100;
                if ismember(traj(j).column, '2')
                    traj(j).dosenum = 1;
                    traj(j).numdoses = 1;
                    traj(j).WPT = [];
                    traj(j).doseints = [];
                    traj(j).prevdose = [];
                elseif ismember(traj(j).column, '3')
                    traj(j).prevdose = 10;
                elseif ismember(traj(j).column, '4')
                    traj(j).prevdose = 20;
                elseif ismember(traj(j).column, '5')
                    traj(j).prevdose = 35;
                elseif ismember(traj(j).column, '6')
                    traj(j).prevdose = 50;
                elseif ismember(traj(j).column, '7')
                    traj(j).prevdose = 75;
                elseif ismember(traj(j).column, '8')
                    traj(j).prevdose = 100;
                elseif ismember(traj(j).column, '9')
                    traj(j).prevdose = 125;
                elseif ismember(traj(j).column, '10')
                    traj(j).prevdose = 150;
                elseif ismember(traj(j).column, '11')    
                    traj(j).prevdose = 300;
                end
            elseif i == 10 %GH1835: 2nd dose varying interval
                traj(j).code = 1835;
                traj(j).date = '12-18-18';
                traj(j).tdose = 60;
                traj(j).celltype = 'BT474';
                traj(j).seed = 2000;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).dose = 35;
                traj(j).prevdose = 35;
                if ismember(traj(j).column, '2')
                    traj(j).WPT = 7;
                    traj(j).doseints = 7;
                elseif ismember(traj(j).column, '3')
                    traj(j).WPT = 6;
                    traj(j).doseints = 6;
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '4')
                    traj(j).WPT = 6;
                    traj(j).doseints = 6;
                elseif ismember(traj(j).column, '5') 
                    traj(j).WPT = 5;
                    traj(j).doseints = 5;
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '6')
                    traj(j).WPT = 5;
                    traj(j).doseints = 5;
                elseif ismember(traj(j).column, '7')
                    traj(j).WPT = 4;
                    traj(j).doseints = 4;
                elseif ismember(traj(j).column, '8')
                    traj(j).WPT = 3;
                    traj(j).doseints = 3;
                elseif ismember(traj(j).column, '9')
                    traj(j).WPT = 2;
                    traj(j).doseints = 2;
                elseif ismember(traj(j).column, '10')
                    traj(j).WPT = 1;
                    traj(j).doseints = 1;
                elseif ismember(traj(j).column, '11')
                    traj(j).WPT = [];
                    traj(j).doseints = [];
                    traj(j).dosenum = 1;
                    traj(j).numdoses = 1;
                    traj(j).prevdose = [];
                end
            elseif i == 11 %GH1902 fourth dose varying dosage
                traj(j).code = 1902;
                traj(j).date = '1-17-19';
                traj(j).tdose = 68.4;
                traj(j).seed = 2000;
                traj(j).dosenum = 4;
                traj(j).prevdose = [75,75,75];
                traj(j).WPT = 4;
                traj(j).numdoses = 4;
                traj(j).doseints = [4,4,4];
                if ismember(traj(j).column, '2')
                    traj(j).dose = 300;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 150;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 125;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 100;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 75;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 35;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 20;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 10;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end
            elseif i == 12 %GH1903 fifth dose varying dosage
                traj(j).code = 1903;
                traj(j).date = '2-13-19';
                traj(j).tdose = 69;
                traj(j).seed = 2000;
                traj(j).dosenum = 5;
                traj(j).prevdose = [75,75,75,75];
                traj(j).WPT = 4;
                traj(j).numdoses = 5;
                traj(j).doseints = [4,4,4,4];
                if ismember(traj(j).column, '2')
                    traj(j).dose = 300;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 150;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 125;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 100;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 75;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 35;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 20;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 10;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end
            elseif i == 13 %GH1827 third dose varying dosage
                traj(j).code = 1827;
                traj(j).celltype = 'BT474';
                traj(j).date = '11-28-18';
                traj(j).tdose = 63.25;
                traj(j).seed = 4000;
                traj(j).dosenum = 3;
                traj(j).prevdose = [35,35];
                traj(j).WPT = 2;
                traj(j).numdoses = 3;
                traj(j).doseints = [2,2];
                if ismember(traj(j).column, '2')
                    traj(j).dose = 120;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 60;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 55;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 42.5;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 35;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 27.5;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 20;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 10;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end
            elseif i == 14 %GH1828 second dose varying dosage
                traj(j).code = 1828;
                traj(j).celltype = 'BT474';                
                traj(j).date = '11-28-18';
                traj(j).tdose = 63.25;
                traj(j).seed = 3000;
                traj(j).dosenum = 2;
                traj(j).prevdose = 35;
                traj(j).WPT = 4;
                traj(j).numdoses = 2;
                traj(j).doseints = 4;
                if ismember(traj(j).column, '2')
                    traj(j).dose = 120;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 60;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 55;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 42.5;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 35;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 27.5;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 20;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 10;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end
            elseif i == 15 %GH1825 first dose varying dosage
                traj(j).code = 1825;
                traj(j).celltype = 'BT474';
                traj(j).date = '10-28-18';
                traj(j).tdose = 78.25;
                traj(j).seed = 2000;
                traj(j).dosenum = 1;
                traj(j).prevdose = [];
                traj(j).WPT = [];
                traj(j).numdoses = 1;
                traj(j).doseints = [];
                if ismember(traj(j).column, '2')
                    traj(j).dose = 120;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 60;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 42.5;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 35;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 27.5;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 20;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 15;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 10;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end
            elseif i == 16 %GH1826 second dose varying dosage
                traj(j).code = 1826;
                traj(j).celltype = 'BT474';                
                traj(j).date = '11-14-18';
                traj(j).tdose = 72;
                traj(j).seed = 4000;
                traj(j).dosenum = 2;
                traj(j).prevdose = 35;
                traj(j).WPT = 2;
                traj(j).numdoses = 2;
                traj(j).doseints = 2;
                if ismember(traj(j).column, '2')
                    traj(j).dose = 120;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 60;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 55;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 42.5;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 35;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 27.5;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 20;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 10;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end
            elseif i == 17 %GH1904 single dose
                traj(j).code = 1904;
                traj(j).date = '2-20-19';
                traj(j).tdose = 63;
                traj(j).seed = 2000;
                traj(j).dosenum = 1;
                traj(j).prevdose = [];
                traj(j).WPT = [];
                traj(j).numdoses = 1;
                traj(j).doseints = [];
                traj(j).celltype = '231';
                if ismember(traj(j).column, '2')
                    traj(j).dose = 500;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 300;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 200;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 150;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 125;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 100;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 75;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 25;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end
            elseif i == 18 %GH1909 single dose
                traj(j).code = 1909;
                traj(j).date = '5-6-19';
                traj(j).tdose = 68.96666667;
                traj(j).seed = 2000;
                traj(j).dosenum = 1;
                traj(j).prevdose = [];
                traj(j).WPT = [];
                traj(j).numdoses = 1;
                traj(j).doseints = [];
                traj(j).celltype = '231';
                if ismember(traj(j).column, '2')
                    traj(j).dose = 1000;
                elseif ismember(traj(j).column, '3')
                    traj(j).dose = 500;
                elseif ismember(traj(j).column, '4')
                    traj(j).dose = 300;
                elseif ismember(traj(j).column, '5')
                    traj(j).dose = 200;
                elseif ismember(traj(j).column, '6')
                    traj(j).dose = 150;
                elseif ismember(traj(j).column, '7')
                    traj(j).dose = 100;
                elseif ismember(traj(j).column, '8')
                    traj(j).dose = 75;
                elseif ismember(traj(j).column, '9')
                    traj(j).dose = 50;
                elseif ismember(traj(j).column, '10')
                    traj(j).dose = 25;
                elseif ismember(traj(j).column, '11')   
                    traj(j).dose = 0;
                end
            elseif i == 19       % GH1910: 2nd dose varying interval
                % What cell type is this?
                traj(j).code = 1910;
                traj(j).dose = 75;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).prevdose = 75;
                traj(j).date = '04-26-19';
                traj(j).seed = 2000;
                % When was the first dose relative to this?? Confused by
                % this whole set up, are both doses monitored or just the
                % second one?
                if ismember(traj(j).column, '2') || ismember(traj(j).column, '3')
                    traj(j).tdose = 93;
                    traj(j).doseintdays = 0;
                elseif ismember(traj(j).column, '4') || ismember(traj(j).column, '5')
                    traj(j).doseintdays = 12;
                    traj(j).tdose = 384;
                elseif ismember(traj(j).column, '6') || ismember(traj(j).column, '7')
                    traj(j).doseintdays = 14;
                    traj(j).tdose = 432;
                elseif ismember(traj(j).column, '8') || ismember(traj(j).column, '9')
                    traj(j).doseintdays = 16;
                    traj(j).tdose = 480;
                elseif ismember(traj(j).column, '10') || ismember(traj(j).column, '11')
                    traj(j).numdoses = 1;
                    traj(j).dosenum = 1;
                end
            elseif i == 20       % GH1911: 2nd dose varying interval
                % Again, what cell type?
                traj(j).code = 1911;
                traj(j).dose = 75;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).prevdose = 75;
                traj(j).tdose = 141;
                traj(j).date = '04-22-19';
                traj(j).seed = 2000;
                if ismember(traj(j).column, '2') || ismember(traj(j).column, '3')
                    traj(j).doseintdays = 2;
                elseif ismember(traj(j).column, '4') || ismember(traj(j).column, '5')
                    traj(j).doseintdays = 4;
                elseif ismember(traj(j).column, '6') || ismember(traj(j).column, '7')
                    traj(j).doseintdays = 6;
                elseif ismember(traj(j).column, '8') || ismember(traj(j).column, '9')
                    traj(j).doseintdays = 8;
                elseif ismember(traj(j).column, '10') || ismember(traj(j).column, '11')
                    traj(j).doseintdays = 10;
                end    
            elseif i == 21       % GH1912: varying num doses
                % What cell type?
                traj(j).code = 1912;
                traj(j).dose = 75;
                traj(j).date = '06-17-19';
                traj(j).seed = 2000;
                if ismember(traj(j).column, '2') || ismember(traj(j).column, '3')
                    traj(j).tdose = 72;
                    traj(j).dosenum = 5;
                    traj(j).numdoses = 5;
                    traj(j).prevdose = [75,75,75,75];
                    traj(j).doseintdays = [13,13,13,13];
                    traj(j).tdose = 406;
                elseif ismember(traj(j).column, '4') || ismember(traj(j).column, '5')
                    traj(j).tdose = 742;
                    traj(j).dosenum = 4;
                    traj(j).numdoses = 4;
                    traj(j).prevdose = [75,75,75];
                    traj(j).doseintdays = [13,13,13];
                elseif ismember(traj(j).column, '6') || ismember(traj(j).column, '7')
                    traj(j).tdose = 1077;
                    traj(j).dosenum = 3;
                    traj(j).numdoses = 3;
                    traj(j).prevdose = [75,75];
                    traj(j).doseintdays = [13,13];
                elseif ismember(traj(j).column, '8') || ismember(traj(j).column, '9')
                    traj(j).tdose = 1413;
                    traj(j).dosenum = 2;
                    traj(j).numdoses = 2;
                    traj(j).prevdose = 75;
                    traj(j).doseintdays = 13;
                elseif ismember(traj(j).column, '10') || ismember(traj(j).column, '11')
                    traj(j).dosenum = 1;
                    traj(j).numdoses = 1;
                end    
            elseif i == 22       % GH1919: varying int
                traj(j).code = 1919;
                traj(j).celltype = 'BT474';
                traj(j).dose = 35;
                traj(j).date = '08-09-19';
                traj(j).seed = 2000;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).prevdose = 35;
                if ismember(traj(j).column, '2') || ismember(traj(j).column, '3')
                    traj(j).tdose = 63;
                    traj(j).doseintdays = 0;
                elseif ismember(traj(j).column, '4') || ismember(traj(j).column, '5')
                    traj(j).tdose = 138;
                    traj(j).doseintdays = 2;
                elseif ismember(traj(j).column, '6') || ismember(traj(j).column, '7')
                    traj(j).tdose = 186;
                    traj(j).doseintdays = 4;
                elseif ismember(traj(j).column, '8') || ismember(traj(j).column, '9')
                    traj(j).tdose = 234;
                    traj(j).doseintdays = 6;
                elseif ismember(traj(j).column, '10') || ismember(traj(j).column, '11')
                    traj(j).tdose = 282;
                    traj(j).numdoses = 8;
                end                 
            elseif i == 23       % GH1920: varying int
                traj(j).code = 1920;
                traj(j).celltype = 'BT474';
                traj(j).dose = 35;
                traj(j).date = '08-09-19';
                traj(j).seed = 2000;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).prevdose = 35;
                if ismember(traj(j).column, '2') || ismember(traj(j).column, '3')
                    traj(j).tdose = 332;
                    traj(j).doseintdays = 10;
                elseif ismember(traj(j).column, '4') || ismember(traj(j).column, '5')
                    traj(j).tdose = 380;
                    traj(j).doseintdays = 12;
                elseif ismember(traj(j).column, '6') || ismember(traj(j).column, '7')
                    traj(j).tdose = 427;
                    traj(j).doseintdays = 14;
                elseif ismember(traj(j).column, '8') || ismember(traj(j).column, '9')
                    traj(j).tdose = 475;
                    traj(j).doseintdays = 16;
                elseif ismember(traj(j).column, '10') || ismember(traj(j).column, '11')
                    traj(j).tdose = 523;
                    traj(j).numdoses = 18;
                end                    
            elseif i == 24       % GH1907: varying int
                %THIS SECTION NEEDS CLARIFICATION
                traj(j).code = 1907;
                traj(j).celltype = '231';
                traj(j).dose = 200;
                traj(j).date = '05-06-19';
                traj(j).seed = 2000;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).prevdose = 200;
                if ismember(traj(j).column, '2') || ismember(traj(j).column, '3')
                    % Is this one a double dose? or just one dose at t=93
                    % hrs? Then we should change numdoses to be =1
                    traj(j).tdose = [69,93];
                    traj(j).doseintdays = 0;
                
                elseif ismember(traj(j).column, '6') || ismember(traj(j).column, '7')
      % I changed this so that if both doses are observed in the data set,
      % that is recorded as two doses
                    traj(j).tdose = [69,141];
                    traj(j).numdoses = 2;
                    traj(j).doseintdays = 2; % Added this 
                end                    
            elseif i == 25       % GH1908: varying int
                % THIS SECTION NEEDS CLARIFICATION
                traj(j).code = 1908;
                traj(j).celltype = '231';
                traj(j).dose = 200;
                traj(j).date = '05-08-19';
                traj(j).seed = 2000;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).prevdose = 200;
                %PLEASE CHECK THAT THIS IS CORRECT
                if ismember(traj(j).column, '2') || ismember(traj(j).column, '3')
                    % are tdoses the start or end of the dose?
                    traj(j).tdose = [69,189];
                    traj(j).doseintdays = 4;
                elseif ismember(traj(j).column, '4') || ismember(traj(j).column, '5')
                    traj(j).tdose = [69,237];
                    traj(j).doseintdays = 6;
                elseif ismember(traj(j).column, '6') || ismember(traj(j).column, '7')
                    traj(j).tdose = [69,285];
                    traj(j).doseintdays = 8;
                elseif ismember(traj(j).column, '8') || ismember(traj(j).column, '9')
                    traj(j).tdose = [93,333];
                    traj(j).doseintdays = 10;
                elseif ismember(traj(j).column, '10') || ismember(traj(j).column, '11')
                    % is this the start or end of the dose?
                    traj(j).tdose = 69;
                    traj(j).numdoses = 1;
                    traj(j).dosenum = 1;
                    traj(j).prevdose = [];
                end                                   
            elseif i == 26       % GH1914: varying int
                %THIS SECTION NEEDS CLARIFICATION
                % Is this plate monitored during both doses or just the
                % second? Please adjust tdose accordingly
                traj(j).code = 1914;
                traj(j).celltype = '231';
                traj(j).dose = 200;
                traj(j).date = '06-19-19';
                traj(j).seed = 2000;
                traj(j).dosenum = 2;
                traj(j).numdoses = 2;
                traj(j).prevdose = 200;
                if ismember(traj(j).column, '2') || ismember(traj(j).column, '3')
                    % were these replated? The tdose and the doseintdays
                    % don't make sens
                    traj(j).tdose = [64,376];
                    traj(j).doseintdays = 12;
                elseif ismember(traj(j).column, '4') || ismember(traj(j).column, '5')
                    traj(j).tdose = [64,424];
                    traj(j).doseintdays = 14;
                elseif ismember(traj(j).column, '6') || ismember(traj(j).column, '7')
                   
                    traj(j).tdose = [64, 472];
                    traj(j).doseintdays = 16;
                end                    
            elseif i == 27       % GH1918: varying num doses
                % THIS SECTION NEEDS CLARIFICATION
                traj(j).code = 1918;
                traj(j).dose = 200;
                traj(j).celltype = '231';
                traj(j).date = '06-17-19';
                traj(j).seed = 4000;
                if ismember(traj(j).column, '2') || ismember(traj(j).column, '3') ...
                        || ismember(traj(j).column, '4') || ismember(traj(j).column, '5')
                    % if multiple doses are monitored in this plate, please
                    % update tdose accordingly. I did my best to do this
                    % but check work
                    traj(j).tdose = [72, 408, 743];
                    traj(j).dosenum = 3;
                    traj(j).numdoses = 3;
                    traj(j).prevdose = [200, 200];
                    traj(j).doseintdays = [13,13];
                elseif ismember(traj(j).column, '6') || ismember(traj(j).column, '7') ...
                        || ismember(traj(j).column, '8') || ismember(traj(j).column, '9')
                    traj(j).tdose = [72, 408];
                    traj(j).dosenum = 2;
                    traj(j).numdoses = 2;
                    traj(j).prevdose = 200;
                    traj(j).doseintdays = 13;
                elseif ismember(traj(j).column, '10') || ismember(traj(j).column, '11')
                    traj(j).tdose = 72;
                    traj(j).dosenum = 1;
                    traj(j).numdoses = 1;
                end                    
        end
        traj(j).accdose = traj(j).dose + sum(traj(j).prevdose);
        if i < 19
            traj(j).doseintdays = traj(j).doseints * 7 - 1;
        end
    end
end
%% random sorting stuff
T = struct2table(traj); % converts traj to a table
sortedT = sortrows(T, 'numdoses'); % sort the table by '(field)'
% check "sortedT"
%% Add color by WPT 
for j = 1:length(traj)
    if ~isempty(traj(j).WPT)
    WPT(j) = traj(j).WPT;
    end
end
colorsets = varycolor(length(unique(WPT))+2);
uniqWPT= unique(WPT);
for i = 1:length(traj)
    traj(i).color = [];
    for j = 1:length(uniqWPT)
        if traj(i).WPT==uniqWPT(j)
            %used to be if flag = 0, none are flagged now
            traj(i).color =colorsets(j,:);
        end
    end
    if isempty(traj(i).WPT)
        traj(i).color = colorsets(end,:); % make untreated control black
    end
end
%% Add color by dose (regardless of how many there are!)
for j = 1:length(traj)
    
    if ~isempty(traj(j).dose)
    dose(j) = traj(j).dose;
    end

end
colorsets = varycolor(length(unique(dose))+1);
uniqdose= unique(dose);
for i = 1:length(traj)
   
    for j = 1:length(uniqdose)
        if traj(i).dose==uniqdose(j)
            %used to be if flag = 0, now none are flagged
            traj(i).color =colorsets(j,:);
        end
    end
    if isempty(traj(i).dose)
        traj(i).color = colorsets(end,:); % make untreated control black
    end
end
%% Save the raw data structure data sets
% This saves the traj structure just containing raw data as trajraw.mat
save('../out/trajraw.mat', 'traj')