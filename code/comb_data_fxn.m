function [new_data_struct] = comb_data_fxn(large_struct, filter_criteria, ntreat, dose, date)
% this is a function that is intended to be used to comb through structure
% of all data and output a data structure containing only the trajectories
% that follow the set of conditions specified or which were taken on that
% date, and finds the mean and standard deviation of those with the same
% condition
traj = large_struct;

switch filter_criteria
    case 'treatmentnum'
    
        uniqWPT = [];
        WPTlist = [];
        for i = 1:length(traj)
            if ~isempty(traj(i).dosenum)
            if traj(i).dosenum==ntreat 
            WPTlist =vertcat(WPTlist, traj(i).WPT);
            end
            end
        end
        uniqWPT= unique(WPTlist);
        % Make a new structure which combines each WPT
         for i = 1:length(uniqWPT)
            traj2(i).Cdox = [];
            traj2(i).Nmat = [];
            traj2(i).nreps = 0;
            traj2(i).tmat = [];
         end
         for i = 1:length(uniqWPT) % number of unique seed numbers
             for j = 1:length(traj)
                 % pull from the same experiment: first treat
                 if contains(traj(j).date, date)  % only want data from this run
                     if traj(j).WPT == uniqWPT(i)
                         if traj(j).dosenum == ntreat
                             traj2(i).nreps = traj2(i).nreps +1;
                             traj2(i).Cdox = traj(j).dose;
                             traj2(i).color = traj(j).color;
                             traj2(i).tmat = horzcat(traj2(i).tmat,traj(j).time);
                             traj2(i).Nmat = horzcat(traj2(i).Nmat, traj(j).rawN);
                             traj2(i).tdose = traj(j).tdose;
                             traj2(i).WPT = traj(j).WPT;
                         end
                     end
                 end
             end
         end
         
         for i = 1:length(traj2)
             if i ==1
                 Nfin = 5.5e4;
             else
                 Nfin = 4e4;
             end
             N = traj2(i).Nmat;
             t = traj2(i).tmat;
             i0 = find(t(:,1)>traj2(i).tdose,1,'first'); % I arbitrarily search for a maximum in the first 200 hours
             iend = find(N(:,1)>=Nfin,1, 'first');
             if ~isempty(iend)
                 tfit = t(i0:iend,:)-t(i0, :);
                 Nfit = N(i0:iend, :);
             end
             if isempty(iend)
                 tfit = t(i0:end, :)-t(i0, :);
                 Nfit = N(i0:end, :);
             end
             
             traj2(i).tfit =round(tfit,0);
             traj2(i).Nfit =Nfit;
             traj2(i).tvec = traj2(i).tfit(:,1);
         end
         for i = 1:length(traj2)
             traj2(i).Nmean = mean(traj2(i).Nfit,2);
             traj2(i).tvec = round(traj2(i).tfit(:,1),0);
             traj2(i).Nstd = std(traj2(i).Nfit,0,2);
         end
         
         
         
         
         
         new_data_struct = traj2;

end