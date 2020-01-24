function[Xtr, ytr] =make_labeled_data(trainImages, trainLabels, digpos, digneg, N)
% This function outputs a matrix Xtr where each column is a 784 x 1 vector
% of pixels from the digits data set
% ytr is the corresponding decision (+1 for digit we want to identify, -1
% for digit we want to dientify from)

% dig pos just denotes which digits we label as having a positive ytr value
% vs negative. i.e. if we want to distinguish between 2 and 4 digpos is 2
% and digneg is 4. Note these are abritrary and should give you the same
% results


A0=[];
trlabs = [];
ytr = [];
ctdig = 0;
for i = 1:N
    % make the column vector A_k
        if ismember(trainLabels(i),digpos)  
            ctdig = ctdig+1;
        M = trainImages(:,:,1,i); % M is square matrix
        l= size(M,1)*size(M,2);
        Acol = reshape(M, [l,1]);
        trlabs=horzcat(trlabs, 1);
        A0 = horzcat(A0,Acol);
        end
        
        
        if ismember(trainLabels(i), digneg) 
            M = trainImages(:,:,1,i); % M is square matrix
            l= size(M,1)*size(M,2);
            Acol = reshape(M, [l,1]);
            trlabs =horzcat(trlabs, 0);% set the new trlabs vector =0 to denote any others
            A0 = horzcat(A0,Acol);   
        end
end
A =im2double(A0);
Xtr = A;
% Make a -1 +1 vector of labels
for i = 1:length(trlabs)
    if trlabs(i)> 0
        ytr(i,1) = 1;
    end
    if trlabs(i) == 0
        ytr(i,1) = -1;
    end
end



end