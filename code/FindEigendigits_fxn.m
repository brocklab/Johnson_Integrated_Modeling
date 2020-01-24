function [m, V] = FindEigendigits_fxn(A)

% Convert A to doubles and let this equal A0 
A0 = im2double(A);

% find mu (mean vector)
% finds the mean along the rows
m = mean(A0,2);

% let X = the A0-m
X = A0-m;
% X is a k x n matrix
% find your large X where
% Sigma = 1/n*sum(X_k)*(X_k)'
% Normally Sigma = X*X'

% In the following comments A = X
% Typically we have:
% Ax=lambda*x
% and we find the eigenvectors of A

% Here we can work with the smaller system and let x = A*v
% A*A'*A*v = mu*A*v
% A*A' =  big cov matrix

% A'*A = small cov matrix

% Find the eigenvectors of the small covariance matrix
smallcov = (X'*X); % this is an n x n where n = # images
[Vsmall, D] = eig(smallcov);
lambdas = diag(D);
[orderedvals, ind] = sort(lambdas, 'descend');

%reorder
Vnew = Vsmall(:,ind); 

% If Vsmall is an eigenvector of A'A then A*Vsmall is an eigenvector of the
% big covariance matrix

% multiply by A0 to get eignevectors of big covariance matrix

%Vbig = A0*Vnew;
Vbig = X*Vnew;
% Want to normalize each eigenvector
k = length(lambdas);
for j = 1:k
    % find the 2-norm of each column
    norms(j) = norm(Vbig(:,j));
    
    V(:,j) = Vbig(:,j)./norms(j);
end





end