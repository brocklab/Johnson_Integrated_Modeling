% This script is intended to be used to refer to when making the classifier
% for our single cell gene expression data based on survival or
% non-survival. It uses the 784 pixels in an image to classify the digit of
% the image by finding the eigenvectors of a number of training images, and
% then projecting new images onto eigenspace and using the coordinates to
% identify the class of the obejct most similar to its coordinates

% 1. Load the data
close all; clear all; clc
%loads in the digits.mat file from the data folder
load ('../data/digits_ML.mat');

whos;

%% Preview an example image from the training set
i = 500;
M=trainImages(:,:,1,i);
figure;
imagesc(M)
title(['Known digit=', num2str(trainLabels(i))])


%% Make the x by k matrix A where x is the total number of pixels in an
% image and k is the number of training images
% So our A for the cell classifier will be x by k where x is the total
% number of genes and k is the number of cells in our training data set

k = 400; % this is the number of training images used
for i = 1:k
    % make the column vector A_k
    M = trainImages(:,:,1,i); % M is square matrix
    x= size(M,1)*size(M,2);
    Acol = reshape(M, [x,1]); % we shouldn't have to do this since already in vector form
    A(:,i) = Acol;
end


%% 2. Find Eigendigits
% Call 
% should take the input matrix A and return a vector m of length x
% containing the mean column vector of A and an (x by k) matrix V
% containing k eigenvectors of the covariance matrix of A

[m, V] = FindEigendigits_fxn(A);
%% Visualize some mean and eigenvector images...

% Won't be doing this for cell classifier since we can't visualize it 
figure;
subplot(1,3,1)
Mfig = reshape(m, [size(M,1), size(M,2)]);
imagesc(Mfig)
title('Mean matrix')
subplot(1,3,2)
v1 = reshape(V(:,1),[size(M,1), size(M,2)]);
imagesc(v1)
title('1st eigenvector')
subplot(1,3,3)
v2 = reshape(V(:,2), [size(M,1), size(M,2)]);
imagesc(v2)
title('2nd eigenvector')


% check that eigenvectors are normalized
test = norm(V(:,1));
%% 3. Experiments
% We are going to want to run these on our pre-treatment (training data) to
% optimize the number of eigenvectors used

% Use different amounts of training data to generate eigenvectors
% Started with k = 5 and go to k = 400
k = 400;
A = uint8.empty([x,0]);
for i = 1:k
    % make the column vector A_k
    M = trainImages(:,:,1,i); % M is square matrix
    x= size(M,1)*size(M,2);
    Acol = reshape(M, [x,1]);
    A(:,i) = Acol;
end

[m, V] = FindEigendigits_fxn(A);

% Display eigenvectors from larger training set
figure;
for i = 1:3
subplot(1,3,i)
V1 = reshape(V(:,i), [sqrt(x), sqrt(x)]);
imagesc(V1)
title(['Eigenvector #', num2str(i),' for k = ', num2str(k)])
end

%% Project data points into eigenspace

% pick a random digit image from test digits to display
for g = 3:5

% Put back into matrix form to make an image
reconimg = reshape(recon, [sqrt(x), sqrt(x)]);% project the digit into eigenspace
M = im2double(testImages(:,:,1,g)); 
x= size(M,1)*size(M,2);
Mcol = reshape(M, [x,1]);
b = (Mcol-m)'*V; %(1 by x) *(x*k) Gives a vector of length k

% Project onto eigenspace
recon = V*b'; % (x by k) * (k *1) Gives x by 1

figure;
subplot(1,2,1)
imagesc(M)
title(['Orignal test image, digit=', num2str(testLabels(g))])
subplot(1,2,2)
imagesc(reconimg)
title('reconstructed image')
end
% See that the images do not look exactly like original digits

%% Try trimming vector and using the top n eigenvectors
figure;
for i = 1:5
testn =[ 10 20 50 100 400];  
n = testn(i);

btrim = b(1:n);

recontest = V(:, 1:n)*btrim'; % (x by k) * (k *1) Gives x by 1
% Put back into matrix form to make an image
reconimgtest = reshape(recontest, [sqrt(x), sqrt(x)]);

% subplot(2,5,i)
% imagesc(M)
% title(['Orignal image, digit=', num2str(testLabels(g))])
subplot(1,5,i)
imagesc(reconimgtest)
title(['Recon n=', num2str(n), ' eigenvectors'])

end

% As we can see, as we increase the number of eigenvectors we get visually
% better and better reconstruction of the projected image.

%% Classify test images using Euclidean distance
% We will do this using the label of S or R from the gene expression
% matrix.

% First vary the number of training images used
 num_trainimgs = [ 20, 50, 100, 400, 1000, 5000];
 n = 20; % number of eigenvectors used (keep this constant)
for r = 1:length(num_trainimgs)
    
k= num_trainimgs(r); % number of training images used
num_testimgs = 100;
ct_true = 0;
A = uint8.empty([x,0]);
V = [];
% find eigenvectors
for i = 1:k
    % make the column vector A_k
    M = trainImages(:,:,1,i); % M is square matrix
    x= size(M,1)*size(M,2);
    Acol = reshape(M, [x,1]);
    A(:,i) = Acol;
end
[m, V] = FindEigendigits_fxn(A);
% Compute omegass for each training image
Omat = zeros([n+1, k]); % Make a matrix  where each column is the w for each training image
% last row of each column is the digit it correspodns to

for i = 1:k
Mtr = im2double(trainImages(:,:,1,i)); 
x= size(M,1)*size(M,2);
Mcol = reshape(Mtr, [x,1]);
b = (Mcol-m)'*V;
omegatr = b(1:n); % find coordinates of each training image in eigenvector space
Omat(1:end-1,i) = omegatr;
Omat(n+1,i) = trainLabels(i);
end

for g = 1:num_testimgs
M = im2double(testImages(:,:,1,g)); 
x= size(M,1)*size(M,2);
Mcol = reshape(M, [x,1]);
b = (Mcol-m)'*V; % this gives you the coordinates of the new image in eigenvector spacee
omega = b(1:n);
% find minimum
Eudist = Omat(1:n,:)-omega';
[val, index] = min(vecnorm(Eudist));
bfdigit = Omat(end, index);
truedigit = testLabels(g);

if bfdigit ==truedigit
    ct_true = ct_true +1;
end
end

acc(r) = ct_true/num_testimgs;
end
%%
figure;
subplot(1,2,1)
plot(num_trainimgs, acc, '*-', 'LineWidth',2)
xlabel(' Number of training images', 'FontSize',14)
ylabel('Accuracy in digit classification','FontSize',14)
title('Accuracy vs. number of training images used','FontSize',14)
% Vary number of eigenvectors used
k = 1000; % keep number of training images constant
A = uint8.empty([x,0]);
V = [];
 num_eigens = [ 5, 10, 20, 50, 100, 400, 1000];
for r = 1:length(num_eigens)
    

num_testimgs = 100;
ct_true = 0;
n = num_eigens(r);
% find eigenvectors
for i = 1:k
    % make the column vector A_k
    M = trainImages(:,:,1,i); % M is square matrix
    x= size(M,1)*size(M,2);
    Acol = reshape(M, [x,1]);
    A(:,i) = Acol;
end
[m, V] = FindEigendigits_fxn(A);
% Compute omegass for each training image
Omat = zeros([n+1, k]); % Make a matrix  where each column is the w for each training image
% last row of each column is the digit it correspodns to

for i = 1:k
Mtr = im2double(trainImages(:,:,1,i)); 
x= size(M,1)*size(M,2);
Mcol = reshape(Mtr, [x,1]);
b = (Mcol-m)'*V;
omegatr = b(1:n); % find coordinates of each training image in eigenvector space
Omat(1:end-1,i) = omegatr;
Omat(n+1,i) = trainLabels(i);
end

for g = 1:num_testimgs
M = im2double(testImages(:,:,1,g)); 
x= size(M,1)*size(M,2);
Mcol = reshape(M, [x,1]);
b = (Mcol-m)'*V; % this gives you the coordinates of the new image in eigenvector spacee
omega = b(1:n);
% find minimum
Eudist = Omat(1:n,:)-omega';
[val, index] = min(vecnorm(Eudist));
bfdigit = Omat(end, index);
truedigit = testLabels(g);

if bfdigit ==truedigit
    ct_true = ct_true +1;
end
end

acc_n(r) = ct_true/num_testimgs;
end

subplot(1,2,2)
plot(num_eigens, acc_n, '*-', 'LineWidth',2)
xlabel(' Number of eigenvectors used', 'FontSize',14)
ylabel('Accuracy in digit classification','FontSize',14)
title('Accuracy vs. number of eigenvectors used','FontSize',14)

%% Vary number of neighbors in k=nearest neighbors
k = 1000; % use a relatively high number of training images 
n = 50; % use the most optimal number of eigenvectors from previous plot
A = uint8.empty([x,0]);
V = [];   
num_testimgs = 100;
ct_true = 0;
% find eigenvectors
for i = 1:k
    % make the column vector A_k
    M = trainImages(:,:,1,i); % M is square matrix
    x= size(M,1)*size(M,2);
    Acol = reshape(M, [x,1]);
    A(:,i) = Acol;
end
[m, V] = hw1FindEigendigits(A);
% Compute omegass for each training image
Omat = zeros([n+1, k]); % Make a matrix  where each column is the w for each training image
% last row of each column is the digit it correspodns to

for i = 1:k
Mtr = im2double(trainImages(:,:,1,i)); 
x= size(M,1)*size(M,2);
Mcol = reshape(Mtr, [x,1]);
b = (Mcol-m)'*V;
omegatr = b(1:n); % find coordinates of each training image in eigenvector space
Omat(1:end-1,i) = omegatr;
Omat(n+1,i) = trainLabels(i);
end
Omat = Omat';

for g = 1:num_testimgs
M = im2double(testImages(:,:,1,g)); 
x= size(M,1)*size(M,2);
Mcol = reshape(M, [x,1]);
b = (Mcol-m)'*V; % this gives you the coordinates of the new image in eigenvector spacee
omega = b(1:n);
Otest(g,:) = omega;
end
%% KNN classifoer
% use knn search to compare omega (test digit) with Omat (training digits)
for k = 1:5
    ks = [1 3 5 10 20];
ktest=ks(k);
ct_true = 0;
ct_truem=0;
ct_truec=0;
[Idx, Dd] = knnsearch(Omat(:,1:n), Otest, 'K', ktest, 'Distance', 'euclidean');
[mIdx, Dm] = knnsearch(Omat(:,1:n), Otest, 'K', ktest, 'Distance', 'minkowski', 'P',5);
[cIdx, Dc] = knnsearch(Omat(:,1:n), Otest, 'K', ktest, 'Distance', 'chebychev');

for i = 1:num_testimgs
    best_dig = Omat(Idx(i,:), end);
    best_digm= Omat(mIdx(i,:),end);
    best_digc = Omat(cIdx(i,:),end);
    fitted_dig =round(mean(best_dig),0);
    if fitted_dig == testLabels(i)
        ct_true = ct_true +1;
    end
    fitted_digm =round(mean(best_digm),0);
    if fitted_digm == testLabels(i)
        ct_truem = ct_truem +1;
    end
    fitted_digc =round(mean(best_digc),0);
    if fitted_digc == testLabels(i)
        ct_truec = ct_truec +1;
    end
end

acc_k(k) = ct_true/num_testimgs
acc_m(k) = ct_truem/num_testimgs
acc_c(k) = ct_truec/num_testimgs
end
figure;
plot(ks, acc_k, '*-', 'LineWidth',2)
hold on
plot(ks, acc_m, '*-', 'LineWidth',2)
plot(ks, acc_c, '*-', 'LineWidth',2)
legend('Euclidean', 'Minkowski', 'Chebychev')
legend boxoff
xlabel(' K neighbors', 'FontSize',14)
ylabel('Accuracy in digit classification','FontSize',14)
title('Accuracy vs. k nearest neighbors','FontSize',14)

%%


figure;
plot(num_eigens, acc_n, '*-', 'LineWidth',2)
xlabel(' Number of eigenvectors used', 'FontSize',14)
ylabel('Accuracy in digit classification','FontSize',14)
title('Accuracy vs. number of eigenvectors used','FontSize',14)