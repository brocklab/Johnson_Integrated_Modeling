#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 16:27:36 2019

@author: kj22643
"""


# -*- coding: utf-8 -*-
"""
Created on Wed May 29 12:29:10 2019

@author: kj22643
"""

%reset

import numpy as np
import pandas as pd
import os
import scanpy as sc
import seaborn as sns
from plotnine import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import plot, show, draw, figure, cm
import matplotlib as plt
import random
os.chdir('/Users/kj22643/Documents/Documents/231_Classifier_Project/code')
from func_file import find_meanvec
from func_file import find_eigvals
from func_file import find_eigvecs
from func_file import project_cells
path = '/Users/kj22643/Documents/Documents/231_Classifier_Project/data'
#path = '/stor/scratch/Brock/231_10X_data/'
os.chdir(path)
sc.settings.figdir = 'KJ_plots'
sc.set_figure_params(dpi_save=300)
sc.settings.verbosity = 3
#%% Load in pre and post data
adata = sc.read('daylin_anndata.h5ad')
adata.obs.head()
#BgL1K
#30hr
#Rel-1
#Rel-2
#%% Assign survivor category in adata.obs
longTreatLins = adata.obs.loc[(adata.obs['sample'].isin(['Rel-1','Rel-2']))&(adata.obs.lineage!='nan'),'lineage'].unique().tolist()

adata.obs.loc[adata.obs.lineage.isin(longTreatLins)==False,'survivor'] = 'sens'
adata.obs.loc[adata.obs.lineage.isin(longTreatLins)==True,'survivor'] = 'res'



sc.pl.umap(adata,color=['survivor'],wspace=0.3,
           save='alltps_res_sens.png')
# %%try to rename the samples by time point
samps= adata.obs['sample'].unique()

timepoint = np.array(['t=0hr', 't=30hr', 't=1344hr'])

adata.obs.loc[adata.obs['sample']==samps[0], 'timepoint']='t=0hr'
adata.obs.loc[adata.obs['sample']==samps[1], 'timepoint']='t=30hr'
adata.obs.loc[adata.obs['sample']==samps[2], 'timepoint']='t=1344hr'
adata.obs.loc[adata.obs['sample']==samps[3], 'timepoint']='t=1344hr'

print(adata.obs['timepoint'].unique())


sc.pl.umap(adata,color = ['timepoint'], palette=['#2c9e2f','#046df7', '#d604f7', '#c91212'], wspace=0.3,
           save = 'TPs_umap.png')
#%% Separately make dataframes for the pre-treatment, intermediate, and post treatment samples
# t=0 hr (pre-treatment), 3182 pre treatment cells
adata_pre = adata[adata.obs['timepoint']=='t=0hr', :]
dfpre = pd.concat([adata_pre.obs['survivor'],
               pd.DataFrame(adata_pre.raw.X,index=adata_pre.obs.index,
                            columns=adata_pre.var_names),],axis=1) 
# t = 30 hr (intermediate timepoint) 5169 int treatment cells
adata_int = adata[adata.obs['timepoint']=='t=30hr', :]
dfint = pd.DataFrame(adata_int.raw.X, index=adata_int.obs.index, columns = adata_int.var_names)

# t=1344 hr (~roughly 8 weeks), 10332 post treatment cells
adata_post = adata[adata.obs['timepoint']=='t=1344hr', :]
dfpost = pd.DataFrame(adata_post.raw.X, index=adata_post.obs.index, columns = adata_post.var_names)

#%% Try making a UMAP of the first sample only

sc.pl.umap(adata_pre,color=['survivor'],wspace=0.3,
           save='pre_treat_res_sens.png')
#%% Play with scanpys PCA
sc.tl.pca(adata_pre, n_comps=50, zero_center=True, svd_solver='auto', random_state=0, return_info=False, use_highly_variable=None, dtype='float32', copy=False, chunked=False, chunk_size=None)
#%%
classvecser= adata_pre.obs['survivor']
classvec = pd.DataFrame(classvecser)

PCs=adata_pre.obsm['X_pca']
PCdf = pd.DataFrame(PCs)
classvec.reset_index(drop=True, inplace=True)
PCdf.reset_index(drop=True, inplace=True)

PC_df=pd.concat([classvec['survivor'],PCdf], axis =1)
#%%
sns.set_style('white')
from matplotlib.pyplot import plot, show, draw, figure, cm
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(6,6))


ax=sns.scatterplot(PC_df[0], PC_df[1], hue= PC_df['survivor'])
ax.set(xlabel ='PC1', ylabel ='PC2') 

ax1=sns.scatterplot(PC_df[1], PC_df[2], hue= PC_df['survivor'])
ax1.set(xlabel ='PC2', ylabel ='PC3') 

ax2=sns.scatterplot(PC_df[2], PC_df[3], hue= PC_df['survivor'])
ax2.set(xlabel ='PC3', ylabel ='PC4') 

ax3=sns.scatterplot(PC_df[0], PC_df[2], hue= PC_df['survivor'])
ax3.set(xlabel ='PC1', ylabel ='PC3') 

ax4=sns.scatterplot(PC_df[0], PC_df[3], hue= PC_df['survivor'])
ax4.set(xlabel ='PC1', ylabel ='PC4')

ax5=sns.scatterplot(PC_df[1], PC_df[3], hue= PC_df['survivor'])
ax5.set(xlabel ='PC2', ylabel ='PC4')
#%% ATTEMPT AT MAKING a 3D scatter plot with PCs 1, 2, & 3

fig = plt.figure(figsize=(15,15))
ax=fig.add_subplot(111,projection='3d')

PC_dfsens = PC_df[PC_df['survivor']=='sens']
PC_dfres = PC_df[PC_df['survivor']=='res']
Xs= np.asarray(PC_dfsens[0])
Ys=np.asarray(PC_dfsens[2])
Zs=np.asarray(PC_dfsens[3])

Xr= np.asarray(PC_dfres[0])
Yr=np.asarray(PC_dfres[1])
Zr=np.asarray(PC_dfres[2])

ax.scatter(Xr, Yr, Zr, c='b', marker='^', alpha = 1)
ax.scatter(Xs, Ys, Zs, c='r', marker='o', alpha = 0.3)


ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')

ax.azim = 100
ax.elev = -50
# NEXT NEED TO FIND OUT HOW TO OUT ARROWS ON THIS

#%% PCA Overview
sc.pl.pca_overview(adata_pre)
#%%
loadings=adata_pre.varm['PCs']

#%%
print(dfpre) # 22192 columns corresponding to 22191 genes
#%% Make series that label the pre-treatment cells as res/sens and label the 
# label the post treatment cells by their sample
labelsdfpre = dfpre['survivor']
print(labelsdfpre)

#%% Make matrices (data frames) of just the cell-gene matrix for the pre treatment and 
# post treatment samples
genematpre = dfpre.loc[:, dfpre.columns !='survivor']

print(genematpre)
# Now genematpre and genemat post are your ncells rows x ngenes columns gene 
# expression matrices.
#%% Now try to emulate your matlab code... 
# Start with just your pre-treatment time point
# In Matlab we have an x by k where x would be all the genes and k are the indivdual 
# cells (so each column is a cell and each row is a gene)

# let's see if we can make that in python and call it Adf
nint = dfint.shape[0]
npost =dfpost.shape[0]
npre = genematpre.shape[0] # this gets the number of rows in the df (number of cells)
# Set your k for your k-fold cross validation to divide up testing and training data sets
kCV=4
ntrain = round(((kCV-1)/kCV)*npre)+1  # start by setting the number of training cells to 1/10th 
ntest = npre-ntrain
#%% Make your full data frames (include both the testing and training data set from the pre-treatment time point)
# Call these Adf
AdfT= genematpre
Adf = AdfT.T
print(Adf)
AintT=dfint
Aint= AintT.T
ApoT= dfpost
Apost= ApoT.T
#%% Make susbets of the full data frame for training and testing
# Going to perform k-fold CV with k=4

indexes = [i for i, _ in enumerate(Adf)]
ordered_ind= random.sample(indexes, ntest*kCV)
train_folds = {}
test_folds = {}
labstest = {}
labstrain = {}
Vfolds = {}
lamfolds = {}
mfolds = {}


#%% Make dictionaries that hold training and testing data set labels. 
for i in range(kCV):
# test out grabbing a subset of columns from a dataframe
    itest=ordered_ind[i*ntest:ntest*(i+1)]
    itrain =[j for j in indexes if j not in itest]
    Atest = Adf.iloc[:,itest]
    Atrain = Adf.iloc[:,itrain]
    labstest[i] = labelsdfpre.iloc[itest]
    labstrain[i] = labelsdfpre.iloc[itrain]
    # Save the testing and training data frames into a data frame that holds all of them.    
    train_folds[i] = Atrain
    test_folds[i] = Atest
    # Find the mean vector, eigenvalues, and eigenvectors
    m = find_meanvec(A=Atrain)
    lams = find_eigvals(A=Atrain, mu=m )
    V = find_eigvecs(A=Atrain, mu = m)
    # Save the training set mean vector, eigenvalues, and eigenvectors
    Vfolds[i]=V
    lamfolds[i]=lams
    mfolds[i]=m
    # find the coordinates of the training and testing cells by projecting onto eigenvectors (V)

#%% Make a dictionary which contains the training matrix, labels, eigenvectors,
    # eigen values, and mean vectors
full_dict = {'trainmat':{}, 'trainlabel':{}, 'eigenvectors':{}, 'eigvals':{}, 'meanvec':{}}
full_dict['trainmat']= Adf
full_dict['trainlabel']=labelsdfpre
mpre = find_meanvec(A=Adf)
lamspre = find_eigvals(A=Adf, mu = mpre)
Vall=find_eigvecs(A=Adf, mu = mpre)
full_dict['eigenvectors']=Vall
full_dict['eigvals']= lamspre
full_dict['meanvec']=mpre

 #%% Merge the dictionaries of testing and training matrices, labels, eigenvectors, eigenvalues, and mean vectorfo
folds_dict = {'trainmats': {},'testmats': {}, 'trainlabels':{}, 'testlabels':{},
              'eigenvectors':{}, 'eigvals':{}, 'meanvecs':{}}

folds_dict['trainmats']= train_folds
folds_dict['testmats']=test_folds
print(folds_dict)



    
    #%%   Start with one example-- last Atrain and Atest
    
    neigs=100
    trainmat = project_cells(Atrain, V, neigs)
    testmat = project_cells(Atest, V,neigs)
    # find distance between columns of testmat and columns of training mat
    #distmat = ordered_dist(trainmat, testmat) 
    # For each testing cell we have two columns- one is the distance to the 
    # nearest to furthest training cell, and the second is the index of the training cell matrix
    # (trainmat) that that corresponds to.
    
    # We can then use this distance and the indices to come up with the class estimates using a few different distance etrics
    # Start with Euclidean distance and nearest neighbor classification


#%% Now you have the Eigenspace, defined by Vnorm, which contains normalized eigenvectors of
# the length of the number of genes, in order of importance.

# Now we need to do the projection of the individual cells (in A) onto the eigenvector space

#Mcol = reshape(M, [x,1]);
#b = (Mcol-m)'*V; %(1 by x) *(x*k) Gives a vector of length k

# Project onto eigenspace
#recon = V*b'; % (x by k) * (k *1) Gives x by 1

# Declare the number of eigenvectors we want to use from the Vnorm matrix 
neigs = 100
#eventually we will loop through the number of eigenvectors and assess the accuracy as afunction of neigs
trainmat = np.zeros((neigs, ntrain))
trainvec= np.zeros((neigs,1))

# First make your Omat which contains the coordinates (the columns of each training image
# and has the corresponding assignment (sens or res) ))
for i in range(ntrain):
    Mcol = Xmat[:,i] # start right after the last training column in Amat 
    McolT= Mcol.T # to get the first testing cell as a row
    # 1xngene x ngene x ntrain = 1 x ntrain
    b=np.matmul(McolT,Vnorm)

    trainvec= b[:neigs]
    trainmat[:,i]= trainvec
#%% Should generate a neigs x ntrain matrix where each column is the coordinates
    # of the cell in eigenspace
print(trainmat)

print(labelsdfpre[:ntrain])
trainlabels = labelsdfpre[:ntrain]
testlabels = labelsdfpre[ntrain:]

#%% Project the testing cells into eigenspace
ntest = ncells-ntrain
Xte = Atest.sub(m, axis=0) # we subtract from each column
Xtest = Xte.as_matrix()
testmat = np.zeros((neigs, ntest))
testvec = np.zeros((neigs,1))
for i in range(ntest):
    Mcol = Xtest[:,i]
    McolT = Mcol.T
    b = np.matmul(McolT, Vnorm)
    testvec = b[:neigs]
    testmat[:,i] = testvec

#%% Compare the coordinates of the testing set with the coordinates of the training set cells
# Now you have a test mat which has class labels corresponding to each column
# For each column of your test mat, you want to find the k training mat columns it is closest to
# make a matrix that stores the ordered indices with the top being the lowest and the bottom the highest
# Euclidean distance
ordered_inds = np.zeros((ntrain, ntest))
dist = np.zeros((ntrain,ntest))
#for i in range(ntest):
for i in range(ntest):
    testvec = testmat[:,i]
    #for j in range(ntrain):
    for j in range(ntrain):
        trainvec = trainmat[:,j]
        # testvec of length neigs and train vec of length neigs
        # find Euclidean distance between the two vectors
        dist_j = np.linalg.norm(testvec-trainvec)
        # fill in your distance vector
        dist[j, i]=dist_j

#%% Now you have a distance matrix where each column is a testing cell
# for each column, we want to output the indices of the training vector distance
# in order from least to greatest
lab_est = [None]*ntest
for i in range(ntest):
    distcol = dist[:,i]
    ind = np.argsort(distcol)
    ordered_inds[:,i] = ind

    
#%% Use the ordered ind to make sense/res matrix
    # Using k=1 nearest neighbors classifier. Need to figure out how to extend this 

for i in range(ntest):
    index = ordered_inds[0,i]
    lab_est[i] = trainlabels[int(index)]
    print(trainlabels[int(index)])
#%%  Make a data frame with your actual and predicted classes
 df = pd.DataFrame({'actual class':testlabels,'predicted class':lab_est})

cnf_matrix = pd.crosstab(df['predicted class'],df['actual class'])

hm = sns.heatmap(cnf_matrix,annot=True,fmt='d',robust=True,
            linewidths=0.1,linecolor='black')   
#%% Calculate some metrics of accuracy from the confusion matrix
TPR=cnf_matrix.iloc[0,0]/ (sum(cnf_matrix.iloc[:,0]))
print(TPR)
TNR = cnf_matrix.iloc[1,1]/(sum(cnf_matrix.iloc[:,1]))
print(TNR)
PPV = cnf_matrix.iloc[0,0]/sum(cnf_matrix.iloc[0,:])
print(PPV)
NPV = cnf_matrix.iloc[1,1]/sum(cnf_matrix.iloc[1,:])
print(NPV)
Acc = (cnf_matrix.iloc[0,0]+ cnf_matrix.iloc[1,1,])/ntest
print(Acc)
prevalance = sum(cnf_matrix.iloc[:,0])/ntest
print(prevalance)
#%% Vary the number of eigenvectors and then compare the accuracy metrics
neigvec = [10, 50, 75, 100, 150, 200, 250,]
TPRi=np.zeros((len(neigvec),1))
TNRi=np.zeros((len(neigvec),1))
PPVi=np.zeros((len(neigvec),1)) 
NPVi=np.zeros((len(neigvec),1)) 
Acci=np.zeros((len(neigvec),1))


#%%
for k in range(len(neigvec)):
    neigi=neigvec[k] # set the number of eigenvectors
    
    trainmat = np.zeros((neigi, ntrain))
   
    trainvec= np.zeros((neigi,1))
    
# First make your Omat which contains the coordinates (the columns of each training image
# and has the corresponding assignment (sens or res) ))
    for i in range(ntrain):
        Mcol = Xmat[:,i] # start right after the last training column in Amat 
        McolT= Mcol.T # to get the first testing cell as a row
        # 1xngene x ngene x ntrain = 1 x ntrain
        b=np.matmul(McolT,Vnorm)

        trainvec= b[:neigi]
        trainmat[:,i]= trainvec
        # We want to save all of these reduced coordinate spaces
        #eigenmatrices.append(trainmat)
    #Project the testing cells into eigenspace
    testmat = np.zeros((neigi, ntest))
    testvec = np.zeros((neigi,1))
    for i in range(ntest):
        Mcol = Xtest[:,i]
        McolT = Mcol.T
        b = np.matmul(McolT, Vnorm)
        testvec = b[:neigi]
        testmat[:,i] = testvec
# Compare the coordinates of the testing set with the coordinates of the training set cells


    dist = np.zeros((ntrain,ntest))
    #for i in range(ntest):
    for i in range(ntest):
        testvec = testmat[:,i]
        #for j in range(ntrain):
        for j in range(ntrain):
            trainvec = trainmat[:,j]
        # testvec of length neigs and train vec of length neigs
        # find Euclidean distance between the two vectors
            dist_j = np.linalg.norm(testvec-trainvec)
        # fill in your distance vector
            dist[j, i]=dist_j
    # Now you have a distance matrix where each column is a testing cell
    # for each column, we want to output the indices of the training vector distance
    # in order from least to greatest
    lab_est = [None]*ntest
    for i in range(ntest):
        distcol = dist[:,i]
        ind = np.argsort(distcol)
        ordered_inds[:,i] = ind
    # Use the ordered ind to make sense/res matrix
    # Using k=1 nearest neighbors classifier. Need to figure out how to extend this 

    for i in range(ntest):
        index = ordered_inds[0,i]
        lab_est[i] = trainlabels[int(index)]

    dfi = pd.DataFrame({'actual class':testlabels,'predicted class':lab_est})

    cnf_matrixi = pd.crosstab(dfi['predicted class'],dfi['actual class'])

    hm = sns.heatmap(cnf_matrixi,annot=True,fmt='d',robust=True,
            linewidths=0.1,linecolor='black')   
# Calculate some metrics of accuracy from the confusion matrix
    TPRi[k] = cnf_matrixi.iloc[0,0]/ (sum(cnf_matrixi.iloc[:,0]))

    TNRi[k] = cnf_matrixi.iloc[1,1]/(sum(cnf_matrixi.iloc[:,1]))

    PPVi[k] = cnf_matrixi.iloc[0,0]/sum(cnf_matrixi.iloc[0,:])

    NPVi[k] = cnf_matrixi.iloc[1,1]/sum(cnf_matrixi.iloc[1,:])

    Acci[k] = (cnf_matrixi.iloc[0,0]+ cnf_matrixi.iloc[1,1,])/ntest

#%% 
print(TPRi)
print(cnf_matrixi)
print(k)
print(cnf_matrixi.iloc[0,0]/ (sum(cnf_matrixi.iloc[:,0])))

#%% Plot the accuracy metrics versus number of eigenvectors
import matplotlib.pyplot as plt


plt.plot(neigvec, TPRi, 'b-', label='TPR')
plt.plot(neigvec, TNRi, 'g-', label='TNR')
plt.plot(neigvec, PPVi, 'r-', label='PPV')
plt.plot(neigvec, NPVi, 'y-', label='NPV')
plt.plot(neigvec, Acci, 'k-', label='Accuracy')
plt.xlabel('Number of eigenvectors')
plt.ylabel('Metrics of Accuracy')
plt.legend()
loc= 'lower right'


#%% BIG GAP FOR Applying the classifier to the post treatment samples!
# Here we are doing it to just the 107 Aziz sample (one of the very late post treatment samples). 
# We will combine this with the other sample, and we should also do the 30 hour time point.
    
    
    
    
    
    
 #%% Project the post treatment cells into eigenspace
dfpost107=dfpost[dfpost['sample'].str.contains("Aziz")]
mpost = Apost.mean(axis =1)
print(npost)
#%%
neigs=100
X107 = Apost.sub(mpost, axis=0) # we subtract from each column
Xpost = X107.as_matrix()
postmat = np.zeros((neigs, npost))
postvec = np.zeros((neigs,1))
for i in range(npost):
    Mcol = Xpost[:,i]
    McolT = Mcol.T
    b = np.matmul(McolT, Vnorm)
    postvec = b[:neigs]
    postmat[:,i] = postvec

#%% Compare the coordinates of the testing set with the coordinates of the training set cells
# Now you have a test mat which has class labels corresponding to each column
# For each column of your test mat, you want to find the k training mat columns it is closest to
# make a matrix that stores the ordered indices with the top being the lowest and the bottom the highest
# Euclidean distance
ordered_inds = np.zeros((ntrain, npost))
dist = np.zeros((ntrain,npost))
#for i in range(ntest):
for i in range(npost):
    postvec = postmat[:,i]
    #for j in range(ntrain):
    for j in range(ntrain):
        trainvec = trainmat[:,j]
        # testvec of length neigs and train vec of length neigs
        # find Euclidean distance between the two vectors
        dist_j = np.linalg.norm(postvec-trainvec)
        # fill in your distance vector
        dist[j, i]=dist_j

#%% Now you have a distance matrix where each column is a testing cell
# for each column, we want to output the indices of the training vector distance
# in order from least to greatest
lab_est_post = [None]*npost
for i in range(npost):
    distcol = dist[:,i]
    ind = np.argsort(distcol)
    ordered_inds[:,i] = ind

    
#%% Use the ordered ind to make sense/res matrix
    # Using k=1 nearest neighbors classifier. Need to figure out how to extend this 

for i in range(npost):
    index = ordered_inds[0,i]
    lab_est_post[i] = trainlabels[int(index)]
#%% Now that you have lab_est post, quantify the proportion.
ct_sens=0
for i in range(npost):
    if lab_est_post[i]=='sens':
        ct_sens+=1

phi_est=ct_sens/npost
print(phi_est)
phi_est0=1-prevalance
print(phi_est0)
    
    
    
    