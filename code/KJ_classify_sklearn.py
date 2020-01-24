"""
Created on Wed Jul  3 12:22:02 2019

@author: kj22643
"""
%reset

# This script performs k-fold cross validation on the cells labeled resistant
# and sensitive by lineage abundance changes.

# This also used to optimize the number of principal components and the number
# of nearest neighbors for the classifier when it's applied to the remainder of 
# the pre=treatment, intermediate, and post-treatment samples


import numpy as np
import pandas as pd
import os
import scanpy as sc
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import plot, show, draw, figure, cm
import matplotlib as plt
import random
from collections import OrderedDict
import copy
import math
import matplotlib.pyplot as plt
from pandas import DataFrame, Series
import plotnine as gg
import scipy as sp
import scipy.stats as stats
import sklearn as sk
import sklearn.model_selection as model_selection
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import StratifiedKFold
import sklearn.feature_selection as feature_selection
import sklearn.linear_model as linear_model
import sklearn.pipeline as pipeline
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cluster import KMeans
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn import svm
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
os.chdir('/Users/kj22643/Documents/Documents/231_Classifier_Project/code')
from func_file import find_mean_AUC
from func_file import find_mean_AUC_SVM



path = '/Users/kj22643/Documents/Documents/231_Classifier_Project/data'
#path = '/stor/scratch/Brock/231_10X_data/'
os.chdir(path)
sc.settings.figdir = 'KJ_plots'
sc.set_figure_params(dpi_save=300)
sc.settings.verbosity = 3

#%% Load in pre and post treatment 231 data
adata = sc.read('post-cell-cycle-regress.h5ad')
adata.obs.head()
# current samples:
#BgL1K
#30hr
#Rel-1
#Rel-2
# We will change these to time points
# Count the number of unique lineages in all the cells
uniquelins = adata.obs['lineage'].unique()
nunique = len(uniquelins)

#%% Identify lineages that have been recalled from the pre-treatment sample
# Make a column labeled recalled lin that can be use to identify the specific lineages of interest

# These correspond to AA161 and AA170
# AA170 is the highly proliferative lineage that ultimately decreases in lineage abundance
# AA161 is a lowly abundant in the pre-treatment but is one of the ones that comes
# to  be highly abundance in both post treatment samples

reslin = ['GTACATTTTCATACCTCTCG']
senslin = ['GTGTGTTGTTTTTGTACGGC']


adata.obs.loc[adata.obs.lineage.isin(reslin)==True,'recalled_lin'] = 'res'
adata.obs.loc[adata.obs.lineage.isin(reslin)==False,'recalled_lin'] = 'na'
adata.obs.loc[adata.obs.lineage.isin(senslin)==True,'recalled_lin'] = 'sens'

print(adata.obs['recalled_lin'])

#%% Add a column to rename the samples by time point
samps= adata.obs['sample'].unique()

timepoint = np.array(['t=0hr', 't=30hr', 't=1344hr'])

adata.obs.loc[adata.obs['sample']==samps[0], 'timepoint']='t=0hr'
adata.obs.loc[adata.obs['sample']==samps[1], 'timepoint']='t=30hr'
adata.obs.loc[adata.obs['sample']==samps[2], 'timepoint']='t=1344hr'
adata.obs.loc[adata.obs['sample']==samps[3], 'timepoint']='t=1344hr'

print(adata.obs['timepoint'].unique())



#%% Separately make dataframes for the pre-treatment, intermediate, and post treatment samples
# t=0 hr (pre-treatment), 3182 pre treatment cells
# We want to keep the info about the lineage so we can potentially
# use it to make evenly divided testing and training data sets
dfall = pd.concat([adata.obs['lineage'],adata.obs['timepoint'], 
               pd.DataFrame(adata.raw.X,index=adata.obs.index,
                            columns=adata.var_names),], axis=1) 
ncells = len(dfall)
print(ncells)
adata_pre = adata[adata.obs['timepoint']=='t=0hr', :]
dfpre = pd.concat([adata_pre.obs['lineage'], adata_pre.obs['recalled_lin'],
               pd.DataFrame(adata_pre.raw.X,index=adata_pre.obs.index,
                            columns=adata_pre.var_names),], axis=1) 
npre = len(dfpre)
print(npre)
# t = 30 hr (intermediate timepoint) 5169 int treatment cells
adata_int = adata[adata.obs['timepoint']=='t=30hr', :]
dfint = pd.concat([adata_int.obs['lineage'], adata_int.obs['recalled_lin'],
                   pd.DataFrame(adata_int.raw.X, index=adata_int.obs.index, 
                                columns = adata_int.var_names),], axis=1)
nint = len(dfint)
print(nint)
# t=1344 hr (~roughly 8 weeks), 10332 post treatment cells
adata_post = adata[adata.obs['timepoint']=='t=1344hr', :]
dfpost = pd.concat([adata_post.obs['lineage'], adata_post.obs['recalled_lin'],
                    pd.DataFrame(adata_post.raw.X, index=adata_post.obs.index, 
                                 columns = adata_post.var_names),],axis =1)
npost = len(dfpost)
print(npost)
#%% Try to add a .obs column that records lineage abundance from the different samples
linAbundpre= adata_pre.obs['lineage'].value_counts()
linAbundint = adata_int.obs['lineage'].value_counts()
linAbundpost = adata_post.obs['lineage'].value_counts()

# Start by adding the linabundpre and lin abund post to the pre-treatment data frame
df1 = pd.DataFrame(linAbundpre)
df1['linabundpre']= df1.lineage
df1=df1.drop(['lineage'], axis=1)
df1['lineage'] = df1.index
df1=df1.drop(index='nan')


df2 = pd.DataFrame(linAbundpost)
df2['linabundpost']= df2.lineage
df2=df2.drop(['lineage'], axis=1)
df2['lineage'] = df2.index
df2=df2.drop(index='nan')
#%% Merge the linage abundance data frames from the pre and post treatment samples into dfpre
dfpre= pd.DataFrame.merge(df1, dfpre, left_on=['lineage'], 
              right_on=['lineage'], how='right')
dfpre = pd.DataFrame.merge(df2, dfpre, left_on=['lineage'],
              right_on=['lineage'], how='right') 
dfpre['linabundpost'] = dfpre['linabundpost'].fillna(0)
dfpre['linabundpre']= dfpre['linabundpre'].fillna(0)
#%% Make a column that is the logfoldchange from post to pre
dfpre['foldchange'] =  (dfpre['linabundpost']-dfpre['linabundpre'])
#dfpre['foldchange'] =  ((dfpre['linabundpost']/npost)-(dfpre['linabundpre']/npre))/(dfpre['linabundpre']/npre)
print(dfpre['foldchange'].unique())
foldchangevec = dfpre['foldchange']
#%% Look at fold change and log fold change for each cell and its correpsonding lineage
dfpre['logfoldchange'] = np.log(dfpre['foldchange'])
dfpre['logfoldchange']= dfpre['logfoldchange'].fillna(0)


print(dfpre['logfoldchange'].unique())

# Figures of logfold change and fold change
plt.figure()
plt.hist(dfpre['logfoldchange'], range = [3, 7], bins = 100)
plt.xlabel('log(foldchange) of lineage abundance')
plt.ylabel('number of cells')

plt.figure()
plt.hist(dfpre['foldchange'],range = [0, 500], bins = 100)
plt.xlabel('foldchange of lineage abundance')
plt.ylabel('number of cells')



#%% Make the survivor column, but don't label it yet. 
#dfpre['survivor'] =np.where(dfpre['linabundpost'] >1000, 'res','sens'
# Want to call cells that have an increase in lineage abundance resistant
dfpre.loc[dfpre.foldchange>0, 'survivor'] = 'res'

dfpre.loc[dfpre.foldchange<-200, 'survivor']= 'sens'
survivorvec = dfpre['survivor']


# Want to call cells that have an signficant in lineage abundance resistant

# Make a dataframe that only contains cells who are given a sens or resistant label

dfsr = dfpre[(dfpre['survivor']=='sens') | (dfpre['survivor'] == 'res')]
nclass = len(dfsr)
print(nclass)
y= pd.factorize(dfsr['survivor'])[0] 
y ^= 1
mu_sr = sum(y)/len(y)
print(mu_sr)
Xsr = dfsr.drop(columns= [ 'lineage', 'linabundpost', 'recalled_lin', 'linabundpre', 'foldchange', 'logfoldchange', 'survivor'])
ntrain = len(dfsr)


y= pd.factorize(dfsr['survivor'])[0] 
y ^= 1


mu_pre = sum(y)/len(y)
# X is your cell gene matrix, y is your class labels
#%% Run PCA on entire pre-treatment data set
X = Xsr
full_dict = {'fullmat':{}, 'prev':{}, 'labels':{}, 'V':{}, 'lambdas':{}, 'varPC':{}}
n_neighbors = 10
n_components = 500
random_state = 0

pca=PCA(copy=True, iterated_power='auto', n_components=500, random_state=0,
  svd_solver='auto', tol=0.0, whiten=False)
V = pca.fit(X)  
varPC= (pca.explained_variance_ratio_)  

lambdas = pca.singular_values_ 
full_dict['full_mat'] = X
full_dict['labels']=y
full_dict['V']= V
full_dict['varPC']= varPC



# X is your cell gene matrix, y is your class labels
#%%
# Split into train/test
kCV = 5
skf = StratifiedKFold(n_splits=kCV, shuffle= True)
Atrain = {}
Atest = {}
ytest = {}
ytrain = {}
proprestest = {}
proprestrain = {}

folds_dict = {'trainmat':{}, 'trainlabel':{}, 'eigenvectors':{}, 'eigvals':{}, 'meanvec':{}}
for i in range(kCV):    
    for train_index, test_index in skf.split(X, y):
        Atrain[i] = X.iloc[train_index, :]
        Atest[i] = X.iloc[test_index, :]
        ytest[i]= y[test_index]
        ytrain[i]= y[train_index]
        proprestest[i] = sum(ytest[i])/len(ytest[i])
        proprestrain[i] = sum(ytrain[i])/len(ytrain[i])
         
# Save all of your stratified folds into a single dictionary. 
folds_dict['trainmat'] = Atrain
folds_dict['trainlabel']= ytrain
folds_dict['testmat'] = Atest
folds_dict['testlabel'] = ytest
folds_dict['prevtest'] = proprestest
folds_dict['prevtrain']= proprestrain



n_classes = len(np.unique(y))
      


# Use a nearest neighbor classifier to evaluate the methods
knn = KNeighborsClassifier(n_neighbors=n_neighbors)

pca=PCA(copy=True, iterated_power='auto', n_components=n_components, random_state=0,
            svd_solver='auto', tol=0.0, whiten=False)

#%% #Compute the ROC curve for each fold


AUC=np.zeros((kCV,1))
acc = np.zeros((kCV,1))
for i in range(kCV):
    X_train = folds_dict['trainmat'][i]
    y_train = folds_dict['trainlabel'][i]
    X_test = folds_dict['testmat'][i]
    y_test = folds_dict['testlabel'][i]
    
    # Build new model based on training data set
    pca.fit(X_train, y_train)
    V_train = pca.fit_transform(X_train)
# Fit a nearest neighbor classifier on the model built on the training data set
    knn.fit(pca.transform(X_train), y_train)
    prob_scores = knn.predict_proba(pca.transform(X_test))
# Compute the nearest neighbor accuracy on the embedded test set
    acc_knn = knn.score(pca.transform(X_test), y_test)
    acc[i]=acc_knn
    fpr, tpr, thresholds = metrics.roc_curve(y_test, prob_scores[:,1], pos_label=1)
    roc_auc = auc(fpr, tpr)
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve from k-fold CV')
    plt.legend(loc="lower right")
    plt.show()

    AUC[i]= roc_auc
  
plt.figure()
plt.plot(np.linspace(0,4, num=kCV), AUC)
plt.ylim([0.5, 1.05])
plt.xlabel('fold')
plt.ylabel('AUC')
plt.title('AUC for each fold')
plt.show()

meanAUC_PCA = AUC.mean()

#%% Test function call
os.chdir('/Users/kj22643/Documents/Documents/231_Classifier_Project/code')
from func_file import find_mean_AUC

meanAUCtest = find_mean_AUC(folds_dict, n_neighbors, n_components)
print(meanAUCtest)
#%% Run loop to vary the number of neighbors 

n_neighborsvec = np.linspace(1,101,num = 11)
meanAUC_neighbs = np.zeros((len(n_neighborsvec),1))

for i in range(len(n_neighborsvec)):
    meanAUC_neighbs[i]=find_mean_AUC(folds_dict, int(n_neighborsvec[i]), n_components)
#%%  
plt.figure()
plt.plot(n_neighborsvec, meanAUC_neighbs)
plt.xlabel('k nearest neighbors')
plt.ylabel('mean AUC')
plt.title('mean AUC as a function of knn')
plt.show()

#%% Run loop to vary the number of components
n_neighbors = 90
n_compsvec = [2, 3, 5, 10, 20, 50, 100, 250, 500]
meanAUC_comps = np.zeros((len(n_compsvec),1))
for i in range(len(n_compsvec)):
    meanAUC_comps[i]=find_mean_AUC(folds_dict, n_neighbors, n_compsvec[i])
#%%
plt.figure()
plt.plot(n_compsvec, meanAUC_comps)
plt.xlabel('Number of principal components')
plt.ylabel('mean AUC')
plt.title('mean AUC as a function of number of PCs')
plt.show()
#%% Set the number of neighbors and number of components to the two optimal values
# and generate the ROC curve
n_neigbors = 100
n_components = 500

X_train = folds_dict['trainmat'][1]
y_train = folds_dict['trainlabel'][1]
X_test = folds_dict['testmat'][1]
y_test = folds_dict['testlabel'][1]
    
    # Build new model based on training data set
pca.fit(X_train, y_train)
V_train = pca.fit_transform(X_train)
# Fit a nearest neighbor classifier on the model built on the training data set
knn.fit(pca.transform(X_train), y_train)
prob_scores = knn.predict_proba(pca.transform(X_test))
# Compute the nearest neighbor accuracy on the embedded test set
acc_knn = knn.score(pca.transform(X_test), y_test)
acc[1]=acc_knn
fpr, tpr, thresholds = metrics.roc_curve(y_test, prob_scores[:,1], pos_label=1)
roc_auc = auc(fpr, tpr)
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',
        lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve k=90, nPCs = 500')
plt.legend(loc="lower right")
plt.show()

AUC_opt= roc_auc

mean_AUC_opt = find_mean_AUC(folds_dict, n_neighbors, n_components)



#%%Repeat and try with svm, varying the C parameter and the type of basis function


clf = svm.SVC(kernel='poly', C=1000, probability = True)
                          

AUC=np.zeros((kCV,1))
acc = np.zeros((kCV,1))
for i in range(kCV):
    X_train = folds_dict['trainmat'][i]
    y_train = folds_dict['trainlabel'][i]
    X_test = folds_dict['testmat'][i]
    y_test = folds_dict['testlabel'][i]
    
    # Build new model based on training data set
    clf.fit(X_train, y_train)
    SVM_prob_scores= clf.predict_proba(X_test)
    
# Compute the nearest neighbor accuracy on the embedded test set

    fpr, tpr, thresholds = metrics.roc_curve(y_test, SVM_prob_scores[:,1], pos_label=1)
    roc_auc = auc(fpr, tpr)
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve from CV for SVM')
    plt.legend(loc="lower right")
    plt.show()

    AUC[i]= roc_auc
  
plt.figure()
plt.plot(np.linspace(0,4, num=kCV), AUC)
plt.ylim([0.5, 1.05])
plt.xlabel('fold')
plt.ylabel('AUC')
plt.title('AUC for each fold from SVM')
plt.show()

meanAUC_SVM = AUC.mean()

#%% Run a loop to vary the C parameters
os.chdir('/Users/kj22643/Documents/Documents/231_Classifier_Project/code')
from func_file import find_mean_AUC_SVM

Cvec = [0.01, 1, 100, 1000, 10000]
meanAUC_C = np.zeros((len(Cvec),1))
for i in range(len(Cvec)):
    meanAUC_C[i]=find_mean_AUC_SVM(folds_dict, Cvec[i], 'rbf')

#%% Plot mean AUC of varying C parameter
plt.figure()
plt.plot(Cvec, meanAUC_C)
plt.ylim([0.9, 1.05])
plt.xlabel('C parameters')
plt.ylabel('mean AUC')
plt.title('mean AUC as a function of SVM C parameters')
plt.show()

print(max(meanAUC_C))
Copt = Cvec[3]
#%% Using the best C, try  different basis functions for SVM
basisvec = ['rbf', 'linear', 'poly', 'sigmoid']
#%%
meanAUC_basis = np.zeros((len(Cvec),1))
for i in range(len(Cvec)):
    meanAUC_basis[i]=find_mean_AUC_SVM(folds_dict, Copt, basisvec[i])
#%% Figure out how to plot these categorical results in a bar graph!
plt.style.use('ggplot')
x_pos = [i for i, _ in enumerate(basisvec)]   
AUCvals = [np.asscalar(meanAUC_basis[0]), np.asscalar(meanAUC_basis[1]), np.asscalar(meanAUC_basis[2]), np.asscalar(meanAUC_basis[3])]

test = [ 3, 4, 5]
print(AUCvals)
#%%
plt.figure()
plt.bar(x_pos, AUCvals)
plt.xlabel('type of basis function')
plt.ylabel('mean AUC')
plt.title('mean AUC as a function of SVM basis function parameters')
plt.xticks(x_pos, basisvec[:4])
plt.ylim([0.9, 0.95])
plt.show()


# Looks like rbf is the best and C=1000 is the best 

#%%
from sklearn import metrics
from sklearn.metrics import roc_curve, auc

fpr, tpr, thresholds = metrics.roc_curve(labstest, prob_scores[:,0], pos_label=1)
roc_auc = auc(fpr, tpr)

plt.style.use('ggplot')
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve from Pre-treatment Sample')
plt.legend(loc="lower right")
plt.show()

# Now you have 0 for sens 1 for res labels (labstest) 
# and the corresponding probabilities (first column of the prob_scores)
# Generate an ROC curve





#%%
# Embed the data set in 2 dimensions using the fitted model
X_embedded = pca.transform(X)
 # Plot the projected points and show the evaluation score
plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c = labels, s=30, cmap='Set1')
plt.title("{}, KNN (k={})\nTest accuracy = {:.2f}".format('PCA',
                                                              n_neighbors,
                                                              acc_knn))
plt.show()



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