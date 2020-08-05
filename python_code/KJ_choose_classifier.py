# This script finds the optimal classifier (PCA+KNN vs. Linear SVM) and 
# justifies the use of the threshold for calling an unknown cell S or R. 

%reset

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
import csv
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cluster import KMeans
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn import svm
from sklearn.svm import LinearSVC      
from sklearn.calibration import CalibratedClassifierCV
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
os.chdir('/Users/kj22643/Documents/Documents/231_Classifier_Project/code')
from func_file import find_mean_AUC
from func_file import find_mean_AUC_SVM


# Set path to the folder containing the scRNAseq data
path = '/Users/kj22643/Documents/Documents/231_Classifier_Project/data'
#path = '/stor/scratch/Brock/231_10X_data/'
os.chdir(path)
sc.settings.figdir = 'KJ_plots'
sc.set_figure_params(dpi_save=300)
sc.settings.verbosity = 3

#%% Load in normalized scRNAseq data from all time points
adata = sc.read('post-cell-cycle-regress.h5ad')
adata.obs.head()
# Note: adata object contains additional column of lineage abundance info. 
# this will be used to identify sensitive and resistant cells in the pre-
# treatment time point. Also contains column for the sample number, which
# corresponds to the time point
uniquelins = adata.obs['lineage'].unique()
nunique = len(uniquelins)

# Name samples by their time point  (0, 7, 10 weeks and pre/post)
samps= adata.obs['sample'].unique()

timepoint = np.array(['t=0 wks', 't=7 wks', 't=10 wks'])
adata.obs.loc[adata.obs['sample']==samps[0], 'timepoint']='t=0 wks'
adata.obs.loc[adata.obs['sample']==samps[1], 'timepoint']='t=7 wks'
adata.obs.loc[adata.obs['sample']==samps[2], 'timepoint']='t=10 wks'

TP = np.array(['pre', 'post'])
adata.obs.loc[adata.obs['sample']==samps[0], 'TP']='pre'
adata.obs.loc[adata.obs['sample']==samps[1], 'TP']='post'
adata.obs.loc[adata.obs['sample']==samps[2], 'TP']='post'

# Plot all data in a umap colored by timepoint
sc.pl.umap(adata,color='timepoint',palette=['#f79d02', '#d604f7', '#00c6c6'])

  
#%% Calculate the pre and post treatment lineage abundances
# lineage counts = the number of cells in a lineage
# lineage abundance = the proportion of a population the cells in that lineage
# take up at a specific time point.
# Calculate the magnitude of change in lineage abundance
# change_i = abund_post_i - abund_pre_i


lincountspre= adata.obs.loc[adata.obs.TP=='pre','lineage'].value_counts()
lincountspost= adata.obs.loc[adata.obs.TP=='post','lineage'].value_counts()

npre = sum(lincountspre)
npost=sum (lincountspost)
linabundpre = lincountspre/npre
linabundpost = lincountspost/npost
# Put into a dictionary to compare the lineage abundances at each time point
dpost = linabundpost.to_dict()
dpre = linabundpre.to_dict()
# dchange = the change in lineage abundance from pre to post of an individual lineage
# key = lineage, val = change in lineage abundance
dchange = {key: dpost[key] - dpre.get(key, 0) for key in dpost.keys()}
# save this dictionary into the anndata object
adata.uns['linabundchange'] = dchange
linabundchange = adata.uns['linabundchange']
# output these lists of lineage abundances to csv files
a_file = open("linabundchange.csv", "w")
writer = csv.writer(a_file)
for key, value in linabundchange.items():
    writer.writerow([key, value])
    
a_file.close()
linabundpre.to_csv("linabundpre.csv")
filename = "linabundpre.csv"
linabundpost.to_csv("linabundpost.csv")
filename = "linabundpost.csv"
#%%Map the lineages in dchange to the cells in the anndata object
# each cell now is labeled by the amount its lineage changes post treatment
adata.obs['linabundchange']= adata.obs['lineage'].map(dchange)
# Set the threshold for calling cells sensitive or resistant, pulling from only
# the lineages with greatest increase (R) or decrease (S) in abundance 
S_threshold = 0.05
R_threshold = 0
# Plot the change in lineage abundance distribution                                       
plt.figure()
plt.hist(adata.obs.loc[adata.obs.timepoint == 't=0 wks','linabundchange'], bins = 100)
plt.plot([R_threshold,R_threshold], [0,900], c='r', alpha = 1)
plt.plot([-S_threshold, -S_threshold],[0,900], c='g', alpha = 1)
plt.xlabel(' Change in lineage abundance')
plt.ylabel('Number of cells in the lineage')
plt.title('Distribution of lineage abundance change')
#plt.legend(loc=1, prop={'size': 15})     
plt.ylim(0, 900)
plt.xlim(-0.5, 0.5)
plt.grid(b=None)
                                    

#%% Make the sensitive and resistant labels within your anndata object

classLabel = np.array(['res', 'sens', 'unknown', 'res_est', 'sens_est'])
adata.obs.loc[adata.obs.timepoint=='t=0 wks','classLabel'] = 'unknown'
adata.obs.loc[adata.obs.timepoint=='t=7 wks','classLabel'] = 'unknown'
adata.obs.loc[adata.obs.timepoint=='t=10 wks','classLabel'] = 'unknown'
adata.obs.loc[(adata.obs.linabundchange>R_threshold)&(adata.obs.timepoint=='t=0 wks'), 'classLabel'] = 'res'
adata.obs.loc[(adata.obs.linabundchange<-S_threshold)&(adata.obs.timepoint=='t=0 wks'), 'classLabel'] = 'sens'
adata.obs.loc[adata.obs.lineage=='nan', 'classLabel'] = 'unknown'
#adata.obs.loc[(adata.obs.linpropchange<-0.5)&(adata.obs.timepoint=='t=0hr'), 'classLabel'] = 'sens'

# plot the labeled cells alongside unknown cells and alone on same axes.
sc.pl.umap(adata,color='classLabel', palette = ['red', 'green', 'gray'])
sc.pl.umap(adata,color='classLabel', palette = ['red', 'green', 'white'])


#%% Now we want to take out only the data that is labeled as res and sens 
# and test the classifier on the labeled dataset
adata_sr= adata[(adata.obs['classLabel']=='res')|(adata.obs['classLabel']=='sens'), :]
dfsr= pd.concat([adata_sr.obs['lineage'],adata_sr.obs['classLabel'],
               pd.DataFrame(adata_sr.raw.X,index=adata_sr.obs.index,
                            columns=adata_sr.var_names),], axis=1) 
# plot only the labeled cells
sc.pl.umap(adata_sr,color='classLabel',palette=['red', 'green']) 

nsr= len(dfsr)
y= pd.factorize(dfsr['classLabel'])[0] 
nres = sum(y)
mu_sr = sum(y)/len(y) # shows how we have unbalanced class sizes
print(mu_sr)
# Xsr is your cell gene matrix, y is your class labels
Xsr = dfsr.drop(columns= [ 'classLabel', 'lineage'])



# %%Assign the parameters for PCA+KNN and Linear SVM and fit the models
# PCA hyperparameters optimized in KJ_find_hyperparameters.py 

# SET UP PCA
n_neighbors = 73
n_components = 500

knn = KNeighborsClassifier(n_neighbors=n_neighbors)

pca=PCA(copy=True, iterated_power='auto', n_components=n_components, random_state=0,
            svd_solver='auto', tol=0.0, whiten=False)
# FIT PCA
pca.fit(Xsr, y)
# the components is a n_components x n-total genes matrix that gives the weight of 
# each gene that goes into making the principal component.
components = pca.components_
compdf = pd.DataFrame(components, columns = adata_sr.var_names)
adata.uns['components'] = compdf
# We cann use the compdf to make the arrow plot of the gene weights in PC space 
V = pca.fit_transform(Xsr) 
PCsr = pca.transform(Xsr)
var_in_PCs= pca.explained_variance_ratio_
cdf_varPCs = var_in_PCs.cumsum()

# SET UP SVM
Copt = 1
clf = svm.LinearSVC( C=Copt)
linear_svc = LinearSVC()

# FIT SVM
clf.fit(Xsr, y)
clf = linear_svc.fit(Xsr, y)
calibrated_svc = CalibratedClassifierCV(base_estimator = linear_svc, cv= "prefit")
calibrated_svc.fit(Xsr, y)
predicted = calibrated_svc.predict(Xsr)
prob = calibrated_svc.predict_proba(Xsr)
weights = clf.coef_
weightsdf = pd.DataFrame(weights, columns = adata_sr.var_names)
adata.uns['weights'] = weightsdf



#%% Split into train/test, fit each classifier to the training data, and test
# accuracy on testing data
kCV = 5
skf = StratifiedKFold(n_splits=kCV, shuffle= True) 
Atrain = {}
Atest = {}
ytest = {}
ytrain = {}
proprestest = {}
proprestrain = {}
X = Xsr

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


#%% #Compute the ROC curve for each fold


AUC_PCA=np.zeros((kCV,1))
AUC_SVM = np.zeros((kCV,1))
acc_PCA = np.zeros((kCV,1))
acc_SVM = np.zeros((kCV,1))
for i in range(kCV):
    X_train = folds_dict['trainmat'][i]
    y_train = folds_dict['trainlabel'][i]
    X_test = folds_dict['testmat'][i]
    y_test = folds_dict['testlabel'][i]
    
    # Build new model based on training data set
    pca.fit(X_train, y_train)
# Fit a nearest neighbor classifier on the model built on the training data set
    knn.fit(pca.transform(X_train), y_train)
    prob_scores = knn.predict_proba(pca.transform(X_test))
# Compute the nearest neighbor accuracy on the embedded test set
    acc_PCAi = knn.score(pca.transform(X_test), y_test)
    acc_PCA[i]=acc_PCAi
    fpr, tpr, thresholds = metrics.roc_curve(y_test, prob_scores[:,1], pos_label=1)
    roc_auc = auc(fpr, tpr)

    AUC_PCA[i]= roc_auc
    
    clf = linear_svc.fit(X_train, y_train)
    calibrated_svc = CalibratedClassifierCV(base_estimator = linear_svc, cv= "prefit")
    calibrated_svc.fit(X_train, y_train)
    predicted = calibrated_svc.predict(X_test)
    prob = calibrated_svc.predict_proba(X_test)
    acc_SVMi = clf.score(X_test, y_test)
    acc_SVM[i] = acc_SVMi
    fpr, tpr, thresholds = metrics.roc_curve(y_test, prob[:,1], pos_label=1)
    roc_auc = auc(fpr, tpr)
    
    AUC_SVM[i] = roc_auc
    
    

 #%% plot AUC for each classifier 
 
 
meanAUC_PCA = AUC_PCA.mean()
meanAUC_SVM = AUC_SVM.mean()
meanacc_PCA = acc_PCA.mean()
meanacc_SVM = acc_SVM.mean()

plt.figure()
plt.plot(np.linspace(1,5, num=kCV), AUC_PCA,  color='darkorange',
        label='PCA+KNN = %0.2f' % meanAUC_PCA)
plt.plot(np.linspace(1,5,num=kCV), AUC_SVM,color='navy',
         label = 'Linear SVM = %0.2f' % meanAUC_SVM, linestyle='--')
plt.ylim([0.5, 1.05])
plt.xlabel('fold')
plt.ylabel('AUC')
plt.title('AUC for each fold')
plt.legend(loc="lower right")
plt.grid(b=None)
plt.show()

plt.figure()
plt.plot(np.linspace(1,5, num=kCV), acc_PCA,  color='darkorange',
         label='PCA+KNN = %0.2f' % meanacc_PCA)
plt.plot(np.linspace(1,5,num=kCV), acc_SVM,color='navy',
         label = 'Linear SVM = %0.2f' % meanacc_SVM, linestyle='--')
plt.ylim([0.5, 1.05])
plt.xlabel('fold')
plt.ylabel('Accuracy')
plt.title('Accuracy in prediction for each fold')
plt.legend(loc="lower right")
plt.grid(b=None)
plt.show()

# Results indicate that the better classifier is linear SVM. Now see what 
# probability thresholds enable accuracy in labeled data set



#%% Generate ROC curve to determine the best threshold for deciding on class

# FOR PCA
n_classes = len(np.unique(y))
      
# Fit a nearest neighbor classifier on the model built on the training data set
knn.fit(pca.transform(Xsr), y)
prob_scores = knn.predict_proba(pca.transform(Xsr))
# Compute the nearest neighbor accuracy on the embedded test set
acc_knn = knn.score(pca.transform(Xsr), y)
fpr, tpr, thresholds = metrics.roc_curve(y, prob_scores[:,1], pos_label=1)
roc_auc = auc(fpr, tpr)

plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC curve labeled data PCA + KNN')
plt.legend(loc="lower right")
plt.grid(b=None)
plt.show()

#%% FOR SVM
clf = linear_svc.fit(Xsr, y)
calibrated_svc = CalibratedClassifierCV(base_estimator = linear_svc, cv= "prefit")
calibrated_svc.fit(Xsr, y)
predicted = calibrated_svc.predict(X)
prob = calibrated_svc.predict_proba(X)


fpr, tpr, thresholds = metrics.roc_curve(y, prob[:,1], pos_label=1)
roc_auc = auc(fpr, tpr)

plt.figure()
lw = 2
plt.plot(fpr, tpr, color='navy',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve labeled data Linear SVM')
plt.legend(loc="lower right")
plt.grid(b=None)
plt.show()


# This tells us that the classes are very separable, and so really any reasonable
# threshold should be able to separate the classes.
# This completes the goal of this script- to determine the classifier and 
# justify the use of any threshold.












# %% Apply SVM to data,
adata_unk= adata[(adata.obs['classLabel']=='unknown'), :]
# Make a data frame with just the raw gene-cell matrix but keep the cell identity
Xunk = pd.DataFrame(adata_unk.raw.X,index=adata_unk.obs.index,
                            columns=adata_unk.var_names)

clf = linear_svc.fit(Xsr, y)
calibrated_svc = CalibratedClassifierCV(base_estimator = linear_svc, cv= "prefit")
calibrated_svc.fit(Xsr, y)
predicted = calibrated_svc.predict(Xunk)
prob = calibrated_svc.predict_proba(Xunk) # output a probability of being in each class



#%% Vary the threshold probability 
thres_prob = 0.9
B = prob[:,1]>thres_prob
y_est=B*1

#%%
#  Try mapping outcomes

# first make an indexed data frame
class_est= pd.DataFrame(y_est, index = adata_unk.obs.index)
class_est.columns=['est']
adata.obs['class_est'] = adata.obs.index.map(class_est.est)
sc.pl.umap(adata,color=['class_est'])
adata.obs.loc[adata.obs.class_est==1, 'classLabel'] = 'res'
adata.obs.loc[adata.obs.class_est==0, 'classLabel'] = 'sens'


#
sc.pl.umap(adata,color=['classLabel'])


adata_pre = adata[adata.obs['timepoint']=='t=0 wks', :]
dfpre = pd.concat([adata_pre.obs['classLabel'],
               pd.DataFrame(adata_pre.raw.X,index=adata_pre.obs.index,
                            columns=adata_pre.var_names),], axis=1) 
lin_list_pre = adata_pre.obs.lineage


#sc.pl.umap(adata_pre,color='classLabel', palette = ['red', 'green']) 
npre = len(dfpre)
print(npre)
ypreb = dfpre.classLabel=='res'
ypre= ypreb*1
phirpre = (sum(ypre))/npre
phispre = 1-phirpre
print(phispre)

adata_post1 = adata[adata.obs['timepoint']=='t=7 wks', :]
dfpost1= pd.concat([adata_post1.obs['classLabel'],
               pd.DataFrame(adata_post1.raw.X,index=adata_post1.obs.index,
                            columns=adata_pre.var_names),], axis=1) 
    
npost1 = len(dfpost1)
print(npost1)
ypost1b = dfpost1.classLabel=='res'
ypost1= ypost1b*1
phirpost1 = (sum(ypost1))/npost1
phispost1 = 1-phirpost1
print(phispost1)


adata_post2 = adata[adata.obs['timepoint']=='t=10 wks', :]
dfpost2= pd.concat([adata_post2.obs['classLabel'],
               pd.DataFrame(adata_post2.raw.X,index=adata_post2.obs.index,
                            columns=adata_post2.var_names),], axis=1) 
    
npost2 = len(dfpost2)
print(npost2)
ypost2b = dfpost1.classLabel=='res'
ypost2= ypost2b*1
phirpost2 = (sum(ypost2))/npost2
phispost2 = 1-phirpost2
print(phispost2)

#%% Look at the weights of each gene! And compare the values of each
weightsdf = weightsdf.abs().sort_values(by= [0], axis =1,ascending = False)
weights1 = weightsdf.iloc[0,:]
indtop = weights1.index[0:50]
sc.pl.matrixplot(adata, indtop, groupby = 'classLabel', swap_axes = True, figsize = [5, 10])
#standard_scale = 'var',

#%%
x = -np.abs(weights[0,:])
ind = np.unravel_index(np.argsort(x, axis=None), x.shape)
weights1 = weightsdf.iloc[0,:]
ordered_weights = weights1.iloc[(-np.abs(weights[1,:]).argsort())]
indtop = ordered_weights.index[0:50]
sc.pl.matrixplot(adata, indtop, groupby = 'classLabel', swap_axes = True)
#%%


ordered_weights=weightsdf.iloc[(-np.abs(weights[0,:])).argsort()]
x = -np.abs(weights[0,:]);
ind = np.unravel_index(np.argsort(x, axis=None), x.shape)


indtop = ordered_comp1.index[0:50]

ordered_comp2=PC2.iloc[(-np.abs(components[1,:])).argsort()]
x2 = -np.abs(components[1,:]);
ind2 = np.unravel_index(np.argsort(x2, axis=None), x2.shape)


indtop = ordered_comp1.index[0:50]
indtop2 = ordered_comp2.index[0:50]
compdf1 = compdf.index==0
compdf2 = compdf.index==1
comps1 = compdf[compdf1]
comps2 = compdf[compdf2]










#%% Fit a nearest neighbor classifier on the model built on the training data set

# Make a nearest neighbor classifier from the down sampled sensitive cells
#knn.fit(pca.transform(Xsmall), ysmall)

knn.fit(pca.transform(Xsr), y)

thres_prob = 0.15
# Apply the pca and knn classifier to the pre, int, and post1&2 treatment samples. 
y_premin = knn.predict(pca.transform(Xpremin))
pre_prob= knn.predict_proba(pca.transform(Xpremin))
B = pre_prob[:,1]>thres_prob
mu_pre = (sum(pre_prob[:,1]) + sum(y))/npre
y_pre = B*1
mu_pre_PCA = (sum(y_pre)+sum(y))/(len(Xpremin) + len(Xsr))
sigmasq_pre_PCA = mu_pre_PCA*(1-mu_pre_PCA)/(len(Xsr)+len(Xpremin))
sigmasq_pre = (mu_pre*(1-mu_pre))/(npre)
mu_pre_k = (sum(y_premin)+sum(y))/npre
#print(mu_pre)
#print(mu_pre_k)
#print(mu_pre_PCA)
#print(sigmasq_pre)

y_int = knn.predict(pca.transform(Xint))
int_prob= knn.predict_proba(pca.transform(Xint))
mu_int = (sum(int_prob[:,1])/nint)
C = int_prob[:,1]>thres_prob
y_intpr = C*1
mu_int_PCA = sum(y_intpr)/(len(y_intpr))
sigmasq_int_PCA = mu_int_PCA*(1-mu_int_PCA)/(len(Xint))
mu_int_k = (sum(y_int))/nint
#print(mu_int)
#rint(mu_int_k)
#print(mu_int_PCA)
#print(sigmasq_int_PCA)

y_post1= knn.predict(pca.transform(Xpost1))
post_prob1= knn.predict_proba(pca.transform(Xpost1))
mu_post1 = sum(post_prob1[:,1])/npost1
D= post_prob1[:,1]>thres_prob
y_postpr1 = D*1
mu_post1_PCA = sum(y_postpr1)/(len(y_postpr1))
sigmasq_post1_PCA = mu_post1_PCA*(1-mu_post1_PCA)/len(Xpost1)
mu_post1_k = sum(y_post1)/npost1
#print(mu_post1)
#print(mu_post1_k)
#print(mu_post1_PCA)
#print(sigmasq_post1_PCA)

y_post2 = knn.predict(pca.transform(Xpost2))
post_prob2= knn.predict_proba(pca.transform(Xpost2))
mu_post2 = sum(post_prob2[:,1])/npost2
E= post_prob2[:,1]>thres_prob
y_postpr2 = E*1
mu_post2_PCA = sum(y_postpr2)/(len(y_postpr2))
sigmasq_post2_PCA = mu_post2_PCA*(1-mu_post2_PCA)/len(Xpost2)
mu_post2_k = sum(y_post2)/npost2
#print(mu_post2)
#print(mu_post2_k)
#print(mu_post2_PCA)
#print(sigmasq_post2_PCA)
#%% 
thres_prob = 0.15
B = pre_prob[:,1]>thres_prob
y_pre = B*1
mu_pre_PCA = (sum(y_pre)+sum(y))/(len(Xpremin) + len(Xsr))

D= post_prob1[:,1]>thres_prob
y_postpr1 = D*1
mu_post1_PCA = sum(y_postpr1)/(len(y_postpr1))

E= post_prob2[:,1]>thres_prob
y_postpr2 = E*1
mu_post2_PCA = sum(y_postpr2)/(len(y_postpr2))

print('phir0=',mu_pre_PCA)
print('phir7wks=', mu_post1_PCA)
print('phi10wks=', mu_post2_PCA)
#%%
phi_est= {'phi_t': [1-mu_pre_PCA, 1-mu_post1_PCA, 1-mu_post2_PCA],
        't': [0, 1176, 1656],
        'ncells': [3157,5262,4900]
        }

dfphi= DataFrame(phi_est, columns= ['phi_t', 't', 'ncells'])

print(dfphi)

dfphi.to_csv("phi_t_est_pyth.csv")





#%%
import os

filename = "phi_t_est_pyth.csv"
path = "/Users/kj22643/Documents/Documents/Grant_dose_optimization/data"
fullpath = os.path.join(path, filename)
dfphi.to_csv("phi_t_est_pyth.csv")

#%% Make data frames for the pre, int, and post-treatment cells in PC space

# Pre-treatment
PCsrdf = pd.DataFrame(PCsr)
PCsrdf['classlabel'] = y
PCsrdf['time'] = 0
PCsrdf.reset_index(drop=True, inplace=True)
PCsrdf['recalled_lin'] = dfsr['recalled_lin']


# Pre-treatment
PCpmdf = pd.DataFrame(PCspremin)
PCpmdf['classlabel'] = y_pre
PCpmdf.reset_index(drop=True, inplace=True)
PCpmdf['time'] = 0
PCpmdf['recalled_lin'] = dfpremin['recalled_lin']


# t= 30 hr
PCintdf = pd.DataFrame(PCsint)
PCintdf['classlabel'] = y_intpr
PCintdf.reset_index(drop=True, inplace=True)
#PCintdf['kclust'] = kclustint
PCintdf['time'] = 30
PCintdf['recalled_lin'] = 'nan'


# t = 1176 hr
PCpost1df = pd.DataFrame(PCspost1)
PCpost1df['classlabel'] = y_postpr1
PCpost1df.reset_index(drop=True, inplace=True)
#PCpostdf['kclust'] = kclustpost
PCpost1df['time'] = 1176
PCpost1df['recalled_lin'] = 'nan'

# t= 1656hr
PCpost2df = pd.DataFrame(PCspost2)
PCpost2df['classlabel'] = y_postpr2
PCpost2df.reset_index(drop=True, inplace=True)
#PCpostdf['kclust'] = kclustpost
PCpost2df['time'] = 1176
PCpost2df['recalled_lin'] = 'nan'

#%% PC1 PC2 and PC3


# First just look at the cells used to build the classifier
PCsrdfs = PCsrdf[PCsrdf['classlabel']==0]
PCsrdfr = PCsrdf[PCsrdf['classlabel']==1]
# Then other breakdowns
PCsrdfls = PCsrdf[PCsrdf['recalled_lin']=='sens']
PCsrdflr = PCsrdf[PCsrdf['recalled_lin']=='res']
PCpmdfs = PCpmdf[PCpmdf['classlabel']==0]
PCpmdfr = PCpmdf[PCpmdf['classlabel']==1]
PCintdfs = PCintdf[PCintdf['classlabel']==0]
PCintdfr = PCintdf[PCintdf['classlabel'] == 1]
PCpost1dfs = PCpost1df[PCpost1df['classlabel']==0]
PCpost1dfr = PCpost1df[PCpost1df['classlabel'] == 1]
PCpost2dfs = PCpost2df[PCpost2df['classlabel']==0]
PCpost2dfr = PCpost2df[PCpost2df['classlabel'] == 1]
#%% WHY DOES RUNNING THIS CHANGE MY PC DATAFRAMES????
# Cells used for classifying
xsrs= np.asarray(PCsrdfs[0])
ysrs=np.asarray(PCsrdfs[1])
zsrs=np.asarray(PCsrdfs[2])
#
xsrr= np.asarray(PCsrdfr[0])
ysrr=np.asarray(PCsrdfr[1])
zsrr=np.asarray(PCsrdfr[2])

# Lineages that we isolated
xsrls= np.asarray(PCsrdfls[0])
ysrls=np.asarray(PCsrdfls[1])
zsrls=np.asarray(PCsrdfls[2])

xsrlr= np.asarray(PCsrdflr[0])
ysrlr=np.asarray(PCsrdflr[1])
zsrlr=np.asarray(PCsrdflr[2])

# Pre-treat samples
xp= np.asarray(PCpmdf[0])
yp=np.asarray(PCpmdf[1])
zp=np.asarray(PCpmdf[2])

#pre-treat sensitive
xps= np.asarray(PCpmdfs[0])
yps=np.asarray(PCpmdfs[1])
zps=np.asarray(PCpmdfs[2])
# pre-treat resistant
xpr= np.asarray(PCpmdfr[0])
ypr=np.asarray(PCpmdfr[1])
zpr=np.asarray(PCpmdfr[2])

# t=30 hr cells
xi= np.asarray(PCintdf[0])
yi=np.asarray(PCintdf[1])
zi=np.asarray(PCintdf[2])
#sens and res labels
xis= np.asarray(PCintdfs[0])
yis=np.asarray(PCintdfs[1])
zis=np.asarray(PCintdfs[2])
xir= np.asarray(PCintdfr[0])
yir=np.asarray(PCintdfr[1])
zir=np.asarray(PCintdfr[2])

# t=1176 hr cells
xpo1= np.asarray(PCpost1df[0])
ypo1=np.asarray(PCpost1df[1])
zpo1=np.asarray(PCpost1df[2])
# sens and res labels
xpos1= np.asarray(PCpost1dfs[0])
ypos1=np.asarray(PCpost1dfs[1])
zpos1=np.asarray(PCpost1dfs[2])
xpor1= np.asarray(PCpost1dfr[0])
ypor1=np.asarray(PCpost1dfr[1])
zpor1=np.asarray(PCpost1dfr[2])

# t=1656 hr cells
xpo2= np.asarray(PCpost2df[0])
ypo2=np.asarray(PCpost2df[1])
zpo2=np.asarray(PCpost2df[2])
# sens and res labels
xpos2= np.asarray(PCpost2dfs[0])
ypos2=np.asarray(PCpost2dfs[1])
zpos2=np.asarray(PCpost2dfs[2])
xpor2= np.asarray(PCpost2dfr[0])
ypor2=np.asarray(PCpost2dfr[1])
zpor2=np.asarray(PCpost2dfr[2])
#%%
fig = plt.figure(figsize=(15,15))
ax=fig.add_subplot(111,projection='3d')

# pre-treat cells for classifying
srs = ax.scatter(xsrs, ysrs, zsrs, c='g', marker='^', alpha = 1, label = 't=0 hr labeled sensitive')
srr = ax.scatter(xsrr, ysrr, zsrr, c='r', marker='^', alpha = 1, label = 't=0 hr labeled resistant')

#pre_cells = ax.scatter(xp, yp, zp, c='b', marker='o', alpha = 0.2, label = 't=0 hr remaining')
# classified pre-treatment cells
#ps = ax.scatter(xps, yps, zps, c='olivedrab', marker='+', alpha = 0.5, label = 't=0 hr est sensitive')
#ps = ax.scatter(xpr, ypr, zpr, c='pink', marker='+', alpha = 0.5, label = 't=0 hr est resistant')

#isolated lineages
#ls_pre = ax.scatter(xsrls, ysrls, zsrls, c='lime', marker='o', alpha = 1, label = 'sensitive lineage AA170')
#lr_pre = ax.scatter(xsrlr, ysrlr, zsrlr, c='fuchsia', marker='o', alpha = 1, label = 'resistant lineage AA161')

#int_cells = ax.scatter(xi, yi, zi, c='grey', marker = 'o', alpha = 0.2, label = 't=30 hr unclassified')
#ints = ax.scatter(xis, yis, zis, c='olivedrab', marker = '+', alpha = 0.5, label = 't=30 hr est sensitive')
#intr = ax.scatter(xir, yir, zir, c='pink', marker = '+', alpha = 0.5, label = 't=30 hr est resistant')

#post_cells = ax.scatter(xpo, ypo, zpo, c='c', marker = 'o', alpha = 0.1, label = 't=1344 hr')
#os = ax.scatter(xpos1, ypos1, zpos1, c='olivedrab', marker = '+', alpha = 0.5, label = 't=1176 hr est sensitive')
por = ax.scatter(xpor1, ypor1, zpor1, c='pink', marker = '+', alpha = 0.5, label = 't=1176 hr est resistant')


#pos2 = ax.scatter(xpos2, ypos2, zpos2, c='olivedrab', marker = '+', alpha = 0.5, label = 't=1656 hr est sensitive')
#por2 = ax.scatter(xpor2, ypor2, zpor2, c='pink', marker = '+', alpha = 0.5, label = 't=1656 hr est resistant')


ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
plt.legend(loc=2, prop={'size': 25})
#plt.title('Pre and post-treatment cells in PC space',fontsize= 20)

#ax.legend([sens_cells, res_cells], ['t=0 hr sens', 't=0 hr res'])

ax.azim = 100
ax.elev = -50

#%% PC1 vs PC2

fig = plt.figure(figsize=(10,10))

srs = plt.scatter(xsrs, ysrs, c='g', marker='^', alpha = 1, label = 't=0 hr labeled sensitive')
srr = plt.scatter(xsrr, ysrr, c='r', marker='^', alpha = 1, label = 't=0 hr labeled resistant')
#pre_cells = plt.scatter(xp, yp, c='b', marker='o', alpha = 0.2, label = 't=0 hr remaining')
#pres = plt.scatter(xps, yps, c='olivedrab', marker='+', alpha = 0.5, label = 't=0 hr est sensitive')
#prer = plt.scatter(xpr, ypr, c='pink', marker='+', alpha = 0.5, label = 't=0 hr est resistant')

#int_cells = plt.scatter(xi, yi, c='grey', marker = 'o', alpha = 0.2, label = 't=30 hr unclassified')
#ints = plt.scatter(xis, yis, c='olivedrab', marker = '+', alpha = 0.5, label = 't=30 hr est sensitive')
#intr = plt.scatter(xir, yir, c='pink', marker = '+', alpha = 0.5, label = 't=30 hr est resistant')
#post_cells = plt.scatter(xpo, ypo, c='c', marker = 'o', alpha = 0.1, label = 't=1344 hr')

#ls_pre = plt.scatter(xsrls, ysrls,  c='lime', marker='o', alpha = 1, label = 'sensitive lineage AA170')
#lr_pre = plt.scatter(xsrlr, ysrlr, c='fuchsia', marker='o', alpha = 1, label = 'resistant lineage AA161')
#pos1 = plt.scatter(xpos1, ypos1,  c='olivedrab', marker = '+', alpha = 0.5, label = 't=1176 hr est sensitive')
#por1 = plt.scatter(xpor1, ypor1,  c='pink', marker = '+', alpha = 0.5, label = 't=1176 hr est resistant')
pos2 = plt.scatter(xpos2, ypos2,  c='olivedrab', marker = '+', alpha = 0.5, label = 't=1656 hr est sensitive')
por2 = plt.scatter(xpor2, ypor2,  c='pink', marker = '+', alpha = 0.5, label = 't=1656 hr est resistant')

plt.xlabel('PC1')
plt.ylabel('PC2')

plt.legend(loc=1, prop={'size': 15})

#%% PC1 vs PC3

fig = plt.figure(figsize=(10,10))

srs = plt.scatter(xsrs, zsrs, c='g', marker='^', alpha = 1, label = 't=0 hr labeled sensitive')
srr = plt.scatter(xsrr, zsrr, c='r', marker='^', alpha = 1, label = 't=0 hr labeled resistant')
pre_cells = plt.scatter(xp, zp, c='b', marker='o', alpha = 0.1, label = 't=0 hr remaining')


plt.xlabel('PC1')
plt.ylabel('PC3')

plt.legend(loc=1, prop={'size': 15})
#%% PC2 vs PC3

fig = plt.figure(figsize=(10,10))

srs = plt.scatter(ysrs, zsrs, c='g', marker='^', alpha = 1, label = 't=0 hr labeled sensitive')
srr = plt.scatter(ysrr, zsrr, c='r', marker='^', alpha = 1, label = 't=0 hr labeled resistant')
pre_cells = plt.scatter(xp, zp, c='b', marker='o', alpha = 0.2, label = 't=0 hr remaining')


plt.xlabel('PC2')
plt.ylabel('PC3')

plt.legend(loc=1, prop={'size': 15})
#%% Make a consensus dataframe

PCalldf = pd.concat([PCsrdf, PCpmdf, PCintdf, PCpostdf], axis=0)
# Use consensus data frame to look at the separation between the classlabels
ydfs = PCalldf[PCalldf['classlabel']==0]
yplots= ydfs[0]
xs = np.random.normal(0, 0.04, size=len(yplots))

plt.figure()
bp = PCalldf.boxplot(column=0, by='classlabel', grid=False)
# Add some random "jitter" to the x-axis
qu = plot(xs, yplots, 'r.', alpha=0.2)
plt.ylabel('PC1')

plt.figure()
bp = PCalldf.boxplot(column=1, by='classlabel', grid=False)

plt.ylabel('PC2')