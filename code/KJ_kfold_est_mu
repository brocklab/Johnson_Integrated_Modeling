#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:22:02 2019

@author: kj22643
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 15:11:10 2019

@author: kj22643
"""

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
adata = sc.read('daylin_anndata.h5ad')
adata.obs.head()
# current samples:
#BgL1K
#30hr
#Rel-1
#Rel-2
# We will change these to time points
#%% Assign survivor category in adata.obs
longTreatLins = adata.obs.loc[(adata.obs['sample'].isin(['Rel-1','Rel-2']))&(adata.obs.lineage!='nan'),'lineage'].unique().tolist()

adata.obs.loc[adata.obs.lineage.isin(longTreatLins)==False,'survivor'] = 'sens'
adata.obs.loc[adata.obs.lineage.isin(longTreatLins)==True,'survivor'] = 'res'

# %%try to rename the samples by time point
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
adata_pre = adata[adata.obs['timepoint']=='t=0hr', :]
dfpre = pd.concat([adata_pre.obs['survivor'], adata_pre.obs['lineage'],
               pd.DataFrame(adata_pre.raw.X,index=adata_pre.obs.index,
                            columns=adata_pre.var_names),],axis=1) 
# t = 30 hr (intermediate timepoint) 5169 int treatment cells
adata_int = adata[adata.obs['timepoint']=='t=30hr', :]
dfint = pd.concat([adata_int.obs['lineage'],
                   pd.DataFrame(adata_int.raw.X, index=adata_int.obs.index, 
                                columns = adata_int.var_names),], axis=1)

# t=1344 hr (~roughly 8 weeks), 10332 post treatment cells
adata_post = adata[adata.obs['timepoint']=='t=1344hr', :]
dfpost = pd.concat([adata_post.obs['lineage'],
                    pd.DataFrame(adata_post.raw.X, index=adata_post.obs.index, 
                                 columns = adata_post.var_names),],axis =1)

#%% Use sklearn to do principle component analysis on the  entire pre-treatment sample
#X = dfpre.loc[:, dfpre.columns !='survivor', dfpre.columns !='lineage']
X = dfpre.drop(columns= ['survivor', 'lineage'])

y= pd.factorize(dfpre['survivor'])[0] 
ncells = len(y)

mu_pre = sum(y)/len(y)
# X is your cell gene matrix, y is your class labels
#%% Set up cross validation where your test set is not contained in your training set at all
# Split into train/test
kCV = 5
# use this function to ensure that the class balance is maintained for each of your test sets
skf = StratifiedKFold(n_splits=kCV, shuffle= True)
Atrain = {}
Atest = {}
ytest = {}
ytrain = {}
mu_true_test = {}
ntest = {}

folds_dict = {'trainmat':{}, 'trainlabel':{}, 'V':{}, 'lambdas':{}, }
for i in range(kCV):    
    for train_index, test_index in skf.split(X, y):
        Atrain[i] = X.iloc[train_index, :]
        Atest[i] = X.iloc[test_index, :]
        ytest[i]= y[test_index]
        ytrain[i]= y[train_index]
        mu_true_test[i] = sum(ytest[i])/len(ytest[i])
        ntest[i]=len(ytest[i])
         
# Save all of your stratified folds into a single dictionary. 
folds_dict['trainmat'] = Atrain
folds_dict['trainlabel']= ytrain
folds_dict['testmat'] = Atest
folds_dict['testlabel'] = ytest
folds_dict['prevtest'] = mu_true_test
folds_dict['ntest'] = ntest




n_classes = len(np.unique(y))
      


# %%Assign the optimal parameters (foudn from KJ_classify_sklearn.py) for building your prediction model 
# p(x|Sj) where Sj is your training set and x is any new cell (in your test set or in future cells)
n_neighbors = 15
n_components = 100
random_state = 0

Copt = 1000
basis = 'rbf'


knn = KNeighborsClassifier(n_neighbors=n_neighbors)

pca=PCA(copy=True, iterated_power='auto', n_components=n_components, random_state=0,
            svd_solver='auto', tol=0.0, whiten=False)

clf = svm.SVC(kernel=basis, C=Copt)

#%% Build a new model for each fold and apply it to the test set to generate mu_hat
# where mu_hat_j is the average value of the test set predictions (sens=0, res =1) from the training set model
y_PCA = {}
mu_PCA = {}
y_SVM = {}
mu_SVM = {}
V_train = {}


for i in range(kCV):
    X_train = folds_dict['trainmat'][i]
    y_train = folds_dict['trainlabel'][i]
    X_test = folds_dict['testmat'][i]
    y_test = folds_dict['testlabel'][i]
    
    # PCA MODEL OUTPUTS FOR EACH FOLD
    pca.fit(X_train, y_train)
    V_train[i] = pca.fit_transform(X_train)
# Fit a nearest neighbor classifier on the model built on the training data set
    knn.fit(pca.transform(X_train), y_train)
    y_PCA[i] = knn.predict(pca.transform(X_test))
# Compute the nearest neighbor accuracy on the embedded test set
    mu_PCA[i] = sum(y_PCA[i])/ folds_dict['ntest'][i]
    
    # SVM MODEL OUTPUTS FOR EACH FOLD
    clf.fit(X_train, y_train)
    y_SVM[i]= clf.predict(X_test)
    mu_SVM[i] = sum(y_SVM[i])/folds_dict['ntest'][i]
    
#%%  Put into folds_dict
folds_dict['V_train'] = V_train
folds_dict['y_PCA']= y_PCA
folds_dict['mu_PCA'] = mu_PCA
folds_dict['y_SVM'] = y_SVM
folds_dict['mu_SVM'] = mu_SVM
   
    
#%%Compare the mu_SVM and mu_PCA test set estimates of the expectation to the known
    # prevalence of the test set
df = pd.DataFrame()
dfprevtest = pd.DataFrame(folds_dict['prevtest'], index=[0])
dfmu_PCA= pd.DataFrame(folds_dict['mu_PCA'], index = [0])
dfmu_SVM = pd.DataFrame(folds_dict['mu_SVM'], index = [0])

npprevtest=np.array(dfprevtest)
npmu_PCA = np.array(dfmu_PCA)
npmu_SVM = np.array(dfmu_SVM)
mu_pre = np.mean(npprevtest)
mu_hat_PCA = np.mean(npmu_PCA)
mu_hat_SVM = np.mean(npmu_SVM)
ntest = folds_dict['ntest'][0]
print(mu_hat_PCA)
print(mu_hat_SVM)
                 
sigmasq_PCA =(1-mu_hat_PCA)*mu_hat_PCA/ntest
sigmasq_SVM = (1-mu_hat_SVM)*mu_hat_SVM/ntest
print(sigmasq_PCA)
print(sigmasq_SVM)


#%% Next step, apply the models to the subsequent time points!

    
    
    
    
    
    
    
    
    
    
    