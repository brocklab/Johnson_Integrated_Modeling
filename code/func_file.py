#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 10:35:07 2019

@author: kj22643
"""
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

def find_mean_AUC(folds_dict, n_neighbors, n_components):

    
    kCV = len(folds_dict['trainmat'])

    pca=PCA(copy=True, iterated_power='auto', n_components=n_components, random_state=0,
            svd_solver='auto', tol=0.0, whiten=False)
# Use a nearest neighbor classifier to evaluate the methods
    knn = KNeighborsClassifier(n_neighbors=n_neighbors)


# #Compute the ROC curve for each fold


    AUCPCA=np.zeros((kCV,1))
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
        

        AUCPCA[i]= roc_auc
  


    meanAUC_PCA = AUCPCA.mean()
    return(meanAUC_PCA)


def find_mean_AUC_SVM(folds_dict, C, basis):

    
    kCV = len(folds_dict['trainmat'])
    clf = svm.SVC(kernel=basis, C=C, probability = True)



# #Compute the ROC curve for each fold


    AUCSVM=np.zeros((kCV,1))
    acc = np.zeros((kCV,1))
    for i in range(kCV):
        X_train = folds_dict['trainmat'][i]
        y_train = folds_dict['trainlabel'][i]
        X_test = folds_dict['testmat'][i]
        y_test = folds_dict['testlabel'][i]
    
    # Build new model based on training data set
        clf.fit(X_train, y_train)
        SVM_prob_scores= clf.predict_proba(X_test)

        fpr, tpr, thresholds = metrics.roc_curve(y_test, SVM_prob_scores[:,1], pos_label=1)
        roc_auc = auc(fpr, tpr)
        

        AUCSVM[i]= roc_auc
  


    meanAUC_SVM = AUCSVM.mean()
    return(meanAUC_SVM)

















# This is the equivalent of your find_eigendigits function in Matlab
def find_meanvec(A):

# It starts by finding the mean vector
    m = A.mean(axis =1)
    return m

def find_eigvals(A,mu):
    X = A.sub(mu, axis=0) # we subtract from each column
    
    # convert the subtracted gene expression matrix to a matrix 
    # Xmat should be ngenes x ntrain
    Xmat= X.as_matrix()
    # transpose it 
    # this should be ntrain x ngenes
    XmatT = Xmat.T
    # since we have less cells than genes, we can use the trick to make a smaller cov matrix
    # which is the (ntrain x ngenes * ngenes x ntrain = ntrain x ntrain square small covariance matrix)
    smallcov = np.matmul(XmatT,Xmat) # Now you have your ntrain x ntrain small covariance matrix

    # should get a vector of ntrain lambdas, and a square matrix ntrain x ntrain of the eigenvectors
    lambdas, Vsmall= LA.eig(smallcov)
    ind_arr = np.argsort(-abs(lambdas))

    # Now your ind_arr should start with the highest eigenvalue and work its way down
    # print the indices to check that the first is the highest lambda and the last is the lowest
 
    ordered_lambdas = lambdas[ind_arr]
    return ordered_lambdas

def find_eigvecs(A, mu):
   
    X = A.sub(mu, axis=0) 
    Xmat= X.as_matrix()
    XmatT = Xmat.T
    smallcov = np.matmul(XmatT,Xmat) 
    lambdas, Vsmall= LA.eig(smallcov)
    ind_arr = np.argsort(-abs(lambdas))
    # Reorder your lambdas and eigenvectors (Vsmall)
    Vnew = Vsmall[ind_arr]
    # Now apply to the big system
    # Since XX'x=mux
    # Let x = Xv
    # XX'Xv = mu Xv
    # XX' = big cov matrix
    # X'X = small cov matrix

    # These are the eigenvectors of the big covariance matrix 
    # ngenes x ntrain x ntrain x ntrain gives ngenes x ntrain
    Vbig = np.matmul(Xmat, Vnew)
    # These vectors of length ntrain are the eigenvectors, in order of importance (eigenvalue)
    # Vbig is now an ngenes by ntrain matrix of eigenvectors (each column =1 eigenvector)

    # Renormalize: Now that you have your ngenes x ntraining cells matrix of 
    #eigenvectors, make sure the big covariance matrix is normalized
    norms = LA.norm(Vbig, axis = 0)
    Vnorm = Vbig/norms # divides each column by i
    return Vnorm

def project_cells(A, V, neigs):
    # Now you have the Eigenspace, defined by Vnorm, which contains normalized eigenvectors of
    # the length of the number of genes, in order of importance.
    
    # Now we need to do the projection of the individual cells (in A) onto the eigenvector space

    #Mcol = reshape(M, [x,1]);
    #b = (Mcol-m)'*V; %(1 by x) *(x*k) Gives a vector of length k

    # Project onto eigenspace
    #recon = V*b'; % (x by k) * (k *1) Gives x by 1

    # Declare the number of eigenvectors we want to use from the Vnorm matrix 
    #eventually we will loop through the number of eigenvectors and assess the accuracy as afunction of neigs
    ntrain = A.shape[1]
    trainmat = np.zeros((neigs, ntrain))
    trainvec= np.zeros((neigs,1))

    # First make your Omat which contains the coordinates (the columns of each training image
    # and has the corresponding assignment (sens or res) ))
    m = A.mean(axis =1)
    X = A.sub(m, axis=0) # we subtract from each column
    Xmat= X.as_matrix()
    for i in range(ntrain):
        Mcol = Xmat[:,i] # start right after the last training column in Amat 
        McolT= Mcol.T # to get the first testing cell as a row
    # 1xngene x ngene x ntrain = 1 x ntrain
        b=np.matmul(McolT,)

        trainvec= b[:neigs]
        trainmat[:,i]= trainvec
    return trainmat