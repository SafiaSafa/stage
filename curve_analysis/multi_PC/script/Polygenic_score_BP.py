# coding: utf-8

# # Polygenic risk score on BP data

# #### Imports

# In[1]:

import sys
sys.path.append("/home/vcabeli/Documents/scripts/python") # Contains "util" package
import util
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.discrete.discrete_model as smd
import scipy.stats
import sklearn.metrics

from scipy.stats import spearmanr

from tqdm import tqdm


def polygen_score(geno, sorted_snps, threshs, test_idces, pheno, betas) :

    pb_cor = np.zeros(len(threshs))
    auc = np.zeros(len(threshs))
    n_samples = len(test_idces)
    n_snps = max(threshs)+1
    sorted_snps = sorted_snps[:n_snps]
    test_cases = [i for i,j in enumerate(test_idces) if pheno['vals'][j]==2]
    test_controls = [i for i,j in enumerate(test_idces) if pheno['vals'][j]==1]

    betas = betas.values[sorted_snps]

    polygen_score_test = np.multiply(np.nan_to_num(geno[test_idces][:,sorted_snps]),
                                     betas)
    polygen_score_test = polygen_score_test.cumsum(axis=1)

    for i, top in enumerate(tqdm(threshs)):
        pb_cor[i] = scipy.stats.pointbiserialr([pheno['vals'][j] for j in test_idces],
                                               np.array(polygen_score_test[:,top]).flatten())[0]
        TPR = np.zeros(n_samples)
        FPR = np.zeros(n_samples)
        for idx, val in enumerate(polygen_score_test[:,top]):
            TPR[idx] = 1.0*np.count_nonzero(polygen_score_test[test_cases,top]>val) / len(test_cases)
            FPR[idx] = 1.0*np.count_nonzero(polygen_score_test[test_controls,top]>val) / len(test_controls)
        auc[i] = sklearn.metrics.auc(x=FPR, y=TPR,reorder=True)

    #return pb_cor
    return auc


def polygen_score_sign(geno, sorted_snps, threshs, test_idces, pheno, betas) :

    pb_cor = np.zeros(len(threshs))
    auc = np.zeros(len(threshs))
    n_samples = len(test_idces)
    n_snps = max(threshs)+1
    sorted_snps = sorted_snps[:n_snps]
    test_cases = [i for i,j in enumerate(test_idces) if pheno['vals'][j]==2]
    test_controls = [i for i,j in enumerate(test_idces) if pheno['vals'][j]==1]

    betas = betas.values[sorted_snps]

    polygen_score_test = np.multiply(np.sign(np.nan_to_num(geno[test_idces][:,sorted_snps])),
                                     betas)
    polygen_score_test = polygen_score_test.cumsum(axis=1)

    for i, top in enumerate(tqdm(threshs)):
        pb_cor[i] = scipy.stats.pointbiserialr([pheno['vals'][j] for j in test_idces],
                                               np.array(polygen_score_test[:,top]).flatten())[0]
        TPR = np.zeros(n_samples)
        FPR = np.zeros(n_samples)
        for idx, val in enumerate(polygen_score_test[:,top]):
            TPR[idx] = 1.0*np.count_nonzero(polygen_score_test[test_cases,top]>val) / len(test_cases)
            FPR[idx] = 1.0*np.count_nonzero(polygen_score_test[test_controls,top]>val) / len(test_controls)
        auc[i] = sklearn.metrics.auc(x=FPR, y=TPR,reorder=True)

    #return pb_cor
    return auc


