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


# #### Load bed bim fam, BP phenotype and covariates

# In[2]:

BP_DIR = "/home/vcabeli/Documents/data/BP_final/"
snp_data, pheno = util.load_data(BP_DIR + "/BP.B37-final")

pheno_sj = pd.read_table("/home/vcabeli/Documents/data/latest_imputed_genotyped_bd/pheno_cleaned.txt",
                         delim_whitespace=True)
pheno['bp_type'] = np.zeros(snp_data.row_count)


# In[3]:

BP1_samples = pheno_sj.IID[pheno_sj.BP1 == 1].values
pheno['bp_type'][[sample_id[1] in BP1_samples for sample_id in pheno['iid']]] = 1

BP2_samples = pheno_sj.IID[pheno_sj.BP2 == 1].values
pheno['bp_type'][[sample_id[1] in BP2_samples for sample_id in pheno['iid']]] = 2


# In[4]:

cov = pd.read_table("/home/vcabeli/Documents/data/BP_final/BP.B37-final.cov",
                    delim_whitespace=True)


# #### Load logistic regression results on complete dataset

# In[5]:

logistic_res = pd.read_table("/home/vcabeli/Documents/analyses_current/assoc/genotyped/class_logistic..assoc.logistic",
                             delim_whitespace=True)
logistic_res = logistic_res[logistic_res.TEST == "ADD"]
logistic_res.reset_index(drop=True, inplace=True)


# In[6]:

assert(np.all(logistic_res.SNP == snp_data.col)) # Same indices correspond to the same SNPs?


# #### Divide data in training/test sets


# In[26]:

i_scores_complete = util.get_i_scores_binary()
G_hw = snp_data.val.copy()


# In[27]:

MAFs = np.nansum(2-G_hw, axis=0, ) / (np.count_nonzero(~np.isnan(G_hw), axis=0) * 2)
G_hw = (2-G_hw - 2*MAFs)/np.sqrt(2*MAFs*(1-MAFs))


# In[7]:
n_iter = 100


# In[18]:

threshs = range(10, 250140, 500)
threshs_clumped = range(10, 80000,300)


res_hw_pval_sign = np.zeros([n_iter, len(threshs)])
res_hw_iscore_sign = np.zeros([n_iter, len(threshs)])

res_hw_pval_sign_clumped = np.zeros([n_iter, len(threshs_clumped)])
res_hw_iscore_sign_clumped = np.zeros([n_iter, len(threshs_clumped)])


for iteration in range(n_iter):

    print "=========================Iteration {}=========================".format(iteration)

    train_idces = np.random.choice(np.arange(snp_data.row_count), size=int(snp_data.row_count*0.5), replace=False)
    test_idces = np.setdiff1d(np.arange(snp_data.row_count), train_idces, assume_unique=True)

    stats_sets = pd.DataFrame(index= ['Training', 'Test'],
                              columns=['n_controls', 'n_cases', 'n_bp1', 'n_bp2'])

    for i, idces in enumerate([train_idces, test_idces]):
        stats_sets.iloc[i]['n_controls'] = np.count_nonzero([pheno['vals'][j] == 1 for j in idces])
        stats_sets.iloc[i]['n_cases'] = np.count_nonzero([pheno['vals'][j] == 2 for j in idces])
        stats_sets.iloc[i]['n_bp1'] = np.count_nonzero(pheno['bp_type'][idces] == 1)
        stats_sets.iloc[i]['n_bp2'] = np.count_nonzero(pheno['bp_type'][idces] == 2)


    # In[8]:

    training_sample_out = "training_samples.keep"

    with open(training_sample_out, 'w') as f:
        for i in train_idces:
            f.write(pheno['iid'][i][0] + "\t" + pheno['iid'][i][1])
            f.write("\n")


    # ## Classical polygenic score : pvalue threshold, logistic $\beta$s
    #
    # Where SNP selection is based on their logitic association test p-value and weights are the $\beta$s.

    # In[9]:

    os.chdir("/home/vcabeli/Documents/analyses_current/polygenic_score/BP")

    cmd_plink = "plink --bfile {} --keep {} --logistic sex beta --covar {} --hide-covar --out {}".format("/home/vcabeli/Documents/data/BP_final/BP.B37-final",
                                                                                                         training_sample_out,
                                                                                                         "/home/vcabeli/Documents/data/BP_final/BP.B37-final.cov",
                                                                                                         "train_log")
    print cmd_plink
    p = subprocess.Popen(cmd_plink, shell=True)
    p.wait()

    train_res = pd.read_table("train_log.assoc.logistic", delim_whitespace=True)
    assert(np.all(logistic_res.SNP == train_res.SNP))



    util.control_idces = [i for i in train_idces if pheno['vals'][i] == 1]
    util.case_idces = [i for i in train_idces if pheno['vals'][i] == 2]

    i_scores_train = util.get_i_scores_binary()


    # In[15]:

    sorted_snps_pval = np.argsort(train_res.P.values)
    sorted_snps_iscore = np.argsort(-i_scores_train)




    # ## Coded genotyped
    #
    # Under Hardy-Weinberg equilibrium, $G_i = \frac{(X_i - 2 f_i)}{\sqrt{2fi(1-fi)}}$ where $X_i$ is the number of minor alleles and $f_i$ the minor allele frequency.

    # In[31]:

    ttest_test = polygen_score_sign(G_hw, sorted_snps_pval,
                                    threshs,
                                    test_idces, pheno,
                                    train_res.BETA)


    # In[32]:

    res_hw_pval_sign[iteration] = ttest_test



    # ## I score threshold

    # In[39]:

    ttest_test = polygen_score_sign(G_hw, sorted_snps_iscore,
                                    threshs,
                                    test_idces, pheno,
                                    train_res.BETA)


    # In[40]:

    res_hw_iscore_sign[iteration] = ttest_test

    # Clump

    os.chdir("/home/vcabeli/Documents/analyses_current/polygenic_score/BP")

    cmd_plink_clump = "plink --bfile {} --keep {} --clump {} --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 500 --out {} ".format("/home/vcabeli/Documents/data/BP_final/BP.B37-final",
                                                                                                         training_sample_out,
                                                                                                         "train_log.assoc.logistic",
                                                                                                         "train_clump_0.2_500")
    print cmd_plink_clump
    p = subprocess.Popen(cmd_plink_clump, shell=True)
    p.wait()

    clumped_res = pd.read_table("train_clump_0.2_500.clumped",
                            delim_whitespace=True)
    clumped_snps = set(clumped_res.SNP)

    clumped_sorted_snps_pval = np.array([i for i in sorted_snps_pval if snp_data.col[i] in clumped_snps])
    clumped_sorted_snps_iscore = np.array([i for i in sorted_snps_iscore if snp_data.col[i] in clumped_snps])

    ttest_test = polygen_score_sign(G_hw, clumped_sorted_snps_pval,
                                threshs_clumped,
                                test_idces, pheno,
                                train_res.BETA)

    res_hw_pval_sign_clumped[iteration] = ttest_test

    ttest_test = polygen_score(G_hw, clumped_sorted_snps_iscore,
                           threshs_clumped,
                           test_idces, pheno,
                           train_res.BETA)

    res_hw_iscore_sign_clumped[iteration] = ttest_test



np.savetxt(X=res_hw_pval_sign, fname="res_hw_pval_sign.txt")
np.savetxt(X=res_hw_iscore_sign, fname="res_hw_iscore_sign.txt")

np.savetxt(X=res_hw_pval_sign_clumped, fname="res_hw_pval_sign_clumped.txt")
np.savetxt(X=res_hw_iscore_sign_clumped, fname="res_hw_iscore_sign_clumped.txt")
