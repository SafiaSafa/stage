import numpy as np
import sklearn.metrics
from tqdm import tqdm


def polygen_score(geno, sorted_snps, threshs, test_idces, pheno, betas) :

    auc = np.zeros(len(threshs))
    n_samples = len(test_idces)
    n_snps = max(threshs)+1
    sorted_snps = sorted_snps[:n_snps]
    test_cases = [i for i,j in enumerate(test_idces) if pheno['vals'][j]==2]
    test_controls = [i for i,j in enumerate(test_idces) if pheno['vals'][j]==1]

    betas = betas.values[sorted_snps]

    polygen_score_test = np.multiply(np.nan_to_num(geno[test_idces][:,sorted_snps]),
                                     np.nan_to_num(betas))
    polygen_score_test = polygen_score_test.cumsum(axis=1)

    for i, top in enumerate(tqdm(threshs)):
        TPR = np.zeros(n_samples)
        FPR = np.zeros(n_samples)
        for idx, val in enumerate(polygen_score_test[:,top]):
            TPR[idx] = 1.0*np.count_nonzero(polygen_score_test[test_cases,top]>val) / len(test_cases)
            FPR[idx] = 1.0*np.count_nonzero(polygen_score_test[test_controls,top]>val) / len(test_controls)
        auc[i] = sklearn.metrics.auc(x=FPR, y=TPR,reorder=True)

    return auc


def polygen_score_sign(geno, sorted_snps, threshs, test_idces, pheno, betas) :

    auc = np.zeros(len(threshs))
    n_samples = len(test_idces)
    n_snps = max(threshs)+1
    sorted_snps = sorted_snps[:n_snps]
    test_cases = [i for i,j in enumerate(test_idces) if pheno['vals'][j]==2]
    #print(len(test_cases))
    test_controls = [i for i,j in enumerate(test_idces) if pheno['vals'][j]==1]
    #print(len(test_controls))
    betas = betas.values[sorted_snps]
    # geno = vecteur de taille m (nb snp) de marqueurs du code genetique [du groupe test] [trie par p-value]
    # np.nan_to_num = Replace nan with zero and inf with large finite numbers.
    polygen_score_test = np.multiply(np.nan_to_num(geno[test_idces][:,sorted_snps]),
                                     np.sign(np.nan_to_num(betas)))
    polygen_score_test = polygen_score_test.cumsum(axis=1)
    #print "val",polygen_score_test[:,140][1]

    # score polygenique calculer sur 70000 snps mais auc tout les 140(=len(threshs)) snps (lisse la courbe)
    # top score polygenique tout les 140 snps
    for i, top in enumerate(tqdm(threshs)):
        TPR = np.zeros(n_samples)
        FPR = np.zeros(n_samples)
        for idx, val in enumerate(polygen_score_test[:,top]):
            #print polygen_score_test[test_cases,top]>val
            #print np.count_nonzero(polygen_score_test[test_cases,top]>val) / len(test_cases)
            TPR[idx] = 1.0*np.count_nonzero(polygen_score_test[test_cases,top]>val) / len(test_cases) ##############PQ??????????
            FPR[idx] = 1.0*np.count_nonzero(polygen_score_test[test_controls,top]>val) / len(test_controls)
        auc[i] = sklearn.metrics.auc(x=FPR, y=TPR,reorder=True)

    return auc