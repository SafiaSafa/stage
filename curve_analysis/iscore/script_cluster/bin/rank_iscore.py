#!/usr/bin/python2
# coding: utf-8

# In[1]:

import sys
sys.path.append("/media/vcabeliNFS2/safia")
import util
from Polygenic_score import *
import argparse
import os
import subprocess
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.metrics
import re

# In[2]:


def clump_sorted_snps(bfile, assoc_file, snp_values, sorted_snps):

    clump_out = assoc_file[:assoc_file.index(".assoc")] + "_clump"
    cmd_plink_clump = "/media/vcabeliNFS2/safia/bin/plink --bfile {} --clump {} --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 500 --threads 8 --out {} ".format(bfile,
            assoc_file,
            clump_out)
    print cmd_plink_clump
    p = subprocess.Popen(cmd_plink_clump, shell=True)
    p.wait()

    clumped_res = pd.read_table(clump_out+".clumped", delim_whitespace=True)
    clumped_snps = set(clumped_res.SNP)

    clumped_sorted_snps = np.array([i for i in sorted_snps if snp_values[i] in clumped_snps])
    return clumped_sorted_snps



# In[3]:

def save_strat(strat_file, strat):
    pca_res = pd.read_table(strat_file, delim_whitespace=True, skiprows=1, header=None)
    PCs_out = strat + ".PCs"
    with open(PCs_out, "a") as f:
        f.write(', '.join(pca_res[0].values))
        f.write('\n')
        f.write(str(pca_res[2].values.tolist())[1:-1])
        f.write('\n')
        f.write(str(pca_res[3].values.tolist())[1:-1])
        f.write('\n')


# In[4]:

def plink_prune(plink_bfile):

    cmd_prune = """/media/vcabeliNFS2/safia/bin/plink --bfile {} --exclude /media/vcabeliNFS2/data/high-LD-regions_37.txt --range --indep-pairwise 50 5 0.2 --allow-extra-chr --threads 8 --out {}""".format(plink_bfile, plink_bfile+"_prune")
    p = subprocess.Popen(cmd_prune, shell=True)
    assert(p.wait() == 0)

    cmd_extract = """/media/vcabeliNFS2/safia/bin/plink --bfile {} --extract {} --allow-extra-chr --chr 1-23 --make-bed --threads 8 --out {}""".format(plink_bfile, plink_bfile + "_prune.prune.in", 
                             plink_bfile + "_pruned")
    p = subprocess.Popen(cmd_extract, shell=True)
    assert(p.wait() == 0)

    print "Wrote {} bed/bim/fam.".format(plink_bfile + "_pruned")



# In[5]:

def write_cov_file(pca_file,nb_pc):

    cmd_sed = "sed -i -e \"s/:/ /g\" -e \"s/\s\+/\t/g\" -e \"s/^\s\+//g\" {}".format(pca_file)
    p = subprocess.Popen(cmd_sed, shell=True)
    assert(p.wait()==0)
    
    cov_file = pca_file[:-9]+".cov"
    a=""
    for i in range(nb_pc):
        a=a+("\tPC{}".format(i+1))
    

    cmd_print = "printf \"FID\tIID{}\n\" > {}".format(a,cov_file)

    p = subprocess.Popen(cmd_print, shell=True)
    assert(p.wait()==0)

    cmd_tail = "tail -n+2 {} | cut -f1-{} >> {}".format(pca_file, 2+nb_pc, cov_file) 

    #-n+2 récup lignes en commançant par le 2ème
    #-f1-2+nb_pc résup les champs 1-nb_pc+2
    
    p = subprocess.Popen(cmd_tail, shell=True)
    assert(p.wait()==0)
    
    print "Wrote {}".format(cov_file)


# In[6]:

def iterate(snp_data, pheno,n_iter,threshs,thresh_LD,nb_pc,repertoire):
    res_hw_pval_sign = np.zeros([n_iter, len(threshs)])
    
    for iteration in range(n_iter):

        print "=========================Iteration {}=========================".format(iteration+1)
        
        # extraction des SNPs indépendantes selon le seuil de LD
        cmd_extract_low_ld_snps = "/media/vcabeliNFS2/safia/bin/plink --bfile {} --extract {} --make-bed --threads 8 --out {}".format("/media/vcabeliNFS2/data/BP/BP.B37-final",
                                                                                             "/media/vcabeliNFS2/data/L2_thresh_{}_BP.extract".format(thresh_LD),
                                                                                             repertoire+"/low_ld")
        print cmd_extract_low_ld_snps
        p = subprocess.Popen(cmd_extract_low_ld_snps, shell=True)
        assert(p.wait()==0)

        plink_prune(repertoire+"/low_ld")


 
        # Covariables
        if(nb_pc!=0):
            cmd_smart_pca2 = "/media/vcabeliNFS2/safia/bin/smartpca.perl -i {} -a {} -b {} -o {} -p {} -e {} -l {} -k {} -t {} -m {}".format(repertoire+"/low_ld_pruned.bed",
                                                                                                           repertoire+"/low_ld_pruned.bim",
                                                                                                           repertoire+"/low_ld_pruned.fam",
                                                                                                           repertoire+"/low_ld_pruned.pca",
                                                                                                           repertoire+"/low_ld_pruned.plot",
                                                                                                           repertoire+"/low_ld_pruned.eval",
                                                                                                           repertoire+"/low_ld_pruned.log",
                                                                                                           nb_pc, 2, 0)
            print cmd_smart_pca2
            p = subprocess.Popen(cmd_smart_pca2, shell=True)
            p.wait()


            write_cov_file(repertoire+"/low_ld_pruned.pca.evec",nb_pc)
            save_strat(repertoire+"/low_ld_pruned.pca.evec", repertoire+"/low_ld{}".format(thresh_LD)) #OUTPUT LINE

        
        # #### Split training / testing datasets
        train_idces = np.random.choice(np.arange(snp_data.row_count), size=int(snp_data.row_count*0.5), replace=False)
        test_idces = np.setdiff1d(np.arange(snp_data.row_count), train_idces, assume_unique=True)




        training_sample_out = repertoire+"/training_samples.keep"

        with open(training_sample_out, 'w') as f:
            for i in train_idces:
                f.write(pheno['iid'][i][0] + "\t" + pheno['iid'][i][1])
                f.write("\n")




        # #### Build training set bed bim fam


        #creation fichier plink (bam,bed,fam) avec les données d'apprentissage seulement
        cmd_keep_plink = "/media/vcabeliNFS2/safia/bin/plink --bfile {} --keep {} --make-bed --threads 8 --out {}".format("/media/vcabeliNFS2/data/BP/BP.B37-final",
                                                                                 training_sample_out,
                                                                                 repertoire+"/training_set")
        print cmd_keep_plink
        p = subprocess.Popen(cmd_keep_plink, shell=True)
        p.wait()





        #test d'association: classification (Rlog)
        if (nb_pc==0):
            cmd_second_gwas = "/media/vcabeliNFS2/safia/bin/plink --bfile {} --logistic sex beta hide-covar --threads 8 --out {}".format(repertoire+"/training_set",repertoire+"/low_ld")
        elif (nb_pc==1):
            cmd_second_gwas = "/media/vcabeliNFS2/safia/bin/plink --bfile {} --logistic sex beta hide-covar --covar {} --covar-name PC1 --threads 8 --out {}".format(repertoire+"/training_set",
                                                                                                           repertoire+"/low_ld_pruned.cov",
                                                                                                           repertoire+"/low_ld")
        else:
            cmd_second_gwas = "/media/vcabeliNFS2/safia/bin/plink --bfile {} --logistic sex beta hide-covar --covar {} --covar-name PC1-PC{} --threads 8 --out {}".format(repertoire+"/training_set",
                                                                                                           repertoire+"/low_ld_pruned.cov",
                                                                                                            nb_pc,
                                                                                                           repertoire+"/low_ld")
        print cmd_second_gwas
        p = subprocess.Popen(cmd_second_gwas, shell=True)
        p.wait()
        low_ld_res = pd.read_table(repertoire+"/low_ld.assoc.logistic", delim_whitespace=True)
        
        # Compute I score on the training set 
        snp_data_iscore, pheno_iscore = util.load_data("training_set")
        i_scores_classification = util.get_i_scores_binary()
        low_ld_res['I_score']=i_scores_classification

        # Classe les SNPs par I scores 
        sorted_snps_low_ld = np.argsort(abs(low_ld_res.I_score))
        # clump les SNPs pour garder les SNPs indépendantes 
        sorted_snps_low_ld = clump_sorted_snps(repertoire+"/training_set", repertoire+"/low_ld.assoc.logistic",
                                               snp_data.col, sorted_snps_low_ld)

        
        # #### Coded genotype (hardy-weinberg)
        # 
        # With sufficiently large cohorts, traning and test sets should have the same MAFs. 
        # The coded genotype is computed on the whole dataset.




        G_hw = snp_data.val.copy()

        MAFs = np.nansum(2-G_hw, axis=0, ) / (np.count_nonzero(~np.isnan(G_hw), axis=0) * 2)
        G_hw = (2-G_hw - 2*MAFs)/np.sqrt(2*MAFs*(1-MAFs))

        # ### Compute PRS
        
        prs_low_ld = polygen_score_sign(G_hw, sorted_snps_low_ld,
                                        threshs,
                                        test_idces, pheno,
                                        low_ld_res.BETA)
        res_hw_pval_sign[iteration] = prs_low_ld
        
        ###### Compute I-score
        #i_scores_classification = util.get_i_scores_binary()

        #df_iscore = pd.DataFrame({'SNP':snp_data.col, 'i_score':i_scores_classification})
        #data frame a 2 colone (pê recupe que la deuxieme)
        #res_hw_iscore[iteration] = df_iscore
        
        ####### Compute LD-score 
        
    return(res_hw_pval_sign)


# In[ ]:

if __name__ == '__main__':

    snp_data, pheno = util.load_data("/media/vcabeliNFS2/data/BP/BP.B37-final")
    n_iter = 100
    threshs = range(10, 95000, 500)
    thresh_LD = 5
    nb_pc_tot=51
    time = range(10,95000,500)
    get_color = lambda : "#" + "".join(np.random.choice(list("02468acef"), size=6))
    step = 10
    ############
    repertoire="/media/vcabeliNFS2/safia/"+str(sys.argv[1])
    os.system("mkdir "+repertoire)
    res_pval=iterate(snp_data, pheno,n_iter,threshs,thresh_LD,int(sys.argv[1]),repertoire)
    os.system("mkdir "+repertoire+"/resultat")
    with open(repertoire+"/resultat/prs_low_ld_{}.res".format(thresh_LD), 'a') as f: #OUTPUT LINE
        df = pd.DataFrame(res_pval)
        df.to_csv(f, header=None,index=None)
