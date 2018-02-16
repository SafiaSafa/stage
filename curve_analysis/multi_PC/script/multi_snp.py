#!/usr/bin/python2
# coding: utf-8

# In[2]:

import sys
import util
from Polygenic_score import *

import os
import subprocess
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.metrics
import re
import sklearn.metrics
from tqdm import tqdm


# In[3]:

def clump_sorted_snps(bfile, assoc_file, snp_values, sorted_snps):

    clump_out = assoc_file[:assoc_file.index(".assoc")] + "_clump"
    cmd_plink_clump = "plink --bfile {} --clump {} --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 500 --out {} ".format(bfile,
            assoc_file,
            clump_out)
    print cmd_plink_clump
    p = subprocess.Popen(cmd_plink_clump, shell=True)
    p.wait()

    clumped_res = pd.read_table(clump_out+".clumped", delim_whitespace=True)
    clumped_snps = set(clumped_res.SNP)

    clumped_sorted_snps = np.array([i for i in sorted_snps if snp_values[i] in clumped_snps])
    return clumped_sorted_snps



# In[4]:

def save_strat(strat_file, strat):
    pca_res = pd.read_table(strat_file, delim_whitespace=True, skiprows=1, header=None)
    PCs_out = strat + ".PCs"
    #with open("res_clumped/"+PCs_out, "a") as f:
    with open(PCs_out, "a") as f:
        f.write(', '.join(str(pca_res[0].values)))
        f.write('\n')
        f.write(str(pca_res[2].values.tolist())[1:-1])
        f.write('\n')
        f.write(str(pca_res[3].values.tolist())[1:-1])
        f.write('\n')


# In[5]:

def plink_prune(plink_bfile):

    cmd_prune = """plink --bfile {}      --exclude /home/vcabeli/Documents/data/high-LD-regions_37.txt           --range --indep-pairwise 50 5 0.2           --allow-extra-chr           --out {}""".format(plink_bfile, plink_bfile+"_prune")
    p = subprocess.Popen(cmd_prune, shell=True)
    assert(p.wait() == 0)

    cmd_extract = """plink --bfile {}           --extract {}           --allow-extra-chr           --chr 1-23           --make-bed           --out {}""".format(plink_bfile, plink_bfile + "_prune.prune.in", 
                             plink_bfile + "_pruned")
    p = subprocess.Popen(cmd_extract, shell=True)
    assert(p.wait() == 0)

    print "Wrote {} bed/bim/fam.".format(plink_bfile + "_pruned")



# In[6]:

def write_cov_file(pca_file,nb_pc=2):

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


# In[18]:

def iterate_bis(snp_data, pheno,n_iter,threshs,thresh_LD,nb_pc,nb_snps):
    res_hw_pval_sign = np.zeros([n_iter, len(threshs)])
    for iteration in range(n_iter):

        print "=========================Iteration {}=========================".format(iteration)
        
        # extraction des SNPs indépendantes selon le seuil de LD
        cmd_extract_low_ld_snps = "plink --bfile {} --extract {} --make-bed --out {}".format("/home/vcabeli/Documents/data/BP_final/BP.B37-final",
                                                                                             "L2_thresh_{}_BP.extract".format(thresh_LD),
                                                                                             "low_ld")
        print cmd_extract_low_ld_snps
        p = subprocess.Popen(cmd_extract_low_ld_snps, shell=True)
        assert(p.wait()==0)

        plink_prune("low_ld")



        # Covariables
        cmd_smart_pca2 = "smartpca.perl -i {} -a {} -b {} -o {} -p {} -e {} -l {} -k {} -t {} -m {}".format("low_ld_pruned.bed",
                                                                                                           "low_ld_pruned.bim",
                                                                                                           "low_ld_pruned.fam",
                                                                                                           "low_ld_pruned.pca",
                                                                                                           "low_ld_pruned.plot",
                                                                                                           "low_ld_pruned.eval",
                                                                                                           "low_ld_pruned.log",
                                                                                                           nb_pc, 2, 0)
        print cmd_smart_pca2
        p = subprocess.Popen(cmd_smart_pca2, shell=True)
        p.wait()


        write_cov_file("low_ld_pruned.pca.evec",2)
        save_strat("low_ld_pruned.pca.evec", "low_ld{}".format(thresh_LD)) #OUTPUT LINE

        # Split training / testing datasets
        indx=np.random.choice(snp_data.shape[0],nb_snps)
        val_pheno={}
        liste=[]
        for i in indx:
            liste.append(pheno['vals'][i])
        val_pheno['vals'] = liste
        copi_data=snp_data[indx,:]
        
        train_idces = np.random.choice(np.arange(copi_data.row_count), size=int(copi_data.row_count*0.5), replace=False)
        print len(train_idces)
        test_idces = np.setdiff1d(np.arange(copi_data.row_count), train_idces, assume_unique=True)
        print len(test_idces)


        training_sample_out = "training_samples.keep"

        with open(training_sample_out, 'w') as f:
            for i in train_idces:
                f.write(pheno['iid'][i][0] + "\t" + pheno['iid'][i][1])
                f.write("\n")


        # #### Coded genotype (hardy-weinberg)
        # 
        # With sufficiently large cohorts, traning and test sets should have the same MAFs. The coded genotype is computed on the whole dataset.



        G_hw = snp_data.val.copy()

        MAFs = np.nansum(2-G_hw, axis=0, ) / (np.count_nonzero(~np.isnan(G_hw), axis=0) * 2)
        G_hw = (2-G_hw - 2*MAFs)/np.sqrt(2*MAFs*(1-MAFs))


        # #### Build training set bed bim fam



        cmd_keep_plink = "plink --bfile {} --keep {} --make-bed --out {}".format("/home/vcabeli/Documents/data/BP_final/BP.B37-final",
                                                                                 training_sample_out,
                                                                                 "training_set")
        print cmd_keep_plink
        p = subprocess.Popen(cmd_keep_plink, shell=True)
        p.wait()




        cmd_second_gwas = "plink --bfile {} --logistic sex beta hide-covar --covar {} --covar-name PC1-PC{} --out {}".format("training_set",
                                                                                                           "low_ld_pruned.cov",
                                                                                                            nb_pc,
                                                                                                           "low_ld")
        print cmd_second_gwas
        p = subprocess.Popen(cmd_second_gwas, shell=True)
        p.wait()
        low_ld_res = pd.read_table("low_ld.assoc.logistic", delim_whitespace=True)



        # ### Compute PRS





        sorted_snps_low_ld = np.argsort(low_ld_res.P)
        sorted_snps_low_ld = clump_sorted_snps("training_set", "low_ld.assoc.logistic",
                                               copi_data.col, sorted_snps_low_ld)


        prs_low_ld = polygen_score_sign(G_hw, sorted_snps_low_ld,
                                        threshs,
                                        test_idces, val_pheno,
                                        low_ld_res.BETA)
        res_hw_pval_sign[iteration] = prs_low_ld
    return(res_hw_pval_sign)


# In[ ]:

snp_data, pheno = util.load_data("/home/vcabeli/Documents/data/BP_final/BP.B37-final")
n_iter = 2
threshs = range(10, 70000, 140)
thresh_LD = 5
nb_pc_tot=5
time = range(10,70000,140)
time_pc = range(2,nb_pc_tot)
col='bgrcmk'

############"
auc_50 = np.zeros([nb_pc_tot-2,n_iter,len(threshs)])
auc_100 = np.zeros([nb_pc_tot-2,n_iter,len(threshs)])
auc_200 = np.zeros([nb_pc_tot-2,n_iter,len(threshs)])
auc_500 = np.zeros([nb_pc_tot-2,n_iter,len(threshs)])
auc_1000= np.zeros([nb_pc_tot-2,n_iter,len(threshs)])
auc_2000= np.zeros([nb_pc_tot-2,n_iter,len(threshs)])
###########PC
for PC in range(2,nb_pc_tot):
    print "========================= {} PCs =========================".format(PC)
    print "========================= {} SNPs =========================".format(50)
    auc_50[PC-2]=iterate_bis(snp_data, pheno,n_iter,threshs,thresh_LD,PC,50)
    print "========================= {} SNPs =========================".format(100)
    auc_100[PC-2]=iterate_bis(snp_data, pheno,n_iter,threshs,thresh_LD,PC,100)
    print "========================= {} SNPs =========================".format(200)
    auc_200[PC-2]=iterate_bis(snp_data, pheno,n_iter,threshs,thresh_LD,PC,200)
    print "========================= {} SNPs =========================".format(500)
    auc_500[PC-2]=iterate_bis(snp_data, pheno,n_iter,threshs,thresh_LD,PC,500)
    print "========================= {} SNPs =========================".format(1000)
    auc_1000[PC-2]=iterate_bis(snp_data, pheno,n_iter,threshs,thresh_LD,PC,1000)
    print "========================= {} SNPs =========================".format(snp_data.shape[0])
    auc_2000[PC-2]=iterate_bis(snp_data, pheno,n_iter,threshs,thresh_LD,PC,snp_data.shape[0])



sns.tsplot(auc_50.mean(axis=(1,2)),time_pc,color=col[0],ci=[0,95], condition="SNP 50")
sns.tsplot(auc_100.mean(axis=(1,2)),time_pc,color=col[1],ci=[0,95], condition="SNP 100")
sns.tsplot(auc_200.mean(axis=(1,2)),time_pc,color=col[2],ci=[0,95], condition="SNP 200")
sns.tsplot(auc_500.mean(axis=(1,2)),time_pc,color=col[3],ci=[0,95], condition="SNP 500")
sns.tsplot(auc_1000.mean(axis=(1,2)),time_pc,color=col[4],ci=[0,95], condition="SNP 1000")
sns.tsplot(auc_2000.mean(axis=(1,2)),time_pc,color=col[5],ci=[0,95], condition="SNP 2000")
sns.plt.ylabel('AUC')
sns.plt.xlabel('N PCs')
sns.plt.savefig("fct_SNP.pdf", format="pdf")

print "Wrote {}.".format("fct_SNP.pdf")
plt.show()


# In[23]:

sns.tsplot(auc_50.mean(axis=(1,2)),time_pc,color=col[0],ci=[0,95], condition="SNP 50")
sns.tsplot(auc_100.mean(axis=(1,2)),time_pc,color=col[1],ci=[0,95], condition="SNP 100")
sns.tsplot(auc_200.mean(axis=(1,2)),time_pc,color=col[2],ci=[0,95], condition="SNP 200")
sns.tsplot(auc_500.mean(axis=(1,2)),time_pc,color=col[3],ci=[0,95], condition="SNP 500")
sns.tsplot(auc_1000.mean(axis=(1,2)),time_pc,color=col[4],ci=[0,95], condition="SNP 1000")
sns.tsplot(auc_2000.mean(axis=(1,2)),time_pc,color=col[5],ci=[0,95], condition="SNP 2000")
sns.plt.ylabel('AUC')
sns.plt.xlabel('N PCs')
sns.plt.savefig("fct_SNP.pdf", format="pdf")

print "Wrote {}.".format("fct_SNP.pdf")
plt.show()
