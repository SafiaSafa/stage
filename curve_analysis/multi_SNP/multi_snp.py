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


# In[218]:

def select_snp(nb_snps):
   snp_gene_table = pd.read_table("/home/vcabeli/Documents/data/BP_final/BP.B37-final.bim",
                                  delim_whitespace=True,names=["chr","snp","distance","position","Allele_1","Allele_2"])


   indx=np.random.choice(snp_gene_table.shape[0],nb_snps)
   snp_extract=snp_gene_table.snp[indx]
   #ecriture du fichier contenant la liste des snp a garder 
   extract_file=str(nb_snps)+"_snps.extract"
   set_file=str(nb_snps)+"_snps_set"
   with open(extract_file, 'w') as f:
       for snp in snp_extract:
               f.write(snp)
               f.write("\n")
   #ecriture fichier bim/bed/fam contenant 50 snps
   cmd_keep_plink = "plink --bfile {} --extract {} --make-bed --out {}".format("/home/vcabeli/Documents/data/BP_final/BP.B37-final",
                                                                                    extract_file,
                                                                                    set_file)
   print cmd_keep_plink
   p = subprocess.Popen(cmd_keep_plink, shell=True)
   p.wait()
   return(set_file)


# In[286]:

def iterate(n_iter,threshs,nb_pc,data):
    res_hw_pval_sign = np.zeros([n_iter, len(threshs)])
    for iteration in range(n_iter):

        print "=========================Iteration {}=========================".format(iteration)
        # extraction des SNPs indépendantes selon le seuil de LD 
        cmd_extract_low_ld_snps = "plink --bfile {} --extract {} --make-bed --out {}".format(data,
                                                                                             "L2_thresh_{}_BP.extract".format(thresh_LD),
                                                                                             "low_ld")
        print cmd_extract_low_ld_snps
        p = subprocess.Popen(cmd_extract_low_ld_snps, shell=True)
        assert(p.wait()==0)

        plink_prune("low_ld")



        # Covariables
        if(nb_pc!=0):
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


            write_cov_file("low_ld_pruned.pca.evec",nb_pc)
            save_strat("low_ld_pruned.pca.evec", "low_ld{}".format(thresh_LD)) #OUTPUT LINE


        # Split training / testing datasets

        snp_data, pheno = util.load_data("/home/vcabeli/Documents/data/BP_final/BP.B37-final")

        train_idces = np.random.choice(np.arange(snp_data.row_count), size=int(snp_data.row_count*0.5), replace=False)

        test_idces = np.setdiff1d(np.arange(snp_data.row_count), train_idces, assume_unique=True)
 



        training_sample_out = "training_samples.keep"

        with open(training_sample_out, 'w') as f:
            for i in train_idces:
                f.write(pheno['iid'][i][0] + "\t" + pheno['iid'][i][1])
                f.write("\n")

        # #### Build training set bed bim fam

        cmd_keep_plink = "plink --bfile {} --keep {} --make-bed --out {}".format("/home/vcabeli/Documents/data/BP_final/BP.B37-final",
                                                                                 training_sample_out,
                                                                                 "training_set")
        print cmd_keep_plink
        p = subprocess.Popen(cmd_keep_plink, shell=True)
        p.wait()
        # Association
        if (nb_pc==0):
            cmd_second_gwas = "plink --bfile {} --logistic sex beta hide-covar --out {}".format("training_set","low_ld")
        elif(nb_pc==1):
            cmd_second_gwas = "plink --bfile {} --logistic sex beta hide-covar --covar {} --covar-name PC1 --out {}".format("training_set",
                                                                                                           "low_ld_pruned.cov",
                                                                                                           "low_ld")
        else:
            cmd_second_gwas = "plink --bfile {} --logistic sex beta hide-covar --covar {} --covar-name PC1-PC{} --out {}".format("training_set",
                                                                                                           "low_ld_pruned.cov",
                                                                                                            nb_pc,
                                                                                                           "low_ld")
        print cmd_second_gwas
        p = subprocess.Popen(cmd_second_gwas, shell=True)
        p.wait()
        low_ld_res = pd.read_table("low_ld.assoc.logistic", delim_whitespace=True)

        # ### Compute PRS


        # Coded genotype (hardy-weinberg)
        G_hw = snp_data.val.copy()

        MAFs = np.nansum(2-G_hw, axis=0, ) / (np.count_nonzero(~np.isnan(G_hw), axis=0) * 2)
        G_hw = (2-G_hw - 2*MAFs)/np.sqrt(2*MAFs*(1-MAFs))



        sorted_snps_low_ld = np.argsort(low_ld_res.P)
        sorted_snps_low_ld = clump_sorted_snps("training_set", "low_ld.assoc.logistic",snp_data.col, sorted_snps_low_ld)

        prs_low_ld = polygen_score_sign(G_hw, sorted_snps_low_ld,
                                        threshs,
                                        test_idces, pheno,
                                        low_ld_res.BETA)
        res_hw_pval_sign[iteration] = prs_low_ld
    return(res_hw_pval_sign)


# In[284]:

data_50=select_snp(50)
data_100=select_snp(100)
data_200=select_snp(200)
data_500=select_snp(500)
data_1000=select_snp(1000)
data_2000=select_snp(2000)
nb_pc_tot=11
threshs = range(10, 95000, 1000)
n_iter=10
time_pc = range(0,nb_pc_tot)
get_color = lambda : "#" + "".join(np.random.choice(list("02468acef"), size=6))
auc_50 = np.zeros([nb_pc_tot,n_iter,len(threshs)])
auc_100 = np.zeros([nb_pc_tot,n_iter,len(threshs)])
auc_200 = np.zeros([nb_pc_tot,n_iter,len(threshs)])
auc_500 = np.zeros([nb_pc_tot,n_iter,len(threshs)])
auc_1000 = np.zeros([nb_pc_tot,n_iter,len(threshs)])
auc_2000 = np.zeros([nb_pc_tot,n_iter,len(threshs)])


# ### Calcul des covariables pour 50 snps

# In[288]:

for PC in range(0,nb_pc_tot):
        print "========================= {} PCs =========================".format(PC)
        auc_50[PC]=iterate(n_iter,threshs,PC,data_50)
        auc_100[PC]=iterate(n_iter,threshs,PC,data_100)
        auc_200[PC]=iterate(n_iter,threshs,PC,data_200)
        auc_500[PC]=iterate(n_iter,threshs,PC,data_500)
        auc_1000[PC]=iterate(n_iter,threshs,PC,data_1000)
        auc_2000[PC]=iterate(n_iter,threshs,PC,data_2000)


# In[291]:

sns.tsplot(auc_50.mean(axis=(1,2)),time_pc,color=get_color(),ci=[0,95], condition="SNP 50")
sns.tsplot(auc_100.mean(axis=(1,2)),time_pc,color=get_color(),ci=[0,95], condition="SNP 100")
sns.tsplot(auc_200.mean(axis=(1,2)),time_pc,color=get_color(),ci=[0,95], condition="SNP 200")
sns.tsplot(auc_500.mean(axis=(1,2)),time_pc,color=get_color(),ci=[0,95], condition="SNP 500")
sns.tsplot(auc_1000.mean(axis=(1,2)),time_pc,color=get_color(),ci=[0,95], condition="SNP 1000")
sns.tsplot(auc_2000.mean(axis=(1,2)),time_pc,color=get_color(),ci=[0,95], condition="SNP 2000")
sns.plt.ylim(0.4,1)
sns.plt.ylabel('AUC')
sns.plt.xlabel('N PCs')
sns.plt.savefig("fct_pc.pdf", format="pdf")

print "Wrote {}.".format("fct_pc.pdf")
plt.show()