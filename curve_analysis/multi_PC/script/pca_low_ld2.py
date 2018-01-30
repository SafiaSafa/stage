# coding: utf-8

# 

# In[1]:

import sys
sys.path.append("/home/vcabeli/Documents/scripts/python") # Contains "util" and "Polygenic_score" packages
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



def save_strat(strat_file, strat):
    pca_res = pd.read_table(strat_file, delim_whitespace=True, skiprows=1, header=None)
    PCs_out = strat + ".PCs"
    with open("res_clumped_raw_genotype/"+PCs_out, "a") as f:
        f.write(', '.join(pca_res[0].values))
        f.write('\n')
        f.write(str(pca_res[2].values.tolist())[1:-1])
        f.write('\n')
        f.write(str(pca_res[3].values.tolist())[1:-1])
        f.write('\n')



# In[2]:

os.chdir("/home/vcabeli/Documents/analyses_current/pca_low_ld/")


# In[3]:

snp_data, pheno = util.load_data("/home/vcabeli/Documents/data/BP_final/BP.B37-final")


# #### Split training / testing datasets

# In[4]:

train_idces = np.random.choice(np.arange(snp_data.row_count), size=int(snp_data.row_count*0.5), replace=False)
test_idces = np.setdiff1d(np.arange(snp_data.row_count), train_idces, assume_unique=True)


# In[5]:

training_sample_out = "training_samples.keep"

with open(training_sample_out, 'w') as f:
    for i in train_idces:
        f.write(pheno['iid'][i][0] + "\t" + pheno['iid'][i][1])
        f.write("\n")


# #### Coded genotype (hardy-weinberg)
# 
# With sufficiently large cohorts, traning and test sets should have the same MAFs. The coded genotype is computed on the whole dataset.

# In[6]:

G_hw = snp_data.val.copy()

MAFs = np.nansum(2-G_hw, axis=0, ) / (np.count_nonzero(~np.isnan(G_hw), axis=0) * 2)
#G_hw = (2-G_hw - 2*MAFs)/np.sqrt(2*MAFs*(1-MAFs))


# #### Build training set bed bim fam

# In[7]:

cmd_keep_plink = "plink --bfile {} --keep {} --make-bed --out {}".format("/home/vcabeli/Documents/data/BP_final/BP.B37-final",
                                                                         training_sample_out,
                                                                         "training_set")
print cmd_keep_plink
p = subprocess.Popen(cmd_keep_plink, shell=True)
p.wait()


# In[8]:

def plink_prune(plink_bfile):

    cmd_prune = """plink --bfile {}      --exclude /home/vcabeli/Documents/data/high-LD-regions_37.txt           --range --indep-pairwise 50 5 0.2           --allow-extra-chr           --out {}""".format(plink_bfile, plink_bfile+"_prune")
    p = subprocess.Popen(cmd_prune, shell=True)
    assert(p.wait() == 0)

    cmd_extract = """plink --bfile {}           --extract {}           --allow-extra-chr           --chr 1-23           --make-bed           --out {}""".format(plink_bfile, plink_bfile + "_prune.prune.in", 
                             plink_bfile + "_pruned")
    p = subprocess.Popen(cmd_extract, shell=True)
    assert(p.wait() == 0)

    print "Wrote {} bed/bim/fam.".format(plink_bfile + "_pruned")



def independent_snps(snps_list):
    with open("non_pruned_list", 'w') as f:
        for snp in snps_list:
            f.write(snp)
            f.write("\n")
    cmd_extract = "plink --bfile {} --extract {} --make-bed --out {}".format("training_set",
                                                                             "non_pruned_list",
                                                                             "snp_list")
    print cmd_extract
    p = subprocess.Popen(cmd_extract, shell=True)
    assert(p.wait()==0)

    plink_prune("snp_list")

    t = pd.read_table("snp_list_pruned.bim", delim_whitespace=True, header=None)
    return t[1].values




def write_cov_file(pca_file):

    cmd_sed = "sed -i -e \"s/:/ /g\" -e \"s/\s\+/\t/g\" -e \"s/^\s\+//g\" {}".format(pca_file)
    p = subprocess.Popen(cmd_sed, shell=True)
    assert(p.wait()==0)
    
    cov_file = pca_file[:-9]+".cov"
    cmd_print = "printf \"FID\tIID\tPC1\tPC2\n\" > {}".format(cov_file)
    p = subprocess.Popen(cmd_print, shell=True)
    assert(p.wait()==0)

    cmd_tail = "tail -n+2 {} | cut -f1-4 >> {}".format(pca_file, cov_file)
    p = subprocess.Popen(cmd_tail, shell=True)
    assert(p.wait()==0)
    
    print "Wrote {}".format(cov_file)
    

def plot_pca(cov_file) :
    
    training_set_pca = pd.read_table(cov_file)
    pheno_df = pd.DataFrame({'IID' : [pheno['iid'][i][1] for i in range(len(pheno['iid']))], 'pheno' : pheno['vals']})
    training_set_pca = training_set_pca.merge(pheno_df, on="IID", how="inner")
    
    plt.scatter(data=training_set_pca[training_set_pca.pheno == 1], x="PC1", y="PC2", 
                label="controls", s=15)
    plt.scatter(data=training_set_pca[training_set_pca.pheno == 2], x="PC1", y="PC2",
                label="cases", s=15)
    plt.legend()
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
    



threshs = range(10, 70000, 300)
thresh_LD = 5



cmd_extract_low_ld_snps = "plink --bfile {} --extract {} --make-bed --out {}".format("training_set",
                                                                                     "L2_thresh_{}_BP.extract".format(thresh_LD),
                                                                                     "low_ld")
print cmd_extract_low_ld_snps
p = subprocess.Popen(cmd_extract_low_ld_snps, shell=True)
assert(p.wait()==0)

plink_prune("low_ld")


# In[16]:

cmd_smart_pca2 = "smartpca.perl -i {} -a {} -b {} -o {} -p {} -e {} -l {} -k {} -t {} -m {}".format("low_ld_pruned.bed",
                                                                                                   "low_ld_pruned.bim",
                                                                                                   "low_ld_pruned.fam",
                                                                                                   "low_ld_pruned.pca",
                                                                                                   "low_ld_pruned.plot",
                                                                                                   "low_ld_pruned.eval",
                                                                                                   "low_ld_pruned.log",
                                                                                                   2, 2, 0)
print cmd_smart_pca2
p = subprocess.Popen(cmd_smart_pca2, shell=True)
p.wait()


# In[17]:

write_cov_file("low_ld_pruned.pca.evec")
save_strat("low_ld_pruned.pca.evec", "low_ld{}".format(thresh_LD)) #OUTPUT LINE
#plot_pca("low_ld_pruned.cov")



# In[18]:

cmd_second_gwas = "plink --bfile {} --logistic sex beta hide-covar --covar {} --covar-name PC1-PC2 --out {}".format("training_set",
                                                                                                   "low_ld_pruned.cov",
                                                                                                   "low_ld")
print cmd_second_gwas
p = subprocess.Popen(cmd_second_gwas, shell=True)
p.wait()

low_ld_res = pd.read_table("low_ld.assoc.logistic", delim_whitespace=True)



# ### Compute PRS

# In[21]:

sorted_snps_low_ld = np.argsort(low_ld_res.P)
sorted_snps_low_ld = clump_sorted_snps("training_set", "low_ld.assoc.logistic",
                                       snp_data.col, sorted_snps_low_ld)


prs_low_ld = polygen_score_sign(G_hw, sorted_snps_low_ld,
                                threshs,
                                test_idces, pheno,
                                low_ld_res.BETA)
with open("res_clumped_raw_genotype/prs_low_ld_{}.res".format(thresh_LD), 'a') as f: #OUTPUT LINE
    np.savetxt(f, prs_low_ld)




### RANDOM

MAFs = pd.read_table("/home/vcabeli/Documents/data/BP_final/BP.B37-final_freqs.frq", delim_whitespace=True)
low_ld_snps = pd.read_table("L2_thresh_{}_BP.extract".format(thresh_LD), header=None)
low_ld_MAFs = MAFs[MAFs.SNP.isin(low_ld_snps.unstack().values)].MAF

low_ld_snps = independent_snps(low_ld_snps[0].values)

random_snps = np.random.choice(MAFs.SNP, len(low_ld_MAFs)*2, replace=False)
random_snps = np.random.choice(independent_snps(random_snps), len(low_ld_snps), replace=False)


with open("random_snps.extract", 'w') as f:
    for snp in random_snps:
        f.write(snp)
        f.write('\n')


cmd_extract_low_ld_snps = "plink --bfile {} --extract {} --make-bed --out {}".format("training_set",
                                                                                     #"L2_thresh_{}_BP.extract".format(thresh_LD),
                                                                                     "random_snps.extract",
                                                                                     "low_ld")
print cmd_extract_low_ld_snps
p = subprocess.Popen(cmd_extract_low_ld_snps, shell=True)
assert(p.wait()==0)

plink_prune("low_ld")


# In[16]:

cmd_smart_pca2 = "smartpca.perl -i {} -a {} -b {} -o {} -p {} -e {} -l {} -k {} -t {} -m {}".format("low_ld_pruned.bed",
                                                                                                   "low_ld_pruned.bim",
                                                                                                   "low_ld_pruned.fam",
                                                                                                   "low_ld_pruned.pca",
                                                                                                   "low_ld_pruned.plot",
                                                                                                   "low_ld_pruned.eval",
                                                                                                   "low_ld_pruned.log",
                                                                                                   2, 2, 0)
print cmd_smart_pca2
p = subprocess.Popen(cmd_smart_pca2, shell=True)
p.wait()


# In[17]:

write_cov_file("low_ld_pruned.pca.evec")
save_strat("low_ld_pruned.pca.evec", "random_{}".format(thresh_LD)) #OUTPUT LINE
#plot_pca("low_ld_pruned.cov")



# In[18]:

cmd_second_gwas = "plink --bfile {} --logistic sex beta hide-covar --covar {} --covar-name PC1-PC2 --out {}".format("training_set",
                                                                                                   "low_ld_pruned.cov",
                                                                                                   "low_ld")
print cmd_second_gwas
p = subprocess.Popen(cmd_second_gwas, shell=True)
p.wait()

low_ld_res = pd.read_table("low_ld.assoc.logistic", delim_whitespace=True)



# ### Compute PRS

# In[21]:

sorted_snps_low_ld = np.argsort(low_ld_res.P)
sorted_snps_low_ld = clump_sorted_snps("training_set", "low_ld.assoc.logistic",
                                       snp_data.col, sorted_snps_low_ld)


prs_low_ld = polygen_score_sign(G_hw, sorted_snps_low_ld,
                                threshs,
                                test_idces, pheno,
                                low_ld_res.BETA)
with open("res_clumped_raw_genotype/prs_random_{}.res".format(thresh_LD), 'a') as f: #OUTPUT LINE
    np.savetxt(f, prs_low_ld)






### RANDOM MATCH

def find_quantile(quantiles, value):
    idx = (np.abs(quantiles-value)).argmin()
    return idx


def sample_match(data, target_data, size, quantiles=100):
    target_quantiles = target_data.quantile(np.linspace(0,1,quantiles)).values
    quantile_map = {row.SNP : find_quantile(target_quantiles, row.MAF) for row in data.itertuples()}
    len_quantiles = np.bincount(quantile_map.values())
    
    prop = [1./(quantiles)/len_quantiles[quantile_map[SNP]] for SNP in data.SNP]
    
    random_sample = np.random.choice(data.SNP, size, False, p=prop)
    return random_sample


MAFs = pd.read_table("/home/vcabeli/Documents/data/BP_final/BP.B37-final_freqs.frq", delim_whitespace=True)
low_ld_snps = pd.read_table("L2_thresh_{}_BP.extract".format(thresh_LD), header=None)
low_ld_MAFs = MAFs[MAFs.SNP.isin(low_ld_snps.unstack().values)].MAF

random_snps = sample_match(MAFs, low_ld_MAFs, len(low_ld_MAFs)*2)

low_ld_snps = independent_snps(low_ld_snps[0].values)
random_snps = np.random.choice(independent_snps(random_snps), len(low_ld_snps), replace=False)




with open("random_snps.extract", 'w') as f:
    for snp in random_snps:
        f.write(snp)
        f.write('\n')


cmd_extract_low_ld_snps = "plink --bfile {} --extract {} --make-bed --out {}".format("training_set",
                                                                                     #"L2_thresh_{}_BP.extract".format(thresh_LD),
                                                                                     "random_snps.extract",
                                                                                     "low_ld")
print cmd_extract_low_ld_snps
p = subprocess.Popen(cmd_extract_low_ld_snps, shell=True)
assert(p.wait()==0)

plink_prune("low_ld")


# In[16]:

cmd_smart_pca2 = "smartpca.perl -i {} -a {} -b {} -o {} -p {} -e {} -l {} -k {} -t {} -m {}".format("low_ld_pruned.bed",
                                                                                                   "low_ld_pruned.bim",
                                                                                                   "low_ld_pruned.fam",
                                                                                                   "low_ld_pruned.pca",
                                                                                                   "low_ld_pruned.plot",
                                                                                                   "low_ld_pruned.eval",
                                                                                                   "low_ld_pruned.log",
                                                                                                   2, 2, 0)
print cmd_smart_pca2
p = subprocess.Popen(cmd_smart_pca2, shell=True)
p.wait()


# In[17]:

write_cov_file("low_ld_pruned.pca.evec")
save_strat("low_ld_pruned.pca.evec", "random_matched_{}".format(thresh_LD)) #OUTPUT LINE
#plot_pca("low_ld_pruned.cov")



# In[18]:

cmd_second_gwas = "plink --bfile {} --logistic sex beta hide-covar --covar {} --covar-name PC1-PC2 --out {}".format("training_set",
                                                                                                   "low_ld_pruned.cov",
                                                                                                   "low_ld")
print cmd_second_gwas
p = subprocess.Popen(cmd_second_gwas, shell=True)
p.wait()

low_ld_res = pd.read_table("low_ld.assoc.logistic", delim_whitespace=True)



# ### Compute PRS

# In[21]:

sorted_snps_low_ld = np.argsort(low_ld_res.P)
sorted_snps_low_ld = clump_sorted_snps("training_set", "low_ld.assoc.logistic",
                                       snp_data.col, sorted_snps_low_ld)


prs_low_ld = polygen_score_sign(G_hw, sorted_snps_low_ld,
                                threshs,
                                test_idces, pheno,
                                low_ld_res.BETA)
with open("res_clumped_raw_genotype/prs_random_matched_{}.res".format(thresh_LD), 'a') as f: #OUTPUT LINE
    np.savetxt(f, prs_low_ld)
