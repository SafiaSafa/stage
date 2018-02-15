#!/usr/bin/env python2
from pysnptools.snpreader.bed import Bed
import pysnptools.util.pheno
import numpy as np
import pandas as pd
import util
import sys

BP_DIR = "/home/vcabeli/Documents/data/BP_final_imputed/0.5_info/merged/plink/bedbimfam"


for chrom in range(1,23):
    print("CHROMOSOME {}".format(chrom))
    snp_data, pheno = util.load_data(BP_DIR + "/BP.B37-imputed.chr{}.merged.1KG_P3.EURMAF0.005.0.5_info".format(chrom))
    i_scores_classification = util.get_i_scores_binary()

    df_iscore = pd.DataFrame({'SNP':snp_data.col, 'i_score':i_scores_classification})
    df_iscore.to_csv("/home/vcabeli/Documents/analyses_current/assoc/imputed/classification/iscore_chr{}".format(chrom), sep=" ", index=False)
