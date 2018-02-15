
BP_DIR="/home/vcabeli/Documents/data/BP_final"

for BP_type in "BP1" "BP2"
do
    gcta64 --bfile ${BP_DIR}/AAO \ #Input PLINK binary PED files, e.g. test.fam, test.bim and test.bed (see PLINK user manual for details).
           --pheno ${BP_DIR}/AAO.pheno_warpedlmm_pheno.txt \ # phenotype file 
           --mlma \ #mixed linear model association analysis;
           --autosome \ #Include SNPs on all of the autosomes in the analysis.
           --keep ${BP_DIR}/${BP_type}_samples.keep \ # list of individuals to be included in the analysis.
           --qcovar ${BP_DIR}/plink.cov \ #Input quantitative covariates from a plain text file. Each quantitative covariate is recognized as a continuous variable.    
           --thread-num 7 \ # number of threads on which the program will be running.
           --out gcta_mlma_warped_${BP_type} #Specify output root filename

    Rscript ../../../plot_manhattan_qq_gcta.R gcta_mlma_warped_${BP_type}.mlma
done


'''
--mlma
This option will initiate an MLM based association analysis including the candidate SNP
y = a + bx + g + e
where y is the phenotype, a is the mean term, b is the additive effect (fixed effect) of the candidate SNP to be tested for association,
x is the SNP genotype indicator variable coded as 0, 1 or 2, g is the polygenic effect (random effect) i.e. the accumulated effect of all
SNPs (as captured by the GRM calculated using all SNPs) and e is the residual. For the ease of computation, the genetic variance, var(g),
is estimated based on the null model i.e. y = a + g + e and then fixed while testing for the association between each SNP and the trait.
This analysis would be similar as that implemented in other software tools such as EMMAX, FaST-LMM and GEMMA. The results will be saved in the *.mlma file.
'''