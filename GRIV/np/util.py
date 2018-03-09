#!/usr/bin/env python2
"""
Vincent Cabeli 02/2017 GBA


IDEAS:
	Compute I score from two matrices? snp_reader.read[case_idces] / [control_idces]
"""

import pysnptools.util.pheno
from pysnptools.snpreader.bed import Bed
import numpy as np
import scipy.stats
import pandas as pd
import multiprocessing
import scipy.sparse

#==============================================================================#
cpus = multiprocessing.cpu_count() #8

# Global vars for multiprocessing i score computation. Set during load_data().
snp_data = None
permute_flag = False

case_idces = None
control_idces = None

cont_idces = None
pheno_vals = None
pheno_vals_total_deviation = None
pheno_vals_mean = None

epistasy_matrix = None
idxg = None

#==============================================================================#

def load_data(snp_file, pheno="fam", pheno_type="binary", covar_file=None,
	      sample_size=1.0, permute=False, pheno_col=-1):
    """
        Loads genotype and phenotype data with the pysnptools library. The user
    can specify either a continuous or case/control phenotype, a fraction of
    the sample size to be randomly drawn, or if the phenotypes are to be
    randomly permuted.

        Args:
            snp_file (str) : The Bed/Bim/Fam common file name.
            pheno_file (optional str) : A file containing iids and phenotypes.
            covar_file (optional str) : A file containing iids and covariables.
            sample_size (optional float) : A float in [0,1], the fraction of the
        whole dataset to be sampled.
            permute (optional boolean) : Is the phenotype to be permuted ?
            pheno_col (optional int) : The column index of the phenotype file
        that contains phenotype values (-1 for last column).

        Returns:
            snp_data, a Snpdata object ,
            pheno, a Pheno object,
    """
    reset_globals()
    # Locate snp data on disk
    snp_reader = Bed(snp_file, count_A1=False)

    # Load phenotype
    if pheno == "fam":
        pheno = pysnptools.util.pheno.loadPhen(snp_file + ".fam")
        if pheno_col == -1: # Default column
        	print("Using last column in .fam file as phenotype.")
    else:
        pheno = pysnptools.util.pheno.loadPhen(pheno)
    pheno['vals'] = [vals[pheno_col] for vals in pheno['vals']]

    # Remove nans in pheno
    pheno_idces = np.where(~ np.isnan(pheno['vals']))[0]
    for key in ['iid', 'vals']:
        pheno[key] = [value for i, value in enumerate(pheno[key]) if i in pheno_idces]
    # Draw random sample ?
    assert sample_size > 0.0 and sample_size <= 1.0, ("Please specify the size",
    	"of the random sample as a ratio of the total size (between 0 and 1)")
    if sample_size < 1.0 :
    	pheno, pheno_idces = draw_sample(pheno, pheno_idces, sample_size)
    # Set global variables for multithreading
    assert pheno_type in ["binary", "continuous"], ("pheno_type {} not understood. Please specifiy either binary or continuous.".format(pheno_type))
    #if permute:
    #	pheno['vals'] = np.random.permutation(pheno['vals'])
    if pheno_type == "binary":
        set_global_idces(pheno)
    else:
        set_global_pheno(pheno)

    ## Load covariates
    #if covar_file is not None:
    #    covar = pysnptools.util.pheno.loadPhen(covar_file)
    #    snp_reader, pheno, covar = srutil.intersect_apply([snp_reader, pheno, covar])
    #    covar = covar['vals']
    #else:
    #    snp_reader, pheno = srutil.intersect_apply([snp_reader, pheno])
    #    covar = None

    # Load SNP data for which there is a valid phenotype
    global snp_data
    if permute:
    	pheno_idces = np.random.permutation(pheno_idces)
    snp_data = snp_reader[pheno_idces, :].read(dtype=np.float32)
    print("Loaded {} samples and {} SNPs.".format(snp_data.row_count, snp_data.col_count))

    assert (((not permute) and np.all(pheno['iid'] == snp_data.iid)) or permute), "The samples are not sorted!"

    return snp_data, pheno


def set_global_pheno(pheno):
    """
        Set global variables for multiprocessing. This function sets pheno_vals,
    a numpy array of phenotype values ordered according to iids. This function is
    called during load_data().

        Args:
            pheno (dict) : A phenotype dictionnary as created by
        pysnptools.pheno.loadPhen

        Returns:
            Nothing.
    """
    print("\tConsidering phenotype as a continuous variable.")
    global cont_idces
    cont_idces = np.s_[0:len(pheno['vals'])]
    global pheno_vals
    pheno_vals = np.array(pheno['vals']).flatten()
    global pheno_vals_total_deviation
    pheno_vals_total_deviation = np.sum((pheno_vals - pheno_vals.mean())**2)
    global pheno_vals_mean
    pheno_vals_mean = pheno_vals.mean()


def set_global_idces(pheno):
    """
        Set global variables for multiprocessing. This function sets case_idces and
    control_idces, which respectively designate Case and Control indices in the list
    of iids. The use of a numpy slice instead of list of indices divides computing
    time by two. If the samples are ordered : controls then cases or cases then
    controls, does the conversion. This function is called during load_data().

        Args:
            pheno (dict) : A phenotype dictionnary as created by
        pysnptools.pheno.loadPhen

        Returns:
            Nothing.
    """
    print("\tConsidering phenotype as a plink binary 1/2.")

    global case_idces
    case_idces = [i for i in range(0, len(pheno['vals'])) if pheno['vals'][i] == 2]
    n_case_idces = len(case_idces)

    if case_idces == range(min(case_idces), max(case_idces) + 1):
        case_idces = np.s_[min(case_idces):(max(case_idces) + 1)]
        n_case_idces = case_idces.stop - case_idces.start

    global control_idces
    control_idces = [i for i in range(0, len(pheno['vals'])) if pheno['vals'][i] == 1]
    n_control_idces = len(control_idces)

    if control_idces == range(min(control_idces), max(control_idces) + 1):
        control_idces = np.s_[min(control_idces):(max(control_idces) + 1)]
        n_control_idces = control_idces.stop - control_idces.start

    print("\tFound {} cases and {} controls.".format(n_case_idces,
                                                     n_control_idces))


def get_permuted_I_scores(snp_file, nPermuts, I_score=None, pheno="fam", 
                          pheno_type="binary", pheno_col=-1):
    """

    """

    assert pheno_type in ["binary", "continuous"], ("pheno_type {} not understood. Please",
                                                    " specifiy either binary or continuous.".format(pheno_type))
    snp_data, pheno_obs = load_data(snp_file, pheno, pheno_type=pheno_type, 
                                    pheno_col = pheno_col)
    i_scores = get_i_scores_binary() if pheno_type == "binary" else get_i_scores_continuous()
    permute_df = pd.DataFrame({'Observed_iscore' : i_scores})

    for i in tqdm(range(nPermuts)):
        snp_data, pheno_permute = load_data(snp_file, pheno, pheno_type=pheno_type, 
                                            pheno_col=pheno_col, permute=True)
        i_scores = get_i_scores_binary() if pheno_type == "binary"  else get_i_scores_continuous()

        permute_df[('Permute_' + str(i))] = i_scores

    return permute_df


def get_i_scores_continuous():
    """"
        Computes and returns I scores for a continuous phenotype.

    Needs the following global objects:
        snp_data (SnpData object) : The object that contains SNP values.
        pheno_vals (numpy ndarray) : A numpy array containing the phenotype
    values.
    Calls i_score_snp_continuous() with $cpus threads (default=8).

        Returns:
            i_scores, a list of I scores in the same order as the snp_data columns.
    """
    assert check_globals(["snp_data", "pheno_vals", "cont_idces"])
    print("Computing I scores of {} SNPs, for {} samples...".format(snp_data.col_count,
                                                                    (len(pheno_vals))))
    pool = multiprocessing.Pool(processes=cpus)
    i_scores = pool.map(i_score_snp_continuous, range(0, snp_data.col_count))
    pool.close()
    print("Done.")

    return np.array(i_scores)


def get_i_scores_binary(sample_size = 1.0):
    """"
        Computes and returns I scores for a binary Case/Control phenotype.

    Needs the following global objects:
        snp_data (SnpData object) : The object that contains SNP values.
        case_idces (numpy slice or list) : The indices of the case samples in
    the snp_data object.
        control_idces (numpy slice or list) : The indices of the control
    samples in the snp_data object.
    Calls i_score_snp_binary() with $cpus threads (default=8).

        Returns:
            i_scores, a list of I scores in the same order as the snp_data columns.
    """
    assert check_globals(["snp_data", "case_idces", "control_idces"])
    print("Computing I scores of {} SNPs, for {} case samples and {} controls...".format(snp_data.col_count,
                                                                                         np.shape(snp_data.val[case_idces])[0],
                                                                                         np.shape(snp_data.val[control_idces])[0]))
    pool = multiprocessing.Pool(processes=cpus)
    i_scores = pool.map(i_score_snp_binary, range(0, snp_data.col_count))
    pool.close()
    print("Done.")

    return np.array(i_scores)


def get_stat_power_binary(sample_size = 1.0):
    """"
        Args:
            snp_data (SnpData object) : The object that contains SNP values.
            pheno (dict) : A dictionnary as returned by the pysnptools loadPhen
        function, with keys 'vals', 'header' and 'iid'.
            slice_iid_cases (tuple) : A numpy slice constructed with np.s_ that
        indicates the case iids.
            slice_iid_controls (tuple) : A numpy slice constructed with np.s_ that
        indicates the controls iids.

        Returns:
            i_scores, a list of I scores in the same order as the snp_data columns.
    """
    assert check_globals(["snp_data", "case_idces", "control_idces"])
    print("Computing statistical power of {} SNPs, for {} case samples and {} controls...".format(snp_data.col_count,
                                                                                                  np.shape(snp_data.val[case_idces])[0],
                                                                                                  np.shape(snp_data.val[control_idces])[0]))
    pool = multiprocessing.Pool(processes=cpus)
    i_scores = pool.map(stat_power_snp_binary, range(0, snp_data.col_count))
    pool.close()
    print("Done.")

    return np.array(i_scores)


def get_i_score_pairwise_interaction(start, stop):
    """
    """
    global epistasy_matrix
    epistasy_matrix = scipy.sparse.dok_matrix((snp_data.col_count, stop-start+1))
    set_idxg(start)

#    for i in range(start, stop+1):
    for i,j in enumerate(np.random.choice(range(snp_data.col_count), size=stop-start)):
        set_idxg(j)
        pool = multiprocessing.Pool(processes=cpus)
        print "\t Iteration {}, SNP {} :".format(idxg, snp_data.col[idxg])
        i_scores_epistasy = pool.map(i_score_snp_pairwise_interaction, range(0, snp_data.col_count))
        epistasy_matrix[:,i-start] = np.reshape(i_scores_epistasy, (len(i_scores_epistasy), 1))
        print "Found {} interactions with I score > 1.".format(epistasy_matrix[:,i-start].count_nonzero())
        pool.close()

    return epistasy_matrix


def i_score_snp_continuous(idx):
    """
        Computes the I score of a single SNP for a continuous phenotype. Global
    variables 'snp_data' and 'pheno_vals' must have been correctly set
    (see load_data() and set_global_pheno()).

        Args:
            idx (int): The index of the SNP to retrieve.

        Returns:
            An I score (float) for the given SNP.
    """
    x = snp_data.val[cont_idces, idx]

    pheno_sum = [(np.count_nonzero(x == geno)**2 *
                  (pheno_vals[x == geno].mean() - pheno_vals_mean)**2)
                 if (np.count_nonzero(x == geno) != 0)
                 else 0
                 for geno in [0, 1, 2]]

    return np.sum(pheno_sum) / pheno_vals_total_deviation


def i_score_snp_binary(idx):
    """
        Computes the I score of a single SNP for a binary phenotype. Global
    variables 'snp_data', 'case_idces' and 'control_idces' must have been
    correctly set (see load_data() and set_global_idces()).

        Args:
            idx (int): The index of the SNP to retrieve.

        Returns:
            An I score (float) for the given SNP.
    """
    t = snp_data.val[case_idces, idx]
    t = t[~np.isnan(t)].astype(np.int)
    x = np.bincount(t, minlength=3)[0:3] #Counts for genotype [0,1,2]

    t = snp_data.val[control_idces, idx]
    t = t[~np.isnan(t)].astype(np.int)
    y = np.bincount(t, minlength=3)[0:3]

    i_score = i_score_2way(x,y)
    return i_score


def i_score_snp_pairwise_interaction(idx1, idx2=None):
    """
        Computes the I score of a single SNP for a binary phenotype. Global
    variables 'snp_data', 'case_idces' and 'control_idces' must have been
    correctly set (see load_data() and set_global_idces()).

        Args:
            idx (int): The index of the SNP to retrieve.

        Returns:
            An I score (float) for the given SNP.
    """
    if idx2 is None:
        idx2 = idxg
    t = snp_data.val[case_idces, [idx1,idx2]]
    t = t[:,0] * 3 + t[:,1] # convert to base 3(2^3) genotype number
    x = np.bincount(t[~np.isnan(t)].astype(np.int), minlength=9)

    t = snp_data.val[control_idces, [idx1,idx2]]
    t = t[:,0] * 3 + t[:,1]
    y = np.bincount(t[~np.isnan(t)].astype(np.int), minlength=9)

    i_score = i_score_2way(x,y)
    return i_score if i_score > 1 else 0;


def i_score_2way(x, y):
    """
        Compute the I score of two sample counts.

        Args:
            x (np.array) : A table that contains the counts of each genotype in
        the first cohort.
            y (np.array) : A table that contains the counts of each genotype in
        the second cohort.
    """
    nx = np.sum(x, dtype=np.float)
    ny = np.sum(y, dtype=np.float)
    i_score = nx * ny * np.sum(np.power((x / nx - y / ny),2)) / (nx + ny)
    return i_score


#def i_score_snp_binary(idx):
#    """
#        Computes the I score of a single SNP for a binary phenotype. Global
#    variables 'snp_data', 'case_idces' and 'control_idces' must have been
#    correctly set (see load_data() and set_global_idces()).
#
#        Args:
#            idx (int): The index of the SNP to retrieve.
#
#        Returns:
#            An I score (float) for the given SNP.
#    """
#    unique, count = np.unique(snp_data.val[case_idces, idx], return_counts=True)
#    table = dict(zip(unique, count))
#    x = [table.get(geno, 0) for geno in [0, 1, 2]]
#
#    unique, count = np.unique(snp_data.val[control_idces, idx], return_counts=True)
#    table = dict(zip(unique, count))
#    y = [table.get(geno, 0) for geno in [0, 1, 2]]
#
#    i_score = i_score_2way(np.array(x, dtype=np.float),
#                           np.array(y, dtype=np.float))
#    return i_score
#
#
#def i_score_2way(x, y):
#    """
#        Compute the I score of two sample counts.
#
#        Args:
#            x (np.array) : A table that contains the counts of each genotype in
#        the first cohort.
#            y (np.array) : A table that contains the counts of each genotype in
#        the second cohort.
#    """
#    nx = np.sum(x)
#    ny = np.sum(y)
#    i_score = nx * ny * np.sum(np.power((x / nx - y / ny),2)) / (nx + ny)
#    return i_score


def stat_power_snp_binary(idx):
    """
        Computes the statistical power on a single SNP for a binary phenotype.
    Global variables 'snp_data', 'case_idces' and 'control_idces' must have
    been correctly set (see load_data() and set_global_idces()).

        Args:
            idx (int): The index of the SNP to retrieve.

        Returns:
            The statistical power (float) of the given SNP.
    """
    unique, count = np.unique(snp_data.val[case_idces, idx], return_counts=True)
    table = dict(zip(unique, count))
    x = [table[geno] if geno in table.keys() else 0 for geno in [0, 1, 2]]

    unique, count = np.unique(snp_data.val[control_idces, idx], return_counts=True)
    table = dict(zip(unique, count))
    y = [table[geno] if geno in table.keys() else 0 for geno in [0, 1, 2]]

    n1 = np.sum(x, dtype=np.float)
    n2 = np.sum(y, dtype=np.float)
    p1 = (x[0] + x[1]) / n1
    p2 = (y[0] + y[1]) / n2

    z = (p1-p2)/(np.sqrt((1/n1 + 1/n2) * ((n1*p1+n2*p2)/(n1+n2)) * (1-(n1*p1+n2*p2)/(n1+n2))))
    stat_pow = 1-scipy.stats.norm.cdf(1.645 - abs(z))

    return stat_pow


def check_globals(globals_to_check):
    """
        Check if the global variables snp_data, case_idces and control_idces are
    correctly instantiated.

	Args:
	    globals_to_check (list) : A list containing the names of the global
	variables that are to be checked.

        Returns:
            A boolean, True if global variables are set, else, False.
    """
    set_globals = True
    if ("snp_data" in globals_to_check):
    	if permute_flag:
    		print "Using permuted snp_data!"
    	if (not type(snp_data) is pysnptools.snpreader.snpdata.SnpData):
        	print "Please set util.snp_data to a SnpData object containing the SNPs informations."
        	set_globals = False
    if ("case_idces" in globals_to_check) and (not type(case_idces) in [slice, list]):
        print "Please set util.case_idces to the list of indices that correspond to case samples."
        set_globals = False
    if ("control_idces" in globals_to_check) and (not type(control_idces) in [slice, list]):
        print "Please set util.control_idces to the list of indices that correspond to control samples."
        set_globals = False
    if ("pheno_vals" in globals_to_check) and (not type(pheno_vals) is np.ndarray):
        print "Please set util.control_idces to the list of indices that correspond to control samples."
        set_globals = False
    if ("cont_idces" in globals_to_check) and (not type(cont_idces) in [slice, list]):
        print "Please set util.cont_idces to the list of indices that should be included in the continuous I score computation."
        set_globals = False
    return set_globals


def draw_sample(pheno, pheno_idces, sample_size):
    """
        Draws a random sample.

        Args:
            pheno (dict) : A dictionnary as returned by the pysnptools loadPhen
        function, with keys 'vals', 'header' and 'iid'.
            pheno_idces (np.ndarray) :
            sample_size (float) : The fraction of the total data to be drawn.

        Returns:
            pheno (dict) : The sampled pheno dictionary.
            pheno_idces (np.ndarray) : The sampled indices.
    """
    idces = np.arange(len(pheno_idces), dtype=np.int)
    shuf_idces = np.random.shuffle(idces)
    sample_idces = np.sort(idces[:np.int(sample_size * len(idces))])
    pheno_idces = pheno_idces[sample_idces]
    for key in ['iid', 'vals']:
            pheno[key] = [value for i, value in enumerate(pheno[key]) if i in sample_idces]

    return pheno, pheno_idces


def reset_globals():
    """
        Resets global variables.

    """
    global snp_data
    global case_idces
    global control_idces
    global pheno_vals
    global pheno_vals_total_deviation
    global pheno_vals_mean

    snp_data = None
    case_idces = None
    control_idces = None
    pheno_vals = None
    pheno_vals_total_deviation = None
    pheno_vals_mean = None


def set_idxg(i):
    global idxg
    idxg=i
