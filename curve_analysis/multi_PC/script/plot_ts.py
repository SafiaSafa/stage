import argparse
import numpy as np
import seaborn as sns

def get_args():
    parser = argparse.ArgumentParser(description="Plots timeseries for polygenic score results.")
    parser.add_argument("-pval", dest='res_pval', type=str, help="p-value results")
    parser.add_argument("-iscore", dest='res_iscore', type=str, help="I score results")
    parser.add_argument("-output", dest='output', type=str, help="Graphics output")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    res_pval = np.loadtxt(args.res_pval)
    res_iscore = np.loadtxt(args.res_iscore)
    #time = range(10, 250140, 500)
    time = range(10, 80000,300)
    for i in range(len(time)):
        res_pval[res_pval[:,i] == 0,i] = np.nan
        res_iscore[res_iscore[:,i] == 0,i] = np.nan
        res_pval[np.isnan(res_pval[:,i]),i] = np.nanmedian(res_pval[:,i])
        res_iscore[np.isnan(res_iscore[:,i]),i] = np.nanmedian(res_iscore[:,i])

    sns.tsplot(np.power(res_pval,1), color='b',  ci=[0,95], time=time, condition="p-value")#, err_style="unit_traces")
    sns.tsplot(np.power(res_iscore,1), color='g', ci=[0,95], time=time, condition="I score")#, err_style="unit_traces")
    sns.plt.ylabel('AUC')
    sns.plt.xlabel('N SNPs')
    sns.plt.savefig(args.output+".pdf", format="pdf")

    print "Wrote {}.".format(args.output+".pdf")

