#!/usr/bin/env python3.6

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import link_library.plot_library as plib
import statsmodels.sandbox.stats.multicomp as sm

desc = """Kai Kammer - 2018-09-02. 
Script to plot p-value distributions from xquest and xtract output and artificial distributions.
Also allows plotting of xtract raw_stat.xls files
"""


parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store",
                    help="Name of the input file")
parser.add_argument('-o', '--outname', action="store", dest="outname", default='plot.png',
                    help="Name for the output figure")
args = parser.parse_args()

def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
    """
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    copied from: https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues

def get_sorted_pvals(n=500):
    # so we get the same random values each time
    np.random.seed(80)
    # using a combination of 4/5 exponential distribution with heavy focus on low values and 1/5 uniform distribution
    # resembles experimental values quite well
    p_vals_exp = np.random.exponential(0.005, round(n-n/5))
    p_vals_uni = np.random.uniform(0.0, 1.0, round(n/5))
    p_vals = np.append(p_vals_exp, p_vals_uni)
    p_vals[p_vals < 0] = 0
    p_vals[p_vals>1] = 1
    p_vals = np.sort(p_vals)
    return p_vals

def get_bonf(p_vals):
    bonf = p_vals * len(p_vals)
    bonf[bonf > 1] = 1
    return bonf

def get_bh_fdr(p_vals, xtract=False):
    n = len(p_vals)
    # computing explicit ranks instead of np.arrange
    # done to handle equal p-values, they should get the same rank, then increase by 1
    # see: https://stats.stackexchange.com/questions/18872/how-to-deal-with-identical-p-values-with-the-benjamini-hochberg-method-for-corre
    ranks = []
    for j, val in enumerate(p_vals):
        if j == 0:
            ranks.append(1)
        else:
            if j > 0 and val == p_vals[j-1]:
                ranks.append(ranks[j-1])
            else:
                ranks.append(ranks[j-1]+1)
    ranks = np.array(ranks)
    bh_s = p_vals * float(n) / ranks
    """ 
    correcting for q-vals which are not sorted anymore; i.e. a q-val with a higher index has a lower value
    xtract does not seem to correct for this
    see also: http://www.biostathandbook.com/multiplecomparisons.html
    relevant:  The largest P value that has P<(i/m)Q is significant, and all of the P values smaller than it
    are also significant, even the ones that aren't less than their Benjamini-Hochberg critical value.       
    """
    if not xtract:
        i = n - 1
        while i > 0:
            if bh_s[i] < bh_s[i - 1]:
                bh_s[i - 1] = bh_s[i]
            i -= 1
    bh_s[bh_s > 1] = 1
    return bh_s#sm.multipletests(p_vals, method='fdr_bh')[1]



def plot_xtract_df(df):
    df = df.sort_values("pvalue")
    df = df.reset_index(drop=True)
    df["fdr_corr"] = get_bh_fdr(df["pvalue"].values, xtract=False)
    df.hist(column=["pvalue", "Bonf", "FDR", "fdr_corr"], density=True, range=[0, 1], sharey=True, grid=False)
    plib.save_fig("xtract_hist")
    df.plot(title="Stats Comparison", kind="line", y=["pvalue", "Bonf", "FDR", "fdr_corr"])
    plib.save_fig("xtract_stats")
    df.plot(title="log2", kind="kde", y=["log2ratio"])
    plib.save_fig("xtract_log2_kde")
    df = df.sort_values("log2ratio")
    df = df.reset_index()
    g = sns.jointplot(x='log2ratio', y='FDR', data=df, kind='hex')
    plib.save_fig("xtract_log2_vs_fdr_hex")
    ax = df.plot(title="log2 vs qval", kind="line", y=["log2ratio", "FDR"], secondary_y=["FDR"])
    ax_min , ax_max = ax.get_xlim()
    ax.hlines(y=[-1,1], xmin=ax_min, xmax=ax_max, linestyles='--', colors='grey')
    plib.save_fig("xtract_log2_and_fdr")
    df = df.sort_values("FDR")
    df = df.reset_index(drop=True)
    print(df["FDR"])
    ax = df.plot(title="qval kde", kind="line", y=["FDR"])

def plot_xquest_df(df):
    df = df.sort_values("FDR")
    df = df.reset_index()
    df.plot(kind="line", y=["ld-Score"], x="FDR")
    plib.save_fig("xquest_fdr_vs_ld")
    # df.plot(kind="line", y=["FDR"])
    df.plot(kind="kde", y=["ld-Score"])
    plib.save_fig("xquest_ld_kde")
    df.hist()
    plib.save_fig("xquest_hist_all")

def plot_sample_df(sample_dict):
    df = pd.DataFrame().from_dict(sample_dict)
    df.plot(title="Sign. Comparison", kind="line")
    plib.save_fig("sample_dist")
    df.hist(density=True, range=[0, 1], sharey=True)
    plib.save_fig("sample_dist_hist")

def plot_xtract_raw_stats(df):
    # FDR=FP/(FP+TP); TPR=TP/(TP+FN) (sensitivity); TNR=TN/(TN+FP) (specificity)
    #qvalue	svalue	TP	FP	TN	FN	FDR	sens	cutoff
    xval = "qvalue"
    df = df.sort_values(xval)
    # crossover is the first qvalue when the svalue (=sensitivity?) becomes one (i.e. no false negatives)
    # xtract uses this as qvalue threshold
    # here we should have no false negatives; our true positive rate is 1
    sens_1_crossover = df[df["svalue"] >= 1].iloc[0][xval]
    # df_norm = (df - df.mean()) / (df.max() - df.min())  # normalize to -1;+1
    # df_norm = (df - df.min()) / (df.max() - df.min())  # normalize to 0;1 intra-group
    df_norm = (df - min(df.min())) / max((df.max()) - min(df.min()))  # normalize to 0;1 inter-group+
    df_norm[xval] = df[xval]
    df_norm["spec"] = df_norm["TN"]/(df_norm["TN"] + df_norm["FP"])
    # note that we have normalized all columns to their min/max values of all columns (i.e. inter-group)
    # this allows to correctly calculate the FDR using more sensible numbers
    # intra-group normalization is good for visualising the contribution of each condition (TP, FP, TN, FN)
    pd_ax = df_norm.plot(x=xval, y=["TP", "FP", "TN", "FN"])

    pd_ax.vlines(x=sens_1_crossover,ymin=0,ymax=1,label="svalue cutoff")
    pd_ax.vlines(x=0.05, ymin=0, ymax=1, label="signi=5%", colors=['purple'])
    # we have to re-plot the legend at the center right position
    pd_ax.legend(loc=7)
    # the following two also work but I prefer the previous
    # pd_ax.axvline(x=sens_1_crossover)
    # pd_ax.plot((sens_1_crossover, sens_1_crossover), (0,1),label="svalue cutoff")
    # df.plot(x=xval, y=["TP", "FP", "TN", "FN"])
    plib.save_fig("xtract_raw_stats")

def main():
    uxid_string = 'uxID'
    uid_string = 'uID'
    tp_string = 'TP'
    n=100
    p_vals = get_sorted_pvals(n)
    # print(p_vals)
    bonf = get_bonf(p_vals)
    bh = get_bh_fdr(p_vals)
    # print(bh)
    bh_xt = get_bh_fdr(p_vals, xtract=True)
    if ".xls" in args.input:
        # the xls files written by xtract are buggy and can only be read this way
        df = pd.read_csv(args.input, delim_whitespace=True)
    else:
        df = pd.read_csv(args.input, engine='python')
    if tp_string in df:
        plot_xtract_raw_stats(df)
    elif uid_string in df:
        plot_xtract_df(df)
    elif uxid_string in df:
        plot_xquest_df(df)
    else:
        print("WARNING: No compatible input found")
    #plot_xtract_df("../files/merged_diff_lfq_xl_all.csv")
    #plot_xquest_df("../files/luci_merge_complete_diff_lfq_2o2.csv")
    #plot_xtract_raw_stats("../files/KK_26S_merged_final_raw_stat.xls")


    plot_dict = {"p_val": p_vals, "bonf": bonf, "fdr_xtract": bh_xt, "fdr_corr": bh}

    #plot_sample_df(plot_dict)
    plt.show()


if __name__ == "__main__":
    main()
