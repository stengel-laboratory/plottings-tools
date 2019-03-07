#!/usr/bin/env python

import pandas as pd
import argparse
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
import os.path
from matplotlib_venn import venn2,venn3
from numpy import sign
import link_library.plot_library as plib

desc = """Kai Kammer - 2018-01-29. 
Script to do a linear regression between datasets by looking at a given regression key. 
"""


parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input csv files separated by spaces")
parser.add_argument('-p', '--plots', action="store", dest="plots", default=['venn'], type=str, nargs='+',
                    help="List of figure types to be plotted separated by spaces. "
                         "Possible keywords: venn (venn diagram), regr (regression), comp (comparison)")
parser.add_argument('-k1', '--key1', action="store", dest="key1", default="uxID", type=str,
                    help="First key used for plotting. "
                         "If plotting in venn mode only this key is used.")
parser.add_argument('-k2', '--key2', action="store", dest="key2", default="ld-Score", type=str,
                    help="Second key used for plotting.")
parser.add_argument('-q', '--quant_mode', action="store_true", dest="quant", default=False,
                    help="Provide this option when working with xtract output")
parser.add_argument('-w', '--overwrite_origin', action="store_true", dest="ov_origin", default=False,
                    help="Overwrite the data origin with the given input files")
args = parser.parse_args()
origin_str = 'origin'
plot_list = ['venn', 'regr', 'comp', 'strip']



def add_origin(df_list):
    # name each dataframe, add the originating dataset
    for n, df in enumerate(df_list): #type : pd.DataFrame
        if origin_str not in df or args.ov_origin:
            df.name = os.path.basename(args.input[n])
            df[origin_str] = df.name


def split_to_origin(df):
    origin_set = set(df[origin_str].tolist())
    sel_list = []
    for origin in sorted(list(origin_set)):
        print("Origin: " + origin)
        sel = df[df[origin_str] == origin]
        sel.name = os.path.basename(origin)
        sel_list.append(sel)
    return sel_list


def plot_regression_lfq_vs_diff(df):
    df_list = split_to_origin(df)
    print("Input list length is {0}".format(len(df_list)))
    if len(df_list) == 2:
        sel1 = df_list[0]
        sel2 = df_list[1]
        key1 = df_list[0].name
        key2 = df_list[1].name
        rename_to_1 = "{0}_{1}".format(args.key2, key1)
        rename_to_2 = "{0}_{1}".format(args.key2, key2)
        sel1 = sel1.sort_values([args.key1], ascending=False)
        sel1 = sel1.rename(columns={args.key2: rename_to_1})
        sel2 = sel2.sort_values([args.key1], ascending=False)
        sel2 = sel2.rename(columns={args.key2: rename_to_2})
        sel1[rename_to_2] = sel2[rename_to_2].values
        sel2[rename_to_1] = sel1[rename_to_1].values
        df_new = pd.concat([sel1, sel2], ignore_index=True, sort=True)
        if args.quant:
            df_new['same_sign'] = df_new.apply(compare_log2ratio, args=(rename_to_1, rename_to_2), axis=1)
            sns.lmplot(x=rename_to_1, y=rename_to_2, hue='same_sign', data=df_new, fit_reg=False, legend_out=False)
            plib.save_fig("regr_log2ratio")
        else:
            df_new['regr_compare'] = df_new.apply(compare_scores, args=(rename_to_1, rename_to_2), axis=1)
            sns.lmplot(x=rename_to_1, y=rename_to_2, hue='regr_compare', data=df_new, fit_reg=False, legend_out=False)
            plib.save_fig("regr_ldscores")
        #sns.regplot(x="ld-Score_lfq", y="ld-Score_diff",data=df_new, color=df_new['score_compare'])
        sns.jointplot(x=rename_to_1, y=rename_to_2, data=df_new, kind="reg")
        # sns.lmplot(x="ld-Score_lfq", y="ld-Score_diff", col=origin_str, hue=origin_str, data=df_new,
        #             col_wrap=2, ci=None, palette="muted", size=4,
        #            scatter_kws={"s": 50, "alpha": 1})


def plot_regression_dist_vs_log2ratio(df):

    df['regulation'] = df.apply(compare_log2ratio_significance, args=("log2ratio",), axis=1)
    df['log2abs'] = abs(df["log2ratio"])*40
    df = df.sort_values(['log2abs'], ascending=True)
    lm = sns.lmplot(x=args.key1, y=args.key2, hue='regulation', data=df, fit_reg=False, legend_out=False, scatter_kws={"s":df["log2abs"]**1.3}) # sorted(df["log2abs"]) does not sort correctly :(
    plib.save_fig("log2_regulation")
    #lm.set_xticklabels(rotation='vertical', labels=df['uID'])
    sel = df[df['regulation'] != "not_signi"]
    #sns.jointplot(x=args.key1, y=args.key2, data=sel, kind="reg")
    #sns.jointplot(x=args.key1, y=args.key2, data=df, kind="reg")


def remove_dups(df, rem_key="ld-Score", sort_key="uxID"):
    df = df.sort_values(rem_key, ascending=False).drop_duplicates(sort_key)
    return df


def plot_venn(df):
    df_list = split_to_origin(df)
    print("Venn: Len of origin is {0}".format(len(df_list)))
    if len(df_list) == 2:
        set1 = set(df_list[0][args.key1].tolist())
        set2 = set(df_list[1][args.key1].tolist())
        venn2([set1, set2], set_labels=(df_list[0].name, df_list[1].name), alpha=0.65)
        plib.save_fig("venn_{0}".format(args.key1))
    elif len(df_list) == 3:
        set1 = set(df_list[0][args.key1].tolist())
        set2 = set(df_list[1][args.key1].tolist())
        set3 = set(df_list[2][args.key1].tolist())
        venn3([set1, set2, set3], alpha=0.65, set_labels=(df_list[0].name,
                                              df_list[1].name,
                                              df_list[2].name))
        plib.save_fig("venn_{0}".format(args.key1))
    else:
        print("ERROR: Too many entries for a venn diagramm: {0}".format(len(df_list)))
        exit(1)

def compare_scores(row, key1, key2):
    if row[key1] == row[key2]:
        val = "equal"
    elif row[key1] > row[key2]:
        val = key1
    else:
        val = key2
    return val

def compare_log2ratio(row, key1, key2):
    if sign(row[key1]) == sign(row[key2]):
        val = "yes"
    else:
        val = "no"
    return val


def compare_log2ratio_significance(row, key):
    if -1 < row[key] < 1:
        val = "not_signi"
    elif row[key] > 1.0:
        val = "up"
    else:
        val = "down"

    return val


def plot_dist(df):
    sns.jointplot(data=df, x="FDR", y="ld-Score")

def plot_strip(df):
    df = df.sort_values([args.key2], ascending=False)
    ax = sns.stripplot(data=df, x=args.key1, y=args.key2, hue=origin_str, size=10, linewidth=1)
    ax.set_xticklabels(rotation='vertical', labels=df[args.key1])
    ax.grid(which='major', axis='both', color='darkgrey')

def plot_facet(df):
    sns.FacetGrid(df, row=origin_str, col=args.key1, margin_titles=True)

def main():
    df_list = [pd.read_csv(f, sep=None, engine='python') for f in args.input]
    add_origin(df_list)
    # df is the final output dataframe concatenated from the modified input files
    df = pd.concat(df_list, sort=True)  # type: pd.DataFrame

    #plot_cat(df)
    if 'venn' in args.plots:
        plot_venn(df)
    if 'regr' in args.plots:
        plot_regression_lfq_vs_diff(df)
    if 'comp' in args.plots and args.quant:
        plot_regression_dist_vs_log2ratio(df)


if __name__ == "__main__":
    main()
