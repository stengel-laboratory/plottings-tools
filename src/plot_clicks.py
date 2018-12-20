#!/usr/bin/env python3.6
import matplotlib
matplotlib.use('Agg')
import argparse
import pandas as pd
import seaborn as sns
import link_library.plot_library as plib
import link_library as ll

desc = """Kai Kammer - 2018-09-17. 
Script to visualize crosslink patterns
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input csv files separated by spaces")
parser.add_argument('-o', '--outname', action="store", dest="outname", default='',
                    help="Append this postfix to the plot file name.")
args = parser.parse_args()

xq_db = ll.xQuestDB()
xt_db = ll.xTractDB()

def get_clicks(df):
    prot_string = xq_db.prot1_string[:-1]
    pos_string = xq_db.pos1_string[:-1]
    vals_list = [[xq_db.prot1_string,xq_db.pos1_string],[xq_db.prot2_string,xq_db.pos2_string]]
    # new_df = pd.DataFrame(df)
    # df_count1 = new_df.groupby(
    #             [xq_db.prot1_string,xq_db.pos1_string]).size().reset_index(name='count')
    # df_count1 = df_count1.rename(
    #     index=str, columns={xq_db.pos1_string: pos_string, xq_db.prot1_string: prot_string})
    # df_count2 = new_df.groupby(
    #             [xq_db.prot2_string,xq_db.pos2_string]).size().reset_index(name='count')
    # df_count2 = df_count2.rename(
    #     index=str, columns={xq_db.pos2_string: pos_string, xq_db.prot2_string: prot_string})
    #
    # df_merge = pd.merge(df_count1, df_count2, on=[prot_string, pos_string]).set_index([prot_string, pos_string]).sum(axis=1).reset_index(name='count')
    # df_merge = df_merge.sort_values([prot_string,'count'],ascending=False).reset_index(drop=True)
    df_merge = ll.get_count_df(df, vals_list=vals_list, merge_vals_list=[prot_string, pos_string],sort_key=prot_string)
    print(df_merge)
    with sns.plotting_context("notebook",font_scale=1.5):
        fg = sns.catplot(data=df_merge, x=pos_string,y='count', col=prot_string, kind='bar', order=df_merge[pos_string])
        plib.facet_grid_vertical_label(fg)
        # ax =sns.barplot(data=df_merge, x=pos_string, y='count', order=df_merge[pos_string])
        # ax.set_xticklabels(ax.get_xticklabels(), rotation=90, size=12)
        plib.save_fig("clicks", out_dir='plots')



def main():
    df_xquest = None
    df_xtract = None
    df_dist = None
    for inp in args.input:
        if ".xls" in inp:
            # the xls files written by xtract are buggy and can only be read this way
            df = pd.read_csv(inp, delim_whitespace=True)
        else:
            df = pd.read_csv(inp, engine='python')
        if xt_db.uxid_string in df:
            df_xtract = df
        elif xq_db.uxid_string in df:
            df_xquest = df
        else:
            print("WARNING: No compatible input found for {0}".format(inp))
    if df_xquest is not None:
        get_clicks(df_xquest)
    else:
        pass


if __name__ == "__main__":
    main()