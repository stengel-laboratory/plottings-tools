#!/usr/bin/env python3.6

import argparse
import pandas as pd
import seaborn as sns
import link_library.plot_library as plib
import link_library as ll

desc = """Kai Kammer - 2019-01-17. 
Script to boxplot log2ratio distributions from xtract output.
Useful for dilution series.
"""


parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store",
                    help="Name of the input file")
parser.add_argument('-o', '--outname', action="store", dest="outname", default='xtract_dilution_series',
                    help="Name for the output figure")
args = parser.parse_args()

xt_db = ll.xTractDB()

def plot_dil(df):
    df = df.sort_values(xt_db.exp_string)
    ax = sns.boxplot(data=df, x=xt_db.exp_string, y=xt_db.log2_string)
    plib.save_fig(args.outname, 'plots')
    print(df.groupby(xt_db.exp_string).mean())


def main():
    if ".xls" in args.input:
        # the xls files written by xtract are buggy and can only be read this way
        df = pd.read_csv(args.input, delim_whitespace=True)
    else:
        df = pd.read_csv(args.input, engine='python')
    if xt_db.uxid_string in df:
        plot_dil(df)
    else:
        print("WARNING: No compatible input found")

if __name__ == "__main__":
    main()