#!/usr/bin/env python3.6

import matplotlib
matplotlib.use('agg')
import os
import argparse
import pandas as pd
from link_library.bag_container_library import plot_bag
from link_library.bag_container_library import process_bag

# TODO: uid level is fine; however for doing the violations on uxid level they would have to be calculated before sum()
# TODO: use regular containers as a control

desc = """Kai Kammer - 2018-09-17. 
Script to plot xTract bag container ms1 areas. All plots are saved to the folder 'plots'
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input csv files separated by spaces")
parser.add_argument('-o', '--outname', action="store", dest="outname", default='plots',
                    help="Name of the plot folder")
parser.add_argument('-l', '--level_ms1', action="store", dest="level", default='uID',
                    help="Level on which the ms1 intensities are summed. Either uID (peptide)"
                         "or uxID (protein).")
parser.add_argument('-f', '--filter', action="store", dest="filter", default="",
                    help="Optionally specify a link type to filter for. Possible values: monolink, xlink")
parser.add_argument('-p', '--plot_type', action="store", dest="plot", default='scatter',
                    help="Type of plot. Possible values: scatter, bar,"
                         " lh (light heavy), rep (replicates), rep_bar, cluster, std (standard deviation)")
args = parser.parse_args()


def main():
    df_list = []  # tuple (dataframe, bag_container)
    uid_string = "b_peptide_uID"
    for inp in args.input:
        if ".xls" in inp:
            # the xls files written by xtract are buggy and can only be read this way
            df = pd.read_csv(inp, engine='python', delimiter='\t', na_values=['-'])
        else:
            df = pd.read_csv(inp, engine='python')
        df.name = os.path.basename(inp)
        df._metadata += ['name']
        if uid_string in df:
            df_list.append(df)
        else:
            print("WARNING: No compatible input found for {0}".format(args.input))
            exit(1)
    bag_cont = process_bag.BagContainer(level=args.level, df_list=df_list, filter=args.filter)
    plotter = plot_bag.PlotMaster(bag_cont)
    if args.plot == 'scatter':
        plotter.plot_scatter()
    elif args.plot == 'bar':
        plotter.plot_bar()
    elif args.plot == 'lh':
        plotter.plot_light_heavy_scatter()
    elif args.plot == 'rep':
        plotter.plot_bio_rep_scatter()
    elif args.plot == 'rep_bar':
        plotter.plot_bio_rep_bar()
    elif args.plot == 'cluster':
        plotter.plot_clustermap()
    elif args.plot == 'std':
        plotter.plot_ms1_area_std()
    else:
        print("WARNING: No compatible plot specified: {0}".format(args.input))
        exit(1)


if __name__ == "__main__":
    main()
