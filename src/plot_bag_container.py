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
parser.add_argument('-p', '--plot_type', action="store", dest="plot", default='scatter',
                    help="Type of plot. Possible values: scatter, bar,"
                         " lh (light heavy), rep (replicates), rep_bar, cluster, std (standard deviation),"
                         " link (ms1 area overview), log2r (log2ratio), dil (dilution series), domain, dist (distance)")
parser.add_argument('-f', '--filter', action="store", dest="filter", default="",
                    help="Optionally specify a link type to filter for. Possible values: monolink, xlink")
parser.add_argument('-e', '--sel_exp', action="store_true", dest="sel_exp", default=False,
                    help="Optionally provide this flag to exclude specific experiments before plotting")
parser.add_argument('-i', '--impute', action="store_true", dest="impute", default=False,
                    help="Optionally provide this flag to impute missing values for log2ratio calculations")
parser.add_argument('-nr', '--norm_replicates', action="store_true", dest="norm_replicates", default=False,
                    help="Optionally normalize replicates to their mean experimental ms1 area")
parser.add_argument('-ne', '--norm_experiments', action="store", dest="norm_experiments", default="yes",
                    help="Optionally select the experiment normalization method (or whether to normalize at all). "
                         "Possible values: yes (default norm method), xt (xTract norm method), no (do not normalize)")
parser.add_argument('-ep', '--experiment_percentage', action="store", dest="experiment_percentage", default="50",
                    help="Optionally specify the (inclusive) percentage of experiments a link has to be found in. "
                         "Only relevant for domain and link type plots"
                         "Possible values: Any value between 0 and 100")
parser.add_argument('-dom', '--domains', action="store", dest="domains", default="",
                    help="Optionally specify a file containing domain ranges to color certain plots.")
parser.add_argument('-dis', '--distance', action="store", dest="distance", default="",
                    help="Optionally specify a file containing domain ranges to color certain plots.")
parser.add_argument('-w', '--whitelist', action="store", dest="whitelist", default="",
                    help="Optionally specify a file containing allowed links (uxids), i.e. a whitelist.")
parser.add_argument('-s', '--sortlist', action="store", dest="sortlist", default="",
                    help="Optionally specify a file containing the order of experiments.")
args = parser.parse_args()


def main():
    df_list = []  # tuple (dataframe, bag_container)
    df_domains = None
    df_dist = None
    df_whitelist = None
    df_sortlist = None
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
    if args.domains:
        df_domains = pd.read_csv(args.domains, engine='python')
    if args.distance:
        df_dist = pd.read_csv(args.distance, engine='python')
    if args.whitelist:
        df_whitelist = pd.read_csv(args.whitelist, engine='python')
    if args.sortlist:
        df_sortlist = pd.read_csv(args.sortlist, engine='python')
    bag_cont = process_bag.BagContainer(level=args.level, df_list=df_list, filter=args.filter, sel_exp=args.sel_exp,
                                        df_domains=df_domains, impute_missing=args.impute, norm_exps=args.norm_experiments,
                                        norm_reps=args.norm_replicates, df_dist=df_dist, whitelist=df_whitelist,
                                        sortlist=df_sortlist)
    plotter = plot_bag.PlotMaster(bag_cont, out_folder=args.outname)
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
    elif args.plot == 'link':
        plotter.plot_link_overview(int(args.experiment_percentage))
    elif args.plot == 'log2r':
        plotter.plot_log2ratio()
    elif args.plot == 'dil':
        plotter.plot_dilution_series()
    elif args.plot == 'domain':
        plotter.plot_domain_overview(int(args.experiment_percentage))
    elif args.plot == 'monoq':
        plotter.plot_mono_vs_xlink_quant()
    elif args.plot == 'dist':
        plotter.plot_dist_vs_quant()
    else:
        print("WARNING: No compatible plot specified: {0}".format(args.input))
        exit(1)


if __name__ == "__main__":
    main()
