#!/usr/bin/env python

import matplotlib

matplotlib.use('agg')
import os
import configargparse
import pandas as pd
import logging
import link_library as ll
from link_library.bag_container_library import process_bag
from link_library.bag_container_library import plot_bag

# TODO: uid level is fine; however for doing the violations on uxid level they would have to be calculated before sum()
# TODO: use regular containers as a control

desc = """Kai Kammer - 2018-09-17.\n
Script to plot xTract bag container ms1 areas. All plots are by default saved to the folder 'plots'.\n
"""


# create a hybrid parser formatter allowing for default arguments and \n new lines in description and epilogue
class Formatter(configargparse.ArgumentDefaultsHelpFormatter, configargparse.RawDescriptionHelpFormatter): pass


# same functionality as hybrid formatter above; additionally recognizes 'R|' in arguments and will allow \n in
# arguments' help texts starting with 'R|'. Taken from: https://stackoverflow.com/a/22157136
class SmartFormatter(configargparse.ArgumentDefaultsHelpFormatter, configargparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return configargparse.HelpFormatter._split_lines(self, text, width)

# helper func to return a dict like this: key1: value1 \n key2: value2 \n etc...
def format_dict(my_dict):
    return "\n".join(": ".join((k, str(v))) for k, v in sorted(my_dict.items()))


valid_plots_dict = {'lh': 'light heavy', 'rep': 'replicates', 'rep_bar': 'replicates bar', 'cluster': 'cluster map',
               'std': 'standard deviation', 'link': 'ms1 area overview', 'log2r': 'log2ratio',
               'dil': 'dilution series', 'domain': 'protein domains', 'domain_sl': 'single links in protein domains',
               'dist': 'distance', 'reaction': 'reaction state (quenched vs. hydrolized monolinks)'}
parser = configargparse.ArgParser(description=desc, formatter_class=SmartFormatter)
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input csv files separated by spaces")
parser.add_argument('-V', '--version', action='version', version='%(prog)s ' + ll.__version__)
parser.add_argument('-cf', '--config_file', is_config_file=True,
                    help='Optionally specify a config file containing your settings')
parser.add_argument('-cw', '--config_write', is_write_out_config_file_arg=True,
                    help='Optionally specify a file to save your current settings to')
parser.add_argument('-o', '--outname', action="store", dest="outname", default='plots',
                    help="Name of the folder to save the plots to")
parser.add_argument('-l', '--level_ms1', action="store", dest="level", default='uID',
                    help="Level on which the ms1 intensities are summed. Either uID (peptide)"
                         "or uxID (protein).")
# choices look ugly with long list; hide them by setting metavar='' and manually print them instead
parser.add_argument('-p', '--plot_type', action="store", dest="plot", default='cluster',
                    choices=valid_plots_dict.keys(), metavar='',
                    help=f"Type of plot. Possible values are\n{format_dict(valid_plots_dict)}")
parser.add_argument('-f', '--filter', action="store", dest="filter", default=None,
                    help="Optionally specify a link type to filter for. Possible values: monolink, xlink, "
                         "intralink (loop link)")
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
parser.add_argument('-v', '--vio_list', action="store", dest='vio_list', default=['lh', 'xt'], type=str, nargs='+',
                    help="List of input possible violation filters separated by spaces. "
                         "Possible values: lh (light/heavy log2 ratio, xt (xTract type violations), none (no filtering")
parser.add_argument('-dom', '--domains', action="store", dest="domains", default=None,
                    help="Optionally specify a file containing domain ranges to color certain plots.")
parser.add_argument('-dis', '--distance', action="store", dest="distance", default=None,
                    help="Optionally specify a file containing domain ranges to color certain plots.")
parser.add_argument('-w', '--whitelist', action="store", dest="whitelist", default=None,
                    help="Optionally specify a whitelist. It is possible to filter on any valid column.")
parser.add_argument('-s', '--sortlist', action="store", dest="sortlist", default=None,
                    help="Optionally specify a file containing the order of experiments.")
parser.add_argument('-m', '--map', action="store", dest="map", default=None,
                    help="Optionally specify a file containing a mapping for values in columns.")
args = parser.parse_args()

log_file = '{0}/plot_bag_container.log'.format(args.outname)
os.makedirs(args.outname, exist_ok=True)
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s',
                    datefmt='%d/%m/%Y %H:%M:%S:')


def main():
    print('\nPlot parameters are written to {0}\n'.format(log_file))
    logging.info(f'Using link library version {ll.__version__}')
    logging.info('Start plot with the following parameters \n' + parser.format_values())
    df_list = []  # tuple (dataframe, bag_container)
    df_domains = None
    df_dist = None
    df_whitelist = None
    df_sortlist = None
    df_map = None
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
    if args.map:
        df_map = pd.read_csv(args.map, engine='python')
    bag_cont = process_bag.BagContainer(level=args.level, df_list=df_list, filter=args.filter, sel_exp=args.sel_exp,
                                        df_domains=df_domains, impute_missing=args.impute,
                                        norm_exps=args.norm_experiments,
                                        norm_reps=args.norm_replicates, df_dist=df_dist, whitelist=df_whitelist,
                                        sortlist=df_sortlist, vio_list=args.vio_list, df_map=df_map)
    plotter = plot_bag.PlotMaster(bag_cont, out_folder=args.outname)
    if args.plot == 'scatter':
        plotter.plot_scatter()
    elif args.plot == 'bar':
        plotter.plot_bar()
    elif args.plot == 'lh':
        plotter.plot_light_heavy_scatter()
    elif args.plot == 'rep':
        plotter.plot_bio_rep_scatter()
        # plotter.plot_bio_rep_ma_scatter()
    elif args.plot == 'rep_bar':
        plotter.plot_bio_rep_bar()
    elif args.plot == 'cluster':
        plotter.plot_clustermap()
    elif args.plot == 'std':
        plotter.plot_ms1_area_std()
    elif args.plot == 'link':
        plotter.plot_link_overview(int(args.experiment_percentage))
        plotter.plot_link_overview(int(args.experiment_percentage), convert_to_log2=False)
    elif args.plot == 'log2r':
        plotter.plot_log2ratio()
        plotter.plot_log2ma()
    elif args.plot == 'dil':
        plotter.plot_dilution_series()
    elif args.plot == 'domain':
        plotter.plot_domain_overview(int(args.experiment_percentage))
    elif args.plot == 'domain_sl':
        plotter.plot_domain_single_link(int(args.experiment_percentage))
        # plotter.plot_domain_single_link(int(args.experiment_percentage), ratio_mean=True)
        plotter.plot_domain_single_link(int(args.experiment_percentage), log2ratio=True)
    elif args.plot == 'monoq':
        plotter.plot_mono_vs_xlink_quant()
    elif args.plot == 'dist':
        plotter.plot_dist_vs_quant_alt()
        plotter.plot_dist_vs_quant_log2_alt()
        plotter.plot_dist_quant_corr()
    elif args.plot == 'reaction':
        plotter.plot_reaction_state()
        plotter.plot_reaction_state_log2ratio()
    else:
        print("WARNING: No compatible plot specified: {0}".format(args.input))
        exit(1)
    logging.info("Plot was successful\n")

if __name__ == "__main__":
    main()
